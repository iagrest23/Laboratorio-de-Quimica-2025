#!/usr/bin/env python
"""
ani_conformer_smiles_clustering.py

Descripción
-----------
Script para:
1. Cargar un dataset ANI (formato HDF5, usando TorchANI ANIDataset).
2. Reagrupar conformeros por fórmula molecular.
3. Calcular matrices de distancia y vectorizarlas para cada conformero.
4. Hacer clustering de conformeros (por fórmula) usando DADApy (ADP).
5. Generar SMILES a partir de (species, coordinates) usando smilestools.
6. "Arreglar" SMILES fallidos dentro de cada cluster (reemplazando "SMILESN'T"
   por un SMILES válido del mismo cluster, si existe).
7. Escribir un nuevo HDF5 con los conformeros (agrupados por número de átomos).
8. Guardar un CSV con metadatos de cada conformero:
   índice global, fórmula molecular, ID de cluster, SMILES original y SMILES tras
   la sanitización por cluster.

Parámetros (CLI)
----------------
--patho : str
    Ruta al archivo HDF5 de entrada con el dataset ANI (ANIDataset).
--group : str
    Identificador del grupo (ej. "004"). Se usa solo para nombrar el CSV de salida
    (conformer_metadata_group_<group>.csv).
--oname : str
    Nombre del archivo HDF5 de salida (nuevo dataset ANI reescrito).

Uso
---
Ejemplo de ejecución:

    python ani_conformer_smiles_clustering.py \
        --patho ani2x_group_004.h5 \
        --group 004 \
        --oname ani2x_SMILES_004.h5

Requisitos
----------
- TorchANI (ANIDataset)
- dadapy
- smilestools (módulo que provee species_coordinates_to_smiles, etc.)
- NumPy, csv, argparse
"""



import numpy as np
from collections import defaultdict
from smilestools import species_coordinates_to_smiles
from smilestools import species_to_formula, distance_matrix, vectorize_distance_matrix, reassign_clusters
from torchani.datasets import ANIDataset
import argparse
import warnings
from dadapy import Data
import seaborn as sns
from collections import Counter
import csv

warnings.filterwarnings("ignore")

class ConformerData:
    """Class to store conformer information including coordinates and energy"""
    def __init__(self, coordinates, energy, species, forces, original_index):
        self.coordinates = coordinates
        self.energy = energy
        self.species = species
        self.forces = forces
        self.original_index = original_index


def parse_input():
    parser = argparse.ArgumentParser(description="Generate SMILES for a given group.")
    
    parser.add_argument("--patho", required=True, type=str,
                        help='Path to dataset')
    parser.add_argument("--group", required=True, type=str,
                        help="Group identifier (e.g., '004')")
    parser.add_argument("--oname", required=True, type=str,
                        help="Output file name (e.g., 'ani2x_SMILES_004.h5')")

    return parser.parse_args()

def refactor_ds_by_formula(ds):
    """
    Returns a dict where each key is a molecular formula string, and its
    associated items are coordinates, species, energies, forces and 
    original index. Uses TorchANI's regroup_by_formula functionality.
    """
    # Regroup the dataset by molecular formula
    regrouped_ds = ds.regroup_by_formula()

    formula_conformers = defaultdict(list)
    
    # Iterate through each formula group in the regrouped dataset
    
    for formula in regrouped_ds.keys():
        
        group_data = regrouped_ds[formula]
        
        coordinates = group_data['coordinates']
        species = group_data['species']
        energies = group_data['energies']
        forces = group_data['forces']
        
        if hasattr(coordinates, 'numpy'):
            coordinates = coordinates.numpy()
        if hasattr(species, 'numpy'):
            species = species.numpy()
        if hasattr(energies, 'numpy'):
            energies = energies.numpy()
        if hasattr(forces, 'numpy'):
            forces = forces.numpy()


        for i in range(len(coordinates)):
            r = coordinates[i]
            z = species[i]
            e = energies[i]
            f = forces[i]
            
            conformer_data = ConformerData(r, e, z, f, i)
            formula_conformers[formula].append(conformer_data)
    
    return formula_conformers

def get_dist_vectors_by_formula(formula_conformers):
    """
    Reads formula_conformers dict  and outputs a dict with 
    the flattened distance matrix of each conformer as items
    of each molecular formula
    """
    dist_vectors = {}
    for formula, conformers in formula_conformers.items():
        vectors = []
        for conf_data in conformers:
            D = distance_matrix(conf_data.coordinates)
            v = vectorize_distance_matrix(D)
            vectors.append(v)
        dist_vectors[formula] = np.vstack(vectors)
    return dist_vectors

def get_clusters_id_for_formula(dist_vectors,formula):
    """
    Get cluster IDs for conformers of a specific molecular formula
    """    
    X = dist_vectors[formula]

    nconfs, _ = np.shape(X)

    if nconfs <=2:
        clusters_id = np.array([0]*nconfs)
    else:
        data = Data(X, verbose=False)
        data.compute_distances()
        data.compute_id_2NN()
        clusters_id = data.compute_clustering_ADP(Z=4, halo=True)
        clusters_id = np.array(clusters_id, dtype=int)

        unique_initial = np.unique(clusters_id)
        n_outliers = np.sum(clusters_id == -1)
        n_clusters = len(unique_initial[unique_initial != -1])

        if n_outliers > 0:
            clusters_id = reassign_clusters(X, clusters_id)
            unique_final = np.unique(clusters_id)
            final_n_clusters = len(unique_final)

    return clusters_id

def get_smiles_for_group_of_conformers(zs,rs):
    """
    Get SMILES for group of conformers (hardcoded to neutral 
    and Huckel on)
    """
    ss = []
    for index, (z, r) in enumerate(zip(zs, rs)):
        try:
            smi = species_coordinates_to_smiles(r, z, 0, True)
        except Exception as e:
            smi = "SMILESN'T"

        ss.append(smi)
    return np.asarray(ss)

def sanitize_ss_by_cs(ss,cs):
    
    """
    Sanitizes SMILES strings by cluster ID.
    
    For each molecule with SMILES == "SMILESN'T", finds other molecules 
    in the same cluster that have valid SMILES and assigns the first 
    valid SMILES found to replace "SMILESN'T".
    
    Args:
        ss (list): List of SMILES strings
        cs (list): List of cluster IDs corresponding to each SMILES
        
    Returns:
        list: Sanitized list of SMILES strings
    """
    
    # Create a copy yo avoid modifying the original list

    sanitized_ss = ss.copy()

    # Create a dictionary to group indices by cluster_id

    cluster_groups = {}

    for i, cluster_id in enumerate(cs):
        if cluster_id not in cluster_groups:
            cluster_groups[cluster_id] = []
        cluster_groups[cluster_id].append(i)

    # Go through each molecule

    for i, smiles in enumerate(ss):
        if smiles == "SMILESN'T":
            cluster_id = cs[i]

            # Find other molecules in the same cluster

            cluster_indices = cluster_groups[cluster_id]

            # Look for a valid SMILES in the same cluster

            valid_smiles = None
            for j in cluster_indices:
                if ss[j] != "SMILESN'T":
                    valid_smiles = ss[j]
                    break
            
            # If we found a valid SMILES, assign it
            if valid_smiles is not None:
                sanitized_ss[i] = valid_smiles

    return sanitized_ss


def extract_conformer_data(molecule):
    """
    Extract conformer data arrays from ConformerData objects
    """
    zs = np.array([conf.species for conf in molecule])
    rs = np.array([conf.coordinates for conf in molecule])
    fs = np.array([conf.forces for conf in molecule])
    es = np.array([conf.energy for conf in molecule])
    return zs, rs, fs, es

def main():
    args = parse_input()
    patho = args.patho
    group = args.group
    oname = args.oname
    ds = ANIDataset(patho)
    formula_conformers = refactor_ds_by_formula(ds)
    dist_vectors = get_dist_vectors_by_formula(formula_conformers)
    new_ds = ANIDataset(oname, grouping="by_num_atoms")
    
    # Open file for writing
    with open(f"conformer_metadata_group_{group}.csv", "w", newline='') as csvfile:
        # Write header
        writer = csv.writer(csvfile)

        writer.writerow(["Index", "Molecular_Formula", "Cluster_ID", "SMILES","Post_Cluster_SMILES"])
        
        c = 0
        for i, formula in enumerate(formula_conformers):
            
            if i == 501:  # Skip missing or unwanted molecular formula
                print(f"Skipping molecular formula index {i}")
                continue

            molecule = formula_conformers[formula]
            zs, rs, fs, es = extract_conformer_data(molecule)
            cs = get_clusters_id_for_formula(dist_vectors, formula)
            fms = np.asarray([formula]*len(zs))
            ss = get_smiles_for_group_of_conformers(zs, rs)
            sss = sanitize_ss_by_cs(ss, cs)
            
            indices = np.asarray([c + j for j in range(len(zs))])

            numpy_conformers = {
                "index": indices,
                "species": zs,
                "coordinates": rs,
                "forces": fs,
                "energies": es
            }

            new_ds.auto_append_conformers(numpy_conformers)

            ## Write metadata to CSV for each conformer
            for j in range(len(zs)):
                
                writer.writerow([c + j, fms[j], cs[j], ss[j], sss[j]])
            
            c = c + len(zs)

if __name__ == "__main__":
    main()

