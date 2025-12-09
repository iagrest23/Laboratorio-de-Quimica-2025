import h5py
from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDetermineBonds
import io
import numpy as np
import os
from collections import defaultdict
from dadapy import Data
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


pt = Chem.GetPeriodicTable()


def mol_to_xyz_block(mol: Chem.Mol, conf_id: int = -1) -> str:
    """
    Convert an RDKit Mol to an XYZ format string.
    Assumes the molecule has 3D coordinates.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformers (no 3D coordinates)")

    conf = mol.GetConformer(conf_id)
    n_atoms = mol.GetNumAtoms()
    lines = [str(n_atoms), ""]

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol()
        lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

    return "\n".join(lines)

def dist(r1, r2):
    d = r1 - r2
    return np.sqrt(np.sum(d**2))

def has_close_contact(mol, threshold=0.1, conf_id=-1):
    conf = mol.GetConformer(conf_id)
    xyz_array = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        current_coordinates = np.array([pos.x, pos.y, pos.z])
        for coordinates in xyz_array:
            distance = dist(coordinates, current_coordinates)
            if distance <= threshold:
                return True
        xyz_array.append(current_coordinates)
    return False


def xyz_to_rdkitmol_openbabel(species: list[int], coordinates) -> Chem.Mol:
    if len(species) != len(coordinates):
        raise ValueError("Mismatched species and coordinates")
    
    from rdkit.Chem import GetPeriodicTable
    xyz_lines = [str(len(species)), "Converted from xyz"]
    pt = Chem.GetPeriodicTable()
    for Z, (x, y, z) in zip(species, coordinates):
        symbol = Z
        xyz_lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    xyz_string = "\n".join(xyz_lines)
   # print(xyz_string)
    # Convert to OBMol
    obmol = pybel.readstring("xyz", xyz_string)
    # Add hydrogens and perceive bond orders BEFORE converting
    obmol.addh()  # Add hydrogens
    obmol.make3D()  # Ensure 3D coordinates
    
    # Convert OBMol to RDKit Mol
    obConversion = ob.OBConversion()
    obConversion.SetOutFormat("mol")
    mol_block = obConversion.WriteString(obmol.OBMol)

    rdkit_mol = Chem.MolFromMolBlock(mol_block, sanitize=False)  # Don't sanitize initially


    if rdkit_mol is None:
        raise ValueError("RDKit failed to parse the mol block")
    
    # Now sanitize the molecule
    try:
        Chem.SanitizeMol(rdkit_mol)
        
    except:
        # If sanitization fails, try without kekulization
        Chem.SanitizeMol(rdkit_mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        

    return rdkit_mol


def species_coordinates_to_smiles(r,z,charge=0,useHueckel=True):
      
    try:
        mol=xyz_to_rdkitmol_openbabel(z,r)
    except Exception as e:
        return "SMILESN'T"

        
    if has_close_contact(mol):
        return "SMILESN'T"
    try:
         rdDetermineBonds.DetermineBonds(mol,charge=0,useHueckel=True)
         

         smiles=Chem.MolToSmiles(mol)
         return smiles
    except Exception as e:
         return "SMILESN'T"

def species_to_formula(species_key):
    counts = Counter(species_key)
    # Convert atomic numbers → symbols, order C, H, then others alphabetically
    def sort_key(sym):
        if sym == "C": return (0, sym)
        if sym == "H": return (1, sym)
        return (2, sym)

    parts = []
    for Z in sorted(counts, key=lambda z: sort_key(pt.GetElementSymbol(z))):
        sym = pt.GetElementSymbol(Z)
        n = counts[Z]
        parts.append(sym if n == 1 else f"{sym}{n}")
    return "".join(parts)


def distance_matrix(r):
    """Calculate distance matrix for a set of coordinates"""
    diff = r[:, None, :] - r[None, :, :]  # shape (n_atoms, n_atoms, 3)
    D = np.sqrt(np.sum(diff**2, axis=-1))  # shape (n_atoms, n_atoms)
    return D


def vectorize_distance_matrix(D):
    """Convert distance matrix to vector using upper triangle"""
    return D[np.triu_indices_from(D, k=1)]

def reassign_clusters(X, clusters_id):
    """
    Reassign Conformers with cluster_id = -1 to their nearest cluster
    based on Euclidian distance to cluster centroids

    """
    clusters_id = clusters_id.copy()

    ## Get Unique cluster IDs (excluding -1)

    unique_clusters_ids = np.unique(clusters_id)
    unique_clusters_ids = unique_clusters_ids[unique_clusters_ids != -1]

    if len(unique_clusters_ids) == 0:
        print("Warning: No valid clusters found")
        return clusters_id

    ## Calculate average distance vector (centroid) for each valid cluster

    cluster_centroids = []
    for cluster_id in unique_clusters_ids:
        # Get mask for current cluster
        mask = clusters_id == cluster_id
        # Get all distance vectors for this cluster
        cluster_distances = X[mask]
        # Calculate centroid
        centroid = np.mean(cluster_distances, axis=0)
        cluster_centroids.append(centroid)

    cluster_centroids = np.array(cluster_centroids)

    # Find outliers (cluster_id = -1) and reassing them

    outlier_indices = np.where(clusters_id == -1)[0]

    if len(outlier_indices) == 0:
        print("No outliers to reassing")
        return clusters_id

    print(f"Reassigning {len(outlier_indices)} outliers to nearest cluster")

    for idx in outlier_indices:
        # Get distance vector for this outlier
        outlier_vector = X[idx]

        # Calculate Euclidean distances to all cluster centroids

        distances_to_centroids = np.linalg.norm(
                cluster_centroids - outlier_vector[np.newaxis, :],
                axis=1
        )

        # Find nearest cluster

        nearest_cluster_idx = np.argmin(distances_to_centroids)
        nearest_cluster_id = unique_clusters_ids[nearest_cluster_idx]

        # Reassign outlier to nearest cluster
        clusters_id[idx] = nearest_cluster_id

    return clusters_id

if __name__ == "__main__":
    # code here runs only when the file is executed directly
    main()

