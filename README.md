# Laboratorio-de-Quimica-2025
Material asociado al proyecto desarrollado en la materia Laboratorio de Química. Tesina de grado de Licenciatura en Ciencias Químicas: Bases de Datos Masivas de Química Cuantica - Caracterizacion Quimioinformatica y Entrenamiento de Redes Neuronales.


# ANI Conformer SMILES & Clustering

Este repositorio contiene el código desarrollado en el marco de una tesina de grado para:

- Reorganizar un dataset ANI (formato HDF5 de TorchANI) por fórmula molecular.
- Calcular descriptores geométricos basados en matrices de distancia.
- Realizar clustering de conformeros para cada fórmula molecular usando DADApy.
- Generar SMILES a partir de coordenadas 3D (OpenBabel + RDKit).
- Sanitizar SMILES fallidos dentro de cada cluster.
- Escribir un nuevo dataset ANI y un archivo CSV con metadatos de cada conformero.

Está pensado para usarse sobre datasets del tipo ANI (por ejemplo, ANI-2x) descargados de Zenodo y manejados a través de `torchani.datasets.ANIDataset`.

---

## Estructura del repositorio

- `ani_conformer_smiles_clustering.py`  
  Script principal.  
  - Lee un archivo HDF5 con un dataset ANI.
  - Reagrupa conformeros por fórmula molecular.
  - Calcula vectores de distancias interatómicas para cada conformero.
  - Aplica clustering (ADP, vía DADApy) dentro de cada fórmula.
  - Genera SMILES para cada conformero (`species_coordinates_to_smiles`).
  - Reemplaza SMILES fallidos (`"SMILESN'T"`) por un SMILES válido del mismo cluster.
  - Escribe:
    - Un nuevo archivo HDF5 con los conformeros (agrupados por número de átomos).
    - Un CSV con metadatos de cada conformero (índice, fórmula, cluster, SMILES).

- `smilestools.py`  
  Módulo de utilidades que provee:
  - Conversión `species + coordinates → RDKit Mol` usando OpenBabel (`xyz_to_rdkitmol_openbabel`).
  - Detección de contactos demasiado cercanos (`has_close_contact`).
  - Generación robusta de SMILES desde coordenadas 3D (`species_coordinates_to_smiles`).
  - Cálculo de fórmula molecular a partir de números atómicos (`species_to_formula`).
  - Cálculo y vectorización de matrices de distancia (`distance_matrix`, `vectorize_distance_matrix`).
  - Reasignación de outliers de clustering al cluster más cercano (`reassign_clusters`).

---

## Requisitos

Se recomienda usar Python 3.9+ y un entorno virtual (conda, mamba o venv).

Dependencias principales:

- [TorchANI](https://github.com/aiqm/torchani) (para `ANIDataset`)
- [NumPy](https://numpy.org/)
- [RDKit](https://www.rdkit.org/)
- [OpenBabel](https://openbabel.org/) con bindings de Python (`openbabel`, `pybel`)
- [DADApy](https://dadapy.readthedocs.io/) (ADP clustering)
- `h5py`
- `matplotlib` y `seaborn` (algunas utilidades del módulo de SMILES)

Ejemplo de instalación (orientativo, según tu entorno):

```bash
conda create -n ani-smiles python=3.10
conda activate ani-smiles

# RDKit (canal conda-forge)
conda install -c conda-forge rdkit

# OpenBabel + pybel
conda install -c conda-forge openbabel

# TorchANI
pip install torchani

# DADApy y otros
pip install dadapy h5py numpy matplotlib seaborn

