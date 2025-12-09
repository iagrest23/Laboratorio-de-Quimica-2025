# Laboratorio de Química 2025  
### Tesina de Grado – Bases de Datos Masivas de Química Cuántica  
**Caracterización quimioinformática y entrenamiento de redes neuronales**

Este repositorio reúne el código desarrollado para la tesina de grado de la Licenciatura en Ciencias Químicas (UBA), cuyo objetivo es:

- Analizar bases de datos cuánticas tipo ANI a gran escala.  
- Realizar caracterización quimioinformática mediante SMILES/SMARTS.  
- Detectar clusters estructurales basados en geometría molecular.  
- Entrenar modelos tipo ANI con TorchANI.  

Incluye herramientas de preprocesamiento, clustering, sanitización de SMILES, análisis de cobertura de grupos funcionales y entrenamiento de redes neuronales.

---

# Contenidos principales del repositorio

## 1. `ani_conformer_smiles_clustering.py`
Script principal del pipeline de procesamiento ANI.

Funciones principales:

- Carga un dataset ANI (HDF5, `torchani.datasets.ANIDataset`).
- Reagrupa conformeros por fórmula molecular.
- Calcula matrices de distancia y las vectoriza.
- Realiza clustering por fórmula usando ADP (DADApy).
- Genera SMILES desde coordenadas (`species_coordinates_to_smiles`).
- Corrige SMILES fallidos dentro del cluster.
- Escribe:
  - Un nuevo archivo ANI reordenado.
  - Un CSV con metadatos de cada conformero (índice, fórmula, cluster, SMILES original y sanitizado).

---

## 2. `smilestools.py`
Módulo de utilidades químico-informáticas:

- Conversión `species + coordinates → RDKit Mol` vía OpenBabel.
- Chequeo de contactos cercanos.
- Generación robusta de SMILES en 3D.
- Cálculo de fórmulas moleculares.
- Cálculo/vectorización de matrices de distancia.
- Reasignación de outliers de clustering.

---

## 3. `train_ani_model.py`
Script para el entrenamiento de un modelo ANI desde cero usando TorchANI.

Incluye:

- Carga del dataset ANI-2x.
- Construcción de un modelo `simple_ani`.
- Optimización con AdamW y scheduler ReduceLROnPlateau.
- Cálculo del RMSE de validación.
- Guardado de checkpoints (`latest_training_state.pt`, `best_model_state.pt`).
- Registro de métricas en TensorBoard.

---

## 4. `smarts_coverage.py`
Herramienta para analizar la cobertura de patrones SMARTS en múltiples listas de SMILES.

Permite:

- Leer un archivo `SMARTS.csv` con columnas `Description` y `SMARTS`.
- Recorrer todos los `.csv` en una carpeta con la columna `smiles`.
- Contar coincidencias SMARTS–SMILES mediante RDKit.
- Exportar un CSV resumen.
- Generar automáticamente un backup del archivo SMARTS original.

---

# Estructura del repositorio

├── ani_conformer_smiles_clustering.py
├── smilestools.py
├── train_ani_model.py
├── smarts_coverage.py
└── README.md


---

# Requisitos

Se recomienda Python 3.9+ y un entorno virtual (conda o mamba).

Dependencias principales:

- [TorchANI](https://github.com/aiqm/torchani) (para ANIDataset)
- [NumPy](https://numpy.org/)
- [RDKit](https://www.rdkit.org/)
- [OpenBabel](https://openbabel.org/) con bindings de Python (openbabel, pybel)
- [DADApy](https://dadapy.readthedocs.io/) (ADP clustering)
- [pandas](https://pandas.pydata.org/)
- [h5py](https://www.h5py.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)

Ejemplo de instalación:

```bash
conda create -n ani-smiles python=3.10
conda activate ani-smiles

conda install -c conda-forge rdkit openbabel
pip install torchani dadapy h5py numpy matplotlib seaborn















