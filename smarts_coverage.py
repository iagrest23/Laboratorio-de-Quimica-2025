#!/usr/bin/env python
"""
smarts_coverage.py

Descripción
-----------
Script para analizar la cobertura de patrones SMARTS en uno o varios
archivos .csv que contienen listas de SMILES. El procedimiento consiste en:

1. Cargar un archivo CSV con definiciones SMARTS 
   (columnas obligatorias: "Description", "SMARTS").
2. Leer todos los archivos .csv dentro de una carpeta especificada,
   cada uno conteniendo una columna llamada 'smiles'.
3. Para cada patrón SMARTS, contar cuántos SMILES de cada archivo presentan
   coincidencia de subestructura (RDKit: Mol.HasSubstructMatch).
4. Guardar una tabla resumen en formato CSV con la cantidad de matches por
   archivo y por SMARTS.
5. Generar una copia de respaldo automática del archivo SMARTS original
   si aún no existe.

Parámetros (CLI)
----------------
--folder : str
    Carpeta que contiene los .csv con la columna 'smiles'.
--smarts : str
    Archivo CSV con definiciones SMARTS. Default: smarts.csv.
--out : str
    Nombre del CSV de salida con la cobertura de SMARTS. Default: smarts_coverage.csv.

Uso
---
Ejemplo de ejecución:

    python smarts_coverage.py \
        --folder datasets_SMILES \
        --smarts functional_groups.csv \
        --out smarts_coverage_results.csv

Requisitos
----------
- RDKit
- pandas
- glob, argparse
- Archivos CSV con columnas válidas

Notas
-----
Este script forma parte del pipeline de análisis químico-informático del 
proyecto "Laboratorio de Química 2025". Permite cuantificar la prevalencia 
de grupos funcionales definidos por SMARTS dentro de grandes colecciones 
de SMILES.
"""


import os
import pandas as pd
from rdkit import Chem
from glob import glob
import shutil

def load_smarts(smarts_file):
    """Load SMARTS definitions from a CSV file."""
    if not os.path.exists(smarts_file):
        raise FileNotFoundError(f"SMARTS file '{smarts_file}' not found.")
    df = pd.read_csv(smarts_file)
    if not {"Description", "SMARTS"}.issubset(df.columns):
        raise ValueError("SMARTS CSV must have 'Description' and 'SMARTS' columns.")
    return df

def count_smarts_matches(smiles_list, smarts):
    """Count number of SMILES that match a SMARTS pattern."""
    pattern = Chem.MolFromSmarts(smarts)
    count = 0
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol and mol.HasSubstructMatch(pattern):
            count += 1
    return count

def analyze_smarts_coverage(folder_path, smarts_file, output_csv="smarts_coverage.csv"):
    """Analyze SMARTS coverage across all CSV files in a folder."""
    # Backup original SMARTS list if not yet done
    backup = os.path.splitext(smarts_file)[0] + "_original.csv"
    if not os.path.exists(backup):
        shutil.copy(smarts_file, backup)
        print(f"💾 Backup of SMARTS file created as: {backup}")

    # Load SMARTS definitions
    smarts_df = load_smarts(smarts_file)
    csv_files = glob(os.path.join(folder_path, "*.csv"))
    if not csv_files:
        raise FileNotFoundError("No .csv files found in the provided folder.")

    # Process each SMARTS pattern
    results = []
    for _, row in smarts_df.iterrows():
        desc, smarts = row["Description"], row["SMARTS"]
        result_row = {"Description": desc, "SMARTS": smarts}
        for file in csv_files:
            df = pd.read_csv(file)
            if "smiles" not in df.columns:
                raise ValueError(f"'smiles' column not found in {file}")
            smiles_list = df["smiles"].dropna().tolist()
            result_row[os.path.basename(file)] = count_smarts_matches(smiles_list, smarts)
        results.append(result_row)

    # Export results
    df_out = pd.DataFrame(results)
    df_out.to_csv(output_csv, index=False)
    print(f"\n✅ SMARTS coverage analysis saved to: {output_csv}\n")
    print(df_out)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyze SMARTS coverage across multiple CSV datasets.")
    parser.add_argument("--folder", required=True, help="Path to folder containing .csv files with a 'smiles' column")
    parser.add_argument("--smarts", default="smarts.csv", help="Path to SMARTS definition file (CSV with Description,SMARTS)")
    parser.add_argument("--out", default="smarts_coverage.csv", help="Output CSV filename")
    args = parser.parse_args()

    analyze_smarts_coverage(args.folder, args.smarts, args.out)

