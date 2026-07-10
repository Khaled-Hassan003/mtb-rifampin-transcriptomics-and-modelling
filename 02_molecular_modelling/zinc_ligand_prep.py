#!/usr/bin/env python3
# zinc_ligand_prep.py
# Automated ZINC Ligand Fetching and RDKit Preparation

import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

# Read Target ZINC IDs
with open("zinc_ids.txt", "r") as f:
    zinc_ids = [line.strip() for line in f if line.strip()]

output_folder = "prepared_ligands"
os.makedirs(output_folder, exist_ok=True)

# Fetch from ZINC Database and Standardize 3D Coordinates
for zid in zinc_ids:
    url = f"https://zinc15.docking.org/substances/{zid}.sdf"
    response = requests.get(url)
    
    if response.status_code == 200:
        sdf_path = f"{output_folder}/{zid}.sdf"
        with open(sdf_path, "wb") as f:
            f.write(response.content)
            
        # Optimize structure using RDKit for simulation embedding
        suppl = Chem.SDMolSupplier(sdf_path)
        for mol in suppl:
            if mol:
                mol_h = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol_h)
                
                writer = Chem.SDWriter(f"{output_folder}/{zid}_optimized.sdf")
                writer.write(mol_h)
                writer.close()
