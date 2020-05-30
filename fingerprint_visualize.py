#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 14:08:49 2020

@author: shibau
"""

from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw
import pubchempy as pcp
from IPython.display import SVG
print(rdBase.rdkitVersion) 
 
taxol = pcp.get_compounds('taxol', 'name')
taxol = taxol[0]
type(taxol) ### pubchempy.Compound
taxol = Chem.MolFromSmiles(taxol.canonical_smiles)
type(taxol) ### rdkit.Chem.rdchem.Mol

### Morganフィンガープリント
bitI_morgan = {}
fp_morgan = AllChem.GetMorganFingerprintAsBitVect(taxol, 2, bitInfo=bitI_morgan)


### RDKitフィンガープリント
bitI_rdkit = {}
fp_rdkit = Chem.RDKFingerprint(taxol, bitInfo=bitI_rdkit)

print(fp_morgan.GetNumBits(),fp_morgan.GetNumOnBits()) ### 2048 86
print(len(bitI_morgan)) ### 86
 
print(len(fp_rdkit), len(bitI_rdkit.keys())) ### (2048, 1444)
 
for key in list(bitI_morgan.keys())[:5]:
    print(bitI_morgan[key])
    
    
### Morgan可視化   
morgan_turples = ((taxol, bit, bitI_morgan) for bit in list(bitI_morgan.keys())[:12])
Draw.DrawMorganBits(morgan_turples, molsPerRow=4, legends=['bit: '+str(x) for x in list(bitI_morgan.keys())[:12]])

### RDKit可視化
rdkit_turples = ((taxol, bit, bitI_rdkit) for bit in list(bitI_rdkit.keys())[:12])
Draw.DrawRDKitBits(rdkit_turples, molsPerRow=4, legends=['bit: '+str(x) for x in list(bitI_rdkit.keys())[:12]])
