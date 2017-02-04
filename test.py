import pandas as pd
import numpy as np
import gzip
#import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdmolops

"""
Read in train and test as Pandas DataFrames
"""
df_train = pd.read_csv("train.csv.gz", nrows=1000, compression='gzip')

# Convert SMILES to molecules
ms = [Chem.MolFromSmiles(s) for s in df_train["smiles"].values]

# Topological fingerprints
# top_fps = [FingerprintMols.FingerprintMol(m) for m in ms]
top_fps = [rdmolops.RDKFingerprint(m, fpSize=256) for m in ms]

# Morgan fingerprints
morgan_fps = [AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=256) for m in ms]

df = pd.DataFrame()
for fp in top_fps:
	unpacked = pd.DataFrame(data=np.array([fp[i] for i in range(len(fp))]).reshape(1, -1))
	df = pd.concat((df, unpacked))	
print df

#df = pd.DataFrame(morgan_fps)
df.to_csv('morgan.csv')
df.to_csv('top.csv')