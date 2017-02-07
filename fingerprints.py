import pandas as pd
import numpy as np
import gzip
# this is for the sketchy molecular fragment stuff
# import os

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
# from rdkit.Chem.AtomPairs import Pairs # atom pairs stuff
# from rdkit.Chem import FragmentCatalog # molecular fragment stuff
# from rdkit import RDConfig

# Number of bits for vectors
NBits = 256

def export_fps(fps, outfile):
	df_list = []
	for fp in fps:
		unpacked = pd.DataFrame(data=np.array([fp[i] for i in range(len(fp))]).reshape(1, -1))
		df_list.append(unpacked) 
	df = pd.concat(df_list)	
	df.to_csv(outfile)
	print df

# Topological fingerprints
def top(ms,outfile):
	top_fps = [rdmolops.RDKFingerprint(m, fpSize=NBits) for m in ms]
	export_fps(top_fps, outfile)

# Morgan fingerprints
def morgan(ms,outfile):
	morgan_fps = [AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=NBits) for m in ms_test]
	export_fps(morgan_fps, outfile)

# MACCS keys
def maccs(ms,outfile):
	maccs_fps = [MACCSkeys.GenMACCSKeys(m) for m in ms]
	export_fps(maccs_fps, outfile)

# Atom pairs and topological torsions
#def pair(ms,outfile):
#	pair_fps = [Pairs.GetAtomPairFingerprintAsBitVect(m) for m in ms]
#	export_fps(pair_fps, outfile) 

# FOR TRAINING
# Read in train as Pandas DataFrame and convert SMILES to molecules
df_train = pd.read_csv("train.csv.gz", nrows=100000, compression='gzip')
ms_train = [Chem.MolFromSmiles(s) for s in df_train["smiles"].values]

# FOR TESTING
# Read in test as Pandas DataFrame and convert SMILES to molecules
df_test = pd.read_csv("test.csv.gz", compression='gzip')
ms_test = [Chem.MolFromSmiles(s) for s in df_test["smiles"].values]

#top(ms,'top.csv')
morgan(ms_test,'morgan_test.csv')
#maccs(ms,'maccs.csv')
# THIS IS VERY SPARSE AND HUMUNGOUS
# pair(ms,'pair.csv')

# Molecular fragments
# fName = os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
# fparams = FragmentCatalog.FragCatParams(1,6,fName)
# fcat = FragmentCatalog.FragCatalog(fparams)
# fpgen = FragmentCatalog.FragFPGenerator()
# frags = [fpgen.GetFPForMol(fpgen,m,fcat) for m in ms] # BROKEN