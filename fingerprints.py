import pandas as pd
import numpy as np
import gzip
# this is for the sketchy molecular fragment stuff
# import os

from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import AllChem
from rdkit.Chem import FragmentCatalog
from rdkit import RDConfig

# Number of bits for vectors
NBits = 256

def export_fps(fps, outfile):
	df = pd.DataFrame()
	for fp in fps:
		unpacked = pd.DataFrame(data=np.array([fp[i] for i in range(len(fp))]).reshape(1, -1))
		df = pd.concat((df, unpacked))	
	df.to_csv(outfile)
	print df

# Topological fingerprints
def top(outfile):
	top_fps = [rdmolops.RDKFingerprint(m, fpSize=NBits) for m in ms]
	export_fps(top_fps, outfile)

# Morgan fingerprints
def morgan(outfile):
	morgan_fps = [AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=NBits) for m in ms]
	export_fps(morgan_fps, outfile)

# MACCS keys
def maccs(outfile):
	maccs_fps = [MACCSkeys.GenMACCSKeys(m) for m in ms]
	export_fps(maccs_fps, outfile)

# Atom pairs and topological torsions
def pair(outfile):
	pair_fps = [Pairs.GetAtomPairFingerprintAsBitVect(m) for m in ms]
	export_fps(pair_fps, outfile) 

# Read in train and test as Pandas DataFrames
df_train = pd.read_csv("train.csv.gz", nrows=10000, compression='gzip')

#print "data done. starting SMILES..."

# Convert SMILES to molecules
ms = [Chem.MolFromSmiles(s) for s in df_train["smiles"].values]

#print "SMILES done. starting top..."
top('top.csv')
#print "top done. starting morgan..."
morgan('morgan.csv')
#print "morgan done. starting maccs..."
maccs('maccs.csv')
#print "maccs done"
# THIS IS VERY SPARSE AND HUMUNGOUS
# pair('pair.csv')

# Molecular fragments
# fName = os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
# fparams = FragmentCatalog.FragCatParams(1,6,fName)
# fcat = FragmentCatalog.FragCatalog(fparams)
# fpgen = FragmentCatalog.FragFPGenerator()
# frags = [fpgen.GetFPForMol(fpgen,m,fcat) for m in ms] # BROKEN