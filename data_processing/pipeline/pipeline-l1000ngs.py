#################################################################
#################################################################
############### L1000 Notebook Generation Server ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, h5py, os
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support3 as S
import L1000Ngs as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
l1000fwd_file = 'rawdata.dir/CD_signatures_LM_42809x978.gctx'
probe_metadata_file = 'rawdata.dir/Probes_L1000_metadata.csv'

##### 2. R Connection #####
r.source('pipeline/scripts/l1000ngs.R')

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Get signatures
#############################################

@follows(mkdir('s1-processed.dir'))

@merge((l1000fwd_file, probe_metadata_file),
       's1-processed.dir/signatures.h5')

def processSignatures(infiles, outfile):

	# Parse GCTx
	gct_data = parse(infiles[0])

	# Get probes from GCT
	gct_probes = [x[2:-1] for x in gct_data.col_metadata_df.index]

	# Get probe metadata
	probe_dataframe = pd.read_csv(infiles[1], index_col='pr_id').reindex(gct_probes).reset_index()
	probe_dataframe.head()

	# Prepare outfile
	f = h5py.File(outfile, 'w')

	# Data
	data_grp = f.create_group('data')
	data_grp.create_dataset('cd', data=gct_data.data_df.T.values, chunks=True, compression="gzip")

	# Gene metadata
	gene_grp = f.create_group('meta/gene')
	for label, colData in probe_dataframe.items():
		gene_grp.create_dataset(label, data=colData.values, dtype=h5py.special_dtype(vlen=str))

	# Sample metadata
	sample_grp = f.create_group('meta/sample')
	for label, colData in gct_data.row_metadata_df.reset_index()[['rid', 'batch', 'cell_id', 'pert_desc', 'pert_dose', 'pert_id', 'pert_time']].items():

		# Get values and data type
		if label in ['pert_time', 'pert_dose']:
			sample_grp.create_dataset(label, data=colData.values)
		else:
			sample_grp.create_dataset(label, data=[x[2:-1].encode('utf-8') for x in colData.values], dtype=h5py.special_dtype(vlen=str))

	# Close
	f.close()

#############################################
########## 2. Get Metadata
#############################################

@transform(processSignatures,
		   suffix('s.h5'),
		   '_metadata.txt')

def getSignatureMetadata(infile, outfile):

	# Read HDF5
	f = h5py.File(infile, 'r')

	# Get data
	metadata_dict = {key: value.value for key, value in f.get('meta/sample').items()}

	# Convert to dataframe
	metadata_dataframe = pd.DataFrame(metadata_dict)

	# Save
	metadata_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 3. Copy Data
#############################################

@transform((processSignatures, getSignatureMetadata),
		   regex(r'.*/(.*)'),
           r'../website/app/static/data/\1')

def copyData(infile, outfile):

	# Copy
	os.system('cp {infile} {outfile}'.format(**locals()))

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
