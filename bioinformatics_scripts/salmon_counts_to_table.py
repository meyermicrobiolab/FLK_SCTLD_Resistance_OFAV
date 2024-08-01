##### I DIDN'T USE THIS SCRIPT TO EXECUTE ON THE COMMAND LINE #######
### I just copied this code into the command line

#$########################################################################
## CODE USED IN THE HPC COMMAND LINE - GENERATE COUNTS FILE   ####
###########################################################################

### NOTE - start a scavenger run before running this code on the HPC

module load python

python

import pandas as pd
import glob

all_files = glob.glob('*.quant/quant.sf') # get the name of every file in a list

first = 1

for filename in all_files:
     if (first == 1): # read in the first file, then make it a dataframe
             master = pd.read_csv(filename, header = 0, sep='\t', usecols=[0,4])
             master.columns = names=["Contig", filename]
			 first = 0
     else: # add subsequent files as additional columns to the master file
             df = pd.read_csv(filename, header = 0, sep='\t', usecols=[0,4])
             df.columns = names=["Contig", filename]
             master = pd.merge(master, df, how = "outer", on="Contig")
			 
master.to_csv('salmon_quant_all_NumReads.tsv', sep = "\t", float_format='%.f') # write out a tab separated file
