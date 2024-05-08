#!/usr/bin/env python
import pandas as pd
import numpy as np
import subprocess
from Bio.Seq import Seq
from random import random
from scipy.stats import itemfreq

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
description="""

generate_non_targeting_gRNAs.py

Given a list of targeting guides, this script generates a list of 
non-targeting guides depending on the number of samples specified 
in the --n_random_samples argument. 

""")


#####################################################################################################
## required args: 

parser.add_argument("--gRNAs", help="""required, .txt file containing the list of targeting guide RNAs""",
                     required=True)

parser.add_argument("-n", "--n_random_samples", type=int,
                    help="""required, number of random samples of non-targeting guides to be 
                    generated from the list of targeting guides given""", required=True)

parser.add_argument("-o", "--output_dir", help="""required, path to the directory where the output of the script
                    will be written into.""", required=True)

args = parser.parse_args()

#######################################################################################################

## optional args:

parser.add_argument("--path_tmp_files", type=str, default=args.output_dir,
                    help="""path to all temporary files created for and during bowtie alignment""",
                    dest="path_tmp_files")
parser.add_argument("--bowtie_index", type=str, default='/data/reddylab/Reference_Data/Genomes/hg19/hg19',
                    help="""path to the directory that contains all the index files required for bowtie 
                    to perform alignment for the particular reference genome. The default is for the hg19
                    human genome but can be changed to any other genome.""", dest="bowtie_index")
parser.add_argument("--runname", type=str, default='NT_guides_generation',
                    help="""name of all the temporary files that will be created using the bowtie
                    alignment""", dest="runname")
parser.add_argument("--btThreshold", type=int, default=31,
                    help="""Threshold for bowtie""", dest='btThreshold')
parser.add_argument("--mFlag", type=int, default=1,
                    help="""Threshold for number of mismatches allowed in alignment""", dest='mFlag')

########################################################################################################


args = parser.parse_args()

global DNA_ALPHABET
DNA_ALPHABET = ['A','T','G','C']

def extract_freq_from_np_uniq(a):
    freqs = dict(zip(DNA_ALPHABET, np.zeros(4)))
    kk,vv = np.unique(a, return_counts=True)
    freqs.update(dict(zip(kk, (1.*vv)/a.shape[0])))
    return [freqs[k] for k in DNA_ALPHABET]


def subsample_alphabet_given_freqs(freqs):
    return np.random.choice(DNA_ALPHABET, p=freqs)


def outputTempBowtieFastq(libraryTable, outputFileName):
    phredString = 'I4!=======44444+++++++' #weighting for how impactful mismatches are along sgRNA sequence 
    with open(outputFileName,'w') as outfile:
        for name, row in libraryTable.iterrows():
            outfile.write('@' + name + '\n')
            outfile.write('CCN' + str(Seq(row['sequence'][1:]).reverse_complement()) + '\n')
            outfile.write('+\n')
            outfile.write(phredString + '\n')

            
### setting arguments to variables

n_random_samples = args.n_random_samples
file_name = args.gRNAs.split('/')[-1].split('.')[0]
fqFile = args.path_tmp_files + '/' + 'temp.fastq'
alignmentList = [(args.btThreshold, args.mFlag, args.bowtie_index, args.runname)]


#### reading the input targeting GRNA file
grnas = pd.read_csv(args.gRNAs, sep=',', names = ["gRNA"])
grnas["gRNA"] = grnas["gRNA"].str.upper()



freqs = np.apply_along_axis(extract_freq_from_np_uniq, 0, 
                            np.vstack([list(t) for t in grnas.gRNA.values.T]))

nt_final_set = set()
alignmentColumns = []
while len(nt_final_set)<n_random_samples:
    random_samples = [''.join(np.apply_along_axis(subsample_alphabet_given_freqs, 0,  freqs)) for ii in range(n_random_samples)]
    negTable = pd.DataFrame(random_samples, 
                            index=['NT_%d' %i  for i in range(1, len(random_samples)+1)], 
                            columns = ['sequence'])
    outputTempBowtieFastq(negTable, fqFile)

    for btThreshold, mflag, bowtieIndex, runname in alignmentList:

        alignedFile = args.path_tmp_files + '/' + args.runname + '_aligned.txt'
        unalignedFile = args.path_tmp_files + '/' + args.runname + '_unaligned.fq'
        maxFile = args.path_tmp_files + '/' + args.runname + '_max.fq'

        bowtieString = '/nfs/software/helmod/apps/Core/bowtie/1.1.1-fasrc01/bowtie -n 3 -l 5 -e ' + \
            str(btThreshold) + \
            ' -m ' + str(mflag) + \
            ' --nomaqround -a --tryhard -p 16 --chunkmbs 256 ' + \
            bowtieIndex + \
            ' --suppress 5,6,7 --un ' + \
            unalignedFile + \
            ' --max ' + maxFile + \
            ' -q ' + \
            fqFile + ' ' + \
            alignedFile

        subprocess.call(bowtieString, shell=True)

        #read unaligned file for negs, and then don't flip boolean of alignmentTable
        with open(unalignedFile) as infile:
            sgsAligning = set()
            for i, line in enumerate(infile):
                if i%4 == 0: #id line
                    sgsAligning.add(line.strip()[1:])

        alignmentColumns.append(negTable.apply(lambda row: row.name in sgsAligning, axis=1))

    alignmentTable = pd.concat(alignmentColumns,axis=1, keys=list(zip(*alignmentList))[3])

    nt_final_set = set(list(nt_final_set) + list(negTable[alignmentTable.values].values.T[0]))

nt_final_df = pd.DataFrame(list(nt_final_set), columns=['nt_sequence'])
nt_final_df.to_csv(args.output_dir + '/' + file_name + '_NT_guides.txt', index=False)