#!/usr/bin/env python

import argparse
import requests
import time
from io import StringIO
api_token = 'bb2472bededc'

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
description="""

getLD_SNPs.py

For a given list of SNPs, expand the list based on Linkage Disequilibrium(LD) using the 
NCBI tool LDLink

""")
 
required = parser.add_argument_group('required arguments')

################################################################
## required args:
required.add_argument("-i", "--input",
	                  help="path to input SNP file", required=True)
required.add_argument("-o", "--output",
	                  help="path to output file", required=True)

#################################################################

args = parser.parse_args()
infile = args.input
outfile = args.output

 ### get the LD SNPs using NCBI LDlink tool

def getLD(snps):
    '''Takes a list of SNP rsIDs or hg19 coordinates and returns all of the SNPs in LD'''
    st_time = time.time()
    api_token = 'bb2472bededc'
    problem_ids = []
    begin = 0
    for i, snp in enumerate(snps):
        payload = {'token' : api_token, # Dict of params for HTTP request
          'var' : snp,
          'pop' : 'EUR',
          'r2_d' : 'r2'}
        r = requests.get('https://ldlink.nci.nih.gov/LDlinkRest/ldproxy', params=payload) # Make HTTP request
        if "traceback" in str(r.content):
            problem_ids.append(snp)
            print("SNP %s was not found! Skipping..." % (snp))
            if i == begin:
                begin += 1
                ld_snps = []
        elif i == begin: # If the first query, instantiate the results df
            ld_snps = pd.read_table(StringIO(str(r.content, 'utf-8'))) # Parse the bytestring into a pandas DF
            ld_snps['Query'] = snp
        else: # If not the first query, add onto the existing results df
            more_snps = pd.read_table(StringIO(str(r.content, 'utf-8')))
            more_snps['Query'] = snp
            ld_snps = pd.concat(objs=[ld_snps, more_snps])
        print(time.ctime())
        print("Completed %d / %d SNPs!" % (i+1, len(snps))) # Progress reporting for sanity
        print("Total: %d linked SNPs" % (len(ld_snps)))
        print("Elapsed: %d minutes %d seconds" % (np.floor((time.time() - st_time)/60), (time.time() - st_time)%60))
        print("")
    return ld_snps, problem_ids
 
def read_snp_lines(filename):
    with open(filename, "r") as f:
        lines = [line.rstrip('\n') for line in f.readlines()]
        return lines

ebi_ibd_uniqueSNP = read_snp_lines(infile)
all_LD, problem_ids = getLD(ebi_ibd_uniqueSNP)
all_LD.to_csv(outfile, sep="\t")