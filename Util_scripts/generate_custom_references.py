#!/usr/bin/env python
import pandas as pd
import pybedtools
from pybedtools import BedTool

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
description="""

generate_custom_references.py

Given a list of genomic coordinates, extract a custom reference fasta file and custom genomic annotation file
from the corresponding reference files. It can be used for any genome build/organism. 

""")

###############################################################################################################################
## required args: 

parser.add_argument("-i", "--input_coordinates", help="""required, .bed file containing the list of genomic coordinates""",
                     required=True)

parser.add_argument("-fa", "--reference_fasta", help="""required, .fna file consisting of the fasta sequences for all chromosomes of a particular genome""", required=True)

parser.add_argument("-gtf", "--genomic_annotation", help="""required, .gtf file with genomic annotations for a given genome""", 
                    required=True)

parser.add_argument("-o", "--output-dir", help="""required, path to the directory where all the output files will be saved""",
                   required=True)

###############################################################################################################################

args = parser.parse_args()

regions = BedTool(args.input_coordinates)
ref_fasta = BedTool(args.reference_fasta)
ref_gtf = BedTool(args.genomic_annotation)

region_fasta = regions.sequence(fi=ref_fasta).save_seqs(args.output_dir + '/coordinates_fasta.fa')
region_gtf = ref_gtf.intersect(regions).saveas(args.output_dir + '/coordinates_gtf.gtf')