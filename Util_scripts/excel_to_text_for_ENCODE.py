
import pandas as pd
import numpy
import xlrd
import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
description="""

excel_to_text_for_ENCODE.py

This script enables the conversion of a metadata objects spreadsheet required for 
data submission into ENCODE portal from the excel form into .txt files for each of 
the sheet present in the spreadsheet. Each sheet represents each metadata object. 
The default sheet options are : [genetic_modification, biosample, library, 
functional_characterization_exp, replicate, file]. The script is WIP as it 
incorporates only for the above mentioned 6 metadata objects for now.

""")

###############################
# required arguments:

parser.add_argument("-i", "--excel-file", help="""required, file path to the metadata objects excel spreadsheet for ENCODE data submission.""", 
                        required=True)
parser.add_argument("-o", "--out-directory", help="""required, path to the directory where the output .txt files should be located.""")
###############################
# optional arguments:

parser.add_argument("--sheet-names", type=str, nargs='*', default=['genetic_modification','biosample','library','functional_characterization_exp','replicate','file','files','functional_characterization_ser'], help="""required, names of the sheets present in the metadata objects excel spreadsheet.""")

args = parser.parse_args()

for sheet in args.sheet_names:
    df = pd.read_excel(args.excel_file, sheet_name = sheet)
    df.fillna('', inplace=True)
    if (sheet=='functional_characterization_exp'):
        sheet = 'functional_characterization_experiment'
    df.to_csv('%s/%s.txt' % (args.out_directory, sheet), sep='\t', index=False)
