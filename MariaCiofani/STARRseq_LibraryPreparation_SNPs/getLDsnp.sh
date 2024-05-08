#!/bin/bash
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
python /data/reddylab/Revathy/collabs/Maria/human-th-ms/notebooks/getLD_SNPs.py -i /data/reddylab/Revathy/collabs/Maria/human-th-ms/data/snp_data/source_snps/ebi_ibd_uniqueSNP.txt -o data/reddylab/Revathy/collabs/Maria/human-th-ms/data/snp_data/source_snps/all_LD_IBD_SNPs.tab
