README.md

## Methods for Library Creation for SNP-Autoimmune Disorders Project

### GWAS Data Acquisition
We downloaded GWAS data from the EBI-GWAS catalog using the EBI API. This catalog is a curated collection of all published genome-wide association studies. Lead variants for each autoimmune disorder (Multiple Sclerosis, Inflammatory Bowel Disorders, Rheumatoid Arthritis, Psoriasis) were identified using the 'Disease/Trait' column from the GWAS catalog data table. 

### Linkage Disequilibrium (LD) Expansion
The lead SNPs were expanded using Linkage Disequilibrium (LD). Initially, the lead SNPs were categorized based on the 'Ancestry' annotation in the GWAS catalog. Given that the catalog included more ancestry subtypes than supported by the ENSEMBL LD calculator tool, we grouped all subtypes into four main ancestry categories: European (EUR), East Asian (EAS), African (ASW), and South Asian (SAS). LD expansion was performed separately for these four lists, and the results were then concatenated. SNPs with an R² ≥ 0.8 were retained. 

All code assocaiated GWAS data acquistion and LD expansion can be found in [SNP_LD_expansion.ipynb](MariaCiofani/STARRseq_LibraryPreparation_SNPs/SNP_LD_expansion.ipynb)

### SNP Overlay with ATAC-Seq Peaks
The expanded list of SNPs (lead SNPs + LD-expanded SNPs) was overlayed with Th1/Th17 ATAC-seq peaks to identify SNPs located in open chromatin regions. This process yielded a total of 2,330 unique variants overlapping with the ATAC-seq peaks. 

All code assocaiated SNP overlay with ATAC-Seq peaks can be found in [linking_snps_loops_genes-all_snps_ancestry.ipynb](MariaCiofani/STARRseq_LibraryPreparation_SNPs/linking_snps_loops_genes-all_snps_ancestry.ipynb)

### Fragment Design
The STARR-seq library was designed with 246 bp fragments, including flanking sequences at both the 5' and 3' ends. Each focal variant was tested for all its alleles at three specific positions within the 246 bp fragment: the 82nd position (1/3rd), the 123rd position (1/2nd), and the 164th position (2/3rd).

### Haplotypes and Multi-Allelic SNPs
We identified a total of 2,330 SNPs that overlapped with the ATAC-seq peak regions. Within this set, we noted that some variants were in close proximity to each other (e.g., 15 SNPs within a 300 bp fragment). This suggests the presence of haplotypes, and the library was designed to test these haplotypes collectively. Also, more than 50% of the SNPs were found to be multi-allelic (with more than two alleles, and in some cases, up to four alleles) based on variant information (allele_string, minor_allele, MAF) from ENSEMBL. To address this, we consulted the 1000 Genomes Project (1KGP) and gnomAD for allele-specific information. It was found that for most multi-allelic variants, the third and fourth alleles had negligible allele frequencies and were disregarded.

Code for this is in given notebook: [allele_frequency.ipynb](MariaCiofani/STARRseq_LibraryPreparation_SNPs/allele_frequency.ipynb)

### Library Creation Strategy
(Done by Dr. Bill Majoros)

Substitution Method (for variants not in 1000 Genomes Project, 10% of cases):
1.	Construct all haplotypes for the combined ancestry group for the variant.
2.	Select the most common haplotype as the background haplotype.
3.	Substitute each annotated allele into the background haplotype.
4.	Ensure at least one annotated allele matches the hg38 reference sequence.

Haplotype Method (for variants present in 1000 Genomes Project, 90% of cases):
1.	Construct all haplotypes for the combined ancestry group for the variant.
2.	Identify haplotypes containing each annotated allele.
3.	Select a representative haplotype for each allele based on its frequency in the ancestry group.
4.	Discard all annotated alleles that is not present in any haplotype.
5.	Verify that the annotated reference allele matches the hg38 sequence.

Code for Library construction (Initial iterations before haplotype mapping and data wrangling after haplotype mapping) along with sanity checks after haplotype mapping is in the notebook: [library_preparation.ipynb](MariaCiofani/STARRseq_LibraryPreparation_SNPs/library_preparation.ipynb)

There were duplications in the library (some of the oligos along with their IDs were duplicated) whose source had to be tracked down. All code associated with this analysis/sanity check of the library is in [library_preparation_recheck.ipynb](MariaCiofani/STARRseq_LibraryPreparation_SNPs/library_preparation_recheck.ipynb)





