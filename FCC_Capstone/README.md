## Contents in notebooks:

1. **[FinalCode_v1](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/FinalCode_v1.ipynb):**

 - FINAL HEATMAPS, PROFILE PLOTS, SUMMARY HEATMAPS, SUMMARY SCORE FILE for scatterplots

2. **[FinalCode_v2](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/FinalCode_v2.ipynb):**

 - Unionset of active enhancer calls compared to common background/negative set. Heatmaps, profile plots, Scatterplots created. For Scatterplots comparing genomic context with Dustin's sequence analysis, average signal value taken across 3 assays for Dustin's sequence analysis score.

 - FeatureMatrix creation for all 3 assays where rows indicate active enhancer calls (from Junke's z-score processed calls) that was used in the heatmaps and columns have genomic scores pulled out from the computeMatrix for each TF/HistoneMark for every 10bp tile (+/-2kb region from center). The matrix has shape - 5000x(400x19) 

 - Create an inactive set where the inactive regions can be used as negative controls with confidence of not possessing any kind of activity - not anchored on ATAC peak calls from the inactive bins of z-score processed enhancer calls from Junke

3. **[ASHG_poster.ipynb](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/ASHG_poster.ipynb):**

 - Notebook was used to generate all plots necessary for the ASHG poster.

4. [Heatmaps_v7](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Heatmaps_v7.ipynb):

 - Versions of Scatterplots comparing TF score between Dustin's sequence models and genomic context analysis score from Profile plots. Also compared different bigwigs of NFE2 to identify the bigwig with better signal

5. [Heatmaps_v6](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Heatmaps_v6.ipynb):

 - Comparison of REST encode .bw files, Heatmap for Tiling MPRA, Heatmaps for repressed hits, Add RNA tpm to the repressed heatmaps (ASTARR) to check the reason for the hole in the middle. Version 1 of scatterplots comparing sequence models to genomic context. 
 
6. [Heatmaps_v3](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Heatmaps_v3.ipynb): 

 - Iterations of Heatmaps with versions of common background sets and summits of peaks for the enhancer call data

7. [Heatmaps_v2](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Heatmaps_v2.ipynb):

 - Previous iterations of heatmaps (separately for all assays - ASTARR, wgSTARR, LentiMPRA, tilingMPRA). Enhancer calls used - previous versions. Split these heatmaps to promoter and distal. GC content matching to match the inactive with the active. Version 1 of common background set. 

8. [Heatmaps_v1](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Heatmaps_v1.ipynb):

 - Version 1 trials of heatmaps, adding new TFs/Histone marks, adding CRISPR data, trying out CRISPRi/k data to see the enrichment patterns, trying random sampling of active enhancer calls instead of using the whole set, separations into promoter and distal. 

9. [20221212_fcc_heatmaps_update.ipynb](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/20221212_fcc_heatmaps_update.ipynb):

 - Initial iteration of profile_plots (multiple iterations present that adjusts for changes made to data or changes in versions of enhancer calls from Junke)

10. [Profile_plots_v2](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Profile_plots_v2.ipynb):

 - Initial versions of profile_plots continued (multiple iterations present that adjusts for changes made to data or changes in versions of enhancer calls from Junke)

11. [ChIPseeker_annotation](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Chipseeker_annotation.ipynb): 

 - All of the ChIPseeker annotation performed on different sets of data at different points.

12. [Common_background_set_v1](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Common_background_set_v1.ipynb): 

 - Initial iterations of creation of common background for all th 3 assays tested

13. [input_vs_inactive peaks comparison](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/input_vs_inactive%20peaks%20comparison.ipynb): 

 - Initial iterations of testing the GC content matching script between the active set and inactive set used

14. [Inactive_Peaks_Data_Wrangling](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Inactive_Peaks_Data_Wrangling.ipynb):

 - Initial iterations of data wrangling for the inactive peaks from Junke's enhancer call pipeline and ATAC-seq data 

15. [Bubble_plots](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Bubble_plots.ipynb):

 - Iterations for creating bubble plots with different sets of assay data

16. [Correlation_Enhancer_Calls.ipynb](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Correlation_Enhancer_Calls.ipynb):
 
 - Explore the correlation using heatmaps, pairplots between the enhancer calls (earlier iterations from Junke's enhancer call pipeline) for ATAC-STARRseq, wgSTARRseq, LentiMPRA and TewheyMPRA data. Also used Jaccard Index to compare the assays

17. [Correlation_biwigs.ipynb](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Correlation_bigwigs.ipynb):

 - Correlation between signal bigwigs of ASTARR and MPRA assays. Used multiBigWigSummary from Deeptools package to plot correlation heatmaps.
 - multiBigwigSummary for ASTARR and MPRA at 3 loci (GATA,MYC&FADS) anchored at all cCREs

18. [Wilcoxon_test.ipynb](https://github.com/RevathyVenukuttan/Reddy_lab/blob/main/FCC_Capstone/Wilcoxon_test.ipynb):

 - Wilcoxon rank sum test between the active and inactive genomic scores for a list of TFs/Histone marks that showed significance (on visual examination) in the heatmaps/profile plots. Purpose of this was to determine statistical significance with a p-value.



## Path to directories/files in HARDAC:

#### 1. Final ComputeMatrix file: 

`/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/10032022/ASTARR_LentiMPRA_wgSTARR.GC_content.random5000.V6.refpoint.mat.gz`

#### 2. Final heatmaps from computeMatrix:

`/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/10032022/ASTARR_LentiMPRA_wgSTARR.GC_content.random5000.V6.png`

#### 3. Final Profile Plots:

ASTARR: `/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/profile_plots/astarr.active_v_inactive.zoomed.V3.pdf`

LentiMPRA: `/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/profile_plots/lentiMPRA.active_v_inactive.zoomed.V3.pdf`

wgSTARR: `/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/profile_plots/wgstarr.active_v_inactive.log2.zoomed.V3.pdf`

#### 4. Final Summary Heatmaps:

`/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/bubble_plots/ASTARR_Lenti_wgSTARR.Promoter_Distal.AUC_summary_heatmaps.V2.pdf`

#### 5. Summary table with scores: 

`/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/bubble_plots/ASTARR_Lenti_wgSTARR.Promoter_Distal.AUC_summary_heatmaps.signal_table.txt`

#### 6. Summary table with scores - including TFs from Dustin's sequence analysis:

`/data/reddylab/Revathy/collabs/Jamborees/fcc_capstone/data/enhancer_call_comparison/bubble_plots/ASTARR_Lenti_wgSTARR.Promoter_Distal.AUC_summary_heatmaps.signal_table.dustin_seq_feature.V2.txt`


