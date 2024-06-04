## Path to directories/files:

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


## Contents in notebooks:

1. [Heatmaps_v8](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Heatmaps_v8.ipynb): 

 - Unionset of active enhancer calls compared to common background/negative set. Heatmaps, profile plots, Scatterplots created. For Scatterplots comparing genomic context with Dustin's sequence analysis, average signal value taken across 3 assays for Dustin's sequence analysis score.

 - FeatureMatrix creation for all 3 assays where rows indicate active enhancer calls (from Junke's z-score processed calls) that was used in the heatmaps and columns have genomic scores pulled out from the computeMatrix for each TF/HistoneMark for every 10bp tile (+/-2kb region from center). The matrix has shape - 5000x(400x19) 

 - Create an inactive set where the inactive regions can be used as negative controls with confidence of not possessing any kind of activity - not anchored on ATAC peak calls from the inactive bins of z-score processed enhancer calls from Junke


2. [Heatmaps_v7](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Heatmaps_v7.ipynb):

Versions of Scatterplots comparing TF score between Dustin's sequence models and genomic context analysis score from Profile plots. Also compared different bigwigs of NFE2 to identify the bigwig with better signal

3. [Heatmaps_v6](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Heatmaps_v6.ipynb):

Comparison of REST encode .bw files, Heatmap for Tiling MPRA, Heatmaps for repressed hits, Add RNA tpm to the repressed heatmaps (ASTARR) to check the reason for the hole in the middle. Version 1 of scatterplots comparing sequence models to genomic context. 

4. [Heatmaps_v5](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Heatmaps_v5.ipynb):

FINAL HEATMAPS, PROFILE PLOTS, SUMMARY HEATMAPS, SUMMARY SCORE FILE for scatterplots + Profile plots used for the ASHG paper (subset from the final profile plots)

5. [Heatmaps_v3](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Heatmaps_v3.ipynb): 

Iterations of Heatmaps with versions of common background sets and summits of peaks for the enhancer call data

6. [Heatmaps_v2](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Heatmaps_v2.ipynb):

Previous iterations of heatmaps (separately for all assays - ASTARR, wgSTARR, LentiMPRA, tilingMPRA). Enhancer calls used - previous versions. Split these heatmaps to promoter and distal. GC content matching to match the inactive with the active. Version 1 of common background set. 

7. [Heatmaps_v1](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Heatmaps_v1.ipynb):

Version 1 trials of heatmaps, adding new TFs/Histone marks, adding CRISPR data, trying out CRISPRi/k data to see the enrichment patterns, trying random sampling of active enhancer calls instead of using the whole set, separations into promoter and distal. 

8. [ChIPseeker_annotation](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Chipseeker_annotation.ipynb): 

All of the ChIPseeker annotation performed on different sets of data at different points.

9. [Profile_plots_v2](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Profile_plots_v2.ipynb):

Initial versions of profile_plots (multiple iterations present that adjusts for changes made to data or changes in versions of enhancer calls from Junke)

10. [Common_background_set_v1](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Common_background_set_v1.ipynb): 

Initial iterations of creation of common background for all th 3 assays tested

11. [input_vs_inactive peaks comparison](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/input_vs_inactive%20peaks%20comparison.ipynb): 

Initial iterations of testing the GC content matching script between the active set and inactive set used

12. [Inactive_Peaks_Data_Wrangling](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Inactive_Peaks_Data_Wrangling.ipynb):

Initial iterations of data wrangling for the inactive peaks from Junke's enhancer call pipeline and ATAC-seq data 

13. [Bubble_plots](https://jupyterhub.genome.duke.edu/user/rv103/notebooks/reddylab_data/Revathy/collabs/Jamborees/fcc_capstone/Bubble_plots.ipynb):

Iterations for creating bubble plots with different sets of assay data