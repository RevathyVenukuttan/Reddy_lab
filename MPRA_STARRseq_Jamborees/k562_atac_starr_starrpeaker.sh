#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
input_dir=/data/reddylab/Revathy/collabs/Jamborees/05172021/
output_dir=/data/reddylab/Revathy/collabs/Jamborees/07192021/peak_calls
file_dir=/data/reddylab/Revathy/collabs/Jamborees/07192021/data
output_files=(\ KS91_K562_hg38_ASTARRseq_Output_rep1.f3q10.sorted.with_umis.chr6.dups_marked.bam \
                KS91_K562_hg38_ASTARRseq_Output_rep2.f3q10.sorted.with_umis.chr6.dups_marked.bam \
                KS91_K562_hg38_ASTARRseq_Output_rep3.f3q10.sorted.with_umis.chr6.dups_marked.bam \
                KS91_K562_hg38_ASTARRseq_Output_rep4.f3q10.sorted.with_umis.chr6.dups_marked.bam \
                KS91_K562_hg38_ASTARRseq_Output_rep5.f3q10.sorted.with_umis.chr6.dups_marked.bam \
                KS91_K562_hg38_ASTARRseq_Output_rep6.f3q10.sorted.with_umis.chr6.dups_marked.bam \
                )
output_file=${output_files[${SLURM_ARRAY_TASK_ID}]}
identifier=$(echo $output_file | cut -d "." -f1)
replicate=$(echo $identifier | cut -d "_" -f6)
input_file=$(/bin/ls /data/reddylab/Revathy/collabs/Jamborees/07192021/data/*_Input_${replicate}.*.bam)
starrpeaker --prefix ${output_dir}/${identifier} \
--chromsize ${input_dir}/GRCh38.chrom.sizes.simple.sorted \
--blacklist ${input_dir}/ENCODE_blacklist_GRCh38_ENCFF419RSJ_merged.bed \
-i ${input_file} \
--output ${output_file} \
--cov ${input_dir}/STARRPeaker_cov_GRCh38_ucsc-gc-5bp.bw \
      ${input_dir}/STARRPeaker_cov_GRCh38_linearfold-folding-energy-100bp.bw \
      ${input_dir}/STARRPeaker_cov_GRCh38_gem-mappability-100mer.bw \
--threshold 0.05
