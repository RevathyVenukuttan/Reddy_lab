#!/bin/sh

sra_ids="SRR10311789 SRR10311790 SRR10311791 SRR10311792"
for sra_id in ${sra_ids}; do
    fastq-dump ${sra_id}
done
