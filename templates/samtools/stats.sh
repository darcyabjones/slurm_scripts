#!/bin/bash --login

set -euo pipefail

if [ -f "08-samtools_stats.log" ]
then
    RESUME_PARAM="--batch-resume 08-samtools_stats.log"
else
    RESUME_PARAM=""
fi

code/slurm_scripts/bin/pt 'samtools stats --ref-seq work/ArME14.fasta {} > output/STAR_ArME14_stats/{be}_samtools_stats.txt' output/STAR_ArME14/*.cram \
| code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 1 \
  --batch-module miniconda3/latest \
  --batch-condaenv ${PWD}/condaenv \
  --partition work \
  --time 4:00:00 \
  --job-name samtools_stats_me14 \
  ${RESUME_PARAM}
