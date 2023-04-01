#!/bin/bash --login

set -euo pipefail

if [ -f "02a-align_ArME14.log" ]
then
    RESUME_PARAM="--batch-resume 02a-align_ArME14.log"
else
    RESUME_PARAM=""
fi

code/slurm_scripts/bin/pt --nparams 2 "code/star.sh work/star_index_ArME14 input/ArME14.fasta.gz {0} {1} {@bp r/_\d+_L\d+_R$/} {@bp r/_R/}" input/fastq/ArME14/*{R1,R2}.fq.gz \
| code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 4 \
  --mem 8G \
  --batch-dry-run \
  --batch-module miniconda3/latest \
  --batch-condaenv ${PWD}/condaenv \
  --partition work \
  --time 4:00:00 \
  --job-name star_align \
  ${RESUME_PARAM}

