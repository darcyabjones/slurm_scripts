#!/bin/bash --login

set -euo pipefail

SCRIPT="$(basename $0)"
LOGNAME="${SCRIPT%.*}.log"

INDIR="input/fastq"
OUTDIR="output/falco"
mkdir -p "${OUTDIR}"

if [ -f "${LOGNAME}" ]
then
    RESUME_PARAM="--batch-resume '${LOGNAME}'"
else
    RESUME_PARAM=""
fi

# Note, I also ran this separately with input/resequencing_fastq

pt 'falco --outdir output/falco/{bee} {}' "${INDIR}/"*.{fq,fastq,fq.gz,fastq.gz} \
| sbatch_jobarray.sh \
  --cpus-per-task 1 \
  --ntasks 8 \
  --batch-pack \
  --batch-dry-run \
  --batch-module miniconda3/latest \
  --batch-condaenv ${PWD}/condaenv \
  --partition work \
  --time 2:00:00 \
  --job-name falco \
  ${RESUME_PARAM}
