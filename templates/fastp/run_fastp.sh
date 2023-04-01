#!/bin/bash --login

set -euo pipefail

SCRIPT="$(basename $0)"
LOGNAME="${SCRIPT%.*}.log"

INDIR="input/fastq"
STATS_OUTDIR="output/fastp"
TRIMMED_OUTDIR="output/trimmed"
mkdir -p "${STATS_OUTDIR}"
mkdir -p "${TRIMMED_OUTDIR}"

if [ -f "${LOGNAME}" ]
then
    RESUME_PARAM="--batch-resume '${LOGNAME}'"
else
    RESUME_PARAM=""
fi

# Note: I also ran this separately using reads in input/resequencing_fastq

pt --nparams 2 "fastp --in1 {0} --in2 {1} --out1 '${TRIMMED_OUTDIR}/{0b}' --out2 '${TRIMMED_OUTDIR}/{1b}' --qualified_quality_phred 5 --length_required 50 --detect_adapter_for_pe --json '${STATS_OUTDIR}/{@bp r/_R/}-fastp.json' --html '${STATS_OUTDIR}/{@bp r/_R/}-fastp.html' -R '{@bp r/_R/}' --thread 3" "${INDIR}"/*{R1.fq,R1.fastq,R1.fq.gz,R1.fastq.gz,R2.fq,R2.fastq,R2.fq.gz,R2.fastq.gz} \
| sbatch_jobarray.sh \
  --cpus-per-task 2 \
  --ntasks 6 \
  --batch-pack \
  --batch-dry-run \
  --batch-module miniconda3/latest \
  --batch-condaenv ${PWD}/condaenv \
  --partition work \
  --time 2:00:00 \
  --job-name fastp \
  ${RESUME_PARAM}


