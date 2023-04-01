#!/bin/bash --login

set -euo pipefail

if [ -f "04-run_bbsplit.log" ]
then
    RESUME_PARAM="--batch-resume 04-run_bbsplit.log"
else
    RESUME_PARAM=""
fi

OUTDIR="${PWD}/output/bbsplit"
mkdir -p ${PWD}/output/bbsplit/ArME14
mkdir -p ${PWD}/output/bbsplit/CaFrontier
mkdir -p ${PWD}/output/bbsplit/stats


code/slurm_scripts/bin/pt --nparams 2 "cd work/bbsplit && bbsplit.sh -eoom -Xmx14g threads=8 build=1 in='{0}' in2='{1}' usejni=t maxindel=200000 ambiguous=best ambiguous2=best basename='${OUTDIR}/%/{@bp}#.fq.gz' refstats='${OUTDIR}/stats/{@bp r/_R/}.refstats' ihist='${OUTDIR}/stats/{@bp r/_R/}.ihist' indelhist='${OUTDIR}/stats/{@bp r/_R/}.indelhist' mhist='${OUTDIR}/stats/{@bp r/_R/}.mhist'" ${PWD}/output/trimmed/*{R1,R2}.fq.gz \
| code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 8 \
  --mem 16G \
  --batch-dry-run \
  --batch-module miniconda3/latest \
  --batch-condaenv ${PWD}/condaenv \
  --partition work \
  --time 24:00:00 \
  --job-name bbsplit \
  ${RESUME_PARAM}

