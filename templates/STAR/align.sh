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

#!/usr/bin/env bash

set -euo pipefail

INDEX="$(realpath "$1")"
REF="$(realpath "$2")"
R1="$(realpath "$3")"
R2=$(realpath "$4")
SAMPLE="$5"
OUTPREFIX="$6"

TMPSAM="${TMPDIR:-/tmp}/$$-star.sam"

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

CODEDIR="${PWD}/code"
WORKDIR="${PWD}/work/star_me14/${OUTPREFIX}"
OUTDIR="${PWD}/output/STAR_ArME14"
STATDIR="${PWD}/output/STAR_ArME14_stats"

mkdir -p "${WORKDIR}" "${OUTDIR}" "${STATDIR}"

cd "${WORKDIR}"

ID=$(zcat "${R1}" \
  | head -n 1 \
  | sed 's/^.*:\([A-Za-z0-9]*\):\([[:digit:]]*\):[[:digit:]]*:[[:digit:]]*:[[:digit:]]* .*$/\1.\2/') || :
PU="${ID}.${SAMPLE}"
LB="${SAMPLE}"
PL="ILLUMINA"

RG="ID:${PU}\tSM:${SAMPLE}\tPL:${PL}\tPU:${PU}\tLB:${LB}"
BNAME="$(basename "${OUTPREFIX}")"

STAR \
    --genomeDir "${INDEX}" \
    --readFilesIn "${R1}" "${R2}" \
    --outFileNamePrefix "${BNAME}." \
    --readFilesCommand zcat \
    --runThreadN "${NCPUS}" \
    --outSAMtype SAM \
    --outMultimapperOrder Random \
    --outSAMstrandField intronMotif \
    --outSJfilterIntronMaxVsReadN 500 2000 \
    --alignIntronMin 20 \
    --alignIntronMax 5000 \
    --peOverlapNbasesMin 20 \
    --outSAMattrRGline "${RG}"

samtools faidx "${REF}"

cat "${BNAME}.Aligned.out.sam" \
  | awk -v strType=2 -f ${CODEDIR}/tagXSstrandedData.awk \
  | samtools fixmate -m - - \
  | samtools sort -@${NCPUS} -T ".samtools_sort$$" --reference "${REF}" -O CRAM,embed_ref -o "${OUTDIR}/${OUTPREFIX}.cram" -

samtools index "${OUTDIR}/${OUTPREFIX}.cram"

samtools view -O BAM "${OUTDIR}/${OUTPREFIX}.cram" \
| bedtools genomecov -bga -split -ibam - \
| sort -k1,1 -k2,2n \
> "genome.bedgraph"

bedGraphToBigWig genome.bedgraph "${REF}.fai" "${OUTDIR}/${OUTPREFIX}.bw"

cp *Log.final.out "${STATDIR}"
cp *Log.out "${STATDIR}"

