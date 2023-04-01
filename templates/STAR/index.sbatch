#!/bin/bash --login
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --account=y95
#SBATCH --time=4:00:00
#SBATCH --job-name=star_build

module load miniconda3/latest
conda activate ./condaenv

set -euo pipefail

zcat input/ArME14.fasta.gz > work/ArME14.fasta

zcat input/ArME14-genes.gff3.gz \
| awk '/^#/ {print} $4 < $5 {print}' \
| awk -F '\t' -v OFS="\t" '$3 == "CDS" {$3="exon"; print} $3 != "CDS" {print}' \
| gt gff3_to_gtf -force -o work/ArME14-genes.gtf

mkdir -p work/star_index_ArME14

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL STAR \
    --runMode genomeGenerate \
    --genomeDir ${PWD}/work/star_index_ArME14 \
    --genomeFastaFiles work/ArME14.fasta \
    --sjdbGTFfile work/ArME14-genes.gtf \
    --genomeSAindexNbases 11 \
    --runThreadN "${NCPUS}"

