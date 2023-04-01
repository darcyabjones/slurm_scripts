#!/bin/bash --login
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --export=NONE
#SBATCH --account=y95
#SBATCH --time=3:00:00
#SBATCH --job-name=merge_stringtie

module load miniconda3/latest
conda activate ./condaenv

set -euo pipefail

ls output/stringtie_ArME14/*.gtf > work/stringtie_ArME14_samples.txt

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL stringtie --merge \
    -G work/ArME14-genes.gff3 \
    -o output/stringtie_ArME14.gtf \
    work/stringtie_ArME14_samples.txt

mkdir -p "output/stringtie_gffcompare"

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL gffcompare \
    -T \
    -r work/ArME14-genes.gff3 \
    -o output/stringtie_gffcompare/ArME14 \
    output/stringtie_ArME14.gtf

