#!/bin/bash --login
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --account=y95
#SBATCH --time=3:00:00
#SBATCH --job-name=count_deseq

module load miniconda3/latest
conda activate ./condaenv

set -euo pipefail

NCPUS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL prepDE.py -i output/stringtie_assembled_ArME14_counts/sample_lst.txt -l 100 -g output/stringtie_assembled_ArME14_counts/gene_count_matrix.csv -t output/stringtie_assembled_ArME14_counts/transcript_count_matrix.csv
srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL prepDE.py -i output/stringtie_assembled_CaFrontier_counts/sample_lst.txt -l 100 -g output/stringtie_assembled_CaFrontier_counts/gene_count_matrix.csv -t output/stringtie_assembled_CaFrontier_counts/transcript_count_matrix.csv

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL prepDE.py -i output/stringtie_genes_ArME14_counts/sample_lst.txt -l 100 -g output/stringtie_genes_ArME14_counts/gene_count_matrix.csv -t output/stringtie_genes_ArME14_counts/transcript_count_matrix.csv
srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL prepDE.py -i output/stringtie_genes_CaFrontier_counts/sample_lst.txt -l 100 -g output/stringtie_genes_CaFrontier_counts/gene_count_matrix.csv -t output/stringtie_genes_CaFrontier_counts/transcript_count_matrix.csv

mkdir -p output/feature_counts

gffread work/ArME14-genes.gff3 -T \
  | awk -F "\t" -v OFS="\t" '{$9=gensub(/gene_id \"([A-Z0-9_]*).*$/, "gene_id \"\\1\";", "g", $9); print}' \
  | awk -F "\t" -v OFS="\t" '$3 == "CDS" {$3="exon"} {print}' \
  | awk -F "\t" -v OFS="\t" '$9 !~ /gene_id/ {gene_id=gensub(/transcript_id \"([A-Z0-9_]*).*;/, "gene_id \"\\1\";", "g", $9); $9=$9" "gene_id} {print}' \
  > work/ArME14-genes2.gtf

awk -v OFS="\t" -F "\t" '$3 == "mRNA" {tid=gensub(/^.*ID=([^;]*).*$/, "\\1", "g", $9); gid=gensub(/^.*Parent=([^;]*).*$/, "\\1", "g", $9); print tid, gid}' work/Carientinum-nuclear.gff3 > work/CaTID2GID.tsv

gffread --keep-genes work/Carientinum-nuclear.gff3 -T \
    | awk -v OFS="\t" -F "\t" '
    BEGIN {while ((getline line< "work/CaTID2GID.tsv") > 0) {split(line, sline, "\t"); tid2gid[sline[1]] = sline[2]; gid2tid[sline[2]] = sline[1];} }
    $3 != "CDS" {
        tid=gensub(/^.*transcript_id\s+"([^"]+)".*$/, "\\1", "g", $9);
        gid=tid2gid[tid];
        $9="transcript_id \""tid"\"; gene_id \""gid"\";";
        if ( gid != "" ) {
            print
        }
    }' \
  > work/Carientinum-nuclear2.gtf

srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL featureCounts -s 2 -p -T 1 -a work/ArME14-genes2.gtf -o output/feature_counts/feature_counts_ArME14.txt work/stringtie_genes_ArME14/*.bam
srun --cpus-per-task "${SLURM_CPUS_PER_TASK:-1}" --export=ALL featureCounts -s 2 -p -T 1 -a work/Carientinum-nuclear2.gtf -o output/feature_counts/feature_counts_CaFrontier.txt work/stringtie_genes_CaFrontier/*.bam

