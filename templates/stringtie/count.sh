#!/bin/bash --login

set -euo pipefail

if [ -f "05a-count_stringtie_assembled_ArME14.log" ]
then
    RESUME_PARAM="--batch-resume 05a-count_stringtie_assembled_ArME14.log"
else
    RESUME_PARAM=""
fi

mkdir -p output/stringtie_assembled_ArME14_counts
mkdir -p work/stringtie_assembled_ArME14


code/slurm_scripts/bin/pt --group "{be r/_.*$/}" "samtools merge -f --reference work/ArME14.fasta -O BAM -o 'work/stringtie_assembled_ArME14/{@:0be r/_.*$/}.bam' {@q} && samtools index 'work/stringtie_assembled_ArME14/{@:0be r/_.*$/}.bam' && stringtie -p 2 -e -B --rf --ref work/ArME14.fasta -G output/stringtie_ArME14.gtf -o 'output/stringtie_assembled_ArME14_counts/{@:0be r/_.*$/}/{@:0be r/_.*$/}.gtf' 'work/stringtie_assembled_ArME14/{@:0be r/_.*$/}.bam'" output/STAR_ArME14/[^Tween]*.cram \
| code/slurm_scripts/bin/sbatch_jobarray.sh \
  --cpus-per-task 2 \
  --batch-dry-run \
  --batch-module miniconda3/latest \
  --batch-condaenv ${PWD}/condaenv \
  --partition work \
  --time 4:00:00 \
  --job-name stringtie \
  ${RESUME_PARAM}

echo -e "ids\tisolate\ttreatment\treplicate" > output/stringtie_assembled_ArME14_counts/samples.tsv
code/slurm_scripts/bin/pt --group "{be r/_.*$/}" "{@:0be r/_.*$/}	{@:0be r/_(L|E|iv)_\d+_.*$/}	{@:0be r/_\d+_.*$/ L/.*_/}	{@:0be r/_.*$/ L/.*_/}" output/STAR_ArME14/[^Tween]*.cram >> output/stringtie_assembled_ArME14_counts/samples.tsv

code/slurm_scripts/bin/pt --group "{be r/_.*$/}" "{@:0be r/_.*$/} output/stringtie_assembled_ArME14_counts/{@:0be r/_.*$/}/{@:0be r/_.*$/}.gtf" output/STAR_ArME14/[^Tween]*.cram > output/stringtie_assembled_ArME14_counts/sample_lst.txt

