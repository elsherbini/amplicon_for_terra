#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH --partition=medium
#SBATCH --job-name=submit_jobs

mkdir -p ./slurm

snakemake --use-conda --conda-prefix=~/snakemake_conda_prefix --cluster-config cluster.yaml \
  --cluster "sbatch --output=slurm/{rulename}.{jobid} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem} --time={cluster.time}" \
  --jobname {rulename}.{jobid} --jobs 500 --group-components filter_trim_se=25 dereplicate_fastq=10 sample_inference=10 --keep-going --rerun-incomplete

