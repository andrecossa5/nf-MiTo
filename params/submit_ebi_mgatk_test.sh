#!/bin/bash
#SBATCH --job-name=nf-test
#SBATCH --mem=1GB
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

source ~/.bashrc
PATH_PIPELINE="/homes/cossa/nf-MiTo"     # git clone of nf-MiTo, on branch v2.0.0 (git pull before submitting)
cd /nfs/research/coorenslab/cossa/nf_mito_test

nextflow run $PATH_PIPELINE/main.nf -params-file $PATH_PIPELINE/params/params_ebi_mgatk_test.json -profile slurm,conda -c $PATH_PIPELINE/config/user.config -with-tower
