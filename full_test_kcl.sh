#!/bin/bash 
#SBATCH --job-name=clipseq_fulltest 
#SBATCH --output=clipseq_fulltest_%j.out 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem=20G
#SBATCH --mail-user=hannah.18.jones@kcl.ac.uk
#SBATCH --time=24:00:00

source /scratch/prj/ppn_rnp_networks/users/hannah.jones/software/mambaforge/bin/activate nf

nextflow run main.nf --telescope_gtf ensembl_GRCh38_rmsk_validation.gtf -c ../../../../shared/nextflow/KCL-CREATE.config -resume -profile test_full,singularity
