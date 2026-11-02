#!/bin/bash
#SBATCH --job-name=clipseq_fulltest
#SBATCH --output=clipseq_fulltest_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --mail-user=hannah.18.jones@kcl.ac.uk
#SBATCH --time=08:00:00


ml Anaconda3 Singularity

source $EBROOTANACONDA3/etc/profile.d/conda.sh

conda activate /camp/home/jonesh1/.conda/envs/nf

nextflow run main.nf -c ../ref_data/crick_config --consensus_peak=false --peakcaller clippy  -profile test_full_te,singularity
