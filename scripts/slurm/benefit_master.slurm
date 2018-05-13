#!/bin/bash
#
# Set compute node type.
#SBATCH -p physical
#
# Set required nodes.
#SBATCH --nodes=1
#
# Set up nodes for job.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#
# Set memory limits (MB).  
#SBATCH --mem=245760
#
# Set walltime. Set in Days-Hours:Minutes:Seconds format
#SBATCH --time=7-00:00:00
#
# Set up job name. 
#SBATCH --job-name=SMP_$1
#
# Set up output file name.
#SBATCH --output=output/slurm_outputs/SMP_$1.out
#
# Set email notifications
#SBATCH --mail-user=davidpw@student.unimelb.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT_90

taxon=$1

module purge

export R_LIB_USER="~/R/lib/"

module load PROJ
module load GDAL
module load R/3.5.0-GCC-4.9.2

cd ~/../../data/cephfs/punim0440

Rscript --vanilla scripts/benefits/Generate_benefit_rasters.R $taxon