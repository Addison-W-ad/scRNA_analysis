#!/bin/bash

#SBATCH --job-name=cr-gex
#SBATCH --array=0-3
#SBATCH --partition=512GB
#SBATCH --time=12:00:00
#SBATCH --output=/Path/logs/cr-gex/2025-05-19/cr-gex_%A_%a.log
#SBATCH --mail-user=s233754@utsouthwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=450G

# This script runs cellranger count on 1 sample, 1 node. 
# since the day0 sample has running, we start from day1.
# Submit from the project directory.
# %A = The master SLURM job array ID (JobID)
# %a = The array index for each sub-task

# Read in sample and input info from the config file
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID \
    'NR>1 && $1 == ArrayTaskID {print $2}' \
    /Path/slurm_array_config_2025-05-19.tsv)

input=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID \
    'NR>1 && $1 == ArrayTaskID {print $3}' \
    /Path/slurm_array_config_2025-05-19.tsv)

REF="/shared/genome_reference/refdata-gex-mm10-2020-A" # Location of reference genome downloaded from 10X
TOKEN="10xcloud_token.json"

cd cellranger_counts

module load cellranger/9.0.0 

srun cellranger count \
--id=${sample} \
--libraries=${input} \
--transcriptome=${REF} \
--tenx-cloud-token-path=${TOKEN} \
--cell-annotation-model=auto \
--create-bam=true \
--chemistry=threeprime \
--localcores=32 \
--localmem=450

echo "Cell Ranger count completed for sample ${sample}" 