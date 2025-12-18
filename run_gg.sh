#!/bin/bash
#SBATCH --job-name=run_GOgetter
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=100GB
#SBATCH --time=48:00:00
#SBATCH --output=run_gg%j.out
#SBATCH --error=run_gg%j.err
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL


module load  BLAST+/2.13.0-gompi-2022a
module load Miniforge3

source activate /home/ely67071/.conda/envs/go_getter

input_dir=$1

#makeblastdb -in "${input_dir}/Araport11_pep_20250411_representative_gene_model" -dbtype prot



bash GOgetter.sh -i ha412_transcript_b3.fasta -D "${input_dir}/ATH_GO_GOSLIM.txt" -d "${input_dir}/Araport11_pep_20250411_representative_gene_model"
