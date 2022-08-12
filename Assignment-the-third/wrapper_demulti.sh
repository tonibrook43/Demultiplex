#!/bin/bash
#SBATCH -A bgmp                ###
#SBATCH --partition=bgmp       ### Partition
#SBATCH --job-name=Wrapper         ### Job Name
#SBATCH --time=24:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --cpus-per-task=1           ### Number of CPU cores per task, same as saying 8 cores

conda activate base 
./demultiplex.py