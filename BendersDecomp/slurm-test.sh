#!/bin/sh
##SBATCH --account=pi-josewalt
#SBATCH --clusters=faculty
#SBATCH --partition=isecc
#SBATCH --qos=isecc
##SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
##SBATCH --mem=120000
#SBATCH --job-name="ccr_test"
#SBATCH --output="ccr_test.out"
#SBATCH --error="ccr_test.err"
##SBATCH --mail-user=caigao@buffalo.edu
##SBATCH --mail-type=ALL
##SBATCH --requeue

##echo "SLURM_JOB_ID="$SLURM_JOB_ID
##echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
##echo "SLURM_NNODES"=$SLURM_NNODES
##echo "SLURMTMPDIR="$SLURMTMPDIR
##echo "SLURM_ARRAYID="$SLURM_ARRAYID
##echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
##echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
##echo "working directory = "$SLURM_SUBMIT_DIR

NPROCS=`srun --nodes=${SLURM_NNODES}$ bash -c 'hostname' |wc -l`
echo "NPROCS="$NPROCS

module load gurobi/9.0.0

make

./bin/main  8 ../InstanceGenerator/ret/n_8_e_1.txt 4

echo "All Done!"

