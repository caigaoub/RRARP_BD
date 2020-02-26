#!/bin/bash
#SBATCH --clusters=mae
#SBATCH --partition=isecc
#SBATCH --qos=isecc
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=125000
#SBATCH --job-name="MSTB_EXP_WEI"
#SBATCH --array=0-299
#SBATCH --output=../data/console/console_%A_%a.out
#SBATCH --error=../data/console/console_%A_%a.err
#SBATCH --mail-user=ningjiwe@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --exclude=cpn-p26-[07-12]
##SBATCH --constraint=CPU-E5645
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_ARRAYID="$SLURM_ARRAYID
echo "SLURM_ARRAY_JOB_ID"=$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID"=$SLURM_ARRAY_TASK_ID
echo "working directory = "$SLURM_SUBMIT_DIR
echo "SLURM_NTASKS_PER_CORE = "$SLURM_NTASKS_PER_CORE

NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo "NPROCS="$NPROCS

module load python/anaconda-4.3.1
module load gurobi/7.0.2
ulimit -s unlimited

echo ""
echo "--> BEGINNING"
echo ""

./main ./configs/config_algo0 local_config${SLURM_ARRAY_TASK_ID}

echo ""
echo "--> ALLDONE"

