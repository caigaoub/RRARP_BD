#!/bin/sh
##SBATCH --account=pi-josewalt
#SBATCH  --clusters=faculty
#SBATCH  --partition=isecc
#SBATCH  --qos=isecc
#SBATCH --time=24:00:00
#SBATCH  --nodes=1
#SBATCH  --ntasks-per-node=12
#SBATCH --mem=120000
#SBATCH  --array=301-450
#SBATCH  --job-name="Al4"
#SBATCH  --output="./ret/console/Al4-%a.out"
#SBATCH  --error="./ret/console/Al4-%a.err"
#SBATCH --mail-user=caigao@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --exclude=cpn-p26-[07-12]
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

# make

#./bin/main 2 1 8 1 ./ret/configs/config_${SLURM_ARRAY_TASK_ID}

# ./bin/main 3 1 8 1 ./ret/configs/config_${SLURM_ARRAY_TASK_ID}

 ./bin/main 4 1 12 1 ./ret/configs/config_${SLURM_ARRAY_TASK_ID}

echo "===>> All Done!"

