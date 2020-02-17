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
##SBATCH --mail-user=josewalt@buffalo.edu
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

echo "start computing now"
PP       = g++
CPPARGS   = -O3 -m64 -std=c++1y -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./dat/
# GRBPATH   = /opt/gurobi900/linux64
GRBPATH =  /util/academic/gurobi/gurobi900/linux64
INCGRB    = -I$(GRBPATH)/include/
CPPLIBGRB = -L$(GRBPATH)/lib/ -lgurobi_c++ -lgurobi90 $(CPPSTDLIB) -lpthread -lm
# INSTANCE	= 7 $(DATPATH)sample_n_25.txt

all:
	$(CPP) $(CPPARGS) $(SRCPATH)DataHandler.cpp \
	$(SRCPATH)PartitionScheme.cpp \
	$(SRCPATH)STEFormulation.cpp \
	$(SRCPATH)BendersCuts.cpp \
	$(SRCPATH)DualFormulation.cpp \
	$(SRCPATH)GlobalMC.cpp \
	-o $(BINPATH)main $(SRCPATH)main.cpp $(INCGRB) $(CPPLIBGRB)

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main

run:
	$(BINPATH)main

./bin/main  8 ../InstanceGenerator/ret/inst_n_10/n_10_e_1.txt 3
echo "All Done!"

