#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00
#SBATCH --job-name=hybrid
cd $SLURM_SUBMIT_DIR

make

# set variable so task placement works as expected
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

mpirun -np $SLURM_NTASKS -genv \
 OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK -genv I_MPI_PIN_DOMAIN=omp \
 ./BSIM