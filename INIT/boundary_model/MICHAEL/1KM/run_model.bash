#!/bin/bash
#SBATCH --constraint=hatteras
#SBATCH -p lowpri
##SBATCH -n 480
#SBATCH --nodes=12
#SBATCH --tasks-per-node=40
#SBATCH -t 72:00:00
#SBATCH --exclude=compute-5-1
#SBATCH --mem-per-cpu=4G
#SBATCH -J job_model
#SBATCH -e boundary_model.%j.err
#SBATCH -o boundary_model.%j.out



module purge
module load compiler/2022.0.2 mpi/2021.5.1 icc/2022.0.2 hdf5/1.12.1-intel netcdf-c/4.8.1-intel netcdf-fortran/4.5.4-intel mvapich2/2.3.7-intel nco

# we need this because the system is a bit messed up (for now)
export SLURM_MPI_TYPE=pmi2

# add our new modules:
hostname > "CONTROL.TXT"
echo $SLURM_JOB_NODELIST >> "CONTROL.TXT"


rm -f OUTPUT/*
echo $PWD


time srun EXEC/boundary_model.exe

rm -f CONTROL.TXT


exit
