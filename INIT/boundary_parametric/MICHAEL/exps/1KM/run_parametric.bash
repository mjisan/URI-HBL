#!/bin/bash 

#SBATCH --constraint=hatteras
#SBATCH -p lowpri
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH -w largemem-5-2
#SBATCH --mem=250GB
#SBATCH -t 12:02:30
#SBATCH -J job_param
#SBATCH -e boundary_parametric.%j.err
#SBATCH -o boundary_parametric.%j.out


module purge
module load compiler/2022.0.2 mpi/2021.5.1 icc/2022.0.2 hdf5/1.12.1-intel netcdf-c/4.8.1-intel netcdf-fortran/4.5.4-intel mvapich2/2.3.7-intel nco


# We need this because the system is a bit messed up (for now)
export SLURM_MPI_TYPE=pmi2

# Add our new modules:
hostname > "CONTROL.TXT"
echo $SLURM_JOB_NODELIST >> "CONTROL.TXT"

#rm -f OUTPUT/*

time srun EXEC/boundary_parametric.exe

rm -f CONTROL.TXT

exit
