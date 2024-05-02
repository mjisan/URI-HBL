#!/bin/bash
#SBATCH --constraint=hatteras
#SBATCH -p lowpri
#SBATCH -t 120:00:00
#SBATCH -w largemem-5-2
#SBATCH -J job_comb
#SBATCH -e parametric.%j.err
#SBATCH -o parametric.%j.out
#SBATCH --mem=250GB


# load nco module
  module load nco

# set name of experiment and resolution
  ROOT=/projects/ees/dhs-crc/mjisan/HBL_V_DEV5/INIT3/boundary_home
  PROJECTS='/projects/ees/dhs-crc/mjisan/HBL_V_DEV5/INIT3/vortex_gen/'
  EXP=MICHAEL
  RES=LANDFALL3

# set path for mppnccombine utility
  mppnccombine=$ROOT/LIBS/mppnccombine/mppnccombine


# combine output diagnostic spatial files
  cd OUTPUT

  diag_files=`ls *.nc.0000`

  diag_list=$(echo $diag_files | rev | cut -c 6- | rev)

  echo ${diag_list}

# make sure $PROJECTS directory exists
  if [ ! -d "$PROJECTS/${EXP}/${RES}/OUTPUT" ]; then
    mkdir -p $PROJECTS/${EXP}/${RES}/OUTPUT
  fi

# move diagnostic spatial files to $PROJECTS/${EXP}/${RES}/OUTPUT
                 mv ${diag_list}.* $PROJECTS/${EXP}/${RES}/OUTPUT

# move diagnostic temporal files to $PROJECTS/${EXP}/${RES}/OUTPUT
         mv diagnostics_temporal.nc $PROJECTS/${EXP}/${RES}/OUTPUT

  cd $PROJECTS/${EXP}/${RES}/OUTPUT

  rm ${diag_list}

# combine diagnostic spatial files 
#  $mppnccombine -v -64 ${diag_list} ${diag_list}.*
  $mppnccombine -v -n4 ${diag_list} ${diag_list}.*

# combine diagnostic spatial and temporal files 

#  ncks -A ${diag_list} diagnostics_temporal.nc 
#  ncks -A diagnostics_temporal.nc ${diag_list}
#  mv diagnostics_temporal.nc ${diag_list}

# remove diagnostic spatial files 
#  rm ${diag_list}.* 

exit

