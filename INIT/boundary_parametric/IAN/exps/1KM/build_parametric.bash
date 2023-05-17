#!/bin/bash -v

#################################################################
# generic compile script for experiments                     
#################################################################

#################################################################
# set environment
#################################################################

  ulimit -s unlimited
  ulimit -a

  alias make="make -j 2"

#################################################################
# set paths
#################################################################


  root=$(echo $PWD | cut -d'/' -f-11)
 
  base=/projects/ees/dhs-crc/mjisan/HBL_V_DEV5/INIT/boundary_home
  src=$base/source
  list_paths=$base/LIBS/fms/list_paths
  mkmf=$base/LIBS/fms/mkmf

  updates=$root/updates
  grids=$PWD/GRIDS

  cppDefs=

  template=$base/templates/intel_default.mk
  executable=$PWD/EXEC/boundary_parametric.exe 

#################################################################
# create Makefile
#################################################################

  cd $root/objs
  rm $root/objs/*

  $list_paths -o $root/objs/pathnames $src/parametric $src/shared $grids

  $mkmf -a $root/objs -m Makefile -t $template -p $executable \
        $updates -c $cppDefs $root/objs/pathnames

#################################################################
# call the main Makefile
#################################################################

  make -f Makefile

  exit
