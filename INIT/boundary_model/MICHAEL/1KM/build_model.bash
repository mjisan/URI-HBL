#!/bin/bash -v

#################################################################
# generic compile script for experiments                     
#################################################################

#################################################################
# set environment
#################################################################

  alias make="make -j 2"
  list_paths="/projects/ees/dhs-crc/mjisan/HBL_V_DEV5/INIT3/LIBS/fms/list_paths"
  mkmf="/projects/ees/dhs-crc/mjisan/HBL_V_DEV5/INIT3/LIBS/fms/mkmf"

#################################################################
# set paths
#################################################################


  root=$(echo $PWD | cut -d'/' -f-11)
 
  base=/projects/ees/dhs-crc/mjisan/HBL_V_DEV5/INIT3/boundary_home
  src=$base/source

  updates=$root/updates
  grids=$PWD/GRIDS

  cppDefs=

  template=$base/templates/intel_default.mk
  executable=$PWD/EXEC/boundary_model.exe

#################################################################
# create Makefile
#################################################################

  cd $root/objs
  rm $root/objs/*

  $list_paths -o $root/objs/pathnames $src/model $src/shared \
               $grids

  $mkmf -a $root/objs -m Makefile -t $template -p $executable \
        $updates -c $cppDefs $root/objs/pathnames

#################################################################
# call the main Makefile
#################################################################

  make -f Makefile

  exit
