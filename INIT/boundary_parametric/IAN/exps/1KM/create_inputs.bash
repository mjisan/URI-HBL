#!/bin/bash -v

#################################################################
# generic compile script for experiments                     
#################################################################

#################################################################
# set environment
#################################################################

  module load nco

#################################################################
# create limits to domain
#################################################################

  xbeg=$(grep -m 1 lonstart GRIDS/storm_domain.h | awk '{print $6}')
  xend=$(grep -m 1 lonend   GRIDS/storm_domain.h | awk '{print $6}')

  ybeg=$(grep -m 1 latstart GRIDS/storm_domain.h | awk '{print $6}')
  yend=$(grep -m 1 latend   GRIDS/storm_domain.h | awk '{print $6}')

  echo ' xbeg = ' $xbeg 
  echo ' xend = ' $xend
  echo ' ybeg = ' $ybeg
  echo ' yend = ' $yend

  BASE=/projects/ees/dhs-crc/mjisan/HBL_V_DEV2/boundary_home/inputs

#################################################################
# create topography file
#################################################################

  topog=$( grep topog_name input.nml |  awk '{print $3}')

  echo '  topog = ' $topog

  if [ $topog = GMTED2010_7arc_sec_north_america.nc ]; then
    topog_file=$BASE/GMTED2010_7arc_sec_north_america.nc 
  fi

  if [ $topog = ETOPO_1arc_min_north_america.nc ]; then
    topog_file=$BASE/ETOPO_1arc_min_north_america.nc
  fi

  if [ $topog = ETOPO_2arc_min.nc ]; then
    topog_file=$BASE/ETOPO_2arc_min.nc
  fi

  ncks -d lon,$xbeg,$xend -d lat,$ybeg,$yend -O $topog_file INPUT/topog_storm_domain.nc 

#################################################################
# create land roughness file
#################################################################

  land_rough=$( grep land_rough_name input.nml |  awk '{print $3}')

  echo '  land_rough = ' $land_rough

  if [ $land_rough = land_roughness_modis_0.5km_usda_tbl.nc ]; then
    land_rough_file=$BASE/land_roughness_modis_0.5km_usda_tbl.nc
  fi

  if [ $land_rough = land_roughness_modis_0.5km_wrf_tbl.nc ]; then
    land_rough_file=$BASE/land_roughness_modis_0.5km_wrf_tbl.nc
  fi

  ncks -d lon,$xbeg,$xend -d lat,$ybeg,$yend -O $land_rough_file INPUT/land_rough_storm_domain.nc 

  exit
