#!/bin/bash

dir="/projects/ees/dhs-crc/mjisan/HBL_V_DEV3/INIT/vortex_gen/FLORENCE/vortex_gen_out/OUTPUT/"
#dir1="/data/ginislab/mjisan/HBL_V_DEV3/INIT2/vortex_gen/IRMA/vortex_gen_out/"
#dir1="/projects/ees/dhs-crc/mjisan/HBL_V_DEV3/INIT2/vortex_gen/FLORENCE/vortex_gen_out_2nd_itr_new/OUTPUT/"
dir1="/projects/ees/dhs-crc/mjisan/HBL_V_DEV3/INIT2/vortex_gen/FLORENCE/vortex_gen_out/OUTPUT/"
file="boundary_model_2018-09-1300.nc"


for file in `cd ${dir};ls -1 ${file}` ;do
   echo $file
NAME=`echo $file`

NAME2=`echo "$NAME" | awk -F'[_.]' '{print $3}'`

echo $NAME2
#echo $NAME

cd /projects/ees/dhs-crc/mjisan/HBL_V_DEV3/INIT2/vortex_gen

rm temp.ncl
rm temp2.ncl
rm vg_plot_cross.ncl

ln -sf $dir/$file boundary_model.nc



ncl create_nc_irma.ncl



ln -sf $dir1/$file boundary_model.nc



ncl create_nc_irma25.ncl
cp matched25.nc matched25.nc


sed 's/res@tiMainString="tdate"/res@tiMainString="'"$NAME2"'"/g' vg_6lines1.ncl >> temp.ncl
sed 's/"filename"/"'"$NAME2"'"/g' temp.ncl >> temp2.ncl

mv temp2.ncl vg_plot_cross.ncl
ncl vg_plot_cross.ncl
#cat vg_6lines1.ncl | sed "s/^.* res@tiMainString =.*$/ res@tiMainString = '$NAME2.nc' /"  > namelist.wps.new
#ncl 'tmain="$NAME2"' 'date="1012"' vg_6lines.ncl

done
