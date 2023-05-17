#!/bin/bash
module purge all
module load compiler/2022.0.2 mpi/2021.5.1 icc/2022.0.2 hdf5/1.12.1-intel netcdf-c/4.8.1-intel netcdf-fortran/4.5.4-intel mvapich2/2.3.7-intel nco



#####################################################################################DIRECTORY SETUP#######################################################################################################################

HBL=/projects/ees/dhs-crc/mjisan/HBL_V_DEV5/INIT
DIR1=$HBL/boundary_parametric/IAN/exps/2KM/INPUT           #Input directory for parametric run
DIR2=$HBL/vortex_gen                                             #Location of scaling code
DIR3=$HBL/boundary_parametric/IAN/exps/2KM                 #Model run directory
DIR4=$HBL/boundary_model/IAN/exps/2KM/INPUT                #Input directory for model run
DIR5=$HBL/boundary_model/IAN/exps/2KM                      #Model root directory
DIR6=$HBL/vortex_gen/IAN/vortex_gen_out/OUTPUT              #Path of saved outputs from initalization run
INPUT_TRACK=track_file_florence                                  #Actual Track file                                                                        
NEW_TRACK=track_file_florence_new                                #Outpur Track file from vortex generation phase

#####################################################################################Track Preprocessing####################################################################################################################

cd $DIR2

rm $NEW_TRACK                                                                                      #Remove any pre-existing track file from the output directory

cd $DIR1

rm -rf text-*                                                                                      #Remove any pre-existing copy of splitted track file


######################################################################################################################################

mv $INPUT_TRACK-sv $INPUT_TRACK                                                                    #use saved track

cp $INPUT_TRACK $INPUT_TRACK-sv                                                                    #Save copy of actual track file

csplit -b '-%d.txt' -f text -- track_file_florence '//' '{*}'                                      #Split track file for each row. 


#for i in $(seq 16 30);do                                                                           #from track file row 56: 2018-09-13 06:00:00
                                                                                                   #from track file row 58: 2018-09-13 18:00:00
i=28
    
########################################################NOTHING NEEDS TO BE CHANGED BELOW THIS LINE###########################################################################################################################
    

    
cd $DIR1

rm track_file_florence
rm text-$i-a.txt
rm text-$i-b.txt
rm text-$i-1.txt

#CHANGE TEXT NAME

myvar=$(expr $i + 1)

cat text-$i.txt >> track_file_florence
cd $DIR2
rm track_file_florence
cp $DIR1/track_file_florence .
cd $DIR1
rm track_file_florence
cat text-$i.txt >> track_file_florence

#CHANGE HOUR AND FILE NAME

HOUR=`awk '{print $4}' track_file_florence`
echo $HOUR

if [ $HOUR -eq 0000 ]
then
sed 's/0000/0600/g' text-$i.txt >> text-$i-a.txt
elif [ $HOUR -eq 0600 ]
then 
sed 's/0600/1200/g' text-$i.txt >> text-$i-a.txt
elif [ $HOUR -eq 1200 ]
then
sed 's/1200/1800/g' text-$i.txt >> text-$i-a.txt
else

GETDATE=`awk 'NR==1 {print $3}' text-$i.txt`
GETDATE1=`awk 'NR==1 {print $3+1}' text-$i.txt`

sed 's/1800/0000/g' text-$i.txt >> text-$i-1.txt
sed "s/$GETDATE/$GETDATE1/g" text-$i-1.txt >> text-$i-a.txt

fi


cat text-$i-a.txt >> track_file_florence

#CHANGE DATE


STARTDATESTRING=`awk '{print $3}' track_file_florence`
STARTDATESTRING1=`awk '{print $3}' text-$i.txt`
STARTHOUR=`awk 'NR==1 {print substr($4,0,2)}' track_file_florence`
dat2=$(date -d "$STARTDATESTRING1" +'%Y-%m-%d'" $STARTHOUR:00:00 +00:00")
dat3=$(date -d "$STARTDATESTRING1" +'%Y-%m-%d'"$STARTHOUR")

cd $DIR3

rm namelist.wps.new
cat input.nml | sed "s/^.* storm_start_time =.*$/ storm_start_time = '$dat2' /"  >> namelist.wps.new
mv namelist.wps.new input.nml
./build_parametric.bash
sbatch run_parametric.bash

sleep 40

cd $DIR4
cp $DIR3/OUTPUT/boundary_parametric.nc .
cd $DIR5
rm namelist.wps.new
cat input.nml1 | sed "s/^.* diag_file =.*$/ diag_file = 'OUTPUT\/boundary_model_$dat3.nc' /"  > namelist.wps.new
mv namelist.wps.new input.nml

./build_model.bash

sbatch  run_model.bash

sleep 1400

sbatch run_combine.bash

sleep 60


cd $DIR1

sed -i -e '2d' track_file_florence
cat text-$myvar.txt >> track_file_florence

cd $DIR2

rm -rf track_file_florence
cp $DIR1/track_file_florence .

ln -sf $DIR6/boundary_model_$dat3.nc boundary_model.nc

ncl scaling_1st_ian.ncl

sed -i -e 's/-2147483647/0000/g' $NEW_TRACK


#done

/usr/bin/sed -i -e '0~2d' $NEW_TRACK

cd $DIR1

rm text-*   #cleanup all the splitted text file


