function [ lat, long ] = lat_long( lat_cent, long_cent, nlat, nlong, stepy, stepx )
%Generates the lat and long at that center point for that number of steps
% Catrina N
% 11-12-2017

%%% Test Case for debuging
% lat_cent = 45.6;
% long_cent = 67.7;
% 
% %number of grid points
% nlat = 301;
% nlong = 301;
% 
% %resolution
% res = 1/24;
lat_cent = abs(lat_cent);
long_cent = abs(long_cent);

rearth = 6371.e3;
d2r = pi/180;

res_lat = stepy/(d2r*rearth);
res_long = stepx/(d2r*rearth*cos(lat_cent*d2r));

%make a vector of the positive and negative steps around zero
latgrid_center = (-floor(nlat/2)):1:(floor(nlat/2)-1);
%multiple the steps by the step size 
latgrid_center = latgrid_center *  res_lat ;
%Now add in what the center lat is
lat = latgrid_center + lat_cent;
%Remove the negatives
lat = abs(lat);
%lat = lat*(-1);

longgrid_center = (-floor(nlong/2)):1:(floor(nlong/2)-1);
longgrid_center = longgrid_center*-1;
longgrid_center = longgrid_center *  res_long ;
long = longgrid_center + long_cent;
long = abs(long);
long = long*(-1);
end

