# Makefile for the "mppnccombine" program
#
# Written by Hans Vahlenkamp (Hans.Vahlenkamp)
# Geophysical Fluid Dynamics Laboratory / NOAA
# Princeton Forrestal Campus
# Last updated: 05/15/08

SOURCES = mppnccombine.c
OBJECTS = mppnccombine.o
CC = cc
CFLAGS = -I${NETCDF}/include
LFLAGS = 
LIBS = -L${NETCDF}/lib -lnetcdf -L/usr/share/Modules/software/RHEL-6.5/hdf5/1.8.11_intel-14.0.3/lib -lhdf5_hl -lhdf5 -lz -limf

mppnccombine: $(OBJECTS)
	$(CC) $(CFLAGS) -o mppnccombine -s $(OBJECTS) $(LFLAGS) $(LIBS)
