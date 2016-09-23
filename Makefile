
INCLUDES = $(wildcard *.o)

ifneq "$(NETCDF)" ""
        INCLUDES += -I$(NETCDF)/include
        LIBS += -L$(NETCDF)/lib
        NCLIB = -lnetcdf
        NCLIBF = -lnetcdff
        ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
                LIBS += $(NCLIBF)
        endif # CHECK FOR NETCDF4
        LIBS += $(NCLIB)
endif


gfortran:
		( $(MAKE) all \
        "FC = gfortran" \
        "FFLAGS = -ffree-form --std=legacy -g -fbacktrace -ffree-line-length-none -fdefault-real-8" \
		"LDFLAGS = ")	

ifort:
		( $(MAKE) all \
        "FC = ifort" \
        "FFLAGS = -autodouble" \
		"LDFLAGS = ")		
		


all: mpasviewer.o

mpasviewer.o: driver.f90 utils.o file_manip.o params.o
	$(FC) $(FFLAGS)  driver.f90 $(INCLUDES) utils.o file_manip.o params.o -o mpasviewer 

file_manip.o: file_manip.f90 params.o
	$(FC) $(FFLAGS) -c file_manip.f90 $(INCLUDES) params.o

params.o: params.f90
	$(FC) $(FFLAGS) -c params.f90

utils.o: utils.f90 file_manip.o params.o
	$(FC) $(FFLAGS) -c utils.f90 file_manip.o params.o $(INCLUDES)

clean:
	rm *.o mpasviewer *.mod 
