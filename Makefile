INCLUDES = $(wildcard *.o)

LIBS = $(shell nc-config --libs)

INCLUDES += -I$(shell nc-config --includedir)

NCLIBF = -lnetcdff
ifneq ($(wildcard $(NETCDF)/lib/libnetcdff.*), ) # CHECK FOR NETCDF4
		  LIBS += $(NCLIBF)
endif # CHECK FOR NETCDF4

all:
	@echo "*********************************************"
	@echo "  Use a target, such as ifort or gfortran"
	@echo "*********************************************"

gfortran:
	cd src; $(MAKE) FC="gfortran" \
	INCLUDES="$(INCLUDES)" \
	LIBS="$(LIBS)" \
	FFLAGS="-ffree-form -fbacktrace -g --std=legacy -ffree-line-length-none -fdefault-real-8" 

ifort:
	cd src; $(MAKE) FC="ifort" \
	INCLUDES="$(INCLUDES)" \
	LIBS="$(LIBS)" \
	FFLAGS="-traceback -g -debug all -autodouble" 
						
clean:
	cd src; $(MAKE) clean
	rm mpasviewer
