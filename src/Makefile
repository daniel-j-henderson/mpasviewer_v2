all: driver.f90 utils.o file_manip.o mesh_rotate.o params.o
	$(FC) $(FFLAGS)  driver.f90 $(INCLUDES) utils.o file_manip.o mesh_rotate.o params.o $(LIBS) -o ../mpasviewer 

file_manip.o: file_manip.f90 params.o
	$(FC) $(FFLAGS) -c file_manip.f90 $(INCLUDES) 

params.o: params.f90	
	$(FC) $(FFLAGS) -c params.f90

utils.o: utils.f90 file_manip.o mesh_rotate.o params.o
	$(FC) $(FFLAGS) -c utils.f90 $(INCLUDES)

mesh_rotate.o: mesh_rotate.f90 params.o
	$(FC) $(FFLAGS) -c mesh_rotate.f90 $(INCLUDES)

clean: 
	rm *.o *.mod mpasviewer
