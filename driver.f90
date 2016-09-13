program main 

   use inputprocessing
   use outputhandler
   use params
   use netcdf
   use mesh_rotate
   implicit none
   
   
   character(len=NF90_MAX_NAME), dimension(20) :: Variables = ' ', cmdVariables = ' '
   character(len=100) :: meshInfoFile = ' ', meshDataFile = ' ', outputfile = 'latlon.output.nc'
   character(len=600) :: commandarg
   character(len=NF90_MAX_NAME) :: var = ' '
   logical :: bothFiles = .false., infoOnly = .false., dataOnly = .false., rotated = .false., file_present
   integer :: l, n
   character (len = NF90_MAX_NAME), dimension(:), allocatable :: varstemp
   real :: original_latitude_degrees = 0.0, original_longitude_degrees = 0.0, new_latitude_degrees = 0.0, new_longitude_degrees = 0.0, birdseye_rotation_counter_clockwise_degrees = 0.0

   type(ncfile) :: ncin1, ncin2, ncout

   
   namelist /interpolator_settings/ Variables, meshInfoFile, meshDataFile, outputfile, gridW, gridH, gridHmin, gridHmax, gridWmin, gridWmax
   namelist /rotate_settings/ original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees

   inquire(file='namelist.input', exist=file_present)
   if (file_present) then
      open(41,file='namelist.input') 
      read(41,interpolator_settings)
      read(41,rotate_settings)
      close(41)
   end if
    
    
    
   call get_command(commandarg, l, ierr)
   
   if(index(commandarg, '-v') > 0) then
      i = index(commandarg, '-v')
      i = i+2
      k=1
      do while(commandarg((i+1):(i+1)) /= '-' .and. commandarg((i+1):(i+2)) /= ' ')
         j=1
         do while(commandarg((i+j):(i+j)) /= ' ')
            var(j:j) = commandarg((i+j):(i+j))
            j=j+1
         end do
         cmdVariables(k) = var
         k = k+1
         var = ' '
         i = i+j
      end do
   end if
   
   if (index(commandarg, '-i') > 0) then
      meshInfoFile = ' '
      i = index(commandarg, '-i')
      i = i+3
      k = 1
      do while (commandarg(i:i) /= ' ')
         meshInfoFile(k:k) = commandarg(i:i)
         k = k+1
         i = i+1
      end do
   end if
   
   if (index(commandarg, '-d') > 0) then
      meshDataFile = ' '
      i = index(commandarg, '-d')
      i = i+3
      k = 1
      do while (commandarg(i:i) /= ' ')
         meshDataFile(k:k) = commandarg(i:i)
         k = k+1
         i = i+1
      end do
   end if
   
   if (index(commandarg, '-o') > 0) then
      outputfile = ' '
      i = index(commandarg, '-o')
      i = i+3
      k = 1
      do while (commandarg(i:i) /= ' ')
         outputfile(k:k) = commandarg(i:i)
         k = k+1
         i = i+1
      end do
   end if
      

    
    ! Determine if we have one or both files, and set the filenames accordingly
   if (len_trim(meshInfoFile) > 0) then
      if (len_trim(meshDataFile) > 0) then
         bothFiles = .true.
         ncin1%filename = trim(meshInfoFile)
         ncin2%filename = trim(meshDataFile)
      else
         infoOnly = .true.
         ncin1%filename = trim(meshInfoFile)
         ncin2%filename = trim(meshInfoFile)
      end if
   else if (len_trim(meshDataFile) > 0) then
      dataOnly = .true.
      ncin1%filename = trim(meshDataFile)
      ncin2%filename = trim(meshDataFile)
   else
      write (0,*) 'You need to provide at least one input file either in the namelist'
      write(0,*) 'variables "meshInfoFile" or "meshDataFile," or you can use the'
      write(0,*) 'command line arguments -i or -d, respectively'
      write(0,*) 'Program Terminated.'
      stop
   end if

   ncout%filename = trim(outputfile)
   
   
   k=size_of(cmdVariables)
   j=size_of(Variables)
   n = j + k
   if (n == 0) then
      write (0,*) 'You need to provide at least one desired variable either in the namelist'
      write(0,*) 'variable "Variables" (e.g. ''var1'', ''var2''...), or you can use the'
      write(0,*) 'command line argument -v'
      write(0,*) 'Program Terminated.'
      stop
   end if
   
   
   !Extract the desired mesh variables from the Variables array
   
    allocate(varstemp(n))
    do i=1,j
      varstemp(i) = trim(Variables(i))
    end do
    do i=j,n
      if (i < n) then
         varstemp(i+1) = trim(cmdVariables(i+1-j))
      end if
    end do
   
   ! edit the open file function to perform the setup() also 
   call open_mpas_file(ncin1%filename, 'NF90_NOWRITE')
   call open_mpas_file(ncin2%filename, 'NF90_NOWRITE')

   !maybe check to make sure they are comparable.

   ! Handle rotation

   ! create type called grid with the elements map, start-lat, end-lat, etc.
   call create_grid_map(grid, rotated, ncin1)
   
   do i=1, n
      if (contains_elem(ncin2%filename, 'VAR', varstemp(i))) then
         call interp_variable_defmode(ncin2, ncout, varstemp(i), grid)
      else if (contains_elem(ncin1%filename, 'VAR', varstemp(i))) then
         call interp_variable_defmode(ncin1, ncout, varstemp(i), grid) 
      else
         write (0,*) "Variable "//trim(varstemp(i))//" does not exist, skipping it."
      end if
   end do

   ierr = nf90_enddef(ncout%ncid)
   if (ierr /= NF90_NOERR) then
      write (0,*) "********************************************"
      write (0,*) "  Error ending define mode"
      write (0,*) "********************************************"
      stop
   end if

   do i=1, n
      if (contains_elem(ncin2%filename, 'VAR', varstemp(i))) then
         call interp_variable_datamode(ncin2, ncout, varstemp(i), grid)
      else if (contains_elem(ncin1%filename, 'VAR', varstemp(i))) then
         call interp_variable_datamode(ncin1, ncout, varstemp(i), grid) 
      end if
   end do
              
   call close_mpas_file(ncin1)
   call close_mpas_file(ncin2)
   call close_mpas_file(ncout)

end program main










    
    
   
   write(*,*) 'Opening Input'
    call open_input(filename, filename2)
   write(*,*) 'Running Setup'
    call setup()
    

    ! if (needs rotated), then call rotate(ncid, MeshX, MeshY, MeshZ, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, &
    ! original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees)
    
    if (needs_rotated(original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees)) then
      call rotate(ncid, MeshX, MeshY, MeshZ, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, &
                    original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees)
        rotated = .true.
    end if
    
    write(*,*) 'Checking for existence of variables, throwing out any for which there is no input data.'
    call check_existence(varstemp)
   write(*,*) 'Interpolating data for variables:', desiredMeshVars
    write(*,*) 'Creating Grid Map'
    call create_grid_map(grid, rotated)
    write(*,*) 'Creating Output File'
    call create_output_file(newFilename)
    
    
   write(*,*) 'Determining Dimensionality'
    call dimensionality(desiredMeshVars)
    write(*,*) 'Putting data into file'
   !For each desired mesh variable, based on its dimension call the appropriate routine
   !to put the data in the output file
   do i = 1, nMeshVars
      if (nDims(i) == 1) then
         call put_data1(meshVarIDs(i), meshDimIDs(i,1), gridVarIDs(i))
      else if (nDims(i) == 2) then
         call put_data2(meshVarIDs(i), meshDimIDs(i,:), gridVarIDs(i))
      else if (nDims(i) == 3) then
         call put_data3(meshVarIDs(i), meshDimIDs(i,:), gridVarIDs(i))
      end if
    end do
    
    call put_latlons(grid(:,:,1), grid(:,:,2))
    
   write(*,*) 'Cleaning Up'
    call clean_up()
    deallocate(desiredMeshVars)
    
    write(*,*) 'Successfully interpolated all your desired variables onto a lat-lon grid'
end program latlondriver
