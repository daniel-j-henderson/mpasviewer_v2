program main 

   !use inputprocessing
   !use outputhandler
   use utils
   use mpas_file_manip
   use params
   use netcdf
   !use mesh_rotate
   implicit none
   
   
   character(len=NF90_MAX_NAME), dimension(20) :: Variables = ' ', cmdVariables = ' '
   character(len=100) :: meshInfoFile = ' ', meshDataFile = ' ', outputfile = 'latlon.output.nc'
   character(len=600) :: commandarg
   character(len=NF90_MAX_NAME) :: iVar = ' '
   logical :: bothFiles = .false., infoOnly = .false., dataOnly = .false., rotated = .false., file_present
   integer :: ierr, iFile, i, j, k, l, n, nx, ny, nFiles, nargs
   character (len = NF90_MAX_NAME), dimension(:), allocatable :: varstemp
   real(kind=RKIND) :: original_latitude_degrees = 0.0, original_longitude_degrees = 0.0, new_latitude_degrees = 0.0, new_longitude_degrees = 0.0, birdseye_rotation_counter_clockwise_degrees = 0.0
   real(kind=RKIND) :: lat_start, lon_start, lat_end, lon_end

   character(len=StrKIND) :: fname, static_file = ''
   character(len=StrKIND), dimension(:), allocatable :: files
   type(ncfile) :: ncin, ncout
   type(interpgrid) :: grid 

   
   namelist /interpolator_settings/ Variables, meshInfoFile, meshDataFile, outputfile, nx, ny, lat_start, lat_end, lon_start, lon_end
!   namelist /rotate_settings/ original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees

   inquire(file='namelist.input', exist=file_present)
   if (file_present) then
      open(41,file='namelist.input') 
      read(41,interpolator_settings)
   !   read(41,rotate_settings)
      close(41)
   end if
    
   grid%nx = nx
   grid%ny = ny
   grid%lat_start = lat_start
   grid%lat_end = lat_end
   grid%lon_start = grid%lon_start
   grid%lon_end = lon_end 
    
   call get_command(commandarg, l, ierr)
   
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

   n = 0
   do i=1, size(Variables)
      if (len(trim(Variables(i))) > 0) n = n+1
   end do

   nargs = command_argument_count()
   do i=1, nargs
      call get_command_argument(i, fname)
      if (trim(fname) == '-v') then
         j = i+1
         call get_command_argument(j, fname)
         do while(fname(1:1) .ne. '-' .and. j <= nargs)
            n = n+1
            j = j+1
            call get_command_argument(j, fname)
         end do
         exit
      end if
   end do
        
   allocate(varstemp(n))

   do i=1, nargs
      call get_command_argument(i, fname)
      if (trim(fname) == '-v') then
         j = i+1
         k = 1
         call get_command_argument(j, fname)
         do while(fname(1:1) .ne. '-' .and. j <= nargs)
            varstemp(k) = fname
            j = j+1
            k = k+1
            call get_command_argument(j, fname)
         end do
         exit
      end if
   end do

   j = 1
   k = k+1
   do while (k <= n)
      varstemp(k) = Variables(j)
      j = j+1
      k = k+1
   end do
    
   
   i=1
   nFiles = 0   
   call get_command_argument(1, fname)
   do while(fname(1:1) .ne. '-')
      if (index(fname, '.nc', back=.true.) .ne. len(trim(fname)) - 2) then
         write (0,*) "All input files must be netcdf files"
         stop
      end if
      nFiles = nFiles + 1
      i = i+1
      call get_command_argument(i, fname)
   end do

   if (nFiles < 1) then
      write (0,*) 'You neet to provide input netcdf files as the first command line arguments'
      stop
   end if

   allocate(files(nFiles))
   
   do i=1, nFiles
      call get_command_argument(i, files(i))
      if (len(trim(static_file)) == 0) then
         if (file_contains_elem(files(i), VAR, 'cellsOnCell')) then
            static_file = files(i)
         end if
      end if
   end do

    
    ! Determine if we have one or both files, and set the filenames accordingly

   ncin%filename = static_file
   if (len(trim(ncin%filename)) == 0) then
      write (0,*) "At least one of the files provided must contain the all the static mesh information fields."
      stop
   end if

   
   if(n == 0) then
      write (0,*) "Must provide the desired variables in the namelist or command line."
      stop
   end if

   
   ! edit the open file function to perform the setup() also 
   call open_mpas_file(ncin, 'NF90_NOWRITE')
   write (0,*) "Opened "//trim(ncin%filename)
   write (0,*) "Calling create_grid"
   call create_grid(grid, ncin)

   write (0,*) "Closing "//trim(ncin%filename)
   call close_mpas_file(ncin)

   do iFile = 1, nFiles
      ncin%filename = files(iFile)

      call open_mpas_file(ncin, 'NF90_NOWRITE')
      write (0,*) "Opened "//trim(ncin%filename)

      ncout%filename = ''
      ncout%filename = trim(outputfile)//trim(ncin%filename)
      call open_mpas_file(ncout, 'CREATE')
      write (0,*) "Created "//trim(ncout%filename)

      !maybe check to make sure they are comparable.

      ! Handle rotation

      ! create type called grid with the elements map, start-lat, end-lat, etc.

      write (0,*) "Calling create_output_from_grid"
      call create_output_from_grid(ncin, ncout, grid)
      

      write (0,*) "Calling define_variables_io"
      call define_variables_io(ncin, ncout, varstemp)

      ierr = nf90_enddef(ncout%ncid)
      if (ierr /= NF90_NOERR) then
         write (0,*) "********************************************"
         write (0,*) "  Error ending define mode"
         write (0,*) "********************************************"
         stop
      end if

      do i=1, n
         if (contains_elem(ncin, VAR, varstemp(i))) then
            call copy_data(ncin, ncout, varstemp(i), grid)
         end if
      end do
                 
      write (0,*) "Closing "//trim(ncin%filename)
      call close_mpas_file(ncin)
         
      write (0,*) "Closing "//trim(ncout%filename)
      call close_mpas_file(ncout)
      write (0,*) "End loop"   
   end do

   write (0,*) "End Program Main"
end program main
