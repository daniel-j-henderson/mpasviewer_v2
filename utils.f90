module utils
   use netcdf
   use mpas_file_manip
   type :: interpgrid
      integer, dimension(:,:), pointer :: cell_map, edge_map, vertex_map
      real (kind=RKIND) :: lat_start, lat_end, lon_start, lon_end
      integer :: nx, ny
   end type interpgrid

   contains

   subroutine create_grid(grid, f)
      implicit none
      type(interpgrid) :: grid
      type(ncfile) :: f

      real(kind=RKIND), dimension(:), pointer :: latElem, lonElem
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell, verticesOnCell
      real(kind=RKIND) :: temp, dist, mindist, lat_pt, lon_pt, dx, dy
      integer :: v, iCell, k, i, j, e, nearest_cell

      if (grid%lat_start == grid%lat_end .or. grid%lon_start == grid%lon_end .or. grid%nx == 0 .or. grid%ny == 0) then
         write (0,*) "Error: Invalid argument ''grid'' in ''create_grid''."
         stop
      end if


      allocate(grid%cell_map(grid%nx, grid%ny))
         
      
      call get_variable_1dREAL(f, 'latCell', latElem)
      call get_variable_1dREAL(f, 'lonCell', lonElem)
      call get_variable_2dINT(f, 'cellsOnCell', cellsOnCell)
      call get_variable_1dINT(f, 'nEdgesOnCell', nEdgesOnCell)

      write (0,*) "  making cell_map"
      nearest_cell = 1
      dx = (grid%lon_end - grid%lon_start) / grid%nx
      dy = (grid%lat_end - grid%lat_start) / grid%ny
      do j=1, grid%ny
         lat_pt = grid%lat_start + (j-1) * dy + dy / 2.0
         lat_pt = lat_pt * PI / 180.0
      do i=1, grid%nx
         lon_pt = grid%lon_start + (i-1) * dx + dx / 2.0
         lon_pt = lon_pt * PI / 180.0
         
         nearest_cell = nearest_cell_path(lat_pt, lon_pt, nearest_cell, &
                                  nEdgesOnCell, cellsOnCell, latElem, lonElem)
         grid%cell_map(i, j) = nearest_cell
      end do
      end do
      
      deallocate(latElem, lonElem, cellsOnCell)
      allocate(grid%vertex_map(grid%nx, grid%ny))
         
      call get_variable_1dREAL(f, 'latVertex', latElem)
      call get_variable_1dREAL(f, 'lonVertex', lonElem)
      call get_variable_2dINT(f, 'verticesOnCell', verticesOnCell)

      write (0,*) "  making vertex_map"
      do j=1, grid%ny
         lat_pt = grid%lat_start + (j-1) * dy + dy / 2.0
      do i=1, grid%nx
         lon_pt = grid%lon_start + (i-1) * dx + dx / 2.0
         iCell = grid%cell_map(i, j)
         v = verticesOnCell(1, iCell)
         mindist = sphere_distance(latElem(v), lonElem(v), lat_pt, lon_pt, 1.0)
         do k = 2, nEdgesOnCell(iCell)  
            temp = verticesOnCell(k, iCell)
            dist = sphere_distance(latElem(temp), lonElem(temp), lat_pt, lon_pt, 1.0)
            if (dist < mindist) then
               mindist = dist
               v = temp
            end if
         end do
         grid%vertex_map(i, j) = v
      end do
      end do
      

      deallocate(latElem, lonElem, verticesOnCell)
      allocate(grid%edge_map(grid%nx, grid%ny))
      
      call get_variable_2dINT(f, 'edgesOnCell', edgesOnCell) 
      call get_variable_1dREAL(f, 'latEdge', latElem)
      call get_variable_1dREAL(f, 'lonEdge', lonElem)
      write(0,*) "  making edge_map"
      do j=1, grid%ny
         lat_pt = grid%lat_start + (j-1) * dy + dy / 2.0
      do i=1, grid%nx
         lon_pt = grid%lon_start + (i-1) * dx + dx / 2.0
         iCell = grid%cell_map(i, j)
         e = edgesOnCell(1, iCell)
         mindist = sphere_distance(latElem(v), lonElem(v), lat_pt, lon_pt, 1.0)
         do k = 2, nEdgesOnCell(iCell)  
            temp = edgesOnCell(k, iCell)
            dist = sphere_distance(latElem(temp), lonElem(temp), lat_pt, lon_pt, 1.0)
            if (dist < mindist) then
               mindist = dist
               e = temp
            end if
         end do
         grid%edge_map(i, j) = e
      end do
      end do

      deallocate(lonElem, latElem, edgesOnCell, nEdgesOnCell)

   end subroutine create_grid

   subroutine create_output_from_grid(fin, fout, grid)
      implicit none
      type(interpgrid), intent(in) :: grid
      type(ncfile), intent(in) :: fin
      type(ncfile), intent(inout) :: fout

      call add_dimension(fout, 'xDim', grid%nx)
      call add_dimension(fout, 'yDim', grid%ny)

      call copy_dimensions(fin, fout) 

   end subroutine create_output_from_grid

   subroutine define_variables_io(fin, fout, vars)
      implicit none
      type(ncfile) :: fin, fout
      character(len=*), dimension(:), intent(inout) :: vars
      integer :: i

      do i=1, size(vars)
         if (.not. fout%contains_elem(VAR, trim(vars(i))) .and. fin%contains_elem(VAR, trim(vars(i)))) then
           call copy_variable_defmode(fin, fout, vars(i))
         end if
      end do 

   end subroutine define_variables_io

   subroutine copy_data(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name
      integer :: ierr, var_id, xtype, ndims

      ierr = nf90_inq_varid(fin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(fin%ncid, var_id, xtype=xtype, ndims=ndims)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      if (xtype == NF90_INT) then
         select case(ndims)
         case(1)
            call copy_data_1dINT(fin, fout, var_name, grid)
         case(2)
            call copy_data_2dINT(fin, fout, var_name, grid)
         case(3)
            call copy_data_3dINT(fin, fout, var_name, grid)
         case default
            write (0,*) "Bad ndims in copy_data"
         end select

      else if (xtype == NF90_REAL .or. xtype == NF90_DOUBLE) then
         select case(ndims)
         case(1)
            call copy_data_1dREAL(fin, fout, var_name, grid)
         case(2)
            call copy_data_2dREAL(fin, fout, var_name, grid)
         case(3)
            call copy_data_3dREAL(fin, fout, var_name, grid)
         case default
            write (0,*) "Bad ndims in copy_data"
         end select
      
      else 
         write (0,*) "Bad xtype in copy_data"
         return
      end if

   end subroutine copy_data

   subroutine copy_data_1dINT(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map
      integer, dimension(:), pointer :: field
      integer, dimension(:,:), pointer :: newfield
      integer :: i, j

      call get_variable_1dINT(fin, var_name, field)
      
      if (size(field) == fin%nCells) then
         map => grid%cell_map
      else if (size(field) == fin%nEdges) then
         map => grid%edge_map 
      else if (size(field) == fin%nVertices) then
         map => grid%vertex_map
      else
         write (0,*) "Variable "//trim(var_name)//" not dimensioned spatially"
         return
      end if

      allocate(newfield(grid%nx, grid%ny))

      do j=1, grid%ny
      do i=1, grid%nx
         newfield(i, j) = field(map(i,j))
      end do
      end do

      call put_variable_2dINT(fout, newfield, var_name)
   end subroutine copy_data_1dINT

   subroutine copy_data_2dINT(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map
      integer, dimension(:,:), pointer :: field
      integer, dimension(:,:,:), pointer :: newfield
      integer :: i, j, t, ierr, var_id, ndims, xtype
      integer, dimension(2) :: dimids = 0, dimlens = 0
      character(len=StrKIND) :: dim_name
      logical :: has_time

      
      ierr = nf90_inq_varid(fin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(fin%ncid, var_id, xtype=xtype, ndims=ndims, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      do i=1, 2
         ierr = nf90_inquire_dimension(fin%ncid, dimids(i), name=dim_name, len=dimlens(i))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimension '//trim(dim_name)//' in '//fin%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'         
         end if
      end do

      if(dimlens(2) == fin%nCells) then
         map => grid%cell_map
         has_time = .false.
      else if (dimlens(1) == fin%ncells) then
         map => grid%cell_map
         has_time = .true.
      else if (dimlens(2) == fin%nEdges) then
         map => grid%edge_map
         has_time = .false.
      else if (dimlens(1) == fin%nEdges) then
         map => grid%edge_map
         has_time = .true.
      else if (dimlens(2) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .false.
      else if (dimlens(1) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .true.
      else
         write (0,*) "Not sure which map to use, copy data mode"
         return
      end if
       
   
      !if (dimlens(1) * dimlens(2) > MAX_CHUNK_SIZE) then
         ! Too big, must handle in slices 
      !   if (has_time) then
      !      do i=1, dimlens(2)
      if(.false.) then         
         write (0,*) "FALSE"
      else
         call get_variable_2dINT(fin, var_name, field)
         if (has_time) then
            allocate(newfield(grid%nx, grid%ny, dimlens(2)))
            do t=1, dimlens(2)
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, t) = field(map(i, j), t)
            end do
            end do   
            end do   
         else
            allocate(newfield(grid%nx, grid%ny, dimlens(1)))
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, :) = field(:, map(i, j))
            end do
            end do   
         end if

         call put_variable_3dINT(fout, newfield, var_name)
      end if

   end subroutine copy_data_2dINT

   subroutine copy_data_3dINT(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map
      integer, dimension(:,:,:), pointer :: field
      integer, dimension(:,:,:,:), pointer :: newfield
      integer :: i, j, ierr, xtype, ndims, var_id
      integer, dimension(3) :: dimids = 0, dimlens = 0
      character(len=StrKIND) :: dim_name
      logical :: has_time 

      
      ierr = nf90_inq_varid(fin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(fin%ncid, var_id, xtype=xtype, ndims=ndims, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      do i=1, 3
         ierr = nf90_inquire_dimension(fin%ncid, dimids(i), name=dim_name, len=dimlens(i))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimension '//trim(dim_name)//' in '//fin%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'         
         end if
      end do

      if(dimlens(3) == fin%nCells) then
         map => grid%cell_map
         has_time = .false.
      else if (dimlens(2) == fin%ncells) then
         map => grid%cell_map
         has_time = .true.
      else if (dimlens(3) == fin%nEdges) then
         map => grid%edge_map
         has_time = .false.
      else if (dimlens(2) == fin%nEdges) then
         map => grid%edge_map
         has_time = .true.
      else if (dimlens(3) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .false.
      else if (dimlens(2) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .true.
      else
         write (0,*) "Not sure which map to use, copy data mode"
         return
      end if
       
   
      !if (dimlens(1) * dimlens(2) > MAX_CHUNK_SIZE) then
         ! Too big, must handle in slices 
      !   if (has_time) then
      !      do i=1, dimlens(2)
      if(.false.) then         
         write (0,*) "FALSE"
      else
         call get_variable_3dINT(fin, var_name, field)
         if (has_time) then
            allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(3)))
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, :, :) = field(:, map(i, j), :)
            end do
            end do   
         else
            allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(2)))
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, :, :) = field(:, :, map(i, j))
            end do
            end do   
         end if

         call put_variable_4dINT(fout, newfield, var_name)
      end if

   end subroutine copy_data_3dINT

   
   subroutine copy_data_1dREAL(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map
      real(kind=RKIND), dimension(:), pointer :: field
      real(kind=RKIND), dimension(:,:), pointer :: newfield
      integer :: i, j

      call get_variable_1dREAL(fin, var_name, field)
      
      if (size(field) == fin%nCells) then
         map => grid%cell_map
      else if (size(field) == fin%nEdges) then
         map => grid%edge_map 
      else if (size(field) == fin%nVertices) then
         map => grid%vertex_map
      else
         write (0,*) "Variable "//trim(var_name)//" not dimensioned spatially"
         return
      end if

      allocate(newfield(grid%nx, grid%ny))

      do j=1, grid%ny
      do i=1, grid%nx
         newfield(i, j) = field(map(i,j))
      end do
      end do

      call put_variable_2dREAL(fout, newfield, var_name)
   end subroutine copy_data_1dREAL

   subroutine copy_data_2dREAL(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map
      real(kind=RKIND), dimension(:,:), pointer :: field
      real(kind=RKIND), dimension(:,:,:), pointer :: newfield
      integer :: i, j, xtype, ierr, ndims, var_id
      integer, dimension(2) :: dimids = 0, dimlens = 0
      character(len=StrKIND) :: dim_name
      logical :: has_time

      
      ierr = nf90_inq_varid(fin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(fin%ncid, var_id, xtype=xtype, ndims=ndims, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      do i=1, 2
         ierr = nf90_inquire_dimension(fin%ncid, dimids(i), name=dim_name, len=dimlens(i))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimension '//trim(dim_name)//' in '//fin%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'         
         end if
      end do

      if(dimlens(2) == fin%nCells) then
         map => grid%cell_map
         has_time = .false.
      else if (dimlens(1) == fin%ncells) then
         map => grid%cell_map
         has_time = .true.
      else if (dimlens(2) == fin%nEdges) then
         map => grid%edge_map
         has_time = .false.
      else if (dimlens(1) == fin%nEdges) then
         map => grid%edge_map
         has_time = .true.
      else if (dimlens(2) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .false.
      else if (dimlens(1) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .true.
      else
         write (0,*) "Not sure which map to use, copy data mode"
         return
      end if
       
   
      !if (dimlens(1) * dimlens(2) > MAX_CHUNK_SIZE) then
         ! Too big, must handle in slices 
      !   if (has_time) then
      !      do i=1, dimlens(2)
      if(.false.) then         
         write (0,*) "FALSE"
      else
         call get_variable_2dREAL(fin, var_name, field)
         if (has_time) then
            allocate(newfield(grid%nx, grid%ny, dimlens(2)))
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, :) = field(map(i, j), :)
            end do   
            end do   
         else
            allocate(newfield(grid%nx, grid%ny, dimlens(1)))
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, :) = field(:, map(i, j))
            end do
            end do   
         end if

         call put_variable_3dREAL(fout, newfield, var_name)
      end if

   end subroutine copy_data_2dREAL

   subroutine copy_data_3dREAL(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map
      real(kind=RKIND), dimension(:,:,:), pointer :: field
      real(kind=RKIND), dimension(:,:,:,:), pointer :: newfield
      integer :: i, j, ierr, ndims, var_id, xtype
      integer, dimension(3) :: dimids = 0, dimlens = 0
      character(len=StrKIND) :: dim_name
      logical :: has_time

      
      ierr = nf90_inq_varid(fin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(fin%ncid, var_id, xtype=xtype, ndims=ndims, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      do i=1, 3
         ierr = nf90_inquire_dimension(fin%ncid, dimids(i), name=dim_name, len=dimlens(i))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimension '//trim(dim_name)//' in '//fin%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'         
         end if
      end do

      if(dimlens(3) == fin%nCells) then
         map => grid%cell_map
         has_time = .false.
      else if (dimlens(2) == fin%ncells) then
         map => grid%cell_map
         has_time = .true.
      else if (dimlens(3) == fin%nEdges) then
         map => grid%edge_map
         has_time = .false.
      else if (dimlens(2) == fin%nEdges) then
         map => grid%edge_map
         has_time = .true.
      else if (dimlens(3) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .false.
      else if (dimlens(2) == fin%nVertices) then
         map => grid%vertex_map
         has_time = .true.
      else
         write (0,*) "Not sure which map to use, copy data mode"
         return
      end if
       
   
      !if (dimlens(1) * dimlens(2) > MAX_CHUNK_SIZE) then
         ! Too big, must handle in slices 
      !   if (has_time) then
      !      do i=1, dimlens(2)
      if(.false.) then         
         write (0,*) "FALSE"
      else
         call get_variable_3dREAL(fin, var_name, field)
         if (has_time) then
            allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(3)))
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, :, :) = field(:, map(i, j), :)
            end do
            end do   
         else
            allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(2)))
            do j=1, grid%ny
            do i=1, grid%nx
               newfield(i, j, :, :) = field(:, :, map(i, j))
            end do
            end do   
         end if

         call put_variable_4dREAL(fout, newfield, var_name)
      end if

   end subroutine copy_data_3dREAL


   !==================================================================================================
    integer function nearest_cell_path(target_lat, target_lon, start_cell, &
                                  nEdgesOnCell, cellsOnCell, latCell, lonCell)
   !==================================================================================================
    implicit none

    real (kind=RKIND), intent(in) :: target_lat, target_lon
    integer, intent(in) :: start_cell
    integer, dimension(:), pointer, intent(in) :: nEdgesOnCell
    integer, dimension(:,:), pointer, intent(in) :: cellsOnCell
    real (kind=RKIND), dimension(:), pointer, intent(in) :: latCell, lonCell

    integer :: i
    integer :: iCell
    integer :: current_cell
    real (kind=RKIND) :: current_distance, d
    real (kind=RKIND) :: nearest_distance

    nearest_cell_path = start_cell
    current_cell = -1

    do while (nearest_cell_path /= current_cell)
       current_cell = nearest_cell_path
       current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                          target_lon, 1.0_RKIND)
       nearest_cell_path = current_cell
       nearest_distance = current_distance
       do i = 1, nEdgesOnCell(current_cell)
          iCell = cellsOnCell(i,current_cell)
          d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_RKIND)
          if (d < nearest_distance) then
             nearest_cell_path = iCell
             nearest_distance = d
          end if
       end do
    end do

    end function nearest_cell_path

   real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                 cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance










end module utils














      
