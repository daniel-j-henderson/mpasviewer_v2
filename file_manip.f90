module mpas_file_manip 
   
   use netcdf
   use utils
   use params
   

   character(len=StrKIND), dimension(2), parameter :: &
      static_vars_1dINT =[character(len=StrKIND) :: 'nEdgesOnCell', &
                           'nEdgesOnEdge']

   character(len=StrKIND), dimension(8), parameter :: &
      static_vars_2dINT =[character(len=StrKIND) :: 'edgesOnCell', &
                           'cellsOnEdge',&
                           'cellsOnCell', &
                           'verticesOnCell',&
                           'verticesOnEdge',& 
                           'edgesOnVertex',&
                           'cellsOnVertex', &
                           'edgesOnEdge']

   character(len=StrKIND), dimension(21), parameter :: &
      static_vars_1dREAL=[character(len=StrKIND) :: 'latCell', 'lonCell',&
                           'xCell', 'yCell', 'zCell', &
                           'meshDensity', &
                           'latEdge', 'lonEdge', &
                           'xEdge', 'yEdge', 'zEdge', &
                           'latVertex', 'lonVertex', &
                           'xVertex', 'yVertex', 'zVertex', &
                           'dvEdge', 'dcEdge', &
                           'angleEdge', &
                           'areaCell', 'areaTriangle']

   character(len=StrKIND), dimension(2), parameter :: &
      static_vars_2dREAL=[character(len=StrKIND) :: 'weightsOnEdge', &
                           'kiteAreasOnVertex']
   

   type :: ncfile
   ! This type contains records about the ncfile so as to avoid repeating netcdf
   ! calls to see if a variable or dimension is present or something like that.
      integer :: ncid=0, ndims=0, nvars=0, natts=0, nCells=0, nEdges=0, nVertices=0
      character(len=StrKIND) :: filename
      character(len=StrKIND), dimension(:), pointer :: dims
      character(len=StrKIND), dimension(:), pointer :: vars
      character(len=StrKIND), dimension(:), pointer :: atts
      
      contains 
      procedure :: set_file_equal
      generic :: assignment(=) => set_file_equal
      procedure :: add_dim_record
      procedure :: add_var_record
      procedure :: add_att_record
      procedure :: contains_elem
      procedure :: is_open
      procedure :: clean

   end type ncfile

   contains

   subroutine clean(this)
      implicit none
      
      class(ncfile) :: this
      
      if (associated(this%dims)) deallocate(this%dims)
      if (associated(this%vars)) deallocate(this%vars)
      if (associated(this%atts)) deallocate(this%atts)

      this%filename = ' '
      this%ncid = 0
      this%ndims = 0
      this%nvars = 0
      this%natts = 0
      this%nCells = 0
      this%nEdges = 0
      this%nVertices = 0

   end subroutine clean

   subroutine set_file_equal(this, f)
      implicit none
      class(ncfile), intent(inout) :: this
      class(ncfile), intent(in) :: f

      this%ncid = f%ncid
      this%ndims = f%ndims
      this%nvars = f%nvars
      this%natts = f%natts
      this%nCells = f%nCells
      this%nEdges = f%nEdges
      this%nVertices = f%nVertices
      this%filename = f%filename
      
      if (associated(this%dims)) deallocate(this%dims)
      if (associated(this%vars)) deallocate(this%vars)
      if (associated(this%atts)) deallocate(this%atts)

      if (associated(f%dims)) then
         allocate(this%dims(size(f%dims)))
         this%dims = f%dims
      end if
      if (associated(f%vars)) then
         allocate(this%vars(size(f%vars)))
         this%vars = f%vars
      end if
      if (associated(f%atts)) then
         allocate(this%atts(size(f%atts)))
         this%atts = f%atts
      end if

   end subroutine set_file_equal






 

   logical function is_open(this)
      implicit none
      class(ncfile) :: this
   
      is_open = .false.
      if (this%ncid .ne. 0) is_open = .true.
      
   end function is_open

   logical function file_contains_elem(filename, type, elem_name)
      implicit none

      character(len=*) :: filename
      integer :: type
      character(len=*) :: elem_name

      integer :: ncid, ierr, el_id

      ierr = nf90_open(trim(filename), NF90_NOWRITE, ncid)
      if (ierr /= NF90_NOERR) then
         write (0,*) "Could not open file "//trim(filename)
         return
      end if

      select case(type)
      case(DIM)
         ierr = nf90_inq_dimid(ncid, trim(elem_name), el_id)
      case(VAR)
         ierr = nf90_inq_varid(ncid, trim(elem_name), el_id)
      case default
      end select
      
      if (ierr == NF90_NOERR) then
         file_contains_elem = .true.
      else
         file_contains_elem = .false.
      end if
   end function file_contains_elem

   logical function contains_elem(this, type, elem_name)
      implicit none

      class(ncfile) :: this
      integer :: type
      character(len=*) :: elem_name

      integer :: i, n
      character(len=StrKIND), dimension(:), pointer :: record

      select case(type)
      case(DIM)
         n = this%ndims
         record => this%dims
      case(VAR)
         n = this%nvars
         record => this%vars
      case(ATT)
         n = this%natts
         record => this%atts
      case default
      end select 

      contains_elem = .false.
      if (n == 0) return      

      do i=1, n
         if (trim(record(i)) == trim(elem_name)) then
            contains_elem = .true.
            exit
         end if
      end do
   end function contains_elem

   subroutine add_dim_record(this, dim_name)
      implicit none
   
      class(ncfile) :: this
      character(len=*) :: dim_name

      this%ndims = this%ndims+1
      this%dims(this%ndims) = trim(dim_name)
   end subroutine add_dim_record

   subroutine add_var_record(this, var_name)
      implicit none

      class(ncfile) :: this
      character(len=*) :: var_name

      this%nvars = this%nvars+1
      this%vars(this%nvars) = trim(var_name)
   end subroutine add_var_record

   subroutine add_att_record(this, att_name)
      implicit none

      class(ncfile) :: this
      character(len=*) :: att_name

      this%natts = this%natts+1
      this%atts(this%natts) = trim(att_name)
   end subroutine add_att_record

   logical function is_static(var_name)
      implicit none
      character(len=*) :: var_name

      integer :: i

      is_static = .false.
      
      do i=1, size(static_vars_1dINT)
         if (trim(static_vars_1dINT(i)) == trim(var_name)) is_static = .true.
      end do

      do i=1, size(static_vars_2dINT)
         if (trim(static_vars_2dINT(i)) == trim(var_name)) is_static = .true.
      end do

      do i=1, size(static_vars_1dREAL)
         if (trim(static_vars_1dREAL(i)) == trim(var_name)) is_static = .true.
      end do
      
      do i=1, size(static_vars_2dREAL)
         if (trim(static_vars_2dREAL(i)) == trim(var_name)) is_static = .true.
      end do

   end function

   ! netcdf file utility wrappers

   subroutine open_mpas_file(f, mode)
   ! If a file already exists, open it in read mode and extract some useful
   ! information about the file. Or create the file.
      type(ncfile), intent(inout) :: f
      character(len=*) :: mode
      integer :: ierr, var_id, dim_id, i, temp
      integer, dimension(:), allocatable :: ids
      character(len=StrKIND) :: elem_name
    
      allocate(f%dims(MAX_NDIMS), f%vars(MAX_NVARS), f%atts(MAX_NATTS))

      select case(mode)
      case('NF90_NOWRITE')

         ierr = nf90_open(trim(f%filename), NF90_NOWRITE, f%ncid)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error opening NetCDF file '//f%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
            stop
         end if

         ierr = nf90_inquire(f%ncid, f%ndims, f%nvars, f%natts)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring about '//f%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if

         if (allocated(ids)) deallocate(ids)
         allocate(ids(f%ndims))
         ierr = nf90_inq_dimids(f%ncid, temp, ids, i)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimids in '//f%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if

         do i=1, f%ndims
            ierr = nf90_inquire_dimension(f%ncid, ids(i), name=elem_name)
            if (ierr /= NF90_NOERR) then
               write(0,*) '*********************************************************************************'
               write(0,*) 'Error inquiring dim name in '//f%filename
               write(0,*) 'ierr = ', ierr
               write(0,*) '*********************************************************************************'
            end if

            f%dims(i) = elem_name
         end do

         if(allocated(ids)) deallocate(ids)
         allocate(ids(f%nvars))
         ierr = nf90_inq_varids(f%ncid, temp, ids)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring varids in '//f%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if

         do i=1, f%nvars
            ierr = nf90_inquire_variable(f%ncid, ids(i), name=elem_name)
            if (ierr /= NF90_NOERR) then
               write(0,*) '*********************************************************************************'
               write(0,*) 'Error inquiring var name in '//f%filename
               write(0,*) 'ierr = ', ierr
               write(0,*) '*********************************************************************************'
            end if

            f%vars(i) = elem_name
         end do


         ierr = nf90_inq_dimid(f%ncid, 'nCells', dim_id)
         if (ierr /= NF90_NOERR) then
            write(0,*) 'nCells not present in '//f%filename
            f%nCells = -1
         else
            ierr = nf90_inquire_dimension(f%ncid, dim_id, len=f%nCells)
            if (ierr /= NF90_NOERR) then
               write(0,*) '*********************************************************************************'
               write(0,*) 'Error inquiring nCells in '//f%filename
               write(0,*) 'ierr = ', ierr
               write(0,*) '*********************************************************************************'
               f%nCells = -1
            end if
         end if

         ierr = nf90_inq_dimid(f%ncid, 'nEdges', dim_id)
         if (ierr /= NF90_NOERR) then
            f%nCells = -1
            write(0,*) 'nEdges not present in '//f%filename
         else
            ierr = nf90_inquire_dimension(f%ncid, dim_id, len=f%nEdges)
            if (ierr /= NF90_NOERR) then
               write(0,*) '*********************************************************************************'
               write(0,*) 'Error inquiring nEdges in '//f%filename
               write(0,*) 'ierr = ', ierr
               write(0,*) '*********************************************************************************'
               f%nEdges = -1
            end if
         end if

         ierr = nf90_inq_dimid(f%ncid, 'nVertices', dim_id)
         if (ierr /= NF90_NOERR) then
            f%nVertices = -1
            write(0,*) 'nVertices not present in '//f%filename
         else 
            ierr = nf90_inquire_dimension(f%ncid, dim_id, len=f%nVertices)
            if (ierr /= NF90_NOERR) then
               write(0,*) '*********************************************************************************'
               write(0,*) 'Error inquiring nVertices in '//f%filename
               write(0,*) 'ierr = ', ierr
               write(0,*) '*********************************************************************************'
               f%nVertices = -1
            end if
         end if


      case('CREATE')
         ierr = nf90_create(f%filename, NF90_CLOBBER, f%ncid)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error creating file '//f%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if
      
      case default
         write (0,*) "Bad open mode"
         stop
      end select

   end subroutine open_mpas_file

   subroutine get_dimension(f, dim_name, field)
   ! More precisely, get dimension length
      implicit none
   
      type(ncfile) :: f
      character(len=*), intent(in) :: dim_name
      integer :: field
      integer :: ierr, id
      
      ierr = nf90_inq_dimid(f%ncid, trim(dim_name), id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring '//trim(dim_name)//' dimid in'//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, id, len=field)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension '//trim(dim_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'         
      end if
   end subroutine get_dimension

   subroutine add_dimension(f, dim_name, field)
   ! Define a dimension in an nc file. Must already be in define mode.
      implicit none

      type(ncfile) :: f
      character(len=*) :: dim_name
      integer :: field, ierr, dimid

      if (.not. f%ndims < MAX_NDIMS) then
         write (0,*) "ERROR: Trying to add too many dimensions to "//f%filename
         return
      end if

      if (f%contains_elem(DIM, dim_name)) then
         write (0,*) "Already contains dimension "//trim(dim_name)//", skipping the add."
         return
      end if

      call f%add_dim_record(dim_name)

      ierr = nf90_def_dim(f%ncid, dim_name, field, dimid)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error defining dimension '//trim(dim_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'         
      end if

   end subroutine add_dimension

   subroutine copy_dimensions(ncin, ncout)
   ! Copy all dimensions, save for n- Cells, Edges, Vertices, from one file to
   ! another.
      type(ncfile) :: ncin, ncout

      integer :: i, field

      do i=1, ncin%ndims
         if (ncout%ndims == MAX_NDIMS) write (0,*)"ERROR: putting too many dims in ncout"
         if (trim(ncin%dims(i)) == 'nCells' .or. trim(ncin%dims(i)) == &
            'nEdges' .or. trim(ncin%dims(i)) == 'nVertices') cycle
         if (ncout%contains_elem(DIM, ncin%dims(i))) cycle
         call get_dimension(ncin, ncin%dims(i), field)
         call add_dimension(ncout, ncin%dims(i), field)
      end do

   end subroutine copy_dimensions

   subroutine copy_variable_defmode(ncin, ncout, var_name)
   ! Copy a variable definition from one file to another. Data part must be done
   ! in data mode.
      implicit none
   
      class(ncfile) :: ncin, ncout
      character(len=*) :: var_name

      integer :: ierr, var_id, xtype, ndims, i
      integer, dimension(:), allocatable :: dimids, newdimids
      character(len=StrKIND), dimension(:), allocatable :: dims

      ierr = nf90_inq_varid(ncin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if 

      ierr = nf90_inquire_variable(ncin%ncid, var_id, xtype=xtype, ndims=ndims)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if


      if (ndims == 0) then
         ierr = nf90_def_var(ncout%ncid, var_name, xtype, varid=var_id)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error defining variable '//trim(var_name)//' in copy def mode'
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
            stop
         end if
         return
      end if
   
      allocate(dimids(ndims), dims(ndims), newdimids(ndims))

      ierr = nf90_inquire_variable(ncin%ncid, var_id, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      do i=1, ndims
         ierr = nf90_inquire_dimension(ncin%ncid, dimids(i), dims(i))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimension in copy def mode'
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if

         if (.not. ncout%contains_elem(DIM, dims(i))) then
            write (0,*) "Trying to copy a variable whose dimensions are not in ncout"
            return
         end if 

         ierr = nf90_inq_dimid(ncout%ncid, dims(i), newdimids(i))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimension in copy def mode'
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if
      end do
      
      ierr = nf90_def_var(ncout%ncid, var_name, xtype, newdimids, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error defining variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

   end subroutine copy_variable_defmode

   subroutine copy_variable_datamode(ncin, ncout, var_name, cell_map, edge_map, vertex_map)
   ! Once in data mode, this function will copy a variable's data from one file
   ! to another, assuming it is dimensioned by nCells, nEdges, or nVertices (and
   ! not time)
      implicit none
      
      class(ncfile) :: ncin, ncout
      character(len=*) :: var_name

      integer :: ierr, i, var_id, ndims, nelems, xtype
      integer, dimension(:), pointer :: map, cell_map, edge_map, vertex_map, dimlens
      integer, dimension(:), pointer :: field_1dINT, newfield_1dINT
      integer, dimension(:,:), pointer :: field_2dINT, newfield_2dINT
      integer, dimension(:,:,:), pointer :: field_3dINT, newfield_3dINT
      real(kind=RKIND), dimension(:), pointer :: field_1dREAL, newfield_1dREAL
      real(kind=RKIND), dimension(:,:), pointer :: field_2dREAL, newfield_2dREAL
      real(kind=RKIND), dimension(:,:,:), pointer :: field_3dREAL, newfield_3dREAL
      character(len=StrKIND) :: field_1dCHAR
      character(len=StrKIND), dimension(:), pointer :: field_2dCHAR
      character(len=StrKIND), dimension(:,:), pointer :: field_3dCHAR
      logical :: has_time

      real(kind=RKIND) :: const

      ierr = nf90_inq_varid(ncin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(ncin%ncid, var_id, xtype=xtype, ndims=ndims)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in copy def mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      if (ndims == 0) then
         ierr = nf90_get_var(ncin%ncid, var_id, const)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error getting variable '//trim(var_name)//' in '//ncin%filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
            stop
         end if
         ierr = nf90_inq_varid(ncout%ncid, trim(var_name), var_id)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring varID of '//trim(var_name)//' in copy def mode'
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
            stop
         end if
         ierr = nf90_put_var(ncout%ncid, var_id, const)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error putting variable '//trim(var_name)//' putting constant var'
            write(0,*) 'ierr = ', ierr
            write(0,*) 'shape(field):',  ncout%nEdges
            write(0,*) '*********************************************************************************'
         end if  
         return
      end if

      allocate(dimlens(ndims))

      has_time = .false.
      if (xtype == NF90_INT) then
         select case(ndims)
         case(1)
            call get_variable_1dINT(ncin, var_name, field_1dINT)
            if(size(field_1dINT) == ncin%nCells) then
               map => cell_map
            else if (size(field_1dINT) == ncin%nEdges) then
               map => edge_map
            else if (size(field_1dINT) == ncin%nVertices) then
               map => vertex_map
            else
            !   write (0,*) "Not sure which map to use, copy data mode"
               call put_variable_1dINT(ncout, field_1dINT, var_name)
            end if

            allocate(newfield_1dINT(size(map)))
            call compact_field_1dINT(field_1dINT, newfield_1dINT, map)
            call put_variable_1dINT(ncout, newfield_1dINT, var_name)
         case(2)
            call get_variable_2dINT(ncin, var_name, field_2dINT)
            dimlens = shape(field_2dINT)
            if(dimlens(2) == ncin%nCells) then
               map => cell_map
               has_time = .false.
            else if (dimlens(1) == ncin%ncells) then
               map => cell_map
               has_time = .true.
            else if (dimlens(2) == ncin%nEdges) then
               map => edge_map
               has_time = .false.
            else if (dimlens(1) == ncin%nEdges) then
               map => edge_map
               has_time = .true.
            else if (dimlens(2) == ncin%nVertices) then
               map => vertex_map
               has_time = .false.
            else if (dimlens(1) == ncin%nVertices) then
               map => vertex_map
               has_time = .true.
            else
               write (0,*) "Not sure which map to use, copy data mode"
            end if

            if (has_time) then
               allocate(newfield_2dINT(size(map), dimlens(2)))
               do i = 1, dimlens(2)
                  field_1dINT => field_2dINT(:,i)
                  newfield_1dINT => newfield_2dINT(:,i)
                  call compact_field_1dINT(field_1dINT, newfield_1dINT, map)
               end do
            else 
               allocate(newfield_2dINT(dimlens(1), size(map)))
               call compact_field_2dINT(field_2dINT, newfield_2dINT, map)
               return
            end if

            call put_variable_2dINT(ncout, newfield_2dINT, var_name)
         case(3)
            call get_variable_3dINT(ncin, var_name, field_3dINT)
            dimlens = shape(field_3dINT)
            if(dimlens(3) == ncin%nCells) then
               map => cell_map
               has_time = .false.
            else if (dimlens(2) == ncin%ncells) then
               map => cell_map
               has_time = .true.
            else if (dimlens(3) == ncin%nEdges) then
               map => edge_map
               has_time = .false.
            else if (dimlens(2) == ncin%nEdges) then
               map => edge_map
               has_time = .true.
            else if (dimlens(3) == ncin%nVertices) then
               map => vertex_map
               has_time = .false.
            else if (dimlens(2) == ncin%nVertices) then
               map => vertex_map
               has_time = .true.
            else
               write (0,*) "Not sure which map to use, copy data mode"
            end if

            if (has_time) then
               allocate(newfield_3dINT(dimlens(1), size(map), dimlens(3)))
               do i = 1, dimlens(3)
                  field_2dINT => field_3dINT(:,:,i)
                  newfield_2dINT => newfield_3dINT(:,:,i)
                  call compact_field_2dINT(field_2dINT, newfield_2dINT, map)
               end do
            else 
               allocate(newfield_3dINT(dimlens(1), dimlens(2), size(map)))
               call compact_field_3dINT(field_3dINT, newfield_3dINT, map)
            end if

            call put_variable_3dINT(ncout, newfield_3dINT, var_name)
         case default 
            write (0,*) "Error in case for copy data mode"
         end select
      else if (xtype == NF90_FLOAT .or. xtype == NF90_DOUBLE) then
         select case(ndims)
         case(1)
            call get_variable_1dREAL(ncin, trim(var_name), field_1dREAL)
            if(size(field_1dREAL) == ncin%nCells) then
               map => cell_map
            else if (size(field_1dREAL) == ncin%nEdges) then
               map => edge_map
            else if (size(field_1dREAL) == ncin%nVertices) then
               map => vertex_map
            else
               ! write (0,*) "Not sure which map to use, copy data mode"
               call put_variable_1dREAL(ncout, field_1dREAL, var_name)
               return
            end if

            allocate(newfield_1dREAL(size(map)))
            call compact_field_1dREAL(field_1dREAL, newfield_1dREAL, map)
            call put_variable_1dREAL(ncout, newfield_1dREAL, trim(var_name))
         case(2)
            call get_variable_2dREAL(ncin, var_name, field_2dREAL)
            dimlens = shape(field_2dREAL)
            has_time = .false.
            if(dimlens(2) == ncin%nCells) then
               map => cell_map
               has_time = .false.
            else if (dimlens(1) == ncin%nCells) then
               map => cell_map
               has_time = .true.
            else if (dimlens(2) == ncin%nEdges) then
               map => edge_map
               has_time = .false.
            else if (dimlens(1) == ncin%nEdges) then  
               map => edge_map
               has_time = .true.
            else if (dimlens(2) == ncin%nVertices) then
               map => vertex_map
               has_time = .false.
            else if (dimlens(1) == ncin%nVertices) then
               map => vertex_map
               has_time = .true.
            else
               write (0,*) "Not sure which map to use, copy data mode"
            end if
            if (has_time) then
               allocate(newfield_2dREAL(size(map), dimlens(2)))
               do i = 1, dimlens(2)
                  field_1dREAL => field_2dREAL(:,i)
                  newfield_1dREAL => newfield_2dREAL(:,i)
                  call compact_field_1dREAL(field_1dREAL, newfield_1dREAL, map)
               end do
            else 
               allocate(newfield_2dREAL(dimlens(1), size(map)))
               call compact_field_2dREAL(field_2dREAL, newfield_2dREAL, map)
            end if
            call put_variable_2dREAL(ncout, newfield_2dREAL, var_name)
         case(3)
            call get_variable_3dREAL(ncin, var_name, field_3dREAL)
            dimlens = shape(field_3dREAL)
            if(dimlens(3) == ncin%nCells) then
               map => cell_map
               has_time = .false.
            else if (dimlens(2) == ncin%nCells) then
               map => cell_map
               has_time = .true.
            else if (dimlens(3) == ncin%nEdges) then
               map => edge_map
               has_time = .false.
            else if (dimlens(2) == ncin%nEdges) then
               map => edge_map
               has_time = .true.
            else if (dimlens(3) == ncin%nVertices) then
               map => vertex_map
               has_time = .false.
            else if (dimlens(2) == ncin%nVertices) then
               map => vertex_map
               has_time = .true.
            else
               write (0,*) "Not sure which map to use, copy data mode"
            end if

            if (has_time) then
               allocate(newfield_3dREAL(dimlens(1), size(map), dimlens(3)))
               do i = 1, dimlens(3)
                  field_2dREAL => field_3dREAL(:,:,i)
                  newfield_2dREAL => newfield_3dREAL(:,:,i)
                  call compact_field_2dREAL(field_2dREAL, newfield_2dREAL, map)
               end do
            else 
               allocate(newfield_3dREAL(dimlens(1), dimlens(2), size(map)))
               call compact_field_3dREAL(field_3dREAL, newfield_3dREAL, map)
            end if
            call put_variable_3dREAL(ncout, newfield_3dREAL, var_name)
         case default 
            write (0,*) "Error in case for copy data mode"
         end select
      else if (xtype == NF90_CHAR) then
         select case(ndims)
         case(1)
            field_1dCHAR = ' ' 
            call get_variable_1dCHAR(ncin, var_name, field_1dCHAR)
            call put_variable_1dCHAR(ncout, field_1dCHAR, var_name)
         case(2)
            field_2dCHAR = ' ' 
            call get_variable_2dCHAR(ncin, var_name, field_2dCHAR)
            call put_variable_2dCHAR(ncout, field_2dCHAR, var_name)
         case(3)
            field_3dCHAR = ' ' 
            call get_variable_3dCHAR(ncin, var_name, field_3dCHAR)
            call put_variable_3dCHAR(ncout, field_3dCHAR, var_name)
         case default
            write (0,*) "Variable "//trim(var_name)//" is not supported"
         end select
      else
         write (0,*) "in copy data mode, neither real nor int"
         write (0,*) "xtype, REAL, DOUBLE, FLOAT, INT", xtype, NF90_REAL, NF90_DOUBLE, NF90_FLOAT, NF90_INT
      end if
   end subroutine copy_variable_datamode  
       


   subroutine create_variable_1dINT(f, var_name, dim_name)
      implicit none
      
      type(ncfile) :: f
      character(len=*) :: dim_name, var_name

      integer :: ierr, dimid, varid

      ierr = nf90_inq_dimid(f%ncid, dim_name, dimid)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimID of '//trim(dim_name)//' in create_variable_1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
      if (.not. f%nvars < MAX_NVARS) then
         write (0,*) "ERROR: Trying to add too many variables to "//f%filename
         return
      end if

      if (f%contains_elem(VAR, var_name)) then
         write (0,*) "Already contains variable "//trim(var_name)//", skipping the add."
         return
      end if

      call f%add_var_record(var_name)

      ierr = nf90_def_var(f%ncid, var_name, NF90_INT, (/dimid/), varid)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error defining variable'//trim(var_name)//'in create_variable_1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
   end subroutine create_variable_1dINT

   subroutine create_static_field_1d(ncout, ncin, var_name)
   ! More accurately, copy static field definition from one file to another
      implicit none
      type(ncfile) :: ncout, ncin
      character(len=*) :: var_name

      integer :: ierr, var_id, n, xtype, temp
      integer, dimension(1) :: dimids
      character(len=StrKIND) :: dim_name
      
      ierr = nf90_inq_varid(ncin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in create_static_field_1d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(ncin%ncid, var_id, xtype=xtype, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable in create_static_field_1d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(ncin%ncid, dimids(1), dim_name, n) 
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension in create_static_field_1d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inq_dimid(ncout%ncid, dim_name, dimids(1))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimID of '//trim(dim_name)//' in create_static_field_1d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      if (.not. ncout%nvars < MAX_NVARS) then
         write (0,*) "ERROR: Trying to add too many variables to "//ncout%filename
         return
      end if

      if (ncout%contains_elem(VAR, var_name)) then
         write (0,*) "Already contains variable "//trim(var_name)//", skipping the add."
         return
      end if

      call ncout%add_var_record(var_name)

      ierr = nf90_def_var(ncout%ncid, var_name, xtype, dimids, temp)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error defining variable '//trim(var_name)//' in create_static_field_1d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if
   
   end subroutine create_static_field_1d

   subroutine create_static_field_2d(ncout, ncin, var_name)
      type(ncfile) :: ncout, ncin
      character(len=StrKIND) :: var_name

      integer :: ierr, var_id, xtype, temp
      integer, dimension(2) :: dimids
      character(len=StrKIND) :: dim_name
      
      ierr = nf90_inq_varid(ncin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in create_static_field_2d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(ncin%ncid, var_id, xtype=xtype, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable in create_static_field_2d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(ncin%ncid, dimids(2), dim_name) 
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension in create_static_field_2d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inq_dimid(ncout%ncid, dim_name, dimids(2))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimID of '//trim(dim_name)//' in create_static_field_2d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_dimension(ncin%ncid, dimids(1), dim_name) 
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension in create_static_field_2d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inq_dimid(ncout%ncid, dim_name, dimids(1))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimID of '//trim(dim_name)//' in create_static_field_2d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      if (.not. ncout%nvars < MAX_NVARS) then
         write (0,*) "ERROR: Trying to add too many variables to "//ncout%filename
         return
      end if

      if (ncout%contains_elem(VAR, var_name)) then
         write (0,*) "Already contains variable "//trim(var_name)//", skipping the add."
         return
      end if

      call ncout%add_var_record(var_name)

      ierr = nf90_def_var(ncout%ncid, var_name, xtype, dimids, temp)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error defining variable in create_static_field_2d'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if
   
   end subroutine create_static_field_2d



   subroutine copyandcompact_static_fields(ncin, ncout, cell_map, edge_map, vertex_map, icell_map, iedge_map, ivertex_map)
   ! Copy the variable definitions, compact the data, then place the data for
   ! all of the static fields defined at the top of this module. Assumes we
   ! start in define mode for ncout.
      type(ncfile) :: ncin, ncout
      integer, dimension(:), pointer :: cell_map, edge_map, vertex_map, icell_map, iedge_map, ivertex_map

      integer, dimension(:), pointer :: map, imap
      integer, dimension(:), pointer :: field_1dINT, newfield_1dINT
      integer, dimension(:,:), pointer :: field_2dINT, newfield_2dINT
      real(kind=RKIND), dimension(:), pointer :: field_1dREAL, newfield_1dREAL
      real(kind=RKIND), dimension(:,:), pointer :: field_2dREAL, newfield_2dREAL
      integer :: id, ierr
      integer, dimension(1) :: dim_ids
      integer, dimension(2) :: dimensions, dims

      do i=1, size(static_vars_1dINT)
         if(.not. ncout%contains_elem(VAR, trim(static_vars_1dINT(i)))) &
            call create_static_field_1d(ncout, ncin, static_vars_1dINT(i))
      end do

      do i=1, size(static_vars_2dINT)
         if(.not. ncout%contains_elem(VAR, trim(static_vars_2dINT(i)))) &
            call create_static_field_2d(ncout, ncin, static_vars_2dINT(i))
      end do

      do i=1, size(static_vars_1dREAL)
         if(.not. ncout%contains_elem(VAR, trim(static_vars_1dREAL(i)))) &
            call create_static_field_1d(ncout, ncin, static_vars_1dREAL(i))
      end do

      do i=1, size(static_vars_2dREAL)
         if(.not. ncout%contains_elem(VAR, trim(static_vars_2dREAL(i)))) &
            call create_static_field_2d(ncout, ncin, static_vars_2dREAL(i))
      end do

      
      ! switch ncout to data mode

      ierr = nf90_enddef(ncout%ncid)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error ending define mode'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      do i=1, size(static_vars_1dINT) 
         call get_variable_1dINT(ncin, static_vars_1dINT(i), field_1dINT)
         if (size(field_1dINT) == ncin%nCells) then
            allocate(newfield_1dINT(ncout%nCells))
            call compact_field_1dINT(field_1dINT, newfield_1dINT, cell_map)
            call put_variable_1dINT(ncout, newfield_1dINT, static_vars_1dINT(i))
         else if (size(field_1dINT) == ncin%nEdges) then
            allocate(newfield_1dINT(ncout%nEdges))
            call compact_field_1dINT(field_1dINT, newfield_1dINT, edge_map)
            if (ncout%nEdges .ne. size(newfield_1dINT)) write (0,*) 'dims dont match, nEdges, size(field):', ncout%nEdges, size(newfield_1dINT)
            call put_variable_1dINT(ncout, newfield_1dINT, static_vars_1dINT(i))
         else if (size(field_1dINT) == ncin%nVertices) then
            allocate(newfield_1dINT(ncout%nVertices))
            call compact_field_1dINT(field_1dINT, newfield_1dINT, vertex_map)
            call put_variable_1dINT(ncout, newfield_1dINT, static_vars_1dINT(i))
         end if 
         deallocate(newfield_1dINT)

      end do

      do i=1, size(static_vars_2dINT)
         call get_variable_2dINT(ncin, static_vars_2dINT(i), field_2dINT)

         dimensions = shape(field_2dINT)
         nullify(map, imap)
         if (dimensions(2) == ncin%nCells) then
            map => cell_map
         else if (dimensions(2) == ncin%nEdges) then
            map => edge_map
         else if (dimensions(2) == ncin%nVertices) then
            map => vertex_map
         else
            write (0,*) "Not sure which map to use"
         end if 

         if (abs(maxval(field_2dINT) - ncin%nCells) < 2) then
            imap => icell_map
         else if (abs(maxval(field_2dINT) - ncin%nEdges) < 2) then
            imap => iedge_map
         else if (abs(maxval(field_2dINT) - ncin%nVertices) < 2) then
            imap => ivertex_map
         else
            write (0,*) "Not sure which imap to use"
         end if

         allocate(newfield_2dINT(dimensions(1), size(map)))
         call compact_field_2dINT(field_2dINT, newfield_2dINT, map)
         call reindex_field_2dINT(newfield_2dINT, imap)
         where(newfield_2dINT == 0) newfield_2dINT = 1
         call put_variable_2dINT(ncout, newfield_2dINT, static_vars_2dINT(i))
         deallocate(newfield_2dINT)
      end do


      do i=1, size(static_vars_1dREAL)
         call get_variable_1dREAL(ncin, static_vars_1dREAL(i), field_1dREAL)
         if (size(field_1dREAL) == ncin%nCells) then
            allocate(newfield_1dREAL(size(cell_map)))
            call compact_field_1dREAL(field_1dREAL, newfield_1dREAL, cell_map)
            call put_variable_1dREAL(ncout, newfield_1dREAL, static_vars_1dREAL(i))
         else if (size(field_1dREAL) == ncin%nEdges) then
            allocate(newfield_1dREAL(size(edge_map)))
            call compact_field_1dREAL(field_1dREAL, newfield_1dREAL, edge_map)
            call put_variable_1dREAL(ncout, newfield_1dREAL, static_vars_1dREAL(i))
         else if (size(field_1dREAL) == ncin%nVertices) then
            allocate(newfield_1dREAL(size(vertex_map)))
            call compact_field_1dREAL(field_1dREAL, newfield_1dREAL, vertex_map)
            call put_variable_1dREAL(ncout, newfield_1dREAL, static_vars_1dREAL(i))
         end if 
         deallocate(newfield_1dREAL)
      end do


      do i=1, size(static_vars_2dREAL)
         call get_variable_2dREAL(ncin, static_vars_2dREAL(i), field_2dREAL)

         dimensions = shape(field_2dREAL)
         if (dimensions(2) == ncin%nCells) then
            map => cell_map
            imap => icell_map
         else if (dimensions(2) == ncin%nEdges) then
            map => edge_map
            imap => iedge_map
         else if (dimensions(2) == ncin%nVertices) then
            map => vertex_map
            imap => ivertex_map
         end if 
         
         allocate(newfield_2dREAL(dimensions(1), size(map)))
         call compact_field_2dREAL(field_2dREAL, newfield_2dREAL, map)

         call put_variable_2dREAL(ncout, newfield_2dREAL, static_vars_2dREAL(i))
         deallocate(newfield_2dREAL)
      end do

   end subroutine copyandcompact_static_fields


   subroutine put_variable_1dREAL(f, field, var_name)
      type(ncfile) :: f 
      real(kind=RKIND), dimension(:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id, n
      integer, dimension(1) :: temp
      
      ierr = nf90_inq_varid(f%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 1dREAL'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=temp)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, temp(1), len=n)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      if(size(field) .ne. n) write (0,*) "size(field), n:", size(field), n, var_name

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error putting variable '//trim(var_name)//' in put 1dREAL'
         write(0,*) 'ierr = ', ierr
         write(0,*) 'shape(field), f%nEdges:', shape(field), f%nEdges
         write(0,*) '*********************************************************************************'
      end if  
   end subroutine put_variable_1dREAL

   
   subroutine put_variable_1dINT(f, field, var_name)
      type(ncfile) :: f
      integer, dimension(:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id, n
      integer, dimension(1) :: temp
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=temp)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, temp(1), len=n)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      if(size(field) .ne. n) write (0,*) "size(field), n:", size(field), n, var_name

      ierr = nf90_put_var(f%ncid, var_id, field, start=(/1/), count=(/size(field)/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error putting variable '//trim(var_name)//' in put 1dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) 'shape(field), f%nEdges:', shape(field), f%nEdges
         write(0,*) '*********************************************************************************'
      end if  

   end subroutine put_variable_1dINT


   subroutine put_variable_2dREAL(f, field, var_name)
      type(ncfile) :: f
      real(kind=RKIND), dimension(:,:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 2dREAL'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error putting variable '//trim(var_name)//' in put 2dREAL'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if  
   end subroutine put_variable_2dREAL

   subroutine put_variable_3dREAL(f, field, var_name)
      type(ncfile) :: f
      real(kind=RKIND), dimension(:,:,:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 3dREAL'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error putting variable '//trim(var_name)//' in put 3dREAL'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if  
   end subroutine put_variable_3dREAL

   subroutine put_variable_2dINT(f, field, var_name)
      type(ncfile) :: f
      integer, dimension(:,:), pointer :: field
      character(len=StrKIND) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 2dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if
      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error putting variable '//trim(var_name)//' in put 2dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if  
   end subroutine put_variable_2dINT

   subroutine put_variable_3dINT(f, field, var_name)
      type(ncfile) :: f
      integer, dimension(:,:,:), pointer :: field
      character(len=StrKIND) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 3dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if
      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error putting variable '//trim(var_name)//' in put 3dINT'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if  
   end subroutine put_variable_3dINT

   subroutine put_variable_1dCHAR(f, field, var_name)
      type(ncfile) :: f
      character(len=StrKIND) :: field
      character(len=StrKIND) :: var_name

      integer :: var_id, n, i
      integer, dimension(2) :: temp

      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 1dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      n = len(trim(field))
      ierr = nf90_put_var(f%ncid, var_id, trim(field), (/1/), (/n/), (/1/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error putting variable '//trim(var_name)//' in put 1dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if  
   end subroutine put_variable_1dCHAR

   subroutine put_variable_2dCHAR(f, field, var_name)
      type(ncfile) :: f
      character(len=StrKIND), dimension(:), pointer :: field
      character(len=StrKIND) :: var_name

      integer :: var_id, n, i
      integer, dimension(2) :: temp
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 2dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=temp)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 2dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, temp(2), len=n)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 2dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if
      
      do i=1, n
         ierr = nf90_put_var(f%ncid, var_id, trim(field(i)))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error putting variable '//trim(var_name)//' in put 2dCHAR'
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if  
      end do
   end subroutine put_variable_2dCHAR

   subroutine put_variable_3dCHAR(f, field, var_name)
      type(ncfile) :: f
      character(len=StrKIND), dimension(:,:), pointer :: field
      character(len=StrKIND) :: var_name

      integer :: var_id, n1, n2, i, j
      integer, dimension(3) :: temp
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 3dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=temp)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 3dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, temp(2), len=n1)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 3dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, temp(3), len=n2)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in put 3dCHAR'
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if
      
      do i=1, n2
      do j=1, n1
         ierr = nf90_put_var(f%ncid, var_id, trim(field(j,i)))
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error putting variable '//trim(var_name)//' in put 3dCHAR'
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if  
      end do
      end do
   end subroutine put_variable_3dCHAR

   subroutine get_variable_1dINT(f, var_name, field)
      implicit none
      
      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      integer, dimension(:), pointer, intent(inout) :: field

      integer :: var_id, n, ierr
      integer, dimension(1) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1/), count = (/n/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
   end subroutine get_variable_1dINT

   subroutine get_variable_2dINT(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      integer, dimension(:,:), pointer, intent(inout) :: field

      integer :: var_id, n1, n2, ierr
      integer, dimension(2) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n1, n2))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1/), count = (/n1, n2/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
   end subroutine get_variable_2dINT

   subroutine get_variable_3dINT(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      integer, dimension(:,:,:), pointer, intent(inout) :: field

      integer :: var_id, n1, n2, n3, ierr
      integer, dimension(3) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(2), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(3), len=n3)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(3), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n1, n2, n3))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1,1/), count = (/n1, n2, n3/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
   end subroutine get_variable_3dINT

   subroutine get_variable_1dREAL(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      real(kind=RKIND), dimension(:), pointer, intent(inout) :: field

      integer :: var_id, n, ierr, xtype
      integer, dimension(1) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids, xtype=xtype)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      if (xtype .ne. NF90_DOUBLE) write (0,*) "var "//trim(var_name)//" is not a double"

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n))
      ierr = nf90_get_var(f%ncid, var_id, field, count = (/n/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
   end subroutine get_variable_1dREAL

   subroutine get_variable_2dREAL(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      real(kind=RKIND), dimension(:,:), pointer, intent(inout) :: field

      integer :: var_id, n1, n2, ierr
      integer, dimension(2) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n1, n2))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1/), count = (/n1, n2/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
   end subroutine get_variable_2dREAL

   subroutine get_variable_3dREAL(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      real(kind=RKIND), dimension(:,:,:), pointer, intent(inout) :: field

      integer :: var_id, n1, n2, n3, ierr
      integer, dimension(3) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(2), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(3), len=n3)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(3), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n1, n2, n3))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1,1/), count = (/n1, n2, n3/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
   end subroutine get_variable_3dREAL

   subroutine get_variable_1dCHAR(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      character(len=StrKIND), intent(out) :: field

      integer :: var_id, ierr, dimid, strlen

      ierr = nf90_inq_dimid(f%ncid, 'StrLen', dimid)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimID of StrLen in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_dimension(f%ncid, dimid, len=strlen)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring StrLen in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_get_var(f%ncid, var_id, field(1:strlen), start=(/1/), count=(/strlen/), stride=(/1/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

   end subroutine get_variable_1dCHAR

   subroutine get_variable_2dCHAR(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      character(len=StrKIND), dimension(:), pointer, intent(inout) :: field

      integer :: var_id, n1, n2, ierr
      integer, dimension(2) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n2))
      field(:) = ' '
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1/), count = (/n1, n2/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      
   end subroutine get_variable_2dCHAR

   subroutine get_variable_3dCHAR(f, var_name, field)
      implicit none

      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      character(len=StrKIND), dimension(:,:), pointer, intent(inout) :: field

      integer :: var_id, n1, n2, n3, ierr
      integer, dimension(3) :: dim_ids

      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring varID of '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if 

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(1), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(2), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(3), len=n3)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring dimension with id ', dim_ids(3), ' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
      end if

      allocate(field(n2, n3))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1,1/), count = (/n1, n2, n3/))
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error getting variable '//trim(var_name)//' in '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if

      
   end subroutine get_variable_3dCHAR

   subroutine get_attribute_REAL(f, att_name, field)
      implicit none
      
      type(ncfile) :: f 
      character(len=*), intent(in) :: att_name
      real (kind=RKIND), intent(out) :: field
      integer :: ierr
      ierr = nf90_get_att(f%ncid, NF90_GLOBAL, trim(att_name), field)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error inquiring attribute '//trim(att_name)//' in file '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
   end subroutine get_attribute_REAL
   
   subroutine copy_attributes(ncin, ncout)
   ! TODO: find a way to copy all attributes from one file to another without
   ! knowing what they are. 
      implicit none
   
      type(ncfile) :: ncin, ncout

      integer :: ierr, i
      character(len=StrKIND) :: att_name

      do i=1,ncin%natts
         ierr = nf90_inq_attname(ncin%ncid, NF90_GLOBAL, i, att_name)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring attribute ', i 
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
         end if
         ierr = nf90_copy_att(ncin%ncid, NF90_GLOBAL, att_name, ncout%ncid, NF90_GLOBAL)
         if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error copying attribute '//trim(att_name) 
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
            stop
         end if
      end do
   end subroutine copy_attributes


   subroutine close_mpas_file(f)
      implicit none
      
      type(ncfile) :: f
      integer :: ierr

      ierr = nf90_close(f%ncid)
      if (ierr /= NF90_NOERR) then
         write(0,*) '*********************************************************************************'
         write(0,*) 'Error while closing NetCDF file '//f%filename
         write(0,*) 'ierr = ', ierr
         write(0,*) '*********************************************************************************'
         stop
      end if
      call f%clean()
   end subroutine

end module mpas_file_manip 
   
