module mpas_file_manip 
   
   use netcdf
   !use utils
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

      this%filename = ''
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
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_open', .false., 'open_mpas_file', f%filename)

         ierr = nf90_inquire(f%ncid, f%ndims, f%nvars, f%natts)
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire', .false., 'open_mpas_file', f%filename)

         if (allocated(ids)) deallocate(ids)
         allocate(ids(f%ndims))
         ierr = nf90_inq_dimids(f%ncid, temp, ids, i)
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_dimids', .false., 'open_mpas_file', f%filename)

         do i=1, f%ndims
            ierr = nf90_inquire_dimension(f%ncid, ids(i), name=elem_name)
            if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .false., 'open_mpas_file', f%filename)

            f%dims(i) = elem_name
         end do

         if(allocated(ids)) deallocate(ids)
         allocate(ids(f%nvars))
         ierr = nf90_inq_varids(f%ncid, temp, ids)
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varids', .false., 'open_mpas_file', f%filename)

         do i=1, f%nvars
            ierr = nf90_inquire_variable(f%ncid, ids(i), name=elem_name)
            if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .false., 'open_mpas_file', f%filename)

            f%vars(i) = elem_name
         end do


         ierr = nf90_inq_dimid(f%ncid, 'nCells', dim_id)
         if (ierr /= NF90_NOERR) then
            write(0,*) 'nCells not present in '//f%filename
            f%nCells = -1
         else
            ierr = nf90_inquire_dimension(f%ncid, dim_id, len=f%nCells)
            if (ierr /= NF90_NOERR) then
               call handle_err(ierr, 'nf90_inquire dimension', .false., 'open_mpas_file', f%filename)
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
               call handle_err(ierr, 'nf90_inquire dimension', .false., 'open_mpas_file', f%filename)
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
               call handle_err(ierr, 'nf90_inquire dimension', .false., 'open_mpas_file', f%filename)
               f%nVertices = -1
            end if
         end if


      case('CREATE')
         ierr = nf90_create(f%filename, NF90_CLOBBER, f%ncid)
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_create', .true., 'open_mpas_file', f%filename)
      
      case default
         write (0,*) "Bad open mode for "//trim(f%filename)
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
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_dimid', .false., 'get_dimension', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, id, len=field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .false., 'get_dimension', f%filename)
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
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_def_dim', .false., 'add_dimension', f%filename)

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
      character(len=StrKIND), intent(inout) :: var_name

      integer :: ierr, var_id, xtype, ndims, i, j
      integer, dimension(:), allocatable :: dimids, newdimids
      character(len=StrKIND), dimension(:), allocatable :: dims
      logical :: is_spatial=.false.

      ierr = nf90_inq_varid(ncin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .true., 'copy_variable_defmode', ncin%filename)

      ierr = nf90_inquire_variable(ncin%ncid, var_id, xtype=xtype, ndims=ndims)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'copy_variable_defmode', ncin%filename)


      if (ndims == 0) then
         ierr = nf90_def_var(ncout%ncid, var_name, xtype, varid=var_id)
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_def_var', .true., 'copy_variable_defmode', ncout%filename)
         return
      end if
   
      allocate(dimids(ndims), dims(ndims), newdimids(ndims+1))

      ierr = nf90_inquire_variable(ncin%ncid, var_id, dimids=dimids)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'copy_variable_defmode', ncin%filename)
      is_spatial = .false.

      ! Assume all variables to be interpolated are actually dimensioned
      ! spatially

      ierr = nf90_inq_dimid(ncout%ncid, 'xDim', newdimids(1))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_dimid', .false., 'copy_variable_defmode', ncout%filename)

      ierr = nf90_inq_dimid(ncout%ncid, 'yDim', newdimids(2))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_dimid', .false., 'copy_variable_defmode', ncout%filename)

      j = 3
      do i=1, ndims
         ierr = nf90_inquire_dimension(ncin%ncid, dimids(i), dims(i))
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .false., 'copy_variable_defmode', ncin%filename)

         if (trim(dims(i)) == 'nCells' .or. trim(dims(i)) == 'nEdges' .or. trim(dims(i)) == 'nVertices') then
            is_spatial = .true.
            cycle
         else
            ierr = nf90_inq_dimid(ncout%ncid, dims(i), newdimids(j))
            if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_dimid', .false., 'copy_variable_defmode', ncout%filename)
            j = j+1
         end if
      end do

      if (.not. is_spatial) then
         write (0,*) "No spatial dimension for "//trim(var_name)//", skipping"
         !var_name = ''
         return
      end if
      
      ierr = nf90_def_var(ncout%ncid, var_name, xtype, newdimids, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_def_var', .true., 'copy_variable_defmode', ncout%filename)

      call ncout%add_var_record(var_name)

   end subroutine copy_variable_defmode

   subroutine put_variable_1dREAL(f, field, var_name)
      type(ncfile) :: f 
      real(kind=RKIND), dimension(:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id

      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_1dREAL', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_1dREAL', f%filename)
   end subroutine put_variable_1dREAL

   
   subroutine put_variable_1dINT(f, field, var_name)
      type(ncfile) :: f
      integer, dimension(:,:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id

      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_1dINT', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_1dINT', f%filename)

   end subroutine put_variable_1dINT


   subroutine put_variable_2dREAL(f, field, var_name)
      type(ncfile) :: f
      real(kind=RKIND), dimension(:,:), pointer :: field
      character(len=*) :: var_name

      integer :: i, var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_2dREAL', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_2dREAL', f%filename)
   end subroutine put_variable_2dREAL

   subroutine put_variable_3dREAL(f, field, var_name)
      type(ncfile) :: f
      real(kind=RKIND), dimension(:,:,:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_3dREAL', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_3dREAL', f%filename)
   end subroutine put_variable_3dREAL

   subroutine put_variable_4dREAL(f, field, var_name)
      type(ncfile) :: f
      real(kind=RKIND), dimension(:,:,:,:), pointer :: field
      character(len=*) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_4dINT', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_4dREAL', f%filename)
   end subroutine put_variable_4dREAL

   subroutine put_variable_2dINT(f, field, var_name)
      type(ncfile) :: f
      integer, dimension(:,:), pointer :: field
      character(len=StrKIND) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_2dINT', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_2dINT', f%filename)
   end subroutine put_variable_2dINT

   subroutine put_variable_3dINT(f, field, var_name)
      type(ncfile) :: f
      integer, dimension(:,:,:), pointer :: field
      character(len=StrKIND) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_3dINT', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_3dINT', f%filename)
   end subroutine put_variable_3dINT

   subroutine put_variable_4dINT(f, field, var_name)
      type(ncfile) :: f
      integer, dimension(:,:,:,:), pointer :: field
      character(len=StrKIND) :: var_name

      integer :: var_id
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .false., 'put_variable_4dINT', f%filename)

      ierr = nf90_put_var(f%ncid, var_id, field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_put_var', .false., 'put_variable_4dINT', f%filename)
   end subroutine put_variable_4dINT

   subroutine get_variable_1dINT(f, var_name, field)
      implicit none
      
      type(ncfile) :: f
      character(len=*), intent(in) :: var_name
      integer, dimension(:), pointer, intent(inout) :: field

      integer :: var_id, n, ierr
      integer, dimension(1) :: dim_ids

      !if (associated(field)) deallocate(field)
      
      ierr = nf90_inq_varid(f%ncid, var_name, var_id)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .true., 'get_variable_1dINT', f%filename)

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'get_variable_1dINT', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_1dINT', f%filename)

      allocate(field(n))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1/), count = (/n/))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_get_var', .true., 'get_variable_1dINT', f%filename)
      
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
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .true., 'get_variable_2dINT', f%filename)

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'get_variable_2dINT', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_2dINT', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_2dINT', f%filename)

      allocate(field(n1, n2))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1/), count = (/n1, n2/))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_get_var', .true., 'get_variable_2dINT', f%filename)
      
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
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .true., 'get_variable_3dINT', f%filename)


      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'get_variable_3dINT', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_3dINT', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_3dINT', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(3), len=n3)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_3dINT', f%filename)

      allocate(field(n1, n2, n3))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1,1/), count = (/n1, n2, n3/))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_get_var', .true., 'get_variable_3dINT', f%filename)
      
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
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .true., 'get_variable_1dREAL', f%filename)

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids, xtype=xtype)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'get_variable_1dREAL', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_1dREAL', f%filename)

      allocate(field(n))
      ierr = nf90_get_var(f%ncid, var_id, field, count = (/n/))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_get_var', .true., 'get_variable_1dREAL', f%filename)
      
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
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .true., 'get_variable_2dREAL', f%filename)

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'get_variable_2dREAL', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_2dREAL', f%filename)


      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_2dREAL', f%filename)

      allocate(field(n1, n2))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1/), count = (/n1, n2/))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_get_var', .true., 'get_variable_2dREAL', f%filename)
      
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
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_varid', .true., 'get_variable_3dREAL', f%filename)

      ierr = nf90_inquire_variable(f%ncid, var_id, dimids=dim_ids)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_variable', .true., 'get_variable_3dREAL', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(1), len=n1)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_3dREAL', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(2), len=n2)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_3dREAL', f%filename)

      ierr = nf90_inquire_dimension(f%ncid, dim_ids(3), len=n3)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .true., 'get_variable_3dREAL', f%filename)

      allocate(field(n1, n2, n3))
      ierr = nf90_get_var(f%ncid, var_id, field, start=(/1,1,1/), count = (/n1, n2, n3/))
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_get_var', .true., 'get_variable_3dREAL', f%filename)
      
   end subroutine get_variable_3dREAL

   subroutine get_attribute_REAL(f, att_name, field)
      implicit none
      
      type(ncfile) :: f 
      character(len=*), intent(in) :: att_name
      real (kind=RKIND), intent(out) :: field
      integer :: ierr
      ierr = nf90_get_att(f%ncid, NF90_GLOBAL, trim(att_name), field)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_get_att', .false., 'get_attribute_real')
   end subroutine get_attribute_REAL
   
   subroutine close_mpas_file(f)
      implicit none
      
      type(ncfile) :: f
      integer :: ierr

      ierr = nf90_close(f%ncid)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_close', .true., 'close_mpas_file')
      call f%clean()
   end subroutine

   subroutine handle_err(ierr, metname, abort, funcname, filename)
      implicit none
      integer :: ierr
      character(len=*) :: metname
      logical, optional :: abort
      character(len=*), optional :: funcname, filename
      

      write(0,*) '*********************************************************************************'
      write(0,*) 'Error in netcdf method '//trim(metname)
      write(0,*) 'ierr = ', ierr
      if (present(funcname)) write(0,*) 'Function: '//trim(funcname)
      if (present(filename)) write(0,*) 'Filename: '//trim(filename)
      if (present(abort)) then
         if(abort) write (0,*) 'Stopping Program'
      end if
      write(0,*) '*********************************************************************************'
      if (present(abort)) then
         if (abort) then
            stop
         end if
      end if

   end subroutine handle_err
      

end module mpas_file_manip 
   
