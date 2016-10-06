module utils
   use netcdf
   use mpas_file_manip

!   integer, parameter :: NN=1, WP=2
!
!   type :: interpgrid
!      integer :: mode = NN
!      integer, dimension(:,:), pointer :: cell_map, edge_map, vertex_map
!      real (kind=RKIND) :: lat_start, lat_end, lon_start, lon_end
!      real (kind=RKIND), dimension(:,:), pointer :: lats, lons
!      integer :: nx, ny
!   end type interpgrid

   contains

   subroutine create_grid(grid, f, orig_lat, orig_lon, new_lat, new_lon, birdseye)
      use mesh_rotate
      implicit none
      type(interpgrid) :: grid
      type(ncfile) :: f
      real(kind=RKIND) :: orig_lat, orig_lon, new_lat, new_lon, birdseye

      real(kind=RKIND) :: dx, dy
      integer :: i, j

      if (grid%lat_start == grid%lat_end .or. grid%lon_start == grid%lon_end .or. grid%nx == 0 .or. grid%ny == 0) then
         write (0,*) "Error: Invalid argument ''grid'' in ''create_grid''."
         stop
      end if


      allocate(grid%lats(grid%nx, grid%ny), grid%lons(grid%nx, grid%ny))
         
      dx = (grid%lon_end - grid%lon_start) / grid%nx
      dy = (grid%lat_end - grid%lat_start) / grid%ny
      do j=1, grid%ny
         grid%lats(:,j) = grid%lat_start + (j-1) * dy + dy / 2.0
      end do
      do i=1, grid%nx
         grid%lons(i,:) = grid%lon_start + (i-1) * dx + dx / 2.0
      end do

      if (needs_rotated(orig_lat, orig_lon, new_lat, new_lon, birdseye)) then
         write (0,*) "needs rotated"
         call rotate(grid, orig_lat, orig_lon, new_lat, new_lon, birdseye) 
      end if

      call create_maps_NN(grid, f)
      if (grid%mode == WP) then
         call calculate_weights_WP(grid, f)
      end if

   end subroutine create_grid

   subroutine calculate_weights_WP(grid, f)
      implicit none
      type(interpgrid) :: grid
      type(ncfile) :: f

      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex
      real(kind=RKIND), dimension(:), pointer :: latCell, lonCell, latEdge, lonEdge, latVertex, lonVertex
      real(kind=RKIND), dimension(:, :), allocatable :: vertCoords
      real(kind=RKIND), dimension(3) :: pt
      integer :: i, j, k, c, v, nVertices

      write (0,*) "  calculating remapping weights"

      ! For cell-based fields
      call get_variable_1dREAL(f, 'latCell', latCell)
      call get_variable_1dREAL(f, 'lonCell', lonCell)
      call get_variable_2dINT(f, 'cellsOnVertex', cellsOnVertex)
      allocate(vertCoords(3,3))

      allocate(grid%cell_weights(3, grid%nx, grid%ny))
      do j=1, grid%ny
      do i=1, grid%nx
         do k=1, 3
            c = cellsOnVertex(k, grid%vertex_map(i,j))
            call con_lx(latCell(c), lonCell(c), 1.0, vertCoords(1, k), vertCoords(2, k), vertCoords(3, k))
         end do
         call con_lx(grid%lats(i,j), grid%lons(i,j), 1.0, pt(1), pt(2), pt(3))
         grid%cell_weights(:,i,j) = mpas_wachspress_coordinates(3,  vertCoords, pt, .true., 1.0)
      end do
      end do
      deallocate(latCell, lonCell, vertCoords)

      ! For vertex-based fields
      call get_variable_1dREAL(f, 'latVertex', latVertex)
      call get_variable_1dREAL(f, 'lonVertex', lonVertex)
      call get_variable_1dINT(f, 'nEdgesOnCell', nEdgesOnCell)
      call get_variable_2dINT(f, 'verticesOnCell', verticesOnCell)

      allocate(grid%vertex_weights(maxval(nEdgesOnCell), grid%nx, grid%ny))
      do j=1, grid%ny
      do i=1, grid%nx
         c = grid%cell_map(i,j)
         nVertices = nEdgesOnCell(c)
         if (allocated(vertCoords)) deallocate(vertCoords)
         allocate(vertCoords(3, nVertices))
         do k=1, nVertices
            v = verticesOnCell(k, c)
            call con_lx(latVertex(v), lonVertex(v), 1.0, vertCoords(1,k), vertCoords(2,k), vertCoords(3,k))
         end do
         call con_lx(grid%lats(i,j), grid%lons(i,j), 1.0, pt(1), pt(2), pt(3))
         grid%vertex_weights(:,i,j) = mpas_wachspress_coordinates(nVertices, vertCoords, pt, .true., 1.0)
      end do
      end do
      deallocate(latVertex, lonVertex, vertCoords, verticesOnCell, nEdgesOnCell)
      
   end subroutine calculate_weights_WP

!***********************************************************************
!
!  function mpas_wachspress_coordinates
!
!> \brief Compute the barycentric Wachspress coordinates for a polygon
!> \author  Phillip Wolfram
!> \date    01/26/2015
!> \details
!>  Computes the barycentric Wachspress coordinates for a polygon with nVertices
!>  points in R3, vertCoords for a particular pointInterp with normalized radius.
!>  Follows Gillette, A., Rand, A., Bajaj, C., 2011.
!>  Error estimates for generalized barycentric interpolation.
!>  Advances in computational mathematics 37 (3), 417–439.
!>  Optimized version of mpas_wachspress_coordinates uses optional cached B_i areas
!------------------------------------------------------------------------
   function mpas_wachspress_coordinates(nVertices, vertCoords, pointInterp, on_a_sphere, sphere_radius) !{{{
      implicit none

      ! input points
      integer, intent(in) :: nVertices
      real (kind=RKIND), dimension(3, nVertices), intent(in) :: vertCoords
      real (kind=RKIND), dimension(3), intent(in) :: pointInterp
      ! output
      real (kind=RKIND), dimension(nVertices) :: mpas_wachspress_coordinates
      ! computational intermediates
      real (kind=RKIND), dimension(nVertices) :: wach       ! The wachpress area-product
      real (kind=RKIND) :: wach_total                       ! The wachpress total weight
      integer :: i, j                                       ! Loop indices
      integer :: im1, i0, ip1                               ! im1 = (i-1), i0 = i, ip1 = (i+1)

      ! triangle areas to compute wachspress coordinate
      real (kind=RKIND), dimension(nVertices) :: areaA
      real (kind=RKIND), dimension(nVertices) :: areaB

      logical, intent(in) :: on_a_sphere
      real(kind=RKIND), intent(in), optional :: sphere_radius
      real(kind=RKIND) :: radiusLocal

      if ( on_a_sphere ) then
         radiusLocal = sphere_radius
      else
         radiusLocal = 1.0_RKIND
      end if

     ! compute areas
     do i = 1, nVertices
        ! compute first area B_i
        ! get vertex indices
        im1 = mod(nVertices + i - 2, nVertices) + 1
        i0  = mod(nVertices + i - 1, nVertices) + 1
        ip1 = mod(nVertices + i    , nVertices) + 1

        ! precompute B_i areas
        ! always the same because B_i independent of xp,yp,zp
        ! (COULD CACHE AND USE RESULT FROM ARRAY FOR FURTHER OPTIMIZATION)
        areaB(i) = mpas_triangle_signed_area(vertCoords(:, im1), vertCoords(:, i0), vertCoords(:, ip1), on_a_sphere, radiusLocal)
     end do

      ! compute areas
      do i = 1, nVertices
         ! compute first area B_i
         ! get vertex indices
         im1 = mod(nVertices + i - 2, nVertices) + 1
         i0  = mod(nVertices + i - 1, nVertices) + 1
         ip1 = mod(nVertices + i    , nVertices) + 1

         ! compute A_ij areas
         ! must be computed each time
         areaA(i0) = mpas_triangle_signed_area(pointInterp, vertCoords(:, i0), vertCoords(:, ip1), on_a_sphere, radiusLocal)

         ! precomputed B_i areas, cached
      end do


      ! for each vertex compute wachpress coordinate
      do i = 1, nVertices
         wach(i) = areaB(i)
         do j = (i + 1), (i + nVertices - 2)
            i0  = mod(nVertices + j - 1, nVertices) + 1
            ! accumulate products for A_ij subareas
            wach(i) = wach(i) * areaA(i0)
         end do
      end do

      ! get summed weights for normalization
      wach_total = 0
      do i = 1, nVertices
         wach_total = wach_total + wach(i)
      end do

      ! compute lambda
      mpas_wachspress_coordinates= 0.0_RKIND
      do i = 1, nVertices
         mpas_wachspress_coordinates(i) = wach(i)/wach_total
      end do

   end function mpas_wachspress_coordinates!}}}


!***********************************************************************
!
!  routine mpas_wachspress_interpolate
!
!> \brief Interpolate using barycentric Wachspress coordinates
!> \author  Phillip Wolfram
!> \date    03/27/2015
!> \details
!>  Interpolate using the barycentric Wachspress coordinates for a polygon with nVertices
!>  having values phi.
!------------------------------------------------------------------------
   real (kind=RKIND) function mpas_wachspress_interpolate(lambda, phi) !{{{
      implicit none

      ! input points
      real (kind=RKIND), dimension(:), intent(in) :: lambda   !< Input: Wachspress coordinate / weight
      real (kind=RKIND), dimension(:), intent(in) :: phi      !< Input: values at lambda weights
      ! output for function
      !real (kind=RKIND), intent(out) :: mpas_wachspress_interpolate

      mpas_wachspress_interpolate = sum(phi * lambda)

   end function mpas_wachspress_interpolate! }}}

!***********************************************************************
!
!  routine mpas_triangle_signed_area
!
!> \brief   Calculates area of a triangle, whether on a sphere or a plane.
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle whether on a sphere or a plane.
!>  Note this does not handle triangles spanning planar periodic meshes because mpas_triangle_signed_area_plane does not!
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_triangle_signed_area(a, b, c, on_a_sphere, radius)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle for which to get the area
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, intent(in) :: on_a_sphere
      real(kind=RKIND), dimension(3) :: normalvec
      real(kind=RKIND), intent(in) :: radius

      if (on_a_sphere) then
         mpas_triangle_signed_area = mpas_triangle_signed_area_sphere(a, b, c, radius)
      else
         normalvec = (/ 0, 0, 1 /)
!         mpas_triangle_signed_area = mpas_triangle_signed_area_plane(a, b, c, normalvec)
      endif
   end function mpas_triangle_signed_area !}}}
   
!***********************************************************************
!
!  routine mpas_triangle_signed_area_plane
!
!> \brief   Calculates signed area of a triangle in a plane
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle in a plane.
!>  Uses cross product.  Signed area will be positive if the vertices are oriented counterclockwise.
!>  Note this does not handle triangles spanning periodic meshes!
!-----------------------------------------------------------------------
!   real(kind=RKIND) function mpas_triangle_signed_area_plane(a, b, c, normalvec)!{{{
!      !-----------------------------------------------------------------
!      ! input variables
!      !-----------------------------------------------------------------
!      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle for which to calculate the area
!      real(kind=RKIND), dimension(3), intent(in) :: normalvec  !< Input: 3d vector indicating the normal direction for the plane for assigning a sign to the area
!      !-----------------------------------------------------------------
!      ! local variables
!      !-----------------------------------------------------------------
!      real(kind=RKIND), dimension(3) :: ab, ac, crossprod, triangleNormal
!
!      ab = b - a
!      ac = c - a
!      call mpas_cross_product_in_r3(ab, ac, crossprod)
!      if (mpas_vec_mag_in_r3(crossprod) == 0.0_RKIND) then
!         mpas_triangle_signed_area_plane = 0.0_RKIND
!      else
!         triangleNormal = crossprod / mpas_vec_mag_in_r3(crossprod)
!         mpas_triangle_signed_area_plane = 0.5_RKIND * (mpas_vec_mag_in_r3(crossprod)) *  &
!              sum(triangleNormal * normalvec)
!      endif
!   end function mpas_triangle_signed_area_plane !}}}


!***********************************************************************
!
!  routine mpas_triangle_signed_area_sphere
!
!> \brief   Calculates area of a triangle on a sphere
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine calculates the area of a triangle on the surface of a sphere.
!>  Uses the spherical analog of Heron's formula.
!>  Copied from mesh generator.  A CCW winding angle is positive.
!-----------------------------------------------------------------------
   real(kind=RKIND) function mpas_triangle_signed_area_sphere(a, b, c, radius)!{{{
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(3), intent(in) :: a, b, c  !< Input: 3d (x,y,z) points forming the triangle in which to calculate the bary weights
      real(kind=RKIND), intent(in) :: radius  !< sphere radius
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND) :: ab, bc, ca, semiperim, tanqe
      real(kind=RKIND), dimension(3) :: ablen, aclen, Dlen

      ab = sphere_distance_x(a(1), a(2), a(3), b(1), b(2), b(3))/radius
      bc = sphere_distance_x(b(1), b(2), b(3), c(1), c(2), c(3))/radius
      ca = sphere_distance_x(c(1), c(2), c(3), a(1), a(2), a(3))/radius
      semiperim = 0.5 * (ab + bc + ca)

      tanqe = sqrt(max(0.0_RKIND,tan(0.5_RKIND * semiperim) * tan(0.5_RKIND * (semiperim - ab)) &
                   * tan(0.5_RKIND * (semiperim - bc)) * tan(0.5_RKIND * (semiperim - ca))))

      mpas_triangle_signed_area_sphere = 4.0_RKIND * radius * radius * atan(tanqe)

      ! computing correct signs (in similar fashion to mpas_sphere_angle)
      ablen(1) = b(1) - a(1)
      ablen(2) = b(2) - a(2)
      ablen(3) = b(3) - a(3)

      aclen(1) = c(1) - a(1)
      aclen(2) = c(2) - a(2)
      aclen(3) = c(3) - a(3)

      dlen(1) =   (ablen(2) * aclen(3)) - (ablen(3) * aclen(2))
      dlen(2) = -((ablen(1) * aclen(3)) - (ablen(3) * aclen(1)))
      dlen(3) =   (ablen(1) * aclen(2)) - (ablen(2) * aclen(1))

      if ((Dlen(1)*a(1) + Dlen(2)*a(2) + Dlen(3)*a(3)) < 0.0) then
        mpas_triangle_signed_area_sphere = -mpas_triangle_signed_area_sphere
      end if

   end function mpas_triangle_signed_area_sphere !}}}

!***********************************************************************
!
!  routine mpas_point_in_polygon
!
!> \brief   Hit test to determine if a point is inside of a polygon
!> \author  Matthew Hoffman
!> \date    13 January 2015
!> \details
!>  This routine determines if a point is inside of a polygon.
!>  This is difficult because floating point arithmetic prevents a precise
!>  determination.  A tolerance is used to allow the point to be within the
!>  the polygon within some tolerance.  This means it is possible for a point
!>  to be identified to be within multiple polygons.  However, it avoids the
!>  situation where a point on a edge could be 'orphaned' - determined
!>  to not belong to *any* polygons.
!-----------------------------------------------------------------------
!   logical function mpas_point_in_polygon(point, polygonVertices, on_a_sphere)!{{{
!      !-----------------------------------------------------------------
!      ! input variables
!      !-----------------------------------------------------------------
!      real(kind=RKIND), dimension(3), intent(in) :: point  !< Input: 3d (x,y,z) point
!      real(kind=RKIND), dimension(:,:), intent(in) :: polygonVertices  !< Input: 3d (x,y,z) points forming the polygon to test, second dimension should be 3
!      logical, intent(in) :: on_a_sphere  !< Input: If on a sphere
!      !-----------------------------------------------------------------
!      ! local variables
!      !-----------------------------------------------------------------
!      real(kind=RKIND), dimension(3) :: normal_vector, crossprod, vec1, vec2
!      integer :: polygonDegree, i
!      integer, dimension(:), allocatable :: vertexNeighborFwd
!      real(kind=RKIND), parameter :: eps = 1.0e-12_RKIND
!
!      if (on_a_sphere) then
!         normal_vector = point
!      else
!         normal_vector = (/ 0.0_RKIND, 0.0_RKIND, 1.0_RKIND /)
!      endif
!
!      polygonDegree = size(polygonVertices, 1)
!      allocate(vertexNeighborFwd(polygonDegree))
!      vertexNeighborFwd  = (/ (i+1, i = 1, polygonDegree) /)
!      vertexNeighborFwd(polygonDegree) = 1
!
!      mpas_point_in_polygon = .true.
!      do i = 1, polygonDegree
!         vec1 = polygonVertices(vertexNeighborFwd(i),:) - polygonVertices(i,:)
!         vec2 = point - polygonVertices(i,:)
!         call mpas_cross_product_in_r3(vec1, vec2, crossprod)
!         if (sum(crossprod * normal_vector) < (0.0_RKIND - eps)) then
!            mpas_point_in_polygon = .false.
!            exit  ! If the point is ouside one of the edges, then we need not look further.
!         endif
!      enddo
!
!      deallocate(vertexNeighborFwd)
!
!   end function mpas_point_in_polygon !}}}
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************


   subroutine create_maps_NN(grid, f)
      implicit none
      type(interpgrid) :: grid
      type(ncfile) :: f

      real(kind=RKIND), dimension(:), pointer :: latElem, lonElem
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell, verticesOnCell
      real(kind=RKIND) :: temp, dist, mindist, lat_pt, lon_pt, dx, dy
      integer :: v, iCell, k, i, j, e, nearest_cell

      allocate(grid%cell_map(grid%nx, grid%ny))

      call get_variable_1dREAL(f, 'latCell', latElem)
      call get_variable_1dREAL(f, 'lonCell', lonElem)
      call get_variable_2dINT(f, 'cellsOnCell', cellsOnCell)
      call get_variable_1dINT(f, 'nEdgesOnCell', nEdgesOnCell)

      write (0,*) "  making cell_map"
      nearest_cell = 1
      do j=1, grid%ny
      do i=1, grid%nx
         nearest_cell = nearest_cell_path(grid%lats(i,j), grid%lons(i,j), nearest_cell, &
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
      do i=1, grid%nx
         iCell = grid%cell_map(i, j)
         v = verticesOnCell(1, iCell)
         mindist = sphere_distance(latElem(v), lonElem(v), grid%lats(i,j), grid%lons(i,j), 1.0)
         do k = 2, nEdgesOnCell(iCell)  
            temp = verticesOnCell(k, iCell)
            dist = sphere_distance(latElem(temp), lonElem(temp), grid%lats(i,j), grid%lons(i,j), 1.0)
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
         mindist = sphere_distance(latElem(v), lonElem(v), grid%lats(i,j), grid%lons(i,j), 1.0)
         do k = 2, nEdgesOnCell(iCell)  
            temp = edgesOnCell(k, iCell)
            dist = sphere_distance(latElem(temp), lonElem(temp), grid%lats(i,j), grid%lons(i,j), 1.0)
            if (dist < mindist) then
               mindist = dist
               e = temp
            end if
         end do
         grid%edge_map(i, j) = e
      end do
      end do

      deallocate(lonElem, latElem, edgesOnCell, nEdgesOnCell)

   end subroutine create_maps_NN

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
      integer :: i, ierr, varid, xdimid, ydimid

      do i=1, size(vars)
         if (.not. fout%contains_elem(VAR, trim(vars(i))) .and. fin%contains_elem(VAR, trim(vars(i)))) then
           call copy_variable_defmode(fin, fout, vars(i))
           call put_att_str(fout, vars(i), 'coordinates', 'lat_pt lon_pt')
         end if
      end do 

      ierr = nf90_inq_dimid(fout%ncid, 'xDim', xdimid)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_dimid', .false., 'define_variables_io')

      ierr = nf90_inq_dimid(fout%ncid, 'yDim', ydimid)
      if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inq_dimid', .false., 'define_variables_io')

      if (.not. fout%contains_elem(VAR, 'lat_pt')) then
         ierr = nf90_def_var(fout%ncid, 'lat_pt', NF90_REAL, (/xdimid, ydimid/), varid)
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_def_var', .false., 'define_variables_io')
         call fout%add_var_record('lat_pt')
         call put_att_str(fout, 'lat_pt', 'units', 'degree_east')
      end if
      if (.not. fout%contains_elem(VAR, 'lon_pt')) then
         ierr = nf90_def_var(fout%ncid, 'lon_pt', NF90_REAL, (/xdimid, ydimid/), varid)
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_def_var', .false., 'define_variables_io')
         call fout%add_var_record('lon_pt')
         call put_att_str(fout, 'lon_pt', 'units', 'degree_north')
      end if
          
   end subroutine define_variables_io

   subroutine copy_data(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name
      integer :: ierr, var_id, xtype, ndims

      integer(kind=8) :: t1, t2, rate

      write (0,*) "Copying "//trim(var_name)
      call system_clock(t1, rate)

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
      call system_clock(t2)
      write (0,*) "Time to copy"//trim(var_name)//":", real(t2-t1) / real(rate)

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
         call put_variable_1dINT(fout, field, var_name)
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
         write (0,*) "Variable "//trim(var_name)//" not dimensioned spatially"
         call put_variable_2dINT(fout, field, var_name)
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
         write (0,*) "Variable "//trim(var_name)//" not dimensioned spatially"
         call put_variable_3dINT(fout, field, var_name)
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

      integer, dimension(:,:), pointer :: map, elOnElem
      real(kind=RKIND), dimension(:,:,:), pointer :: weights
      real(kind=RKIND), dimension(:), pointer :: field, vals
      real(kind=RKIND), dimension(:,:), pointer :: newfield
      integer :: i, j, k, elem
      integer, dimension(3) :: lens

      call get_variable_1dREAL(fin, var_name, field)
      
      if (size(field) == fin%nCells) then
         if(grid%mode == NN) then
            map => grid%cell_map
         else if (grid%mode == WP) then
            weights => grid%cell_weights
            map => grid%vertex_map
            call get_variable_2dINT(fin, 'cellsOnVertex', elOnElem)
         end if
      else if (size(field) == fin%nEdges) then
         map => grid%edge_map 
         weights => grid%edge_weights
      else if (size(field) == fin%nVertices) then
         if(grid%mode == NN) then
            map => grid%vertex_map
         else if (grid%mode == WP) then
            weights => grid%vertex_weights
            map => grid%cell_map
            call get_variable_2dINT(fin, 'verticesOnCell', elOnElem)
         end if
      else
         write (0,*) "Variable "//trim(var_name)//" not dimensioned spatially"
         call put_variable_1dREAL(fout, field, var_name)
         return
      end if

      allocate(newfield(grid%nx, grid%ny))

      if (grid%mode == NN) then
         do j=1, grid%ny
         do i=1, grid%nx
            newfield(i, j) = field(map(i,j))
         end do
         end do
      else if (grid%mode == WP) then
         lens = shape(weights)
         allocate(vals(lens(1)))
         do j=1, grid%ny
         do i=1, grid%nx
            elem = map(i,j)
            do k=1, lens(1)
               vals(k) = field(elOnElem(k, elem))
            end do 
            newfield(i, j) = mpas_wachspress_interpolate(weights(:,i,j), vals)
         end do
         end do
      end if

      call put_variable_2dREAL(fout, newfield, var_name)
   end subroutine copy_data_1dREAL

   subroutine copy_data_2dREAL(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map, elOnElem
      real(kind=RKIND), dimension(:,:,:), pointer :: weights
      real(kind=RKIND), dimension(:,:), pointer :: field
      real(kind=RKIND), dimension(:), pointer :: vals
      real(kind=RKIND), dimension(:,:,:), pointer :: newfield
      integer :: i, j, k, t, xtype, ierr, ndims, var_id, elem
      integer :: iSlice, nSlices, slice_dim_len
      integer, dimension(3) :: lens
      integer, dimension(2) :: dimids = 0, dimlens = 0
      character(len=StrKIND) :: dim_name
      logical :: has_time

      
      ierr = nf90_inq_varid(fin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         call handle_err(ierr, 'nf90_inq_varid', .false., 'copy_data_2dREAL', fin%filename) 
         return
      end if

      ierr = nf90_inquire_variable(fin%ncid, var_id, xtype=xtype, ndims=ndims, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         call handle_err(ierr, 'nf90_inquire_variable', .false., 'copy_data_2dREAL', fin%filename) 
         return
      end if

      do i=1, 2
         ierr = nf90_inquire_dimension(fin%ncid, dimids(i), name=dim_name, len=dimlens(i))
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .false., 'copy_data_2dREAL', fin%filename) 
      end do

      
      ! Determine what maps and remapping weights need to be used based on the
      ! variable's dimensions
      if (dimlens(1) == fin%nCells .or. dimlens(2) == fin%nCells) then
         if(grid%mode == NN) then
            map => grid%cell_map
         else if (grid%mode == WP) then
            weights => grid%cell_weights
            map => grid%vertex_map
            call get_variable_2dINT(fin, 'cellsOnVertex', elOnElem)
         end if
      else if (dimlens(1) == fin%nEdges .or. dimlens(2) == fin%nEdges) then
         map => grid%edge_map 
         weights => grid%edge_weights
      else if (dimlens(1) == fin%nVertices .or. dimlens(2) == fin%nVertices) then
         if(grid%mode == NN) then
            map => grid%vertex_map
         else if (grid%mode == WP) then
            weights => grid%vertex_weights
            map => grid%cell_map
            call get_variable_2dINT(fin, 'verticesOnCell', elOnElem)
         end if
      else
         ! If the variable is not dimensioned spatially, simply copy it over
         ! as-is
         call put_variable_2dREAL(fout, field, var_name)
         return
      end if


      ! If the least-rapidly-varying dimension of a field is not the spatial
      ! dimension for that field, then it must be the time dimension
      ! (NC_UNLIMITED)
      if(dimlens(1) == fin%nCells  .or. dimlens(1) == fin%nEdges .or. dimlens(1) == fin%nVertices) then
         has_time = .true.
         slice_dim_len = dimlens(2)
      else 
         slice_dim_len = dimlens(1)
         has_time = .false.
      end if
      

      ! If the field is getting so large that we don't want to hold the whole
      ! thing in memory, then we can interpolate in slices (each slice being
      ! along the dimensions that are not the spatial dimension
      if (dimlens(1) * dimlens(2) > MAX_CHUNK_SIZE) then
         nSlices = slice_dim_len
      else
         nSlices = 1
      end if
      

      ! Perform interpolation in 'slices' of the non-spatial dimension. If the
      ! field is small enough to do in 1 slice, then nSlices == 1
      do iSlice = 1, nSlices
         call get_variable_2dREAL(fin, var_name, field, start=(/merge(1, iSlice, has_time), merge(iSlice, 1, has_time)/), &
                                  cnt=(/merge(dimlens(1), slice_dim_len/nSlices, has_time), merge(slice_dim_len/nSlices, dimlens(2), has_time)/))
         select case (grid%mode)
         case (NN)
            if (has_time) then
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(2)))
               do j=1, grid%ny
               do i=1, grid%nx
                  newfield(i, j, iSlice:iSlice+slice_dim_len/nSlices-1) = field(map(i, j), :)
               end do   
               end do   
            else
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(1)))
               do j=1, grid%ny
               do i=1, grid%nx
                  newfield(i, j, iSlice:iSlice+slice_dim_len/nSlices-1) = field(:, map(i, j))
               end do
               end do   
            end if
         case (WP)
            lens = shape(weights)
            allocate(vals(lens(1)))
            if (has_time) then
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(2)))
               do j=1, grid%ny
               do i=1, grid%nx
                  elem = map(i,j)
                  do t=1, size(field(1,:)) 
                  do k=1, lens(1)
                     vals(k) = field(elOnElem(k, elem), t)
                  end do 
                  newfield(i, j, iSlice+t-1) = mpas_wachspress_interpolate(weights(:,i,j), vals)
                  end do
               end do
               end do
            else 
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(1)))
               do j=1, grid%ny
               do i=1, grid%nx
                  elem = map(i,j)
                  do t=1, size(field(:,1)) 
                  do k=1, lens(1)
                     vals(k) = field(t,elOnElem(k, elem))
                  end do 
                  newfield(i, j, iSlice+t-1) = mpas_wachspress_interpolate(weights(:,i,j), vals)
                  end do
               end do
               end do
            end if
         end select

      end do

      call put_variable_3dREAL(fout, newfield, var_name)

      deallocate(field, newfield)

   end subroutine copy_data_2dREAL

   subroutine copy_data_3dREAL(fin, fout, var_name, grid)
      implicit none
      type(ncfile) :: fin, fout
      type(interpgrid) :: grid
      character(len=*) :: var_name

      integer, dimension(:,:), pointer :: map, elOnElem
      real(kind=RKIND), dimension(:,:,:), pointer :: weights
      real(kind=RKIND), dimension(:,:,:), pointer :: field
      real(kind=RKIND), dimension(:,:,:,:), pointer :: newfield
      real(kind=RKIND), dimension(:), pointer :: vals
      integer, dimension(3) :: lens
      integer :: i, j, k, u, t, elem, ierr, ndims, var_id, xtype
      integer :: iSlice, jSlice, niSlices, njSlices, islice_dim_len, jslice_dim_len
      integer, dimension(3) :: dimids = 0, dimlens = 0
      character(len=StrKIND) :: dim_name
      logical :: has_time


      integer(kind=8) :: t1, t2, t3, rate
      real(kind=8) :: read_time


      ! Get variable info from input file 
      ierr = nf90_inq_varid(fin%ncid, trim(var_name), var_id)
      if (ierr /= NF90_NOERR) then
         call handle_err(ierr, 'nf90_inq_varid', .false., 'copy_data_3dREAL', fin%filename) 
         return
      end if

      ierr = nf90_inquire_variable(fin%ncid, var_id, xtype=xtype, ndims=ndims, dimids=dimids)
      if (ierr /= NF90_NOERR) then
         call handle_err(ierr, 'nf90_inquire_variable', .false., 'copy_data_3dREAL', fin%filename) 
         return
      end if

      do i=1, 3
         ierr = nf90_inquire_dimension(fin%ncid, dimids(i), name=dim_name, len=dimlens(i))
         if (ierr /= NF90_NOERR) call handle_err(ierr, 'nf90_inquire_dimension', .false., 'copy_data_3dREAL', fin%filename) 
      end do

      ! Determine what maps and remapping weights need to be used based on the
      ! variable's dimensions
      if (dimlens(3) == fin%nCells .or. dimlens(2) == fin%nCells) then
         if(grid%mode == NN) then
            map => grid%cell_map
         else if (grid%mode == WP) then
            weights => grid%cell_weights
            map => grid%vertex_map
            call get_variable_2dINT(fin, 'cellsOnVertex', elOnElem)
         end if
      else if (dimlens(3) == fin%nEdges .or. dimlens(2) == fin%nEdges) then
         map => grid%edge_map 
         weights => grid%edge_weights
      else if (dimlens(3) == fin%nVertices .or. dimlens(2) == fin%nVertices) then
         if(grid%mode == NN) then
            map => grid%vertex_map
         else if (grid%mode == WP) then
            weights => grid%vertex_weights
            map => grid%cell_map
            call get_variable_2dINT(fin, 'verticesOnCell', elOnElem)
         end if
      else
         ! If the variable is not dimensioned spatially, simply copy it over
         ! as-is
         call put_variable_3dREAL(fout, field, var_name)
         return
      end if


      ! If the least-rapidly-varying dimension of a field is not the spatial
      ! dimension for that field, then it must be the time dimension
      ! (NC_UNLIMITED)
      if(dimlens(2) == fin%nCells  .or. dimlens(2) == fin%nEdges .or. dimlens(2) == fin%nVertices) then
         has_time = .true.
         islice_dim_len = dimlens(3)
      else 
         has_time = .false.
         islice_dim_len = dimlens(2)
      end if
      jslice_dim_len = dimlens(1)


      ! If the field is getting so large that we don't want to hold the whole
      ! thing in memory, then we can interpolate in slices (each slice being
      ! along the dimensions that are not the spatial dimension
      if (dimlens(1) * dimlens(2) * dimlens(3) > MAX_CHUNK_SIZE) then
         niSlices = islice_dim_len
         njSlices = jslice_dim_len
      else
         njSlices = 1
         niSlices = 1
      end if

      
      ! For each slice of the field, interpolate it based on a) 'has_time' and
      ! b) 'mode'
      nullify(field)

      call system_clock(t1, rate)
      read_time = 0.0
      do iSlice = 1, niSlices
      do jSlice = 1, njSlices
         
         call system_clock(t2)
         call get_variable_3dREAL(fin, var_name, field, start=(/jSlice, merge(1, iSlice, has_time), merge(iSlice, 1, has_time)/), &
                                  cnt=(/jslice_dim_len/njSlices, merge(dimlens(2), islice_dim_len/niSlices, has_time), merge(islice_dim_len/niSlices, dimlens(3), has_time)/))
         call system_clock(t3)
         read_time = read_time + real(t3-t2) / real(rate)

         select case (grid%mode)

         ! Nearest-Neighbor interpolation mode
         case (NN)
            ! If has_time, then the spatial dimension is the
            ! second-to-least-rapidly-varying dimension (by MPAS convention),
            ! otherwise it is the least-rapidly varying.
            if (has_time) then
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(3)))
               do j=1, grid%ny
               do i=1, grid%nx
                  newfield(i, j, jSlice:jSlice+jslice_dim_len/njSlices-1, iSlice:iSlice+islice_dim_len/niSlices-1) = field(:, map(i, j), :)
               end do   
               end do   
            else
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(2)))
               do j=1, grid%ny
               do i=1, grid%nx
                  newfield(i, j, jSlice:jSlice+jslice_dim_len/njSlices-1, iSlice:iSlice+islice_dim_len/niSlices-1) = field(:, :, map(i, j))
               end do
               end do   
            end if


         ! Wachspress weighted interpolation mode
         case (WP)
            lens = shape(weights)
            if(.not. associated(vals)) allocate(vals(lens(1)))
            if (has_time) then
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(3)))
               ! For every point to be interpolated to...
               !write (0,*) shape(field)
               do j=1, grid%ny
               do i=1, grid%nx
                  elem = map(i,j)
                  ! And for every time/vert-level/soil-level/whatever...
                  do t=1, size(field(1, 1, :))
                  do u=1, size(field(:, 1, 1)) 
                     ! Determine the values that must be weighted, and call the
                     ! weighting function.
                     ! e.g. for a cell-based field, 'elem' is the nearest
                     ! vertex, elOnElem is cells on vertex 
                     do k=1, lens(1)
                        vals(k) = field(u, elOnElem(k, elem), t)
                     end do 
                     newfield(i, j, jSlice+u-1, iSlice+t-1) = mpas_wachspress_interpolate(weights(:,i,j), vals)
                  end do
                  end do
               end do
               end do
            else 
               ! Same procedure as above, but using a different dimension of the
               ! field. 
               if(.not. associated(newfield)) allocate(newfield(grid%nx, grid%ny, dimlens(1), dimlens(2)))
               do j=1, grid%ny
               do i=1, grid%nx
                  elem = map(i,j)
                  do t=1, size(field(1, :, 1)) 
                  do u=1, size(field(:, 1, 1)) 
                     do k=1, lens(1)
                        vals(k) = field(u,t,elOnElem(k, elem))
                     end do 
                     newfield(i, j, iSlice+u-1, jSlice+t-1) = mpas_wachspress_interpolate(weights(:,i,j), vals)
                  end do
                  end do
               end do
               end do
            end if
         end select

      
      ! End slices
      end do
      end do

      call system_clock(t2)
      write (0,*) "  interpolation time:", real(t2-t1) / real(rate)
      write (0,*) "  read time:", read_time

      call put_variable_4dREAL(fout, newfield, var_name)

      call system_clock(t3)
      write (0,*) "  write time:", real(t3-t2) / real(rate)

      deallocate(field, newfield)

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

   real(kind=RKIND) function sphere_distance_x(x1, y1, z1, x2, y2, z2)
      implicit none
      real(kind=RKIND), intent(in) :: x1, y1, z1, x2, y2, z2
      real(kind=RKIND) :: radius

      radius = sqrt(x1**2 + y1**2 + z1**2)
      sphere_distance_x = acos(x1*x2 + y1*y2 + z1*z2) / radius
   end function sphere_distance_x

   subroutine con_lx(lat, lon, radius, x, y, z)
      implicit none

      real (kind=RKIND), intent(in) :: radius, lat, lon
      real (kind=RKIND), intent(out) :: x, y, z

      z = radius * sin(lat)
      x = radius * cos(lon) * cos(lat)
      y = radius * sin(lon) * cos(lat)
   end subroutine con_lx

end module utils
     
