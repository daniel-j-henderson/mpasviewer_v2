module params
	use netcdf
	public
		
		!This module is used to hold a bunch of variables that are (or at least were at one point) useful globally. 
		
		
		character (len=100) :: newFilename, filename, filename2 ! Names of NetCDF files

		integer :: i, j, k, moved
		integer, parameter :: MAX_VARIABLES = 20 !Max number of desired variables you can look at at one time
		integer, parameter :: NUM_VALID_MPAS_DIMS = 7, MAX_MPAS_DIMENSION = 3 !Number of supported MPAS dimensions and maximum number of dimensions per variable
		integer :: ierr       ! Return error code from NetCDF calls
		integer :: ncid, ncidNew, ncid2      ! Handle to the NetCDF files
		integer :: xID, yID, zID, nCoCID, CoCID, cellmin, terID, nCellsID, maxEdgesID, tdimID, xdimID, ydimID, TimeID, VoCID, nVertID, lonsID, latsID     ! ID of the variables needed
		integer :: xVertID, yVertID, zVertID, nEdgesID, xEdgeID, yEdgeID, zEdgeID, EoVID, nslID, nvlID, nvlP1ID
		integer, dimension(NUM_VALID_MPAS_DIMS-1) :: gridDimIDs !x, y, t, vert, soil, vertP1
		integer, dimension(NUM_VALID_MPAS_DIMS) :: meshDimIDRef, dimSizes !cells, vertices, edges, time, vert, soil, vertP1
		integer :: nCells, maxEdges, nearestOne, nearestVert, elapsedTime, nMeshVars, secondEntry, nVertices, nEdges, nearestEdge, nVertLevels = 0, nSoilLevels = 0, nVertLevelsP1 = 0
		character (len = NF90_MAX_NAME), dimension(:), allocatable :: desiredMeshVars !Array of variable names
	
		integer, dimension(:), allocatable :: nCellsOnCell, gridVarIDs, meshVarIDs, nDims, meshVarType !Mesh Parameters to be extracted
		integer, dimension(:,:), allocatable :: cellsOnCell, meshDimIDs, verticesOnCell, edgesOnVertex !Mesh Parameters to be extracted

	
		real :: r2d, radius, latSpacing, lonSpacing, x_search, y_search, z_search, d, dn
		real, allocatable, dimension(:,:,:) :: grid
		real, dimension(:), allocatable :: MeshX, MeshY, MeshZ, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, latCell, lonCell !Mesh Parameters to be extracted
	
		!Grid Parameters
		integer :: gridH = 100, gridW = 200
		real :: pi = 3.141592653589793, gridHmin = -90.0, gridHmax = 90.0, gridWmin = 0.0, gridWmax = 360.0 !let's do these in degrees
		
		
		contains
		!Sets up the grid and its parameters (converting to radians and such)
		subroutine setup()
		r2d = 180.0 / pi
		gridHmin = gridHmin / r2d
		gridHmax = gridHmax / r2d
		gridWmin = gridWmin / r2d
		gridWmax = gridWmax / r2d
		latSpacing = (gridHmax - gridHmin) / (gridH - 1)
		lonSpacing = (gridWmax - gridWmin) / (gridW - 1)
		
		ierr = nf90_get_att(ncid, NF90_GLOBAL, 'sphere_radius', radius)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error, there is no sphere_radius attribute in your info file'//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		end subroutine setup
		
		function size_of(stringArray)
			implicit none
		
			character(len=*), dimension(:), intent(in) :: stringArray
			integer :: size_of
			integer :: x, c
		
			c = 0
			do x=1,size(stringArray)
				if (len_trim(stringArray(x)) > 0) then
					c = c + 1
				end if
			end do

			size_of = c
			return
		end function size_of
	
		! Returns true if the dimension id refers to a spatial variable, else false
		function is_spatial(ID)
			implicit none
			integer, intent(in) :: ID
			logical :: is_spatial
			if (ID == nCellsID .or. ID == nVertID .or. ID == nEdgesID) then
				is_spatial = .true.
			else
				is_spatial = .false.
			end if
			return
		end function is_spatial
		
end module params
