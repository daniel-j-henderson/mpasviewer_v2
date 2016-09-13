!module a
module inputprocessing
	
	use netcdf
	use params
	
	
	contains
	
	
	subroutine open_input(filename, filename2)
		implicit none
		character(len=*), intent(in) :: filename
		character(len=*), optional, intent(in) :: filename2
		integer :: l=0
                real :: start, finish
		
		
		ierr = nf90_open(filename, NF90_NOWRITE, ncid)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error opening NetCDF file '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		if (filename2 /= filename)then
			ierr = nf90_open(filename2, NF90_NOWRITE, ncid2)
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error opening NetCDF file '//filename2
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if
		else 
			ncid2 = ncid
		end if
	
		!IDs section
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!$OMP SECTIONS
		!$OMP SECTION
		ierr = nf90_inq_varid(ncid, 'xCell', xID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of xCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'yCell', yID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of yCell in'//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'zCell', zID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of zCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'latCell', latsID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of zCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'lonCell', lonsID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of zCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'nEdgesOnCell', nCoCID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of nEdgesOnCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'cellsOnCell', CoCID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of cellsOnCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		!$OMP SECTION
		ierr = nf90_inq_varid(ncid, 'verticesOnCell', VoCID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID of verticesOnCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inq_varid(ncid, 'xVertex', xVertID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID xVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inq_varid(ncid, 'yVertex', yVertID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID yVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'zVertex', zVertID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID zVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inq_varid(ncid, 'xEdge', xEdgeID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID xEdge in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		!$OMP SECTION
		ierr = nf90_inq_varid(ncid, 'yEdge', yEdgeID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID yEdge in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'zEdge', zEdgeID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID zEdge in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_varid(ncid, 'edgesOnVertex', EoVID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring varID edgesOnVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_inq_dimid(ncid, 'nCells', nCellsID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nCells dimid in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		end if
		
		ierr = nf90_inq_dimid(ncid, 'nVertLevels', nvlID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nVertLevels dimid in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		end if
		!$OMP SECTION
		ierr = nf90_inq_dimid(ncid, 'nSoilLevels', nslID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nSoilLevels dimid in'//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
		end if
        ierr = nf90_inq_dimid(ncid, 'nVertLevelsP1', nvlP1ID)
        if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring nVertLevelsP1 dimid in'//filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
        end if
		
		ierr = nf90_inq_dimid(ncid, 'maxEdges', maxEdgesID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring dimID maxEdges in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inq_dimid(ncid2, 'Time', TimeID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring time dim ID in '//filename2
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inq_dimid(ncid, 'nVertices', nVertID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nVertices dimid in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inq_dimid(ncid, 'nEdges', nEdgesID)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nEdges dimid in'//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		!$OMP END SECTIONS
	
		!Inquire Section
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ierr = nf90_inquire_dimension(ncid, nCellsID, len=nCells)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring dimension nCells in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		end if
	
		ierr = nf90_inquire_dimension(ncid, nvlID, len=nVertLevels)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring dimension nVertLevels in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
		end if
	
		ierr = nf90_inquire_dimension(ncid, nslID, len=nSoilLevels)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring dimension nSoilLevels in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
		end if
		ierr = nf90_inquire_dimension(ncid, nvlP1ID, len=nVertLevelsP1)
        if (ierr /= NF90_NOERR) then
            write(0,*) '*********************************************************************************'
            write(0,*) 'Error inquiring dimension nSoilLevels in '//filename
            write(0,*) 'ierr = ', ierr
            write(0,*) '*********************************************************************************'
        end if
		ierr = nf90_inquire_dimension(ncid, maxEdgesID, len=maxEdges)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring dimension maxEdges in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inquire_dimension(ncid2, TimeID, len=elapsedTime)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring dimension Time in '//filename2
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inquire_dimension(ncid, nVertID, len=nVertices)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nVertices dimension in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		ierr = nf90_inquire_dimension(ncid, nEdgesID, len=nEdges)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring dimension nEdges in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
	
		
		
		
		!Get Mesh Info Data Section
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    call cpu_time(start)
		!$OMP SECTIONS
		!$OMP SECTION
		allocate(MeshX(nCells))
		ierr = nf90_get_var(ncid, xID, MeshX, count = (/nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable xCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		allocate(MeshY(nCells))
		ierr = nf90_get_var(ncid, yID, MeshY, count = (/nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable yCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		!$OMP SECTION
		allocate(MeshZ(nCells))
		ierr = nf90_get_var(ncid, zID, MeshZ, count = (/nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable zCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		allocate(latCell(nCells), lonCell(nCells))
		ierr = nf90_get_var(ncid, latsID, latCell, count = (/nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable zCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		ierr = nf90_get_var(ncid, lonsID, lonCell, count = (/nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable zCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		allocate(nCellsOnCell(nCells))
		ierr = nf90_get_var(ncid, nCoCID, nCellsOnCell, count = (/nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable nEdgesOnCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		!$OMP SECTION
		allocate(cellsOnCell(maxEdges, nCells))
		ierr = nf90_get_var(ncid, CoCID, cellsOnCell, count = (/maxEdges, nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable cellsOnCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		allocate(verticesOnCell(maxEdges, nCells))
		ierr = nf90_get_var(ncid, VoCID, verticesOnCell, count = (/maxEdges, nCells/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable verticesOnCell in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if

		!$OMP SECTION
		allocate(xVertex(nVertices))
		ierr = nf90_get_var(ncid, xVertID, xVertex, count = (/nVertices/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable xVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if

		allocate(yVertex(nVertices))
		ierr = nf90_get_var(ncid, yVertID, yVertex, count = (/nVertices/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable yVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		!$OMP SECTION
		allocate(zVertex(nVertices))
		ierr = nf90_get_var(ncid, zVertID, zVertex, count = (/nVertices/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable zVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
	

		allocate(xEdge(nEdges))
		ierr = nf90_get_var(ncid, xEdgeID, xEdge, count = (/nEdges/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable xEdge in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
	
	
		!$OMP SECTION
		allocate(yEdge(nEdges))
		ierr = nf90_get_var(ncid, yEdgeID, yEdge, count = (/nEdges/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable yEdge in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
	
	
	
		!$OMP SECTION
		allocate(zEdge(nEdges))
		ierr = nf90_get_var(ncid, zEdgeID, zEdge, count = (/nEdges/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable zEdge in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
	
		!$OMP SECTION
		allocate(edgesOnVertex(3, nVertices))
		ierr = nf90_get_var(ncid, EoVID, edgesOnVertex, count = (/3, nVertices/))
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error getting variable edgesOnVertex in '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		!$OMP END SECTIONS
		call cpu_time(finish)
        write (*,*) 'Time for retrieving Mesh Info:', finish - start
		
		ierr = nf90_inq_dimid(ncid2, 'nCells', l)
		if (ierr /= NF90_NOERR) then
			nCellsID = l
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nCells dimid in'//filename2
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		else
			nCellsID = l
		end if
		
		ierr = nf90_inq_dimid(ncid2, 'nVertices', l)
		if (ierr /= NF90_NOERR) then
			
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nVertices dimid in'//filename2
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		else
			nVertID = l
		end if
		
		ierr = nf90_inq_dimid(ncid2, 'nEdges', l)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nEdges dimid in'//filename2
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		else
			nEdgesID = l

		end if
		
		ierr = nf90_inq_dimid(ncid2, 'nVertLevels', l)
		if (ierr /= NF90_NOERR) then
			
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nVertLevels dimid in'//filename2
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		else
			nvlID = l
		end if
		
        ierr = nf90_inq_dimid(ncid2, 'nVertLevelsP1', l)
        if (ierr /= NF90_NOERR) then
            
            write(0,*)'*********************************************************************************'
            write(0,*) 'Non-fatal error inquiring nVertLevelsP1 dimid in'//filename2
            write(0,*) 'ierr = ', ierr
            write(0,*)'*********************************************************************************'
            
        else
            nvlP1ID = l
        end if

		ierr = nf90_inq_dimid(ncid2, 'nSoilLevels', l)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error inquiring nSoilLevels dimid in'//filename2
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			
		else 
			nslID = l
		end if
		
		
		
		meshDimIDRef = (/nCellsID, nVertID, nEdgesID, TimeID, nvlID, nslID, nvlP1ID/)
		dimSizes = (/nCells, nVertices, nEdges, elapsedTime, nVertLevels, nSoilLevels, nVertLevelsP1/)
	end subroutine open_input
		
		
	
		
		
		
		

		
	subroutine create_grid_map(grid, rotated)
		implicit none
		real, dimension(:,:,:), allocatable, intent(out) :: grid
		logical, intent(in) :: rotated
                real :: start, finish
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!	Create the grid map containing lon, lat, nearest Mesh Cell, vertex, and edge information. It is a lookup table for the rest of our stuff. 
		allocate(grid(gridW, gridH, 5)) !Note // grid(gridW, gridH, :) = (/lon, lat, nearestCell, nearestVertex, nearestEdge/)
		cellmin = 1
        call cpu_time(start)
        !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(j, x_search, y_search, z_search, moved, k, d, dn, nearestOne, nearestVert, nearestEdge) FIRSTPRIVATE(cellmin)
		do i = 1, gridW
            do j = 1, gridH
				grid(i, j, 1) = gridWmin + (i - 1) * lonSpacing !longitude
				grid(i, j, 2) = gridHmin + (j - 1) * latSpacing	!latitude
				x_search = radius * cos(grid(i, j, 1)) * cos(grid(i, j, 2))
				y_search = radius * sin(grid(i, j, 1)) * cos(grid(i, j, 2))
				z_search = radius * sin(grid(i, j, 2))
				moved = 1
				do while (moved /= 0)
			   
					moved = 0
					d = sqrt((MeshX(cellmin) - x_search)**2.0 + (MeshY(cellmin) - y_search)**2.0 + (MeshZ(cellmin) - z_search)**2.0)
					do k=1, nCellsOnCell(cellmin) 
						dn = sqrt((MeshX(cellsOnCell(k,cellmin)) - x_search)**2.0 + (MeshY(cellsOnCell(k,cellmin)) - y_search)**2.0 + (MeshZ(cellsOnCell(k,cellmin)) - z_search)**2.0)    
						if (dn < d) then     
							d = dn
							nearestOne = cellsOnCell(k,cellmin)
							moved = 1
						end if
					
					end do
				
					if (moved == 1) then
						cellmin = nearestOne
					end if
				
				end do
				grid(i,j,3) = cellmin
				if (rotated) then
					grid(i,j,1) = lonCell(cellmin)
					grid(i,j,2) = latCell(cellmin)
				end if
				! Create grid map for nearest Mesh Vertex. Used for vertex-based fields
				d = sqrt((xVertex(verticesOnCell(1, cellmin)) - x_search)**2.0 + (yVertex(verticesOnCell(1, cellmin)) - y_search)**2.0 + (zVertex(verticesOnCell(1, cellmin)) - z_search)**2.0)
				nearestVert = verticesOnCell(1, cellmin)
				do k = 2, maxEdges
					dn = sqrt((xVertex(verticesOnCell(k, cellmin)) - x_search)**2.0 + (yVertex(verticesOnCell(k, cellmin)) - y_search)**2.0 + (zVertex(verticesOnCell(k, cellmin)) - z_search)**2.0)
					if (dn < d) then
						d = dn
						nearestVert = verticesOnCell(k, cellmin)
					end if
				end do
				grid(i, j, 4) = nearestVert
				
				! Create grid map for nearest mesh Edge. Used for Edge-based Fields
				d = sqrt((xEdge(edgesOnVertex(1, nearestVert)) - x_search)**2 + (yEdge(edgesOnVertex(1, nearestVert)) - y_search)**2 + (zEdge(edgesOnVertex(1, nearestVert)) - z_search)**2)
				nearestEdge = edgesOnVertex(1, nearestVert)
				do k = 2, 3
					dn = sqrt((xEdge(edgesOnVertex(k, nearestVert)) - x_search)**2 + (yEdge(edgesOnVertex(k, nearestVert)) - y_search)**2 + (zEdge(edgesOnVertex(k, nearestVert)) - z_search)**2)
					if (dn < d) then 
						d = dn
						nearestEdge = edgesOnVertex(k, nearestVert)
					end if
				end do
				grid(i,j,5) = nearestEdge
			end do
		end do
                !$OMP END DO
                call cpu_time(finish)
                write (*,*) 'Time to create Grid Map:', finish - start
	deallocate(MeshX, MeshY, MeshZ, nCellsOnCell, cellsOnCell, verticesOnCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, edgesOnVertex)
	
	end subroutine create_grid_map
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
end module inputprocessing
	
	
	
