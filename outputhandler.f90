

module outputhandler
	use netcdf
	use params
	

	
	public
	integer :: nvertlID, nsoillID, nvertlP1ID
	contains
		
	subroutine create_output_file(newFilename)
	
		implicit none
		character(len=*), intent(in) :: newFilename
		
		ierr = nf90_create(newFilename, NF90_CLOBBER, ncidNew)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error creating NetCDF file '//newFilename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if

		ierr = nf90_def_dim(ncidNew, 'timeDim', NF90_UNLIMITED, tdimID)
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error defining dimension TIMEDim of unlimited length in file '//newFilename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   stop
		end if

		ierr = nf90_def_dim(ncidNew, 'nVertLevels', nVertLevels, nvertlID)
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error defining dimension nVertLevels in file '//newFilename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   !stop
		end if
		
		ierr = nf90_def_dim(ncidNew, 'nSoilLevels', nSoilLevels, nsoillID)
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error defining dimension nSoilLevels in file '//newFilename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   !stop
		end if

        ierr = nf90_def_dim(ncidNew, 'nVertLevelsP1', nVertLevelsP1, nvertlP1ID)
        if (ierr /= NF90_NOERR) then
           write(0,*) '*********************************************************************************'
           write(0,*) 'Error defining dimension nVertLevelsP1in file '//newFilename
           write(0,*) 'ierr = ', ierr
           write(0,*) '*********************************************************************************'
           !stop
        end if

		ierr = nf90_def_dim(ncidNew, 'xDim', gridW, xdimID)
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error defining dimension xDim of unlimited length in file '//newFilename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   stop
		end if
	
		ierr = nf90_def_dim(ncidNew, 'yDim', gridH, ydimID)
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error defining dimension xDim of unlimited length in file '//newFilename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   stop
		end if
	
		gridDimIDs = (/xdimID, ydimID, tdimID, nvertlID, nsoillID, nvertlP1ID/) !reference of all the dimIDs we just created. Will grow as more dimensions become supported
	end subroutine create_output_file
	
	
	subroutine check_existence(variableList)
		implicit none
		
		character(len=*), dimension(:), intent(in) :: variableList
		character(len=NF90_MAX_NAME), dimension(MAX_VARIABLES) :: varList
		integer, dimension(MAX_VARIABLES) :: IDholder
		integer :: l, temp
		l = 0
		do i=1,size_of(variableList)
			ierr = nf90_inq_varid(ncid2, trim(variableList(i)), temp)
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error getting variable ID of ', trim(variableList(i)), 'from file '//filename2
				write(0,*) 'Skipping this variable.'
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
			else 
				l = l+1
				varList(l) = trim(variableList(i))
				IDholder(l) = temp
			end if 
		end do
		if (l == 0) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'None of your desired variables are in the input data file'
			write(0,*) '*********************************************************************************'
		else
			allocate(desiredMeshVars(l), meshVarIDs(l))
			desiredMeshVars = varList(1:l)
			meshVarIDs = IDholder(1:l)
			nMeshVars = l
		end if
	
	end subroutine check_existence
	
	
	
	
	
	
	subroutine dimensionality(desiredMeshVars)
		implicit none
		character(len=NF90_MAX_NAME), dimension(:), intent(in) :: desiredMeshVars
		integer, dimension(:), allocatable :: newDimIDs
		
		allocate(meshVarType(nMeshVars))
		allocate(meshDimIDs(nMeshVars, MAX_MPAS_DIMENSION))
		allocate(gridVarIDs(nMeshVars))
		allocate(nDims(nMeshVars))
		
		do i=1, nMeshVars 
			
			ierr = nf90_inquire_variable(ncid2, meshVarIDs(i), xtype=meshVarType(i), ndims=nDims(i), dimids=meshDimIDs(i,:))
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error inquiring variable', trim(desiredMeshVars(i)), 'from file '//filename2
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if 
			
			!Sets up the dimensions of the new variable to be created in the output file, puts in newDimIDs
			call setup_dimensions(meshDimIDs(i,:), nDims(i), newDimIDs)

			ierr = nf90_def_var(ncidNew, 'grid'//trim(desiredMeshVars(i)), meshVarType(i), newDimIDs, gridVarIDs(i))
			if (ierr /= NF90_NOERR) then
			   write(0,*) '*********************************************************************************'
			   write(0,*) 'Error defining variable', 'grid'//trim(desiredMeshVars(i)), 'in file '//newFilename
			   write(0,*) 'ierr = ', ierr
			   write(0,*) '*********************************************************************************'
			   stop
			end if
			
			ierr = nf90_put_att(ncidNew, gridVarIDs(i), 'coordinates', 'gridLons gridLats')
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error putting new attribute in'//newFilename
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if 
		
		end do
		
		ierr = nf90_def_var(ncidNew, 'gridLons', NF90_REAL, (/xdimID, ydimID/), gridVarIDs(nMeshVars + 1))
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error defining variable gridLats in file '//newFilename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   stop
		end if
		
		ierr = nf90_def_var(ncidNew, 'gridLats', NF90_REAL, (/xdimID, ydimID/), gridVarIDs(nMeshVars + 2))
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error defining variable gridLats in file '//newFilename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   stop
		end if
	
		ierr = nf90_put_att(ncidNew, gridVarIDs(nMeshVars+1), 'units', 'degree_north')
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error putting new attribute in'//newFilename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if 
			
		ierr = nf90_put_att(ncidNew, gridVarIDs(nMeshVars+2), 'units', 'degree_east')
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error putting new attribute in'//newFilename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
			
		ierr = nf90_enddef(ncidNew)
		if (ierr /= NF90_NOERR) then
		   write(0,*) '*********************************************************************************'
		   write(0,*) 'Error ending def in file '//filename
		   write(0,*) 'ierr = ', ierr
		   write(0,*) '*********************************************************************************'
		   stop
		end if
	end subroutine dimensionality
	
	
	
	
	! Takes the dimensions of a mesh variable and translates to the corresponding dimensions on the grid,
	! keeping time on the rightmost dimension if it is there at all
	! sourceIDs is the array of dimIDs for the input data variable, num is number of dims, new DimIDs is corresponding grid dimensions (output)
	subroutine setup_dimensions(sourceIDs, num, newDimIDs)
		implicit none
		
		integer, dimension(MAX_MPAS_DIMENSION), intent(in) :: sourceIDs
		integer, intent(in) :: num
		integer, dimension(:), allocatable, intent(out) :: newDimIDs
		integer :: l
		logical :: found = .false.
		allocate(newDimIDs(num+1))
		l=3
		do k=1,num
			if(sourceIDs(k) == TimeID) then
				newDimIDs(num + 1) = tdimID
			else if (is_spatial(sourceIDs(k))) then
				
				newDimIDs(1) = gridDimIDs(1)
				newDimIDs(2) = gridDimIDs(2)
				
			else
				do j=1,NUM_VALID_MPAS_DIMS
					if (sourceIDs(k) == meshDimIDRef(j)) then
						newDimIDs(l) = gridDimIDs(j-1)
						l=l+1
						found = .true.
					end if
				end do
				if (.not. found) then
					write(0,*) '*********************************************************************************'
					write(*,*) 'One of your variables has a dimension which is not supported.'
					write(*,*) 'Program Terminated.'
					write(0,*) '*********************************************************************************'
					stop
				end if
				found = .false.
			end if
		end do
	
	end subroutine setup_dimensions
	
	
	
	
	
	
	
	
	
	subroutine put_data1(desiredVarID, varDimID, gridVarID)
		implicit none
		integer, intent(in) :: desiredVarID, varDimID, gridVarID
		real, dimension(:), allocatable :: meshData
		real, dimension(gridW, gridH) :: putData
		integer :: j, k
		if (varDimID == nCellsID) then
			allocate(meshData(nCells))
			ierr = nf90_get_var(ncid2, desiredVarID, meshData(:))
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if	
			
			do j=1, gridW
				do k=1, gridH
					putData(j,k) = meshData(grid(j, k, 3))
				end do
			end do
			
			ierr = nf90_put_var(ncidNew, gridVarID, putData, count=(/gridW, gridH/))
				if (ierr /= NF90_NOERR) then
					write(0,*) '*********************************************************************************'
					write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
					write(0,*) 'ierr = ', ierr
					write(0,*) '*********************************************************************************'
					stop
				else 
					write(0,*) '*********************************************************************************'
					write(0,*) 'Successfully put a 1d variable in '//newFilename
					write(0,*) '*********************************************************************************'
				end if
			deallocate(meshData)
		
		else if (varDimID == nVertID) then
			allocate(meshData(nVertices))
			ierr = nf90_get_var(ncid2, desiredVarID, meshData(:))
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if	
			
			do j=1, gridW
				do k=1, gridH
					putData(j,k) = meshData(grid(j, k, 4))
				end do
			end do
			
			ierr = nf90_put_var(ncidNew, gridVarID, putData, count=(/gridW, gridH/))
				if (ierr /= NF90_NOERR) then
					write(0,*) '*********************************************************************************'
					write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
					write(0,*) 'ierr = ', ierr
					write(0,*) '*********************************************************************************'
					stop
				else 
					write(0,*) '*********************************************************************************'
					write(0,*) 'Successfully put a 1d variable in '//newFilename
					write(0,*) '*********************************************************************************'
				end if
			deallocate(meshData)
		
		else if (varDimID == nEdgesID) then
			allocate(meshData(nEdges))
			ierr = nf90_get_var(ncid2, desiredVarID, meshData(:))
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if	
			
			do j=1, gridW
				do k=1, gridH
					putData(j,k) = meshData(grid(j, k, 5))
				end do
			end do
			
			ierr = nf90_put_var(ncidNew, gridVarID, putData, count=(/gridW, gridH/))
				if (ierr /= NF90_NOERR) then
					write(0,*) '*********************************************************************************'
					write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
					write(0,*) 'ierr = ', ierr
					write(0,*) '*********************************************************************************'
					stop
				else 
					write(0,*) '*********************************************************************************'
					write(0,*) 'Successfully put a 1d variable in '//newFilename
					write(0,*) '*********************************************************************************'
				end if
			deallocate(meshData)

		else 
			write(*,*) 'One of your chosen 1-D variables is not of a supported spatial dimension'
		
		end if
		
		
	end subroutine put_data1
	
	
	
	
	
	
	subroutine put_data2(desiredVarID, varDimIDs, gridVarID)
		implicit none
			integer, intent(in) :: desiredVarID, gridVarID
			integer, dimension(2), intent(in) :: varDimIDs
			integer :: spatialDim, i, j, k, spatialDimSize, otherDim                      
			real, dimension(:,:,:), allocatable :: putData
			real, dimension(:,:), allocatable :: meshData
			logical :: tooBig = .false., left = .false.
     	    		integer(kind=8) :: bigness, x, y
     	    		integer, dimension(2) :: dimensionSizes
     	    
     	    
			if ((varDimIDs(1) == nCellsID) .or. (varDimIDs(2) == nCellsID)) then
					spatialDim = 3
					spatialDimSize = nCells
					if ((varDimIDs(2) == TimeID) .or. (varDimIDs(1) == TimeID)) then
						otherDim = elapsedTime
						if (varDimIDs(1) == TimeID) then
							left = .true.
						end if
					!case(nvlID)
					else if ((varDimIDs(2) == nvlID) .or. (varDimIDs(1) == nvlID)) then
						otherDim = nVertLevels
						if (varDimIDs(1) == nvlID) then
							left = .true.
						end if
					!case(nslID)
					else if ((varDimIDs(2) == nslID) .or. (varDimIDs(1) == nslID)) then
						otherDim = nSoilLevels
						if (varDimIDs(1) == nslID) then
							left = .true.
						end if
                    else if ((varDimIDs(2) == nvlP1ID) .or. (varDimIDs(1) == nvlP1ID)) then
                        otherDim = nVertLevelsP1
                        if (varDimIDs(1) == nvlP1ID) then
                            left = .true.
                        end if
					else
						print *, 'Invalid dimension of 2-d variable'
					end if
			else if ((varDimIDs(1) == nVertID) .or. (varDimIDs(2) == nVertID)) then
					spatialDimSize = nVertices
					spatialDim = 4
					if ((varDimIDs(2) == TimeID) .or. (varDimIDs(1) == TimeID)) then
						otherDim = elapsedTime
						if (varDimIDs(1) == TimeID) then
							left = .true.
						end if
						
					!case(nvlID)
					else if ((varDimIDs(2) == nvlID) .or. (varDimIDs(1) == nvlID)) then
						otherDim = nVertLevels
						if (varDimIDs(1) == nvlID) then
							left = .true.
						end if
						
					!case(nslID)
					else if ((varDimIDs(2) == nslID) .or. (varDimIDs(1) == nslID)) then
						otherDim = nSoilLevels
						if (varDimIDs(1) == nslID) then
							left = .true.
						end if
						
                    else if ((varDimIDs(2) == nvlP1ID) .or. (varDimIDs(1) == nvlP1ID)) then
                        otherDim = nVertLevelsP1
                        if (varDimIDs(1) == nvlP1ID) then
                            left = .true.
                        end if
					else
						print *, 'Invalid dimension of 2-d variable'
					end if
			else if ((varDimIDs(1) == nEdgesID) .or. (varDimIDs(2) == nEdgesID)) then
					spatialDim = 5
					spatialDimSize = nEdges
					if ((varDimIDs(2) == TimeID) .or. (varDimIDs(1) == TimeID)) then
						otherDim = elapsedTime
						if (varDimIDs(1) == TimeID) then
							left = .true.
						end if
						
					!case(nvlID)
					else if ((varDimIDs(2) == nvlID) .or. (varDimIDs(1) == nvlID)) then
						otherDim = nVertLevels
						if (varDimIDs(1) == nvlID) then
							left = .true.
						end if
					!case(nslID)
					else if ((varDimIDs(2) == nslID) .or. (varDimIDs(1) == nslID)) then
						otherDim = nSoilLevels
						if (varDimIDs(1) == nslID) then
							left = .true.
						end if
					
                    !case(nvlP1ID)
                    else if ((varDimIDs(2) == nvlP1ID) .or. (varDimIDs(1) == nvlP1ID)) then
                        otherDim = nVertLevelsP1
                        if (varDimIDs(1) == nvlP1ID) then
                            left = .true.
                        end if
                    else
						print *, 'Invalid dimension of 2-d variable'
					end if
			else 
				print *, 'No spatial dimension in one of your 2-d variables, could not put any data.'
			end if
			

			
			
			if(left) then
				dimensionSizes = (/otherDim, spatialDimSize/)
			else
				dimensionSizes = (/spatialDimSize, otherDim/)
			end if
                        x = spatialDimSize
                        y = otherDim
                        bigness = x*y
			if (bigness > 800000000) then
				tooBig = .true.
				if (.not. left) then
					dimensionSizes = (/spatialDimSize, 1/)
				else
					dimensionSizes = (/1, spatialDimSize/)
				end if
				write (*,*) 'Your data for your 2-d variable is so large that we are going to handle it in 1-d slices.'    
				allocate(meshData(dimensionSizes(1), dimensionSizes(2)))
				allocate(putData(gridW, gridH, 1))
			else
				allocate(meshData(dimensionSizes(1), dimensionSizes(2)))
				allocate(putData(gridW, gridH, otherDim))
			end if
			

			if (tooBig) then
				do i=1, otherDim
					if (.not. left) then
					ierr = nf90_get_var(ncid2, desiredVarID, meshData(:,:), (/1, i/), (/dimensionSizes(1), 1/))
					if (ierr /= NF90_NOERR) then
						write(0,*) '*********************************************************************************'
						write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
						write(0,*) 'ierr = ', ierr
						write(0,*) '*********************************************************************************'
						stop
					end if
					
					do j=1, gridW
						do k=1, gridH
							putData(j,k,1) = meshData(grid(j, k, spatialDim),1)
						end do
					end do
						
					else
						ierr = nf90_get_var(ncid2, desiredVarID, meshData, (/i, 1/), (/1, spatialDimSize/))
						if (ierr /= NF90_NOERR) then
							write(0,*) '*********************************************************************************'
							write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
							write(0,*) 'ierr = ', ierr
							write(0,*) '*********************************************************************************'
							stop
						end if
						do j=1, gridW
							do k=1, gridH
								putData(j,k,1) = meshData(1, grid(j, k, spatialDim))
							end do
						end do
					end if
					ierr = nf90_put_var(ncidNew, gridVarID, putData(:,:,1), (/1, 1, i/), count=(/gridW, gridH, 1/))
					if (ierr /= NF90_NOERR) then
						write(0,*) '*********************************************************************************'
						write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
						write(0,*) 'at slice:', i
						write(0,*) 'ierr = ', ierr
						write(0,*) '*********************************************************************************'
						stop
					else 
						write(0,*) '*********************************************************************************'
						write(0,*) 'Successfully put slice number', i, 'of a 2d variable in '//newFilename
						write(0,*) '*********************************************************************************'
					end if
                                        
				end do
				
			else !If not too big
				ierr = nf90_get_var(ncid2, desiredVarID, meshData(:,:))
				if (ierr /= NF90_NOERR) then
					write(0,*) '*********************************************************************************'
					write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
					write(0,*) 'ierr = ', ierr
					write(0,*) '*********************************************************************************'
					stop
				end if
				
				if(.not. left) then
					do j=1, gridW
						do k=1, gridH
							putData(j,k,:) = meshData(grid(j, k, spatialDim), :)
						end do
					end do
				else
					do j=1, gridW
						do k=1, gridH
							putData(j,k,:) = meshData(:, grid(j, k, spatialDim))
						end do
					end do
				end if
			
			
			
				ierr = nf90_put_var(ncidNew, gridVarID, putData, count=(/gridW, gridH, otherDim/))
				if (ierr /= NF90_NOERR) then
					write(0,*) '*********************************************************************************'
					write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
					write(0,*) 'ierr = ', ierr
					write(0,*) '*********************************************************************************'
					stop
				else 
					write(0,*) '*********************************************************************************'
					write(0,*) 'Successfully put a 2d variable in '//newFilename
					write(0,*) '*********************************************************************************'
				end if	
                        end if
			deallocate(meshData)
			deallocate(putData)
				
	end subroutine put_data2
		
		
		
	
	
	
	
	subroutine put_data3(desiredVarID, varDimIDs, gridVarID)
		implicit none 
		
		integer, intent(in) :: desiredVarID, gridVarID
		integer, dimension(3), intent(in) :: varDimIDs
		
		integer :: spatialDim, dim2, dim3, j, k, m, n, l, spatialDimIndex
		integer, dimension(3,2) :: whichDim
		integer, dimension(2) :: dataDim
		real, dimension(:,:,:,:), allocatable :: putData
		real, dimension(:,:,:), allocatable :: meshData
		logical :: tooBig = .false.
        integer(kind=8) :: bigness, x, y, z
		
 		l=0
 		do j=1,3
 			do k=1, NUM_VALID_MPAS_DIMS
 				if (varDimIDs(j) == meshDimIDRef(k)) then
 					l=l+1
 					if (is_spatial(varDimIDs(j))) then
 						spatialDimIndex = j
 						spatialDim = k+2
 					else
 						whichDim(l,1) = gridDimIDs(k-1)
 					end if
 					whichDim(l,2) = dimSizes(k)
 				end if
 			end do
 		end do
 		
 		if (l < 3) then
 			write(*,*) '**************************************************'
 			write(*,*) 'Your 3-d Variable does not have 3 supported dimensions'
 			write(*,*) '**************************************************'
 		end if

		x = whichDim(1, 2)
		y = whichDim(2, 2)
		z = whichDim(3, 2)
		
		!Bigness is the number of elements in the array which would hold the data from the MPAS file for this variable.
		!If Bigness is more than 800M, then the array will take up over 3.2 GB of memory, and we'd get kicked
		!off the Yellowstone login node pretty soon, so we'll handle the data in chunks.
        bigness = x*y*z
 		if (bigness > 800000000) then
 			tooBig = .true.
            write (*,*) 'Your data for your 3-d variable is so large that we are going to handle it in 1-d slices.'    
 		end if

 		! Collapse down to get rid of spatial dimension size
 		l=0
 		do j=1,3
 			if (j /= spatialDimIndex) then
 				l=l+1
 				dataDim(l) = whichDim(j,2)
 			end if
 		end do
 
 		if (tooBig) then
 			allocate(putData(gridW, gridH, 1, 1))
 			if (spatialDimIndex == 1) then
 				allocate(meshData(whichDim(1,2), 1, 1))
 				do j=1, dataDim(1)
					do k=1, dataDim(2)
						ierr = nf90_get_var(ncid2, desiredVarID, meshData(:,1,1), (/1, j, k/), (/whichDim(1,2), 1, 1/))
						if (ierr /= NF90_NOERR) then
							write(0,*) '*********************************************************************************'
							write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
							write(0,*) 'ierr = ', ierr
							write(0,*) '*********************************************************************************'
							stop
						end if
						
						do m=1, gridW
							do n=1, gridH
								putData(m,n,1,1) = meshData(grid(m, n, spatialDim),1,1)
							end do
						end do
						
						ierr = nf90_put_var(ncidNew, gridVarID, putData, (/1, 1, j, k/), (/gridW, gridH, 1, 1/))
						if (ierr /= NF90_NOERR) then
							write(0,*) '*********************************************************************************'
							write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
							write(0,*) 'at slice (',j,',',k,')'
							write(0,*) 'ierr = ', ierr
							write(0,*) '*********************************************************************************'
							stop
						else 
							write(0,*) '*********************************************************************************'
							write(0,*) 'Successfully put a 3d variable in '//newFilename
							write(0,*) 'at slice (',j,',',k,')'
							write(0,*) '*********************************************************************************'
						end if
					end do
				end do
				deallocate(meshData)
				
			else if (spatialDimIndex == 2) then
				allocate(meshData(1, whichDim(2,2), 1))
 				do j=1, dataDim(1)
					do k=1, dataDim(2)
						ierr = nf90_get_var(ncid2, desiredVarID, meshData(1,:,1), (/j, 1, k/), (/1, whichDim(2,2), 1/))
						if (ierr /= NF90_NOERR) then
							write(0,*) '*********************************************************************************'
							write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
							write(0,*) 'ierr = ', ierr
							write(0,*) '*********************************************************************************'
							stop
						end if
						
						do m=1, gridW
							do n=1, gridH
								putData(m,n,1,1) = meshData(1,grid(m, n, spatialDim),1)
							end do
						end do
						
						ierr = nf90_put_var(ncidNew, gridVarID, putData, (/1, 1, j, k/), (/gridW, gridH, 1, 1/))
						if (ierr /= NF90_NOERR) then
							write(0,*) '*********************************************************************************'
							write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
							write(0,*) 'at slice (',j,',',k,')'
							write(0,*) 'ierr = ', ierr
							write(0,*) '*********************************************************************************'
							stop
						else 
							write(0,*) '*********************************************************************************'
							write(0,*) 'Successfully put a 3d variable in '//newFilename
							write(0,*) 'at slice (',j,',',k,')'
							write(0,*) '*********************************************************************************'
						end if
					end do
				end do
				deallocate(meshData)
				
			else if (spatialDimIndex == 3) then
				allocate(meshData(1, 1, whichDim(3,2)))
 				do j=1, dataDim(1)
					do k=1, dataDim(2)
						ierr = nf90_get_var(ncid2, desiredVarID, meshData(1,1,:), (/j, k, 1/), (/1, 1, whichDim(3,2)/))
						if (ierr /= NF90_NOERR) then
							write(0,*) '*********************************************************************************'
							write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
							write(0,*) 'ierr = ', ierr
							write(0,*) '*********************************************************************************'
							stop
						end if
						
						do m=1, gridW
							do n=1, gridH
								putData(m,n,1,1) = meshData(1,1,grid(m, n, spatialDim))
							end do
						end do
						
						ierr = nf90_put_var(ncidNew, gridVarID, putData, (/1, 1, j, k/), (/gridW, gridH, 1, 1/))
						if (ierr /= NF90_NOERR) then
							write(0,*) '*********************************************************************************'
							write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
							write(0,*) 'at slice (',j,',',k,')'
							write(0,*) 'ierr = ', ierr
							write(0,*) '*********************************************************************************'
							stop
						else 
							write(0,*) '*********************************************************************************'
							write(0,*) 'Successfully put a 3d variable in '//newFilename
							write(0,*) 'at slice (',j,',',k,')'
							write(0,*) '*********************************************************************************'
						end if
					end do
				end do
				deallocate(meshData)
			else
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error with spatial dim index in put_var3'
				write(0,*) '*********************************************************************************'
			end if
		else
			allocate(putData(gridW, gridH, dataDim(1), dataDim(2)))
			allocate(meshData(whichDim(1,2), whichDim(2,2), whichDim(3,2)))
            ierr = nf90_get_var(ncid2, desiredVarID, meshData(:,:,:))
            if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error getting data from mesh variable with ID', desiredVarID, 'in'//filename
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if
			! Populate the putData array before copying it into the file
			do j=1, gridW
				do k=1, gridH
					do l=1, dataDim(1)
						if (spatialDimIndex == 1) then
						putData(j,k,l,:) = meshData(grid(j, k, spatialDim),l,:)
						else if (spatialDimIndex == 2) then
						putData(j,k,l,:) = meshData(l,grid(j, k, spatialDim),:)
						else if (spatialDimIndex == 3) then
						putData(j,k,l,:) = meshData(l,:,grid(j, k, spatialDim))
						end if
					end do
				end do
			end do
		
		
			ierr = nf90_put_var(ncidNew, gridVarID, putData, count=(/gridW, gridH, dataDim(1), dataDim(2)/))
				if (ierr /= NF90_NOERR) then
					write(0,*) '*********************************************************************************'
					write(0,*) 'Error putting grid variable data with grid varID', gridVarID, 'in'//filename
					write(0,*) 'ierr = ', ierr
					write(0,*) '*********************************************************************************'
					stop
				else 
					write(0,*) '*********************************************************************************'
					write(0,*) 'Successfully put a 3d variable in '//newFilename
					write(0,*) '*********************************************************************************'
				end if	
			deallocate(meshData)
		end if
		
		deallocate(putData)		
 
 	end subroutine put_data3
 
 
	subroutine put_latlons(lons, lats)
		implicit none
		real, dimension(gridW, gridH), intent(in) :: lats, lons
		
		
		
		
		ierr = nf90_put_var(ncidNew, gridVarIDs(nMeshVars+2), lats*r2d)
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error putting grid lats in'//newFilename
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if 
		
			
			
		ierr = nf90_put_var(ncidNew, gridVarIDs(nMeshVars+1), lons*r2d)
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error putting grid lons in'//newFilename
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if 
		
		
		
	end subroutine put_latlons
	
	
	subroutine clean_up()
		ierr = nf90_close(ncid)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error while closing NetCDF file '//filename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		if (ncid2 /= ncid) then
			ierr = nf90_close(ncid2)
			if (ierr /= NF90_NOERR) then
				write(0,*) '*********************************************************************************'
				write(0,*) 'Error while closing NetCDF file '//filename2
				write(0,*) 'ierr = ', ierr
				write(0,*) '*********************************************************************************'
				stop
			end if
		end if
	
		ierr = nf90_close(ncidNew)
		if (ierr /= NF90_NOERR) then
			write(0,*) '*********************************************************************************'
			write(0,*) 'Error while closing NetCDF file '//newFilename
			write(0,*) 'ierr = ', ierr
			write(0,*) '*********************************************************************************'
			stop
		end if
		
		deallocate(meshVarType, MeshDimIDs, gridVarIDs, nDims)
		
	end subroutine clean_up
end module outputhandler
				
			
			
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
