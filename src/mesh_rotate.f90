!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE MESH_ROTATE
!
! Rotate a mesh on a unit sphere from a reference latitude,longitude point 
! to a destination latitude,longitude point and optionally rotate the mesh 
! counter-cloclwise around the destination point. This simply changes the 
! values for {x,y,z}{Cell,Vertex,Edge} that create_grid_map will use.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mesh_rotate

  use netcdf
  use params
  implicit none

  integer, parameter :: RKIND = 8
  real(kind=RKIND) :: pii, omega


  



contains
	
	
   function needs_rotated(original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees)
      implicit none
      
      real(kind=RKIND), intent(in) :: original_latitude_degrees
      real(kind=RKIND), intent(in) :: original_longitude_degrees
      real(kind=RKIND), intent(in) :: new_latitude_degrees
      real(kind=RKIND), intent(in) :: new_longitude_degrees
      real(kind=RKIND), intent(in) :: birdseye_rotation_counter_clockwise_degrees   
      logical :: needs_rotated
      needs_rotated = .false.
      
      if (original_latitude_degrees /= 0.0) then
           needs_rotated = .true.
           return
      else if (original_longitude_degrees /= 0.0) then
           needs_rotated = .true.
           return
      else if (new_latitude_degrees /= 0.0) then
           needs_rotated = .true.
           return
      else if (new_longitude_degrees /= 0.0) then
           needs_rotated = .true.
           return
      else if (birdseye_rotation_counter_clockwise_degrees /= 0.0) then
           needs_rotated = .true.
           return
      end if
      
      return 
    end function needs_rotated
      
      
      
      
      
      
   subroutine rotate(ncid, xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, &
    				 original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees)

      implicit none

      integer, intent(in) :: ncid

      real(kind=RKIND), intent(in) :: original_latitude_degrees
      real(kind=RKIND), intent(in) :: original_longitude_degrees
      real(kind=RKIND), intent(in) :: new_latitude_degrees
      real(kind=RKIND), intent(in) :: new_longitude_degrees
      real(kind=RKIND), intent(in) :: birdseye_rotation_counter_clockwise_degrees   

      integer :: i, varid

      real(kind=RKIND), dimension(:), intent(inout) :: xCell
      real(kind=RKIND), dimension(:), intent(inout) :: yCell
      real(kind=RKIND), dimension(:), intent(inout) :: zCell
      real(kind=RKIND), dimension(:), allocatable :: latCell
      real(kind=RKIND), dimension(:), allocatable :: lonCell

      real(kind=RKIND), dimension(:), intent(inout) :: xVertex
      real(kind=RKIND), dimension(:), intent(inout) :: yVertex
      real(kind=RKIND), dimension(:), intent(inout) :: zVertex
      real(kind=RKIND), dimension(:), allocatable :: latVertex
      real(kind=RKIND), dimension(:), allocatable :: lonVertex

      real(kind=RKIND), dimension(:), intent(inout) :: xEdge
      real(kind=RKIND), dimension(:), intent(inout) :: yEdge
      real(kind=RKIND), dimension(:), intent(inout) :: zEdge
      real(kind=RKIND), dimension(:), allocatable :: latEdge
      real(kind=RKIND), dimension(:), allocatable :: lonEdge


      real (kind=RKIND) :: original_latitude_radians, original_longitude_radians, new_latitude_radians, new_longitude_radians
      real (kind=RKIND) :: thetaLat, thetaLon, thetaBirdsEye
      real (kind=RKIND) :: x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator
      real (kind=RKIND) :: uCrossProduct, vCrossProduct, wCrossProduct
      real (kind=RKIND) :: xNew, yNew, zNew

      real (kind=RKIND) :: v
      real (kind=RKIND) :: ax, ay, az
      real (kind=RKIND) :: bx, by, bz
      real (kind=RKIND) :: cx, cy, cz

      allocate(latCell(1:nCells))
      allocate(lonCell(1:nCells))
      
      allocate(latVertex(1:nVertices))
      allocate(lonVertex(1:nVertices))

      
      allocate(latEdge(1:nEdges))
      allocate(lonEdge(1:nEdges))
	  
	  pii = 2.*asin(1.0)
      omega = 2.0*pii / 86400.0




      original_latitude_radians = degreesToRadians(original_latitude_degrees)
      original_longitude_radians = degreesToRadians(original_longitude_degrees)
      new_latitude_radians = degreesToRadians(new_latitude_degrees)
      new_longitude_radians = degreesToRadians(new_longitude_degrees)

      thetaLat = new_latitude_radians - original_latitude_radians
      thetaLon = new_longitude_radians - original_longitude_radians
      thetaBirdsEye = degreesToRadians(birdseye_rotation_counter_clockwise_degrees)

      ! create the unit vector <x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator>
      call convert_lx(x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator, 1.0, 0.0, original_longitude_radians)

      ! create the unit vector <xNew, yNew, zNew>
      call convert_lx(xNew, yNew, zNew, 1.0_RKIND, new_latitude_radians, new_longitude_radians)

      ! create the unit vector <uCrossProduct, vCrossProduct, wCrossProduct> by using a right-angle cross-product of two unit vectors in the perpendicular plane
      call cross_product(x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator, &
                         0.0_RKIND, 0.0_RKIND, 1.0_RKIND, &
                         uCrossProduct, vCrossProduct, wCrossProduct)

      
         call executeRotation(xCell, yCell, zCell, thetaLat, thetaLon, thetaBirdsEye, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)
         call convert_xl(xCell, yCell, zCell, latCell, lonCell)
      

     
         call executeRotation(xVertex, yVertex, zVertex, thetaLat, thetaLon, thetaBirdsEye, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)
         call convert_xl(xVertex, yVertex, zVertex, latVertex, lonVertex)
     

      
         call executeRotation(xEdge, yEdge, zEdge, thetaLat, thetaLon, thetaBirdsEye, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)
         call convert_xl(xEdge, yEdge, zEdge, latEdge, lonEdge)
         
        


      


   end subroutine rotate

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Note that the since <u,v,w> must be a unit vector (where u**2 + v**2 + w**2 = 1),
   !  <uCrossProduct, vCrossProduct, wCrossProduct> and <xNew,yNew,zNew> must be unit vectors as well
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine executeRotation(x, y, z, thetaLat, thetaLon, thetaBirdsEye, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)

      implicit none

         real (kind=RKIND), dimension(:), intent(inout) :: x, y, z
         real (kind=RKIND), intent(in) :: thetaLat, thetaLon, thetaBirdsEye
         real (kind=RKIND), intent(in) :: uCrossProduct, vCrossProduct, wCrossProduct
         real (kind=RKIND), intent(in) :: xNew, yNew, zNew
         real (kind=RKIND) ::  u, v, w

         ! latitude rotation (rotate around cross product of xyz corresponding to original point's longitude at the equator and the z axis unit vector)
         u = uCrossProduct
         v = vCrossProduct
         w = wCrossProduct
        
         call rotate_about_vector(x, y, z, thetaLat, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND, u, v, w)

         ! longitude rotation (rotate around z axis)
         u = 0.0
         v = 0.0
         w = 1.0
         call rotate_about_vector(x, y, z, thetaLon, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND, u, v, w)

         ! bird's eye rotation (rotate around vector from origin to geolocation)
         u = xNew
         v = yNew
         w = zNew
         call rotate_about_vector(x, y, z, thetaBirdsEye, 0.0_RKIND, 0.0_RKIND, 0.0_RKIND, u, v, w)

   end subroutine executeRotation


   real function degreesToRadians(degAngle)

      implicit none

      real(kind=RKIND) :: degAngle
      degreesToRadians = degAngle * 2 * pii / 360 
   end function degreesToRadians






   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CONVERT_XL
   !
   ! Convert (x, y, z) to a (lat, lon) location on a sphere with
   !    radius sqrt(x^2 + y^2 + z^2).
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine convert_xl(x, y, z, lat, lon)
   
      implicit none
   
      real (kind=RKIND), dimension(:), intent(in) :: x, y, z
      real (kind=RKIND), dimension(:), intent(out) :: lat, lon
   
      real (kind=RKIND), dimension(size(x)) :: dl
      real (kind=RKIND) :: clat, eps
      parameter (eps=1.e-10)
      
      integer :: i
   
      dl = sqrt(x*x + y*y + z*z)
      lat = asin(z/dl)
   
   !  check for being close to either pole
   do i=1, size(x)
      if (abs(x(i)) > eps) then
   
         if (abs(y(i)) > eps) then
   
            lon(i) = atan(abs(y(i)/x(i)))
   
            if ((x(i) <= 0.) .and. (y(i) >= 0.)) then
               lon(i) = pii-lon(i)
            else if ((x(i) <= 0.) .and. (y(i) < 0.)) then
               lon(i) = lon(i)+pii
            else if ((x(i) >= 0.) .and. (y(i) <= 0.)) then
               lon(i) = 2*pii-lon(i)
            end if
   
         else ! we're either on longitude 0 or 180
   
            if (x(i) > 0) then
               lon(i) = 0.
            else
               lon(i) = pii
            end if
   
         end if
   
      else if (abs(y(i)) > eps) then
   
         if (y(i) > 0) then
            lon(i) = pii/2.
         else
            lon(i) = 3.*pii/2.
         end if
   
      else  ! we are at a pole
   
         lon(i) = 0.
   
      end if
	end do
   end subroutine convert_xl


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CONVERT_LX
   !
   ! Convert (lat,lon) to an (x, y, z) location on a sphere with specified radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine convert_lx(x, y, z, radius, lat, lon)

      implicit none

      real (kind=RKIND), intent(in) :: radius
      real (kind=RKIND), intent(out) :: x, y, z
      real (kind=RKIND), intent(in) :: lat, lon
      
      integer :: i
      
	  
      z = radius * sin(lat)
      x = radius * cos(lon) * cos(lat)
      y = radius * sin(lon) * cos(lat)
   end subroutine convert_lx


   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ROTATE_ABOUT_VECTOR
   !
   ! Rotates the point (x,y,z) through an angle theta about the line through (a,b,c) 
   ! with direction vector <u,v,w>. Note that the uvw must describe a unit vector (where u**2 + v**2 + w**2 = 1) 
   !
   ! Reference: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine rotate_about_vector(x, y, z, theta, a, b, c, u, v, w)

      implicit none

      real (kind=RKIND), dimension(:), intent(inout) :: x, y, z
      real (kind=RKIND), intent(in) :: theta, a, b, c, u, v, w
      real (kind=RKIND) :: xp, yp, zp

      real (kind=RKIND) :: vw2, uw2, uv2
      real (kind=RKIND) :: m
      
      integer :: i
      


      vw2 = v**2.0 + w**2.0
      uw2 = u**2.0 + w**2.0
      uv2 = u**2.0 + v**2.0
      m = sqrt(u**2.0 + v**2.0 + w**2.0)
	  
	  do i=1, size(x)
      xp = (a*vw2 + u*(-b*v-c*w+u*x(i)+v*y(i)+w*z(i)) + ((x(i)-a)*vw2+u*(b*v+c*w-v*y(i)-w*z(i)))*cos(theta) + m*(-c*v+b*w-w*y(i)+v*z(i))*sin(theta))/m**2.0
      yp = (b*uw2 + v*(-a*u-c*w+u*x(i)+v*y(i)+w*z(i)) + ((y(i)-b)*uw2+v*(a*u+c*w-u*x(i)-w*z(i)))*cos(theta) + m*( c*u-a*w+w*x(i)-u*z(i))*sin(theta))/m**2.0
      zp = (c*uv2 + w*(-a*u-b*v+u*x(i)+v*y(i)+w*z(i)) + ((z(i)-c)*uv2+w*(a*u+b*v-u*x(i)-v*y(i)))*cos(theta) + m*(-b*u+a*v-v*x(i)+u*y(i))*sin(theta))/m**2.0
      
      ! alternate calculation
      !xp = (a*vw2 - u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(theta)) + x*cos(theta) + (-c*v+b*w-w*y+v*z)*sin(theta)
      !yp = (b*uw2 - v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(theta)) + y*cos(theta) + ( c*u-a*w+w*x-u*z)*sin(theta)
      !zp = (c*uv2 - w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(theta)) + z*cos(theta) + (-b*u+a*v-v*x+u*y)*sin(theta)

      x(i) = xp
      y(i) = yp
      z(i) = zp
      end do


   end subroutine rotate_about_vector


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! SUBROUTINE CROSS_PRODUCT
   !
   ! Computes C = A x B
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   subroutine cross_product(Ax, Ay, Az, &
                            Bx, By, Bz, &
                            Cx, Cy, Cz)

      real (kind=RKIND), intent(in)  :: Ax, Ay, Az
      real (kind=RKIND), intent(in)  :: Bx, By, Bz
      real (kind=RKIND), intent(out) :: Cx, Cy, Cz

      Cx = (Ay * Bz) - (Az * By)
      Cy = (Az * Bx) - (Ax * Bz)
      Cz = (Ax * By) - (Ay * Bx)

   end subroutine cross_product                                 

end module
