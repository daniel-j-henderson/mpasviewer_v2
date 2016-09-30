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
      
      
      
      
      
      
   subroutine rotate(grid, original_latitude_degrees, original_longitude_degrees, new_latitude_degrees, new_longitude_degrees, birdseye_rotation_counter_clockwise_degrees)

      implicit none

      real(kind=RKIND), intent(in) :: original_latitude_degrees
      real(kind=RKIND), intent(in) :: original_longitude_degrees
      real(kind=RKIND), intent(in) :: new_latitude_degrees
      real(kind=RKIND), intent(in) :: new_longitude_degrees
      real(kind=RKIND), intent(in) :: birdseye_rotation_counter_clockwise_degrees   

      integer :: i, j, varid

      type(interpgrid), intent(inout) :: grid

      real (kind=RKIND) :: original_latitude_radians, original_longitude_radians, new_latitude_radians, new_longitude_radians
      real (kind=RKIND) :: thetaLat, thetaLon, thetaBirdsEye
      real (kind=RKIND) :: x0LongitudeAtEquator, y0LongitudeAtEquator, z0LongitudeAtEquator
      real (kind=RKIND) :: uCrossProduct, vCrossProduct, wCrossProduct
      real (kind=RKIND) :: xNew, yNew, zNew

      real (kind=RKIND), dimension(:,:), allocatable :: x_pt, y_pt, z_pt

      real (kind=RKIND) :: v
      real (kind=RKIND) :: ax, ay, az
      real (kind=RKIND) :: bx, by, bz
      real (kind=RKIND) :: cx, cy, cz

     pii = 2.*asin(1.0)
      omega = 2.0*pii / 86400.0


      allocate(x_pt(grid%nx, grid%ny), y_pt(grid%nx, grid%ny), z_pt(grid%nx, grid%ny))
      do j=1, grid%ny
      do i=1, grid%nx
         call convert_lx(x_pt(i, j), y_pt(i, j), z_pt(i, j), 1.0_RKIND, grid%lats(i, j), grid%lons(i, j))
      end do
      end do

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

      
         call executeRotation(x_pt, y_pt, z_pt, thetaLat, thetaLon, thetaBirdsEye, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)
         call convert_xl(x_pt, y_pt, z_pt, grid%lats, grid%lons)
      

   end subroutine rotate

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Note that the since <u,v,w> must be a unit vector (where u**2 + v**2 + w**2 = 1),
   !  <uCrossProduct, vCrossProduct, wCrossProduct> and <xNew,yNew,zNew> must be unit vectors as well
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine executeRotation(x, y, z, thetaLat, thetaLon, thetaBirdsEye, uCrossProduct, vCrossProduct, wCrossProduct, xNew, yNew, zNew)

      implicit none

         real (kind=RKIND), dimension(:,:), intent(inout) :: x, y, z
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
   
      real (kind=RKIND), dimension(:,:), intent(in) :: x, y, z
      real (kind=RKIND), dimension(:,:), intent(out) :: lat, lon
   
      real (kind=RKIND), dimension(:,:), allocatable :: dl
      real (kind=RKIND) :: clat, eps
      parameter (eps=1.e-10)
      
      integer :: i,j
      integer, dimension(2) :: lens

      lens = shape(x)
   
      allocate(dl, mold=x)
      dl = sqrt(x*x + y*y + z*z)
      lat = asin(z/dl)
   
   !  check for being close to either pole
      do j=1, lens(2)
      do i=1, lens(1)
         if (abs(x(i,j)) > eps) then
            if (abs(y(i,j)) > eps) then
               lon(i,j) = atan(abs(y(i,j)/x(i,j)))
      
               if ((x(i,j) <= 0.) .and. (y(i,j) >= 0.)) then
                  lon(i,j) = pii-lon(i,j)
               else if ((x(i,j) <= 0.) .and. (y(i,j) < 0.)) then
                  lon(i,j) = lon(i,j)+pii
               else if ((x(i,j) >= 0.) .and. (y(i,j) <= 0.)) then
                  lon(i,j) = 2*pii-lon(i,j)
               end if
      
            else ! we're either on longitude 0 or 180
      
               if (x(i,j) > 0) then
                  lon(i,j) = 0.
               else
                  lon(i,j) = pii
               end if
      
            end if
      
         else if (abs(y(i,j)) > eps) then
      
            if (y(i,j) > 0) then
               lon(i,j) = pii/2.
            else
               lon(i,j) = 3.*pii/2.
            end if
      
         else  ! we are at a pole
      
            lon(i,j) = 0.
      
         end if
      end do
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

      real (kind=RKIND), dimension(:,:), intent(inout) :: x, y, z
      real (kind=RKIND), intent(in) :: theta, a, b, c, u, v, w
      real (kind=RKIND) :: xp, yp, zp

      real (kind=RKIND) :: vw2, uw2, uv2
      real (kind=RKIND) :: m
      
      integer :: i, j
      integer, dimension(2) :: lens
      
      lens = shape(x)

      vw2 = v**2.0 + w**2.0
      uw2 = u**2.0 + w**2.0
      uv2 = u**2.0 + v**2.0
      m = sqrt(u**2.0 + v**2.0 + w**2.0)
     
      do j=1, lens(2)
      do i=1, lens(1)
         xp = (a*vw2 + u*(-b*v-c*w+u*x(i,j)+v*y(i,j)+w*z(i,j)) + ((x(i,j)-a)*vw2+u*(b*v+c*w-v*y(i,j)-w*z(i,j)))*cos(theta) + m*(-c*v+b*w-w*y(i,j)+v*z(i,j))*sin(theta))/m**2.0
         yp = (b*uw2 + v*(-a*u-c*w+u*x(i,j)+v*y(i,j)+w*z(i,j)) + ((y(i,j)-b)*uw2+v*(a*u+c*w-u*x(i,j)-w*z(i,j)))*cos(theta) + m*( c*u-a*w+w*x(i,j)-u*z(i,j))*sin(theta))/m**2.0
         zp = (c*uv2 + w*(-a*u-b*v+u*x(i,j)+v*y(i,j)+w*z(i,j)) + ((z(i,j)-c)*uv2+w*(a*u+b*v-u*x(i,j)-v*y(i,j)))*cos(theta) + m*(-b*u+a*v-v*x(i,j)+u*y(i,j))*sin(theta))/m**2.0
      
      ! alternate calculation
      !xp = (a*vw2 - u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(theta)) + x*cos(theta) + (-c*v+b*w-w*y+v*z)*sin(theta)
      !yp = (b*uw2 - v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(theta)) + y*cos(theta) + ( c*u-a*w+w*x-u*z)*sin(theta)
      !zp = (c*uv2 - w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(theta)) + z*cos(theta) + (-b*u+a*v-v*x+u*y)*sin(theta)

         x(i,j) = xp
         y(i,j) = yp
         z(i,j) = zp
      end do
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
