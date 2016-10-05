module params
   integer, parameter :: RKIND=8, StrKIND=100, MAX_NDIMS=30, MAX_NVARS=200, MAX_NATTS=75, DIM=1, VAR=2, ATT=3, NVARS_MAX_NAMELIST=20
   real(kind=RKIND), parameter :: PI = 3.14159265359
   integer, parameter :: NN=1, WP=2

   type :: interpgrid
      integer :: mode = NN
      integer, dimension(:,:), pointer :: cell_map, edge_map, vertex_map
      real (kind=RKIND), dimension(:,:,:), pointer :: cell_weights, edge_weights, vertex_weights
      real (kind=RKIND) :: lat_start, lat_end, lon_start, lon_end
      real (kind=RKIND), dimension(:,:), pointer :: lats, lons
      integer :: nx, ny
   end type interpgrid
end module params
