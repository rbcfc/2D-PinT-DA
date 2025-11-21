module shallow_water_utils

contains

subroutine compute_nx_1D(nx,ny,nx_1D,indices_1D)
integer :: i,j,ir
integer :: nx, ny
integer :: nx_1D
integer,dimension(:,:,:),allocatable :: indices_1D

  allocate(indices_1D(nx,ny,3))
  ir = 0
  do j=1,ny
  do i=1,nx
    ir = ir + 1
    indices_1D(i,j,1) = ir ! eta
  enddo
  enddo

  do j=1,ny
  do i=2,nx
    ir = ir + 1
    indices_1D(i,j,2) = ir ! u
  enddo
  enddo

  do j=2,ny
  do i=1,nx
    ir = ir + 1
    indices_1D(i,j,3) = ir ! v
  enddo
  enddo

  nx_1D = ir
  print *,'nx 1D = ',nx_1D
end subroutine compute_nx_1D

subroutine trans_2Dto1D(un,vn,etan,xn,nx,ny,nx_1D,indices_1D)
  integer :: nx,ny
  real,dimension(1:nx+1,0:ny+1) :: un
  real,dimension(0:nx+1,1:ny+1) :: vn
  real,dimension(0:nx+1,0:ny+1) :: etan
  real,dimension(nx_1D) :: xn
  integer :: nx_1D
  integer,dimension(:,:,:),allocatable :: indices_1D
  integer :: i,j

  do j=1,ny
  do i=1,nx
    xn(indices_1D(i,j,1)) = etan(i,j)
  enddo
enddo

  do j=1,ny
  do i=2,nx
    xn(indices_1D(i,j,2)) = un(i,j)
  enddo
enddo

do j=2,ny
do i=1,nx
  xn(indices_1D(i,j,3)) = vn(i,j)
enddo
enddo

end subroutine trans_2Dto1D

subroutine trans_1Dto2D(xn,un,vn,etan,nx,ny,nx_1D,indices_1D)
  integer :: nx, ny
  real,dimension(1:nx+1,0:ny+1) :: un
  real,dimension(0:nx+1,1:ny+1) :: vn
  real,dimension(0:nx+1,0:ny+1) :: etan
  real,dimension(nx_1D) :: xn
  integer :: nx_1D
  integer,dimension(:,:,:),allocatable :: indices_1D
  integer :: i,j

  do j=1,ny
  do i=1,nx
     etan(i,j) = xn(indices_1D(i,j,1))
  enddo
enddo

  etan(0,:) =etan(1,:)
  etan(nx+1,:) = etan(nx,:)
  etan(:,0) = etan(:,1)
  etan(:,ny+1) = etan(:,ny)

  do j=1,ny
  do i=2,nx
    un(i,j) = xn(indices_1D(i,j,2))
  enddo
 enddo

 un(1,:)=0.
 un(nx+1,:)=0.
 un(:,0) = -un(:,1)
 un(:,ny+1) = -un(:,ny)

  do j=2,ny
  do i=1,nx
    vn(i,j) = xn(indices_1D(i,j,3))
  enddo
 enddo

  vn(:,1)=0.
  vn(:,ny+1)=0.
  vn(0,:) = -vn(1,:)
  vn(nx+1,:) = -vn(nx,:)

end subroutine trans_1Dto2D

end module shallow_water_utils
