!> Shallow water_step module

!> @brief
!! Compute the right hand side of the shallow water equations
module shallow_water_step

contains

!> @param[in] un,vn Velocity components
!> @param[in] etan Free surface elevation
!> @param[out] rhs_u,rhs_v,rhs_eta Right hand sides

subroutine rhs_sw(un,vn,etan,rhs_u,rhs_v,rhs_eta,nx,ny,g,H,coriolis,dx,dy)
integer :: nx, ny
real,dimension(1:nx+1,0:ny+1) :: un        !> x-velocity component
real,dimension(0:nx+1,1:ny+1) :: vn        !> y-velocity component
real,dimension(0:nx+1,0:ny+1) :: etan      !> free surface elevation

real,dimension(1:nx+1,0:ny+1) :: rhs_u
real,dimension(0:nx+1,1:ny+1) :: rhs_v
real,dimension(0:nx+1,0:ny+1) :: rhs_eta

real :: dx, dy
real :: g,H,coriolis

integer :: i,j

real, parameter :: beta = 0.00

do j=1,ny
do i=2,nx
   rhs_u(i,j) = -g * (etan(i,j)-etan(i-1,j))/dx + &
     (coriolis+j*beta)*0.25*(vn(i-1,j)+vn(i,j)+vn(i-1,j+1)+vn(i,j+1))
enddo
enddo

do j=2,ny
do i=1,nx
   rhs_v(i,j) = -g * (etan(i,j)-etan(i,j-1))/dy - &
   (coriolis+j*beta)*0.25*(un(i,j-1)+un(i+1,j-1)+un(i,j)+un(i+1,j))
enddo
enddo

do j=1,ny
do i=1,nx
   rhs_eta(i,j) = -H * (un(i+1,j)-un(i,j))/dx -H * (vn(i,j+1)-vn(i,j))/dy
enddo
enddo

end subroutine rhs_sw
end module shallow_water_step
