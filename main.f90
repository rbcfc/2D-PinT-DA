!> @mainpage Main program
!> @authors
!> Rishabh Bhatt, Laurent Debreu, Arthur Vidard

program main

use shallow_water
use parareal_utils
use eigenvalues
use output_netcdf
use shallow_water_adjoint
use inexact_version
  
integer :: i
real, dimension(:), allocatable :: result,tmp,obs,B,x0,x1,x2,test,test2,r
real, dimension(:,:), allocatable :: diffsol, cg_solution
real :: t1,t2,t3,t4
real :: t5,t6,t7,t8
real :: rnorm,xout,eps1
integer :: it, nitermax
real :: e_quad, e_cf

double precision :: dnrm2
call initialise_sw()

call CreatenetCDF()
!call compute_matrix_sw()
!print *,'Matrix computed'

!Computing the true solution using sequential fine solver
allocate(Truesolution(FineGrid%nx_1D,0:N_time_windows))

Truesolution(:,0) = initial_xn

call OutPutNetCDF()

write(*,*)
call my_cpu_time(t3)
do i=1,N_time_windows

  print *,'I = ',i,maxval(FineGrid%eta)
  call integrate_sw(etan=FineGrid%eta,un=Finegrid%u,vn=FineGrid%v,nt=Nfine)

  call trans_2Dto1D(FineGrid%u,FineGrid%v,FineGrid%eta,Truesolution(:,i),FineGrid%nx,FineGrid%ny, &
                FineGrid%nx_1D,FineGrid%indices_1D)
  call OutPutNetCDF()
end do
call my_cpu_time(t4)
t4=t4-t3
print *,'Timing non parareal = ',t4
write(*,*)

call CloseNetCDF()


allocate(diffsol,mold=FineGrid%eta)
diffsol = FineGrid%eta
open(10,file='output',form='formatted')
do i=1,FineGrid%nx
   !write(10,*) FineGrid%xr(i),FineGrid%eta(i)
end do
close(10)

write(*,*) '##############################'

call initialise_parareal(FineGrid,CoarseGrid)
print *,'start parareal'

!call Parareal(FineGrid,CoarseGrid,initial_xn,result,maxk=41,iter_out=it)
 
! call trans_1Dto2D(result,FineGrid%u,FineGrid%v,FineGrid%eta,FineGrid%nx, &
!               FineGrid%ny,FineGrid%nx_1D,FineGrid%indices_1D)
! diffsol=diffsol - FineGrid%eta
! print *,'Error = ',sqrt(DOT_PRODUCT(reshape(diffsol,(/size(diffsol)/)),reshape(diffsol(/size(diffsol)/)))


allocate(obs,B,x0,x1,x2,mold=initial_xn)

obs = initial_xn
do i=1,N_time_windows
  call integrate_sw(xn=obs,nt=Nfine)
end do

do i=1,FineGrid%nx_1D
  ! print *,'obs = ',i,obs(i),initial_xn(i)
end do

call Bvector(FineGrid%nx_1D,obs,B)

do i=1,FineGrid%nx_1D
   !print *,'b = ',i,obs(i),B(i)
end do

print *,'B vector computed'

x0 = 0.

!call Parareal(FineGrid,CoarseGrid,-B,result,maxk=41,iter_out=it)

print *,'Begin conjgrad'
call set_verbose(.TRUE.)
!call my_cpu_time(t1)

if (regularisation) then
  print *, 'Using regularisation constant, alpharegul = ', alpharegul
else
  print *, 'No regularisation used'
end if

write (*,*) '############ Test 1 - Run exact CG with exact matrix ###########'
write (*,*)
!call conjgrad(Mmatrix,B,x0,FineGrid%nx_1D,eps=1.D-2,observation=obs,exact_matrix=Mmatrix,forward_matrix=Fmatrix)

!call my_cpu_time(t2)
!print *,'timing CG = ',t2-t1 
!print *,'difference from initial state = ', dnrm2(FineGrid%nx_1D,x0-initial_xn,1)


!allocate(test2(FineGrid%nx_1D))
!e_cf = 0.5*(dnrm2(FineGrid%nx_1D,test2 - obs,1))**2
!if (regularisation) then
!  e_cf = e_cf + 0.5*alpharegul*(dnrm2(FineGrid%nx_1D,x0,1))**2
!end if
!print *, "exact cost function value ", e_cf
!deallocate(test2)

!allocate(test2(FineGrid%nx_1D))
!call Mmatrix(FineGrid%nx_1D,x0,test2)
!xout = dnrm2(FineGrid%nx_1D,test2 - B,1)
!print *, "exact grad value ", xout
!deallocate(test2)

!allocate(test(FineGrid%nx_1D))
!call Mmatrix(FineGrid%nx_1D,x0,test)
!e_quad = 0.5*DOT_PRODUCT(x0,test) - DOT_PRODUCT(B,x0)
!print *, "exact quadratic value ", e_quad
!deallocate(test)

x1 = 0.

!call my_cpu_time(t5)
write (*,*)
write (*,*) '############ Test 2 - Run exact CG with matrix from parareal ###########'
write (*,*)

call conjgrad(Mmatrix_parareal,B,x1,FineGrid%nx_1D,eps=1.D0,observation=obs,exact_matrix=Mmatrix,forward_matrix=Fmatrix)

!call my_cpu_time(t6)
!print *,'timing CG = ',t6-t5
!stop


!!!!!!!! Find eps for inexact cg !!!!!!!!!!!

!allocate(r(FineGrid%nx_1D),test(FineGrid%nx_1D))

!call Mmatrix(FineGrid%nx_1D,x1,test)
!r = test - B
!call matrixnorminv(Mmatrix,r,FineGrid%nx_1D,rnorm)
!print *, 'rnorm = ', rnorm

!call matrixnorminv(Mmatrix,B,FineGrid%nx_1D,xout)
!print *, 'bnorm = ', xout

!call get_tol(rnorm,xout,eps1)

!print *, "inexact cg eps = ",eps1
!!!!!!!!!!!!!!!!!!!!!!

!deallocate(r)
!eps1 = 9.0799496506838787E-002

! alpharegul = 5
eps1 = 1.0753282746828472E-002
print *, 'inexact cg eps = ', eps1

x2 = 0.
nitermax=20

call my_cpu_time(t7)
write (*,*)
write (*,*) '############ Test 3 - Run inexact CG with parareal ###########'
!call inexact_conjgrad(B,x2,FineGrid%nx_1D,eps=eps1,nitermax=nitermax,observation=obs)

call my_cpu_time(t8)
print *, 'timing inexact CG = ',t8-t7

!call give_output(cg_solution)


contains

subroutine give_output(sol)

real, dimension(:,:), allocatable :: sol


call change_output(.true.)

call initialise_sw(x0)
call CreateNetCDF()
allocate(sol(FineGrid%nx_1D,0:N_time_windows))

initial_xn = x0
sol(:,0) = initial_xn 

call OutPutNetCDF()

do i=1,N_time_windows
  print *,'I = ',i,maxval(FineGrid%eta)
  call integrate_sw(etan=FineGrid%eta,un=FineGrid%u,vn=FineGrid%v,nt=Nfine)

  call trans_2Dto1D(FineGrid%u,FineGrid%v,FineGrid%eta,sol(:,i),FineGrid%nx,FineGrid%ny, &
                FineGrid%nx_1D,FineGrid%indices_1D)
  call outPutNetCDF()
end do

call CloseNetCDF()
deallocate(sol)

end subroutine give_output  




end program main
