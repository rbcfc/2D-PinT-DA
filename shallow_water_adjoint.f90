module shallow_water_adjoint
use shallow_water
use SHALLOW_WATER_STEP_DIFF

contains
SUBROUTINE INTEGRATE_SW_adjoint(xn, xnb, etan, etanb, un, unb, vn, vnb, nt)
  USE EIGENVALUES

  REAL, DIMENSION(curgrid%nx_1d), OPTIONAL :: xn
  REAL, DIMENSION(curgrid%nx_1d), OPTIONAL :: xnb
  REAL, DIMENSION(0:curgrid%nx+1,0:curgrid%ny+1), OPTIONAL :: etan
  REAL, DIMENSION(0:curgrid%nx+1,0:curgrid%ny+1), OPTIONAL :: etanb
  REAL, DIMENSION(1:curgrid%nx+1,0:curgrid%ny+1), OPTIONAL :: un
  REAL, DIMENSION(1:curgrid%nx+1,0:curgrid%ny+1), OPTIONAL :: unb
  REAL, DIMENSION(0:curgrid%nx+1,1:curgrid%ny+1), OPTIONAL :: vn
  REAL, DIMENSION(0:curgrid%nx+1,1:curgrid%ny+1), OPTIONAL :: vnb
! number of time steps
  INTEGER :: nt
  INTEGER :: nb
  REAL, DIMENSION(curgrid%nx_1d) :: rhs_xn, xcur
  REAL, DIMENSION(curgrid%nx_1d) :: rhs_xnb, xcurb,rhs_a
  REAL, PARAMETER :: eps=1.d-8
  INTEGER :: nx, ny, nx_1d
  REAL :: dt

  nx = curgrid%nx
  ny = curgrid%ny
  dt = curgrid%dt
  nx_1d = curgrid%nx_1d

  IF (PRESENT(etanb)) THEN
   call trans_2Dto1D(unb,vnb,etanb,xcurb,nx,ny,curgrid%nx_1D,curgrid%indices_1D)
  ELSE
    xcurb = xnb
  END IF
  rhs_xnb = 0.0
  ! call set_verbose(.TRUE.)

  DO nb=nt,1,-1
    CALL GMRES(SHALLOW_WATER_THETA_MATRIX_adjoint, nx_1d, rhs_xnb, &
&            xcurb, itr_max=max_iterations, mr=max(10,nx_1D/100), tol_abs=get_solver_prec(), &
&            tol_rel=get_solver_prec())
    xcurb = rhs_xnb
    rhs_xnb = dt*(1.-theta)*rhs_xnb
    CALL SHALLOW_WATER_MATRIX_adjoint(nx_1d, rhs_xnb,rhs_a)
    xcurb = xcurb + rhs_a
  END DO

  IF (PRESENT(etanb)) THEN
    call trans_1Dto2D(xcurb,unb,vnb,etanb,nx,ny,curgrid%nx_1D,curgrid%indices_1D)
  ELSE
    xnb = xcurb
  END IF
END SUBROUTINE INTEGRATE_SW_adjoint

subroutine shallow_water_matrix_adjoint(n,rhs_xnb,xnb,k)
  integer :: n
  real,dimension(n) :: xnb
  real,dimension(n) :: rhs_xnb
  integer, optional :: k

  real,dimension(1:curgrid%nx+1,0:curgrid%ny+1) :: unb
  real,dimension(0:curgrid%nx+1,1:curgrid%ny+1) :: vnb
  real,dimension(0:curgrid%nx+1,0:curgrid%ny+1) :: etanb
  real, dimension(0:curgrid%nx+1,0:curgrid%ny+1) :: rhs_etab
  real, dimension(1:curgrid%nx+1,0:curgrid%ny+1) :: rhs_ub
  real, dimension(0:curgrid%nx+1,1:curgrid%ny+1) :: rhs_vb

  call trans_1Dto2D(rhs_xnb,rhs_ub,rhs_vb,rhs_etab,curgrid%nx,curgrid%ny,curgrid%nx_1D,curgrid%indices_1D)
  call rhs_sw_b(unb,vnb,etanb,rhs_ub,rhs_vb,rhs_etab,curgrid%nx,curgrid%ny,g,H,coriolis,curgrid%dx,curgrid%dy)
  call trans_2Dto1D(unb,vnb,etanb,xnb,curgrid%nx,curgrid%ny,curgrid%nx_1D,curgrid%indices_1D)

end subroutine shallow_water_matrix_adjoint

subroutine shallow_water_theta_matrix_adjoint(n,rhs_xnb,xnb,k)
  integer :: n
  real,dimension(n) :: xnb
  real,dimension(n) :: rhs_xnb
  integer, optional :: k

  real,dimension(n) :: rhs_a,rhs_b

  xnb = rhs_xnb

  rhs_b = -curgrid%dt*theta*rhs_xnb

  call shallow_water_matrix_adjoint(n,rhs_b,rhs_a)

  xnb = xnb + rhs_a

end subroutine shallow_water_theta_matrix_adjoint
end module shallow_water_adjoint
