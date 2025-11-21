!> Eigenvalues module

!> @brief
!! A series of linear algebra computational routines 


module eigenvalues

  abstract interface
      subroutine matrix_vector(n,x,y,k)
        integer :: n
        integer, optional :: k
        real,dimension(n) :: x,y
      end subroutine matrix_vector
  end interface


  logical :: verbose = .FALSE.
  logical :: verbose_cg = .FALSE.

  real,dimension(:),allocatable :: work_dgesvd

  real,dimension(:),allocatable :: pold_conjgrad

contains

  subroutine set_verbose(lverbose)
    logical :: lverbose

    verbose_cg  = lverbose
  end subroutine set_verbose


  subroutine compute_eigenvalues(matrix_v,n,max_eigenvalue_sw,svd,nb_vectors, &
                 eigenvectors,resid_in,resid_out,SM)
    !include 'debug.h'
    integer :: n
    real,dimension(:) :: max_eigenvalue_sw
    procedure(matrix_vector) :: matrix_v
    logical :: svd
    integer, optional :: nb_vectors
    real,dimension(:,:),optional :: eigenvectors
    real,dimension(n), optional :: resid_in,resid_out
    logical, optional :: SM

    integer :: ldv
    real,dimension(:,:),allocatable :: v, d,z
    real,dimension(:),allocatable :: workl, workd,resid,ax,workev,dr,di
    logical,dimension(:),allocatable :: select
    integer,dimension(11) :: iparam
    integer,dimension(14) :: ipntr

    integer :: ncv,nev,ido,info
    real :: tol,sigmar,sigmai
    character(len=2) :: which
    character(len=1) :: bmat,HOWMNY
    integer :: lworkl
    logical :: rvec
    integer :: i
    integer :: iloop

    ido = 0
    bmat = 'I'
    if (present(SM)) then
      if (SM) THEN
        which = 'SM'
      else
        which = 'LM'
      endif
    else
      which = 'LM'
    endif
    if (svd) then
      nev = min(nb_vectors,n-1)
    else
      if (present(nb_vectors)) then
        nev = min(nb_vectors,n-1)
      else
        nev = min(40,n-1)
      endif
    endif
    tol = 0.
    !tol=1.D-2
    info = 0
    ncv = max(20,min(2*nev+1,n))


    allocate(v(n,ncv),workl(3*ncv*ncv+6*ncv))
    allocate(workd(3*n),resid(n),ax(n))

    lworkl=size(workl)
    ldv=size(v,1)

    iparam=0
    iparam(1) = 1 ! ishift
    iparam(3) = 1000 ! maxitr
    iparam(4) = 1 ! NB
    iparam(7) = 1 ! mode 1


    allocate(select(ncv))


    ! ndigit = -3
    ! logfil = 6
    ! msgets = 0
    ! msaitr = 0
    ! msapps = 0
    ! msaupd = 0
    ! msaup2 = 0
    ! mseigt = 0
    ! mseupd = 0

    if (present(resid_in)) then
      info = 1
      resid = resid_in
    endif

    iloop = 0
    do while (.TRUE.)
     ! print *,'loop'
     ! if (mod(iloop,20) == 0) print *,'iloop = ',iloop
      iloop = iloop+1
      if (svd) THEN
        ! the matrix is symmetric
        call dsaupd ( ido, bmat, n, which, nev, tol, resid,      &
                          ncv, v, n, iparam, ipntr, workd,    &
                          workl, lworkl, info )
      else
        call dnaupd ( ido, bmat, n, which, nev, tol, resid,      &
                      ncv, v, n, iparam, ipntr, workd,    &
                      workl, lworkl, info )
      endif

    if (ido == -1) then
      call matrix_v(n,workd(ipntr(1)), workd(ipntr(2)))
    elseif (ido == 1) then
        workd(ipntr(3):ipntr(3)+n-1)=workd(ipntr(1):ipntr(1)+n-1)
        call matrix_v(n,workd(ipntr(3)), workd(ipntr(2)))
    elseif (ido == 99) then
      exit
    else
      print *,'ido = ',ido
      stop
    endif
  enddo
    print *,'iparam(5) = ',iparam(5)
  if (info/=0) then
    print *,'Caution : Eigenvalues not found (info = ',info,')'
    print *,'Try to increase Nev (Current value = ',Nev,')'
    stop
  endif


  allocate(workev(3*ncv))
  allocate(DR(nev+1),DI(nev+1))
  DR=0.
  DI=0.
  allocate(z(n,nev))

  if (svd) THEN
    rvec = .TRUE.
    HOWMNY='A'
    print *,'ldv = ',ldv
    call dseupd ( rvec, HOWMNY, select, DR(1:nev), z, n,                   &
            sigmar, bmat, n, which, nev, tol,                       &
            resid, ncv, v, ldv, iparam, ipntr, workd,               &
            workl, lworkl, info )
    do i=0,nb_vectors-1
      eigenvectors(:,i+1) = z(:,nev-i)
    enddo
  else
    if (present(nb_vectors)) then
      rvec=.TRUE.
    else
    rvec = .FALSE.
    endif
    HOWMNY  = 'A'
  call dneupd ( rvec, HOWMNY, select, DR, DI, z, n,                 &
          sigmar, sigmai, workev, bmat, n, which, nev, tol,         &
          resid, ncv, v, ldv, iparam, ipntr, workd,                 &
          workl, lworkl, info )
          if (present(nb_vectors)) then
          do i=0,nb_vectors-1
            eigenvectors(:,i+1) = z(:,nev-i)
          enddo
        endif
  endif
  if (info/=0) then
    print *,'Caution (eupd) : Eigenvalues not found (info = ',info,')'
    stop
  endif

  do i=0,nev-1
    print *,'eigen = ',i,sqrt(DR(nev-i)**2+DI(nev-i)**2),DR(nev-i),DI(nev-i)
  enddo

  if (present(nb_vectors)) then
    do i=0,nb_vectors-1
      max_eigenvalue_sw(i+1) = sqrt(DR(nev-i)**2+DI(nev-i)**2)
   enddo
 else
   max_eigenvalue_sw(1) = sqrt(DR(nev)**2+DI(nev)**2)
endif

if (present(resid_out)) then
  resid_out = resid
endif

end subroutine compute_eigenvalues

subroutine conjgrad(matrix_v,b,x,n,eps,observation,exact_matrix,forward_matrix)


  procedure(matrix_vector) :: matrix_v
  procedure(matrix_vector), optional :: exact_matrix, forward_matrix 
  integer :: n
  real,dimension(n) :: b,x,obs
  real :: eps

  real,dimension(n) :: r,Ax,Ap,p,Apold,pold,ptmp
  integer :: iteration

  real :: rsold,rsnew,alpha,tmp,cf,res,tmp2,grad,quad
  real, dimension(:), allocatable :: Aptmp, cftmp

  logical :: exists
  double precision :: dnrm2

  real, dimension(n), optional :: observation
  real :: equad

  call matrix_v(n,x,Ax)

  r = b -Ax

  if (verbose_cg .eqv. .true.) then
    print *,'r = ',DOT_PRODUCT(r,r)
    print *,'bnorm = ',DOT_PRODUCT(b,b),sqrt(DOT_PRODUCT(b,b))
  end if
  p = r
  rsold = DOT_PRODUCT(r,r)

  equad = -46.491068930860784

  iteration = 1
  do while (.TRUE.)

    print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print *, 'CG iteration number ', iteration
    write (*,*)
    call matrix_v(n,p,Ap)

    alpha = rsold / DOT_PRODUCT(p,Ap)
    x = x + alpha * p
    r = r - alpha * Ap
    rsnew = DOT_PRODUCT(r,r)

    if (verbose_cg .eqv. .true.) then
      print *,'ratio = ',rsnew/rsold
      print *,'iter = ',iteration,rsnew,'r2-norm: ',sqrt(rsnew)
    end if
   
    if(PRESENT(observation)) then
      if(PRESENT(forward_matrix)) then
        allocate (cftmp(n))
        call forward_matrix(n,x,cftmp)
        cf = 0.5*dnrm2(n,cftmp - observation,1)
        cf = cf + 0.5*5*(dnrm2(n,x,1))**2
        print *, 'cost function value : ', cf

        deallocate(cftmp)
      end if
    end if 

    !call matrixnorminv(exact_matrix,r,n,tmp)
    !print *, 'A-1 norm of r: ',tmp
    
    if (PRESENT(exact_matrix)) then
      allocate(Aptmp(n))
      call exact_matrix(n,x,Aptmp)

      quad = 0.5*DOT_PRODUCT(x,APtmp) - DOT_PRODUCT(b,x)
      print *, 'quadratic value : ', quad

      tmp =  SQRT(2*(quad - equad))
      print *, 'A-1 norm of r', tmp
      deallocate(Aptmp)
    end if

    inquire(file="ecg_norm.csv",exist = exists)
    if(exists) then
      open(2,file='ecg_norm.csv',form='formatted',status = 'old', position = 'append',action='write')
    else
      open(2,file='ecg_norm.csv',form='formatted',status = 'new', position = 'append',action='write')
    end if

    write(2,'(F15.10,F15.10,F15.10,F15.10)') sqrt(rsnew), cf, quad, tmp
    close(2)

    if (sqrt(rsnew) < eps) exit

    p = r + (rsnew/rsold) * p
    rsold = rsnew

    iteration = iteration + 1
    !print *, 'interation counter changes to iter = ', iteration
  enddo

  if (verbose_cg .eqv. .true.) then
    print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print *,'lastiter = ',iteration,rsnew,sqrt(rsnew)
  end if

end subroutine conjgrad

subroutine gmres( matrix_v, n, x, rhs, itr_max, mr, tol_abs, &
  tol_rel)

!*****************************************************************************80
!
!! MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
!
!  Discussion:
!
!    The linear system A*X=B is solved iteratively.
!
!    The matrix A is assumed to be stored in sparse triplet form.  Only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations
!    to take.  0 < MR <= N.
!
!    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  procedure(matrix_vector) :: matrix_v
  integer ( kind = 4 ) mr
  integer ( kind = 4 ) n
  real ( kind = 8 ) av
  real ( kind = 8 ) c(1:mr)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(1:mr+1)
  real ( kind = 8 ) h(1:mr+1,1:mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j

  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(1:n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ) rhs(1:n)
  real ( kind = 8 ) s(1:mr)
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  real ( kind = 8 ) v(1:n,1:mr+1)
  real ( kind = 8 ) x(1:n)
  real ( kind = 8 ) y(1:mr+1)

  itr_used = 0

  if ( n < mr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
    write ( *, '(a)' ) '  N < MR.'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  MR = ', mr
    stop
  end if

  do itr = 1, itr_max

    call matrix_v(n,x,r)

    r(1:n) = rhs(1:n) - r(1:n)

    rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )

    if ( verbose ) then
      write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( rho <= 1.D-40 ) then
      exit
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call matrix_v(n,v(1:n,k), v(1:n,k+1))

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - h(j,k) * v(1:n,j)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( av + delta * h(k+1,k) == av ) then

        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do

        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then

        y(1:k+1) = h(1:k+1,k)

        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y(1:k+1) )
        end do

        h(1:k+1,k) = y(1:k+1)

      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g(1:k+1) )
      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .or. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .or. rho <= tol_abs ) then
      exit
    end if

  end do

  if (itr>= itr_max) THEN
    print *,'Caution : GMRES performed maximum number of iterations'
  endif
 !print *,'itr = ',itr_used,itr
  if ( verbose ) then
    write ( *, '(a)'       ) ' '
    write ( *, '(a)'       ) 'MGMRES_ST:'
    write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end
subroutine mult_givens ( c, s, k, g )

!*****************************************************************************80
!
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer ( kind = 4 ) K, indicates the location of the first
!    vector entry.
!
!    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer ( kind = 4 ) k

  real ( kind = 8 ) c
  real ( kind = 8 ) g(1:k+1)
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end

subroutine  matrixnorminv(matrix_v,x,n,xout,eps)

integer :: n
real,dimension(n) :: x
real,optional :: eps
procedure(matrix_vector) :: matrix_v

real, intent(out) :: xout
real,dimension(n) :: y,x0
real :: eps1

y = x
x0 = 0.

if (present(eps)) THEN
  eps1=eps
else
  eps1=1.D-1
endif
!call set_verbose(.false.)
!call conjgrad(matrix_v,y,x0,n,eps=eps1)
call cg_for_norm(matrix_v,y,x0,n,eps=eps1)
!call set_verbose(.true.)

xout =sqrt(DOT_PRODUCT(x,x0))

end subroutine matrixnorminv

! inverse of a dense matrix
subroutine matrix_inv(n,A,Ainv)
integer :: n
real,dimension(n,n) :: A, Ainv

integer :: info
integer,dimension(n) :: ipiv
real,dimension(n) :: work

Ainv = A
!print *,'A = ',A

! LU Factotization
call DGETRF(n,n,Ainv,n,ipiv,info)
if (info /= 0 ) then
print *,'info DGETRF = ',info
stop
endif

call DGETRI(n,Ainv,n,ipiv,work,n,info)
if (info /= 0 ) then
print *,'info DGETRI = ',info
stop
endif

end subroutine matrix_inv

! pseudo inverse of a matrix
subroutine matrix_pinv_old(n,nlim,A,Ainv,U,Sdiag)
integer :: n,nlim
real,dimension(n,nlim) :: A, Acopy
real,dimension(nlim,n) :: Ainv

integer :: info

real,dimension(nlim,nlim) :: S
real,dimension(nlim) :: Sdiag
real,dimension(n,nlim) :: U
real,dimension(nlim,nlim) :: VT
integer :: lwork,i
real,dimension(1) :: work_tmp

Acopy = A


lwork=-1
call DGESVD('S','S',n,nlim,Acopy,n,Sdiag,U,n,VT,nlim,work_tmp,lwork,info)

!lwork = 5*n! minimal size
!lwork = lwork * 5
!print *,'WORK = ',lwork,nint(work_tmp(1))
!allocate(work(lwork))

lwork=nint(work_tmp(1))
if (.not.allocated(work_dgesvd)) then
  allocate(work_dgesvd(lwork))
else
  if (size(work_dgesvd) < lwork) then
    deallocate(work_dgesvd)
    allocate(work_dgesvd(lwork))
  endif
endif
call DGESVD('S','S',n,nlim,Acopy,n,Sdiag,U,n,VT,nlim,work_dgesvd,lwork,info)
!print *,'S = ',Sdiag
if (info /=0) then
   print *,'DGESVD info = ',info
   stop
endif

S=0.
do i=1,nlim
!  print *,'ieigen = ',n,nlim,i,Sdiag(i)
  if (abs(Sdiag(i)) >= 1.D-12) THEN
    S(i,i)=1./Sdiag(i)
  endif
enddo
!stop


Ainv = MATMUL(TRANSPOSE(VT),MATMUL(S,TRANSPOSE(U)))
!stop
!print *,'A = ',A
end subroutine matrix_pinv_old

! pseudo inverse of a matrix
subroutine matrix_pinv(n,nlim,A,Ainv)
integer :: n,nlim
real,dimension(n,nlim) :: A, Acopy
real,dimension(nlim,n) :: Ainv

integer :: info
real,dimension(:),allocatable :: work
real,dimension(nlim,nlim) :: S
real,dimension(nlim) :: Sdiag
real,dimension(n,nlim) :: U
real,dimension(nlim,nlim) :: VT
integer :: lwork,i

Acopy = A

lwork = 5*n! minimal size
lwork = lwork * 5
allocate(work(lwork))

call DGESVD('S','S',n,nlim,Acopy,n,Sdiag,U,n,VT,nlim,work,lwork,info)
!print *,'S = ',Sdiag
if (info /=0) then
   print *,'DGESVD info = ',info
   stop
endif

S=0.
do i=1,nlim
 ! print *,'ieigen = ',n,nlim,i,Sdiag(i)
  if (abs(Sdiag(i)) >= 1.D-20) THEN
    S(i,i)=1./Sdiag(i)
  endif
enddo
!stop


Ainv = MATMUL(TRANSPOSE(VT),MATMUL(S,TRANSPOSE(U)))
!stop
!print *,'A = ',A
end subroutine matrix_pinv



subroutine cg_for_norm(matrix_v,b,x,n,eps)

  procedure(matrix_vector) :: matrix_v
  integer :: n
  real,dimension(n) :: b,x
  real :: eps

  real,dimension(n) :: r,Ax,Ap,p,Apold,pold,ptmp,Aptmp
  integer :: iteration

  real :: rsold,rsnew,alpha,tmp

  logical :: exists

  call matrix_v(n,x,Ax)

  r = b -Ax

  p = r
  rsold = DOT_PRODUCT(r,r)


  iteration = 1
  do while (.TRUE.)
    call matrix_v(n,p,Ap)

    alpha = rsold / DOT_PRODUCT(p,Ap)
    x = x + alpha * p
    r = r - alpha * Ap
    rsnew = DOT_PRODUCT(r,r)

    if (sqrt(rsnew) < eps) exit

    p = r + (rsnew/rsold) * p
    rsold = rsnew

    iteration = iteration + 1
    !print *, 'iteration counter changes to iter = ', iteration, 'rnorm = ', sqrt(rsnew)
  enddo
  print *, 'cg iterations to calculate norm = ', iteration

end subroutine cg_for_norm



end module eigenvalues
