!> Parareal module
!> @brief module to provide all necessary subroutines for running parareal in data assimilation

module parareal_utils

use shallow_water
use shallow_water_adjoint
use eigenvalues

integer, parameter :: N_time_windows = 40   !< Number of time windows
integer :: Nfine = 20                       !< Number of fine time steps per time window
integer :: q = 5                           !< Number of coarse time steps per time window

real, parameter :: coarse_solver_prec = 1.d-4   !< Coarse solver precision
logical, parameter :: use_coarse_solver = .true.

logical :: verbose_p = .true.

!################### krylov_enhanced ######################
logical, parameter :: Krylov_Enhanced = .true.  !< Use the Krylov_Enhanced parareal version
! if (debug = .TRUE.) we check that the image of the subspace S are correctly computed
! (Krylov Enhanced case only)
logical, parameter :: debug = .false.
!#########################################################

!################## regularisation #######################
logical, parameter :: regularisation = .true. 
real :: alpharegul = 5 
!#########################################################

! The maximum number of vectors in the subspace S (Krylov Enhanced) is given by decal*N_time_windows
integer, parameter :: decal = 30  !< the shift, used for keeping only latest vectors 

logical, parameter :: incomplete_projection = .true.

real, dimension(:,:,:), allocatable :: lambdank       !< array for storing solution from all the parareal iterations
real, dimension(:,:,:), allocatable :: Flambdank      !< array for storing all fine solver evaluations   
real, dimension(:,:), allocatable :: U_matrix_GS      !< array of orthogonal vectors (basis for subspace S, can call U) 
                                                      !< obtained from the Gram-Schmidt process
real, dimension(:,:), allocatable :: U_matrix_test 

integer :: nb_vectors_limgs   !< number of linearly independent vectors obtained from Gram-Schmidt process

! The max number of possible parareal iterations = N, for krylov enhanced during each iteration the maximum number of vectors we can
! take for the subspace is N + 1(from initial state)
! integer, parameter :: max_nb_vectors = N*(N+1)
! decal is used for keeping only the latest vectors at a particular iteration.
integer, parameter :: max_nb_vectors = (decal+1)*(N_time_windows)   !< maximum number of vectors that can be taken to form the
!< subspace S for the Krylov enhanced version

real, dimension(max_nb_vectors) :: svd_eigenvalues, gain

integer :: nb_vectors     !< number of vectors to take for forming the subspace 

integer :: minvector, nb_times 
integer :: bwindows   !< beginning (index) of time window
integer :: ewindows   !< end (index) of time window

integer :: p_para
integer :: max_p_para   !< maximum parareal iterations
integer :: nx_F

real, dimension(:), allocatable :: lambda_in, lambda_out
!$omp threadprivate(lambda_in,lambda_out)

type(grid), pointer :: Finesolver, Coarsesolver

contains

!> @brief Initialise parareal. Done by defining coarse and fine solvers, allocating memory to compute iterations.
!> @param[in] F fine solver
!> @param[in] G coarse solver
subroutine initialise_parareal(F,G)

  type(grid), pointer :: F, G           ! Fine(F) and Coarse(G) solvers
  F%nt_time_windows = Nfine             ! No. of fine time steps per time window
  G%nt_time_windows = q                 ! No. of coarse time steps per time window
  G%dt = F%dt * real(Nfine)/real(q)     ! Coarse time step

print *, 'size = ', F%nx_1D
print *, 'Fine step size = ', F%dt, 'Coarse step size = ', G%dt

  allocate(lambdank(F%nx_1D,0:N_time_windows,0:N_time_windows+1))   ! allocating maximum parareal iterations (at most N)
  allocate(Flambdank(F%nx_1D,-1:N_time_windows,0:N_time_windows+1)) ! same as above (the index can start from -1)
  allocate(U_matrix_GS(F%nx_1D,max_nb_vectors))
  if (debug) then
    allocate(U_matrix_test(F%nx_1D,max_nb_vectors))
  end if

!$OMP PARALLEL
  allocate(lambda_in(F%nx_1D),lambda_out(F%nx_1D))    ! input - lambda_in, output - lambda_out
!$OMP END PARALLEL

max_p_para = 0
nb_vectors_limgs = 0

end subroutine initialise_parareal


! To print more information in parareal subroutine or not
subroutine set_pverbose(option)
  logical :: option

  verbose_p = option

end subroutine set_pverbose



!> @brief performs the parareal iterations

!> @param[in] F fine solver
!> @param[in] G coarse solver
!> @param[in] lambda0 initial condition
!> @param[out] lambdaout solution
!> @param[in] maxk maximum parareal iterations
!> @param[in] eps parareal tolerance
subroutine Parareal(F,G,lambda0,lambdaout,maxk,eps,iter_out,full_sol,guess)

  type(grid), pointer :: F,G
  real, dimension(F%nx_1D) :: lambda0
  real, dimension(F%nx_1D) :: lambda1,lambda2,lambda3    
  real, dimension(:), allocatable :: lambdaout
  real, dimension(:,:), allocatable, optional :: full_sol
  real, dimension(F%nx_1D), optional :: guess
  integer, optional :: maxk    
  real, optional :: eps       
  integer, intent(out) :: iter_out

  integer :: maxind   ! maximum parareal index
  integer :: p,nw     ! parareal iteration and time window index
  real, dimension(size(lambda0)) :: y1, y2   ! arrays for storing the results of the correction procedure

  real, dimension(size(lambda0),max_nb_vectors) :: Flambdank_s  ! image of basis vectors in S by applying F
  real, dimension(max_nb_vectors) :: C,C2,C3
  real, dimension(max_nb_vectors,max_nb_vectors) :: s,sinv,spinv   ! The subspace S 
  real, dimension(F%nx_1D, max_nb_vectors) :: Fmatrix, lambdamatrix     ! matrix of fine evaluations and matrix of initial
  !conditions
  real, dimension(F%nx_1D,max_nb_vectors) :: Fmatrix_rli, lambdamatrix_rli  ! same matrix with removed linearly independent vectors 
  ! as seen after Gram Schmidt

  real :: err
  real :: t1,t2,t3,t4,t5,t6,t7,t8
  real :: t9,t10,t11,t12,t13,t14
  real :: t15,t16,t17
  real :: t18,t19,t20
  double precision :: dnrm2

  logical :: exists

! Warning: below it is not a comment, it is an instruction compiled only when openmp is turned on (hence the $)
! get the number of threads
!$ integer :: omp_get_thread_num

  integer :: i, j, k

  real, dimension(F%nx_1D) :: max_eigenvalue_sw
  logical :: svd = .true.
  real, dimension(F%nx_1D,100) :: eigenvectors

  Finesolver => F
  Coarsesolver => G

  t17 = 0.
  call my_cpu_time(t1)
  t5 = 0.

  if (present(guess)) then
    maxind = 2
  else
    if (present(maxk)) then
      maxind = maxk+1
    else 
      maxind = N_time_windows 
    end if
  end if

  !if (present(maxk)) then
    !maxind = maxk 
  !else
    !maxind = N_time_windows
  !end if

  !write (*,*) 'maxind = ', maxind

  ! define the size of lambdank --> allocate(lambdank(size(lambda0),0:N_time_windows,maxind+1))
  ! done in the initialise_parareal module
  t20 = 0.

  ! Initialisation by coarse solver
  if (verbose_p .eqv. .true.) then
    print *,'Parareal coarse integration, N = ',N_time_windows
  end if

  p = 0   ! iteration 0

  if (present(guess)) then
    lambdank(:,0,p) = guess
  else
    lambdank(:,0,p) = lambda0
  end if

  !lambdank(:,0,p) = lambda0
  curgrid => G
  ! put the initial state in the subspace
  Flambdank(:,-1,p) = lambda0

  call my_cpu_time(t18)
  call set_solver_prec(coarse_solver_prec)

  ! Find initial configuration (Gx0, G^2x0,...,G^Nx0) over the time windows, nw -> time window no.
  nb_vectors_limgs = 0

  do nw=0,N_time_windows-1
    call integrate_over_time_windows(lambdank(:,nw,p),lambdank(:,nw+1,p)) ! storing in lambdank(:,:,0)
    Flambdank(:,nw,p) = lambdank(:,nw+1,p)  ! also storing in Flambdank
  end do

  call set_solver_prec(1.e-12)
  call my_cpu_time(t19)
  t20 = t20 + t19 - t18

  if (verbose_p .eqv. .true.) then
    print *, 'End coarse initialisation. Timing = ', t20
  end if

  nb_vectors_limgs = 0

  t9 = 0.

! index -- (lambda0 size, time window, parareal iteration)
! #### parareal iterations start #####
  do p=1,maxind-1
    if (verbose_p .eqv. .true.) then
      print *,' Parareal iteration ',p
    end if
    if (present(guess)) then
      lambdank(:,0,p) = guess
      Flambdank(:,-1,p) = guess
    else
      lambdank(:,0,p) = lambda0
      Flambdank(:,-1,p) = lambda0
    end if

    !lambdank(:,0,p) = lambda0
    curgrid => F
    ! put the initial state in the subspace
    !Flambdank(:,-1,p) = lambda0

    ! parallel integration of the fine grid over time windows
    call my_cpu_time(t11)

! Fine integrations in parallel
  curgrid => F


!$OMP PARALLEL
!$OMP DO

    do nw=0,N_time_windows-1
      if (p>=1) then
         call integrate_over_time_windows(lambdank(:,nw,p-1),lambda_out)    ! lambda_out here is F(lambdank(:,nw,p-1)) 
         Flambdank(:,nw,p-1) = lambda_out
      else
        ! mathematically the same step as above, but implementation wise this slight modification improves efficiency
        ! by using fewer gmres iterations for the implicit integration scheme.  
        !!!!!!!!!!!!!!!!!!!!!! changing the way we use less gmres iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call integrate_over_time_windows(lambdank(:,nw,p-1) - lambdank(:,nw,p-2), Flambdank(:,nw,p-1))
         Flambdank(:,nw,p-1) = Flambdank(:,nw,p-1) + Flambdank(:,nw,p-2)  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !call integrate_over_time_windows(lambdank(:,nw,p-1), Flambdank(:,nw,p-1))
      end if
    end do

!$OMP END DO
!$OMP END PARALLEL

    call my_cpu_time(t12)
      t9=t9+t12-t11
    if (verbose_p .eqv. .true.) then
      print *,'Fine integration timing = ',t12-t11
    end if

    if (Krylov_Enhanced) then

      p_para = p-1

      minvector = max(p_para - decal,0)   ! if parareal iterations exceed decal, then we remove (p_para - decal)
      ! number from the current number of parareal iterations; otherwise we take out nothing.

      nb_times = p_para - minvector + 1   ! number of times means the number of parareal iterations we want to include
      ! for determining the total number of vectors in the subspace S and also the size of S. For each parareal iteration 
      ! the size of the subspace changes depending on the value of nb_times

      nx_F = F%nx_1D
      bwindows = 0
      ewindows = N_time_windows - 1
      nb_vectors = (ewindows - bwindows +1)*nb_times    ! total number of vectors in S

      ! lambdamatrix is the S subspace
      lambdamatrix(:,1:nb_vectors) = RESHAPE(lambdank(:,bwindows:ewindows,minvector:p_para),(/nx_F,nb_vectors/))

      ! lambdamatrix is the F(S) subspace
      Fmatrix(:,1:nb_vectors) = RESHAPE(Flambdank(:,bwindows:ewindows,minvector:p_para),(/nx_F,nb_vectors/))

      call my_cpu_time(t18)
      ! construct an orthogonal basis of S (using CGS2 Gram-Schmidt algorithm)
      nb_vectors_limgs = 1
      ! print *, 'Norm U = ',1,dnrm2(nx_F,lambdamatrix(:,1),1)

      U_matrix_GS(:,nb_vectors_limgs) = lambdamatrix(:,nb_vectors_limgs)/dnrm2(nx_F,lambdamatrix(:,1),1)

      lambdamatrix_rli(:,1) = lambdamatrix(:,1)
      Fmatrix_rli(:,1) = Fmatrix(:,1)

      do i=2,nb_vectors
        C2(1:nb_vectors_limgs) = MATMUL(TRANSPOSE(U_matrix_GS(:,1:nb_vectors_limgs)), lambdamatrix(:,i))
        lambda2 = lambdamatrix(:,i) - MATMUL(U_matrix_GS(:,1:nb_vectors_limgs),C2(1:nb_vectors_limgs))

        !!!!!!!!!!!!!!!!!!! second pass of the Gram Schmidt process!!!!!!!!!!!!!!!!!!

        C3(1:nb_vectors_limgs) = MATMUL(TRANSPOSE(U_matrix_GS(:,1:nb_vectors_limgs)),lambda2)
        lambda2 = lambda2 - MATMUL(U_matrix_GS(:,1:nb_vectors_limgs),C3(1:nb_vectors_limgs))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! print *, 'vectors = ', nb_vectors_limgs+1, dnrm2(nx_F,lambda2,1)
        if (dnrm2(nx_F,lambda2,1)>1.D-10) then  ! remove linearly dependent vectors
          nb_vectors_limgs = nb_vectors_limgs + 1
          U_matrix_GS(:,nb_vectors_limgs) = lambda2/dnrm2(nx_F,lambda2,1)
          lambdamatrix_rli(:,nb_vectors_limgs) = lambdamatrix(:,i)
          Fmatrix_rli(:,nb_vectors_limgs) = Fmatrix(:,i)
        end if
      end do


      if (debug) then
        ! Test the orthogonality of the resulting vectors
        t15=0
        t14=0
        do j=1,nb_vectors_limgs
          do i=1,nb_vectors_limgs
            !if ((i/=j) .and. (DOT_PRODUCT(U_matrix_GS(:,i),U_matrix_GS(:,j))>1.D-10)) then 
            if ((i/=j) .and. (DOT_PRODUCT(U_matrix_GS(:,i),U_matrix_GS(:,j))<1.D0)) then 
              !print *, 'error non orthognal basis vectors = ',i,j,DOT_PRODUCT(U_matrix_GS(:,i), U_matrix_GS(:,j))
              !stop
              t14 = abs(DOT_PRODUCT(U_matrix_GS(:,i),U_matrix_GS(:,j)))
              print *, 'values = ',i,j,t14
              if (t14 .ge. t15) then
                t15 = t14
              end if

            end if
          end do
          print *, 'max value',j,'=',t15
          inquire(file="para_error.csv",exist = exists)
          if(exists) then
            open(4,file='para_error.csv',form='formatted',status = 'old', position = 'append',action='write')
          else
            open(4,file='para_error.csv',form='formatted',status = 'new', position = 'append',action='write')
          end if

          write(4,'(F25.18)') t15  
          close(4)
        end do
      end if

      if (verbose_p .eqv. .true.) then
        print *, 'nb_vectors_limgs = ', nb_vectors_limgs, 'total number of vectors = ', nb_vectors
      end if
      
      ! Now compute spinv = (U^TS)^{-1}

      ! U^T S
      s(1:nb_vectors_limgs,1:nb_vectors_limgs) = &
        MATMUL(TRANSPOSE(U_matrix_GS(:,1:nb_vectors_limgs)),lambdamatrix_rli(:,1:nb_vectors_limgs))
      ! spinv = (U^TS)^{-1}
      spinv = 0.
      call matrix_pinv(nb_vectors_limgs,nb_vectors_limgs,s(1:nb_vectors_limgs,1:nb_vectors_limgs), &
        spinv(1:nb_vectors_limgs,1:nb_vectors_limgs))

      if (debug) then
        ! check the inversion (U^TS)^{-1} (U^TS) = Identity

        s(1:nb_vectors_limgs,1:nb_vectors_limgs) = &
          MATMUL(spinv(1:nb_vectors_limgs,1:nb_vectors_limgs),s(1:nb_vectors_limgs,1:nb_vectors_limgs))
      do i=1,nb_vectors_limgs
        do j=1,nb_vectors_limgs
          if (i==j) then
            if (abs(s(i,j)-1) > 1.D-12) then
              print *, 'bug in inversion of (U^TS) i=j : ',i,abs(s(i,j)-1)
            end if
          else
            if (abs(s(i,j)) > 1.D-12) then
              print *, 'bug in inversion of (U^TS) i/=j ',i,j,abs(s(i,j))
            end if
          end if
        end do
      end do

      ! Check that U = S (U^T S)^{-1}
      U_matrix_test(:,1:nb_vectors_limgs) = &
        MATMUL(lambdamatrix_rli(:,1:nb_vectors_limgs),spinv(1:nb_vectors_limgs,1:nb_vectors_limgs))
      do i=1,nb_vectors_limgs
        lambda1 = U_matrix_test(:,i) - U_matrix_GS(:,i)
        ! print *, 'error = ',i,DOT_PRODUCT(lambda1,lambda1)
        if (DOT_PRODUCT(lambda1,lambda1) > 1.D-12) then
          print *, 'The test U = S (U^T S)^{-1} failed ',i,DOT_PRODUCT(lambda1,lambda1)
          stop
        end if
      end do
    end if

    ! compute FU = (FS) (U^TS)^{-1} = (FS) spinv
    if (debug) then
      print *, 'check the correct computation of FU'
    end if

    do i=1,nb_vectors_limgs
      Flambdank_s(:,i) = MATMUL(Fmatrix_rli(:,1:nb_vectors_limgs),spinv(1:nb_vectors_limgs,i))
      if (debug) then
        ! check with a direct integration of the orthogonal vector
        curgrid => F
        ! first check that Fmatrix = F*lambdamatrix
        call integrate_over_time_windows(lambdamatrix_rli(:,i),lambda1)
        lambda2 = lambda1 - Fmatrix_rli(:,i)
        print *, 'vector00 ',i,'error = ',sqrt(DOT_PRODUCT(lambda1,lambda1)), &
          sqrt(DOT_PRODUCT(lambda2,lambda2))/sqrt(DOT_PRODUCT(lambda1,lambda1))

        call integrate_over_time_windows(U_matrix_GS(:,i),lambda1)
        lambda2 = lambda1 - Flambdank_s(:,i)
        print *, 'vector01 ',i,'error = ',sqrt(DOT_PRODUCT(lambda2,lambda2))/sqrt(DOT_PRODUCT(lambda1,lambda1))
        ! Flambdank_s(:,i) = lambda1
      end if
    end do

    call my_cpu_time(t19)
    t5=t5+t19-t18
    if (debug) then
      print *, 'Fin check'
    end if
  end if   ! end if Krylov enhanced



  gain = 0.
  curgrid => G
    do nw=0,N_time_windows-1
      if (nw< p-1) then
        lambdank(:,nw+1,p) = lambdank(:,nw+1,p-1)
        cycle
      end if
      if (Krylov_Enhanced) then
        if (incomplete_projection) then
          lambda1 = lambdank(:,nw,p) - lambdank(:,nw,p-1)
        else
          lambda1 = lambdank(:,nw,p)
        end if

        ! Compute C = U^TX
        C(1:nb_vectors_limgs) = MATMUL(TRANSPOSE(U_matrix_GS(:,1:nb_vectors_limgs)),lambda1)

        !!!!!!! Projection = 0, we don't use Krylov enhanced!!!!!!!
        !C=0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Compute FPX = (FU) (U^TX)
        lambdank(:,nw+1,p) = MATMUL(Flambdank_s(:,1:nb_vectors_limgs),C(1:nb_vectors_limgs))

        if (use_coarse_solver) then   ! updated coarse solver from Krylov_Enhanced
          ! Compute (I-P)*X = X - U U^T X = X - U C

          lambda2 = lambda1

          ! U C

          lambda3 = MATMUL(U_matrix_GS(:,1:nb_vectors_limgs),C(1:nb_vectors_limgs))

          ! For information check the projection against U

          do i =1,nb_vectors_limgs
            ! print *, 'Proj on basis vector ',i,' = ', DOT_PRODUCT(lambda3,U_matrix_GS(:,i))
          end do

          lambda1 = lambda1 - lambda3

          ! print *, 'projection = ',dnrm2(nx_F,lambda1,1)/dnrm2(nx_F,lambda2,1),dnrm2(nx_F,lambda2,1),dnrm2(nx_F,lambda1,1)

          call my_cpu_time(t18)

          call set_solver_prec(coarse_solver_prec)

          curgrid=>G
          ! compute G (I-P)*X

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !if (p .lt. 2) then
            call integrate_over_time_windows(lambda1,y1)
          !else
            !y1=0
          !end if
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call set_solver_prec(1.e-12)
          call my_cpu_time(t19)
          t20 = t20+t19-t18

          ! Compute the coarse correction = FPX + G(I-P)X
          lambdank(:,nw+1,p) = lambdank(:,nw+1,p) + y1
        end if

        if (incomplete_projection) then
          ! add the correction

          lambdank(:,nw+1,p) = lambdank(:,nw+1,p) + Flambdank(:,nw,p-1)
        end if
      else
          call integrate_over_time_windows(lambdank(:,nw,p-1),y1)     ! y1 = G(lambda^{p-1})
          call integrate_over_time_windows(lambdank(:,nw,p),y2)       ! y2 = G(lambda^{p})
          lambdank(:,nw+1,p) = Flambdank(:,nw,p-1) - y1 + y2          
      end if
    end do

    ! using reshape and size to get the error vector
    err = DOT_PRODUCT(RESHAPE(lambdank(:,:,p) - lambdank(:,:,p-1),(/size(lambdank(:,:,p))/)), &
                        RESHAPE(lambdank(:,:,p) - lambdank(:,:,p-1),(/size(lambdank(:,:,p))/)))

    err = sqrt(err)
    print *,'Iteration = ',p,'err = ',err

    !inquire(file="para_error.csv",exist = exists)
    !if(exists) then
       !open(4,file='para_error.csv',form='formatted',status = 'old', position = 'append',action='write')
    !else
       !open(4,file='para_error.csv',form='formatted',status = 'new', position = 'append',action='write')
    !end if

    !write(4,'(F20.14)') err  
    !close(4)


    if (present(eps)) then
      if (err<=eps) exit
    end if
    err = DOT_PRODUCT(RESHAPE(lambdank(:,N_time_windows,p) - lambdank(:,N_time_windows,p-1),&
                (/size(lambdank(:,N_time_windows,p))/)), &
            RESHAPE(lambdank(:,N_time_windows,p) - lambdank(:,N_time_windows,p-1),&
                (/size(lambdank(:,N_time_windows,p))/)))

    err = sqrt(err)
    if (verbose_p .eqv. .true.) then
      print *,'Iteration = ',p,'err2 = ',err
    end if

  end do

  call my_cpu_time(t2)
  if (verbose_p .eqv. .true.) then
    write(*,*) '##################################################'
    print *,'Parareal timing = ',t2-t1
    if (Krylov_Enhanced) then
      print *,'Krylov_Enhanced timing = ',t5
    end if

    print *,'Fine Integration cost = ',t9
    print *,'Coarse Integration cost = ',t20
  end if

  curgrid => F
  if (.not. allocated(lambdaout)) then
    allocate(lambdaout(F%nx_1D))
  end if
  
  if (p==maxind) p = p-1   ! after loop is finished with maxind the counter adds 1 to p 

  if (present(full_sol)) then
    if (.not. allocated(full_sol)) then
      allocate(full_sol(F%nx_1D,p))
      full_sol = lambdank(:,N_time_windows,0:p)
    end if
  end if
  

  lambdaout = lambdank(:,N_time_windows,p)
  iter_out = p

  !inquire(file="ecg_norm.csv",exist = exists)
  !if(exists) then
     !open(1,file='ecg_norm.csv',form='formatted',status = 'old', position = 'append',action='write')
  !else
     !open(1,file='ecg_norm.csv',form='formatted',status = 'new', position = 'append',action='write')
  !end if

  !write(1,"(I2)") iter_out
  !close(1)

  write (*,*) 'Total parareal iterations', iter_out  
  if (verbose_p .eqv. .true.) then
    write(*,*) '##################################################'
  end if 
  write(*,*)
end subroutine Parareal



!> @brief returns the solution y after integrating over one time window of a fine or coarse grid with initial condition x
subroutine integrate_over_time_windows(x,y) 

  ! provide x as the input, returns y as the solution after integrating over curgrid%nt_time_windows with initial condition x
  real, dimension(curgrid%nx_1D) :: x,y

  real :: w1,w2
  !call my_cpu_time(w1)
    y=x

    ! subroutine integrate_sw has some optional inputs and the only required input is nt.
    ! provide 2D arrays or a 1D array (y here)
    call integrate_sw(y,nt=curgrid%nt_time_windows)

  !call my_cpu_time(w2)
  ! print *,'Timing to integrate = ', w2-w1

end subroutine integrate_over_time_windows



!> @brief the exact matrix for the minimisation. Available as a matrix vector product (i.e. y = (F^N)^T*(F^N)*x)
subroutine Mmatrix(nx,x,y,k)

  integer :: nx
  real, dimension(nx) :: x,y,y1

  integer :: i
  integer, optional :: k

  y=x
  ! forward computation FN*x
  do i=0,N_time_windows-1
    ! print *,'I = ',I
    call integrate_sw(xn=y,nt=curgrid%nt_time_windows)
  end do

  ! backward computation FN^T*(FN*x)
  do i=N_time_windows-1,0,-1
    ! print *,'I2 = ',I
    call INTEGRATE_SW_ADJOINT(xnb=y,nt=curgrid%nt_time_windows)
  end do
  
  if (regularisation) then
    call identity_matrix(nx,x,y1)
    y = y + alpharegul*y1
  end if

end subroutine Mmatrix


!> @brief approximation of exact matrix obtained from parareal. Available as matrix vector product (i.e. y = (F^N)^TP(k)*x); k is
! the parareal iteration index
subroutine Mmatrix_parareal(nx,x,y,k)

  integer :: nx
  real, dimension(nx) :: x,y,y1
  ! input - x,  result - y

  real, dimension(:), allocatable :: result
  integer :: i, it_out
  integer, optional :: k

  y=x
  !do i=0,N_time_windows-1
    ! print *,'I = ',I
    ! call integrate_sw(xn=y,nt=curgrid%nt_time_windows)
  !end do
  !print *,'x = ',sum(abs(x))
  !print *,'y = ',sum(abs(y))
  
  if (sum(abs(x)) == 0.) then
    ! if the starting vector is zero the result is zero
    y = 0.
  else
    call set_pverbose(.false.)
    if (PRESENT(k)) then
      call Parareal(FineGrid,CoarseGrid,x,result,maxk=k,iter_out=it_out)
    else
      call Parareal(FineGrid,CoarseGrid,x,result,eps=1.D-1,iter_out=it_out)
    end if
    y = result
  end if
  
  call set_pverbose(.true.)

  do i=N_time_windows-1,0,-1
    ! print *,'I2 = ', I
    call INTEGRATE_SW_ADJOINT(xnb=y,nt=curgrid%nt_time_windows)
  end do
  
  if (regularisation) then
    call identity_matrix(nx,x,y1)
    y = y + alpharegul*y1
  end if

end subroutine Mmatrix_parareal

subroutine Fmatrix(nx,x,y,k)

  integer :: nx
  real, dimension(nx) :: x,y,y1

  integer :: i
  integer, optional :: k

  y=x
  ! forward computation FN*x
  do i=0,N_time_windows-1
    ! print *,'I = ',I
    call integrate_sw(xn=y,nt=curgrid%nt_time_windows)
  end do

end subroutine Fmatrix

subroutine Ematrix(nx,x,y,k)

  integer :: nx
  integer, optional :: k
  real, dimension(nx) :: x,y
  real, dimension(nx) :: xout,xout1
  
  call set_pverbose(.false.)
  call Mmatrix(nx,x,xout,k)
  call Mmatrix_parareal(nx,x,xout1,k)

  y = xout1 - xout

end subroutine


!>  RHS vector (i.e. y = (F^N)^T*x); x is the true observation
subroutine Bvector(nx,x,y)

  integer :: nx
  real, dimension(nx) :: x,y

  integer :: i
  
  y=x
  do i=N_time_windows-1,0,-1
    call INTEGRATE_SW_ADJOINT(xnb=y,nt=curgrid%nt_time_windows)
  end do

end subroutine Bvector


!Define the identity matrix, x -> Ix
subroutine identity_matrix(nx,x,y)
  
  integer, dimension(nx,nx) :: id_mat
  real, dimension(nx) :: x,y
  integer :: i, nx
  
  ! penalising only the velocities
  id_mat = 0
  forall(i=curgrid%nx*curgrid%ny+1:nx) id_mat(i,i) = 1
  
  y = MATMUL(id_mat,x)
  
end subroutine identity_matrix



end module parareal_utils
