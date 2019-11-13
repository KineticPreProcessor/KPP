!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Adjoint of SDIRK - Singly-Diagonally-Implicit Runge-Kutta method       !
!            * Sdirk 2a, 2b: L-stable, 2 stages, order 2                  !
!            * Sdirk 3a:     L-stable, 3 stages, order 2, adj-invariant   !
!            * Sdirk 4a, 4b: L-stable, 5 stages, order 4                  !
!  By default the code employs the KPP sparse linear algebra routines     !
!  Compile with -DFULL_ALGEBRA to use full linear algebra (LAPACK)        !
!                                                                         !
!    (C)  Adrian Sandu, July 2005                                         !
!    Virginia Polytechnic Institute and State University                  !
!    Contact: sandu@cs.vt.edu                                             !
!    Revised by Philipp Miehe and Adrian Sandu, May 2006                  !
!    This implementation is part of KPP - the Kinetic PreProcessor        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE KPP_ROOT_Integrator

  USE KPP_ROOT_Precision
  USE KPP_ROOT_Global, ONLY: FIX, RCONST, TIME
  USE KPP_ROOT_Parameters, ONLY: NVAR, NSPEC, NFIX, LU_NONZERO
  USE KPP_ROOT_JacobianSP, ONLY: LU_DIAG
  USE KPP_ROOT_Jacobian, ONLY: Jac_SP_Vec, JacTR_SP_Vec
  USE KPP_ROOT_LinearAlgebra, ONLY: KppDecomp, KppSolve,    &
               KppSolveTR, Set2zero, WLAMCH, WCOPY, WAXPY, WSCAL, WADD
  
  IMPLICIT NONE
  PUBLIC
  SAVE
  
!~~~>  Statistics on the work performed by the SDIRK method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4,  &
           Nrej=5, Ndec=6, Nsol=7, Nsng=8,               &
           Ntexit=1, Nhexit=2, Nhnew=3
                 
CONTAINS

SUBROUTINE INTEGRATE_ADJ( NADJ, Y, Lambda, TIN, TOUT, &
           ATOL_adj, RTOL_adj,                        &
           ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U )

   USE KPP_ROOT_Parameters
   USE KPP_ROOT_Global
   IMPLICIT NONE

!~~~> Y - Concentrations
   KPP_REAL :: Y(NVAR)
!~~~> NADJ - No. of cost functionals for which adjoints
!                are evaluated simultaneously
!            If single cost functional is considered (like in
!                most applications) simply set NADJ = 1      
   INTEGER :: NADJ
!~~~> Lambda - Sensitivities w.r.t. concentrations
!     Note: Lambda (1:NVAR,j) contains sensitivities of
!           the j-th cost functional w.r.t. Y(1:NVAR), j=1...NADJ
   KPP_REAL  :: Lambda(NVAR,NADJ)
!~~~> Tolerances for adjoint calculations
!     (used for full continuous adjoint, and for controlling
!      iterations when used to solve the discrete adjoint)   
   KPP_REAL, INTENT(IN)  :: ATOL_adj(NVAR,NADJ), RTOL_adj(NVAR,NADJ)
   KPP_REAL, INTENT(IN) :: TIN  ! Start Time
   KPP_REAL, INTENT(IN) :: TOUT ! End Time
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   KPP_REAL, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   KPP_REAL, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: Ierr_U

   INTEGER, SAVE :: Ntotal = 0
   KPP_REAL :: RCNTRL(20), RSTATUS(20), T1, T2
   INTEGER       :: ICNTRL(20), ISTATUS(20), Ierr

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp

   !~~~> fine-tune the integrator:
   ICNTRL(5) = 8    ! Max no. of Newton iterations
   ICNTRL(7) = 1    ! Adjoint solution by: 0=Newton, 1=direct
   ICNTRL(8) = 1    ! Save fwd LU factorization: 0 = do *not* save, 1 = save

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
   
   T1 = TIN; T2 = TOUT
   CALL SDIRKADJ( NVAR, NADJ, T1, T2, Y, Lambda,            &
                  RTOL, ATOL, ATOL_adj, RTOL_adj,           &
                  RCNTRL, ICNTRL, RSTATUS, ISTATUS, Ierr )

   !~~~> Debug option: number of steps
   ! Ntotal = Ntotal + ISTATUS(Nstp)
   ! WRITE(6,777) ISTATUS(Nstp),Ntotal,VAR(ind_O3),VAR(ind_NO2)
   ! 777 FORMAT('NSTEPS=',I5,' (',I5,')  O3=',E24.14,'  NO2=',E24.14)    

   IF (Ierr < 0) THEN
        PRINT *,'SDIRK: Unsuccessful exit at T=',TIN,' (Ierr=',Ierr,')'
   ENDIF
   
   ! if optional parameters are given for output they to return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(Ierr_U))    Ierr_U       = Ierr

   END SUBROUTINE INTEGRATE_ADJ


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SDIRKADJ(N, NADJ, Tinitial, Tfinal, Y, Lambda,   &
                       RelTol, AbsTol, RelTol_adj, AbsTol_adj, &
                       RCNTRL, ICNTRL, RSTATUS, ISTATUS, Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using a Singly-Diagonally-Implicit
!    Runge-Kutta (SDIRK) method.
!
!    This implementation is based on the book and the code Sdirk4:
!
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.
!    This code is based on the SDIRK4 routine in the above book.
!
!    Methods:
!            * Sdirk 2a, 2b: L-stable, 2 stages, order 2                  
!            * Sdirk 3a:     L-stable, 3 stages, order 2, adjoint-invariant   
!            * Sdirk 4a, 4b: L-stable, 5 stages, order 4                  
!
!    (C)  Adrian Sandu, July 2005
!    Virginia Polytechnic Institute and State University
!    Contact: sandu@cs.vt.edu
!    Revised by Philipp Miehe and Adrian Sandu, May 2006                  
!    This implementation is part of KPP - the Kinetic PreProcessor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!-     Y(NVAR)    = vector of initial conditions (at T=Tinitial)
!-    [Tinitial,Tfinal]  = time range of integration
!     (if Tinitial>Tfinal the integration is performed backwards in time)
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE ode_Fun( T, Y, Ydot ) = ODE function,
!                       returns Ydot = Y' = F(T,Y)
!- SUBROUTINE ode_Fun( T, Y, Ydot ) = Jacobian of the ODE function,
!                       returns Jcb = dF/dY
!-    ICNTRL(1:20)    = integer inputs parameters
!-    RCNTRL(1:20)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!-    Y(NVAR)         -> vector of final states (at T->Tfinal)
!-    ISTATUS(1:20)   -> integer output parameters
!-    RSTATUS(1:20)   -> real output parameters
!-    Ierr            -> job status upon return
!                        success (positive value) or
!                        failure (negative value)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    Note: For input parameters equal to zero the default values of the
!          corresponding variables are used.
!~~~>  
!    ICNTRL(1) = not used
!
!    ICNTRL(2) = 0: AbsTol, RelTol are NVAR-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3) = Method
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0 the default value of 1500 is used
!        Note: use a conservative estimate, since the checkpoint
!              buffers are allocated to hold Max_no_steps
!
!    ICNTRL(5)  -> maximum number of Newton iterations
!        For ICNTRL(5)=0 the default value of 8 is used
!
!    ICNTRL(6)  -> starting values of Newton iterations:
!        ICNTRL(6)=0 : starting values are interpolated (the default)
!        ICNTRL(6)=1 : starting values are zero
!
!    ICNTRL(7)  -> method to solve ADJ equations:
!        ICNTRL(7)=0 : modified Newton re-using LU (the default)
!        ICNTRL(7)=1 : direct solution (additional one LU factorization per stage)
!
!    ICNTRL(8)  -> checkpointing the LU factorization at each step:
!        ICNTRL(8)=0 : do *not* save LU factorization (the default)
!        ICNTRL(8)=1 : save LU factorization
!        Note: if ICNTRL(7)=1 the LU factorization is *not* saved
!
!~~~>  Real parameters
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!                  It is strongly recommended to keep Hmin = ZERO
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                 (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!                  than the predicted value  (default=0.9)
!    RCNTRL(8)  -> ThetaMin. If Newton convergence rate smaller
!                  than ThetaMin the Jacobian is not recomputed;
!                  (default=0.001)
!    RCNTRL(9)  -> NewtonTol, stopping criterion for Newton's method
!                  (default=0.03)
!    RCNTRL(10) -> Qmin
!    RCNTRL(11) -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
!                  step size is kept constant and the LU factorization
!                  reused (default Qmin=1, Qmax=1.2)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call to Rosenbrock adds the current no. of fcn calls
!      to previous value of ISTATUS(1), and similar for the other params.
!      Set ISTATUS(1:10) = 0 before call to avoid this accumulation.
!
!    ISTATUS(1) = No. of function calls
!    ISTATUS(2) = No. of jacobian calls
!    ISTATUS(3) = No. of steps
!    ISTATUS(4) = No. of accepted steps
!    ISTATUS(5) = No. of rejected steps (except at the beginning)
!    ISTATUS(6) = No. of LU decompositions
!    ISTATUS(7) = No. of forward/backward substitutions
!    ISTATUS(8) = No. of singular matrix decompositions
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                   computed Y upon return
!    RSTATUS(2)  -> Hexit,last accepted step before return
!    RSTATUS(3)  -> Hnew, last predicted step before return
!        For multiple restarts, use Hnew as Hstart in the following run
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE

! Arguments      
      INTEGER, INTENT(IN)          :: N, NADJ, ICNTRL(20)
      KPP_REAL, INTENT(INOUT) :: Y(NVAR), Lambda(NVAR,NADJ)
      KPP_REAL, INTENT(IN)    :: Tinitial, Tfinal, &
                    RelTol(NVAR), AbsTol(NVAR), RCNTRL(20), &
                    RelTol_adj(NVAR,NADJ), AbsTol_adj(NVAR,NADJ)
      INTEGER, INTENT(OUT)         :: Ierr
      INTEGER, INTENT(INOUT)       :: ISTATUS(20) 
      KPP_REAL, INTENT(OUT)   :: RSTATUS(20)
       
!~~~>  SDIRK method coefficients, up to 5 stages
      INTEGER, PARAMETER :: Smax = 5
      INTEGER, PARAMETER :: S2A=1, S2B=2, S3A=3, S4A=4, S4B=5
      KPP_REAL :: rkGamma, rkA(Smax,Smax), rkB(Smax), rkC(Smax), &
                       rkD(Smax),  rkE(Smax), rkBhat(Smax), rkELO,    &
                       rkAlpha(Smax,Smax), rkTheta(Smax,Smax)
      INTEGER :: sdMethod, rkS ! The number of stages

!~~~>  Checkpoints in memory buffers
      INTEGER :: stack_ptr = 0 ! last written entry in checkpoint
      KPP_REAL, DIMENSION(:),     POINTER :: chk_H, chk_T
      KPP_REAL, DIMENSION(:,:),   POINTER :: chk_Y
      KPP_REAL, DIMENSION(:,:,:), POINTER :: chk_Z
      INTEGER,       DIMENSION(:,:),   POINTER :: chk_P
#ifdef FULL_ALGEBRA
      KPP_REAL, DIMENSION(:,:,:), POINTER :: chk_J
#else
      KPP_REAL, DIMENSION(:,:),   POINTER :: chk_J
#endif
! Local variables      
      KPP_REAL :: Hmin, Hmax, Hstart, Roundoff,    &
                       FacMin, Facmax, FacSafe, FacRej, &
                       ThetaMin, NewtonTol, Qmin, Qmax
      LOGICAL       :: SaveLU, DirectADJ                 
      INTEGER       :: ITOL, NewtonMaxit, Max_no_steps, i
      KPP_REAL, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
      KPP_REAL, PARAMETER :: DeltaMin = 1.0d-5

       
!~~~>  Initialize statistics
      ISTATUS(1:20) = 0
      RSTATUS(1:20) = ZERO
      Ierr          = 0

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:NVAR) and RelTol(1:NVAR)
      IF (ICNTRL(2) == 0) THEN
         ITOL = 1
      ELSE
         ITOL = 0
      END IF

!~~~> ICNTRL(3) - method selection       
      SELECT CASE (ICNTRL(3))
      CASE (0,1)
         CALL Sdirk2a
      CASE (2)
         CALL Sdirk2b
      CASE (3)
         CALL Sdirk3a
      CASE (4)
         CALL Sdirk4a
      CASE (5)
         CALL Sdirk4b
      CASE DEFAULT
         CALL Sdirk2a
      END SELECT
      
!~~~>   The maximum number of time steps admitted
      IF (ICNTRL(4) == 0) THEN
         Max_no_steps = 200000
      ELSEIF (ICNTRL(4) > 0) THEN
         Max_no_steps = ICNTRL(4)
      ELSE
         PRINT * ,'User-selected ICNTRL(4)=',ICNTRL(4)
         CALL SDIRK_ErrorMsg(-1,Tinitial,ZERO,Ierr)
   END IF
 
!~~~> The maximum number of Newton iterations admitted
      IF(ICNTRL(5) == 0)THEN
         NewtonMaxit=8
      ELSE
         NewtonMaxit=ICNTRL(5)
         IF(NewtonMaxit <= 0)THEN
             PRINT * ,'User-selected ICNTRL(5)=',ICNTRL(5)
             CALL SDIRK_ErrorMsg(-2,Tinitial,ZERO,Ierr)
         END IF
      END IF

!~~~> Solve ADJ equations directly or by Newton iterations 
      DirectADJ = (ICNTRL(7) == 1)
 
!~~~> Save or not the forward LU factorization
      SaveLU = (ICNTRL(8) /= 0) .AND. (.NOT.DirectADJ)

!~~~>  Unit roundoff (1+Roundoff>1)
      Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
      IF (RCNTRL(1) == ZERO) THEN
         Hmin = ZERO
      ELSEIF (RCNTRL(1) > ZERO) THEN
         Hmin = RCNTRL(1)
      ELSE
         PRINT * , 'User-selected RCNTRL(1)=', RCNTRL(1)
         CALL SDIRK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Upper bound on the step size: (positive value)
      IF (RCNTRL(2) == ZERO) THEN
         Hmax = ABS(Tfinal-Tinitial)
      ELSEIF (RCNTRL(2) > ZERO) THEN
         Hmax = MIN(ABS(RCNTRL(2)),ABS(Tfinal-Tinitial))
      ELSE
         PRINT * , 'User-selected RCNTRL(2)=', RCNTRL(2)
         CALL SDIRK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Starting step size: (positive value)
      IF (RCNTRL(3) == ZERO) THEN
         Hstart = MAX(Hmin,Roundoff)
      ELSEIF (RCNTRL(3) > ZERO) THEN
         Hstart = MIN(ABS(RCNTRL(3)),ABS(Tfinal-Tinitial))
      ELSE
         PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
         CALL SDIRK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax
      IF (RCNTRL(4) == ZERO) THEN
         FacMin = 0.2_dp
      ELSEIF (RCNTRL(4) > ZERO) THEN
         FacMin = RCNTRL(4)
      ELSE
         PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
      IF (RCNTRL(5) == ZERO) THEN
         FacMax = 10.0_dp
      ELSEIF (RCNTRL(5) > ZERO) THEN
         FacMax = RCNTRL(5)
      ELSE
         PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
      IF (RCNTRL(6) == ZERO) THEN
         FacRej = 0.1_dp
      ELSEIF (RCNTRL(6) > ZERO) THEN
         FacRej = RCNTRL(6)
      ELSE
         PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
      IF (RCNTRL(7) == ZERO) THEN
         FacSafe = 0.9_dp
      ELSEIF (RCNTRL(7) > ZERO) THEN
         FacSafe = RCNTRL(7)
      ELSE
         PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF

!~~~> ThetaMin: decides whether the Jacobian should be recomputed
      IF(RCNTRL(8) == 0.D0)THEN
         ThetaMin = 1.0d-3
      ELSE
         ThetaMin = RCNTRL(8)
      END IF

!~~~> Stopping criterion for Newton's method
      IF(RCNTRL(9) == ZERO)THEN
         NewtonTol = 3.0d-2
      ELSE
         NewtonTol = RCNTRL(9)
      END IF

!~~~> Qmin, Qmax: IF Qmin < Hnew/Hold < Qmax, STEP SIZE = CONST.
      IF(RCNTRL(10) == ZERO)THEN
         Qmin=ONE
      ELSE
         Qmin=RCNTRL(10)
      END IF
      IF(RCNTRL(11) == ZERO)THEN
         Qmax=1.2D0
      ELSE
         Qmax=RCNTRL(11)
      END IF

!~~~>  Check if tolerances are reasonable
      IF (ITOL == 0) THEN
         IF (AbsTol(1) <= ZERO .OR. RelTol(1) <= 10.D0*Roundoff) THEN
            PRINT * , ' Scalar AbsTol = ',AbsTol(1)
            PRINT * , ' Scalar RelTol = ',RelTol(1)
            CALL SDIRK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
         END IF
      ELSE
         DO i=1,N
            IF (AbsTol(i) <= 0.D0.OR.RelTol(i) <= 10.D0*Roundoff) THEN
              PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
              PRINT * , ' RelTol(',i,') = ',RelTol(i)
              CALL SDIRK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
            END IF
         END DO
      END IF
    
    IF (Ierr < 0) RETURN
    
!~~~>  Allocate memory buffers    
    CALL SDIRK_AllocBuffers

!~~~>  Call forward integration    
    CALL SDIRK_FwdInt( N, Tinitial, Tfinal, Y, Ierr )

!~~~>  Call adjoint integration    
    CALL SDIRK_DadjInt( N, NADJ, Lambda, Ierr )

!~~~>  Free memory buffers    
    CALL SDIRK_FreeBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CONTAINS  !  Procedures internal to SDIRKADJ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SDIRK_FwdInt( N,Tinitial,Tfinal,Y,Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      USE KPP_ROOT_Parameters
      IMPLICIT NONE

!~~~> Arguments:      
      INTEGER :: N
      KPP_REAL, INTENT(INOUT) :: Y(NVAR)
      KPP_REAL, INTENT(IN) :: Tinitial, Tfinal
      INTEGER, INTENT(OUT) :: Ierr
      
!~~~> Local variables:      
      KPP_REAL :: Z(NVAR,Smax), G(NVAR), TMP(NVAR),        &   
                       NewtonRate, SCAL(NVAR), RHS(NVAR),       &
                       T, H, Theta, Hratio, NewtonPredictedErr, &
                       Qnewton, Err, Fac, Hnew,Tdirection,      &
                       NewtonIncrement, NewtonIncrementOld
      INTEGER :: j, IER, istage, NewtonIter, IP(NVAR)
      LOGICAL :: Reject, FirstStep, SkipJac, SkipLU, NewtonDone
      
#ifdef FULL_ALGEBRA      
      KPP_REAL, DIMENSION(NVAR,NVAR) :: FJAC, E
#else      
      KPP_REAL, DIMENSION(LU_NONZERO):: FJAC, E
#endif      


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initializations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      T = Tinitial
      Tdirection = SIGN(ONE,Tfinal-Tinitial)
      H = MAX(ABS(Hmin),ABS(Hstart))
      IF (ABS(H) <= 10.D0*Roundoff) H=1.0D-6
      H=MIN(ABS(H),Hmax)
      H=SIGN(H,Tdirection)
      SkipLU =.FALSE.
      SkipJac = .FALSE.
      Reject=.FALSE.
      FirstStep=.TRUE.

      CALL SDIRK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( (Tfinal-T)*Tdirection - Roundoff > ZERO )


!~~~>  Compute E = 1/(h*gamma)-Jac and its LU decomposition
      IF ( .NOT.SkipLU ) THEN ! This time around skip the Jac update and LU
         CALL SDIRK_PrepareMatrix ( H, T, Y, FJAC, &
                   SkipJac, SkipLU, E, IP, Reject, IER )
         IF (IER /= 0) THEN
             CALL SDIRK_ErrorMsg(-8,T,H,Ierr); RETURN
         END IF
      END IF      

      IF (ISTATUS(Nstp) > Max_no_steps) THEN
             CALL SDIRK_ErrorMsg(-6,T,H,Ierr); RETURN
      END IF   
      IF ( (T+0.1d0*H == T) .OR. (ABS(H) <= Roundoff) ) THEN
             CALL SDIRK_ErrorMsg(-7,T,H,Ierr); RETURN
      END IF   

stages:DO istage = 1, rkS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Simplified Newton iterations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~>  Starting values for Newton iterations
       CALL Set2zero(N,Z(1,istage))
       
!~~~>   Prepare the loop-independent part of the right-hand side
       CALL Set2zero(N,G)
       IF (istage > 1) THEN
           DO j = 1, istage-1
               ! Gj(:) = sum_j Theta(i,j)*Zj(:) = H * sum_j A(i,j)*Fun(Zj(:))
               CALL WAXPY(N,rkTheta(istage,j),Z(1,j),1,G,1)
               ! Zi(:) = sum_j Alpha(i,j)*Zj(:)
               CALL WAXPY(N,rkAlpha(istage,j),Z(1,j),1,Z(1,istage),1)
           END DO
       END IF

       !~~~>  Initializations for Newton iteration
       NewtonDone = .FALSE.
       Fac = 0.5d0 ! Step reduction factor if too many iterations
            
NewtonLoop:DO NewtonIter = 1, NewtonMaxit

!~~~>   Prepare the loop-dependent part of the right-hand side
 	    CALL WADD(N,Y,Z(1,istage),TMP)         	! TMP <- Y + Zi
            CALL FUN_CHEM(T+rkC(istage)*H,TMP,RHS)	! RHS <- Fun(Y+Zi)
            ISTATUS(Nfun) = ISTATUS(Nfun) + 1
!            RHS(1:N) = G(1:N) - Z(1:N,istage) + (H*rkGamma)*RHS(1:N)
	    CALL WSCAL(N, H*rkGamma, RHS, 1)
	    CALL WAXPY (N, -ONE, Z(1,istage), 1, RHS, 1)
            CALL WAXPY (N, ONE, G,1, RHS,1)

!~~~>   Solve the linear system
            CALL SDIRK_Solve ( 'N', H, N, E, IP, IER, RHS )
            
!~~~>   Check convergence of Newton iterations
            CALL SDIRK_ErrorNorm(N, RHS, SCAL, NewtonIncrement)
            IF ( NewtonIter == 1 ) THEN
                Theta      = ABS(ThetaMin)
                NewtonRate = 2.0d0 
            ELSE
                Theta = NewtonIncrement/NewtonIncrementOld
                IF (Theta < 0.99d0) THEN
                    NewtonRate = Theta/(ONE-Theta)
                    ! Predict error at the end of Newton process 
                    NewtonPredictedErr = NewtonIncrement &
                               *Theta**(NewtonMaxit-NewtonIter)/(ONE-Theta)
                    IF (NewtonPredictedErr >= NewtonTol) THEN
                      ! Non-convergence of Newton: predicted error too large
                      Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                      Fac = 0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                      EXIT NewtonLoop
                    END IF
                ELSE ! Non-convergence of Newton: Theta too large
                    EXIT NewtonLoop
                END IF
            END IF
            NewtonIncrementOld = NewtonIncrement
            ! Update solution: Z(:) <-- Z(:)+RHS(:)
            CALL WAXPY(N,ONE,RHS,1,Z(1,istage),1) 
            
            ! Check error in Newton iterations
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
            IF (NewtonDone) EXIT NewtonLoop
            
            END DO NewtonLoop
            
            IF (.NOT.NewtonDone) THEN
                 !CALL RK_ErrorMsg(-12,T,H,Ierr);
                 H = Fac*H; Reject=.TRUE.
                 SkipJac = .TRUE.; SkipLU = .FALSE.
                 CYCLE Tloop
            END IF

!~~~>  End of implified Newton iterations

 
   END DO stages


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Error estimation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      CALL Set2zero(N,TMP)
      DO i = 1,rkS
         IF (rkE(i)/=ZERO) CALL WAXPY(N,rkE(i),Z(1,i),1,TMP,1)
      END DO  

      CALL SDIRK_Solve( 'N', H, N, E, IP, IER, TMP )
      CALL SDIRK_ErrorNorm(N, TMP, SCAL, Err)

!~~~> Computation of new step size Hnew
      Fac  = FacSafe*(Err)**(-ONE/rkELO)
      Fac  = MAX(FacMin,MIN(FacMax,Fac))
      Hnew = H*Fac

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Accept/Reject step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
accept: IF ( Err < ONE ) THEN !~~~> Step is accepted

         FirstStep=.FALSE.
         ISTATUS(Nacc) = ISTATUS(Nacc) + 1

!~~~> Checkpoint solution
         CALL SDIRK_Push( T, H, Y, Z, E, IP )

!~~~> Update time and solution
         T  =  T + H
         ! Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
         DO i = 1,rkS 
            IF (rkD(i)/=ZERO) CALL WAXPY(N,rkD(i),Z(1,i),1,Y,1)
         END DO  
       
!~~~> Update scaling coefficients
         CALL SDIRK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)

!~~~> Next time step
         Hnew = Tdirection*MIN(ABS(Hnew),Hmax)
         ! Last T and H
         RSTATUS(Ntexit) = T
         RSTATUS(Nhexit) = H
         RSTATUS(Nhnew)  = Hnew
         ! No step increase after a rejection
         IF (Reject) Hnew = Tdirection*MIN(ABS(Hnew),ABS(H))
         Reject = .FALSE.
         IF ((T+Hnew/Qmin-Tfinal)*Tdirection > ZERO) THEN
            H = Tfinal-T
         ELSE
            Hratio=Hnew/H
            ! If step not changed too much keep Jacobian and reuse LU
            SkipLU = ( (Theta <= ThetaMin) .AND. (Hratio >= Qmin) &
                                     .AND. (Hratio <= Qmax) )
            IF (.NOT.SkipLU) H = Hnew
         END IF
         ! If convergence is fast enough, do not update Jacobian
!         SkipJac = (Theta <= ThetaMin)
         SkipJac = .FALSE.

      ELSE accept !~~~> Step is rejected

         IF (FirstStep .OR. Reject) THEN
             H = FacRej*H
         ELSE
             H = Hnew
         END IF
         Reject  = .TRUE.
         SkipJac = .TRUE.
         SkipLU  = .FALSE. 
         IF (ISTATUS(Nacc) >= 1) ISTATUS(Nrej) = ISTATUS(Nrej) + 1
         
      END IF accept
      
      END DO Tloop

      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE SDIRK_FwdInt



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SDIRK_DadjInt( N, NADJ, Lambda, Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      USE KPP_ROOT_Parameters
      IMPLICIT NONE

!~~~> Arguments:      
      INTEGER, INTENT(IN) :: N, NADJ
      KPP_REAL, INTENT(INOUT) :: Lambda(NVAR,NADJ)
      INTEGER, INTENT(OUT) :: Ierr
      
!~~~> Local variables:      
      KPP_REAL :: Y(NVAR)
      KPP_REAL :: Z(NVAR,Smax), U(NVAR,NADJ,Smax),   &
                       TMP(NVAR), G(NVAR),                &   
                       NewtonRate, SCAL(NVAR), DU(NVAR),  &
                       T, H, Theta, NewtonPredictedErr,   &
                       NewtonIncrement, NewtonIncrementOld
      INTEGER :: j, IER, istage, iadj, NewtonIter, &
                 IP(NVAR), IP_adj(NVAR)
      LOGICAL :: Reject, SkipJac, SkipLU, NewtonDone
      
#ifdef FULL_ALGEBRA      
      KPP_REAL, DIMENSION(NVAR,NVAR) :: E, Jac, E_adj
#else      
      KPP_REAL, DIMENSION(LU_NONZERO):: E, Jac, E_adj
#endif      


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( stack_ptr > 0 )
        
   !~~~>  Recover checkpoints for stage values and vectors
      CALL SDIRK_Pop( T, H, Y, Z, E, IP )

!~~~>  Compute E = 1/(h*gamma)-Jac and its LU decomposition
      IF (.NOT.SaveLU) THEN
          SkipJac = .FALSE.; SkipLU = .FALSE.
          CALL SDIRK_PrepareMatrix ( H, T, Y, Jac, &
                   SkipJac, SkipLU, E, IP, Reject, IER )
         IF (IER /= 0) THEN
             CALL SDIRK_ErrorMsg(-8,T,H,Ierr); RETURN
         END IF
      END IF      

stages:DO istage = rkS, 1, -1

!~~~>  Jacobian of the current stage solution
       TMP(1:N) = Y(1:N) + Z(1:N,istage)
       CALL JAC_CHEM(T+rkC(istage)*H,TMP,Jac)
       ISTATUS(Njac) = ISTATUS(Njac) + 1
       
       IF (DirectADJ) THEN 
#ifdef FULL_ALGEBRA
         E_adj(1:N,1:N) = -Jac(1:N,1:N)
         DO j=1,N
            E_adj(j,j) = E_adj(j,j) + ONE/(H*rkGamma)
         END DO
         CALL DGETRF( N, N, E_adj, N, IP_adj, IER )
#else
         E_adj(1:LU_NONZERO) = -Jac(1:LU_NONZERO)
         DO i = 1,NVAR
            j = LU_DIAG(i); E_adj(j) = E_adj(j) + ONE/(H*rkGamma)
         END DO
         CALL KppDecomp ( E_adj, IER)
#endif
         ISTATUS(Ndec) = ISTATUS(Ndec) + 1
         IF (IER /= 0)  THEN
             PRINT*,'At stage ',istage,' the matrix used in adjoint', &
                      ' computation is singular'
             CALL SDIRK_ErrorMsg(-8,T,H,Ierr); RETURN
         END IF
       END IF

adj:   DO iadj = 1, NADJ
       
!~~~> Update scaling coefficients
       CALL SDIRK_ErrorScale(ITOL, AbsTol_adj(1:NVAR,iadj), &
             RelTol_adj(1:NVAR,iadj), Lambda(1:NVAR,iadj), SCAL)
      
!~~~>   Prepare the loop-independent part of the right-hand side
!       G(:) = H*Jac^T*( B(i)*Lambda + sum_j A(j,i)*Uj(:) )
       G(1:N) = rkB(istage)*Lambda(1:N,iadj)
       IF (istage < rkS) THEN
           DO j = istage+1, rkS
               CALL WAXPY(N,rkA(j,istage),U(1,iadj,j),1,G,1)
           END DO
       END IF
#ifdef FULL_ALGEBRA  
       TMP = MATMUL(TRANSPOSE(Jac),G)    ! DZ <- Jac(Y+Z)*Y_tlm
#else      
       CALL JacTR_SP_Vec ( Jac, G, TMP )    
#endif      
       G(1:N) = H*TMP(1:N)

DirADJ:IF (DirectADJ) THEN 

            CALL SDIRK_Solve ( 'T', H, N, E_adj, IP_adj, IER, G )
            U(1:N,iadj,istage) = G(1:N)
      
       ELSE DirADJ

            !~~~>  Initializations for Newton iteration
            CALL Set2zero(N,U(1,iadj,istage))
            NewtonDone = .FALSE.
            
NewtonLoop:DO NewtonIter = 1, NewtonMaxit

!~~~>   Prepare the loop-dependent part of the right-hand side
#ifdef FULL_ALGEBRA  
            TMP = MATMUL(TRANSPOSE(Jac),U(1:N,iadj,istage))    
#else      
            CALL JacTR_SP_Vec ( Jac, U(1:N,iadj,istage), TMP )    
#endif      
            DU(1:N) = U(1:N,iadj,istage) - (H*rkGamma)*TMP(1:N) - G(1:N)

!~~~>   Solve the linear system
            CALL SDIRK_Solve ( 'T', H, N, E, IP, IER, DU )
            
!~~~>   Check convergence of Newton iterations
            
            CALL SDIRK_ErrorNorm(N, DU, SCAL, NewtonIncrement)
            IF ( NewtonIter == 1 ) THEN
                Theta      = ABS(ThetaMin)
                NewtonRate = 2.0d0 
            ELSE
                Theta = NewtonIncrement/NewtonIncrementOld
                IF (Theta < 0.99d0) THEN
                    NewtonRate = Theta/(ONE-Theta)
                     ! Predict error at the end of Newton process 
                    NewtonPredictedErr = NewtonIncrement &
                               *Theta**(NewtonMaxit-NewtonIter)/(ONE-Theta)
                     ! Non-convergence of Newton: predicted error too large
                    IF (NewtonPredictedErr >= NewtonTol) EXIT NewtonLoop
                ELSE ! Non-convergence of Newton: Theta too large
                    EXIT NewtonLoop
                END IF
            END IF
            NewtonIncrementOld = NewtonIncrement
            ! Update solution
            U(1:N,iadj,istage) = U(1:N,iadj,istage) - DU(1:N)
            
            ! Check error in Newton iterations
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
            ! AbsTol is often inappropriate for adjoints -
            !    we do at least 4 Newton iterations to ensure convergence
            !    of all adjoint components
            IF ((NewtonIter>=4) .AND. NewtonDone) EXIT NewtonLoop
            
            END DO NewtonLoop
            
            !~~~> If Newton iterations fail employ the direct solution
            IF (.NOT.NewtonDone) THEN                 
               PRINT*,'Problems with Newton Adjoint!!!'
#ifdef FULL_ALGEBRA
               E_adj(1:N,1:N) = -Jac(1:N,1:N)
               DO j=1,N
                  E_adj(j,j) = E_adj(j,j) + ONE/(H*rkGamma)
               END DO
               CALL DGETRF( N, N, E_adj, N, IP_adj, IER )
#else
               E_adj(1:LU_NONZERO) = -Jac(1:LU_NONZERO)
               DO i = 1,NVAR
                  j = LU_DIAG(i); E_adj(j) = E_adj(j) + ONE/(H*rkGamma)
               END DO
               CALL KppDecomp ( E_adj, IER)
#endif
               ISTATUS(Ndec) = ISTATUS(Ndec) + 1
               IF (IER /= 0)  THEN
                   PRINT*,'At stage ',istage,' the matrix used in adjoint', &
                            ' computation is singular'
                   CALL SDIRK_ErrorMsg(-8,T,H,Ierr); RETURN
               END IF
                  CALL SDIRK_Solve ( 'T', H, N, E_adj, IP_adj, IER, G )
                  U(1:N,iadj,istage) = G(1:N)
                 
            END IF

!~~~>  End of simplified Newton iterations
      
       END IF DirADJ

   END DO adj
 
   END DO stages

!~~~> Update adjoint solution
         ! Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
         DO istage = 1,rkS 
             DO iadj = 1,NADJ
                 Lambda(1:N,iadj) = Lambda(1:N,iadj) + U(1:N,iadj,istage)
                 !CALL WAXPY(N,ONE,U(1:N,iadj,istage),1,Lambda(1,iadj),1)
             END DO  
         END DO  

      END DO Tloop

      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE SDIRK_DadjInt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_AllocBuffers
!~~~>  Allocate buffer space for checkpointing
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       INTEGER :: i
   
       ALLOCATE( chk_H(Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer H'; STOP
       END IF   
       ALLOCATE( chk_T(Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer T'; STOP
       END IF   
       ALLOCATE( chk_Y(NVAR,Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer Y'; STOP
       END IF   
       ALLOCATE( chk_Z(NVAR,rkS,Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer K'; STOP
       END IF   
       IF (SaveLU) THEN
#ifdef FULL_ALGEBRA
          ALLOCATE( chk_J(NVAR,NVAR,Max_no_steps),  STAT=i )
#else
          ALLOCATE( chk_J(LU_NONZERO,Max_no_steps), STAT=i )
#endif
          IF (i/=0) THEN
             PRINT*,'Failed allocation of buffer J'; STOP
          END IF   
          ALLOCATE( chk_P(NVAR,Max_no_steps), STAT=i )
          IF (i/=0) THEN
             PRINT*,'Failed allocation of buffer P'; STOP
          END IF   
       END IF   
 
     END SUBROUTINE SDIRK_AllocBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE SDIRK_FreeBuffers
!~~~>  Dallocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       INTEGER :: i
   
       DEALLOCATE( chk_H, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer H'; STOP
       END IF   
       DEALLOCATE( chk_T, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer T'; STOP
       END IF   
       DEALLOCATE( chk_Y, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer Y'; STOP
       END IF   
       DEALLOCATE( chk_Z, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer K'; STOP
       END IF   
       IF (SaveLU) THEN
          DEALLOCATE( chk_J, STAT=i )
          IF (i/=0) THEN
             PRINT*,'Failed deallocation of buffer J'; STOP
          END IF   
          DEALLOCATE( chk_P, STAT=i )
          IF (i/=0) THEN
             PRINT*,'Failed deallocation of buffer P'; STOP
          END IF   
       END IF   
 
     END SUBROUTINE SDIRK_FreeBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE SDIRK_Push( T, H, Y, Z, E, P )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       KPP_REAL :: T, H, Y(NVAR), Z(NVAR,Smax) 
       INTEGER       :: P(NVAR)
#ifdef FULL_ALGEBRA
       KPP_REAL :: E(NVAR,NVAR)
#else       
       KPP_REAL :: E(LU_NONZERO)
#endif
   
       stack_ptr = stack_ptr + 1
       IF ( stack_ptr > Max_no_steps ) THEN
         PRINT*,'Push failed: buffer overflow'
         STOP
       END IF  
       chk_H( stack_ptr ) = H
       chk_T( stack_ptr ) = T
       chk_Y(1:NVAR,stack_ptr) = Y(1:NVAR)
       chk_Z(1:NVAR,1:rkS,stack_ptr) = Z(1:NVAR,1:rkS)
       IF (SaveLU) THEN 
#ifdef FULL_ALGEBRA
          chk_J(1:NVAR,1:NVAR,stack_ptr) = E(1:NVAR,1:NVAR)
          chk_P(1:NVAR,stack_ptr)        = P(1:NVAR)
#else       
          chk_J(1:LU_NONZERO,stack_ptr)  = E(1:LU_NONZERO)
#endif
       END IF
  
      END SUBROUTINE SDIRK_Push
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_Pop( T, H, Y, Z, E, P )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       KPP_REAL :: T, H, Y(NVAR), Z(NVAR,Smax)
       INTEGER       :: P(NVAR)
#ifdef FULL_ALGEBRA
       KPP_REAL :: E(NVAR,NVAR)
#else       
       KPP_REAL :: E(LU_NONZERO)
#endif
   
       IF ( stack_ptr <= 0 ) THEN
         PRINT*,'Pop failed: empty buffer'
         STOP
       END IF  
       H = chk_H( stack_ptr )
       T = chk_T( stack_ptr )
       Y(1:NVAR) = chk_Y(1:NVAR,stack_ptr)
       Z(1:NVAR,1:rkS) = chk_Z(1:NVAR,1:rkS,stack_ptr)
       IF (SaveLU) THEN
#ifdef FULL_ALGEBRA
          E(1:NVAR,1:NVAR) = chk_J(1:NVAR,1:NVAR,stack_ptr)
          P(1:NVAR)        = chk_P(1:NVAR,stack_ptr)
#else       
          E(1:LU_NONZERO)  = chk_J(1:LU_NONZERO,stack_ptr)
#endif
       END IF

       stack_ptr = stack_ptr - 1
  
      END SUBROUTINE SDIRK_Pop

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER :: i, ITOL
      KPP_REAL :: AbsTol(NVAR), RelTol(NVAR), &
                       Y(NVAR), SCAL(NVAR)
      IF (ITOL == 0) THEN
        DO i=1,NVAR
          SCAL(i) = ONE / ( AbsTol(1)+RelTol(1)*ABS(Y(i)) )
        END DO
      ELSE
        DO i=1,NVAR
          SCAL(i) = ONE / ( AbsTol(i)+RelTol(i)*ABS(Y(i)) )
        END DO
      END IF
      END SUBROUTINE SDIRK_ErrorScale
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_ErrorNorm(N, Y, SCAL, Err)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      
      INTEGER :: i, N
      KPP_REAL :: Y(N), SCAL(N), Err      
      Err = ZERO
      DO i=1,N
           Err = Err+(Y(i)*SCAL(i))**2
      END DO
      Err = MAX( SQRT(Err/DBLE(N)), 1.0d-10 )
!
      END SUBROUTINE SDIRK_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
SUBROUTINE SDIRK_ErrorMsg(Code,T,H,Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   KPP_REAL, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: Ierr

   Ierr = Code
   PRINT * , &
     'Forced exit from SDIRK due to the following error:'

   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)
      PRINT * , '--> Improper value for maximal no of Newton iterations'
    CASE (-3)
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5)
      PRINT * , '--> Improper tolerance values'
    CASE (-6)
      PRINT * , '--> No of steps exceeds maximum bound', max_no_steps
    CASE (-7)
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE SDIRK_ErrorMsg
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_PrepareMatrix ( H, T, Y, FJAC, &
                   SkipJac, SkipLU, E, IP, Reject, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Compute the matrix E = 1/(H*GAMMA)*Jac, and its decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
      IMPLICIT NONE
      
      KPP_REAL, INTENT(INOUT) :: H
      KPP_REAL, INTENT(IN)    :: T, Y(NVAR)
      LOGICAL, INTENT(INOUT)       :: SkipJac,SkipLU,Reject
      INTEGER, INTENT(OUT)         :: ISING, IP(NVAR)
#ifdef FULL_ALGEBRA
      KPP_REAL, INTENT(INOUT) :: FJAC(NVAR,NVAR)
      KPP_REAL, INTENT(OUT)   :: E(NVAR,NVAR)
#else
      KPP_REAL, INTENT(INOUT) :: FJAC(LU_NONZERO)
      KPP_REAL, INTENT(OUT)   :: E(LU_NONZERO)
#endif
      KPP_REAL                :: HGammaInv
      INTEGER                      :: i, j, ConsecutiveSng

      ConsecutiveSng = 0
      ISING = 1
      
Hloop: DO WHILE (ISING /= 0)
      
      HGammaInv = ONE/(H*rkGamma)

!~~~>  Compute the Jacobian
!      IF (SkipJac) THEN
!          SkipJac = .FALSE.
!      ELSE
      IF (.NOT. SkipJac) THEN
          CALL JAC_CHEM( T, Y, FJAC )
          ISTATUS(Njac) = ISTATUS(Njac) + 1
      END IF  
      
#ifdef FULL_ALGEBRA
      DO j=1,NVAR
         DO i=1,NVAR
            E(i,j) = -FJAC(i,j)
         END DO
         E(j,j) = E(j,j)+HGammaInv
      END DO
      CALL DGETRF( NVAR, NVAR, E, NVAR, IP, ISING )
#else
      DO i = 1,LU_NONZERO
            E(i) = -FJAC(i)
      END DO
      DO i = 1,NVAR
          j = LU_DIAG(i); E(j) = E(j) + HGammaInv
      END DO
      CALL KppDecomp ( E, ISING)
      IP(1) = 1
#endif
      ISTATUS(Ndec) = ISTATUS(Ndec) + 1

      IF (ISING /= 0) THEN
          WRITE (6,*) ' MATRIX IS SINGULAR, ISING=',ISING,' T=',T,' H=',H
          ISTATUS(Nsng) = ISTATUS(Nsng) + 1; ConsecutiveSng = ConsecutiveSng + 1
          IF (ConsecutiveSng >= 6) RETURN ! Failure
          H = 0.5d0*H
          SkipJac = .FALSE.
          SkipLU  = .FALSE.
          Reject  = .TRUE.
      END IF
      
      END DO Hloop

      END SUBROUTINE SDIRK_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_Solve ( Transp, H, N, E, IP, ISING, RHS )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Solves the system (H*Gamma-Jac)*x = R
!      using the LU decomposition of E = I - 1/(H*Gamma)*Jac
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER, INTENT(IN)       :: N, IP(N), ISING
      CHARACTER, INTENT(IN)     :: Transp
      KPP_REAL, INTENT(IN) :: H
#ifdef FULL_ALGEBRA
      KPP_REAL, INTENT(IN) :: E(NVAR,NVAR)
#else
      KPP_REAL, INTENT(IN) :: E(LU_NONZERO)
#endif
      KPP_REAL, INTENT(INOUT) :: RHS(N)
      KPP_REAL                :: HGammaInv
      
      HGammaInv = ONE/(H*rkGamma)
      CALL WSCAL(N,HGammaInv,RHS,1)
      SELECT CASE (TRANSP)
      CASE ('N')
#ifdef FULL_ALGEBRA  
         CALL DGETRS( 'N', N, 1, E, N, IP, RHS, N, ISING )
#else
         CALL KppSolve(E, RHS)
#endif
      CASE ('T')
#ifdef FULL_ALGEBRA  
         CALL DGETRS( 'T', N, 1, E, N, IP, RHS, N, ISING )
#else
         CALL KppSolveTR(E, RHS, RHS)
#endif
      CASE DEFAULT
         PRINT*,'Error in SDIRK_Solve. Unknown Transp argument:',Transp
         STOP
      END SELECT
      ISTATUS(Nsol) = ISTATUS(Nsol) + 1
 
      END SUBROUTINE SDIRK_Solve


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk4a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S4A
! Number of stages
      rkS = 5

! Method coefficients
      rkGamma = .2666666666666666666666666666666667d0

      rkA(1,1) = .2666666666666666666666666666666667d0
      rkA(2,1) = .5000000000000000000000000000000000d0
      rkA(2,2) = .2666666666666666666666666666666667d0
      rkA(3,1) = .3541539528432732316227461858529820d0
      rkA(3,2) = -.5415395284327323162274618585298197d-1
      rkA(3,3) = .2666666666666666666666666666666667d0
      rkA(4,1) = .8515494131138652076337791881433756d-1
      rkA(4,2) = -.6484332287891555171683963466229754d-1
      rkA(4,3) = .7915325296404206392428857585141242d-1
      rkA(4,4) = .2666666666666666666666666666666667d0
      rkA(5,1) = 2.100115700566932777970612055999074d0
      rkA(5,2) = -.7677800284445976813343102185062276d0
      rkA(5,3) = 2.399816361080026398094746205273880d0
      rkA(5,4) = -2.998818699869028161397714709433394d0
      rkA(5,5) = .2666666666666666666666666666666667d0

      rkB(1)   = 2.100115700566932777970612055999074d0
      rkB(2)   = -.7677800284445976813343102185062276d0
      rkB(3)   = 2.399816361080026398094746205273880d0
      rkB(4)   = -2.998818699869028161397714709433394d0
      rkB(5)   = .2666666666666666666666666666666667d0

      rkBhat(1)= 2.885264204387193942183851612883390d0
      rkBhat(2)= -.1458793482962771337341223443218041d0
      rkBhat(3)= 2.390008682465139866479830743628554d0
      rkBhat(4)= -4.129393538556056674929560012190140d0
      rkBhat(5)= 0.d0

      rkC(1)   = .2666666666666666666666666666666667d0
      rkC(2)   = .7666666666666666666666666666666667d0
      rkC(3)   = .5666666666666666666666666666666667d0
      rkC(4)   = .3661315380631796996374935266701191d0
      rkC(5)   = 1.d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.d0
      rkD(2)   = 0.d0
      rkD(3)   = 0.d0
      rkD(4)   = 0.d0
      rkD(5)   = 1.d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   = -.6804000050475287124787034884002302d0
      rkE(2)   = 1.558961944525217193393931795738823d0
      rkE(3)   = -13.55893003128907927748632408763868d0
      rkE(4)   = 15.48522576958521253098585004571302d0
      rkE(5)   = 1.d0

! Local order of Err estimate
      rkElo    = 4

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = 1.875000000000000000000000000000000d0
      rkTheta(3,1) = 1.708847304091539528432732316227462d0
      rkTheta(3,2) = -.2030773231622746185852981969486824d0
      rkTheta(4,1) = .2680325578937783958847157206823118d0
      rkTheta(4,2) = -.1828840955527181631794050728644549d0
      rkTheta(4,3) = .2968246986151577397160821594427966d0
      rkTheta(5,1) = .9096171815241460655379433581446771d0
      rkTheta(5,2) = -3.108254967778352416114774430509465d0
      rkTheta(5,3) = 12.33727431701306195581826123274001d0
      rkTheta(5,4) = -11.24557012450885560524143016037523d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = 2.875000000000000000000000000000000d0
      rkAlpha(3,1) = .8500000000000000000000000000000000d0
      rkAlpha(3,2) = .4434782608695652173913043478260870d0
      rkAlpha(4,1) = .7352046091658870564637910527807370d0
      rkAlpha(4,2) = -.9525565003057343527941920657462074d-1
      rkAlpha(4,3) = .4290111305453813852259481840631738d0
      rkAlpha(5,1) = -16.10898993405067684831655675112808d0
      rkAlpha(5,2) = 6.559571569643355712998131800797873d0
      rkAlpha(5,3) = -15.90772144271326504260996815012482d0
      rkAlpha(5,4) = 25.34908987169226073668861694892683d0
               
!~~~> Coefficients for continuous solution
!          rkD(1,1)= 24.74416644927758d0
!          rkD(1,2)= -4.325375951824688d0
!          rkD(1,3)= 41.39683763286316d0
!          rkD(1,4)= -61.04144619901784d0
!          rkD(1,5)= -3.391332232917013d0
!          rkD(2,1)= -51.98245719616925d0
!          rkD(2,2)= 10.52501981094525d0
!          rkD(2,3)= -154.2067922191855d0
!          rkD(2,4)= 214.3082125319825d0
!          rkD(2,5)= 14.71166018088679d0
!          rkD(3,1)= 33.14347947522142d0
!          rkD(3,2)= -19.72986789558523d0
!          rkD(3,3)= 230.4878502285804d0
!          rkD(3,4)= -287.6629744338197d0
!          rkD(3,5)= -18.99932366302254d0
!          rkD(4,1)= -5.905188728329743d0
!          rkD(4,2)= 13.53022403646467d0
!          rkD(4,3)= -117.6778956422581d0
!          rkD(4,4)= 134.3962081008550d0
!          rkD(4,5)= 8.678995715052762d0
!
!         DO i=1,4  ! CONTi <-- Sum_j rkD(i,j)*Zj
!           CALL Set2zero(N,CONT(1,i))
!           DO j = 1,rkS
!             CALL WAXPY(N,rkD(i,j),Z(1,j),1,CONT(1,i),1)
!           END DO
!         END DO
          
          rkELO = 4.0d0
          
      END SUBROUTINE Sdirk4a

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk4b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S4B
! Number of stages
      rkS = 5

! Method coefficients
      rkGamma = .25d0

      rkA(1,1) = 0.25d0
      rkA(2,1) = 0.5d00
      rkA(2,2) = 0.25d0
      rkA(3,1) = 0.34d0
      rkA(3,2) =-0.40d-1
      rkA(3,3) = 0.25d0
      rkA(4,1) = 0.2727941176470588235294117647058824d0
      rkA(4,2) =-0.5036764705882352941176470588235294d-1
      rkA(4,3) = 0.2757352941176470588235294117647059d-1
      rkA(4,4) = 0.25d0
      rkA(5,1) = 1.041666666666666666666666666666667d0
      rkA(5,2) =-1.020833333333333333333333333333333d0
      rkA(5,3) = 7.812500000000000000000000000000000d0
      rkA(5,4) =-7.083333333333333333333333333333333d0
      rkA(5,5) = 0.25d0

      rkB(1)   =  1.041666666666666666666666666666667d0
      rkB(2)   = -1.020833333333333333333333333333333d0
      rkB(3)   =  7.812500000000000000000000000000000d0
      rkB(4)   = -7.083333333333333333333333333333333d0
      rkB(5)   =  0.250000000000000000000000000000000d0

      rkBhat(1)=  1.069791666666666666666666666666667d0
      rkBhat(2)= -0.894270833333333333333333333333333d0
      rkBhat(3)=  7.695312500000000000000000000000000d0
      rkBhat(4)= -7.083333333333333333333333333333333d0
      rkBhat(5)=  0.212500000000000000000000000000000d0

      rkC(1)   = 0.25d0
      rkC(2)   = 0.75d0
      rkC(3)   = 0.55d0
      rkC(4)   = 0.50d0
      rkC(5)   = 1.00d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.0d0
      rkD(2)   = 0.0d0
      rkD(3)   = 0.0d0
      rkD(4)   = 0.0d0
      rkD(5)   = 1.0d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   =  0.5750d0
      rkE(2)   =  0.2125d0
      rkE(3)   = -4.6875d0
      rkE(4)   =  4.2500d0
      rkE(5)   =  0.1500d0

! Local order of Err estimate
      rkElo    = 4

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = 2.d0
      rkTheta(3,1) = 1.680000000000000000000000000000000d0
      rkTheta(3,2) = -.1600000000000000000000000000000000d0
      rkTheta(4,1) = 1.308823529411764705882352941176471d0
      rkTheta(4,2) = -.1838235294117647058823529411764706d0
      rkTheta(4,3) = 0.1102941176470588235294117647058824d0
      rkTheta(5,1) = -3.083333333333333333333333333333333d0
      rkTheta(5,2) = -4.291666666666666666666666666666667d0
      rkTheta(5,3) =  34.37500000000000000000000000000000d0
      rkTheta(5,4) = -28.33333333333333333333333333333333d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = 3.
      rkAlpha(3,1) = .8800000000000000000000000000000000d0
      rkAlpha(3,2) = .4400000000000000000000000000000000d0
      rkAlpha(4,1) = .1666666666666666666666666666666667d0
      rkAlpha(4,2) = -.8333333333333333333333333333333333d-1
      rkAlpha(4,3) = .9469696969696969696969696969696970d0
      rkAlpha(5,1) = -6.d0
      rkAlpha(5,2) = 9.d0
      rkAlpha(5,3) = -56.81818181818181818181818181818182d0
      rkAlpha(5,4) = 54.d0
      
      END SUBROUTINE Sdirk4b

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk2a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S2A
! Number of stages
      rkS = 2

! Method coefficients
      rkGamma = .2928932188134524755991556378951510d0

      rkA(1,1) = .2928932188134524755991556378951510d0
      rkA(2,1) = .7071067811865475244008443621048490d0
      rkA(2,2) = .2928932188134524755991556378951510d0

      rkB(1)   = .7071067811865475244008443621048490d0
      rkB(2)   = .2928932188134524755991556378951510d0

      rkBhat(1)= .6666666666666666666666666666666667d0
      rkBhat(2)= .3333333333333333333333333333333333d0

      rkC(1)   = 0.292893218813452475599155637895151d0
      rkC(2)   = 1.0d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.0d0
      rkD(2)   = 1.0d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   =  0.4714045207910316829338962414032326d0
      rkE(2)   = -0.1380711874576983496005629080698993d0

! Local order of Err estimate
      rkElo    = 2

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = 2.414213562373095048801688724209698d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = 3.414213562373095048801688724209698d0
          
      END SUBROUTINE Sdirk2a

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk2b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S2B
! Number of stages
      rkS      = 2

! Method coefficients
      rkGamma  = 1.707106781186547524400844362104849d0

      rkA(1,1) = 1.707106781186547524400844362104849d0
      rkA(2,1) = -.707106781186547524400844362104849d0
      rkA(2,2) = 1.707106781186547524400844362104849d0

      rkB(1)   = -.707106781186547524400844362104849d0
      rkB(2)   = 1.707106781186547524400844362104849d0

      rkBhat(1)= .6666666666666666666666666666666667d0
      rkBhat(2)= .3333333333333333333333333333333333d0

      rkC(1)   = 1.707106781186547524400844362104849d0
      rkC(2)   = 1.0d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.0d0
      rkD(2)   = 1.0d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   = -.4714045207910316829338962414032326d0
      rkE(2)   =  .8047378541243650162672295747365659d0

! Local order of Err estimate
      rkElo    = 2

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = -.414213562373095048801688724209698d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = .5857864376269049511983112757903019d0
      
      END SUBROUTINE Sdirk2b


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk3a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S3A
! Number of stages
      rkS = 3

! Method coefficients
      rkGamma = .2113248654051871177454256097490213d0

      rkA(1,1) = .2113248654051871177454256097490213d0
      rkA(2,1) = .2113248654051871177454256097490213d0
      rkA(2,2) = .2113248654051871177454256097490213d0
      rkA(3,1) = .2113248654051871177454256097490213d0
      rkA(3,2) = .5773502691896257645091487805019573d0
      rkA(3,3) = .2113248654051871177454256097490213d0

      rkB(1)   = .2113248654051871177454256097490213d0
      rkB(2)   = .5773502691896257645091487805019573d0
      rkB(3)   = .2113248654051871177454256097490213d0

      rkBhat(1)= .2113248654051871177454256097490213d0
      rkBhat(2)= .6477918909913548037576239837516312d0
      rkBhat(3)= .1408832436034580784969504064993475d0

      rkC(1)   = .2113248654051871177454256097490213d0
      rkC(2)   = .4226497308103742354908512194980427d0
      rkC(3)   = 1.d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.d0
      rkD(2)   = 0.d0
      rkD(3)   = 1.d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   =  0.9106836025229590978424821138352906d0
      rkE(2)   = -1.244016935856292431175815447168624d0
      rkE(3)   =  0.3333333333333333333333333333333333d0

! Local order of Err estimate
      rkElo    = 2

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) =  1.0d0
      rkTheta(3,1) = -1.732050807568877293527446341505872d0
      rkTheta(3,2) =  2.732050807568877293527446341505872d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) =   2.0d0
      rkAlpha(3,1) = -12.92820323027550917410978536602349d0
      rkAlpha(3,2) =   8.83012701892219323381861585376468d0

      END SUBROUTINE Sdirk3a

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   END SUBROUTINE SDIRKADJ ! AND ALL ITS INTERNAL PROCEDURES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE FUN_CHEM( T, Y, P )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      USE KPP_ROOT_Parameters, ONLY: NVAR
      USE KPP_ROOT_Global, ONLY: TIME, FIX, RCONST
      USE KPP_ROOT_Function, ONLY: Fun
      USE KPP_ROOT_Rates, ONLY: Update_SUN, Update_RCONST, Update_PHOTO

      KPP_REAL :: T, Told
      KPP_REAL :: Y(NVAR), P(NVAR)
      
      Told = TIME
      TIME = T
      CALL Update_SUN()
      CALL Update_RCONST()
      
      CALL Fun( Y, FIX, RCONST, P )
      
      TIME = Told
      
      END SUBROUTINE FUN_CHEM


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE JAC_CHEM( T, Y, JV )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      USE KPP_ROOT_Parameters, ONLY: NVAR, LU_NONZERO
      USE KPP_ROOT_Global, ONLY: TIME, FIX, RCONST
      USE KPP_ROOT_Jacobian, ONLY: Jac_SP,LU_IROW,LU_ICOL
      USE KPP_ROOT_Rates, ONLY: Update_SUN, Update_RCONST, Update_PHOTO
  
      KPP_REAL ::  T, Told
      KPP_REAL ::  Y(NVAR)
#ifdef FULL_ALGEBRA
      KPP_REAL :: JS(LU_NONZERO), JV(NVAR,NVAR)
      INTEGER :: i, j
#else
      KPP_REAL :: JV(LU_NONZERO)
#endif
 
      Told = TIME
      TIME = T
      CALL Update_SUN()
      CALL Update_RCONST()

#ifdef FULL_ALGEBRA
      CALL Jac_SP(Y, FIX, RCONST, JS)
      DO j=1,NVAR
         DO j=1,NVAR
          JV(i,j) = 0.0d0
         END DO
      END DO
      DO i=1,LU_NONZERO
         JV(LU_IROW(i),LU_ICOL(i)) = JS(i)
      END DO
#else
      CALL Jac_SP(Y, FIX, RCONST, JV)
#endif
      TIME = Told

      END SUBROUTINE JAC_CHEM

END MODULE KPP_ROOT_Integrator


