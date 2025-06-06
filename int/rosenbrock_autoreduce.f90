!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Rosenbrock - Implementation of several Rosenbrock methods:             !
!               * Ros2                                                    !
!               * Ros3                                                    !
!               * Ros4                                                    !
!               * Rodas3                                                  !
!               * Rodas4                                                  !
!  By default the code employs the KPP sparse linear algebra routines     !
!  Compile with -DFULL_ALGEBRA to use full linear algebra (LAPACK)        !
!                                                                         !
!    (C)  Adrian Sandu, August 2004                                       !
!    Virginia Polytechnic Institute and State University                  !
!    Contact: sandu@cs.vt.edu                                             !
!    Revised by Philipp Miehe and Adrian Sandu, May 2006                  !
!                                                                         !
!    Revised by Mike Long and Haipeng Lin to add auto-reduce fun.         !
!    Harvard University, Atmospheric Chemistry Modeling Group             !
!    Contact: hplin@seas.harvard.edu                April 2022            !
!                                                                         !
!    This implementation is part of KPP - the Kinetic PreProcessor        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE KPP_ROOT_Integrator

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_Global
  IMPLICIT NONE
  PUBLIC
  SAVE

!~~~> Flags to determine if we should call the UPDATE_* routines from within
!~~~> the integrator.  If using KPP in an external model, you might want to
!~~~> disable these calls (via ICNTRL(15)) to avoid excess computations.
  LOGICAL, PRIVATE :: Do_Update_RCONST
  LOGICAL, PRIVATE :: Do_Update_PHOTO
  LOGICAL, PRIVATE :: Do_Update_SUN

!~~~>  Statistics on the work performed by the Rosenbrock method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                        Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                        Ntexit=1, Nhexit=2, Nhnew = 3,  &
                        NARthr=4

CONTAINS

SUBROUTINE INTEGRATE( TIN,       TOUT,      ICNTRL_U, RCNTRL_U,  &
                      ISTATUS_U, RSTATUS_U, IERR_U              )

   USE KPP_ROOT_Util, ONLY : Integrator_Update_Options

   IMPLICIT NONE

   KPP_REAL, INTENT(IN) :: TIN  ! Start Time
   KPP_REAL, INTENT(IN) :: TOUT ! End Time
   !~~~> Optional input parameters and statistics
   INTEGER,  INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   KPP_REAL, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   KPP_REAL, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: IERR_U

   KPP_REAL :: RCNTRL(20), RSTATUS(20)
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR

   INTEGER, SAVE :: Ntotal = 0

   !~~~> Zero input and output arrays for safety's sake
   ICNTRL     = 0
   RCNTRL     = 0.0_dp
   ISTATUS    = 0
   RSTATUS    = 0.0_dp

   !~~~> fine-tune the integrator:
   ICNTRL(15) = 5       ! Call Update_SUN and Update_RCONST from w/in the int.

   !~~~> if optional parameters are given, and if they are /= 0,
   !     then use them to overwrite default settings
   IF ( PRESENT( ICNTRL_U ) ) THEN
      WHERE( ICNTRL_U /= 0 ) ICNTRL = ICNTRL_U
   ENDIF
   IF ( PRESENT( RCNTRL_U ) ) THEN
      WHERE( RCNTRL_U > 0 ) RCNTRL = RCNTRL_U
   ENDIF

   !~~~> Determine the settings of the Do_Update_* flags, which determine
   !~~~> whether or not we need to call Update_* routines in the integrator
   !~~~> (or not, if we are calling them from a higher-level)
   ! ICNTRL(15) = -1 ! Do not call Update_* functions within the integrator
   !            =  0 ! Status quo
   !            =  1 ! Call Update_RCONST from within the integrator
   !            =  2 ! Call Update_PHOTO from within the integrator
   !            =  3 ! Call Update_RCONST and Update_PHOTO from w/in the int.
   !            =  4 ! Call Update_SUN from within the integrator
   !            =  5 ! Call Update_SUN and Update_RCONST from within the int.
   !            =  6 ! Call Update_SUN and Update_PHOTO from within the int.
   !            =  7 ! Call Update_SUN, Update_PHOTO, Update_RCONST w/in int.
   CALL Integrator_Update_Options( ICNTRL(15),          &
                                   Do_Update_RCONST,    &
                                   Do_Update_PHOTO,     &
                                   Do_Update_Sun       )

   !~~~> In order to remove the prior EQUIVALENCE statements (which
   !~~~> are not thread-safe), we now have declared VAR and FIX as
   !~~~> threadprivate pointer variables that can point to C.
   VAR => C(1:NVAR )
   FIX => C(NVAR+1:NSPEC)

   !~~~> Call the integrator
   CALL Rosenbrock( NVAR,   VAR,    TIN,     TOUT,    ATOL, RTOL,  &
                    RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR        )

   !~~~> Free pointers
   VAR => NULL()
   FIX => NULL()

   !~~~> Debug option: show number of steps
   !Ntotal = Ntotal + ISTATUS(Nstp)
   !PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

   STEPMIN = RSTATUS(Nhexit)

   !~~~> if optional parameters are given for output
   !~~~> use them to store information in them
   IF ( PRESENT( ISTATUS_U ) ) ISTATUS_U = ISTATUS
   IF ( PRESENT( RSTATUS_U ) ) RSTATUS_U = RSTATUS
   IF ( PRESENT( IERR_U    ) ) IERR_U    = IERR

END SUBROUTINE INTEGRATE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Rosenbrock(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using a Rosenbrock method defined by:
!
!     G = 1/(H*gamma(1)) - Jac(t0,Y0)
!     T_i = t0 + Alpha(i)*H
!     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
!     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
!         gamma(i)*dF/dT(t0, Y0)
!     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j
!
!    For details on Rosenbrock methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.
!    The codes contained in the book inspired this implementation.
!
!    (C)  Adrian Sandu, August 2004
!    Virginia Polytechnic Institute and State University
!    Contact: sandu@cs.vt.edu
!    Revised by Philipp Miehe and Adrian Sandu, May 2006
!    This implementation is part of KPP - the Kinetic PreProcessor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!-     Y(N)    = vector of initial conditions (at T=Tstart)
!-    [Tstart,Tend]  = time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE  Fun( T, Y, Ydot ) = ODE function,
!                       returns Ydot = Y' = F(T,Y)
!- SUBROUTINE  Jac( T, Y, Jcb ) = Jacobian of the ODE function,
!                       returns Jcb = dFun/dY
!-    ICNTRL(1:20)    = integer inputs parameters
!-    RCNTRL(1:20)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!-    Y(N)    -> vector of final states (at T->Tend)
!-    ISTATUS(1:20)   -> integer output parameters
!-    RSTATUS(1:20)   -> real output parameters
!-    IERR            -> job status upon return
!                        success (positive value) or
!                        failure (negative value)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1) = 1: F = F(y)   Independent of T (AUTONOMOUS)
!              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)
!
!    ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3)  -> selection of a particular Rosenbrock method
!        = 0 :    Rodas3 (default)
!        = 1 :    Ros2
!        = 2 :    Ros3
!        = 3 :    Ros4
!        = 4 :    Rodas3
!        = 5 :    Rodas4
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0) the default value of 200000 is used
!
!    ICNTRL(12)  -> use auto-reduce solver? set threshold in RCNTRL(12)
!    ICNTRL(13)  -> ... append slow species when auto-reducing?
!    ICNTRL(14) -> choose a target species instead for determining threshold?
!                  if yes, specify idx. then RCNTRL(12) is obsolete.
!
!    ICNTRL(15) -> Toggles calling of Update_* functions w/in the integrator
!        = -1 :  Do not call Update_* functions within the integrator
!        =  0 :  Status quo
!        =  1 :  Call Update_RCONST from within the integrator
!        =  2 :  Call Update_PHOTO from within the integrator
!        =  3 :  Call Update_RCONST and Update_PHOTO from w/in the int.
!        =  4 :  Call Update_SUN from within the integrator
!        =  5 :  Call Update_SUN and Update_RCONST from within the int.
!        =  6 :  Call Update_SUN and Update_PHOTO from within the int.
!        =  7 :  Call Update_SUN, Update_PHOTO, Update_RCONST w/in the int.
!
!    ICNTRL(16) -> 
!        = 0 : allow negative concentrations (default)
!        = 1 : set negative concentrations to zero
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                          (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!         than the predicted value  (default=0.9)
!
!    RCNTRL(12) -> threshold for auto-reduction (req. ICNTRL(12)) (default=100)
!    RCNTRL(14) -> AR threshold ratio (default=0.01)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!    OUTPUT ARGUMENTS:
!    -----------------
!
!    T           -> T value for which the solution has been computed
!                     (after successful return T=Tend).
!
!    Y(N)        -> Numerical solution at T
!
!    IDID        -> Reports on successfulness upon return:
!                    = 1 for success
!                    < 0 for error (value equals error code)
!
!    ISTATUS(1)  -> No. of function calls
!    ISTATUS(2)  -> No. of jacobian calls
!    ISTATUS(3)  -> No. of steps
!    ISTATUS(4)  -> No. of accepted steps
!    ISTATUS(5)  -> No. of rejected steps (except at very beginning)
!    ISTATUS(6)  -> No. of LU decompositions
!    ISTATUS(7)  -> No. of forward/backward substitutions
!    ISTATUS(8)  -> No. of singular matrix decompositions
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                     computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!                   For multiple restarts, use Hnew as Hstart
!                     in the subsequent run
!    RSTATUS(4)  -> ARthr, last auto-reduction threshold determined
!                   only if AR is on (ICNTRL(12)) and key spc (ICNTRL(14))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_LinearAlgebra
  IMPLICIT NONE

!~~~>  Arguments
   INTEGER,       INTENT(IN)    :: N
   KPP_REAL, INTENT(INOUT) :: Y(N)
   KPP_REAL, INTENT(IN)    :: Tstart,Tend
   KPP_REAL, INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   KPP_REAL, INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   KPP_REAL, INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   KPP_REAL :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   KPP_REAL :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   KPP_REAL :: Hmin, Hmax, Hstart
   KPP_REAL :: Texit, Redux_Threshold
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol, Autoreduce, Autoreduce_Append
   INTEGER       :: AR_target_spc
   KPP_REAL :: AR_thr_ratio
!~~~>   Parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   KPP_REAL, PARAMETER :: DeltaMin = 1.0E-5_dp

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:4) = ZERO

!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3)
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~> Auto-reduction toggle
   Autoreduce    = .false.
   IF (ICNTRL(12) == 1) Autoreduce = .true.

   Autoreduce_Append = ICNTRL(13) == 1
!~~~> Target species (if zero, uses the regular threshold)
   AR_target_spc = ICNTRL(14)

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
!~~~> Auto-reduction threshold
    Redux_threshold = 1.d2
    IF (RCNTRL(12) > ZERO) THEN
       Redux_Threshold = RCNTRL(12)
    ELSEIF (RCNTRL(12) < ZERO) THEN
       Autoreduce = .false.
    ENDIF
!~~~> Auto-reduction threshold ratio (only if ICNTRL(14) is not zero)
    AR_thr_ratio = RCNTRL(14)
!~~~>  CALL Auto-reducing Rosenbrock method
    IF ( Autoreduce .and. .not. Autoreduce_Append ) THEN
         ! ros_yIntegrator is the aggressively micro-optimized revision by Haipeng Lin.
         ! ros_cIntegrator is the original auto-reduce implementation by Mike Long.
         CALL ros_yIntegrator(Y, Tstart, Tend, Texit,   &
         AbsTol, RelTol,                          &
         !  Integration parameters
         Autonomous, VectorTol, Max_no_steps,     &
         Roundoff, Hmin, Hmax, Hstart,            &
         FacMin, FacMax, FacRej, FacSafe,         &
         ! Autoreduce threshold
         redux_threshold, AR_target_spc,          &
         AR_thr_ratio,                            &
         !  Error indicator
         IERR)
    ENDIF

!~~~> CALL Auto-reducing Rosenbrock method capable of append
! (this version is less aggressively optimized.)
    IF ( Autoreduce .and. Autoreduce_Append ) THEN
         ! ros_yIntegratorA is the append version of the AR integrator. It has less
         ! optimizations because of the need to update Prod/Loss.
         CALL ros_yIntegratorA(Y, Tstart, Tend, Texit,   &
         AbsTol, RelTol,                          &
         !  Integration parameters
         Autonomous, VectorTol, Max_no_steps,     &
         Roundoff, Hmin, Hmax, Hstart,            &
         FacMin, FacMax, FacRej, FacSafe,         &
         ! Autoreduce threshold
         redux_threshold, AR_target_spc,          &
         AR_thr_ratio,                            &
         !  Error indicator
         IERR)
    ENDIF
    
!~~~>  CALL Normal Rosenbrock method
    IF ( .not. Autoreduce .or. IERR .eq. -99 ) &
         CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
         AbsTol, RelTol,                          &
         !  Integration parameters
         Autonomous, VectorTol, Max_no_steps,     &
         Roundoff, Hmin, Hmax, Hstart,            &
         FacMin, FacMax, FacRej, FacSafe,         &
         !  Error indicator
         IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   KPP_REAL, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:'

   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5)
      PRINT * , '--> Improper tolerance values'
    CASE (-6)
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7)
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   KPP_REAL, INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   KPP_REAL, INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   KPP_REAL, INTENT(OUT) ::  T
!~~~> Input: tolerances
   KPP_REAL, INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   KPP_REAL, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   KPP_REAL, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   KPP_REAL :: Ynew(N), Fcn0(N), Fcn(N)
   KPP_REAL :: K(N*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA
   KPP_REAL :: Jac0(N,N), Ghimj(N,N)
#else
   KPP_REAL :: Jac0(LU_NONZERO), Ghimj(LU_NONZERO)
#endif
   KPP_REAL :: H, Hnew, HC, HG, Fac, Tau
   KPP_REAL :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   KPP_REAL, PARAMETER :: DeltaMin = 1.0E-5_dp
!~~~>  Locally called functions
!    KPP_REAL WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~>  Initial preparations
   DO_SLV  = .true.
   DO_FUN  = .true.
   DO_JVS  = .true.
   
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FunTemplate( T, Y, Fcn0 )
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate( T, Y, Jac0 )
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !slim: CALL WCOPY(N,Fcn0,1,Fcn,1)
         Fcn(1:N) = Fcn0(1:N)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         !slim: CALL WCOPY(N,Y,1,Ynew,1)
         Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
           CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate( Tau, Ynew, Fcn )
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       !slim: CALL WCOPY(N,Fcn,1,K(ioffset+1),1)
       K(ioffset+1:ioffset+N) = Fcn(1:N)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))

   END DO Stage


!~~~>  Compute the new solution
   !slim: CALL WCOPY(N,Y,1,Ynew,1)
   Ynew(1:N) = Y(1:N)
   DO j=1,ros_S
         CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   !slim: CALL WSCAL(N,ZERO,Yerr,1)
   Yerr(1:N) = ZERO
   DO j=1,ros_S
        CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      IF (ICNTRL(16) == 1) THEN
        ! new value is non-negative:
        Y = MAX(Ynew,ZERO)
      ELSE
        !slim: CALL WCOPY(N,Ynew,1,Y,1)
        Y(1:N) = Ynew(1:N)
      ENDIF      
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
   END IF ! Err <= 1

   END DO UntilAccepted

   END DO TimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_Integrator

 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_yIntegrator (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Autoreduce threshold
        threshold, AR_target_spc, AR_thr_ratio,  &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Alternative micro-optimized implementation, hplin, 4/10/22
! which does not resize arrays to rNVAR, instead always keeping to full N size
! and skipping using DO_SLV
!
! All compiled assembly code was verified just short of linking to a proper BLAS.

   USE KPP_ROOT_Global,  ONLY : cNONZERO, rNVAR
   use KPP_ROOT_Monitor, ONLY : SPC_NAMES
   USE KPP_ROOT_JacobianSP

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   KPP_REAL, INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   KPP_REAL, INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   KPP_REAL, INTENT(OUT) ::  T
!~~~> Input: tolerances
   KPP_REAL, INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   KPP_REAL, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   KPP_REAL, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Autoreduction threshold
   KPP_REAL, INTENT(IN) :: threshold
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   KPP_REAL :: Ynew(N), Fcn0(N), Fcn(N), Prod(N), Loss(N), LossY(N)
   KPP_REAL :: K(NVAR*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA
   KPP_REAL :: Jac0(N,N), Ghimj(N,N)
#else
   KPP_REAL :: Jac0(LU_NONZERO), Ghimj(LU_NONZERO)
   KPP_REAL :: cGhimj(LU_NONZERO) ! not known at this point what cNONZERO will be
#endif
   KPP_REAL :: H, Hnew, HC, HG, Fac, Tau
   KPP_REAL :: Err, Yerr(N), Yerrsub(NVAR)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular, Reduced
!~~~>  Local parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   KPP_REAL, PARAMETER :: DeltaMin = 1.0E-5_dp
   INTEGER :: SPC
   KPP_REAL :: alpha_factor ! hplin 4/10/22
!      Inline local parameters for AR.
   INTEGER :: II, III, idx, nrmv, s
   KPP_REAL, INTENT(IN) :: AR_thr_ratio
   KPP_REAL :: AR_thr
   INTEGER, INTENT(IN) :: AR_target_spc
!~~~>  Initial preparations
   DO_SLV  = .true.
   DO_FUN  = .true.
   DO_JVS  = .true.
   Reduced = .false.

   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

   ! reset K - hplin. this is a multiplier that gets applied to ros_A and ros_M
   K = 0.0_dp
   Ghimj = 0.0_dp

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   ENDIF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   ENDIF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

   ! ... 0.40% ... 5.2 %

!~~~>   Compute the function at current time
! this is necessary anyway for the rest of the time loop. do not optimize out.
   IF (T .eq. Tstart) THEN
      CALL FunSplitF(T,Y,Fcn0,Prod,Loss,LossY) ! always calculates PL.
   ELSE  ! ELSE or ne seems reasonably close in performance.
      CALL FunSplitN(T,Y,Fcn0)
   ENDIF
   ! The above Prod, Loss were updated at every TimeLoop, so they no longer
   ! reflect initial condition. The copy operation is extremely expensive,
   ! so we instead run it once at the beginning of the timeloop and use a IF
   ! to switch it to a wholy different codepath...
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

   ! ... 1.02% ... 13.00%

!~~~>  Parse species for reduced computation
   if (.not. reduced) then
      ! Inline the entire reduction operation here.
      iSPC_MAP = 0
      NRMV     = 0
      S        = 1
      AR_thr   = threshold

      ! Target species?
      if(AR_target_spc .gt. 0) then
          AR_thr = AR_thr_ratio * max(LossY(AR_target_spc), Prod(AR_target_spc))           ! Lin et al., 2022 in prep.
          RSTATUS(NARthr) = AR_thr
      endif

      ! Checks should be kept out of tight inner loops.
      IF(keepActive) THEN
       DO i=1,NVAR
         ! Short-circuiting using SKIP is very important here.
         if (.not. keepSpcActive(i) .and. &
             abs(LossY(i)).lt.AR_thr .and. abs(Prod(i)).lt.AR_thr) then ! per Shen et al., 2020
            NRMV=NRMV+1
            ! RMV(NRMV) = i ! not needed unless in append version.
            DO_SLV(i) = .false.
            ! DO_FUN(i) = .false.
            cycle
         endif
         SPC_MAP(S)  = i ! Add to full spc map.
         iSPC_MAP(i) = S
         S=S+1
       ENDDO
      ENDIF

      IF (.not. keepActive) THEN
        DO i=1,NVAR
         ! Short-circuiting using SKIP is very important here.
         if (abs(LossY(i)).lt.AR_thr .and. abs(Prod(i)).lt.AR_thr) then ! per Shen et al., 2020
            NRMV=NRMV+1
            ! RMV(NRMV) = i ! not needed unless in append version.
            DO_SLV(i) = .false.
            ! DO_FUN(i) = .false.
            cycle
         endif
         SPC_MAP(S)  = i ! Add to full spc map.
         iSPC_MAP(i) = S
         S=S+1
       ENDDO
      ENDIF

      rNVAR    = NVAR-NRMV ! Number of active species in the reduced mechanism
      II  = 1
      III = 1
      idx = 0
      ! This loop can be unrolled into two branches.
      ScanFirstNonZero: DO i = 1,LU_NONZERO
         IF ((DO_SLV(LU_IROW(i))).and.(DO_SLV(LU_ICOL(i)))) THEN
            idx = 1
            cLU_IROW(1) = iSPC_MAP(LU_IROW(i))
            cLU_ICOL(1) = iSPC_MAP(LU_ICOL(i))
            JVS_MAP(1)  = i
            EXIT ScanFirstNonZero
         ENDIF
         DO_JVS(i) = .false.
      ENDDO ScanFirstNonZero
      DO i = i+1,LU_NONZERO ! There is no escape for looping through LU_NONZERO here.
         IF ((DO_SLV(LU_IROW(i))).and.(DO_SLV(LU_ICOL(i)))) THEN
            idx=idx+1 ! counter for the number of non-zero elements in the reduced Jacobian
            cLU_IROW(idx) = iSPC_MAP(LU_IROW(i))
            cLU_ICOL(idx) = iSPC_MAP(LU_ICOL(i))
            JVS_MAP(idx)  = i

            IF (cLU_IROW(idx).ne.cLU_IROW(idx-1)) THEN
               II=II+1
               cLU_CROW(II) = idx
            ENDIF
            IF (cLU_IROW(idx).eq.cLU_ICOL(idx)) THEN
               III=III+1
               cLU_DIAG(III) = idx
            ENDIF
            CYCLE
         ENDIF
         DO_JVS(i) = .false.
      ENDDO

     cNONZERO = idx
     cLU_CROW(1)       = 1 ! 1st index = 1
     cLU_DIAG(1)       = 1 ! 1st index = 1
     cLU_CROW(rNVAR+1) = cNONZERO+1
     cLU_DIAG(rNVAR+1) = cLU_DIAG(rNVAR)+1

     ! this loop approximately 0.6% ... 7.8%
     reduced = .true.
   endif
   !return

   ! ... 1.86% ... 26% // 1.77% ...

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

   ! 1.86 ~ 1.90%

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate(T,Y,Jac0) ! Reacts to DO_JVS()
   ISTATUS(Njac) = ISTATUS(Njac) + 1
!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_cPrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,cGhimj,Pivot,Singular) ! calculate Ghimj as cGhimj(cNONZERO).
   ! cGhimj all elem are rewritten so no need to zero out - warning - hplin 4/10/22
   !
   ! this step has one skip addressing step, -Jac0(JVS_MAP(1:cNONZERO)). slow.
   ! ros_cDecomp is continuous.

   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   ENDIF

   ! map cNONZERO back to full space. this step is very slow.
   DO i = 1, cNONZERO
      Ghimj(JVS_MAP(i)) = cGhimj(i)
   ENDDO
   ! Ghimj(JVS_MAP(1:cNONZERO)) = cGhimj
   ! The above implementation seems to make -fcheck=bounds unhappy
   ! and cause a segfault at -O3. The asm jumps seemed weird, so we
   ! just explicitly write out the loop.

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1) ! note this is full space

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         call WCOPY(N,Fcn0,1,Fcn,1)
         ! Fcn(1:N) = Fcn0(1:N)
         ! istage>1 and a new function evaluation is needed at the current istage
         ! K = 0.0_dp ! is this fix needed? hplin 14:04 -- not. 3 hours wiser later
       ELSEIF ( ros_NewF(istage) ) THEN
         call WCOPY(N,Y,1,Ynew,1)
         ! Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
            ! In full vector space. Just use WAXPY as normal
            ! other entries in K are set to 1 previously.
            ! the rest are filled by Fcn.

            !          N  alpha                             x                y .... Y <- Y + a*X.
            ! i.e. Ynew <- Ynew + ros_A(..) * K. for !DO_FUN, we want K === 0
            ! if there are entries in x that are zero naturally they do not carry into Ynew.
            !
            ! otherwise Ynew will be wrongly updated

            ! order of operations here is from K(N*(j-1)+1:N*j+1), total of N elem.
            ! K for !DO_FUN should not be updated.
            ! in fact, K is reasonably sparse here.

            ! full version:
            ! CALL WAXPY(N, ros_A((istage-1)*(istage-2)/2+j), K(N*(j-1)+1), 1, Ynew, 1)

            ! only rNVAR version - maybe loops need to be unrolled: (15:39)
            alpha_factor = ros_A((istage-1)*(istage-2)/2+j)
            DO i = 1,rNVAR
               Ynew(SPC_MAP(i)) = Ynew(SPC_MAP(i)) + alpha_factor * K(N*(j-1)+SPC_MAP(i))
            ENDDO
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunSplitN(Tau,Ynew,Fcn)
         ! this step reacts to DO_FUN.
         ! Fcn is updated thru Fun_Split(Ynew,..,..,P,D), Fcn <- P - D*Ynew
         ! P and D are only updated iff. DO_FUN is set to true. otherwise they are zero
         ! the purpose is for Fcn to be zero for !DO_FUNs so K is set to 0 for them, so delta is 0.
         ! this looks okay then. 4/10/22 13:36 hplin
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       ! K(ioffset+1:ioffset+rNVAR) = Fcn(SPC_MAP(1:rNVAR))

      ! now operate on full space for all of below.
      ! full version:
      ! K(ioffset+1:ioffset+N) = Fcn(1:N)

      ! faster version: (this copy also feels expensive...)
      ! unroll j = 1 stage iff. istage-1>1, otherwise this is skipped.
      IF(istage .gt. 1) THEN
         HC = ros_C((istage-1)*(istage-2)/2+1)/(Direction*H)
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = Fcn(SPC_MAP(i)) + HC * K(SPC_MAP(i))
         ENDDO
      ENDIF

      IF(istage .eq. 1) THEN
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = Fcn(SPC_MAP(i))
         ENDDO
      ENDIF

      DO j = 2, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)

         ! full version
         ! K(ioffset+1:ioffset+1+N) <- K(ioffset+1:ioffset+1+N) + HC*K(N*(j-1)+1)
         ! CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
         ! K also updated here...
         ! write(6,*) "istage,kupd2",istage,K(ioffset+1:ioffset+1+N)

         ! faster version:
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = K(ioffset+SPC_MAP(i)) + HC * K(N*(j-1)+SPC_MAP(i))
         ENDDO
         ! CALL zWAXPY(N,HC,K(N*(j-1)+1),K(ioffset+1),SPC_MAP)
         ! loop unrolling is consistently slower here. 18:58
      ENDDO

      IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)

         ! full version: CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
         ! faster version:
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = K(ioffset+SPC_MAP(i)) + HG * dFdT(SPC_MAP(i))
         ENDDO
      ENDIF

      ! CALL ros_cSolve(Ghimj(1:cNONZERO), Pivot, K(ioffset+1), JVS_MAP, SPC_MAP)
      ! this means that ros_cSolve now does not need to be mapped back to full space here
      ! avoiding several full maps.
      !
      ! K is passed to ros_Solve -> KppSolve, which responds to DO_SLV.
      ! If DO_SLV is 0, K is not updated at that point
      ! Note other terms may still depend on other JVS(LU_NONZERO) or K(N) terms ...

      ! this fix is generally necessary. this means that before K gets to
      ! ros_Solve, the non-DO_SLV terms need to be zeroed out.
      !
      ! of course this is inefficient, so care should be taken to only update K
      ! where necessary without doing this final sweep.
      ! each one of these sweeps costs approximately 2% of total time.
      ! DO j = 1, N
      !    IF (.not. DO_SLV(j)) K(ioffset+j) = 0.0_dp
      ! ENDDO
      ! ... 2.90% (39%)
      CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))
      ! ... 3.10% (40%)

      ! note in ros_cSolve -> kppSolve the back-substitution for !DO_SLV==!DO_FUN are
      ! not resolved so Ynew remains the same.

   ENDDO Stage

   ! ... 4.44%
! roll the new solution and error estimation into one loop.
   Ynew(:) = Y(:)
   Yerr(:) = ZERO
   DO j = 1, ros_S
      DO i = 1,rNVAR
         Ynew(SPC_MAP(i)) = Ynew(SPC_MAP(i)) + ros_M(j) * K(N*(j-1)+SPC_MAP(i))
         Yerr(SPC_MAP(i)) = Yerr(SPC_MAP(i)) + ros_E(j) * K(N*(j-1)+SPC_MAP(i))
      ENDDO
   ENDDO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
   ! ... 4.73% (~0.3%). i will leave this alone for now

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      CALL WCOPY(N,Ynew,1,Y,1)
      !Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      ENDIF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew

      ! write(6,*) "STEP ACCEPTED"

      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
      ! write(6,*) "STEP REJECTED"
   END IF ! Err <= 1

   END DO UntilAccepted

   ! ... 4.70%
   END DO TimeLoop

   ! 1st order calculation for removed species per Shen et al. (2020) Eq. 4
   ! -- currently, DO_FUN() selects
   ! -- DO_FUN loops over 1,NVAR. Only needs to loop over NVAR-rNVAR
   !    but the structure doesn't exist. Maybe worth considering 
   !    for efficiency purposes.
   DO i=1,N
      IF (.not. DO_SLV(i)) THEN
           call autoreduce_1stOrder(i,Y(i),Prod(i),Loss(i),Tstart,Tend)
      ENDIF
   ENDDO

!~~~> Successful exit
   IERR = 1  !~~~> The integration was successful

 END SUBROUTINE ros_yIntegrator
 ! Aggressively micro-optimized by hplin
 ! 4/10/22

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_yIntegratorA (Y, Tstart, Tend, T,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Autoreduce threshold
        threshold, AR_target_spc, AR_thr_ratio,  &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Alternative micro-optimized implementation, hplin, 4/10/22
! which does not resize arrays to rNVAR, instead always keeping to full N size
! and skipping using DO_SLV
!
! All compiled assembly code was verified just short of linking to a proper BLAS.

   USE KPP_ROOT_Global,  ONLY : cNONZERO, rNVAR
   use KPP_ROOT_Monitor, ONLY : SPC_NAMES
   USE KPP_ROOT_JacobianSP

  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   KPP_REAL, INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval
   KPP_REAL, INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   KPP_REAL, INTENT(OUT) ::  T
!~~~> Input: tolerances
   KPP_REAL, INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   KPP_REAL, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   KPP_REAL, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Autoreduction threshold
   KPP_REAL, INTENT(IN) :: threshold
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables
   KPP_REAL :: Ynew(N), Fcn0(N), Fcn(N), Prod(N), Loss(N), LossY(N), Prd0(N), Los0(N)
   KPP_REAL :: K(NVAR*ros_S), dFdT(N)
#ifdef FULL_ALGEBRA
   KPP_REAL :: Jac0(N,N), Ghimj(N,N)
#else
   KPP_REAL :: Jac0(LU_NONZERO), Ghimj(LU_NONZERO)
   KPP_REAL :: cGhimj(LU_NONZERO) ! not known at this point what cNONZERO will be
#endif
   KPP_REAL :: H, Hnew, HC, HG, Fac, Tau
   KPP_REAL :: Err, Yerr(N), Yerrsub(NVAR)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular, Reduced
!~~~>  Local parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   KPP_REAL, PARAMETER :: DeltaMin = 1.0E-5_dp
   INTEGER :: SPC
   KPP_REAL :: alpha_factor ! hplin 4/10/22
!      Inline local parameters for AR.
   INTEGER :: II, III, idx, nrmv, s
   KPP_REAL, INTENT(IN) :: AR_thr_ratio
   KPP_REAL :: AR_thr
   INTEGER, INTENT(IN) :: AR_target_spc
!~~~>  Initial preparations
   DO_SLV  = .true.
   DO_FUN  = .true.
   DO_JVS  = .true.
   Reduced = .false.
   RMV     = 0

   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

   ! reset K - hplin. this is a multiplier that gets applied to ros_A and ros_M
   K = 0.0_dp
   Ghimj = 0.0_dp

!~~~> Time loop begins below

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   ENDIF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   ENDIF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
! this is necessary anyway for the rest of the time loop. do not optimize out.
   CALL FunSplitF(T,Y,Fcn0,Prod,Loss,LossY) ! always calculates PL.
   ! The above Prod, Loss were updated at every TimeLoop, so they no longer
   ! reflect initial condition. The copy operation is extremely expensive,
   ! so we instead run it once at the beginning of the timeloop and use a IF
   ! to switch it to a wholy different codepath...
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

   ! ... 1.02% ... 13.00%

!~~~>  Parse species for reduced computation
   if (.not. reduced) then
      Prd0 = Prod
      Los0 = Loss ! Save original prod loss vector for 1st order approx. do not optimize out.

      ! Inline the entire reduction operation here.
      iSPC_MAP = 0
      NRMV     = 0
      S        = 1
      AR_thr   = threshold

      ! Target species?
      if(AR_target_spc .gt. 0) then
          AR_thr = AR_thr_ratio * max(LossY(AR_target_spc), Prod(AR_target_spc))           ! Lin et al., 2022 in prep.
          RSTATUS(NARthr) = AR_thr
      endif

      ! Checks should be kept out of tight inner loops.
      IF(keepActive) THEN
       DO i=1,NVAR
         ! Short-circuiting using SKIP is very important here.
         if (.not. keepSpcActive(i) .and. &
             abs(LossY(i)).lt.AR_thr .and. abs(Prod(i)).lt.AR_thr) then ! per Shen et al., 2020
            NRMV=NRMV+1
            RMV(NRMV) = i
            DO_SLV(i) = .false.
            ! DO_FUN(i) = .false.
            cycle
         endif
         SPC_MAP(S)  = i ! Add to full spc map.
         iSPC_MAP(i) = S
         S=S+1
       ENDDO
      ENDIF

      IF (.not. keepActive) THEN
        DO i=1,NVAR
         ! Short-circuiting using SKIP is very important here.
         if (abs(LossY(i)).lt.AR_thr .and. abs(Prod(i)).lt.AR_thr) then ! per Shen et al., 2020
            NRMV=NRMV+1
            RMV(NRMV) = i
            DO_SLV(i) = .false.
            ! DO_FUN(i) = .false.
            cycle
         endif
         SPC_MAP(S)  = i ! Add to full spc map.
         iSPC_MAP(i) = S
         S=S+1
       ENDDO
      ENDIF

      rNVAR    = NVAR-NRMV ! Number of active species in the reduced mechanism
      II  = 1
      III = 1
      idx = 0
      ! This loop can be unrolled into two branches.
      ScanFirstNonZero: DO i = 1,LU_NONZERO
         IF ((DO_SLV(LU_IROW(i))).and.(DO_SLV(LU_ICOL(i)))) THEN
            idx = 1
            cLU_IROW(1) = iSPC_MAP(LU_IROW(i))
            cLU_ICOL(1) = iSPC_MAP(LU_ICOL(i))
            JVS_MAP(1)  = i
            EXIT ScanFirstNonZero
         ENDIF
         DO_JVS(i) = .false.
      ENDDO ScanFirstNonZero
      DO i = i+1,LU_NONZERO ! There is no escape for looping through LU_NONZERO here.
         IF ((DO_SLV(LU_IROW(i))).and.(DO_SLV(LU_ICOL(i)))) THEN
            idx=idx+1 ! counter for the number of non-zero elements in the reduced Jacobian
            cLU_IROW(idx) = iSPC_MAP(LU_IROW(i))
            cLU_ICOL(idx) = iSPC_MAP(LU_ICOL(i))
            JVS_MAP(idx)  = i

            IF (cLU_IROW(idx).ne.cLU_IROW(idx-1)) THEN
               II=II+1
               cLU_CROW(II) = idx
            ENDIF
            IF (cLU_IROW(idx).eq.cLU_ICOL(idx)) THEN
               III=III+1
               cLU_DIAG(III) = idx
            ENDIF
            CYCLE
         ENDIF
         DO_JVS(i) = .false.
      ENDDO

     cNONZERO = idx
     cLU_CROW(1)       = 1 ! 1st index = 1
     cLU_DIAG(1)       = 1 ! 1st index = 1
     cLU_CROW(rNVAR+1) = cNONZERO+1
     cLU_DIAG(rNVAR+1) = cLU_DIAG(rNVAR)+1

     ! this loop approximately 0.6% ... 7.8%
     reduced = .true.
   endif ! not reduced

   if (reduced) then
      ! Scan prod/loss for condition change for append functionality.
      ! Note because this requires internal Prod/Loss update, yIntegrator will fall through
      ! here. (hplin, 4/10/22)
      ! RMV is filled at most to NRMV, so read to that. Note that NRMV is only here because
      ! we inlined Reduce(), so this cannot be ported to ros_cIntegrator. (hplin, 4/12/22)
      DO i=1,NRMV
         SPC = RMV(i)
         if (SPC .eq. 0) cycle ! Species is already appended
         if (abs(LossY(SPC)).gt.threshold .or. abs(Prod(SPC)).gt.threshold) then
            CALL APPEND(SPC)
            RMV(i) = 0
         endif
      ENDDO
   endif ! reduced. append functionality.

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

   ! 1.86 ~ 1.90%

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate(T,Y,Jac0) ! Reacts to DO_JVS()
   ISTATUS(Njac) = ISTATUS(Njac) + 1
!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

   CALL ros_cPrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,cGhimj,Pivot,Singular) ! calculate Ghimj as cGhimj(cNONZERO).
   ! cGhimj all elem are rewritten so no need to zero out - warning - hplin 4/10/22
   !
   ! this step has one skip addressing step, -Jac0(JVS_MAP(1:cNONZERO)). slow.
   ! ros_cDecomp is continuous.

   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   ENDIF

   ! map cNONZERO back to full space. this step is very slow.
   DO i = 1, cNONZERO
      Ghimj(JVS_MAP(i)) = cGhimj(i)
   ENDDO
   ! Ghimj(JVS_MAP(1:cNONZERO)) = cGhimj
   ! The above implementation seems to make -fcheck=bounds unhappy
   ! and cause a segfault at -O3. The asm jumps seemed weird, so we
   ! just explicitly write out the loop.

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1) ! note this is full space

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         call WCOPY(N,Fcn0,1,Fcn,1)
         ! Fcn(1:N) = Fcn0(1:N)
         ! istage>1 and a new function evaluation is needed at the current istage
         ! K = 0.0_dp ! is this fix needed? hplin 14:04 -- not. 3 hours wiser later
       ELSEIF ( ros_NewF(istage) ) THEN
         call WCOPY(N,Y,1,Ynew,1)
         ! Ynew(1:N) = Y(1:N)
         DO j = 1, istage-1
            ! In full vector space. Just use WAXPY as normal
            ! other entries in K are set to 1 previously.
            ! the rest are filled by Fcn.

            !          N  alpha                             x                y .... Y <- Y + a*X.
            ! i.e. Ynew <- Ynew + ros_A(..) * K. for !DO_FUN, we want K === 0
            ! if there are entries in x that are zero naturally they do not carry into Ynew.
            !
            ! otherwise Ynew will be wrongly updated

            ! order of operations here is from K(N*(j-1)+1:N*j+1), total of N elem.
            ! K for !DO_FUN should not be updated.
            ! in fact, K is reasonably sparse here.

            ! full version:
            ! CALL WAXPY(N, ros_A((istage-1)*(istage-2)/2+j), K(N*(j-1)+1), 1, Ynew, 1)

            ! only rNVAR version - maybe loops need to be unrolled: (15:39)
            alpha_factor = ros_A((istage-1)*(istage-2)/2+j)
            DO i = 1,rNVAR
               Ynew(SPC_MAP(i)) = Ynew(SPC_MAP(i)) + alpha_factor * K(N*(j-1)+SPC_MAP(i))
            ENDDO
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunSplitN(Tau,Ynew,Fcn)
         ! this step reacts to DO_FUN.
         ! Fcn is updated thru Fun_Split(Ynew,..,..,P,D), Fcn <- P - D*Ynew
         ! P and D are only updated iff. DO_FUN is set to true. otherwise they are zero
         ! the purpose is for Fcn to be zero for !DO_FUNs so K is set to 0 for them, so delta is 0.
         ! this looks okay then. 4/10/22 13:36 hplin
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       ! K(ioffset+1:ioffset+rNVAR) = Fcn(SPC_MAP(1:rNVAR))

      ! now operate on full space for all of below.
      ! full version:
      ! K(ioffset+1:ioffset+N) = Fcn(1:N)

      ! faster version: (this copy also feels expensive...)
      ! unroll j = 1 stage iff. istage-1>1, otherwise this is skipped.
      IF(istage .gt. 1) THEN
         HC = ros_C((istage-1)*(istage-2)/2+1)/(Direction*H)
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = Fcn(SPC_MAP(i)) + HC * K(SPC_MAP(i))
         ENDDO
      ENDIF

      IF(istage .eq. 1) THEN
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = Fcn(SPC_MAP(i))
         ENDDO
      ENDIF

      DO j = 2, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)

         ! full version
         ! K(ioffset+1:ioffset+1+N) <- K(ioffset+1:ioffset+1+N) + HC*K(N*(j-1)+1)
         ! CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
         ! K also updated here...
         ! write(6,*) "istage,kupd2",istage,K(ioffset+1:ioffset+1+N)

         ! faster version:
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = K(ioffset+SPC_MAP(i)) + HC * K(N*(j-1)+SPC_MAP(i))
         ENDDO
         ! CALL zWAXPY(N,HC,K(N*(j-1)+1),K(ioffset+1),SPC_MAP)
         ! loop unrolling is consistently slower here. 18:58
      ENDDO

      IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)

         ! full version: CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
         ! faster version:
         DO i = 1,rNVAR
            K(ioffset+SPC_MAP(i)) = K(ioffset+SPC_MAP(i)) + HG * dFdT(SPC_MAP(i))
         ENDDO
      ENDIF

      ! CALL ros_cSolve(Ghimj(1:cNONZERO), Pivot, K(ioffset+1), JVS_MAP, SPC_MAP)
      ! this means that ros_cSolve now does not need to be mapped back to full space here
      ! avoiding several full maps.
      !
      ! K is passed to ros_Solve -> KppSolve, which responds to DO_SLV.
      ! If DO_SLV is 0, K is not updated at that point
      ! Note other terms may still depend on other JVS(LU_NONZERO) or K(N) terms ...

      ! this fix is generally necessary. this means that before K gets to
      ! ros_Solve, the non-DO_SLV terms need to be zeroed out.
      !
      ! of course this is inefficient, so care should be taken to only update K
      ! where necessary without doing this final sweep.
      ! each one of these sweeps costs approximately 2% of total time.
      ! DO j = 1, N
      !    IF (.not. DO_SLV(j)) K(ioffset+j) = 0.0_dp
      ! ENDDO
      ! ... 2.90% (39%)
      CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))
      ! ... 3.10% (40%)

      ! note in ros_cSolve -> kppSolve the back-substitution for !DO_SLV==!DO_FUN are
      ! not resolved so Ynew remains the same.

   ENDDO Stage

! roll the new solution and error estimation into one loop.
   Ynew(:) = Y(:)
   Yerr(:) = ZERO
   DO j = 1, ros_S
      DO i = 1,rNVAR
         Ynew(SPC_MAP(i)) = Ynew(SPC_MAP(i)) + ros_M(j) * K(N*(j-1)+SPC_MAP(i))
         Yerr(SPC_MAP(i)) = Yerr(SPC_MAP(i)) + ros_E(j) * K(N*(j-1)+SPC_MAP(i))
      ENDDO
   ENDDO
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
   ! ... 4.73% (~0.3%). i will leave this alone for now

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      CALL WCOPY(N,Ynew,1,Y,1)
      !Y(1:N) = Ynew(1:N)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      ENDIF
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew

      ! write(6,*) "STEP ACCEPTED"

      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
      ! write(6,*) "STEP REJECTED"
   END IF ! Err <= 1

   END DO UntilAccepted

   ! ... 4.70%
   END DO TimeLoop

   ! 1st order calculation for removed species per Shen et al. (2020) Eq. 4
   ! -- currently, DO_FUN() selects
   ! -- DO_FUN loops over 1,NVAR. Only needs to loop over NVAR-rNVAR
   !    but the structure doesn't exist. Maybe worth considering 
   !    for efficiency purposes.
   DO i=1,N
      IF (.not. DO_SLV(i)) THEN
           call autoreduce_1stOrder(i,Y(i),Prd0(i),Los0(i),Tstart,Tend)
      ENDIF
   ENDDO

!~~~> Successful exit
   IERR = 1  !~~~> The integration was successful

 END SUBROUTINE ros_yIntegratorA

 SUBROUTINE AutoReduce_1stOrder(i,Y,P,k,Ti,Tf)
   KPP_REAL, INTENT(INOUT) :: Y
   KPP_REAL, INTENT(IN)    :: P,k,Ti,Tf
   INTEGER, INTENT(IN)          :: i
   KPP_REAL                :: term

   if (k .le. 1.d-30) return
   if (Y .le. 1.d-30) return
   term = P/k
   Y = term+(Y-term)*exp(-k*(Tf-Ti))
 END SUBROUTINE AutoReduce_1stOrder

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  KPP_REAL FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   KPP_REAL, INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   KPP_REAL :: Err, Scale, Ymax
   INTEGER  :: i
   KPP_REAL, PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   KPP_REAL, INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N)
!~~~> Output arguments
   KPP_REAL, INTENT(OUT) :: dFdT(N)
!~~~> Local variables
   KPP_REAL :: Delta
   KPP_REAL, PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FunTemplate( T+Delta, Y, dFdT )
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL WAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_cPrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(IN) ::  Jac0(N,N)
#else
   KPP_REAL, INTENT(IN) ::  Jac0(LU_NONZERO)
#endif   
   KPP_REAL, INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(OUT) :: Ghimj(N,N)
#else
   KPP_REAL, INTENT(OUT) :: Ghimj(cNONZERO)
#endif   
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   KPP_REAL, INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   KPP_REAL :: ghinv
   KPP_REAL, PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,rNVAR
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:cNONZERO) = -Jac0(JVS_MAP(1:cNONZERO))
     ghinv = ONE/(Direction*H*gam)
     DO i=1,rNVAR
       Ghimj(cLU_DIAG(i)) = Ghimj(cLU_DIAG(i))+ghinv
     END DO
#endif   
!~~~>    Compute LU decomposition
     CALL ros_cDecomp( Ghimj, Pivot, ISING )
     IF (ISING == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ISING .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ISING

   END DO ! WHILE Singular

 END SUBROUTINE ros_cPrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_cDecomp( A, Pivot, ISING )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(INOUT) :: A(N,N)
#else   
   KPP_REAL, INTENT(INOUT) :: A(cNONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ISING

#ifdef FULL_ALGEBRA
   CALL  DGETRF( N, N, A, N, Pivot, ISING )
#else   
   CALL cKppDecomp ( A, ISING )
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

 END SUBROUTINE ros_cDecomp
 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_cSolve( A, Pivot, b, map1, map2 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(IN) :: A(N,N)
   INTEGER :: ISING
#else   
   KPP_REAL, INTENT(IN) :: A(cNONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   KPP_REAL, INTENT(INOUT) :: b(rNVAR)
   INTEGER, INTENT(IN)          :: map1(LU_NONZERO), map2(NVAR)
   KPP_REAL                :: btmp(N), Atmp(LU_NONZERO)

#ifdef FULL_ALGEBRA
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ISING )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ISING=",ISING
   END IF  
#else   

   Atmp = 0.d0
   Btmp = 0.d0
   Atmp(map1(1:cNONZERO)) = A
   btmp(map2(1:rNVAR))    = b
!   call cWCOPY(cNONZERO,LU_NONZERO,A,1,Atmp,1,map1)
!   call cWCOPY(rNVAR,NVAR,B,1,Btmp,1,map2)
   CALL KppSolve( Atmp, btmp )
   b = btmp(map2(1:rNVAR))
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_cSolve

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(IN) ::  Jac0(N,N)
#else
   KPP_REAL, INTENT(IN) ::  Jac0(LU_NONZERO)
#endif
   KPP_REAL, INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(OUT) :: Ghimj(N,N)
#else
   KPP_REAL, INTENT(OUT) :: Ghimj(LU_NONZERO)
#endif
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   KPP_REAL, INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ising, Nconsecutive
   KPP_REAL :: ghinv
   KPP_REAL, PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
#ifdef FULL_ALGEBRA
     !slim: CALL WCOPY(N*N,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(N*N,(-ONE),Ghimj,1)
     Ghimj = -Jac0
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(i,i) = Ghimj(i,i)+ghinv
     END DO
#else
     !slim: CALL WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     !slim: CALL WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     Ghimj(1:LU_NONZERO) = -Jac0(1:LU_NONZERO)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,N
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO
#endif
!~~~>    Compute LU decomposition
     CALL ros_Decomp( Ghimj, Pivot, ising )
     IF (ising == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ising .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ising

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ising )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(INOUT) :: A(N,N)
#else
   KPP_REAL, INTENT(INOUT) :: A(LU_NONZERO)
#endif
!~~~> Output variables
   INTEGER, INTENT(OUT) :: Pivot(N), ising

#ifdef FULL_ALGEBRA
   CALL  DGETRF( N, N, A, N, Pivot, ising )
#else
   CALL KppDecomp ( A, ising )
   Pivot(1) = 1
#endif
   ISTATUS(Ndec) = ISTATUS(Ndec) + 1

  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, Pivot, b )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution
!  (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
#ifdef FULL_ALGEBRA
   KPP_REAL, INTENT(IN) :: A(N,N)
   INTEGER :: ising
#else
   KPP_REAL, INTENT(IN) :: A(LU_NONZERO)
#endif
   INTEGER, INTENT(IN) :: Pivot(N)
!~~~> InOut variables
   KPP_REAL, INTENT(INOUT) :: b(N)

#ifdef FULL_ALGEBRA
   CALL  DGETRS( 'N', N , 1, A, N, Pivot, b, N, ising )
   IF ( Info < 0 ) THEN
      PRINT*,"Error in DGETRS. ising=",ising
   END IF
#else
   CALL KppSolve( A, b )
#endif

   ISTATUS(Nsol) = ISTATUS(Nsol) + 1

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;


!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   End of the set of internal Rosenbrock subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
END SUBROUTINE Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE FunTemplate( T, Y, Ydot, P_VAR, D_VAR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE function call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 USE KPP_ROOT_Parameters, ONLY : NVAR, LU_NONZERO
 USE KPP_ROOT_Global,     ONLY : FIX, RCONST, TIME
 USE KPP_ROOT_Function,   ONLY : Fun, Fun_SPLIT
 USE KPP_ROOT_Rates,      ONLY : Update_SUN, Update_RCONST

!~~~> Input variables
   KPP_REAL :: T, Y(NVAR)
!~~~> Output variables
   KPP_REAL :: Ydot(NVAR)
   KPP_REAL, OPTIONAL :: P_VAR(NVAR), D_VAR(NVAR)
!~~~> Local variables
   KPP_REAL :: Told, P(NVAR), D(NVAR)

   Told = TIME
   TIME = T
   IF ( Do_Update_SUN    ) CALL Update_SUN()
   IF ( Do_Update_RCONST ) CALL Update_RCONST(Y)
   CALL KPP_FUN_OR_FUN_SPLIT
   TIME = Told

   IF (Present(P_VAR)) P_VAR=P
   IF (Present(D_VAR)) D_VAR=D

END SUBROUTINE FunTemplate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE FunSplitF( T, Y, Ydot, P_VAR, D_VAR, DY_VAR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE function call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!  This version does not react to DO_FUN even within autoreduce.
!  It also does not have any optional arguments.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 USE KPP_ROOT_Parameters, ONLY: NVAR, LU_NONZERO
 USE KPP_ROOT_Global, ONLY: FIX, RCONST, TIME
 USE KPP_ROOT_Function, ONLY: Fun_SPLITF
!~~~> Input variables
   KPP_REAL :: T, Y(NVAR)
!~~~> Output variables
   KPP_REAL :: Ydot(NVAR)
   KPP_REAL :: P_VAR(NVAR), D_VAR(NVAR), DY_VAR(NVAR)
!~~~> Local variables
   KPP_REAL :: Told, P(NVAR), D(NVAR)
   P    = 0.d0
   D    = 0.d0
   Told = TIME
   TIME = T
   CALL Fun_SPLITF( Y, FIX, RCONST, P, D )
   DY_VAR = D*y ! this can be used later.
   Ydot = P - DY_VAR
   TIME = Told

   P_VAR=P
   D_VAR=D

 END SUBROUTINE FunSplitF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE FunSplitN( T, Y, Ydot)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE function call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!  This version does not have any optional arguments.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 USE KPP_ROOT_Parameters, ONLY: NVAR, LU_NONZERO
 USE KPP_ROOT_Global, ONLY: FIX, RCONST, TIME
 USE KPP_ROOT_Function, ONLY: Fun_SPLITF
!~~~> Input variables
   KPP_REAL :: T, Y(NVAR)
!~~~> Output variables
   KPP_REAL :: Ydot(NVAR)
   KPP_REAL :: P_VAR(NVAR), D_VAR(NVAR)
!~~~> Local variables
   KPP_REAL :: Told, P(NVAR), D(NVAR)

   P    = 0.d0
   D    = 0.d0
   Told = TIME
   TIME = T
   CALL Fun_SPLITF( Y, FIX, RCONST, P, D )
   Ydot = P - D*y
   TIME = Told

 END SUBROUTINE FunSplitN

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE JacTemplate( T, Y, Jcb )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE Jacobian call.
!  Updates the rate coefficients (and possibly the fixed species) at each call
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 USE KPP_ROOT_Parameters,    ONLY : NVAR, LU_NONZERO
 USE KPP_ROOT_Global,        ONLY : FIX, RCONST, TIME
 USE KPP_ROOT_Jacobian,      ONLY : Jac_SP, LU_IROW, LU_ICOL
 USE KPP_ROOT_LinearAlgebra
 USE KPP_ROOT_Rates,         ONLY : Update_SUN, Update_RCONST
!~~~> Input variables
    KPP_REAL :: T, Y(NVAR)
!~~~> Output variables
#ifdef FULL_ALGEBRA
    KPP_REAL :: JV(LU_NONZERO), Jcb(NVAR,NVAR)
#else
    KPP_REAL :: Jcb(LU_NONZERO)
#endif
!~~~> Local variables
    KPP_REAL :: Told
#ifdef FULL_ALGEBRA
    INTEGER :: i, j
#endif

    Told = TIME
    TIME = T
    IF ( Do_Update_SUN    ) CALL Update_SUN()
    IF ( Do_Update_RCONST ) CALL Update_RCONST(Y)
#ifdef FULL_ALGEBRA
    CALL Jac_SP(Y, FIX, RCONST, JV)
    DO j=1,NVAR
      DO i=1,NVAR
         Jcb(i,j) = 0.0_dp
      END DO
    END DO
    DO i=1,LU_NONZERO
       Jcb(LU_IROW(i),LU_ICOL(i)) = JV(i)
    END DO
#else
    CALL Jac_SP( Y, FIX, RCONST, Jcb )
#endif
    TIME = Told

END SUBROUTINE JacTemplate

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE cKppDecomp( JVS, IER )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Sparse LU factorization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER  :: IER
      KPP_REAL :: JVS(cNONZERO), W(rNVAR), a
      INTEGER  :: k, kk, j, jj

      a = 0. ! mz_rs_20050606
      IER = 0
      DO k=1,rNVAR
        ! mz_rs_20050606: don't check if real value == 0
        ! IF ( JVS( LU_DIAG(k) ) .EQ. 0. ) THEN
        IF ( ABS(JVS(cLU_DIAG(k))) < TINY(a) ) THEN
            IER = k
            RETURN
        END IF
        DO kk = cLU_CROW(k), cLU_CROW(k+1)-1
              W( cLU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = cLU_CROW(k), cLU_DIAG(k)-1
            j = cLU_ICOL(kk)
            a = -W(j) / JVS( cLU_DIAG(j) )
            W(j) = -a
            DO jj = cLU_DIAG(j)+1, cLU_CROW(j+1)-1
               W( cLU_ICOL(jj) ) = W( cLU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = cLU_CROW(k), cLU_CROW(k+1)-1
            JVS(kk) = W( cLU_ICOL(kk) )
         END DO
      END DO
      
END SUBROUTINE cKppDecomp

SUBROUTINE APPEND(IDX)
  USE KPP_ROOT_JacobianSP
  ! Reactivate a deactivated species
  INTEGER, INTENT(IN) :: IDX ! Index of deactivated KPP species to append
  INTEGER :: I
  ! set the do_* logicals
  DO_SLV(IDX) = .true.
  DO_FUN(IDX) = .true.
  ! increment rNVAR
  rNVAR    = rNVAR+1 
  ! append SPC_MAP & iSPC_MAP
  SPC_MAP(rNVAR) = IDX ! From AR to full species
  iSPC_MAP(IDX)  = rNVAR ! From full to AR species
  ! -- the following requires scanning LU_NONZERO elements
  DO I = 1, LU_NONZERO
     IF (LU_IROW(i).eq.IDX .and. DO_SLV(LU_ICOL(i))) THEN ! TERM IS ACTIVE
        cNONZERO = cNONZERO+1 ! Add a non-zero term
        ! append cLU_IROW
        ! append cLU_ICOL
        ! append JVS_MAP
        cLU_IROW(cNONZERO) = iSPC_MAP(LU_IROW(I))
        cLU_ICOL(cNONZERO) = iSPC_MAP(LU_ICOL(I))
        JVS_MAP(cNONZERO)  = I
        DO_JVS(I)          = .true.
        ! append cLU_CROW
        ! append cLU_DIAG
        IF (cLU_IROW(cNONZERO).ne.cLU_IROW(cNONZERO-1)) THEN
           cLU_CROW(rNVAR) = cNONZERO
        ENDIF
        IF (cLU_IROW(cNONZERO).eq.cLU_ICOL(cNONZERO)) THEN
           cLU_DIAG(rNVAR) = cNONZERO
        ENDIF
     ENDIF
  ENDDO
  cLU_CROW(rNVAR+1) = cNONZERO+1
  cLU_DIAG(rNVAR+1) = cLU_DIAG(rNVAR)+1
END SUBROUTINE APPEND


END MODULE KPP_ROOT_Integrator
