!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  RungeKuttaTLM - Tangent Linear Model of Fully Implicit 3-stage RK:     !
!          * Radau-2A   quadrature (order 5)                              !
!          * Radau-1A   quadrature (order 5)                              !
!          * Lobatto-3C quadrature (order 4)                              !
!          * Gauss      quadrature (order 6)                              !
!  By default the code employs the KPP sparse linear algebra routines     !
!  Compile with -DFULL_ALGEBRA to use full linear algebra (LAPACK)        !
!                                                                         !
!    (C)  Adrian Sandu, August 2005                                       !
!    Virginia Polytechnic Institute and State University                  !
!    Contact: sandu@cs.vt.edu                                             !
!    Revised by Philipp Miehe and Adrian Sandu, May 2006                  !
!    This implementation is part of KPP - the Kinetic PreProcessor        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE KPP_ROOT_Integrator

  USE KPP_ROOT_Precision
  USE KPP_ROOT_Parameters, ONLY: NVAR, NSPEC, NFIX, LU_NONZERO
  USE KPP_ROOT_Global, ONLY: FIX, RCONST, TIME
  USE KPP_ROOT_Jacobian
  USE KPP_ROOT_LinearAlgebra

  IMPLICIT NONE
  PUBLIC
  SAVE

!~~~>  Statistics on the work performed by the Runge-Kutta method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
    Nrej=5, Ndec=6, Nsol=7, Nsng=8, Ntexit=1, Nhacc=2, Nhnew=3

CONTAINS

  ! **************************************************************************

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INTEGRATE_TLM( NTLM, Y, Y_tlm, TIN, TOUT, ATOL_tlm, RTOL_tlm, &
       ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    USE KPP_ROOT_Parameters, ONLY: NVAR
    USE KPP_ROOT_Global,     ONLY: ATOL,RTOL,VAR

    IMPLICIT NONE

!~~~> Y - Concentrations
    KPP_REAL :: Y(NVAR)
!~~~> NTLM - No. of sensitivity coefficients
    INTEGER NTLM
!~~~> Y_tlm - Sensitivities of concentrations
!     Note: Y_tlm (1:NVAR,j) contains sensitivities of
!               Y(1:NVAR) w.r.t. the j-th parameter, j=1...NTLM
    KPP_REAL :: Y_tlm(NVAR,NTLM)
    KPP_REAL :: TIN  ! TIN - Start Time
    KPP_REAL :: TOUT ! TOUT - End Time
    KPP_REAL, INTENT(IN),  OPTIONAL :: RTOL_tlm(NVAR,NTLM),ATOL_tlm(NVAR,NTLM)
    INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
    KPP_REAL, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
    INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
    KPP_REAL, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
    INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

    INTEGER :: IERR

    KPP_REAL :: RCNTRL(20), RSTATUS(20), T1, T2
    INTEGER :: ICNTRL(20), ISTATUS(20)
    INTEGER, SAVE :: Ntotal = 0

    ICNTRL(1:20) = 0
    RCNTRL(1:20) = 0.0_dp

    !~~~> fine-tune the integrator:
    ICNTRL(2)  = 0   ! Tolerances: 0=vector, 1=scalar
    ICNTRL(5)  = 8   ! Max no. of Newton iterations
    ICNTRL(6)  = 0   ! Starting values for Newton are: 0=interpolated, 1=zero
    ICNTRL(7)  = 1   ! TLM solution: 0=iterations, 1=direct
    ICNTRL(9)  = 0   ! TLM Newton iteration error: 0=fwd Newton, 1=TLM Newton error est
    ICNTRL(10) = 1   ! FWD error estimation: 0=classic, 1=SDIRK
    ICNTRL(11) = 1   ! Step size ontroller: 1=Gustaffson, 2=classic
    ICNTRL(12) = 0   ! Trunc. error estimate: 0=fwd only, 1=fwd and TLM

    !~~~> if optional parameters are given, and if they are >0,
    !     then use them to overwrite default settings
    IF (PRESENT(ICNTRL_U)) THEN
      WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
    END IF
    IF (PRESENT(RCNTRL_U)) THEN
      WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
    END IF

    T1 = TIN; T2 = TOUT
    CALL RungeKuttaTLM(  NVAR, NTLM, Y, Y_tlm, T1, T2, RTOL, ATOL, &
                           RTOL_tlm, ATOL_tlm,                     &
			   RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR  )

    Ntotal = Ntotal + ISTATUS(Nstp)
    PRINT*,'NSTEPS=',ISTATUS(Nstp),' (',Ntotal,')','  O3=', VAR(ind_O3)

    ! if optional parameters are given for output
    ! use them to store information in them
    IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
    IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
    IF (PRESENT(IERR_U)) IERR_U = IERR

    IF (IERR < 0) THEN
      PRINT *,'Runge-Kutta-TLM: Unsuccessful exit at T=', TIN,' (IERR=',IERR,')'
    ENDIF

  END SUBROUTINE INTEGRATE_TLM


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE RungeKuttaTLM(  N, NTLM, Y, Y_tlm, T, Tend,RelTol,AbsTol,    &
                             RelTol_tlm,AbsTol_tlm,RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  This implementation is based on the book and the code Radau5:
!
!         E. HAIRER AND G. WANNER
!         "SOLVING ORDINARY DIFFERENTIAL EQUATIONS II. 
!              STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS."
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
!         SPRINGER-VERLAG (1991)
!
!         UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!         CH-1211 GENEVE 24, SWITZERLAND
!         E-MAIL:  HAIRER@DIVSUN.UNIGE.CH,  WANNER@DIVSUN.UNIGE.CH
!
!   Methods:
!          * Radau-2A   quadrature (order 5)                              
!          * Radau-1A   quadrature (order 5)                              
!          * Lobatto-3C quadrature (order 4)                              
!          * Gauss      quadrature (order 6)                              
!                                                                         
!   (C)  Adrian Sandu, August 2005                                       
!   Virginia Polytechnic Institute and State University                  
!   Contact: sandu@cs.vt.edu                                             
!   Revised by Philipp Miehe and Adrian Sandu, May 2006                  
!   This implementation is part of KPP - the Kinetic PreProcessor        
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!       ----------------
!
!    Note: For input parameters equal to zero the default values of the
!          corresponding variables are used.
!
!     N           Dimension of the system
!     T           Initial time value
!
!     Tend        Final T value (Tend-T may be positive or negative)
!
!     Y(N)        Initial values for Y
!
!     RelTol,AbsTol   Relative and absolute error tolerances. 
!          for ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
!                        = 1: AbsTol, RelTol are scalars
!
!~~~>  Integer input parameters:
!  
!    ICNTRL(1) = not used
!
!    ICNTRL(2) = 0: AbsTol, RelTol are NVAR-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3) = RK method selection       
!              = 1:  Radau-2A    (the default)
!              = 2:  Lobatto-3C
!              = 3:  Gauss
!	       = 4:  Radau-1A
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0 the default value of 100000 is used
!
!    ICNTRL(5)  -> maximum number of Newton iterations
!        For ICNTRL(5)=0 the default value of 8 is used
!
!    ICNTRL(6)  -> starting values of Newton iterations:
!        ICNTRL(6)=0 : starting values are obtained from 
!                      the extrapolated collocation solution
!                      (the default)
!        ICNTRL(6)=1 : starting values are zero
!
!    ICNTRL(7)  -> method to solve TLM equations:
!        ICNTRL(7)=0 : modified Newton re-using LU (the default)
!        ICNTRL(7)=1 : direct solution (additional one *big* LU factorization)
!
!    ICNTRL(9) -> switch for TLM Newton iteration error estimation strategy
!		ICNTRL(9) = 0: base number of iterations as forward solution
!		ICNTRL(9) = 1: use RTOL_tlm and ATOL_tlm to calculate
!				error estimation for TLM at Newton stages
!
!    ICNTRL(10) -> switch for error estimation strategy
!		ICNTRL(10) = 0: one additional stage: at c=0, 
!				see Hairer (default)
!		ICNTRL(10) = 1: two additional stages: one at c=0 
!				and one SDIRK at c=1, stiffly accurate
!
!    ICNTRL(11) -> switch for step size strategy
!              ICNTRL(11)=1:  mod. predictive controller (Gustafsson, default)
!              ICNTRL(11)=2:  classical step size control
!              the choice 1 seems to produce safer results;
!              for simple problems, the choice 2 produces
!              often slightly faster runs
!
!    ICNTRL(12) -> switch for TLM truncation error control
!		ICNTRL(12) = 0: TLM error is not used
!		ICNTRL(12) = 1: TLM error is computed and used
!
!
!~~~>  Real input parameters:
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!                  (highly recommended to keep Hmin = ZERO, the default)
!
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!
!    RCNTRL(3)  -> Hstart, the starting step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                 (default=0.1)
!
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!                  than the predicted value  (default=0.9)
!
!    RCNTRL(8)  -> ThetaMin. If Newton convergence rate smaller
!                  than ThetaMin the Jacobian is not recomputed;
!                  (default=0.001)
!
!    RCNTRL(9)  -> NewtonTol, stopping criterion for Newton's method
!                  (default=0.03)
!
!    RCNTRL(10) -> Qmin
!
!    RCNTRL(11) -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
!                  step size is kept constant and the LU factorization
!                  reused (default Qmin=1, Qmax=1.2)
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
!    IERR        -> Reports on successfulness upon return:
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
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IMPLICIT NONE
      
      INTEGER       :: N, NTLM
      KPP_REAL :: Y(N),Y_tlm(N,NTLM)
      KPP_REAL :: AbsTol(N),RelTol(N),RCNTRL(20),RSTATUS(20)
      KPP_REAL :: AbsTol_tlm(N,NTLM),RelTol_tlm(N,NTLM)
      INTEGER :: ICNTRL(20), ISTATUS(20)
      LOGICAL :: StartNewton, Gustafsson, SdirkError, TLMNewtonEst, &
                       TLMtruncErr, TLMDirect
      INTEGER :: IERR, ITOL
      KPP_REAL :: T,Tend

      !~~~> Control arguments
      INTEGER :: Max_no_steps, NewtonMaxit, rkMethod
      KPP_REAL :: Hmin,Hmax,Hstart,Qmin,Qmax
      KPP_REAL :: Roundoff, ThetaMin, NewtonTol
      KPP_REAL :: FacSafe,FacMin,FacMax,FacRej
      ! Runge-Kutta method parameters
      INTEGER, PARAMETER :: R2A=1, R1A=2, L3C=3, GAU=4, L3A=5
      KPP_REAL :: rkT(3,3), rkTinv(3,3), rkTinvAinv(3,3), rkAinvT(3,3), &
                       rkA(3,3), rkB(3),  rkC(3), rkD(0:3), rkE(0:3), 	&
	   	       rkBgam(0:4), rkBhat(0:4), rkTheta(0:3),          &
                       rkGamma,  rkAlpha, rkBeta, rkELO
      !~~~> Local variables
      INTEGER :: i
      KPP_REAL, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        SETTING THE PARAMETERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IERR = 0
      ISTATUS(1:20) = 0
      RSTATUS(1:20) = ZERO
       
!~~~> ICNTRL(1) - autonomous system - not used       
!~~~> ITOL: 1 for vector and 0 for scalar AbsTol/RelTol
      IF (ICNTRL(2) == 0) THEN
         ITOL = 1
      ELSE
         ITOL = 0
      END IF
!~~~> Error control selection  
      IF (ICNTRL(10) == 0) THEN 
         SdirkError = .FALSE.
      ELSE
         SdirkError = .TRUE.
      END IF      
!~~~> Method selection  
      SELECT CASE (ICNTRL(3))     
      CASE (0,1)
         CALL Radau2A_Coefficients
      CASE (2)
         CALL Lobatto3C_Coefficients
      CASE (3)
         CALL Gauss_Coefficients
      CASE (4)
         CALL Radau1A_Coefficients
      CASE DEFAULT
         WRITE(6,*) 'ICNTRL(3)=',ICNTRL(3)
         CALL RK_ErrorMsg(-13,T,ZERO,IERR)
      END SELECT
!~~~> Max_no_steps: the maximal number of time steps
      IF (ICNTRL(4) == 0) THEN
         Max_no_steps = 200000
      ELSE
         Max_no_steps=ICNTRL(4)
         IF (Max_no_steps <= 0) THEN
            WRITE(6,*) 'ICNTRL(4)=',ICNTRL(4)
            CALL RK_ErrorMsg(-1,T,ZERO,IERR)
         END IF
      END IF
!~~~> NewtonMaxit    maximal number of Newton iterations
      IF (ICNTRL(5) == 0) THEN
         NewtonMaxit = 8
      ELSE
         NewtonMaxit=ICNTRL(5)
         IF (NewtonMaxit <= 0) THEN
            WRITE(6,*) 'ICNTRL(5)=',ICNTRL(5)
            CALL RK_ErrorMsg(-2,T,ZERO,IERR)
          END IF
      END IF
!~~~> StartNewton:  Use extrapolation for starting values of Newton iterations
      IF (ICNTRL(6) == 0) THEN
         StartNewton = .TRUE.
      ELSE
         StartNewton = .FALSE.
      END IF      
!~~~> Solve TLM equations directly or by Newton iterations 
      IF (ICNTRL(7) == 0) THEN
         TLMDirect = .FALSE.
      ELSE
         TLMDirect = .TRUE.
      END IF      
!~~~> Newton iteration error control selection  
      IF (ICNTRL(9) == 0) THEN 
         TLMNewtonEst = .FALSE.
      ELSE
         TLMNewtonEst = .TRUE.
      END IF      
!~~~> Gustafsson: step size controller
      IF(ICNTRL(11) == 0)THEN
         Gustafsson = .TRUE.
      ELSE
         Gustafsson = .FALSE.
      END IF
!~~~> TLM truncation error control selection  
      IF (ICNTRL(12) == 0) THEN 
         TLMtruncErr = .FALSE.
      ELSE
         TLMtruncErr = .TRUE.
      END IF      

!~~~> Roundoff: smallest number s.t. 1.0 + Roundoff > 1.0
      Roundoff=WLAMCH('E');

!~~~> Hmin = minimal step size
      IF (RCNTRL(1) == ZERO) THEN
         Hmin = ZERO
      ELSE
         Hmin = MIN(ABS(RCNTRL(1)),ABS(Tend-T))
      END IF
!~~~> Hmax = maximal step size
      IF (RCNTRL(2) == ZERO) THEN
         Hmax = ABS(Tend-T)
      ELSE
         Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-T))
      END IF
!~~~> Hstart = starting step size
      IF (RCNTRL(3) == ZERO) THEN
         Hstart = ZERO
      ELSE
         Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-T))
      END IF
!~~~> FacMin: lower bound on step decrease factor
      IF(RCNTRL(4) == ZERO)THEN
         FacMin = 0.2d0
      ELSE
         FacMin = RCNTRL(4)
      END IF
!~~~> FacMax: upper bound on step increase factor
      IF(RCNTRL(5) == ZERO)THEN
         FacMax = 8.D0
      ELSE
         FacMax = RCNTRL(5)
      END IF
!~~~> FacRej: step decrease factor after 2 consecutive rejections
      IF(RCNTRL(6) == ZERO)THEN
         FacRej = 0.1d0
      ELSE
         FacRej = RCNTRL(6)
      END IF
!~~~> FacSafe:  by which the new step is slightly smaller
!               than the predicted value
      IF (RCNTRL(7) == ZERO) THEN
         FacSafe=0.9d0
      ELSE
         FacSafe=RCNTRL(7)
      END IF
      IF ( (FacMax < ONE) .OR. (FacMin > ONE) .OR. &
           (FacSafe <= 1.0d-3) .OR. (FacSafe >= ONE) ) THEN
            WRITE(6,*)'RCNTRL(4:7)=',RCNTRL(4:7)
            CALL RK_ErrorMsg(-4,T,ZERO,IERR)
      END IF

!~~~> ThetaMin:  decides whether the Jacobian should be recomputed
      IF (RCNTRL(8) == ZERO) THEN
         ThetaMin = 1.0d-3
      ELSE
         ThetaMin=RCNTRL(8)
         IF (ThetaMin <= 0.0d0 .OR. ThetaMin >= 1.0d0) THEN
            WRITE(6,*) 'RCNTRL(8)=', RCNTRL(8)
            CALL RK_ErrorMsg(-5,T,ZERO,IERR)
         END IF
      END IF
!~~~> NewtonTol:  stopping crierion for Newton's method
      IF (RCNTRL(9) == ZERO) THEN
         NewtonTol = 3.0d-2
      ELSE
         NewtonTol = RCNTRL(9)
         IF (NewtonTol <= Roundoff) THEN
            WRITE(6,*) 'RCNTRL(9)=',RCNTRL(9)
            CALL RK_ErrorMsg(-6,T,ZERO,IERR)
         END IF
      END IF
!~~~> Qmin AND Qmax: IF Qmin < Hnew/Hold < Qmax then step size = const.
      IF (RCNTRL(10) == ZERO) THEN
         Qmin=1.D0
      ELSE
         Qmin=RCNTRL(10)
      END IF
      IF (RCNTRL(11) == ZERO) THEN
         Qmax=1.2D0
      ELSE
         Qmax=RCNTRL(11)
      END IF
      IF (Qmin > ONE .OR. Qmax < ONE) THEN
         WRITE(6,*) 'RCNTRL(10:11)=',Qmin,Qmax
         CALL RK_ErrorMsg(-7,T,ZERO,IERR)
      END IF
!~~~> Check if tolerances are reasonable
      IF (ITOL == 0) THEN
          IF (AbsTol(1) <= ZERO.OR.RelTol(1) <= 10.d0*Roundoff) THEN
              WRITE (6,*) 'AbsTol/RelTol=',AbsTol,RelTol 
              CALL RK_ErrorMsg(-8,T,ZERO,IERR)
          END IF
      ELSE
          DO i=1,N
          IF (AbsTol(i) <= ZERO.OR.RelTol(i) <= 10.d0*Roundoff) THEN
              WRITE (6,*) 'AbsTol/RelTol(',i,')=',AbsTol(i),RelTol(i)
              CALL RK_ErrorMsg(-8,T,ZERO,IERR)
          END IF
          END DO
      END IF

!~~~> Parameters are wrong
      IF (IERR < 0) RETURN

!~~~> Call the core method
      CALL RK_IntegratorTLM( N,NTLM,T,Tend,Y,Y_tlm,IERR )

   CONTAINS ! Internal procedures to RungeKutta



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_IntegratorTLM( N,NTLM,T,Tend,Y,Y_tlm,IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IMPLICIT NONE
!~~~> Arguments
      INTEGER,  INTENT(IN)         :: N, NTLM
      KPP_REAL, INTENT(IN)    :: Tend
      KPP_REAL, INTENT(INOUT) :: T, Y(N), Y_tlm(NVAR,NTLM)
      INTEGER,  INTENT(OUT)        :: IERR

!~~~> Local variables
#ifdef FULL_ALGEBRA
      KPP_REAL, DIMENSION(NVAR,NVAR)    ::   &
                             FJAC, E1, Jac1, Jac2, Jac3
      COMPLEX(kind=dp), DIMENSION(NVAR,NVAR) :: E2   
      KPP_REAL, DIMENSION(3*NVAR,3*NVAR)  :: Jbig, Ebig
      KPP_REAL, DIMENSION(3*NVAR)         :: Zbig
      INTEGER,       DIMENSION(3*NVAR)         :: IPbig
#else
      KPP_REAL, DIMENSION(LU_NONZERO)   ::   &
                             FJAC, E1, Jac1, Jac2, Jac3
      COMPLEX(kind=dp), DIMENSION(LU_NONZERO)  :: E2 
!  Next 3 commented lines for sparse big linear algebra:
!      KPP_REAL, DIMENSION(3,3,LU_NONZERO) :: Jbig
!      KPP_REAL, DIMENSION(3,NVAR)         :: Zbig
!      INTEGER, DIMENSION(3,NVAR)               :: IPbig
      KPP_REAL, DIMENSION(3*NVAR,3*NVAR)  :: Jbig, Ebig
      KPP_REAL, DIMENSION(3*NVAR)         :: Zbig
      INTEGER,       DIMENSION(3*NVAR)         :: IPbig
#endif                
      KPP_REAL, DIMENSION(NVAR) :: Z1,Z2,Z3,Z4,SCAL,DZ1,DZ2,DZ3,DZ4, &
                 G,TMP,SCAL_tlm
      KPP_REAL, DIMENSION(NVAR,NTLM) :: Z1_tlm,Z2_tlm,Z3_tlm
      KPP_REAL  :: CONT(NVAR,3), Tdirection,  H, Hacc, Hnew, Hold, Fac, &
                 FacGus, Theta, Err, ErrOld, NewtonRate, NewtonIncrement,    &
                 Hratio, Qnewton, NewtonPredictedErr,NewtonIncrementOld,    &
		 ThetaTLM, ThetaSD
      INTEGER :: i,j,IP1(NVAR),IP2(NVAR),NewtonIter, ISING, Nconsecutive, itlm
      INTEGER :: saveNiter, NewtonIterTLM, info
      LOGICAL :: Reject, FirstStep, SkipJac, NewtonDone, SkipLU
      
            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initial setting
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Tdirection = SIGN(ONE,Tend-T)
      H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , Hmax )
      IF (ABS(H) <= 10.d0*Roundoff) H = 1.0d-6
      H = SIGN(H,Tdirection)
      Hold      = H
      Reject    = .FALSE.
      FirstStep = .TRUE.
      SkipJac   = .FALSE.
      SkipLU    = .FALSE.
      IF ((T+H*1.0001D0-Tend)*Tdirection >= ZERO) THEN
         H = Tend-T
      END IF
      Nconsecutive = 0
      CALL RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( (Tend-T)*Tdirection - Roundoff > ZERO )

      IF ( .NOT.SkipLU ) THEN ! This time around skip the Jac update and LU
        !~~~> Compute the Jacobian matrix
        IF ( .NOT.SkipJac ) THEN
          CALL JAC_CHEM(T,Y,FJAC)
          ISTATUS(Njac) = ISTATUS(Njac) + 1
        END IF
        !~~~> Compute the matrices E1 and E2 and their decompositions
        CALL RK_Decomp(N,H,FJAC,E1,IP1,E2,IP2,ISING)
        IF (ISING /= 0) THEN
          ISTATUS(Nsng) = ISTATUS(Nsng) + 1; Nconsecutive = Nconsecutive + 1
          IF (Nconsecutive >= 5) THEN
            CALL RK_ErrorMsg(-12,T,H,IERR); RETURN
          END IF
          H=H*0.5d0; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
          CYCLE Tloop
        ELSE
          Nconsecutive = 0    
        END IF   
      END IF ! SkipLU
   
      ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      IF (ISTATUS(Nstp) > Max_no_steps) THEN
        PRINT*,'Max number of time steps is ',Max_no_steps
        CALL RK_ErrorMsg(-9,T,H,IERR); RETURN
      END IF
      IF (0.1D0*ABS(H) <= ABS(T)*Roundoff)  THEN
        CALL RK_ErrorMsg(-10,T,H,IERR); RETURN
      END IF
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Loop for the simplified Newton iterations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
      !~~~>  Starting values for Newton iteration
      IF ( FirstStep .OR. (.NOT.StartNewton) ) THEN
         CALL Set2zero(N,Z1)
         CALL Set2zero(N,Z2)
         CALL Set2zero(N,Z3)
      ELSE
         ! Evaluate quadratic polynomial
         CALL RK_Interpolate('eval',N,H,Hold,Z1,Z2,Z3,CONT)
      END IF
      
      !~~~>  Initializations for Newton iteration
      NewtonDone = .FALSE.
      Fac = 0.5d0 ! Step reduction if too many iterations
      
NewtonLoop:DO  NewtonIter = 1, NewtonMaxit  
 
            !~~~> Prepare the right-hand side
            CALL RK_PrepareRHS(N,T,H,Y,Z1,Z2,Z3,DZ1,DZ2,DZ3)
            
            !~~~> Solve the linear systems
            CALL RK_Solve( N,H,E1,IP1,E2,IP2,DZ1,DZ2,DZ3,ISING )
            
            NewtonIncrement = SQRT( ( RK_ErrorNorm(N,SCAL,DZ1)**2 + &
                                RK_ErrorNorm(N,SCAL,DZ2)**2 +       &
                                RK_ErrorNorm(N,SCAL,DZ3)**2 )/3.0d0 )
            
            IF ( NewtonIter == 1 ) THEN
                Theta      = ABS(ThetaMin)
                NewtonRate = 2.0d0 
            ELSE
                Theta = NewtonIncrement/NewtonIncrementOld
                IF (Theta < 0.99d0) THEN
                    NewtonRate = Theta/(ONE-Theta)
                ELSE ! Non-convergence of Newton: Theta too large
                    EXIT NewtonLoop
                END IF
                IF ( NewtonIter < NewtonMaxit ) THEN 
                  ! Predict error at the end of Newton process 
                  NewtonPredictedErr = NewtonIncrement &
                               *Theta**(NewtonMaxit-NewtonIter)/(ONE-Theta)
                  IF (NewtonPredictedErr >= NewtonTol) THEN
                    ! Non-convergence of Newton: predicted error too large
                    Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                    Fac=0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                    EXIT NewtonLoop
                  END IF
                END IF
            END IF

            NewtonIncrementOld = MAX(NewtonIncrement,Roundoff) 
            ! Update solution
            CALL WAXPY(N,-ONE,DZ1,1,Z1,1) ! Z1 <- Z1 - DZ1
            CALL WAXPY(N,-ONE,DZ2,1,Z2,1) ! Z2 <- Z2 - DZ2
            CALL WAXPY(N,-ONE,DZ3,1,Z3,1) ! Z3 <- Z3 - DZ3
            
            ! Check error in Newton iterations
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
            IF (NewtonDone) THEN
	      ! Tune error in TLM variables by defining the minimal number of Newton iterations.
	      saveNiter = NewtonIter+1
	      EXIT NewtonLoop
	    END IF
   	    IF (NewtonIter == NewtonMaxit) THEN
		PRINT*, 'Slow or no convergence in Newton Iteration: Max no. of', &
		 	'Newton iterations reached'
	    END IF
            
      END DO NewtonLoop
            
      IF (.NOT.NewtonDone) THEN
           !CALL RK_ErrorMsg(-12,T,H,IERR);
           H = Fac*H; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
           CYCLE Tloop
      END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> SDIRK Stage
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IF (SdirkError) THEN

!~~~>  Starting values for Newton iterations
       Z4(1:N) = Z3(1:N)
       
!~~~>   Prepare the loop-independent part of the right-hand side
       CALL FUN_CHEM(T,Y,DZ4)
       ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       
!       G = H*rkBgam(0)*DZ4 + rkTheta(1)*Z1 + rkTheta(2)*Z2 + rkTheta(3)*Z3
       CALL Set2Zero(N, G)
       CALL WAXPY(N,rkBgam(0)*H, DZ4,1,G,1) 
       CALL WAXPY(N,rkTheta(1),Z1,1,G,1)
       CALL WAXPY(N,rkTheta(2),Z2,1,G,1)
       CALL WAXPY(N,rkTheta(3),Z3,1,G,1)

       !~~~>  Initializations for Newton iteration
       NewtonDone = .FALSE.
       Fac = 0.5d0 ! Step reduction factor if too many iterations
            
SDNewtonLoop:DO NewtonIter = 1, NewtonMaxit

!~~~>   Prepare the loop-dependent part of the right-hand side
	    CALL WADD(N,Y,Z4,TMP)         ! TMP <- Y + Z4
	    CALL FUN_CHEM(T+H,TMP,DZ4) 	  ! DZ4 <- Fun(Y+Z4)         
            ISTATUS(Nfun) = ISTATUS(Nfun) + 1
!            DZ4(1:N) = (G(1:N)-Z4(1:N))*(rkGamma/H) + DZ4(1:N)
	    CALL WAXPY (N, -ONE*rkGamma/H, Z4, 1, DZ4, 1)
            CALL WAXPY (N, rkGamma/H, G,1, DZ4,1)

!~~~>   Solve the linear system
#ifdef FULL_ALGEBRA  
            CALL DGETRS( 'N', N, 1, E1, N, IP1, DZ4, N, ISING )
#else
            CALL KppSolve(E1, DZ4)
#endif
            
!~~~>   Check convergence of Newton iterations
            NewtonIncrement = RK_ErrorNorm(N,SCAL,DZ4)
            IF ( NewtonIter == 1 ) THEN
                ThetaSD      = ABS(ThetaMin)
                NewtonRate = 2.0d0 
            ELSE
                ThetaSD = NewtonIncrement/NewtonIncrementOld
                IF (ThetaSD < 0.99d0) THEN
                    NewtonRate = ThetaSD/(ONE-ThetaSD)
                    ! Predict error at the end of Newton process 
                    NewtonPredictedErr = NewtonIncrement &
                               *ThetaSD**(NewtonMaxit-NewtonIter)/(ONE-ThetaSD)
                    IF (NewtonPredictedErr >= NewtonTol) THEN
                      ! Non-convergence of Newton: predicted error too large
                      Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                      Fac = 0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                      EXIT SDNewtonLoop
                    END IF
                ELSE ! Non-convergence of Newton: Theta too large
                    print*,'Theta too large: ',ThetaSD
                    EXIT SDNewtonLoop
                END IF
            END IF
            NewtonIncrementOld = NewtonIncrement
            ! Update solution: Z4 <-- Z4 + DZ4
            CALL WAXPY(N,ONE,DZ4,1,Z4,1) 
            
            ! Check error in Newton iterations
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
            IF (NewtonDone) EXIT SDNewtonLoop
            
            END DO SDNewtonLoop
            
            IF (.NOT.NewtonDone) THEN
                H = Fac*H; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
                CYCLE Tloop
            END IF

!~~~>  End of implified SDIRK Newton iterations
      END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Error estimation, forward solution
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (SdirkError) THEN
!         DZ4(1:N) =  rkD(1)*Z1 + rkD(2)*Z2 + rkD(3)*Z3 - Z4    
	 CALL Set2Zero(N, DZ4)
	 IF (rkD(1) /= ZERO) CALL WAXPY(N, rkD(1), Z1, 1, DZ4, 1)
	 IF (rkD(2) /= ZERO) CALL WAXPY(N, rkD(2), Z2, 1, DZ4, 1)
	 IF (rkD(3) /= ZERO) CALL WAXPY(N, rkD(3), Z3, 1, DZ4, 1)
	 CALL WAXPY(N, -ONE, Z4, 1, DZ4, 1)
         Err = RK_ErrorNorm(N,SCAL,DZ4)    
      ELSE
        CALL  RK_ErrorEstimate(N,H,Y,T, &
               E1,IP1,Z1,Z2,Z3,SCAL,Err,FirstStep,Reject)
      END IF
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> If error small enough, compute TLM solution
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (Err < ONE) THEN

      !~~~>  Jacobian values
      CALL WADD(N,Y,Z1,TMP)              ! TMP  <- Y + Z1
      CALL JAC_CHEM(T+rkC(1)*H,TMP,Jac1) ! Jac1 <- Jac(Y+Z1)  
      CALL WADD(N,Y,Z2,TMP)              ! TMP  <- Y + Z2
      CALL JAC_CHEM(T+rkC(2)*H,TMP,Jac2) ! Jac2 <- Jac(Y+Z2)  
      CALL WADD(N,Y,Z3,TMP)              ! TMP  <- Y + Z3
      CALL JAC_CHEM(T+rkC(3)*H,TMP,Jac3) ! Jac3 <- Jac(Y+Z3)  
       
TLMDIR: IF (TLMDirect) THEN 
       
#ifdef FULL_ALGEBRA
      !~~~>  Construct the full Jacobian
      DO j=1,N
        DO i=1,N
          Jbig(i,j)         = -H*rkA(1,1)*Jac1(i,j)
          Jbig(N+i,j)       = -H*rkA(2,1)*Jac1(i,j)
          Jbig(2*N+i,j)     = -H*rkA(3,1)*Jac1(i,j)
          Jbig(i,N+j)       = -H*rkA(1,2)*Jac2(i,j)
          Jbig(N+i,N+j)     = -H*rkA(2,2)*Jac2(i,j)
          Jbig(2*N+i,N+j)   = -H*rkA(3,2)*Jac2(i,j)
          Jbig(i,2*N+j)     = -H*rkA(1,3)*Jac3(i,j)
          Jbig(N+i,2*N+j)   = -H*rkA(2,3)*Jac3(i,j)
          Jbig(2*N+i,2*N+j) = -H*rkA(3,3)*Jac3(i,j)
        END DO  
      END DO
      Ebig = -Jbig
      DO i=1,3*NVAR
          Jbig(i,i) = ONE + Jbig(i,i)
      END DO
      !~~~>  Solve the big system
      ! CALL DGETRF(3*NVAR,3*NVAR,Jbig,3*NVAR,IPbig,j) 
      CALL WGEFA(3*N,Jbig,IPbig,info) 
      IF (info /= 0) THEN
        PRINT*,'Big big guy is singular'; STOP
      END IF  
      DO itlm = 1, NTLM      
        DO j=1,NVAR
            Zbig(j)        = Y_tlm(j,itlm)
            Zbig(NVAR+j)   = Y_tlm(j,itlm)
            Zbig(2*NVAR+j) = Y_tlm(j,itlm)
        END DO
        Zbig = MATMUL(Ebig,Zbig)
        !CALL DGETRS ('N',3*NVAR,1,Jbig,3*NVAR,IPbig,Zbig,3*NVAR,ISING) 
        CALL WGESL('N',3*N,Jbig,IPbig,Zbig)
        DO j=1,NVAR
            Z1_tlm(j,itlm) = Zbig(j)
            Z2_tlm(j,itlm) = Zbig(NVAR+j)
            Z3_tlm(j,itlm) = Zbig(2*NVAR+j)
        END DO
      END DO  
#else       
!  Commented code for sparse big linear algebra:
!      !~~~>  Construct the big Jacobian
!      DO j=1,LU_NONZERO
!        DO i=1,3
!          Jbig(i,1,j) = -H*rkA(i,1)*Jac1(j)
!          Jbig(i,2,j) = -H*rkA(i,2)*Jac2(j)
!          Jbig(i,3,j) = -H*rkA(i,3)*Jac3(j)
!        END DO
!      END DO
!      DO j=1,NVAR
!        DO i=1,3
!          Jbig(i,i,LU_DIAG(j)) = ONE + Jbig(i,i,LU_DIAG(j))
!        END DO
!      END DO
!      !~~~>  Solve the big system
!      CALL KppDecompBig( Jbig, IPbig, j )
!   Use full big linear algebra:
      Jbig(1:3*N,1:3*N) = 0.0d0
      DO i=1,LU_NONZERO
	  Jbig(LU_irow(i),LU_icol(i))         = -H*rkA(1,1)*Jac1(i)
	  Jbig(LU_irow(i),N+LU_icol(i))       = -H*rkA(1,2)*Jac2(i)
	  Jbig(LU_irow(i),2*N+LU_icol(i))     = -H*rkA(1,3)*Jac3(i)
	  Jbig(N+LU_irow(i),LU_icol(i))       = -H*rkA(2,1)*Jac1(i)
	  Jbig(N+LU_irow(i),N+LU_icol(i))     = -H*rkA(2,2)*Jac2(i)
	  Jbig(N+LU_irow(i),2*N+LU_icol(i))   = -H*rkA(2,3)*Jac3(i)
	  Jbig(2*N+LU_irow(i),LU_icol(i))     = -H*rkA(3,1)*Jac1(i)
	  Jbig(2*N+LU_irow(i),N+LU_icol(i))   = -H*rkA(3,2)*Jac2(i)
	  Jbig(2*N+LU_irow(i),2*N+LU_icol(i)) = -H*rkA(3,3)*Jac3(i)
      END DO
      Ebig = -Jbig
      DO i=1, 3*N
	 Jbig(i,i) = ONE + Jbig(i,i)
      END DO
      ! CALL DGETRF(3*N,3*N,Jbig,3*N,IPbig,info) 
      CALL WGEFA(3*N,Jbig,IPbig,info) 
      IF (info /= 0) THEN
        PRINT*,'Big guy is singular'; STOP
      END IF  

      DO itlm = 1, NTLM      
!  Commented code for sparse big linear algebra:
!       ! Compute RHS
!        CALL RK_PrepareRHS_TLMdirect(N,H,Jac1,Jac2,Jac3,Y_tlm(1,itlm),Zbig)
!       ! Solve the system
!        CALL KppSolveBig( Jbig, IPbig, Zbig )
!        DO j=1,NVAR
!            Z1_tlm(j,itlm) = Zbig(1,j)
!            Z2_tlm(j,itlm) = Zbig(2,j)
!            Z3_tlm(j,itlm) = Zbig(3,j)
!        END DO
        ! Compute RHS
        CALL RK_PrepareRHS_TLMdirect(N,H,Jac1,Jac2,Jac3,Y_tlm(1,itlm),Zbig)
        ! Solve the system
	! CALL DGETRS('N',3*N,1,Jbig,3*N,IPbig,Zbig,3*N,ISING) 
        CALL WGESL('N',3*N,Jbig,IPbig,Zbig)
        Z1_tlm(1:NVAR,itlm) = Zbig(1:NVAR)
        Z2_tlm(1:NVAR,itlm) = Zbig(NVAR+1:2*NVAR)
        Z3_tlm(1:NVAR,itlm) = Zbig(2*NVAR+1:3*NVAR)
      END DO  
#endif

    ELSE TLMDIR
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Loop for Newton iterations, TLM variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
Tlm:DO itlm = 1, NTLM      
            
      !~~~>  Starting values for Newton iteration
      CALL Set2zero(N,Z1_tlm(1,itlm))
      CALL Set2zero(N,Z2_tlm(1,itlm))
      CALL Set2zero(N,Z3_tlm(1,itlm))
      
      !~~~>  Initializations for Newton iteration
      IF (TLMNewtonEst) THEN
         NewtonDone = .FALSE.
         Fac = 0.5d0 ! Step reduction if too many iterations

         CALL RK_ErrorScale(N,ITOL,AbsTol_tlm(1,itlm),RelTol_tlm(1,itlm), &
              Y_tlm(1,itlm),SCAL_tlm)
      END IF
      
NewtonLoopTLM:DO  NewtonIterTLM = 1, NewtonMaxit  
 
            !~~~> Prepare the right-hand side
            CALL RK_PrepareRHS_TLM(N,H,Jac1,Jac2,Jac3,Y_tlm(1,itlm), &
                      Z1_tlm(1,itlm),Z2_tlm(1,itlm),Z3_tlm(1,itlm),  &
                      DZ1,DZ2,DZ3) 
                                 
            !~~~> Solve the linear systems
            CALL RK_Solve( N,H,E1,IP1,E2,IP2,DZ1,DZ2,DZ3,ISING )
            
	    IF (TLMNewtonEst) THEN
!~~~>   Check convergence of Newton iterations
            NewtonIncrement = SQRT( ( RK_ErrorNorm(N,SCAL_tlm,DZ1)**2 + &
                                RK_ErrorNorm(N,SCAL_tlm,DZ2)**2 +       &
                                RK_ErrorNorm(N,SCAL_tlm,DZ3)**2 )/3.0d0 )
            IF ( NewtonIterTLM == 1 ) THEN
                ThetaTLM      = ABS(ThetaMin)
                NewtonRate = 2.0d0 
            ELSE
                ThetaTLM = NewtonIncrement/NewtonIncrementOld
                IF (ThetaTLM < 0.99d0) THEN
                    NewtonRate = ThetaTLM/(ONE-ThetaTLM)
                ELSE ! Non-convergence of Newton: ThetaTLM too large
                    EXIT NewtonLoopTLM
                END IF
                IF ( NewtonIterTLM < NewtonMaxit ) THEN 
                  ! Predict error at the end of Newton process 
                  NewtonPredictedErr = NewtonIncrement &
                               *ThetaTLM**(NewtonMaxit-NewtonIterTLM)/(ONE-ThetaTLM)
                  IF (NewtonPredictedErr >= NewtonTol) THEN
                    ! Non-convergence of Newton: predicted error too large
                    Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                    Fac=0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIterTLM))
                    EXIT NewtonLoopTLM
                  END IF
                END IF
            END IF

            NewtonIncrementOld = MAX(NewtonIncrement,Roundoff) 
	    END IF !(TLMNewtonEst)
            ! Update solution
            CALL WAXPY(N,-ONE,DZ1,1,Z1_tlm(1,itlm),1) ! Z1 <- Z1 - DZ1
            CALL WAXPY(N,-ONE,DZ2,1,Z2_tlm(1,itlm),1) ! Z2 <- Z2 - DZ2
            CALL WAXPY(N,-ONE,DZ3,1,Z3_tlm(1,itlm),1) ! Z3 <- Z3 - DZ3
            
            ! Check error in Newton iterations
            IF (TLMNewtonEst) THEN
  	       NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
	       IF (NewtonDone) EXIT NewtonLoopTLM
	    ELSE
  	       ! Minimum number of iterations same as FWD iterations
               IF (NewtonIterTLM>=saveNiter) EXIT NewtonLoopTLM
	    END IF
            
      END DO NewtonLoopTLM
            
      IF ((TLMNewtonEst) .AND. (.NOT.NewtonDone)) THEN
           !CALL RK_ErrorMsg(-12,T,H,IERR);
           H = Fac*H; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
           CYCLE Tloop
      END IF
              
  END DO Tlm    

  END IF TLMDIR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Error estimation, TLM solution
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (TLMtruncErr) THEN
        CALL RK_ErrorEstimate_tlm(N,NTLM,T,H,FJAC,Y, Y_tlm, &
               E1,IP1,Z1_tlm,Z2_tlm,Z3_tlm,Err,FirstStep,Reject)
      END IF

  END IF ! (Err<ONE)
            

!~~~> Computation of new step size Hnew
      Fac  = Err**(-ONE/rkELO)*   &
             MIN(FacSafe,(ONE+2*NewtonMaxit)/(NewtonIter+2*NewtonMaxit))
      Fac  = MIN(FacMax,MAX(FacMin,Fac))
      Hnew = Fac*H

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Accept/reject step 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
accept:IF (Err < ONE) THEN !~~~> STEP IS ACCEPTED

         FirstStep=.FALSE.
         ISTATUS(Nacc) = ISTATUS(Nacc) + 1
         IF (Gustafsson) THEN
            !~~~> Predictive controller of Gustafsson
            IF (ISTATUS(Nacc) > 1) THEN
               FacGus=FacSafe*(H/Hacc)*(Err**2/ErrOld)**(-0.25d0)
               FacGus=MIN(FacMax,MAX(FacMin,FacGus))
               Fac=MIN(Fac,FacGus)
               Hnew = Fac*H
            END IF
            Hacc=H
            ErrOld=MAX(1.0d-2,Err)
         END IF
         Hold = H
         T = T+H 
         ! Update solution: Y <- Y + sum(d_i Z_i)
         IF (rkD(1)/=ZERO) CALL WAXPY(N,rkD(1),Z1,1,Y,1)
         IF (rkD(2)/=ZERO) CALL WAXPY(N,rkD(2),Z2,1,Y,1)
         IF (rkD(3)/=ZERO) CALL WAXPY(N,rkD(3),Z3,1,Y,1)
         ! Update TLM solution: Y <- Y + sum(d_i*Z_i_tlm)
         DO itlm = 1,NTLM
            IF (rkD(1)/=ZERO) CALL WAXPY(N,rkD(1),Z1_tlm(1,itlm),1,Y_tlm(1,itlm),1)
            IF (rkD(2)/=ZERO) CALL WAXPY(N,rkD(2),Z2_tlm(1,itlm),1,Y_tlm(1,itlm),1)
            IF (rkD(3)/=ZERO) CALL WAXPY(N,rkD(3),Z3_tlm(1,itlm),1,Y_tlm(1,itlm),1)
         END DO   
         ! Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3
         IF (StartNewton) CALL RK_Interpolate('make',N,H,Hold,Z1,Z2,Z3,CONT)
         CALL RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)
         RSTATUS(Ntexit) = T
         RSTATUS(Nhnew)  = Hnew
         RSTATUS(Nhacc)  = H
         Hnew = Tdirection*MIN( MAX(ABS(Hnew),Hmin) , Hmax )
         IF (Reject) Hnew = Tdirection*MIN(ABS(Hnew),ABS(H))
         Reject = .FALSE.
         IF ((T+Hnew/Qmin-Tend)*Tdirection >=  ZERO) THEN
            H = Tend-T
         ELSE
            Hratio=Hnew/H
            ! Reuse the LU decomposition
            SkipLU = (Theta<=ThetaMin) .AND. (Hratio>=Qmin) .AND. (Hratio<=Qmax)
            ! For TLM: do not skip LU (decrease TLM error)
  	    SkipLU = .FALSE.
            IF (.NOT.SkipLU) H=Hnew
         END IF
         ! If convergence is fast enough, do not update Jacobian
!         SkipJac = (Theta <= ThetaMin)
         SkipJac  = .FALSE.	! For TLM: do not skip jac

      ELSE accept !~~~> Step is rejected

         IF (FirstStep .OR. Reject) THEN
             H = FacRej*H
         ELSE
             H = Hnew
         END IF
         Reject   = .TRUE.
         SkipJac  = .TRUE.
         SkipLU   = .FALSE. 
         IF (ISTATUS(Nacc) >= 1) ISTATUS(Nrej) = ISTATUS(Nrej) + 1

      END IF accept
      
    END DO Tloop
    
    ! Successful exit
    IERR = 1  

 END SUBROUTINE RK_IntegratorTLM



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE RK_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   KPP_REAL, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from RungeKutta due to the following error:'


   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)
      PRINT * , '--> Improper value for maximal no of Newton iterations'
    CASE (-3)
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)
      PRINT * , '--> Improper values for FacMin/FacMax/FacSafe/FacRej'
    CASE (-5)
      PRINT * , '--> Improper value for ThetaMin'
    CASE (-6)
      PRINT * , '--> Newton stopping tolerance too small'
    CASE (-7)
      PRINT * , '--> Improper values for Qmin, Qmax'
    CASE (-8)
      PRINT * , '--> Tolerances are too small'
    CASE (-9)
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-10)
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-11)
      PRINT * , '--> Matrix is repeatedly singular'
    CASE (-12)
      PRINT * , '--> Non-convergence of Newton iterations'
    CASE (-13)
      PRINT * , '--> Requested RK method not implemented'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   WRITE(6,FMT="(5X,'T=',E12.5,'  H=',E12.5)") T, H 

 END SUBROUTINE RK_ErrorMsg


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: N, ITOL
   KPP_REAL, INTENT(IN) :: AbsTol(*), RelTol(*), Y(N)
   KPP_REAL, INTENT(OUT) :: SCAL(N)
   INTEGER :: i
   
   IF (ITOL==0) THEN
       DO i=1,N
          SCAL(i)= ONE/(AbsTol(1)+RelTol(1)*ABS(Y(i)))
       END DO
   ELSE
       DO i=1,N
          SCAL(i)=ONE/(AbsTol(i)+RelTol(i)*ABS(Y(i)))
       END DO
   END IF
      
 END SUBROUTINE RK_ErrorScale


!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SUBROUTINE RK_Transform(N,Tr,Z1,Z2,Z3,W1,W2,W3)
!!~~~>                 W <-- Tr x Z
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      IMPLICIT NONE
!      INTEGER :: N, i
!      KPP_REAL :: Tr(3,3),Z1(N),Z2(N),Z3(N),W1(N),W2(N),W3(N)
!      KPP_REAL :: x1, x2, x3
!      DO i=1,N
!          x1 = Z1(i); x2 = Z2(i); x3 = Z3(i)
!          W1(i) = Tr(1,1)*x1 + Tr(1,2)*x2 + Tr(1,3)*x3
!          W2(i) = Tr(2,1)*x1 + Tr(2,2)*x2 + Tr(2,3)*x3
!          W3(i) = Tr(3,1)*x1 + Tr(3,2)*x2 + Tr(3,3)*x3
!      END DO
!  END SUBROUTINE RK_Transform
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE RK_Interpolate(action,N,H,Hold,Z1,Z2,Z3,CONT)
!~~~>   Constructs or evaluates a quadratic polynomial
!         that interpolates the Z solution at current step
!         and provides starting values for the next step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER, INTENT(IN) :: N
      KPP_REAL, INTENT(IN) :: H,Hold
      KPP_REAL, INTENT(INOUT) :: Z1(N),Z2(N),Z3(N),CONT(N,3)
      CHARACTER(LEN=4), INTENT(IN) :: action
      KPP_REAL :: r, x1, x2, x3, den
      INTEGER :: i
       
      SELECT CASE (action) 
      CASE ('make')
         ! Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3
         den = (rkC(3)-rkC(2))*(rkC(2)-rkC(1))*(rkC(1)-rkC(3))
         DO i=1,N
             CONT(i,1)=(-rkC(3)**2*rkC(2)*Z1(i)+Z3(i)*rkC(2)*rkC(1)**2 &
                        +rkC(2)**2*rkC(3)*Z1(i)-rkC(2)**2*rkC(1)*Z3(i) &
                        +rkC(3)**2*rkC(1)*Z2(i)-Z2(i)*rkC(3)*rkC(1)**2)&
                        /den-Z3(i)
             CONT(i,2)= -( rkC(1)**2*(Z3(i)-Z2(i)) + rkC(2)**2*(Z1(i)  &
                          -Z3(i)) +rkC(3)**2*(Z2(i)-Z1(i)) )/den
             CONT(i,3)= ( rkC(1)*(Z3(i)-Z2(i)) + rkC(2)*(Z1(i)-Z3(i))  &
                           +rkC(3)*(Z2(i)-Z1(i)) )/den
         END DO
      CASE ('eval')
          ! Evaluate quadratic polynomial
          r = H/Hold
         x1 = ONE + rkC(1)*r
         x2 = ONE + rkC(2)*r
         x3 = ONE + rkC(3)*r
         DO i=1,N
            Z1(i) = CONT(i,1)+x1*(CONT(i,2)+x1*CONT(i,3))
            Z2(i) = CONT(i,1)+x2*(CONT(i,2)+x2*CONT(i,3))
            Z3(i) = CONT(i,1)+x3*(CONT(i,2)+x3*CONT(i,3))
         END DO
       END SELECT   
  END SUBROUTINE RK_Interpolate


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_PrepareRHS(N,T,H,Y,Z1,Z2,Z3,R1,R2,R3)
!~~~> Prepare the right-hand side for Newton iterations
!     R = Z - hA x F
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N
      KPP_REAL, INTENT(IN) :: T, H
      KPP_REAL, INTENT(IN), DIMENSION(N) :: Y,Z1,Z2,Z3
      KPP_REAL, INTENT(INOUT), DIMENSION(N) :: R1,R2,R3
      KPP_REAL, DIMENSION(N) :: F, TMP
      
      CALL WCOPY(N,Z1,1,R1,1) ! R1 <- Z1
      CALL WCOPY(N,Z2,1,R2,1) ! R2 <- Z2
      CALL WCOPY(N,Z3,1,R3,1) ! R3 <- Z3

      CALL WADD(N,Y,Z1,TMP)              ! TMP <- Y + Z1
      CALL FUN_CHEM(T+rkC(1)*H,TMP,F)    ! F1 <- Fun(Y+Z1)         
      CALL WAXPY(N,-H*rkA(1,1),F,1,R1,1) ! R1 <- R1 - h*A_11*F1
      CALL WAXPY(N,-H*rkA(2,1),F,1,R2,1) ! R2 <- R2 - h*A_21*F1
      CALL WAXPY(N,-H*rkA(3,1),F,1,R3,1) ! R3 <- R3 - h*A_31*F1

      CALL WADD(N,Y,Z2,TMP)              ! TMP <- Y + Z2
      CALL FUN_CHEM(T+rkC(2)*H,TMP,F)    ! F2 <- Fun(Y+Z2)        
      CALL WAXPY(N,-H*rkA(1,2),F,1,R1,1) ! R1 <- R1 - h*A_12*F2
      CALL WAXPY(N,-H*rkA(2,2),F,1,R2,1) ! R2 <- R2 - h*A_22*F2
      CALL WAXPY(N,-H*rkA(3,2),F,1,R3,1) ! R3 <- R3 - h*A_32*F2

      CALL WADD(N,Y,Z3,TMP)              ! TMP <- Y + Z3
      CALL FUN_CHEM(T+rkC(3)*H,TMP,F)    ! F3 <- Fun(Y+Z3)     
      CALL WAXPY(N,-H*rkA(1,3),F,1,R1,1) ! R1 <- R1 - h*A_13*F3
      CALL WAXPY(N,-H*rkA(2,3),F,1,R2,1) ! R2 <- R2 - h*A_23*F3
      CALL WAXPY(N,-H*rkA(3,3),F,1,R3,1) ! R3 <- R3 - h*A_33*F3
            
  END SUBROUTINE RK_PrepareRHS
  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_PrepareRHS_TLM(N,H,Jac1,Jac2,Jac3,Y_tlm, &
                      Z1_tlm,Z2_tlm,Z3_tlm,R1,R2,R3)
!~~~> Prepare the right-hand side for Newton iterations
!     R = Z_tlm - hA x Jac*Z_tlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N
      KPP_REAL, INTENT(IN) :: H
      KPP_REAL, INTENT(IN), DIMENSION(N) :: Y_tlm,Z1_tlm,Z2_tlm,Z3_tlm
      KPP_REAL, INTENT(INOUT), DIMENSION(N) :: R1,R2,R3
#ifdef FULL_ALGEBRA      
      KPP_REAL, INTENT(IN), DIMENSION(NVAR,NVAR)  :: Jac1, Jac2, Jac3
#else      
      KPP_REAL, INTENT(IN), DIMENSION(LU_NONZERO) :: Jac1, Jac2, Jac3
#endif
      KPP_REAL, DIMENSION(N) :: F, TMP

      CALL WCOPY(N,Z1_tlm,1,R1,1) ! R1 <- Z1_tlm
      CALL WCOPY(N,Z2_tlm,1,R2,1) ! R2 <- Z2_tlm
      CALL WCOPY(N,Z3_tlm,1,R3,1) ! R3 <- Z3_tlm

      CALL WADD(N,Y_tlm,Z1_tlm,TMP)      ! TMP <- Y + Z1
#ifdef FULL_ALGEBRA  
      F = MATMUL(Jac1,TMP)    
#else      
      CALL Jac_SP_Vec ( Jac1, TMP, F )   ! F1 <- Jac(Y+Z1)*(Y_tlm+Z1_tlm)  
#endif      
      CALL WAXPY(N,-H*rkA(1,1),F,1,R1,1) ! R1 <- R1 - h*A_11*F1
      CALL WAXPY(N,-H*rkA(2,1),F,1,R2,1) ! R2 <- R2 - h*A_21*F1
      CALL WAXPY(N,-H*rkA(3,1),F,1,R3,1) ! R3 <- R3 - h*A_31*F1

      CALL WADD(N,Y_tlm,Z2_tlm,TMP)      ! TMP <- Y + Z2
#ifdef FULL_ALGEBRA      
      F = MATMUL(Jac2,TMP)    
#else      
      CALL Jac_SP_Vec ( Jac2, TMP, F )   ! F2 <- Jac(Y+Z2)*(Y_tlm+Z2_tlm)  
#endif      
      CALL WAXPY(N,-H*rkA(1,2),F,1,R1,1) ! R1 <- R1 - h*A_12*F2
      CALL WAXPY(N,-H*rkA(2,2),F,1,R2,1) ! R2 <- R2 - h*A_22*F2
      CALL WAXPY(N,-H*rkA(3,2),F,1,R3,1) ! R3 <- R3 - h*A_32*F2

      CALL WADD(N,Y_tlm,Z3_tlm,TMP)      ! TMP <- Y + Z3
#ifdef FULL_ALGEBRA      
      F = MATMUL(Jac3,TMP)    
#else      
      CALL Jac_SP_Vec ( Jac3, TMP, F )   ! F3 <- Jac(Y+Z3)*(Y_tlm+Z3_tlm)  
#endif      
      CALL WAXPY(N,-H*rkA(1,3),F,1,R1,1) ! R1 <- R1 - h*A_13*F3
      CALL WAXPY(N,-H*rkA(2,3),F,1,R2,1) ! R2 <- R2 - h*A_23*F3
      CALL WAXPY(N,-H*rkA(3,3),F,1,R3,1) ! R3 <- R3 - h*A_33*F3
            
  END SUBROUTINE RK_PrepareRHS_TLM
  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_PrepareRHS_TLMdirect(N,H,Jac1,Jac2,Jac3,Y_tlm,Zbig)
!~~~> Prepare the right-hand side for direct solution
!     Z =  hA x Jac*Y_tlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N
      KPP_REAL, INTENT(IN) :: H
      KPP_REAL, INTENT(IN), DIMENSION(N) :: Y_tlm
! Commented code for sparse big linear algebra:  
!      KPP_REAL, INTENT(OUT), DIMENSION(3,N) :: Zbig
      KPP_REAL, INTENT(OUT), DIMENSION(3*N) :: Zbig
#ifdef FULL_ALGEBRA      
      KPP_REAL, INTENT(IN), DIMENSION(NVAR,NVAR)  :: Jac1, Jac2, Jac3
#else      
      KPP_REAL, INTENT(IN), DIMENSION(LU_NONZERO) :: Jac1, Jac2, Jac3
#endif
      KPP_REAL, DIMENSION(N) :: F, TMP

#ifdef FULL_ALGEBRA  
      F = MATMUL(Jac1,Y_tlm)    
#else      
      CALL Jac_SP_Vec ( Jac1, Y_tlm, F )   ! F1 <- Jac(Y+Z1)*(Y_tlm)  
#endif      
! Commented code for sparse big linear algebra:  
!      Zbig(1,1:N) = H*rkA(1,1)*F(1:N)
!      Zbig(2,1:N) = H*rkA(2,1)*F(1:N)
!      Zbig(3,1:N) = H*rkA(3,1)*F(1:N)
      Zbig(1:N)       = H*rkA(1,1)*F(1:N)
      Zbig(N+1:2*N)   = H*rkA(2,1)*F(1:N)
      Zbig(2*N+1:3*N) = H*rkA(3,1)*F(1:N)

#ifdef FULL_ALGEBRA      
      F = MATMUL(Jac2,Y_tlm)    
#else      
      CALL Jac_SP_Vec ( Jac2, Y_tlm, F )   ! F2 <- Jac(Y+Z2)*(Y_tlm)  
#endif    
! Commented code for sparse big linear algebra:  
!      Zbig(1,1:N) = Zbig(1,1:N) + H*rkA(1,2)*F(1:N)
!      Zbig(2,1:N) = Zbig(2,1:N) + H*rkA(2,2)*F(1:N)
!      Zbig(3,1:N) = Zbig(3,1:N) + H*rkA(3,2)*F(1:N)
      Zbig(1:N)       = Zbig(1:N)       + H*rkA(1,2)*F(1:N)
      Zbig(N+1:2*N)   = Zbig(N+1:2*N)   + H*rkA(2,2)*F(1:N)
      Zbig(2*N+1:3*N) = Zbig(2*N+1:3*N) + H*rkA(3,2)*F(1:N)

#ifdef FULL_ALGEBRA      
      F = MATMUL(Jac3,Y_tlm)    
#else      
      CALL Jac_SP_Vec ( Jac3, Y_tlm, F )   ! F3 <- Jac(Y+Z3)*(Y_tlm)  
#endif      
! Commented code for sparse big linear algebra:  
!      Zbig(1,1:N) = Zbig(1,1:N) + H*rkA(1,3)*F(1:N)
!      Zbig(2,1:N) = Zbig(2,1:N) + H*rkA(2,3)*F(1:N)
!      Zbig(3,1:N) = Zbig(3,1:N) + H*rkA(3,3)*F(1:N)
      Zbig(1:N)       = Zbig(1:N)       + H*rkA(1,3)*F(1:N)
      Zbig(N+1:2*N)   = Zbig(N+1:2*N)   + H*rkA(2,3)*F(1:N)
      Zbig(2*N+1:3*N) = Zbig(2*N+1:3*N) + H*rkA(3,3)*F(1:N)
            
  END SUBROUTINE RK_PrepareRHS_TLMdirect

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_Decomp(N,H,FJAC,E1,IP1,E2,IP2,ISING)
   !~~~> Compute the matrices E1 and E2 and their decompositions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N
      KPP_REAL, INTENT(IN) :: H
#ifdef FULL_ALGEBRA      
      KPP_REAL, INTENT(IN) :: FJAC(NVAR,NVAR)
      KPP_REAL, INTENT(OUT) ::E1(NVAR,NVAR)
      COMPLEX(kind=dp), INTENT(OUT) :: E2(N,N)
#else      
      KPP_REAL, INTENT(IN) :: FJAC(LU_NONZERO)
      KPP_REAL, INTENT(OUT) :: E1(LU_NONZERO)
      COMPLEX(kind=dp),INTENT(OUT) :: E2(LU_NONZERO)
#endif      
      INTEGER, INTENT(OUT) :: IP1(N), IP2(N), ISING

      INTEGER :: i, j
      KPP_REAL    :: Alpha, Beta, Gamma
      
      Gamma = rkGamma/H
      Alpha = rkAlpha/H
      Beta  = rkBeta /H

#ifdef FULL_ALGEBRA      
      DO j=1,N
         DO  i=1,N
            E1(i,j)=-FJAC(i,j)
         END DO
         E1(j,j)=E1(j,j)+Gamma
      END DO
      CALL DGETRF(N,N,E1,N,IP1,ISING) 
#else      
      DO i=1,LU_NONZERO
         E1(i)=-FJAC(i)
      END DO
      DO i=1,NVAR
         j=LU_DIAG(i); E1(j)=E1(j)+Gamma
      END DO
      CALL KppDecomp(E1,ISING)
#endif      
      
      IF (ISING /= 0) THEN
         ISTATUS(Ndec) = ISTATUS(Ndec) + 1
         RETURN
      END IF
     
#ifdef FULL_ALGEBRA      
      DO j=1,N
        DO i=1,N
          E2(i,j) = DCMPLX( -FJAC(i,j), ZERO )
        END DO
        E2(j,j) = E2(j,j) + CMPLX( Alpha, Beta )
      END DO
      CALL ZGETRF(N,N,E2,N,IP2,ISING)     
#else  
      DO i=1,LU_NONZERO
         E2(i) = DCMPLX( -FJAC(i), ZERO )
      END DO
      DO i=1,NVAR
         j = LU_DIAG(i); E2(j)=E2(j) + CMPLX( Alpha, Beta )
      END DO
      CALL KppDecompCmplx(E2,ISING)    
#endif      
      ISTATUS(Ndec) = ISTATUS(Ndec) + 1
      
   END SUBROUTINE RK_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_Solve(N,H,E1,IP1,E2,IP2,R1,R2,R3,ISING)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N,IP1(NVAR),IP2(NVAR)
      INTEGER, INTENT(OUT) :: ISING
#ifdef FULL_ALGEBRA      
      KPP_REAL, INTENT(IN) :: E1(NVAR,NVAR)
      COMPLEX(kind=dp), INTENT(IN) :: E2(NVAR,NVAR)
#else      
      KPP_REAL, INTENT(IN) :: E1(LU_NONZERO)
      COMPLEX(kind=dp), INTENT(IN) :: E2(LU_NONZERO)
#endif      
      KPP_REAL, INTENT(IN) :: H
      KPP_REAL, INTENT(INOUT) :: R1(N),R2(N),R3(N)

      KPP_REAL    :: x1, x2, x3
      COMPLEX(kind=dp) :: BC(N)
      INTEGER :: i
!      
     ! Z <- h^{-1) T^{-1) A^{-1) x Z
      DO i=1,N
          x1 = R1(i)/H; x2 = R2(i)/H; x3 = R3(i)/H
          R1(i) = rkTinvAinv(1,1)*x1 + rkTinvAinv(1,2)*x2 + rkTinvAinv(1,3)*x3
          R2(i) = rkTinvAinv(2,1)*x1 + rkTinvAinv(2,2)*x2 + rkTinvAinv(2,3)*x3
          R3(i) = rkTinvAinv(3,1)*x1 + rkTinvAinv(3,2)*x2 + rkTinvAinv(3,3)*x3
      END DO

#ifdef FULL_ALGEBRA      
      CALL DGETRS ('N',N,1,E1,N,IP1,R1,N,ISING) 
#else      
      CALL KppSolve (E1,R1)
#endif      
!      
      DO i=1,N
        BC(i) = DCMPLX(R2(i),R3(i))
      END DO
#ifdef FULL_ALGEBRA      
      CALL ZGETRS ('N',N,1,E2,N,IP2,BC,N,ISING) 
#else      
      CALL KppSolveCmplx (E2,BC)
#endif      
      DO i=1,N
        R2(i) = DBLE( BC(i) )
        R3(i) = AIMAG( BC(i) )
      END DO

      ! Z <- T x Z
      DO i=1,N
          x1 = R1(i); x2 = R2(i); x3 = R3(i)
          R1(i) = rkT(1,1)*x1 + rkT(1,2)*x2 + rkT(1,3)*x3
          R2(i) = rkT(2,1)*x1 + rkT(2,2)*x2 + rkT(2,3)*x3
          R3(i) = rkT(3,1)*x1 + rkT(3,2)*x2 + rkT(3,3)*x3
      END DO

      ISTATUS(Nsol) = ISTATUS(Nsol) + 1

   END SUBROUTINE RK_Solve


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_ErrorEstimate(N,H,Y,T,      &
               E1,IP1,Z1,Z2,Z3,SCAL,Err,     &
               FirstStep,Reject)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N, IP1(N)
#ifdef FULL_ALGEBRA      
      KPP_REAL, INTENT(IN) :: E1(NVAR,NVAR)
      INTEGER :: ISING
#else      
      KPP_REAL, INTENT(IN) :: E1(LU_NONZERO)
#endif      
      KPP_REAL, INTENT(IN) :: T,H,SCAL(N),Z1(N),Z2(N),Z3(N),Y(N)
      LOGICAL,INTENT(IN) :: FirstStep,Reject
      KPP_REAL, INTENT(INOUT) :: Err

      KPP_REAL :: F1(N),F2(N),F0(N),TMP(N)
      INTEGER :: i
      KPP_REAL :: HEE1,HEE2,HEE3

      HEE1  = rkE(1)/H
      HEE2  = rkE(2)/H
      HEE3  = rkE(3)/H

      CALL FUN_CHEM(T,Y,F0)
      ISTATUS(Nfun) = ISTATUS(Nfun) + 1

      DO  i=1,N
         F2(i)  = HEE1*Z1(i)+HEE2*Z2(i)+HEE3*Z3(i)
         TMP(i) = rkE(0)*F0(i) + F2(i)
      END DO

#ifdef FULL_ALGEBRA      
      CALL DGETRS ('N',N,1,E1,N,IP1,TMP,N,ISING) 
      IF ((rkMethod==R1A).OR.(rkMethod==GAU).OR.(rkMethod==L3A)) THEN
            CALL DGETRS ('N',N,1,E1,N,IP1,TMP,N,ISING)
      END IF	    
      IF (rkMethod==GAU) THEN  
            CALL DGETRS ('N',N,1,E1,N,IP1,TMP,N,ISING)
      ENDIF	    
#else      
      CALL KppSolve (E1, TMP)
      IF ((rkMethod==R1A).OR.(rkMethod==GAU).OR.(rkMethod==L3A)) THEN
            CALL KppSolve (E1,TMP)
      END IF	    
      IF (rkMethod==GAU) THEN
            CALL KppSolve (E1,TMP)
      END IF	    
#endif      

      Err = RK_ErrorNorm(N,SCAL,TMP)
!
      IF (Err < ONE) RETURN
firej:IF (FirstStep.OR.Reject) THEN
          DO i=1,N
             TMP(i)=Y(i)+TMP(i)
          END DO
          CALL FUN_CHEM(T,TMP,F1)    
          ISTATUS(Nfun) = ISTATUS(Nfun) + 1
          DO i=1,N
             TMP(i)=F1(i)+F2(i)
          END DO

#ifdef FULL_ALGEBRA      
          CALL DGETRS ('N',N,1,E1,N,IP1,TMP,N,ISING) 
#else      
          CALL KppSolve (E1, TMP)
#endif      
          Err = RK_ErrorNorm(N,SCAL,TMP)
       END IF firej
 
   END SUBROUTINE RK_ErrorEstimate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_ErrorEstimate_tlm(N,NTLM,T,H,FJAC,Y,Y_tlm, &
               E1,IP1,Z1_tlm,Z2_tlm,Z3_tlm,FWD_Err,     &
               FirstStep,Reject)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N,NTLM,IP1(N)
#ifdef FULL_ALGEBRA      
      KPP_REAL, INTENT(IN) :: FJAC(N,N),  E1(N,N)
      KPP_REAL :: J0(N,N)
      INTEGER  :: ISING
#else      
      KPP_REAL, INTENT(IN) :: FJAC(LU_NONZERO), E1(LU_NONZERO)
      KPP_REAL :: J0(LU_NONZERO)
#endif      
      KPP_REAL, INTENT(IN) :: T,H, Z1_tlm(N,NTLM),Z2_tlm(N,NTLM),Z3_tlm(N,NTLM), &
      				Y_tlm(N,NTLM), Y(N)
      LOGICAL, INTENT(IN) :: FirstStep, Reject
      KPP_REAL, INTENT(INOUT) :: FWD_Err

      INTEGER :: itlm
      KPP_REAL :: HEE1,HEE2,HEE3, SCAL_tlm(N), Err, TMP(N), TMP2(N), JY_tlm(N)
	
	HEE1  = rkE(1)/H
	HEE2  = rkE(2)/H
	HEE3  = rkE(3)/H
	
      DO itlm=1,NTLM
        CALL RK_ErrorScale(N,ITOL,AbsTol_tlm(1,itlm),RelTol_tlm(1,itlm), &
              Y_tlm(1,itlm),SCAL_tlm)

 	CALL JAC_CHEM(T,Y,J0)
	ISTATUS(Njac) = ISTATUS(Njac) + 1
	CALL JAC_SP_Vec(J0,Y_tlm(1,itlm),JY_tlm)

	DO  i=1,N
          TMP2(i)  = HEE1*Z1_tlm(i,itlm)+HEE2*Z2_tlm(i,itlm)+HEE3*Z3_tlm(i,itlm)
          TMP(i) = rkE(0)*JY_tlm(i) + TMP2(i)
	END DO

#ifdef FULL_ALGEBRA      
	CALL DGETRS ('N',N,1,E1,N,IP1,TMP,N,ISING) 
!	IF ((ICNTRL(3)==3).OR.(ICNTRL(3)==4)) 
!               CALL DGETRS('N',N,1,E1,N,IP1,TMP,N,ISING)
!	IF (ICNTRL(3)==3) CALL DGETRS ('N',N,1,E1,N,IP1,TMP,N,ISING)
#else      
	CALL KppSolve (E1, TMP)
!	IF ((ICNTRL(3)==3).OR.(ICNTRL(3)==4)) THEN
!          CALL KppSolve (E1,TMP)
!	END IF  
!	IF (ICNTRL(3)==3) THEN
!          CALL KppSolve (E1,TMP)
!	END IF  
#endif      

        Err = RK_ErrorNorm(N,SCAL_tlm,TMP)
!
	FWD_Err = MAX(FWD_Err, Err)
      END DO
      
   END SUBROUTINE RK_ErrorEstimate_tlm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   KPP_REAL FUNCTION RK_ErrorNorm(N,SCAL,DY)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N
      KPP_REAL, INTENT(IN) :: SCAL(N),DY(N)
      INTEGER :: i

      RK_ErrorNorm = ZERO
      DO i=1,N
          RK_ErrorNorm = RK_ErrorNorm + (DY(i)*SCAL(i))**2
      END DO
      RK_ErrorNorm = MAX( SQRT(RK_ErrorNorm/N), 1.0d-10 )
 
   END FUNCTION RK_ErrorNorm
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Radau2A_Coefficients
!    The coefficients of the 3-stage Radau-2A method
!    (given to ~30 accurate digits)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
! The coefficients of the Radau2A method
      KPP_REAL :: b0

!      b0 = 1.0d0
      IF (SdirkError) THEN
        b0 = 0.2d-1
      ELSE
        b0 = 0.5d-1
      END IF

! The coefficients of the Radau2A method
      rkMethod = R2A

      rkA(1,1) =  1.968154772236604258683861429918299d-1
      rkA(1,2) = -6.55354258501983881085227825696087d-2
      rkA(1,3) =  2.377097434822015242040823210718965d-2
      rkA(2,1) =  3.944243147390872769974116714584975d-1
      rkA(2,2) =  2.920734116652284630205027458970589d-1
      rkA(2,3) = -4.154875212599793019818600988496743d-2
      rkA(3,1) =  3.764030627004672750500754423692808d-1
      rkA(3,2) =  5.124858261884216138388134465196080d-1
      rkA(3,3) =  1.111111111111111111111111111111111d-1

      rkB(1) = 3.764030627004672750500754423692808d-1
      rkB(2) = 5.124858261884216138388134465196080d-1
      rkB(3) = 1.111111111111111111111111111111111d-1

      rkC(1) = 1.550510257216821901802715925294109d-1
      rkC(2) = 6.449489742783178098197284074705891d-1
      rkC(3) = 1.0d0
      
      ! New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
      rkD(1) = 0.0d0
      rkD(2) = 0.0d0
      rkD(3) = 1.0d0

      ! Classical error estimator: 
      ! H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) = 1.0d0*b0
      rkE(1) = -10.04880939982741556246032950764708d0*b0
      rkE(2) = 1.382142733160748895793662840980412d0*b0
      rkE(3) = -.3333333333333333333333333333333333d0*b0

      ! Sdirk error estimator
      rkBgam(0) = b0
      rkBgam(1) = .3764030627004672750500754423692807d0-1.558078204724922382431975370686279d0*b0
      rkBgam(2) = .8914115380582557157653087040196118d0*b0+.5124858261884216138388134465196077d0
      rkBgam(3) = -.1637777184845662566367174924883037d0-.3333333333333333333333333333333333d0*b0
      rkBgam(4) = .2748888295956773677478286035994148d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -1.520677486405081647234271944611547d0-10.04880939982741556246032950764708d0*b0
      rkTheta(2) = 2.070455145596436382729929151810376d0+1.382142733160748895793662840980413d0*b0
      rkTheta(3) = -.3333333333333333333333333333333333d0*b0-.3744441479783868387391430179970741d0

      ! Local order of error estimator 
      IF (b0==0.0d0) THEN
        rkELO  = 6.0d0
      ELSE	
        rkELO  = 4.0d0
      END IF	

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 3.637834252744495732208418513577775d0
      rkAlpha = 2.681082873627752133895790743211112d0
      rkBeta  = 3.050430199247410569426377624787569d0

      rkT(1,1) =  9.443876248897524148749007950641664d-2
      rkT(1,2) = -1.412552950209542084279903838077973d-1
      rkT(1,3) = -3.00291941051474244918611170890539d-2
      rkT(2,1) =  2.502131229653333113765090675125018d-1
      rkT(2,2) =  2.041293522937999319959908102983381d-1
      rkT(2,3) =  3.829421127572619377954382335998733d-1
      rkT(3,1) =  1.0d0
      rkT(3,2) =  1.0d0
      rkT(3,3) =  0.0d0

      rkTinv(1,1) =  4.178718591551904727346462658512057d0
      rkTinv(1,2) =  3.27682820761062387082533272429617d-1
      rkTinv(1,3) =  5.233764454994495480399309159089876d-1
      rkTinv(2,1) = -4.178718591551904727346462658512057d0
      rkTinv(2,2) = -3.27682820761062387082533272429617d-1
      rkTinv(2,3) =  4.766235545005504519600690840910124d-1
      rkTinv(3,1) = -5.02872634945786875951247343139544d-1
      rkTinv(3,2) =  2.571926949855605429186785353601676d0
      rkTinv(3,3) = -5.960392048282249249688219110993024d-1

      rkTinvAinv(1,1) =  1.520148562492775501049204957366528d+1
      rkTinvAinv(1,2) =  1.192055789400527921212348994770778d0
      rkTinvAinv(1,3) =  1.903956760517560343018332287285119d0
      rkTinvAinv(2,1) = -9.669512977505946748632625374449567d0
      rkTinvAinv(2,2) = -8.724028436822336183071773193986487d0
      rkTinvAinv(2,3) =  3.096043239482439656981667712714881d0
      rkTinvAinv(3,1) = -1.409513259499574544876303981551774d+1
      rkTinvAinv(3,2) =  5.895975725255405108079130152868952d0
      rkTinvAinv(3,3) = -1.441236197545344702389881889085515d-1

      rkAinvT(1,1) = .3435525649691961614912493915818282d0
      rkAinvT(1,2) = -.4703191128473198422370558694426832d0
      rkAinvT(1,3) = .3503786597113668965366406634269080d0
      rkAinvT(2,1) = .9102338692094599309122768354288852d0
      rkAinvT(2,2) = 1.715425895757991796035292755937326d0
      rkAinvT(2,3) = .4040171993145015239277111187301784d0
      rkAinvT(3,1) = 3.637834252744495732208418513577775d0
      rkAinvT(3,2) = 2.681082873627752133895790743211112d0
      rkAinvT(3,3) = -3.050430199247410569426377624787569d0

  END SUBROUTINE Radau2A_Coefficients

    

    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Lobatto3C_Coefficients
!    The coefficients of the 3-stage Lobatto-3C method
!    (given to ~30 accurate digits)
!    The parameter b0 can be chosen to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
      KPP_REAL :: b0

      rkMethod = L3C

!      b0 = 1.0d0
      IF (SdirkError) THEN
        b0 = 0.2d0
      ELSE
        b0 = 0.5d0
      END IF
! The coefficients of the Lobatto3C method

      rkA(1,1) =  .1666666666666666666666666666666667d0
      rkA(1,2) = -.3333333333333333333333333333333333d0
      rkA(1,3) =  .1666666666666666666666666666666667d0
      rkA(2,1) =  .1666666666666666666666666666666667d0
      rkA(2,2) =  .4166666666666666666666666666666667d0
      rkA(2,3) = -.8333333333333333333333333333333333d-1
      rkA(3,1) =  .1666666666666666666666666666666667d0
      rkA(3,2) =  .6666666666666666666666666666666667d0
      rkA(3,3) =  .1666666666666666666666666666666667d0

      rkB(1) = .1666666666666666666666666666666667d0
      rkB(2) = .6666666666666666666666666666666667d0
      rkB(3) = .1666666666666666666666666666666667d0

      rkC(1) = 0.0d0
      rkC(2) = 0.5d0
      rkC(3) = 1.0d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) = .16666666666666666666666666666666667d0-b0
      rkBhat(2) = .66666666666666666666666666666666667d0
      rkBhat(3) = .16666666666666666666666666666666667d0
      
      ! New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
      rkD(1) = 0.0d0
      rkD(2) = 0.0d0
      rkD(3) = 1.0d0

      ! Classical error estimator: 
      !   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) =   .3808338772072650364017425226487022*b0
      rkE(1) = -1.142501631621795109205227567946107*b0
      rkE(2) = -1.523335508829060145606970090594809*b0
      rkE(3) =   .3808338772072650364017425226487022*b0

      ! Sdirk error estimator
      rkBgam(0) = b0
      rkBgam(1) = .1666666666666666666666666666666667d0-1.d0*b0
      rkBgam(2) = .6666666666666666666666666666666667d0
      rkBgam(3) = -.2141672105405983697350758559820354d0
      rkBgam(4) = .3808338772072650364017425226487021d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -3.d0*b0-.3808338772072650364017425226487021d0
      rkTheta(2) = -4.d0*b0+1.523335508829060145606970090594808d0
      rkTheta(3) = -.142501631621795109205227567946106d0+b0

      ! Local order of error estimator 
      IF (b0==0.0d0) THEN
        rkELO  = 5.0d0
      ELSE	
        rkELO  = 4.0d0
      END IF	

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 2.625816818958466716011888933765284d0
      rkAlpha = 1.687091590520766641994055533117359d0
      rkBeta  = 2.508731754924880510838743672432351d0

      rkT(1,1) = 1.d0
      rkT(1,2) = 1.d0
      rkT(1,3) = 0.d0
      rkT(2,1) = .4554100411010284672111720348287483d0
      rkT(2,2) = -.6027050205505142336055860174143743d0
      rkT(2,3) = -.4309321229203225731070721341350346d0
      rkT(3,1) = 2.195823345445647152832799205549709d0
      rkT(3,2) = -1.097911672722823576416399602774855d0
      rkT(3,3) = .7850032632435902184104551358922130d0

      rkTinv(1,1) = .4205559181381766909344950150991349d0
      rkTinv(1,2) = .3488903392193734304046467270632057d0
      rkTinv(1,3) = .1915253879645878102698098373933487d0
      rkTinv(2,1) = .5794440818618233090655049849008650d0
      rkTinv(2,2) = -.3488903392193734304046467270632057d0
      rkTinv(2,3) = -.1915253879645878102698098373933487d0
      rkTinv(3,1) = -.3659705575742745254721332009249516d0
      rkTinv(3,2) = -1.463882230297098101888532803699806d0
      rkTinv(3,3) = .4702733607340189781407813565524989d0

      rkTinvAinv(1,1) = 1.104302803159744452668648155627548d0
      rkTinvAinv(1,2) = .916122120694355522658740710823143d0
      rkTinvAinv(1,3) = .5029105849749601702795812241441172d0
      rkTinvAinv(2,1) = 1.895697196840255547331351844372453d0
      rkTinvAinv(2,2) = 3.083877879305644477341259289176857d0
      rkTinvAinv(2,3) = -1.502910584974960170279581224144117d0
      rkTinvAinv(3,1) = .8362439183082935036129145574774502d0
      rkTinvAinv(3,2) = -3.344975673233174014451658229909802d0
      rkTinvAinv(3,3) = .312908409479233358005944466882642d0

      rkAinvT(1,1) = 2.625816818958466716011888933765282d0
      rkAinvT(1,2) = 1.687091590520766641994055533117358d0
      rkAinvT(1,3) = -2.508731754924880510838743672432351d0
      rkAinvT(2,1) = 1.195823345445647152832799205549710d0
      rkAinvT(2,2) = -2.097911672722823576416399602774855d0
      rkAinvT(2,3) = .7850032632435902184104551358922130d0
      rkAinvT(3,1) = 5.765829871932827589653709477334136d0
      rkAinvT(3,2) = .1170850640335862051731452613329320d0
      rkAinvT(3,3) = 4.078738281412060947659653944216779d0

  END SUBROUTINE Lobatto3C_Coefficients

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Gauss_Coefficients
!    The coefficients of the 3-stage Gauss method
!    (given to ~30 accurate digits)
!    The parameter b3 can be chosen by the user
!    to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
      KPP_REAL :: b0
! The coefficients of the Gauss method


      rkMethod = GAU
      
!      b0 = 4.0d0
      b0 = 0.1d0
      
! The coefficients of the Gauss method

      rkA(1,1) =  .1388888888888888888888888888888889d0
      rkA(1,2) = -.359766675249389034563954710966045d-1
      rkA(1,3) =  .97894440153083260495800422294756d-2
      rkA(2,1) =  .3002631949808645924380249472131556d0
      rkA(2,2) =  .2222222222222222222222222222222222d0
      rkA(2,3) = -.224854172030868146602471694353778d-1
      rkA(3,1) =  .2679883337624694517281977355483022d0
      rkA(3,2) =  .4804211119693833479008399155410489d0
      rkA(3,3) =  .1388888888888888888888888888888889d0

      rkB(1) = .2777777777777777777777777777777778d0
      rkB(2) = .4444444444444444444444444444444444d0
      rkB(3) = .2777777777777777777777777777777778d0

      rkC(1) = .1127016653792583114820734600217600d0
      rkC(2) = .5000000000000000000000000000000000d0
      rkC(3) = .8872983346207416885179265399782400d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) =-1.4788305577012361475298775666303999d0*b0 &
                  +.27777777777777777777777777777777778d0
      rkBhat(2) =  .44444444444444444444444444444444444d0 &
                  +.66666666666666666666666666666666667d0*b0
      rkBhat(3) = -.18783610896543051913678910003626672d0*b0 &
                  +.27777777777777777777777777777777778d0

      ! New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
      rkD(1) = .1666666666666666666666666666666667d1
      rkD(2) = -.1333333333333333333333333333333333d1
      rkD(3) = .1666666666666666666666666666666667d1

      ! Classical error estimator: 
      !   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) = .2153144231161121782447335303806954d0*b0
      rkE(1) = -2.825278112319014084275808340593191d0*b0
      rkE(2) = .2870858974881495709929780405075939d0*b0
      rkE(3) = -.4558086256248162565397206448274867d-1*b0

      ! Sdirk error estimator
      rkBgam(0) = 0.d0
      rkBgam(1) = .2373339543355109188382583162660537d0
      rkBgam(2) = .5879873931885192299409334646982414d0
      rkBgam(3) = -.4063577064014232702392531134499046d-1
      rkBgam(4) = .2153144231161121782447335303806955d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -2.594040933093095272574031876464493d0
      rkTheta(2) = 1.824611539036311947589425112250199d0
      rkTheta(3) = .1856563166634371860478043996459493d0

      ! ELO = local order of classical error estimator 
      rkELO = 4.0d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 4.644370709252171185822941421408064d0
      rkAlpha = 3.677814645373914407088529289295970d0
      rkBeta  = 3.508761919567443321903661209182446d0

      rkT(1,1) =  .7215185205520017032081769924397664d-1
      rkT(1,2) = -.8224123057363067064866206597516454d-1
      rkT(1,3) = -.6012073861930850173085948921439054d-1
      rkT(2,1) =  .1188325787412778070708888193730294d0
      rkT(2,2) =  .5306509074206139504614411373957448d-1
      rkT(2,3) =  .3162050511322915732224862926182701d0
      rkT(3,1) = 1.d0
      rkT(3,2) = 1.d0
      rkT(3,3) = 0.d0

      rkTinv(1,1) =  5.991698084937800775649580743981285d0
      rkTinv(1,2) =  1.139214295155735444567002236934009d0
      rkTinv(1,3) =   .4323121137838583855696375901180497d0
      rkTinv(2,1) = -5.991698084937800775649580743981285d0
      rkTinv(2,2) = -1.139214295155735444567002236934009d0
      rkTinv(2,3) =   .5676878862161416144303624098819503d0
      rkTinv(3,1) = -1.246213273586231410815571640493082d0
      rkTinv(3,2) =  2.925559646192313662599230367054972d0
      rkTinv(3,3) =  -.2577352012734324923468722836888244d0

      rkTinvAinv(1,1) =  27.82766708436744962047620566703329d0
      rkTinvAinv(1,2) =   5.290933503982655311815946575100597d0
      rkTinvAinv(1,3) =   2.007817718512643701322151051660114d0
      rkTinvAinv(2,1) = -17.66368928942422710690385180065675d0
      rkTinvAinv(2,2) = -14.45491129892587782538830044147713d0
      rkTinvAinv(2,3) =   2.992182281487356298677848948339886d0
      rkTinvAinv(3,1) = -25.60678350282974256072419392007303d0
      rkTinvAinv(3,2) =   6.762434375611708328910623303779923d0
      rkTinvAinv(3,3) =   1.043979339483109825041215970036771d0
      
      rkAinvT(1,1) = .3350999483034677402618981153470483d0
      rkAinvT(1,2) = -.5134173605009692329246186488441294d0
      rkAinvT(1,3) = .6745196507033116204327635673208923d-1
      rkAinvT(2,1) = .5519025480108928886873752035738885d0
      rkAinvT(2,2) = 1.304651810077110066076640761092008d0
      rkAinvT(2,3) = .9767507983414134987545585703726984d0
      rkAinvT(3,1) = 4.644370709252171185822941421408064d0
      rkAinvT(3,2) = 3.677814645373914407088529289295970d0
      rkAinvT(3,3) = -3.508761919567443321903661209182446d0
      
  END SUBROUTINE Gauss_Coefficients



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Radau1A_Coefficients
!    The coefficients of the 3-stage Gauss method
!    (given to ~30 accurate digits)
!    The parameter b3 can be chosen by the user
!    to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
!      KPP_REAL :: b0 = 0.3d0
      KPP_REAL :: b0 = 0.1d0

! The coefficients of the Radau1A method

      rkMethod = R1A

      rkA(1,1) =  .1111111111111111111111111111111111d0
      rkA(1,2) = -.1916383190435098943442935597058829d0
      rkA(1,3) =  .8052720793239878323318244859477174d-1
      rkA(2,1) =  .1111111111111111111111111111111111d0
      rkA(2,2) =  .2920734116652284630205027458970589d0
      rkA(2,3) = -.481334970546573839513422644787591d-1
      rkA(3,1) =  .1111111111111111111111111111111111d0
      rkA(3,2) =  .5370223859435462728402311533676479d0
      rkA(3,3) =  .1968154772236604258683861429918299d0

      rkB(1) = .1111111111111111111111111111111111d0
      rkB(2) = .5124858261884216138388134465196080d0
      rkB(3) = .3764030627004672750500754423692808d0

      rkC(1) = 0.d0
      rkC(2) = .3550510257216821901802715925294109d0
      rkC(3) = .8449489742783178098197284074705891d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) = .11111111111111111111111111111111111d0-b0
      rkBhat(2) = .51248582618842161383881344651960810d0
      rkBhat(3) = .37640306270046727505007544236928079d0

      ! New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
      rkD(1) = .3333333333333333333333333333333333d0
      rkD(2) = -.8914115380582557157653087040196127d0
      rkD(3) = .1558078204724922382431975370686279d1

      ! Classical error estimator: 
      ! H* Sum (b_j-bhat_j) f(Z_j) = H*E(0)*F(0) + Sum E_j Z_j
      rkE(0) =   .2748888295956773677478286035994148d0*b0
      rkE(1) = -1.374444147978386838739143017997074d0*b0
      rkE(2) = -1.335337922441686804550326197041126d0*b0
      rkE(3) =   .235782604058977333559011782643466d0*b0

      ! Sdirk error estimator
      rkBgam(0) = 0.0d0
      rkBgam(1) = .1948150124588532186183490991130616d-1
      rkBgam(2) = .7575249005733381398986810981093584d0
      rkBgam(3) = -.518952314149008295083446116200793d-1
      rkBgam(4) = .2748888295956773677478286035994148d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -1.224370034375505083904362087063351d0
      rkTheta(2) = .9340045331532641409047527962010133d0
      rkTheta(3) = .4656990124352088397561234800640929d0

      ! ELO = local order of classical error estimator 
      rkELO = 4.0d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 3.637834252744495732208418513577775d0
      rkAlpha = 2.681082873627752133895790743211112d0
      rkBeta  = 3.050430199247410569426377624787569d0

      rkT(1,1) =  .424293819848497965354371036408369d0
      rkT(1,2) = -.3235571519651980681202894497035503d0
      rkT(1,3) = -.522137786846287839586599927945048d0
      rkT(2,1) =  .57594609499806128896291585429339d-1
      rkT(2,2) =  .3148663231849760131614374283783d-2
      rkT(2,3) =  .452429247674359778577728510381731d0
      rkT(3,1) = 1.d0
      rkT(3,2) = 1.d0
      rkT(3,3) = 0.d0

      rkTinv(1,1) = 1.233523612685027760114769983066164d0
      rkTinv(1,2) = 1.423580134265707095505388133369554d0
      rkTinv(1,3) = .3946330125758354736049045150429624d0
      rkTinv(2,1) = -1.233523612685027760114769983066164d0
      rkTinv(2,2) = -1.423580134265707095505388133369554d0
      rkTinv(2,3) = .6053669874241645263950954849570376d0
      rkTinv(3,1) = -.1484438963257383124456490049673414d0
      rkTinv(3,2) = 2.038974794939896109682070471785315d0
      rkTinv(3,3) = -.544501292892686735299355831692542d-1

      rkTinvAinv(1,1) =  4.487354449794728738538663081025420d0
      rkTinvAinv(1,2) =  5.178748573958397475446442544234494d0
      rkTinvAinv(1,3) =  1.435609490412123627047824222335563d0
      rkTinvAinv(2,1) = -2.854361287939276673073807031221493d0
      rkTinvAinv(2,2) = -1.003648660720543859000994063139137d+1
      rkTinvAinv(2,3) =  1.789135380979465422050817815017383d0
      rkTinvAinv(3,1) = -4.160768067752685525282947313530352d0
      rkTinvAinv(3,2) =  1.124128569859216916690209918405860d0
      rkTinvAinv(3,3) =  1.700644430961823796581896350418417d0

      rkAinvT(1,1) = 1.543510591072668287198054583233180d0
      rkAinvT(1,2) = -2.460228411937788329157493833295004d0
      rkAinvT(1,3) = -.412906170450356277003910443520499d0
      rkAinvT(2,1) = .209519643211838264029272585946993d0
      rkAinvT(2,2) = 1.388545667194387164417459732995766d0
      rkAinvT(2,3) = 1.20339553005832004974976023130002d0
      rkAinvT(3,1) = 3.637834252744495732208418513577775d0
      rkAinvT(3,2) = 2.681082873627752133895790743211112d0
      rkAinvT(3,3) = -3.050430199247410569426377624787569d0

  END SUBROUTINE Radau1A_Coefficients

  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE RungeKuttaTLM ! and all its internal procedures
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE FUN_CHEM(T, V, FCT)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    USE KPP_ROOT_Parameters
    USE KPP_ROOT_Global
    USE KPP_ROOT_Function, ONLY: Fun
    USE KPP_ROOT_Rates, ONLY: Update_SUN, Update_RCONST, Update_PHOTO

    IMPLICIT NONE

    KPP_REAL, INTENT(IN) :: V(NVAR), T
    KPP_REAL, INTENT(INOUT) :: FCT(NVAR)
    KPP_REAL :: Told

    Told = TIME
    TIME = T
    CALL Update_SUN()
    CALL Update_RCONST()
    CALL Update_PHOTO()
    TIME = Told
    
    CALL Fun(V, FIX, RCONST, FCT)
    
  END SUBROUTINE FUN_CHEM


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE JAC_CHEM (T, V, JF)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    USE KPP_ROOT_Parameters
    USE KPP_ROOT_Global
    USE KPP_ROOT_JacobianSP
    USE KPP_ROOT_Jacobian, ONLY: Jac_SP
    USE KPP_ROOT_Rates, ONLY: Update_SUN, Update_RCONST, Update_PHOTO

    IMPLICIT NONE

    KPP_REAL, INTENT(IN) :: V(NVAR), T 
#ifdef FULL_ALGEBRA    
    KPP_REAL, INTENT(INOUT) :: JF(NVAR,NVAR)
    KPP_REAL :: JV(LU_NONZERO)
    INTEGER       :: i, j 
#else
    KPP_REAL, INTENT(INOUT) :: JF(LU_NONZERO)
#endif   

    KPP_REAL :: Told

    Told = TIME
    TIME = T
    CALL Update_SUN()
    CALL Update_RCONST()
    CALL Update_PHOTO()
    TIME = Told
    
#ifdef FULL_ALGEBRA    
    CALL Jac_SP(V, FIX, RCONST, JV)
    DO j=1,NVAR
      DO i=1,NVAR
         JF(i,j) = 0.0d0
      END DO
    END DO
    DO i=1,LU_NONZERO
       JF(LU_IROW(i),LU_ICOL(i)) = JV(i)
    END DO
#else
    CALL Jac_SP(V, FIX, RCONST, JF) 
#endif   

  END SUBROUTINE JAC_CHEM

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
END MODULE KPP_ROOT_Integrator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

