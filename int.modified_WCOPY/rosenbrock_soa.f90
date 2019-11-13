MODULE KPP_ROOT_Integrator
 
 USE KPP_ROOT_Precision
 USE KPP_ROOT_Parameters
 USE KPP_ROOT_Global
 USE KPP_ROOT_Function
 USE KPP_ROOT_Jacobian
 USE KPP_ROOT_Hessian
 USE KPP_ROOT_LinearAlgebra
 USE KPP_ROOT_Rates
 
   IMPLICIT NONE
   PUBLIC
   SAVE
!~~~>  Statistics on the work performed by the Rosenbrock method
   INTEGER :: Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng
   INTEGER, PARAMETER :: ifun=11,  ijac=12, istp=13,  &
                iacc=14, irej=15,  idec=16, isol=17,  &
                isng=18, itexit=11,ihexit=12
!~~~>  Checkpoints in memory
   INTEGER, PARAMETER :: bufsize = 1500
   INTEGER :: stack_ptr = 0 ! last written entry
   KPP_REAL, DIMENSION(:), POINTER :: buf_H, buf_T
   KPP_REAL, DIMENSION(:,:), POINTER :: buf_Y, buf_K
   KPP_REAL, DIMENSION(:,:), POINTER :: buf_Y_tlm, buf_K_tlm

CONTAINS

 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE INTEGRATE_SOA( NSOA, Y, Y_tlm, Y_adj, Y_soa, TIN, TOUT, &
       AtolAdj, RtolAdj, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   IMPLICIT NONE        
    
!~~~> NSOA - No. of vectors to multiply SOA with
   INTEGER :: NSOA
!~~~> Y -   Forward model variables
   KPP_REAL, INTENT(INOUT)  :: Y(NVAR)
!~~~> Y_adj -   Tangent linear variables
   KPP_REAL, INTENT(INOUT)  :: Y_tlm(NVAR,NSOA)
!~~~> Y_adj   - First order adjoint
   KPP_REAL, INTENT(INOUT)  :: Y_adj(NVAR)
!~~~> Y_soa - Second order adjoint
   KPP_REAL, INTENT(INOUT)  :: Y_soa(NVAR,NSOA)
!~~~> 
   KPP_REAL, INTENT(IN)         :: TIN  ! TIN  - Start Time
   KPP_REAL, INTENT(IN)         :: TOUT ! TOUT - End Time
!~~~> Optional input parameters and statistics
   INTEGER,  INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   KPP_REAL, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   KPP_REAL, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)

   INTEGER N_stp, N_acc, N_rej, N_sng, IERR
   SAVE N_stp, N_acc, N_rej, N_sng
   INTEGER i
   KPP_REAL :: RCNTRL(20), RSTATUS(20)
   KPP_REAL :: AtolAdj(NVAR), RtolAdj(NVAR)
   INTEGER :: ICNTRL(20), ISTATUS(20)

   ICNTRL(1:20)  = 0
   RCNTRL(1:20)  = 0.0_dp
   ISTATUS(1:20) = 0
   RSTATUS(1:20) = 0.0_dp
   
   ICNTRL(1) = 0       ! 0 = non-autonomous, 1 = autonomous
   ICNTRL(2) = 1       ! 0 = scalar, 1 = vector tolerances
   RCNTRL(3) = STEPMIN ! starting step
   ICNTRL(4) = 5       ! choice of the method for forward and adjoint integration

! Tighter tolerances, especially atol, are needed for the full continuous adjoint
!  (Atol on sensitivities is different than on concentrations)
!   CADJ_ATOL(1:NVAR) = 1.0d-5
!   CADJ_RTOL(1:NVAR) = 1.0d-4

   ! if optional parameters are given, and if they are >=0, 
   !       then they overwrite default settings
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) >= 0) ICNTRL(:) = ICNTRL_U(:)
   ENDIF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) >= 0) RCNTRL(:) = RCNTRL_U(:)
   ENDIF

   
   CALL RosenbrockSOA(NSOA,                   &
         Y, Y_tlm, Y_adj, Y_soa,                &
         TIN,TOUT,                              &
         ATOL,RTOL,                             &
         Fun_Template,Jac_Template,Hess_Template,  &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)

             
!   N_stp = N_stp + ICNTRL(istp)
!   N_acc = N_acc + ICNTRL(iacc)
!   N_rej = N_rej + ICNTRL(irej)
!   N_sng = N_sng + ICNTRL(isng)
!   PRINT*,'Step=',N_stp,' Acc=',N_acc,' Rej=',N_rej, &
!        ' Singular=',N_sng

   IF (IERR < 0) THEN
     print *,'RosenbrockSOA: Unsucessful step at T=', &
         TIN,' (IERR=',IERR,')'
   ENDIF

   STEPMIN = RCNTRL(ihexit)
   ! if optional parameters are given for output they to return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)

END SUBROUTINE INTEGRATE_SOA



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE RosenbrockSOA( NSOA,            &
           Y, Y_tlm, Y_adj, Y_soa,         &
           Tstart,Tend,                    &
           AbsTol,RelTol,                  &
           ode_Fun,ode_Jac , ode_Hess,     &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   
!    ADJ = Adjoint of the Tangent Linear Model of a RosenbrockSOA Method
!
!    Solves the system y'=F(t,y) using a RosenbrockSOA method defined by:
!
!     G = 1/(H*gamma(1)) - ode_Jac(t0,Y0)
!     T_i = t0 + Alpha(i)*H
!     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
!     G * K_i = ode_Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
!         gamma(i)*dF/dT(t0, Y0)
!     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j 
!
!    For details on RosenbrockSOA methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.  
!    The codes contained in the book inspired this implementation.       
!
!    (C)  Adrian Sandu, August 2004
!    Virginia Polytechnic Institute and State University    
!    Contact: sandu@cs.vt.edu
!    This implementation is part of KPP - the Kinetic PreProcessor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    
!~~~>   INPUT ARGUMENTS: 
!    
!-     Y(NVAR)     -> vector of initial conditions (at T=Tstart)
!      NSOA        -> dimension of linearized system, 
!                     i.e. the number of sensitivity coefficients
!-     Y_adj(NVAR) -> vector of initial sensitivity conditions (at T=Tstart)
!-    [Tstart,Tend]  = time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)  
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE ode_Fun( T, Y, Ydot ) = ODE function, 
!                       returns Ydot = Y' = F(T,Y) 
!- SUBROUTINE ode_Fun( T, Y, Ydot ) = Jacobian of the ODE function,
!                       returns Jcb = dF/dY 
!-    ICNTRL(1:10)    = integer inputs parameters
!-    RCNTRL(1:10)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:  
!     
!-    Y(NVAR)         -> vector of final states (at T->Tend)
!-    Y_adj(NVAR)     -> vector of final sensitivities (at T=Tend)
!-    ICNTRL(11:20)   -> integer output parameters
!-    RCNTRL(11:20)   -> real output parameters
!-    IERR       -> job status upon return
!       - succes (positive value) or failure (negative value) -
!           =  1 : Success
!           = -1 : Improper value for maximal no of steps
!           = -2 : Selected RosenbrockSOA method not implemented
!           = -3 : Hmin/Hmax/Hstart must be positive
!           = -4 : FacMin/FacMax/FacRej must be positive
!           = -5 : Improper tolerance values
!           = -6 : No of steps exceeds maximum bound
!           = -7 : Step size too small
!           = -8 : Matrix is repeatedly singular
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1)   = 1: F = F(y)   Independent of T (AUTONOMOUS)
!              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)
!    ICNTRL(2)   = 0: AbsTol, RelTol are NVAR-dimensional vectors
!              = 1:  AbsTol, RelTol are scalars
!    ICNTRL(3)  -> maximum number of integration steps
!        For ICNTRL(3)=0) the default value of 100000 is used
!
!    ICNTRL(4)  -> selection of a particular Rosenbrock method
!        = 0 :  default method is Rodas3
!        = 1 :  method is  Ros2
!        = 2 :  method is  Ros3 
!        = 3 :  method is  Ros4 
!        = 4 :  method is  Rodas3
!        = 5:   method is  Rodas4
!
!    ICNTRL(5) -> Type of adjoint algorithm
!         = 0 : default is discrete adjoint ( of method ICNTRL(4) )
!         = 1 : no adjoint       
!         = 2 : discrete adjoint ( of method ICNTRL(4) )
!         = 3 : fully adaptive continuous adjoint ( with method ICNTRL(6) )
!         = 4 : simplified continuous adjoint ( with method ICNTRL(6) )
!
!    ICNTRL(6)  -> selection of a particular Rosenbrock method for the
!                continuous adjoint integration - for cts adjoint it
!                can be different than the forward method ICNTRL(4)
!         Note 1: to avoid interpolation errors (which can be huge!) 
!                   it is recommended to use only ICNTRL(6) = 1 or 4
!         Note 2: the performance of the full continuous adjoint
!                   strongly depends on the forward solution accuracy Abs/RelTol
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO 
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!          
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMin,upper bound on step increase factor (default=6)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!            (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller 
!         than the predicted value  (default=0.9)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call to RosenbrockSOA adds the corrent no. of fcn calls
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
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    For multiple restarts, use Hexit as Hstart in the following run 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
   
!~~~>  Arguments   
   INTEGER, INTENT(IN)     :: NSOA
   KPP_REAL, INTENT(INOUT) :: Y(NVAR)
   KPP_REAL, INTENT(INOUT) :: Y_tlm(NVAR,NSOA)
   KPP_REAL, INTENT(INOUT) :: Y_adj(NVAR)
   KPP_REAL, INTENT(INOUT) :: Y_soa(NVAR,NSOA)
   KPP_REAL, INTENT(IN)    :: Tstart, Tend
   KPP_REAL, INTENT(IN)    :: AbsTol(NVAR),RelTol(NVAR)
   INTEGER, INTENT(IN)     :: ICNTRL(20)
   KPP_REAL, INTENT(IN)    :: RCNTRL(20)
   INTEGER, INTENT(INOUT)  :: ISTATUS(20)
   KPP_REAL, INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)    :: IERR
!~~~>  The method parameters   
   INTEGER, PARAMETER :: Smax = 6
   INTEGER  :: Method, ros_S
   KPP_REAL, DIMENSION(Smax) :: ros_M, ros_E, ros_Alpha, ros_Gamma
   KPP_REAL, DIMENSION(Smax*(Smax-1)/2) :: ros_A, ros_C
   KPP_REAL :: ros_ELO
   LOGICAL, DIMENSION(Smax) :: ros_NewF
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables     
   KPP_REAL :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   KPP_REAL :: Hmin, Hmax, Hstart, Hexit
   KPP_REAL :: T
   INTEGER :: i, UplimTol, Max_no_steps
   LOGICAL :: Autonomous, VectorTol
!~~~>   Parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0
   KPP_REAL, PARAMETER :: DeltaMin = 1.0d-5
!~~~>   Functions
   EXTERNAL ode_Fun, ode_Jac, ode_Hess

!~~~>  Initialize statistics
   Nfun = ISTATUS(ifun)
   Njac = ISTATUS(ijac)
   Nstp = ISTATUS(istp)
   Nacc = ISTATUS(iacc)
   Nrej = ISTATUS(irej)
   Ndec = ISTATUS(idec)
   Nsol = ISTATUS(isol)
   Nsng = ISTATUS(isng)
   
!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  
!               the code uses AbsTol(1) and RelTol(1)
!      For Vector tolerances (ICNTRL(2) == 0) 
!               the code uses AbsTol(1:NVAR) and RelTol(1:NVAR)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
         UplimTol  = NVAR
   ELSE 
      VectorTol = .FALSE.
         UplimTol  = 1
   END IF
   
!~~~>   The maximum number of steps admitted
   IF (ICNTRL(3) == 0) THEN
      Max_no_steps = bufsize - 1
   ELSEIF (Max_no_steps > 0) THEN
      Max_no_steps=ICNTRL(3)
   ELSE 
      PRINT * ,'User-selected max no. of steps: ICNTRL(3)=',ICNTRL(3)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN      
   END IF

!~~~>  The particular Rosenbrock method chosen
   IF (ICNTRL(4) == 0) THEN
      Method = 5
   ELSEIF ( (ICNTRL(4) >= 1).AND.(ICNTRL(4) <= 5) ) THEN
      Method = ICNTRL(4)
   ELSE  
      PRINT * , 'User-selected Rosenbrock method: ICNTRL(4)=', Method
      CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN      
   END IF

 
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
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax 
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2d0
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE  
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0d0
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE  
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1d0
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE  
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9d0
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE  
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.d0*Roundoff) &
         .OR. (RelTol(i) >= 1.0d0) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
     
 
!~~~>   Initialize the particular RosenbrockSOA method
   SELECT CASE (Method)
     CASE (1)
       CALL Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   & 
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   & 
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   & 
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, & 
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, & 
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR) 
       RETURN     
   END SELECT
 
!~~~>  Allocate checkpoint space or open checkpoint files
   CALL ros_AllocateDBuffers( ros_S )
   
!~~~>  Forward Rosenbrock and TLM integration   
   CALL ros_TlmInt (NSOA, Y, Y_tlm,              &
        Tstart, Tend, T,                         &
        AbsTol, RelTol,                          &
        ode_Fun, ode_Jac, ode_Hess,              &
!~~~> Rosenbrock method coefficients     
        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )

   PRINT*,'FORWARD STATISTICS'
   PRINT*,'Step=',Nstp,' Acc=',Nacc,             &
        ' Rej=',Nrej, ' Singular=',Nsng
   Nstp = 0
   Nacc = 0
   Nrej = 0
   Nsng = 0

!~~~>  If Forward integration failed return   
   IF (IERR<0) RETURN

!~~~>  Backward ADJ and SOADJ Rosenbrock integration   
   CALL ros_SoaInt (                             &
        NSOA, Y_adj, Y_soa,                      &
        Tstart, Tend, T,                         &
        AbsTol, RelTol,                          &
        ode_Fun, ode_Jac, ode_Hess,              &
!~~~> RosenbrockSOA method coefficients     
        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )


   PRINT*,'ADJOINT STATISTICS'
   PRINT*,'Step=',Nstp,' Acc=',Nacc,             &
        ' Rej=',Nrej, ' Singular=',Nsng

!~~~>  Free checkpoint space or close checkpoint files
   CALL ros_FreeDBuffers

!~~~>  Collect run statistics
   ISTATUS(ifun) = Nfun
   ISTATUS(ijac) = Njac
   ISTATUS(istp) = Nstp
   ISTATUS(iacc) = Nacc
   ISTATUS(irej) = Nrej
   ISTATUS(idec) = Ndec
   ISTATUS(isol) = Nsol
   ISTATUS(isng) = Nsng
!~~~> Last T and H
   RSTATUS(itexit) = T
   RSTATUS(ihexit) = Hexit    

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 END SUBROUTINE RosenbrockSOA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
   KPP_REAL, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from RosenbrockSOA due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)    
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected RosenbrockSOA method not implemented'
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
    CASE (-9)    
      PRINT * , '--> Improper type of adjoint selected'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg
   
     
     
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_TlmInt (NSOA, Y, Y_tlm,          &
        Tstart, Tend, T,                         &
        AbsTol, RelTol,                          &
        ode_Fun, ode_Jac, ode_Hess,              &
!~~~> Rosenbrock method coefficients     
        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
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
   KPP_REAL, INTENT(INOUT) :: Y(NVAR)
!~~~> Input: Number of sensitivity coefficients
   INTEGER, INTENT(IN) :: NSOA 
!~~~> Input: the initial sensitivites at Tstart; Output: the sensitivities at T   
   KPP_REAL, INTENT(INOUT) :: Y_tlm(NVAR,NSOA)
!~~~> Input: integration interval   
   KPP_REAL, INTENT(IN) :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   KPP_REAL, INTENT(OUT) ::  T      
!~~~> Input: tolerances      
   KPP_REAL, INTENT(IN) ::  AbsTol(NVAR), RelTol(NVAR)
!~~~> Input: ode function and its Jacobian      
   EXTERNAL ode_Fun, ode_Jac, ode_Hess
!~~~> Input: The Rosenbrock method parameters   
   INTEGER, INTENT(IN) ::  ros_S
   KPP_REAL, INTENT(IN) :: ros_M(ros_S), ros_E(ros_S),  & 
       ros_Alpha(ros_S), ros_A(ros_S*(ros_S-1)/2),      &
       ros_Gamma(ros_S), ros_C(ros_S*(ros_S-1)/2), ros_ELO
   LOGICAL, INTENT(IN) :: ros_NewF(ros_S)
!~~~> Input: integration parameters   
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   KPP_REAL, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   KPP_REAL, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe 
!~~~> Output: last accepted step   
   KPP_REAL, INTENT(OUT) :: Hexit 
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables        
   KPP_REAL :: Ystage(NVAR*ros_S), Fcn0(NVAR), Fcn(NVAR) 
   KPP_REAL :: K(NVAR*ros_S), Tmp(NVAR)
   KPP_REAL :: Ystage_tlm(NVAR*ros_S,NSOA), Fcn0_tlm(NVAR,NSOA), Fcn_tlm(NVAR,NSOA) 
   KPP_REAL :: K_tlm(NVAR*ros_S,NSOA)
   KPP_REAL :: Hes0(NHESS)
   KPP_REAL :: dFdT(NVAR), dJdT(LU_NONZERO)
   KPP_REAL :: Jac0(LU_NONZERO), Jac(LU_NONZERO), Ghimj(LU_NONZERO)
   KPP_REAL :: H, Hnew, HC, HG, Fac, Tau 
   KPP_REAL :: Err, Yerr(NVAR), Ynew(NVAR), Ynew_tlm(NVAR,NSOA)
   INTEGER :: Pivot(NVAR), Direction, ioffset, joffset, j, istage, mtlm
   LOGICAL :: RejectLastH, RejectMoreH, Singular
!~~~>  Local parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
   KPP_REAL, PARAMETER :: DeltaMin = 1.0d-5
!~~~>  Locally called functions
!   KPP_REAL WLAMCH
!   EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
!~~~>  Initial preparations
   T = Tstart
   Hexit = 0.0_dp
   H = MIN(Hstart,Hmax) 
   IF (ABS(H) <= 10.D0*Roundoff) H = DeltaMin
   
   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF               

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.
   
!~~~> Time loop begins below 

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) ) 
      
   IF ( Nstp > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1d0*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   
!~~~>  Limit H if necessary to avoid going beyond Tend   
   Hexit = H
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL ode_Fun(T,Y,Fcn0)
  
!~~~>   Compute the Jacobian at current time
   CALL ode_Jac(T,Y,Jac0)
  
!~~~>   Compute the Hessian at current time
   CALL ode_Hess(T,Y,Hes0)
   
!~~~>   Compute the TLM function at current time
   DO mtlm = 1, NSOA
      CALL Jac_SP_Vec ( Jac0, Y_tlm(1:NVAR,mtlm), Fcn0_tlm(1:NVAR,mtlm) )
   END DO  
   
!~~~>  Compute the function and Jacobian derivatives with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, ode_Fun, dFdT )
      CALL ros_JacTimeDerivative ( T, Roundoff, Y, &
                Jac0, ode_Jac, dJdT )
   END IF
 
!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO  
   
   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1),&
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S
      
      ! Current istage vector is K(ioffset+1:ioffset+NVAR:ioffset+1+NVAR-1)
       ioffset = NVAR*(istage-1)
         
      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         CALL WCOPY(NVAR,Y,1,Ystage(ioffset+1:ioffset+NVAR),1)
         DO mtlm=1,NSOA                    
            CALL WCOPY(NVAR,Y_tlm(1:NVAR,mtlm),1,Ystage_tlm(ioffset+1:ioffset+NVAR,mtlm),1)
         END DO                    
         CALL WCOPY(NVAR,Fcn0,1,Fcn,1)
         CALL WCOPY(NVAR*NSOA,Fcn0_tlm,1,Fcn_tlm,1)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL WCOPY(NVAR,Y,1,Ystage(ioffset+1:ioffset+NVAR),1)
         DO mtlm=1,NSOA                    
            CALL WCOPY(NVAR,Y_tlm(1:NVAR,mtlm),1,Ystage_tlm(ioffset+1:ioffset+NVAR,mtlm),1)
         END DO                    
         DO j = 1, istage-1
           joffset = NVAR*(j-1)
           CALL WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j),    &
                     K(joffset+1:joffset+NVAR),1,Ystage(ioffset+1:ioffset+NVAR),1) 
           DO mtlm=1,NSOA                    
              CALL WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
                     K_tlm(joffset+1:joffset+NVAR,mtlm),1,Ystage_tlm(ioffset+1:ioffset+NVAR,mtlm),1) 
           END DO                    
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL ode_Fun(Tau,Ystage(ioffset+1:ioffset+NVAR),Fcn)
         CALL ode_Jac(Tau,Ystage(ioffset+1:ioffset+NVAR),Jac)
         DO mtlm=1,NSOA 
           CALL Jac_SP_Vec ( Jac, Ystage_tlm(ioffset+1:ioffset+NVAR,mtlm), Fcn_tlm(1:NVAR,mtlm) )
         END DO              
       END IF ! if istage == 1 elseif ros_NewF(istage)
       
       CALL WCOPY(NVAR,Fcn,1,K(ioffset+1:ioffset+NVAR),1)
       DO mtlm=1,NSOA   
          CALL WCOPY(NVAR,Fcn_tlm(1:NVAR,mtlm),1,K_tlm(ioffset+1:ioffset+NVAR,mtlm),1)
       END DO                
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL WAXPY(NVAR,HC,K(NVAR*(j-1)+1:NVAR*j),1,K(ioffset+1:ioffset+NVAR),1)
         DO mtlm=1,NSOA 
           CALL WAXPY(NVAR,HC,K_tlm(NVAR*(j-1)+1:NVAR*j,mtlm),1,K_tlm(ioffset+1:ioffset+NVAR,mtlm),1)
         END DO              
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL WAXPY(NVAR,HG,dFdT,1,K(ioffset+1:ioffset+NVAR),1)
         DO mtlm=1,NSOA 
           CALL Jac_SP_Vec ( dJdT, Ystage_tlm(ioffset+1:ioffset+NVAR,mtlm), Fcn_tlm(1:NVAR,mtlm) )
           CALL WAXPY(NVAR,HG,Fcn_tlm(1:NVAR,mtlm),1,K_tlm(ioffset+1:ioffset+NVAR,mtlm),1)
         END DO              
       END IF
       CALL ros_Solve('N', Ghimj, Pivot, K(ioffset+1:ioffset+NVAR))
       DO mtlm=1,NSOA   
         CALL Hess_Vec ( Hes0,  K(ioffset+1:ioffset+NVAR), Y_tlm(1:NVAR,mtlm), Tmp )
         CALL WAXPY(NVAR,ONE,Tmp,1,K_tlm(ioffset+1:ioffset+NVAR,mtlm),1)
         CALL ros_Solve('N', Ghimj, Pivot, K_tlm(ioffset+1:ioffset+NVAR,mtlm))
       END DO                
      
   END DO Stage     
            

!~~~>  Compute the new solution 
   CALL WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
      CALL WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1:NVAR*j),1,Ynew,1)
   END DO
   DO mtlm=1,NSOA       
     CALL WCOPY(NVAR,Y_tlm(1:NVAR,mtlm),1,Ynew_tlm(1:NVAR,mtlm),1)
     DO j=1,ros_S
       joffset = NVAR*(j-1)
       CALL WAXPY(NVAR,ros_M(j),K_tlm(joffset+1:joffset+NVAR,mtlm),1,Ynew_tlm(1:NVAR,mtlm),1)
     END DO
   END DO

!~~~>  Compute the error estimation 
   CALL WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S     
        CALL WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1:NVAR*j),1,Yerr,1)
   END DO 
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac  

!~~~>  Check the error magnitude and adjust step size
   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      Nacc = Nacc+1
      !~~~>  Checkpoints for stage values and vectors
      CALL ros_DPush( ros_S, NSOA, T, H, Ystage, K, Ystage_tlm, K_tlm )
      !~~~>  Accept new solution, etc.
      CALL WCOPY(NVAR,Ynew,1,Y,1)
      CALL WCOPY(NVAR*NSOA,Ynew_tlm,1,Y_tlm,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H) 
      END IF   
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
      IF (Nacc >= 1) THEN
         Nrej = Nrej+1
      END IF    
   END IF ! Err <= 1

   END DO UntilAccepted 

   END DO TimeLoop 
   
!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE ros_TlmInt
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_SoaInt (                         &
        NSOA, Lambda, Sigma,                     &
        Tstart, Tend, T,                         &
        AbsTol, RelTol,                          &
        ode_Fun, ode_Jac, ode_Hess,              &
!~~~> RosenbrockSOA method coefficients     
        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockSOA method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
   INTEGER, INTENT(IN)     :: NSOA
!~~~> First order adjoint   
   KPP_REAL, INTENT(INOUT) :: Lambda(NVAR)
!~~~> Second order adjoint   
   KPP_REAL, INTENT(INOUT) :: Sigma(NVAR,NSOA)
!~~~> Input: integration interval   
   KPP_REAL, INTENT(IN) :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   KPP_REAL, INTENT(OUT) ::  T      
!~~~> Input: tolerances      
   KPP_REAL, INTENT(IN) ::  AbsTol(NVAR), RelTol(NVAR)
!~~~> Input: ode function and its Jacobian      
   EXTERNAL ode_Fun, ode_Jac, ode_Hess
!~~~> Input: The RosenbrockSOA method parameters   
   INTEGER, INTENT(IN) ::  ros_S
   KPP_REAL, INTENT(IN) :: ros_M(ros_S), ros_E(ros_S),  & 
       ros_Alpha(ros_S), ros_A(ros_S*(ros_S-1)/2),      &
       ros_Gamma(ros_S), ros_C(ros_S*(ros_S-1)/2), ros_ELO
   LOGICAL, INTENT(IN) :: ros_NewF(ros_S)
!~~~> Input: integration parameters   
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   KPP_REAL, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   KPP_REAL, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe 
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables        
   !KPP_REAL :: Ystage_adj(NVAR,NSOA)
   !KPP_REAL :: dFdT(NVAR)
   KPP_REAL :: Ystage(NVAR*ros_S), K(NVAR*ros_S)
   KPP_REAL :: Ystage_tlm(NVAR*ros_S,NSOA), K_tlm(NVAR*ros_S,NSOA)
   KPP_REAL :: U(NVAR*ros_S), V(NVAR*ros_S)
   KPP_REAL :: W(NVAR*ros_S,NSOA), Z(NVAR*ros_S,NSOA)
   KPP_REAL :: Jac(LU_NONZERO), dJdT(LU_NONZERO), Ghimj(LU_NONZERO)
   KPP_REAL :: Hes0(NHESS), Hes1(NHESS), dHdT(NHESS)
   KPP_REAL :: Tmp(NVAR), Tmp2(NVAR)
   KPP_REAL :: H, HC, HA, Tau 
   INTEGER :: Pivot(NVAR), Direction, i, ioffset, joffset
   INTEGER :: msoa, j, istage
!~~~>  Local parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
   KPP_REAL, PARAMETER :: DeltaMin = 1.0d-5
!~~~>  Locally called functions
!    KPP_REAL WLAMCH
!   EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
   
   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF               

!~~~> Time loop begins below 
TimeLoop: DO WHILE ( stack_ptr > 0 )
        
   !~~~>  Recover checkpoints for stage values and vectors
   CALL ros_DPop( ros_S, NSOA, T, H, Ystage, K, Ystage_tlm, K_tlm )

   Nstp = Nstp+1

!~~~>    Compute LU decomposition 
   CALL ode_Jac(T,Ystage(1:NVAR),Ghimj)
   CALL WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
   Tau = ONE/(Direction*H*ros_Gamma(1))
   DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+Tau
   END DO
   CALL ros_Decomp( Ghimj, Pivot, j )
            
!~~~>   Compute Hessian at the beginning of the interval
   CALL ode_Hess(T,Ystage(1),Hes0)
   
!~~~>   Compute the stages
Stage: DO istage = ros_S, 1, -1
      
      !~~~> Current istage offset. 
       ioffset = NVAR*(istage-1)
      
      !~~~> Compute U
       CALL WCOPY(NVAR,Lambda,1,U(ioffset+1:ioffset+NVAR),1)
       CALL WSCAL(NVAR,ros_M(istage),U(ioffset+1:ioffset+NVAR),1)
       DO j = istage+1, ros_S
         joffset = NVAR*(j-1)
         HA = ros_A((j-1)*(j-2)/2+istage)
         HC = ros_C((j-1)*(j-2)/2+istage)/(Direction*H)
         CALL WAXPY(NVAR,HA,V(joffset+1:joffset+NVAR),1,U(ioffset+1:ioffset+NVAR),1) 
         CALL WAXPY(NVAR,HC,U(joffset+1:joffset+NVAR),1,U(ioffset+1:ioffset+NVAR),1) 
       END DO
       CALL ros_Solve('T', Ghimj, Pivot, U(ioffset+1:ioffset+NVAR))
      !~~~> Compute W
       DO msoa = 1, NSOA
         CALL WCOPY(NVAR,Sigma(1:NVAR,msoa),1,W(ioffset+1:ioffset+NVAR,msoa),1)
         CALL WSCAL(NVAR,ros_M(istage),W(ioffset+1:ioffset+NVAR,msoa),1)
       END DO
       DO j = istage+1, ros_S
         joffset = NVAR*(j-1)
         HA = ros_A((j-1)*(j-2)/2+istage)
         HC = ros_C((j-1)*(j-2)/2+istage)/(Direction*H)
         DO msoa = 1, NSOA
           CALL WAXPY(NVAR,HA,   &
               Z(joffset+1:joffset+NVAR,msoa),1,W(ioffset+1:ioffset+NVAR,msoa),1) 
           CALL WAXPY(NVAR,HC,   &
               W(joffset+1:joffset+NVAR,msoa),1,W(ioffset+1:ioffset+NVAR,msoa),1) 
         END DO
       END DO
       DO msoa = 1, NSOA
         CALL HessTR_Vec( Hes0, U(ioffset+1:ioffset+NVAR), Ystage_tlm(ioffset+1:ioffset+NVAR,msoa), Tmp )
         CALL WAXPY(NVAR,ONE,Tmp,1,W(ioffset+1:ioffset+NVAR,msoa),1)
         CALL ros_Solve('T', Ghimj, Pivot, W(ioffset+1:ioffset+NVAR,msoa))
       END DO       
      !~~~> Compute V 
       Tau = T + ros_Alpha(istage)*Direction*H
       CALL ode_Jac(Tau,Ystage(ioffset+1:ioffset+NVAR),Jac)
       CALL JacTR_SP_Vec(Jac,U(ioffset+1:ioffset+NVAR),V(ioffset+1:ioffset+NVAR)) 
      !~~~> Compute Z 
       CALL ode_Hess(T,Ystage(ioffset+1:ioffset+NVAR),Hes1)
       DO msoa = 1, NSOA
         CALL JacTR_SP_Vec(Jac,W(ioffset+1:ioffset+NVAR,msoa),Z(ioffset+1:ioffset+NVAR,msoa)) 
         CALL HessTR_Vec( Hes1, U(ioffset+1:ioffset+NVAR), Ystage_tlm(ioffset+1:ioffset+NVAR,msoa), Tmp )
         CALL WAXPY(NVAR,ONE,Tmp,1,Z(ioffset+1:ioffset+NVAR,msoa),1)
       END DO  
             
   END DO Stage     

   IF (.NOT.Autonomous) THEN
!~~~>  Compute the Jacobian derivative with respect to T. 
!      Last "Jac" computed for stage 1
      CALL ros_JacTimeDerivative ( T, Roundoff, Ystage(1),        &
                Jac, ode_Jac, dJdT )
!~~~>  Compute the Hessian derivative with respect to T. 
!      Last "Jac" computed for stage 1
      CALL ros_HesTimeDerivative ( T, Roundoff, Ystage(1),        &
                Hes0, ode_Hess, dHdT )
   END IF

!~~~>  Compute the new solution 
   
      !~~~>  Compute Lambda 
      DO istage=1,ros_S
         ioffset = NVAR*(istage-1)
         ! Add V_i
         CALL WAXPY(NVAR,ONE,V(ioffset+1:ioffset+NVAR),1,Lambda,1)
         ! Add (H0xK_i)^T * U_i
         CALL HessTR_Vec ( Hes0, U(ioffset+1:ioffset+NVAR), K(ioffset+1:ioffset+NVAR), Tmp )
         CALL WAXPY(NVAR,ONE,Tmp,1,Lambda,1)
      END DO
     ! Add H * dJac_dT_0^T * \sum(gamma_i U_i)
     ! Tmp holds sum gamma_i U_i
      IF (.NOT.Autonomous) THEN
        Tmp(1:NVAR) = ZERO
        DO istage = 1, ros_S
          ioffset = NVAR*(istage-1)
          CALL WAXPY(NVAR,ros_Gamma(istage),U(ioffset+1:ioffset+NVAR),1,Tmp,1)
        END DO  
        CALL JacTR_SP_Vec(dJdT,Tmp,Tmp2) 
        CALL WAXPY(NVAR,H,Tmp2,1,Lambda,1)
      END IF ! .NOT.Autonomous
      
      !~~~>  Compute Sigma 
    DO msoa = 1, NSOA
    
      DO istage=1,ros_S
         ioffset = NVAR*(istage-1)
         ! Add Z_i
         CALL WAXPY(NVAR,ONE,Z(ioffset+1:ioffset+NVAR,msoa),1,Sigma(1:NVAR,msoa),1)
         ! Add (Hess_0 x K_i)^T * W_i
         CALL HessTR_Vec ( Hes0, W(ioffset+1:ioffset+NVAR,msoa), K(ioffset+1:ioffset+NVAR), Tmp )
         CALL WAXPY(NVAR,ONE,Tmp,1,Sigma(1:NVAR,msoa),1)
         ! Add (Hess_0 x K_tlm_i)^T * U_i
         CALL HessTR_Vec ( Hes0, U(ioffset+1:ioffset+NVAR), K_tlm(ioffset+1:ioffset+NVAR,msoa), Tmp )
         CALL WAXPY(NVAR,ONE,Tmp,1,Sigma(1:NVAR,msoa),1)
      END DO  
      
      !~~~> Add high derivative terms
      DO istage=1,ros_S
         ioffset = NVAR*(istage-1)
         CALL ros_HighDerivative ( T, Roundoff, Ystage(1), Hes0, K(ioffset+1:ioffset+NVAR), &
                 U(ioffset+1:ioffset+NVAR), Ystage_tlm(1:NVAR,msoa), ode_Hess, Tmp)
         CALL WAXPY(NVAR,ONE,Tmp,1,Sigma(1:NVAR,msoa),1)
      END DO  

      IF (.NOT.Autonomous) THEN
        ! Add H * dJac_dT_0^T * \sum(gamma_i W_i)
        ! Tmp holds sum gamma_i W_i
        Tmp(1:NVAR) = ZERO
        DO istage = 1, ros_S
          ioffset = NVAR*(istage-1)
          CALL WAXPY(NVAR,ros_Gamma(istage),W(ioffset+1:ioffset+NVAR,msoa),1,Tmp,1)
        END DO  
        CALL JacTR_SP_Vec(dJdT,Tmp,Tmp2) 
        CALL WAXPY(NVAR,H,Tmp2,1,Sigma(1:NVAR,msoa),1)
        ! Add H * ( dHess_dT_0 x Y_tlm_0)^T * \sum(gamma_i U_i)
        ! Tmp holds sum gamma_i U_i
        Tmp(1:NVAR) = ZERO
        DO istage = 1, ros_S
          ioffset = NVAR*(istage-1)
          CALL WAXPY(NVAR,ros_Gamma(istage),U(ioffset+1:ioffset+NVAR),1,Tmp,1)
        END DO  
        CALL HessTR_Vec ( dHdT, Tmp, Ystage_tlm(ioffset+1:ioffset+NVAR,msoa), Tmp2 )
        CALL WAXPY(NVAR,H,Tmp2,1,Sigma(1:NVAR,msoa),1)
      END IF ! .NOT.Autonomous
      
    END DO ! msoa


   END DO TimeLoop 
   
!~~~> Save last state
   
!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE ros_SoaInt
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
   
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  KPP_REAL FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, & 
               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

! Input arguments   
   KPP_REAL, INTENT(IN) :: Y(NVAR), Ynew(NVAR),    &
          Yerr(NVAR), AbsTol(NVAR), RelTol(NVAR)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables     
   KPP_REAL :: Err, Scale, Ymax   
   INTEGER  :: i
   KPP_REAL, PARAMETER :: ZERO = 0.0d0
   
   Err = ZERO
   DO i=1,NVAR
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/NVAR)

   ros_ErrorNorm = MAX(Err,1.0d-10)
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END FUNCTION ros_ErrorNorm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, ode_Fun, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

!~~~> Input arguments   
   KPP_REAL, INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR) 
   EXTERNAL ode_Fun
!~~~> Output arguments   
   KPP_REAL, INTENT(OUT) :: dFdT(NVAR)   
!~~~> Local variables     
   KPP_REAL :: Delta  
   KPP_REAL, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
   
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL ode_Fun(T+Delta,Y,dFdT)
   CALL WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL WSCAL(NVAR,(ONE/Delta),dFdT,1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE ros_FunTimeDerivative
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_JacTimeDerivative ( T, Roundoff, Y, &
                Jac0, ode_Jac, dJdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the Jacobian by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

!~~~> Input arguments   
   KPP_REAL, INTENT(IN) :: T, Roundoff, Y(NVAR), Jac0(LU_NONZERO) 
   EXTERNAL ode_Jac
!~~~> Output arguments   
   KPP_REAL, INTENT(OUT) :: dJdT(LU_NONZERO)   
!~~~> Local variables     
   KPP_REAL Delta  
   KPP_REAL, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
   
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL ode_Jac( T+Delta, Y, dJdT )   
   CALL WAXPY(LU_NONZERO,(-ONE),Jac0,1,dJdT,1)
   CALL WSCAL(LU_NONZERO,(ONE/Delta),dJdT,1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE ros_JacTimeDerivative
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_HesTimeDerivative ( T, Roundoff, Y, Hes0, ode_Hess, dHdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the Hessian by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

!~~~> Input arguments   
   KPP_REAL, INTENT(IN) :: T, Roundoff, Y(NVAR), Hes0(NHESS) 
   EXTERNAL ode_Hess
!~~~> Output arguments   
   KPP_REAL, INTENT(OUT) :: dHdT(NHESS)   
!~~~> Local variables     
   KPP_REAL Delta  
   KPP_REAL, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
   
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL ode_Hess( T+Delta, Y, dHdT )   
   CALL WAXPY(NHESS,(-ONE),Hes0,1,dHdT,1)
   CALL WSCAL(NHESS,(ONE/Delta),dHdT,1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE ros_HesTimeDerivative
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_HighDerivative ( T, Roundoff, Y, Hes0, K, U, Y_tlm, &
                                  ode_Hess, Term)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> High order derivative by finite differences:
!     d/dy { (Hes0 x K_i)^T * U_i } * Y_tlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

!~~~> Input arguments   
   KPP_REAL, INTENT(IN) :: T, Roundoff, Y(NVAR), Hes0(NHESS) 
   KPP_REAL, INTENT(IN) :: K(NVAR), U(NVAR), Y_tlm(NVAR) 
   EXTERNAL ode_Hess
!~~~> Output arguments   
   KPP_REAL, INTENT(OUT) :: Term(NVAR)   
!~~~> Local variables     
   KPP_REAL :: Delta, Y1(NVAR), Hes1(NHESS), Tmp(NVAR)  
   KPP_REAL, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
   
   CALL HessTR_Vec ( Hes0, U, K, Tmp )
   
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   Y1(1:NVAR) = Y(1:NVAR) + Delta*Y_tlm(1:NVAR)
   CALL ode_Hess( T, Y1, Hes1 )   
   ! Add (Hess_0 x K_i)^T * U_i
   CALL HessTR_Vec ( Hes1, U, K, Term )

   CALL WAXPY(NVAR,(-ONE),Tmp,1,Term,1)
   CALL WSCAL(NVAR,(ONE/Delta),Term,1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE ros_HighDerivative
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
   KPP_REAL, INTENT(IN) ::  gam, Jac0(LU_NONZERO)
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments   
   KPP_REAL, INTENT(OUT) :: Ghimj(LU_NONZERO)
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(NVAR)
!~~~> Inout arguments   
   KPP_REAL, INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables     
   INTEGER  :: i, ising, Nconsecutive
   KPP_REAL ::  ghinv
   KPP_REAL, PARAMETER :: ONE  = 1.0d0, HALF = 0.5d0
   
   Nconsecutive = 0
   Singular = .TRUE.
   
   DO WHILE (Singular)
   
!~~~>    Construct Ghimj = 1/(H*ham) - Jac0
     CALL WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO
!~~~>    Compute LU decomposition 
     CALL ros_Decomp( Ghimj, Pivot, ising )
     IF (ising == 0) THEN
!~~~>    If successful done 
        Singular = .FALSE. 
     ELSE ! ising .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        Nsng = Nsng+1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE. 
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decomps
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ising 
         
   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix

  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Decomp( A, Pivot, ising )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the LU decomposition   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Inout variables     
   KPP_REAL, INTENT(INOUT) :: A(LU_NONZERO)
!~~~> Output variables     
   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   
   CALL KppDecomp ( A, ising )
!~~~> Note: for a full matrix use Lapack:
!     CALL  DGETRF( NVAR, NVAR, A, NVAR, Pivot, ising ) 
   Pivot(1) = 1
    
   Ndec = Ndec + 1

  END SUBROUTINE ros_Decomp
 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( C, A, Pivot, b )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomp)   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   IMPLICIT NONE
!~~~> Input variables 
   CHARACTER, INTENT(IN) :: C    
   KPP_REAL, INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)
!~~~> InOut variables     
   KPP_REAL, INTENT(INOUT) :: b(NVAR)
   
   SELECT CASE (C)
     CASE ('N')
         CALL KppSolve( A, b )
     CASE ('T')
         CALL KppSolveTR( A, b, b )
     CASE DEFAULT
         PRINT*,'Unknown C = (',C,') in ros_Solve'
         STOP
   END SELECT
!~~~> Note: for a full matrix use Lapack:
!     NRHS = 1
!     CALL  DGETRS( C, NVAR , NRHS, A, NVAR, Pivot, b, NVAR, INFO )
     
   Nsol = Nsol+1

  END SUBROUTINE ros_Solve
  
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  
   INTEGER, PARAMETER :: S = 2
   INTEGER, INTENT(OUT) ::  ros_S
   KPP_REAL, DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   KPP_REAL, DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   KPP_REAL, INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   KPP_REAL :: g
   
    g = 1.0d0 + 1.0d0/SQRT(2.0d0)
   
!~~~> Name of the method
    ros_Name = 'ROS-2'   
!~~~> Number of stages
    ros_S = S
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
   
    ros_A(1) = (1.d0)/g
    ros_C(1) = (-2.d0)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.d0)/(2.d0*g)
    ros_M(2)= (1.d0)/(2.d0*g)
! E_i = Coefficients for error estimator    
    ros_E(1) = 1.d0/(2.d0*g)
    ros_E(2) = 1.d0/(2.d0*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0d0
    ros_Alpha(2) = 1.0d0 
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = g
    ros_Gamma(2) =-g
   
 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  IMPLICIT NONE
  
   INTEGER, PARAMETER :: S = 3
   INTEGER, INTENT(OUT) ::  ros_S
   KPP_REAL, DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   KPP_REAL, DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   KPP_REAL, INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   
!~~~> Name of the method
   ros_Name = 'ROS-3'   
!~~~> Number of stages
   ros_S = S
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
   ros_A(1)= 1.d0
   ros_A(2)= 1.d0
   ros_A(3)= 0.d0

   ros_C(1) = -0.10156171083877702091975600115545d+01
   ros_C(2) =  0.40759956452537699824805835358067d+01
   ros_C(3) =  0.92076794298330791242156818474003d+01
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1d+01
   ros_M(2) =  0.61697947043828245592553615689730d+01
   ros_M(3) = -0.42772256543218573326238373806514d+00
! E_i = Coefficients for error estimator    
   ros_E(1) =  0.5d+00
   ros_E(2) = -0.29079558716805469821718236208017d+01
   ros_E(3) =  0.22354069897811569627360909276199d+00
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0d+00
   ros_Alpha(2)= 0.43586652150845899941601945119356d+00
   ros_Alpha(3)= 0.43586652150845899941601945119356d+00
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1)= 0.43586652150845899941601945119356d+00
   ros_Gamma(2)= 0.24291996454816804366592249683314d+00
   ros_Gamma(3)= 0.21851380027664058511513169485832d+01

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3 
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  
   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   KPP_REAL, DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   KPP_REAL, DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   KPP_REAL, INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   
   
!~~~> Name of the method
   ros_Name = 'ROS-4'   
!~~~> Number of stages
   ros_S = S
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
   ros_A(1) = 0.2000000000000000d+01
   ros_A(2) = 0.1867943637803922d+01
   ros_A(3) = 0.2344449711399156d+00
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0D0

   ros_C(1) =-0.7137615036412310d+01
   ros_C(2) = 0.2580708087951457d+01
   ros_C(3) = 0.6515950076447975d+00
   ros_C(4) =-0.2137148994382534d+01
   ros_C(5) =-0.3214669691237626d+00
   ros_C(6) =-0.6949742501781779d+00
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735d+01
   ros_M(2) = 0.2870493262186792d+00
   ros_M(3) = 0.4353179431840180d+00
   ros_M(4) = 0.1093502252409163d+01
!~~~> E_i  = Coefficients for error estimator    
   ros_E(1) =-0.2815431932141155d+00
   ros_E(2) =-0.7276199124938920d-01
   ros_E(3) =-0.1082196201495311d+00
   ros_E(4) =-0.1093502252409163d+01
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.D0
   ros_Alpha(2) = 0.1145640000000000d+01
   ros_Alpha(3) = 0.6552168638155900d+00
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1) = 0.5728200000000000d+00
   ros_Gamma(2) =-0.1769193891319233d+01
   ros_Gamma(3) = 0.7592633437920482d+00
   ros_Gamma(4) =-0.1049021087100450d+00

  END SUBROUTINE Ros4
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  
   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   KPP_REAL, DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   KPP_REAL, DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   KPP_REAL, INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   
   
!~~~> Name of the method
   ros_Name = 'RODAS-3'   
!~~~> Number of stages
   ros_S = S
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
 
   ros_A(1) = 0.0d+00
   ros_A(2) = 2.0d+00
   ros_A(3) = 0.0d+00
   ros_A(4) = 2.0d+00
   ros_A(5) = 0.0d+00
   ros_A(6) = 1.0d+00

   ros_C(1) = 4.0d+00
   ros_C(2) = 1.0d+00
   ros_C(3) =-1.0d+00
   ros_C(4) = 1.0d+00
   ros_C(5) =-1.0d+00 
   ros_C(6) =-(8.0d+00/3.0d+00) 
         
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0d+00
   ros_M(2) = 0.0d+00
   ros_M(3) = 1.0d+00
   ros_M(4) = 1.0d+00
!~~~> E_i  = Coefficients for error estimator    
   ros_E(1) = 0.0d+00
   ros_E(2) = 0.0d+00
   ros_E(3) = 0.0d+00
   ros_E(4) = 1.0d+00
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0d+00    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0d+00
   ros_Alpha(2) = 0.0d+00
   ros_Alpha(3) = 1.0d+00
   ros_Alpha(4) = 1.0d+00
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1) = 0.5d+00
   ros_Gamma(2) = 1.5d+00
   ros_Gamma(3) = 0.0d+00
   ros_Gamma(4) = 0.0d+00

  END SUBROUTINE Rodas3
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
             ros_Gamma,ros_NewF,ros_ELO,ros_Name)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  
   INTEGER, PARAMETER :: S=6
   INTEGER, INTENT(OUT) ::  ros_S
   KPP_REAL, DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   KPP_REAL, DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   KPP_REAL, INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   

!~~~> Name of the method
    ros_Name = 'RODAS-4'   
!~~~> Number of stages
    ros_S = S

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000d0
    ros_Alpha(2) = 0.386d0
    ros_Alpha(3) = 0.210d0 
    ros_Alpha(4) = 0.630d0
    ros_Alpha(5) = 1.000d0
    ros_Alpha(6) = 1.000d0
        
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = 0.2500000000000000d+00
    ros_Gamma(2) =-0.1043000000000000d+00
    ros_Gamma(3) = 0.1035000000000000d+00
    ros_Gamma(4) =-0.3620000000000023d-01
    ros_Gamma(5) = 0.0d0
    ros_Gamma(6) = 0.0d0

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
    ros_A(1) = 0.1544000000000000d+01
    ros_A(2) = 0.9466785280815826d+00
    ros_A(3) = 0.2557011698983284d+00
    ros_A(4) = 0.3314825187068521d+01
    ros_A(5) = 0.2896124015972201d+01
    ros_A(6) = 0.9986419139977817d+00
    ros_A(7) = 0.1221224509226641d+01
    ros_A(8) = 0.6019134481288629d+01
    ros_A(9) = 0.1253708332932087d+02
    ros_A(10) =-0.6878860361058950d+00
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0d+00

    ros_C(1) =-0.5668800000000000d+01
    ros_C(2) =-0.2430093356833875d+01
    ros_C(3) =-0.2063599157091915d+00
    ros_C(4) =-0.1073529058151375d+00
    ros_C(5) =-0.9594562251023355d+01
    ros_C(6) =-0.2047028614809616d+02
    ros_C(7) = 0.7496443313967647d+01
    ros_C(8) =-0.1024680431464352d+02
    ros_C(9) =-0.3399990352819905d+02
    ros_C(10) = 0.1170890893206160d+02
    ros_C(11) = 0.8083246795921522d+01
    ros_C(12) =-0.7981132988064893d+01
    ros_C(13) =-0.3152159432874371d+02
    ros_C(14) = 0.1631930543123136d+02
    ros_C(15) =-0.6058818238834054d+01

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0d+00
    ros_M(6) = 1.0d+00

!~~~> E_i  = Coefficients for error estimator    
    ros_E(1) = 0.0d+00
    ros_E(2) = 0.0d+00
    ros_E(3) = 0.0d+00
    ros_E(4) = 0.0d+00
    ros_E(5) = 0.0d+00
    ros_E(6) = 1.0d+00

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
    ros_ELO = 4.0d0
     
  END SUBROUTINE Rodas4
  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Fun_Template( T, Y, Ydot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE function call.
!  Updates the rate coefficients (and possibly the fixed species) at each call    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Input variables     
   KPP_REAL T, Y(NVAR)
!~~~> Output variables     
   KPP_REAL Ydot(NVAR)
!~~~> Local variables
   KPP_REAL Told     

   Told = TIME
   TIME = T
   CALL Update_SUN()
   CALL Update_RCONST()
   CALL Fun( Y, FIX, RCONST, Ydot )
   TIME = Told
     
   Nfun = Nfun+1
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE Fun_Template
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Jac_Template( T, Y, Jcb )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE Jacobian call.
!  Updates the rate coefficients (and possibly the fixed species) at each call    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
!~~~> Input variables     
    KPP_REAL T, Y(NVAR)
!~~~> Output variables     
    KPP_REAL Jcb(LU_NONZERO)
!~~~> Local variables
    KPP_REAL Told     

    Told = TIME
    TIME = T   
    CALL Update_SUN()
    CALL Update_RCONST()
    CALL Jac_SP( Y, FIX, RCONST, Jcb )
    TIME = Told
     
    Njac = Njac+1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE Jac_Template                                      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Hess_Template( T, Y, Hes )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the ODE Hessian call.
!  Updates the rate coefficients (and possibly the fixed species) at each call    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Input variables     
    KPP_REAL T, Y(NVAR)
!~~~> Output variables     
    KPP_REAL Hes(NHESS)
!~~~> Local variables
    KPP_REAL Told     

    Told = TIME
    TIME = T   
    CALL Update_SUN()
    CALL Update_RCONST()
    CALL Hessian( Y, FIX, RCONST, Hes )
    TIME = Told

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE Hess_Template                                      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_AllocateDBuffers( S )
!~~~>  Allocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INTEGER :: i, S
   
   ALLOCATE( buf_H(bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer H'; STOP
   END IF   
   ALLOCATE( buf_T(bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer T'; STOP
   END IF   
   ALLOCATE( buf_Y(NVAR*S,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer Y'; STOP
   END IF   
   ALLOCATE( buf_K(NVAR*S,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer K'; STOP
   END IF   
   ALLOCATE( buf_Y_tlm(NVAR*S,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer Y_tlm'; STOP
   END IF   
   ALLOCATE( buf_K_tlm(NVAR*S,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer K_tlm'; STOP
   END IF   
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 END SUBROUTINE ros_AllocateDBuffers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_FreeDBuffers
!~~~>  Dallocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INTEGER :: i
   
   DEALLOCATE( buf_H, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer H'; STOP
   END IF   
   DEALLOCATE( buf_T, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer T'; STOP
   END IF   
   DEALLOCATE( buf_Y, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer Y'; STOP
   END IF   
   DEALLOCATE( buf_K, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer K'; STOP
   END IF   
   DEALLOCATE( buf_Y_tlm, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer Y_tlm'; STOP
   END IF   
   DEALLOCATE( buf_K_tlm, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer K_tlm'; STOP
   END IF   
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 END SUBROUTINE ros_FreeDBuffers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_DPush( S, NSOA, T, H, Ystage, K, Ystage_tlm, K_tlm )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   INTEGER, INTENT(IN) :: S    ! no of stages
   INTEGER, INTENT(IN) :: NSOA ! no of second order adjoints
   KPP_REAL :: T, H, Ystage(NVAR*S), K(NVAR*S)
   KPP_REAL :: Ystage_tlm(NVAR*S,NSOA), K_tlm(NVAR*S,NSOA) 
   
   stack_ptr = stack_ptr + 1
   IF ( stack_ptr > bufsize ) THEN
     PRINT*,'Push failed: buffer overflow'
     STOP
   END IF  
   buf_H( stack_ptr ) = H
   buf_T( stack_ptr ) = T
   CALL WCOPY(NVAR*S,Ystage,1,buf_Y(1:NVAR*S,stack_ptr),1)
   CALL WCOPY(NVAR*S,K,1,buf_K(1:NVAR*S,stack_ptr),1)
   CALL WCOPY(NVAR*S*NSOA,Ystage_tlm,1,buf_Y_tlm(1:NVAR*S*NSOA,stack_ptr),1)
   CALL WCOPY(NVAR*S*NSOA,K_tlm,1,buf_K_tlm(1:NVAR*S*NSOA,stack_ptr),1)
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 END SUBROUTINE ros_DPush
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_DPop( S, NSOA, T, H, Ystage, K, Ystage_tlm, K_tlm  )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   INTEGER, INTENT(IN) :: S    ! no of stages
   INTEGER, INTENT(IN) :: NSOA ! no of second order adjoints
   KPP_REAL, INTENT(OUT) :: T, H, Ystage(NVAR*S), K(NVAR*S)
   KPP_REAL, INTENT(OUT) :: Ystage_tlm(NVAR*S,NSOA), K_tlm(NVAR*S,NSOA) 
   
   IF ( stack_ptr <= 0 ) THEN
     PRINT*,'Pop failed: empty buffer'
     STOP
   END IF  
   H = buf_H( stack_ptr )
   T = buf_T( stack_ptr )
   CALL WCOPY(NVAR*S,buf_Y(1:NVAR*S,stack_ptr),1,Ystage,1)
   CALL WCOPY(NVAR*S,buf_K(1:NVAR*S,stack_ptr),1,K,1)
   CALL WCOPY(NVAR*S*NSOA,buf_Y_tlm(1:NVAR*S*NSOA,stack_ptr),1,Ystage_tlm,1)
   CALL WCOPY(NVAR*S*NSOA,buf_K_tlm(1:NVAR*S*NSOA,stack_ptr),1,K_tlm,1)

   stack_ptr = stack_ptr - 1
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 END SUBROUTINE ros_DPop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE KPP_ROOT_Integrator

