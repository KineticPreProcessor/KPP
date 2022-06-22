MODULE KPP_ROOT_Integrator

  USE KPP_ROOT_Parameters, ONLY: NVAR, NFIX, NSPEC, LU_NONZERO
  USE KPP_ROOT_Global
  IMPLICIT NONE
  PUBLIC
  SAVE

  INTEGER, PARAMETER :: Nhnew = 3
  
CONTAINS

SUBROUTINE INTEGRATE( TIN, TOUT, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

   IMPLICIT NONE

   KPP_REAL, INTENT(IN) :: TIN  ! Start Time
   KPP_REAL, INTENT(IN) :: TOUT ! End Time
   !~~~> Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   KPP_REAL, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   KPP_REAL, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

   KPP_REAL :: RCNTRL(20), RSTATUS(20)
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR

   INTEGER, SAVE :: Ntotal = 0
   
   !~~~> Zero input and output arrays for safety's sake
   ICNTRL     = 0
   RCNTRL     = 0.0_dp
   ISTATUS    = 0
   RSTATUS    = 0.0_dp

   !~~~> fine-tune the integrator:
   ICNTRL(1)  = 0       ! 0 - non-autonomous, 1 - autonomous
   ICNTRL(2)  = 0       ! 0 - vector tolerances, 1 - scalars
   ICNTRL(15) = 5       ! Call Update_SUN and Update_RCONST from w/in the int. 

   !~~~> if optional parameters are given, and if they are /= 0,
   !     then use them to overwrite default settings
   IF ( PRESENT( ICNTRL_U ) ) THEN
      WHERE( ICNTRL_U /= 0 ) ICNTRL = ICNTRL_U
   ENDIF
   IF ( PRESENT( RCNTRL_U ) ) THEN
      WHERE( RCNTRL_U > 0 ) RCNTRL = RCNTRL_U
   ENDIF

   !~~~> In order to remove the prior EQUIVALENCE statements (which
   !~~~> are not thread-safe), we now have declared VAR and FIX as
   !~~~> threadprivate pointer variables that can point to C.
   VAR => C(1:NVAR )
   FIX => C(NVAR+1:NSPEC)

   !~~~> Call the integrator
   CALL Exponential(NVAR,VAR,NFIX,FIX,TIN,TOUT,IERR)

   !~~~> Free pointers
   VAR => NULL()
   FIX => NULL()

   !~~~> if optional parameters are given for output
   !~~~> use them to store information in them
   IF ( PRESENT( ISTATUS_U ) ) ISTATUS_U = ISTATUS
   IF ( PRESENT( RSTATUS_U ) ) RSTATUS_U = RSTATUS
   IF ( PRESENT( IERR_U    ) ) IERR_U    = IERR
   
END SUBROUTINE INTEGRATE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Exponential(N,Y,NF,FIX,Tstart,Tend, IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using a an explicit method defined by
!    Eq. XXX in Brasseur & Jacob, 2016  
!
!    This implementation is part of KPP - the Kinetic PreProcessor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!-     Y(N)    = vector of initial conditions (at T=Tstart)
!-    [Tstart,Tend]  = time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)
!-    RelTol, AbsTol = user precribed accuracy
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_LinearAlgebra
  USE KPP_ROOT_Function
  IMPLICIT NONE

!~~~>  Arguments
   INTEGER,       INTENT(IN)    :: N, NF
   KPP_REAL, INTENT(INOUT) :: Y(N), FIX(NF)
   KPP_REAL, INTENT(IN)    :: Tstart,Tend
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Local variables
   KPP_REAL            :: Ydot(N), P(N), D(N)
!~~~>   Parameters
   KPP_REAL, PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   KPP_REAL, PARAMETER :: DeltaMin = 1.0E-5_dp

   CALL KPP_FUN_OR_FUN_SPLIT
   
   Y=Y*dexp(-1._dp*D*(Tstart-Tend))+(1._dp-dexp(-1._dp*D*(Tstart-Tend))*(P/D))

   ! Can be rewritten as a linear algebra operation for possible optimization. - MSL
   
 END SUBROUTINE Exponential
END MODULE KPP_ROOT_Integrator
