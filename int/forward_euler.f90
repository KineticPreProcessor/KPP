!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!    Forward Euler method for non-stiff ODEs                              !
!                                                                         !
!    Version History:                                                     !
!    Apr. 30, 2022 - Initial Version - M.S.Long                           !
!                                                                         !
!    This implementation is part of KPP - the Kinetic PreProcessor        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
MODULE KPP_ROOT_Integrator

  USE KPP_ROOT_Precision, ONLY: dp
  USE KPP_ROOT_Parameters, ONLY: NVAR, NFIX, NSPEC
  USE KPP_ROOT_Global
  IMPLICIT NONE
  PUBLIC
  SAVE

  ! description of the error numbers IERR
  CHARACTER(LEN=50), PARAMETER, DIMENSION(1) :: IERR_NAMES = (/ &
    'dummy value                                       ' /)

CONTAINS

  ! **************************************************************************

  SUBROUTINE INTEGRATE( TIN, TOUT )
    IMPLICIT NONE
    
    KPP_REAL, INTENT(IN) :: TIN  ! Start Time
    KPP_REAL, INTENT(IN) :: TOUT ! End Time

    CALL FORWARDEULER(NVAR,VAR,TIN,TOUT)
  END SUBROUTINE INTEGRATE

  SUBROUTINE FORWARDEULER(N,Y,Tstart,Tend)


    ! Get P-kY
    CALL FunTEMPLATE( Y, dYdt )

    ! Integrate
    Ynew = Y+dYdt*(Tend-Tstart)

    ! Check for negatives

  END SUBROUTINE FORWARDEULER

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE FunTemplate( Y, Ydot )
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !  Template for the ODE function call.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USE KPP_ROOT_Parameters, ONLY: NVAR
    USE KPP_ROOT_Global,     ONLY: FIX, RCONST
    USE KPP_ROOT_Function,   ONLY: Fun
    
    !~~~> Input variables
    KPP_REAL :: Y(NVAR)
    !~~~> Output variables
    KPP_REAL :: Ydot(NVAR)

    CALL Fun( Y, FIX, RCONST, Ydot )

  END SUBROUTINE FunTemplate

END MODULE KPP_ROOT_Integrator
