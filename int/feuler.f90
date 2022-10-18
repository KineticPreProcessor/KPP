!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!    Forward Euler method for non-stiff ODEs                              !
!                                                                         !
!    Version History:                                                     !
!    Apr. 30, 2022 - Initial Version - M.S.Long                           !
!                                                                         !
!    This implementation is part of KPP - the Kinetic PreProcessor        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
MODULE KPP_ROOT_Integrator

  USE KPP_ROOT_Precision,  ONLY: dp
  USE KPP_ROOT_Parameters, ONLY: NVAR, NFIX, NSPEC
  USE KPP_ROOT_Global
  USE KPP_ROOT_Monitor

  IMPLICIT NONE
  PUBLIC

  SAVE

!~~~> Flags to determine if we should call the UPDATE_* routines from within
!~~~> the integrator.  If using KPP in an external model, you might want to
!~~~> disable these calls (via ICNTRL(15)) to avoid excess computations.
  LOGICAL, PRIVATE :: Do_Update_RCONST
  LOGICAL, PRIVATE :: Do_Update_PHOTO
  LOGICAL, PRIVATE :: Do_Update_SUN

!~~~> description of the error numbers IERR
  CHARACTER(LEN=50), PARAMETER, DIMENSION(1) :: IERR_NAMES = (/ &
    'dummy value                                       ' /)

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE Integrate( TIN,       TOUT,      ICNTRL_U , RCNTRL_U,           &
                        ISTATUS_U, RSTATUS_U, IERR_U                        )
    ! ICNTRL(16)
    ! 0 -> do nothing.
    ! 1 -> set negative values to zero
    ! 2 -> return with error code
    ! 3 -> stop at negative
    !
    ! ICNTRL(17) = verbose error output

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! User interface routine to the KPP forward Euler integrator
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    USE KPP_ROOT_Util, ONLY : Integrator_Update_Options

    !~~~> Inputs
    KPP_REAL, INTENT(IN)            :: TIN           ! Start Time
    KPP_REAL, INTENT(IN)            :: TOUT          ! End Time
    INTEGER,  INTENT(IN),  OPTIONAL :: ICNTRL_U(20)  ! Input options
    KPP_REAL, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)  ! Input options

    !~~~> Outputs
    INTEGER,  INTENT(IN),  OPTIONAL :: ISTATUS_U(20) ! Returned status values
    KPP_REAL, INTENT(OUT), OPTIONAL :: RSTATUS_U(20) ! Returned status values
    INTEGER,  INTENT(OUT), OPTIONAL :: IERR_U        ! Error code

    ! Local variables
    INTEGER                         :: ICNTRL(20)
    KPP_REAL                        :: RSTATUS(20)
    INTEGER                         :: IERR

    !~~~> Zero input and output arrays for safety's sake
    ICNTRL     = 0
    RSTATUS    = 0.0_dp

    !~~~> fine-tune the integrator:
    ICNTRL(1 ) = 0       ! Verbose error output
    ICNTRL(2 ) = 0       ! Stop upon negative results
    ICNTRL(15) = 5       ! Call Update_SUN and Update_RCONST from w/in the int.

   !~~~> if optional parameters are given, and if they are /= 0,
   !     then use them to overwrite default settings
   IF ( PRESENT( ICNTRL_U ) ) THEN
      WHERE( ICNTRL_U /= 0 ) ICNTRL = ICNTRL_U
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
   VAR => C(1:NVAR)
   FIX => C(NVAR+1:NSPEC)

   !~~~> Call the integrator
   CALL ForwardEuler( NVAR, C(1:NVAR), TIN, TOUT, ICNTRL, IERR )

   !~~~> Free pointers
   VAR => NULL()
   FIX => NULL()

   !~~~> Return error status (NOTE: ISTATUS_U does nothing)
   IF ( PRESENT( IERR_U    ) ) IERR_U    = IERR
   IF ( PRESENT( RSTATUS_U ) ) RSTATUS_U = RSTATUS

  END SUBROUTINE Integrate

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE ForwardEuler( N, Y, Tstart, Tend, ICNTRL, IERR )
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !  Forward Euler integrator
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! External functions
    USE KPP_ROOT_Rates, ONLY : Update_SUN, Update_RCONST

    ! Input: Arguments
    INTEGER,  INTENT(IN)    :: N            ! Dimension for Y(N)
    KPP_REAL, INTENT(INOUT) :: Y(N)         ! Initial condition
    KPP_REAL, INTENT(IN)    :: Tstart       ! Starting time
    KPP_REAL, INTENT(IN)    :: Tend         ! Ending time
    INTEGER,  INTENT(IN)    :: ICNTRL(20)   ! KPP integrator options

    ! Output:
    INTEGER, INTENT(OUT)    :: IERR         ! Error status

    ! Local variables
    INTEGER  :: i
    LOGICAL  :: hasNegative
    KPP_REAL :: Ynew(N), dYdt(N)

    ! Update rates before integration if desired
    IF ( Do_Update_SUN    ) CALL Update_SUN()
    IF ( Do_Update_RCONST ) CALL Update_RCONST()

    ! Get P-kY
    CALL FunTemplate( Y, dYdt )

    !~~~> Do the integration
    Ynew = Y + ( dYdt * ( Tend - Tstart ) )

    ! Check for negatives
    IF (ICNTRL(16) .gt. 0) THEN ! Don't perform DO() loop if you don't care
       DO i=1,N
          IF (Ynew(i) .lt. 0._dp) THEN
             IF (ICNTRL(17) /= 0) THEN
                write(*,*) trim(SPC_NAMES(i)), " is negative: ", Ynew(i)
             ENDIF
             IERR = -9
             IF (ICNTRL(16) == 1) THEN
                Ynew(i) = 0._dp
             ELSE IF (ICNTRL(16) == 2) THEN
                write(*,*) '(ICNTRL(16) = 2) Negative value. Returning.'
                RETURN
             ELSE IF (ICNTRL(16) == 3) THEN
                write(*,*) '(ICNTRL(16) = 3) Negative value. Stopping.'
                STOP
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    ! Return updated concentrations
    Y = Ynew

    ! Succesful exit
    IERR = 1

  END SUBROUTINE ForwardEuler

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE FunTemplate( Y, Ydot )
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !  Template for the ODE function call.
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USE KPP_ROOT_Parameters, ONLY : NVAR
    USE KPP_ROOT_Global,     ONLY : FIX, RCONST
    USE KPP_ROOT_Function,   ONLY : Fun

    !~~~> Input variables
    KPP_REAL :: Y(NVAR)

    !~~~> Output variables
    KPP_REAL :: Ydot(NVAR)

    !~~~> Compute equation rates and time derivative of variable species
    CALL Fun( Y, FIX, RCONST, Ydot )

  END SUBROUTINE FunTemplate

END MODULE KPP_ROOT_Integrator
