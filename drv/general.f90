PROGRAM KPP_ROOT_Driver

!~~~> Global variables and routines
      USE KPP_ROOT_Model
      USE KPP_ROOT_Initialize, ONLY: Initialize

!~~~> Local variables
      KPP_REAL :: T, DVAL(NSPEC)
      KPP_REAL :: RSTATE(20)
      INTEGER :: i

!~~~> Control (in) arguments for the integration
!~~~> These are set to zero by default, which will invoke default behavior
      INTEGER :: ICNTRL(20)
      KPP_REAL :: RCNTRL(20)

      ICNTRL  = 0
      RCNTRL  = 0.d0

!~~~> Initialization
      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO

!~~~> Set default species concentrations
      CALL Initialize()

!~~~> Open log file for write
      CALL InitSaveData()

!~~~> Time loop
      T = TSTART
kron: DO WHILE (T < TEND)

        TIME = T

!~~~> Write values of monitored species at each iteration
        CALL GetMass( C, DVAL )
        WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,       &
                   ( TRIM(SPC_NAMES(MONITOR(i))),           &
                     C(MONITOR(i))/CFACTOR, i=1,NMONITOR )
        CALL SaveData()

!~~~> Update sunlight intensity and reaction rates
!~~~> before calling the integrator.
        CALL Update_SUN()
        CALL Update_RCONST()

!~~~> Call the integrator
        CALL INTEGRATE( TIN       = T,        &
                        TOUT      = T+DT,     &
                        ICNTRL_U  = ICNTRL,   &
                        RCNTRL_U  = RCNTRL,   &
                        RSTATUS_U = RSTATE   )
        T = RSTATE(1)

      END DO kron
!~~~> End Time loop

!~~~> Write final values of monitored species and close log file
      CALL GetMass( C, DVAL )
      WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,     &
               ( TRIM(SPC_NAMES(MONITOR(i))),           &
                 C(MONITOR(i))/CFACTOR, i=1,NMONITOR )
      TIME = T
      CALL SaveData()
      CALL CloseSaveData()

991   FORMAT(F6.1,'%. T=',E9.3,2X,200(A,'=',E11.4,'; '))

END PROGRAM KPP_ROOT_Driver
