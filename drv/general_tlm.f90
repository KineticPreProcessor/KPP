!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Driver for the tangent linear model (TLM)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM KPP_ROOT_TLM_Driver

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize

      KPP_REAL :: T, DVAL(NSPEC)
      INTEGER :: i, j, ind_1 = 1, ind_2 = 2
      ! INTEGER :: ind_1 = ind_NO2, ind_2 = ind_O3

!~~~>  NTLM = Number of sensitivity coefficients to compute
!      Note: this value is set for sensitivities w.r.t. all initial values
!            the setting may have to be changed for other applications
      INTEGER, PARAMETER :: NTLM = NVAR
      KPP_REAL, DIMENSION(NVAR,NTLM) ::  Y_tlm, ATOL_tlm, RTOL_tlm

  
!~~~>  Control (in) and status (out) vectors for the integrator
      INTEGER,  DIMENSION(20) :: ICNTRL, ISTATUS
      KPP_REAL, DIMENSION(20) :: RCNTRL, RSTATUS

      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

!~~~> Tolerances for calculating concentrations       
      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO
      
!~~~> Tolerances for calculating sensitivities 
!     are used for controlling sensitivity truncation error
!     and for solving the linear sensitivity equations by iterations  
!     Note: Sensitivities typically span many orders of magnitude
!           and a careful tuning of ATOL_tlm may be necessary     
      DO j=1,NTLM
        DO i=1,NVAR
          RTOL_tlm(i,j) = 1.0d-4
          ATOL_tlm(i,j) = 1.0d-3
        END DO
      END DO
     
      CALL Initialize()
!~~~> Note: the initial values below are for sensitivities 
!           w.r.t. initial values;
!           they may have to be changed for other applications
      Y_tlm(1:NVAR,1:NTLM) = 0.0d0
      DO j=1,NTLM
        Y_tlm(j,j) = 1.0d0
      END DO

!~~~> Default control options
      ICNTRL(1:20) = 0
      RCNTRL(1:20) = 0.0d0       

!~~~> Begin time loop

      CALL InitSaveData()

      T = TSTART
      
kron: DO WHILE (T < TEND)
       
        TIME = T
        CALL GetMass( C, DVAL )
        WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,      &
                    (TRIM(SPC_NAMES(MONITOR(i))),          &
                      C(MONITOR(i))/CFACTOR, i=1,NMONITOR),&
                    (TRIM(SMASS(i)),DVAL(i)/CFACTOR, i=1,NMASS)

        CALL SaveData()
        CALL Update_SUN() 
        CALL Update_RCONST()

        CALL INTEGRATE_TLM( NTLM, VAR, Y_tlm, T, T+DT,     &
                            ATOL_tlm, RTOL_tlm,            &
                            ICNTRL, RCNTRL, ISTATUS, RSTATUS )

        T = T+DT

      END DO kron

      CALL GetMass( C, DVAL )
      WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,        &
                  (TRIM(SPC_NAMES(MONITOR(i))),            &
                    C(MONITOR(i))/CFACTOR, i=1,NMONITOR),  &
                  (TRIM(SMASS(i)),DVAL(i)/CFACTOR, i=1,NMASS)

      CALL SaveData()

!~~~> End time loop ~~~~~~~~~~

      OPEN(20, FILE='KPP_ROOT_TLM_results.m')
      WRITE(6,*) '**************************************************'
      WRITE(6,*) ' Concentrations and Sensitivities at final time '
      WRITE(6,*) ' were written in the file KPP_ROOT_TLM_results.m'
      WRITE(6,*) '**************************************************'
      DO j=1,NTLM
        WRITE(20,993) ( Y_tlm(i,j), i=1,NVAR )          
      END DO

      CALL CloseSaveData()

      WRITE(6,995) TRIM(SPC_NAMES(ind_1)),TRIM(SPC_NAMES(ind_1)), &
                   Y_tlm(ind_1,1)
      WRITE(6,995) TRIM(SPC_NAMES(ind_2)),TRIM(SPC_NAMES(ind_2)), &
                   Y_tlm(ind_2,2)
      WRITE(6,995) TRIM(SPC_NAMES(ind_1)),TRIM(SPC_NAMES(ind_2)), &
                   Y_tlm(ind_2,1)
      WRITE(6,995) TRIM(SPC_NAMES(ind_2)),TRIM(SPC_NAMES(ind_1)), &
                   Y_tlm(ind_1,2)
                  
                 
 991  FORMAT(F6.1,'%. T=',E10.3,3X,20(A,'=',E10.4,';',1X))
 993  FORMAT(1000(E24.16,2X))
 995  FORMAT('TLM: d[',A,'](tf)/d[',A,'](t0)=',E14.7)
 
 
      !~~~> The entire matrix of sensitivities
      WRITE(6,996) ( 'd ',TRIM(SPC_NAMES(i)), i=1,NVAR ) 
      DO j=1,NTLM
        WRITE(6,997) TRIM(SPC_NAMES(j)),( Y_tlm(i,j), i=1,NVAR )          
      END DO
 996  FORMAT(12X,100('  ',A2,A6,4X))
 997  FORMAT('d/d',A6,' = ',100(E12.5,2X))
 

END PROGRAM KPP_ROOT_TLM_Driver

