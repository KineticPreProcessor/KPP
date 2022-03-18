!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Driver for the Adjoint (ADJ) model
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM KPP_ROOT_ADJ_Driver

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize

      KPP_REAL :: T, DVAL(NSPEC)
      INTEGER :: i, j, ind_1 = 1, ind_2 = 2
      ! INTEGER :: ind_1 = ind_NO2, ind_2 = ind_O3

!~~>  NADJ = Number of functionals for which sensitivities are computed
!     Note: the setting below is for sensitivities of all final concentrations;
!           the setting may have to be changed for other applications
      INTEGER, PARAMETER :: NADJ = NVAR
      KPP_REAL, DIMENSION(NVAR,NADJ) :: Y_adj, ATOL_adj, RTOL_adj

!~~>  Control (in) and status (out) arguments for the integration
      KPP_REAL, DIMENSION(20) :: RCNTRL, RSTATUS
      INTEGER,       DIMENSION(20) :: ICNTRL, ISTATUS
  
      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

!~~~> Tolerances for calculating concentrations       
      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO
      
!~~~> Tolerances for calculating adjoints 
!     are used for controlling adjoint truncation error
!     and for solving the linear adjoint equations by iterations  
!     Note: Adjoints typically span many orders of magnitude
!           and a careful tuning of ATOL_adj may be necessary     
      DO j=1,NADJ
        DO i=1,NVAR
          RTOL_adj(i,j) = 1.0d-4
          ATOL_adj(i,j) = 1.0d-10
        END DO
      END DO
     
      CALL Initialize()
      
!~~~>  The adjoint values at the final time
      Y_adj(1:NVAR,1:NADJ) = 0.0d0
      DO j=1,NADJ
        Y_adj(j,j) = 1.0d0
      END DO

!~~~> Default control options
      ICNTRL(1:20) = 0
      RCNTRL(1:20) = 0.0d0       

!~~~> Begin time loop

      TIME = TSTART
      CALL InitSaveData()

      T = TSTART

      CALL GetMass( C, DVAL )
      WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,        &
                  (TRIM(SPC_NAMES(MONITOR(i))),            &
                    C(MONITOR(i))/CFACTOR, i=1,NMONITOR),  &
                  (TRIM(SMASS(i)),DVAL(i)/CFACTOR, i=1,NMASS)

      TIME = T
      CALL SaveData()

      CALL INTEGRATE_ADJ( NADJ, VAR, Y_adj, T, TEND,       &
                          ATOL_adj, RTOL_adj,              &
                          ICNTRL, RCNTRL, ISTATUS, RSTATUS )


      CALL GetMass( C, DVAL )
      WRITE(6,991) (TEND-TSTART)/(TEND-TSTART)*100, TEND,  &
                  (TRIM(SPC_NAMES(MONITOR(i))),            &
                    C(MONITOR(i))/CFACTOR, i=1,NMONITOR),  &
                  (TRIM(SMASS(i)),DVAL(i)/CFACTOR, i=1,NMASS)

      TIME = T
      CALL SaveData()

!~~~> End time loop ~~~~~~~~~~

      OPEN(20, FILE='KPP_ROOT_ADJ_results.m')
      WRITE(6,*) '**************************************************'
      WRITE(6,*) ' Concentrations and Sensitivities at final time '
      WRITE(6,*) ' were written in the file KPP_ROOT_ADJ_results.m'
      WRITE(6,*) '**************************************************'
      DO j=1,NADJ
        WRITE(20,993) ( Y_adj(i,j), i=1,NVAR )          
      END DO

      WRITE(6,995) TRIM(SPC_NAMES(ind_1)),TRIM(SPC_NAMES(ind_1)), &
                   Y_adj(ind_1,1)
      WRITE(6,995) TRIM(SPC_NAMES(ind_2)),TRIM(SPC_NAMES(ind_2)), &
                   Y_adj(ind_2,2)
      WRITE(6,995) TRIM(SPC_NAMES(ind_2)),TRIM(SPC_NAMES(ind_1)), &
                   Y_adj(ind_1,2)
      WRITE(6,995) TRIM(SPC_NAMES(ind_1)),TRIM(SPC_NAMES(ind_2)), &
                   Y_adj(ind_2,1)

      CALL CloseSaveData()
      
 991  FORMAT(F6.1,'%. T=',E10.3,3X,20(A,'=',E10.4,';',1X))
 993  FORMAT(1000(E24.16,2X))
 995  FORMAT('ADJ: d[',A,'](tf)/d[',A,'](t0)=',E14.7)
 
      !~~~> The entire matrix of sensitivities
      WRITE(6,996) ( 'd ',TRIM(SPC_NAMES(i)), i=1,NVAR ) 
      DO j=1,NADJ
        WRITE(6,997) TRIM(SPC_NAMES(j)),( Y_adj(j,i), i=1,NVAR )          
      END DO
 996  FORMAT(12X,100('  ',A2,A6,4X))
 997  FORMAT('d/d',A6,' = ',100(E12.5,2X))

END PROGRAM KPP_ROOT_ADJ_Driver

