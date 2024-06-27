! Time-stamp: <2024-06-27 12:29:53 sander>

! Author:
! Rolf Sander, MPICH, Mainz, 2022

! setenv KPP_HOME ../../
! $KPP_HOME/bin/kpp mcm.kpp
! gmake -f Makefile_mcm clean
! gmake -f Makefile_mcm EXTERNAL_RATES_F90=constants_mcm.f90
! ./mcm.exe

!*****************************************************************************

MODULE sbrfunc ! subroutines and functions

  USE constants_mcm
  USE mcm_Global
  USE mcm_Rates
  USE mcm_Integrator, ONLY: integrate
  USE mcm_Initialize, ONLY: initialize
  IMPLICIT NONE
  REAL, PARAMETER :: pi = 3.14159265

CONTAINS

  SUBROUTINE print_data
    
    WRITE (10, "(F8.0,20ES11.3E2)") time, &
      C(IND_O3)/M,  C(IND_HO2)/M, H2O/M,         C(IND_CH4)/M, C(IND_NO)/M,  &
      C(IND_NO2)/M, C(IND_NO3)/M, C(IND_C5H8)/M, C(IND_OH)/M

    WRITE (20, "(F8.0,F5.0,20ES11.3E2)") time, zenith*180./pi, &
      J(1), J(2), J(3), J(4), J(5), J(6), J(7), J(8)
    
  END SUBROUTINE print_data
  
END MODULE sbrfunc

!*****************************************************************************

PROGRAM driver

  USE sbrfunc
  REAL, PARAMETER :: MAX_SZA = 89.5 *pi / 180. ! max solar zenith angle
  INTEGER  :: ICNTRL(20)
  KPP_REAL :: RCNTRL(20)

  tstart = 0.             ! time start
  tend = tstart + 86400.  ! time end
  temp = 298.             ! temperature
  dt   = 1200.            ! results every dt seconds
  zenith = MAX_SZA        ! solar zenith angle [radians] starts at midnight
  
  ICNTRL(:)  = 0
  ICNTRL(16) = 1 ! set negative concentrations to zero 
  ! ICNTRL(18) = 1 ! DAE QSSA

  RCNTRL(:)  = 0.
  RCNTRL(15) = 1.0 ! b
  RCNTRL(16) = 1.7 ! k

  ! define tolerances for the integrator:
  rtol(:) = 1e-2
  atol(:) = 1e-4

  ! open output files:
  OPEN(10, FILE="mixrat.dat", status="UNKNOWN") ! mixing ratios
  OPEN(20, FILE="jval.dat",   status="UNKNOWN") ! j-values
  
  CALL initialize() ! set up the default initial values
  M   = 2.5E19        ! convert mixing ratio to mcl/cm3
  O2  = 0.21 * M
  N2  = 0.78 * M
  H2O = 1E-2 * M
  C(IND_O3)   = 3.0E-08 * M
  C(IND_NO2)  = 1.0E-10 * M
  C(IND_CH4)  = 1.8E-06 * M
  C(IND_C5H8) = 1.0E-09 * M

  time = tstart
  WRITE (10,"(A)") &
    "time O3 HO2 H2O CH4 NO NO2 NO3 C5H8 OH"
  WRITE (20,"(A)") &
    "time sza J_O3_O1D J_O3_O3P J_H2O2 J_NO2 J_NO3_NO J_NO3_NO2 J_HONO J_HNO3"

  CALL print_data()

  ! time loop:
  DO WHILE (time .LT. tend)
    CALL update_rconst()                  ! update the rate constants
    CALL integrate(time, time+dt, ICNTRL, RCNTRL) ! call the integrator
    time = time + dt
    zenith = MIN(MAX_SZA, ABS(2.*pi*time/86400.-pi)) ! diurnal cycle (at night: 0.95*90 deg)
  CALL print_data()
  ENDDO

  ! close output files:
  CLOSE(10)
  CLOSE(20)

END PROGRAM driver
