    MODULE EXAMPLE1

! DEMONSTRATION PROGRAM FOR THE DVODE_F90 PACKAGE.

! The following is a simple example problem, with the coding
! needed for its solution by DVODE_F90. The problem is from
! chemical kinetics, and consists of the following three rate
! equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.

! The following coding solves this problem with DVODE_F90,
! using a user supplied Jacobian and printing results at
! t = .4, 4.,...,4.d10. It uses ITOL = 2 and ATOL much smaller
! for y2 than y1 or y3 because y2 has much smaller values. At
! the end of the run, statistical quantities of interest are
! printed. (See optional output in the full DVODE description
! below.) Output is written to the file example1.dat.

    CONTAINS

      SUBROUTINE FEX(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)

        YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
        YDOT(3) = 3.E7*Y(2)*Y(2)
        YDOT(2) = -YDOT(1) - YDOT(3)
        RETURN
      END SUBROUTINE FEX

      SUBROUTINE JEX(NEQ,T,Y,ML,MU,PD,NRPD)
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NRPD
        DOUBLE PRECISION PD, T, Y
        DIMENSION Y(NEQ), PD(NRPD,NEQ)

        PD(1,1) = -.04D0
        PD(1,2) = 1.D4*Y(3)
        PD(1,3) = 1.D4*Y(2)
        PD(2,1) = .04D0
        PD(2,3) = -PD(1,3)
        PD(3,2) = 6.E7*Y(2)
        PD(2,2) = -PD(1,2) - PD(3,2)
        RETURN
      END SUBROUTINE JEX

    END MODULE EXAMPLE1

!******************************************************************

    PROGRAM RUNEXAMPLE1

      USE DVODE_F90_M
      USE EXAMPLE1

      IMPLICIT NONE
      DOUBLE PRECISION ATOL, RTOL, T, TOUT, Y, RSTATS
      INTEGER NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, I
      DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31)

      TYPE (VODE_OPTS) :: OPTIONS

      OPEN (UNIT=6,FILE='example1.dat')
      IERROR = 0
      NEQ = 3
      Y(1) = 1.0D0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      T = 0.0D0
      TOUT = 0.4D0
      RTOL = 1.D-4
      ATOL(1) = 1.D-8
      ATOL(2) = 1.D-14
      ATOL(3) = 1.D-6
      ITASK = 1
      ISTATE = 1
!     OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR=RTOL, &
!       USER_SUPPLIED_JACOBIAN=.TRUE.)
      OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,      &
        RELERR=RTOL,USER_SUPPLIED_JACOBIAN=.TRUE.)
      DO IOUT = 1, 12
        CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEX)
        CALL GET_STATS(RSTATS,ISTATS)
        WRITE (6,90003) T, Y(1), Y(2), Y(3)
        DO I = 1, NEQ
          IF (Y(I)<0.0D0) IERROR = 1
        END DO
        IF (ISTATE<0) THEN
          WRITE (6,90004) ISTATE
          STOP
        END IF
        TOUT = TOUT*10.0D0
      END DO
      WRITE (6,90000) ISTATS(11), ISTATS(12), ISTATS(13), ISTATS(19), &
        ISTATS(20), ISTATS(21), ISTATS(22)
      IF (IERROR==1) THEN
        WRITE (6,90001)
      ELSE
        WRITE (6,90002)
      END IF
90000 FORMAT (/'  No. steps =',I4,'   No. f-s =',I4,'  No. J-s =',I4, &
        '   No. LU-s =',I4/'  No. nonlinear iterations =', &
        I4/'  No. nonlinear convergence failures =', &
        I4/'  No. error test failures =',I4/)
90001 FORMAT (/' An error occurred.')
90002 FORMAT (/' No errors occurred.')
90003 FORMAT (' At t =',D12.4,'   y =',3D14.6)
90004 FORMAT (///' Error halt: ISTATE =',I3)
      STOP
    END PROGRAM RUNEXAMPLE1
