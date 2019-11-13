    MODULE EXAMPLE2

! DEMONSTRATION PROGRAM FOR THE DVODE_F90 PACKAGE.

! The following is a modification of the previous example
! program to illustrate root finding. The problem is from
! chemical kinetics, and consists of the following three
! rate equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
! In addition, we want to find the values of t, y1, y2,
! and y3 at which:
!   (1) y1 reaches the value 1.d-4, and
!   (2) y3 reaches the value 1.d-2.

! The following coding solves this problem with DVODE_F90
! using an internally generated dense Jacobian and
! printing results at t = .4, 4., ..., 4.d10, and at the
! computed roots. It uses ITOL = 2 and ATOL much smaller
! for y2 than y1 or y3 because y2 has much smaller values.
! At the end of the run, statistical quantities of interest
! are printed (see optional outputs in the full description
! below). Output is written to the file example2.dat.

    CONTAINS

      SUBROUTINE FEX(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(3), YDOT(3)

        YDOT(1) = -0.04D0*Y(1) + 1.0D4*Y(2)*Y(3)
        YDOT(3) = 3.0D7*Y(2)*Y(2)
        YDOT(2) = -YDOT(1) - YDOT(3)
        RETURN
      END SUBROUTINE FEX

      SUBROUTINE GEX(NEQ,T,Y,NG,GOUT)
        IMPLICIT NONE
        INTEGER NEQ, NG
        DOUBLE PRECISION T, Y, GOUT
        DIMENSION Y(3), GOUT(2)

        GOUT(1) = Y(1) - 1.0D-4
        GOUT(2) = Y(3) - 1.0D-2
        RETURN
      END SUBROUTINE GEX

    END MODULE EXAMPLE2

!******************************************************************

    PROGRAM RUNEXAMPLE2

      USE DVODE_F90_M
      USE EXAMPLE2

      IMPLICIT NONE
      INTEGER ITASK, ISTATE, NG, NEQ, IOUT, JROOT, ISTATS, IERROR, I
      DOUBLE PRECISION ATOL, RTOL, RSTATS, T, TOUT, Y
      DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31), JROOT(2)

      TYPE (VODE_OPTS) :: OPTIONS

      OPEN (UNIT=6,FILE='example2.dat')
      IERROR = 0
      NEQ = 3
      Y(1) = 1.0D0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      T = 0.0D0
      TOUT = 0.4D0
      RTOL = 1.0D-4
      ATOL(1) = 1.0D-8
      ATOL(2) = 1.0D-10
      ATOL(3) = 1.0D-8
      ITASK = 1
      ISTATE = 1
      NG = 2
!     OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR_VECTOR=ATOL, &
!       NEVENTS=NG)
      OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,             &
        ABSERR_VECTOR=ATOL,NEVENTS=NG)
      DO 20 IOUT = 1, 12
10      CONTINUE
        CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=GEX)
        CALL GET_STATS(RSTATS,ISTATS,NG,JROOT)
        WRITE (6,90000) T, Y(1), Y(2), Y(3)
        DO I = 1, NEQ
          IF (Y(I)<0.0D0) IERROR = 1
        END DO
90000   FORMAT (' At t =',D12.4,'   Y =',3D14.6)
        IF (ISTATE<0) GO TO 30
        IF (ISTATE==2) GO TO 20
        WRITE (6,90001) JROOT(1), JROOT(2)
90001   FORMAT (5X,' The above line is a root, JROOT =',2I5)
        ISTATE = 2
        GO TO 10
20    TOUT = TOUT*10.0D0
      WRITE (6,90002) ISTATS(11), ISTATS(12), ISTATS(13), ISTATS(10)
      IF (IERROR==1) THEN
        WRITE (6,90003)
      ELSE
        WRITE (6,90004)
      END IF
90002 FORMAT (/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4, &
        '  No. g-s =',I4/)
      STOP
30    WRITE (6,90005) ISTATE
90003 FORMAT (/' An error occurred.')
90004 FORMAT (/' No errors occurred.')
90005 FORMAT (///' Error halt.. ISTATE =',I3)
      STOP
    END PROGRAM RUNEXAMPLE2
