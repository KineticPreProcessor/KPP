! DVODE_F90 demonstration program
! Enforces nonnegativity on the Class F problems

! Toronto Stiff Test Set

! This program solves each of the thirty problems for several error
! tolerances using each of the dense, banded, and sparse solution
! options and with the numerical Jacobians formed using each of the
! original VODE algorithm and Doug Salane's JACSP. The ODEs are
! solved both in scaled and unscaled form. Scalar error tolerances
! are used for all problems even though this is not appropriate in
! some cases. The program takes 10-15 minutes to run because several
! thousand integrations are performed in all. Details of each
! integration are written to the output file 'stiffoptions.dat'.
! Summaries of the integration statistics may be found by searching
! on the string 'Individual'. There are several blocks of such summaries.
! When you run the program, the last line in the output file should
! read: 'No errors occurred.' Please contact one of the authors if
! any errors do occur.

! Reference:
! ALGORITHM 648, COLLECTED ALGORITHMS FROM ACM.
! TRANSACTIONS ON MATHEMATICAL SOFTWARE,
! VOL. 13, NO. 1, P. 28

! Caution:
! The original test routines are altered in this version.
! Some subroutine arguments, COMMON variables, and local
! variables have been moved to the following global header.
! You should obtain the original test suite from netlib
! if you plan to use the routines for other purposes.

    MODULE STIFFSET

! Note:
! In this version ID and IID are defined in the main
! program demostiff (not in the Toronto routines).
! The other parameters are obtained by calling IVALU
! with a modified argument list.

      IMPLICIT NONE
      INTEGER IWT, N, ID, IID, NFCN, NJAC, NLUD
      DOUBLE PRECISION W, DY
      DIMENSION W(20), DY(400)

    CONTAINS

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
!       Define the system derivatives.
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)
        CALL FCN(T,Y,YDOT)
        RETURN
      END SUBROUTINE DERIVS

      SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!       Load the Jacobian as a dense matrix.
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD, I, J
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(NEQ), PD(NROWPD,NEQ)
        CALL PDERV(T,Y)
        DO J = 1, NEQ
          DO I = 1, NEQ
            PD(I,J) = DY(I+(J-1)*NROWPD)
          END DO
        END DO
        RETURN
      END SUBROUTINE JACD

      SUBROUTINE JACB(NEQ,T,Y,ML,MU,PD,NROWPD)
!       Load the Jacobian as a banded matrix.
        IMPLICIT NONE
        INTEGER NEQ, ML, MU, NROWPD, I, J, K
        DOUBLE PRECISION T, Y, PD
        DIMENSION Y(NEQ), PD(NROWPD,NEQ)
        CALL PDERV(T,Y)
        DO J = 1, NEQ
           DO I = 1, NEQ
              K = I - J + MU + 1
              PD(K,J) = DY(I+(J-1)*NEQ)
           END DO
        END DO
      END SUBROUTINE JACB

      SUBROUTINE JACS(NEQ,T,Y,IA,JA,NZ,P)
!       Load the Jacobian as a sparse matrix.
        IMPLICIT NONE
        INTEGER NEQ, IA, JA, NZ, ML, MU, COL, ROW, I, NROWPD
!       INTEGER NZSAVE
        DOUBLE PRECISION T, Y, P
        DIMENSION Y(*), IA(*), JA(*), P(*)
        IF (NZ<=0) THEN
          NZ = NEQ*NEQ
          RETURN
        END IF
        ML = NEQ - 1
        MU = NEQ - 1
        NROWPD = NEQ
        CALL JACD(NEQ,T,Y,ML,MU,P,NROWPD)
        IA(1) = 1
        DO I = 1, NEQ
          IA(I+1) = IA(I) + NEQ
        END DO
        I = 0
        DO COL = 1, NEQ
          I = NEQ*(COL-1)
          DO ROW = 1, NEQ
            I = I + 1
            JA(I) = ROW
          END DO
        END DO
        RETURN
      END SUBROUTINE JACS

      SUBROUTINE IVALU(XSTART,XEND,HBEGIN,HMAX,Y)

!      ROUTINE TO PROVIDE THE INITIAL VALUES REQUIRED TO SPECIFY
!      THE MATHEMATICAL PROBLEM AS WELL AS VARIOUS PROBLEM
!      PARAMETERS REQUIRED BY THE TESTING PACKAGE. THE APPROPRIATE
!      SCALING VECTOR IS ALSO INITIALISED IN CASE THIS OPTION IS
!      SELECTED.

!      PARAMETERS (OUTPUT)
!         N      - DIMENSION OF THE PROBLEM
!         XSTART - INITIAL VALUE OF THE INDEPENDENT VARIABLE
!         XEND   - FINAL VALUE OF THE INDEPENDENT VARIABLE
!         HBEGIN - APPROPRIATE STARTING STEPSIZE
!         Y      - VECTOR OF INITIAL CONDITIONS FOR THE DEPENDENT
!                  VARIABLES
!         W      - VECTOR OF WEIGHTS USED TO SCALE THE PROBLEM IF
!                  THIS OPTION IS SELECTED.

!      PARAMETER  (INPUT)
!         IWT    - FLAG TO INDICATE IF SCALED OPTION IS SELECTED
!         ID     - FLAG IDENTIFYING WHICH EQUATION IS BEING SOLVED

        IMPLICIT NONE

!     .. Scalar Arguments ..
        DOUBLE PRECISION HBEGIN, HMAX, XEND, XSTART
!     .. Array Arguments ..
        DOUBLE PRECISION Y(20)
!     .. Local Scalars ..
        DOUBLE PRECISION XS
        INTEGER I, IOUT, ITMP
!     .. Data statements ..
        DATA XS/0.D0/

!     .. Executable Statements ..

        XSTART = XS
!     GOTO (40, 80, 120, 160, 20, 20, 20, 20, 20, 20, 200, 220, 220,    &
!     220, 220, 20, 20, 20, 20, 20, 360, 400, 400, 400, 400, 20, 20, 20,&
!     20, 20, 540, 580, 600, 640, 660, 680, 20, 20, 20, 20, 700, 740,   &
!     760, 780, 800, 20, 20, 20, 20, 20, 840, 860, 880, 900, 920) ID
        GOTO (20,30,40,50,10,10,10,10,10,10,60,70,70,70,70,10,10,10,10,10, &
          130,140,140,140,140,10,10,10,10,10,200,210,220,230,240,250,10,10,10, &
          10,260,270,280,290,300,10,10,10,10,10,310,320,330,340,350) ID
10      IOUT = 6
        WRITE (IOUT,FMT=90000) ID
        STOP

!     PROBLEM CLASS A - LINEAR WITH REAL EIGENVALUES

20      CONTINUE
!     PROBLEM A1
        N = 4
        W(1) = 0.100D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        XEND = 20.D0
        HBEGIN = 1.0D-2
        HMAX = 20.D0
        DO I = 1, N
          Y(I) = 1.0D0
        END DO
        GOTO 360

30      CONTINUE
!     PROBLEM A2
        N = 9
        W(1) = 0.100D+00
        W(2) = 0.200D+00
        W(3) = 0.300D+00
        W(4) = 0.400D+00
        W(5) = 0.500D+00
        W(6) = 0.600D+00
        W(7) = 0.700D+00
        W(8) = 0.800D+00
        W(9) = 0.900D+00
        XEND = 120.D0
        HBEGIN = 5.D-4
        HMAX = 120.D0
        DO I = 1, N
          Y(I) = 0.D0
        END DO
        GOTO 360

40      CONTINUE
!     PROBLEM A3
        N = 4
        W(1) = 0.100D+01
        W(2) = 0.100D+01
        W(3) = 0.782D+01
        W(4) = 0.100D+01
        HBEGIN = 1.D-5
        XEND = 20.D0
        HMAX = 20.D0
        DO I = 1, N
          Y(I) = 1.D0
        END DO
        GOTO 360

50      CONTINUE
!     PROBLEM A4
        N = 10
        W(1) = 0.100D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        W(5) = 0.100D+01
        W(6) = 0.100D+01
        W(7) = 0.100D+01
        W(8) = 0.100D+01
        W(9) = 0.100D+01
        W(10) = 0.100D+01
        XEND = 1.D0
        HBEGIN = 1.D-5
        HMAX = 1.D0
        DO I = 1, N
          Y(I) = 1.D0
        END DO
        GOTO 360

!     PROBLEM CLASS B - LINEAR WITH NON-REAL EIGENVALUES

60      CONTINUE
!     PROBLEM B1
        N = 4
        W(1) = 0.100D+01
        W(2) = 0.859D+01
        W(3) = 0.100D+01
        W(4) = 0.322D+02
        XEND = 20.D0
        HBEGIN = 7.D-3
        HMAX = 20.D0
        Y(1) = 1.D0
        Y(2) = 0.D0
        Y(3) = 1.D0
        Y(4) = 0.D0
        GOTO 360

70      CONTINUE
!     PROBLEM B2, B3, B4, B5
        N = 6
        ITMP = IID - 1
        GOTO (80,90,100,110) ITMP
80      CONTINUE
        W(1) = 0.100D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        W(5) = 0.100D+01
        W(6) = 0.100D+01
        GOTO 120
90      CONTINUE
        W(1) = 0.100D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        W(5) = 0.100D+01
        W(6) = 0.100D+01
        GOTO 120
100     CONTINUE
        W(1) = 0.112D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        W(5) = 0.100D+01
        W(6) = 0.100D+01
        GOTO 120
110     CONTINUE
        W(1) = 0.131D+01
        W(2) = 0.112D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        W(5) = 0.100D+01
        W(6) = 0.100D+01
120     CONTINUE
        XEND = 20.D0
        HBEGIN = 1.D-2
        HMAX = 20.D0
        DO I = 1, N
          Y(I) = 1.D0
        END DO
        GOTO 360

!     PROBLEM CLASS C - NON-LINEAR COUPLING FROM
!                       STEADY STATE TO TRANSIENT

130     CONTINUE
!     PROBLEM C1
        N = 4
        W(1) = 0.102D+01
        W(2) = 0.103D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        XEND = 20.D0
        HBEGIN = 1.D-2
        HMAX = 20.D0
        DO I = 1, N
          Y(I) = 1.D0
        END DO
        GOTO 360

140     CONTINUE
!     PROBLEM C2, C3, C4, C5
        N = 4
        ITMP = IID - 1
        GOTO (150,160,170,180) ITMP
150     CONTINUE
        W(1) = 0.200D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        GOTO 190
160     CONTINUE
        W(1) = 0.200D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        W(4) = 0.100D+01
        GOTO 190
170     CONTINUE
        W(1) = 0.200D+01
        W(2) = 0.400D+01
        W(3) = 0.200D+02
        W(4) = 0.420D+03
        GOTO 190
180     CONTINUE
        W(1) = 0.200D+01
        W(2) = 0.800D+01
        W(3) = 0.136D+03
        W(4) = 0.371D+05
190     CONTINUE
        XEND = 20.D0
        HBEGIN = 1.D-2
        HMAX = 20.D0
        DO I = 1, N
          Y(I) = 1.D0
        END DO
        GOTO 360

!     PROBLEM CLASS D - NON-LINEAR WITH REAL EIGENVALUES

200     CONTINUE
!     PROBLEM D1
        N = 3
        W(1) = 0.223D+02
        W(2) = 0.271D+02
        W(3) = 0.400D+03
        XEND = 400.D0
        HBEGIN = 1.7D-2
        HMAX = 400.D0
        DO I = 1, N
          Y(I) = 0.D0
        END DO
        GOTO 360

210     CONTINUE
!     PROBLEM D2
        N = 3
        W(1) = 0.100D+01
        W(2) = 0.365D+00
        W(3) = 0.285D+02
        XEND = 40.D0
        HBEGIN = 1.D-5
        HMAX = 40.D0
        Y(1) = 1.D0
        Y(2) = 0.D0
        Y(3) = 0.D0
        GOTO 360

220     CONTINUE
!     PROBLEM D3
        N = 4
        W(1) = 0.100D+01
        W(2) = 0.100D+01
        W(3) = 0.360D+00
        W(4) = 0.485D+00
        XEND = 20.D0
        HBEGIN = 2.5D-5
        HMAX = 20.D0
        DO I = 1, 2
          Y(I) = 1.D0
          Y(I+2) = 0.D0
        END DO
        GOTO 360

230     CONTINUE
!     PROBLEM D4
        N = 3
        W(1) = 0.100D+01
        W(2) = 0.142D+01
        W(3) = 0.371D-05
        XEND = 50.D0
        HBEGIN = 2.9D-4
        HMAX = 50.D0
        Y(1) = 1.D0
        Y(2) = 1.D0
        Y(3) = 0.D0
        GOTO 360

240     CONTINUE
!     PROBLEM D5
        N = 2
        W(1) = 0.992D+00
        W(2) = 0.984D+00
        XEND = 1.D2
        HBEGIN = 1.D-4
        HMAX = 1.D2
        Y(1) = 0.D0
        Y(2) = 0.D0
        GOTO 360

250     CONTINUE
!     PROBLEM D6
        N = 3
        W(1) = 0.100D+01
        W(2) = 0.148D+00
        W(3) = 0.577D-07
        XEND = 1.D0
        HBEGIN = 3.3D-8
        HMAX = 1.D0
        Y(1) = 1.D0
        Y(2) = 0.D0
        Y(3) = 0.D0
        GOTO 360

!     PROBLEM CLASS E - NON-LINEAR WITH NON-REAL EIGENVALUES

260     CONTINUE
!     PROBLEM E1
        N = 4
        W(1) = 0.100D-07
        W(2) = 0.223D-06
        W(3) = 0.132D-04
        W(4) = 0.171D-02
        XEND = 1.D0
        HBEGIN = 6.8D-3
        HMAX = 1.D0
        DO I = 1, N
          Y(I) = 0.D0
        END DO
        GOTO 360

270     CONTINUE
!     PROBLEM E2
        N = 2
        W(1) = 0.202D+01
        W(2) = 0.764D+01
        XEND = 1.D1
        HBEGIN = 1.D-3
        HMAX = 1.D1
        Y(1) = 2.D0
        Y(2) = 0.D0
        GOTO 360

280     CONTINUE
!     PROBLEM E3
        N = 3
        W(1) = 0.163D+01
        W(2) = 0.160D+01
        W(3) = 0.263D+02
        XEND = 5.D2
        HBEGIN = .2D-1
        HMAX = 5.D2
        Y(1) = 1.D0
        Y(2) = 1.D0
        Y(3) = 0.D0
        GOTO 360

290     CONTINUE
!     PROBLEM E4
        N = 4
        W(1) = 0.288D+02
        W(2) = 0.295D+02
        W(3) = 0.155D+02
        W(4) = 0.163D+02
        XEND = 1.D3
        HBEGIN = 1.D-3
        HMAX = 1.D3
        Y(1) = 0.D0
        Y(2) = -2.D0
        Y(3) = -1.D0
        Y(4) = -1.D0
        GOTO 360

300     CONTINUE
!     PROBLEM E5
        N = 4
        W(1) = 0.176D-02
        W(2) = 0.146D-09
        W(3) = 0.827D-11
        W(4) = 0.138D-09
        XEND = 1.D3
        HBEGIN = 5.D-5
        HMAX = 1.D3
        Y(1) = 1.76D-3
        DO I = 2, N
          Y(I) = 0.D0
        END DO
        GOTO 360

!     PROBLEM CLASS F - CHEMICAL KINETICS EQUATIONS

310     CONTINUE
!     PROBLEM F1
        N = 4
        W(1) = 0.121D+04
        W(2) = 0.835D-01
        W(3) = 0.121D+04
        W(4) = 0.100D+00
        HMAX = 1.D3
        HBEGIN = 1.D-4
        XEND = 1.D3
        Y(1) = 761.D0
        Y(2) = 0.D0
        Y(3) = 600.D0
        Y(4) = .1D0
        GOTO 360

320     CONTINUE
!     PROBLEM F2
        N = 2
        W(1) = 0.100D+01
        W(2) = 0.253D-02
        HMAX = 240.D0
        HBEGIN = 1.D-2
        XEND = 240.D0
        Y(1) = 1.0D0
        Y(2) = 0.D0
        GOTO 360

330     CONTINUE
!     PROBLEM F3
        N = 5
        W(1) = 0.400D-05
        W(2) = 0.100D-05
        W(3) = 0.374D-08
        W(4) = 0.765D-06
        W(5) = 0.324D-05
        HBEGIN = 1.D-6
        HBEGIN = 1.D-10
        HMAX = 100.D0
        XEND = 100.D0
        Y(1) = 4.D-6
        Y(2) = 1.D-6
        Y(3) = 0.0D0
        Y(4) = 0.0D0
        Y(5) = 0.0D0
        GOTO 360

340     CONTINUE
!     PROBLEM F4
        N = 3
        W(1) = 0.118D+06
        W(2) = 0.177D+04
        W(3) = 0.313D+05
        HBEGIN = 1.D-3
!       HMAX = 50.D0
        HMAX = 1.D0
        XEND = 300.D0
        Y(1) = 4.D0
        Y(2) = 1.1D0
        Y(3) = 4.D0
        GOTO 360

350     CONTINUE
!     PROBLEM F5
        N = 4
        W(1) = 0.336D-06
        W(2) = 0.826D-02
        W(3) = 0.619D-02
        W(4) = 0.955D-05
        HBEGIN = 1.D-7
        HMAX = 100.D0
        XEND = 100.D0
        Y(1) = 3.365D-7
        Y(2) = 8.261D-3
        Y(3) = 1.642D-3
        Y(4) = 9.380D-6
360     CONTINUE
        IF (IWT<0) GOTO 370
        DO I = 1, N
          Y(I) = Y(I)/W(I)
        END DO
370     CONTINUE
        RETURN

90000   FORMAT (' AN INVALID INTERNAL PROBLEM ID OF ',I4, &
          ' WAS FOUND BY THE IVALU ROUTINE',' RUN TERMINATED. CHECK THE DATA.' &
          )
      END SUBROUTINE IVALU

      SUBROUTINE EVALU(Y)

!     ROUTINE TO PROVIDE THE 'TRUE' SOLUTION OF THE DIFFERENTIAL
!     EQUATION EVALUATED AT THE ENDPOINT OF THE INTEGRATION.

!     1986 REVISION:  SOME VERY SMALL CONSTANTS HAVE BEEN RECAST IN THE
!     (NOT SO SMALL CONST)/(1.E38) TO AVOID COMPILE-TIME UNDERFLOW ERROR
!     IT IS ASSUMED 1E+38 WON'T OVERFLOW.
!     PARAMETER  (OUTPUT)
!        Y      - THE TRUE SOLUTION VECTOR EVALUATED AT THE ENDPOINT

!     PARAMETERS (INPUT)
!        N      - DIMENSION OF THE PROBLEM
!        W      - VECTOR OF WEIGHTS USED TO SCALE THE PROBLEM
!                 IF THIS OPTION IS SELECTED
!        IWT    - FLAG USED TO SIGNAL WHEN THE SCALED PROBLEM IS
!                 BEING SOLVED
!        ID     - FLAG USED TO INDICATE WHICH EQUATION IS BEING
!                 SOLVED

        IMPLICIT NONE

!     .. Parameters ..
        DOUBLE PRECISION TENE38
        PARAMETER (TENE38=1.D38)
!     .. Array Arguments ..
        DOUBLE PRECISION Y(20)
!     .. Local Scalars ..
        INTEGER I

!     .. Executable Statements ..

!     GOTO (20, 40, 60, 80, 620, 620, 620, 620, 620, 620, 100, 120, 140,&
!     160, 180, 620, 620, 620, 620, 620, 200, 220, 240, 260, 280, 620,  &
!     620, 620, 620, 620, 300, 320, 340, 360, 380, 400, 620, 620, 620,  &
!     620, 420, 440, 460, 480, 500, 620, 620, 620, 620, 620, 520, 540,  &
!     560, 580, 600, 620, 620, 620, 620, 620) ID
        GOTO (10,20,30,40,310,310,310,310,310,310,50,60,70,80,90,310,310,310, &
          310,310,100,110,120,130,140,310,310,310,310,310,150,160,170,180,190, &
          200,310,310,310,310,210,220,230,240,250,310,310,310,310,310,260,270, &
          280,290,300,310,310,310,310,310) ID
        GOTO 310

!     PROBLEM CLASS A

!     PROBLEM A1
10      Y(1) = 4.539992969929191D-05
        Y(2) = 2.061153036149920D-09
        Y(3) = 2.823006338263857D-18/TENE38
        Y(4) = 5.235792540515384D-14/TENE38
        GOTO 310

!     PROBLEM A2
20      Y(1) = 9.999912552999704D-02
        Y(2) = 1.999982511586291D-01
        Y(3) = 2.999975543202422D-01
        Y(4) = 3.999971057541257D-01
        Y(5) = 4.999969509963023D-01
        Y(6) = 5.999971057569546D-01
        Y(7) = 6.999975543256127D-01
        Y(8) = 7.999982511659962D-01
        Y(9) = 8.999991255386128D-01
        GOTO 310

!     PROBLEM A3
30      Y(1) = -1.353352661867235D-03
        Y(2) = 1.368526917891521D-02
        Y(3) = 1.503725348455117D+00
        Y(4) = 1.353352832366099D-01
        GOTO 310

!     PROBLEM A4
40      Y(1) = 3.678794411714325D-01
        Y(2) = 1.265870722340194D-14
        Y(3) = 1.911533219339204D-04/TENE38
        Y(4) = 2.277441666729596D-17/TENE38
        Y(5) = 0.0D0
        Y(6) = 0.0D0
        Y(7) = 0.0D0
        Y(8) = 0.0D0
        Y(9) = 0.0D0
        Y(10) = 0.0D0
        GOTO 310

!     PROBLEM CLASS B

!     PROBLEM B1
50      Y(1) = 1.004166730990124D-09
        Y(2) = 1.800023280346500D-08
        Y(3) = 0.0D0
        Y(4) = -6.042962877027475D-03/TENE38/TENE38
        GOTO 310

!     PROBLEM B2
60      Y(1) = 6.181330838820067D-31
        Y(2) = 8.963657877626303D-31
        Y(3) = 2.738406773453261D-27
        Y(4) = 2.061153063164016D-09
        Y(5) = 4.539992973654118D-05
        Y(6) = 1.353352832365270D-01
        GOTO 310

!     PROBLEM B3
70      Y(1) = -1.076790816984970D-28
        Y(2) = 5.455007683862160D-28
        Y(3) = 2.738539964946867D-27
        Y(4) = 2.061153071123456D-09
        Y(5) = 4.539992974611305D-05
        Y(6) = 1.353352832365675D-01
        GOTO 310

!     PROBLEM B4
80      Y(1) = 1.331242472678293D-22
        Y(2) = -2.325916064237926D-22
        Y(3) = 1.517853928534857D-35
        Y(4) = 2.061152428936651D-09
        Y(5) = 4.539992963392291D-05
        Y(6) = 1.353352832363442D-01
        GOTO 310

!     PROBLEM B5
90      Y(1) = -3.100634584292190D-14
        Y(2) = 3.862788998076547D-14
        Y(3) = 1.804851385304217D-35
        Y(4) = 2.061153622425655D-09
        Y(5) = 4.539992976246673D-05
        Y(6) = 1.353352832366126D-01
        GOTO 310

!     PROBLEM CLASS C

!     PROBLEM C1
100     Y(1) = 4.003223925456179D-04
        Y(2) = 4.001600000000000D-04
        Y(3) = 4.000000000000000D-04
        Y(4) = 2.000000000000000D-02
        GOTO 310

!     PROBLEM C2
110     Y(1) = 1.999999997938994D+00
        Y(2) = 3.999999990839974D-02
        Y(3) = 4.001599991537078D-02
        Y(4) = 4.003201271914461D-02
        GOTO 310

!     PROBLEM C3
120     Y(1) = 1.999999997939167D+00
        Y(2) = 3.999999990840744D-01
        Y(3) = 4.159999990793773D-01
        Y(4) = 4.333055990159567D-01
        GOTO 310

!     PROBLEM C4
130     Y(1) = 1.999999997938846D+00
        Y(2) = 3.999999990839318D+00
        Y(3) = 1.999999991637941D+01
        Y(4) = 4.199999965390368D+02
        GOTO 310

!     PROBLEM C5
140     Y(1) = 1.999999997938846D+00
        Y(2) = 7.999999981678634D+00
        Y(3) = 1.359999993817714D+02
        Y(4) = 3.712799965967762D+04
        GOTO 310

!     PROBLEM CLASS D

!     PROBLEM D1
150     Y(1) = 2.224222010616901D+01
        Y(2) = 2.711071334484136D+01
        Y(3) = 3.999999999999999D+02
        GOTO 310

!     PROBLEM D2
160     Y(1) = 7.158270687193941D-01
        Y(2) = 9.185534764557338D-02
        Y(3) = 2.841637457458413D+01
        GOTO 310

!     PROBLEM D3
170     Y(1) = 6.397604446889910D-01
        Y(2) = 5.630850708287990D-03
        Y(3) = 3.602395553110090D-01
        Y(4) = 3.170647969903515D-01
        GOTO 310

!     PROBLEM D4
180     Y(1) = 5.976546980673215D-01
        Y(2) = 1.402343408546138D+00
        Y(3) = -1.893386540441913D-06
        GOTO 310

!     PROBLEM D5
190     Y(1) = -9.916420698713913D-01
        Y(2) = 9.833363588544478D-01
        GOTO 310

!     PROBLEM D6
200     Y(1) = 8.523995440749948D-01
        Y(2) = 1.476003981941319D-01
        Y(3) = 5.773087333950041D-08
        GOTO 310

!     PROBLEM CLASS E

!     PROBLEM E1
210     Y(1) = 1.000000000000012D-08
        Y(2) = -1.625323873316817D-19
        Y(3) = 2.025953375595861D-17
        Y(4) = -1.853149807630002D-15
        GOTO 310

!     PROBLEM E2
220     Y(1) = -1.158701266031984D+00
        Y(2) = 4.304698089780476D-01
        GOTO 310

!     PROBLEM E3
230     Y(1) = 4.253052197643089D-03
        Y(2) = 5.317019548450387D-03
        Y(3) = 2.627647748753926D+01
        GOTO 310

!     PROBLEM E4
240     Y(1) = 1.999999977523654D+01
        Y(2) = -2.000000022476345D+01
        Y(3) = -2.247634567084293D-07
        Y(4) = 2.247634567084293D-07
        GOTO 310

!     PROBLEM E5
250     Y(1) = 1.618076919919600D-03
        Y(2) = 1.382236955418478D-10
        Y(3) = 8.251573436034144D-12
        Y(4) = 1.299721221058136D-10
        GOTO 310

!     PROBLEM CLASS F

!     PROBLEM F1
260     Y(1) = 1.211129474696585D+03
        Y(2) = 1.271123619113051D-05
        Y(3) = 1.208637804660361D+03
        Y(4) = 3.241981171933418D-04
        GOTO 310

!     PROBLEM F2
270     Y(1) = 3.912699122292088D-01
        Y(2) = 1.329964166084866D-03
        GOTO 310

!     PROBLEM F3
280     Y(1) = 3.235910070806680D-13
        Y(2) = 2.360679774997897D-07
        Y(3) = 7.639319089351045D-14
        Y(4) = 7.639319461070194D-07
        Y(5) = 3.236067653908783D-06
        GOTO 310

!     PROBLEM F4
290     Y(1) = 4.418303324022590D+00
        Y(2) = 1.290244712916425D+00
        Y(3) = 3.019282584050490D+00
        GOTO 310

!     PROBLEM F5
300     Y(1) = 1.713564284690712D-07
        Y(2) = 3.713563071160676D-03
        Y(3) = 6.189271785267793D-03
        Y(4) = 9.545143571530929D-06
310     CONTINUE
        IF (IWT<0) GOTO 320
        DO I = 1, N
          Y(I) = Y(I)/W(I)
        END DO
320     CONTINUE
        RETURN
      END SUBROUTINE EVALU

      SUBROUTINE FCN(X,Y,YP)

!     ROUTINE TO EVALUATE THE DERIVATIVE F(X,Y) CORRESPONDING TO
!     THE DIFFERENTIAL EQUATION:
!                    DY/DX = F(X,Y) .
!     THE ROUTINE STORES THE VECTOR OF DERIVATIVES IN YP(*). THE
!     PARTICULAR EQUATION BEING INTEGRATED IS INDICATED BY THE
!     VALUE OF THE FLAG ID WHICH IS PASSED THROUGH COMMON. THE
!     DIFFERENTIAL EQUATION IS SCALED BY THE WEIGHT VECTOR W(*)
!     IF THIS OPTION HAS BEEN SELECTED (IF SO IT IS SIGNALLED
!     BY THE FLAG IWT).

        IMPLICIT NONE

!     .. Scalar Arguments ..
        DOUBLE PRECISION X
!     .. Array Arguments ..
        DOUBLE PRECISION Y(20), YP(20)
!     .. Local Scalars ..
        DOUBLE PRECISION F, Q, S, SUM, T, TEMP, XTEMP
        INTEGER I
!     .. Local Arrays ..
        DOUBLE PRECISION BPARM(4), CPARM(4), VECT1(4), VECT2(4), YTEMP(20)
!     .. Data statements ..
        DATA BPARM/3.D0, 8.D0, 25.D0, 1.D2/
        DATA CPARM/1.D-1, 1.D0, 1.D1, 2.D1/

!     .. Executable Statements ..

        NFCN = NFCN + 1
        IF (IWT<0) GOTO 10
        DO I = 1, N
          YTEMP(I) = Y(I)
          Y(I) = Y(I)*W(I)
        END DO
10      CONTINUE
        GOTO (20,30,40,50,260,260,260,260,260,260,60,70,70,70,70,260,260,260, &
          260,260,80,90,90,90,90,260,260,260,260,260,100,110,120,130,140,150, &
          260,260,260,260,160,170,180,190,200,260,260,260,260,260,210,220,230, &
          240,250) ID
        GOTO 260

!     PROBLEM CLASS A - LINEAR WITH REAL EIGENVALUES

!     PROBLEM A1
20      YP(1) = -.5D0*Y(1)
        YP(2) = -1.D0*Y(2)
        YP(3) = -1.D2*Y(3)
        YP(4) = -9.D1*Y(4)
        GOTO 260

!     PROBLEM A2
30      YP(1) = -1.8D3*Y(1) + 9.D2*Y(2)
        DO I = 2, 8
          YP(I) = Y(I-1) - 2.D0*Y(I) + Y(I+1)
        END DO
        YP(9) = 1.D3*Y(8) - 2.D3*Y(9) + 1.D3
        GOTO 260

!     PROBLEM A3
40      YP(1) = -1.D4*Y(1) + 1.D2*Y(2) - 1.D1*Y(3) + 1.D0*Y(4)
        YP(2) = -1.D3*Y(2) + 1.D1*Y(3) - 1.D1*Y(4)
        YP(3) = -1.D0*Y(3) + 1.D1*Y(4)
        YP(4) = -1.D-1*Y(4)
        GOTO 260

!     PROBLEM A4
50      DO I = 1, 10
          YP(I) = -REAL(I)**5*Y(I)
        END DO
        GOTO 260

!     PROBLEM CLASS B - LINEAR WITH NON-REAL EIGENVALUES

!     PROBLEM B1
60      YP(1) = -Y(1) + Y(2)
        YP(2) = -1.D2*Y(1) - Y(2)
        YP(3) = -1.D2*Y(3) + Y(4)
        YP(4) = -1.D4*Y(3) - 1.D2*Y(4)
        GOTO 260

!     PROBLEMS B2, B3, B4, B5
70      YP(1) = -1.D1*Y(1) + BPARM(IID-1)*Y(2)
        YP(2) = -BPARM(IID-1)*Y(1) - 1.D1*Y(2)
        YP(3) = -4.D0*Y(3)
        YP(4) = -1.D0*Y(4)
        YP(5) = -.5D0*Y(5)
        YP(6) = -.1D0*Y(6)
        GOTO 260

!     PROBLEM CLASS C - NON-LINEAR COUPLING FROM
!                       STEADY STATE TO TRANSIENT

!     PROBLEM C1
80      YP(1) = -Y(1) + (Y(2)*Y(2)+Y(3)*Y(3)+Y(4)*Y(4))
        YP(2) = -1.D1*Y(2) + 1.D1*(Y(3)*Y(3)+Y(4)*Y(4))
        YP(3) = -4.D1*Y(3) + 4.D1*Y(4)*Y(4)
        YP(4) = -1.D2*Y(4) + 2.D0
        GOTO 260

!     PROBLEMS C2, C3, C4, C5
90      YP(1) = -Y(1) + 2.D0
        YP(2) = -1.D1*Y(2) + CPARM(IID-1)*Y(1)*Y(1)
        YP(3) = -4.D1*Y(3) + (Y(1)*Y(1)+Y(2)*Y(2))*CPARM(IID-1)*4.D0
        YP(4) = (Y(1)*Y(1)+Y(2)*Y(2)+Y(3)*Y(3))*CPARM(IID-1)*1.D1 - 1.D2*Y(4)
        GOTO 260

!     PROBLEM CLASS D - NON-LINEAR WITH REAL EIGENVALUES

!     PROBLEM D1
100     YP(1) = .2D0*Y(2) - .2D0*Y(1)
        YP(2) = 1.D1*Y(1) - (6.D1-.125D0*Y(3))*Y(2) + .125D0*Y(3)
        YP(3) = 1.D0
        GOTO 260

!     PROBLEM D2
110     YP(1) = -.04D0*Y(1) + .01D0*Y(2)*Y(3)
        YP(2) = 4.D2*Y(1) - 1.D2*Y(2)*Y(3) - 3.D3*Y(2)**2
        YP(3) = 3.D1*Y(2)**2
        GOTO 260

!     PROBLEM D3
120     YP(1) = Y(3) - 1.D2*Y(1)*Y(2)
        YP(3) = -YP(1)
        YP(4) = -Y(4) + 1.D4*Y(2)**2
        YP(2) = YP(1) - YP(4) + Y(4) - 1.D4*Y(2)**2
        GOTO 260

!     PROBLEM D4
130     YP(1) = -.013D0*Y(1) - 1.D3*Y(1)*Y(3)
        YP(2) = -2.5D3*Y(2)*Y(3)
        YP(3) = YP(1) + YP(2)
        GOTO 260

!     PROBLEM D5
140     XTEMP = .01D0 + Y(1) + Y(2)
        YP(1) = .01D0 - XTEMP*(1.D0+(Y(1)+1.D3)*(Y(1)+1.D0))
        YP(2) = .01D0 - XTEMP*(1.D0+Y(2)**2)
        GOTO 260

!     PROBLEM D6
150     YP(1) = -Y(1) + 1.D8*Y(3)*(1.D0-Y(1))
        YP(2) = -1.D1*Y(2) + 3.D7*Y(3)*(1.D0-Y(2))
        YP(3) = -YP(1) - YP(2)
        GOTO 260

!     PROBLEM CLASS E - NON-LINEAR WITH NON-REAL EIGENVALUES

!     PROBLEM E1
160     YP(1) = Y(2)
        YP(2) = Y(3)
        YP(3) = Y(4)
        YP(4) = (Y(1)**2-SIN(Y(1))-1.D8)*Y(1) + (Y(2)*Y(3)/(Y(1)**2+1.D0)-4.D6 &
          )*Y(2) + (1.D0-6.D4)*Y(3) + (1.D1*EXP(-Y(4)**2)-4.D2)*Y(4) + 1.D0
        GOTO 260

!     PROBLEM E2
170     YP(1) = Y(2)
        YP(2) = 5.D0*Y(2) - 5.D0*Y(1)*Y(1)*Y(2) - Y(1)
        GOTO 260

!     PROBLEM E3
180     YP(1) = -55.D0*Y(1) - Y(3)*Y(1) + 65.D0*Y(2)
        YP(2) = .785D-1*Y(1) - .785D-1*Y(2)
        YP(3) = .1D0*Y(1)
        GOTO 260

!     PROBLEM E4
190     SUM = Y(1) + Y(2) + Y(3) + Y(4)
        DO I = 1, 4
          VECT2(I) = -Y(I) + .5D0*SUM
        END DO
        VECT1(1) = .5D0*(VECT2(1)**2-VECT2(2)**2)
        VECT1(2) = VECT2(1)*VECT2(2)
        VECT1(3) = VECT2(3)**2
        VECT1(4) = VECT2(4)**2
        TEMP = -1.D1*VECT2(1) - 1.D1*VECT2(2)
        VECT2(2) = 1.D1*VECT2(1) - 1.D1*VECT2(2)
        VECT2(1) = TEMP
        VECT2(3) = 1.D3*VECT2(3)
        VECT2(4) = 1.D-2*VECT2(4)
        SUM = 0.D0
        DO I = 1, 4
          SUM = SUM + VECT1(I) - VECT2(I)
        END DO
        DO I = 1, 4
          YP(I) = VECT2(I) - VECT1(I) + .5D0*SUM
        END DO
        GOTO 260

!     PROBLEM E5
200     XTEMP = -7.89D-10*Y(1)
        YP(1) = XTEMP - 1.1D7*Y(1)*Y(3)
        YP(2) = -XTEMP - 1.13D9*Y(2)*Y(3)
        YP(4) = 1.1D7*Y(1)*Y(3) - 1.13D3*Y(4)
        YP(3) = YP(2) - YP(4)
        GOTO 260

!     PROBLEM CLASS F - CHEMICAL KINETICS EQUATIONS

!     PROBLEM F1
210     TEMP = 6.D-3*EXP(20.7D0-1.5D4/Y(1))
        YP(1) = 1.3D0*(Y(3)-Y(1)) + 1.04D4*TEMP*Y(2)
        YP(2) = 1.88D3*(Y(4)-Y(2)*(1.D0+TEMP))
        YP(3) = 1752.D0 - 269.D0*Y(3) + 267.D0*Y(1)
        YP(4) = .1D0 + 320.D0*Y(2) - 321.D0*Y(4)
        GOTO 260

!     PROBLEM F2
220     YP(1) = -Y(1) - Y(1)*Y(2) + 294.D0*Y(2)
        YP(2) = Y(1)*(1.D0-Y(2))/98.D0 - 3.D0*Y(2)
        GOTO 260

!     PROBLEM F3
230     YP(1) = -1.0D7*Y(2)*Y(1) + 1.D1*Y(3)
        YP(2) = -1.0D7*Y(2)*Y(1) - 1.D7*Y(2)*Y(5) + 1.D1*Y(3) + 1.D1*Y(4)
        YP(3) = 1.0D7*Y(2)*Y(1) - 1.001D4*Y(3) + 1.D-3*Y(4)
        YP(4) = 1.D4*Y(3) - 1.0001D1*Y(4) + 1.D7*Y(2)*Y(5)
        YP(5) = 1.D1*Y(4) - 1.D7*Y(2)*Y(5)
        GOTO 260

!     PROBLEM F4
240     S = 77.27D0
        T = 0.161D0
        Q = 8.375D-6
        F = 1.D0
        YP(1) = S*(Y(2)-Y(1)*Y(2)+Y(1)-Q*Y(1)*Y(1))
        YP(2) = (-Y(2)-Y(1)*Y(2)+F*Y(3))/S
        YP(3) = T*(Y(1)-Y(3))
        GOTO 260

!     PROBLEM F5
250     YP(1) = -3.D11*Y(1)*Y(2) + 1.2D8*Y(4) - 9.D11*Y(1)*Y(3)
        YP(2) = -3.D11*Y(1)*Y(2) + 2.D7*Y(4)
        YP(3) = -9.D11*Y(1)*Y(3) + 1.D8*Y(4)
        YP(4) = 3.D11*Y(1)*Y(2) - 1.2D8*Y(4) + 9.D11*Y(1)*Y(3)
260     CONTINUE
        IF (IWT<0) GOTO 270
        DO I = 1, N
          YP(I) = YP(I)/W(I)
          Y(I) = YTEMP(I)
        END DO
270     CONTINUE
        RETURN
      END SUBROUTINE FCN

      SUBROUTINE PDERV(X,Y)

!     ROUTINE TO EVALUATE THE JACOBIAN MATRIX OF PARTIAL DERIVATIVES
!     CORRESPONDING TO THE DIFFERENTIAL EQUATION:
!                   DY/DX = F(X,Y).
!     THE N**2 ELEMENTS OF THE ARRAY DY(*) ARE ASSIGNED THE VALUE OF
!     THE JACOBIAN MATRIX WITH ELEMENT I+(J-1)*N BEING ASSIGNED THE
!     VALUE OF DF(I)/DY(J). THE PARTICULAR EQUATION BEING INTEGRATED
!     IS INDICATED BY THE VALUE OF THE FLAG ID WHICH IS PASSED THROUGH
!     COMMON. IF A SCALED DIFFERENTIAL EQUATION IS BEING SOLVED (AS
!     SIGNALLED IWT) THE ELEMENTS OF THE JACOBIAN ARE SCALED ACCORDING-
!     LY BY THE WEIGHT VECTOR W(*).

        IMPLICIT NONE

!     .. Scalar Arguments ..
        DOUBLE PRECISION X
!     .. Array Arguments ..
        DOUBLE PRECISION Y(20)
!     .. Local Scalars ..
        DOUBLE PRECISION F, Q, S, SUM, T, TEMP, XTEMP1, XTEMP2, XTEMP3
        INTEGER I, ITMP, J, L
!     .. Local Arrays ..
        DOUBLE PRECISION BPARM(4), CPARM(4), VECT2(4), YTEMP(20)
!     .. Data statements ..
        DATA BPARM/3.D0, 8.D0, 25.D0, 1.D2/
        DATA CPARM/1.D-1, 1.D0, 1.D1, 2.D1/

!     .. Executable Statements ..

        NJAC = NJAC + 1
        IF (IWT<0) GOTO 10
        DO I = 1, N
          YTEMP(I) = Y(I)
          Y(I) = Y(I)*W(I)
        END DO
10      CONTINUE
        GOTO (20,30,40,50,260,260,260,260,260,260,60,70,70,70,70,260,260,260, &
          260,260,80,90,90,90,90,260,260,260,260,260,100,110,120,130,140,150, &
          260,260,260,260,160,170,180,190,200,260,260,260,260,260,210,220,230, &
          240,250) ID
        GOTO 260

!     PROBLEM CLASS A - LINEAR WITH REAL EIGENVALUES

!     PROBLEM A1
20      DO I = 1, 16
          DY(I) = 0.D0
        END DO
        DY(1) = -.5D0
        DY(6) = -1.D0
        DY(11) = -1.D2
        DY(16) = -9.D1
        GOTO 260

!     PROBLEM A2
30      DO I = 1, 81
          DY(I) = 0.D0
        END DO
        DO I = 2, 62, 10
          DY(I) = 1.D0
          DY(I+9) = -2.D0
          DY(I+18) = 1.D0
        END DO
        DY(1) = -1.8D3
        DY(10) = 9.D2
        DY(72) = 1.D3
        DY(81) = -2.D3
        GOTO 260

!     PROBLEM A3
40      DO I = 1, 16
          DY(I) = 0.D0
        END DO
        DY(1) = -1.D4
        DY(5) = 1.D2
        DY(6) = -1.D3
        DY(9) = -1.D1
        DY(10) = 1.D1
        DY(11) = -1.D0
        DY(13) = 1.D0
        DY(14) = -1.D1
        DY(15) = 1.D1
        DY(16) = -1.D-1
        GOTO 260

!     PROBLEM A4
50      DO I = 1, 100
          DY(I) = 0.D0
        END DO
        DO I = 1, 10
          DY((I-1)*10+I) = -REAL(I)**5
        END DO
        GOTO 260

!     PROBLEM CLASS B - LINEAR WITH NON-REAL EIGENVALUES

!     PROBLEM B1
60      DO I = 1, 16
          DY(I) = 0.D0
        END DO
        DY(1) = -1.D0
        DY(2) = -1.D2
        DY(5) = 1.D0
        DY(6) = -1.D0
        DY(11) = -1.D2
        DY(12) = -1.D4
        DY(15) = 1.D0
        DY(16) = -1.D2
        GOTO 260

!     PROBLEMS B2, B3, B4, B5
70      DO I = 1, 36
          DY(I) = 0.D0
        END DO
        DY(1) = -1.D1
        DY(2) = -BPARM(IID-1)
        DY(7) = BPARM(IID-1)
        DY(8) = -1.D1
        DY(15) = -4.D0
        DY(22) = -1.D0
        DY(29) = -.5D0
        DY(36) = -.1D0
        GOTO 260

!     PROBLEM CLASS C - NON-LINEAR COUPLING FROM
!                       STEADY STATE TO TRANSIENT

!     PROBLEM C1
80      DO I = 1, 16
          DY(I) = 0.D0
        END DO
        DY(1) = -1.D0
        DY(5) = 2.D0*Y(2)
        DY(6) = -1.D1
        DY(9) = 2.D0*Y(3)
        DY(10) = 2.D1*Y(3)
        DY(11) = -4.D1
        DY(13) = 2.D0*Y(4)
        DY(14) = 2.D1*Y(4)
        DY(15) = 8.D1*Y(4)
        DY(16) = -1.D2
        GOTO 260

!     PROBLEMS C2, C3, C4, C5
90      DO I = 1, 16
          DY(I) = 0.D0
        END DO
        DY(1) = -1.D0
        DY(2) = 2.D0*Y(1)*CPARM(IID-1)
        DY(3) = 8.D0*Y(1)*CPARM(IID-1)
        DY(4) = 2.D1*Y(1)*CPARM(IID-1)
        DY(6) = -1.D1
        DY(7) = 8.D0*Y(2)*CPARM(IID-1)
        DY(8) = 2.D1*Y(2)*CPARM(IID-1)
        DY(11) = -4.D1
        DY(12) = 2.D1*Y(3)*CPARM(IID-1)
        DY(16) = -1.D2
        GOTO 260

!     PROBLEM CLASS D - NON-LINEAR WITH REAL EIGENVALUES

!     PROBLEM D1
100     DY(1) = -.2D0
        DY(2) = 1.D1
        DY(3) = 0.D0
        DY(4) = .2D0
        DY(5) = -6.D1 + .125D0*Y(3)
        DY(6) = 0.D0
        DY(7) = 0.D0
        DY(8) = .125D0*Y(2) + .125D0
        DY(9) = 0.D0
        GOTO 260

!     PROBLEM D2
110     DY(1) = -4.D-2
        DY(2) = 4.D2
        DY(3) = 0.D0
        DY(4) = 1.D-2*Y(3)
        DY(5) = -1.D2*Y(3) - 6.D3*Y(2)
        DY(6) = 6.D1*Y(2)
        DY(7) = .1D-1*Y(2)
        DY(8) = -1.D2*Y(2)
        DY(9) = 0.D0
        GOTO 260

!     PROBLEM D3
120     DY(1) = -1.D2*Y(2)
        DY(2) = DY(1)
        DY(3) = -DY(1)
        DY(4) = 0.D0
        DY(5) = -1.D2*Y(1)
        DY(7) = -DY(5)
        DY(8) = 2.D4*Y(2)
        DY(6) = DY(5) - DY(8)
        DY(6) = DY(6) - 2.D4*Y(2)
        DY(9) = 1.D0
        DY(10) = 1.D0
        DY(11) = -1.D0
        DY(12) = 0.D0
        DY(13) = 0.D0
        DY(14) = 2.D0
        DY(15) = 0.D0
        DY(16) = -1.D0
        GOTO 260

!     PROBLEM D4
130     DY(1) = -.013D0 - 1.D3*Y(3)
        DY(2) = 0.D0
        DY(4) = 0.D0
        DY(5) = -2.5D3*Y(3)
        DY(7) = -1.D3*Y(1)
        DY(8) = -2.5D3*Y(2)
        DO I = 3, 9, 3
          DY(I) = DY(I-1) + DY(I-2)
        END DO
        GOTO 260

!     PROBLEM D5
140     XTEMP1 = Y(1) + 1.D3
        XTEMP2 = Y(1) + 1.D0
        XTEMP3 = .01D0 + Y(1) + Y(2)
        DY(2) = -(1.D0+Y(2)**2)
        DY(3) = -(1.D0+XTEMP1*XTEMP2)
        DY(1) = -(-DY(3)+XTEMP3*(XTEMP1+XTEMP2))
        DY(4) = -(2.D0*XTEMP3*Y(2)-DY(2))
        GOTO 260

!     PROBLEM D6
150     DY(1) = -1.D0 - 1.D8*Y(3)
        DY(2) = 0.D0
        DY(4) = 0.D0
        DY(5) = -1.D1 - 3.D7*Y(3)
        DY(7) = 1.D8*(1.D0-Y(1))
        DY(8) = 3.D7*(1.D0-Y(2))
        DO I = 3, 9, 3
          DY(I) = -DY(I-2) - DY(I-1)
        END DO
        GOTO 260

!     PROBLEM CLASS E - NON-LINEAR WITH NON-REAL EIGENVALUES

!     PROBLEM E1
160     DO I = 1, 16
          DY(I) = 0.D0
        END DO
        DY(5) = 1.D0
        DY(10) = 1.D0
        DY(15) = 1.D0
        XTEMP1 = Y(1)
        XTEMP2 = Y(2)/(XTEMP1**2+1.D0)**2
        DY(4) = 3.D0*XTEMP1**2 - XTEMP1*COS(XTEMP1) - SIN(XTEMP1) - 1.D8 - &
          2.D0*XTEMP1*Y(2)*Y(3)*XTEMP2
        DY(8) = 2.D0*Y(3)*Y(2)/(1.D0+Y(1)**2) - 4.D6
        DY(12) = Y(2)*Y(2)/(1.D0+Y(1)**2) + 1.D0 - 6.D4
        DY(16) = 1.D1*EXP(-Y(4)**2)*(1.D0-2.D0*Y(4)**2) - 4.D2
        GOTO 260

!     PROBLEM E2
170     DY(1) = 0.D0
        DY(2) = -1.D1*Y(1)*Y(2) - 1.D0
        DY(3) = 1.D0
        DY(4) = 5.D0 - 5.D0*Y(1)*Y(1)
        GOTO 260

!     PROBLEM E3
180     DY(1) = -55.D0 - Y(3)
        DY(2) = .785D-1
        DY(3) = 0.1D0
        DY(4) = 65.D0
        DY(5) = -.785D-1
        DY(6) = 0.D0
        DY(7) = -Y(1)
        DY(8) = 0.D0
        DY(9) = 0.D0
        GOTO 260

!     PROBLEM E4
190     SUM = Y(1) + Y(2) + Y(3) + Y(4)
        DO I = 1, 4
          VECT2(I) = -Y(I) + .5D0*SUM
        END DO
        DO I = 1, 16
          DY(I) = 0.D0
        END DO
        DY(1) = VECT2(1) + 1.D1
        DY(2) = VECT2(2) - 1.D1
        DY(5) = -DY(2)
        DY(6) = DY(1)
        DY(11) = 2.D0*VECT2(3) - 1.D3
        DY(16) = 2.D0*VECT2(4) - 1.D-2
        DO I = 1, 4
          SUM = 0.D0
          DO J = 1, 4
            L = I + (J-1)*4
            SUM = SUM + DY(L)
          END DO
          DO J = 1, 4
            L = I + (J-1)*4
            DY(L) = -DY(L) + .5D0*SUM
          END DO
        END DO
        DO J = 1, 4
          SUM = 0.D0
          DO I = 1, 4
            L = I + (J-1)*4
            SUM = SUM + DY(L)
          END DO
          DO I = 1, 4
            L = I + (J-1)*4
            DY(L) = -DY(L) + .5D0*SUM
          END DO
        END DO
        GOTO 260

!     PROBLEM E5
200     DY(1) = -7.89D-10 - 1.1D7*Y(3)
        DY(2) = 7.89D-10
        DY(4) = 1.1D7*Y(3)
        DY(5) = 0.D0
        DY(6) = -1.13D9*Y(3)
        DY(8) = 0.D0
        DY(9) = -1.1D7*Y(1)
        DY(10) = -1.13D9*Y(2)
        DY(12) = -DY(9)
        DY(13) = 0.D0
        DY(14) = 0.D0
        DY(16) = -1.13D3
        DO I = 3, 15, 4
          DY(I) = DY(I-1) - DY(I+1)
        END DO
        GOTO 260

!     PROBLEM CLASS F - CHEMICAL KINETICS EQUATIONS

!     PROBLEM F1
210     TEMP = 90.D0*EXP(20.7D0-1.5D4/Y(1))/Y(1)**2
        DY(1) = -1.3D0 + 1.04D4*TEMP*Y(2)
        DY(2) = -1.88D3*Y(2)*TEMP
        DY(3) = 267.D0
        DY(4) = 0.D0
        TEMP = 6.D-3*EXP(20.7D0-1.5D4/Y(1))
        DY(5) = 1.04D4*TEMP
        DY(6) = -1.88D3*(1.D0+TEMP)
        DY(7) = 0.D0
        DY(8) = 320.D0
        DY(9) = 1.3D0
        DY(10) = 0.D0
        DY(11) = -269.D0
        DY(12) = 0.0D0
        DY(13) = 0.0D0
        DY(14) = 1.88D3
        DY(15) = 0.0D0
        DY(16) = -321.0D0
        GOTO 260

!     PROBLEM F2
220     DY(1) = -1.D0 - Y(2)
        DY(2) = (1.D0-Y(2))/98.D0
        DY(3) = -Y(1) + 294.D0
        DY(4) = -Y(1)/98.D0 - 3.D0
        GOTO 260

!     PROBLEM F3
230     DY(1) = -1.D7*Y(2)
        DY(2) = -1.D7*Y(2)
        DY(3) = 1.D7*Y(2)
        DY(4) = 0.0D0
        DY(5) = 0.0D0
        DY(6) = -1.D7*Y(1)
        DY(7) = -1.D7*Y(1) - 1.D7*Y(5)
        DY(8) = 1.D7*Y(1)
        DY(9) = 1.D7*Y(5)
        DY(10) = -1.D7*Y(5)
        DY(11) = 1.D1
        DY(12) = 1.D1
        DY(13) = -1.001D4
        DY(14) = 1.D4
        DY(15) = 0.0D0
        DY(16) = 0.0D0
        DY(17) = 1.D1
        DY(18) = 1.D-3
        DY(19) = -1.0001D1
        DY(20) = 1.D1
        DY(21) = 0.0D0
        DY(22) = -1.D7*Y(2)
        DY(23) = 0.0D0
        DY(24) = 1.D7*Y(2)
        DY(25) = -1.0D7*Y(2)
        GOTO 260

!     PROBLEM F4
240     S = 77.27D0
        T = 0.161D0
        Q = 8.375D-6
        F = 1.D0
        DY(1) = S*(-Y(2)+1.D0-2.D0*Q*Y(1))
        DY(2) = -Y(2)/S
        DY(3) = T
        DY(4) = S*(1.D0-Y(1))
        DY(5) = (-1.D0-Y(1))/S
        DY(6) = 0.D0
        DY(7) = 0.D0
        DY(8) = F/S
        DY(9) = -T
        GOTO 260

!     PROBLEM F5
250     DY(1) = -3.D11*Y(2) - 9.D11*Y(3)
        DY(2) = -3.D11*Y(2)
        DY(3) = -9.D11*Y(3)
        DY(4) = 3.D11*Y(2) + 9.D11*Y(3)
        DY(5) = -3.D11*Y(1)
        DY(6) = -3.D11*Y(1)
        DY(7) = 0.0D0
        DY(8) = 3.D11*Y(1)
        DY(9) = -9.D11*Y(1)
        DY(10) = 0.0D0
        DY(11) = -9.D11*Y(1)
        DY(12) = 9.D11*Y(1)
        DY(13) = 1.2D8
        DY(14) = 2.D7
        DY(15) = 1.D8
        DY(16) = -1.2D8
260     CONTINUE
        IF (IWT<0) GOTO 270
        DO I = 1, N
          Y(I) = YTEMP(I)
          DO J = 1, N
            ITMP = I + (J-1)*N
            DY(ITMP) = DY(ITMP)*W(J)/W(I)
          END DO
        END DO
270     CONTINUE
        RETURN
      END SUBROUTINE PDERV

    END MODULE STIFFSET

!******************************************************************

    PROGRAM DEMOSTIFF

      USE STIFFSET
      USE DVODE_F90_M

      IMPLICIT NONE
      INTEGER ITASK, ISTATE, ISTATS, NEQ, I, CLASS, PROBLEM, MYID, ITEST, ISTR, &
        IJAC, JTYPE, CONSTANTJ, ITOL, IADIM, IA, JADIM, JA, ROW, COL, ML, MU,   &
        IJACSP , NUM_S, NUM_J, NUM_D, Total_Solutions, ISCALE, ITIME, COMPS,    &
        METHOD, WHICH_PROBS, WHICH_TOLS, COUNTS
      DOUBLE PRECISION RSTATS, T, TOUT, HBEGIN, HBOUND, TBEGIN, TEND, Y, EPS,   &
        YINIT, YFINAL, RELERR_TOLERANCES, ABSERR_TOLERANCES, AERROR, ERRMAX,    &
        ERSTATS, LBOUNDS, UBOUNDS, ERSTATS2
      REAL DVTIME, DVTIME1, DVTIME2, EXECUTION_TIMES 
      DIMENSION Y(20), RSTATS(22), ISTATS(31), YINIT(20), YFINAL(20), MYID(55)
      DIMENSION RELERR_TOLERANCES(20), ABSERR_TOLERANCES(20), AERROR(20),       &
        IA(21), JA(400), ERSTATS(12,55,3), NUM_S(12,55,3), NUM_J(12,55,3),      &
        NUM_D(12,55,3), EXECUTION_TIMES(12), COMPS(20), LBOUNDS(20),            &
        UBOUNDS(20), ERSTATS2(2,12,55,6), WHICH_PROBS(55), WHICH_TOLS(12),      &
        COUNTS(2,12,55,6,3)
      LOGICAL J_IS_CONSTANT, ANALYTIC_JACOBIAN, ERRORS, SUPPLY_STRUCTURE,       &
        TOLVEC, USE_JACSP, BOUNDS, USEW, USEHBEGIN
      INTRINSIC ABS, MAX, MAXVAL, EPSILON
      TYPE (VODE_OPTS) :: OPTIONS
      DATA MYID/1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 11, 12, 13, 14, 15,  0, 0, 0, 0,  &
        0, 21, 22, 23, 24, 25, 0, 0, 0, 0, 0, 31, 32, 33, 34, 35, 36, 0, 0, 0,  &
        0, 41, 42, 43, 44, 45, 0, 0, 0, 0, 0, 51, 52, 53, 54, 55/

      OPEN (UNIT=6,FILE='stiffoptions.dat')

!     Initialize the cumulative solutions total.
      Total_Solutions = 0

!     Use vector tolerances (same solution in either case)?
      TOLVEC = .FALSE.
      
!     Initialize the flag to indicate whether any errors occurred.
      ERRORS = .FALSE.

!     Execution times per tolerance.
      EXECUTION_TIMES(1:12) = 0.0E0
      ERSTATS2(1:2,1:12,1:55,1:6) = 0.0D0
      WHICH_PROBS(1:55) = 0
      WHICH_TOLS(1:12) = 0
      COUNTS(1:2,1:12,1:55,1:6,1:3) = 0

      WRITE(6,90032)

!     Enter the scale equations loop.
!        1 = the ODEs will be scaled
!        2 = the ODEs will not be scaled
      DO ISCALE = 1 , 2
      USEW = .TRUE.
      IF (ISCALE == 2) USEW = .FALSE.

!     Enter the numerical Jacobian formation loop.
!        1 = use Doug's JACSP for Jacobian
!        2 = use original VODE Jacobian algorithm
      DO IJACSP = 1, 2
      USE_JACSP = .TRUE.
      IF (IJACSP == 2) USE_JACSP = .FALSE. 

!     Initialize the max error array for this value of IJACSP.
      ERSTATS(1:12,1:15,1:3) = 0.0D0

!     Initialize the stats arrays.
      NUM_S(1:12,1:15,1:3) = 0
      NUM_J(1:12,1:15,1:3) = 0
      NUM_D(1:12,1:15,1:3) = 0

!     Enter the sparsity structure loop.
!        1 = numerically determined sparse pointer arrays
!        2 = supplied pointer arrays
!     Make this a one-trip loop (1 to 1 or 2 to 2) for stats to be meaningful:
      DO ISTR = 2, 2
      SUPPLY_STRUCTURE = .FALSE.
      IF (ISTR == 2) SUPPLY_STRUCTURE = .TRUE.

!     Enter the error tolerance loop.
!        Error tolerance = 10^(-ITOL)
      DO ITOL = 5, 12, 1

      WHICH_TOLS(ITOL) = 1

!     Enter the linear algebra type loop.
!        1 = dense Jacobian
!        2 = sparse Jacobian
!        3 = banded Jacobian
      DO JTYPE = 1, 3

!     Enter the constant Jacobian loop.
!        1 = Jacobian is not constant
!        2 = Jacobian is constant for classes 0 and 1
!     Make this a one-trip loop (1 to 1 or 2 to 2) for stats to be meaningful:
      DO CONSTANTJ = 1, 1

!     Enter the numerical or analytical Jacobian loop.
!        1 = numerical Jacobian
!        2 = analytical Jacobian
!     Make this a one-trip loop (1 to 1 or 2 to 2) for stats to be meaningful:
      DO IJAC = 1, 1
      ANALYTIC_JACOBIAN = .TRUE.
      IF (IJAC == 1) ANALYTIC_JACOBIAN = .FALSE.

      IF (JTYPE==1 .AND. IJACSP==1) METHOD = 1
      IF (JTYPE==2 .AND. IJACSP==1) METHOD = 2
      IF (JTYPE==3 .AND. IJACSP==1) METHOD = 3
      IF (JTYPE==1 .AND. IJACSP==2) METHOD = 4
      IF (JTYPE==2 .AND. IJACSP==2) METHOD = 5
      IF (JTYPE==3 .AND. IJACSP==2) METHOD = 6

!     Enter the problem test loop.
      DO ITEST = 1, 55

!       JACS supplies the full matrix in the analytic Jacobian case.
        IF (IJAC == 2 .AND. ISTR == 1) GOTO 20
        ID = MYID(ITEST)
        IF (ID == 0) GOTO 20
        WHICH_PROBS(ITEST) = 1
        IID = MOD(ID,10)
        WRITE (6,90010)
        CLASS = ID/10
        PROBLEM = ID - 10*CLASS
        J_IS_CONSTANT = .FALSE.
        IF ((CLASS==0 .OR. CLASS==1) .AND. CONSTANTJ==2) J_IS_CONSTANT = .TRUE.

!       Output header.
        IF (CLASS == 0) THEN
          WRITE (6,90000) PROBLEM
        ELSE IF (CLASS == 1) THEN
          WRITE (6,90001) PROBLEM
        ELSE IF (CLASS == 2) THEN
          WRITE (6,90002) PROBLEM
        ELSE IF (CLASS == 3) THEN
          WRITE (6,90003) PROBLEM
        ELSE IF (CLASS == 4) THEN
          WRITE (6,90004) PROBLEM
        ELSE IF (CLASS == 5) THEN
          WRITE (6,90005) PROBLEM
        END IF
        WRITE(6,*)    ' Results for ITOL = ', ITOL, ' follow.'
        PRINT *,      ' Results for ITOL = ', ITOL, ' follow.'
        IF (IJACSP == 1) THEN
           WRITE(6,*) ' Results for JACSP Jacobian follow.'
           PRINT *,   ' Results for JACSP Jacobian follow.'
        ELSE
           WRITE(6,*) ' Results for VODE Jacobian follow.'
           PRINT *,   ' Results for VODE Jacobian follow.'
        END IF
        IF (ISCALE == 1) THEN
           WRITE(6,*) ' The ODEs will be scaled.'
           PRINT *,   ' The ODEs will be scaled.'
        ELSE
           WRITE(6,*) ' The ODEs will not be scaled.'
           PRINT *,   ' The ODEs will not be scaled.'
        END IF
        IF (ISTR == 1) THEN
           WRITE(6,*) ' Results for numerically determined sparse pointers follow.'
           PRINT *,   ' Results for numerically determined sparse pointers follow.'
        ELSE
           WRITE(6,*) ' Results for user supplied sparse pointers follow.'
           PRINT *,   ' Results for user supplied sparse pointers follow.'
        END IF
        IF (JTYPE == 1) THEN
           WRITE(6,*) ' Results for dense Jacobian option follow.'
           PRINT *,   ' Results for dense Jacobian option follow.'
        END IF
        IF (JTYPE == 2) THEN
           WRITE(6,*) ' Results for sparse Jacobian option follow.'
           PRINT *,   ' Results for sparse Jacobian option follow.'
        END IF
        IF (JTYPE == 3) THEN
           WRITE(6,*) ' Results for banded Jacobian option follow.'
           PRINT *,   ' Results for banded Jacobian option follow.'
        END IF
        IF (CONSTANTJ == 1) THEN
           WRITE(6,*) ' Results for non-constant Jacobian option follow.'
           PRINT *,   ' Results for non-constant Jacobian option follow.'
        ELSE
           WRITE(6,*) ' Results for constant Jacobian option follow.'
           PRINT *,   ' Results for constant Jacobian option follow.'
        END IF
        IF (IJAC == 1) THEN
           WRITE(6,*) ' Results for numerical Jacobian option follow.'
           PRINT *,   ' Results for numerical Jacobian option follow.'
        ELSE
           WRITE(6,*) ' Results for analytical Jacobian option follow.'
           PRINT *,   ' Results for analytical Jacobian option follow.'
        END IF

!       Initialize the problem.
!       Scale the ODEs?
!         USEW = .TRUE.
          IWT = -1
          IF (USEW) IWT = 1
!       Get the problem information.
          CALL IVALU(TBEGIN,TEND,HBEGIN,HBOUND,YINIT)
!       Use the IVALU starting step size?
          USEHBEGIN = .TRUE.
          IF (.NOT. USEHBEGIN) HBEGIN = 0.0D0
!       Define the number of ODEs.
          NEQ = N
!       Use the full matrix for the banded solution.
          ML = NEQ - 1
          MU = NEQ - 1
!       Initial and final integration times.
          T = TBEGIN
          TOUT = TEND
!       Initial solution.
          Y(1:NEQ) = YINIT(1:NEQ)
!       Error tolerances.
          EPS = 1.0D0 / 10.0D0**ITOL
          RELERR_TOLERANCES(1:NEQ) = EPS
          ABSERR_TOLERANCES(1:NEQ) = EPS
          WRITE (6,90007) ID,TBEGIN,TEND,HBEGIN,HBOUND,IWT,N,EPS,Y(1:NEQ)
          BOUNDS = .FALSE.
!         Force nonnegativity for the chemical kinetics problems in Class F.
          IF (ITEST >= 51) THEN
             BOUNDS = .TRUE.
             COMPS(1) = 1
             COMPS(2) = 2
             COMPS(3) = 3
             DO I = 1, NEQ
                COMPS(I) = I
             END DO
             LBOUNDS(1:NEQ) = 0.0D0
             UBOUNDS(1:NEQ) = 1.0D50
          END IF
        IF (BOUNDS) WRITE(6,90040)
!       Initialize the integration flags.
          ITASK = 1
          ISTATE = 1

!       Use the full matrix for sparsity structure.
        IF (JTYPE == 2 .AND. SUPPLY_STRUCTURE) THEN
           IADIM = NEQ + 1
           JADIM = NEQ * NEQ
           IA(1) = 1
           DO I = 1, NEQ
             IA(I+1) = IA(I) + NEQ
           END DO
           I = 0
           DO COL = 1, NEQ
             I = NEQ*(COL-1)
             DO ROW = 1, NEQ
               I = I + 1
               JA(I) = ROW
             END DO
           END DO
        END IF

        CALL CPU_TIME(DVTIME1)

!       Set the integration options.
        IF (JTYPE == 1) THEN
!          Dense solution ...
           IF (TOLVEC) THEN
!             Supply vectors of tolerances.
              IF (BOUNDS) THEN
                 OPTIONS = SET_OPTS(DENSE_J=.TRUE.,                          &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR_VECTOR=RELERR_TOLERANCES(1:NEQ),     &
                                 ABSERR_VECTOR=ABSERR_TOLERANCES(1:NEQ),     &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 CONSTRAINED=COMPS(1:NEQ),                   &
                                 CLOWER=LBOUNDS(1:NEQ),CUPPER=UBOUNDS(1:NEQ))
              ELSE
                 OPTIONS = SET_OPTS(DENSE_J=.TRUE.,                          &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR_VECTOR=RELERR_TOLERANCES(1:NEQ),     &
                                 ABSERR_VECTOR=ABSERR_TOLERANCES(1:NEQ),     &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP)
              END IF
           ELSE
!             Supply scalar tolerances.
              IF (BOUNDS) THEN
                 OPTIONS = SET_OPTS(DENSE_J=.TRUE.,                          &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR=EPS,ABSERR=EPS,                      &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 CONSTRAINED=COMPS(1:NEQ),                   &
                                 CLOWER=LBOUNDS(1:NEQ),CUPPER=UBOUNDS(1:NEQ))
              ELSE
                 OPTIONS = SET_OPTS(DENSE_J=.TRUE.,                          &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR=EPS,ABSERR=EPS,                      &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP)
              END IF
           END IF
        END IF

        IF (JTYPE == 2) THEN
!          Sparse solution ...
           IF (TOLVEC) THEN
!             Supply vectors of tolerances.
              IF (BOUNDS) THEN
                 OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,                         &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR_VECTOR=RELERR_TOLERANCES(1:NEQ),     &
                                 ABSERR_VECTOR=ABSERR_TOLERANCES(1:NEQ),     &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 USER_SUPPLIED_SPARSITY=SUPPLY_STRUCTURE,    &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 MA28_RPS=.TRUE.,CONSTRAINED=COMPS(1:NEQ),   &
                                 CLOWER=LBOUNDS(1:NEQ),CUPPER=UBOUNDS(1:NEQ))
              ELSE
                 OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,                         &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR_VECTOR=RELERR_TOLERANCES(1:NEQ),     &
                                 ABSERR_VECTOR=ABSERR_TOLERANCES(1:NEQ),     &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 USER_SUPPLIED_SPARSITY=SUPPLY_STRUCTURE,    &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 MA28_RPS=.TRUE.)
              END IF
           ELSE
!             Supply scalar tolerances.
              IF (BOUNDS) THEN
                 OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,                         &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR=EPS,ABSERR=EPS,                      &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 USER_SUPPLIED_SPARSITY=SUPPLY_STRUCTURE,    &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 MA28_RPS=.TRUE.,CONSTRAINED=COMPS(1:NEQ),   &
                                 CLOWER=LBOUNDS(1:NEQ),CUPPER=UBOUNDS(1:NEQ))
              ELSE
                 OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,                         &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR=EPS,ABSERR=EPS,                      &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 USER_SUPPLIED_SPARSITY=SUPPLY_STRUCTURE,    &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 MA28_RPS=.TRUE.)
              END IF
           END IF
        END IF

        IF (JTYPE == 3) THEN
!          Banded solution ...
           IF (TOLVEC) THEN
!             Supply vectors of tolerances.
              IF (BOUNDS) THEN
                 OPTIONS = SET_OPTS(BANDED_J=.TRUE.,                         &
                                 LOWER_BANDWIDTH=ML,UPPER_BANDWIDTH=MU,      &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR_VECTOR=RELERR_TOLERANCES(1:NEQ),     &
                                 ABSERR_VECTOR=ABSERR_TOLERANCES(1:NEQ),     &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 CONSTRAINED=COMPS(1:NEQ),                   &
                                 CLOWER=LBOUNDS(1:NEQ),CUPPER=UBOUNDS(1:NEQ))
              ELSE
                 OPTIONS = SET_OPTS(BANDED_J=.TRUE.,                         &
                                 LOWER_BANDWIDTH=ML,UPPER_BANDWIDTH=MU,      &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR_VECTOR=RELERR_TOLERANCES(1:NEQ),     &
                                 ABSERR_VECTOR=ABSERR_TOLERANCES(1:NEQ),     &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP)
              END IF
           ELSE
!             Supply scalar tolerances.
              IF (BOUNDS) THEN
                 OPTIONS = SET_OPTS(BANDED_J=.TRUE.,                         &
                                 LOWER_BANDWIDTH=ML,UPPER_BANDWIDTH=MU,      &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR=EPS,ABSERR=EPS,                      &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP,                &
                                 CONSTRAINED=COMPS(1:NEQ),                   &
                                 CLOWER=LBOUNDS(1:NEQ),CUPPER=UBOUNDS(1:NEQ))
              ELSE
                 OPTIONS = SET_OPTS(BANDED_J=.TRUE.,                         &
                                 LOWER_BANDWIDTH=ML,UPPER_BANDWIDTH=MU,      &
                                 USER_SUPPLIED_JACOBIAN=ANALYTIC_JACOBIAN,   &
                                 RELERR=EPS,ABSERR=EPS,                      &
                                 MXSTEP=100000,H0=HBEGIN,HMAX=HBOUND,        &
                                 CONSTANT_JACOBIAN=J_IS_CONSTANT,            &
                                 JACOBIAN_BY_JACSP=USE_JACSP)
              END IF
           END IF
        END IF

!       Supply the sparse structure arrays directly.
        IF (JTYPE==2 .AND. SUPPLY_STRUCTURE) THEN
           CALL USERSETS_IAJA(IA(1:IADIM),IADIM,JA(1:JADIM),JADIM)
        END IF

!       Perform the integration.
        IF (JTYPE == 1) THEN
           CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JACD)
        END IF
        IF (JTYPE == 2) THEN
           CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JACS)
        END IF
        IF (JTYPE == 3) THEN
           CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JACB)
        END IF

        CALL CPU_TIME(DVTIME2)
        DVTIME = DVTIME2 - DVTIME1
        EXECUTION_TIMES(ITOL) = EXECUTION_TIMES(ITOL) + DVTIME

        Total_Solutions = Total_Solutions + 1

!       Check if an error occurred.
        IF (ISTATE < 0) THEN
          WRITE (6,90006) ISTATE
          ERRORS = .TRUE.
        END IF

!       Gather the integration statistics.
        CALL GET_STATS(RSTATS,ISTATS)
        WRITE (6,90009) ISTATS(11), ISTATS(12), ISTATS(13)
        NUM_S(ITOL,ITEST,JTYPE) = ISTATS(11)
        NUM_D(ITOL,ITEST,JTYPE) = ISTATS(12)
        NUM_J(ITOL,ITEST,JTYPE) = ISTATS(13)

!       Calculate the error in the final solution.
        CALL EVALU(YFINAL)
        DO I = 1, NEQ
          AERROR(I) = ABS(Y(I)-YFINAL(I))
        END DO
        ERRMAX = MAXVAL(AERROR(1:NEQ))
        ERSTATS(ITOL,ITEST,JTYPE) = ERRMAX
        ERSTATS2(ISCALE,ITOL,ITEST,METHOD) = ERRMAX
        COUNTS(ISCALE,ITOL,ITEST,METHOD,1) = ISTATS(11)
        COUNTS(ISCALE,ITOL,ITEST,METHOD,2) = ISTATS(12)
        COUNTS(ISCALE,ITOL,ITEST,METHOD,3) = ISTATS(13)
        WRITE (6,90008) (I,Y(I),YFINAL(I),AERROR(I),I=1,NEQ)

!       Release the work arrays and determine how much storage was required.
        CALL RELEASE_ARRAYS

20    CONTINUE
      END DO ! End of ITEST Loop
      END DO ! End of IJAC Loop
      END DO ! End of JTYPE Loop
      END DO ! End of CONSTANTJ Loop
      END DO ! End of ITOL Loop
      END DO ! End of ISTR Loop

!     Print the maximum errors for this value of IJACSP.
      DO ITOL = 5, 12, 1
         WRITE(6,90016) ITOL
         WRITE(6,90020)
         IF (IJACSP == 1) WRITE(6,90023)
         IF (IJACSP == 2) WRITE(6,90022)
         IF (ISCALE == 1) WRITE(6,90030)
         IF (ISCALE == 2) WRITE(6,90031)
         WRITE(6,90021)
         WRITE(6,90015)
         DO ITEST = 1, 55
            IF (WHICH_PROBS(ITEST) == 1) THEN
               ID = MYID(ITEST)
               IF (ID /= 0) WRITE(6,90014) ITOL, ITEST,         &
                            (ERSTATS(ITOL,ITEST,JTYPE),JTYPE=1,3)
            END IF
         END DO
      END DO 

!     Print the integration stats for this value of IJACSP.
      DO ITOL = 5, 12, 1
         WRITE(6,90016) ITOL
         WRITE(6,90024)
         IF (IJACSP == 1) WRITE(6,90023)
         IF (IJACSP == 2) WRITE(6,90022)
         IF (ISCALE == 1) WRITE(6,90030)
         IF (ISCALE == 2) WRITE(6,90031)
         WRITE(6,90029)
         DO ITEST = 1, 55
            IF (WHICH_PROBS(ITEST) == 1) THEN
               ID = MYID(ITEST)
               IF (ID /= 0) THEN
                  WRITE(6,90026) ITEST, (NUM_S(ITOL,ITEST,JTYPE),JTYPE=1,3)
                  WRITE(6,90027)        (NUM_J(ITOL,ITEST,JTYPE),JTYPE=1,3)
                  WRITE(6,90028)        (NUM_D(ITOL,ITEST,JTYPE),JTYPE=1,3)
                  WRITE(6,*) ' '
               END IF
            END IF
         END DO
      END DO 

      END DO ! End of IJACSP Loop
      END DO ! End of ISCALE Loop

      WRITE(6,90045)

      DO ITEST = 1, 55
         IF (WHICH_PROBS(ITEST) == 1) THEN
            ID = MYID(ITEST)
            IID = MOD(ID,10)
            WRITE (6,90010)
            CLASS = ID/10
            PROBLEM = ID - 10*CLASS
            IF (CLASS == 0) THEN
              WRITE (6,90000) PROBLEM
            ELSE IF (CLASS == 1) THEN
              WRITE (6,90001) PROBLEM
            ELSE IF (CLASS == 2) THEN
              WRITE (6,90002) PROBLEM
            ELSE IF (CLASS == 3) THEN
              WRITE (6,90003) PROBLEM
            ELSE IF (CLASS == 4) THEN
              WRITE (6,90004) PROBLEM
            ELSE IF (CLASS == 5) THEN
              WRITE (6,90005) PROBLEM
            END IF
            DO ITOL = 5, 12, 1
               IF (WHICH_TOLS(ITOL) == 1) THEN
                  DO ISCALE = 1,  2
                     IF (ISCALE == 1) WRITE(6,90035) ITEST, ITOL
                     IF (ISCALE == 2) WRITE(6,90036) ITEST, ITOL
                     WRITE(6,90037) ERSTATS2(ISCALE,ITOL,ITEST,1), &
                                    ERSTATS2(ISCALE,ITOL,ITEST,4)
                     WRITE(6,90038) ERSTATS2(ISCALE,ITOL,ITEST,2), &
                                    ERSTATS2(ISCALE,ITOL,ITEST,5)
                     WRITE(6,90039) ERSTATS2(ISCALE,ITOL,ITEST,3), &
                                    ERSTATS2(ISCALE,ITOL,ITEST,6)
                     WRITE(6,90044)
                     WRITE(6,90041) COUNTS(ISCALE,ITOL,ITEST,1,1), &
                                    COUNTS(ISCALE,ITOL,ITEST,4,1), &
                                    COUNTS(ISCALE,ITOL,ITEST,1,2), &
                                    COUNTS(ISCALE,ITOL,ITEST,4,2), &
                                    COUNTS(ISCALE,ITOL,ITEST,1,3), &
                                    COUNTS(ISCALE,ITOL,ITEST,4,3)
                     WRITE(6,90042) COUNTS(ISCALE,ITOL,ITEST,2,1), &
                                    COUNTS(ISCALE,ITOL,ITEST,5,1), &
                                    COUNTS(ISCALE,ITOL,ITEST,2,2), &
                                    COUNTS(ISCALE,ITOL,ITEST,5,2), &
                                    COUNTS(ISCALE,ITOL,ITEST,2,3), &
                                    COUNTS(ISCALE,ITOL,ITEST,5,3)
                     WRITE(6,90043) COUNTS(ISCALE,ITOL,ITEST,3,1), &
                                    COUNTS(ISCALE,ITOL,ITEST,6,1), &
                                    COUNTS(ISCALE,ITOL,ITEST,3,2), &
                                    COUNTS(ISCALE,ITOL,ITEST,6,2), &
                                    COUNTS(ISCALE,ITOL,ITEST,3,3), &
                                    COUNTS(ISCALE,ITOL,ITEST,6,3)
                     WRITE (6,90010)
                  END DO
               END IF
            END DO
         END IF
      END DO

!     Total number of solutions performed.
      WRITE(6,90025) Total_Solutions
      WRITE (6,90010)

!     Execution Times.
      WRITE(6,90033)
      DO ITOL = 5, 12, 1
         ITIME = INT(EXECUTION_TIMES(ITOL))
         WRITE(6,90034) ITOL, ITIME
      END DO
      WRITE (6,90010)

!     Where to search in output file.
      WRITE(6,90032)

!     Indicate whether any DVODE_F90 made any error returns.
      IF (ERRORS) THEN
         WRITE(6,90011)
      ELSE
         WRITE(6,90012)
      END IF

!     Format statements follow.
90000 FORMAT (' Class/Problem = A',I1)
90001 FORMAT (' Class/Problem = B',I1)
90002 FORMAT (' Class/Problem = C',I1)
90003 FORMAT (' Class/Problem = D',I1)
90004 FORMAT (' Class/Problem = E',I1)
90005 FORMAT (' Class/Problem = F',I1)
90006 FORMAT (' An error occurred in VODE_F90. ISTATE = ',I3)
90007 FORMAT (' Problem ID       = ',I3,/,   ' Initial time     = ',D15.5,/, &
              ' Final time       = ',D15.5,/,' Initial stepsize = ',D15.5,/, &
              ' Maximum stepsize = ',D15.5,/,' IWT flag         = ',I3,/,    &
              ' Number of ODEs   = ',I3,/,   ' Error tolerance  = ',D15.5,/, &
              ' Initial solution = ',/,(D15.5))
90008 FORMAT (' Computed and reference solutions and absolute' &
              ' errors follow:',/,(I3,3D15.5))
90009 FORMAT (' Steps = ',I10,' f-s = ',I10,' J-s = ',I10)
90010 FORMAT (' ____________________________________________________________')
90011 FORMAT (/,' Errors occurred. (Search on ISTATE.)')
90012 FORMAT (/,' No errors occurred.')
90013 FORMAT (' For ITOL = ', I3,' Dense  max. observed error = ',D15.5)
90014 FORMAT (I5,6X,I5,3D15.5)
90015 FORMAT (/,'Tolerance',3X,'Problem',4X,'Dense',10X,'Sparse',9X,'Banded')
90016 FORMAT (//,20X,' Results for ITOL = ',I3)
90017 FORMAT (' For ITOL = ', I3,' Sparse max. observed error = ',D15.5)
90018 FORMAT (' For ITOL = ', I3,' Banded max. observed error = ',D15.5)
90019 FORMAT (' The dense and banded results are not identical.')
90020 FORMAT (/,' Maximum Errors for All Problems Follow:')
90021 FORMAT (/,10X,' Maximum Errors for Individual Problems Follow:')
90022 FORMAT (' The original VODE Jacobian was used.')
90023 FORMAT (' The JACSP Jacobian was used.')
90024 FORMAT (/,' Integration Statistics Follow:')
90025 FORMAT (/,' Total number of solutions performed = ',I5)
90026 FORMAT (I5,'   Steps: ',3I10)
90027 FORMAT (5X,'   Jacs:  ',3I10)
90028 FORMAT (5X,'   Ders:  ',3I10)
90029 FORMAT (/,' Problem',12X,'Dense',4X,'Sparse',4X,'Banded')
90030 FORMAT (' The ODEs were scaled.')
90031 FORMAT (' The ODEs were not scaled.')
90032 FORMAT (//,' To locate the summary results, search on the phrase',/ &
                 ' Summary Statistics Follow in the output file',         &
                 ' (stiffoptions.dat).',//)
90033 FORMAT (/,' Tolerance', 10X, 'Total Execution Time')
90034 FORMAT ((4X,I5,10X,I10))
90035 FORMAT (/,' Problem No. = ',I3,/,' ITOL = ',I3,/,' Scaled ODEs', &
              13X,'Absolute Errors')
90036 FORMAT (' Problem No. = ',I3,/,' ITOL = ',I3,/,' Unscaled ODEs', &
              18X,'Absolute Errors')
90037 FORMAT (' Dense : JACSP - VODE: ',2D15.5)
90038 FORMAT (' Sparse: JACSP - VODE: ',2D15.5)
90039 FORMAT (' Banded: JACSP - VODE: ',2D15.5)
90040 FORMAT (' Nonnegativity will be enforced for this problem.')
90041 FORMAT (' Dense  SDJ: ',6I7)
90042 FORMAT (' Sparse SDJ: ',6I7)
90043 FORMAT (' Banded SDJ: ',6I7)
90044 FORMAT (19X,'J',6X,'V',6X,'J',6X,'V',6X,'J',6X,'V')
90045 FORMAT (///,20X,'Summary Statistics Follow')

      STOP
    END PROGRAM DEMOSTIFF
