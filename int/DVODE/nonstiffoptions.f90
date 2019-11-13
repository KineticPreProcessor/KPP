! DVODE_F90 demonstration program

! Toronto Nonstiff Test Set

! Reference:
! ALGORITHM 648, COLLECTED ALGORITHMS FROM ACM.
! TRANSACTIONS ON MATHEMATICAL SOFTWARE,
! VOL. 13, NO. 1, P. 28

! Caution:
! This version of the test set is intended for use only
! in conjunction with this demo program for DVODE_F90.
! The original test routines are altered in this version.
! Some subroutine arguments, COMMON variables, and local
! variables have been moved to the following global header.
! You should obtain the original test suite from netlib
! if you plan to use the routines for another solver.

!******************************************************************

    MODULE NONSTIFFSET

! Note:
! In this version ID is defined in the main program
! demostiff. The other parameters are obtained by
! calling IVALU with a modified argument list.

      IMPLICIT NONE
      INTEGER IWT, N, ID, NFCN
      DOUBLE PRECISION W
      DIMENSION W(51)

    CONTAINS

      SUBROUTINE DERIVS(NEQ,T,Y,YDOT)
        IMPLICIT NONE
        INTEGER NEQ
        DOUBLE PRECISION T, Y, YDOT
        DIMENSION Y(NEQ), YDOT(NEQ)

        CALL FCN(T,Y,YDOT)
        RETURN
      END SUBROUTINE DERIVS

      SUBROUTINE IVALU(XSTART,XEND,HBEGIN,HMAX,Y)

!      ROUTINE TO PROVIDE THE INITIAL VALUES REQUIRED TO SPECIFY
!      THE MATHEMATICAL PROBLEM AS WELL AS VARIOUS PROBLEM
!      PARAMETERS REQUIRED BY THE TESTING PACKAGE. THE APPROPRIATE
!      SCALING VECTOR IS ALSO INITIALISED IN CASE THIS OPTION IS
!      SELECTED.

!      PARAMETERS (OUTPUT)
!      N      - DIMENSION OF THE PROBLEM
!      XSTART - INITIAL VALUE OF THE INDEPENDENT VARIABLE
!      XEND   - FINAL VALUE OF THE INDEPENDENT VARIABLE
!      HBEGIN - APPROPRIATE STARTING STEPSIZE
!      Y      - VECTOR OF INITIAL CONDITIONS FOR THE DEPENDENT
!               VARIABLES
!      WT     - VECTOR OF WEIGHTS USED TO SCALE THE PROBLEM IF
!               THIS OPTION IS SELECTED.

!      PARAMETER  (INPUT)
!      IWT    - FLAG TO INDICATE IF SCALED OPTION IS SELESTED
!      ID     - FLAG IDENTIFYING WHICH EQUATION IS BEING SOLVED

        IMPLICIT NONE

!     .. Scalar Arguments ..
        DOUBLE PRECISION HBEGIN, HMAX, XEND, XSTART
!     INTEGER ID, IWT, N
!     .. Array Arguments ..
        DOUBLE PRECISION Y(51)
!     DOUBLE PRECISION W(51)
!     .. Local Scalars ..
        DOUBLE PRECISION E, HB, HM, XE, XS
        INTEGER I, IOUT
!     .. Data statements ..
        DATA HM, HB, XS, XE/20.0D0, 1.0D0, 0.0D0, 20.0D0/

!     .. Executable Statements ..

        HMAX = HM
        HBEGIN = HB
        XSTART = XS
        XEND = XE
!     GOTO (40, 60, 80, 100, 120, 20, 20, 20, 20, 20, 140, 160, 180,    &
!     200, 220, 20, 20, 20, 20, 20, 240, 280, 320, 360, 400, 20, 20, 20,&
!     20, 20, 420, 420, 420, 420, 420, 20, 20, 20, 20, 20, 540, 560,    &
!     580, 600, 620, 20, 20, 20, 20, 20, 640, 660, 680, 700, 720) ID
        GO TO (20,30,40,50,60,10,10,10,10,10,70,80,90,100,110,10,10,10,10,10, &
          120,130,140,150,160,10,10,10,10,10,170,170,170,170,170,10,10,10,10, &
          10,230,240,250,260,270,10,10,10,10,10,280,290,300,310,320) ID
10      IOUT = 6
        WRITE (IOUT,FMT=90000) ID
        STOP

!     PROBLEM CLASS A

!     1:
20      CONTINUE
!     PROBLEM A1
        N = 1
        W(1) = 0.100D+01
        Y(1) = 1.0D0
        GO TO 330
!     2:
30      CONTINUE
!     PROBLEM A2
        N = 1
        W(1) = 0.100D+01
        Y(1) = 1.0D0
        GO TO 330
!     3:
40      CONTINUE
!     PROBLEM A3
        N = 1
        W(1) = 0.271D+01
        Y(1) = 1.0D0
        GO TO 330
!     4:
50      CONTINUE
!     PROBLEM A4
        N = 1
        W(1) = 0.177D+02
        Y(1) = 1.0D0
        GO TO 330
!     5:
60      CONTINUE
!     PROBLEM A5
        N = 1
        W(1) = 0.620D+01
        Y(1) = 4.0D0
        GO TO 330

!     PROBLEM CLASS B

!     11:
70      CONTINUE
!     PROBLEM B1
        N = 2
        W(1) = 0.425D+01
        W(2) = 0.300D+01
        Y(1) = 1.0D0
        Y(2) = 3.0D0
        GO TO 330
!     12:
80      CONTINUE
!     PROBLEM B2
        N = 3
        W(1) = 0.200D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        Y(1) = 2.0D0
        Y(2) = 0.0D0
        Y(3) = 1.0D0
        GO TO 330
!     13:
90      CONTINUE
!     PROBLEM B3
        N = 3
        W(1) = 0.100D+01
        W(2) = 0.519D+00
        W(3) = 0.947D+00
        Y(1) = 1.0D0
        Y(2) = 0.0D0
        Y(3) = 0.0D0
        GO TO 330
!     14:
100     CONTINUE
!     PROBLEM B4
        N = 3
        W(1) = 0.300D+01
        W(2) = 0.220D+01
        W(3) = 0.100D+01
        Y(1) = 3.0D0
        Y(2) = 0.0D0
        Y(3) = 0.0D0
        GO TO 330
!     15:
110     CONTINUE
!     PROBLEM B5
        N = 3
        W(1) = 0.100D+01
        W(2) = 0.100D+01
        W(3) = 0.100D+01
        Y(1) = 0.0D0
        Y(2) = 1.0D0
        Y(3) = 1.0D0
        GO TO 330

!     PROBLEM CLASS C

!     21:
120     CONTINUE
!     PROBLEM C1
        N = 10
        W(1) = 0.100D+01
        W(2) = 0.368D+00
        W(3) = 0.271D+00
        W(4) = 0.224D+00
        W(5) = 0.195D+00
        W(6) = 0.175D+00
        W(7) = 0.161D+00
        W(8) = 0.149D+00
        W(9) = 0.139D+00
        W(10) = 0.998D+00
        Y(1) = 1.0D0
        DO I = 2, N
          Y(I) = 0.0D0
        END DO
        GO TO 330
!     22:
130     CONTINUE
!     PROBLEM C2
        N = 10
        W(1) = 0.100D+01
        W(2) = 0.250D+00
        W(3) = 0.148D+00
        W(4) = 0.105D+00
        W(5) = 0.818D-01
        W(6) = 0.669D-01
        W(7) = 0.566D-01
        W(8) = 0.491D-01
        W(9) = 0.433D-01
        W(10) = 0.100D+01
        Y(1) = 1.0D0
        DO I = 2, N
          Y(I) = 0.0D0
        END DO
        GO TO 330
!     23:
140     CONTINUE
!     PROBLEM C3
        N = 10
        W(1) = 0.100D+01
        W(2) = 0.204D+00
        W(3) = 0.955D-01
        W(4) = 0.553D-01
        W(5) = 0.359D-01
        W(6) = 0.252D-01
        W(7) = 0.184D-01
        W(8) = 0.133D-01
        W(9) = 0.874D-02
        W(10) = 0.435D-02
        Y(1) = 1.0D0
        DO I = 2, N
          Y(I) = 0.0D0
        END DO
        GO TO 330
!     24:
150     CONTINUE
!     PROBLEM C4
        N = 51
        W(1) = 0.100D+01
        W(2) = 0.204D+00
        W(3) = 0.955D-01
        W(4) = 0.553D-01
        W(5) = 0.359D-01
        W(6) = 0.252D-01
        W(7) = 0.186D-01
        W(8) = 0.143D-01
        W(9) = 0.113D-01
        W(10) = 0.918D-02
        W(11) = 0.760D-02
        W(12) = 0.622D-02
        W(13) = 0.494D-02
        W(14) = 0.380D-02
        W(15) = 0.284D-02
        W(16) = 0.207D-02
        W(17) = 0.146D-02
        W(18) = 0.101D-02
        W(19) = 0.678D-03
        W(20) = 0.444D-03
        W(21) = 0.283D-03
        W(22) = 0.177D-03
        W(23) = 0.107D-03
        W(24) = 0.637D-04
        W(25) = 0.370D-04
        W(26) = 0.210D-04
        W(27) = 0.116D-04
        W(28) = 0.631D-05
        W(29) = 0.335D-05
        W(30) = 0.174D-05
        W(31) = 0.884D-06
        W(32) = 0.440D-06
        W(33) = 0.215D-06
        W(34) = 0.103D-06
        W(35) = 0.481D-07
        W(36) = 0.221D-07
        W(37) = 0.996D-08
        W(38) = 0.440D-08
        W(39) = 0.191D-08
        W(40) = 0.814D-09
        W(41) = 0.340D-09
        W(42) = 0.140D-09
        W(43) = 0.564D-10
        W(44) = 0.224D-10
        W(45) = 0.871D-11
        W(46) = 0.334D-11
        W(47) = 0.126D-11
        W(48) = 0.465D-12
        W(49) = 0.169D-12
        W(50) = 0.600D-13
        W(51) = 0.189D-13
        Y(1) = 1.0D0
        DO I = 2, N
          Y(I) = 0.0D0
        END DO
        GO TO 330
!     25:
160     CONTINUE
!     PROBLEM C5
        N = 30
        W(1) = 0.545D+01
        W(2) = 0.471D+01
        W(3) = 0.203D+01
        W(4) = 0.664D+01
        W(5) = 0.834D+01
        W(6) = 0.346D+01
        W(7) = 0.113D+02
        W(8) = 0.172D+02
        W(9) = 0.748D+01
        W(10) = 0.302D+02
        W(11) = 0.411D+01
        W(12) = 0.144D+01
        W(13) = 0.244D+02
        W(14) = 0.284D+02
        W(15) = 0.154D+02
        W(16) = 0.764D+00
        W(17) = 0.661D+00
        W(18) = 0.284D+00
        W(19) = 0.588D+00
        W(20) = 0.366D+00
        W(21) = 0.169D+00
        W(22) = 0.388D+00
        W(23) = 0.190D+00
        W(24) = 0.877D-01
        W(25) = 0.413D-01
        W(26) = 0.289D+00
        W(27) = 0.119D+00
        W(28) = 0.177D+00
        W(29) = 0.246D+00
        W(30) = 0.319D-01
        Y(1) = 3.42947415189D0
        Y(2) = 3.35386959711D0
        Y(3) = 1.35494901715D0
        Y(4) = 6.64145542550D0
        Y(5) = 5.97156957878D0
        Y(6) = 2.18231499728D0
        Y(7) = 11.2630437207D0
        Y(8) = 14.6952576794D0
        Y(9) = 6.27960525067D0
        Y(10) = -30.1552268759D0
        Y(11) = 1.65699966404D0
        Y(12) = 1.43785752721D0
        Y(13) = -21.1238353380D0
        Y(14) = 28.4465098142D0
        Y(15) = 15.3882659679D0
        Y(16) = -.557160570446D0
        Y(17) = .505696783289D0
        Y(18) = .230578543901D0
        Y(19) = -.415570776342D0
        Y(20) = .365682722812D0
        Y(21) = .169143213293D0
        Y(22) = -.325325669158D0
        Y(23) = .189706021964D0
        Y(24) = .0877265322780D0
        Y(25) = -.0240476254170D0
        Y(26) = -.287659532608D0
        Y(27) = -.117219543175D0
        Y(28) = -.176860753121D0
        Y(29) = -.216393453025D0
        Y(30) = -.0148647893090D0
        GO TO 330

!     PROBLEM CLASS D

170     CONTINUE
!     PROBLEM D1, D2, D3, D4, D5
        E = .2D0*(REAL(ID)-3.D1) - .1D0
        N = 4
        IF (ID/=31) GO TO 180
        W(1) = 0.110D+01
        W(2) = 0.995D+00
        W(3) = 0.101D+01
        W(4) = 0.111D+01
180     IF (ID/=32) GO TO 190
        W(1) = 0.130D+01
        W(2) = 0.954D+00
        W(3) = 0.105D+01
        W(4) = 0.136D+01
190     IF (ID/=33) GO TO 200
        W(1) = 0.150D+01
        W(2) = 0.866D+00
        W(3) = 0.115D+01
        W(4) = 0.173D+01
200     IF (ID/=34) GO TO 210
        W(1) = 0.170D+01
        W(2) = 0.714D+00
        W(3) = 0.140D+01
        W(4) = 0.238D+01
210     IF (ID/=35) GO TO 220
        W(1) = 0.190D+01
        W(2) = 0.436D+00
        W(3) = 0.229D+01
        W(4) = 0.436D+01
220     CONTINUE
        Y(1) = 1.0D0 - E
        Y(2) = 0.0D0
        Y(3) = 0.0D0
        Y(4) = SQRT((1.0D0+E)/(1.0D0-E))
        GO TO 330

!     PROBLEM CLASS E

!     41:
230     CONTINUE
!     PROBLEM E1
        N = 2
        W(1) = 0.679D+00
        W(2) = 0.478D+00
        E = .79788456080286536D0
        Y(1) = E*.84147098480789651D0
        Y(2) = E*.11956681346419146D0
        GO TO 330
!     42:
240     CONTINUE
!     PROBLEM E2
        N = 2
        W(1) = 0.201D+01
        W(2) = 0.268D+01
        Y(1) = 2.0D0
        Y(2) = 0.0D0
        GO TO 330
!     43:
250     CONTINUE
!     PROBLEM E3
        N = 2
        W(1) = 0.116D+01
        W(2) = 0.128D+01
        Y(1) = 0.0D0
        Y(2) = 0.0D0
        GO TO 330
!     44:
260     CONTINUE
!     PROBLEM E4
        N = 2
        W(1) = 0.340D+02
        W(2) = 0.277D+00
        Y(1) = 3.D1
        Y(2) = 0.D0
        GO TO 330
!     45:
270     CONTINUE
!     PROBLEM E5
        N = 2
        W(1) = 0.141D+02
        W(2) = 0.240D+01
        Y(1) = 0.0D0
        Y(2) = 0.0D0
        GO TO 330

!     PROBLEM CLASS F

!     51:
280     CONTINUE
!     PROBLEM F1
        N = 2
        W(1) = 0.129D+02
        W(2) = 0.384D+02
        Y(1) = 0.0D0
        Y(2) = 0.0D0
        HMAX = 1.0D0
        GO TO 330
!     52:
290     CONTINUE
!     PROBLEM F2
        N = 1
        W(1) = 0.110D+03
        Y(1) = 110.0D0
        HMAX = 1.0D0
        GO TO 330
!     53:
300     CONTINUE
!     PROBLEM F3
        N = 2
        W(1) = 0.131D+01
        W(2) = 0.737D+00
        Y(1) = 0.0D0
        Y(2) = 0.0D0
        HMAX = 1.0D0
        HBEGIN = 0.9D0
        GO TO 330
!     54:
310     CONTINUE
!     PROBLEM F4
        N = 1
        W(1) = 0.152D+01
        Y(1) = 1.0D0
        HMAX = 10.0D0
        GO TO 330
!     55:
320     CONTINUE
!     PROBLEM F5
        N = 1
        W(1) = 0.100D+01
        Y(1) = 1.0D0
        HMAX = 20.0D0
330     CONTINUE
        IF (IWT<0.) GO TO 340
        DO I = 1, N
          Y(I) = Y(I)/W(I)
        END DO
340     CONTINUE
        RETURN

90000   FORMAT ('AN INVALID INTERNAL PROBLEM ID OF ',I4, &
          ' WAS FOUND BY THE IVALU ROUTINE',/ &
          ' RUN TERMINATED. CHECK THE DATA.')
      END SUBROUTINE IVALU

      SUBROUTINE EVALU(Y)

!     ROUTINE TO PROVIDE THE 'TRUE' SOLUTION OF THE DIFFERENTIAL
!     EQUATION EVALUATED AT THE ENDPOINT OF THE INTEGRATION.

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

!     .. Scalar Arguments ..
!     INTEGER ID, IWT, N
!     .. Array Arguments ..
        DOUBLE PRECISION Y(51)
!     DOUBLE PRECISION W(51)
!     .. Local Scalars ..
        INTEGER I

!     .. Executable Statements ..

!     GOTO (20, 40, 60, 80, 100, 620, 620, 620, 620, 620, 120, 140, 160,&
!     180, 200, 620, 620, 620, 620, 620, 220, 240, 260, 280, 300, 620,  &
!     620, 620, 620, 620, 320, 340, 360, 380, 400, 620, 620, 620, 620,  &
!     620, 420, 440, 460, 480, 500, 620, 620, 620, 620, 620, 520, 540,  &
!     560, 580, 600) ID
        GO TO (10,20,30,40,50,310,310,310,310,310,60,70,80,90,100,310,310,310, &
          310,310,110,120,130,140,150,310,310,310,310,310,160,170,180,190,200, &
          310,310,310,310,310,210,220,230,240,250,310,310,310,310,310,260,270, &
          280,290,300) ID
        GO TO 310

!     PROBLEM CLASS A

!     1:
!        PROBLEM A1
10      Y(1) = 2.061153353012535D-09
        GO TO 310
!     2:
!        PROBLEM A2
20      Y(1) = 2.182178902359887D-01
        GO TO 310
!     3:
!        PROBLEM A3
30      Y(1) = 2.491650271850414D+00
        GO TO 310
!     4:
!        PROBLEM A4
40      Y(1) = 1.773016648131483D+01
        GO TO 310
!     5:
!        PROBLEM A5
50      Y(1) = -7.887826688964196D-01
        GO TO 310

!     PROBLEM CLASS B

!     11:
!        PROBLEM B1
60      Y(1) = 6.761876008576667D-01
        Y(2) = 1.860816099640036D-01
        GO TO 310
!     12:
!        PROBLEM B2
70      Y(1) = 1.000000001030576D+00
        Y(2) = 1.000000000000000D+00
        Y(3) = 9.999999989694235D-01
        GO TO 310
!     13:
!        PROBLEM B3
80      Y(1) = 2.061153488557776D-09
        Y(2) = 5.257228022048349D-02
        Y(3) = 9.474277177183630D-01
        GO TO 310
!     14:
!        PROBLEM B4
90      Y(1) = 9.826950928005993D-01
        Y(2) = 2.198447081694832D+00
        Y(3) = 9.129452507276399D-01
        GO TO 310
!     15:
!        PROBLEM B5
100     Y(1) = -9.396570798729192D-01
        Y(2) = -3.421177754000779D-01
        Y(3) = 7.414126596199957D-01
        GO TO 310

!     PROBLEM CLASS C

!     21:
!        PROBLEM C1
110     Y(1) = 2.061153622240064D-09
        Y(2) = 4.122307244619555D-08
        Y(3) = 4.122307244716968D-07
        Y(4) = 2.748204829855288D-06
        Y(5) = 1.374102414941961D-05
        Y(6) = 5.496409659803266D-05
        Y(7) = 1.832136553274552D-04
        Y(8) = 5.234675866508716D-04
        Y(9) = 1.308668966628220D-03
        Y(10) = 9.979127409508656D-01
        GO TO 310
!     22:
!        PROBLEM C2
120     Y(1) = 2.061153577984930D-09
        Y(2) = 2.061153573736588D-09
        Y(3) = 2.061153569488245D-09
        Y(4) = 2.061153565239902D-09
        Y(5) = 2.061153560991560D-09
        Y(6) = 2.061153556743217D-09
        Y(7) = 2.061153552494874D-09
        Y(8) = 2.061153548246532D-09
        Y(9) = 2.061153543998189D-09
        Y(10) = 9.999999814496180D-01
        GO TO 310
!     23:
!        PROBLEM C3
130     Y(1) = 2.948119211022058D-03
        Y(2) = 5.635380154844266D-03
        Y(3) = 7.829072515926013D-03
        Y(4) = 9.348257908594937D-03
        Y(5) = 1.007943610301970D-02
        Y(6) = 9.982674171429909D-03
        Y(7) = 9.088693332766085D-03
        Y(8) = 7.489115195185912D-03
        Y(9) = 5.322964130953349D-03
        Y(10) = 2.762434379029886D-03
        GO TO 310
!     24:
!        PROBLEM C4
140     Y(1) = 3.124111453721466D-03
        Y(2) = 6.015416842150318D-03
        Y(3) = 8.470021834842650D-03
        Y(4) = 1.033682931733337D-02
        Y(5) = 1.153249572873923D-02
        Y(6) = 1.204549525737964D-02
        Y(7) = 1.192957068015293D-02
        Y(8) = 1.128883207111195D-02
        Y(9) = 1.025804501391024D-02
        Y(10) = 8.982017581934167D-03
        Y(11) = 7.597500902492453D-03
        Y(12) = 6.219920556824985D-03
        Y(13) = 4.935916341009131D-03
        Y(14) = 3.801432544256119D-03
        Y(15) = 2.844213677587894D-03
        Y(16) = 2.069123394222672D-03
        Y(17) = 1.464687282843915D-03
        Y(18) = 1.009545263941126D-03
        Y(19) = 6.779354330227017D-04
        Y(20) = 4.437815269118510D-04
        Y(21) = 2.833264542938954D-04
        Y(22) = 1.765005798796805D-04
        Y(23) = 1.073342592697238D-04
        Y(24) = 6.374497601777217D-05
        Y(25) = 3.698645309704183D-05
        Y(26) = 2.097466832643746D-05
        Y(27) = 1.162956710412555D-05
        Y(28) = 6.306710405783322D-06
        Y(29) = 3.346286430868515D-06
        Y(30) = 1.737760074184334D-06
        Y(31) = 8.835366904275847D-07
        Y(32) = 4.399520411127637D-07
        Y(33) = 2.146181897152360D-07
        Y(34) = 1.025981211654928D-07
        Y(35) = 4.807864068784215D-08
        Y(36) = 2.209175152474847D-08
        Y(37) = 9.956251263138180D-09
        Y(38) = 4.402193653748924D-09
        Y(39) = 1.910149382204028D-09
        Y(40) = 8.135892921473050D-10
        Y(41) = 3.402477118549235D-10
        Y(42) = 1.397485617545782D-10
        Y(43) = 5.638575303049199D-11
        Y(44) = 2.235459707956947D-11
        Y(45) = 8.710498036398032D-12
        Y(46) = 3.336554275346643D-12
        Y(47) = 1.256679567784939D-12
        Y(48) = 4.654359053128788D-13
        Y(49) = 1.693559145599857D-13
        Y(50) = 5.996593816663054D-14
        Y(51) = 1.891330702629865D-14
        GO TO 310
!     25:
!        PROBLEM C5
150     Y(1) = -4.792730224323733D+00
        Y(2) = -2.420550725448973D+00
        Y(3) = -9.212509306014886D-01
        Y(4) = -4.217310404035213D+00
        Y(5) = 7.356202947498970D+00
        Y(6) = 3.223785985421212D+00
        Y(7) = 4.035559443262270D+00
        Y(8) = 1.719865528670555D+01
        Y(9) = 7.478910794233703D+00
        Y(10) = -2.998759326324844D+01
        Y(11) = -4.107310937550929D+00
        Y(12) = -9.277008321754407D-01
        Y(13) = -2.442125302518482D+01
        Y(14) = 2.381459045746554D+01
        Y(15) = 1.492096306951359D+01
        Y(16) = 3.499208963063806D-01
        Y(17) = -5.748487687912825D-01
        Y(18) = -2.551694020879149D-01
        Y(19) = -5.237040978903326D-01
        Y(20) = -2.493000463579661D-01
        Y(21) = -8.045341642044464D-02
        Y(22) = -3.875289237334110D-01
        Y(23) = 5.648603288767891D-02
        Y(24) = 3.023606472143342D-02
        Y(25) = 4.133856546712445D-02
        Y(26) = -2.862393029841379D-01
        Y(27) = -1.183032405136207D-01
        Y(28) = -1.511986457359206D-01
        Y(29) = -2.460068894318766D-01
        Y(30) = -3.189687411323877D-02
        GO TO 310

!     PROBLEM CLASS D

!     31:
!        PROBLEM D1
160     Y(1) = 2.198835352008397D-01
        Y(2) = 9.427076846341813D-01
        Y(3) = -9.787659841058176D-01
        Y(4) = 3.287977990962036D-01
        GO TO 310
!     32:
!        PROBLEM D2
170     Y(1) = -1.777027357140412D-01
        Y(2) = 9.467784719905892D-01
        Y(3) = -1.030294163192969D+00
        Y(4) = 1.211074890053952D-01
        GO TO 310
!     33:
!        PROBLEM D3
180     Y(1) = -5.780432953035361D-01
        Y(2) = 8.633840009194193D-01
        Y(3) = -9.595083730380727D-01
        Y(4) = -6.504915126712089D-02
        GO TO 310
!     34:
!        PROBLEM D4
190     Y(1) = -9.538990293416394D-01
        Y(2) = 6.907409024219432D-01
        Y(3) = -8.212674270877433D-01
        Y(4) = -1.539574259125825D-01
        GO TO 310
!     35:
!        PROBLEM D5
200     Y(1) = -1.295266250987574D+00
        Y(2) = 4.003938963792321D-01
        Y(3) = -6.775390924707566D-01
        Y(4) = -1.270838154278686D-01
        GO TO 310

!     PROBLEM CLASS E

!     41:
!        PROBLEM E1
210     Y(1) = 1.456723600728308D-01
        Y(2) = -9.883500195574063D-02
        GO TO 310
!     42:
!        PROBLEM E2
220     Y(1) = 2.008149762174948D+00
        Y(2) = -4.250887527320057D-02
        GO TO 310
!     43:
!        PROBLEM E3
230     Y(1) = -1.004178858647128D-01
        Y(2) = 2.411400132095954D-01
        GO TO 310
!     44:
!        PROBLEM E4
240     Y(1) = 3.395091444646555D+01
        Y(2) = 2.767822659672869D-01
        GO TO 310
!     45:
!        PROBLEM E5
250     Y(1) = 1.411797390542629D+01
        Y(2) = 2.400000000000002D+00
        GO TO 310

!     PROBLEM CLASS F

!     51:
!        PROBLEM F1
260     Y(1) = -1.294460621213470D1
        Y(2) = -2.208575158908672D-15
        GO TO 310
!     52:
!        PROBLEM F2
270     Y(1) = 70.03731057008607D0
        GO TO 310
!     53:
!        PROBLEM F3
280     Y(1) = -3.726957553088175D-1
        Y(2) = -6.230137949234190D-1
        GO TO 310
!     54:
!        PROBLEM F4
290     Y(1) = 9.815017249707434D-11
        GO TO 310
!     55:
!        PROBLEM F5
300     Y(1) = 1.0D0
310     CONTINUE
        IF (IWT<0) GO TO 320
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
        DOUBLE PRECISION Y(51), YP(51)
!     .. Local Scalars ..
        DOUBLE PRECISION C1, C2, D, EX, K2, M0, P, TEMP
        INTEGER I, I3, I3M2, ITEMP, J, L, LL, MM
!     .. Local Arrays ..
        DOUBLE PRECISION M(5), Q(5,5), R(5), YTEMP(51)
!     .. Data statements ..
!     THE FOLLOWING DATA IS FOR PROBLEM C5 AND DEFINES THE MASSES
!     OF THE 5 OUTER PLANETS ETC. IN SOLAR UNITS.
!     K2 IS THE GRAVITATIONAL CONSTANT.
!     THE NEXT DATA IS FOR PROBLEMS F1 AND F5.
!     C1 IS PI**2 + 0.1**2 AND C2 IS SUM I**(4/3) FOR I=1 TO 19.
        DATA M0/1.00000597682D0/, M/.954786104043D-3, .285583733151D-3, &
          .437273164546D-4, .517759138449D-4, .277777777778D-5/
        DATA K2/2.95912208286D0/
        DATA EX, C1, C2/.33333333333333333D0, 9.879604401089358D0, &
          438.4461015267790D0/

!     .. Executable Statements ..

        NFCN = NFCN + 1
        IF (IWT<0) GO TO 10
        DO I = 1, N
          YTEMP(I) = Y(I)
          Y(I) = Y(I)*W(I)
        END DO
10      CONTINUE
        GO TO (20,30,40,50,60,330,330,330,330,330,70,80,90,100,110,330,330, &
          330,330,330,120,130,140,150,160,330,330,330,330,330,190,190,190,190, &
          190,330,330,330,330,330,200,210,220,230,240,330,330,330,330,330,250, &
          270,290,300,320) ID
        GO TO 330

!     PROBLEM CLASS A

!     1:
!        PROBLEM A1
20      YP(1) = -Y(1)
        GO TO 330
!     2:
!        PROBLEM A2
30      YP(1) = -.5D0*Y(1)*Y(1)*Y(1)
        GO TO 330
!     3:
!        PROBLEM A3
40      YP(1) = Y(1)*COS(X)
        GO TO 330
!     4:
!        PROBLEM A4
50      YP(1) = (1.0D0-Y(1)/20.0D0)*Y(1)/4.D0
        GO TO 330
!     5:
!        PROBLEM A5
60      YP(1) = (Y(1)-X)/(Y(1)+X)
        GO TO 330

!     PROBLEM CLASS B

!     11:
!        PROBLEM B1
70      D = Y(1) - Y(1)*Y(2)
        YP(1) = D + D
        YP(2) = -(Y(2)-Y(1)*Y(2))
        GO TO 330
!     12:
!        PROBLEM B2
80      YP(1) = -Y(1) + Y(2)
        YP(3) = Y(2) - Y(3)
        YP(2) = -YP(1) - YP(3)
        GO TO 330
!     13:
!        PROBLEM B3
90      D = Y(2)*Y(2)
        YP(1) = -Y(1)
        YP(2) = Y(1) - D
        YP(3) = D
        GO TO 330
!     14:
!        PROBLEM B4
100     D = SQRT(Y(1)*Y(1)+Y(2)*Y(2))
        YP(1) = -Y(2) - Y(1)*Y(3)/D
        YP(2) = Y(1) - Y(2)*Y(3)/D
        YP(3) = Y(1)/D
        GO TO 330
!     15:
!        PROBLEM B5
110     YP(1) = Y(2)*Y(3)
        YP(2) = -Y(1)*Y(3)
        YP(3) = -.51D0*Y(1)*Y(2)
        GO TO 330

!     PROBLEM CLASS C

!     21:
!        PROBLEM C1
120     YP(1) = -Y(1)
        DO I = 2, 9
          YP(I) = Y(I-1) - Y(I)
        END DO
        YP(10) = Y(9)
        GO TO 330
!     22:
!        PROBLEM C2
130     YP(1) = -Y(1)
        DO I = 2, 9
          YP(I) = REAL(I-1)*Y(I-1) - REAL(I)*Y(I)
        END DO
        YP(10) = 9.0D0*Y(9)
        GO TO 330
!     23:
!        PROBLEM C3
140     YP(1) = -2.0D0*Y(1) + Y(2)
        DO I = 2, 9
          YP(I) = Y(I-1) - 2.0D0*Y(I) + Y(I+1)
        END DO
        YP(10) = Y(9) - 2.0D0*Y(10)
        GO TO 330
!     24:
!        PROBLEM C4
150     YP(1) = -2.0D0*Y(1) + Y(2)
        DO I = 2, 50
          YP(I) = Y(I-1) - 2.0D0*Y(I) + Y(I+1)
        END DO
        YP(51) = Y(50) - 2.0D0*Y(51)
        GO TO 330
!     25:
!        PROBLEM C5
160     I = 0
        DO L = 3, 15, 3
          I = I + 1
          P = Y(L-2)**2 + Y(L-1)**2 + Y(L)**2
          R(I) = 1.0D0/(P*SQRT(P))
          J = 0
          DO LL = 3, 15, 3
            J = J + 1
            IF (LL/=L) GO TO 170
            GO TO 180
!              THEN
170         P = (Y(L-2)-Y(LL-2))**2 + (Y(L-1)-Y(LL-1))**2 + (Y(L)-Y(LL))**2
            Q(I,J) = 1.0D0/(P*SQRT(P))
            Q(J,I) = Q(I,J)
180         CONTINUE
          END DO
        END DO
        I3 = 0
        DO I = 1, 5
          I3 = I3 + 3
          I3M2 = I3 - 2
          DO LL = I3M2, I3
            MM = LL - I3
            YP(LL) = Y(LL+15)
            P = 0.0D0
            DO J = 1, 5
              MM = MM + 3
              IF (J/=I) P = P + M(J)*(Y(MM)*(Q(I,J)-R(J))-Y(LL)*Q(I,J))
            END DO
            YP(LL+15) = K2*(-(M0+M(I))*Y(LL)*R(I)+P)
          END DO
        END DO
        GO TO 330

!     PROBLEM CLASS D

!     31:32:33:34:35:
!        PROBLEMS D1, D2, D3, D4, D5
190     YP(1) = Y(3)
        YP(2) = Y(4)
        D = Y(1)*Y(1) + Y(2)*Y(2)
        D = SQRT(D*D*D)
        YP(3) = -Y(1)/D
        YP(4) = -Y(2)/D
        GO TO 330

!     PROBLEM CLASS E

!     41:
!        PROBLEM E1
200     YP(1) = Y(2)
        YP(2) = -(Y(2)/(X+1.0D0)+(1.0D0-.25D0/(X+1.0D0)**2)*Y(1))
        GO TO 330
!     42:
!        PROBLEM E2
210     YP(1) = Y(2)
        YP(2) = (1.0D0-Y(1)*Y(1))*Y(2) - Y(1)
        GO TO 330
!     43:
!        PROBLEM E3
220     YP(1) = Y(2)
        YP(2) = Y(1)**3/6.0D0 - Y(1) + 2.0D0*SIN(2.78535D0*X)
        GO TO 330
!     44:
!        PROBLEM E4
230     YP(1) = Y(2)
        YP(2) = .032D0 - .4D0*Y(2)*Y(2)
        GO TO 330
!     45:
!        PROBLEM E5
240     YP(1) = Y(2)
        YP(2) = SQRT(1.0D0+Y(2)*Y(2))/(25.0D0-X)
        GO TO 330

!     PROBLEM CLASS F

!     51:
!        PROBLEM F1
250     YP(1) = Y(2)
        YP(2) = .2D0*Y(2) - C1*Y(1)
        ITEMP = INT(X)
        IF ((ITEMP/2)*2==ITEMP) GO TO 260
        YP(2) = YP(2) - 1.0D0
        GO TO 330
260     YP(2) = YP(2) + 1.0D0
        GO TO 330
!     52:
!        PROBLEM F2
270     ITEMP = INT(X)
        IF ((ITEMP/2)*2==ITEMP) GO TO 280
        YP(1) = 55.0D0 - .50D0*Y(1)
        GO TO 330
280     YP(1) = 55.0D0 - 1.50D0*Y(1)
        GO TO 330
!     53:
!        PROBLEM F3
290     YP(1) = Y(2)
        YP(2) = .01D0*Y(2)*(1.0D0-Y(1)**2) - Y(1) - ABS(SIN( &
          3.1415926535897932D0*X))
        GO TO 330
!     54:
!        PROBLEM F4
300     IF (X>10.) GO TO 310
        TEMP = X - 5.0D0
        YP(1) = -2.0D0/21.0D0 - 120.0D0*TEMP/(1.0D0+4.0D0*TEMP**2)**16
        GO TO 330
310     YP(1) = -2.0D0*Y(1)
        GO TO 330
!     55:
!        PROBLEM F5
320     YP(1) = Y(1)*(4.0D0/(3.0D0*C2))*(SIGN(ABS(X- &
          1.0D0)**EX,X-1.0D0)+SIGN(ABS(X-2.0D0)**EX,X-2.0D0)+SIGN(ABS(X- &
          3.0D0)**EX,X-3.0D0)+SIGN(ABS(X-4.0D0)**EX,X-4.0D0)+SIGN(ABS(X- &
          5.0D0)**EX,X-5.0D0)+SIGN(ABS(X-6.0D0)**EX,X-6.0D0)+SIGN(ABS(X- &
          7.0D0)**EX,X-7.0D0)+SIGN(ABS(X-8.0D0)**EX,X-8.0D0)+SIGN(ABS(X- &
          9.0D0)**EX,X-9.0D0)+SIGN(ABS(X-10.0D0)**EX,X-10.0D0)+SIGN(ABS(X- &
          11.0D0)**EX,X-11.0D0)+SIGN(ABS(X-12.0D0)**EX,X-12.0D0)+SIGN(ABS(X- &
          13.0D0)**EX,X-13.0D0)+SIGN(ABS(X-14.0D0)**EX,X-14.0D0)+SIGN(ABS(X- &
          15.0D0)**EX,X-15.0D0)+SIGN(ABS(X-16.0D0)**EX,X-16.0D0)+SIGN(ABS(X- &
          17.0D0)**EX,X-17.0D0)+SIGN(ABS(X-18.0D0)**EX,X-18.0D0)+SIGN(ABS(X- &
          19.0D0)**EX,X-19.0D0))
330     CONTINUE
        IF (IWT<0) GO TO 340
        DO I = 1, N
          YP(I) = YP(I)/W(I)
          Y(I) = YTEMP(I)
        END DO
340     CONTINUE
        RETURN
      END SUBROUTINE FCN

    END MODULE NONSTIFFSET

!******************************************************************

    PROGRAM DEMONONSTIFF

      USE NONSTIFFSET
      USE DVODE_F90_M

      IMPLICIT NONE
      INTEGER ITASK, ISTATE, ISTATS, NEQ, I, CLASS, PROBLEM, MYID, ITEST, ITOL
      DOUBLE PRECISION RSTATS, T, TOUT, HBEGIN, HBOUND, TBEGIN, TEND, Y, EPS, &
        YINIT, YFINAL, RELERR_TOLERANCES, ABSERR_TOLERANCES, AERROR
      LOGICAL USEW, USEHBEGIN, ERRORS
      DIMENSION Y(51), RSTATS(22), ISTATS(31), YINIT(51), YFINAL(51), MYID(55)
      DIMENSION RELERR_TOLERANCES(51), ABSERR_TOLERANCES(51), AERROR(51)
      TYPE (VODE_OPTS) :: OPTIONS
      DATA MYID/1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 11, 12, 13, 14, 15, 0, 0, 0, 0, &
        0, 21, 22, 23, 24, 25, 0, 0, 0, 0, 0, 31, 32, 33, 34, 35, 0, 0, 0, 0, &
        0, 41, 42, 43, 44, 45, 0, 0, 0, 0, 0, 51, 52, 53, 54, 55/

      OPEN (UNIT=6,FILE='nonstiffoptions.dat')

      ERRORS = .FALSE.

      DO ITOL = 4, 10, 2
      WRITE(6,*) 'Results for ITOL = ', ITOL, ' follow.'
      DO ITEST = 1, 55
        ID = MYID(ITEST)
        IF (ID==0) GO TO 20
        WRITE (6,90010)
        CLASS = ID/10
        PROBLEM = ID - 10*CLASS
        IF (CLASS==0) THEN
          WRITE (6,90000) PROBLEM
        ELSE IF (CLASS==1) THEN
          WRITE (6,90001) PROBLEM
        ELSE IF (CLASS==2) THEN
          WRITE (6,90002) PROBLEM
        ELSE IF (CLASS==3) THEN
          WRITE (6,90003) PROBLEM
        ELSE IF (CLASS==4) THEN
          WRITE (6,90004) PROBLEM
        ELSE IF (CLASS==5) THEN
          WRITE (6,90005) PROBLEM
        END IF
!     Scale the odes?
        USEW = .TRUE.
!     Use the IVALU starting step size?
        USEHBEGIN = .TRUE.
        IWT = -1
        IF (USEW) IWT = 1
        CALL IVALU(TBEGIN,TEND,HBEGIN,HBOUND,YINIT)
        IF (.NOT. USEHBEGIN) HBEGIN = 0.0D0
        NEQ = N
        T = TBEGIN
        TOUT = TEND
        Y(1:NEQ) = YINIT(1:NEQ)
        EPS = 1.0D0 / 10.0D0**ITOL
        RELERR_TOLERANCES(1:NEQ) = EPS
        ABSERR_TOLERANCES(1:NEQ) = EPS
        WRITE (6,90007) ID, TBEGIN, TEND, HBEGIN, HBOUND, IWT, N, EPS, Y(1:NEQ)
        ITASK = 1
        ISTATE = 1
        OPTIONS = SET_OPTS(RELERR_VECTOR=RELERR_TOLERANCES(1:NEQ), &
          ABSERR_VECTOR=ABSERR_TOLERANCES(1:NEQ),MXSTEP=100000,H0=HBEGIN, &
          HMAX=HBOUND)
10      CONTINUE
        CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
        IF (ISTATE<0) THEN
          WRITE (6,90006) ISTATE
          ERRORS = .TRUE.
!         STOP
        END IF
        CALL GET_STATS(RSTATS,ISTATS)
        WRITE (6,90009) ISTATS(11), ISTATS(12)
        IF (TOUT<TEND) GO TO 10
        CALL EVALU(YFINAL)
        DO I = 1, NEQ
          AERROR(I) = ABS(Y(I)-YFINAL(I))
        END DO
        WRITE (6,90008) (I,Y(I),YFINAL(I),AERROR(I),I=1,NEQ)
20    END DO ! End of ITEST Loop
      END DO ! End of ITOL Loop

      IF (ERRORS) THEN
         WRITE(6,90011)
      ELSE
         WRITE(6,90012)
      END IF

!     Format statements for this problem:
90000 FORMAT (' Class/Problem = A',I1)
90001 FORMAT (' Class/Problem = B',I1)
90002 FORMAT (' Class/Problem = C',I1)
90003 FORMAT (' Class/Problem = D',I1)
90004 FORMAT (' Class/Problem = E',I1)
90005 FORMAT (' Class/Problem = F',I1)
90006 FORMAT (' An error occurred in VODE_F90. ISTATE = ',I3)
90007 FORMAT (' Problem ID       = ',I3,/,' Initial time     = ',D15.5,/, &
        ' Final time       = ',D15.5,/,' Initial stepsize = ',D15.5,/, &
        ' Maximum stepsize = ',D15.5,/,' IWT flag         = ',I3,/, &
        ' Number of odes   = ',I3,/,' Error tolerance  = ',D15.5,/, &
        ' Initial solution = ',/,(D15.5))
90008 FORMAT (' Computed and reference solutions and absolute', &
        ' errors follow:',/,(I3,3D15.5))
90009 FORMAT (' Steps = ',I10,' f-s = ',I10)
90010 FORMAT (' _________________________________________')
90011 FORMAT (' Errors occurred.')
90012 FORMAT (' No errors occurred.')
      STOP

    END PROGRAM DEMONONSTIFF
