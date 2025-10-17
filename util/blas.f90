!--------------------------------------------------------------
!
! BLAS/LAPACK-like subroutines used by the integration algorithms
! It is recommended to replace them by calls to the optimized
!      BLAS/LAPACK library for your machine
!
!  (C) Adrian Sandu, Aug. 2004
!      Virginia Polytechnic Institute and State University
!--------------------------------------------------------------

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: The following BLAS functions have been removed, %%%
!%%% as they now have been replaced by pure F90 code       %%%
!%%% in the various integrator modules;                    %%%
!%%%                                                       %%%
!%%% (1) WCOPY                                             %%%
!%%% (2) WAXPY                                             %%%
!%%% (3) WSCAL                                             %%%
!%%% (4) WLAMCH                                            %%%
!%%% (5) WLAMCH_ADD                                        %%%
!%%% (6) SET2ZERO                                          %%%
!%%% (7) WADD                                              %%%
!%%%                                                       %%%
!%%% @yantosca, 17 Oct 2025                                %%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!--------------------------------------------------------------
      KPP_REAL FUNCTION WDOT (N, DX, incX, DY, incY) 
!--------------------------------------------------------------
!     dot produce: wdot = x(1:N)*y(1:N) 
!     only for incX=incY=1
!     after BLAS
!     replace this by the function from the optimized BLAS implementation:
!         CALL SDOT(N,X,1,Y,1) or  CALL DDOT(N,X,1,Y,1)
!--------------------------------------------------------------
!      USE messy_mecca_kpp_Precision
!--------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: N, incX, incY
      KPP_REAL :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 
!                                                                       
!     Code for unequal or nonpositive increments.                       
!                                                                       
    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        WDOT = WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 
!                                                                       
!     Code for both increments equal to 1.                              
!                                                                       
!     Clean-up loop so remaining vector length is a multiple of 5.      
!                                                                       
   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         WDOT = WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          WDOT = WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) + DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 
!                                                                       
!     Code for equal, positive, non-unit increments.                    
!                                                                       
   60 NS = N*incX 
      DO i = 1,NS,incX 
        WDOT = WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION WDOT                                          

!--------------------------------------------------------------
      SUBROUTINE WGEFA(N,A,Ipvt,info)
!--------------------------------------------------------------
!     WGEFA FACTORS THE MATRIX A (N,N) BY
!           GAUSS ELIMINATION WITH PARTIAL PIVOTING
!     LINPACK - LIKE 
!--------------------------------------------------------------
!
      INTEGER       :: N,Ipvt(N),info
      KPP_REAL :: A(N,N)
      KPP_REAL :: t, dmax, da
      INTEGER       :: j,k,l
      KPP_REAL, PARAMETER :: ZERO = 0.0, ONE = 1.0

      info = 0

size: IF (n > 1) THEN
      
col:  DO k = 1, n-1

!        find l = pivot index
!        l = idamax(n-k+1,A(k,k),1) + k - 1
         l = k; dmax = abs(A(k,k))
         DO j = k+1,n
            da = ABS(A(j,k))
            IF (da > dmax) THEN
              l = j; dmax = da
            END IF
         END DO
         Ipvt(k) = l

!        zero pivot implies this column already triangularized
         IF (ABS(A(l,k)) < TINY(ZERO)) THEN
            info = k
            return
         ELSE   
            IF (l /= k) THEN
               t = A(l,k); A(l,k) = A(k,k); A(k,k) = t
            END IF
            t = -ONE/A(k,k)
            CALL WSCAL(n-k,t,A(k+1,k),1)
            DO j = k+1, n
               t = A(l,j)
               IF (l /= k) THEN
                  A(l,j) = A(k,j); A(k,j) = t
               END IF
               CALL WAXPY(n-k,t,A(k+1,k),1,A(k+1,j),1)
            END DO         
         END IF
         
       END DO col
       
      END IF size
      
      Ipvt(N) = N
      IF (ABS(A(N,N)) == ZERO) info = N
      
      END SUBROUTINE WGEFA


!--------------------------------------------------------------
      SUBROUTINE WGESL(Trans,N,A,Ipvt,b)
!--------------------------------------------------------------
!     WGESL solves the system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by WGEFA.
!
!     Trans      = 'N'   to solve  A*x = b ,
!                = 'T'   to solve  transpose(A)*x = b
!     LINPACK - LIKE 
!--------------------------------------------------------------

      INTEGER       :: N,Ipvt(N)
      CHARACTER     :: trans
      KPP_REAL :: A(N,N),b(N)
      KPP_REAL :: t
      INTEGER       :: k,kb,l

      
      SELECT CASE (Trans)

      CASE ('n','N')  !  Solve  A * x = b

!        first solve  L*y = b
         IF (n >= 2) THEN
          DO k = 1, n-1
            l = Ipvt(k)
            t = b(l)
            IF (l /= k) THEN
               b(l) = b(k)
               b(k) = t
            END IF
            CALL WAXPY(n-k,t,a(k+1,k),1,b(k+1),1)
          END DO
         END IF
!        now solve  U*x = y
         DO kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            CALL WAXPY(k-1,t,a(1,k),1,b(1),1)
         END DO
      
      CASE ('t','T')  !  Solve transpose(A) * x = b

!        first solve  trans(U)*y = b
         DO k = 1, n
            t = WDOT(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
         END DO
!        now solve trans(L)*x = y
         IF (n >= 2) THEN
         DO kb = 1, n-1
            k = n - kb
            b(k) = b(k) + WDOT(n-k,a(k+1,k),1,b(k+1),1)
            l = Ipvt(k)
            IF (l /= k) THEN
               t = b(l); b(l) = b(k); b(k) = t
            END IF
         END DO
         END IF
   
      END SELECT

      END SUBROUTINE WGESL
