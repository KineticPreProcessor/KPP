
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppDecomp( JVS, IER )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Sparse LU factorization
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER  :: IER
      KPP_REAL :: JVS(LU_NONZERO), W(NVAR), a
      INTEGER  :: k, kk, j, jj

      a = 0. ! mz_rs_20050606
      IER = 0
      DO k=1,NVAR
        ! mz_rs_20050606: don't check if real value == 0
        ! IF ( JVS( LU_DIAG(k) ) .EQ. 0. ) THEN
        IF ( ABS(JVS(LU_DIAG(k))) < TINY(a) ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
      
END SUBROUTINE KppDecomp


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppDecompCmplx( JVS, IER )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Sparse LU factorization, complex
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER        :: IER
      DOUBLE COMPLEX :: JVS(LU_NONZERO), W(NVAR), a
      KPP_REAL  :: b = 0.0
      INTEGER        :: k, kk, j, jj

      IER = 0
      DO k=1,NVAR
        IF ( ABS(JVS(LU_DIAG(k))) < TINY(b) ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
      
END SUBROUTINE KppDecompCmplx


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppDecompCmplxR( JVSR, JVSI, IER )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Sparse LU factorization, complex
!   (Real and Imaginary parts are used instead of complex data type)     
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER       :: IER
      KPP_REAL :: JVSR(LU_NONZERO), JVSI(LU_NONZERO) 
      KPP_REAL :: WR(NVAR), WI(NVAR), ar, ai, den
      INTEGER       :: k, kk, j, jj

      IER = 0
      ar  = 0.0
      DO k=1,NVAR
        IF (  ( ABS(JVSR(LU_DIAG(k))) < TINY(ar) ) .AND. &
              ( ABS(JVSI(LU_DIAG(k))) < TINY(ar) ) )  THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              WR( LU_ICOL(kk) ) = JVSR(kk)
              WI( LU_ICOL(kk) ) = JVSI(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            den = JVSR(LU_DIAG(j))**2 + JVSI(LU_DIAG(j))**2
            ar = -(WR(j)*JVSR(LU_DIAG(j)) + WI(j)*JVSI(LU_DIAG(j)))/den
            ai = -(WI(j)*JVSR(LU_DIAG(j)) - WR(j)*JVSI(LU_DIAG(j)))/den
            WR(j) = -ar
            WI(j) = -ai
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               WR( LU_ICOL(jj) ) = WR( LU_ICOL(jj) ) + ar*JVSR(jj) - ai*JVSI(jj)
               WI( LU_ICOL(jj) ) = WI( LU_ICOL(jj) ) + ar*JVSI(jj) + ai*JVSR(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVSR(kk) = WR( LU_ICOL(kk) )
            JVSI(kk) = WI( LU_ICOL(kk) )
         END DO
      END DO

END SUBROUTINE KppDecompCmplxR


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveIndirect( JVS, X )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Sparse solve subroutine using indirect addressing
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER  :: i, j
      KPP_REAL :: JVS(LU_NONZERO), X(NVAR), sum

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO  
      END DO

      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
      
END SUBROUTINE KppSolveIndirect


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveTRIndirect( JVS, X )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Complex sparse solve transpose subroutine using indirect addressing
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER       :: i, j
      KPP_REAL :: JVS(LU_NONZERO), X(NVAR)

      DO i=1,NVAR
        X(i) = X(i)/JVS(LU_DIAG(i))
	! subtract all nonzero elements in row i of JVS from X
        DO j=LU_DIAG(i)+1,LU_CROW(i+1)-1
	  X(LU_ICOL(j)) = X(LU_ICOL(j))-JVS(j)*X(i)
	END DO
      END DO

      DO i=NVAR, 1, -1
	! subtract all nonzero elements in row i of JVS from X
        DO j=LU_CROW(i),LU_DIAG(i)-1
	  X(LU_ICOL(j)) = X(LU_ICOL(j))-JVS(j)*X(i)
	END DO
      END DO
      
END SUBROUTINE KppSolveTRIndirect


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveCmplx( JVS, X )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Complex sparse solve subroutine using indirect addressing
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER        :: i, j
      DOUBLE COMPLEX :: JVS(LU_NONZERO), X(NVAR), sum

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO  
      END DO

      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
      
END SUBROUTINE KppSolveCmplx

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveCmplxR( JVSR, JVSI, XR, XI )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Complex sparse solve subroutine using indirect addressing
!   (Real and Imaginary parts are used instead of complex data type)     
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER       ::  i, j
      KPP_REAL ::  JVSR(LU_NONZERO), JVSI(LU_NONZERO), XR(NVAR), XI(NVAR), sumr, sumi, den

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             XR(i) = XR(i) - (JVSR(j)*XR(LU_ICOL(j)) - JVSI(j)*XI(LU_ICOL(j)))
             XI(i) = XI(i) - (JVSR(j)*XI(LU_ICOL(j)) + JVSI(j)*XR(LU_ICOL(j)))
         END DO  
      END DO

      DO i=NVAR,1,-1
        sumr = XR(i); sumi = XI(i)
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
            sumr = sumr - (JVSR(j)*XR(LU_ICOL(j)) - JVSI(j)*XI(LU_ICOL(j)))
            sumi = sumi - (JVSR(j)*XI(LU_ICOL(j)) + JVSI(j)*XR(LU_ICOL(j)))
        END DO
        den   = JVSR(LU_DIAG(i))**2 + JVSI(LU_DIAG(i))**2
        XR(i) = (sumr*JVSR(LU_DIAG(i)) + sumi*JVSI(LU_DIAG(i)))/den
        XI(i) = (sumi*JVSR(LU_DIAG(i)) - sumr*JVSI(LU_DIAG(i)))/den
      END DO
      
END SUBROUTINE KppSolveCmplxR


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveTRCmplx( JVS, X )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Complex sparse solve transpose subroutine using indirect addressing
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER        :: i, j
      DOUBLE COMPLEX :: JVS(LU_NONZERO), X(NVAR)

      DO i=1,NVAR
        X(i) = X(i)/JVS(LU_DIAG(i))
	! subtract all nonzero elements in row i of JVS from X
        DO j=LU_DIAG(i)+1,LU_CROW(i+1)-1
	  X(LU_ICOL(j)) = X(LU_ICOL(j))-JVS(j)*X(i)
	END DO
      END DO

      DO i=NVAR, 1, -1
	! subtract all nonzero elements in row i of JVS from X
        DO j=LU_CROW(i),LU_DIAG(i)-1
	  X(LU_ICOL(j)) = X(LU_ICOL(j))-JVS(j)*X(i)
	END DO
      END DO
      
END SUBROUTINE KppSolveTRCmplx


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE KppSolveTRCmplxR( JVSR, JVSI, XR, XI )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Complex sparse solve transpose subroutine using indirect addressing
!   (Real and Imaginary parts are used instead of complex data type)     
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE KPP_ROOT_Parameters
  USE KPP_ROOT_JacobianSP

      INTEGER       ::  i, j
      KPP_REAL ::  JVSR(LU_NONZERO), JVSI(LU_NONZERO), XR(NVAR), XI(NVAR), den

      DO i=1,NVAR
        den   = JVSR(LU_DIAG(i))**2 + JVSI(LU_DIAG(i))**2
        XR(i) = (XR(i)*JVSR(LU_DIAG(i)) + XI(i)*JVSI(LU_DIAG(i)))/den
        XI(i) = (XI(i)*JVSR(LU_DIAG(i)) - XR(i)*JVSI(LU_DIAG(i)))/den
	! subtract all nonzero elements in row i of JVS from X
        DO j=LU_DIAG(i)+1,LU_CROW(i+1)-1
	  XR(LU_ICOL(j)) = XR(LU_ICOL(j))-(JVSR(j)*XR(i) - JVSI(j)*XI(i))
	  XI(LU_ICOL(j)) = XI(LU_ICOL(j))-(JVSI(j)*XR(i) + JVSR(j)*XI(i))
	END DO
      END DO

      DO i=NVAR, 1, -1
	! subtract all nonzero elements in row i of JVS from X
        DO j=LU_CROW(i),LU_DIAG(i)-1
	  XR(LU_ICOL(j)) = XR(LU_ICOL(j))-(JVSR(j)*XR(i) - JVSI(j)*XI(i))
	  XI(LU_ICOL(j)) = XI(LU_ICOL(j))-(JVSI(j)*XR(i) + JVSR(j)*XI(i))
	END DO
      END DO
      
END SUBROUTINE KppSolveTRCmplxR


!
! Next few commented subroutines perform sparse big linear algebra
!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!SUBROUTINE KppDecompBig( JVS, IP, IER )
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!        Sparse LU factorization
!!        for the Runge Kutta (3n)x(3n) linear system
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  USE KPP_ROOT_Parameters
!  USE KPP_ROOT_JacobianSP
!
!      INTEGER  :: IP3(3), IER, IP(3,NVAR)
!      KPP_REAL :: JVS(3,3,LU_NONZERO), W(3,3,NVAR), a(3,3), E(3,3)
!      INTEGER  :: k, kk, j, jj
!
!      a = 0.0d0
!      IER = 0
!      DO k=1,NVAR
!        DO kk = LU_CROW(k), LU_CROW(k+1)-1
!              W( 1:3,1:3,LU_ICOL(kk) ) = JVS(1:3,1:3,kk)
!        END DO
!        DO kk = LU_CROW(k), LU_DIAG(k)-1
!            j = LU_ICOL(kk)
!            E(1:3,1:3) = JVS( 1:3,1:3,LU_DIAG(j) )
!            ! CALL DGETRF(3,3,E,3,IP3,IER) 
!            CALL FAC3(E,IP3,IER)
!            IF ( IER /= 0 )  RETURN
!            ! a = W(j) / JVS( LU_DIAG(j) )
!            a(1:3,1:3) = W( 1:3,1:3,j )
!            ! CALL DGETRS ('N',3,3,E,3,IP3,a,3,IER) 
!            CALL SOL3('N',E,IP3,a(1,1))
!            CALL SOL3('N',E,IP3,a(1,2))
!            CALL SOL3('N',E,IP3,a(1,3))
!            W(1:3,1:3,j) = a(1:3,1:3)
!            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
!               W( 1:3,1:3,LU_ICOL(jj) ) = W( 1:3,1:3,LU_ICOL(jj) ) &
!                        - MATMUL( a(1:3,1:3) , JVS(1:3,1:3,jj) )
!            END DO
!         END DO
!         DO kk = LU_CROW(k), LU_CROW(k+1)-1
!            JVS(1:3,1:3,kk) = W( 1:3,1:3,LU_ICOL(kk) )
!         END DO
!      END DO
!
!      DO k=1,NVAR
!         ! CALL WGEFA(JVS(1,1,LU_DIAG(k)),3,3,IP(1,k),IER)
!         ! CALL DGETRF(3,3,JVS(1,1,LU_DIAG(k)),3,IP(1,k),IER)
!         CALL FAC3(JVS(1,1,LU_DIAG(k)),IP(1,k),IER)
!         IF ( IER /= 0 )  RETURN
!      END DO 
!      
!END SUBROUTINE KppDecompBig
!
!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!SUBROUTINE KppSolveBig( JVS, IP, X )
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!        Sparse solve subroutine using indirect addressing
!!        for the Runge Kutta (3n)x(3n) linear system
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  USE KPP_ROOT_Parameters
!  USE KPP_ROOT_JacobianSP
!
!      INTEGER  :: i, j, k, m, IP3(3), IP(3,NVAR), IER
!      KPP_REAL :: JVS(3,3,LU_NONZERO), X(3,NVAR), sum(3)
!
!      DO i=1,NVAR
!        DO j = LU_CROW(i), LU_DIAG(i)-1 
!          !X(1:3,i) = X(1:3,i) - MATMUL(JVS(1:3,1:3,j),X(1:3,LU_ICOL(j)));
!          DO k=1,3
!            DO m=1,3
!	       X(k,i) = X(k,i) - JVS(k,m,j)*X(m,LU_ICOL(j))
!            END DO
!          END DO
!        END DO  
!      END DO
!
!      DO i=NVAR,1,-1
!        sum(1:3) = X(1:3,i);
!        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
!          !sum(1:3) = sum(1:3) - MATMUL(JVS(1:3,1:3,j),X(1:3,LU_ICOL(j)));
!          DO k=1,3
!            DO m=1,3
!	       sum(k) = sum(k) - JVS(k,m,j)*X(m,LU_ICOL(j))
!            END DO
!          END DO
!        END DO
!        ! X(i) = sum/JVS(LU_DIAG(i));
!        ! CALL DGETRS ('N',3,1,JVS(1:3,1:3,LU_DIAG(i)),3,IP(1,i),sum,3,0) 
!        ! CALL WGESL('N',JVS(1,1,LU_DIAG(i)),3,3,IP(1,i),sum)
!        CALL SOL3('N',JVS(1,1,LU_DIAG(i)),IP(1,i),sum)
!        X(1:3,i) = sum(1:3)
!      END DO
!      
!END SUBROUTINE KppSolveBig
!
!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!SUBROUTINE KppSolveBigTR( JVS, IP, X )
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!        Big sparse transpose solve using indirect addressing
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  USE KPP_ROOT_Parameters
!  USE KPP_ROOT_JacobianSP
!
!      INTEGER       :: i, j, k, m, IP(3,NVAR)
!      KPP_REAL :: JVS(3,3,LU_NONZERO), X(3,NVAR)
!
!      DO i=1,NVAR
!        ! X(i) = X(i)/JVS(LU_DIAG(i))
!        CALL SOL3('T',JVS(1,1,LU_DIAG(i)),IP(1,i),X(1,i))
!        DO j=LU_DIAG(i)+1,LU_CROW(i+1)-1
!	  !X(1:3,LU_ICOL(j)) = X(1:3,LU_ICOL(j)) &
!          !    - MATMUL( TRANSPOSE(JVS(1:3,1:3,j)), X(1:3,i) )
!          DO k=1,3
!            DO m=1,3
!	       X(k,LU_ICOL(j)) = X(k,LU_ICOL(j)) - JVS(m,k,j)*X(m,i)
!            END DO
!          END DO
!	END DO
!      END DO
!
!      DO i=NVAR, 1, -1
!        DO j=LU_CROW(i),LU_DIAG(i)-1
!	  !X(1:3,LU_ICOL(j)) = X(1:3,LU_ICOL(j)) &
!          !   - MATMUL( TRANSPOSE(JVS(1:3,1:3,j)), X(1:3,i) )
!          DO k=1,3
!            DO m=1,3
!	       X(k,LU_ICOL(j)) = X(k,LU_ICOL(j)) - JVS(m,k,j)*X(m,i)
!            END DO
!          END DO
!	END DO
!      END DO
!      
!END SUBROUTINE KppSolveBigTR
!
!
!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!SUBROUTINE FAC3(A,IPVT,INFO)
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!     FAC3 FACTORS THE MATRIX A (3,3) BY
!!           GAUSS ELIMINATION WITH PARTIAL PIVOTING
!!     LINPACK - LIKE 
!!
!!     Remove comments to perform pivoting
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!
!      KPP_REAL :: A(3,3)
!      INTEGER       :: IPVT(3),INFO
!!      INTEGER       :: L
!!      KPP_REAL :: t, dmax, da, TMP(3)
!      KPP_REAL, PARAMETER :: ZERO = 0.0, ONE = 1.0
!
!      info = 0
!!      t = TINY(da)
!!      
!!      da = ABS(A(1,1)); L = 1
!!      IF ( ABS(A(2,1))>da ) THEN
!!        da = ABS(A(2,1)); L = 2
!!        IF ( ABS(A(3,1))>da ) THEN
!!          L = 3
!!        END IF  
!!      END IF  
!!      IPVT(1)  = L
!!      IF (L /=1 ) THEN
!!         TMP(1:3) = A(L,1:3)
!!         A(L,1:3) = A(1,1:3)
!!         A(1,1:3) = TMP(1:3)
!!      END IF
!!      IF (ABS(A(1,1)) < t) THEN
!!         info = 1
!!         return
!!      END IF   
!!
!      A(2,1) = A(2,1)/A(1,1)
!      A(2,2) = A(2,2) - A(2,1)*A(1,2)
!      A(2,3) = A(2,3) - A(2,1)*A(1,3)
!      A(3,1) = A(3,1)/A(1,1)
!      A(3,2) = A(3,2) - A(3,1)*A(1,2)
!      A(3,3) = A(3,3) - A(3,1)*A(1,3)
!      
!!      IPVT(2)  = 2
!!      IF (ABS(A(3,2))>ABS(A(2,2))) THEN
!!         IPVT(2)  = 3
!!         TMP(2:3) = A(3,2:3)
!!         A(3,2:3) = A(2,2:3)
!!         A(2,2:3) = TMP(2:3)
!!      END IF
!!      IF (ABS(A(2,2)) < t) THEN
!!         info = 1
!!         return
!!      END IF   
!!      
!      A(3,2)   = A(3,2)/A(2,2)
!      A(3,3)   = A(3,3) - A(3,2)*A(2,3)
!      IPVT(3)  = 3
!      
!END SUBROUTINE FAC3
!
!
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!SUBROUTINE SOL3(Trans,A,IPVT,b)
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!     SOL3 solves the system 3x3
!!     A * x = b  or  trans(a) * x = b
!!     using the factors computed by WGEFA.
!!
!!     Trans      = 'N'   to solve  A*x = b ,
!!                = 'T'   to solve  transpose(A)*x = b
!!     LINPACK - LIKE 
!!
!!     Remove comments to use pivoting
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!      CHARACTER     :: Trans
!      KPP_REAL :: a(3,3),b(3)
!      INTEGER       :: IPVT(3)
!!      INTEGER       :: L
!!      KPP_REAL :: TMP
!      
!      SELECT CASE (Trans)
!
!      CASE ('n','N')  !  Solve  A * x = b
!
!!     Solve  L*y = b
!!         L = IPVT(1)
!!         IF (L /= 1) THEN
!!            TMP = B(1); B(1) = B(L); B(L) = TMP
!!         END IF
!         b(2) = b(2)-A(2,1)*b(1)
!         b(3) = b(3)-A(3,1)*b(1)
!         
!!         L = IPVT(2)
!!         IF (L /= 2) THEN
!!            TMP = B(2); B(2) = B(L); B(L) = TMP
!!         END IF
!         b(3) = b(3)-A(3,2)*b(2)
!
!!     Solve  U*x = y
!         b(3) = b(3)/A(3,3)
!         b(2) = (b(2)-A(2,3)*b(3))/A(2,2)
!         b(1) = (b(1)-A(1,3)*b(3)-A(1,2)*b(2))/A(1,1)
!      
!      
!      CASE ('t','T')  !  Solve transpose(A) * x = b
!
!!      Solve transpose(U)*y = b
!         b(1) = b(1)/A(1,1)
!         b(2) = (b(2)-A(1,2)*b(1))/A(2,2)
!         b(3) = (b(3)-A(1,3)*b(1)-A(2,3)*b(2))/A(3,3)
!
!!      Solve transpose(L)*x = y
!         b(2) = b(2)-A(3,2)*b(3)
!!         L = ipvt(2)
!!         IF (L /= 2) THEN
!!            TMP = B(2); B(2) = B(L); B(L) = TMP
!!         END IF
!         b(1) = b(1)-A(3,1)*b(3)-A(2,1)*b(2)
!!         L = ipvt(1)
!!         IF (L /= 1) THEN
!!            TMP = B(1); B(1) = B(L); B(L) = TMP
!!         END IF
!   
!      END SELECT
!
!END SUBROUTINE SOL3

