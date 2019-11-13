!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Driver to test Hessian Singular Vectors
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM KPP_ROOT_Driver

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize

      KPP_REAL :: T, DVAL(NSPEC), DV(NVAR)
      KPP_REAL :: Hs(NVAR,NVAR), SOA(NVAR)
      INTEGER :: i, j, ind_1 = 1, ind_2 = 2, ind_COST
      LOGICAL :: LHessian
      INTEGER, PARAMETER :: NSOA = 1
      KPP_REAL :: Y_tlm(NVAR,NSOA)
      KPP_REAL :: Y_adj(NVAR,1)
      KPP_REAL :: Y_soa(NVAR,NSOA)
 
      integer          lwork, n
      parameter        (n = NVAR)
      
      ! Number of desired eigenvectors
      INTEGER, PARAMETER :: NEIG = 4
      double complex   zwork(2000), alpha(NEIG), beta(NEIG), &
                       eivec(n,NEIG), tmp(n), residu(n)

      integer          kmax, jmax, jmin, maxstep, method, m, l, maxnmv,  &
                       order, testspace
      double precision tol, lock, dznrm2
      logical          wanted
      double complex   target
      real             elapse
      real             etime, tarray(2)

      target = (0.0,0.0)
      tol = 1.d-9
      kmax = NEIG           ! Number of wanted solutions
      jmin = 10; jmax = 20  ! min/max size of search space
      
      maxstep = 30  ! Max number of Jacobi-Davidson iterations
      lock = 1.d-9
!...     order =  0: nearest to target
!...     order = -1: smallest real part
!...     order =  1: largest real part
!...     order = -2: smallest complex part
!...     order =  2: largest complex part
      order = 1

      method = 2  ! 1=gmres(m), 2=cgstab(l)
      m = 30      ! for gmres(m):
      l = 2       ! for cgstab(l):

      maxnmv = 20 !...     maximum number of matvecs in cgstab or gmres
      
      testspace = 3
!...     Testspace 1: w = "Standard Petrov" * v (Section 3.1.1)
!...     Testspace 2: w = "Standard 'variable' Petrov" * v (Section 3.1.2)
!...     Testspace 3: w = "Harmonic Petrov" * v (Section 3.5.1)

      if ( method .eq. 1 ) then
         lwork =  4 +  m  + 5*jmax + 3*kmax
      else
         lwork = 10 + 6*l + 5*jmax + 3*kmax
      end if
      wanted = .true.

      call jdqz(alpha, beta, eivec, wanted, n, target, tol, &
          kmax, jmax, jmin, &
          method, m, l, maxnmv, maxstep, &
          lock, order, testspace, zwork, lwork )

      elapse = etime( tarray )

      print*,'Number of converged eigenvalues = ',kmax
      
!...     Compute the norms of the residuals:
      do j = 1, kmax
         call amul  ( n, eivec(1,j), residu )
         call zscal ( n, beta(j), residu, 1)
         call bmul  ( n, eivec(1,j), tmp )
         call zaxpy( n, -alpha(j), tmp, 1, residu, 1 )
         print '("lambda(",i2,"): (",1p,e11.4,",",e11.4, &
             " )")', j,alpha(j)/beta(j)
         print '(a30,d13.6)', '||beta Ax - alpha Bx||:', &
               dznrm2( n, residu, 1 )
      end do
      write(*,10) tarray(1), elapse

   10 format(1x,'END JDQZ AFTER ',f6.2,' SEC. CPU-TIME AND ', f6.2, &
            ' SEC. ELAPSED TIME' )

      open(120,file='KPP_ROOT_vec.m')
      do i=1,NVAR
        write(120,120) (eivec(i,j), j=1,NEIG)
      end do
      close(120)	

      open(120,file='KPP_ROOT_val.m')
      do i=1,NEIG
        write(120,120) ALPHA(i), BETA(i), ALPHA(i)/BETA(i)
      end do
      close(120)	
      
120   format(10000(E14.7,1X))

END PROGRAM KPP_ROOT_Driver


!==========================================================================
!     Dummy preconditioner subroutine
!==========================================================================
      subroutine PRECON( neq, q )
!...............................................
!...     Subroutine to compute q = K^-1 q
!...............................................
      integer        neq, i
      double complex q(neq)
      end subroutine PRECON


!==========================================================================
!     matrix vector multiplication subroutine
!
!
!     Compute the matrix vector multiplication y<---A*x
!     where A = adjoint(Model)*tlm(Model)
!
!     Note: we multiply the real and imaginaray parts separately
!           x = (xr,xi)
!           y = (yr,yi)
!           We run the code twice such that:
!           yr <-- A*xr and yi <-- A*xi
!
!==========================================================================
      subroutine AMUL (n, u, v)
  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize
      integer           n, j
      Double Complex u(n), v(n)
      Double precision xr(n), yr(n),  xi(n), yi(n), three, two 
      parameter         (three = 3.0D+0, two = 2.0D+0)
      INTEGER, PARAMETER :: NSOA = 1
      KPP_REAL :: Y_tlm(NVAR,NSOA)
      KPP_REAL :: Y_adj(NVAR,1)
      KPP_REAL :: Y_soa(NVAR,NSOA)
      integer, save :: indexA = 0
      
      indexA = indexA + 1
      print*,'AMUL #',indexA

! Real and imaginary parts of complex inputs
      xr(1:n) = DBLE( u(1:n) )
      xi(1:n) = DIMAG( u(1:n) )


      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO
!~~~>  The cost function is 0.5*VAR(ind_COST, tF)**2
      ind_COST = ind_O3

! First do a FWD and a TLM integration
      CALL Initialize()
      Y_tlm(1:NVAR,1) = xr(1:NVAR)
      Y_adj(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)

! Set Lambda(t_final) = Y_tlm(t_final) and do backward integration
      CALL Initialize()
      Y_adj(1:NVAR,1) = Y_tlm(1:NVAR)
      Y_tlm(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      yr(1:NVAR) = Y_adj(1:NVAR)

! First do a FWD and a TLM integration
      CALL Initialize()
      Y_tlm(1:NVAR,1) = xi(1:NVAR)
      Y_adj(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)

! Set Lambda(t_final) = Y_tlm(t_final) and do backward integration
      CALL Initialize()
      Y_adj(1:NVAR,1) = Y_tlm(1:NVAR)
      Y_tlm(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      yi(1:NVAR) = Y_adj(1:NVAR)

! Convert double outputs to complex
      DO i=1,n
        v(i) = DCMPLX( yr(i), yi(i) )
      END DO	


      print*,'   (END AMUL #',indexA,')'
      
      end  subroutine AMUL

!==========================================================================
!     matrix vector multiplication subroutine
!
!     Compute the matrix vector multiplication y<---B*x
!     where B = Hessian, and y = SOA
!
!     Note: we multiply the real and imaginaray parts separately
!           x = (xr,xi)
!           y = (yr,yi)
!           We run the code twice such that:
!           yr <-- B*xr and yi <-- B*xi
!==========================================================================

      subroutine BMUL (n, u, v)
  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize
      integer           n, j
      Double Complex u(n), v(n)
      Double precision xr(n), yr(n), xi(n), yi(n)
      INTEGER, PARAMETER :: NSOA = 1
      KPP_REAL :: Y_tlm(NVAR,NSOA)
      KPP_REAL :: Y_adj(NVAR,1)
      KPP_REAL :: Y_soa(NVAR,NSOA)
      KPP_REAL :: Final(NVAR), FinalTlm(NVAR)
      integer, save :: indexB = 0
      
      indexB = indexB + 1
      print*,'BMUL #',indexB

! Real and imaginary parts of complex inputs
      xr(1:n) = DBLE( u(1:n) )
      xi(1:n) = DIMAG( u(1:n) )

      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO
!~~~>  The cost function is 0.5*VAR(ind_COST, tF)**2
      ind_COST = ind_O3

! Forward & TLM integration
      CALL Initialize()
      Y_tlm(1:NVAR,1) = xr
      Y_adj(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      Final(1:NVAR)    = VAR(1:NVAR)			    
      FinalTlm(1:NVAR) = Y_tlm(1:NVAR)			    
			    
! Adjoint and SOA integration, real part
      CALL Initialize()
      Y_adj(1:NVAR,1) = 0.0d0; Y_adj(ind_COST,1) = Final(ind_COST)
      Y_soa(1:NVAR,1) = 0.0d0; Y_soa(ind_COST,1) = FinalTlm(ind_COST)
      Y_tlm(1:NVAR,1) = xr
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
			    
      yr(1:NVAR) = Y_soa(1:NVAR,1)			    


! Forward & TLM integration
      CALL Initialize()
      Y_tlm(1:NVAR,1) = xi
      Y_adj(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      Final(1:NVAR)    = VAR(1:NVAR)			    
      FinalTlm(1:NVAR) = Y_tlm(1:NVAR)			    

! Adjoint and SOA integration, imaginary part
      CALL Initialize()
      Y_adj(1:NVAR,1) = 0.0d0; Y_adj(ind_COST,1) = Final(ind_COST)
      Y_soa(1:NVAR,1) = 0.0d0; Y_soa(ind_COST,1) = FinalTlm(ind_COST)
      Y_tlm(1:NVAR,1) = xi
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
			    
      yi(1:NVAR) = Y_soa(1:NVAR,1)			    

! Convert double outputs to complex
      DO i=1,n
        v(i) = DCMPLX( yr(i), yi(i) )
      END DO	


      print*,'   (END BMUL #',indexB,')'
      
      end subroutine BMUL 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


