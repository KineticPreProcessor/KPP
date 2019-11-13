!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Driver to test the Second Order Adjoint (SOA) model
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM KPP_ROOT_Driver

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize

      KPP_REAL :: T, DVAL(NSPEC), DV(NVAR)
      KPP_REAL :: ADJ1(NVAR), ADJ2(NVAR)
      KPP_REAL :: FinalRef(NVAR), FinalRef_tlm(NVAR)
      KPP_REAL :: FinalPert(NVAR), FinalPert_tlm(NVAR)
      KPP_REAL :: Hs(NVAR,NVAR), AdjTlm(NVAR,NVAR), SOA(NVAR)
      INTEGER :: i, j, ind_1 = 1, ind_2 = 2, ind_COST
      LOGICAL :: LHessian, LAdjTlm
 
!~~~>  Number of second order adjoints
!      i.e., number of vectors U_i s.t. Sigma_i = Hess*U_i
!      1 <= i <= NSOA

      INTEGER, PARAMETER :: NSOA = 1

      KPP_REAL :: Y_tlm(NVAR,NSOA)
      KPP_REAL :: Y_adj(NVAR,1)
      KPP_REAL :: Y_soa(NVAR,NSOA)
  
      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO
     
     
!~~~>  The cost function is 0.5*VAR(ind_COST, tF)**2
      ind_COST = ind_O3
     
     
!~~~>  Compute the full Hessian NVARxNVAR (or not)
      LHessian = .TRUE.
     
     
!~~~>  Compute the full Adjoint*Tlm NVARxNVAR (or not)
      LAdjTlm = .TRUE.
     
!~~~~~~~~~~~~~~~~~~~~~~~~
      PRINT*,'FIRST RUN, Reference'
      CALL Initialize()
      DV = 0.0d0
      DV(ind_NO2) =  0.02*VAR(ind_NO2)
      DV(ind_NO)  =  0.04*VAR(ind_NO)
!      DV(ind_PAN) = -0.2*VAR(ind_PAN)
!      DV(ind_O3)  =  0.1*VAR(ind_O3)
!      DO i=1,NVAR
!        DV(i) = rand()*VAR(i)
!      END DO	
      Y_tlm(1:NVAR,1) = DV
      Y_adj(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1:NSOA) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      FinalRef     = VAR
      FinalRef_tlm(1:NVAR) = Y_tlm(1:NVAR,1)
!~~~~~~~~~~~~~~~~~~~~~~~~   
     
!~~~~~~~~~~~~~~~~~~~~~~~~
      PRINT*,'FIRST RUN, Perturbed'
      CALL Initialize()
      VAR = VAR + DV
      Y_tlm(1:NVAR,1) = DV
      Y_adj(1:NVAR,1) = 0.0d0
      Y_soa(1:NVAR,1:NSOA) = 0.0d0
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      FinalPert     = VAR
      FinalPert_tlm(1:NVAR) = Y_tlm(1:NVAR,1)
!~~~~~~~~~~~~~~~~~~~~~~~~
    
!~~~~~~~~~~~~~~~~~~~~~~~~
      PRINT*,'PERTURBED RUN'
      CALL Initialize()
      VAR = VAR + DV
      Y_tlm(1:NVAR,1) = DV
      Y_adj(1:NVAR,1) = 0.0d0
      Y_adj(ind_COST,1) = FinalPert(ind_COST)
      Y_soa(1:NVAR,1) = 0.0d0
      Y_soa(ind_COST,1) = FinalPert_tlm(ind_COST)
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      F2 = VAR(ind_COST)**2/2
      ADJ2 = Y_adj(:,1)
      PRINT*,'PERTURBED FUN = ',F2
!~~~~~~~~~~~~~~~~~~~~~~~~      
      
!~~~~~~~~~~~~~~~~~~~~~~~~
      PRINT*,'REFERENCE RUN'
      CALL Initialize()
      Y_tlm(1:NVAR,1) = DV
      Y_adj(1:NVAR,1) = 0.0d0
      Y_adj(ind_COST,1) = FinalRef(ind_COST)
      Y_soa(1:NVAR,1) = 0.0d0
      Y_soa(ind_COST,1) = FinalRef_tlm(ind_COST)
      CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
      F1 = VAR(ind_COST)**2/2
      ADJ1 = Y_adj(:,1)
      SOA  = Y_soa(:,1)
      PRINT*,'REFERENCE FUN = ',F1
!~~~~~~~~~~~~~~~~~~~~~~~~   
      
!~~~~~~~~~~~~~~~~~~~~~~~~
      IF (LHessian) THEN
      DO i = 1,NVAR
        PRINT*,'.... FULL HESSIAN column ',i,' of ',NVAR
        CALL Initialize()
        Y_tlm(1:NVAR,1) = 0.0d0
        Y_tlm(i,1)      = 1.0d0
        Y_adj(1:NVAR,1) = 0.0d0
        Y_adj(ind_COST,1) = FinalRef(ind_COST)
        Y_soa(1:NVAR,1) = 0.0d0
        Y_soa(ind_COST,1) = FinalRef_tlm(ind_COST)
        CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
        CALL Initialize()
        Y_soa(1:NVAR,1) = 0.0d0
        Y_soa(ind_COST,1) = Y_tlm(ind_COST,1)
        Y_tlm(1:NVAR,1) = 0.0d0
        Y_tlm(i,1)      = 1.0d0
        Y_adj(1:NVAR,1) = 0.0d0
        Y_adj(ind_COST,1) = FinalRef(ind_COST)
        CALL INTEGRATE_soa(Nsoa, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
	Hs(1:NVAR,i) = Y_soa(1:NVAR,1)
      END DO
      ! Write Hessian in file
      OPEN(150,FILE='KPP_ROOT_Hess.m')
      DO i=1,NVAR
        WRITE(150,224) (Hs(i,j),j=1,NVAR)
      END DO  
      CLOSE(150)      
      END IF
!~~~~~~~~~~~~~~~~~~~~~~~~
    
!~~~~~~~~~~~~~~~~~~~~~~~~
      IF (LAdjTlm) THEN
      DO i = 1,NVAR
        PRINT*,'.... FULL ADJ*TLM column ',i,' of ',NVAR
        CALL Initialize()
        Y_tlm(1:NVAR,1) = 0.0d0
        Y_tlm(i,1)      = 1.0d0
        Y_adj(1:NVAR,1) = 0.0d0
        Y_adj(ind_COST,1) = FinalRef(ind_COST)
        Y_soa(1:NVAR,1) = 0.0d0
        Y_soa(ind_COST,1) = FinalRef_tlm(ind_COST)
        CALL INTEGRATE_SOA( NSOA, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
        CALL Initialize()
        Y_adj(1:NVAR,1) = Y_tlm(1:NVAR)
        Y_soa(1:NVAR,1) = 0.0d0
        Y_tlm(1:NVAR,1) = 0.0d0
        CALL INTEGRATE_soa(Nsoa, VAR, Y_tlm, Y_adj, Y_soa, TSTART, TEND, &
                            ATOL, RTOL)
	AdjTlm(1:NVAR,i) = Y_adj(1:NVAR,1)
      END DO
      ! Write Hessian in file
      OPEN(150,FILE='KPP_ROOT_AdjTlm.m')
      DO i=1,NVAR
        WRITE(150,224) (AdjTlm(i,j),j=1,NVAR)
      END DO  
      CLOSE(150)      
      END IF
!~~~~~~~~~~~~~~~~~~~~~~~~

      WRITE(6,*) 'FD versus SOA:'      
      DO i=1,NVAR
        WRITE(6,222) ADJ2(i)-ADJ1(i), SOA(i),           &
	  ABS(ADJ2(i)-ADJ1(i)-SOA(i))/                  &
          MAX(ABS(ADJ2(i)-ADJ1(i)),ABS(SOA(i)),1.d-16), &
	  TRIM(SPC_NAMES(i))
      END DO  
222   FORMAT('FD=',E12.5,'   soa=',E12.5,'  Diff=',E12.5,'  (',A,')')      

      WRITE(6,FMT="(A14,1X,E14.7)") 'F2-F1 = ',F2-F1
      WRITE(6,FMT="(A14,1X,E14.7)") 'First = ',DOT_PRODUCT(ADJ1,DV)      
      WRITE(6,FMT="(A14,1X,E14.7)") 'F2-F1-First = ',F2-F1-DOT_PRODUCT(ADJ1,DV)      
      WRITE(6,FMT="(A14,1X,E14.7)") 'Second = ',DOT_PRODUCT(SOA,DV)/2      

      IF (LHessian) THEN
      WRITE(6,*) 'The Hessian:'      
      DO i=1,NVAR
        WRITE(6,223) (Hs(i,j),j=1,NVAR)
      END DO  
223   FORMAT(100(E12.5,1X))
      END IF


      OPEN(150,FILE='KPP_ROOT_Cost.m')
      WRITE(150,FMT="(A14,1X,E14.7)") 'F2-F1 = ',F2-F1
      WRITE(150,FMT="(A14,1X,E14.7)") 'First = ',DOT_PRODUCT(ADJ1,DV)      
      WRITE(150,FMT="(A14,1X,E14.7)") 'F2-F1-First = ',F2-F1-DOT_PRODUCT(ADJ1,DV)      
      WRITE(150,FMT="(A14,1X,E14.7)") 'Second = ',DOT_PRODUCT(SOA,DV)/2      
      CLOSE(150)      

      OPEN(150,FILE='KPP_ROOT_FD.m')
      DO i=1,NVAR
        WRITE(150,224) ADJ2(i)-ADJ1(i), SOA(i),        &
	  ABS(ADJ2(i)-ADJ1(i)-SOA(i))/                 &
          MAX(ABS(ADJ2(i)-ADJ1(i)),ABS(SOA(i)),1.d-16)
      END DO  
      CLOSE(150)      

      OPEN(150,FILE='KPP_ROOT_Sol.m')
      DO i=1,NVAR
        WRITE(150,224) FinalRef(i)
      END DO  
      CLOSE(150)      

      OPEN(150,FILE='KPP_ROOT_Adj.m')
      DO i=1,NVAR
        WRITE(150,224) ADJ1(i)
      END DO  
      CLOSE(150)      

224   FORMAT(10000(E24.14,1X))      


END PROGRAM KPP_ROOT_Driver


