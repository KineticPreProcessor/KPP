 #define MAX(a,b) ( ((a) >= (b)) ?(a):(b)  )
 #define MIN(b,c) ( ((b) <  (c)) ?(b):(c)  )
 #define ABS(x)   ( ((x) >=  0 ) ?(x):(-x) )
 #define SQRT(d)  ( pow((d),0.5)  )
 /* SIGN transfer function */
 #define SIGN(x,y) (((y) >=  0 ) ?(ABS(x)):(-ABS(x)) ) 

/*~~> Numerical constants */
 #define  ZERO     (KPP_REAL)0.0
 #define  ONE      (KPP_REAL)1.0

/*~~~> Statistics on the work performed by the Runge-Kutta method */
 #define Nfun	0
 #define Njac	1
 #define Nstp	2
 #define Nacc	3
 #define Nrej	4
 #define Ndec	5
 #define Nsol	6
 #define Nsng	7
 #define Ntexit	0
 #define Nhacc	1
 #define Nhnew	2

/*~~~> Runge-Kutta method parameters */
 #define RKmax 3

 #define R2A 1
 #define R1A 2
 #define L3C 3
 #define GAU 4
 #define L3A 5

 int rkMethod,
     SdirkError;
 KPP_REAL rkT[RKmax][RKmax],
	rkTinv[RKmax][RKmax],
       	rkTinvAinv[RKmax][RKmax],
   	rkAinvT[RKmax][RKmax],
   	rkA[RKmax+1][RKmax+1],
   	rkB[RKmax+1],
      	rkC[RKmax+1],
       	rkD[RKmax+1],
       	rkE[RKmax+1],
   	rkBgam[RKmax+2],
       	rkBhat[RKmax+2],
       	rkTheta[RKmax+1],
       	rkF[RKmax+2],
   	rkGamma,
   	rkAlpha,
   	rkBeta,
	rkELO;
/*~~~> Function headers */
// void INTEGRATE(KPP_REAL TIN, KPP_REAL TOUT, int ICNTRL_U[], KPP_REAL RCNTRL_U[],
//		int ISTATUS_U[], KPP_REAL RSTATUS_U[], int IERR_U); 
 void INTEGRATE(KPP_REAL TIN, KPP_REAL TOUT); 
 void RungeKutta(int N, KPP_REAL T, KPP_REAL Tend, KPP_REAL Y[],
     	       KPP_REAL RelTol[], KPP_REAL AbsTol[], KPP_REAL RCNTRL[],
	       int ICNTRL[], KPP_REAL RSTATUS[], int ISTATUS[], int* IERR);
 void RK_Integrator(int N, KPP_REAL*  T, KPP_REAL  Tend, KPP_REAL Y[],
     		   KPP_REAL AbsTol[], KPP_REAL RelTol[], int ITOL,
     		   int ISTATUS[], KPP_REAL RSTATUS[], KPP_REAL Hmin,
		   KPP_REAL Hmax, KPP_REAL Hstart, KPP_REAL Roundoff,
		   int Max_no_steps, int NewtonMaxit, int StartNewton,
		   int Gustafsson, KPP_REAL ThetaMin, KPP_REAL NewtonTol,
		   KPP_REAL FacSafe, KPP_REAL FacMax, KPP_REAL FacMin,
     		   KPP_REAL FacRej, KPP_REAL Qmin, KPP_REAL Qmax, int* IERR);
 void RK_ErrorMsg(int Code, KPP_REAL T, KPP_REAL H, int* IERR);
 void RK_ErrorScale(int N, int ITOL, KPP_REAL AbsTol[], KPP_REAL RelTol[], 
		KPP_REAL Y[], KPP_REAL SCAL[]);
 /*void RK_Transform(int N, KPP_REAL Tr[][RKmax], KPP_REAL Z1[], KPP_REAL Z2[],
                KPP_REAL Z3[], KPP_REAL W1[], KPP_REAL W2[], KPP_REAL W3[]);*/
 void RK_Interpolate(char action[], int N, KPP_REAL H, KPP_REAL Hold,
	 KPP_REAL Z1[], KPP_REAL Z2[], KPP_REAL Z3[], KPP_REAL CONT[][RKmax]);
 void RK_PrepareRHS(int N, KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL FO[],
        KPP_REAL Z1[], KPP_REAL Z2[], KPP_REAL Z3[], KPP_REAL R1[], KPP_REAL R2[],
        KPP_REAL R3[]);
 void RK_Decomp(int N, KPP_REAL H, KPP_REAL FJAC[], KPP_REAL E1[],
                int IP1[], KPP_REAL E2R[], KPP_REAL E2I[], 
		int IP2[], int* ISING, int ISTATUS[]);
 void RK_Solve(int N, KPP_REAL H, KPP_REAL E1[], int IP1[], KPP_REAL E2R[],
	KPP_REAL E2I[], int IP2[], KPP_REAL R1[], KPP_REAL R2[], KPP_REAL R3[],
	int ISTATUS[]);
 void RK_ErrorEstimate(int N, KPP_REAL H, KPP_REAL T, KPP_REAL Y[],
        KPP_REAL FO[], KPP_REAL E1[], int IP1[], KPP_REAL Z1[],
        KPP_REAL Z2[], KPP_REAL Z3[], KPP_REAL SCAL[], KPP_REAL* Err,
        int FirstStep, int Reject, int ISTATUS[]);
 void Radau2A_Coefficients();
 void Lobatto3C_Coefficients ();
 void Gauss_Coefficients();
 void Radau1A_Coefficients();
 void Lobatto3A_Coefficients();
 void FUN_CHEM(KPP_REAL T, KPP_REAL V[], KPP_REAL FCT[]);
 void JAC_CHEM(KPP_REAL T, KPP_REAL V[], KPP_REAL JF[]);
 KPP_REAL RK_ErrorNorm(int N, KPP_REAL SCAL[], KPP_REAL DY[]);
 void Fun(KPP_REAL Y[], KPP_REAL FIX[], KPP_REAL RCONST[], KPP_REAL Ydot[]);
 void Jac_SP(KPP_REAL Y[], KPP_REAL FIX[], KPP_REAL RCONST[], KPP_REAL Ydot[]);
 void WCOPY(int N, KPP_REAL X[], int incX, KPP_REAL Y[], int incY);
 void WADD(int N, KPP_REAL Y[], KPP_REAL Z[], KPP_REAL TMP[]);
 void WAXPY(int N, KPP_REAL Alpha, KPP_REAL X[], int incX, KPP_REAL Y[], int incY );
 KPP_REAL WLAMCH( char C );
 void KppSolveCmplxR(KPP_REAL JVSR[], KPP_REAL JVSI[], KPP_REAL XR[], KPP_REAL XI[]);
 int KppDecompCmplxR(KPP_REAL *JVSR, KPP_REAL *JVSI);
 void Set2Zero(int N, KPP_REAL A[]);
 int KppDecomp( KPP_REAL A[] );
 void KppSolve ( KPP_REAL A[], KPP_REAL b[] );
 void Update_SUN();
 void Update_RCONST();
 void Update_PHOTO();
 
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//void INTEGRATE(KPP_REAL TIN, KPP_REAL TOUT, int ICNTRL_U[], KPP_REAL RCNTRL_U[],
//	       int ISTATUS_U[], KPP_REAL RSTATUS_U[], int IERR_U )
void INTEGRATE(KPP_REAL TIN, KPP_REAL TOUT )
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* RungeKutta - Fully Implicit 3-stage Runge-Kutta methods based on:
          * Radau-2A   quadrature (order 5)                        
          * Radau-1A   quadrature (order 5)                              
          * Lobatto-3C quadrature (order 4)                              
          * Gauss      quadrature (order 6)                              
  By default the code employs the KPP sparse linear algebra routines     
  Compile with -DFULL_ALGEBRA to use full linear algebra (LAPACK)        
                                                                         
    (C)  Adrian Sandu, August 2005                                       
    Virginia Polytechnic Institute and State University                  
    Contact: sandu@cs.vt.edu                                             
    Revised by Philipp Miehe and Adrian Sandu, May 2006  
    F90 to C translation by Tinting Jiang and Don Jacob, July 2006                
    This implementation is part of KPP - the Kinetic PreProcessor        
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
   int IERR;
   KPP_REAL RCNTRL[20],
  	  RSTATUS[20],
	  T1,
	  T2;
   int ICNTRL[20],
       ISTATUS[20];
   static int Ntotal = 0; /* for printing the number of steps */
   
   int i;
   for ( i = 0; i < 20; i++ ) {
     RCNTRL[i] = ZERO;
     ICNTRL[i] = 0;
   } /* for */

   /*~~~> fine-tune the integrator: */
   ICNTRL[1]  = 0;    /* 0=vector tolerances, 1=scalar tolerances */
   ICNTRL[4]  = 8;    /* Max no. of Newton iterations */
   ICNTRL[5]  = 0;    /* Starting values for Newton are interpolated(0) or zero(1) */
   ICNTRL[9]  = 1;    /* 0 - classic or 1 - SDIRK error estimation */
   ICNTRL[10] = 0;    /* Gustaffson(0) or classic(1) controller */ 
   
//   /*~~~> if optional parameters are given, and if they are >0,
//          then use them to overwrite default settings */	
//   if (ICNTRL_U != NULL) {
//	for ( i = 0; i < 20; i++) {
//	   if (ICNTRL_U[i] > 0) {
//		  ICNTRL[i] = ICNTRL_U[i];
//	   } /* end if */
//	} /* end for */
//   } /* end if */
//   if (RCNTRL_U != NULL) {
//	for (i = 0; i < 20; i++) {
//   	   if (RCNTRL_U[i] > 0) {
//		RCNTRL[i] = RCNTRL_U[i];
//	   } /* end if */
//	} /* end for */
//   } /* end if */

   T1 = TIN;
   /*printf("T1=%f\n", T1);*/
   T2 = TOUT;
   /*printf("T2=%f\n", T2);*/
   RungeKutta(NVAR, T1, T2, VAR, RTOL, ATOL, RCNTRL,ICNTRL,RSTATUS,ISTATUS, &IERR);
   
   Ntotal += ISTATUS[Nstp];
   printf("NSTEPS=%d (%d)  O3=%E ", ISTATUS[Nstp], Ntotal, VAR[ind_O3]);

//   /* if optional parameters are given for output
//    * use them to store information in them */
//   if (ISTATUS_U != NULL) {
//	for (i = 0; i < 20; i++) {
//	   ISTATUS_U[i] = ISTATUS[i];
//	} /* end for */
//    } /* end if */
//   if (RSTATUS_U != NULL) {
//	for (i = 0; i < 20; i++) {
//	   RSTATUS_U[i] = RSTATUS[i];
//	} /* end for */
//   } /* end if */
//   /*if (IERR_U != NULL) */
//	IERR_U = IERR;

   if (IERR < 0) {
	printf("Runge-Kutta: Unsuccessful exit at T=%f(IERR=%d)", TIN, IERR);
   } /* end if */
   
} /* end INTEGRATE */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void RungeKutta(int N, KPP_REAL T, KPP_REAL Tend, KPP_REAL Y[],
     	       KPP_REAL RelTol[], KPP_REAL AbsTol[], KPP_REAL RCNTRL[],
	       int ICNTRL[], KPP_REAL RSTATUS[], int ISTATUS[], int* IERR)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This implementation is based on the book and the code Radau5:

     	E. HAIRER AND G. WANNER
     	"SOLVING ORDINARY DIFFERENTIAL EQUATIONS II.
       	     STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS."
     	SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
        SPRINGER-VERLAG (1991)

	UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
	CH-1211 GENEVE 24, SWITZERLAND
    	E-MAIL:  HAIRER@DIVSUN.UNIGE.CH,  WANNER@DIVSUN.UNIGE.CH
   
    Methods:
	 * Radau-2A	quadrature (order 5)
      	 * Radau-1A	quadrature (order 5)
      	 * Lobatto-3C	quadrature (order 4)
    	 * Gauss	quadrature (order 6)

    (C)  Adrian Sandu, August 2005
    Virginia Polytechnic Institute and State University
    Contact: sandu@cs.vt.edu
    This implementation is part of KPP - the Kinetic PreProcessor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  *~~~>   INPUT ARGUMENTS:
	  ----------------
    Note: For input parameters equal to zero the default values of the
	  corresponding variables are used.

     N		Dimension of the system
     T		Initial time value

     Tend	Final T value (Tend-T may be positive or negative)

     Y(N)	Initial values for Y

     RelTol, AbsTol   Relative and absolute error tolerances.
           for ICNTRL[1] = 0: AbsTol, RelTol are N-dimensional vectors
                         = 1: AbsTol, RelTol are scalars
  *~~~>   Integer input parameters:

     ICNTRL[0] = not used
     
     ICNTRL[1] = 0: AbsTol, RelTol are NVAR-dimensional vectors
               = 1: AbsTol, RelTol are scalars

     ICNTRL[2] = RK method selection
	     = 1: Radau-2A	(the default)
	     = 2: Radau-1A
	     = 3: Lobatto-3C
	     = 4: Gauss
	     = 5: Lobatto-3A (not yet implemented)

     ICNTRL[3] -> maximum number of integration steps
	For ICNTRL[3] = 0 the default value of 10000 is used

     ICNTRL[4] -> maximum number of Newton iterations
	For ICNTRL[4] = 0 the default value of 8 is used

     ICNTRL[5] -> starting values of Newton iterations:
	ICNTRL[5] = 0: starting values are obtained from
		     the extrapolated collocation solution
		     (the default)
   	ICNTRL[5] = 1: starting values are zero

     ICNTRL[9] -> switch for error estimation strategy
	ICNTRL[9] = 0: one additional stage at c = 0,
		     see Hairer (default)
	ICNTRL[9] = 1: two additional stages at c = 0,
		     and SDIRK at c = 1, stiffly accurate

     ICNTRL[10] -> switch for step size strategy
	ICNTRL[10] = 0: mod. predictive controller (Gustafsson, default)
	ICNTRL[10] = 1: classical step size control
 	the choice 1 seems to produce safer results;
	for simple problems, the choice 2 produces
	often slightly faster runs

  *~~~>   Real input parameters:
     RCNTRL[0]  -> Hmin, lower bound for the integration step size
		(highly recommended to keep Hmin = ZERO, the default)

     RCNTRL[1]  -> Hmax, upper bound for the integration step size

     RCNTRL[2]  -> Hstart, the starting step size

     RCNTRL[3]  -> FacMin, lower bound on step decrease factor 
		(default=0.2)

     RCNTRL[4]  -> FacMax, upper bound on step increase factor
		(default=6)

     RCNTRL[5]  -> FacRej, step decrease factor after multiple rejections
		(default=0.1)

     RCNTRL[6]  -> FacSafe, by which the new step is slightly smaller
		than the predicted value (default=0.9)

     RCNTRL[7]  -> ThetaMin. If Newton convergence rate smaller
		than ThetaMin the Jacobian is not recomputed;
		(default=0.001)

     RCNTRL[8]  -> NewtonTol, stopping criterion for Newton's method
		(default=0.03)

     RCNTRL[9]  -> Qmin

     RCNTRL[10] -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
		 step size is kept constant and the LU factorization
		 reused (default Qmin=1, Qmax=1.2)	
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  *~~~>     OUTPUT ARGUMENTS:
            -----------------
     T	      -> T value for which the solution has been computed
		 (after successful return T=Tend).

     Y(N)  -> Numerical solution at T

  Note: each call to RungeKutta adds the current no. of fcn calls
	to previous value of ISTATUS[0], and similar to other params.
	Set ISTATUS[1:9] = 0 before call to avoid this accumulation.

     ISTATUS[0] -> No. of function calls
     ISTATUS[1] -> No. of Jacobian calls
     ISTATUS[2] -> No. of steps
     ISTATUS[3] -> No. of accepted steps
     ISTATUS[4] -> No. of rejected steps (except at very beginning)
     ISTATUS[5] -> No. of LU decompositions
     ISTATUS[6] -> No. of forward/backward substitutions
     ISTATUS[7] -> No. of singular matrix decompositions

     RSTATUS[0] -> Texit, the time corresponding to the 
		 computed Y upon return
     RSTATUS[1] -> Hexit, last accepted step before exit
     RSTATUS[2] -> Hnew, last predicted step (not yet taken)
		 For multiple restarts, use Hnew as Hstart
  		 in the subsequent run

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  *~~~>    RETURN VALUE (int):

     IERR     -> Reports on successfulness upon return:
	= 1 for success
	< 0 for error (value equals error code)
	= -1  : Improper value for maximal no of steps
	= -2  : Improper value for maximal no of Newton iterations
	= -3  : Hmin/Hmax/Hstart must be positive
    	= -4  : Improper values for FacMin/FacMax/FacSafe/FacRej
    	= -5  : Improper value for ThetaMin
    	= -6  : Newton stopping tolerance too small
    	= -7  : Improper values for Qmin, Qmax
    	= -8  : Tolerances are too small
    	= -9  : No of steps exceeds maximum bound
    	= -10 : Step size too small
	= -11 : Matrix is repeatedly singular
   	= -12 : Non-convergence of Newton iterations
    	= -13 : Requested RK method not implemented
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
{
  /*~~~>  Control arguments */
  /*printf("Starting RungeKutta\n");*/
   int Max_no_steps;
   KPP_REAL Hmin,
  	  Hmax,
	  Hstart,
	  Qmin,
	  Qmax,
   	  Roundoff,
	  ThetaMin,
	  NewtonTol,
   	  FacSafe,
          FacMin,
	  FacMax,
	  FacRej;
  /*~~~> Local variables */
   int NewtonMaxit,
       ITOL,
       i,
       StartNewton,
       Gustafsson; 

   *IERR = 0;
   for (i = 0; i < 20; i++) {
	ISTATUS[i] = 0;
	RSTATUS[i] = ZERO;
   } /* end for */
   
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  /*~~~> ICNTRL[0] - autonomous system - not used */
  /*~~~> ITOL: 1 for vector and 0 for scalar AbsTol/RelTol */
   if (ICNTRL[1] == 0) {
	/*printf("Entering if ICNTRL[1] == 0\n");*/
	ITOL = 1;
   }
   else {
	ITOL = 0;
   } /*end if */
  /*~~~> Error control selection */
   if (ICNTRL[9] == 0) {
	SdirkError = 0;
   }
   else {
	/*printf("Entering if ICNTRL[9] == 1\n");*/    
	SdirkError = 1;
   } /* end if */
  /*~~~> Method Selection */
  /*printf("Starting case ICNTRL[2]\n");*/
   switch (ICNTRL[2]) {
	case 0:
	case 1:	Radau2A_Coefficients(); /*printf("case 0 or 1\n");*/
		break;
	case 2: Lobatto3C_Coefficients(); printf("case 2\n");
		break;
	case 3: Gauss_Coefficients(); printf("case 3\n");
		break;
	case 4: Radau1A_Coefficients(); printf("case 4\n");
		break;
	case 5: Lobatto3A_Coefficients(); printf("case 5\n");
		break;
	default:
		printf("\n ICNTRL[2]=%d\n", ICNTRL[2]);
		RK_ErrorMsg(-13, T, ZERO, IERR);
   } /* end switch */
 
 /*~~~> Max_no_steps: the maximal number of time steps */
   if (ICNTRL[3] == 0) {
	Max_no_steps = 200000;
	/*printf("Max_no_steps = %d\n", Max_no_steps);*/
   }
   else {
	Max_no_steps = ICNTRL[3];
	if (Max_no_steps <= 0)
	{
		printf("\n ICNTRL[3]=%d\n", ICNTRL[3]);
		RK_ErrorMsg(-1, T, ZERO, IERR);
	} /* end if */
   } /* end if */

 /*~~~> NewtonMaxit: maximal number of Newton iterations */
   if (ICNTRL[4] == 0)
	NewtonMaxit = 8;
   else 
   {
	NewtonMaxit = ICNTRL[4];
	if (NewtonMaxit <= 0) {
		printf("\n ICNTRL[4]=%d\n", ICNTRL[4]);
		RK_ErrorMsg(-2, T, ZERO, IERR);
	}
	/*printf("NewtonMaxit = %d\n", NewtonMaxit);*/
   } /* end if */

 /*~~~> StartNewton: Use extrapolation for starting values of Newton iterations */
   if (ICNTRL[5] == 0) {
	StartNewton = 1;
	/*printf("StartNewton = %d\n", StartNewton);*/
   }
   else {
	StartNewton = 0;
   } /* end if */

 /*~~~> Gustafsson: step size controller */
   if (ICNTRL[10] == 0) {
	Gustafsson = 1;
	/*printf("Gustafsson = %d\n", Gustafsson);*/
   }
   else {
	Gustafsson = 0;
   } /* end if */
  
 /*~~~> Roundoff: smallest number s.t. 1.0 + Roundoff > 1.0 */
   Roundoff = WLAMCH('E');
   /*printf("Roundoff=%f\n", Roundoff);*/
 /*~~~> Hmin = minimal step size */
   if (RCNTRL[0] == ZERO) {
	/*printf("Entering if RCNTRL[0]=%f == ZERO\n", RCNTRL[0]);*/   
	Hmin = ZERO;
   }
   else {
	Hmin = MIN(ABS(RCNTRL[0]), ABS(Tend-T));
   } /* end if */
 /*~~~> Hmax = maximal step size */
   if (RCNTRL[1] == ZERO) {
	/*printf("Entering if RCNTRL[1]=%f == ZERO\n", RCNTRL[1]);*/   
      	Hmax = ABS(Tend-T);
   }
   else {
      	Hmax = MIN(ABS(RCNTRL[1]),ABS(Tend-T));
   } /* end if */
 /*~~~> Hstart = starting step size */
   if (RCNTRL[2] == ZERO) {
	/*printf("Entering if RCNTRL[2]=%f == ZERO\n", RCNTRL[2]);*/   
    	Hstart = ZERO;
   }
   else {
	Hstart = MIN(ABS(RCNTRL[2]),ABS(Tend-T));
   } /* end if */
 /*~~~> FacMin: lower bound on step decrease factor */
   if (RCNTRL[3] == ZERO) {
	/*printf("Entering if RCNTRL[3]=%f == ZERO\n", RCNTRL[3]);*/   
	FacMin = (KPP_REAL)0.2;
   }
   else {
   	FacMin = RCNTRL[3];
   } /* end if */
 /*~~~> FacMax: upper bound on step increase factor */
   if (RCNTRL[4] == ZERO) {
	/*printf("Entering if RCNTRL[4]=%f == ZERO\n", RCNTRL[4]);*/   
    	FacMax = (KPP_REAL)8.0;
   }
   else {
    	FacMax = RCNTRL[4];
   } /* end if */
 /*~~~> FacRej: step decrease factor after 2 consecutive rejections */
   if (RCNTRL[5] == ZERO) {
	/*printf("Entering if RCNTRL[5]=%f == ZERO\n", RCNTRL[5]);*/
	FacRej = (KPP_REAL)0.1;
   }
   else {
     	FacRej = RCNTRL[5];
   } /* end if */
 /*~~~> FacSafe: by which the new step is slightly smaller
		 than the predicted value */
   if (RCNTRL[6] == ZERO) {
	/*printf("Entering if RCNTRL[6]=%f == ZERO\n", RCNTRL[6]);*/   
    	FacSafe = (KPP_REAL)0.9;
   }
   else {
 	FacSafe = RCNTRL[6];
   } /* end if */
   if ((FacMax < ONE) || (FacMin > ONE) || (FacSafe <= 1.0e-03) || (FacSafe >= ONE)) {
	printf("\n RCNTRL[3]=%f, RCNTRL[4]=%f, RCNTRL[5]=%f, RCNTRL[6]=%f\n",
		RCNTRL[3], RCNTRL[4], RCNTRL[5], RCNTRL[6]);
	RK_ErrorMsg(-4, T, ZERO, IERR);
   } /* end if */

 /*~~~> ThetaMin: decides whether the Jacobian should be recomputed */
   if (RCNTRL[7] == ZERO)
  	ThetaMin = (KPP_REAL)1.0e-03;
   else {
	ThetaMin = RCNTRL[7];
	if (ThetaMin <= (KPP_REAL)0.0 || ThetaMin >= (KPP_REAL)1.0) {
		printf("\n RCNTRL[7]=%f\n", RCNTRL[7]);
		RK_ErrorMsg(-5, T, ZERO, IERR);
	}
   } /* end if */
 /*~~~> NewtonTol: stopping criterion for Newton's method */
   if (RCNTRL[8] == ZERO)
	NewtonTol = (KPP_REAL)3.0e-02;
   else {
	NewtonTol = RCNTRL[8];
  	if (NewtonTol <= Roundoff) {
		printf("\n RCNTRL[8]=%f\n", RCNTRL[8]);
		RK_ErrorMsg(-6, T, ZERO, IERR);
	}
   } /* end if */
 /*~~~> Qmin AND Qmax: IF Qmin < Hnew/Hold < Qmax then step size = const. */
   if (RCNTRL[9] == ZERO) {
	Qmin = ONE;
   }
   else {
	Qmin = RCNTRL[9];
   } /* end if */
   if (RCNTRL[10] == ZERO) {
	Qmax = (KPP_REAL)1.2;
   }
   else {
	Qmax = RCNTRL[10];
   } /* end if */
   if (Qmin > ONE || Qmax < ONE) {
	printf("\n RCNTRL[9]=%f\n", Qmin);
	printf("\n RCNTRL[10]=%f\n", Qmax);
	RK_ErrorMsg(-7, T, ZERO, IERR);
   } /* end if */
 /*~~~> Check if tolerances are reasonable */
   if (ITOL == 0) {
	if ( AbsTol[0] <= ZERO || RelTol[0] <= ( ((KPP_REAL)10.0)*Roundoff) ) {
	   printf("\n AbsTol=%f\n", AbsTol[0]);
	   printf("\n RelTol=%f\n", RelTol[0]);
	   RK_ErrorMsg(-8, T, ZERO, IERR);
	}
   }
   else {
	for (i = 0; i < N; i++) {
	   if ( (AbsTol[i] <= ZERO) || (RelTol[i] <= ((KPP_REAL)10.0)*Roundoff) ) {
                printf("\n AbsTol[%d] = %f\n", i, AbsTol[i]);
		printf("\n AbsTol[%d] = %f\n", i, RelTol[i]);
		RK_ErrorMsg(-8, T, ZERO, IERR);
	   } /* end if */
	} /* end for */
   } /* end if */

 /*~~~> Parameters are wrong */
   if (*IERR < 0)
	return;

 /*~~~> Call the core method */	
   RK_Integrator(N, &T, Tend, Y, AbsTol, RelTol, ITOL, ISTATUS, RSTATUS,
		   Hmin, Hmax, Hstart, Roundoff, Max_no_steps, NewtonMaxit,
		   StartNewton, Gustafsson, ThetaMin, NewtonTol,
		   FacSafe, FacMax, FacMin, FacRej, Qmin, Qmax, IERR);
  /*printf("*IERR = %d\n", *IERR);*/
  /*printf("Ending RungeKutta\n");*/
} /* RungeKutta */
   
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void RK_Integrator( int N, 
	/*~~~> Input: integration interval */        
     	KPP_REAL*  T, KPP_REAL  Tend, KPP_REAL Y[],
     	KPP_REAL AbsTol[], KPP_REAL RelTol[], int ITOL,
  /*~~~> Input: the initial condition at T; output: the solution at Tend */     
     	int ISTATUS[], KPP_REAL RSTATUS[], KPP_REAL Hmin, KPP_REAL Hmax, 
	KPP_REAL Hstart, KPP_REAL Roundoff, int Max_no_steps, int NewtonMaxit,
       	int StartNewton, int Gustafsson, KPP_REAL ThetaMin,
       	KPP_REAL NewtonTol, KPP_REAL FacSafe, KPP_REAL FacMax, KPP_REAL FacMin,
     	KPP_REAL FacRej, KPP_REAL Qmin, KPP_REAL Qmax, int* IERR)
{  
   /*printf("Starting RK_Integrator\n");*/
   KPP_REAL FJAC[LU_NONZERO],
  	  E1[LU_NONZERO],
	  E2R[LU_NONZERO],
	  E2I[LU_NONZERO];
   KPP_REAL Z1[NVAR],
  	  Z2[NVAR],
	  Z3[NVAR],
	  Z4[NVAR],
	  SCAL[NVAR],
	  DZ1[NVAR],
	  DZ2[NVAR],
	  DZ3[NVAR],
	  DZ4[NVAR],
	  G[NVAR],
	  TMP[NVAR],
	  FO[NVAR];   
   KPP_REAL CONT[NVAR][RKmax],
  	  Tdirection,
	  H,
	  Hacc,
	  Hnew,
	  Hold,
	  Fac,
	  FacGus, 
   	  Theta,
	  Err,
	  ErrOld,
	  NewtonRate,
	  NewtonIncrement,
	  Hratio,
	  Qnewton, 
	  NewtonPredictedErr,
	  NewtonIncrementOld,
	  ThetaSD;
   int IP1[NVAR],
       IP2[NVAR],
       NewtonIter,
       Nconsecutive,
       i;
   int ISING;
   int Reject,
       FirstStep,
       SkipJac,
       NewtonDone,
       SkipLU;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
  /*~~~>  INITIAL setting */
   /*printf("ONE=%f Tend-*T=%f\n", ONE, Tend-*T);*/
   Tdirection = SIGN(ONE, Tend-*T);
   /*printf("Tdirection=%f\n", Tdirection);*/
   /*printf("Hmin=%g Hstart=%g \n", Hmin, Hstart);*/
   /*printf("ABS(Hmin)=%g ABS(Hstart)=%g \n", ABS(Hmin), ABS(Hstart));*/
   H = MIN( MAX(ABS(Hmin), ABS(Hstart)), Hmax ); 
   if (ABS(H) <= ((KPP_REAL)10.0*Roundoff)) {
	/*printf("Entering if ABS(H)=%g <= 10.0*Roundoff=%g\n", ABS(H),(KPP_REAL)10.0*Roundoff);*/
        H = (KPP_REAL)(1.0e-06);
   } /* end if */
   H = SIGN(H, Tdirection);
   Hold = H;
   /*printf("Hold=%f\n", Hold);*/ 
   Reject = 0;
   FirstStep = 1;
   SkipJac = 0;
   SkipLU = 0;
   if ((*T+H*((KPP_REAL)1.0001)-Tend)*Tdirection >= ZERO) {   
	H = Tend - *T;
   } /* end if */
   Nconsecutive = 0;
   RK_ErrorScale(N, ITOL, AbsTol, RelTol, Y, SCAL);
   /*for(i=0; i<NVAR; i++) {
        printf("AbsTol=%g RelTol=%g Y=%g SCAL=%g \n", AbsTol[i], RelTol[i], Y[i], SCAL[i] );
   }*/
   
 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /*~~~> Time loop begins */ 
  /* while Tloop */
Tloop:   while ( (Tend-*T)*Tdirection - Roundoff > ZERO ) {
      /*printf("Starting Tloop: (Tend-*T)*Tdirection - Roundoff=%g > ZERO\n", (Tend-*T)*Tdirection-Roundoff);*/
      /*if ( Reject == 0 ) { */
	FUN_CHEM(*T,Y,FO);
	ISTATUS[Nfun]++;
      /* } * end if */ 
      if ( SkipLU == 0 ) { /* This time around skip the Jac update and LU */
	/*~~~> Compute the Jacobian matrix */
	/*printf("Entering if SkipLU == 0\n");*/
	if ( SkipJac == 0 ) {
	   /*printf("Entering if SkipJac == 0\n");*/
	   JAC_CHEM(*T,Y,FJAC);
	   ISTATUS[Njac]++;
	} /* end if */
	/*~~~> Compute the matrices E1 and E2 and their decompositions */
	RK_Decomp(N,H,FJAC,E1,IP1,E2R,E2I,IP2,&ISING,ISTATUS);
	/*printf("ISING=%d\n", ISING);*/
	if ( ISING != 0 ) {
	   /*printf("Entering if ISING != 0\n");*/
	   ISTATUS[Nsng]++;
	   Nconsecutive++;
	   if (Nconsecutive >= 5) {
		RK_ErrorMsg(-12,*T,H,IERR);
	   }
	   H = H * ((KPP_REAL)0.5);
	   Reject = 1;
	   SkipJac = 1;
	   SkipLU = 0;
	   goto Tloop;
	}
	else {
	   /*printf("Entering if ISING == 0\n");*/
	   Nconsecutive = 0;
	} /* end if */
      } /* end if !SkipLU */

      /*printf("NSTEPS=%d\n", ISTATUS[Nstp]);*/
      ISTATUS[Nstp]++;
      /*printf("NSTEPS=%d\n", ISTATUS[Nstp]);*/
      if (ISTATUS[Nstp] > Max_no_steps) {
	printf("\n Max number of time steps is = %d", Max_no_steps);
	RK_ErrorMsg(-9,*T,H,IERR);
      } /* end if */
      if (((KPP_REAL)0.1)*ABS(H) <= ABS(*T)*Roundoff)
	RK_ErrorMsg(-10,*T,H,IERR);
 
 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 /*~~~>  Loop for the simplified Newton iterations */
 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

      /*~~~> Starting values for Newton iteration */  
      if ( (FirstStep == 1) || (StartNewton == 0) ) {
	/*printf("Entering if FirstStep or Not StartNewton\n");*/
	Set2Zero(N,Z1);
	Set2Zero(N,Z2);
	Set2Zero(N,Z3);
      }
      else
      {
    	/* Evaluate quadratic polynomial */
   	RK_Interpolate("eval",N,H,Hold,Z1,Z2,Z3,CONT);
      } /* end if */

      /*~~~> Initializations for Newton iteration */
      NewtonDone = 0;
      Fac = (KPP_REAL)0.5; /* Step reduction if too many iterations */
  
      /* for NewtonLoop */
      for (NewtonIter = 1; NewtonIter <= NewtonMaxit; NewtonIter++)
      {
	/*printf("Starting NewtonIter loop: NewtonIter=%d\n", NewtonIter);*/
   	/* Prepare the right-hand side */
      	RK_PrepareRHS(N,*T,H,Y,FO,Z1,Z2,Z3,DZ1,DZ2,DZ3);

	/* Solve the linear systems */
  	RK_Solve(N,H,E1,IP1,E2R,E2I,IP2,DZ1,DZ2,DZ3,ISTATUS);

	NewtonIncrement = 
		SQRT((RK_ErrorNorm(N,SCAL,DZ1) * RK_ErrorNorm(N,SCAL,DZ1)
		    + RK_ErrorNorm(N,SCAL,DZ2) * RK_ErrorNorm(N,SCAL,DZ2) 
		    + RK_ErrorNorm(N,SCAL,DZ3) * RK_ErrorNorm(N,SCAL,DZ3))
	            / (KPP_REAL)3.0 );
	/*printf( "NewtonIncrement=%g \n", NewtonIncrement );*/
	if (NewtonIter == 1) {
	   /*printf("Entering NewtonIter=%d == 1\n", NewtonIter);*/
	   Theta = ABS(ThetaMin);
	   NewtonRate = (KPP_REAL)2.0;
	}
	else {
	   /*printf("Entering else from NewtonIter=%d == 1\n", NewtonIter);*/
	   Theta = NewtonIncrement / NewtonIncrementOld;
	   /*printf("Theta=%f\n", Theta);*/
	   if (Theta < (KPP_REAL)0.99) {
		/*printf( "Entering if Theta=%g < 0.99\n", Theta );*/
		NewtonRate = Theta / (ONE-Theta);
	   }
	   else { /* Non-convergence of Newton: Theta too large */
       	      break; /* EXIT NewtonLoop */
	   } /* end if */
	   if (NewtonIter < NewtonMaxit) {
	      /*printf("Entering if NewtonIter=%d < NewtonMaxit=%d\n", NewtonIter, NewtonMaxit);*/
	      /* Predict error at the end of Newton process */
 	      NewtonPredictedErr = NewtonIncrement
                             *pow(Theta,(NewtonMaxit-NewtonIter))/(ONE-Theta);
	      /*printf("NewtonRate=%g NewtonPredictedErr=%g \n", NewtonRate, NewtonPredictedErr);*/
	      if (NewtonPredictedErr >= NewtonTol) {
		/*printf( "Entering if NewtonPredictedErr=%g >= NewtonTol%g \n", NewtonPredictedErr, NewtonTol );*/
	      	/* Non-convergence of Newton: predicted error too large */
	        Qnewton = MIN((KPP_REAL)10.0, NewtonPredictedErr/NewtonTol);
	        Fac = (KPP_REAL)0.8*pow(Qnewton,(-ONE/(1+NewtonMaxit-NewtonIter)));
	        break; /* EXIT NewtonLoop */
	      } /* end if */
	   } /* end if */
   	} /* end if */
	    
   	NewtonIncrementOld = MAX(NewtonIncrement, Roundoff);
	/*printf("NewtonIncrementOld=%f\n", NewtonIncrementOld);*/
  	/*~~~> Update solution */
   	WAXPY(N,-ONE,DZ1,1,Z1,1);	/* Z1 <- Z1 - DZ1 */
   	WAXPY(N,-ONE,DZ2,1,Z2,1);	/* Z2 <- Z2 - DZ2 */
   	WAXPY(N,-ONE,DZ3,1,Z3,1);	/* Z3 <- Z3 - DZ3 */
        /*for(i=0; i<N; i++) 
	   printf("DZ1[%d]=%g Z1[%d]=%g\n", i, DZ1[i], i, Z1[i]);*/
  	/*~~~> Check error in Newton iterations */
	/*printf( "NewtonRate=%g NewtonIncrement=%g NewtonTol=%g \n", NewtonRate, NewtonIncrement, NewtonTol );*/
	NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol);
	/*printf( "NewtonDone=%d \n", NewtonDone );*/
   	if (NewtonDone == 1)
	   break; /* Exit NewtonLoop */
   	if (NewtonIter == NewtonMaxit)
   	{
	   printf("\n Slow or no convergence in Newton Iteration:");
	   printf(" Max no. of Newton iterations reached");
   	} /* end if */
      }/* end for NewtonLoop */

      if ( NewtonDone == 0)
      {
	/*printf( "Entering if NewtonDone == 0\n" );*/
	/*RK_ErrorMsg(-12,*T,H,IERR); */
	H = Fac*H;
	Reject  = 1;
	SkipJac = 1;
	SkipLU  = 0;
	goto Tloop;
      } /* end if */

 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 /*~~~> SDIRK Stage */
 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (SdirkError == 1)
      {
	/*printf("Entering if SdirkError=%d==1\n", SdirkError);*/
      	/*~~~> Starting values for Newton iterations */
   	for (i = 0; i < N; i++)
	{
	   Z4[i] = Z3[i];
 	}
      	/*~~~> Prepare the loop-independent part of the right-hand side */
      	/* G = H*rkBgam(0)*FO + rkTheta(1)*Z1 
		+ rkTheta(2)*Z2 + rkTheta(3)*Z3; */
      	Set2Zero(N,G);
      	if (rkMethod != L3A)
	   WAXPY(N,rkBgam[0]*H, FO,1,G,1);
      	WAXPY(N,rkTheta[0],Z1,1,G,1);  
      	WAXPY(N,rkTheta[1],Z2,1,G,1);  
      	WAXPY(N,rkTheta[2],Z3,1,G,1);  

      	/*~~~> Initializations for Newton iteration */
      	NewtonDone = 0;
      	Fac = (KPP_REAL)0.5; /* Step reduction factor if too many iterations */

      	/* for SDNewtonLoop */
      	for (NewtonIter = 1; NewtonIter <= NewtonMaxit; NewtonIter++)
      	{
	   /*printf("Starting NewtonIter loop: NewtonIter=%d\n", NewtonIter);*/
	   /*~~~> Prepare the loop-dependent part of right-hand side */
	   WADD(N,Y,Z4,TMP);	/* TMP <- Y + Z4 */
	   FUN_CHEM(*T+H,TMP,DZ4);	/* DZ4 <- Fun(Y+Z4) */
	   ISTATUS[Nfun]++;
	   /* DZ4[1,N] = (D[1, N]-Z4[1,N])*(rkGamma/H) + DZ4[1,N]; */
	   WAXPY(N,-ONE*rkGamma/H,Z4,1,DZ4,1);
	   WAXPY(N,rkGamma/H,G,1,DZ4,1);

	   /*~~~> Solve the linear system */
	   KppSolve(E1, DZ4);
	   /*~~~> Note: for a full matrix use Lapack:
		DGETRS('N', 5, 1, E1, N, IP1, DZ4, 5, ISING) */

	   /*~~~> Check convergence of Newton iterations */
	   NewtonIncrement = RK_ErrorNorm(N,SCAL,DZ4);
	   /*printf( "NewtonIncrement=%g \n", NewtonIncrement );*/
	   if (NewtonIter == 1) {
	      /*printf("Entering NewtonIter=%d == 1\n", NewtonIter);*/
	      ThetaSD = ABS(ThetaMin);
              NewtonRate = (KPP_REAL)2.0;
	   }
	   else {
	      /*printf("Entering else from NewtonIter=%d == 1\n", NewtonIter);*/
	      ThetaSD = NewtonIncrement / NewtonIncrementOld;
	      if (ThetaSD < (KPP_REAL)0.99) {
		/*printf( "Entering if Theta=%g < 0.99\n", Theta );*/
		NewtonRate = ThetaSD / (ONE-ThetaSD);
		/* Predict error at the end of Newton process */
		NewtonPredictedErr = NewtonIncrement
			*pow(ThetaSD,(NewtonMaxit-NewtonIter))/(ONE-ThetaSD);
	        /*printf("NewtonRate=%g NewtonPredictedErr=%g \n", NewtonRate, NewtonPredictedErr);*/
		if (NewtonPredictedErr >= NewtonTol) {
		   /*printf( "Entering if NewtonPredictedErr=%g >= NewtonTol%g \n", NewtonPredictedErr, NewtonTol );*/
		   /* Non-convergence of Newton: predicted error too large */
		   /* printf("\n Error too large", NewtonPredictedErr); */
		   Qnewton = MIN((KPP_REAL)10.0,NewtonPredictedErr/NewtonTol);
		   Fac = (KPP_REAL)0.8*pow(Qnewton,(-ONE/(1+NewtonMaxit-NewtonIter)));
		   break; /* EXIT SDNewtonLoop */
		} /* end if */
	      }
	      /* Non-convergence of Newton: Theta too large */
	      else {
		/* prinf("\n Theta too large", ThetaSD); */
		break; /* EXIT SDNewtonLoop */
	      } /* end if */
	   } /* end if */
     	   NewtonIncrementOld = NewtonIncrement;
	   /*printf("NewtonIncrementOld=%f\n", NewtonIncrementOld);*/
	   /* Update solution: Z4 <-- Z4 + DZ4; */
	   WAXPY(N,ONE,DZ4,1,Z4,1);

	   /* Check error in Newton iterations */
	   NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol);
	   /*printf( "NewtonDone=%d \n", NewtonDone );*/
	   if (NewtonDone == 1)
		break;	/* EXIT SDNewtonLoop */
      	} /* end for SDNewtonLoop */
      
      	if ( NewtonDone == 0 ) {
	   /*printf( "Entering if NewtonDone == 0\n" );*/
	   H       = Fac*H;
	   Reject  = 1;
	   SkipJac = 1;
	   SkipLU  = 0;
	   goto Tloop;
      	} /* end if */
      } /* end if */
      /*~~~> End of implified SDIRK Newton iterations */

 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 /*~~~> Error estimation */
 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (SdirkError == 1) {
      	Set2Zero(N, DZ4);
      	if (rkMethod == L3A) {
	   for (i = 0; i < N; i++)
		DZ4[i] = H*rkF[0]*FO[i];
	   if (rkF[0] != ZERO)
		WAXPY(N, rkF[0], Z1, 1, DZ4, 1);
	   if (rkF[1] != ZERO)
		WAXPY(N, rkF[1], Z2, 1, DZ4, 1);
	   if (rkF[2] != ZERO)
		WAXPY(N, rkF[2], Z3, 1, DZ4, 1);
	   for (i = 0; i < N; i++)
		 TMP[i] = Y[i] + Z4[i];
	   FUN_CHEM(*T+H, TMP, DZ1);
	   WAXPY(N, H*rkBgam[4], DZ1, 1, DZ4, 1);
      	}
      	else {
	   /* DZ4(1,N) = rkD(1)*Z1 + rkD(2)*Z2 + rkD(3)*Z3 - Z4; */
	   if (rkD[0] != ZERO)
		WAXPY(N, rkD[0], Z1, 1, DZ4, 1);
	   if (rkD[1] != ZERO)
		WAXPY(N, rkD[1], Z2, 1, DZ4, 1);
	   if (rkD[2] != ZERO)
		WAXPY(N, rkD[2], Z3, 1, DZ4, 1);
	   WAXPY(N, -ONE, Z4, 1, DZ4, 1);
      	} /* end if */
      	Err = RK_ErrorNorm(N,SCAL,DZ4);
      }
      else
      {
      	RK_ErrorEstimate(N,H,*T,Y,FO,E1,IP1,Z1,Z2,Z3,SCAL,&Err,
			FirstStep,Reject,ISTATUS);
      } /* end if */

      /*~~~> Computation of new step size Hnew */
      Fac = pow(Err, (-ONE/rkELO))
	    *MIN(FacSafe,(ONE+2*NewtonMaxit)/(NewtonIter+2*NewtonMaxit));
      Fac = MIN(FacMax,MAX(FacMin,Fac));
      Hnew = Fac*H;

 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 /*~~~> Accept/reject step */
 /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      // accept:
      if (Err < ONE) { /*~~~> STEP IS ACCEPTED */
      	FirstStep = 0;
      	ISTATUS[Nacc]++;
      	if (Gustafsson == 1) {
	   /*~~~> Predictive controller of Gustafsson */
	   if (ISTATUS[Nacc] > 1) {
	      FacGus = FacSafe*(H/Hacc)*pow(Err*Err/ErrOld,(KPP_REAL)(-0.25));
	      FacGus = MIN(FacMax,MAX(FacMin,FacGus));
	      Fac = MIN(Fac,FacGus);	
	      Hnew = Fac*H;
	   } /* end if */
	   Hacc = H;
	   ErrOld = MAX((KPP_REAL)1.0e-02,Err);
      	} /* end if */
        Hold = H;
      	*T = *T + H;
      	/* Update solution: Y <- Y + sum(d_i Z-i) */
      	if (rkD[0] != ZERO)
	   WAXPY(N, rkD[0], Z1, 1, Y, 1);
      	if (rkD[1] != ZERO)
	   WAXPY(N, rkD[1], Z2, 1, Y, 1);
      	if (rkD[2] != ZERO)
	   WAXPY(N, rkD[2], Z3, 1, Y, 1);
      	/* Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3 */
      	if (StartNewton == 1)
	   RK_Interpolate("make",N,H,Hold,Z1,Z2,Z3,CONT);
      	RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL);
	RSTATUS[Ntexit] = *T;
	RSTATUS[Nhnew] = Hnew;
	RSTATUS[Nhacc] = H;
      	Hnew = Tdirection*MIN( MAX(ABS(Hnew),Hmin), Hmax );
      	if (Reject == 1)
	   Hnew = Tdirection*MIN(ABS(Hnew),ABS(H));
      	Reject = 0;
      	if ((*T+Hnew/Qmin-Tend)*Tdirection >= ZERO)
	   H = Tend - *T;
      	else {
	   Hratio = Hnew/H;
	   /* Reuse the LU decomposition */
	   SkipLU = (Theta <= ThetaMin) && (Hratio >= Qmin) && (Hratio <= Qmax);
	   if ( SkipLU == 0)
	   	H = Hnew;
        } /* end if */
        /* If convergence is fast enough, do not update Jacobian */
      	/* SkipJac = (Theta <= ThetaMin); */
      	SkipJac = 0;
      }

      /*~~~> Step is rejected */
      else {
      	if ((FirstStep == 1)  || (Reject == 1))
	   H = FacRej*H;
        else
	   H = Hnew;
        Reject  = 1;
        SkipJac = 1; /* Skip if rejected - Jac is independent of H */
        SkipLU  = 0;
      	if (ISTATUS[Nacc] >= 1)
	   ISTATUS[Nrej]++;
      } /* end if accept */
   } /* while: time Tloop */
 
   /*~~~> Successful exit */
   *IERR = 1;

} /* RK_Integrator */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Handles all error messages and returns IERR = error Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void RK_ErrorMsg(int Code, KPP_REAL T, KPP_REAL H, int* IERR)
{
   Code = *IERR;
   printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
   printf("\nForced exit from RungeKutta due to the following error:\n");

   switch (Code) {
   case -1:
      printf("--> Improper value for maximal no of steps");
      break;
   case -2:
      printf("--> Improper value for maximal no of Newton iterations");
      break;
   case -3:
      printf("--> Hmin/Hmax/Hstart must be positive");
      break;
   case -4:
      printf("--> Improper values for FacMin/FacMax/FacSafe/FacRej");
      break;
   case -5:
      printf("--> Improper value for ThetaMin");
      break;
   case -6:
      printf("--> Newton stopping tolerance too small");
      break;
   case -7:
      printf("--> Improper values for Qmin, Qmax");
      break;
   case -8:
      printf("--> Tolerances are too small");
      break;
   case -9:
      printf("--> No of steps exceeds maximum bound");
      break;
   case -10:
      printf("--> Step size too small: (T + 10*H = T) or H < Roundoff");
      break;
   case -11:
      printf("--> Matrix is repeatedly singular");
      break;
   case -12:
      printf("--> Non-convergence of Newton iterations");
      break;
   case -13:
      printf("--> Requested RK method not implemented");
      break;
   default:
      printf("Unknown Error code: %d \n", Code);
   } /* end switch */ 

   printf("\n     T=%e,  H =%e", T, H);
   printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

} /* RK_ErrorMsg */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*	Handles all error messages and returns SCAL */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void RK_ErrorScale(int N, 
   /*~~~> Input arguments: */
	int ITOL, KPP_REAL AbsTol[], KPP_REAL RelTol[], KPP_REAL Y[],
   /*~~~> Output arguments: */
	KPP_REAL SCAL[])
{
   /*printf("Starting RK_ErrorScale\n");*/
   int i;
   if (ITOL == 0) {
	for (i = 0; i < N; i++) {
	   SCAL[i] = ONE/(AbsTol[0]+RelTol[0]*ABS(Y[i]));
	} /* end for loop*/
   }
   else {
	for (i = 0; i< N; i++) {
	   SCAL[i] = ONE/(AbsTol[i]+RelTol[i]*ABS(Y[i]));
	} /* end for loop */
   } /* end if */
   /*printf("Ending RK_ErrorScale\n");*/
} /* RK_ErrorScale */  

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void RK_Transform(int N, KPP_REAL Tr[][3], KPP_REAL Z1[], KPP_REAL Z2[],
		KPP_REAL Z3[], KPP_REAL W1[], KPP_REAL W2[], KPP_REAL W3[])
-->	W <-- Tr x Z 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
   int i;
   KPP_REAL x1,x2,x3;
   for (i = 0; i < N; i++)
   {
	x1 = Z1[i];
	x2 = Z2[i]; 
	x3 = Z3[i]; 
	W1[i] = Tr[0][0]*x1 + Tr[0][1]*x2 + Tr[0][2]*x3;
	W2[i] = Tr[1][0]*x1 + Tr[1][1]*x2 + Tr[1][2]*x3;
	W1[i] = Tr[2][0]*x1 + Tr[2][1]*x2 + Tr[2][2]*x3;
   }  end for loop
}  end RK_Transform */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*--> Constructs or evaluates a quadratic polynomial that
      interpolates the Z solution at current step and
      provides starting values for the next step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void RK_Interpolate(char action[], int N, KPP_REAL H, KPP_REAL Hold, KPP_REAL Z1[],
	KPP_REAL Z2[], KPP_REAL Z3[], KPP_REAL CONT[][RKmax])
{
   int i;
   KPP_REAL r,x1,x2,x3,den;
 
   /* Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3 */
   if (action == "make") {
	den = (rkC[2]-rkC[1])*(rkC[1]-rkC[0])*(rkC[0]-rkC[2]);
	for (i = 0; i < N; i++)
	{
	   CONT[i][0] = (-rkC[2]*rkC[2]*rkC[1]*Z1[i]
			+Z3[i]*rkC[1]*rkC[0]*rkC[0]
	   		+rkC[1]*rkC[1]*rkC[2]*Z1[i]
	   		-rkC[1]*rkC[1]*rkC[0]*Z3[i]
	   		+rkC[2]*rkC[2]*rkC[0]*Z2[i]
	   		-Z2[i]*rkC[2]*rkC[0]*rkC[0])/den-Z3[i];
	   CONT[i][1] = -( rkC[0]*rkC[0]*(Z3[i]-Z2[i])
	   		+ rkC[1]*rkC[1]*(Z1[i]-Z3[i])
	   		+ rkC[2]*rkC[2]*(Z2[i]-Z1[i]) )/den;
	   CONT[i][2] = ( rkC[0]*(Z3[i]-Z2[i])
	   		+ rkC[1]*(Z1[i]-Z3[i])
	   		+ rkC[2]*(Z2[i]-Z1[i]) )/den;
	} /* end for loop */
   }
   /* Evaluate quadratic polynomial */
   else if (action == "eval") {
   	r  = H / Hold;
	x1 = ONE + rkC[0]*r;
	x2 = ONE + rkC[1]*r;
	x3 = ONE + rkC[2]*r;
	for (i = 0; i < N; i++)
	{
	   Z1[i] = CONT[i][0]+x1*(CONT[i][1]+x1*CONT[i][2]);	
	   Z2[i] = CONT[i][0]+x2*(CONT[i][1]+x2*CONT[i][2]);	
	   Z3[i] = CONT[i][0]+x3*(CONT[i][1]+x3*CONT[i][2]);
	} /* end for loop */
   } /* end if */
} /* RK_Interpolate */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~> Prepare the right-hand side for Newton iterations
       R = Z - hA x F
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void RK_PrepareRHS(int N, KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL FO[],
	KPP_REAL Z1[], KPP_REAL Z2[], KPP_REAL Z3[], KPP_REAL R1[], KPP_REAL R2[],
	KPP_REAL R3[])
{
   KPP_REAL TMP[N], F[N];
   WCOPY(N,Z1,1,R1,1); /* R1 <- Z1 */	
   WCOPY(N,Z2,1,R2,1); /* R2 <- Z2 */	
   WCOPY(N,Z3,1,R3,1); /* R3 <- Z3 */	
 
   if (rkMethod == L3A)
   {
	WAXPY(N,-H*rkA[0][0],FO,1,R1,1); /* R1 <- R1 - h*A_10*FO */
	WAXPY(N,-H*rkA[1][0],FO,1,R2,1); /* R2 <- R2 - h*A_20*FO */
	WAXPY(N,-H*rkA[2][0],FO,1,R3,1); /* R3 <- R3 - h*A_30*FO */
   } /* end if */

   WADD(N,Y,Z1,TMP);			/* TMP <- Y + Z1 */
   FUN_CHEM(T+rkC[0]*H,TMP,F);		/* F1 <- Fun(Y+Z1) */
   WAXPY(N,-H*rkA[0][0],F,1,R1,1); 	/* R1 <- R1 - h*A_11*F1 */
   WAXPY(N,-H*rkA[1][0],F,1,R2,1); 	/* R2 <- R2 - h*A_21*F1 */
   WAXPY(N,-H*rkA[2][0],F,1,R3,1); 	/* R3 <- R3 - h*A_31*F1 */

   WADD(N,Y,Z2,TMP);			/* TMP <- Y + Z2 */
   FUN_CHEM(T+rkC[1]*H,TMP,F);		/* F2 <- Fun(Y+Z2) */
   WAXPY(N,-H*rkA[0][1],F,1,R1,1); 	/* R1 <- R1 - h*A_12*F2 */
   WAXPY(N,-H*rkA[1][1],F,1,R2,1); 	/* R2 <- R2 - h*A_22*F2 */
   WAXPY(N,-H*rkA[2][1],F,1,R3,1); 	/* R3 <- R3 - h*A_32*F2 */

   WADD(N,Y,Z3,TMP);			/* TMP <- Y + Z3 */
   FUN_CHEM(T+rkC[2]*H,TMP,F);		/* F3 <- Fun(Y+Z3) */
   WAXPY(N,-H*rkA[0][2],F,1,R1,1); 	/* R1 <- R1 - h*A_13*F3 */
   WAXPY(N,-H*rkA[1][2],F,1,R2,1); 	/* R2 <- R2 - h*A_23*F3 */
   WAXPY(N,-H*rkA[2][2],F,1,R3,1); 	/* R3 <- R3 - h*A_33*F3 */
} /* RK_PrepareRHS */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~> Compute the matrices E1 and E2 and their decompositions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void RK_Decomp(int N, KPP_REAL H, KPP_REAL FJAC[], KPP_REAL E1[],
	int IP1[], KPP_REAL E2R[], KPP_REAL E2I[], int IP2[], 
	int* ISING, int ISTATUS[]) 
{
   /*printf("Starting RK_Decomp\n");*/
   KPP_REAL Alpha, Beta, Gamma;
   int i,j;

   Gamma = rkGamma / H;
   Alpha = rkAlpha / H;
   Beta  = rkBeta / H;

   for (i = 0; i < LU_NONZERO; i++) {
	E1[i] = -FJAC[i];
   }
   for (i = 0; i < NVAR; i++) {
	j = LU_DIAG[i];
	E1[j] = E1[j] + Gamma;
   }
   *ISING = KppDecomp(E1);
   /*~~~> Note: for a full matrix use Lapack:
   for (j = 0; j < N; j++) {
	for (i = 0; i < N; i++)
	{  
	   E1[i,j] = -FJAC[i][j];
	}
	E1[i][j] = E1[i][j]+Gamma;
   }
   DGETRF(N, N, E1, N, IP1, ISING); */

   if (*ISING != 0) {
	/*printf("Matrix is singular, ISING=%d\n", *ISING);*/
	ISTATUS[Ndec]++;
	return;
   }

   for (i = 0; i < LU_NONZERO; i++) {
	E2R[i] = (KPP_REAL)(-FJAC[i]);
	E2I[i] = ZERO;
   }
   for (i = 0; i < NVAR; i++) {
	j = LU_DIAG[i];
	E2R[j] = E2R[j] + Alpha;
	E2I[j] = E2I[j] + Beta;
   }
   *ISING = KppDecompCmplxR(E2R, E2I);
   /*printf("Matrix is singular, ISING=%d\n", *ISING);*/
   /*~~~> Note: for a full matrix use Lapack:
   for (j = 0; j < N; j++) {
      	for (i = 0; i < N; i++) {
	   E2[i,j] = DCMPLX( -FJAC[i,j], ZERO );
	}
	E2[j,j] = E2[j,j] + CMPLX( Alpha, Beta );
   }
   ZGETRF(N,N,E2,N,IP2,ISING); */
   ISTATUS[Ndec]++;
} /*RK_Decomp */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/ 
void RK_Solve(int N, KPP_REAL H, KPP_REAL E1[], int IP1[], KPP_REAL E2R[],
	KPP_REAL E2I[], int IP2[],KPP_REAL R1[], KPP_REAL R2[], KPP_REAL R3[],
	int ISTATUS[])
{
   KPP_REAL x1, x2, x3;
   KPP_REAL BCR[N], BCI[N];
   int i;

   /* Z <- h^{-1) T^{-1) A^{-1) x Z */
   for (i = 0; i < N; i++)
   {
	x1 = R1[i]/H;
	x2 = R2[i]/H;
	x3 = R3[i]/H;
        R1[i] = rkTinvAinv[0][0]*x1 + rkTinvAinv[0][1]*x2 + rkTinvAinv[0][2]*x3;
	R2[i] = rkTinvAinv[1][0]*x1 + rkTinvAinv[1][1]*x2 + rkTinvAinv[1][2]*x3;
	R3[i] = rkTinvAinv[2][0]*x1 + rkTinvAinv[2][1]*x2 + rkTinvAinv[2][2]*x3;
   } 
   KppSolve(E1,R1);
   /*~~~> Note: for a full matrix use Lapack:
   DGETRS('N',5,1,E1,N,IP1,R1,5,0); */
   for(i = 0; i < N; i++)
   {
        BCR[i] = (KPP_REAL)(R2[i]);
        BCI[i] = (KPP_REAL)(R3[i]);
   }
   KppSolveCmplxR(E2R,E2I,BCR,BCI);
   /*~~~> Note: for a full matrix use Lapack:
   ZGETRS ('N',N,1,E2,N,IP2,BC,N,0); */
   for (i = 0; i < N; i++)
   {
	R2[i] = (KPP_REAL)( BCR[i] );
        R3[i] = (KPP_REAL)( BCI[i] );
   }

   /* Z <- T x Z */
   for (i = 0; i < N; i++)
   {
	x1 = R1[i];
	x2 = R2[i];
	x3 = R3[i];
        R1[i] = rkT[0][0]*x1 + rkT[0][1]*x2 + rkT[0][2]*x3;
	R2[i] = rkT[1][0]*x1 + rkT[1][1]*x2 + rkT[1][2]*x3;
	R3[i] = rkT[2][0]*x1 + rkT[2][1]*x2 + rkT[2][2]*x3;
   } 
   ISTATUS[Nsol]++;
} /* RK_Solve */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void RK_ErrorEstimate(int N, KPP_REAL H, KPP_REAL T, KPP_REAL Y[], 
	KPP_REAL FO[], KPP_REAL E1[], int IP1[], KPP_REAL Z1[], 
	KPP_REAL Z2[], KPP_REAL Z3[], KPP_REAL SCAL[], KPP_REAL* Err,
	int FirstStep, int Reject, int ISTATUS[])
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
   KPP_REAL F1[N],F2[N],TMP[N];
   KPP_REAL HrkE1,HrkE2,HrkE3;
   int i;

   HrkE1  = rkE[0]/H;
   HrkE2  = rkE[1]/H;
   HrkE3  = rkE[2]/H;

   for (i = 0; i < N; i++)
   {
	F2[i]  = HrkE1*Z1[i]+HrkE2*Z2[i]+HrkE3*Z3[i];
	TMP[i] = rkE[0]*FO[i] + F2[i];
   }
   
   KppSolve(E1, TMP);
   if ((rkMethod == R1A) || (rkMethod == GAU) || (rkMethod == L3A))
	KppSolve(E1,TMP);
   if (rkMethod == GAU)
	KppSolve(E1,TMP);
   /*~~~> Note: for a full matrix use Lapack:
	DGETRS ('N',N,1,E1,N,IP1,TMP,N,0);
      	if ((rkMethod == R1A) || (rkMethod == GAU) || (rkMethod == L3A))
	   DGETRS('N',N,1,E1,N,IP1,TMP,N,0);
   	if (rkMethod == GAU)
	   DGETRS('N',N,1,E1,N,IP1,TMP,N,0); */

   *Err = RK_ErrorNorm(N,SCAL,TMP);

   if (*Err < ONE)
	return;
   /* firej */
   if ((FirstStep == 1) || (Reject == 1)) {
	for (i = 0; i < N; i++) {
	   TMP[i] = Y[i] + TMP[i];
        } /* end for */
        FUN_CHEM(T,TMP,F1);
        ISTATUS[Nfun]++;
	for (i = 0; i < N; i++) {
	   TMP[i] = F1[i] + F2[i];
	} /* end for */

        KppSolve(E1, TMP);
	/*~~~> Note: for a full matrix use Lapack:
	   DGETRS ('N',N,1,E1,N,IP1,TMP,N,0); */
        *Err = RK_ErrorNorm(N,SCAL,TMP);
   } /* end if firej */
} /* RK_ErrorEstimate */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   KPP_REAL RK_ErrorNorm(int N, KPP_REAL SCAL[], KPP_REAL DY[])
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
   int i;
   KPP_REAL RK_ErrorNorm = ZERO;
   for (i = 0; i < N; i++) {
	RK_ErrorNorm = RK_ErrorNorm + (DY[i]*SCAL[i]) * (DY[i]*SCAL[i]);
   }    
   RK_ErrorNorm = MAX( SQRT(RK_ErrorNorm/N), (KPP_REAL)1.0e-10 );
   return RK_ErrorNorm;
} /* RK_ErrorNorm */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	The coefficients of the 3-stage Radau-2A method
	(given to ~30 accurate digits)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Radau2A_Coefficients()
{
   /*printf("Starting Radau2A\n");*/
   /* The coefficients of the Radau2A method */
   KPP_REAL b0;
   /* b0 = (KPP_REAL)1.0; */
   if (SdirkError == 1) {
   	b0 = (KPP_REAL)(0.2e-01);
	/*printf("b0 = %f\n", b0);*/
   }
   else {
	b0 = (KPP_REAL)(0.5e-01);
	/*printf("b0 = %f\n", b0);*/
   } /*end if */
   /* The coefficients of the Radau2A method */
   rkMethod = R2A;

   rkA[0][0] = (KPP_REAL)(1.968154772236604258683861429918299e-01);
   rkA[0][1] = (KPP_REAL)(-6.55354258501983881085227825696087e-02);
   rkA[0][2] = (KPP_REAL)(2.377097434822015242040823210718965e-02);
   rkA[1][0] = (KPP_REAL)(3.944243147390872769974116714584975e-01);
   rkA[1][1] = (KPP_REAL)(2.920734116652284630205027458970589e-01);
   rkA[1][2] = (KPP_REAL)(-4.154875212599793019818600988496743e-02);
   rkA[2][0] = (KPP_REAL)(3.764030627004672750500754423692808e-01);
   rkA[2][1] = (KPP_REAL)(5.124858261884216138388134465196080e-01);
   rkA[2][2] = (KPP_REAL)(1.111111111111111111111111111111111e-01);

   rkB[0] = (KPP_REAL)(3.764030627004672750500754423692808e-01);
   rkB[1] = (KPP_REAL)(5.124858261884216138388134465196080e-01);
   rkB[2] = (KPP_REAL)(1111111111111111111111111111111111e-01);

   rkC[0] = (KPP_REAL)(1.550510257216821901802715925294109e-01);
   rkC[1] = (KPP_REAL)(6.449489742783178098197284074705891e-01);
   rkC[2] = ONE;

   /* New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j */
   rkD[0] = ZERO;
   rkD[1] = ZERO;
   rkD[2] = ONE;
   
   /* Classical error estimator: */
   /* H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */
   rkE[0] = ONE*b0;
   rkE[1] = (KPP_REAL)(-10.04880939982741556246032950764708)*b0;
   rkE[2] = (KPP_REAL)1.382142733160748895793662840980412*b0;
   rkE[3] = (KPP_REAL)(-.3333333333333333333333333333333333)*b0;

   /* Sdirk error estimator */
   rkBgam[0] = b0;
   rkBgam[1] = (KPP_REAL).3764030627004672750500754423692807
		-(KPP_REAL)1.558078204724922382431975370686279*b0;
   rkBgam[2] = (KPP_REAL).8914115380582557157653087040196118*b0
		+(KPP_REAL).5124858261884216138388134465196077;
   rkBgam[3] = (KPP_REAL)(-.1637777184845662566367174924883037)
		-(KPP_REAL).3333333333333333333333333333333333*b0;
   rkBgam[4] = (KPP_REAL).2748888295956773677478286035994148;

   /* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta[0] = (KPP_REAL)(-1.520677486405081647234271944611547)
		-(KPP_REAL)10.04880939982741556246032950764708*b0;
   rkTheta[1] = (KPP_REAL)2.070455145596436382729929151810376
		+(KPP_REAL)1.382142733160748895793662840980413*b0;
   rkTheta[2] = (KPP_REAL)(-.3333333333333333333333333333333333)*b0
		-(KPP_REAL).3744441479783868387391430179970741;

   /* Local order of error estimator */
   if ( b0 == ZERO)
   	rkELO  = (KPP_REAL)6.0;
   else
        rkELO  = (KPP_REAL)4.0;

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ~~~> Diagonalize the RK matrix:
    	 rkTinv * inv(rkA) * rkT =
    	           |  rkGamma      0           0     |
      	           |      0      rkAlpha   -rkBeta   |
      	           |      0      rkBeta     rkAlpha  |
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = (KPP_REAL)3.637834252744495732208418513577775;
   rkAlpha = (KPP_REAL)2.681082873627752133895790743211112;
   rkBeta  = (KPP_REAL)3.050430199247410569426377624787569;

   rkT[0][0] = (KPP_REAL)(9.443876248897524148749007950641664e-02);
   rkT[0][1] = (KPP_REAL)(-1.412552950209542084279903838077973e-01);
   rkT[0][2] = (KPP_REAL)(-3.00291941051474244918611170890539e-02);
   rkT[1][0] = (KPP_REAL)(2.502131229653333113765090675125018e-01);
   rkT[1][1] = (KPP_REAL)(2.041293522937999319959908102983381e-01);
   rkT[1][2] = (KPP_REAL)(3.829421127572619377954382335998733e-01);
   rkT[2][0] =  ONE;
   rkT[2][1] =  ONE;
   rkT[2][2] =  ZERO;

   rkTinv[0][0] = (KPP_REAL)4.178718591551904727346462658512057;
   rkTinv[0][1] = (KPP_REAL)(3.27682820761062387082533272429617e-01);
   rkTinv[0][2] = (KPP_REAL)(5.233764454994495480399309159089876e-01);
   rkTinv[1][0] = (KPP_REAL)(-4.178718591551904727346462658512057);
   rkTinv[1][1] = (KPP_REAL)(-3.27682820761062387082533272429617e-01);
   rkTinv[1][2] = (KPP_REAL)(4.766235545005504519600690840910124e-01);
   rkTinv[2][0] = (KPP_REAL)(-5.02872634945786875951247343139544e-01);
   rkTinv[2][1] = (KPP_REAL)2.571926949855605429186785353601676;
   rkTinv[2][2] = (KPP_REAL)(-5.960392048282249249688219110993024e-01);

   rkTinvAinv[0][0] = (KPP_REAL)(1.520148562492775501049204957366528e+01);
   rkTinvAinv[0][1] = (KPP_REAL)1.192055789400527921212348994770778;
   rkTinvAinv[0][2] = (KPP_REAL)1.903956760517560343018332287285119;
   rkTinvAinv[1][0] = (KPP_REAL)(-9.669512977505946748632625374449567);
   rkTinvAinv[1][1] = (KPP_REAL)(-8.724028436822336183071773193986487);
   rkTinvAinv[1][2] = (KPP_REAL)3.096043239482439656981667712714881;
   rkTinvAinv[2][0] = (KPP_REAL)(-1.409513259499574544876303981551774e+01);
   rkTinvAinv[2][1] = (KPP_REAL)5.895975725255405108079130152868952;
   rkTinvAinv[2][2] = (KPP_REAL)(-1.441236197545344702389881889085515e-01);

   rkAinvT[0][0] = (KPP_REAL).3435525649691961614912493915818282;
   rkAinvT[0][1] = (KPP_REAL)(-.4703191128473198422370558694426832);
   rkAinvT[0][2] = (KPP_REAL).3503786597113668965366406634269080;
   rkAinvT[1][0] = (KPP_REAL).9102338692094599309122768354288852;
   rkAinvT[1][1] = (KPP_REAL)1.715425895757991796035292755937326;
   rkAinvT[1][2] = (KPP_REAL).4040171993145015239277111187301784;
   rkAinvT[2][0] = (KPP_REAL)3.637834252744495732208418513577775;
   rkAinvT[2][1] = (KPP_REAL)2.681082873627752133895790743211112;
   rkAinvT[2][2] = (KPP_REAL)(-3.050430199247410569426377624787569);
   /*printf("Ending Radau2A");*/
} /* Radau2A_Coefficients */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	The coefficients of the 3-stage Lobatto-3C method
	(given to ~30 accurate digits)
	The parameter b0 can be chosen to tune the error estimator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Lobatto3C_Coefficients ()
{
   KPP_REAL b0;	
   rkMethod = L3C;
   /* b0 = 1.0d0 */
   if (SdirkError == 1)
	b0 = (KPP_REAL)0.2;
   else
	b0 = (KPP_REAL)0.5;
    
   /* The coefficients of the Lobatto3C method */
   rkA[0][0] = (KPP_REAL).1666666666666666666666666666666667;
   rkA[0][1] = (KPP_REAL)(-.3333333333333333333333333333333333);
   rkA[0][2] = (KPP_REAL).1666666666666666666666666666666667;
   rkA[1][0] = (KPP_REAL).1666666666666666666666666666666667;
   rkA[1][1] = (KPP_REAL).4166666666666666666666666666666667;
   rkA[1][2] = (KPP_REAL)(-.8333333333333333333333333333333333e-01);
   rkA[2][0] = (KPP_REAL).1666666666666666666666666666666667;
   rkA[2][1] = (KPP_REAL).6666666666666666666666666666666667;
   rkA[2][2] = (KPP_REAL).1666666666666666666666666666666667;

   rkB[0] = (KPP_REAL).1666666666666666666666666666666667;
   rkB[1] = (KPP_REAL).6666666666666666666666666666666667;
   rkB[2] = (KPP_REAL).1666666666666666666666666666666667;

   rkC[0] = ZERO;
   rkC[1] = (KPP_REAL)0.5;
   rkC[2] = ONE;

   /* Classical error estimator, embedded solution: */
   rkBhat[0] = b0;
   rkBhat[1] = (KPP_REAL).16666666666666666666666666666666667-b0;
   rkBhat[2] = (KPP_REAL).66666666666666666666666666666666667;
   rkBhat[3] = (KPP_REAL).16666666666666666666666666666666667;

   /* New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j */
   rkD[0] = ZERO;
   rkD[1] = ZERO;
   rkD[2] = ONE;

   /* Classical error estimator: */
   /* H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */
   rkE[0] = (KPP_REAL).3808338772072650364017425226487022*b0;
   rkE[1] = (KPP_REAL)(-1.142501631621795109205227567946107)*b0;
   rkE[2] = (KPP_REAL)(-1.523335508829060145606970090594809)*b0;
   rkE[3] = (KPP_REAL).3808338772072650364017425226487022*b0;

   /* Sdirk error estimator */
   rkBgam[0] = b0;
   rkBgam[1] = (KPP_REAL).1666666666666666666666666666666667-1.0*b0;
   rkBgam[2] = (KPP_REAL).6666666666666666666666666666666667;
   rkBgam[3] = (KPP_REAL)(-.2141672105405983697350758559820354);
   rkBgam[4] = (KPP_REAL).3808338772072650364017425226487021;

   /* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta[0] = -3.0*b0-(KPP_REAL).3808338772072650364017425226487021;
   rkTheta[1] = -4.0*b0+(KPP_REAL)1.523335508829060145606970090594808;
   rkTheta[2] = (KPP_REAL)(-.142501631621795109205227567946106)+b0;

   /* Local order of error estimator */
   if (b0 == ZERO)
	rkELO  = 5.0;
   else
	rkELO  = 4.0;

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ~~~> Diagonalize the RK matrix:
    rkTinv * inv(rkA) * rkT =
              |  rkGamma      0           0     |
              |      0      rkAlpha   -rkBeta   |
              |      0      rkBeta     rkAlpha  |
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = (KPP_REAL)2.625816818958466716011888933765284;
   rkAlpha = (KPP_REAL)1.687091590520766641994055533117359;
   rkBeta  = (KPP_REAL)2.508731754924880510838743672432351;

   rkT[0][1] = ONE;
   rkT[0][1] = ONE;
   rkT[0][2] = ZERO;
   rkT[1][0] = (KPP_REAL).4554100411010284672111720348287483;
   rkT[1][1] = (KPP_REAL)(-.6027050205505142336055860174143743);
   rkT[1][2] = (KPP_REAL)(-.4309321229203225731070721341350346);
   rkT[2][0] = (KPP_REAL)2.195823345445647152832799205549709;
   rkT[2][1] = (KPP_REAL)(-1.097911672722823576416399602774855);
   rkT[2][2] = (KPP_REAL).7850032632435902184104551358922130;

   rkTinv[0][0] = (KPP_REAL).4205559181381766909344950150991349;
   rkTinv[0][1] = (KPP_REAL).3488903392193734304046467270632057;
   rkTinv[0][2] = (KPP_REAL).1915253879645878102698098373933487;
   rkTinv[1][0] = (KPP_REAL).5794440818618233090655049849008650;
   rkTinv[1][1] = (KPP_REAL)(-.3488903392193734304046467270632057);
   rkTinv[1][2] = (KPP_REAL)(-.1915253879645878102698098373933487);
   rkTinv[2][0] = (KPP_REAL)(-.3659705575742745254721332009249516);
   rkTinv[2][1] = (KPP_REAL)(-1.463882230297098101888532803699806);
   rkTinv[2][2] = (KPP_REAL).4702733607340189781407813565524989;

   rkTinvAinv[0][0] = (KPP_REAL)1.104302803159744452668648155627548;
   rkTinvAinv[0][1] = (KPP_REAL).916122120694355522658740710823143;
   rkTinvAinv[0][2] = (KPP_REAL).5029105849749601702795812241441172;
   rkTinvAinv[1][0] = (KPP_REAL)1.895697196840255547331351844372453;
   rkTinvAinv[1][1] = (KPP_REAL)3.083877879305644477341259289176857;
   rkTinvAinv[1][2] = (KPP_REAL)(-1.502910584974960170279581224144117);
   rkTinvAinv[2][0] = (KPP_REAL).8362439183082935036129145574774502;
   rkTinvAinv[2][1] = (KPP_REAL)(-3.344975673233174014451658229909802);
   rkTinvAinv[2][2] = (KPP_REAL).312908409479233358005944466882642;

   rkAinvT[0][0] = (KPP_REAL)2.625816818958466716011888933765282;
   rkAinvT[0][1] = (KPP_REAL)1.687091590520766641994055533117358;
   rkAinvT[0][2] = (KPP_REAL)(-2.508731754924880510838743672432351);
   rkAinvT[1][0] = (KPP_REAL)1.195823345445647152832799205549710;
   rkAinvT[1][1] = (KPP_REAL)(-2.097911672722823576416399602774855);
   rkAinvT[1][2] = (KPP_REAL).7850032632435902184104551358922130;
   rkAinvT[2][0] = (KPP_REAL)5.765829871932827589653709477334136;
   rkAinvT[2][1] = (KPP_REAL).1170850640335862051731452613329320;
   rkAinvT[2][2] = (KPP_REAL)4.078738281412060947659653944216779;
} /* Lobatto3C_Coefficients */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	The coefficients of the 3-stage Gauss method
	(given to ~30 accurate digits)
	The parameter b3 can be chosen by the user
	to tune the error estimator
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Gauss_Coefficients()
{
   rkMethod = GAU;
   
   /* The coefficients of the Gauss method */
   KPP_REAL b0;
   /*b0 = (KPP_REAL)4.0; */
   b0 = (KPP_REAL)0.1;

   /* The coefficients of the Gauss method */
   rkA[0][0] = (KPP_REAL).1388888888888888888888888888888889;
   rkA[0][1] = (KPP_REAL)(-.359766675249389034563954710966045e-01);
   rkA[0][2] = (KPP_REAL)(.97894440153083260495800422294756e-02);
   rkA[1][0] = (KPP_REAL).3002631949808645924380249472131556;
   rkA[1][1] = (KPP_REAL).2222222222222222222222222222222222;
   rkA[1][2] = (KPP_REAL)(-.224854172030868146602471694353778e-01);
   rkA[2][0] = (KPP_REAL).2679883337624694517281977355483022;
   rkA[2][1] = (KPP_REAL).4804211119693833479008399155410489;
   rkA[2][2] = (KPP_REAL).1388888888888888888888888888888889;

   rkB[0] = (KPP_REAL).2777777777777777777777777777777778;
   rkB[1] = (KPP_REAL).4444444444444444444444444444444444;
   rkB[2] = (KPP_REAL).2777777777777777777777777777777778;

   rkC[0] = (KPP_REAL).1127016653792583114820734600217600;
   rkC[1] = (KPP_REAL).5000000000000000000000000000000000;
   rkC[2] = (KPP_REAL).8872983346207416885179265399782400;

   /* Classical error estimator, embedded solution: */
   rkBhat[0] = b0;
   rkBhat[1] = (KPP_REAL)(-1.4788305577012361475298775666303999)*b0
   		 +(KPP_REAL).27777777777777777777777777777777778;
   rkBhat[2] = (KPP_REAL).44444444444444444444444444444444444
    	         +(KPP_REAL).66666666666666666666666666666666667*b0;
   rkBhat[3] = (KPP_REAL)(-.18783610896543051913678910003626672)*b0
  	         +(KPP_REAL).27777777777777777777777777777777778;

   /* New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j */
   rkD[0] = (KPP_REAL)(.1666666666666666666666666666666667e+01);
   rkD[1] = (KPP_REAL)(-.1333333333333333333333333333333333e+01);
   rkD[2] = (KPP_REAL)(.1666666666666666666666666666666667e+01);

   /* Classical error estimator: */
   /* H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */
   rkE[0] = (KPP_REAL).2153144231161121782447335303806954*b0;
   rkE[1] = (KPP_REAL)(-2.825278112319014084275808340593191)*b0;
   rkE[2] = (KPP_REAL).2870858974881495709929780405075939*b0;
   rkE[3] = (KPP_REAL)(-.4558086256248162565397206448274867e-01)*b0;

   /* Sdirk error estimator */
   rkBgam[0] = ZERO;
   rkBgam[1] = (KPP_REAL).2373339543355109188382583162660537;
   rkBgam[2] = (KPP_REAL).5879873931885192299409334646982414;
   rkBgam[3] = (KPP_REAL)(-.4063577064014232702392531134499046e-01);
   rkBgam[4] = (KPP_REAL).2153144231161121782447335303806955;

   /* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta[0] = (KPP_REAL)(-2.594040933093095272574031876464493);
   rkTheta[1] = (KPP_REAL)1.824611539036311947589425112250199;
   rkTheta[2] = (KPP_REAL).1856563166634371860478043996459493;

   /* ELO = local order of classical error estimator */
   rkELO = (KPP_REAL)4.0;

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ~~~> Diagonalize the RK matrix:
    rkTinv * inv(rkA) * rkT =
              |  rkGamma      0           0     |
              |      0      rkAlpha   -rkBeta   |
              |      0      rkBeta     rkAlpha  |
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = (KPP_REAL)4.644370709252171185822941421408064;
   rkAlpha = (KPP_REAL)3.677814645373914407088529289295970;
   rkBeta  = (KPP_REAL)3.508761919567443321903661209182446;

   rkT[0][0] = (KPP_REAL)(.7215185205520017032081769924397664e-01);
   rkT[0][1] = (KPP_REAL)(-.8224123057363067064866206597516454e-01);
   rkT[0][2] = (KPP_REAL)(-.6012073861930850173085948921439054e-01);
   rkT[1][0] = (KPP_REAL).1188325787412778070708888193730294;
   rkT[1][1] = (KPP_REAL)(.5306509074206139504614411373957448e-01);
   rkT[1][2] = (KPP_REAL).3162050511322915732224862926182701;
   rkT[2][0] = ONE;
   rkT[2][1] = ONE;
   rkT[2][2] = ZERO;

   rkTinv[0][0] = (KPP_REAL)5.991698084937800775649580743981285;
   rkTinv[0][1] = (KPP_REAL)1.139214295155735444567002236934009;
   rkTinv[0][2] = (KPP_REAL).4323121137838583855696375901180497;
   rkTinv[1][0] = (KPP_REAL)(-5.991698084937800775649580743981285);
   rkTinv[1][1] = (KPP_REAL)(-1.139214295155735444567002236934009);
   rkTinv[1][2] = (KPP_REAL).5676878862161416144303624098819503;
   rkTinv[2][0] = (KPP_REAL)(-1.246213273586231410815571640493082);
   rkTinv[2][1] = (KPP_REAL)2.925559646192313662599230367054972;
   rkTinv[2][2] = (KPP_REAL)(-.2577352012734324923468722836888244);

   rkTinvAinv[0][0] = (KPP_REAL)27.82766708436744962047620566703329;
   rkTinvAinv[0][1] = (KPP_REAL)5.290933503982655311815946575100597;
   rkTinvAinv[0][2] = (KPP_REAL)2.007817718512643701322151051660114;
   rkTinvAinv[1][0] = (KPP_REAL)(-17.66368928942422710690385180065675);
   rkTinvAinv[1][1] = (KPP_REAL)(-14.45491129892587782538830044147713);
   rkTinvAinv[1][2] = (KPP_REAL)2.992182281487356298677848948339886;
   rkTinvAinv[2][0] = (KPP_REAL)(-25.60678350282974256072419392007303);
   rkTinvAinv[2][1] = (KPP_REAL)6.762434375611708328910623303779923;
   rkTinvAinv[2][2] = (KPP_REAL)1.043979339483109825041215970036771;

   rkAinvT[0][0] = (KPP_REAL).3350999483034677402618981153470483;
   rkAinvT[0][1] = (KPP_REAL)(-.5134173605009692329246186488441294);
   rkAinvT[0][2] = (KPP_REAL)(.6745196507033116204327635673208923e-01);
   rkAinvT[1][0] = (KPP_REAL).5519025480108928886873752035738885;
   rkAinvT[1][1] = (KPP_REAL)1.304651810077110066076640761092008;
   rkAinvT[1][2] = (KPP_REAL).9767507983414134987545585703726984;
   rkAinvT[2][0] = (KPP_REAL)4.644370709252171185822941421408064;
   rkAinvT[2][1] = (KPP_REAL)3.677814645373914407088529289295970;
   rkAinvT[2][2] = (KPP_REAL)(-3.508761919567443321903661209182446);
} /* Gauss_Coefficients */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	The coefficients of the 3-stage Gauss method
	(given to ~30 accurate digits)
	The parameter b3 can be chosen by the user
	to tune the error estimator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Radau1A_Coefficients()
{
   /* The coefficients of the Radau1A method */
   KPP_REAL b0;
   /*b0 = (KPP_REAL)0.3; */
   b0 = (KPP_REAL)0.1;

   /* The coefficients of the Radau1A method */
   rkMethod = R1A;

   rkA[0][0] = (KPP_REAL) .1111111111111111111111111111111111;
   rkA[0][1] = (KPP_REAL)(-.1916383190435098943442935597058829);
   rkA[0][2] = (KPP_REAL)(.8052720793239878323318244859477174e-01);
   rkA[1][0] = (KPP_REAL).1111111111111111111111111111111111;
   rkA[1][1] = (KPP_REAL).2920734116652284630205027458970589;
   rkA[1][2] = (KPP_REAL)(-.481334970546573839513422644787591e-01);
   rkA[2][0] = (KPP_REAL).1111111111111111111111111111111111;
   rkA[2][1] = (KPP_REAL).5370223859435462728402311533676479;
   rkA[2][2] = (KPP_REAL).1968154772236604258683861429918299;

   rkB[0] = (KPP_REAL).1111111111111111111111111111111111;
   rkB[1] = (KPP_REAL).5124858261884216138388134465196080;
   rkB[2] = (KPP_REAL).3764030627004672750500754423692808;

   rkC[0] = ZERO;
   rkC[1] = (KPP_REAL).3550510257216821901802715925294109;
   rkC[2] = (KPP_REAL).8449489742783178098197284074705891;

   /* Classical error estimator, embedded solution: */
   rkBhat[0] = b0;
   rkBhat[1] = (KPP_REAL).11111111111111111111111111111111111-b0;
   rkBhat[2] = (KPP_REAL).51248582618842161383881344651960810;
   rkBhat[3] = (KPP_REAL).37640306270046727505007544236928079;

   /* New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j */
   rkD[0] = (KPP_REAL).3333333333333333333333333333333333;
   rkD[1] = (KPP_REAL)(-.8914115380582557157653087040196127);
   rkD[2] = (KPP_REAL)(1.558078204724922382431975370686279);

   /* Classical error estimator: */
   /* H* Sum (b_j-bhat_j) f(Z_j) = H*E(0)*F(0) + Sum E_j Z_j */
   rkE[0] = (KPP_REAL).2748888295956773677478286035994148*b0;
   rkE[1] = (KPP_REAL)(-1.374444147978386838739143017997074)*b0;
   rkE[2] = (KPP_REAL)(-1.335337922441686804550326197041126)*b0;
   rkE[3] = (KPP_REAL).235782604058977333559011782643466*b0;

   /* Sdirk error estimator */
   rkBgam[0] = ZERO;
   rkBgam[1] = (KPP_REAL)(.1948150124588532186183490991130616e-01);
   rkBgam[2] = (KPP_REAL).7575249005733381398986810981093584;
   rkBgam[3] = (KPP_REAL)(-.518952314149008295083446116200793e-01);
   rkBgam[4] = (KPP_REAL).2748888295956773677478286035994148;

   /* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta[0] = (KPP_REAL)(-1.224370034375505083904362087063351);
   rkTheta[1] = (KPP_REAL).9340045331532641409047527962010133;
   rkTheta[2] = (KPP_REAL).4656990124352088397561234800640929;

   /* ELO = local order of classical error estimator */
   rkELO = (KPP_REAL)4.0;

   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ~~~> Diagonalize the RK matrix:
    rkTinv * inv(rkA) * rkT =
              |  rkGamma      0           0     |
              |      0      rkAlpha   -rkBeta   |
              |      0      rkBeta     rkAlpha  |
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = (KPP_REAL)3.637834252744495732208418513577775;
   rkAlpha = (KPP_REAL)2.681082873627752133895790743211112;
   rkBeta  = (KPP_REAL)3.050430199247410569426377624787569;

   rkT[0][0] = (KPP_REAL).424293819848497965354371036408369;
   rkT[0][1] = (KPP_REAL)(-.3235571519651980681202894497035503);
   rkT[0][2] = (KPP_REAL)(-.522137786846287839586599927945048);
   rkT[1][0] = (KPP_REAL)(.57594609499806128896291585429339e-01);
   rkT[1][1] = (KPP_REAL)(.3148663231849760131614374283783e-02);
   rkT[1][2] = (KPP_REAL).452429247674359778577728510381731;
   rkT[2][0] = ONE;
   rkT[2][1] = ONE;
   rkT[2][2] = ZERO;

   rkTinv[0][0] = (KPP_REAL)1.233523612685027760114769983066164;
   rkTinv[0][1] = (KPP_REAL)1.423580134265707095505388133369554;
   rkTinv[0][2] = (KPP_REAL).3946330125758354736049045150429624;
   rkTinv[1][0] = (KPP_REAL)(-1.233523612685027760114769983066164);
   rkTinv[1][1] = (KPP_REAL)(-1.423580134265707095505388133369554);
   rkTinv[1][2] = (KPP_REAL).6053669874241645263950954849570376;
   rkTinv[2][0] = (KPP_REAL)(-.1484438963257383124456490049673414);
   rkTinv[2][1] = (KPP_REAL)2.038974794939896109682070471785315;
   rkTinv[2][2] = (KPP_REAL)(-.544501292892686735299355831692542e-01);

   rkTinvAinv[0][0] = (KPP_REAL)4.487354449794728738538663081025420;
   rkTinvAinv[0][1] = (KPP_REAL)5.178748573958397475446442544234494;
   rkTinvAinv[0][2] = (KPP_REAL)1.435609490412123627047824222335563;
   rkTinvAinv[1][0] = (KPP_REAL)(-2.854361287939276673073807031221493);
   rkTinvAinv[1][1] = (KPP_REAL)(-1.003648660720543859000994063139137e+01);
   rkTinvAinv[1][2] = (KPP_REAL)1.789135380979465422050817815017383;
   rkTinvAinv[2][0] = (KPP_REAL)(-4.160768067752685525282947313530352);
   rkTinvAinv[2][1] = (KPP_REAL)1.124128569859216916690209918405860;
   rkTinvAinv[2][2] = (KPP_REAL)1.700644430961823796581896350418417;

   rkAinvT[0][0] = (KPP_REAL)1.543510591072668287198054583233180;
   rkAinvT[0][1] = (KPP_REAL)(-2.460228411937788329157493833295004);
   rkAinvT[0][2] = (KPP_REAL)(-.412906170450356277003910443520499);
   rkAinvT[1][0] = (KPP_REAL).209519643211838264029272585946993;
   rkAinvT[1][1] = (KPP_REAL)1.388545667194387164417459732995766;
   rkAinvT[1][2] = (KPP_REAL)1.20339553005832004974976023130002;
   rkAinvT[2][0] = (KPP_REAL)3.637834252744495732208418513577775;
   rkAinvT[2][1] = (KPP_REAL)2.681082873627752133895790743211112;
   rkAinvT[2][2] = (KPP_REAL)(-3.050430199247410569426377624787569);
} /* Radau1A_Coefficients */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	The coefficients of the 4-stage Lobatto-3A method
	(given to ~30 accurate digits)
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Lobatto3A_Coefficients()
{
   /* The coefficients of the Lobatto-3A method */
   rkMethod = L3A;

   rkA[0][0] = ZERO;
   rkA[0][1] = ZERO;
   rkA[0][2] = ZERO;
   rkA[0][3] = ZERO;
   rkA[1][0] = (KPP_REAL).11030056647916491413674311390609397;
   rkA[1][1] = (KPP_REAL).1896994335208350858632568860939060;
   rkA[1][2] = (KPP_REAL)(-.339073642291438837776604807792215e-01);
   rkA[1][3] = (KPP_REAL)(.1030056647916491413674311390609397e-01);
   rkA[2][0] = (KPP_REAL)(.73032766854168419196590219427239365e-01);
   rkA[2][1] = (KPP_REAL).4505740308958105504443271474458881;
   rkA[2][2] = (KPP_REAL).2269672331458315808034097805727606;
   rkA[2][3] = (KPP_REAL)(-.2696723314583158080340978057276063e-01);
   rkA[3][0] = (KPP_REAL)(.83333333333333333333333333333333333e-01);
   rkA[3][1] = (KPP_REAL).4166666666666666666666666666666667;
   rkA[3][2] = (KPP_REAL).4166666666666666666666666666666667;
   rkA[3][3] = (KPP_REAL)(.8333333333333333333333333333333333e-01);

   rkB[0] = (KPP_REAL)(.83333333333333333333333333333333333e-01);
   rkB[1] = (KPP_REAL).4166666666666666666666666666666667;
   rkB[2] = (KPP_REAL).4166666666666666666666666666666667;
   rkB[3] = (KPP_REAL)(.8333333333333333333333333333333333e-01);

   rkC[0] = ZERO;
   rkC[1] = (KPP_REAL).2763932022500210303590826331268724;
   rkC[2] = (KPP_REAL).7236067977499789696409173668731276;
   rkC[3] = ONE;

   /* New solution: H*Sum B_j*f(Z_j) = Sum D_j*Z_j */
   rkD[0] = ZERO;
   rkD[1] = ZERO;
   rkD[2] = ZERO;
   rkD[3] = ONE;

   /* Classical error estimator, embedded solution: */
   rkBhat[0] = (KPP_REAL)(.90909090909090909090909090909090909e-01);
   rkBhat[1] = (KPP_REAL).39972675774621371442114262372173276;
   rkBhat[2] = (KPP_REAL).43360657558711961891219070961160058;
   rkBhat[3] = (KPP_REAL)(.15151515151515151515151515151515152e-01);

   /* Classical error estimator: */
   /* H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j */

   rkE[0] = (KPP_REAL)(.1957403846510110711315759367097231e-01);
   rkE[1] = (KPP_REAL)(-.1986820345632580910316020806676438);
   rkE[2] = (KPP_REAL).1660586371214229125096727578826900;
   rkE[3] = (KPP_REAL)(-.9787019232550553556578796835486154e-01);

   /* Sdirk error estimator: */
   rkF[0] = ZERO;
   rkF[1] = (KPP_REAL)(-.66535815876916686607437314126436349);
   rkF[2] = (KPP_REAL)1.7419302743497277572980407931678409;
   rkF[3] = (KPP_REAL)(-1.2918865386966730694684011822841728);

   /* ELO = local order of classical error estimator */
   rkELO = (KPP_REAL)4.0;

   /* Sdirk error estimator: */
   rkBgam[0] = (KPP_REAL)(.2950472755430528877214995073815946e-01);
   rkBgam[1] = (KPP_REAL).5370310883226113978352873633882769;
   rkBgam[2] = (KPP_REAL).2963022450107219354980459699450564;
   rkBgam[3] = (KPP_REAL)(-.7815248400375080035021681445218837e-01);
   rkBgam[4] = (KPP_REAL).2153144231161121782447335303806956;

   /* H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j */
   rkTheta[0] = (KPP_REAL)0.0;
   rkTheta[1] = (KPP_REAL)(-.6653581587691668660743731412643631);
   rkTheta[2] = (KPP_REAL)1.741930274349727757298040793167842;
   rkTheta[3] = (KPP_REAL)(-.291886538696673069468401182284174);


   /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ~~~> Diagonalize the RK matrix:
    rkTinv * inv(rkA) * rkT =
              |  rkGamma      0           0     |
              |      0      rkAlpha   -rkBeta   |
              |      0      rkBeta     rkAlpha  |
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   rkGamma = (KPP_REAL)4.644370709252171185822941421408063;
   rkAlpha = (KPP_REAL)3.677814645373914407088529289295968;
   rkBeta  = (KPP_REAL)3.508761919567443321903661209182446;

   rkT[0][0] = (KPP_REAL)(.5303036326129938105898786144870856e-01);
   rkT[0][1] = (KPP_REAL)(-.7776129960563076320631956091016914e-01);
   rkT[0][2] = (KPP_REAL)(.6043307469475508514468017399717112e-02);
   rkT[1][0] = (KPP_REAL).2637242522173698467283726114649606;
   rkT[1][1] = (KPP_REAL).2193839918662961493126393244533346;
   rkT[1][2] = (KPP_REAL).3198765142300936188514264752235344;
   rkT[2][0] = ONE;
   rkT[2][1] = ZERO;
   rkT[2][2] = ZERO;

   rkTinv[0][0] = (KPP_REAL)7.695032983257654470769069079238553;
   rkTinv[0][1] = (KPP_REAL)(-.1453793830957233720334601186354032);
   rkTinv[0][2] = (KPP_REAL).6302696746849084900422461036874826;
   rkTinv[1][0] = (KPP_REAL)(-7.695032983257654470769069079238553);
   rkTinv[1][1] = (KPP_REAL).1453793830957233720334601186354032;
   rkTinv[1][2] = (KPP_REAL).3697303253150915099577538963125174;
   rkTinv[2][0] = (KPP_REAL)(-1.066660885401270392058552736086173);
   rkTinv[2][1] = (KPP_REAL)3.146358406832537460764521760668932;
   rkTinv[2][2] = (KPP_REAL)(-.7732056038202974770406168510664738);

   rkTinvAinv[0][0] = (KPP_REAL)35.73858579417120341641749040405149;
   rkTinvAinv[0][1] = (KPP_REAL)(-.675195748578927863668368190236025);
   rkTinvAinv[0][2] = (KPP_REAL)2.927206016036483646751158874041632;
   rkTinvAinv[1][0] = (KPP_REAL)(-24.55824590667225493437162206039511);
   rkTinvAinv[1][1] = (KPP_REAL)(-10.50514413892002061837750015342036);
   rkTinvAinv[1][2] = (KPP_REAL)4.072793983963516353248841125958369;
   rkTinvAinv[2][0] = (KPP_REAL)(-30.92301972744621647251975054630589);
   rkTinvAinv[2][1] = (KPP_REAL)12.08182467154052413351908559269928;
   rkTinvAinv[2][2] = (KPP_REAL)(-1.546411207640594954081233702132946);

   rkAinvT[0][0] = (KPP_REAL).2462926658317812882584158369803835;
   rkAinvT[0][1] = (KPP_REAL)(-.2647871194157644619747121197289574);
   rkAinvT[0][2] = (KPP_REAL).2950720515900466654896406799284586;
   rkAinvT[1][0] = (KPP_REAL)1.224833192317784474576995878738004;
   rkAinvT[1][1] = (KPP_REAL)1.929224190340981580557006261869763;
   rkAinvT[1][2] = (KPP_REAL).4066803323234419988910915619080306;
   rkAinvT[2][0] = (KPP_REAL)4.644370709252171185822941421408064;
   rkAinvT[2][1] = (KPP_REAL)3.677814645373914407088529289295968;
   rkAinvT[2][2] = (KPP_REAL)(-3.508761919567443321903661209182446);
} /* Lobatto3A_Coefficients */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	END SUBROUTINE RungeKutta ! and all its internal procedures
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void FUN_CHEM(KPP_REAL T, KPP_REAL V[], KPP_REAL FCT[])
{
   KPP_REAL Told;

   Told = TIME;
   TIME = T;
   Update_SUN();
   Update_RCONST();
   Update_PHOTO();
   TIME = Told;

   Fun(V, FIX, RCONST, FCT);

} /* FUN_CHEM */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void JAC_CHEM (KPP_REAL T, KPP_REAL V[], KPP_REAL JF[])
{
   KPP_REAL Told;

   Told = TIME;
   TIME = T;
   Update_SUN();
   Update_RCONST();
   Update_PHOTO();
   TIME = Told;
   Jac_SP(V, FIX, RCONST, JF);
} /* JAC_CHEM */
