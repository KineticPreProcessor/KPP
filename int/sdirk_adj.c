
#define MAX(a,b) ( ((a) >= (b)) ? (a):(b)  )
#define MIN(b,c) ( ((b) < (c))  ? (b):(c)  )
#define ABS(x)	 ( ((x) >= 0 )  ? (x):(-x) )
#define SIGN(x,y)( ( (x*y) >= 0 ) ?(x):(-x) )
#define SQRT(d)  ( pow((d),0.5)  )

/* Numerical Constants */
#define ZERO	    (KPP_REAL)0.0
#define ONE	    (KPP_REAL)1.0
#define HALF        (KPP_REAL)0.5
#define DeltaMin    (KPP_REAL)1.0e-5
enum boolean { FALSE=0, TRUE=1 };
  
/*~~~>  Statistics on the work performed by the SDIRK method */
enum statistics { Nfun=1, Njac=2, Nstp=3, Nacc=4, Nrej=5, Ndec=6, Nsol=7, 
		  Nsng=8, Ntexit=1, Nhexit=2, Nhnew=3 };

/*~~~>  SDIRK method coefficients, up to 5 stages */
int Smax = 5;
enum ros_Params { S2A=1, S2B=2, S3A=3, S4A=4, S4B=5 };
KPP_REAL rkGamma, rkA[5][5], rkB[5], rkC[5], rkD[5], rkE[5], 
  rkBhat[5], rkELO, rkAlpha[5][5], rkTheta[5][5];
int sdMethod, rkS; /* The number of stages */

/*~~~>  Checkpoints in memory buffers */
int stack_ptr = -1; /* last written entry in checkpoint */
KPP_REAL *chk_H, *chk_T;
KPP_REAL **chk_Y;
KPP_REAL ***chk_Z;
int **chk_P;
#ifdef FULL_ALGEBRA
  KPP_REAL ***chk_J;
#else
  KPP_REAL **chk_J;
#endif

/* Function Headers */
void INTEGRATE_ADJ( int NADJ, KPP_REAL Y[], KPP_REAL Lambda[][NVAR], 
		    KPP_REAL TIN, KPP_REAL TOUT, KPP_REAL ATOL_adj[][NVAR], 
		    KPP_REAL RTOL_adj[][NVAR], int ICNTRL_U[], 
		    KPP_REAL RCNTRL_U[], int ISTATUS_U[], 
		    KPP_REAL RSTATUS_U[] );
int SDIRKADJ( int N, int NADJ, KPP_REAL Tinitial, KPP_REAL Tfinal, KPP_REAL Y[],
	      KPP_REAL Lambda[][NVAR], KPP_REAL RelTol[], KPP_REAL AbsTol[], 
	      KPP_REAL RelTol_adj[][NVAR], KPP_REAL AbsTol_adj[][NVAR], 
	      KPP_REAL RCNTRL[], int ICNTRL[], KPP_REAL RSTATUS[], 
	      int ISTATUS[]);
int SDIRK_FwdInt( int N, KPP_REAL Tinitial, KPP_REAL Tfinal, KPP_REAL Y[], 
		  KPP_REAL Hmax, KPP_REAL Hmin, KPP_REAL Hstart, 
		  KPP_REAL Roundoff, KPP_REAL AbsTol[], KPP_REAL RelTol[], 
		  int ISTATUS[], KPP_REAL RSTATUS[], int Max_no_steps, 
		  int NewtonMaxit, KPP_REAL NewtonTol, KPP_REAL ThetaMin, 
		  KPP_REAL FacSafe, KPP_REAL FacMin, KPP_REAL FacMax, 
		  KPP_REAL FacRej, KPP_REAL Qmin, KPP_REAL Qmax, int ITOL, 
		  int SaveLU );
int SDIRK_DadjInt( int N, int NADJ, KPP_REAL Lambda[][NVAR], int SaveLU, 
		   int ISTATUS[], int ITOL, KPP_REAL AbsTol_adj[][NVAR], 
		   KPP_REAL RelTol_adj[][NVAR], int NewtonMaxit,
		   KPP_REAL ThetaMin, KPP_REAL NewtonTol, int DirectADJ );
void SDIRK_AllocBuffers( int Max_no_steps, int rkS, int SaveLU );
void SDIRK_FreeBuffers( int Max_no_steps, int SaveLU );
void SDIRK_Push( KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL Z[][NVAR], 
		 KPP_REAL E[], int P[], int Max_no_steps, int SaveLU );
void SDIRK_Pop( KPP_REAL* T, KPP_REAL* H, KPP_REAL* Y, KPP_REAL* Z, 
		KPP_REAL* E, int* P, int SaveLU );
void SDIRK_ErrorScale( int N, int ITOL, KPP_REAL AbsTol[], KPP_REAL RelTol[],
		       KPP_REAL Y[],KPP_REAL SCAL[]);
KPP_REAL SDIRK_ErrorNorm( int N, KPP_REAL Y[], KPP_REAL SCAL[] );
int SDIRK_ErrorMsg( int code, KPP_REAL T, KPP_REAL H );
void SDIRK_PrepareMatrix( KPP_REAL H, KPP_REAL T, KPP_REAL Y[], KPP_REAL FJAC[],
			  int SkipJac, int SkipLU, KPP_REAL E[], int IP[], 
			  int Reject, int ISING, int ISTATUS[] );
void SDIRK_Solve ( char Transp, KPP_REAL H, int N, KPP_REAL E[], 
		   int IP[], int ISING, KPP_REAL RHS[], int ISTATUS[] );
void Sdirk4a();
void Sdirk4b();
void Sdirk2a();
void Sdirk2b();
void Sdirk3a();
void FUN_CHEM( KPP_REAL T, KPP_REAL Y[], KPP_REAL P[] );
void JAC_CHEM( KPP_REAL T, KPP_REAL Y[], KPP_REAL JV[] );
KPP_REAL WLAMCH( char C );
void Set2Zero( int N, KPP_REAL Y[] );
void WAXPY( int N, KPP_REAL Alpha, KPP_REAL X[], int incX, KPP_REAL Y[], 
	    int incY );
void WADD( int N, KPP_REAL Y[], KPP_REAL Z[], KPP_REAL TMP[] );
void WSCAL( int N, KPP_REAL Alpha, KPP_REAL X[], int incX );
void JacTR_SP_Vec( KPP_REAL Jac[], KPP_REAL Fcn[], KPP_REAL K[] );
void KppSolve( KPP_REAL A[], KPP_REAL b[] );
void KppSolveTR( KPP_REAL JVS[], KPP_REAL X[], KPP_REAL XX[] );
int KppDecomp( KPP_REAL A[] );
void Update_SUN();
void Update_RCONST();
void Fun( KPP_REAL Y[], KPP_REAL FIX[], KPP_REAL RCONST[], KPP_REAL Ydot[] );
void Jac_SP( KPP_REAL Y[], KPP_REAL FIX[], KPP_REAL RCONST[], KPP_REAL Ydot[] );
void Set2Zero( int N, KPP_REAL Y[] );

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void INTEGRATE_ADJ( int NADJ, KPP_REAL Y[], KPP_REAL Lambda[][NVAR], 
		    KPP_REAL TIN, KPP_REAL TOUT, KPP_REAL ATOL_adj[][NVAR], 
		    KPP_REAL RTOL_adj[][NVAR], int ICNTRL_U[], 
		    KPP_REAL RCNTRL_U[], int ISTATUS_U[], 
		    KPP_REAL RSTATUS_U[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  //int Ntotal = 0;
  KPP_REAL RCNTRL[20], RSTATUS[20], T1, T2;
  int ICNTRL[20], ISTATUS[20], Ierr, i;

  for(i=0; i<20; i++) {
    ICNTRL[i] = 0;
    RCNTRL[i] = (KPP_REAL)0.0;
    ISTATUS[i] = 0;
    RSTATUS[i] = (KPP_REAL)0.0;
  }

  /*~~~> fine-tune the integrator: */
  ICNTRL[4] = 8;    /* Max no. of Newton iterations */
  ICNTRL[6] = 1;   /* Adjoint solution by: 0=Newton, 1=direct */
  ICNTRL[7] = 1;    /* Save fwd LU factorization: 0 = do *not* save, 1 = save */

  /* If optional parameters are given, and if they are >0, 
     then they overwrite default settings. */ 
  //if(ICNTRL_U != NULL) { /* Check to see if ICNTRL_U is not NULL */
  //  for(i=0; i<20; i++) {
  //	if(ICNTRL_U[i] > 0)
  //	  ICNTRL[i] = ICNTRL_U[i];
  //  }
  //}
  //
  //if(RCNTRL_U != NULL) { /* Check to see if RCNTRL_U is not NULL */
  //  for(i=0; i<20; i++) {
  //	if(RCNTRL_U[i] > 0)
  //	  RCNTRL[i] = RCNTRL_U[i];
  //  }
  //}
  
  T1 = TIN; 
  T2 = TOUT;
  Ierr = SDIRKADJ( NVAR, NADJ, T1, T2, Y, Lambda, RTOL, ATOL, ATOL_adj, 
		   RTOL_adj, RCNTRL, ICNTRL, RSTATUS, ISTATUS );

  /*~~~> Debug option: number of steps */
  // Ntotal = Ntotal + ISTATUS(Nstp)
  // printf( "NSTEP=%d Ntotal=%d O3=%e NO2=%e\n", ISTATUS(Nstp), Ntotal, 
  //          VAR(ind_O3), VAR(ind_NO2) );    

  if (Ierr < 0)
    printf("SDIRK: Unsuccessful exit at T=%f (Ierr=%d)\n", TIN, Ierr );
   
  /* if optional parameters are given for output they to return information */
  //if(ISTATUS_U != NULL) { /* Check to see if ISTATUS_U is not NULL */
  //	for(i=0; i<20; i++)
  //	  ISTATUS_U[i] = ISTATUS[i];
  //} 
  //
  //if(RSTATUS_U != NULL) { /* Check to see if RSTATUS_U is not NULL */
  //	for(i=0; i<20; i++)
  //	  RSTATUS_U[i] = RSTATUS[i];
  //}
  //
  //Ierr_U = Ierr;
  
} /* END of SUBROUTINE INTEGRATE_ADJ */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int SDIRKADJ( int N, int NADJ, KPP_REAL Tinitial, KPP_REAL Tfinal, KPP_REAL Y[],
	      KPP_REAL Lambda[][NVAR], KPP_REAL RelTol[], KPP_REAL AbsTol[], 
	      KPP_REAL RelTol_adj[][NVAR], KPP_REAL AbsTol_adj[][NVAR], 
	      KPP_REAL RCNTRL[], int ICNTRL[], KPP_REAL RSTATUS[], 
	      int ISTATUS[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Solves the system y'=F(t,y) using a Singly-Diagonally-Implicit
    Runge-Kutta (SDIRK) method.

    This implementation is based on the book and the code Sdirk4:

      E. Hairer and G. Wanner
      "Solving ODEs II. Stiff and differential-algebraic problems".
      Springer series in computational mathematics, Springer-Verlag, 1996.
    This code is based on the SDIRK4 routine in the above book.

    Methods:
            * Sdirk 2a, 2b: L-stable, 2 stages, order 2                  
            * Sdirk 3a:     L-stable, 3 stages, order 2, adjoint-invariant   
            * Sdirk 4a, 4b: L-stable, 5 stages, order 4                  

    (C)  Adrian Sandu, July 2005
    Virginia Polytechnic Institute and State University
    Contact: sandu@cs.vt.edu
    Revised by Philipp Miehe and Adrian Sandu, May 2006     
    Translation F90 to C by Paul Eller, May 2007
    This implementation is part of KPP - the Kinetic PreProcessor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>   INPUT ARGUMENTS:

-     Y[NVAR]    = vector of initial conditions (at T=Tinitial)
-    [Tinitial,Tfinal]  = time range of integration
     (if Tinitial>Tfinal the integration is performed backwards in time)
-    RelTol, AbsTol = user precribed accuracy
- SUBROUTINE ode_Fun( T, Y, Ydot ) = ODE function,
                       returns Ydot = Y' = F(T,Y)
- SUBROUTINE ode_Fun( T, Y, Ydot ) = Jacobian of the ODE function,
                       returns Jcb = dF/dY
-    ICNTRL[1:20]    = integer inputs parameters
-    RCNTRL[1:20]    = real inputs parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>     OUTPUT ARGUMENTS:

-    Y[NVAR]         -> vector of final states (at T->Tfinal)
-    ISTATUS[0:19]   -> integer output parameters
-    RSTATUS[0:19]   -> real output parameters
-    Ierr            -> job status upon return
                        success (positive value) or
                        failure (negative value)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>     INPUT PARAMETERS:

    Note: For input parameters equal to zero the default values of the
       corresponding variables are used.

    Note: For input parameters equal to zero the default values of the
          corresponding variables are used.
~~~>  
    ICNTRL[0] = not used

    ICNTRL[1] = 0: AbsTol, RelTol are NVAR-dimensional vectors
              = 1: AbsTol, RelTol are scalars

    ICNTRL[2] = Method

    ICNTRL[3]  -> maximum number of integration steps
        For ICNTRL[3]=0 the default value of 1500 is used
        Note: use a conservative estimate, since the checkpoint
              buffers are allocated to hold Max_no_steps

    ICNTRL[4]  -> maximum number of Newton iterations
        For ICNTRL[4]=0 the default value of 8 is used

    ICNTRL[5]  -> starting values of Newton iterations:
        ICNTRL[5]=0 : starting values are interpolated (the default)
        ICNTRL[5]=1 : starting values are zero

    ICNTRL[6]  -> method to solve ADJ equations:
        ICNTRL[6]=0 : modified Newton re-using LU (the default)
        ICNTRL[6]=1 : direct solution(additional one LU factorization per stage)

    ICNTRL[7]  -> checkpointing the LU factorization at each step:
        ICNTRL[7]=0 : do *not* save LU factorization (the default)
        ICNTRL[7]=1 : save LU factorization
        Note: if ICNTRL[6]=1 the LU factorization is *not* saved

~~~>  Real parameters

    RCNTRL[0]  -> Hmin, lower bound for the integration step size
                  It is strongly recommended to keep Hmin = ZERO
    RCNTRL[1]  -> Hmax, upper bound for the integration step size
    RCNTRL[2]  -> Hstart, starting value for the integration step size

    RCNTRL[3]  -> FacMin, lower bound on step decrease factor (default=0.2)
    RCNTRL[4]  -> FacMax, upper bound on step increase factor (default=6)
    RCNTRL[5]  -> FacRej, step decrease factor after multiple rejections
                 (default=0.1)
    RCNTRL[6]  -> FacSafe, by which the new step is slightly smaller
                  than the predicted value  (default=0.9)
    RCNTRL[7]  -> ThetaMin. If Newton convergence rate smaller
                  than ThetaMin the Jacobian is not recomputed;
                  (default=0.001)
    RCNTRL[8]  -> NewtonTol, stopping criterion for Newton's method
                  (default=0.03)
    RCNTRL[9] -> Qmin
    RCNTRL[10] -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
                  step size is kept constant and the LU factorization
                  reused (default Qmin=1, Qmax=1.2)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>     OUTPUT PARAMETERS:

    Note: each call to Rosenbrock adds the current no. of fcn calls
      to previous value of ISTATUS[0], and similar for the other params.
      Set ISTATUS[1:10] = 0 before call to avoid this accumulation.

    ISTATUS[0] = No. of function calls
    ISTATUS[1] = No. of jacobian calls
    ISTATUS[2] = No. of steps
    ISTATUS[3] = No. of accepted steps
    ISTATUS[4] = No. of rejected steps (except at the beginning)
    ISTATUS[5] = No. of LU decompositions
    ISTATUS[6] = No. of forward/backward substitutions
    ISTATUS[7] = No. of singular matrix decompositions

    RSTATUS[0]  -> Texit, the time corresponding to the
                   computed Y upon return
    RSTATUS[1]  -> Hexit,last accepted step before return
    RSTATUS[2]  -> Hnew, last predicted step before return
        For multiple restarts, use Hnew as Hstart in the following run

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* Local variables */     
  KPP_REAL Hmin=0.0, Hmax=0.0, Hstart=0.0, Roundoff, FacMin=0.0, FacMax=0.0, 
    FacSafe=0.0, FacRej=0.0, ThetaMin, NewtonTol,Qmin, Qmax;
  int SaveLU, DirectADJ; /* Boolean variables */                 
  int ITOL, NewtonMaxit, Max_no_steps=0, i, Ierr=0;

  stack_ptr = -1;

/*~~~>  Initialize statistics */
  for(i=0; i<20; i++) {
    ISTATUS[i] = 0;
    RSTATUS[i] = ZERO;
  }

/*~~~>  For Scalar tolerances (ICNTRL[1] != 0)  the code uses AbsTol[0] 
            and RelTol[0]
	For Vector tolerances (ICNTRL[1] == 0) the code uses AbsTol[0:NVAR-1] 
            and RelTol[0:NVAR-1]*/
  if (ICNTRL[1]==0)
    ITOL = 1;
  else
    ITOL = 0;

/*~~~> ICNTRL(3) - method selection */
  switch (ICNTRL[2]) {
    case 0:
    case 1:
      Sdirk2a();
      break;
    case 2:
      Sdirk2b();
      break;
    case 3:
      Sdirk3a();
      break;
    case 4:
      Sdirk4a();
      break;
    case 5:
      Sdirk4b();
      break;
    default:
      Sdirk2a();
  }
      
/*~~~>   The maximum number of time steps admitted */
  if (ICNTRL[3] == 0)
    Max_no_steps = 200000;
  else if (ICNTRL[3] > 0)
    Max_no_steps = ICNTRL[3];
  else {
    printf("User-selected ICNTRL(4)=%d", ICNTRL[3]);
    Ierr = SDIRK_ErrorMsg(-1,Tinitial,ZERO);
  }

/*~~~>The maximum number of Newton iterations admitted */
  if(ICNTRL[4]==0)
    NewtonMaxit = 8;
  else {
    NewtonMaxit=ICNTRL[4];
    if(NewtonMaxit <=0) {
      printf("User-selected ICNTRL(5)=%d", ICNTRL[4] );
      Ierr = SDIRK_ErrorMsg(-2,Tinitial,ZERO);
    }
  }

/*~~~> Solve ADJ equations directly or by Newton iterations */
  DirectADJ = (ICNTRL[6] == 1);
 
/*~~~> Save or not the forward LU factorization */
  SaveLU = ((ICNTRL[7] != 0) && (DirectADJ == 0));

/*~~~>  Unit roundoff (1+Roundoff>1) */
  Roundoff = WLAMCH('E');

/*~~~>  Lower bound on the step size: (positive value) */
  if (RCNTRL[0] == ZERO)
    Hmin = ZERO;
  else if (RCNTRL[0] > ZERO)
    Hmin = RCNTRL[0];
  else {
    printf("User-selected RCNTRL[0]=%f", RCNTRL[0]);
    Ierr = SDIRK_ErrorMsg(-3,Tinitial,ZERO);
  }
   
/*~~~>  Upper bound on the step size: (positive value) */
  if (RCNTRL[1] == ZERO)
    Hmax = ABS(Tfinal-Tinitial);
  else if (RCNTRL[1] > ZERO)
    Hmax = MIN( ABS(RCNTRL[1]), ABS(Tfinal-Tinitial) );
  else {
    printf("User-selected RCNTRL[1]=%f", RCNTRL[1]);
    Ierr = SDIRK_ErrorMsg(-3,Tinitial,ZERO);
  }

/*~~~>  Starting step size: (positive value) */
  if (RCNTRL[2] == ZERO)
    Hstart = MAX( Hmin, Roundoff);
  else if (RCNTRL[2] > ZERO)
    Hstart = MIN( ABS(RCNTRL[2]), ABS(Tfinal-Tinitial) );
  else {
    printf("User-selected Hstart: RCNTRL[2]=%f", RCNTRL[2]);
    Ierr = SDIRK_ErrorMsg(-3,Tinitial,ZERO);
  }

/*~~~>  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax */
  if (RCNTRL[3] == ZERO)
    FacMin = (KPP_REAL)0.2;
  else if (RCNTRL[3] > ZERO)
    FacMin = RCNTRL[3];
  else {
    printf("User-selected FacMin: RCNTRL[3]=%f", RCNTRL[3]);
    Ierr = SDIRK_ErrorMsg(-4,Tinitial,ZERO);
  }
   
  if (RCNTRL[4] == ZERO)
    FacMax = (KPP_REAL)10.0;
  else if (RCNTRL[4] > ZERO)
    FacMax = RCNTRL[4];
  else {
    printf("User-selected FacMax: RCNTRL[4]=%f", RCNTRL[4]);
    Ierr = SDIRK_ErrorMsg(-4,Tinitial,ZERO);
  }

/*~~~>   FacRej: Factor to decrease step after 2 succesive rejections */
  if (RCNTRL[5] == ZERO)
    FacRej = (KPP_REAL)0.1;
  else if (RCNTRL[5] > ZERO)
    FacRej = RCNTRL[5];
  else {
    printf("User-selected FacRej: RCNTRL[5]=%f", RCNTRL[5]);
    Ierr = SDIRK_ErrorMsg(-4,Tinitial,ZERO);
  }

/* ~~~>   FacSafe: Safety Factor in the computation of new step size  */
  if (RCNTRL[6] == ZERO)
    FacSafe = (KPP_REAL)0.9;
  else if (RCNTRL[6] > ZERO)
    FacSafe = RCNTRL[6];
  else {
    printf("User-selected FacSafe: RCNTRL[6]=%f", RCNTRL[6]);
    Ierr = SDIRK_ErrorMsg(-4,Tinitial,ZERO);
  }

/*~~~> ThetaMin: decides whether the Jacobian should be recomputed */
  if (RCNTRL[7] == ZERO)
    ThetaMin = (KPP_REAL)1.0e-03;
  else
    ThetaMin = RCNTRL[7];

/*~~~> Stopping criterion for Newton's method */
  if (RCNTRL[8] == ZERO)
    NewtonTol = (KPP_REAL)3.0e-02;
  else
    NewtonTol = RCNTRL[8];

/* ~~~> Qmin, Qmax: IF Qmin < Hnew/Hold < Qmax, STEP SIZE = CONST. */
  if (RCNTRL[9] == ZERO)
    Qmin = ONE;
  else
    Qmin = RCNTRL[9];

  if (RCNTRL[10] == ZERO)
    Qmax = (KPP_REAL)1.2;
  else
    Qmax = RCNTRL [10];

/* ~~~>  Check if tolerances are reasonable */
  if (ITOL == 0) {
    if ((AbsTol[0]<=ZERO || RelTol[0])<=(((KPP_REAL)10.0)*Roundoff))
      Ierr = SDIRK_ErrorMsg(-5,Tinitial,ZERO);
  }
  else {
    for (i = 0; i < N; i++) {
      if ((AbsTol[i]<=ZERO)||(RelTol[i]<=((KPP_REAL)10.0)*Roundoff))
	Ierr = SDIRK_ErrorMsg(-5,Tinitial,ZERO);
    }
  }

  if (Ierr < 0)
    return Ierr;
    
/*~~~>  Allocate memory buffers */
  SDIRK_AllocBuffers(Max_no_steps, rkS, SaveLU);

/*~~~>  Call forward integration */
  Ierr = SDIRK_FwdInt( N, Tinitial, Tfinal, Y, Hmax, Hmin, Hstart, Roundoff, 
		       AbsTol, RelTol, ISTATUS, RSTATUS, Max_no_steps, 
		       NewtonMaxit, NewtonTol, ThetaMin, FacSafe, FacMin, 
		       FacMax, FacRej, Qmin, Qmax, ITOL, SaveLU );

/*~~~>  Call adjoint integration */  
  Ierr = SDIRK_DadjInt( N, NADJ, Lambda, SaveLU, ISTATUS, ITOL, AbsTol_adj, 
			RelTol_adj, NewtonMaxit, ThetaMin, NewtonTol, 
			DirectADJ );

/*~~~>  Free memory buffers */ 
  SDIRK_FreeBuffers(Max_no_steps, SaveLU);

  return Ierr;

} /* End of main SDIRK_ADJ */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int SDIRK_FwdInt( int N, KPP_REAL Tinitial, KPP_REAL Tfinal, KPP_REAL Y[], 
		  KPP_REAL Hmax, KPP_REAL Hmin, KPP_REAL Hstart, 
		  KPP_REAL Roundoff, KPP_REAL AbsTol[], KPP_REAL RelTol[], 
		  int ISTATUS[], KPP_REAL RSTATUS[], int Max_no_steps, 
		  int NewtonMaxit, KPP_REAL NewtonTol, KPP_REAL ThetaMin, 
		  KPP_REAL FacSafe, KPP_REAL FacMin, KPP_REAL FacMax, 
		  KPP_REAL FacRej, KPP_REAL Qmin, KPP_REAL Qmax, int ITOL, 
		  int SaveLU ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      
/*~~~> Local variables: */   
  KPP_REAL Z[Smax][NVAR], G[NVAR], TMP[NVAR], NewtonRate, SCAL[NVAR], RHS[NVAR],
    T, H, Theta=0.0, Hratio, NewtonPredictedErr, Qnewton, Err, Fac, Hnew, 
    Tdirection, NewtonIncrement, NewtonIncrementOld=0.0;
  int i, j, IER=0, istage, NewtonIter, IP[NVAR];

  /*Boolean Variables*/
  int Reject, FirstStep, SkipJac, SkipLU, NewtonDone, CycleTloop;
      
#ifdef FULL_ALGEBRA      
  KPP_REAL FJAC[NVAR][NVAR], E[NVAR][NVAR];
#else      
  KPP_REAL FJAC[LU_NONZERO], E[LU_NONZERO];
#endif

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~>   Initializations              */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  T = Tinitial;
  Tdirection = SIGN(ONE, Tfinal-Tinitial);
  H = MAX(ABS(Hmin), ABS(Hstart));
  if (ABS(H) <= ((KPP_REAL)10.0 * Roundoff))
    H = (KPP_REAL)(1.0e-06);
  H = MIN(ABS(H), Hmax);
  H = SIGN(H, Tdirection);
  SkipLU = 0;
  SkipJac = 0;
  Reject = 0;
  FirstStep = 1;
  CycleTloop = 0;

  SDIRK_ErrorScale(N, ITOL, AbsTol, RelTol, Y, SCAL);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~>  Time loop begins                */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  while((Tfinal-T)*Tdirection - Roundoff > ZERO) {  /* Tloop */

  /*~~~>  Compute E = 1/(h*gamma)-Jac and its LU decomposition */
    if (SkipLU == 0) { /* This time around skip the Jac update and LU */
      SDIRK_PrepareMatrix( H, T, Y, FJAC, SkipJac, SkipLU, E, IP, 
	                   Reject, IER, ISTATUS);
      if (IER != 0)
	return SDIRK_ErrorMsg(-8, T, H);
    }

    if (ISTATUS[Nstp] > Max_no_steps)
      return SDIRK_ErrorMsg(-6, T, H);

    if ((T + ((KPP_REAL)0.1) * H == T) || (ABS(H) <= Roundoff)) {
      return SDIRK_ErrorMsg(-7, T, H);
    }

    for (istage=0; istage < rkS; istage++) {  /*stages*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~>  Simplified Newton iterations          */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*~~~>  Starting values for Newton iterations */
      Set2Zero(N, &Z[istage][0]);

    /*~~~>   Prepare the loop-independent part of the right-hand side */
      Set2Zero(N, G);
      if (istage > 0) {
	for (j=0; j < istage; j++) {
	/* Gj(:) = sum_j Theta(i,j)*Zj(:) = H * sum_j A(i,j)*Fun(Zj(:)) */
          WAXPY(N, rkTheta[j][istage], &Z[j][0], 1, G, 1);
	/* Zi(:) = sum_j Alpha(i,j)*Zj(:) */
	  WAXPY(N, rkAlpha[j][istage],&Z[j][0], 1, &Z[istage][0], 1);
	}
      }

    /*~~~>  Initializations for Newton iteration */
      NewtonDone = 0; /* false */
      Fac = (KPP_REAL)0.5; /* Step reduction factor */

      for (NewtonIter=0; NewtonIter<NewtonMaxit; NewtonIter++ ) { /*NewtonLoop*/
      /*~~~>  Prepare the loop-dependent part of the right-hand side */
	WADD(N, Y, &Z[istage][0], TMP); /* TMP <- Y + Zi */
        FUN_CHEM(T+rkC[istage]*H, TMP, RHS); /* RHS <- Run(Y+Zi) */
	ISTATUS[Nfun]++;
      /* RHS[0:N-1] = G[0:N-1] - Z[istage][0:N-1] + (H*rkGamma)*RHS[1:N] */
	WSCAL(N, H*rkGamma, RHS, 1);
	WAXPY(N, -ONE, &Z[istage][0], 1, RHS, 1);
	WAXPY(N, ONE, G, 1, RHS, 1 );

      /*~~~>   Solve the linear system  */
	SDIRK_Solve('N', H, N, E, IP, IER, RHS, ISTATUS);
            
      /*~~~>   Check convergence of Newton iterations */
	NewtonIncrement = SDIRK_ErrorNorm(N, RHS, SCAL);
	if (NewtonIter == 0) {
	  Theta = ABS(ThetaMin);
	  NewtonRate = (KPP_REAL)2.0;
       	}
	else {
	  Theta = NewtonIncrement/NewtonIncrementOld;
   	  if (Theta < (KPP_REAL)0.99) {
	    NewtonRate = Theta/(ONE-Theta);
	  /* Predict error at the end of Newton process */
	    NewtonPredictedErr = (NewtonIncrement*pow(Theta,
				 (NewtonMaxit - (NewtonIter+1))/(ONE - Theta)));
	    if(NewtonPredictedErr >= NewtonTol) {
	    /* Non-convergence of Newton: predicted error too large*/
	      Qnewton = MIN((KPP_REAL)10.0, NewtonPredictedErr/NewtonTol);
	      Fac = (KPP_REAL)0.8 * pow(Qnewton, (-ONE / (1 + NewtonMaxit - 
							  NewtonIter + 1)));
	      break;
	    }
	  }
	  else  /* Non-convergence of Newton: Theta too large */ {
	    break;
	  }
	}

	NewtonIncrementOld = NewtonIncrement;

      /* Update solution: Z(:) <-- Z(:)+RHS(:) */
	WAXPY(N, ONE, RHS, 1, &Z[istage][0], 1);

      /* Check error in Newton iterations */
	NewtonDone=(NewtonRate*NewtonIncrement<=NewtonTol);
	if (NewtonDone == 1)
	  break;
      } /* end NewtonLoop for */
                      
      if(NewtonDone == 0) {
	H = Fac*H;
	Reject = 1;
	SkipJac = 1;
	SkipLU = 0;
        CycleTloop = 1;
      } /* end if */

      if (CycleTloop == 1) {
	CycleTloop=0;
	break;
      }
    /* End of implified Newton iterations */
    } /* end stages for */

    if (CycleTloop==0) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~>  Error estimation      */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      ISTATUS[Nstp]++;
      Set2Zero(N, TMP);

      for (i=0; i<rkS; i++) {
        if (rkE[i] != ZERO)
	  WAXPY(N, rkE[i], &Z[i][0], 1, TMP, 1);
      }

      SDIRK_Solve('N', H, N, E, IP, IER, TMP, ISTATUS);
      Err = SDIRK_ErrorNorm(N, TMP, SCAL);

    /*~~~~> Computation of new step size Hnew */
      Fac = FacSafe * pow((Err), (-ONE/rkELO));
      Fac = MAX(FacMin, MIN(FacMax, Fac));
      Hnew = H*Fac;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~>  Accept/Reject step    */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (Err < ONE) { /*~~~> Step is accepted */
        FirstStep = 0; /* false */
        ISTATUS[Nacc]++;

      /* Checkpoint solution */
        SDIRK_Push( T, H, Y, Z, E, IP, Max_no_steps, SaveLU );

      /*~~~> Update time and solution */
        T = T + H;
      /* Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:) */
        for(i=0; i<rkS; i++) {
	  if(rkD[i] != ZERO)
	    WAXPY(N, rkD[i], &Z[i][0], 1, Y, 1);
        }

      /*~~~> Update scaling coefficients */
        SDIRK_ErrorScale(N, ITOL, AbsTol, RelTol, Y, SCAL);

      /*~~~> Next time step */
        Hnew = Tdirection*MIN(ABS(Hnew), Hmax);

      /* Last T and H */
        RSTATUS[Ntexit] = T;
        RSTATUS[Nhexit] = H;
        RSTATUS[Nhnew] = Hnew;

      /* No step increase after a rejection */
        if (Reject==1)
	  Hnew = Tdirection*MIN(ABS(Hnew), ABS(H));
        Reject = 0; /* false */
        if ((T+Hnew/Qmin-Tfinal)*Tdirection > ZERO)
	  H = Tfinal-T;
        else {
	  Hratio = Hnew/H;
        /* If step not changed too much keep Jacobian and reuse LU */
          SkipLU = ((Theta <= ThetaMin) && (Hratio >= Qmin) 
		    && (Hratio <= Qmax));
	  if (SkipLU==0)
	    H = Hnew;
        }

        /* If convergence is fast enough, do not update Jacobian */
        /* SkipJac = (Theta <= ThetaMin); */
        SkipJac = 0;
      }
      else { /*~~~> Step is rejected */
        if ((FirstStep==1) || (Reject==1))
	  H = FacRej * H;
        else  
          H = Hnew;

        Reject = 1;
        SkipJac = 1;
        SkipLU = 0;
        if (ISTATUS[Nacc] >=1)
	  ISTATUS[Nrej]++;

      }
    } /* end CycleTloop if */
  } /* end Tloop */

  /* Successful return */
  return 1;

} /* end SDIRK_FwdInt */
         
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int SDIRK_DadjInt( int N, int NADJ, KPP_REAL Lambda[][NVAR], int SaveLU, 
		   int ISTATUS[], int ITOL, KPP_REAL AbsTol_adj[][NVAR], 
		   KPP_REAL RelTol_adj[][NVAR], int NewtonMaxit, 
		   KPP_REAL ThetaMin, KPP_REAL NewtonTol, int DirectADJ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables: */
  KPP_REAL Y[NVAR];
  KPP_REAL Z[Smax][NVAR], U[Smax][NADJ][NVAR], TMP[NVAR], G[NVAR], NewtonRate, 
    SCAL[NVAR], DU[NVAR], T, H, Theta, NewtonPredictedErr, NewtonIncrement, 
    NewtonIncrementOld=0.0;
  int i, j, IER=0, istage, iadj, NewtonIter, IP[NVAR], IP_adj[NVAR];
  int Reject=0, SkipJac, SkipLU, NewtonDone; /* Boolean */
      
#ifdef FULL_ALGEBRA      
  KPP_REAL E[NVAR][NVAR], Jac[NVAR][NVAR], E_adj[NVAR][NVAR];
#else      
  KPP_REAL E[LU_NONZERO], Jac[LU_NONZERO], E_adj[LU_NONZERO];
#endif      

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~>  Time loop begins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  while ( stack_ptr > -1 ) { /* Tloop */
        
  /*~~~>  Recover checkpoints for stage values and vectors */
    SDIRK_Pop( &T, &H, &Y[0], &Z[0][0], &E[0], &IP[0], SaveLU );

  /*~~~>  Compute E = 1/(h*gamma)-Jac and its LU decomposition */
    if (!SaveLU) {
      SkipJac = FALSE; 
      SkipLU = FALSE;
      SDIRK_PrepareMatrix ( H, T, Y, Jac, SkipJac, SkipLU, E, IP, Reject, IER, 
			    ISTATUS );
      if (IER != 0)
	return SDIRK_ErrorMsg(-8,T,H); 
    }  

    for( istage = rkS-1; istage >=0; istage--) { /* Stages Loop */
    /*~~~>  Jacobian of the current stage solution */
      for(i=0; i<N; i++)
	TMP[i] = Y[i] + Z[istage][i];
      JAC_CHEM(T+rkC[istage]*H,TMP,Jac);
      ISTATUS[Njac]++;
       
      if (DirectADJ) {
#ifdef FULL_ALGEBRA
	for(i=0; i<N; i++) {
	  for(j=0; j<N; j++)
	    E_adj[i][j] = -Jac[i][j];
	}
	for(i=0; i<N; i++)
	  E_adj[i][i] = E_adj[i][i] + ONE/(H*rkGamma);

	DGETRF( N, N, E_adj, N, IP_adj, IER );
#else
	for(i=0; i<LU_NONZERO; i++)
	  E_adj[i] = -Jac[i];
	for(i=0; i<NVAR; i++) {
	  j = LU_DIAG[i];
	  E_adj[j] = E_adj[j] + ONE/(H*rkGamma);
	}
        IER = KppDecomp (E_adj);
#endif
	ISTATUS[Ndec]++;
        if (IER != 0) {
	  printf("At stage %d the matrix used in adjoint computation is " 
		 "singular\n", istage);
          return SDIRK_ErrorMsg(-8,T,H);
	}
      }

      for(iadj=0; iadj<NADJ; iadj++) { /* adj loop */
       
      /*~~~> Update scaling coefficients */
	for(i=0; i<NVAR; i++)
	  SDIRK_ErrorScale(N, ITOL, &AbsTol_adj[iadj][i], &RelTol_adj[iadj][i], 
			   &Lambda[iadj][i], SCAL);
      
      /*~~~>   Prepare the loop-independent part of the right-hand side
	       G(:) = H*Jac^T*( B(i)*Lambda + sum_j A(j,i)*Uj(:) ) */
	for(i=0; i<N; i++)
	  G[i] = rkB[istage]*Lambda[iadj][i];
	if (istage+1 < rkS) {
	  for (j=istage+1; j<rkS; j++)
	    WAXPY(N,rkA[istage][j],&U[j][iadj][0],1,G,1);
	}
#ifdef FULL_ALGEBRA  
	TMP = MATMUL(TRANSPOSE(Jac),G); /* DZ <- Jac(Y+Z)*Y_tlm */
#else      
	JacTR_SP_Vec ( Jac, G, TMP );    
#endif     
	for(i=0; i<N; i++)
	  G[i] = H*TMP[i];
	
	if (DirectADJ) {
	  SDIRK_Solve ( 'T', H, N, E_adj, IP_adj, IER, G, ISTATUS );
	  for(i=0; i<N; i++)
	    U[istage][iadj][i] = G[i];
	} else {
	  /*~~~>  Initializations for Newton iteration */
	  Set2Zero(N,&U[istage][iadj][0]);
	  NewtonDone = FALSE;
            
	  /* Newton Loop */
	  for( NewtonIter=0; NewtonIter<NewtonMaxit; NewtonIter++) {

	  /*~~~>   Prepare the loop-dependent part of the right-hand side */
#ifdef FULL_ALGEBRA  
	    for(i=0; i<N; i++)
	      TMP = MATMUL(TRANSPOSE(Jac),U[istage][iadj][i]);    
#else      
	    for(i=0; i<N; i++)
	      JacTR_SP_Vec ( Jac, &U[istage][iadj][i], TMP );    
#endif      
	    for(i=0; i<N; i++)
	      DU[i] = U[istage][iadj][i] - (H*rkGamma)*TMP[i] - G[i];

	  /*~~~>   Solve the linear system */
            SDIRK_Solve ( 'T', H, N, E, IP, IER, DU, ISTATUS );
            
	  /*~~~>   Check convergence of Newton iterations */
            NewtonIncrement = SDIRK_ErrorNorm(N, DU, SCAL);
            if ( NewtonIter == 0 ) {
	      Theta = ABS(ThetaMin);
	      NewtonRate = (KPP_REAL)2.0;
	    }
            else {
	      Theta = NewtonIncrement/NewtonIncrementOld;
	      if (Theta < (KPP_REAL)0.99) {
		NewtonRate = Theta/(ONE-Theta);
	      /* Predict error at the end of Newton process */ 
		NewtonPredictedErr = NewtonIncrement*
		  pow(Theta,(NewtonMaxit-NewtonIter)) / (ONE-Theta);
	      /* Non-convergence of Newton: predicted error too large */
                if (NewtonPredictedErr >= NewtonTol)
		  break; /* Exit NewtonLoop */
	      } else { /* Non-convergence of Newton: Theta too large */
		break;
	      }
            }
            NewtonIncrementOld = NewtonIncrement;
	  /* Update solution */
	    for(i=0; i<N; i++)
	      U[istage][iadj][i] -= DU[i];
            
          /* Check error in Newton iterations */
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol);
          /* AbsTol is often inappropriate for adjoints -
              we do at least 4 Newton iterations to ensure convergence
              of all adjoint components */
            if ((NewtonIter>=4) && NewtonDone)
	      break; /* Exit Newton Loop */ 
	  }
            
	/*~~~> If Newton iterations fail employ the direct solution */
	  if (!NewtonDone) {
	    printf("Problems with Newton Adjoint!!!\n");
#ifdef FULL_ALGEBRA
	    for(i=0; i<N; i++)
	      E_adj[i][j] = -Jac[i][j];
	    for(j=0; j<N; j++)
	      E_adj[j][j] += ONE/(H*rkGamma);
	    DGETRF( N, N, E_adj, N, IP_adj, IER );
#else
	    for(i=0; i<LU_NONZERO; i++)
	      E_adj[i] = -Jac[i];
	    for(i=0; i<NVAR; i++) {
	      j = LU_DIAG[i];
	      E_adj[j] += ONE/(H*rkGamma);
	    }
	    IER = KppDecomp ( E_adj);
#endif
	    ISTATUS[Ndec]++;
	    if (IER != 0) {
	      printf("At stage %d the matrix used in adjoint computation is "
		     "singular", istage);
	      return SDIRK_ErrorMsg(-8,T,H);
	    }
	    SDIRK_Solve ( 'T', H, N, E_adj, IP_adj, IER, G, ISTATUS );
	    for(i=0; i<N; i++)
	      U[istage][iadj][i] = G[i];  
	  }

	/*~~~>  End of implified Newton iterations */
	} /* End of DirADJ */

      } /* End of adj */
 
    } /* End of stages */

    /*~~~> Update adjoint solution
           Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:) */
    for (istage=0; istage<rkS; istage++) {
      for (iadj=0; iadj<NADJ; iadj++) {
	for(i=0; i<N; i++)
	  Lambda[iadj][i] += U[istage][iadj][i];
      }
    }
  } /* End of Tloop */

  /* Successful return */
  return 1;
  
} /* End of SDIRK_DadjInt */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void SDIRK_AllocBuffers(int Max_no_steps, int rkS, int SaveLU) {
/*~~~>  Allocate buffer space for checkpointing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int i,j;

  chk_H = (KPP_REAL*) malloc(Max_no_steps * sizeof(KPP_REAL));
  if (chk_H == NULL) {
    printf("Failed allocation of buffer H\n");
    exit(0);
  }

  chk_T = (KPP_REAL*) malloc(Max_no_steps * sizeof(KPP_REAL));
  if (chk_T == NULL) {
    printf("Failed allocation of buffer T\n");
    exit(0);
  }

  chk_Y = (KPP_REAL**) malloc(Max_no_steps * sizeof(KPP_REAL*));
  if (chk_Y == NULL) {
    printf("Failed allocation of buffer Y\n");
    exit(0);
  }
  for(i=0; i<Max_no_steps; i++) {
    chk_Y[i] = (KPP_REAL*) malloc(NVAR * sizeof(KPP_REAL));
    if (chk_Y[i] == NULL) {
      printf("Failed allocation of buffer Y\n");
      exit(0);
    }
  }
   
  chk_Z = (KPP_REAL***) malloc(Max_no_steps * sizeof(KPP_REAL**));
  if (chk_Z == NULL) {
    printf("Failed allocation of buffer Z\n");
    exit(0);
  }
  for(i=0; i<Max_no_steps; i++) {
    chk_Z[i] = (KPP_REAL**) malloc(rkS * sizeof(KPP_REAL*));
    if (chk_Z[i] == NULL) {
      printf("Failed allocation of buffer Z\n");
      exit(0);
    }
    for(j=0; j<rkS; j++) {
      chk_Z[i][j] = (KPP_REAL*) malloc(NVAR * sizeof(KPP_REAL));
      if (chk_Z[i][j] == NULL) {
	printf("Failed allocation of buffer Z\n");
	exit(0);
      }
    }
  }

  if (SaveLU) {
#ifdef FULL_ALGEBRA
    chk_J = (KPP_REAL***) malloc(Max_no_steps * sizeof(KPP_REAL**));
    if (chk_J == NULL) {
      printf("Failed allocation of buffer J\n");
      exit(0);
    }
    for(i=0; i<Max_no_steps; i++) {
      chk_J[i] = (KPP_REAL**) malloc(NVAR * sizeof(KPP_REAL*));
      if (chk_J[i] == NULL) {
	printf("Failed allocation of buffer J\n");
	exit(0);
      }
      for(j=0; j<NVAR; j++) {
	chk_J[i][j] = (KPP_REAL*) malloc(NVAR * sizeof(KPP_REAL));
	if (chk_J[i][j] == NULL) {
	  printf("Failed allocation of buffer J\n");
	  exit(0);
	}
      }
    }
#else
    chk_J = (KPP_REAL**) malloc(Max_no_steps * sizeof(KPP_REAL*));
    if (chk_J == NULL) {
      printf("Failed allocation of buffer J\n");
      exit(0);
    }
    for(i=0; i<Max_no_steps; i++) {
      chk_J[i] = (KPP_REAL*) malloc(LU_NONZERO * sizeof(KPP_REAL));
      if (chk_J[i] == NULL) {
	printf("Failed allocation of buffer J\n");
	exit(0);
      }
    }
#endif
 
    chk_P = (int**) malloc(Max_no_steps * sizeof(int*));
    if (chk_P == NULL) {
      printf("Failed allocation of buffer P\n");
      exit(0);
    }
    for(i=0; i<Max_no_steps; i++) {
      chk_P[i] = (int*) malloc(NVAR * sizeof(int));
      if (chk_P[i] == NULL) {
	printf("Failed allocation of buffer P\n");
	exit(0);
      }
    }
  }
} /* End of SDIRK_AllocBuffers */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void SDIRK_FreeBuffers(int Max_no_steps, int SaveLU) {
/*~~~>  Deallocate buffer space for discrete adjoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  int i,j;

  free(chk_H);
  free(chk_T);

  for(i=0; i<Max_no_steps; i++)
    free(chk_Y[i]);
  free(chk_Y);

  for(i=0; i<Max_no_steps; i++) {
    for(j=0; j<rkS; j++) {
      free(chk_Z[i][j]);
    }
    free(chk_Z[i]);
  }
  free(chk_Z);

  if(SaveLU) {
#ifdef FULL_ALGEBRA
    for(i=0; i<Max_no_steps; i++) {
      for(j=0; j<rkS; j++) {
	free(chk_J[i][j]);
      }
      free(chk_J[i]);
    }
    free(chk_Z);
#else
    for(i=0; i<Max_no_steps; i++)
      free(chk_J[i]);
    free(chk_J);
#endif

    for(i=0; i<Max_no_steps; i++)
      free(chk_P[i]);
    free(chk_P);
  }
} /* End of SDIRK_FreeBuffers */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void SDIRK_Push( KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL Z[][NVAR], 
		 KPP_REAL E[], int P[], int Max_no_steps, int SaveLU ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> Saves the next trajectory snapshot for discrete adjoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  int i,j;

  stack_ptr++;
  if( stack_ptr > Max_no_steps ) {
    printf( "Push failed: buffer overflow");
    exit(0);
  }

  chk_H[ stack_ptr ] = H;
  chk_T[ stack_ptr ] = T;
  for(i=0; i<NVAR; i++) {
    chk_Y[stack_ptr][i] = Y[i];
    for(j=0; j<rkS; j++)
      chk_Z[stack_ptr][j][i] = Z[j][i];
  }

  if (SaveLU) {
#ifdef FULL_ALGEBRA
    for(i=0; i<NVAR; i++) {
      for(j=0; j<NVAR; j++)
	chk_J[stack_ptr][i][j] = E[i][j];
      chk_P[stack_ptr][i] = P[i];
    }
#else   
    for(i=0; i<LU_NONZERO; i++)
      chk_J[stack_ptr][i] = E[i];
#endif
  }
  
} /* End of SDIRK_Push */
  
   
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void SDIRK_Pop( KPP_REAL* T, KPP_REAL* H, KPP_REAL* Y, KPP_REAL* Z, KPP_REAL* E,
		int* P, int SaveLU ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> Retrieves the next trajectory snapshot for discrete adjoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
  int i,j;

  if ( stack_ptr < 0 ) {
    printf( "Pop failed: empty buffer\n" );
    exit(0);
  }

  *H = chk_H[ stack_ptr ];
  *T = chk_T[ stack_ptr ];
  for(i=0; i<NVAR; i++) {
    Y[i] = chk_Y[stack_ptr][i];
    for(j=0; j<rkS; j++)
      Z[(j * NVAR) + i] = chk_Z[stack_ptr][j][i];
  }

  if (SaveLU) {
#ifdef FULL_ALGEBRA
    for(i=0; i<NVAR; i++) {
      for(j=0; j<NVAR; j++)
	E[(j*NVAR)+i] = chk_J[stack_ptr][j][i];
      P[i] = chk_P[stack_ptr][i];
    }
#else
    for(i=0; i<LU_NONZERO; i++)
      E[i] = chk_J[stack_ptr][i];
#endif
  }

  stack_ptr--;
  
} /* End of SDIRK_Pop */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void SDIRK_ErrorScale( int N, int ITOL, KPP_REAL AbsTol[], KPP_REAL RelTol[],
		       KPP_REAL Y[], KPP_REAL SCAL[]) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  int i;
  if (ITOL == 0){
    for (i = 0; i < NVAR; i++)
      SCAL[i] = ONE / (AbsTol[0]+RelTol[0]*ABS(Y[i]) );
  }
  else {
    for (i = 0; i < NVAR; i++)
      SCAL[i] = ONE / (AbsTol[i]+RelTol[i]*ABS(Y[i]) );
  }

} /*  End of SDIRK_ErrorScale */
      
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
KPP_REAL SDIRK_ErrorNorm( int N, KPP_REAL Y[], KPP_REAL SCAL[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  int i;
  KPP_REAL Err = ZERO;

  for (i = 0; i < N; i++)
    Err = Err + pow( (Y[i]*SCAL[i]), 2);
  Err = MAX( SQRT(Err/(KPP_REAL)N), (KPP_REAL)1.0e-10);

  return Err;

} /*  End of SDIRK_ErrorNorm  */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int SDIRK_ErrorMsg(int code, KPP_REAL T, KPP_REAL H) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Handles all error messages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

  printf("\nForced exit from Sdirk due to the following error:\n");

  switch (code) {
    case -1:
      printf("--> Improper value for maximal no of steps\n");
      break;
    case -2:
      printf("--> Selected Rosenbrock method not implemented\n");
      break;
    case -3:
      printf("--> Hmin/Hmax/Hstart must be positive\n");
      break;
    case -4:
      printf("--> FacMin/FacMax/FacRej must be positive\n");
      break;
    case -5:
      printf("--> Improper tolerance values\n");
      break;
    case -6:
      printf("--> No of steps exceeds maximum bound\n");
      break;
    case -7:
      printf("--> Step size too small T + 10*H = T or H < Roundoff\n");
      break;
    case -8:
      printf("--> Matrix is repeatedly singular\n");
      break;
    default: /* causing an error */
      printf("Unknown Error code: %d\n", code);
  }

  printf("\nTime = %f and H = %f\n", T, H );
  return code;

} /*  end SDIRK_ErrorMsg   */
      
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void SDIRK_PrepareMatrix ( KPP_REAL H, KPP_REAL T, KPP_REAL Y[], 
			   KPP_REAL FJAC[], int SkipJac, int SkipLU, 
			   KPP_REAL E[], int IP[], int Reject, int ISING, 
			   int ISTATUS[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~>  Compute the matrix E = 1/(H*GAMMA)*Jac, and its decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  KPP_REAL HGammaInv;
  int i, j, ConsecutiveSng = 0;

  ISING = 1;

  while (ISING != 0) {
    HGammaInv = ONE/(H*rkGamma);

  /*~~~>  Compute the Jacobian */
    if (SkipJac==0) {
      JAC_CHEM(T,Y,FJAC);
      ISTATUS[Njac]++;
    }

#ifdef FULL_ALGEBRA
    for(j=0; j<NVAR; j++) {
      for(i=0; i<NVAR; i++)
	E[j][i] = -FJAC[j][i];
      E[j][j] = E[j][j] + HGammaInv;
    }
    DGETRF(NVAR, NVAR, E, NVAR, IP, ISING);
#else
    for(i=0; i<LU_NONZERO; i++)
      E[i] = -FJAC[i];
    for(i=0; i<NVAR; i++) {
      j = LU_DIAG[i];
      E[j]=E[j] + HGammaInv;
    }

    ISING = KppDecomp(E);
    IP[0] = 1;
#endif
      
    ISTATUS[Ndec]++;

    if (ISING != 0) {
      printf("MATRIX IS SINGULAR, ISING=%d T=%e H=%e\n", ISING, T, H);
      ISTATUS[Nsng]++;
      ConsecutiveSng++;
      if (ConsecutiveSng >= 6) 
	return; /* Failure */
      H = (KPP_REAL)(0.5)*H;
      SkipJac = 0; /* False */
      SkipLU  = 0; /* False */
      Reject  = 1; /* True */
    }
  }
} /* End of SDIRK_PrepareMatrix */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void SDIRK_Solve ( char Transp, KPP_REAL H, int N, KPP_REAL E[], int IP[], 
		   int ISING, KPP_REAL RHS[], int ISTATUS[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~>  Solves the system (H*Gamma-Jac)*x = R
      using the LU decomposition of E = I - 1/(H*Gamma)*Jac
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  KPP_REAL HGammaInv;
      
  HGammaInv = ONE/(H*rkGamma);
  WSCAL(N,HGammaInv,RHS,1);
  switch (Transp) {
    case 'N':
#ifdef FULL_ALGEBRA  
      DGETRS( 'N', N, 1, E, N, IP, RHS, N, ISING );
#else
      KppSolve(E, RHS);
#endif
      break;
    case 'T':
#ifdef FULL_ALGEBRA  
      DGETRS( 'T', N, 1, E, N, IP, RHS, N, ISING );
#else
      KppSolveTR(E, RHS, RHS);
#endif
      break;
    default:
      printf("Error in SDIRK_Solve. Unknown Transp argument: %c\n", Transp);
      exit(0);
  }
  ISTATUS[Nsol]++;
} /* End of SDIRK_Solve */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Sdirk4a()
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
  sdMethod = S4A;

/* Number of stages */
  rkS = 5;

/* Method Coefficients */
  rkGamma = (KPP_REAL)0.2666666666666666666666666666666667;

  rkA[0][0] = (KPP_REAL)0.2666666666666666666666666666666667;
  rkA[0][1] = (KPP_REAL)0.5000000000000000000000000000000000;
  rkA[1][1] = (KPP_REAL)0.2666666666666666666666666666666667;
  rkA[0][2] = (KPP_REAL)0.3541539528432732316227461858529820;
  rkA[1][2] = (KPP_REAL)(-0.5415395284327323162274618585298197e-01);
  rkA[2][2] = (KPP_REAL)0.2666666666666666666666666666666667;
  rkA[0][3] = (KPP_REAL)0.8515494131138652076337791881433756e-01;
  rkA[1][3] = (KPP_REAL)(-0.6484332287891555171683963466229754e-01);
  rkA[2][3] = (KPP_REAL)0.7915325296404206392428857585141242e-01;
  rkA[3][3] = (KPP_REAL)0.2666666666666666666666666666666667;
  rkA[0][4] = (KPP_REAL)2.100115700566932777970612055999074;
  rkA[1][4] = (KPP_REAL)(-0.7677800284445976813343102185062276);
  rkA[2][4] = (KPP_REAL)2.399816361080026398094746205273880;
  rkA[3][4] = (KPP_REAL)(-2.998818699869028161397714709433394);
  rkA[4][4] = (KPP_REAL)0.2666666666666666666666666666666667;
  rkB[0] = (KPP_REAL)2.100115700566932777970612055999074;
  rkB[1] = (KPP_REAL)(-0.7677800284445976813343102185062276);
  rkB[2] = (KPP_REAL)2.399816361080026398094746205273880;
  rkB[3] = (KPP_REAL)(-2.998818699869028161397714709433394);
  rkB[4] = (KPP_REAL)0.2666666666666666666666666666666667;

  rkBhat[0] = (KPP_REAL)2.885264204387193942183851612883390;
  rkBhat[1] = (KPP_REAL)(-0.1458793482962771337341223443218041);
  rkBhat[2] = (KPP_REAL)2.390008682465139866479830743628554;
  rkBhat[3] = (KPP_REAL)(-4.129393538556056674929560012190140);
  rkBhat[4] = ZERO;

  rkC[0] = (KPP_REAL)0.2666666666666666666666666666666667;
  rkC[1] = (KPP_REAL)0.7666666666666666666666666666666667;
  rkC[2] = (KPP_REAL)0.5666666666666666666666666666666667;
  rkC[3] = (KPP_REAL)0.3661315380631796996374935266701191;
  rkC[4] = ONE;

/* Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i} */
  rkD[0] = ZERO;
  rkD[1] = ZERO;
  rkD[2] = ZERO;
  rkD[3] = ZERO;
  rkD[4] = ONE;

/* Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i} */
  rkE[0] = (KPP_REAL)(-0.6804000050475287124787034884002302);
  rkE[1] = (KPP_REAL)(1.558961944525217193393931795738823);
  rkE[2] = (KPP_REAL)(-13.55893003128907927748632408763868);
  rkE[3] = (KPP_REAL)(15.48522576958521253098585004571302);
  rkE[4] = ONE;

/* Local order of Err estimate */
  rkELO = 4;

/* h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j} */
  rkTheta[0][1] = (KPP_REAL)1.875000000000000000000000000000000;
  rkTheta[0][2] = (KPP_REAL)1.708847304091539528432732316227462;
  rkTheta[1][2] = (KPP_REAL)(-0.2030773231622746185852981969486824);
  rkTheta[0][3] = (KPP_REAL)0.2680325578937783958847157206823118;
  rkTheta[1][3] = (KPP_REAL)(-0.1828840955527181631794050728644549);
  rkTheta[2][3] = (KPP_REAL)0.2968246986151577397160821594427966;
  rkTheta[0][4] = (KPP_REAL)0.9096171815241460655379433581446771;
  rkTheta[1][4] = (KPP_REAL)(-3.108254967778352416114774430509465);
  rkTheta[2][4] = (KPP_REAL)12.33727431701306195581826123274001;
  rkTheta[3][4] = (KPP_REAL)(-11.24557012450885560524143016037523);

/* Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j} */
  rkAlpha[0][1] = (KPP_REAL)2.875000000000000000000000000000000;
  rkAlpha[0][2] = (KPP_REAL)0.8500000000000000000000000000000000;
  rkAlpha[1][2] = (KPP_REAL)0.4434782608695652173913043478260870;
  rkAlpha[0][3] = (KPP_REAL)0.7352046091658870564637910527807370;
  rkAlpha[1][3] = (KPP_REAL)(-0.9525565003057343527941920657462074e-01);
  rkAlpha[2][3] = (KPP_REAL)0.4290111305453813852259481840631738;
  rkAlpha[0][4] = (KPP_REAL)(-16.10898993405067684831655675112808);
  rkAlpha[1][4] = (KPP_REAL)6.559571569643355712998131800797873;
  rkAlpha[2][4] = (KPP_REAL)(-15.90772144271326504260996815012482);
  rkAlpha[3][4] = (KPP_REAL)25.34908987169226073668861694892683;

  rkELO = (KPP_REAL)4.0;

} /* end Sdirk4a */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Sdirk4b() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  sdMethod = S4B;

/* Number of stages */
  rkS = 5;

/* Method coefficients */
  rkGamma = (KPP_REAL)0.25;

  rkA[0][0] = (KPP_REAL)0.25;
  rkA[0][1] = (KPP_REAL)0.5;
  rkA[1][1] = (KPP_REAL)0.25;
  rkA[0][2] = (KPP_REAL)0.34;
  rkA[1][2] = (KPP_REAL)(-0.40e-01);
  rkA[2][2] = (KPP_REAL)0.25;
  rkA[0][3] = (KPP_REAL)0.2727941176470588235294117647058824;
  rkA[1][3] = (KPP_REAL)(-0.5036764705882352941176470588235294e-01);
  rkA[2][3] = (KPP_REAL)0.2757352941176470588235294117647059e-01;
  rkA[3][3] = (KPP_REAL)0.25;
  rkA[0][4] = (KPP_REAL)1.041666666666666666666666666666667;
  rkA[1][4] = (KPP_REAL)(-1.020833333333333333333333333333333);
  rkA[2][4] = (KPP_REAL)7.812500000000000000000000000000000;
  rkA[3][4] = (KPP_REAL)(-7.083333333333333333333333333333333);
  rkA[4][4] = (KPP_REAL)0.25;

  rkB[0] = (KPP_REAL)1.041666666666666666666666666666667;
  rkB[1] = (KPP_REAL)(-1.020833333333333333333333333333333);
  rkB[2] = (KPP_REAL)7.812500000000000000000000000000000;
  rkB[3] = (KPP_REAL)(-7.083333333333333333333333333333333);
  rkB[4] = (KPP_REAL)0.250000000000000000000000000000000;

  rkBhat[0] = (KPP_REAL)1.069791666666666666666666666666667;
  rkBhat[1] = (KPP_REAL)(-0.894270833333333333333333333333333);
  rkBhat[2] = (KPP_REAL)7.695312500000000000000000000000000;
  rkBhat[3] = (KPP_REAL)(-7.083333333333333333333333333333333);
  rkBhat[4] = (KPP_REAL)0.212500000000000000000000000000000;

  rkC[0] = (KPP_REAL)0.25;
  rkC[1] = (KPP_REAL)0.75;
  rkC[2] = (KPP_REAL)0.55;
  rkC[3] = (KPP_REAL)0.5;
  rkC[4] = ONE;

/* Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i} */
  rkD[0] = ZERO;
  rkD[1] = ZERO;
  rkD[2] = ZERO;
  rkD[3] = ZERO;
  rkD[4] = ONE;

/* Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i} */
  rkE[0] = (KPP_REAL)0.5750;
  rkE[1] = (KPP_REAL)0.2125;
  rkE[2] = (KPP_REAL)(-4.6875);
  rkE[3] = (KPP_REAL)4.2500;
  rkE[4] = (KPP_REAL)0.1500;

/* Local order of Err estimate */
  rkELO = 4;

/* h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j} */
  rkTheta[0][1] = (KPP_REAL)2.0;
  rkTheta[0][2] = (KPP_REAL)1.680000000000000000000000000000000;
  rkTheta[1][2] = (KPP_REAL)(-0.1600000000000000000000000000000000);
  rkTheta[0][3] = (KPP_REAL)1.308823529411764705882352941176471;
  rkTheta[1][3] = (KPP_REAL)(-0.1838235294117647058823529411764706);
  rkTheta[2][3] = (KPP_REAL)0.1102941176470588235294117647058824;
  rkTheta[0][4] = (KPP_REAL)(-3.083333333333333333333333333333333);
  rkTheta[1][4] = (KPP_REAL)(-4.291666666666666666666666666666667);
  rkTheta[2][4] = (KPP_REAL)34.37500000000000000000000000000000;
  rkTheta[3][4] = (KPP_REAL)(-28.3333333333333333333333333333);

/* Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j} */
  rkAlpha[0][1] = (KPP_REAL)3.0;
  rkAlpha[0][2] = (KPP_REAL)0.8800000000000000000000000000000000;
  rkAlpha[1][2] = (KPP_REAL)0.4400000000000000000000000000000000;
  rkAlpha[0][3] = (KPP_REAL)0.1666666666666666666666666666666667;
  rkAlpha[1][3] = (KPP_REAL)(-0.8333333333333333333333333333333333e-01);
  rkAlpha[2][3] = (KPP_REAL)0.9469696969696969696969696969696970;
  rkAlpha[0][4] = (KPP_REAL)(-6.0);
  rkAlpha[1][4] = (KPP_REAL)9.0;
  rkAlpha[2][4] = (KPP_REAL)(-56.81818181818181818181818181818182);
  rkAlpha[3][4] = (KPP_REAL)54.0;

} /* end Sdirk4b */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Sdirk2a() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  sdMethod = S2A;

/* ~~~> Number of stages */
  rkS = 2;

/* ~~~> Method coefficients */
  rkGamma = (KPP_REAL)0.2928932188134524755991556378951510;
  rkA[0][0] = (KPP_REAL)0.2928932188134524755991556378951510;
  rkA[0][1] = (KPP_REAL)0.7071067811865475244008443621048490;
  rkA[1][1] = (KPP_REAL)0.2928932188134524755991556378951510;
  rkB[0] = (KPP_REAL)0.7071067811865475244008443621048490;
  rkB[1] = (KPP_REAL)0.2928932188134524755991556378951510;
  rkBhat[0] = (KPP_REAL)0.6666666666666666666666666666666667;
  rkBhat[1] = (KPP_REAL)0.3333333333333333333333333333333333;
  rkC[0] = (KPP_REAL)0.292893218813452475599155637895151;
  rkC[1] = ONE;

/* ~~~> Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i} */
  rkD[0] = ZERO;
  rkD[1] = ONE;

/* ~~~> Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i} */
  rkE[0] = (KPP_REAL)0.4714045207910316829338962414032326;
  rkE[1] = (KPP_REAL)(-0.1380711874576983496005629080698993);

/* ~~~> Local order of Err estimate */
  rkELO = 2;

/* ~~~> h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j} */
  rkTheta[0][1] = (KPP_REAL)2.414213562373095048801688724209698;

/* ~~~> Starting value for Newton iterations */
  rkAlpha[0][1] = (KPP_REAL)3.414213562373095048801688724209698;

} /* end Sdirk2a */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Sdirk2b() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  sdMethod = S2B;

/* ~~~> Number of stages */
  rkS = 2;

/* ~~~> Method coefficients */
  rkGamma = (KPP_REAL)1.707106781186547524400844362104849;
  rkA[0][0] = (KPP_REAL)1.707106781186547524400844362104849;
  rkA[0][1] = (KPP_REAL)(-0.707106781186547524400844362104849);
  rkA[1][1] = (KPP_REAL)1.707106781186547524400844362104849;
  rkB[0] = (KPP_REAL)(-0.707106781186547524400844362104849);
  rkB[1] = (KPP_REAL)1.707106781186547524400844362104849;
  rkBhat[0] = (KPP_REAL)0.6666666666666666666666666666666667;
  rkBhat[1] = (KPP_REAL)0.3333333333333333333333333333333333;
  rkC[0] = (KPP_REAL)1.707106781186547524400844362104849;
  rkC[1] = ONE;

/* ~~~> Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i} */
  rkD[0] = ZERO;
  rkD[1] = ONE;

/* ~~~> Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i} */
  rkE[0] = (KPP_REAL)(-0.4714045207910316829338962414032326);
  rkE[1] = (KPP_REAL)0.8047378541243650162672295747365659;

/* ~~~> Local order of Err estimate */
  rkELO = 2;

/* ~~~> h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j} */
  rkTheta[0][1] = (KPP_REAL)(-0.414213562373095048801688724209698);

/* ~~~> Starting value for Newton iterations */
  rkAlpha[0][1] = (KPP_REAL)0.5857864376269049511983112757903019;

} /* end Sdirk2b */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Sdirk3a() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  sdMethod = S3A;

/* ~~~> Number of stages */
  rkS = 3;

/* ~~~> Method coefficients */
  rkGamma = (KPP_REAL)0.2113248654051871177454256097490213;
  rkA[0][0] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkA[0][1] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkA[1][1] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkA[0][2] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkA[1][2] = (KPP_REAL)0.5773502691896257645091487805019573;
  rkA[2][2] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkB[0] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkB[1] = (KPP_REAL)0.5773502691896257645091487805019573;
  rkB[2] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkBhat[0]= (KPP_REAL)0.2113248654051871177454256097490213;
  rkBhat[1]= (KPP_REAL)0.6477918909913548037576239837516312;
  rkBhat[2]= (KPP_REAL)0.1408832436034580784969504064993475;
  rkC[0] = (KPP_REAL)0.2113248654051871177454256097490213;
  rkC[1] = (KPP_REAL)0.4226497308103742354908512194980427;
  rkC[2] = ONE;

/* ~~~> Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i} */
  rkD[0] = ZERO;
  rkD[1] = ZERO;
  rkD[2] = ONE;

/* ~~~> Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i} */
  rkE[0] = (KPP_REAL)0.9106836025229590978424821138352906;
  rkE[1] = (KPP_REAL)(-1.244016935856292431175815447168624);
  rkE[2] = (KPP_REAL)0.3333333333333333333333333333333333;

/* ~~~> Local order of Err estimate    */
  rkELO    = 2;

/* ~~~> h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j} */
  rkTheta[0][1] =  ONE;
  rkTheta[0][2] = (KPP_REAL)(-1.732050807568877293527446341505872);
  rkTheta[1][2] = (KPP_REAL)2.732050807568877293527446341505872;

/* ~~~> Starting value for Newton iterations */
  rkAlpha[0][1] = (KPP_REAL)2.0;
  rkAlpha[0][2] = (KPP_REAL)(-12.92820323027550917410978536602349);
  rkAlpha[1][2] = (KPP_REAL)8.83012701892219323381861585376468;

} /* end Sdirk3a */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void FUN_CHEM(KPP_REAL T, KPP_REAL Y[], KPP_REAL P[])
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{

  KPP_REAL Told;

  Told = TIME;
  TIME = T;
  Update_SUN();
  Update_RCONST();
  Fun( Y, FIX, RCONST, P );
  TIME = Told;

}  /*  end FUN_CHEM  */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void JAC_CHEM(KPP_REAL T, KPP_REAL Y[], KPP_REAL JV[]) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  KPP_REAL Told;

#ifdef FULL_ALGEBRA
  KPP_REAL JS[LU_NONZERO];
  int i,j;
#endif

  Told = TIME;
  TIME = T;
  Update_SUN();
  Update_RCONST();

#ifdef FULL_ALGEBRA
  Jac_SP( Y, FIX, RCONST, JS);

  for(j=0; j<NVAR; j++) {
    for(i=0; i<NVAR; i++)
      JV[j][i] = (KPP_REAL)0.0;
  } /* end for */

  for(i=0; i<LU_NONZERO; i++)
    JV[LU_ICOL[i]][LU_IROW[i]] = JS[i];
#else
  Jac_SP(Y, FIX, RCONST, JV);
#endif

  TIME = Told;

} /* end JAC_CHEM  */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
