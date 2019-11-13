#define MAX(a,b) ( ((a) >= (b)) ? (a):(b)  )
#define MIN(b,c) ( ((b) < (c))  ? (b):(c)  )
#define ABS(x)	 ( ((x) >= 0 )  ? (x):(-x) )
#define SQRT(d)  ( pow((d),0.5) )

/* Numerical Constants */
#define ZERO	    (KPP_REAL)0.0
#define ONE	    (KPP_REAL)1.0
#define HALF        (KPP_REAL)0.5
#define DeltaMin    (KPP_REAL)1.0e-5
enum boolean { FALSE=0, TRUE=1 };

/*  Statistics on the work performed by the Rosenbrock method */
enum statistics { Nfun=1, Njac=2, Nstp=3, Nacc=4, Nrej=5, Ndec=6, Nsol=7, 
		  Nsng=8, Ntexit=1, Nhexit=2, Nhnew=3 };

/*~~~>  Parameters of the Rosenbrock method, up to 6 stages */
int ros_S, rosMethod;
enum ros_Params { RS2=1, RS3=2, RS4=3, RD3=4, RD4=5 };
KPP_REAL ros_A[15], ros_C[15], ros_M[6], ros_E[6], ros_Alpha[6], ros_Gamma[6], 
       ros_ELO;
int ros_NewF[6]; /* Holds Boolean values */
char ros_Name[12]; /* Length 12 */

/*~~~>  Types of Adjoints Implemented */
enum adjoint { Adj_none=1, Adj_discrete=2, Adj_continuous=3, 
	       Adj_simple_continuous=4 };

/*~~~>  Checkpoints in memory */
int bufsize = 200000;
int stack_ptr; /* last written entry */
KPP_REAL *chk_H, *chk_T;
KPP_REAL **chk_Y, **chk_K, **chk_J; /* 2D arrays */
KPP_REAL **chk_dY, **chk_d2Y; /* 2D arrays */

/* Function Headers */
void INTEGRATE_ADJ(int NADJ, KPP_REAL Y[], KPP_REAL Lambda[][NVAR], 
		   KPP_REAL TIN, KPP_REAL TOUT, KPP_REAL ATOL_adj[][NVAR], 
		   KPP_REAL RTOL_adj[][NVAR], int ICNTRL_U[], 
		   KPP_REAL RCNTRL_U[], int ISTATUS_U[], KPP_REAL RSTATUS_U[]);
int RosenbrockADJ( KPP_REAL Y[], int NADJ, KPP_REAL Lambda[][NVAR],
		   KPP_REAL Tstart, KPP_REAL Tend, KPP_REAL AbsTol[], 
		   KPP_REAL RelTol[], KPP_REAL AbsTol_adj[][NVAR], 
		   KPP_REAL RelTol_adj[][NVAR], KPP_REAL RCNTRL[], 
		   int ICNTRL[], KPP_REAL RSTATUS[], int ISTATUS[] );
void ros_AllocateDBuffers( int S, int SaveLU );
void ros_FreeDBuffers( int SaveLU );
void ros_AllocateCBuffers();
void ros_FreeCBuffers();
void ros_DPush( int S, KPP_REAL T, KPP_REAL H, KPP_REAL Ystage[], 
		KPP_REAL K[], KPP_REAL E[], int P[], int SaveLU );
void ros_DPop( int S, KPP_REAL* T, KPP_REAL* H, KPP_REAL* Ystage, 
	       KPP_REAL* K, KPP_REAL* E, int* P, int SaveLU );
void ros_CPush( KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL dY[], 
		KPP_REAL d2Y[] );
void ros_CPop( KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL dY[], 
	       KPP_REAL d2Y[] );
int ros_ErrorMsg( int Code, KPP_REAL T, KPP_REAL H);
int ros_FwdInt (KPP_REAL Y[], KPP_REAL Tstart, KPP_REAL Tend, KPP_REAL T, 
		KPP_REAL AbsTol[], KPP_REAL RelTol[], int AdjointType, 
		KPP_REAL Hmin, KPP_REAL Hstart, KPP_REAL Hmax, 
		KPP_REAL Roundoff, int ISTATUS[], int Max_no_steps, 
		KPP_REAL RSTATUS[], int Autonomous, int VectorTol, 
		KPP_REAL FacMax, KPP_REAL FacMin, KPP_REAL FacSafe, 
		KPP_REAL FacRej, int SaveLU);
int ros_DadjInt ( int NADJ, KPP_REAL Lambda[][NVAR], KPP_REAL Tstart, 
		  KPP_REAL Tend, KPP_REAL T, int SaveLU, int ISTATUS[], 
		  KPP_REAL Roundoff, int Autonomous);
int ros_CadjInt ( int NADJ, KPP_REAL Y[][NVAR], KPP_REAL Tstart, KPP_REAL Tend,
		  KPP_REAL T, KPP_REAL AbsTol_adj[][NVAR], 
		  KPP_REAL RelTol_adj[][NVAR], KPP_REAL RSTATUS[], 
		  KPP_REAL Hmin, KPP_REAL Hmax, KPP_REAL Hstart, 
		  KPP_REAL Roundoff, int Max_no_steps, int Autonomous, 
		  int VectorTol, KPP_REAL FacMax, KPP_REAL FacMin, 
		  KPP_REAL FacSafe, KPP_REAL FacRej, int ISTATUS[] );
int ros_SimpleCadjInt ( int NADJ, KPP_REAL Y[][NVAR], KPP_REAL Tstart, 
			KPP_REAL Tend, KPP_REAL T, int ISTATUS[], 
			int Autonomous,	KPP_REAL Roundoff );
KPP_REAL ros_ErrorNorm ( KPP_REAL Y[], KPP_REAL Ynew[], KPP_REAL Yerr[], 
		       KPP_REAL AbsTol[], KPP_REAL RelTol[], int VectorTol );
void ros_FunTimeDerivative ( KPP_REAL T, KPP_REAL Roundoff, KPP_REAL Y[], 
			     KPP_REAL Fcn0[], KPP_REAL dFdT[], int ISTATUS[] );
void ros_JacTimeDerivative ( KPP_REAL T, KPP_REAL Roundoff, KPP_REAL Y[], 
			     KPP_REAL Jac0[], KPP_REAL dJdT[], int ISTATUS[] );
int ros_PrepareMatrix ( KPP_REAL H, int Direction, KPP_REAL gam, 
			 KPP_REAL Jac0[], KPP_REAL Ghimj[], int Pivot[], 
			 int ISTATUS[] );
void ros_Decomp( KPP_REAL A[], int Pivot[], int* ising, int ISTATUS[] );
void ros_Solve( char How, KPP_REAL A[], int Pivot[], KPP_REAL b[], 
		int ISTATUS[] );
void ros_cadj_Y( KPP_REAL T, KPP_REAL Y[] );
void ros_Hermite3( KPP_REAL a, KPP_REAL b, KPP_REAL T, KPP_REAL Ya[], 
		   KPP_REAL Yb[], KPP_REAL Ja[], KPP_REAL Jb[], KPP_REAL Y[] );
void ros_Hermite5( KPP_REAL a, KPP_REAL b, KPP_REAL T, KPP_REAL Ya[], 
		   KPP_REAL Yb[], KPP_REAL Ja[], KPP_REAL Jb[], KPP_REAL Ha[], 
		   KPP_REAL Hb[], KPP_REAL Y[] );
void Ros2();
void Ros3();
void Ros4();
void Rodas3();
void Rodas4();
void JacTemplate( KPP_REAL T, KPP_REAL Y[], KPP_REAL Jcb[] );
void HessTemplate( KPP_REAL T, KPP_REAL Y[], KPP_REAL Hes[] );
void FunTemplate( KPP_REAL T, KPP_REAL Y[], KPP_REAL Fun [] ); 
void WSCAL( int N, KPP_REAL Alpha, KPP_REAL X[], int incX );
void WAXPY( int N, KPP_REAL Alpha, KPP_REAL X[], int incX, KPP_REAL Y[], 
	    int incY );
void WCOPY( int N, KPP_REAL X[], int incX, KPP_REAL Y[], int incY );
KPP_REAL WLAMCH( char C );
void Update_SUN();
void Update_RCONST();
void Fun( KPP_REAL Y[], KPP_REAL FIX[], KPP_REAL RCONST[], KPP_REAL Ydot[] );
void Jac_SP( KPP_REAL Y[], KPP_REAL FIX[], KPP_REAL RCONST[], KPP_REAL Ydot[]);
void Jac_SP_Vec( KPP_REAL Jac[], KPP_REAL Fcn[], KPP_REAL K[] );
void JacTR_SP_Vec( KPP_REAL Jac[], KPP_REAL Fcn[], KPP_REAL K[] );
void HessTR_Vec(KPP_REAL Hess[], KPP_REAL U1[], KPP_REAL U2[], KPP_REAL HTU[]);
void KppSolve( KPP_REAL A[], KPP_REAL b[] );
int KppDecomp( KPP_REAL A[] );
void KppSolveTR( KPP_REAL JVS[], KPP_REAL X[], KPP_REAL XX[] );
void Hessian( KPP_REAL V[], KPP_REAL F[], KPP_REAL RCT[], KPP_REAL Hess[] );

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void INTEGRATE_ADJ( int NADJ, KPP_REAL Y[], KPP_REAL Lambda[][NVAR], 
		    KPP_REAL TIN, KPP_REAL TOUT, KPP_REAL ATOL_adj[][NVAR], 
		    KPP_REAL RTOL_adj[][NVAR], int ICNTRL_U[], 
		    KPP_REAL RCNTRL_U[], int ISTATUS_U[], 
		    KPP_REAL RSTATUS_U[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*  Local Variables */
  KPP_REAL RCNTRL[20], RSTATUS[20];
  int ICNTRL[20], ISTATUS[20], IERR, i;

  for( i = 0; i < 20; i++ ) {
    ICNTRL[i]  = 0;
    RCNTRL[i]  = ZERO;
    ISTATUS[i] = 0;
    RSTATUS[i] = ZERO;
  }

/*~~~> fine-tune the integrator:
  ICNTRL(1) = 0       ! 0 = non-autonomous, 1 = autonomous
  ICNTRL(2) = 1       ! 0 = scalar, 1 = vector tolerances
  RCNTRL(3) = STEPMIN ! starting step
  ICNTRL(3) = 5       ! choice of the method for forward integration
  ICNTRL(6) = 1       ! choice of the method for continuous adjoint
  ICNTRL(7) = 2       ! 1=none, 2=discrete, 3=full continuous, 
                        4=simplified continuous adjoint
  ICNTRL(8) = 1       ! Save fwd LU factorization: 0=*don't* save, 1=save */

/* if optional parameters are given, and if they are >=0, then they overwrite 
   default settings */
//  if(ICNTRL_U != NULL) {
//    for(i=0; i<20; i++)
//	if(ICNTRL_U[i] > 0)
//	  ICNTRL[i] = ICNTRL_U[i];
//    } /* end for */
//  } /* end if */

//  if(RCNTRL_U != NULL) {
//    for(i=0; i<20; i++)
//	if(RCNTRL_U[i] > 0)
//	  RCNTRL[i] = RCNTRL_U[i];
//    } /* end for */
//  } /* end if */

  IERR = RosenbrockADJ( Y, NADJ, Lambda, TIN, TOUT, ATOL, RTOL, ATOL_adj, 
			RTOL_adj, RCNTRL, ICNTRL, RSTATUS, ISTATUS );

  if (IERR < 0)
    printf( "RosenbrockADJ: Unsucessful step at T=%f (IERR=%d)", TIN, IERR );

  STEPMIN = RSTATUS[Nhexit];

/* if optional parameters are given for output 
         copy to them to return information */
//  if(ISTATUS_U != NULL)
//    for(i=0; i<20; i++)
//	ISTATUS_U[i] = ISTATUS[i];
//  }

//  if(RSTATUS_U != NULL)
//    for(i=0; i<20; i++)
//	RSTATUS_U[i] = RSTATUS[i];
//  }

} /* End of INTEGRATE_ADJ */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int RosenbrockADJ( KPP_REAL Y[], int NADJ, KPP_REAL Lambda[][NVAR],
		   KPP_REAL Tstart, KPP_REAL Tend, KPP_REAL AbsTol[], 
		   KPP_REAL RelTol[], KPP_REAL AbsTol_adj[][NVAR],
		   KPP_REAL RelTol_adj[][NVAR], KPP_REAL RCNTRL[],
		   int ICNTRL[], KPP_REAL RSTATUS[], int ISTATUS[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ADJ = Adjoint of the Tangent Linear Model of a Rosenbrock Method

    Solves the system y'=F(t,y) using a RosenbrockADJ method defined by:

     G = 1/(H*gamma(1)) - Jac(t0,Y0)
     T_i = t0 + Alpha(i)*H
     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
         gamma(i)*dF/dT(t0, Y0)
     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j 

    For details on RosenbrockADJ methods and their implementation consult:
      E. Hairer and G. Wanner
      "Solving ODEs II. Stiff and differential-algebraic problems".
      Springer series in computational mathematics, Springer-Verlag, 1996.  
    The codes contained in the book inspired this implementation.

    (C)  Adrian Sandu, August 2004
    Virginia Polytechnic Institute and State University
    Contact: sandu@cs.vt.edu
    Revised by Philipp Miehe and Adrian Sandu, May 2006
    Translation F90 to C by Paul Eller, April 2007
    This implementation is part of KPP - the Kinetic PreProcessor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>   INPUT ARGUMENTS:

-     Y[NVAR]    = vector of initial conditions (at T=Tstart)
      NADJ       -> dimension of linearized system,
                   i.e. the number of sensitivity coefficients
-     Lambda[NVAR][NADJ] -> vector of initial sensitivity conditions 
                            (at T=Tstart)
-    [Tstart,Tend]  = time range of integration
     (if Tstart>Tend the integration is performed backwards in time)
-    RelTol, AbsTol = user precribed accuracy
- void Fun( T, Y, Ydot ) = ODE function,
                       returns Ydot = Y' = F(T,Y)
- void Jac( T, Y, Jcb ) = Jacobian of the ODE function,
                       returns Jcb = dF/dY
-    ICNTRL[0:9]    = integer inputs parameters
-    RCNTRL[0:9]    = real inputs parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>     OUTPUT ARGUMENTS:

-    Y[NVAR]    -> vector of final states (at T->Tend)
-    Lambda[NVAR][NADJ] -> vector of final sensitivities (at T=Tend)
-    ICNTRL[9:18]   -> integer output parameters
-    RCNTRL[9:18]   -> real output parameters
-    IERR       -> job status upon return
       - succes (positive value) or failure (negative value) -
           =  1 : Success
           = -1 : Improper value for maximal no of steps
           = -2 : Selected RosenbrockADJ method not implemented
           = -3 : Hmin/Hmax/Hstart must be positive
           = -4 : FacMin/FacMax/FacRej must be positive
           = -5 : Improper tolerance values
           = -6 : No of steps exceeds maximum bound
           = -7 : Step size too small
           = -8 : Matrix is repeatedly singular
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>     INPUT PARAMETERS:

    Note: For input parameters equal to zero the default values of the
       corresponding variables are used.

    ICNTRL[0]   = 1: F = F(y)   Independent of T (AUTONOMOUS)
              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)

    ICNTRL[1]   = 0: AbsTol, RelTol are NVAR-dimensional vectors
              = 1:  AbsTol, RelTol are scalars

    ICNTRL[2]  -> selection of a particular Rosenbrock method
        = 0 :  default method is Rodas3
        = 1 :  method is  Ros2
        = 2 :  method is  Ros3 
        = 3 :  method is  Ros4 
        = 4 :  method is  Rodas3
        = 5:   method is  Rodas4

    ICNTRL[3]  -> maximum number of integration steps
        For ICNTRL[4]=0) the default value of BUFSIZE is used

    ICNTRL[5]  -> selection of a particular Rosenbrock method for the
                continuous adjoint integration - for cts adjoint it
                can be different than the forward method ICNTRL(3)
         Note 1: to avoid interpolation errors (which can be huge!) 
                   it is recommended to use only ICNTRL[6] = 2 or 4
         Note 2: the performance of the full continuous adjoint
                   strongly depends on the forward solution accuracy Abs/RelTol

    ICNTRL[6] -> Type of adjoint algorithm
         = 0 : default is discrete adjoint ( of method ICNTRL[3] )
         = 1 : no adjoint
         = 2 : discrete adjoint ( of method ICNTRL[3] )
         = 3 : fully adaptive continuous adjoint ( with method ICNTRL[6] )
         = 4 : simplified continuous adjoint ( with method ICNTRL[6] )

    ICNTRL[7]  -> checkpointing the LU factorization at each step:
        ICNTRL[7]=0 : do *not* save LU factorization (the default)
        ICNTRL[7]=1 : save LU factorization
        Note: if ICNTRL[7]=1 the LU factorization is *not* saved

~~~>  Real input parameters:

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

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~>     OUTPUT PARAMETERS:

    Note: each call to RosenbrockADJ adds the corrent no. of fcn calls
      to previous value of ISTATUS(1), and similar for the other params.
      Set ISTATUS[0:9] = 0 before call to avoid this accumulation.

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
    RSTATUS[1]  -> Hexit, last accepted step before exit
    For multiple restarts, use Hexit as Hstart in the following run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~>  Local variables */
  KPP_REAL Roundoff, FacMin, FacMax, FacRej, FacSafe;
  KPP_REAL Hmin, Hmax, Hstart;
  KPP_REAL Texit=0.0;
  int i, UplimTol, Max_no_steps=0, IERR;
  int AdjointType=0, CadjMethod=0;
  int Autonomous, VectorTol, SaveLU; /* Holds boolean values */

  stack_ptr = -1;

/*~~~>  Initialize statistics */
  for(i=0; i<20; i++) {
    ISTATUS[i] = 0;
    RSTATUS[i] = ZERO;
  }

/*~~~>  Autonomous or time dependent ODE. Default is time dependent. */
  Autonomous = !(ICNTRL[0] == 0);

/*~~~> For Scalar tolerances(ICNTRL[1] != 0) the code uses AbsTol[0] and 
         RelTol[0]
   For Vector tolerances(ICNTRL[1] == 0) the code uses AbsTol[1:NVAR] and 
     RelTol[1:NVAR] */

  if (ICNTRL[1] == 0) {
    VectorTol = TRUE;
    UplimTol  = NVAR;
  }
  else { 
    VectorTol = FALSE;
    UplimTol  = 1;
  }

/*~~~>   Initialize the particular Rosenbrock method selected */
  switch( ICNTRL[2] ) {
    case 0:
    case 4:
      Rodas3();
      break;
    case 1:
      Ros2();
      break;
    case 2:
      Ros3();
      break;
    case 3:
      Ros4();
      break;
    case 5:
      Rodas4();
      break;
    default:
      printf( "Unknown Rosenbrock method: ICNTRL[2]=%d", ICNTRL[2] );
      return ros_ErrorMsg(-2, Tstart, ZERO);
  } /* End switch */

/*~~~>   The maximum number of steps admitted */
  if (ICNTRL[3] == 0)
    Max_no_steps = bufsize - 1;
  else if (Max_no_steps > 0)
    Max_no_steps = ICNTRL[3];
  else {
    printf("User-selected max no. of steps: ICNTRL[3]=%d",ICNTRL[3] );
    return ros_ErrorMsg(-1,Tstart,ZERO);
  }

/*~~~>The particular Rosenbrock method chosen for integrating the cts adjoint*/
  if (ICNTRL[5] == 0)
    CadjMethod = 4;
  else if ( (ICNTRL[5] >= 1) && (ICNTRL[5] <= 5) )
    CadjMethod = ICNTRL[5];
  else {
    printf( "Unknown CADJ Rosenbrock method: ICNTRL[5]=%d", CadjMethod );
    return ros_ErrorMsg(-2,Tstart,ZERO);
  }

/*~~~>  Discrete or continuous adjoint formulation */
  if ( ICNTRL[6] == 0 )
    AdjointType = Adj_discrete;
  else if ( (ICNTRL[6] >= 1) && (ICNTRL[6] <= 4) )
    AdjointType = ICNTRL[6];
  else {
    printf( "User-selected adjoint type: ICNTRL[6]=%d", AdjointType );
    return ros_ErrorMsg(-9,Tstart,ZERO);
  }

/*~~~> Save or not the forward LU factorization */
  SaveLU = (ICNTRL[7] != 0);

/*~~~>  Unit roundoff (1+Roundoff>1)  */
  Roundoff = WLAMCH('E');

/*~~~>  Lower bound on the step size: (positive value) */
  if (RCNTRL[0] == ZERO)
    Hmin = ZERO;
  else if (RCNTRL[0] > ZERO)
    Hmin = RCNTRL[0];
  else {
    printf( "User-selected Hmin: RCNTRL[0]=%f", RCNTRL[0] );
    return ros_ErrorMsg(-3,Tstart,ZERO);
  }

/*~~~>  Upper bound on the step size: (positive value) */
  if (RCNTRL[1] == ZERO)
    Hmax = ABS(Tend-Tstart);
  else if (RCNTRL[1] > ZERO)
    Hmax = MIN(ABS(RCNTRL[1]),ABS(Tend-Tstart));
  else {
    printf( "User-selected Hmax: RCNTRL[1]=%f", RCNTRL[1] );
    return ros_ErrorMsg(-3,Tstart,ZERO);
  }

/*~~~>  Starting step size: (positive value) */
  if (RCNTRL[2] == ZERO) {
    Hstart = MAX(Hmin,DeltaMin);
  }
  else if (RCNTRL[2] > ZERO)
    Hstart = MIN(ABS(RCNTRL[2]),ABS(Tend-Tstart));
  else {
    printf( "User-selected Hstart: RCNTRL[2]=%f", RCNTRL[2] );
    return ros_ErrorMsg(-3,Tstart,ZERO);
  }

/*~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax */
  if (RCNTRL[3] == ZERO)
    FacMin = (KPP_REAL)0.2;
  else if (RCNTRL[3] > ZERO)
    FacMin = RCNTRL[3];
  else {
    printf( "User-selected FacMin: RCNTRL[3]=%f", RCNTRL[3] );
    return ros_ErrorMsg(-4,Tstart,ZERO);
  }
  if (RCNTRL[4] == ZERO)
    FacMax = (KPP_REAL)6.0;
  else if (RCNTRL[4] > ZERO)
    FacMax = RCNTRL[4];
  else {
    printf( "User-selected FacMax: RCNTRL[4]=%f", RCNTRL[4] );
    return ros_ErrorMsg(-4,Tstart,ZERO);
  }

/*~~~>   FacRej: Factor to decrease step after 2 succesive rejections */
  if (RCNTRL[5] == ZERO)
    FacRej = (KPP_REAL)0.1;
  else if (RCNTRL[5] > ZERO)
    FacRej = RCNTRL[5];
  else {
    printf( "User-selected FacRej: RCNTRL[5]=%f", RCNTRL[5] );
    return ros_ErrorMsg(-4,Tstart,ZERO);
  }

/*~~~>   FacSafe: Safety Factor in the computation of new step size */
  if (RCNTRL[6] == ZERO)
    FacSafe = (KPP_REAL)0.9;
  else if (RCNTRL[6] > ZERO)
    FacSafe = RCNTRL[6];
  else {
    printf( "User-selected FacSafe: RCNTRL[6]=%f", RCNTRL[6] );
    return ros_ErrorMsg(-4,Tstart,ZERO);
  }

/*~~~>  Check if tolerances are reasonable */
  for(i=0; i < UplimTol; i++) {
    if ( (AbsTol[i] <= ZERO) || (RelTol[i] <= (KPP_REAL)10.0*Roundoff)
	 || (RelTol[i] >= (KPP_REAL)1.0) ) {
      printf( " AbsTol[%d] = %f", i, AbsTol[i] );
      printf( " RelTol[%d] = %f", i, RelTol[i] );
      return ros_ErrorMsg(-5,Tstart,ZERO);
    }
  }

/*~~~>  Allocate checkpoint space or open checkpoint files */
  if (AdjointType == Adj_discrete) {
    ros_AllocateDBuffers( ros_S, SaveLU );
  }
  else if ( (AdjointType == Adj_continuous) || 
	    (AdjointType == Adj_simple_continuous) ) {
    ros_AllocateCBuffers();
  }

/*~~~>  CALL Forward Rosenbrock method */
  IERR = ros_FwdInt(Y, Tstart, Tend, Texit, AbsTol, RelTol, AdjointType, Hmin,
		    Hstart, Hmax, Roundoff, ISTATUS, Max_no_steps, 
		    RSTATUS, Autonomous, VectorTol, FacMax, FacMin,
		    FacSafe, FacRej, SaveLU);

  printf( "\n\nFORWARD STATISTICS\n" );
  printf( "Step=%d Acc=%d Rej=%d Singular=%d\n\n", Nstp, Nacc, Nrej, Nsng );

/*~~~>  If Forward integration failed return */
  if (IERR<0)
    return IERR;

/*~~~>   Initialize the particular Rosenbrock method for continuous adjoint */
  if ( (AdjointType == Adj_continuous) || 
       (AdjointType == Adj_simple_continuous) ) {
    switch (CadjMethod) {
      case 1:
	Ros2();
	break;
      case 2:
	Ros3();
	break;
      case 3:
	Ros4();
	break;
      case 4:
	Rodas3();
	break;
      case 5:
	Rodas4();
	break;
      default:
	printf( "Unknown Rosenbrock method: ICNTRL[2]=%d", ICNTRL[2] );
	return ros_ErrorMsg(-2,Tstart,ZERO);
    }
  } /* End switch */

  switch( AdjointType ) {
    case Adj_discrete:
      IERR = ros_DadjInt (NADJ, Lambda, Tstart, Tend, Texit, SaveLU, ISTATUS, 
			  Roundoff, Autonomous );
      break;
    case Adj_continuous:
      IERR = ros_CadjInt (NADJ, Lambda, Tend, Tstart, Texit, AbsTol_adj, 
			  RelTol_adj, RSTATUS, Hmin, Hmax, Hstart, Roundoff,
			  Max_no_steps, Autonomous, VectorTol, FacMax, FacMin,
			  FacSafe, FacRej, ISTATUS);
      break;
    case Adj_simple_continuous:
      IERR = ros_SimpleCadjInt (NADJ, Lambda, Tstart, Tend, Texit, ISTATUS, 
				Autonomous, Roundoff);
  } /* End switch for AdjointType */

  printf( "ADJOINT STATISTICS\n" );
  printf( "Step=%d Acc=%d Rej=%d Singular=%d\n",Nstp,Nacc,Nrej,Nsng );

/*~~~>  Free checkpoint space or close checkpoint files */
  if (AdjointType == Adj_discrete)
    ros_FreeDBuffers( SaveLU );
  else if ( (AdjointType == Adj_continuous) || 
	    (AdjointType == Adj_simple_continuous) )
    ros_FreeCBuffers();

  return IERR;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_AllocateDBuffers( int S, int SaveLU ) {
/*~~~>  Allocate buffer space for discrete adjoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
  int i;
    
  chk_H = (KPP_REAL*) malloc(bufsize * sizeof(KPP_REAL));
  if (chk_H == NULL) {
    printf("Failed allocation of buffer H");
    exit(0);
  }
    
  chk_T = (KPP_REAL*) malloc(bufsize * sizeof(KPP_REAL));
  if (chk_T == NULL) {
    printf("Failed allocation of buffer T");
    exit(0);
  }
  
  chk_Y = (KPP_REAL**) malloc(bufsize * sizeof(KPP_REAL*));
  if (chk_Y == NULL) {
    printf("Failed allocation of buffer Y");
    exit(0);
  }
  for(i=0; i<bufsize; i++) {
    chk_Y[i] = (KPP_REAL*) malloc(NVAR * S * sizeof(KPP_REAL));
    if (chk_Y[i] == NULL) {
      printf("Failed allocation of buffer Y");
      exit(0);
    }
  }
    
  chk_K = (KPP_REAL**) malloc(bufsize * sizeof(KPP_REAL*));
  if (chk_K == NULL) {
    printf("Failed allocation of buffer K");
    exit(0);
  }
  for(i=0; i<bufsize; i++) {
    chk_K[i] = (KPP_REAL*) malloc(NVAR * S * sizeof(KPP_REAL));
    if (chk_K == NULL) {
      printf("Failed allocation of buffer K");
      exit(0);
    }
  }
  
  if (SaveLU) {
    chk_J = (KPP_REAL**) malloc(bufsize * sizeof(KPP_REAL*));
    if (chk_J == NULL) {
      printf( "Failed allocation of buffer J");
      exit(0);
    }
#ifdef FULL_ALGEBRA
    for(i=0; i<bufsize; i++) {
      chk_J[i] = (KPP_REAL*) malloc(NVAR * NVAR * sizeof(KPP_REAL));
      if (chk_J == NULL) {
	printf( "Failed allocation of buffer J");
	exit(0);
      }
    }
#else
    for(i=0; i<bufsize; i++) {
      chk_J[i] = (KPP_REAL*) malloc(LU_NONZERO * sizeof(KPP_REAL));
      if (chk_J == NULL) {
	printf( "Failed allocation of buffer J");
	exit(0);
      }
    }
#endif
  } 

} /* End of ros_AllocateDBuffers */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_FreeDBuffers( int SaveLU ) {
/*~~~>  Deallocate buffer space for discrete adjoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  int i;

  free(chk_H);
  free(chk_T);

  for(i=0; i<bufsize; i++) {
    free(chk_Y[i]);
    free(chk_K[i]);
  }
  free(chk_Y);
  free(chk_K);

  if (SaveLU) {
    for(i=0; i<bufsize; i++)
      free(chk_J[i]);
    free(chk_J);
  }

} /* End of ros_FreeDBuffers */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_AllocateCBuffers() {
/*~~~>  Allocate buffer space for continuous adjoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  int i;

  chk_H = (KPP_REAL*) malloc(bufsize * sizeof(KPP_REAL));
  if (chk_H == NULL) {
    printf( "Failed allocation of buffer H");
    exit(0);
  }

  chk_T = (KPP_REAL*) malloc(bufsize * sizeof(KPP_REAL));
  if (chk_T == NULL) {
    printf( "Failed allocation of buffer T");
    exit(0);
  }

  chk_Y = (KPP_REAL**) malloc(sizeof(KPP_REAL*) * bufsize);
  if (chk_Y == NULL) {
    printf( "Failed allocation of buffer Y");
    exit(0);
  }
  for(i=0; i<bufsize; i++) {
    chk_Y[i] = (KPP_REAL*) malloc( sizeof(KPP_REAL)* NVAR );
    if (chk_Y == NULL) {
      printf( "Failed allocation of buffer Y");
      exit(0);
    }
  }

  chk_dY = (KPP_REAL**) malloc(sizeof(KPP_REAL*) * bufsize);
  if (chk_dY == NULL) {
    printf( "Failed allocation of buffer dY");
    exit(0);
  }
  for(i=0; i<bufsize; i++) {
    chk_dY[i] = (KPP_REAL*) malloc( sizeof(KPP_REAL) * NVAR);
    if (chk_dY == NULL) {
      printf( "Failed allocation of buffer dY");
      exit(0);
    }
  }

  chk_d2Y = (KPP_REAL**) malloc(sizeof(KPP_REAL*) * bufsize);
  if (chk_d2Y == NULL) {
    printf( "Failed allocation of buffer d2Y");
    exit(0);
  }
  for(i=0; i<bufsize; i++) {
    chk_d2Y[i] = (KPP_REAL*) malloc( sizeof(KPP_REAL) * NVAR);
    if (chk_d2Y == NULL) {
      printf( "Failed allocation of buffer d2Y");
      exit(0);
    }
  }
} /* End of ros_AllocateCBuffers */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_FreeCBuffers() {
/*~~~>  Dallocate buffer space for continuous adjoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  int i;

  free(chk_H);
  free(chk_T);

  for(i=0; i<bufsize; i++) {
    free(chk_Y[i]);
    free(chk_dY[i]);
    free(chk_d2Y[i]);
  }

  free(chk_Y);
  free(chk_dY);
  free(chk_d2Y);

} /* End of ros_FreeCBuffers */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_DPush( int S, KPP_REAL T, KPP_REAL H, KPP_REAL Ystage[], 
		KPP_REAL K[], KPP_REAL E[], int P[], int SaveLU ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> Saves the next trajectory snapshot for discrete adjoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  int i;

  stack_ptr = stack_ptr + 1;
  if ( stack_ptr >= bufsize ) {
    printf( "Push failed: buffer overflow" );
    exit(0);
  }

  chk_H[ stack_ptr ] = H;
  chk_T[ stack_ptr ] = T;
  for(i=0; i<NVAR*S; i++) {
    chk_Y[stack_ptr][i] = Ystage[i];
    chk_K[stack_ptr][i] = K[i];
  }

  if (SaveLU) {
#ifdef FULL_ALGEBRA
    int j;
    for(j=0; j<NVAR; j++) {
      for(i=0; i<NVAR; i++)
	chk_J[stack_ptr][i][j] = E[i][j];
      chk_P[stack_ptr][j] = P[j];
    }
#else
    for(i=0; i<LU_NONZERO; i++)
      chk_J[stack_ptr][i] = E[i];
#endif
  }

} /* End of ros_DPush */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_DPop( int S, KPP_REAL* T, KPP_REAL* H, KPP_REAL* Ystage,
	       KPP_REAL* K, KPP_REAL* E, int* P, int SaveLU ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> Retrieves the next trajectory snapshot for discrete adjoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
  int i;
  
  if ( stack_ptr < 0 ) {
    printf( "Pop failed: empty buffer" );
    exit(0);
  }

  *H = chk_H[ stack_ptr ];
  *T = chk_T[ stack_ptr ];

  for(i=0; i<NVAR*S; i++) {
    Ystage[i] = chk_Y[stack_ptr][i];
    K[i]      = chk_K[stack_ptr][i];
  }

  if (SaveLU) {
#ifdef FULL_ALGEBRA
    int j;
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

} /* End of ros_DPop */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_CPush( KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL dY[], 
		KPP_REAL d2Y[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> Saves the next trajectory snapshot for discrete adjoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
  int i;
  
  stack_ptr++;
  if ( stack_ptr > bufsize ) {
    printf( "Push failed: buffer overflow" );
    exit(0);
  }
  chk_H[ stack_ptr ] = H;
  chk_T[ stack_ptr ] = T;

  for(i = 0; i< NVAR; i++ ) {
    chk_Y[stack_ptr][i] = Y[i];
    chk_dY[stack_ptr][i] = dY[i];
    chk_d2Y[stack_ptr][i] = d2Y[i];
  }
} /* End of ros_CPush */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_CPop( KPP_REAL T, KPP_REAL H, KPP_REAL Y[], KPP_REAL dY[], 
	       KPP_REAL d2Y[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> Retrieves the next trajectory snapshot for discrete adjoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
  int i;

  if ( stack_ptr <= 0 ) {
    printf( "Pop failed: empty buffer" );
    exit(0);
  }
  H = chk_H[ stack_ptr ];
  T = chk_T[ stack_ptr ];

  for(i=0; i<NVAR; i++) {
    Y[i] = chk_Y[stack_ptr][i];
    dY[i] = chk_dY[stack_ptr][i];
    d2Y[i] = chk_d2Y[stack_ptr][i];
  }
  stack_ptr--;
} /* End of ros_CPop */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int ros_ErrorMsg( int Code, KPP_REAL T, KPP_REAL H) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Handles all error messages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  int IERR = Code;
  printf( "Forced exit from RosenbrockADJ due to the following error:");
  
  switch (Code) {
    case -1:
      printf( "--> Improper value for maximal no of steps" );
      break;
    case -2:
      printf( "--> Selected RosenbrockADJ method not implemented" );
      break;
    case -3:
      printf( "--> Hmin/Hmax/Hstart must be positive" );
      break;
    case -4:
      printf( "--> FacMin/FacMax/FacRej must be positive" );
      break;
    case -5:
      printf( "--> Improper tolerance values" );
      break;
    case -6:
      printf( "--> No of steps exceeds maximum buffer bound" );
      break;
    case -7:
      printf( "--> Step size too small: T + 10*H = T or H < Roundoff" );
      break;
    case -8:
      printf( "--> Matrix is repeatedly singular" );
      break;
    case -9:
      printf( "--> Improper type of adjoint selected" );
      break;
    default:
      printf( "Unknown Error code: %d", Code );
  } /* End of switch */

  return IERR;
} /* End of ros_ErrorMsg */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int ros_FwdInt ( KPP_REAL Y[], KPP_REAL Tstart, KPP_REAL Tend, KPP_REAL T, 
		 KPP_REAL AbsTol[], KPP_REAL RelTol[], int AdjointType, 
		 KPP_REAL Hmin, KPP_REAL Hstart, KPP_REAL Hmax, 
		 KPP_REAL Roundoff, int ISTATUS[], int Max_no_steps, 
		 KPP_REAL RSTATUS[], int Autonomous, int VectorTol, 
		 KPP_REAL FacMax, KPP_REAL FacMin, KPP_REAL FacSafe, 
		 KPP_REAL FacRej, int SaveLU ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Template for the implementation of a generic RosenbrockADJ method
      defined by ros_S (no of stages)
      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
~~~> Y - Input: the initial condition at Tstart; Output: the solution at T
~~~> Tstart, Tend - Input: integration interval
~~~> T - Output: time at which the solution is returned (T=Tend if success)
~~~> AbsTol, RelTol - Input: tolerances
~~~> IERR - Output: Error indicator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* ~~~~ Local variables */
  KPP_REAL Ynew[NVAR], Fcn0[NVAR], Fcn[NVAR];
  KPP_REAL K[NVAR*ros_S], dFdT[NVAR];
  KPP_REAL *Ystage = NULL; /* Array pointer */
#ifdef FULL_ALGEBRA
  KPP_REAL Jac0[NVAR][NVAR],  Ghimj[NVAR][NVAR];
#else
  KPP_REAL Jac0[LU_NONZERO], Ghimj[LU_NONZERO];
#endif
  KPP_REAL H, Hnew, HC, HG, Fac, Tau;
  KPP_REAL Err, Yerr[NVAR];
  int Pivot[NVAR], Direction, ioffset, i, j=0, istage;
  int RejectLastH, RejectMoreH, Singular; /* Boolean Values */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~>  Allocate stage vector buffer if needed */
  if (AdjointType == Adj_discrete) {
    Ystage = (KPP_REAL*) malloc(NVAR*ros_S*sizeof(KPP_REAL));
    /* Uninitialized Ystage may lead to NaN on some compilers */
    if (Ystage == NULL) {
      printf( "Allocation of Ystage failed" );
      exit(0);
    }
  }

/*~~~>  Initial preparations */
  T = Tstart;
  RSTATUS[Nhexit] = ZERO;
  H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) );
  if (ABS(H) <= ((KPP_REAL)10.0)*Roundoff)
    H = DeltaMin;

  if (Tend >= Tstart)
    Direction = 1;
  else
    Direction = -1;

  H = Direction*H;
  
  RejectLastH = FALSE;
  RejectMoreH = FALSE;

/*~~~> Time loop begins below */
  while ( ((Direction > 0) && ((T-Tend)+Roundoff <= ZERO)) || 
	  ((Direction < 0) && ((Tend-T)+Roundoff <= ZERO)) ) { /* TimeLoop */

    if ( ISTATUS[Nstp] > Max_no_steps )  /* Too many steps */
      return ros_ErrorMsg(-6,T,H);

    if ( ((T+((KPP_REAL)0.1)*H) == T) || (H <= Roundoff) ) /* Step size 
							    too small */
      return ros_ErrorMsg(-7,T,H);

/*~~~>  Limit H if necessary to avoid going beyond Tend */
    RSTATUS[Nhexit] = H;
    H = MIN(H,ABS(Tend-T));

/*~~~>   Compute the function at current time */
    FunTemplate(T,Y,Fcn0);
    ISTATUS[Nfun] = ISTATUS[Nfun] + 1;

/*~~~>  Compute the function derivative with respect to T */
    if (!Autonomous)
      ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT, ISTATUS );

/*~~~>   Compute the Jacobian at current time */
    JacTemplate(T,Y,Jac0);
    ISTATUS[Njac] = ISTATUS[Njac] + 1;

/*~~~>  Repeat step calculation until current step accepted */
    do {  /* UntilAccepted */

      Singular = ros_PrepareMatrix ( H,Direction,ros_Gamma[0],Jac0,Ghimj,Pivot,
				     ISTATUS );

      if (Singular) /* More than 5 consecutive failed decompositions */
	return ros_ErrorMsg(-8,T,H);

/*~~~>   Compute the stages */
      for( istage = 0; istage < ros_S; istage++ ) { /* Stage */

	/* Current istage offset. Current istage vector is 
	   K(ioffset+1:ioffset+NVAR) */
	ioffset = NVAR*istage;

	/*For the 1st istage the function has been computed previously*/
	if ( istage == 0 ) {
	  WCOPY(NVAR,Fcn0,1,Fcn,1);
	  if (AdjointType == Adj_discrete) { /* Save stage solution */
	    for(i=0; i<NVAR; i++)
	      Ystage[i] = Y[i];
	    WCOPY(NVAR,Y,1,Ynew,1);
	  }
	}
	/* istage>0 and a new function evaluation is needed at the 
	   current istage */
	else if ( ros_NewF[istage] ) {
	  WCOPY(NVAR,Y,1,Ynew,1);
	  for ( j = 0; j < istage; j++ ) {
	    WAXPY( NVAR,ros_A[(istage)*(istage-1)/2+j], 
		   &K[NVAR*j],1,Ynew,1 );
	  }
	  Tau = T + ros_Alpha[istage]*Direction*H;
	  FunTemplate(Tau,Ynew,Fcn);
	  ISTATUS[Nfun] = ISTATUS[Nfun] + 1;
	} /* if istage == 1 elseif ros_NewF[istage] */

	/* Save stage solution every time even if ynew is not updated */
	if ( ( istage > 0 ) && (AdjointType == Adj_discrete) ) {
	  for(i=0; i<NVAR; i++)
	    Ystage[ioffset+i] = Ynew[i];
	}
	WCOPY(NVAR,Fcn,1,&K[ioffset],1);
	for( j = 0; j < istage; j++ ) {
	  HC = ros_C[(istage)*(istage-1)/2+j]/(Direction*H);
	  WAXPY(NVAR,HC,&K[NVAR*j],1,&K[ioffset],1);
	}
	if (( !Autonomous) && (ros_Gamma[istage] != ZERO)) {
	  HG = Direction*H*ros_Gamma[istage];
	  WAXPY(NVAR,HG,dFdT,1,&K[ioffset],1);
	}
	ros_Solve('N', Ghimj, Pivot, &K[ioffset], ISTATUS);
      } /* End of Stage loop */

/*~~~>  Compute the new solution */
      WCOPY(NVAR,Y,1,Ynew,1);
      for( j=0; j<ros_S; j++ )
	WAXPY(NVAR,ros_M[j],&K[NVAR*j],1,Ynew,1);

/*~~~>  Compute the error estimation */ 
      WSCAL(NVAR,ZERO,Yerr,1);
      for( j=0; j<ros_S; j++ )
	WAXPY(NVAR,ros_E[j],&K[NVAR*j],1,Yerr,1);
      Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol );

/*~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax */
      Fac  = MIN(FacMax,MAX(FacMin,FacSafe/pow(Err,(ONE/ros_ELO))));
      Hnew = H*Fac;

/*~~~>  Check the error magnitude and adjust step size */
      ISTATUS[Nstp] = ISTATUS[Nstp] + 1;
      if ( (Err <= ONE) || (H <= Hmin) ) {  /*~~~> Accept step */
	ISTATUS[Nacc]++;
	if (AdjointType == Adj_discrete) { /* Save current state */
	  ros_DPush( ros_S, T, H, Ystage, K, Ghimj, Pivot, SaveLU );
	}
	else if ( (AdjointType == Adj_continuous) || 
		  (AdjointType == Adj_simple_continuous) ) {
#ifdef FULL_ALGEBRA
	  K = MATMUL(Jac0,Fcn0);
#else
	  Jac_SP_Vec( Jac0, Fcn0, &K[0] );
#endif
	  if ( !Autonomous)
	    WAXPY(NVAR,ONE,dFdT,1,&K[0],1);
	  ros_CPush( T, H, Y, Fcn0, &K[0] );
	}
	WCOPY(NVAR,Ynew,1,Y,1);
	T = T + Direction*H;
	Hnew = MAX(Hmin,MIN(Hnew,Hmax));
	if (RejectLastH) { /* No step size increase after a 
			     rejected step */
	  Hnew = MIN(Hnew,H);
	}
	RSTATUS[Nhexit] = H;
	RSTATUS[Nhnew]  = Hnew;
	RSTATUS[Ntexit] = T;
	RejectLastH = FALSE;
	RejectMoreH = FALSE;
	H = Hnew;
	break; /* UntilAccepted - EXIT THE LOOP: WHILE STEP NOT ACCEPTED */
      }

      else { /*~~~> Reject step */
	if (RejectMoreH)
	  Hnew = H*FacRej;
        RejectMoreH = RejectLastH;
        RejectLastH = TRUE;
        H = Hnew;
        if (ISTATUS[Nacc] >= 1)
	  ISTATUS[Nrej]++;
      } /* End if else - Err <= 1 */
    
    } while(1); /* End of UntilAccepted do loop */
  } /* End of TimeLoop */
   
/*~~~> Save last state: only needed for continuous adjoint */
  if ( (AdjointType == Adj_continuous) || 
       (AdjointType == Adj_simple_continuous) ) {
    FunTemplate(T,Y,Fcn0);
    ISTATUS[Nfun]++;
    JacTemplate(T,Y,Jac0);
    ISTATUS[Njac]++;
#ifdef FULL_ALGEBRA
    K = MATMUL(Jac0,Fcn0);
#else
    Jac_SP_Vec( Jac0, Fcn0, &K[0] );
#endif
    if (!Autonomous) {
      ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT, ISTATUS );
      WAXPY(NVAR,ONE,dFdT,1,&K[0],1);
    }
    ros_CPush( T, H, Y, Fcn0, &K[0] );
/*~~~> Deallocate stage buffer: only needed for discrete adjoint */
  }
  else if (AdjointType == Adj_discrete) {
    free(Ystage);
  }

/*~~~> Succesful exit */
  return 1;  /*~~~> The integration was successful */
} /* End of ros_FwdInt */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int ros_DadjInt ( int NADJ, KPP_REAL Lambda[][NVAR], KPP_REAL Tstart, 
		  KPP_REAL Tend, KPP_REAL T, int SaveLU, int ISTATUS[], 
		  KPP_REAL Roundoff, int Autonomous) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Template for the implementation of a generic RosenbrockSOA method
      defined by ros_S (no of stages)
      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~> NADJ - Input: the initial condition at Tstart; Output: the solution at T
!~~~> Lambda[NADJ][NVAR] -  First order adjoint
!~~~> Tstart, Tend - Input: integration interval
!~~~> T - Output: time at which the solution is returned (T=Tend if success)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~ Local variables */
  KPP_REAL Ystage[NVAR*ros_S], K[NVAR*ros_S];
  KPP_REAL U[NADJ][NVAR*ros_S], V[NADJ][NVAR*ros_S];
#ifdef FULL_ALGEBRA
  KPP_REAL Jac[NVAR][NVAR], dJdT[NVAR][NVAR], Ghimj[NVAR][NVAR];
#else
  KPP_REAL Jac[LU_NONZERO], dJdT[LU_NONZERO], Ghimj[LU_NONZERO];
#endif
  KPP_REAL Hes0[NHESS];
  KPP_REAL Tmp[NVAR], Tmp2[NVAR];
  KPP_REAL H=0.0, HC, HA, Tau;
  int Pivot[NVAR], Direction;
  int i, j, m, istage, istart, jstart;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  if (Tend  >=  Tstart)
    Direction = 1;
  else
    Direction = -1;
  
  /*~~~> Time loop begins below */
  while ( stack_ptr > 0 ) { /* TimeLoop */

    /*~~~>  Recover checkpoints for stage values and vectors */
    ros_DPop( ros_S, &T, &H, &Ystage[0], &K[0], &Ghimj[0], &Pivot[0], SaveLU );

    /*~~~>    Compute LU decomposition */
    if (!SaveLU) {
      JacTemplate(T,&Ystage[0],Ghimj);
      ISTATUS[Njac] = ISTATUS[Njac] + 1;
      Tau = ONE/(Direction*H*ros_Gamma[0]);
#ifdef FULL_ALGEBRA
      for(j=0; j<NVAR; j++) {
	for(i=0; i<NVAR; i++)
	  Ghimj[i][j] = -Ghimj[i][j];
      }
      for(i=0; i<NVAR; i++)
	Ghimj[i][i] = Ghimj[i][i]+Tau;
#else
      WSCAL(LU_NONZERO,(-ONE),Ghimj,1);
      for (i=0; i<NVAR; i++)
	Ghimj[LU_DIAG[i]] = Ghimj[LU_DIAG[i]]+Tau;
#endif
      ros_Decomp(Ghimj, Pivot, &j, ISTATUS);
    }

/*~~~>   Compute Hessian at the beginning of the interval */
    HessTemplate(T,&Ystage[0],Hes0);

/*~~~>   Compute the stages */
    for (istage = ros_S - 1; istage >= 0; istage--) { /* Stage loop */

      /*~~~> Current istage first entry */
      istart = NVAR*istage;
	 
      /*~~~> Compute U */
      for (m = 0; m<NADJ; m++) {
	WCOPY(NVAR,&Lambda[m][0],1,&U[m][istart],1);
	WSCAL(NVAR,ros_M[istage],&U[m][istart],1);
      } /* m=0:NADJ-1 */
      for (j = istage+1; j < ros_S; j++) {
	jstart = NVAR*j;
	HA = ros_A[j*(j-1)/2+istage];
	HC = ros_C[j*(j-1)/2+istage]/(Direction*H);
	for ( m = 0; m < NADJ; m++ ) {
	  WAXPY(NVAR,HA,&V[m][jstart],1,&U[m][istart],1);
	  WAXPY(NVAR,HC,&U[m][jstart],1,&U[m][istart],1);
	} /* m=0:NADJ-1 */
      }
      for ( m = 0; m < NADJ; m++ )
	ros_Solve('T', Ghimj, Pivot, &U[m][istart], ISTATUS); /* m=1:NADJ-1 */

      /*~~~> Compute V */
      Tau = T + ros_Alpha[istage]*Direction*H;
      JacTemplate(Tau,&Ystage[istart],Jac);
      ISTATUS[Njac]++;
      for ( m = 0; m < NADJ; m++ ) {
#ifdef FULL_ALGEBRA
	for (i=istart; i < istart+NVAR-1; i++ )
	  V[[m][i] = MATMUL(TRANSPOSE(Jac),U[m][i]];
#else
	JacTR_SP_Vec(Jac,&U[m][istart],&V[m][istart]); 
#endif
      } /* m=0:NADJ-1 */
    } /*End of Stage loop */

    if (!Autonomous)
/*~~~>  Compute the Jacobian derivative with respect to T.
        Last "Jac" computed for stage 1 */
      ros_JacTimeDerivative ( T, Roundoff, &Ystage[0], Jac, dJdT, ISTATUS );

/*~~~>  Compute the new solution */
    /*~~~>  Compute Lambda */ 
    for( istage = 0; istage < ros_S; istage++ ) {
      istart = NVAR*istage;
      for (m = 0; m < NADJ; m++) {
	/* Add V_i */
	WAXPY(NVAR,ONE,&V[m][istart],1,&Lambda[m][0],1);
	/* Add (H0xK_i)^T * U_i */
	HessTR_Vec ( Hes0, &U[m][istart], &K[istart], Tmp );
	WAXPY(NVAR,ONE,Tmp,1,&Lambda[m][0],1);
      } /* m=0:NADJ-1 */
    }

    /* Add H * dJac_dT_0^T * \sum(gamma_i U_i) */
    /* Tmp holds sum gamma_i U_i */
    if (!Autonomous) {
      for( m = 0; m < NADJ; m++ ) {
	for(i=0; i<NVAR; i++)
	  Tmp[i] = ZERO;
	for( istage = 0; istage < ros_S; istage++ ) {
	  istart = NVAR*istage;
	  WAXPY(NVAR,ros_Gamma[istage],&U[m][istart],1,Tmp,1);
	}
#ifdef FULL_ALGEBRA
	Tmp2 = MATMUL(TRANSPOSE(dJdT),Tmp);
#else
	JacTR_SP_Vec(dJdT,Tmp,Tmp2);
#endif
	WAXPY(NVAR,H,Tmp2,1,&Lambda[m][0],1);
      } /* m=0:NADJ-1 */
    } /* .NOT.Autonomous */
  } /* End of TimeLoop */

  /*~~~> Save last state */
  /*~~~> Succesful exit */
  return 1;  /*~~~> The integration was successful */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
} /* End of ros_DadjInt */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int ros_CadjInt ( int NADJ, KPP_REAL Y[][NVAR], KPP_REAL Tstart, KPP_REAL Tend,
		  KPP_REAL T, KPP_REAL AbsTol_adj[][NVAR], 
		  KPP_REAL RelTol_adj[][NVAR], KPP_REAL RSTATUS[], 
		  KPP_REAL Hmin, KPP_REAL Hmax, KPP_REAL Hstart, 
		  KPP_REAL Roundoff, int Max_no_steps, int Autonomous, 
		  int VectorTol, KPP_REAL FacMax, KPP_REAL FacMin, 
		  KPP_REAL FacSafe, KPP_REAL FacRej, int ISTATUS[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Template for the implementation of a generic RosenbrockADJ method
      defined by ros_S (no of stages)
      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
~~~> NADJ, Y[NADJ][NVAR] - Input: the initial condition at Tstart;
                           Output: the solution at T
~~~> Tstart, Tend - Input: integration interval 
~~~> AbsTol_adj[NADJ][NVAR], RelTol_adj[NADJ][NVAR] - Input: adjoint tolerances
~~~> T - Output: time at which the solution is returned (T=Tend if success)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~ Local variables */
  KPP_REAL Y0[NVAR];
  KPP_REAL Ynew[NADJ][NVAR], Fcn0[NADJ][NVAR], Fcn[NADJ][NVAR];
  KPP_REAL K[NADJ][NVAR*ros_S], dFdT[NADJ][NVAR];
#ifdef FULL_ALGEBRA
  KPP_REAL Jac0[NVAR][NVAR], Ghimj[NVAR][NVAR], Jac[NVAR][NVAR], 
         dJdT[NVAR][NVAR];
#else
  KPP_REAL Jac0[LU_NONZERO], Ghimj[LU_NONZERO], Jac[LU_NONZERO], 
         dJdT[LU_NONZERO];
#endif
  KPP_REAL H, Hnew, HC, HG, Fac, Tau; 
  KPP_REAL Err, Yerr[NADJ][NVAR];
  int Pivot[NVAR], Direction, ioffset, j, istage, iadj;
  int RejectLastH, RejectMoreH, Singular; /* Boolean values */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~>  Initial preparations */
  T = Tstart;
  RSTATUS[Nhexit] = (KPP_REAL)0.0;
  H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) );
  if (ABS(H) <= ((KPP_REAL)10.0)*Roundoff) 
    H = DeltaMin;
  
  if (Tend  >=  Tstart)
    Direction = 1;
  else
    Direction = -1;
  
  H = Direction*H;
  RejectLastH = FALSE;
  RejectMoreH = FALSE;

/*~~~> Time loop begins below */
  while ( ((Direction > 0) && ((T-Tend)+Roundoff <= ZERO))
	  || ((Direction < 0) && ((Tend-T)+Roundoff <= ZERO)) ) { /*TimeLoop*/

    if ( ISTATUS[Nstp] > Max_no_steps )  /* Too many steps */
      return ros_ErrorMsg(-6,T,H);

    /* Step size too small */
    if ( ((T+((KPP_REAL)0.1)*H) == T) || (H <= Roundoff) )
      return ros_ErrorMsg(-7,T,H);

/*~~~>  Limit H if necessary to avoid going beyond Tend */
    RSTATUS[Nhexit] = H;
    H = MIN(H,ABS(Tend-T));

/*~~~>   Interpolate forward solution */
    ros_cadj_Y( T, Y0 );
/*~~~>   Compute the Jacobian at current time */
    JacTemplate(T, Y0, Jac0);
    ISTATUS[Njac]++;

/*~~~>  Compute the function derivative with respect to T */
    if (!Autonomous) {
      ros_JacTimeDerivative ( T, Roundoff, Y0, Jac0, dJdT, ISTATUS );
      for (iadj = 0; iadj < NADJ; iadj++) {
#ifdef FULL_ALGEBRA
	for (i=0; i<NVAR; i++)
	  dFdT[iadj][i] = MATMUL(TRANSPOSE(dJdT),Y[iadj][i]);
#else
	JacTR_SP_Vec(dJdT,&Y[iadj][0],&dFdT[iadj][0]);
#endif
	WSCAL(NVAR,(-ONE),&dFdT[iadj][0],1);
      } /* End for loop */
    } /* End if */

/*~~~>  Ydot = -J^T*Y */
#ifdef FULL_ALGEBRA
    int i;
    for(i=0; i<NVAR; i++) {
      for(j=0; j<NVAR; j++)
	Jac0[i][j] = -Jac0[i][j];
    }
#else
    WSCAL(LU_NONZERO,(-ONE),Jac0,1);
#endif

    for ( iadj = 0; iadj < NADJ; iadj++ ) {
#ifdef FULL_ALGEBRA
      int i;
      for(i=0; i<NVAR; i++)
	Fcn0[iadj][i] = MATMUL(TRANSPOSE(Jac0),Y[iadj][i]);
#else
      JacTR_SP_Vec(Jac0,&Y[iadj][0],&Fcn0[iadj][0]);
#endif
    }

/*~~~>  Repeat step calculation until current step accepted */
    do { /* UntilAccepted */

      Singular = ros_PrepareMatrix(H,Direction,ros_Gamma[0], Jac0,Ghimj,Pivot,
				   ISTATUS);
      if (Singular) /* More than 5 consecutive failed decompositions */
	return ros_ErrorMsg(-8,T,H);

/*~~~>   Compute the stages */
      for ( istage = 0; istage < ros_S; istage++ ) { /* Stage loop */

	/* Current istage offset. Current istage vector 
	   is K[ioffset+1:ioffset+NVAR] */
	ioffset = NVAR*(istage-1);

	/* For the 1st istage the function has been computed previously */
	if ( istage == 0 ) {
	  for ( iadj = 0; iadj < NADJ; iadj++ )
	    WCOPY(NVAR,&Fcn0[iadj][0],1,&Fcn[iadj][0],1);

	  /* istage>0 and a new function evaluation is needed at 
	     the current istage */
	}
	else if ( ros_NewF[istage] ) {
	  WCOPY(NVAR*NADJ,&Y[0][0],1,&Ynew[0][0],1);
	  for (j = 0; j < istage-1; j++) {
	    for ( iadj = 0; iadj < NADJ; iadj++ )
	      WAXPY( NVAR,ros_A[(istage-1)*(istage-2)/2+j],
		     &K[iadj][NVAR*(j-1)+1],1,&Ynew[iadj][0],1);
	  } /* End for loop */
	  Tau = T + ros_Alpha[istage]*Direction*H;
	  ros_cadj_Y( Tau, Y0 );
	  JacTemplate(Tau, Y0, Jac);
	  ISTATUS[Njac]++;

#ifdef FULL_ALGEBRA
	  for(i=0; i<NVAR; i++) {
	    for(j=0; j<NVAR; j++)
	      Jac[i][j] = -Jac[i][j];
	  }
#else
	  WSCAL(LU_NONZERO,(-ONE),Jac,1);
#endif

	  for ( iadj = 0; iadj < NADJ; iadj++ ) {
#ifdef FULL_ALGEBRA
	    for(i=0; i<NVAR; i++)
	      Fcn[iadj][i] = MATMUL(TRANSPOSE(Jac),Ynew[iadj][i]);
#else
	    JacTR_SP_Vec(Jac,&Ynew[iadj][0],&Fcn[iadj][0]);
#endif
	  } /* End for loop */
	} /*  if istage == 1 elseif ros_NewF(istage) */

	for ( iadj = 0; iadj < NADJ; iadj++ )
	  WCOPY(NVAR,&Fcn[iadj][0],1,&K[iadj][ioffset+1],1);
	for ( j = 0; j < istage-1; j++ ) {
	  HC = ros_C[(istage-1)*(istage-2)/2+j]/(Direction*H);
	  for ( iadj = 0; iadj < NADJ; iadj++ )
	    WAXPY(NVAR,HC,&K[iadj][NVAR*(j-1)+1],1,&K[iadj][ioffset+1],1);
	} /* End for loop */
	if ((!Autonomous) && (ros_Gamma[istage] != ZERO)) {
	  HG = Direction*H*ros_Gamma[istage];
	  for ( iadj = 0; iadj < NADJ; iadj++ )
	    WAXPY(NVAR,HG,&dFdT[iadj][0],1,&K[iadj][ioffset+1],1);
	} /* End if */
	for ( iadj = 0; iadj < NADJ; iadj++ )
	  ros_Solve('T', Ghimj, Pivot, &K[iadj][ioffset+1], ISTATUS);

      } /* End of Stage loop */

/*~~~>  Compute the new solution */
      for ( iadj = 0; iadj < NADJ; iadj++ ) {
	WCOPY(NVAR,&Y[iadj][0],1,&Ynew[iadj][0],1);
	for ( j=0; j<ros_S; j++ )
	  WAXPY(NVAR,ros_M[j],&K[iadj][NVAR*(j-1)+1],1,&Ynew[iadj][0],1);
      } /* End for loop */

/*~~~>  Compute the error estimation */
      WSCAL(NVAR*NADJ,ZERO,&Yerr[0][0],1);
      for ( j=0; j<ros_S; j++ ) {
	for ( iadj = 0; iadj < NADJ; iadj++ )
	  WAXPY(NVAR,ros_E[j],&K[iadj][NVAR*(j-1)+1],1,&Yerr[iadj][0],1);
      } /* End for loop */

/*~~~> Max error among all adjoint components */
      iadj = 1;
      Err = ros_ErrorNorm ( &Y[iadj][0], &Ynew[iadj][0], &Yerr[iadj][0],
			    &AbsTol_adj[iadj][0], &RelTol_adj[iadj][0], 
			    VectorTol );

/*~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax */
      Fac  = MIN(FacMax,MAX(FacMin,FacSafe/pow(Err,(ONE/ros_ELO))));
      Hnew = H*Fac;

/*~~~>  Check the error magnitude and adjust step size */
      /* ISTATUS[Nstp] = ISTATUS[Nstp] + 1 */
      if ( (Err <= ONE) || (H <= Hmin) ) {  /*~~~> Accept step */
	ISTATUS[Nacc] = ISTATUS[Nacc] + 1;
	WCOPY(NVAR*NADJ,&Ynew[0][0],1,&Y[0][0],1);
	T = T + Direction*H;
	Hnew = MAX(Hmin,MIN(Hnew,Hmax));
	if (RejectLastH) /* No step size increase after a rejected step */
	  Hnew = MIN(Hnew,H);
	RSTATUS[Nhexit] = H;
	RSTATUS[Nhnew] = Hnew;
	RSTATUS[Ntexit] = T;
	RejectLastH = FALSE;
	RejectMoreH = FALSE;
	H = Hnew;
	break; /* UntilAccepted - EXIT THE LOOP: WHILE STEP NOT ACCEPTED */
      }
      else {         /*~~~> Reject step */
	if (RejectMoreH)
	  Hnew = H*FacRej;
	RejectMoreH = RejectLastH;
	RejectLastH = TRUE;
	H = Hnew;
	if (ISTATUS[Nacc] >= 1)
	  ISTATUS[Nrej]++;
      } /* Err <= 1 */
      
    } while(1); /* End of UntilAccepted do loop */

  } /* End of TimeLoop */

  /*~~~> Succesful exit */
  return 1;  /*~~~> The integration was successful */
} /* End of ros_CadjInt */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int ros_SimpleCadjInt ( int NADJ, KPP_REAL Y[][NVAR], KPP_REAL Tstart,
			 KPP_REAL Tend, KPP_REAL T, int ISTATUS[], 
			int Autonomous, KPP_REAL Roundoff ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Template for the implementation of a generic RosenbrockADJ method 
      defined by ros_S (no of stages)
      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
~~~> NADJ, Y[NADJ][NVAR] - Input: the initial condition at Tstart; 
                           Output: the solution at T
~~~> Tstart, Tend - Input: integration interval 
~~~> T - Output: time at which the solution is returned (T=Tend if success) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~ Local variables */
   KPP_REAL Y0[NVAR];
   KPP_REAL Ynew[NADJ][NVAR], Fcn0[NADJ][NVAR], Fcn[NADJ][NVAR];
   KPP_REAL K[NADJ][NVAR*ros_S], dFdT[NADJ][NVAR];
#ifdef FULL_ALGEBRA
   KPP_REAL Jac0[NVAR][NVAR], Ghimj[NVAR][NVAR], Jac[NVAR][NVAR], 
          dJdT[NVAR][NVAR];
#else
   KPP_REAL Jac0[LU_NONZERO], Ghimj[LU_NONZERO], Jac[LU_NONZERO], 
          dJdT[LU_NONZERO];
#endif
   KPP_REAL H, HC, HG, Tau; 
   KPP_REAL ghinv;
   int Pivot[NVAR], Direction, ioffset, i, j, istage, iadj;
   int istack;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~>  INITIAL PREPARATIONS */
  if (Tend  >=  Tstart)
    Direction = -1;
  else
    Direction = 1;
 
/*~~~> Time loop begins below */
  for( istack = stack_ptr; istack >= 1; istack-- ) { /* TimeLoop */
    T = chk_T[istack];
    H = chk_H[istack-1];
    for(i=0; i<NVAR; i++)
      Y0[i] = chk_Y[istack][i];

/*~~~>   Compute the Jacobian at current time */
    JacTemplate(T, Y0, Jac0);
    ISTATUS[Njac] = ISTATUS[Njac] + 1;

/*~~~>  Compute the function derivative with respect to T */
    if (!Autonomous) {
      ros_JacTimeDerivative ( T, Roundoff, Y0, Jac0, dJdT, ISTATUS );
      for ( iadj = 0; iadj < NADJ; iadj++ ) {
#ifdef FULL_ALGEBRA
	for( i=0; i<NVAR; i++)
	  dFdT[iadj][i] = MATMUL(TRANSPOSE(dJdT),Y[iadj][i]);
#else
	JacTR_SP_Vec(dJdT,&Y[iadj][0],&dFdT[iadj][0]);
#endif
	WSCAL(NVAR,(-ONE),&dFdT[iadj][0],1);
      }
    }

/*~~~>  Ydot = -J^T*Y */
#ifdef FULL_ALGEBRA
    for(i=0; i<NVAR; i++) {
      for(j=0; j<NVAR; j++)
	Jac0[i][j] = -Jac0[i][j];
    }
#else
    WSCAL(LU_NONZERO,(-ONE),Jac0,1);
#endif
    
    for(iadj=0; iadj<NADJ; iadj++) {
#ifdef FULL_ALGEBRA
      for(i=0; i<NVAR; i++)
	Fcn0[iadj][i] = MATMUL(TRANSPOSE(Jac0),Y[iadj][i]);
#else
      JacTR_SP_Vec(Jac0,&Y[iadj][0],&Fcn0[iadj][0]);
#endif
    }
    
/*~~~>    Construct Ghimj = 1/(H*ham) - Jac0 */
    ghinv = ONE/(Direction*H*ros_Gamma[0]);
#ifdef FULL_ALGEBRA
    for(i=0; i<NVAR; i++) {
      for(j=0; j<NVAR; j++)
	Ghimj[i][j] = -Jac0[i][j];
    }
    for(i=0; i<NVAR; i++)
      Ghimj[i][i] = Ghimj[i][i]+ghinv;
#else
    WCOPY(LU_NONZERO,Jac0,1,Ghimj,1);
    WSCAL(LU_NONZERO,(-ONE),Ghimj,1);
    for(i=0; i<NVAR; i++)
      Ghimj[LU_DIAG[i]] = Ghimj[LU_DIAG[i]]+ghinv;
#endif

/*~~~>    Compute LU decomposition */
    ros_Decomp( Ghimj, Pivot, &j, ISTATUS );
    if (j != 0) {
      ros_ErrorMsg(-8,T,H);
      printf( "The matrix is singular !");
      exit(0);
    }

/*~~~>   Compute the stages */
    for(istage=0; istage<ros_S; istage++) { /* Stage */
      /* Current istage offset. Current istage vector 
	 is K(ioffset+1:ioffset+NVAR) */
      ioffset = NVAR*istage;

      /* For the 1st istage the function has been computed previously */
      if ( istage == 0 ) {
	for(iadj=0; iadj<NADJ; iadj++)
	  WCOPY(NVAR,&Fcn0[iadj][0],1,&Fcn[iadj][0],1);
      }
      /* istage>=1 and a new function evaluation is needed 
	 at the current istage */
      else if ( ros_NewF[istage] ) {
	WCOPY(NVAR*NADJ,&Y[0][0],1,&Ynew[0][0],1);
	for(j=0; j<istage; j++) {
	  for(iadj=0; iadj<NADJ; iadj++)
	    WAXPY(NVAR,ros_A[istage*(istage-1)/2+j], &K[iadj][NVAR*j],1,
		  &Ynew[iadj][0],1);
	}
	
	Tau = T + ros_Alpha[istage]*Direction*H;
	for(i=0; i<NVAR; i++)
	  ros_Hermite3( chk_T[istack-1], chk_T[istack], Tau,
			&chk_Y[istack-1][i], &chk_Y[istack][i],
			&chk_dY[istack-1][i], &chk_dY[istack][i], Y0 );
	JacTemplate(Tau, Y0, Jac);
	ISTATUS[Njac]++;

#ifdef FULL_ALGEBRA
	for(i=0; i<NVAR; i++) {
	  for(j=0; j<NVAR; j++)
	    Jac[i][j] = -Jac[i][j];
	}
#else
	WSCAL(LU_NONZERO,(-ONE),Jac,1);
#endif

	for(iadj=0; iadj<NADJ; iadj++) {
#ifdef FULL_ALGEBRA
	  for(i=0; i<NVAR; i++)
	    Fcn[iadj][i] = MATMUL(TRANSPOSE(Jac),Ynew[iadj][i]);
#else
	  JacTR_SP_Vec(Jac,&Ynew[iadj][0],&Fcn[iadj][0]);
#endif
	}
      } /* if istage == 1 elseif ros_NewF(istage) */

      for(iadj=0; iadj<NADJ; iadj++)
	WCOPY(NVAR,&Fcn[iadj][0],1,&K[iadj][ioffset],1);
      for(j=0; j<istage-1; j++) {
	HC = ros_C[istage*(istage-1)/2+j]/(Direction*H);
	for(iadj=0; iadj<NADJ; iadj++)
	  WAXPY(NVAR,HC,&K[iadj][NVAR*j],1,&K[iadj][ioffset],1);
      }
      if((!Autonomous) && (ros_Gamma[istage] != ZERO)) {
	HG = Direction*H*ros_Gamma[istage];
	for(iadj=0; iadj<NADJ; iadj++)
	  WAXPY(NVAR,HG,&dFdT[iadj][0],1,&K[iadj][ioffset],1);
      }
      for(iadj=0; iadj<NADJ; iadj++)
	ros_Solve('T', Ghimj, Pivot, &K[iadj][ioffset], ISTATUS);
    } /* End of Stage loop */
    
/*~~~>  Compute the new solution */
    for(iadj=0; iadj<NADJ; iadj++) {
      for(j=0; j<ros_S; j++)
	WAXPY(NVAR,ros_M[j],&K[iadj][NVAR*j],1,&Y[iadj][0],1);
    }
  } /* End of TimeLoop */

/*~~~> Succesful exit */
  return 1;  /*~~~> The integration was successful */

} /* End of ros_SimpleCadjInt */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
KPP_REAL ros_ErrorNorm ( KPP_REAL Y[], KPP_REAL Ynew[], KPP_REAL Yerr[],
		       KPP_REAL AbsTol[], KPP_REAL RelTol[], int VectorTol ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> Computes the "scaled norm" of the error vector Yerr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* Local variables */
  KPP_REAL Err, Scale, Ymax;
  int i;

  Err = ZERO;
  for(i=0; i<NVAR; i++) {
    Ymax = MAX(ABS(Y[i]),ABS(Ynew[i]));
    if (VectorTol)
      Scale = AbsTol[i]+RelTol[i]*Ymax;
    else
      Scale = AbsTol[0]+RelTol[0]*Ymax;

    Err = Err+pow((Yerr[i]/Scale),2);
  }
  Err  = SQRT(Err/NVAR);

  return MAX(Err,(KPP_REAL)1.0e-10);
} /* End of ros_ErrorNorm */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_FunTimeDerivative ( KPP_REAL T, KPP_REAL Roundoff, KPP_REAL Y[],
			     KPP_REAL Fcn0[], KPP_REAL dFdT[], int ISTATUS[]) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> The time partial derivative of the function by finite differences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  KPP_REAL Delta;
  
  Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T));
  FunTemplate(T+Delta,Y,dFdT);
  ISTATUS[Nfun]++;
  WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1);
  WSCAL(NVAR,(ONE/Delta),dFdT,1);

} /* End of ros_FunTimeDerivative */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_JacTimeDerivative ( KPP_REAL T, KPP_REAL Roundoff, KPP_REAL Y[],
			     KPP_REAL Jac0[], KPP_REAL dJdT[], int ISTATUS[]) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~> The time partial derivative of the Jacobian by finite differences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  KPP_REAL Delta;

  Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T));
  JacTemplate(T+Delta,Y,dJdT);
  ISTATUS[Njac]++;
#ifdef FULL_ALGEBRA
  WAXPY(NVAR*NVAR,(-ONE),Jac0,1,dJdT,1);
  WSCAL(NVAR*NVAR,(ONE/Delta),dJdT,1);
#else
  WAXPY(LU_NONZERO,(-ONE),Jac0,1,dJdT,1);
  WSCAL(LU_NONZERO,(ONE/Delta),dJdT,1);
#endif
} /* End of ros_JacTimeDerivative */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int ros_PrepareMatrix ( KPP_REAL H, int Direction, KPP_REAL gam,
			KPP_REAL Jac0[], KPP_REAL Ghimj[], int Pivot[], 
			int ISTATUS[] ) {
/* --- --- --- --- --- --- --- --- --- --- --- --- ---
  Prepares the LHS matrix for stage calculations
  1.  Construct Ghimj = 1/(H*gam) - Jac0
      "(Gamma H) Inverse Minus Jacobian"
  2.  Repeat LU decomposition of Ghimj until successful.
       -half the step size if LU decomposition fails and retry
       -exit after 5 consecutive fails
 --- --- --- --- --- --- --- --- --- --- --- --- --- */

/*~~~> Local variables */
  int i, ising, Nconsecutive;
  int Singular; /* Boolean value */
  KPP_REAL ghinv;

  Nconsecutive = 0;
  Singular = TRUE;

  while (Singular) {

/*~~~>    Construct Ghimj = 1/(H*gam) - Jac0 */
#ifdef FULL_ALGEBRA
    WCOPY(NVAR*NVAR,Jac0,1,Ghimj,1);
    WSCAL(NVAR*NVAR,(-ONE),Ghimj,1);
    ghinv = ONE/(Direction*H*gam);
    for(i=0; i<NVAR; i++)
      Ghimj[i][i] = Ghimj[i][i]+ghinv;
#else
    WCOPY(LU_NONZERO,Jac0,1,Ghimj,1);
    WSCAL(LU_NONZERO,(-ONE),Ghimj,1);
    ghinv = ONE/(Direction*H*gam);
    for(i=0; i<NVAR; i++)
      Ghimj[LU_DIAG[i]] = Ghimj[LU_DIAG[i]]+ghinv;
#endif

/*~~~>    Compute LU decomposition */
    ros_Decomp( Ghimj, Pivot, &ising, ISTATUS );
    if (ising == 0)
/*~~~>    If successful done */
      Singular = FALSE;
    
    else { /* ising != 0 */
/*~~~>    If unsuccessful half the step size; 
          if 5 consecutive fails then return */
      ISTATUS[Nsng]++;
      Nconsecutive++;
      Singular = TRUE;
      printf( "Warning: LU Decomposition returned ising = %d", ising );
      if (Nconsecutive <= 5) /*Less than 5 consecutive failed decompositions*/
	H = H*HALF;
      else  /* More than 5 consecutive failed decompositions */
	return Singular;

    }  /* End of ising */
  }  /* while Singular */

  return Singular;
} /* End of ros_PrepareMatrix */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_Decomp( KPP_REAL A[], int Pivot[], int* ising, int ISTATUS[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Template for the LU decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef FULL_ALGEBRA
  DGETRF( NVAR, NVAR, A, NVAR, Pivot, ising );
#else
  *ising = KppDecomp ( A );
  Pivot[0] = 1;
#endif
  ISTATUS[Ndec]++;

} /* End of ros_Decomp */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_Solve( char How, KPP_REAL A[], int Pivot[], KPP_REAL b[], 
		int ISTATUS[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Template for the forward/backward substitution 
  (using pre-computed LU decomposition)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  switch (How) {
    case 'N':
#ifdef FULL_ALGEBRA    
      DGETRS( 'N', NVAR , 1, A, NVAR, Pivot, b, NVAR, 0 );
#else
      KppSolve( A, b );
#endif
      break;
    case 'T':
#ifdef FULL_ALGEBRA
      DGETRS( 'T', NVAR , 1, A, NVAR, Pivot, b, NVAR, 0 );
#else
      KppSolveTR( A, b, b );
#endif
      break;
    default:
      printf( "Error: unknown argument in ros_Solve: How=%d", How );
      exit(0);
  }
  ISTATUS[Nsol]++;

} /* End of ros_Solve */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_cadj_Y( KPP_REAL T, KPP_REAL Y[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Finds the solution Y at T by interpolating the stored forward trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  int i,j;

  if( (T < chk_T[0]) || (T> chk_T[stack_ptr]) ) {
    printf( "Cannot locate solution at T = %f", T );
    printf( "Stored trajectory is between Tstart = %f", chk_T[0] );
    printf( "    and Tend = %f", chk_T[stack_ptr] );
    exit(0);
  }
  for(i=0; i<stack_ptr-1; i++) {
    if( (T >= chk_T[i]) &&(T <= chk_T[i+1]) )
      exit(0);
  }

  for(j=0; j<NVAR; j++ )
    ros_Hermite3( chk_T[i], chk_T[i+1], T, &chk_Y[i][j], &chk_Y[i+1][j],
		  &chk_dY[i][j],  &chk_dY[i+1][j], Y );

} /* End of ros_cadj_Y */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_Hermite3( KPP_REAL a, KPP_REAL b, KPP_REAL T, KPP_REAL Ya[],
		   KPP_REAL Yb[], KPP_REAL Ja[], KPP_REAL Jb[], KPP_REAL Y[]) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Template for Hermite interpolation of order 5 on the interval [a,b]
   P = c(1) + c(2)*(x-a) + ... + c(4)*(x-a)^3
   P[a,b] = [Ya,Yb], P'[a,b] = [Ja,Jb]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  KPP_REAL Tau, amb[3], C[4][NVAR];
  int i, j;

  amb[0] = ((KPP_REAL)1.0)/(a-b);
  for(i=1; i<3; i++)
    amb[i] = amb[i-1]*amb[0];

/* c(1) = ya; */
  WCOPY(NVAR,Ya,1,&C[0][0],1);
/* c(2) = ja; */
  WCOPY(NVAR,Ja,1,&C[1][0],1);
/* c(3) = 2/(a-b)*ja + 1/(a-b)*jb - 3/(a - b)^2*ya + 3/(a - b)^2*yb; */
  WCOPY(NVAR,Ya,1,&C[2][0],1);
  WSCAL(NVAR,-3.0*amb[1],&C[2][0],1);
  WAXPY(NVAR,3.0*amb[1],Yb,1,&C[2][0],1);
  WAXPY(NVAR,2.0*amb[0],Ja,1,&C[2][0],1);
  WAXPY(NVAR,amb[0],Jb,1,&C[2][0],1);
/* c(4) =  1/(a-b)^2*ja + 1/(a-b)^2*jb - 2/(a-b)^3*ya + 2/(a-b)^3*yb */
  WCOPY(NVAR,Ya,1,&C[3][0],1);
  WSCAL(NVAR,-2.0*amb[2],&C[3][0],1);
  WAXPY(NVAR,2.0*amb[2],Yb,1,&C[3][0],1);
  WAXPY(NVAR,amb[1],Ja,1,&C[3][0],1);
  WAXPY(NVAR,amb[1],Jb,1,&C[3][0],1);
   
  Tau = T - a;
  WCOPY(NVAR,&C[3][0],1,Y,1);
  WSCAL(NVAR,pow(Tau,3),Y,1);
  for(j=2; j>=0; j--)
    WAXPY(NVAR,pow(Tau,(j-1)),&C[j][0],1,Y,1);

} /* End of ros_Hermite3 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void ros_Hermite5( KPP_REAL a, KPP_REAL b, KPP_REAL T, KPP_REAL Ya[], 
		   KPP_REAL Yb[], KPP_REAL Ja[], KPP_REAL Jb[], KPP_REAL Ha[], 
		   KPP_REAL Hb[], KPP_REAL Y[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Template for Hermite interpolation of order 5 on the interval [a,b]
 P = c(1) + c(2)*(x-a) + ... + c(6)*(x-a)^5
 P[a,b] = [Ya,Yb], P'[a,b] = [Ja,Jb], P"[a,b] = [Ha,Hb]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  KPP_REAL Tau, amb[5], C[6][NVAR];
  int i, j;

  amb[0] = ((KPP_REAL)1.0)/(a-b);
  for(i=1; i<5; i++)
    amb[i] = amb[i-1]*amb[0];

/* c(1) = ya; */
  WCOPY(NVAR,Ya,1,&C[0][0],1);
/* c(2) = ja; */
  WCOPY(NVAR,Ja,1,&C[1][0],1);
/* c(3) = ha/2; */
  WCOPY(NVAR,Ha,1,&C[2][0],1);
  WSCAL(NVAR,HALF,&C[2][0],1);

/* c(4) = 10*amb(3)*ya - 10*amb(3)*yb - 6*amb(2)*ja - 4*amb(2)*jb  
          + 1.5*amb(1)*ha - 0.5*amb(1)*hb ; */
  WCOPY(NVAR,Ya,1,&C[3][0],1);
  WSCAL(NVAR,10.0*amb[2],&C[3][0],1);
  WAXPY(NVAR,-10.0*amb[2],Yb,1,&C[3][0],1);
  WAXPY(NVAR,-6.0*amb[1],Ja,1,&C[3][0],1);
  WAXPY(NVAR,-4.0*amb[1],Jb,1,&C[3][0],1);
  WAXPY(NVAR, 1.5*amb[0],Ha,1,&C[3][0],1);
  WAXPY(NVAR,-0.5*amb[0],Hb,1,&C[3][0],1);

/* c(5) =   15*amb(4)*ya - 15*amb(4)*yb - 8.*amb(3)*ja - 7*amb(3)*jb 
            + 1.5*amb(2)*ha - 1*amb(2)*hb ; */
  WCOPY(NVAR,Ya,1,&C[4][0],1);
  WSCAL(NVAR, 15.0*amb[3],&C[4][0],1);
  WAXPY(NVAR,-15.0*amb[3],Yb,1,&C[4][0],1);
  WAXPY(NVAR,-8.0*amb[2],Ja,1,&C[4][0],1);
  WAXPY(NVAR,-7.0*amb[2],Jb,1,&C[4][0],1);
  WAXPY(NVAR,1.5*amb[1],Ha,1,&C[4][0],1);
  WAXPY(NVAR,-amb[1],Hb,1,&C[4][0],1);

/* c(6) =   6*amb(5)*ya - 6*amb(5)*yb - 3.*amb(4)*ja - 3.*amb(4)*jb 
            + 0.5*amb(3)*ha -0.5*amb(3)*hb ; */
  WCOPY(NVAR,Ya,1,&C[5][0],1);
  WSCAL(NVAR, 6.0*amb[4],&C[5][0],1);
  WAXPY(NVAR,-6.0*amb[4],Yb,1,&C[5][0],1);
  WAXPY(NVAR,-3.0*amb[3],Ja,1,&C[5][0],1);
  WAXPY(NVAR,-3.0*amb[3],Jb,1,&C[5][0],1);
  WAXPY(NVAR, 0.5*amb[2],Ha,1,&C[5][0],1);
  WAXPY(NVAR,-0.5*amb[2],Hb,1,&C[5][0],1);

  Tau = T - a;
  WCOPY(NVAR,&C[5][0],1,Y,1);
  for(j=4; j>=0; j--) {
    WSCAL(NVAR,Tau,Y,1);
    WAXPY(NVAR,ONE,&C[j][0],1,Y,1);
  }

} /* End of ros_Hermite5 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Ros2() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 --- AN L-STABLE METHOD, 2 stages, order 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  KPP_REAL g;

  g = (KPP_REAL)1.0 + ((KPP_REAL)1.0)/SQRT((KPP_REAL)2.0);

  rosMethod = RS2;
/*~~~> Name of the method */
  strcpy(ros_Name, "ROS-2");
/*~~~> Number of stages */
  ros_S = 2;

/*~~~> The coefficient matrices A and C are strictly lower triangular.
  The lower triangular (subdiagonal) elements are stored in row-wise order:
  A[0][1] = ros_A[0], A[0][2]=ros_A[1], A[1][2]=ros_A[2], etc.
  The general mapping formula is:
      A[i][j] = ros_A[ (i-1)*(i-2)/2 + j ]
      C[i][j] = ros_C[ (i-1)*(i-2)/2 + j ] */

  ros_A[0] = ((KPP_REAL)1.0)/g;
  ros_C[0] = ((KPP_REAL)-2.0)/g;

/*~~~> Does the stage i require a new function evaluation (ros_NewF[i]=TRUE)
       or does it re-use the function evaluation from stage i-1 
       (ros_NewF[i]=FALSE) */
  ros_NewF[0] = TRUE;
  ros_NewF[1] = TRUE;

/*~~~> M_i = Coefficients for new step solution */
  ros_M[0]= ((KPP_REAL)3.0)/((KPP_REAL)2.0*g);
  ros_M[1]= ((KPP_REAL)1.0)/((KPP_REAL)2.0*g);

/* E_i = Coefficients for error estimator */
  ros_E[0] = ((KPP_REAL)1.0)/((KPP_REAL)2.0*g);
  ros_E[1] = ((KPP_REAL)1.0)/((KPP_REAL)2.0*g);

/*~~~> ros_ELO = estimator of local order - the minimum between the
       main and the embedded scheme orders plus one */
  ros_ELO = (KPP_REAL)2.0;

/*~~~> Y_stage_i ~ Y( T + H*Alpha_i ) */
  ros_Alpha[0] = (KPP_REAL)0.0;
  ros_Alpha[1] = (KPP_REAL)1.0;

/*~~~> Gamma_i = \sum_j  gamma_{i,j} */
  ros_Gamma[0] = g;
  ros_Gamma[1] =-g;

} /* End of Ros2 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Ros3() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  rosMethod = RS3;
/*~~~> Name of the method */
  strcpy(ros_Name, "ROS-3");
/*~~~> Number of stages */
  ros_S = 3;

/*~~~> The coefficient matrices A and C are strictly lower triangular.
   The lower triangular (subdiagonal) elements are stored in row-wise order:
   A[0][1] = ros_A[0], A[0][2]=ros_A[1], A[1][2]=ros_A[2], etc.
   The general mapping formula is:
       A[i][j] = ros_A[ (i-1)*(i-2)/2 + j ]
       C[i][j] = ros_C[ (i-1)*(i-2)/2 + j ] */

  ros_A[0]= (KPP_REAL)1.0;
  ros_A[1]= (KPP_REAL)1.0;
  ros_A[2]= (KPP_REAL)0.0;
  ros_C[0] = (KPP_REAL)-0.10156171083877702091975600115545e01;
  ros_C[1] = (KPP_REAL) 0.40759956452537699824805835358067e01;
  ros_C[2] = (KPP_REAL) 0.92076794298330791242156818474003e01;

/*~~~> Does the stage i require a new function evaluation (ros_NewF[i]=TRUE)
       or does it re-use the function evaluation from stage i-1 
       (ros_NewF[i]=FALSE) */
  ros_NewF[0] = TRUE;
  ros_NewF[1] = TRUE;
  ros_NewF[2] = FALSE;
/*~~~> M_i = Coefficients for new step solution */
  ros_M[0] = (KPP_REAL) 0.1e01;
  ros_M[1] = (KPP_REAL) 0.61697947043828245592553615689730e01;
  ros_M[2] = (KPP_REAL)-0.42772256543218573326238373806514;
/* E_i = Coefficients for error estimator */
  ros_E[0] = (KPP_REAL) 0.5;
  ros_E[1] = (KPP_REAL)-0.29079558716805469821718236208017e01;
  ros_E[2] = (KPP_REAL) 0.22354069897811569627360909276199;

/*~~~> ros_ELO = estimator of local order - the minimum between the
       main and the embedded scheme orders plus 1 */
  ros_ELO = (KPP_REAL)3.0;
/*~~~> Y_stage_i ~ Y( T + H*Alpha_i ) */
  ros_Alpha[0]= (KPP_REAL)0.0;
  ros_Alpha[1]= (KPP_REAL)0.43586652150845899941601945119356;
  ros_Alpha[2]= (KPP_REAL)0.43586652150845899941601945119356;
/*~~~> Gamma_i = \sum_j  gamma_{i,j} */
  ros_Gamma[0]= (KPP_REAL)0.43586652150845899941601945119356;
  ros_Gamma[1]= (KPP_REAL)0.24291996454816804366592249683314;
  ros_Gamma[2]= (KPP_REAL)0.21851380027664058511513169485832e01;

} /* End of Ros3 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Ros4() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3

      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
      SPRINGER-VERLAG (1990)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  rosMethod = RS4;
/*~~~> Name of the method */
  strcpy(ros_Name, "ROS-4");
/*~~~> Number of stages */
  ros_S = 4;

/*~~~> The coefficient matrices A and C are strictly lower triangular.
   The lower triangular (subdiagonal) elements are stored in row-wise order:
   A[0][1] = ros_A[0], A[0][2]=ros_A[1], A[1][2]=ros_A[2], etc.
   The general mapping formula is:
       A[i][j] = ros_A[ (i-1)*(i-2)/2 + j ]
       C[i][j] = ros_C[ (i-1)*(i-2)/2 + j ] */

  ros_A[0] = (KPP_REAL)0.2000000000000000e01;
  ros_A[1] = (KPP_REAL)0.1867943637803922e01;
  ros_A[2] = (KPP_REAL)0.2344449711399156;
  ros_A[3] = ros_A[1];
  ros_A[4] = ros_A[2];
  ros_A[5] = (KPP_REAL)0.0;

  ros_C[0] = (KPP_REAL)-0.7137615036412310e01;
  ros_C[1] = (KPP_REAL) 0.2580708087951457e01;
  ros_C[2] = (KPP_REAL) 0.6515950076447975;
  ros_C[3] = (KPP_REAL)-0.2137148994382534e01;
  ros_C[4] = (KPP_REAL)-0.3214669691237626;
  ros_C[5] = (KPP_REAL)-0.6949742501781779;

/*~~~> Does the stage i require a new function evaluation (ros_NewF[i]=TRUE)
       or does it re-use the function evaluation from stage i-1 
       (ros_NewF[i]=FALSE) */
  ros_NewF[0] = TRUE;
  ros_NewF[1] = TRUE;
  ros_NewF[2] = TRUE;
  ros_NewF[3] = FALSE;
/*~~~> M_i = Coefficients for new step solution */
  ros_M[0] = (KPP_REAL)0.2255570073418735e01;
  ros_M[1] = (KPP_REAL)0.2870493262186792;
  ros_M[2] = (KPP_REAL)0.4353179431840180;
  ros_M[3] = (KPP_REAL)0.1093502252409163e01;
/*~~~> E_i  = Coefficients for error estimator */
  ros_E[0] = (KPP_REAL)-0.2815431932141155;
  ros_E[1] = (KPP_REAL)-0.7276199124938920e-01;
  ros_E[2] = (KPP_REAL)-0.1082196201495311;
  ros_E[3] = (KPP_REAL)-0.1093502252409163e01;

/*~~~> ros_ELO  = estimator of local order - the minimum between the
       main and the embedded scheme orders plus 1 */
  ros_ELO = (KPP_REAL)4.0;
/*~~~> Y_stage_i ~ Y( T + H*Alpha_i ) */
  ros_Alpha[0] = (KPP_REAL)0.0;
  ros_Alpha[1] = (KPP_REAL)0.1145640000000000e01;
  ros_Alpha[2] = (KPP_REAL)0.6552168638155900;
  ros_Alpha[3] = ros_Alpha[2];
/*~~~> Gamma_i = \sum_j  gamma_{i,j} */ 
  ros_Gamma[0] = (KPP_REAL) 0.5728200000000000;
  ros_Gamma[1] = (KPP_REAL)-0.1769193891319233e01;
  ros_Gamma[2] = (KPP_REAL) 0.7592633437920482;
  ros_Gamma[3] = (KPP_REAL)-0.1049021087100450;

} /* End of Ros4 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Rodas3() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 --- A STIFFLY-STABLE METHOD, 4 stages, order 3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  rosMethod = RD3;
/*~~~> Name of the method */
  strcpy(ros_Name, "RODAS-3");
/*~~~> Number of stages */
  ros_S = 4;

/*~~~> The coefficient matrices A and C are strictly lower triangular.
   The lower triangular (subdiagonal) elements are stored in row-wise order:
   A[0][1] = ros_A[0], A[0][2]=ros_A[1], A[1][2]=ros_A[2], etc.
   The general mapping formula is:
       A[i][j] = ros_A[ (i-1)*(i-2)/2 + j ]
       C[i][j] = ros_C[ (i-1)*(i-2)/2 + j ] */

  ros_A[0] = (KPP_REAL)0.0;
  ros_A[1] = (KPP_REAL)2.0;
  ros_A[2] = (KPP_REAL)0.0;
  ros_A[3] = (KPP_REAL)2.0;
  ros_A[4] = (KPP_REAL)0.0;
  ros_A[5] = (KPP_REAL)1.0;

  ros_C[0] = (KPP_REAL) 4.0;
  ros_C[1] = (KPP_REAL) 1.0;
  ros_C[2] = (KPP_REAL)-1.0;
  ros_C[3] = (KPP_REAL) 1.0;
  ros_C[4] = (KPP_REAL)-1.0;
  ros_C[5] = -(((KPP_REAL)8.0)/((KPP_REAL)3.0));

/*~~~> Does the stage i require a new function evaluation (ros_NewF[i]=TRUE)
       or does it re-use the function evaluation from stage i-1 
       (ros_NewF[i]=FALSE) */
  ros_NewF[0] = TRUE;
  ros_NewF[1] = FALSE;
  ros_NewF[2] = TRUE;
  ros_NewF[3] = TRUE;
/*~~~> M_i = Coefficients for new step solution */
  ros_M[0] = (KPP_REAL)2.0;
  ros_M[1] = (KPP_REAL)0.0;
  ros_M[2] = (KPP_REAL)1.0;
  ros_M[3] = (KPP_REAL)1.0;
/*~~~> E_i  = Coefficients for error estimator */
  ros_E[0] = (KPP_REAL)0.0;
  ros_E[1] = (KPP_REAL)0.0;
  ros_E[2] = (KPP_REAL)0.0;
  ros_E[3] = (KPP_REAL)1.0;

/*~~~> ros_ELO  = estimator of local order - the minimum between the
    main and the embedded scheme orders plus 1 */
  ros_ELO  = (KPP_REAL)3.0;
/*~~~> Y_stage_i ~ Y( T + H*Alpha_i ) */
  ros_Alpha[0] = (KPP_REAL)0.0;
  ros_Alpha[1] = (KPP_REAL)0.0;
  ros_Alpha[2] = (KPP_REAL)1.0;
  ros_Alpha[3] = (KPP_REAL)1.0;
/*~~~> Gamma_i = \sum_j  gamma_{i,j} */
  ros_Gamma[0] = (KPP_REAL)0.5;
  ros_Gamma[1] = (KPP_REAL)1.5;
  ros_Gamma[2] = (KPP_REAL)0.0;
  ros_Gamma[3] = (KPP_REAL)0.0;

} /* End of Rodas3 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Rodas4() {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES

      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
      SPRINGER-VERLAG (1996)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  rosMethod = RD4;
/*~~~> Name of the method */
  strcpy(ros_Name, "RODAS-4");
/*~~~> Number of stages */
  ros_S = 6;

/*~~~> Y_stage_i ~ Y( T + H*Alpha_i ) */
  ros_Alpha[0] = (KPP_REAL)0.000;
  ros_Alpha[1] = (KPP_REAL)0.386;
  ros_Alpha[2] = (KPP_REAL)0.210;
  ros_Alpha[3] = (KPP_REAL)0.630;
  ros_Alpha[4] = (KPP_REAL)1.000;
  ros_Alpha[5] = (KPP_REAL)1.000;

/*~~~> Gamma_i = \sum_j  gamma_{i,j} */
  ros_Gamma[0] = (KPP_REAL) 0.2500000000000000;
  ros_Gamma[1] = (KPP_REAL)-0.1043000000000000;
  ros_Gamma[2] = (KPP_REAL) 0.1035000000000000;
  ros_Gamma[3] = (KPP_REAL)-0.3620000000000023e-01;
  ros_Gamma[4] = (KPP_REAL) 0.0;
  ros_Gamma[5] = (KPP_REAL) 0.0;

/*~~~> The coefficient matrices A and C are strictly lower triangular.
   The lower triangular (subdiagonal) elements are stored in row-wise order:
   A[0][1] = ros_A[0], A[0][2]=ros_A[1], A[1][2]=ros_A[2], etc.
   The general mapping formula is:  A[i][j] = ros_A[ (i-1)*(i-2)/2 + j ]
                  C[i][j] = ros_C[ (i-1)*(i-2)/2 + j ] */

  ros_A[0] = (KPP_REAL) 0.1544000000000000e01;
  ros_A[1] = (KPP_REAL) 0.9466785280815826;
  ros_A[2] = (KPP_REAL) 0.2557011698983284;
  ros_A[3] = (KPP_REAL) 0.3314825187068521e01;
  ros_A[4] = (KPP_REAL) 0.2896124015972201e01;
  ros_A[5] = (KPP_REAL) 0.9986419139977817;
  ros_A[6] = (KPP_REAL) 0.1221224509226641e01;
  ros_A[7] = (KPP_REAL) 0.6019134481288629e01;
  ros_A[8] = (KPP_REAL) 0.1253708332932087e02;
  ros_A[9] = (KPP_REAL)-0.6878860361058950;
  ros_A[10] = ros_A[6];
  ros_A[11] = ros_A[7];
  ros_A[12] = ros_A[8];
  ros_A[13] = ros_A[9];
  ros_A[14] = (KPP_REAL)1.0;
  
  ros_C[0]  = (KPP_REAL)-0.5668800000000000e01;
  ros_C[1]  = (KPP_REAL)-0.2430093356833875e01;
  ros_C[2]  = (KPP_REAL)-0.2063599157091915;
  ros_C[3]  = (KPP_REAL)-0.1073529058151375;
  ros_C[4]  = (KPP_REAL)-0.9594562251023355e01;
  ros_C[5]  = (KPP_REAL)-0.2047028614809616e02;
  ros_C[6]  = (KPP_REAL) 0.7496443313967647e01;
  ros_C[7]  = (KPP_REAL)-0.1024680431464352e02;
  ros_C[8]  = (KPP_REAL)-0.3399990352819905e02;
  ros_C[9]  = (KPP_REAL) 0.1170890893206160e02;
  ros_C[10] = (KPP_REAL) 0.8083246795921522e01;
  ros_C[11] = (KPP_REAL)-0.7981132988064893e01;
  ros_C[12] = (KPP_REAL)-0.3152159432874371e02;
  ros_C[13] = (KPP_REAL) 0.1631930543123136e02;
  ros_C[14] = (KPP_REAL)-0.6058818238834054e01;

/*~~~> M_i = Coefficients for new step solution */
  ros_M[0] = ros_A[6];
  ros_M[1] = ros_A[7];
  ros_M[2] = ros_A[8];
  ros_M[3] = ros_A[9];
  ros_M[4] = (KPP_REAL)1.0;
  ros_M[5] = (KPP_REAL)1.0;

/*~~~> E_i  = Coefficients for error estimator */
  ros_E[0] = (KPP_REAL)0.0;
  ros_E[1] = (KPP_REAL)0.0;
  ros_E[2] = (KPP_REAL)0.0;
  ros_E[3] = (KPP_REAL)0.0;
  ros_E[4] = (KPP_REAL)0.0;
  ros_E[5] = (KPP_REAL)1.0;

/*~~~> Does the stage i require a new function evaluation (ros_NewF[i]=TRUE)
       or does it re-use the function evaluation from stage i-1 
       (ros_NewF[i]=FALSE) */
  ros_NewF[0] = TRUE;
  ros_NewF[1] = TRUE;
  ros_NewF[2] = TRUE;
  ros_NewF[3] = TRUE;
  ros_NewF[4] = TRUE;
  ros_NewF[5] = TRUE;

/*~~~> ros_ELO  = estimator of local order - the minimum between the
        main and the embedded scheme orders plus 1 */
  ros_ELO = (KPP_REAL)4.0;

} /* End of Rodas4 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void FunTemplate( KPP_REAL T, KPP_REAL Y[], KPP_REAL Ydot[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Template for the ODE function call.
  Updates the rate coefficients (and possibly the fixed species) at each call
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  KPP_REAL Told;

  Told = TIME;
  TIME = T;
  Update_SUN();
  Update_RCONST();
  Fun( Y, FIX, RCONST, Ydot );
  TIME = Told;

} /* End of FunTemplate */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void JacTemplate( KPP_REAL T, KPP_REAL Y[], KPP_REAL Jcb[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Template for the ODE Jacobian call.
  Updates the rate coefficients (and possibly the fixed species) at each call
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  KPP_REAL Told;
#ifdef FULL_ALGEBRA
  KPP_REAL JV[LU_NONZERO];
  int i, j;
#endif

  Told = TIME;
  TIME = T;
  Update_SUN();
  Update_RCONST();
#ifdef FULL_ALGEBRA
  Jac_SP(Y, FIX, RCONST, JV);
  for(j=0; j<NVAR; j++) {
    for(i=0; i<NVAR; i++)
      Jcb[i][j] = (KPP_REAL)0.0;
  }
  for(i=0; i<LU_NONZERO; i++)
    Jcb[LU_ICOL[i]][LU_IROW[i]] = JV[i];
#else
  Jac_SP( Y, FIX, RCONST, Jcb );
#endif
  TIME = Told;
} /* End of JacTemplate */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void HessTemplate( KPP_REAL T, KPP_REAL Y[], KPP_REAL Hes[] ) {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Template for the ODE Hessian call.
  Updates the rate coefficients (and possibly the fixed species) at each call
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~> Local variables */
  KPP_REAL Told;

  Told = TIME;
  TIME = T;
  Update_SUN();
  Update_RCONST();
  Hessian( Y, FIX, RCONST, Hes );
  TIME = Told;

} /* End of HessTemplate */

/* End of INTEGRATE function                                                 */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
