/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Driver for the Adjoint (ADJ) model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

int  InitSaveData();
void Initialize();
int  SaveData();
int  CloseSaveData();
int  GenerateMatlab( char * prefix );
void GetMass( KPP_REAL CL[], KPP_REAL Mass[] );
void INTEGRATE_ADJ( int NADJ, KPP_REAL Y[], KPP_REAL Lambda[][NVAR], 
		    KPP_REAL TIN, KPP_REAL TOUT, KPP_REAL ATOL_adj[][NVAR], 
		    KPP_REAL RTOL_adj[][NVAR], int ICNTRL_U[], 
		    KPP_REAL RCNTRL_U[], int ISTATUS_U[], 
		    KPP_REAL RSTATUS_U[] );

int main() {

  KPP_REAL T, DVAL[NSPEC];
  int i, j, ind_1 = ind_O1D, ind_2 = ind_O3;

  /*~~>  NADJ = Number of functionals for which sensitivities are computed
     Note: the setting below is for sensitivities of all final concentrations;
     the setting may have to be changed for other applications */
  int NADJ = NVAR;
  KPP_REAL Y_adj[NADJ][NVAR], ATOL_adj[NADJ][NVAR], RTOL_adj[NADJ][NVAR];

/*~~>  Control (in) and status (out) arguments for the integration */
  KPP_REAL RCNTRL[20], RSTATUS[20];
  int ICNTRL[20], ISTATUS[20];
  
  STEPMIN = (double)0.0;
  STEPMAX = (double)0.0;

/*~~~> Tolerances for calculating concentrations */       
  for( i=0; i<NVAR; i++) {
    RTOL[i] = 1.0e-4;
    ATOL[i] = 1.0e-3;
  }
      
/*~~~> Tolerances for calculating adjoints 
       are used for controlling adjoint truncation error
       and for solving the linear adjoint equations by iterations  
       Note: Adjoints typically span many orders of magnitude
       and a careful tuning of ATOL_adj may be necessary */    
  for(i=0; i<NADJ; i++) {
    for(j=0; j<NVAR; j++) {
      RTOL_adj[i][j] = 1.0e-4;
      ATOL_adj[i][j] = 1.0e-10;
    }
  }
     
  Initialize();
      
/*~~~>  The adjoint values at the final time */
  for(i=0; i<NADJ; i++) {
    for(j=0; j<NVAR; j++)
      Y_adj[i][j] = (double)0.0;
  }
  for(j=0; j<NADJ; j++)
    Y_adj[j][j] = (double)1.0;

/*~~~> Default control options */
  for(i=0; i<20; i++) {
    ICNTRL[i] = 0;
    RCNTRL[i] = (double)0.0;
    ISTATUS[i] = 0;
    RSTATUS[i] = (double)0.0;
  }     

/*~~~> Begin time loop */

  TIME = TSTART;
  InitSaveData();

  T = TSTART;

  GetMass( C, DVAL );
  printf("\n%6.1f%% %7.2f   ", (TIME-TSTART)/(TEND-TSTART)*100, TIME/3600 );
  for( i = 0; i < NMONITOR; i++ ) 
    printf( "%9.3e  ", C[ MONITOR[i] ]/CFACTOR );
  for( i = 0; i < NMASS; i++ ) 
    printf( "%9.3e  ", DVAL[i]/CFACTOR );

  TIME = T;
  SaveData();

  INTEGRATE_ADJ( NADJ, VAR, Y_adj, T, TEND, ATOL_adj, RTOL_adj, ICNTRL, 
		 RCNTRL, ISTATUS, RSTATUS );

  GetMass( C, DVAL );

  printf("\n%6.1f%% %7.2f   ", (TEND-TSTART)/(TEND-TSTART)*100, TIME/3600 );
  for( i = 0; i < NMONITOR; i++ ) 
    printf( "%9.3e  ", C[ MONITOR[i] ]/CFACTOR );
  for( i = 0; i < NMASS; i++ ) 
    printf( "%9.3e   ", DVAL[i]/CFACTOR );

  TIME = T;
  SaveData();

/*~~~> End time loop ~~~~~~~~~~*/

  printf( "\n\n****************************************************\n" );
  printf( " Concentrations and Sensitivities at final time\n" );
  printf( " were written in the file KPP_ROOT_ADJ_results.m\n");
  printf( "****************************************************\n");

  FILE *out;
  out = fopen("KPP_ROOT_ADJ_results.m", "w");
  if(out == NULL) {
    printf("Unable to open file KPP_ROOT_ADJ_results.m\n");
    exit(1);
  }

  for(j=0; j<NADJ; j++) {
    for(i=0; i<NVAR; i++)
      fprintf( out, "%7.6e   ", Y_adj[j][i] );
    fprintf( out, "\n" );
  }

  fclose(out);

  printf( "\nADJ: d[%s](tf)/d[%s](t0) = %7.6e\n", SPC_NAMES[ind_1], 
	  SPC_NAMES[ind_1], Y_adj[ind_1][ind_1] );
  printf( "ADJ: d[%s](tf)/d[%s](t0) = %7.6e\n", SPC_NAMES[ind_2], 
	  SPC_NAMES[ind_2], Y_adj[ind_2][ind_2] );
  printf( "ADJ: d[%s](tf)/d[%s](t0) = %7.6e\n", SPC_NAMES[ind_2], 
	  SPC_NAMES[ind_1], Y_adj[ind_2][ind_1] );
  printf( "ADJ: d[%s](tf)/d[%s](t0) =  %7.6e\n\n", SPC_NAMES[ind_1], 
	  SPC_NAMES[ind_2], Y_adj[ind_1][ind_2] );

  CloseSaveData();
 
/*~~~> The entire matrix of sensitivities */
  printf("   ");
  for(i=0; i<NVAR; i++)
    printf( "%15s", SPC_NAMES[i] );
  printf("\n");
  for(j=0; j<NVAR; j++) {
    printf("%5s = ", SPC_NAMES[j]);
    for(i=0; i<NADJ; i++)
      printf( " %7.6e   ",  Y_adj[i][j] );
    printf("\n");
  }

  return 0;
}
