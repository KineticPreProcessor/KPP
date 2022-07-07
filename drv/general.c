//=========================================================================
// general.c -- KPP "box model" driver program for C
//=========================================================================

// Function prototypes
int  InitSaveData();
void Initialize();
int  SaveData();
int  CloseSaveData();
int  GenerateMatlab( char * prefix );
void GetMass( double CL[], double Mass[] );
void INTEGRATE( double TIN, double TOUT );

int main()
{
  KPP_REAL dval[NSPEC];
  int i;


  //=========================================================================
  // Initialization
  //=========================================================================

  // Time variables
  TSTART = 3600.0 * 12.0;
  TEND   = TSTART + 3600.0 * 24.0 * 5.0;
  DT     = 3600.0;

  // Initial temperature
  TEMP   = 236.21;

  // Initial species concentrations
  Initialize();

  // Set relative & absolute tolerances (use same values as for Fortran90)
  for( i = 0; i < NVAR; i++ ) {
    RTOL[i] = 1.0e-4;
    ATOL[i] = 1.0e-3;
  }

  // Set min & max iteration step bounds
  STEPMIN = 0.0;
  STEPMAX = 900.0;

  //=========================================================================
  // Time loop
  //=========================================================================

  // Open MONITOR output file
  InitSaveData();

  // Print MONITOR data to stdout
  printf("\n%7s %7s   ", "done[%]", "Time[h]");
  for( i=0; i<NMONITOR; i++ ) printf( "%8s  ",   SPC_NAMES[MONITOR[i]] );
  for( i=0; i<NMASS;    i++ ) printf( "(%6s)  ", SMASS[i]              );

  // Loop over time
  TIME = TSTART;
  while (TIME <= TEND) {

    // Print species information
    GetMass( C, dval );
    printf("\n%6.1f%% %7.2f   ",
	   (TIME-TSTART)/(TEND-TSTART)*100.0, TIME/3600.0 );
    for( i = 0; i < NMONITOR; i++ )
      printf( "%9.3e  ", C[ MONITOR[i] ]/CFACTOR );
    for( i = 0; i < NMASS; i++ )
      printf( "%9.3e  ", dval[i]/CFACTOR );

    // Write to MONITOR output file
    SaveData();

    // Do the forward integration
    INTEGRATE( TIME, TIME+DT );
    TIME += DT;
  }

  //=========================================================================
  // Cleanup
  //=========================================================================

  // Close the MONITOR output file
  printf("\n");
  CloseSaveData();

  return 0;
}
