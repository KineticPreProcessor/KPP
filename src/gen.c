/******************************************************************************

  KPP - The Kinetic PreProcessor
        Builds simulation code for chemical kinetic systems

  Copyright (C) 1995-1996 Valeriu Damian and Adrian Sandu
  Copyright (C) 1997-2005 Adrian Sandu

  KPP is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation (http://www.gnu.org/copyleft/gpl.html); either version 2 of the
  License, or (at your option) any later version.

  KPP is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, consult http://www.gnu.org/copyleft/gpl.html or
  write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
  Boston, MA  02111-1307,  USA.

  Adrian Sandu
  Computer Science Department
  Virginia Polytechnic Institute and State University
  Blacksburg, VA 24060
  E-mail: sandu@cs.vt.edu

******************************************************************************/


#include <string.h>
#include <math.h>
#include "gdata.h"
#include "code.h"
#include "scan.h"

#define MAX_MONITOR 8
enum strutypes { PLAIN, LU };

int **structB;
int **structJ;
int **LUstructJ;

ICODE InlineCode[ INLINE_OPT ];

int NSPEC, NVAR, NVARP1, NVARACT, NFIX, NREACT;
int NVARST, NFIXST;
int C_DEFAULT, C;
int DC;
int ARP, JVRP, NJVRP, CROW_JVRP, IROW_JVRP, ICOL_JVRP;
int V, F, VAR, FIX;
int RCONST, RCT;
int Vdot, P_VAR, D_VAR, Aout;
int StoichNum;
int KR, A, BV, BR, IV, RR;
int JV, UV, JUV, JTUV, JVS;
int JR, UR, JUR, JRS;
int U1, U2, HU, HTU;
int X, XX, NTMPB;
int D2A, NTMPD2A, NHESS, HESS, IHESS_I, IHESS_J, IHESS_K;
int DDMTYPE;
int STOICM, NSTOICM, IROW_STOICM, ICOL_STOICM, CCOL_STOICM, CNEQN;
int IROW, ICOL, CROW, DIAG;
int LU_IROW, LU_ICOL, LU_CROW, LU_DIAG, CNVAR;
int LOOKAT, NLOOKAT, MONITOR, NMONITOR;
int NMASS, SMASS;
int SPC_NAMES, EQN_NAMES, FAM_NAMES;
int EQN_TAGS;
int NONZERO, LU_NONZERO;
int TIME, SUN, TEMP;
int TSTART, TEND, DT;
int ATOL, RTOL, STEPMIN, STEPMAX, CFACTOR;
int V_USER, CL;
int NMLCV, NMLCF, SCT, PROPENSITY, VOLUME, IRCT;
int FAM,NFAM;

int Jac_NZ, LU_Jac_NZ, nzr;

/* for autoreduce functionality - msl, hplin, 12/8/21 */
int DO_SLV, DO_JVS, DO_FUN, cLU_IROW, cLU_ICOL, cLU_DIAG, cLU_CROW;
int SPC_MAP, iSPC_MAP, RMV, JVS_MAP, RNVAR, cNONZERO, keepActive, keepSpcActive;

NODE *sum, *prod;
int real;
int nlookat;
int nmoni;
int ntrans;
int nmass;
char * CommonName;

int Hess_NZ, *iHess_i, *iHess_j, *iHess_k;
int nnz_stoicm;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
char * ascii(int x)
{
static char s[40];

  sprintf(s, "%d", x);
  return s;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
char * ascid(double x)
{
static char s[40];

  sprintf(s, "%12.6e", x);
  if (useDouble && (useLang==F77_LANG))
    s[strlen(s)-4] = 'd';
  if (useDouble && (useLang==F90_LANG))
    strncat( s, "_dp", 4 );
  return s;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
NODE * RConst( int n )
{
  switch( kr[n].type ) {
    case NUMBER:    return Const( kr[n].val.f );
    case PHOTO:
    case EXPRESION: return Elm( RCT, n );
  }
  return 0;
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void InitGen()
{
int i,j;

  NSPEC   = DefConst( "NSPEC",   INT, "Number of chemical species" );
  NVAR    = DefConst( "NVAR",    INT, "Number of Variable species" );
  NVARACT = DefConst( "NVARACT", INT, "Number of Active species" );
  NFIX    = DefConst( "NFIX",    INT, "Number of Fixed species" );
  NREACT  = DefConst( "NREACT",  INT, "Number of reactions" );
  NVARST  = DefConst( "NVARST",  INT, "Starting of variables in conc. vect." );
  NFIXST  = DefConst( "NFIXST",  INT, "Starting of fixed in conc. vect." );
  NFAM    = DefConst( "NFAM",    INT, "Number of Prod/Loss Families");
  NONZERO = DefConst( "NONZERO", INT, "Number of nonzero entries in Jacobian" );
  LU_NONZERO = DefConst( "LU_NONZERO", INT, "Number of nonzero entries in LU factoriz. of Jacobian" );
  CNVAR   = DefConst( "CNVAR",   INT, "(NVAR+1) Number of elements in compressed row format" );
  CNEQN   = DefConst( "CNEQN",   INT, "(NREACT+1) Number stoicm elements in compressed col format" );

  /* PI      = DefConst( "PI",     real, "Value of pi" ); */

//============================================================================
// MODIFICATION by Bob Yantosca (25 Apr 2022)
//
// For Fortran-90, declare VAR and FIX as POINTER variables, which will
// allow us to remove the non-threadsafe EQUIVALENCE statements.
//
  if ( useLang == F90_LANG ) {
    VAR = DefvElmP( "VAR", real,
                    "Concentrations of variable species (global)" );
    FIX = DefvElmP( "FIX", real,
                    "Concentrations of fixed species (global)" );
  } else {
    VAR = DefvElm( "VAR", real, -NVAR,
                   "Concentrations of variable species (global)" );
    FIX = DefvElm( "FIX", real, -NFIX,
                   "Concentrations of fixed species (global)" );
}
//============================================================================
  V = DefvElm( "V", real, -NVAR, "Concentrations of variable species (local)" );
  F = DefvElm( "F", real, -NFIX, "Concentrations of fixed species (local)" );

  V_USER = DefvElm( "V_USER", real, -NVAR, "Concentration of variable species in USER's order" );

  FAM    = DefvElm( "FAM", real, -NFAM, "Accumulated user-defined prod/loss families." );

  RCONST = DefvElm( "RCONST", real, -NREACT, "Rate constants (global)" );
  RCT    = DefvElm( "RCT",    real, -NREACT, "Rate constants (local)" );

  Vdot = DefvElm( "Vdot", real, -NVAR, "Time derivative of variable species concentrations" );
  P_VAR = DefvElm( "P_VAR", real, -NVAR, "Production term" );
  D_VAR = DefvElm( "D_VAR", real, -NVAR, "Destruction term" );
  StoichNum = DefmElm( "StoichNum", real, -NVAR, -NREACT, "Stoichiometric numbers" );


  JVS   = DefvElm( "JVS", real, -LU_NONZERO, "sparse Jacobian of variables" );

  JV    = DefmElm( "JV",  real, -NVAR, -NVAR, "full Jacobian of variables" );

  UV    = DefvElm( "UV",  real, -NVAR, "User vector for variables" );
  JUV   = DefvElm( "JUV", real, -NVAR, "Jacobian times user vector" );
  JTUV  = DefvElm( "JTUV",real, -NVAR, "Jacobian transposed times user vector" );

  X     = DefvElm( "X",  real, -NVAR, "Vector for variables" );
  XX    = DefvElm( "XX", real, -NVAR, "Vector for output variables" );

  TIME  = DefElm( "TIME", real, "Current integration time");
  SUN   = DefElm( "SUN", real, "Sunlight intensity between [0,1]");
  TEMP  = DefElm( "TEMP", real, "Temperature");
  TSTART = DefElm( "TSTART", real, "Integration start time");
  TEND   = DefElm( "TEND", real, "Integration end time");
  DT     = DefElm( "DT", real, "Integration step");

  A  = DefvElm( "A", real, -NREACT, "Rate for each equation" );

  ARP  = DefvElm( "ARP", real, -NREACT, "Reactant product in each equation" );
  NJVRP    = DefConst( "NJVRP",   INT, "Length of sparse Jacobian JVRP" );
  JVRP = DefvElm( "JVRP", real, -NJVRP, "d ARP(1:NREACT)/d VAR (1:NVAR)" );
  CROW_JVRP= DefvElm( "CROW_JVRP", INT, -CNEQN, "Beginning of rows in JVRP" );
  ICOL_JVRP= DefvElm( "ICOL_JVRP", INT, -NJVRP, "Column indices in JVRP" );
  IROW_JVRP= DefvElm( "IROW_JVRP", INT, -NJVRP, "Row indices in JVRP" );

  NTMPB  = DefConst( "NTMPB",   INT, "Length of Temporary Array B" );
  BV = DefvElm( "B", real, -NTMPB, "Temporary array" );

  NSTOICM      = DefConst("NSTOICM", INT, "Length of Sparse Stoichiometric Matrix" );
  STOICM       = DefvElm( "STOICM", real, -NSTOICM, "Stoichiometric Matrix in compressed column format" );
  IROW_STOICM  = DefvElm( "IROW_STOICM", INT, -NSTOICM, "Row indices in STOICM" );
  ICOL_STOICM  = DefvElm( "ICOL_STOICM", INT, -NSTOICM, "Column indices in STOICM" );
  CCOL_STOICM  = DefvElm( "CCOL_STOICM", INT, -CNEQN, "Beginning of columns in STOICM" );

  DDMTYPE      = DefElm( "DDMTYPE", INT, "DDM sensitivity w.r.t.: 0=init.val., 1=params" );

  NTMPD2A= DefConst( "NTMPD2A",   INT, "Length of Temporary Array D2A" );
  D2A    = DefvElm( "D2A", real, -NTMPD2A, "Second derivatives of equation rates" );
  NHESS = DefConst( "NHESS", INT, "Length of Sparse Hessian" );
  HESS  = DefvElm( "HESS", real, -NHESS, "Hessian of Var (i.e. the 3-tensor d Jac / d Var)" );
  IHESS_I  = DefvElm( "IHESS_I", INT, -NHESS, "Index i of Hessian element d^2 f_i/dv_j.dv_k" );
  IHESS_J  = DefvElm( "IHESS_J", INT, -NHESS, "Index j of Hessian element d^2 f_i/dv_j.dv_k" );
  IHESS_K  = DefvElm( "IHESS_K", INT, -NHESS, "Index k of Hessian element d^2 f_i/dv_j.dv_k" );
  U1  = DefvElm( "U1",  real, -NVAR, "User vector" );
  U2  = DefvElm( "U2",  real, -NVAR, "User vector" );
  HU  = DefvElm( "HU", real, -NVAR, "Hessian times user vectors: (Hess x U2) * U1 = [d (Jac*U1)/d Var] * U2" );
  HTU = DefvElm( "HTU", real, -NVAR, "Transposed Hessian times user vectors: (Hess x U2)^T * U1 = [d (Jac^T*U1)/d Var] * U2 " );

  KR = DefeElm( "KR", 0 );

  IROW  = DefvElm( "IROW", INT, -NONZERO, "Row indexes of the Jacobian of variables" );
  ICOL  = DefvElm( "ICOL", INT, -NONZERO, "Column indexes of the Jacobian of variables" );
  CROW  = DefvElm( "CROW", INT, -CNVAR, "Compressed row indexes of the Jacobian of variables" );
  DIAG  = DefvElm( "DIAG", INT, -CNVAR, "Diagonal indexes of the Jacobian of variables" );
  LU_IROW  = DefvElm( "LU_IROW", INT, -LU_NONZERO, "Row indexes of the LU Jacobian of variables" );
  LU_ICOL  = DefvElm( "LU_ICOL", INT, -LU_NONZERO, "Column indexes of the LU Jacobian of variables" );
  LU_CROW  = DefvElm( "LU_CROW", INT, -CNVAR, "Compressed row indexes of the LU Jacobian of variables" );
  LU_DIAG  = DefvElm( "LU_DIAG", INT, -CNVAR, "Diagonal indexes of the LU Jacobian of variables" );

  IV = DefeElm( "IV", 0 );

  C_DEFAULT = DefvElm( "C_DEFAULT", real, -NSPEC, "Default concentration for all species" );
//============================================================================
// MODIFICATION by Bob Yantosca (25 Apr 2002)
// For Fortran-90, declare C with the TARGET attribute.  This will allow
// VAR and FIX to point to C, which will let us remove the non-threadsafe
// EQUIVALENCE statements.
//
  if ( useLang == F90_LANG ) {
    C = DefvElmT( "C", real, -NSPEC, "Concentration of all species" );
  } else {
    C = DefvElm( "C", real, -NSPEC, "Concentration of all species" );
  }
//============================================================================
  CL = DefvElm( "CL", real, -NSPEC, "Concentration of all species (local)" );
  ATOL = DefvElm( "ATOL", real, -NVAR, "Absolute tolerance" );
  RTOL = DefvElm( "RTOL", real, -NVAR, "Relative tolerance" );

  STEPMIN  = DefElm( "STEPMIN", real, "Lower bound for integration step");
  STEPMAX  = DefElm( "STEPMAX", real, "Upper bound for integration step");

  NLOOKAT = DefConst( "NLOOKAT", INT, "Number of species to look at" );
  LOOKAT  = DefvElm( "LOOKAT", INT, -NLOOKAT, "Indexes of species to look at" );

  NMONITOR = DefConst( "NMONITOR", INT, "Number of species to monitor" );
  MONITOR  = DefvElm( "MONITOR", INT, -NMONITOR, "Indexes of species to monitor" );

  NMASS  = DefConst( "NMASS", INT, "Number of atoms to check mass balance" );
  SMASS  = DefvElm( "SMASS", STRING, -NMASS, "Names of atoms for mass balance" );

  EQN_TAGS    = DefvElm( "EQN_TAGS", STRING, -NREACT, "Equation tags" );
  EQN_NAMES  = DefvElm( "EQN_NAMES", DOUBLESTRING, -NREACT, "Equation names" );
  SPC_NAMES  = DefvElm( "SPC_NAMES", STRING, -NSPEC, "Names of chemical species" );
  FAM_NAMES  = DefvElm( "FAM_NAMES", STRING, -NFAM, "Names of chemical familes" );

  CFACTOR  = DefElm( "CFACTOR", real, "Conversion factor for concentration units");

  /* Autoreduction structures */
  NVARP1        = DefConst( "NVAR+1",  INT, "Number of Variable species + 1" );
  DO_JVS        = DefvElm("DO_JVS", LOGICAL, -LU_NONZERO, "");
  DO_SLV        = DefvElm("DO_SLV", LOGICAL, -NVARP1    , "");
  DO_FUN        = DefvElm("DO_FUN", LOGICAL, -NVAR      , "");
  cLU_IROW      = DefvElm("cLU_IROW", INT, -LU_NONZERO  , "");
  cLU_ICOL      = DefvElm("cLU_ICOL", INT, -LU_NONZERO  , "");
  cLU_CROW      = DefvElm("cLU_CROW", INT, -NVARP1  , "");
  cLU_DIAG      = DefvElm("cLU_DIAG", INT, -NVARP1  , "");
  JVS_MAP       = DefvElm("JVS_MAP", INT, -LU_NONZERO, "");
  SPC_MAP       = DefvElm("SPC_MAP", INT, -NVAR, "");
  iSPC_MAP      = DefvElm("iSPC_MAP", INT, -NVAR, "");
  RMV           = DefvElm("RMV", INT, -NVAR, "");
  RNVAR         = DefElm("rNVAR", INT, "");
  cNONZERO      = DefElm("cNONZERO", INT, "");
  keepSpcActive = DefvElm("KEEPSPCACTIVE", LOGICAL, -NVAR, "");
  keepActive    = DefElm("KEEPACTIVE", LOGICAL,"");

  Aout = DefvElmO( "Aout", real, -NREACT,
                   "Optional argument to return equation rate constants" );

  /* Elements of Stochastic simulation*/
  NMLCV = DefvElm( "NmlcV", INT, -NVAR, "No. molecules of variable species" );
  NMLCF = DefvElm( "NmlcF", INT, -NFIX, "No. molecules of fixed species" );
  SCT   = DefvElm( "SCT",  real, -NREACT, "Stochastic rate constants" );
  PROPENSITY  = DefvElm( "Prop",  real, -NREACT, "Propensity vector" );
  VOLUME = DefElm( "Volume", real, "Volume of the reaction container" );
  IRCT  = DefElm( "IRCT", INT, "Index of chemical reaction" );

  for ( i=0; i<EqnNr; i++ )
    for ( j=0; j<SpcNr; j++ )
      structB[i][j] = ( Stoich_Left[j][i] != 0 ) ? 1 : 0;

  /* Constant values are useful to declare vectors of this size */
  if (useDeclareValues) {
    varTable[ NSPEC ]   -> value  = max(SpcNr,1);
    varTable[ NVAR  ]   -> value  = max(VarNr,1);
    varTable[ NVARACT ] -> value  = max(VarActiveNr,1);
    varTable[ NFIX   ]  -> value  = max(FixNr,1);
    varTable[ NREACT ]  -> value  = max(EqnNr,1);
    varTable[ NVARST ]  -> value  = Index(0);
    varTable[ NFIXST ]  -> value  = Index(VarNr);
  }
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int NonZero( int stru, int start, int end,
             int *row, int *col, int *crow, int *diag )
{
int nElm;
int i,j;

  nElm = 0;
  for (i = 0; i < end-start; i++) {
    crow[i] = Index(nElm);
    for (j = 0; j < end-start; j++) {
      if( (i == j) || ( (stru) ? LUstructJ[i+start][j+start]
                               : structJ[i+start][j+start] ) ) {
        row[nElm] = Index(i);
        col[nElm] = Index(j);
        nElm++;
      }
      if( i == j ) {
        diag[i] = Index(nElm-1);
      }
    }
  }
  crow[i] = Index(nElm);
  diag[i] = Index(nElm);
  return nElm;
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateGData()
{
  if ( (useLang != C_LANG)&&(useLang != MATLAB_LANG) ) return;

  UseFile( driverFile );

  NewLines(1);

  GlobalDeclare( C );
  C_Inline("%s * %s = & %s[%d];", C_types[real],
            varTable[VAR]->name, varTable[C]->name, 0 );
  C_Inline("%s * %s = & %s[%d];", C_types[real],
            varTable[FIX]->name, varTable[C]->name, VarNr );


  GlobalDeclare( RCONST );
  GlobalDeclare( TIME );
  GlobalDeclare( SUN );
  GlobalDeclare( TEMP );
  GlobalDeclare( TSTART );
  GlobalDeclare( TEND );
  GlobalDeclare( DT );
  GlobalDeclare( ATOL );
  GlobalDeclare( RTOL );
  GlobalDeclare( STEPMIN );
  GlobalDeclare( STEPMAX );
  GlobalDeclare( CFACTOR );
  if (useStochastic)
      GlobalDeclare( VOLUME );

  MATLAB_Inline("  %s_Parameters;",rootFileName);
  MATLAB_Inline("  %s_Global_defs;",rootFileName);
  MATLAB_Inline("  %s_Sparse;",rootFileName);
  MATLAB_Inline("  %s_Monitor;",rootFileName);
  if (useJacSparse )
    MATLAB_Inline("  %s_JacobianSP;",rootFileName);
  if (useHessian )
    MATLAB_Inline("  %s_HessianSP;",rootFileName);
  if (useStoicmat )
    MATLAB_Inline("  %s_StoichiomSP;",rootFileName);

  NewLines(1);

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateMonitorData()
{
int i;
int  *lookat;
int  *moni;
char *snames[MAX_SPECIES];
int  *trans;
char *smass[MAX_ATOMS];
char *seqn[MAX_EQN];
char *sfam[MAX_FAMILIES];
char *bufeqn, *p;
int dim;


  /* Allocate local data structures */
  dim    = SpcNr+2;
  lookat = AllocIntegerVector( dim, "lookat in GenerateMonitorData");
  moni   = AllocIntegerVector( dim, "moni in GenerateMonitorData");
  trans  = AllocIntegerVector( dim, "trans in GenerateMonitorData");

  UseFile( monitorFile );

  F77_Inline("%6sBLOCK DATA MONITOR_DATA\n", " " );
  F77_Inline("%6sINCLUDE '%s_Parameters.h'", " ",rootFileName);
  F77_Inline("%6sINCLUDE '%s_Global.h'", " ",rootFileName);
  F77_Inline("%6sINTEGER i", " " );

  NewLines(1);

  for (i = 0; i < SpcNr; i++) {
      snames[i] = SpeciesTable[Code[i]].name;
  }
  InitDeclare( SPC_NAMES, SpcNr, (void*)snames );

  nlookat = 0;
  for (i = 0; i < SpcNr; i++)
    if ( SpeciesTable[Code[i]].lookat ) {
      lookat[nlookat] = Index(i);
      nlookat++;
    }

  if (useDeclareValues)
        varTable[ NLOOKAT ]  -> value  = max(nlookat,1);
  InitDeclare( LOOKAT, nlookat, (void*)lookat );

  nmoni = 0;
  for (i = 0; i < SpcNr; i++)
    if ( SpeciesTable[Code[i]].moni ) {
      moni[nmoni] = Index(i);
      nmoni++;
    }

  if( nmoni > MAX_MONITOR ) {
    Warning( "%d species to monitorize. Too many, keeping %d.",
             nmoni, MAX_MONITOR );
    nmoni = MAX_MONITOR;
  }

  if (useDeclareValues)
        varTable[ NMONITOR ] -> value  = max(nmoni,1);
  InitDeclare( MONITOR, nmoni, (void*)moni );

  ntrans = 0;
  for (i = 0; i < SpcNr; i++)
    if ( SpeciesTable[Code[i]].trans ) {
      trans[ntrans] = Index(i);
      ntrans++;
    }

  nmass = 0;
  for (i = 0; i < AtomNr; i++)
    if ( AtomTable[i].masscheck ) {
      smass[nmass] = AtomTable[i].name;
      nmass++;
    }
  if (useDeclareValues)
        varTable[ NMASS ]    -> value  = max(nmass,1);
  InitDeclare( SMASS, nmass, (void*)smass );

  if ( (bufeqn = (char*)malloc(MAX_EQNLEN*EqnNr+2))==NULL ) {
    FatalError(-30,"GenerateMonitorData: Cannot allocate bufeqn (%d chars)",
                    MAX_EQNLEN*EqnNr);
  }

  p = bufeqn;
  for (i = 0; i < EqnNr; i++) {
    EqnString(i, p);
    seqn[i] = p;
    p += MAX_EQNLEN;
  }
  InitDeclare( EQN_NAMES, EqnNr, (void*)seqn );

  free( bufeqn );

  if (useEqntags==1) {
    for (i = 0; i < EqnNr; i++) {
      seqn[i] = kr[i].label;
    }
    InitDeclare( EQN_TAGS, EqnNr, (void*)seqn );
  }

  for (i = 0; i < FamilyNr; i++) {
    sfam[i] = FamilyTable[ i ].name;
  }
  InitDeclare( FAM_NAMES, FamilyNr, (void*)sfam );

  NewLines(1);
  WriteComment("INLINED global variables");

  switch( useLang ) {
    case C_LANG:   bprintf( InlineCode[ C_DATA ].code );
                 break;
    case F77_LANG: bprintf( InlineCode[ F77_DATA ].code );
                 break;
    case F90_LANG: bprintf( InlineCode[ F90_DATA ].code );
                 break;
    case MATLAB_LANG: bprintf( InlineCode[ MATLAB_DATA ].code );
                 break;
  }
  FlushBuf();

  NewLines(1);
  WriteComment("End INLINED global variables");
  NewLines(1);

  F77_Inline( "%6sEND\n\n", " " );

  /* Free local data structures */
  free(lookat); free(moni); free(trans);

}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateJacobianSparseData()
{
int* irow;
int* icol;
int* crow;
int* diag;
int dim;

  if( !useJacSparse ) return;

  /* Allocate local arrays */
  dim=MAX_SPECIES;
  irow = AllocIntegerVector( dim*dim, "irow in GenerateJacobianSparseData" );
  icol = AllocIntegerVector( dim*dim, "icol in GenerateJacobianSparseData" );
  crow = AllocIntegerVector( dim,     "crow in GenerateJacobianSparseData" );
  diag = AllocIntegerVector( dim,     "diag in GenerateJacobianSparseData" );

  UseFile( sparse_jacFile );

  NewLines(1);
  WriteComment("Sparse Jacobian Data");
  NewLines(1);

  F77_Inline("%6sBLOCK DATA JACOBIAN_SPARSE_DATA\n", " " );
  F77_Inline("%6sINCLUDE '%s_Sparse.h'", " ",rootFileName);
  F77_Inline("%6sINTEGER i"," ");


  Jac_NZ = NonZero( PLAIN, 0, VarNr, irow, icol, crow, diag );
  LU_Jac_NZ = NonZero( LU, 0, VarNr, irow, icol, crow, diag );
  if (useDeclareValues) {
     varTable[NONZERO] -> value = Jac_NZ;
     varTable[LU_NONZERO] -> value = LU_Jac_NZ;
  }

  switch (useJacobian) {
  case JAC_ROW:
    Jac_NZ = NonZero( PLAIN, 0, VarNr, irow, icol, crow, diag );
    InitDeclare( IROW, Jac_NZ, (void*)irow );
    InitDeclare( ICOL, Jac_NZ, (void*)icol );
    InitDeclare( CROW, VarNr+1, (void*)crow );
    InitDeclare( DIAG, VarNr+1, (void*)diag );
    break;
  case JAC_LU_ROW:
    LU_Jac_NZ = NonZero( LU, 0, VarNr, irow, icol, crow, diag );
    InitDeclare( LU_IROW, LU_Jac_NZ, (void*)irow );
    InitDeclare( LU_ICOL, LU_Jac_NZ, (void*)icol );
    InitDeclare( LU_CROW, VarNr+1, (void*)crow );
    InitDeclare( LU_DIAG, VarNr+1, (void*)diag );
  }
  NewLines(1);
  F77_Inline( "%6sEND\n\n", " " );

  /* Free local arrays */
  free(irow); free(icol); free(crow); free(diag);
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateJacobianSparseHeader()
{
  UseFile( sparse_dataFile );

  CommonName = "SDATA";

  NewLines(1);
  WriteComment(" ----------> Sparse Jacobian Data");
  NewLines(1);

  switch (useJacobian) {
  case JAC_ROW:
    ExternDeclare( IROW );
    ExternDeclare( ICOL );
    ExternDeclare( CROW );
    ExternDeclare( DIAG );
     break;
  case JAC_LU_ROW:
    ExternDeclare( LU_IROW );
    ExternDeclare( LU_ICOL );
    ExternDeclare( LU_CROW );
    ExternDeclare( LU_DIAG );
  }

  NewLines(1);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*  mz_rs_20160201+ */
void GenerateStoichNum()
{
  int i, j;
  int CalcStoichNum;

  if( VarNr == 0 ) return;
  if (useLang != MATLAB_LANG)  /* Matlab generates an additional file per function */
    UseFile( functionFile );
  CalcStoichNum = DefFnc( "CalcStoichNum", 1, "calculate stoichiometric numbers");
  FunctionBegin( CalcStoichNum, StoichNum );
  F90_Inline("  StoichNum(:,:) = 0.");
  for (i = 0; i < VarNr; i++) {
    for (j = 0; j < EqnNr; j++) {
      if ( Stoich[i][j] != 0 )
        Assign( Elm( StoichNum, i, j ), Const( Stoich[i][j] ));
    }
  }
  FunctionEnd( CalcStoichNum );
  FreeVariable( CalcStoichNum );
}
/*  mz_rs_20160201- */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateFun( int z_useAggregate )
{
int i, j, k;
int used;
int l, m;
int F_VAR, FSPLIT_VAR;

  if( VarNr == 0 ) return;

  if (useLang != MATLAB_LANG)  /* Matlab generates an additional file per function */
       UseFile( functionFile );

  // Note: When changing the FunctionBegin declarations below,
  // the number of arguments minus one must be changed in DefFnc here
  // as well. (hplin, 7/6/22)
  F_VAR = DefFnc( "Fun", 5,
                  "time derivatives of variables - Aggregate form");

  FSPLIT_VAR = DefFnc( "Fun_SPLIT", 7,
                       "time derivatives of variables - Split form");


  // We have added the capability to return equation rates and the
  // time derivative of variable species from Fun via optional arguments
  // Aout and VdotOut (when z_useAggregate=1)
  //   -- Bob Yantosca (03 May 2022)
  //
  // Vdotout functionality can be accomplished using Vdot (hplin, 7/6/22)
  if( z_useAggregate ) {
    FunctionBegin( F_VAR, V, F, RCT, Vdot, Aout );
  }
  else {
    FunctionBegin( FSPLIT_VAR, V, F, RCT, Vdot, P_VAR, D_VAR, Aout );
  }

  if ( (useLang==MATLAB_LANG)&&(!z_useAggregate) )
     printf("\nWarning: in the function definition move P_VAR to output vars\n");

  if ( useLang!=F90_LANG ) { /* A is a module variable in F90 */
    NewLines(1);
    WriteComment("Local variables");
    Declare( A );
    if ( useLang==F77_LANG ) WriteOMPThreadPrivate("A");
  }
  if (useLang==MATLAB_LANG) {
    MATLAB_Inline("A=zeros(1,length(RCT));");
    MATLAB_Inline("Vdot=zeros(1,length(V));");
  }
  NewLines(1);
  WriteComment("Computation of equation rates");

  for(j=0; j<EqnNr; j++) {
    used = 0;
    if( z_useAggregate ) {
      for (i = 0; i < VarNr; i++)
        if ( Stoich[i][j] != 0 ) {
          used = 1;
          break;
        }
    } else {
      for (i = 0; i < VarNr; i++)
        if ( Stoich_Right[i][j] != 0 ) {
          used = 1;
          break;
        }
    }

    if ( used ) {
      prod = RConst( j );
      for (i = 0; i < VarNr; i++)
        for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
          prod = Mul( prod, Elm( V, i ) );
      for ( ; i < SpcNr; i++)
        for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
          prod = Mul( prod, Elm( F, i - VarNr ) );
      Assign( Elm( A, j ), prod );
    }
  }

  // Add code to return equation rates via optional argument Aout
  //   -- Bob Yantosca (29 Apr 2022)
  NewLines(1);
  if ( useLang == F90_LANG ) {
    fprintf(functionFile, "  !### Use Aout to return equation rates\n");
    fprintf(functionFile, "  IF ( PRESENT( Aout ) ) Aout = A\n");
  }

  if( z_useAggregate ) {
    NewLines(1);
    WriteComment("Aggregate function");

    for (i = 0; i < VarNr; i++) {
      sum = Const(0);
      for (j = 0; j < EqnNr; j++) {
        sum = Add( sum, Mul( Const( Stoich[i][j] ), Elm( A, j ) ) );
      }
      if ( doAutoReduce ) F90_Inline("IF (DO_FUN(%d)) &",i+1);
      Assign( Elm( Vdot, i ), sum );
    }

  } else {

    NewLines(1);
    WriteComment("Production function");

    for (i = 0; i < VarNr; i++) {
      sum = Const(0);
      for (j = 0; j < EqnNr; j++)
        sum = Add( sum, Mul( Const( Stoich_Right[i][j] ), Elm( A, j ) ) );
      if( doAutoReduce ) F90_Inline("IF (DO_FUN(%d)) &",i+1);
      Assign( Elm( P_VAR, i ), sum );
    }

    NewLines(1);
    WriteComment("Destruction function");

    for (i = 0; i < VarNr; i++) {
      sum = Const(0);
      for(j=0; j<EqnNr; j++) {
        if ( Stoich_Left[i][j] == 0 ) continue;
        prod = Mul( RConst( j ), Const( Stoich_Left[i][j] ) );
        for (l = 0; l < VarNr; l++) {
          m=(int)Stoich_Left[l][j] - (l==i);
          for (k = 1; k <= m; k++ )
            prod = Mul( prod, Elm( V, l ) );
        }
        for ( ; l < SpcNr; l++)
          for (k = 1; k <= (int)Stoich_Left[l][j]; k++ )
            prod = Mul( prod, Elm( F, l - VarNr  ) );
        sum = Add( sum, prod );
      }
      if( doAutoReduce ) F90_Inline("IF (DO_FUN(%d)) &",i+1);
      Assign( Elm( D_VAR, i ), sum );
    }

    // Add code to calculate Vdot also for split function:
    NewLines(1);
    if ( useLang == F90_LANG ) {
      // For F90, we can do operations on arrays:
      fprintf(functionFile, "  Vdot = P_VAR - D_VAR*V\n");
    } else if ( useLang == C_LANG ) {
      // For C, we have to explicitly loop over elements:
      fprintf(functionFile, "  for( int n=0; n<NVAR; n++ ) {\n");
      fprintf(functionFile, "    Vdot[n] = P_VAR[n] - D_VAR[n]*V[n];\n");
      fprintf(functionFile, "  }\n");
    }
  }

  if( !z_useAggregate )
    MATLAB_Inline("\n   P_VAR = P_VAR(:);\n   D_VAR = D_VAR(:);\n");

  if( z_useAggregate )
    FunctionEnd( F_VAR );
  else
    FunctionEnd( FSPLIT_VAR );

  FreeVariable( F_VAR );
  FreeVariable( FSPLIT_VAR );

  /** hplin 4/10/22 add FUN_Split2 which overrides DO_FUN. only used for 1st order-autoreduce **/
  if(!z_useAggregate && doAutoReduce) {
    /* add a copy of FSPLIT that does not react to DO_FUN() */
    fprintf(functionFile, "! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    fprintf(functionFile, "!\n");
    fprintf(functionFile, "! Fun_SPLITF - time derivatives of variables - Split form\n");
    fprintf(functionFile, "!  same as Fun_Split, but does not react to DO_FUN.\n");
    fprintf(functionFile, "!   Arguments :\n");
    fprintf(functionFile, "!      V         - Concentrations of variable species (local)\n");
    fprintf(functionFile, "!      F         - Concentrations of fixed species (local)\n");
    fprintf(functionFile, "!      RCT       - Rate constants (local)\n");
    fprintf(functionFile, "!      P_VAR     - Production term\n");
    fprintf(functionFile, "!      D_VAR     - Destruction term\n");
    fprintf(functionFile, "!      Aout      - Array to return rxn rates for diagnostics (OPTIONAL)\n");
    fprintf(functionFile, "! Haipeng Lin (hplin) - Apr 10 2022 KPP 3.0.0-AR\n");
    fprintf(functionFile, "!\n");
    fprintf(functionFile, "! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
    fprintf(functionFile, "SUBROUTINE Fun_SPLITF ( V, F, RCT, P_VAR, D_VAR, Aout )\n\n");
    fprintf(functionFile, "! V - Concentrations of variable species (local)\n");
    fprintf(functionFile, "  REAL(kind=dp) :: V(NVAR)\n");
    fprintf(functionFile, "! F - Concentrations of fixed species (local)\n");
    fprintf(functionFile, "  REAL(kind=dp) :: F(NFIX)\n");
    fprintf(functionFile, "! RCT - Rate constants (local)\n");
    fprintf(functionFile, "  REAL(kind=dp) :: RCT(NREACT)\n");
    fprintf(functionFile, "! P_VAR - Production term\n");
    fprintf(functionFile, "  REAL(kind=dp) :: P_VAR(NVAR)\n");
    fprintf(functionFile, "! D_VAR - Destruction term\n");
    fprintf(functionFile, "  REAL(kind=dp) :: D_VAR(NVAR)\n");
    fprintf(functionFile, "!### Aout - Array for returning KPP reaction rates for diagnostics\n");
    fprintf(functionFile, "  REAL(kind=dp), OPTIONAL :: Aout(NREACT)\n");

    NewLines(1);
    WriteComment("Computation of equation rates");

    for(j=0; j<EqnNr; j++) {
      used = 0;
      for (i = 0; i < VarNr; i++) {
        if ( Stoich_Right[i][j] != 0 ) {
          used = 1;
          break;
        }
      }

      if ( used ) {
        prod = RConst( j );
        for (i = 0; i < VarNr; i++)
          for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
            prod = Mul( prod, Elm( V, i ) );
        for ( ; i < SpcNr; i++)
          for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
            prod = Mul( prod, Elm( F, i - VarNr ) );
        Assign( Elm( A, j ), prod );
      }
    }

    if ( useLang == F90_LANG ) {
      fprintf(functionFile, "\n\n!### KPP 2.3.0_gc, Bob Yantosca (11 Feb 2021)\n");
      fprintf(functionFile, "!### Use Aout to return reaction rates\n");
      fprintf(functionFile, "  IF ( PRESENT( Aout ) ) Aout = A\n\n");
    }
    NewLines(1);
    WriteComment("Production function");

    for (i = 0; i < VarNr; i++) {
      sum = Const(0);
      for (j = 0; j < EqnNr; j++)
        sum = Add( sum, Mul( Const( Stoich_Right[i][j] ), Elm( A, j ) ) );
      /* if( doAutoReduce ) F90_Inline("IF (DO_FUN(%d)) &",i+1); */ /* purpose of existence */
      Assign( Elm( P_VAR, i ), sum );
    }

    NewLines(1);
    WriteComment("Destruction function");

    for (i = 0; i < VarNr; i++) {
      sum = Const(0);
      for(j=0; j<EqnNr; j++) {
        if ( Stoich_Left[i][j] == 0 ) continue;
        prod = Mul( RConst( j ), Const( Stoich_Left[i][j] ) );
        for (l = 0; l < VarNr; l++) {
          m=(int)Stoich_Left[l][j] - (l==i);
          for (k = 1; k <= m; k++ )
            prod = Mul( prod, Elm( V, l ) );
        }
        for ( ; l < SpcNr; l++)
          for (k = 1; k <= (int)Stoich_Left[l][j]; k++ )
            prod = Mul( prod, Elm( F, l - VarNr  ) );
        sum = Add( sum, prod );
      }
      Assign( Elm( D_VAR, i ), sum );
    }

    fprintf(functionFile, "END SUBROUTINE Fun_SPLITF\n");
  } // end if doAutoReduce
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateStochastic()
{
int i, j, k, m, n, jnr;
int used;
int F_VAR;

  if( VarNr == 0 ) return;

  if (useLang != MATLAB_LANG)  /* Matlab generates an additional file per function */
       UseFile( stochasticFile );

  /* ~~~~~~~> 1. PROPENSITY FUNCTION  ~~~~~~~~~~~~ */
  F_VAR = DefFnc( "Propensity", 4, "Propensity function");
  FunctionBegin( F_VAR, NMLCV, NMLCF, SCT, PROPENSITY );

  if ( (useLang==MATLAB_LANG)&&(!useAggregate) )
     printf("\nWarning: in the function definition move P_VAR to output vars\n");

  NewLines(1);

  for(j=0; j<EqnNr; j++) {
    used = 0;
    for (i = 0; i < VarNr; i++)
      if ( Stoich[i][j] != 0 ) {
        used = 1;
        break;
      }
    if ( used ) {
      prod = Elm( SCT, j );
      for (i = 0; i < VarNr; i++)
        for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
          if (k==1)
             prod = Mul( prod, Elm( NMLCV, i ) );
          else
             prod = Mul( prod, Add( Elm( NMLCV, i ), Const(-k+1) ) );

      for ( ; i < SpcNr; i++)
        for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
          if (k==1)
             prod = Mul( prod, Elm( NMLCF, i - VarNr ) );
          else
             prod = Mul( prod, Add( Elm( NMLCF, i - VarNr ), Const(-k+1) ) );
      Assign( Elm( PROPENSITY, j ), prod );
    } /* if used */
  } /* for j */

  MATLAB_Inline("\n   Prop = Prop(:);\n");
  FunctionEnd( F_VAR );
  FreeVariable( F_VAR );

  /* ~~~~~~~> 2. RATE CONVERSION ~~~~~~~~~~~~ */
  F_VAR = DefFnc( "StochasticRates", 3, "Convert deterministic rates to stochastic");
  FunctionBegin( F_VAR, RCT, VOLUME, SCT );
  WriteComment("No. of molecules = Concentration x Volume");
  WriteComment("For a reaction with k reactants:");
  WriteComment("         RCT [ (molec/Volume)^(1-k) * sec^(-1) ]");
  WriteComment("         SCT [ (molec)^(1-k) * sec^(-1) ] = RCT*Volume^(k-1)");
  WriteComment("For p molecules of the same type:  SCT = SCT/(p!)");

  NewLines(1);

  for(j=0; j<EqnNr; j++) {
      prod = Elm( RCT, j );
      m = 0;
      for (i = 0; i < SpcNr; i++)
        m += (int)Stoich_Left[i][j];
      for ( i=2 ; i <= m; i++)
         prod = Mul( prod, Elm(VOLUME, 1) );
      n = 1;
      for (i = 0; i < SpcNr; i++)
        for (k = 2; k <= (int)Stoich_Left[i][j]; k++ )
          n *= k;
      prod = Div( prod, Const( n ) );
      Assign( Elm( SCT, j ), prod );
  } /* for j */

  MATLAB_Inline("\n   SCT = SCT(:);\n");
  FunctionEnd( F_VAR );
  FreeVariable( F_VAR );

  /* ~~~~~~~> 3. THE CHANGE IN NUMBER OF MOLECULES */
  if (useLang == MATLAB_LANG) {
    F_VAR = DefFnc( "MoleculeChange", 3, "Change in the number of molecules");
    FunctionBegin( F_VAR, IRCT, NMLCV, NMLCV );
  } else {
    F_VAR = DefFnc( "MoleculeChange", 2, "Change in the number of molecules");
    FunctionBegin( F_VAR, IRCT, NMLCV );
  }

  NewLines(1);

  F90_Inline("\n SELECT CASE (IRCT)\n");
  C_Inline  ("\n switch (IRCT) { \n");
  MATLAB_Inline("\n switch (IRCT) \n");
  for(j=0; j<EqnNr; j++) {
      jnr = j+1;
      if (j==0)
        F77_Inline("\n%6sIF (IRCT.EQ.%d) THEN"," ",jnr);
      else
        F77_Inline("\n%6sELSEIF (IRCT.EQ.%d) THEN "," ",jnr);
      F90_Inline("\n  CASE (%d) ",jnr);
      C_Inline("\n  case %d: ",jnr);
      MATLAB_Inline("\n case %d, ",jnr);
      for (i = 0; i < VarNr; i++) {
        if ( Stoich_Left[i][j] || Stoich_Right[i][j] )
          Assign( Elm( NMLCV, i ), Add(Elm( NMLCV, i ),
               Const(Stoich_Right[i][j]-Stoich_Left[i][j])) );
      } /* for i */
      C_Inline("  break;",j);
  } /* for j */
  F77_Inline("\n%6sEND IF ! n\n"," ");
  F90_Inline("\n END SELECT\n");
  C_Inline("\n } /* switch (IRCT) */ \n");
  MATLAB_Inline("\n end\n");

  FunctionEnd( F_VAR );
  FreeVariable( F_VAR );

}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateReactantProd()
{
int i, j, k;
int F_STOIC;

  if( VarNr == 0 ) return;

  UseFile( stoichiomFile );

  F_STOIC = DefFnc( "ReactantProd",3, "Reactant Products in each equation");

  FunctionBegin( F_STOIC, V, F, ARP );

  NewLines(1);
  WriteComment("Reactant Products in each equation are useful in the");
  WriteComment("    stoichiometric formulation of mass action law");

  for(j=0; j<EqnNr; j++) {

      prod = Const( 1 );
      for (i = 0; i < VarNr; i++)
        for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
          prod = Mul( prod, Elm( V, i ) );
      for ( ; i < SpcNr; i++)
        for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
          prod = Mul( prod, Elm( F, i - VarNr  ) );
      Assign( Elm( ARP, j ), prod );

  } /* for j EqnNr */

  FunctionEnd(  F_STOIC );
  FreeVariable( F_STOIC );
}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateJacReactantProd()
{
int i, j, k, l, m, JVRP_NZ, newrow;
int F_STOIC;
int crow_JVRP[EqnNr], icol_JVRP[EqnNr*VarNr];
int irow_JVRP[EqnNr*VarNr];

  if( VarNr == 0 ) return;

  if (useDeclareValues) {
    JVRP_NZ = -1;
    for ( i=0; i<EqnNr; i++ )
      for ( j=0; j<VarNr; j++ )
         if ( Stoich_Left[j][i] != 0 ) JVRP_NZ++;
    varTable[ NJVRP ]  -> value  = JVRP_NZ + 1;
  }

  UseFile( stoichiomFile );

  F_STOIC = DefFnc( "JacReactantProd",3, "Jacobian of Reactant Products vector");

  FunctionBegin( F_STOIC, V, F, JVRP );


  NewLines(1);
  WriteComment("Reactant Products in each equation are useful in the");
  WriteComment("   stoichiometric formulation of mass action law");
  WriteComment("Below we compute the Jacobian of the Reactant Products vector");
  WriteComment("   w.r.t. variable species: d ARP(1:NREACT) / d Var(1:NVAR)");

  NewLines(1);

  JVRP_NZ = -1;
  for ( i=0; i<EqnNr; i++ ) {
    newrow = 0;
    crow_JVRP[i] = JVRP_NZ+1;
    for ( j=0; j<VarNr; j++ ) {
      if ( Stoich_Left[j][i] != 0 ) {
        JVRP_NZ++;
        icol_JVRP[JVRP_NZ] = j;
        irow_JVRP[JVRP_NZ] = i;
        if ( newrow == 0 ) { /* a new row begins here */
          crow_JVRP[i] = JVRP_NZ;
          newrow = 1;
        }
        prod = Const( Stoich_Left[j][i] ) ;
        for (l = 0; l < VarNr; l++) {
          m = (int)Stoich_Left[l][i] - (l==j);
          for (k = 1; k <= m; k++ )
            prod = Mul( prod, Elm( V, l ) );
        }
        for ( ; l < SpcNr; l++)
          for (k = 1; k <= (int)Stoich_Left[l][i]; k++ )
            prod = Mul( prod, Elm( F, l - VarNr ) );
        /* Comment the B */
        WriteComment("JVRP(%d) = dARP(%d)/dV(%d)",Index(JVRP_NZ),Index(i),Index(j));
        Assign( Elm( JVRP, JVRP_NZ ), prod );
      }
    }
  }
  crow_JVRP[EqnNr] = JVRP_NZ + 1;

  FunctionEnd(  F_STOIC );
  FreeVariable( F_STOIC );


  UseFile( sparse_stoicmFile );
  NewLines(1);
  WriteComment("Row-compressed sparse data for the Jacobian of reaction products JVRP");
  F77_Inline("%6sBLOCK DATA JVRP_SPARSE_DATA\n", " " );
  F77_Inline("%6sINCLUDE '%s_Sparse.h'", " ", rootFileName);
  F77_Inline("%6sINTEGER i", " ");
  if( (useLang==F77_LANG)||(useLang==F90_LANG) ) {
        for (k=0; k<JVRP_NZ+1; k++) {
           irow_JVRP[k]++;
           icol_JVRP[k]++;
           }
        for (k=0; k<EqnNr+1; k++)
           crow_JVRP[k]++;
  }
  InitDeclare( CROW_JVRP, EqnNr+1, (void*)crow_JVRP );
  InitDeclare( ICOL_JVRP, JVRP_NZ + 1, (void*)icol_JVRP );
  InitDeclare( IROW_JVRP, JVRP_NZ + 1, (void*)irow_JVRP );
  NewLines(1);
  F77_Inline( "%6sEND\n\n", " " );
  NewLines(1);


  UseFile( param_headerFile );
  CommonName = "GDATA";
  NewLines(1);
  DeclareConstant( NJVRP,   ascii( JVRP_NZ + 1 ) );

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateJac()
{
int i,j,k,l,m;
int nElm, nonzeros_B;
int Jac_SP, Jac;

  if( VarNr == 0 ) return;
  if (useJacobian == JAC_OFF) return;

  if (useLang != MATLAB_LANG)  /* Matlab generates an additional file per function */
       UseFile( jacobianFile );

  Jac_SP  = DefFnc( "Jac_SP", 4,
                  "the Jacobian of Variables in sparse matrix representation");
  Jac     = DefFnc( "Jac", 4, "the Jacobian of Variables");

  if( useJacSparse )
    FunctionBegin( Jac_SP, V, F, RCT, JVS );
  else
    FunctionBegin( Jac, V, F, RCT, JV );

  if (useLang == MATLAB_LANG) {
    switch (useJacobian) {
    case JAC_ROW:
      ExternDeclare(IROW); ExternDeclare(ICOL);
      break;
    case JAC_LU_ROW:
      ExternDeclare(LU_IROW); ExternDeclare(LU_ICOL);
      break;
    }
  }

  /* Each nonzero entry of B now counts its rank */
  nonzeros_B = 0;
  for ( i=0; i<EqnNr; i++ )
    for ( j=0; j<SpcNr; j++ )
       if ( structB[i][j] != 0 ) {
         nonzeros_B++;
         structB[i][j] = nonzeros_B;
         }

  if ( (useLang==C_LANG)||(useLang==F77_LANG)||(useLang==F90_LANG) ) {
    NewLines(1);
    WriteComment("Local variables");
    varTable[ NTMPB ] -> value = nonzeros_B;
    Declare( BV );
  }
  if (useLang == MATLAB_LANG) {
    MATLAB_Inline("B=zeros(1,%d);",nonzeros_B);
    MATLAB_Inline("JVS=zeros(1,length(LU_IROW));");
  }

  NewLines(1);

  for ( i=0; i<EqnNr; i++ ) {
    for ( j=0; j<VarNr; j++ ) {
      if ( Stoich_Left[j][i] != 0 ) {
        prod = Mul( RConst( i ), Const( Stoich_Left[j][i] ) );
        for (l = 0; l < VarNr; l++) {
          m = (int)Stoich_Left[l][i] - (l==j);
          for (k = 1; k <= m; k++ )
            prod = Mul( prod, Elm( V, l ) );
        }
        for ( ; l < SpcNr; l++)
          for (k = 1; k <= (int)Stoich_Left[l][i]; k++ )
            prod = Mul( prod, Elm( F, l - VarNr ) );
        /* Comment the B */
        WriteComment("B(%d) = dA(%d)/dV(%d)",Index(structB[i][j]-1),Index(i),Index(j));
        Assign( Elm( BV, structB[i][j]-1 ), prod );
      }
    }
  }

  nElm = 0;
  NewLines(1);
  WriteComment("Construct the Jacobian terms from B's");

  if ( useJacSparse ) {
  for (i = 0; i < VarNr; i++) {
    for (j = 0; j < VarNr; j++) {
      if( LUstructJ[i][j] ) {
        sum = Const(0);
        for (k = 0; k < EqnNr; k++) {
          if( Stoich[i][k]*structB[k][j] != 0 )
            sum = Add( sum, Mul( Const( Stoich[i][k] ), Elm( BV, structB[k][j]-1 ) ) );
        }
        /* Comment the B */
         if( doAutoReduce ) F90_Inline("IF (DO_JVS(%d)) &",nElm+1);
         WriteComment("JVS(%d) = Jac_FULL(%d,%d)",
                  Index(nElm),Index(i),Index(j));
         Assign( Elm( JVS, nElm ), sum );
         nElm++;
      } else {
        if( i == j ) {
          Assign( Elm( JVS, nElm ), Const(0) );
          nElm++;
        }
      }
    }
  }
  } else { /* full Jacobian */
  for (i = 0; i < VarNr; i++) {
    for (j = 0; j < VarNr; j++) {
      if( structJ[i][j] ) {
        sum = Const(0);
        for (k = 0; k < EqnNr; k++) {
          if( Stoich[i][k]*structB[k][j] != 0 )
            sum = Add( sum, Mul( Const( Stoich[i][k] ), Elm( BV, structB[k][j]-1 ) ) );
        }
        Assign( Elm( JV, i, j ), sum );
      }
    }
  }
  }

  if (useLang == MATLAB_LANG) {
    switch (useJacobian) {
    case JAC_ROW:
      MATLAB_Inline("\n   JVS = sparse(IROW,ICOL,JVS,%d,%d);\n",VarNr,VarNr);
      break;
    case JAC_LU_ROW:
      MATLAB_Inline("\n   JVS = sparse(LU_IROW,LU_ICOL,JVS,%d,%d);\n",VarNr,VarNr);
      break;
    }
  }

  if( useJacSparse )
    FunctionEnd( Jac_SP );
  else
    FunctionEnd( Jac );

  FreeVariable( Jac_SP );
  FreeVariable( Jac );

}






/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateHessian()
/* Unlike Hess, this function deffers the sparse Data structure generation */
{
int i, j, k;
int m, i1, i2, nElm;
int F_Hess, F_Hess_VEC, F_HessTR_VEC;
int *coeff_j, *coeff_i1, *coeff_i2;
int Djv_isElm;

  if ( VarNr == 0 ) return;

  if (useLang != MATLAB_LANG)  /* Matlab generates an additional file per function */
      UseFile( hessianFile );

/*  Calculate the number of nonzero terms of the form d^2 A(j)/ ( d v(i1) d v(i2) )*/
  nElm = 0;
  for(j=0; j<EqnNr; j++)
    for (i1 = 0; i1 < VarNr; i1++)
      for (i2 = i1; i2 < VarNr; i2++)
        if (i1==i2) {
          if (Stoich_Left[i1][j]>=2)
            nElm++;
        } else {  /* i1 != i2 */
          if ( (Stoich_Left[i1][j]>=1)&&(Stoich_Left[i2][j]>=1) )
            nElm++;
        }

/* Allocate temporary index arrays */
  coeff_j  = AllocIntegerVector(nElm, "coeff_j  in GenerateHess");
  coeff_i1 = AllocIntegerVector(nElm, "coeff_i1 in GenerateHess");
  coeff_i2 = AllocIntegerVector(nElm, "coeff_i2 in GenerateHess");

/*  Fill in temporary index arrays */
  nElm = 0;
  for(j=0; j<EqnNr; j++)
    for (i1 = 0; i1 < VarNr; i1++)
      for (i2 = i1; i2 < VarNr; i2++)
        if (i1==i2) {
          if (Stoich_Left[i1][j]>=2) {
            coeff_j[nElm] = j; coeff_i1[nElm] = i1; coeff_i2[nElm] = i2;
            nElm++;
            }
        } else {  /* i1 != i2 */
          if ( (Stoich_Left[i1][j]>=1)&&(Stoich_Left[i2][j]>=1) ) {
            coeff_j[nElm] = j; coeff_i1[nElm] = i1; coeff_i2[nElm] = i2;
            nElm++;
            }
        }
/*  Number of nonzero terms of the form d^2 f(i)/ ( d v(i1) d v(i2) ) */
    Hess_NZ = 0;
    for (i = 0; i < VarNr; i++)
      for (i1 = 0; i1 < VarNr; i1++)
         for (i2 = i1; i2 < VarNr; i2++) {
            Djv_isElm = 0;
            for (j = 0; j < EqnNr; j++)
              if ( Stoich[i][j] != 0 )
                for (k = 0; k < nElm; k++)
                  if ( (coeff_j[k]==j) && (coeff_i1[k]==i1)
                                      && (coeff_i2[k]==i2) ) {
                    Djv_isElm = 1;
                  }
            if (Djv_isElm == 1) Hess_NZ++ ;
    }  /* for i, i1, i2 */
  if (useDeclareValues)
       varTable[ NHESS ] -> value = max( Hess_NZ, 1 );

/* Allocate temporary index arrays */
  iHess_i = AllocIntegerVector(Hess_NZ, "iHess_i in GenerateHess");
  iHess_j = AllocIntegerVector(Hess_NZ, "iHess_j in GenerateHess");
  iHess_k = AllocIntegerVector(Hess_NZ, "iHess_k in GenerateHess");

  F_Hess  = DefFnc( "Hessian", 4, "function for Hessian (Jac derivative w.r.t. variables)");
  FunctionBegin( F_Hess, V, F, RCT, HESS );

  WriteComment("--------------------------------------------------------");
  WriteComment("Note: HESS is represented in coordinate sparse format: ");
  WriteComment("      HESS(m) = d^2 f_i / dv_j dv_k = d Jac_{i,j} / dv_k");
  WriteComment("      where i = IHESS_I(m), j = IHESS_J(m), k = IHESS_K(m).");
  WriteComment("--------------------------------------------------------");
  WriteComment("Note: d^2 f_i / dv_j dv_k = d^2 f_i / dv_k dv_j, ");
  WriteComment("      therefore only the terms d^2 f_i / dv_j dv_k");
  WriteComment("      with j <= k are computed and stored in HESS.");
  WriteComment("--------------------------------------------------------");

  if ( (useLang==C_LANG)||(useLang==F77_LANG)||(useLang==F90_LANG) ) {
    NewLines(1);
    WriteComment("Local variables");
    varTable[ NTMPD2A ] -> value = max( nElm, 1 );
    Declare( D2A );
  }

  NewLines(1);
  WriteComment("Computation of the second derivatives of equation rates");

/*  Generate d^2 A(j)/ ( d v(i1) d v(i2) )*/
  nElm = 0;
  for(j=0; j<EqnNr; j++)
    for (i1 = 0; i1 < VarNr; i1++)
      for (i2 = i1; i2 < VarNr; i2++) {

        if (i1==i2) {

         if (Stoich_Left[i1][j]>=2) {
            prod = RConst( j );
            for (i = 0; i < i1; i++)
              for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
                prod = Mul( prod, Elm( V, i ) );
            prod = Mul( prod, Const( Stoich_Left[i1][j] ) );
            prod = Mul( prod, Const( Stoich_Left[i1][j]-1 ) );
            for (k = 1; k <= (int)Stoich_Left[i1][j]-2; k++ )
              prod = Mul( prod, Elm( V, i1 ) );
            for (i = i1+1; i < VarNr; i++)
              for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
                prod = Mul( prod, Elm( V, i ) );
            for ( ; i < SpcNr; i++)
              for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
                prod = Mul( prod, Elm( F, i - VarNr ) );
            /* Comment the D2A */
            WriteComment("D2A(%d) = d^2 A(%d)/{dV(%d)dV(%d)}",Index(nElm),Index(j),Index(i1),Index(i2));
            Assign( Elm( D2A, nElm ), prod );
            nElm++;
          } /* if (Stoich_Left[i1][j]>=2) */

         } else {  /* i1 != i2 */
            if ( (Stoich_Left[i1][j]>=1)&&(Stoich_Left[i2][j]>=1) ) {
               prod = RConst( j );
               for (i = 0; i < i1; i++)
                 for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
                   prod = Mul( prod, Elm( V, i ) );
               prod = Mul( prod, Const( Stoich_Left[i1][j] ) );
               for (k = 1; k <= (int)Stoich_Left[i1][j]-1; k++ )
                 prod = Mul( prod, Elm( V, i1 ) );
               for (i = i1+1; i < i2; i++)
                 for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
                   prod = Mul( prod, Elm( V, i ) );
               prod = Mul( prod, Const( Stoich_Left[i2][j] ) );
               for (k = 1; k <= (int)Stoich_Left[i2][j]-1; k++ )
                 prod = Mul( prod, Elm( V, i2 ) );
               for (i = i2+1; i < VarNr; i++)
                 for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
                   prod = Mul( prod, Elm( V, i ) );
               for ( ; i < SpcNr; i++)
                 for (k = 1; k <= (int)Stoich_Left[i][j]; k++ )
                   prod = Mul( prod, Elm( F, i - VarNr ) );
            /* Comment the D2A */
               WriteComment("D2A(%d) = d^2 A(%d) / dV(%d)dV(%d)",
                                 Index(nElm),Index(j),Index(i1),Index(i2));
            Assign( Elm( D2A, nElm ), prod );
            nElm++;
            } /* if ( (Stoich_Left[i1][j]>=1)&&(Stoich_Left[i2][j]>=1) )  */
         }  /* if i1==i2 */

  } /* for j, i1, i2 */

    NewLines(1);
    WriteComment("Computation of the Jacobian derivative");

/*  Generate d^2 f(i)/ ( d v(i1) d v(i2) ) */
    Hess_NZ = 0;
    for (i = 0; i < VarNr; i++)
      for (i1 = 0; i1 < VarNr; i1++)
         for (i2 = i1; i2 < VarNr; i2++) {
            sum = Const(0);
            Djv_isElm = 0;
            for (j = 0; j < EqnNr; j++)
              if ( Stoich[i][j] != 0 )
                for (k = 0; k < nElm; k++)
                  if ( (coeff_j[k]==j) && (coeff_i1[k]==i1)
                                      && (coeff_i2[k]==i2) ) {
                    sum = Add( sum,
                            Mul( Const( Stoich[i][j] ), Elm( D2A, k ) ) );
                    Djv_isElm = 1;
                  }
            if (Djv_isElm == 1) {
               WriteComment("HESS(%d) = d^2 Vdot(%d)/{dV(%d)dV(%d)} = d^2 Vdot(%d)/{dV(%d)dV(%d)}",
                       Index(Hess_NZ),Index(i),Index(i1),Index(i2),Index(i),Index(i2),Index(i1));
               Assign( Elm( HESS, Hess_NZ ), sum );
               iHess_i[ Hess_NZ ] = i;
               iHess_j[ Hess_NZ ] = i1;
               iHess_k[ Hess_NZ ] = i2;
               Hess_NZ++;
            } else {
	      // Avoid memory leaks (Killian Murphy, 06 Jul 2022)
	      free(sum->elm);
	      free(sum);
            }

    }  /* for i, i1, i2 */


/* free temporary index arrays */
  free(coeff_j);  free(coeff_i1);  free(coeff_i2);

  MATLAB_Inline("\n   HESS = HESS(:);");

  FunctionEnd( F_Hess );

  FreeVariable( F_Hess );


  F_HessTR_VEC  = DefFnc( "HessTR_Vec", 4, "Hessian transposed times user vectors");
  FunctionBegin( F_HessTR_VEC, HESS, U1, U2, HTU );
  WriteComment("Compute the vector HTU =(Hess x U2)^T * U1 = d (Jac^T*U1)/d Var * U2 ");

  for (i=0; i<VarNr; i++) {
      sum = Const(0);
        for (k=0; k<Hess_NZ; k++) {
          if (iHess_j[k]==i)
             sum = Add( sum,
                   Mul( Elm( HESS, k ),
                   Mul( Elm( U1, iHess_i[k] ), Elm( U2, iHess_k[k] ) ) ) );
          else if (iHess_k[k]==i)
             sum = Add( sum,
                   Mul( Elm( HESS, k ),
                   Mul( Elm( U1, iHess_i[k] ), Elm( U2, iHess_j[k] ) ) ) );
        }
      Assign( Elm( HTU, i ), sum );
  }

  MATLAB_Inline("\n   HTU = HTU(:);");

  FunctionEnd(  F_HessTR_VEC );
  FreeVariable( F_HessTR_VEC );


  F_Hess_VEC  = DefFnc( "Hess_Vec", 4, "Hessian times user vectors");
  FunctionBegin( F_HessTR_VEC, HESS, U1, U2, HU );
  WriteComment("Compute the vector HU =(Hess x U2) * U1 = d (Jac*U1)/d Var * U2 ");

  for (i=0; i<VarNr; i++) {
      sum = Const(0);
      for (m=0; m<Hess_NZ; m++) {
          if (iHess_i[m]==i) {
             j = iHess_j[m];
             k = iHess_k[m];
             if (j==k) {
               sum = Add( sum,
                     Mul( Elm( HESS, m ),
                     Mul( Elm( U1, k ), Elm( U2, k ) ) ) );
             }
             else {
               /* This is the optimized code. It can lead to problems when
                          splitting for continuation lines. Therefore we
                          use the non-optimized code below the comment.
               sum = Add( sum,
                     Mul( Elm( HESS, m ),
                     Add( Mul( Elm( U1, j ), Elm( U2, k ) ),
                          Mul( Elm( U1, k ), Elm( U2, j ) ) ) ) ); */
               sum = Add( sum,
                     Mul( Elm( HESS, m ),
                     Mul( Elm( U1, j ), Elm( U2, k ) )  ) );
               sum = Add( sum,
                     Mul( Elm( HESS, m ),
                          Mul( Elm( U1, k ), Elm( U2, j ) ) ) );
             }
          }
      }
      Assign( Elm( HU, i ), sum );
  }

  MATLAB_Inline("\n   HU = HU(:);");

  FunctionEnd(  F_Hess_VEC );
  FreeVariable( F_Hess_VEC );
}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateHessianSparseData()
{
int k;


  UseFile( sparse_hessFile );

  NewLines(1);


  WriteComment("Hessian Sparse Data");
  WriteComment("");

  F77_Inline("%6sBLOCK DATA HESSIAN_SPARSE_DATA\n", " " );
  F77_Inline("%6sINCLUDE '%s_Sparse.h'", " ", rootFileName);
  F77_Inline("%6sINTEGER i", " ");

  if( (useLang==F77_LANG)||(useLang==F90_LANG)||(useLang==MATLAB_LANG) ) {
        for (k=0; k<Hess_NZ; k++) {
           iHess_i[k]++; iHess_j[k]++; iHess_k[k]++;
        }
  }

  InitDeclare( IHESS_I, Hess_NZ, (void*)iHess_i );
  InitDeclare( IHESS_J, Hess_NZ, (void*)iHess_j );
  InitDeclare( IHESS_K, Hess_NZ, (void*)iHess_k );

  if( (useLang==F77_LANG)||(useLang==F90_LANG) ) {
        for (k=0; k<Hess_NZ; k++) {
           iHess_i[k]--; iHess_j[k]--; iHess_k[k]--;
        }
  }

  NewLines(1);
  F77_Inline( "%6sEND\n\n", " " );

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateHessianSparseHeader()
{
  UseFile( sparse_dataFile );

  CommonName = "HESSDATA";

  NewLines(1);
  WriteComment(" ----------> Sparse Hessian Data");
  NewLines(1);

  ExternDeclare( IHESS_I );
  ExternDeclare( IHESS_J );
  ExternDeclare( IHESS_K );

  NewLines(1);
}





/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateStoicmSparseData()
{
int i,j,k, nnz_stoicm;
int *irow_stoicm;
int *ccol_stoicm;
int *icol_stoicm;
double *stoicm;

/* Compute the sparsity structure and allocate data structure vectors */
  nnz_stoicm = 0;
  for (j=0; j<EqnNr; j++)
    for (i=0; i<VarNr; i++)
      if ( Stoich[i][j] != 0.0 )
         nnz_stoicm++;
  if ( (irow_stoicm=(int*)calloc(nnz_stoicm+2,sizeof(int)) ) == NULL )
     FatalError(-30,"GenerateStoicmSparseData: Cannot allocate irow_stoicm");
  if ( (ccol_stoicm=(int*)calloc(EqnNr+2,sizeof(int)) ) == NULL )
     FatalError(-30,"GenerateStoicmSparseData: Cannot allocate ccol_stoicm");
  if ( (icol_stoicm=(int*)calloc(nnz_stoicm+2,sizeof(int)) ) == NULL )
     FatalError(-30,"GenerateStoicmSparseData: Cannot allocate icol_stoicm");
  if ( (stoicm=(double*)calloc(nnz_stoicm+2,sizeof(double)) ) == NULL )
     FatalError(-30,"GenerateStoicmSparseData: Cannot allocate stoicm");

  UseFile( sparse_stoicmFile );

  nnz_stoicm = 0;
  for (j=0; j<EqnNr; j++) {
    ccol_stoicm[ j ] =  nnz_stoicm;
    for (i=0; i<VarNr; i++) {
      if ( Stoich[i][j] != 0 ) {
         irow_stoicm[ nnz_stoicm ] = i;
         icol_stoicm[ nnz_stoicm ] = j;
         stoicm[ nnz_stoicm ] = Stoich[i][j];
         nnz_stoicm++;
      }
    }
  }
  ccol_stoicm[ EqnNr ] =  nnz_stoicm;

  if( (useLang==F77_LANG)||(useLang==F90_LANG) ) {
        for (k=0; k<nnz_stoicm; k++) {
           irow_stoicm[k]++; icol_stoicm[k]++;
        }
        for (k=0; k<=EqnNr; k++) {
           ccol_stoicm[k]++;
        }
  }


  NewLines(1);
  WriteComment(" Stoichiometric Matrix in Compressed Column Sparse Format");
  NewLines(1);
  F77_Inline("%6sBLOCK DATA STOICM_MATRIX\n", " " );
  F77_Inline("%6sINCLUDE '%s_Sparse.h'", " ", rootFileName);
  F77_Inline("%6sINTEGER i", " ");
  InitDeclare( CCOL_STOICM, EqnNr+1, (void*)ccol_stoicm );
  InitDeclare( IROW_STOICM, nnz_stoicm, (void*)irow_stoicm );
  InitDeclare( ICOL_STOICM, nnz_stoicm, (void*)icol_stoicm );
  InitDeclare( STOICM, nnz_stoicm, (void*)stoicm );
  NewLines(1);
  F77_Inline( "%6sEND\n\n", " " );


  UseFile( param_headerFile );
  CommonName = "GDATA";
  NewLines(1);
  DeclareConstant( NSTOICM, ascii( max( nnz_stoicm, 1 ) ) );

/* Free data structure vectors */
  free(irow_stoicm); free(ccol_stoicm); free(icol_stoicm); free(stoicm);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateStoicmSparseHeader()
{
  UseFile( sparse_dataFile );

  NewLines(1);
  WriteComment(" ----------> Sparse Stoichiometric Matrix");
  NewLines(1);
  CommonName = "STOICM_VALUES";
  ExternDeclare( STOICM );
  CommonName = "STOICM_DATA";
  ExternDeclare( IROW_STOICM );
  ExternDeclare( CCOL_STOICM );
  ExternDeclare( ICOL_STOICM );
  NewLines(1);

  NewLines(1);
  WriteComment(" ----------> Sparse Data for Jacobian of Reactant Products");
  NewLines(1);
  CommonName = "JVRP";
  ExternDeclare( ICOL_JVRP );
  ExternDeclare( IROW_JVRP );
  ExternDeclare( CROW_JVRP );
  NewLines(1);

}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateJacVect()
{
int i, j, nElm;
int Jac_VEC;
int Jac_SP_VEC;

  if( useLang == MATLAB_LANG ) return;

  if( VarNr == 0 ) return;

  UseFile( jacobianFile );
  Jac_VEC = DefFnc( "Jac_Vec", 3,
       "function for sparse multiplication: square Jacobian times vector");
  Jac_SP_VEC = DefFnc( "Jac_SP_Vec", 3,
        "function for sparse multiplication: sparse Jacobian times vector");

  if ( useJacSparse ) {
    FunctionBegin( Jac_SP_VEC, JVS, UV, JUV );
    nElm = 0;
    for( i = 0; i < VarNr; i++) {
      sum = Const(0);
      for( j = 0; j < VarNr; j++ )
        if( LUstructJ[i][j] ) {
          if( structJ[i][j] != 0 )
            sum = Add( sum, Mul( Elm( JVS, nElm ), Elm( UV, j ) ) );
          nElm++;
          }
      Assign( Elm( JUV, i ), sum );
    }
    FunctionEnd( Jac_SP_VEC );
    }

  else {
    FunctionBegin( Jac_VEC, JV, UV, JUV );
    for( i = 0; i < VarNr; i++) {
      sum = Const(0);
      for( j = 0; j < VarNr; j++ )
        if( structJ[i][j] != 0 )
          sum = Add( sum, Mul( Elm( JV, i, j ), Elm( UV, j ) ) );
      Assign( Elm( JUV, i ), sum );
    }
    FunctionEnd( Jac_VEC );
  }

  FreeVariable( Jac_VEC );
  FreeVariable( Jac_SP_VEC );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateJacTRVect()
{
int i, j, nElm;
int JacTR_VEC;
int JacTR_SP_VEC;
int **TmpStruct;

  if( useLang == MATLAB_LANG ) return;

  if ( VarNr == 0 ) return;

  UseFile( jacobianFile );

  JacTR_VEC = DefFnc( "JacTR_Vec", 3,
       "sparse multiplication: square Jacobian transposed times vector");
  JacTR_SP_VEC = DefFnc( "JacTR_SP_Vec", 3,
        "sparse multiplication: sparse Jacobian transposed times vector");

  if ( useJacSparse ) {

  /* The temporary array of structure */
  TmpStruct = AllocIntegerMatrix( VarNr, VarNr, "TmpStruct in GenerateJacTRVect" );

    nElm = 0;
    for( i = 0; i < VarNr; i++)
      for( j = 0; j < VarNr; j++ )
        if( LUstructJ[i][j] ) {
          TmpStruct[i][j] = nElm;
          nElm++;
        }

    FunctionBegin( JacTR_SP_VEC, JVS, UV, JTUV );
    nElm = 0;
    for( i = 0; i < VarNr; i++) {
      sum = Const(0);
      for( j = 0; j < VarNr; j++ )
         if( structJ[j][i] != 0 )
           sum = Add( sum, Mul( Elm( JVS, TmpStruct[j][i] ), Elm( UV, j ) ) );
      Assign( Elm( JTUV, i ), sum );
    }
    FunctionEnd( JacTR_SP_VEC );

  /* Free the temporary array of structure */
  FreeIntegerMatrix( TmpStruct, VarNr, VarNr );

    } /* useJacSparse*/

  else {
    FunctionBegin( JacTR_VEC, JV, UV, JTUV );
    for( i = 0; i < VarNr; i++) {
      sum = Const(0);
      for( j = 0; j < VarNr; j++ )
        if( structJ[j][i] != 0 )
          sum = Add( sum, Mul( Elm( JV, j, i ), Elm( UV, j ) ) );
      Assign( Elm( JTUV, i ), sum );
    }
    FunctionEnd( JacTR_VEC );
  }

  FreeVariable( JacTR_VEC );
  FreeVariable( JacTR_SP_VEC );
}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateSparseUtil()
{
int SUTIL;

  if ( useLang == MATLAB_LANG ) return;

  UseFile( linalgFile );

  SUTIL = DefFnc( "SPARSE_UTIL", 0, "SPARSE utility functions");
  CommentFunctionBegin( SUTIL );

  IncludeCode( "%s/util/sutil", Home );

  CommentFunctionEnd( SUTIL );
  FreeVariable( SUTIL );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateBlas()
{
int BLAS;

  if ( useLang == MATLAB_LANG ) return;

  UseFile( linalgFile );

  BLAS = DefFnc( "BLAS_UTIL", 0, "BLAS-LIKE utility functions");
  CommentFunctionBegin( BLAS );

  IncludeCode( "%s/util/blas", Home );

  CommentFunctionEnd( BLAS );
  FreeVariable( BLAS );
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateDFunDRcoeff()
{

  UseFile( stoichiomFile );

  NewLines(1);
  WriteComment("Begin Derivative w.r.t. Rate Coefficients");
  NewLines(1);

  IncludeCode( "%s/util/dFun_dRcoeff", Home );

  NewLines(1);
  WriteComment("End Derivative w.r.t. Rate Coefficients");
  NewLines(1);

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateDJacDRcoeff()
{

  UseFile( stoichiomFile );

  NewLines(1);
  WriteComment("Begin Jacobian Derivative w.r.t. Rate Coefficients");
  NewLines(1);

  IncludeCode( "%s/util/dJac_dRcoeff", Home );

  NewLines(1);
  WriteComment("End Jacobian Derivative w.r.t. Rate Coefficients");
  NewLines(1);

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateSolve()
{
int i, j;
int SOLVE;
int *irow;
int *icol;
int *crow;
int *diag;
int ibgn, iend;
int useLangOld;
int dim;

  if( useLang == MATLAB_LANG ) return;

  /* Allocate local arrays for dimension dim */
  dim  = VarNr+2;
  irow = AllocIntegerVector( dim*dim, "irow in GenerateSolve" );
  icol = AllocIntegerVector( dim*dim, "icol in GenerateSolve" );
  crow = AllocIntegerVector( dim,     "crow in GenerateSolve" );
  diag = AllocIntegerVector( dim,     "diag in GenerateSolve" );

  useLangOld = useLang;
  useLang = C_LANG;
  NonZero( LU, 0, VarNr, irow, icol, crow, diag );
  useLang = useLangOld;

  UseFile( linalgFile );

  SOLVE = DefFnc( "KppSolve", 2, "sparse back substitution");
  FunctionBegin( SOLVE, JVS, X );

  for( i = 0; i < VarNr; i++) {
    ibgn = crow[i];
    iend = diag[i];
    if( ibgn <= iend ) {
      sum = Elm( X, i );
      if ( ibgn < iend ) {
              if( doAutoReduce ) F90_Inline("IF (DO_SLV(%d)) &",i+1);
        for( j = ibgn; j < iend; j++ ) {
          sum = Sub( sum, Mul( Elm( JVS, j ), Elm( X, icol[j] ) ) );
        }
        Assign( Elm( X, i ), sum );
      } else {
	// Avoid memory leaks (Killian Murphy, 07 Jul 2022)
	free(sum->elm);
	free(sum);
      }
    }
  }

  for( i = VarNr-1; i >=0; i--) {
    ibgn = diag[i] + 1;
    iend = crow[i+1];
    sum = Elm( X, i );
    if( doAutoReduce ) F90_Inline("IF (DO_SLV(%d)) &",i+1);
    for( j = ibgn; j < iend; j++ ) {
      sum = Sub( sum, Mul( Elm( JVS, j ), Elm( X, icol[j] ) ) );
    }
    sum = Div( sum, Elm( JVS, diag[i] ) );
    Assign( Elm( X, i ), sum );
  }

  FunctionEnd( SOLVE );
  FreeVariable( SOLVE );

  /* Free Local Arrays */
  free(irow);
  free(icol);
  free(crow);
  free(diag);
}




/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateTRSolve()
{
int i, j;
int SOLVETR;
int *irow;
int *icol;
int *crow;
int *diag;
int ibgn, iend;
int useLangOld;
int **pos;
int dim;

  if( useLang == MATLAB_LANG ) return;

  /* Allocate local arrays for dimension dim */
  dim  = VarNr+2;
  irow = AllocIntegerVector( dim*dim, "irow in GenerateTRSolve" );
  icol = AllocIntegerVector( dim*dim, "icol in GenerateTRSolve" );
  crow = AllocIntegerVector( dim,     "crow in GenerateTRSolve" );
  diag = AllocIntegerVector( dim,     "diag in GenerateTRSolve" );
  pos  = AllocIntegerMatrix( dim+1, dim+1, "pos in GenerateTRSolve");

  useLangOld = useLang;
  useLang = C_LANG;
  NonZero( LU, 0, VarNr, irow, icol, crow, diag );
  useLang = useLangOld;

  UseFile( linalgFile );

  SOLVETR = DefFnc( "KppSolveTR", 3, "sparse, transposed back substitution");
  FunctionBegin( SOLVETR, JVS, X, XX );
  for( i = 0; i < VarNr; i++) {
   for( j = 0; j < VarNr; j++)
    pos[i][j]=-1;
  }
  for( i = 0; i < VarNr; i++) {
    ibgn = crow[i];
    iend = diag[i];
    if( ibgn <= iend ) {
      if ( ibgn < iend ) {
        for( j = ibgn; j < iend; j++ )
          pos[icol[j]][i]=j;
      }
    }
  }

  for( i = VarNr-1; i >=0; i--) {
    ibgn = diag[i] + 1;
    iend = crow[i+1];
    for( j = ibgn; j < iend; j++ )
      pos[icol[j]][i]=j;
    pos[i][i]=diag[i];
  }

  for( i = 0; i<VarNr; i++) {
    sum = Elm( X, i );
    for (j=0; j<i; j++ ){
     if(pos[i][j] >= 0) {
      sum=Sub( sum, Mul ( Elm(JVS,pos[i][j] ), Elm( XX, j ) ) );
     }
    }
    sum=Div( sum, Elm(JVS, diag[i] ) );
    Assign( Elm( XX, i ), sum );
  }
  for( i = VarNr-1; i >=0; i--) {
    sum = Elm( XX, i );
    for (j=i+1; j<VarNr; j++) {
     if(pos[i][j] >= 0) {
      sum=Sub( sum, Mul ( Elm(JVS,pos[i][j] ), Elm( XX, j ) ) );
     }
    }
    Assign( Elm( XX, i ), sum );
   }

  FunctionEnd( SOLVETR );
  FreeVariable( SOLVETR );
  /* Free Local Arrays */
  free(irow); free(icol); free(crow); free(diag);
  FreeIntegerMatrix(pos, dim+1, dim+1);
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateRateLaws()
{

  UseFile( rateFile );

  NewLines(1);
  WriteComment("Begin Rate Law Functions from KPP_HOME/util/UserRateLaws");
  NewLines(1);
  IncludeCode( "%s/util/UserRateLaws", Home );
  NewLines(1);
  WriteComment("End Rate Law Functions from KPP_HOME/util/UserRateLaws");
  NewLines(1);

  NewLines(1);
  WriteComment("Begin INLINED Rate Law Functions");
  NewLines(1);

  switch( useLang ) {
    case C_LANG:  bprintf( InlineCode[ C_RATES ].code );
                 break;
    case F77_LANG: bprintf( InlineCode[ F77_RATES ].code );
                 break;
    case F90_LANG: bprintf( InlineCode[ F90_RATES ].code );
                 break;
    case MATLAB_LANG: bprintf( InlineCode[ MATLAB_RATES ].code );
                 break;
  }
  FlushBuf();

  NewLines(1);
  WriteComment("End INLINED Rate Law Functions");
  NewLines(1);


}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateUpdateSun()
{
int UPDATE_SUN;

  UseFile( rateFile );

  UPDATE_SUN = DefFnc( "Update_SUN", 0, "update SUN light using TIME");
  CommentFunctionBegin( UPDATE_SUN );

  IncludeCode( "%s/util/UpdateSun", Home );

  CommentFunctionEnd( UPDATE_SUN );
  FreeVariable( UPDATE_SUN );
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateUpdateRconst()
{
int i;
int UPDATE_RCONST;

  UseFile( rateFile );

  UPDATE_RCONST = DefFnc( "Update_RCONST", 0, "function to update rate constants");

  FunctionBegin( UPDATE_RCONST );
  F77_Inline("      INCLUDE '%s_Global.h'", rootFileName);
  MATLAB_Inline("global SUN TEMP RCONST");

  if ( useLang==F77_LANG )
      IncludeCode( "%s/util/UserRateLaws_FcnHeader", Home );

  NewLines(1);

  NewLines(1);
  WriteComment("Begin INLINED RCONST");
  NewLines(1);

  switch( useLang ) {
    case C_LANG:  bprintf( InlineCode[ C_RCONST ].code );
                 break;
    case F77_LANG: bprintf( InlineCode[ F77_RCONST ].code );
                 break;
    case F90_LANG: bprintf( InlineCode[ F90_RCONST ].code );
                 break;
    case MATLAB_LANG: bprintf( InlineCode[ MATLAB_RCONST ].code );
                 break;
  }
  FlushBuf();

  NewLines(1);
  WriteComment("End INLINED RCONST");
  NewLines(1);

  for( i = 0; i < EqnNr; i++) {
    /*  mz_rs_20220701+ */
    if (useEqntags==1) {
      F90_Inline("! <%s>:", kr[i].label);
    }
    /*  mz_rs_20220701- */
    if( kr[i].type == EXPRESION )
      Assign( Elm( RCONST, i ), Elm( KR, kr[i].val.st ) );
    if( kr[i].type == PHOTO )
      Assign( Elm( RCONST, i ), Elm( KR, kr[i].val.st ) );
    /* mz_rs_20050117+ */
    if ( kr[i].type == NUMBER ) {
      F90_Inline("! RCONST(%d) = constant rate coefficient", i+1);
    }
    /* mz_rs_20050117- */
  }

  MATLAB_Inline("   RCONST = RCONST(:);");

  FunctionEnd( UPDATE_RCONST );
  FreeVariable( UPDATE_RCONST );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateUpdatePhoto()
{
int i;
int UPDATE_PHOTO;

  UseFile( rateFile );

  UPDATE_PHOTO = DefFnc( "Update_PHOTO", 0, "function to update photolytical rate constants");

  FunctionBegin( UPDATE_PHOTO );
  F77_Inline("      INCLUDE '%s_Global.h'", rootFileName);
  /*  mz_rs_20220212+ */
  /* global is already used in the rates module, don't use it twice */
  /* F90_Inline("   USE %s_Global", rootFileName); */
  /*  mz_rs_20220212- */
  MATLAB_Inline("global SUN TEMP RCONST");

  NewLines(1);
  WriteComment("Begin INLINED RCONST");
  NewLines(1);

  switch( useLang ) {
    case C_LANG:  bprintf( InlineCode[ C_RCONST ].code );
                 break;
    case F77_LANG: bprintf( InlineCode[ F77_RCONST ].code );
                 break;
    case F90_LANG: bprintf( InlineCode[ F90_RCONST ].code );
                 break;
    case MATLAB_LANG: bprintf( InlineCode[ MATLAB_RCONST ].code );
                 break;
  }
  FlushBuf();

  NewLines(1);
  WriteComment("End INLINED RCONST");
  NewLines(1);

  for( i = 0; i < EqnNr; i++) {
    if( kr[i].type == PHOTO )
      Assign( Elm( RCONST, i ), Elm( KR, kr[i].val.st ) );
  }

  FunctionEnd( UPDATE_PHOTO );
  FreeVariable( UPDATE_PHOTO );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateIntegrator()
{
int TIN, TOUT, INTEGRATE;

  UseFile( integratorFile );

  TIN = DefElm( "TIN", real, "Start Time for Integration");
  TOUT = DefElm( "TOUT", real, "End Time for Integration");
  INTEGRATE = DefFnc( "INTEGRATE", 2, "Integrator routine");
  CommentFunctionBegin( INTEGRATE, TIN, TOUT );

  if( strchr( integrator, '/' ) )
    IncludeCode( integrator );
  else
    IncludeCode( "%s/int/%s", Home, integrator );

  CommentFunctionEnd( INTEGRATE );
  FreeVariable( INTEGRATE );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateDriver()
{
int MAIN;

  UseFile( driverFile );

  MAIN = DefFnc( "MAIN", 0, "Main program - driver routine");
  CommentFunctionBegin( MAIN );

  if( strchr( driver, '/' ) )
    IncludeCode( driver );
  else
    IncludeCode( "%s/drv/%s", Home, driver );

  CommentFunctionEnd( MAIN );
  FreeVariable( MAIN );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateUtil()
{
int UTIL;

  UseFile( utilFile );
  NewLines(1);
  WriteComment("User INLINED Utility Functions");

  switch( useLang ) {
    case C_LANG:  bprintf( InlineCode[ C_UTIL ].code );
                 break;
    case F77_LANG:  bprintf( InlineCode[ F77_UTIL ].code );
                 break;
    case F90_LANG:  bprintf( InlineCode[ F90_UTIL ].code );
                 break;
    case MATLAB_LANG:bprintf( InlineCode[ MATLAB_UTIL ].code );
                 break;
  }
  FlushBuf();

  NewLines(1);
  WriteComment("End INLINED Utility Functions");
  NewLines(1);

  WriteComment("Utility Functions from KPP_HOME/util/util");
  UTIL = DefFnc( "UTIL", 0, "Utility functions");
  CommentFunctionBegin( UTIL);

  IncludeCode( "%s/util/util", Home );

  if ((useLang == F90_LANG) && (useEqntags==1)) {
    IncludeCode( "%s/util/tag2num", Home );
  }

  WriteComment("End Utility Functions from KPP_HOME/util/util");
  CommentFunctionEnd( UTIL );
  FreeVariable( UTIL );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateParamHeader()
{
int spc;
int i;
char name[MAX_SPNAME];

int j,dummy_species;

/*  ---------->  First declaration of constants */
  UseFile( param_headerFile );

  NewLines(1);
  DeclareConstant( NSPEC,   ascii( max(SpcNr, 1) ) );
  DeclareConstant( NVAR,    ascii( max(VarNr, 1) ) );
  DeclareConstant( NFAM,    ascii( max(FamilyNr,1)  ) );
  DeclareConstant( NVARACT, ascii( max(VarActiveNr, 1) ) );
  DeclareConstant( NFIX,    ascii( max(FixNr, 1) ) );
  DeclareConstant( NREACT,  ascii( max(EqnNr, 1) ) );
  DeclareConstant( NVARST,  ascii( VarStartNr ) );
  DeclareConstant( NFIXST,  ascii( FixStartNr ) );
  DeclareConstant( NONZERO, ascii( max(Jac_NZ, 1) ) );
  DeclareConstant( LU_NONZERO, ascii( max(LU_Jac_NZ, 1) ) );
  DeclareConstant( CNVAR,   ascii( VarNr+1 ) );
  if ( useStoicmat ) {
        DeclareConstant( CNEQN,   ascii( EqnNr+1 ) );
  }
  if ( useHessian ) {
        DeclareConstant( NHESS,   ascii( max(Hess_NZ, 1) ) );
  }

  DeclareConstant( NLOOKAT,  ascii( nlookat ) );
  DeclareConstant( NMONITOR,  ascii( nmoni ) );
  DeclareConstant( NMASS, ascii( nmass ) );

  NewLines(1);
  WriteComment("Index declaration for variable species in C and VAR");
  WriteComment("  VAR(ind_spc) = C(ind_spc)");
  NewLines(1);
  for( i = 0; i < VarNr; i++) {
    sprintf( name, "%s", "ind_" );
    strncat( name, SpeciesTable[ Code[i] ].name,
             strlen(SpeciesTable[ Code[i] ].name)+1 );
    spc = DefConst( name, INT, 0 );
    DeclareConstant( spc, ascii( Index(i) ) );
    FreeVariable( spc );
  }

  NewLines(1);
  WriteComment("Index declaration for fixed species in C");
  WriteComment("  C(ind_spc)");
  NewLines(1);
  for( i = 0; i < FixNr; i++) {
    sprintf( name, "%s", "ind_" );
    strncat( name, SpeciesTable[ Code[i + VarNr] ].name,
             strlen(SpeciesTable[ Code[i + VarNr] ].name)+1 );
    spc = DefConst( name, INT, 0 );
    DeclareConstant( spc, ascii( Index(i+VarNr) ) );
    FreeVariable( spc );
  }

  if (useDummyindex==1) {
    NewLines(1);
    WriteComment("Index declaration for dummy species");
    NewLines(1);
    for( i = 0; i < MAX_SPECIES; i++) {
      if (SpeciesTable[i].type == 0) continue;
      dummy_species = 1;
      for( j = 0; j < MAX_SPECIES; j++)
        if (Code[j] == i) dummy_species = 0;
      if (dummy_species) {
        sprintf( name, "%s", "ind_" );
        strncat( name, SpeciesTable[i].name, strlen(SpeciesTable[i].name)+1 );
        spc = DefConst( name, INT, 0 );
        DeclareConstant( spc, ascii( 0 ) );
        FreeVariable( spc );
      }
    }
  }

  NewLines(1);
  WriteComment("Index declaration for fixed species in FIX");
  WriteComment("   FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)");
  NewLines(1);
  for( i = 0; i < FixNr; i++) {
    sprintf( name, "%s", "indf_" );
    strncat( name, SpeciesTable[ Code[i + VarNr] ].name,
             strlen( SpeciesTable[ Code[i + VarNr] ].name)+1 );
    spc = DefConst( name, INT, 0 );
    DeclareConstant( spc, ascii( Index(i) ) );
    FreeVariable( spc );
  }
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateGlobalHeader()
{
  int useFortran;

//===========================================================================
// MODIFICATION: Bob Yantosca (26 Mar 2021)
//
// Modify code to inline the F77/F90 THREADPRIVATE declarations, and also
// declare extra arrays and scalars.  These declarations begin with a
// comment character and will be ignored unless the Fortran-90 code is
// compiled with OpenMP.
//===========================================================================

  /*** Define a flag to denote if we are using F90 or F77 ***/
  if ( useLang == F90_LANG || useLang == F77_LANG ) { useFortran = 1; }
  else                                              { useFortran = 0; }

  UseFile( global_dataFile );

  CommonName = "GDATA";

  /*** Write comment header ***/
  NewLines(1);
  WriteComment("Declaration of global variables");
  if ( useFortran ) {
    NewLines(1);
    WriteComment(
      "~~~ If you are using KPP within an OpenMP parallel environment,");
    WriteComment(
      "~~~ then these variables must be declared THREADPRIVATE.  This means");
    WriteComment(
      "~~~ that the compiler will make a private copy of these variables");
    WriteComment(
      "~~~ (in stack memory) for each execution thread.  At the end of ");
    WriteComment(
      "~~~ the OpenMP parallel loop, these variables will be finalized,");
    WriteComment(
      "~~~ and their memory deallocated.");
    WriteComment("~~~");
    WriteComment(
      "~~~ NOTE: Because the OpenMP commands all begin with a comment");
    WriteComment(
      "~~~ character, they will be ignored unless the code is compiled");
    WriteComment(
      "~~~ with OpenMP parallelization turned on.");
  }
  NewLines(1);

  /*** Declare C, concentration array ***/
  ExternDeclare( C );
  if ( useFortran ) { WriteOMPThreadPrivate("C"); }

  /*** Declare VAR and FIX for F90 (these are now pointers!) ***/
  if ( useLang == F90_LANG ) {
    ExternDeclare( VAR );
    WriteOMPThreadPrivate("VAR");

    ExternDeclare( FIX );
    WriteOMPThreadPrivate("FIX");
  }

  /*** Declare VAR and fix for F77                             ***/
  /*** We need to keep the EQUIVALENCE statement for F77 only! ***/
  if ( useLang == F77_LANG ) {
    Declare( VAR );
    WriteOMPThreadPrivate("VAR");
    WriteComment("VAR, FIX are chunks of array C");
    F77_Inline("!      EQUIVALENCE( %s(%d),%s(1) )",
               varTable[C]->name, 1, varTable[VAR]->name );

    Declare( FIX );
    WriteOMPThreadPrivate("FIX");
    if ( FixNr > 0 ) {
      F77_Inline("!      EQUIVALENCE( %s(%d),%s(1) )",
                 varTable[C]->name, VarNr+1, varTable[FIX]->name );
    }
  }

  /*** Declare VAR and FIX for MatLab ***/
  if ( useLang == MATLAB_LANG ) {
    ExternDeclare( VAR );
    ExternDeclare( FIX );
  }

  /*** Declare VAR and FIX for C with the extern property ***/
  if ( useLang == C_LANG ) {
    C_Inline("  extern %s * %s;", C_types[real], varTable[VAR]->name );
    C_Inline("  extern %s * %s;", C_types[real], varTable[FIX]->name );
  }

  /*** Declare all other threadprivate variables ***/
  ExternDeclare( RCONST );
  if ( useFortran ) { WriteOMPThreadPrivate("RCONST"); }

  ExternDeclare( TIME );
  if ( useFortran ) { WriteOMPThreadPrivate("TIME"); }

  ExternDeclare( SUN );
  if ( useFortran ) { WriteOMPThreadPrivate("SUN"); }

  ExternDeclare( TEMP );
  if ( useFortran ) { WriteOMPThreadPrivate("TEMP"); }

  /*** Declare non-threadprivate variables ***/
  NewLines(1);
  if ( useFortran ) {
    WriteComment(
     "~~~ If you are using KPP within an OpenMP parallel environment,");
    WriteComment(
     "~~~ these variables DO NOT need to be declared THREADPRIVATE.");
    NewLines(1);
  }

  ExternDeclare( TSTART );
  ExternDeclare( TEND );
  ExternDeclare( DT );
  ExternDeclare( ATOL );
  ExternDeclare( RTOL );
  ExternDeclare( STEPMIN );
  ExternDeclare( STEPMAX );
  if (doAutoReduce) {
    ExternDeclare(DO_JVS);
    ExternDeclare(DO_SLV);
    ExternDeclare(DO_FUN);
    ExternDeclare(cLU_IROW);
    ExternDeclare(cLU_ICOL);
    ExternDeclare(cLU_CROW);
    ExternDeclare(cLU_DIAG);
    ExternDeclare(JVS_MAP);
    ExternDeclare(SPC_MAP);
    ExternDeclare(iSPC_MAP);
    ExternDeclare(RMV);
    ExternDeclare(RNVAR);
    ExternDeclare(cNONZERO);
    ExternDeclare(keepSpcActive);
    ExternDeclare(keepActive);
    /* hplin 10/19/21 */
    WriteOMPThreadPrivate(" DO_JVS, DO_SLV, DO_FUN, cLU_IROW, cLU_ICOL, cLU_CROW");
    WriteOMPThreadPrivate(" cLU_DIAG, JVS_MAP, SPC_MAP, iSPC_MAP, RMV, rNVAR, cNONZERO, KEEPACTIVE "); /* msl_20160419 */
  }
  ExternDeclare( CFACTOR );
  if (useStochastic)
      ExternDeclare( VOLUME );

  CommonName = "INTGDATA";
  if ( useHessian ) {
     ExternDeclare( DDMTYPE );
  }


  if ( (useLang == C_LANG) || (useLang == F77_LANG) ) {
     CommonName = "INTGDATA";
     ExternDeclare( LOOKAT );
     ExternDeclare( MONITOR );
     CommonName = "CHARGDATA";
     ExternDeclare( SPC_NAMES );
     ExternDeclare( SMASS );
     ExternDeclare( EQN_NAMES );
     ExternDeclare( EQN_TAGS );
     ExternDeclare( FAM_NAMES );
  }

  NewLines(1);
  WriteComment("INLINED global variable declarations");

  switch( useLang ) {
    case C_LANG:  bprintf( InlineCode[ C_GLOBAL ].code );
                 break;
    case F77_LANG: bprintf( InlineCode[ F77_GLOBAL ].code );
                 break;
    case F90_LANG: bprintf( InlineCode[ F90_GLOBAL ].code );
                 break;
    case MATLAB_LANG: bprintf( InlineCode[ MATLAB_GLOBAL ].code );
                 break;
  }
  FlushBuf();

  NewLines(1);
  WriteComment("INLINED global variable declarations");
  NewLines(1);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void WriteSpec( int i, int j )
{
char buf[100];

  if( Reactive[j] )
    sprintf( buf, "%s (r)", SpeciesTable[ Code[j] ].name );
  else
    sprintf( buf, "%s (n)", SpeciesTable[ Code[j] ].name );
  WriteAll("%4d = %-23s", 1 + i, buf );
  FlushBuf();
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int EqnStr( int eq, char * buf, float** mat )
{
int spc, first;

/* bugfix if stoichiometric factor is not an integer */
char s[40];
char tmpStr[MAX_EQNLEN];

  first = 1;
  *buf = 0;
  for( spc = 0; spc < SpcNr; spc++ )
    if( mat[spc][eq] != 0 ) {
      if( ((mat[spc][eq] == 1)||(mat[spc][eq] == -1)) ) {
        s[0] = '\0'; // necessary to remove old contents
      } else {
        /* real */
        /*  mz_rs_20050130+ */
        /* sprintf(s, "%g", mat[spc][eq]); */
        /* remove the minus sign with fabs(), it will be re-inserted later */
        sprintf(s, "%g", mat[spc][eq]?mat[spc][eq]:(-mat[spc][eq]));
        /*  mz_rs_20050130- */
        /* remove trailing zeroes */
        /* -- Commented out my MSL - Oct 26, 2016: Unable to manage integers ending in '0'*/
        /*for (n= strlen(s) - 1; n >= 0; n--)
          if (s[n] != '0') break; */
        s[strlen(s)]= '\0';
        strcat( s, " " );
      }

      if( first ) {
        if( mat[spc][eq] > 0 ) {
          sprintf( tmpStr, "%s", s );
          strncat( buf, tmpStr, strlen(tmpStr)+1 );
        } else {
          sprintf( tmpStr, "%s",   buf              );
          strncat( tmpStr, "- ",   3                );
          strncat( tmpStr, s,      strlen(s)+1      );
          strncpy( buf,    tmpStr, strlen(tmpStr)+1 );
        }
        first = 0;
      } else {
        if( mat[spc][eq] > 0 ) {
          sprintf( tmpStr, "%s",   buf              );
          strncat( tmpStr, " + ",  4                );
          strncat( tmpStr, s,      strlen(s)+1      );
          strncpy( buf,    tmpStr, strlen(tmpStr)+1 );
        } else {
          sprintf( tmpStr, "%s",   buf              );
          strncat( tmpStr, " - ",  4                );
          strncat( tmpStr, s,      strlen(s)+1      );
          strncpy( buf,    tmpStr, strlen(tmpStr)+1 );
        }
      }
      sprintf( tmpStr, "%s", buf );
      strncat( tmpStr, SpeciesTable[ Code[spc] ].name,
               strlen(  SpeciesTable[ Code[spc] ].name )+1 );
      strncpy( buf, tmpStr, strlen(tmpStr)+1 );
    }

  return strlen(buf);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int EqnString( int eq, char * buf )
{
static int lhs = 0;
static int rhs = 0;
int len_rhsbuf;

int i, l;
char lhsbuf[MAX_EQNLEN], rhsbuf[MAX_EQNLEN];

  if(lhs == 0) for( i = 0; i < EqnNr; i++ ) {
                 l = EqnStr( i, lhsbuf, Stoich_Left);
                 lhs = (lhs > l) ? lhs : l;
               }

  if(rhs == 0) for( i = 0; i < EqnNr; i++ ) {
                 l = EqnStr( i, rhsbuf, Stoich_Right);
                 rhs = (rhs > l) ? rhs : l;
               }

  EqnStr( eq, lhsbuf, Stoich_Left);
  len_rhsbuf=EqnStr( eq, rhsbuf, Stoich_Right);

  if(lhs+5+len_rhsbuf>100) // 100 = len of EQN_NAMES string; 5 = len of " --> "
    {
      // truncate the list of products in kpp_monitor file:
      rhsbuf[100-5-lhs-8]=0; // 8 = len of "... etc."
      sprintf(buf, "%*s --> ", lhs, lhsbuf);
      strcat(buf, rhsbuf);
      strcat(buf, "... etc.");
    }
  else
    {
      sprintf(buf, "%*s --> %-*s", lhs, lhsbuf, rhs, rhsbuf);
    }

  return strlen(buf);
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateLog()
{
int i, dn;

  UseFile( logFile );

  WriteAll("### Options -------------------------------------------\n");
  NewLines(1);
  if( useDeclareValues ) { WriteAll("#DECLARE      - ON\n");                 }
     else                { WriteAll("#DECLARE      - OFF\n");                }
  if( useDouble )        { WriteAll("#DOUBLE       - ON\n");                 }
     else                { WriteAll("#DOUBLE       - OFF\n");                }
                           WriteAll("#DRIVER       - %s\n", driver);
  if( useDummyindex )    { WriteAll("#DUMMYINDEX   - ON\n");                 }
     else                { WriteAll("#DUMMYINDEX   - OFF\n");                }
  if( useEqntags )       { WriteAll("#EQNTAGS      - ON\n");                 }
     else                { WriteAll("#EQNTAGS      - OFF\n");                }
  if( useAggregate )     { WriteAll("#FUNCTION     - AGGREGATE\n");          }
    else                 { WriteAll("#FUNCTION     - SPLIT\n");              }
  if( useHessian )       { WriteAll("#HESSIAN      - ON\n");                 }
     else                { WriteAll("#HESSIAN      - OFF\n");                }
                           WriteAll("#INTEGRATOR   - %s\n", integrator);
  switch ( useJacobian ) {
     case JAC_OFF:         WriteAll("#JACOBIAN     - OFF\n");           break;
     case JAC_FULL:        WriteAll("#JACOBIAN     - FULL\n");          break;
     case JAC_LU_ROW:      WriteAll("#JACOBIAN     - SPARSE_LU_ROW\n"); break;
     case JAC_ROW:         WriteAll("#JACOBIAN     - SPARSE_ROW\n");    break;
  }
  switch ( useLang ) {
     case C_LANG:          WriteAll("#LANGUAGE     - C\n");             break;
     case F90_LANG:        WriteAll("#LANGUAGE     - Fortran90\n");     break;
     case MATLAB_LANG:     WriteAll("#LANGUAGE     - matlab\n");        break;
  }
  if( useMex )           { WriteAll("#MEX          - ON\n");                 }
     else                { WriteAll("#MEX          - OFF\n");                }
                           WriteAll("#MINVERSION   - %s\n", minKppVersion);
  if( useReorder )       { WriteAll("#REORDER      - ON\n");                 }
    else                 { WriteAll("#REORDER      - OFF\n");                }
  if( useStochastic )    { WriteAll("#STOCHASTIC   - ON\n");                 }
     else                { WriteAll("#STOCHASTIC   - OFF\n");                }
  if( useStoicmat )      { WriteAll("#STOICMAT     - ON\n");                 }
     else                { WriteAll("#STOICMAT     - OFF\n");                }
  if( upperCaseF90 )     { WriteAll("#UPPERCASEF90 - ON\n");                 }
     else                { WriteAll("#UPPERCASEF90 - OFF\n");                }

  NewLines(1);
  WriteAll("### Parameters ----------------------------------------\n");
  NewLines(1);

  VarStartNr = Index(0);
  FixStartNr = Index(VarNr);

  DeclareConstant( NSPEC,   ascii( SpcNr ) );
  DeclareConstant( NVAR,    ascii( max( VarNr, 1 ) ) );
  DeclareConstant( NVARACT, ascii( max( VarActiveNr, 1 ) ) );
  DeclareConstant( NFIX,    ascii( max( FixNr, 1 ) ) );
  DeclareConstant( NREACT,  ascii( EqnNr ) );
  DeclareConstant( NVARST,  ascii( VarStartNr ) );
  DeclareConstant( NFIXST,  ascii( FixStartNr ) );

  NewLines(1);
  WriteAll("### Species -------------------------------------------\n");

  NewLines(1);
  WriteAll("Variable species\n");

  dn = VarNr/3 + 1;
  for( i = 0; i < dn; i++ ) {
             if( i < VarNr ) WriteSpec( i, i );
    i += dn; if( i < VarNr ) WriteSpec( i, i );
    i += dn; if( i < VarNr ) WriteSpec( i, i );
    i -= 2*dn; WriteAll("\n");
  }


  NewLines(1);
  WriteAll("Fixed species\n");

  dn = FixNr/3 + 1;
  for( i = 0; i < dn; i++ ) {
             if( i < FixNr ) WriteSpec( i, i + VarNr );
    i += dn; if( i < FixNr ) WriteSpec( i, i + VarNr );
    i += dn; if( i < FixNr ) WriteSpec( i, i + VarNr );
    i -= 2*dn; WriteAll("\n");
  }

  NewLines(1);
  WriteAll("### Subroutines ---------------------------------------\n");
  NewLines(1);
  FlushBuf();
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateInitialize()
{
int i;
int I, X;
int INITVAL;

  if ( (useLang==C_LANG)||(useLang==F77_LANG)||(useLang==F90_LANG) )
      UseFile( initFile );

  INITVAL    = DefFnc( "Initialize",    0, "function to initialize concentrations");
  FunctionBegin( INITVAL );
  F77_Inline("      INCLUDE '%s_Global.h'", rootFileName);
  F90_Inline("  USE %s_Global\n", rootFileName);
  MATLAB_Inline("global CFACTOR VAR FIX NVAR NFIX", rootFileName);

  F77_Inline("      INCLUDE '%s_Parameters.h'", rootFileName);
  F90_Inline("  USE %s_Parameters\n", rootFileName);

  I = DefElm( "i", INT, 0);
  X = DefElm( "x", real, 0);
  Declare( I );
  Declare( X );

  NewLines(1);
  WriteComment("~~~ Define scale factor for units");
  WriteAssign( varTable[CFACTOR]->name , ascid( (double)cfactor ) );
  NewLines(1);

  //=========================================================================
  // MODIFICATION by Bob Yantosca (25 Apr 2022)
  //
  // For F90, assign values to C directly (since VAR and FIX now point to C
  // instead of the other way around).  Otherwise, preserve prior code.
  //=========================================================================
  if ( useLang == F90_LANG ) {

    //
    // Fortran-90
    //
    WriteComment("~~~ Zero C array");
    if ( useDouble )
      F90_Inline( "  C = 0.0_dp" );
    else
      F90_Inline( "  C = 0.0" );
    NewLines(1);

    WriteComment("~~~ Set initial species concentrations");
    for( i = 0; i < SpcNr; i++) {
      if( *SpeciesTable[ Code[i] ].ival == 0 ) continue;
      Assign( Elm( C, i ),
              Mul( Elm( IV, SpeciesTable[ Code[i] ].ival ),
              Elm( CFACTOR ) ) );
    }

    NewLines(1);

  } else {

    //
    // Fortran-77, C, and Matlab
    //
    Assign( Elm( X ), Mul( Elm( IV, varDefault ), Elm( CFACTOR ) ) );
    C_Inline("  for( i = 0; i < NVAR; i++ )" );
    F77_Inline("      DO i = 1, NVAR" );
    MATLAB_Inline("   for i = 1:NVAR" );
    ident++;
    Assign( Elm( VAR, -I ), Elm( X ) );
    ident--;
    F77_Inline("      END DO" );
    MATLAB_Inline("   end" );

    NewLines(1);
    Assign( Elm( X ), Mul( Elm( IV, fixDefault ), Elm( CFACTOR ) ) );
    C_Inline("  for( i = 0; i < NFIX; i++ )" );
    F77_Inline("      DO i = 1, NFIX" );
    MATLAB_Inline("   for i = 1:NFIX" );
    ident++;
    Assign( Elm( FIX, -I ), Elm( X ) );
    ident--;
    F77_Inline("      END DO" );
    MATLAB_Inline("   end" );

    NewLines(1);

    for( i = 0; i < VarNr; i++) {
      if( *SpeciesTable[ Code[i] ].ival == 0 ) continue;
      Assign( Elm( VAR, i ),
              Mul( Elm( IV, SpeciesTable[ Code[i] ].ival ),
              Elm( CFACTOR ) ) );
    }

    for( i = 0; i < FixNr; i++) {
      if( *SpeciesTable[ Code[i + VarNr] ].ival == 0 ) continue;
      Assign( Elm( FIX, i ),
              Mul( Elm( IV, SpeciesTable[ Code[i + VarNr] ].ival ),
              Elm( CFACTOR ) ) );
    }
  }

  WriteComment("constant rate coefficients");
  for( i = 0; i < EqnNr; i++) {
    if ( kr[i].type == NUMBER )
       Assign( Elm( RCONST, i ), Const( kr[i].val.f ) );
  }
  WriteComment("END constant rate coefficients");

  NewLines(1);
  WriteComment("INLINED initializations");

  switch( useLang ) {
    case C_LANG: bprintf( InlineCode[ C_INIT ].code );
                 break;
    case F77_LANG:  bprintf( InlineCode[ F77_INIT ].code );
                 break;
    case F90_LANG:  bprintf( InlineCode[ F90_INIT ].code );
                 break;
    case MATLAB_LANG: bprintf( InlineCode[ MATLAB_INIT ].code );
                 break;
  }
  FlushBuf();

  NewLines(1);
  WriteComment("End INLINED initializations");
  NewLines(1);

  MATLAB_Inline("   VAR = VAR(:);\n   FIX = FIX(:);\n" );

  FreeVariable( X );
  FreeVariable( I );
  FunctionEnd( INITVAL );
  FreeVariable( INITVAL );
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateShuffle_user2kpp()
{
int i,k;
int Shuffle_user2kpp;

  UseFile( utilFile );

  Shuffle_user2kpp    = DefFnc( "Shuffle_user2kpp", 2, "function to copy concentrations from USER to KPP");
  FunctionBegin( Shuffle_user2kpp, V_USER, V );

  k = 0;
  for( i = 1; i < SpcNr; i++) {
    if( ReverseCode[i] < 0 ) {
      if( SpeciesTable[i].type == VAR_SPC ) k++;
      continue;
    }
    switch( SpeciesTable[i].type ) {
      case VAR_SPC:
                     if( k < initNr ) {
                       Assign( Elm( V, ReverseCode[i] ), Elm( V_USER, k++ ) );
                       break;
                     }
      case FIX_SPC:
      case DUMMY_SPC:
      default:       break;
    }
  }

  FunctionEnd( Shuffle_user2kpp );
  FreeVariable( Shuffle_user2kpp );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateShuffle_kpp2user()
{
int i,k;
int Shuffle_kpp2user;

  UseFile( utilFile );

  Shuffle_kpp2user    = DefFnc( "Shuffle_kpp2user", 2, "function to restore concentrations from KPP to USER");
  FunctionBegin( Shuffle_kpp2user, V, V_USER );

  k = 0;
  for( i = 0; i < SpcNr; i++) {
    if( ReverseCode[i] < 0 ) {
      if( SpeciesTable[i].type == VAR_SPC ) k++;
      continue;
    }
    switch( SpeciesTable[i].type ) {
      case VAR_SPC:
                     if( k < initNr )
                       Assign( Elm( V_USER, k++ ), Elm( V, ReverseCode[i] ) );
                     break;
      case FIX_SPC:
      case DUMMY_SPC:
      default:       break;
    }
  }

  FunctionEnd( Shuffle_kpp2user );
  FreeVariable( Shuffle_kpp2user );
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateGetMass()
{
int i;
int atm, spc;
int GETMASS, MASS;
SPECIES_DEF *sp;
int numass;

  UseFile( utilFile );

  nmass = 0;
  for( atm = 0; atm < AtomNr; atm++ )
    if( AtomTable[atm].masscheck ) nmass++;
  if( nmass == 0 ) nmass = 1;

  MASS  = DefvElm( "Mass", real, nmass, "value of mass balance" );
  GETMASS = DefFnc( "GetMass", 2, "compute total mass of selected atoms");
  FunctionBegin( GETMASS, CL, MASS);

  numass = 0;
  for( atm = 0; atm < AtomNr; atm++ ) {
    if( AtomTable[atm].masscheck ) {
      sum = Const( 0 );
      for( spc = 0; spc < SpcNr; spc++ ) {
        sp = &SpeciesTable[ Code[spc] ];
        for( i = 0; i < sp->nratoms; i++ ) {
          if( sp->atoms[i].code == atm ) {
            sum = Add( sum, Mul( Const( sp->atoms[i].nr ),
                                 Elm( CL, spc ) ) );
          }
        }
      }
      Assign( Elm( MASS, numass ), sum );
      numass++;
    }
  }

  FunctionEnd( GETMASS );
  FreeVariable( GETMASS );
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateMakefile()
{
char buf[100];

  if ( useLang == MATLAB_LANG ) return;

  sprintf( buf, "Makefile_%s", rootFileName );
  makeFile = fopen(buf, "w");
  if( makeFile == 0 ) {
    FatalError(3,"%s: Can't create file", buf );
  }

  UseFile( makeFile );

  IncludeCode( "%s/util/Makefile", Home );

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateMex()
{
char buf[100], suffix[5];

  if (useLang == MATLAB_LANG) return;
  if (useMex == 0) return;

  // Because this CASE statement is used to create the _Mex* files,
  // use the f90Suffix variable to use the proper .f90 or .F90 extension,
  // depending on the value of the #UPPERCASEF90 switch.
  //   -- Bob Yantosca (22 Apr 2022)
  switch( useLang ) {
    case F77_LANG: sprintf( suffix, "%s", "f"       ); break;
    case F90_LANG: sprintf( suffix, "%s", f90Suffix ); break;
    case C_LANG:   sprintf( suffix, "%s", "c"       ); break;
    default:
      printf("\nCannot create mex files for language %d\n", useLang);
      exit(1);
  }

  sprintf( buf, "%s_mex_Fun.%s", rootFileName, suffix );
  mex_funFile = fopen(buf, "w");
  if( mex_funFile == 0 ) {
    FatalError(3,"%s: Can't create file", buf );
  }
  UseFile( mex_funFile );
  IncludeCode( "%s/util/Mex_Fun", Home );

  if (useJacSparse) {
    sprintf( buf, "%s_mex_Jac_SP.%s", rootFileName, suffix );
    mex_jacFile = fopen(buf, "w");
    if( mex_jacFile == 0 ) {
      FatalError(3,"%s: Can't create file", buf );
    }
    UseFile( mex_jacFile );
    IncludeCode( "%s/util/Mex_Jac_SP", Home );
  }

  if (useHessian) {
    sprintf( buf, "%s_mex_Hessian.%s", rootFileName, suffix );
    mex_hessFile = fopen(buf, "w");
    if( mex_hessFile == 0 ) {
      FatalError(3,"%s: Can't create file", buf );
    }
    UseFile( mex_hessFile );
    IncludeCode( "%s/util/Mex_Hessian", Home );
  }

}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateMatlabTemplates()
{
char buf[200];

  if (useLang != MATLAB_LANG) return;


  sprintf( buf, "%s_Fun_Chem.m", rootFileName );
  mex_funFile = fopen(buf, "w");
  if( mex_funFile == 0 ) {
    FatalError(3,"%s: Can't create file", buf );
  }
  UseFile( mex_funFile );
  IncludeCode( "%s/util/Template_Fun_Chem", Home );

  sprintf( buf, "%s_Update_SUN.m", rootFileName );
  mex_funFile = fopen(buf, "w");
  if( mex_funFile == 0 ) {
    FatalError(3,"%s: Can't create file", buf );
  }
  UseFile( mex_funFile );
  IncludeCode( "%s/util/UpdateSun", Home );

  if (useJacSparse) {
    sprintf( buf, "%s_Jac_Chem.m", rootFileName );
    mex_jacFile = fopen(buf, "w");
    if( mex_jacFile == 0 ) {
      FatalError(3,"%s: Can't create file", buf );
    }
    UseFile( mex_jacFile );
    IncludeCode( "%s/util/Template_Jac_Chem", Home );
  }

  if (useHessian) {
  }

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void GenerateF90Modules(char where)
{
char buf[200];

if (useLang != F90_LANG) return;

switch (where) {
case 'h':

  sprintf( buf, "%s_Precision.%s", rootFileName, f90Suffix );
  sparse_dataFile = fopen(buf, "w");
  if( sparse_dataFile == 0 ) {
    FatalError(3,"%s: Can't create file", buf );
  }
  UseFile( sparse_dataFile );
    F90_Inline("\nMODULE %s_Precision\n", rootFileName );
    F90_Inline("!");
    F90_Inline("! Definition of different levels of accuracy");
    F90_Inline("! for REAL variables using KIND parameterization");
    F90_Inline("!");
    F90_Inline("! KPP SP - Single precision kind");
    F90_Inline("  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)");
    F90_Inline("! KPP DP - Double precision kind");
    F90_Inline("  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)");
    F90_Inline("! KPP QP - Quadruple precision kind");
    F90_Inline("  INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(18,400)");
    F90_Inline("\nEND MODULE %s_Precision\n\n", rootFileName );

  UseFile( initFile );
    F90_Inline("MODULE %s_Initialize\n", rootFileName );
    F90_Inline("  USE %s_Parameters, ONLY: dp, NVAR, NFIX", rootFileName);
    F90_Inline("  IMPLICIT NONE\n", rootFileName );
    F90_Inline("CONTAINS\n\n");

  UseFile( param_headerFile );
    F90_Inline("MODULE %s_Parameters\n", rootFileName );
    F90_Inline("  USE %s_Precision", rootFileName );
    F90_Inline("  PUBLIC\n  SAVE\n");

  UseFile( global_dataFile );
    F90_Inline("MODULE %s_Global\n", rootFileName );
    if ( useDeclareValues ) {
      F90_Inline("  USE %s_Precision", rootFileName );
    }
    else {
      if ( useDouble ) {
              F90_Inline("  USE %s_Parameters, ONLY: dp, NSPEC, NVAR, NFIX, NREACT, LU_NONZERO", rootFileName);
      }
      else {
              F90_Inline("  USE %s_Parameters, ONLY: sp, NSPEC, NVAR, NFIX, NREACT, LU_NONZERO", rootFileName);
      }
    }

    F90_Inline("  PUBLIC\n  SAVE\n");

  UseFile( functionFile );
    F90_Inline("MODULE %s_Function\n", rootFileName );
    if ( useDeclareValues ) {
      F90_Inline("  USE %s_Precision", rootFileName );
    }
    else {
      if( doAutoReduce ) F90_Inline("  USE %s_Global, only : DO_FUN", rootFileName );
      F90_Inline("  USE %s_Parameters", rootFileName );
    }
    F90_Inline("  IMPLICIT NONE\n", rootFileName );
    Declare( A ); /*  mz_rs_20050117 */
    if ( useLang == F90_LANG || useLang == F77_LANG ) WriteOMPThreadPrivate(" A "); /* msl_20160419 */
    F90_Inline("\nCONTAINS\n\n");

  UseFile( rateFile );
    F90_Inline("MODULE %s_Rates\n", rootFileName );
    if ( useDeclareValues )
      F90_Inline("  USE %s_Precision", rootFileName );
    else
      F90_Inline("  USE %s_Parameters", rootFileName );
    F90_Inline("  USE %s_Global", rootFileName );
    F90_Inline("  IMPLICIT NONE", rootFileName );
    F90_Inline("  INTEGER, PARAMETER :: ASSOC = 1, DISSOC = 2", rootFileName );
    NewLines(1);
    IncludeCode( "%s/util/UserRateLawsInterfaces", Home );
    F90_Inline("\nCONTAINS\n\n");

  if ( useStochastic ) {
    UseFile(stochasticFile);
    F90_Inline("MODULE %s_Stochastic\n", rootFileName);
    if ( useDeclareValues )
      F90_Inline("  USE %s_Precision", rootFileName );
    else
      F90_Inline("  USE %s_Parameters, ONLY: NVAR, NFIX, NREACT", rootFileName );
    F90_Inline("  PUBLIC\n  SAVE\n");
    F90_Inline("\nCONTAINS\n\n");
  }

  if ( useJacSparse ) {
    UseFile(sparse_jacFile);
    F90_Inline("MODULE %s_JacobianSP\n", rootFileName);
    F90_Inline("  PUBLIC\n  SAVE\n");
  }

  UseFile( jacobianFile );
    F90_Inline("MODULE %s_Jacobian\n", rootFileName );
    if ( useDeclareValues ) {
      F90_Inline("  USE %s_Precision", rootFileName );
    }
    else {
      if( doAutoReduce ) F90_Inline("  USE %s_Global, ONLY: DO_JVS", rootFileName);
      F90_Inline("  USE %s_Parameters", rootFileName );
    }
    if ( useJacSparse )
      F90_Inline("  USE %s_JacobianSP\n", rootFileName);
    F90_Inline("  IMPLICIT NONE", rootFileName );
    F90_Inline("\nCONTAINS\n\n");

  if ( useStoicmat ) {
    UseFile(sparse_stoicmFile);
    F90_Inline("MODULE %s_StoichiomSP\n", rootFileName);
    F90_Inline("  USE %s_Precision", rootFileName);
    F90_Inline("  PUBLIC\n  SAVE\n");

    UseFile( stoichiomFile );
    F90_Inline("MODULE %s_Stoichiom\n", rootFileName);
    if ( useDeclareValues )
      F90_Inline("  USE %s_Precision", rootFileName );
    else
      F90_Inline("  USE %s_Parameters", rootFileName );
    F90_Inline("  USE %s_StoichiomSP\n", rootFileName);
    F90_Inline("  IMPLICIT NONE", rootFileName );
    F90_Inline("\nCONTAINS\n\n");
  }

  if ( useHessian ) {
    UseFile(sparse_hessFile);
    F90_Inline("MODULE %s_HessianSP\n", rootFileName);
    F90_Inline("  PUBLIC\n  SAVE\n");

    UseFile( hessianFile );
    F90_Inline("MODULE %s_Hessian\n", rootFileName);
    if ( useDeclareValues )
      F90_Inline("  USE %s_Precision", rootFileName );
    else
      F90_Inline("  USE %s_Parameters", rootFileName );
    F90_Inline("  USE %s_HessianSP\n", rootFileName);
    F90_Inline("  IMPLICIT NONE", rootFileName );
    F90_Inline("\nCONTAINS\n\n");
  }

   UseFile( monitorFile );
     F90_Inline("MODULE %s_Monitor", rootFileName);

   UseFile( linalgFile );
    F90_Inline("MODULE %s_LinearAlgebra\n", rootFileName);
    if( doAutoReduce ) F90_Inline("  USE %s_Global, ONLY: DO_SLV", rootFileName);
    F90_Inline("  USE %s_Parameters", rootFileName );
    if ( useJacSparse )
      F90_Inline("  USE %s_JacobianSP\n", rootFileName);
    F90_Inline("  IMPLICIT NONE", rootFileName );
    F90_Inline("\nCONTAINS\n\n");

   UseFile( utilFile );
    F90_Inline("MODULE %s_Util\n", rootFileName);
    F90_Inline("  USE %s_Parameters", rootFileName );
    F90_Inline("  IMPLICIT NONE", rootFileName );
    F90_Inline("\nCONTAINS\n\n");

    /* Here we define the model module which aggregates everything */
    /* put module rootFileName_Model into separate file */
    /* (reusing "sparse_dataFile" as done above for _Precision file) */
    //sprintf( buf, "%s_Model.f90", rootFileName );
    sprintf( buf, "%s_Model.%s", rootFileName, f90Suffix );
    sparse_dataFile = fopen(buf, "w");
    if( sparse_dataFile == 0 ) {
      FatalError(3,"%s: Can't create file", buf );
    }
    UseFile( sparse_dataFile );
    F90_Inline("MODULE %s_Model\n", rootFileName);
    F90_Inline("!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    F90_Inline("!  Completely defines the model %s",  rootFileName);
    F90_Inline("!    by using all the associated modules");
    F90_Inline("!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    F90_Inline("\n  USE %s_Precision", rootFileName );
    F90_Inline("  USE %s_Parameters", rootFileName );
    F90_Inline("  USE %s_Global", rootFileName );
    F90_Inline("  USE %s_Function", rootFileName );
    F90_Inline("  USE %s_Integrator", rootFileName );
    F90_Inline("  USE %s_Rates", rootFileName );
    if ( useStochastic )
       F90_Inline("  USE %s_Stochastic", rootFileName );
    if ( useJacobian )
       F90_Inline("  USE %s_Jacobian", rootFileName );
    if ( useHessian )
       F90_Inline("  USE %s_Hessian", rootFileName);
    if ( useStoicmat )
       F90_Inline("  USE %s_Stoichiom", rootFileName);
    F90_Inline("  USE %s_LinearAlgebra", rootFileName);
    F90_Inline("  USE %s_Monitor", rootFileName);
    F90_Inline("  USE %s_Util", rootFileName);
    F90_Inline("\nEND MODULE %s_Model\n", rootFileName);

  break;

case 't':

  UseFile( initFile );
  F90_Inline("\nEND MODULE %s_Initialize\n", rootFileName );

  UseFile( param_headerFile );
  F90_Inline("\nEND MODULE %s_Parameters\n", rootFileName );

  UseFile( global_dataFile );
  F90_Inline("\nEND MODULE %s_Global\n", rootFileName );

  UseFile( functionFile );
  F90_Inline("\nEND MODULE %s_Function\n", rootFileName );

  UseFile( rateFile );
  F90_Inline("\nEND MODULE %s_Rates\n", rootFileName );

  if ( useStochastic ) {
    UseFile(stochasticFile);
    F90_Inline("\nEND MODULE %s_Stochastic\n", rootFileName);
  }

  if ( useJacSparse ) {
    UseFile(sparse_jacFile);
    F90_Inline("\nEND MODULE %s_JacobianSP\n", rootFileName);
  }

  UseFile( jacobianFile );
  F90_Inline("\nEND MODULE %s_Jacobian\n", rootFileName );

  if ( useStoicmat ) {
    UseFile(sparse_stoicmFile);
    F90_Inline("\nEND MODULE %s_StoichiomSP\n", rootFileName);

    UseFile( stoichiomFile );
    F90_Inline("\nEND MODULE %s_Stoichiom\n", rootFileName);
  }

  if ( useHessian ) {
    UseFile(sparse_hessFile);
    F90_Inline("\nEND MODULE %s_HessianSP\n", rootFileName);

    UseFile( hessianFile );
    F90_Inline("\nEND MODULE %s_Hessian\n", rootFileName );
  }

   UseFile(monitorFile);
     F90_Inline("\nEND MODULE %s_Monitor", rootFileName);

   UseFile( linalgFile );
    F90_Inline("\nEND MODULE %s_LinearAlgebra\n", rootFileName);

   UseFile( utilFile );
    F90_Inline("\nEND MODULE %s_Util\n", rootFileName);

  break;

default:
  printf("\n Unrecognized option '%c' in GenerateF90Modules\n", where);
  break;
}
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void Generate()
{
int i;
int n;

  VarStartNr = 0;
  FixStartNr = VarNr;

  real = useDouble ? DOUBLE : REAL;

  n = MAX_OUTBUF;
  for( i = 1; i < INLINE_OPT; i++ )
    if( InlineCode[i].maxlen > n )
      n = InlineCode[i].maxlen;

  outBuf = (char*)malloc( n );
  outBuffer = outBuf;

  switch( useLang ) {
    case F77_LANG: Use_F( rootFileName );
                 break;
    case F90_LANG: Use_F90( rootFileName );
                 break;
    case C_LANG: Use_C( rootFileName );
                 break;
    case MATLAB_LANG: Use_MATLAB( rootFileName );
                 break;
    default: printf("\n Language no '%d' unknown\n",useLang );
  }
  printf("\nKPP is initializing the code generation.");
  InitGen();

  if ( useLang == F90_LANG )
      GenerateF90Modules('h');

  GenerateLog();

  printf("\nKPP is generating the monitor data:");
  printf("\n    - %s_Monitor",rootFileName);
  GenerateMonitorData();

  printf("\nKPP is generating the utility data:");
  printf("\n    - %s_Util",rootFileName);
  GenerateUtil();

  printf("\nKPP is generating the global declarations:");
  printf("\n    - %s_Main",rootFileName);
  GenerateGData();

  printf("\nKPP is generating the ODE function:");
  printf("\n    - %s_Function",rootFileName);
  GenerateFun(1); /* setting useAggregate=1, generate SUBROUTINE Fun*/
  GenerateFun(0); /* setting useAggregate=0, generate SUBROUTINE FUN_SPLIT */
  GenerateStoichNum();

  if ( useStochastic ) {
    printf("\nKPP is generating the Stochastic description:");
    printf("\n    - %s_Function",rootFileName);
    GenerateStochastic();
  }

  if ( useJacobian ) {
    printf("\nKPP is generating the ODE Jacobian:");
    printf("\n    - %s_Jacobian\n    - %s_JacobianSP",rootFileName,rootFileName);
    GenerateJacobianSparseData();
    GenerateJac();
    if ( (useLang == F77_LANG)||(useLang == F90_LANG)||(useLang == C_LANG) ) {
       GenerateJacVect();
       GenerateJacTRVect();
       if( useJacSparse ) {
          printf("\nKPP is generating the linear algebra routines:");
          printf("\n    - %s_LinearAlgebra",rootFileName);
          GenerateSparseUtil();
          GenerateSolve();
          GenerateTRSolve();
       }
    }
 }

 GenerateBlas();

  if( useHessian ) {
    printf("\nKPP is generating the Hessian:");
    printf("\n    - %s_Hessian\n    - %s_HessianSP",rootFileName,rootFileName);
    GenerateHessian();
    GenerateHessianSparseData();
  }

  printf("\nKPP is generating the utility functions:");
  printf("\n    - %s_Util",rootFileName);

  GenerateInitialize();

  GenerateShuffle_user2kpp();
  GenerateShuffle_kpp2user();

  printf("\nKPP is generating the rate laws:");
  printf("\n    - %s_Rates",rootFileName);

  GenerateRateLaws();
  GenerateUpdateSun();
  GenerateUpdateRconst();
  GenerateUpdatePhoto();
  GenerateGetMass();

  printf("\nKPP is generating the parameters:");
  printf("\n    - %s_Parameters",rootFileName);

  GenerateParamHeader();

  printf("\nKPP is generating the global data:");
  printf("\n    - %s_Global",rootFileName);

  GenerateGlobalHeader();

  if ( (useLang == F77_LANG)||(useLang == C_LANG)||(useLang == MATLAB_LANG) ) {
    printf("\nKPP is generating the sparsity data:");
    if( useJacSparse ) {
        GenerateJacobianSparseHeader();
        printf("\n    - %s_JacobianSP",rootFileName);
        }
    if( useHessian ) {
        GenerateHessianSparseHeader();
        printf("\n    - %s_HessianSP",rootFileName);
        }
    }

  if ( useStoicmat ) {
    printf("\nKPP is generating the stoichiometric description files:");
    printf("\n    - %s_Stoichiom\n    - %s_StoichiomSP",rootFileName,rootFileName);
    GenerateReactantProd();
    GenerateJacReactantProd();
    GenerateStoicmSparseData();
    if ( (useLang == F77_LANG)||(useLang == C_LANG)||(useLang == MATLAB_LANG) )
        GenerateStoicmSparseHeader();
    GenerateDFunDRcoeff();
    GenerateDJacDRcoeff();
  }

  printf("\nKPP is generating the driver from %s.%s:", driver, f90Suffix);
  printf("\n    - %s_Main",rootFileName);

  if ( (useLang == F77_LANG)||(useLang == F90_LANG)||(useLang == C_LANG) )
    GenerateIntegrator();

  if( strcmp( driver, "none" ) != 0 )
    GenerateDriver();

  if ( (useLang == F77_LANG)||(useLang == F90_LANG)||(useLang == C_LANG) )
      GenerateMakefile();

  if ( useLang == F90_LANG )
      GenerateF90Modules('t');

  if ( useLang == MATLAB_LANG )
      GenerateMatlabTemplates();

  if ( (useLang == F77_LANG)||(useLang == F90_LANG)||(useLang == C_LANG) )
      GenerateMex();

  if( initFile )          fclose( initFile );
  if( driverFile )        fclose( driverFile );
  if( functionFile )      fclose( functionFile );
  if( global_dataFile )   fclose( global_dataFile );
  if( hessianFile )       fclose( hessianFile );
  if( integratorFile )    fclose( integratorFile );
  if( jacobianFile )      fclose( jacobianFile );
  if( linalgFile )        fclose( linalgFile );
  if( logFile )           fclose( logFile );
  if( makeFile )          fclose( makeFile );
  if( monitorFile )       fclose( monitorFile );
  if( mex_funFile )       fclose( mex_funFile );
  if( mex_jacFile )       fclose( mex_jacFile );
  if( mex_hessFile )      fclose( mex_hessFile );
  if( param_headerFile )  fclose( param_headerFile );
  if( rateFile )          fclose( rateFile );
  if( sparse_dataFile )   fclose( sparse_dataFile );
  if( sparse_jacFile )    fclose( sparse_jacFile );
  if( sparse_hessFile )   fclose( sparse_hessFile );
  if( sparse_stoicmFile ) fclose( sparse_stoicmFile );
  if( stoichiomFile )     fclose( stoichiomFile );
  if( utilFile )          fclose( utilFile );
  if( stochasticFile )    fclose( stochasticFile );

}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int* AllocIntegerVector(int n, char* message)
{
int* vec;
if ( ( vec=(int*)calloc(n,sizeof(int)) ) == NULL )
   FatalError(-30,"%s: Cannot allocate vector.",message);
return vec;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Allocates a matrix of integers */
int** AllocIntegerMatrix(int m, int n, char* message)
{
int** mat;
int i;
if ( (mat = (int**)calloc(m,sizeof(int*)))==NULL ) {
    FatalError(-30,"%s: Cannot allocate matrix.", message);
    }
for (i=0; i<m; i++)
  if ( (mat[i] = (int*)calloc(n,sizeof(int)))==NULL ) {
    FatalError(-30,"%s: Cannot allocate matrix[%d].", message, i);
    }
return mat;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Frees the memory allocated by AllocIntegerMatrix */
void FreeIntegerMatrix(int** mat, int m, int n)
{
int i;
for (i=0; i<m; i++)
  free(mat[i]);
free(mat);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* i = C-type index; returns language-appropriate index */
int Index( int i )
{
        switch( useLang ) {
           case C_LANG:
             return i;
           case F77_LANG:
             return i+1;
           case F90_LANG:
             return i+1;
           case MATLAB_LANG:
             return i+1;
           default: printf("\n Unknown language no %d\n",useLang);
             exit(1);
        }
}
