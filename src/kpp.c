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


#include <stdlib.h>
#include <string.h>
#include "gdata.h"
#include "scan.h"

char *eqFileName;
char *rootFileName = "ff";
char Home[ MAX_PATH ] = ""; 

short int linStru[ MAX_SPECIES ];
short int colStru[ MAX_SPECIES ];
short int bestStru[ MAX_SPECIES ];
short int *Stru;

enum stru_criteria { UNSORT, LINSORT, COLSORT, BESTSORT };

void EqCopy( EQ_VECT e1, EQ_VECT e2 )
{
int i;

  for( i = 0; i < EqnNr; i++ ) e2[i] = e1[i];
}

int NoSort( const void *p1, const void *p2 )
{
  return -1;
}

int CodeCmp( const void *p1, const void *p2 )
{
CODE *c1, *c2;

  c1 = (CODE*)p1;  
  c2 = (CODE*)p2;  
  
  if ( *c1 < *c2 ) return -1;
  if ( *c1 > *c2 ) return 1;
  return 0;
}

int CodeRCmp( const void *p1, const void *p2 )
{
int rc1, rc2;
CODE *c1, *c2;

  c1 = (CODE*)p1;  
  c2 = (CODE*)p2;  

  rc1 = Reactive[ ReverseCode[ *c1 ] ]; 
  rc2 = Reactive[ ReverseCode[ *c2 ] ];
  if ( rc1 > rc2 ) return -1; 
  if ( rc1 < rc2 ) return 1;
  if ( *c1 < *c2 ) return -1;
  if ( *c1 > *c2 ) return 1;
  return 0;
}

int CodeSSCmp( const void *p1, const void *p2 )
{
  return -CodeRCmp(p1,p2);
}

int CodeSCmp( const void *p1, const void *p2 )
{
CODE *c1, *c2;
short int sc1, sc2;

  c1 = (CODE*)p1;  
  c2 = (CODE*)p2;  
  
  sc1 = Stru[ ReverseCode[ *c1 ] ];
  sc2 = Stru[ ReverseCode[ *c2 ] ];
  
  if ( sc1 > sc2 ) return 1; 
  if ( sc1 < sc2 ) return -1;
  if ( *c1 < *c2 ) return 1;
  if ( *c1 > *c2 ) return -1;
  return 0;
}

void UpdateStructJ()
{
int i,j,k;

  for ( i=0; i<VarNr; i++ )
    for ( j=0; j<VarNr; j++ )
      structJ[i][j]=(i==j)?1:0;
              
  for (i = 0; i < VarNr; i++)
    for (j = 0; j < VarNr; j++)
      for (k = 0; k < EqnNr; k++)
        if ( Stoich[i][k]*((Stoich_Left[j][k])?1:0) != 0.0 )
          structJ[i][j]=1;

  for ( i=0; i<VarNr; i++ )
    for ( j=0; j<VarNr; j++ )
      LUstructJ[i][j]=structJ[i][j];
              
}

int ComputeLUStructJ()
{
int i,j,k;
int nu,nl;

  for (j = 0; j < VarNr-1; j++) {
    for (i = j+1; i < VarNr; i++) {
      if( LUstructJ[i][j] ) {
        for (k = j; k < VarNr; k++) 
	   if ( LUstructJ[j][k] != 0 )
              LUstructJ[i][k] = 1; 
      }   
    }
  }  

  
  nu = 0; nl = 0;
  for (i = 0; i < VarNr; i++) 
    for (j = 0; j < VarNr; j++)
      if( LUstructJ[i][j] ) {
        if(i > j) nl++;
        if(i <= j) nu++;     
      }          

  return nu+nl;
}

int LUnonZero()
{
CODE v[MAX_SPECIES];
CODE *var;
int i,j;

  var = v;
  if( Stru != bestStru ) {
    for( i=0; i<VarNr; i++ )
      var[i] = Code[i];
    qsort( (void*)var, VarNr, sizeof(CODE), CodeSCmp );
  } else {
    var = bestStru;
  }
          
  for (i = 0; i < VarNr; i++)
    for (j = 0; j < VarNr; j++)
      LUstructJ[i][j] = structJ[ ReverseCode[var[i]] ][ ReverseCode[var[j]] ];  

  return ComputeLUStructJ();
}

void LinColSparsity()
{
int i,j,k;
int nlin, ncol;

  for ( i=0; i<VarNr; i++ )
    for ( j=0; j<VarNr; j++ )
      structJ[i][j]=(i==j)?1:0;
              
  for (i = 0; i < VarNr; i++)
    for (j = 0; j < VarNr; j++)
      for (k = 0; k < EqnNr; k++)
        if ( Stoich[i][k]*((Stoich_Left[j][k])?1:0) != 0.0 )
          structJ[i][j]=1;

  for ( i=0; i<VarNr; i++ ) {
    linStru[i] = 0;
    for (j = 0; j < VarNr; j++)
      linStru[i] += structJ[i][j];
  }

  for ( i=0; i<VarNr; i++ ) {
    colStru[i] = 0;
    for (j = 0; j < VarNr; j++)
      colStru[i] += structJ[j][i];
    colStru[i] *= linStru[i];  
  }

  Stru = linStru;
  nlin = LUnonZero();
  Stru = colStru;
  ncol = LUnonZero();
  if( nlin <= ncol ) { 
    Stru = linStru;
    LUnonZero();
  } 
}

void BestSparsity()
{
int i,j,k;
int cnz, lnz;
int best, crt;
int best_i;
int tmp;
int s;

  UpdateStructJ();
      
  for ( i=0; i<VarNr; i++ )
    bestStru[i] = Code[i];

  for ( s=0; s<VarNr-1; s++ ) {
    best = MAX_SPECIES*MAX_SPECIES; best_i = 0;
    for ( i=s; i<VarNr; i++ ) {
      cnz = 0;lnz = 0;
      for (j = s; j < VarNr; j++) {
        cnz += (LUstructJ[i][j]?1:0);
        lnz += (LUstructJ[j][i]?1:0);
      }
      crt = (cnz-1)*(lnz-1);
      if( crt < best ) {
        best = crt;
        best_i = i;
      }    
    }
    for ( i=0; i<VarNr; i++ ) {
      tmp = LUstructJ[s][i];
      LUstructJ[s][i] = LUstructJ[best_i][i];
      LUstructJ[best_i][i] = tmp;      
    }   
    for ( i=0; i<VarNr; i++ ) {
      tmp = LUstructJ[i][s];
      LUstructJ[i][s] = LUstructJ[i][best_i];
      LUstructJ[i][best_i] = tmp;      
    }   
    tmp = bestStru[s];
    bestStru[s] = bestStru[best_i];
    bestStru[best_i] = tmp;
    
    for (i = s+1; i < VarNr; i++) {
      if( LUstructJ[i][s] ) {
        for (k = s; k < VarNr; k++) 
          LUstructJ[i][k] += LUstructJ[s][k];
      }   
    }
  }

  Stru = bestStru;
}

void ReorderSpecies( int criteria )
{
CODE *var;
CODE *fix;
CODE *dummy;
EQ_VECT *tmpStoich_Left;
EQ_VECT *tmpStoich_Right;
EQ_VECT *tmpStoich;
CODE *tmpCode;
CODE *tmpReact;
int i, k;
int new;
int (*cmpVar)(const void *, const void *);
int (*cmpFix)(const void *, const void *);
int dummyNr;

  cmpVar = CodeCmp;
  cmpFix = CodeCmp;

  switch( criteria ) {
    case UNSORT:   cmpVar = useJacobian ? CodeRCmp : CodeCmp;
                   break;  
    case LINSORT:  cmpVar = useJacobian ? CodeSCmp : CodeCmp;
                   Stru = linStru;  
                   break;  
    case COLSORT:  cmpVar = useJacobian ? CodeSCmp : CodeCmp;
                   Stru = colStru;
                   break;  
    case BESTSORT: cmpVar = useJacobian ? NoSort : CodeCmp;
                   break;  
  }

  VarNr       = 0;
  VarActiveNr = 0;
  FixNr       = 0;
  dummyNr     = 0;
  plNr        = 0;

  var = (CODE*)malloc( SpcNr * sizeof(CODE) );
  fix = (CODE*)malloc( SpcNr * sizeof(CODE) );
  dummy = (CODE*)malloc( 5 * sizeof(CODE) );
  tmpStoich_Left = (EQ_VECT*)malloc( SpcNr * sizeof(EQ_VECT) );
  tmpStoich_Right = (EQ_VECT*)malloc( SpcNr * sizeof(EQ_VECT) );
  tmpStoich = (EQ_VECT*)malloc( SpcNr * sizeof(EQ_VECT) );
  tmpCode = (CODE*)malloc( SpcNr * sizeof(CODE) );
  tmpReact = (CODE*)malloc( SpcNr * sizeof(CODE) );

  for( i = 0; i < SpcNr; i++ ) {
    switch( SpeciesTable[ Code[i] ].type ) {
      case VAR_SPC:  var[ VarNr++ ] = Code[ i ];
		     break;
      case FIX_SPC:  fix[ FixNr++ ] = Code[ i ]; 
		     break;
      case DUMMY_SPC:dummy[ dummyNr++ ] = Code[ i ]; 
		     break;
    }
  }

  /*msl*/
  if (doFlux == 1) { /* determine number of PL species and map their index to CODE */
    for( i = 0; i < VarNr; i++ ) {
      /* index mapping */
      if ( SpeciesTable[ Code[i] ].flux ) SpeciesTable[ Code[i] ].flux = Index(plNr++);
    }
  }

  if( Stru != bestStru ) {
    qsort( (void*)var, VarNr, sizeof(CODE), cmpVar );
  } else {
    for( i = 0; i < SpcNr; i++ )
      var[i] = bestStru[i];
  }  
  qsort( (void*)fix, FixNr, sizeof(CODE), cmpFix );

  for( i = 0; i < SpcNr; i++ ) {
    EqCopy( Stoich_Left[i], tmpStoich_Left[i] );
    EqCopy( Stoich_Right[i], tmpStoich_Right[i] );
    EqCopy( Stoich[i], tmpStoich[i] );
    tmpCode[i] = Code[i];
    tmpReact[i] = Reactive[i];
  }

  SpcNr -= dummyNr;
  dummyNr = 0;

  k = 0;
  for( i = 0; i < VarNr; i++ ) {
    new = ReverseCode[ var[i] ];
    EqCopy( tmpStoich_Left[ new ], Stoich_Left[ k ] );
    EqCopy( tmpStoich_Right[ new ], Stoich_Right[ k ] );
    EqCopy( tmpStoich[ new ], Stoich[ k ] );
    Code[ k ] = tmpCode[ new ];
    Reactive[ k ] = tmpReact[ new ];
    if( Reactive[ k ] ) VarActiveNr++; 
    k++;
  }
  for( i = 0; i < FixNr; i++ ) {
    new = ReverseCode[ fix[i] ];
    EqCopy( tmpStoich_Left[ new ], Stoich_Left[ k ] );
    EqCopy( tmpStoich_Right[ new ], Stoich_Right[ k ] );
    EqCopy( tmpStoich[ new ], Stoich[ k ] );
    Code[ k ] = tmpCode[ new ];
    Reactive[ k ] = tmpReact[ new ];
    k++;
  }
  for( i = 0; i < dummyNr; i++ ) {
    new = ReverseCode[ dummy[i] ];
    EqCopy( tmpStoich_Left[ new ], Stoich_Left[ k ] );
    EqCopy( tmpStoich_Right[ new ], Stoich_Right[ k ] );
    EqCopy( tmpStoich[ new ], Stoich[ k ] );
    Code[ k ] = tmpCode[ new ];
    Reactive[ k ] = tmpReact[ new ];
    k++;
  }


  for( i = 0; i < SpcNr+dummyNr; i++ )
    ReverseCode[ Code[i] ] = i;

  free( tmpReact );
  free( tmpCode );
  free( tmpStoich );
  free( tmpStoich_Right );
  free( tmpStoich_Left );
  free( dummy );
  free( fix );
  free( var );   

  fflush(stdout);
}

/* Allocate Internal Arrays */
void  AllocInternalArrays( void )
{
int i;

if ( (Stoich_Left =(float**)calloc(MAX_SPECIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Stoich_Left.\n");

for (i=0; i<MAX_SPECIES; i++)    
    if ( (Stoich_Left[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Stoich_Left[%d]",i,MAX_SPECIES);
    }

if ( (Stoich_Right = (float**)calloc(MAX_SPECIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Stoich_Right.\n");

for (i=0; i<MAX_SPECIES; i++)    
    if ( (Stoich_Right[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Stoich_Right[%d].",i);
    }

if ( (Stoich = (float**)calloc(MAX_SPECIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Stoich.\n");

for (i=0; i<MAX_SPECIES; i++)    
    if ( (Stoich[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Stoich[%d].",i);
    }

if ( (Loss_Coeff =(float**)calloc(MAX_FAMILIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Loss_Coeff.\n");

 for (i=0; i<MAX_FAMILIES; i++)    
  if ( (Loss_Coeff[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Loss_Coeff[%d]",i,MAX_FAMILIES);
    }
if ( (Prod_Coeff = (float**)calloc(MAX_FAMILIES,sizeof(float*)))==NULL ) 
    FatalError(-30,"Cannot allocate Prod_Coeff.\n");

for (i=0; i<MAX_FAMILIES; i++)    
    if ( (Prod_Coeff[i] = (float*)calloc(MAX_EQN,sizeof(float)))==NULL ) {
        FatalError(-30,"Cannot allocate Prod_Coeff[%d].",i);
    }

}


/* Allocate Structure Arrays */
void  AllocStructArrays( void )
{
int i;


if ( (structB = (int**)calloc(EqnNr,sizeof(int*)))==NULL ) 
    FatalError(-30, "Cannot allocate structB.");

for (i=0; i<EqnNr; i++)    
    if ( (structB[i] =(int*) calloc(SpcNr,sizeof(int)))==NULL )
        FatalError(-30, "Cannot allocate structB[%d].\n",i);
    
if ( (structJ = (int**)calloc(SpcNr,sizeof(int*)))==NULL ) 
    FatalError(-30, "Cannot allocate structJ.");

for (i=0; i<SpcNr; i++)    
    if ( (structJ[i] =(int*) calloc(SpcNr,sizeof(int)))==NULL ) 
        FatalError(-30, "Cannot allocate structJ[%d].\n",i);
    
if ( (LUstructJ = (int**)calloc(SpcNr,sizeof(int*)))==NULL ) 
    FatalError(-30, "Cannot allocate LUstructJ.");

for (i=0; i<SpcNr; i++)    
    if ( (LUstructJ[i] = (int*)calloc(SpcNr,sizeof(int)))==NULL ) 
        FatalError(-30, "Cannot allocate LUstructJ[%d].\n",i);

}

/*******************************************************************/                    
void Postprocess( char * root )
{
char buf[ 200 ];

  if ( useLang == MATLAB_LANG ) {
 /*  Add rate function definitions as internal functions to the Update_RCONST file*/
    sprintf( buf, "cat %s_Update_RCONST.m %s_Rates.m > tmp; mv tmp %s_Update_RCONST.m;", 
        root, root, root ); 
    system( buf );	 
  }
} 

/*******************************************************************/

int KppVersionIsTooOld() {
  //
  // Determine if the current KPP version is older than the minimum
  // version number (specified by #MINVERSION in the equation file).
  //
  // Return value: 0 = Version is the same or later than min. version
  //               1 = Version is too old
  //
  //  -- Bob Yantosca (27 Apr 2022)
  //
  int v0, v1, v2, r0, r1, r2;

  // If no minimum version is specified, then exit w/ success.
  // Use strncmp to restrict the search to the 1st 4 characters
  // (indexed from 0) in order to make the test robust enough
  // for running in the Docker container for C-I testing.
  // We have found that just using strcmp( minKppVersion, "none" ),
  // or if ( minKppVersion == "none" ) would cause the C-I tests to
  // fail when running on Azure dev pipelines.
  //  -- Lucas Estrada and Bob Yantosca (28 Apr 2022)
  if ( strncmp( minKppVersion, "none", 3 ) == 0 ) return 0;

  // Read versions into major, minor, patch numbers
  sscanf( KPP_VERSION,   "%d.%d.%d", &v0, &v1, &v2 );  // Current version
  sscanf( minKppVersion, "%d.%d.%d", &r0, &r1, &r2 );  // Required version

  // Compare versions
  if ( v0 > r0 ) return 0;      // Major version # is higher, so not too old
  if ( v0 == r0 ) {             // Otherwise need to test minor version #
    if ( v1 > r1 ) return 0;    // Minor version # is higher, so not too old
    if ( v1 == r1 ) {           // Otherwise test patch version #
      if ( v2 >= r2 ) return 0; // Patch version # is OK, so not too old
    }
  }
  nError += 1;                  // Update error count
  return 1;                     // Otherwise version is too old
}

/*******************************************************************/                    
int main( int argc, char * argv[] )
{
int status;
char name[ 200 ];
char *p;
  
  AllocInternalArrays();

  p = getenv("KPP_HOME");
  if( p ) strcpy( Home, p );

  switch( argc ) {
    case 3: eqFileName = argv[1];
            rootFileName = argv[2];
            break;
    case 2: eqFileName = argv[1];
            strcpy( name, eqFileName );
            p = name + strlen(name);
            while( p > name ) {
              if( *p == '.') { 
                *p = '\0';
                break;
              }
              p--;  
            } 
	    rootFileName = name;
	    break;
    default: FatalError(1,"\nUsage :"
		          "\n        kpp <equations file> [output file]\n");
  }

  printf("\nThis is KPP-%s.\n", KPP_VERSION);

  printf("\nKPP is parsing the equation file.");
  status = ParseEquationFile( argv[1] );

  if( status ) FatalError(2,"%d errors and %d warnings encountered.", 
                           nError, nWarning ); 

  // Stop if this version of KPP is older than #MINVERSION
  // (or continue along if #MINVERSION is not supplied)
  //   -- Bob Yantosca (27 Apr 2022)
  status = KppVersionIsTooOld();
  if ( status != 0 ) {
    sprintf( name, "You are using KPP version %s,\nbut your mechanism requires version %s or later.\n", KPP_VERSION, minKppVersion);
    FatalError( -50, name );
  }

  /* Allocate some internal data structures */
  AllocStructArrays();

  printf("\nKPP is computing Jacobian sparsity structure.");
  ReorderSpecies( UNSORT );
  if (useReorder==1){
    BestSparsity(); 
    ReorderSpecies( BESTSORT );
    }
  UpdateStructJ();
  ComputeLUStructJ();

  if( initNr == -1 ) initNr = VarNr;

  printf("\nKPP is starting the code generation.");
  Generate( rootFileName );
  
  printf("\nKPP is starting the code post-processing.");
  Postprocess( rootFileName );
  
  printf("\n\nKPP has successfully created the model \"%s\".\n\n",rootFileName);

  if( nError ) exit(4);
  if( nWarning ) exit(5);
  
  exit(0);
}
