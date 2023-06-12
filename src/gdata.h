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

// Version numbers must be synchronized in CHANGELOG.md, src/gdata.h,
// and docs/source/conf.py
#define KPP_VERSION "3.0.2"

#ifndef _GDATA_H_
#define _GDATA_H_

#include <stdio.h>
#include <string.h>

// - Many limits can be changed here by adjusting the MAX_* constants
// - To increase the max size of inlined code (F90_GLOBAL etc.),
//   change MAX_INLINE in scan.h.
//
//   NOTES:
//   ------
//   (1) Note: MAX_EQN or MAX_SPECIES over 1023 causes a seg fault in CI build
//         -- Lucas Estrada, 10/13/2021
//
//   (2) MacOS has a hard limit of 65332 bytes for stack memory.  To make
//       sure that you are using this max amount of stack memory, add
//       "ulimit -s 65532" in your .bashrc or .bash_aliases script.  We must
//       also set smaller limits for MAX_EQN and MAX_SPECIES here so that we
//       do not exceed the avaialble stack memory (which will result in the
//       infamous "Segmentation fault 11" error).  If you are stll having
//       problems on MacOS then consider reducing MAX_EQN and MAX_SPECIES
//       to smaller values than are listed below.
//         -- Bob Yantosca (03 May 2022)
#ifdef MACOS
#define MAX_EQN        2000     // Max number of equations (MacOS only)
#define MAX_SPECIES    1000     // Max number of species   (MacOS only)
#else
#define MAX_EQN       11000     // Max number of equations
#define MAX_SPECIES    6000     // Max number of species
#endif
#define MAX_SPNAME       30     // Max char length of species name
#define MAX_IVAL         40     // Max char length of species ID ?
#define MAX_EQNTAG       32     // Max length of equation ID in eqn file
#define MAX_K          1000     // Max length of rate expression in eqn file
#define MAX_ATOMS        10     // Max number of atoms
#define MAX_ATNAME       10     // Max char length of atom name
#define MAX_ATNR        250     // Max number of atom tables
#define MAX_PATH        300     // Max char length of directory paths
#define MAX_FILES        20     // Max number of files to open
#define MAX_FAMILIES    300     // Max number of family definitions
#define MAX_MEMBERS     150     // Max number of family members
#define MAX_EQNLEN      300     // Max char length of equations

#define NO_CODE 	-1
#define max( x, y ) (x) > (y) ? (x) : (y)
#define min( x, y ) (x) < (y) ? (x) : (y)

#define IncName(x)   FileName((x),"MODELS","models","")
#define ModelName(x) FileName((x),"MODELS","models",".def")
#define IntegName(x) FileName((x),"INTEG","int",".def")

enum krtypes { NUMBER, EXPRESION, PHOTO };
enum table_modes { F_TEXT, FC_TEXT, C_TEXT, S_TEXT };
enum lang { NO_LANG, C_LANG, F77_LANG, F90_LANG, MATLAB_LANG };
enum inl_code {
  F77_GLOBAL,  F77_INIT,    F77_DATA,     F77_UTIL,      F77_RATES,
  F77_RCONST,  F90_GLOBAL,  F90_INIT,     F90_DATA,      F90_UTIL,
  F90_RATES,   F90_RCONST,  C_GLOBAL,     C_INIT,        C_DATA,
  C_UTIL,      C_RATES,     C_RCONST,     MATLAB_GLOBAL, MATLAB_INIT,
  MATLAB_DATA, MATLAB_UTIL, MATLAB_RATES, MATLAB_RCONST,
  INLINE_OPT
};

enum jacobian_format { JAC_OFF, JAC_FULL, JAC_LU_ROW, JAC_ROW };


typedef short int CODE;
typedef float EQ_VECT[ MAX_EQN ];

typedef struct {
                 char name[ MAX_ATNAME ];
                 char check;
                 char masscheck;
               } ATOM_DEF;

typedef struct {
                 unsigned char code;
                 unsigned char nr;
               } ATOM;

typedef struct {
		 char type;
		 char lookat;
		 char moni;
		 char trans;
                 short int nratoms;
		 char name[ MAX_SPNAME ];
                 char ival[ MAX_IVAL ];
                 ATOM atoms[ MAX_ATOMS ];
                 int flux; /* msl_290416 */
	       } SPECIES_DEF;

typedef struct {
                 char name[ MAX_SPNAME ];
                 char ival[ MAX_IVAL ];
                 int code;
                 unsigned char nr;
                 float coeff;
                 char  coeff_str;
               } MEMBER;

typedef struct {
		 char type;
                 short int nrmembers;
		 char name[ MAX_SPNAME ];
                 char ival[ MAX_IVAL ];
                 MEMBER members[ MAX_MEMBERS ];
	       } FAMILY_DEF;

typedef struct {
		 char type;
		 union {
		   char st[ MAX_K ];
		   float f;
		 } val;
                 char label[ MAX_EQNTAG ];
	       } KREACT;

typedef struct {
		 char * code;
		 int maxlen;
	       } ICODE;


extern int SpeciesNr;
extern int FamilyNr;
extern int EqnNr;
extern int SpcNr;
extern int AtomNr;
extern int VarNr;
extern int VarActiveNr;
extern int FixNr;
extern int plNr;
extern int VarStartNr;
extern int FixStartNr;
extern int Hess_NZ;
extern int LU_Jac_NZ;
extern int Jac_NZ;

extern int generateSD;

extern int initNr;
extern int xNr;
extern int yNr;
extern int zNr;

extern int falseSpcNr;

extern int useAggregate;
extern int useJacobian;
extern int useJacSparse;
extern int useHessian;
extern int useStoicmat;
extern int useDouble;
extern int useReorder;
extern int useMex;
extern int useDummyindex;
extern int useEqntags;
extern int useLang;
extern int useStochastic;
extern int doFlux;
extern int doAutoReduce;
extern int upperCaseF90;
extern char f90Suffix[4];
extern char minKppVersion[30]; // size [30] must be the same as in scanner.c

/* if useValues=1 KPP replaces parameters like NVAR etc.
       by their values in vector/matrix declarations */
extern int useDeclareValues;

extern char Home[ MAX_PATH ];
extern char integrator[ MAX_PATH ];
extern char driver[ MAX_PATH ];
extern char runArgs[  MAX_PATH ];

extern char *eqFileName;
extern char *rootFileName;

extern ATOM_DEF AtomTable[ MAX_ATNR ];
extern SPECIES_DEF SpeciesTable[ MAX_SPECIES ];
extern FAMILY_DEF FamilyTable[ MAX_FAMILIES ];
extern KREACT 	kr	 [ MAX_EQN ];
extern CODE 	ReverseCode[ MAX_SPECIES ];
extern CODE 	Code	 [ MAX_SPECIES ];
extern float** 	Stoich_Left;
extern float** 	Stoich;
extern float**  Stoich_Right;
extern float**  Prod_Coeff;
extern float**  Loss_Coeff;
extern int 	Reactive [ MAX_SPECIES ];

extern CODE *Prod_Spc[ MAX_EQN ];
extern CODE *Loss_Spc[ MAX_EQN ];

extern int **structB;
extern int **structJ;
extern int **LUstructJ;

extern ICODE InlineCode[ INLINE_OPT ];

extern char *fileList[ MAX_FILES ];
extern int fileNr;

extern char varDefault[ MAX_IVAL ];
extern char radDefault[ MAX_IVAL ];
extern char fixDefault[ MAX_IVAL ];
extern double cfactor;

// Prototypes for functions in scanner.c
void CmdFunction( char *cmd );
void CmdJacobian( char *cmd );
void CmdHessian( char *cmd );
void CmdDeclareValues( char *cmd );
void CmdDouble( char *cmd );
void CmdReorder( char *cmd );
void CmdMex( char *cmd );
void CmdDummyindex( char *cmd );
void CmdEqntags( char *cmd );
void CmdUse( char *cmd );
void CmdLanguage( char *cmd );
void CmdIntegrator( char *cmd );
void CmdDriver( char *cmd );
void CmdRun( char *cmd );
void CmdStochastic( char *cmd );
void CmdFlux( char *cmd );
void Generate();

char * FileName( char *name, char* env, char *dir, char *ext );

int*  AllocIntegerVector( int n, char* message );
int** AllocIntegerMatrix( int m, int n, char* message );
void FreeIntegerMatrix ( int** mat, int m, int n );
int Index( int i );

// Add function prototpyes flagged as missing by the gfortran compiler,
// in order to remove -Wimplicit-function-declaration warnings.
//   -- Bob Yantosca (27 Apr 2022)
//void FatalError( int status, char *fmt, ... );
int KppVersionIsTooOld();

#endif
