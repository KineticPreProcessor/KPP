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


#include "gdata.h"
#include "code.h"
#include <string.h>
#include <stdio.h>

#define MAX_LINE 300   // NOTE: Same as MAX_PATH in gdata.h

/* Setting LEN=32 avoids problems with long species names and long 
   equation tags. A consistent change in F90_DeclareData (see below) 
   is probably also necessary.                                        */
char *F90_types[] = { "",                   /* VOID */ 
                      "INTEGER",            /* INT */
                      "REAL(kind=sp)",      /* FLOAT */
                      "REAL(kind=dp)",      /* DOUBLE */
                      "CHARACTER(LEN=32)",  /* STRING */
                      "CHARACTER(LEN=100)", /* DOUBLESTRING */
                      "LOGICAL"             /* LOGICAL */
                  };

/*************************************************************************************************/
void F90_WriteElm( NODE * n )
{
ELEMENT *elm;
char * name;
char maxi[20];
char maxj[20];

  elm = n->elm;
  name = varTable[ elm->var ]->name;

  switch( n->type ) {
    case CONST: bprintf("%g", elm->val.cnst);
		break;
    case ELM:   bprintf("%s", name);
		break;
    case VELM:  if( elm->val.idx.i >= 0 ) sprintf( maxi, "%d", elm->val.idx.i+1 );
                  else sprintf( maxi, "%s", varTable[ -elm->val.idx.i ]->name );
                bprintf("%s(%s)", name, maxi );
		break;
    case MELM:  if( elm->val.idx.i >= 0 ) sprintf( maxi, "%d", elm->val.idx.i+1 );
                  else sprintf( maxi, "%s", varTable[ -elm->val.idx.i ]->name );
                if( elm->val.idx.j >= 0 ) sprintf( maxj, "%d", elm->val.idx.j+1 );
                  else sprintf( maxj, "%s", varTable[ -elm->val.idx.j ]->name );
                bprintf("%s(%s,%s)", name, maxi, maxj );
		break;
    case EELM:  bprintf("(%s)", elm->val.expr );
		break;
  }
}

/*************************************************************************************************/
void F90_WriteSymbol( int op )
{
  switch( op ) {
    case ADD:   bprintf("+");
                AllowBreak();
    		break;
    case SUB:   bprintf("-");
                AllowBreak();
    		break;
    case MUL:   bprintf("*");
                AllowBreak();
    		break;
    case DIV:   bprintf("/");
                AllowBreak();
    		break;
    case POW:   bprintf("power");
                break;
    case O_PAREN: bprintf("(");
                AllowBreak();
                break;
    case C_PAREN: bprintf(")");
                break;
    case NONE:
                break;
  }
}

/*************************************************************************************************/
void F90_WriteAssign( char *ls, char *rs )
{
int start;
int linelg;
int i, j;
int jfound;
char c;
int first;
int crtident;

/* if MAX_NO_OF_LINES is too small, KPP will split lines incorrectly */
int number_of_lines = 1, MAX_NO_OF_LINES = 2000;

/*  Operator Mapping: 0xaa = '*' | 0xab = '+' | 0xac = ','
                      0xad = '-' | 0xae ='.' | 0xaf = '/' */
/* char op_mult=0xaa, op_plus=0xab, op_minus=0xad, op_dot=0xae, op_div=0xaf; */
char op_plus='+', op_minus='-';  //, op_dot='.', op_div='/'. op_mult='*'; 

  crtident = 2 + ident * 2;
  bprintf("%*s%s = ", crtident, "", ls);
  start = strlen( ls ) + 2;
  linelg = 120 - crtident - start - 1; /* F90 max line length is 132 */

  first = 1;
  while( strlen(rs) > linelg ) {
    jfound = 0;
    if ( number_of_lines >= MAX_NO_OF_LINES ) {
     /* If a new line needs to be started.
          Note: the approach below will create erroneous code if the +/- is within a subexpression, e.g. for
          A*(B+C) one cannot start a new continuation line by splitting at the + sign */
     for( j=linelg; j>5; j-- ) /* split row here if +, -, or comma */
       if ( ( rs[j] == op_plus )||( rs[j] == op_minus )||( rs[j]==',' ) ) {
        jfound = 1; i=j; break;
	}
    }
    if ( ( number_of_lines < MAX_NO_OF_LINES )||( !jfound ) ) {
     for( i=linelg; i>10; i-- ) /* split row here if operator or comma */
       if ( ( rs[i] & 0x80 )||( rs[i]==',' ) ) {
        break;
	}
     if( i <= 10 ) {
       printf("\n Warning: double-check continuation lines for:\n   %s = %s\n",ls,rs);
       i = linelg;
     }
    }
    while ( rs[i-1] & 0x80 ) i--; /* put all operators on the next row */
    while ( rs[i] == ',' ) i++;   /* put commas on the current row */

    c = rs[i];
    rs[i] = 0;

    if ( first ) { /* first line in a split row */
      bprintf("%s", rs );
      linelg++;
      first = 0;
    } else {/* continuation line in a split row - but not last line*/
      bprintf("&\n     %*s&%s", start, "", rs );
      if ( jfound ) {
         bprintf("\n%*s%s = %s", crtident, "", ls, ls);
	 number_of_lines = 1;
	 }
    }
    rs[i] = c;
    rs += i;  /* jump to the first not-yet-written character */
    number_of_lines++;
  }

  if ( number_of_lines > MAX_NO_OF_LINES ) {
     printf("\n Warning: %d continuation lines for %s = ...",number_of_lines,ls);
     }

  if ( first ) bprintf("%s\n", rs );  /* non-split row */
          else bprintf("&\n     %*s&%s\n", start, "", rs ); /* last line in a split row */


  FlushBuf();
}


/*************************************************************************************************/
void F90_WriteComment( char *fmt, ... )
{
Va_list args;
int n;
char buf[ MAX_LINE ];

  Va_start( args, fmt );
  vsprintf( buf, fmt, args );
  va_end( args );
  /* remove trailing spaces */
  /* taken from http://www.cs.bath.ac.uk/~pjw/NOTES/ansi_c/ch10-idioms.pdf */
  for (n= strlen(buf) - 1; n >= 0; n--)
    if (buf[n] != ' ') break;
  buf[n + 1]= '\0';
  bprintf( "! %s\n", buf );
  FlushBuf();
}

/*************************************************************************************************/
void F90_WriteOMPThreadPrivate( char *fmt, ... )
{
Va_list args;
int n;
char buf[ MAX_LINE ];

  Va_start( args, fmt );
  vsprintf( buf, fmt, args );
  va_end( args );
  /* remove trailing spaces */
  /* taken from http://www.cs.bath.ac.uk/~pjw/NOTES/ansi_c/ch10-idioms.pdf */
  for (n= strlen(buf) - 1; n >= 0; n--)
    if (buf[n] != ' ') break;
  buf[n + 1]= '\0';
  bprintf( "  !$OMP THREADPRIVATE( %s )\n", buf );
  FlushBuf();
}

/*************************************************************************************************/
char * F90_Decl( int v )
{
static char buf[120];
VARIABLE *var;
char *baseType;
char maxi[20];
char maxj[20];

  var = varTable[ v ];
  baseType = F90_types[ var->baseType ];

  *buf = 0;

  switch( var->type ) {
    case ELM:
                sprintf( buf, "%s :: %s", baseType, var->name );
		break;
    case VELM:
                if( var->maxi > 0 ) sprintf( maxi, "%d", var->maxi );
                /*  else sprintf( maxi, "%s", varTable[ -var->maxi ]->name); */
                /*sprintf( buf, "%s, DIMENSION(%s) :: %s", baseType, maxi, var->name );*/
                if( var->maxi == 0 ) sprintf( maxi, "%d", 1 );
                /*  else sprintf( maxi, "%s", varTable[ -var->maxi ]->name);  */
                if ( var->maxi < 0 ) {
		    if (varTable[ -var->maxi ]->value < 0)
		      sprintf( maxi, "%s", varTable[ -var->maxi ]->name );
		    else
		      sprintf( maxi, "%d", (varTable[-var->maxi]->value)==0?
		           1:varTable[-var->maxi]->value );
		}

		//=============================================================
		// MODIFICATION by Bob Yantosca (28 Apr 2022)
		//
		// Modify the IF block so that F90 variables can be declared
		// with either the POINTER, TARGET, or OPTIONAL attribute.
		//
		if ( var->attr == ATTR_F90_PTR )
		  sprintf( buf, "%s, POINTER :: %s(:)", baseType,
			                                var->name );
		else if ( var->attr == ATTR_F90_TGT )
		  sprintf( buf, "%s, TARGET :: %s(%s)", baseType,
			                                var->name,
			                                maxi );
		else if ( var->attr == ATTR_F90_OPT )
		  sprintf( buf, "%s, OPTIONAL :: %s(%s)", baseType,
			                                  var->name,
			                                  maxi );
		else if ( var->attr == ATTR_F90_OPT_PTR )
		  sprintf( buf, "%s, OPTIONAL, POINTER :: %s(%s)", baseType,
			                                           var->name,
			                                           maxi );
		else
		  sprintf( buf, "%s :: %s(%s)", baseType, var->name, maxi );
		//=============================================================
 		break;

    case MELM:
                if( var->maxi > 0 ) sprintf( maxi, "%d", var->maxi );
                else {
		 if (varTable[ -var->maxi ]->value < 0)
		    sprintf( maxi, "%s", varTable[ -var->maxi ]->name );
		 else
		    sprintf( maxi, "%d", (varTable[-var->maxi]->value)==0?
		           1:varTable[-var->maxi]->value );
	        }
                if( var->maxj > 0 ) sprintf( maxj, "%d", var->maxj );
                else {
		 if (varTable[ -var->maxj ]->value < 0)
		    sprintf( maxj, "%s", varTable[ -var->maxj ]->name );
		 else
		    sprintf( maxj, "%d", (varTable[-var->maxj]->value)==0?
		           1:varTable[-var->maxj]->value );
	        }
                sprintf( buf, "%s :: %s(%s,%s)",
                         baseType, var->name, maxi, maxj );
		break;

    default:
                printf( "Can not declare type %d\n", var->type );
                break;
  }
  return buf;
}

/*************************************************************************************************/
char * F90_DeclareData( int v, void * values, int n)
{
int i;
int nlines;
static char buf[120];
VARIABLE *var;
int * ival;
double * dval;
char ** cval;
char *baseType;
char maxi[20];
char maxj[20];
int maxCols = MAX_COLS;

 int i_from, i_to;
 int isplit;
 int splitsize;
 int maxi_mod;
 int maxi_div;

 char mynumber[30];

  var = varTable[ v ];
  ival = (int*) values;
  dval = (double *) values;
  cval = (char **) values;

  nlines = 1;
  var -> maxi = max( n, 1 );

  baseType = F90_types[ var->baseType ];

  *buf = 0;

  switch( var->type ) {
    case ELM:
	    bprintf( "  %s :: %s = ", baseType, var->name );
		switch ( var->baseType ) {
		  case INT: bprintf( "%d", *ival ); break;
		  case DOUBLE: bprintf( "%f", *dval); break;
		  case REAL: bprintf( "%lg", *dval ); break;
		  case STRING: bprintf( "'%3s'", *cval ); break;
		}
		break;

    case VELM:
      /* define maxCols here already and choose suitable splitsize */
      switch( var -> baseType ) {
      case INT:          maxCols =12; break;
      case DOUBLE:       maxCols = 5; break;
      case REAL:         maxCols = 5; break;
      case STRING:       maxCols = 3; break;
      case DOUBLESTRING: maxCols = 1; break;
      }
      splitsize = 30 * maxCols; /* elements = lines * columns */
      maxi_mod = var->maxi % splitsize;
      maxi_div = var->maxi / splitsize;
      /* correction if var->maxi is a multiple of splitsize */
      if ( (maxi_div>0) && (maxi_mod==0) ) {
        maxi_mod = splitsize;
        maxi_div--;
      }
      for ( isplit=0; isplit <= maxi_div; isplit++ ) {
        if( var->maxi > 0 ) sprintf( maxi, "%d", var->maxi );
        else sprintf( maxi, "%s", varTable[ -var->maxi ]->name );
        if( (var->maxi == 0) ||
            ((var->maxi < 0) && (varTable[ -var->maxi ]->maxi == 0)) )
          strcat( maxi, "+1");
	bprintf( "  %s, " , baseType);
        if( n>0 ) bprintf( "PARAMETER, " ); /* if values are assigned now */
	//====================================================================
	// MODIFICATION by Bob Yantosca (28 Apr 2022)
	//
	// Add the POINTER, TARGET, or OPTIONAL attributes to F90 variables.
	//
	if ( var->attr == ATTR_F90_PTR ) {
          bprintf( ", POINTER :: %s(:)", var->name ) ;
	  if( n<=0 ) break;
	}
	if ( var->attr == ATTR_F90_TGT ) {
          bprintf( ", TARGET :: %s(%s)", var->name, maxi ) ;
	  if( n<=0 ) break;
        }
	if ( var->attr == ATTR_F90_OPT ) {
          bprintf( ", OPTIONAL :: %s(%s)", var->name, maxi ) ;
	  if( n<=0 ) break;
        }
	if ( var->attr == ATTR_F90_OPT_PTR ) {
          bprintf( ", OPTIONAL, POINTER :: %s(%s)", var->name, maxi ) ;
	  if( n<=0 ) break;
        }
	//====================================================================
        if ( maxi_div==0 ) { /* define array in one piece */
          bprintf( "DIMENSION(%s) :: %s",
                   maxi, var->name) ;
        } else {/* define partial arrays */
          if ( isplit==maxi_div ) { /* last part has size maxi_mod */
            bprintf( "DIMENSION(%d) :: %s_%d",
                     maxi_mod, var->name, isplit) ;
          } else { /* all other parts have size splitsize */
            bprintf( "DIMENSION(%d) :: %s_%d",
                     splitsize, var->name, isplit) ;
          }
        }
        if( n<=0 ) break;

        /* now list values */
        bprintf( " = (/ &\n     " );
        /*   if the array is defined in one piece, then the for loop will
                 go from 0 to n. Otherwise, there will be partial arrays from
                 i_from to i_to which are of size splitsize except for the
                 last one which is usually smaller and contains the rest */
        i_from = isplit * splitsize;
        i_to   = min(i_from+splitsize,n);
        for ( i=i_from; i < i_to; i++ ) {
          switch( var -> baseType ) {
          case INT:
            bprintf( "%3d", ival[i] ); break;
          case DOUBLE:
            /* bprintf( "%4f", dval[i] ); maxCols = 5; break; */
	    sprintf(mynumber, "%12.6e_dp",dval[i]);
            /* mynumber[ strlen(mynumber)-4 ] = 'd'; */
            bprintf( "  %s", mynumber ); break;
          case REAL:
            bprintf( "%12.6e", dval[i] ); break;
          case STRING:
            /* Setting length to 32 avoids problems with long species names and long equation tags. */
            /* A consistent change in "char *F90_types" (see above) is probably also necessary */
            bprintf( "'%-32s'", cval[i] ); break;
          case DOUBLESTRING:
            bprintf( "'%-100.100s'", cval[i] ); break;
          }
          if( i < i_to-1 ) {
            bprintf( "," );
            if( (i+1) % maxCols == 0 ) {
              if (maxCols == 1 ) {
                bprintf( " & ! index %d\n     ", i+1 ); }
              else {
                bprintf( " & ! index %d - %d\n     ", i-maxCols+2, i+1 ); }
              nlines++;
            }
          }
        }
        bprintf( " /) ! index up to %d\n", i_to );
        FlushBuf();
      }

      /* combine the partial arrays */
      if ( maxi_div != 0 ) {
        bprintf( "  %s, PARAMETER, DIMENSION(%s) :: %s = (/&\n    ",
                 baseType, maxi, var->name) ;
        for ( isplit=0; isplit <= maxi_div; isplit++ ) {
          bprintf( "%s_%d", var->name, isplit) ;
          if( isplit < maxi_div ) { /* more parts will follow */
            bprintf( ", " );
            /* line break after 5 variables */
            if( (isplit+1) % 5 == 0 ) bprintf( "&\n    " );
          } else { /* after last part */
            bprintf( " /)\n" );
          }
        }
      }

      break;

    case MELM:  if( var->maxi > 0 ) sprintf( maxi, "%d", var->maxi );
                  else sprintf( maxi, "%s", varTable[ -var->maxi ]->name );
                if( (var->maxi == 0) ||
                    ((var->maxi < 0) && (varTable[ -var->maxi ]->maxi == 0)) )
                  strcat( maxi, "+1");
                if( var->maxj > 0 ) sprintf( maxj, "%d", var->maxj );
                  else sprintf( maxj, "%s", varTable[ -var->maxj ]->name );
                if( (var->maxj == 0) ||
                    ((var->maxj < 0 ) && (varTable[ -var->maxj ]->maxi == 0)) )
                  strcat( maxj, "+1");
                sprintf( buf, "%s, DIMENSION(%s,%s) :: %s\n",	/* changed here */
                         baseType, maxi, maxj,var->name );
		break;

    default:
                printf( "Can not declare type %d", var->type );
                break;
  }
  return buf;
}

/*************************************************************************************************/
void F90_Declare( int v )
{
  if( varTable[ v ]->comment ) {
    F90_WriteComment( "%s - %s",
                    varTable[ v ]->name, varTable[ v ]->comment );
  }
  bprintf("  %s\n", F90_Decl(v) );

  FlushBuf();
}

/*************************************************************************************************/
void F90_ExternDeclare( int v )
{
  F90_Declare( v );
}

/*************************************************************************************************/
void F90_GlobalDeclare( int v )
{
}

/*************************************************************************************************/
void F90_DeclareConstant( int v, char *val )
{
VARIABLE *var;
int ival;
char dummy_val[100];           /* used just to avoid strange behaviour of
                                  sscanf when compiled with gcc */

  strcpy(dummy_val,val);val = dummy_val;

  var = varTable[ v ];

  if( sscanf(val, "%d", &ival) == 1 )
    if( ival == 0 ) var->maxi = 0;
               else var->maxi = 1;
  else
    var->maxi = -1;

  if( var->comment )
    F90_WriteComment( "%s - %s", var->name, var->comment );

  switch( var->type ) {
    case CONST: bprintf("  %s, PARAMETER :: %s = %s \n",
                       F90_types[ var->baseType ], var->name, val );
                break;
    default:
                printf( "Invalid constant %d", var->type );
                break;
  }

  FlushBuf();
}


/*************************************************************************************************/
void F90_WriteVecData( VARIABLE * var, int min, int max, int split )
{
char buf[80];

  if( split )
    sprintf( buf, "%6sdata( %s(i), i = %d, %d ) / &\n%5s",
                " ", var->name, min, max, " " );
  else
    sprintf( buf, "%6sdata %s / &\n%5s",
                    " ", var->name, " " );

  FlushThisBuf( buf );
  bprintf( " / \n\n" );
  FlushBuf();
}

/*************************************************************************************************/
void F90_DeclareDataOld( int v, int * values, int n )
{
int i;
int nlines, min, max;
int split;
VARIABLE *var;
int * ival;
double * dval;
char **cval;
int maxCols = MAX_COLS;
char dsbuf[55];

  var = varTable[ v ];
  ival = (int*) values;
  dval = (double*) values;
  cval = (char**) values;

  nlines = 1;
  min = max = 1;
  split = 0;

  switch( var->type ) {
     case VELM: if( n <= 0 ) break;
    	       for( i = 0; i < n; i++ ) {
                 switch( var->baseType ) {
                   case INT: bprintf( "%3d",  ival[i] ); maxCols=12; break;
                   case DOUBLE:
                   case REAL:bprintf( "%5lg", dval[i] ); maxCols=8; break;
                   case STRING:bprintf( "'%s'", cval[i] ); maxCols=5; break;
                   case DOUBLESTRING:
		        strncpy( dsbuf, cval[i], 54 ); dsbuf[54]='\0';
		        bprintf( "'%48s'", dsbuf ); maxCols=1; break;
                 }
                 if( ( (i+1) % 12 == 0 ) && ( nlines > MAX_LINES ) ) {
                     split = 1; nlines = 1;
                     F90_WriteVecData( var, min, max, split );
                     min = max + 1;
                 }
                 else {
                   if( i < n-1 ) bprintf( "," );
                   if( (i+1) % maxCols == 0 ) {
                     bprintf( "\n%5s", " " );
                     nlines++;
                   }
                 }
                 max ++;
               }
               F90_WriteVecData( var, min, max-1, split );
               break;

    case ELM:  bprintf( "%6sdata %s / ", " ", var->name );
               switch( var->baseType ) {
                 case INT: bprintf( "%d",  *ival ); break;
                 case DOUBLE:
                 case REAL:bprintf( "%lg", *dval ); break;
                 case STRING:bprintf( "'%s'", *cval ); break;
                 case DOUBLESTRING:
		        strncpy( dsbuf, *cval, 54 ); dsbuf[54]='\0';
		        bprintf( "'%s'", dsbuf ); maxCols=1; break;
               }
               bprintf( " / \n" );
               FlushBuf();
               break;
    default:
               printf( "\n Function not defined !\n" );
               break;
  }
}

/*************************************************************************************************/
void F90_InitDeclare( int v, int n, void * values )
{
VARIABLE * var;

  var = varTable[ v ];
  var->maxi = max( n, 1 );

  NewLines(1);
  F90_DeclareData( v, values, n );
}

/*************************************************************************************************/
void F90_FunctionStart( int f, int *vars )
{
int i;
int v;
char * name;
int narg;

  name = varTable[ f ]->name;
  narg = varTable[ f ]->maxi;

  bprintf("SUBROUTINE %s ( ", name );
  for( i = 0; i < narg-1; i++ ) {
    v = vars[ i ];
    bprintf("%s, ", varTable[ v ]->name );
  }
  if( narg >= 1 ) {
    v = vars[ narg-1 ];
    bprintf("%s ", varTable[ v ]->name );
  }
  bprintf(")\n");

  FlushBuf();
}

/*************************************************************************************************/
void F90_FunctionPrototipe( int f, ... )
{
char * name;

  name = varTable[ f ]->name;

  bprintf("      EXTERNAL %s\n", name );

  FlushBuf();
}

/*************************************************************************************************/
void F90_FunctionBegin( int f, ... )
{
Va_list args;
int i;
int vars[20];
int narg;

  narg = varTable[ f ]->maxi;

  Va_start( args, f );
  for( i = 0; i < narg; i++ )
    vars[ i ] = va_arg( args, int );
  va_end( args );

  CommentFncBegin( f, vars );
  F90_FunctionStart( f, vars );
  NewLines(1);

  FlushBuf();

  for( i = 0; i < narg; i++ )
    F90_Declare( vars[ i ] );

  bprintf("\n");
  FlushBuf();

  LogFunctionComment( f, vars );
}

/*************************************************************************************************/
void F90_FunctionEnd( int f )
{
  bprintf("      \nEND SUBROUTINE %s\n\n", varTable[ f ]->name );

  FlushBuf();

  CommentFunctionEnd( f );
}

/*************************************************************************************************/
void F90_Inline( char *fmt, ... )
{
va_list args;
char buf[ MAX_K ];

  if( useLang != F90_LANG ) return;

  va_start( args, fmt );
  vsprintf( buf, fmt, args );
  va_end( args );
  bprintf( "%s\n", buf );

  FlushBuf();

}

/*************************************************************************************************/
void Use_F90( char* rootFileName )
{
  // Temporary string variable
  char buf[ MAX_K ];

  // Use functions for writing code to the F90 format files
  WriteElm              = F90_WriteElm;
  WriteSymbol           = F90_WriteSymbol;
  WriteAssign           = F90_WriteAssign;
  WriteComment          = F90_WriteComment;
  WriteOMPThreadPrivate = F90_WriteOMPThreadPrivate;
  DeclareConstant       = F90_DeclareConstant;
  Declare               = F90_Declare;
  ExternDeclare         = F90_ExternDeclare;
  GlobalDeclare         = F90_GlobalDeclare;
  InitDeclare           = F90_InitDeclare;
  FunctionStart         = F90_FunctionStart;
  FunctionPrototipe     = F90_FunctionPrototipe;
  FunctionBegin         = F90_FunctionBegin;
  FunctionEnd           = F90_FunctionEnd;

  // Parameters
  sprintf( buf, "_Parameters.%s", f90Suffix );
  OpenFile( &param_headerFile, rootFileName, buf, "Parameter Module File" );

  // Initialize
  sprintf( buf, "_Initialize.%s", f90Suffix );
  OpenFile( &initFile, rootFileName, buf, "Initialization File" );

  // Main
  if( strcmp( driver, "none" ) != 0 ) {
    sprintf( buf, "_Main.%s", f90Suffix );
    OpenFile( &driverFile, rootFileName, buf, "Main Program File" );
  }

  // Integrator
  sprintf( buf, "_Integrator.%s", f90Suffix );
  OpenFile( &integratorFile, rootFileName, buf,
	    "Numerical Integrator (Time-Stepping) File" );

  // LinearAlgebra
  sprintf( buf, "_LinearAlgebra.%s", f90Suffix );
  OpenFile( &linalgFile, rootFileName, buf,
	    "Linear Algebra Data and Routines File" );

  // Function
  sprintf( buf, "_Function.%s", f90Suffix );
  OpenFile( &functionFile, rootFileName, buf,
	    "The ODE Function of Chemical Model File" );

  // Jacobian
  sprintf( buf, "_Jacobian.%s", f90Suffix );
  OpenFile( &jacobianFile, rootFileName, buf,
	    "The ODE Jacobian of Chemical Model File" );

  // _Rates
  sprintf( buf, "_Rates.%s", f90Suffix );
  OpenFile( &rateFile, rootFileName, buf, "The Reaction Rates File" );

  // _Stochastic
  if ( useStochastic ) {
    sprintf( buf, "_Stochastic.%s", f90Suffix );
    OpenFile( &stochasticFile, rootFileName, buf,
	      "The Stochastic Chemical Model File" );
  }

  // _Stoichiom and _StoichiomSP
  if ( useStoicmat ) {
    sprintf( buf, "_Stoichiom.%s", f90Suffix );
    OpenFile( &stoichiomFile, rootFileName, buf,
	      "The Stoichiometric Chemical Model File" );

    sprintf( buf, "_StoichiomSP.%s", f90Suffix );
    OpenFile( &sparse_stoicmFile, rootFileName, buf,
	      "Sparse Stoichiometric Data Structures File" );
  }

  // _Util
  sprintf( buf, "_Util.%s", f90Suffix );
  OpenFile( &utilFile, rootFileName, buf, "Auxiliary Routines File" );

  // _Global
  sprintf( buf, "_Global.%s", f90Suffix );
  OpenFile( &global_dataFile, rootFileName, buf, "Global Data Module File" );

  // _JacobianSP
  if ( useJacSparse ) {
    sprintf( buf, "_JacobianSP.%s", f90Suffix );
    OpenFile( &sparse_jacFile, rootFileName, buf,
	      "Sparse Jacobian Data Structures File" );
  }

  // _Hessian and _HessianSP
  if ( useHessian ) {
    sprintf( buf, "_Hessian.%s", f90Suffix );
    OpenFile( &hessianFile, rootFileName, buf, "Hessian File" );

    sprintf( buf, "_HessianSP.%s", f90Suffix );
    OpenFile( &sparse_hessFile, rootFileName, buf,
	      "Sparse Hessian Data Structures File" );
  }

  // .log
  OpenFile( &logFile, rootFileName, ".log",
	    "Log File with Human-Readable Information" );

  // _Monitor
  sprintf( buf, "_Monitor.%s", f90Suffix );
  OpenFile( &monitorFile, rootFileName, buf, "Utility Data Module File" );
}
