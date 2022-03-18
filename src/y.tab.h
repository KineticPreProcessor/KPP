/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    JACOBIAN = 258,
    DOUBLE = 259,
    FUNCTION = 260,
    DEFVAR = 261,
    DEFRAD = 262,
    DEFFIX = 263,
    SETVAR = 264,
    SETRAD = 265,
    SETFIX = 266,
    HESSIAN = 267,
    STOICMAT = 268,
    STOCHASTIC = 269,
    DECLARE = 270,
    INITVALUES = 271,
    EQUATIONS = 272,
    FAMILIES = 273,
    LUMP = 274,
    INIEQUAL = 275,
    EQNEQUAL = 276,
    EQNCOLON = 277,
    LMPCOLON = 278,
    LMPPLUS = 279,
    SPCPLUS = 280,
    SPCEQUAL = 281,
    FAMCOLON = 282,
    ATOMDECL = 283,
    CHECK = 284,
    CHECKALL = 285,
    REORDER = 286,
    MEX = 287,
    DUMMYINDEX = 288,
    EQNTAGS = 289,
    LOOKAT = 290,
    LOOKATALL = 291,
    TRANSPORT = 292,
    TRANSPORTALL = 293,
    MONITOR = 294,
    USES = 295,
    SPARSEDATA = 296,
    WRITE_ATM = 297,
    WRITE_SPC = 298,
    WRITE_MAT = 299,
    WRITE_OPT = 300,
    INITIALIZE = 301,
    XGRID = 302,
    YGRID = 303,
    ZGRID = 304,
    USE = 305,
    LANGUAGE = 306,
    INTFILE = 307,
    DRIVER = 308,
    RUN = 309,
    INLINE = 310,
    ENDINLINE = 311,
    PARAMETER = 312,
    SPCSPC = 313,
    INISPC = 314,
    INIVALUE = 315,
    EQNSPC = 316,
    EQNSIGN = 317,
    EQNCOEF = 318,
    RATE = 319,
    LMPSPC = 320,
    SPCNR = 321,
    ATOMID = 322,
    LKTID = 323,
    MNIID = 324,
    INLCTX = 325,
    INCODE = 326,
    SSPID = 327,
    EQNLESS = 328,
    EQNTAG = 329,
    EQNGREATER = 330,
    TPTID = 331,
    USEID = 332,
    FLUX = 333,
    PLSPC = 334
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 72 "scan.y" /* yacc.c:1909  */

  char str[80];

#line 138 "y.tab.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
