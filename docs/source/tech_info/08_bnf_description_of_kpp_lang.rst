.. _bnf-description:

###################################
BNF description of the KPP language
###################################

Following is the BNF-like specification of the KPP language:

.. code-block:: text

 program ::=                module | module program

 module ::=                 section | command |inline_code

 section ::=                #ATOMS atom_definition_list                           |
                            #CHECK atom_list                                      |
                            #DEFFIX species_definition_list                       |
                            #DEFVAR species_definition_list                       |
                            #EQUATIONS equation_list                              |
                            #FAMILIES family_list                                 |
                            #INITVALUES initvalues_list                           |
                            #LOOKAT species_list atom_list                        |
                            #LUMP lump_list                                       |
                            #MONITOR species_list atom_list                       |
                            #SETFIX species_list_plus                             |
                            #SETVAR species_list_plus                             |
                            #TRANSPORT species_list

 command ::=                #CHECKALL                                             |
                            #DECLARE [ SYMBOL | VALUE ]                           |
                            #DOUBLE [ ON | OFF ]                                  |
                            #DRIVER driver_name                                   |
                            #DUMMYINDEX [ ON | OFF ]                              |
                            #EQNTAGS [ ON | OFF ]                                 |
                            #FUNCTION [ AGGREGATE | SPLIT ]                       |
                            #HESSIAN [ ON | OFF ]                                 |
                            #INCLUDE file_name                                    |
                            #INTEGRATOR integrator_name                           |
                            #INTFILE integrator_name                              |
                            #JACOBIAN [ OFF | FULL | SPARSE_LU_ROW | SPARSE_ROW ] |
                            #LANGUAGE[ Fortran90 | Fortran77 | C | Matlab ]       |
                            #LOOKATALL                                            |
                            #MEX [ ON | OFF ]                                     |
                            #MINVERSION minimum_version_number                    |
                            #MODEL model_name                                     |
                            #REORDER [ ON | OFF ]                                 |
                            #STOCHASTIC [ ON | OFF ]                              |
                            #STOICHMAT [ ON | OFF ]                               |
                            #TRANSPORTALL [ ON | OFF ]                            |
                            #UPPERCASEF90 [ ON | OFF ]

 inline_code ::=            #INLINE inline_type
                            inline_code
		            #ENDINLINE

 atom_count ::=             integer atom_name                                     |
                            atom_name

 atom_definition_list :=    atom_definition                                       |
                            atom_definition_list

 atom_list ::=              atom_name;                                            |
                            atom_name; atom_list

 equation ::=               <equation_tag> expression = expression : rate;        |
                            expression =  expression : rate;

 equation_list ::=          equation                                              |
                            equation equation_list

 equation_tag ::=           Alphanumeric expression, also including the
                            underscore. In scan.l it is defined as
                            "[a-zA-Z_0-0]+".

 expression ::=             term                                                  |
                            term + expression                                     |
                            term - expression

 initvalues_assignment :=   species_name_plus = program_expression;               |
                            CFACTOR = program_expression

 initvalues_list ::=        initvalues_assignment                                 |
                            initvalues_assignment initvalues_list

 inline_type ::=            F90_RATES    | F90_RCONST    | F90_GLOBAL             |
                            F90_INIT     | F90_DATA      | F90_UTIL               |
                            F77_RATES    | F77_RCONST    | F77_GLOBAL             |
                            F77_INIT     | F77_DATA      | F77_UTIL               |
                            C_RATES      | C_RCONST      | C_GLOBAL               |
                            C_INIT       | C_DATA        | C_UTIL                 |
                            MATLAB_RATES | MATLAB_RCONST | MATLAB_GLOBAL          |
                            MATLAB_INIT  | MATLAB_DATA   | MATLAB_UTIL

 lump ::=                   lump_sum : species_name;

 lump_list ::=              lump                                                  |
                            lump lump_list

 lump_sum ::=               species_name                                          |
                            species_name + lump_sum

 rate ::=                   number                                                |
                            program_expression

 species_composition ::=    atom_count                                            |
                            atom_count + species_composition                      |
                            IGNORE

 species_definition ::=     species_name = species_composition;

 species_definition_list := species_definition                                    |
                            species_definition species_definition_list

 species_list ::=           species_name;                                         |
                            species_name; species_list

 species_list_plus ::=      species_name_plus;                                    |
                            species_name_plus; species_list_plus

 species_name ::=           Alphanumeric expression, also including the
                            underscore, starting with a letter.  In
                            scan.l it is defined as "[a-zA-Z_][a-ZA-Z_0-9]*".
                            Its maximum length is 32.

 species_name_plus ::=      species_name                                          |
                            VAR_SPEC                                              |
                            FIX_SPEC                                              |
			    ALL_SPEC

 term ::=                   number species_name                                   |
                            species_name                                          |
                            PROD                                                  |
                            hv
