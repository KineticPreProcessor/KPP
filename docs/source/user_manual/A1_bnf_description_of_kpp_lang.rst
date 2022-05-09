.. _bnf-description:

###################################
BNF Description of the KPP Language
###################################

Following is the BNF-like specification of the language:

.. code-block:: text

    program ::=  module | module program

    module  ::=  section | command |inline_code

    section ::=  #ATOMS atom_definition_list                                  |
                 #CHECK atom_list                                             |
                 #DEFFIX species_definition_list                              |
                 #DEFVAR species_definition_list                              |
                 #EQUATIONS equation_list                                     |
                 #FAMILIES family_list                                        |
                 #INITVALUES initvalues_list                                  |
                 #LOOKAT species_list atom_list                               |
                 #LUMP lump_list                                              |
                 #MONITOR species_list atom_list                              |
                 #SETFIX species_list_plus                                    |
                 #SETVAR species_list_plus                                    |
                 #TRANSPORT species_list

    command ::=  #CHECKALL                                                    |
                 #DECLARE [ SYMBOL | VALUE ]                                  |
                 #DOUBLE [ ON | OFF ]                                         |
                 #DRIVER driver_name                                          |
                 #DUMMYINDEX [ ON | OFF ]                                     |
                 #EQNTAGS [ ON | OFF ]                                        |
                 #FUNCTION [ AGGREGATE | SPLIT ]                              |
                 #HESSIAN [ ON | OFF ]                                        |
                 #INCLUDE file_name                                           |
                 #INTEGRATOR integrator_name                                  |
                 #INTFILE integrator_name                                     |
                 #JACOBIAN [ OFF | FULL | SPARSE_LU_ROW | SPARSE_ROW ]        |
                 #LANGUAGE[ Fortran90 | Fortran77 | C | Matlab ]              |
                 #LOOKATALL                                                   |
                 #MEX [ ON | OFF ]                                            |
                 #MINVERSION minimum_version_number                           |
                 #MODEL model_name                                            |
                 #REORDER [ ON | OFF ]                                        |
                 #STOCHASTIC [ ON | OFF ]                                     |
                 #STOICHMAT [ ON | OFF ]                                      |
                 #TRANSPORTALL [ ON | OFF]                                    |
                 #UPPERCASEF90 [ ON | OFF ]

 inline_code ::= #INLINE inline_type
                 inline_code
		 #ENDINLINE

 atom_count ::=               integer atom_name                              |
                              atom_name

 atom_definition_list :=      atom_definition                                |
                              atom_definition_list                               |
 atom_list ::=                atom_name;                                     |
                              atom_name; atom_list

 equation ::=                 <equation_tag> expression = expression : rate; |
                              expression =  expression : rate;

 equation_list ::=            equation                                       |
                              equation equation_list

 equation_tag ::=             Alphanumeric expression, also including the
                              underscore. In scan.l it is defined as .   

 expression ::=          | term                    |   |
                           | term :math:`+`          |  |
                           | expression              |           |
                           | term :math:`-`          |           |
                           | expression              |           |
 initvalues_assignment   | species_name_plus =     |  |
 ::=                       | program_expression;     |           |
                           | = program_expression;   |           |
 initvalues_list ::=     | initvalues_assignment   |  |
                           | initvalues_assignment   |           |
                           | initvalues_list         |           |
 inline_type ::=         |         |  |

 lump ::=                | lump_sum :              |           |
                           | species_name;           |           |
 lump_list ::=           | lump                    |  |
                           | lump lump_list        |           |
 lump_sum ::=            | species_name            |  |
                           | species_name :math:`+`  |           |
                           | lump_sum                |           |
 rate ::=                | number                  |  |
                           | program_expression      |           |
 species_composition ::= | atom_count              |  |
                           | atom_count :math:`+`    |  |
                           | species_composition     |           |
                           |                           |           |
 species_definition ::=  | species_name =          |           |
                           | species_composition;    |           |
 species_definition_list | species_definition      |  |
 ::=                       |                           |           |
                           | species_definition      |           |
                           | species_definition_list |           |
 species_list ::=        | species_name;           |  |
                           | species_name;           |           |
                           | species_list            |           |
 species_list_plus ::=   | species_name_plus;      |  |
                           | species_name_plus;      |           |
                           | species_list_plus       |           |
 species_name ::=        | Alphanumeric expression,  |           |
                           | also including the        |           |
                           | underscore, starting with |           |
                           | a letter.                 |           |
                           | In , it is defined as .   |           |
                           | Its maximum length is 15. |           |
 species_name_plus ::=   | species_name            |  |
                           |                           |  |
 term ::=                | number species_name   |  |
                           | species_name            |  |

.. |image| image:: Figures/kpp-logo.pdf
   :width: 70.0%
