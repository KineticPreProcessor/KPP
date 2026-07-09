.. _kpp-internal-modules:

####################
KPP internal modules
####################

.. _scanner-parser:

==================
Scanner and parser
==================

This module is responsible for reading the kinetic description files and
extracting the information necessary in the code generation phase. We
make use of the flex and bison generic tools in implementing our own
scanner and parser. Using these tools, this module gathers information
from the input files and fills in the following data structures in
memory:

-  The atom list

-  The species list

-  The left hand side matrix of coefficients

-  The right hand side matrix of coefficients

-  The equation rates

-  The option list

Error checking is performed at each step in the scanner and the parser.
For each syntax error the exact line and input file, along with an
appropriate error message are produced. Some other errors like mass
balance, and equation duplicates, are tested at the end of this phase.

.. _species-reordering:

==================
Species reordering
==================

When parsing the input files, the species list is updated as soon as a
new species is encountered in a chemical equation. Therefore the
ordering of the species is the order in which they appear in the
equation description section. This is not a useful order for subsequent
operations. The species have to be first sorted such that all variable
species and all fixed species are put together. Then if a sparsity
structure of the Jacobian is required, it might be better to reorder the
species in such a way that the factorization of the Jacobian will
preserve the sparsity. This reordering is done using a Markovitz type
algorithm.

.. _expression-trees:

============================
Expression trees computation
============================

This is the core of the preprocessor. This module generates the
production/destruction functions, the Jacobian and all the data
structure nedeed by these functions. It builds a language-independent
structure of each function and statement in the target source file.
Instead of using an intermediate format for this as some other compilers
do, KPP generates the intermediate format for just one statement at a
time. The vast majority of the statements in the target source file are
assignments. The expression tree for each assignment is incrementally
built by scanning the coefficient matrices and the rate constant vector.
At the end, these expression trees are simplified. Similar approaches are
applied to function declaration and prototypes, data declaration and
initialization.

.. _code-generation:

===============
Code generation
===============

There are basically two modules, each dealing with the syntax
particularities of the target language. For example, the C module
includes a function that generates a valid C assignment when given an
expression tree. Similarly there are functions for data declaration,
initializations, comments, function prototypes, etc. Each of these
functions produce the code into an output buffer. A language-specific
routine reads from this buffer and splits the statements into lines to
improve readability of the generated code.
