.. _System Requirements:

###################
System requirements
###################

:program:`KPP-for-GEOS-Chem` requires the following packages:

#. A C-language compiler (such as :program:`gcc`, from the `GNU Compiler
   Collection <https://gcc.gnu.org/>`__)
#. The :program:`flex` lexical parser library.  This is often installed on
   many computer systems by default.  It can also be easily installed
   with a package manager such as :program:`spack`.
#. Python (2.7.5 or greater).  This is required in order
   to post-process source code created by KPP-for-GEOS-Chem to include
   modifications for the OH reactivity diagnostic.

Depending on your setup, you might have to load these packages with the
:command:`module load` or :command:`spack load` commands. Ask your sysadmin
for more information.
