[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html) [![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/KPP/blob/GC_updates/LICENSE.txt)

# KPP

This is the repository for the Kinetic PreProcessor (KPP).

We are in the process of updating this repository.  We will post
updated documentation to our ReadTheDocs site:
https://kpp.readthedocs.io in the near future.


## Versions
KPP - symbolic chemistry Kinetics PreProcessor:

- Current: Version 2.2.3 (in the **main** branch)
- In development: Version 3.0.0

## License

  KPP is distributed under GPL, the [general public license](http://www.gnu.org/copyleft/gpl.html)
  
    (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
    (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
        with contributions from:
        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany

## Getting Started
   Read [user's manual](https://github.com/KineticPreProcessor/KPP/blob/master/doc/kpp_UserManual.pdf)

## Installation

1. Make sure that FLEX (public domain lexical analizer) is installed
   on your machine. Type "flex --version" to test this.

2. Note down the exact path name where the FLEX library is installed. The
   library is called:
	libfl.a or libfl.sh 

3. Define the KPP_HOME environment variable to point to the complete 
   path location of KPP. If, for example, KPP is installed in $HOME/kpp:

   - with C shell (or tcsh) edit the file .cshrc (or .tcshrc) in your
     home directory and add:
	setenv KPP_HOME $HOME/kpp
	set path=( $path $HOME/kpp/bin )
     Execute 'source .cshrc' (or 'source .tcshrc') to make sure these
     changes are in effect.

   - with bash shell edit the file .bashrc in your home directory and add:
	export KPP_HOME=$HOME/kpp
	export PATH=$PATH:$HOME/kpp/bin

3. In KPP_HOME directory edit: 
	Makefile.defs 
   and follow the instructions included to specify the compiler, 
   the location of the FLEX library, etc.
 
4. In KPP_HOME directory build the sources using:
	make

## Cleanup 

1. Delete the KPP object files with:
	make clean

2. Delete the whole distribution (including the KPP binaries) with:
	make distclean

## Troubleshooting
If you have any problems please send the detailed error report and the machine
environment to:

	geos-chem-support@g.harvard.edu

## Branches of development

### main

This is the default branch, which corresponds to the unmodified KPP version 2.2.3_01 as published by A. Sandu et al.  For more information, [see the KPP website](https://people.cs.vt.edu/~asandu/Software/Kpp/).

* NOTE: In keeping with Github best practices, we have renamed "master" to "main".

### dev

This is the development branch for KPP 3.0.0.

### GC_updates

This branch includes modifications to KPP that were made specifically for the GEOS-Chem model.  For more information, please see the [README.md file in the GC_updates branch](https://github.com/KineticPreProcessor/KPP/blob/GC_updates/README.md).

### F77

This branch incorporates fixes from Josue Bock for generating chemical mechanism solver files in Fortran-77 format.

### mistra

This branch contains modifications to KPP that were made specifically for the mistra project.

## Documentation:

* [KPP user manual](https://github.com/KineticPreProcessor/KPP/blob/main/doc/kpp_UserManual.pdf)
* [KPP references](https://people.cs.vt.edu/~asandu/Software/Kpp/docsforkpp.htm)
