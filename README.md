[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/KPP/blob/GC_updates/LICENSE.txt) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4552707.svg)](https://doi.org/10.5281/zenodo.4552707) [![ReadTheDocs](https://img.shields.io/readthedocs/kpp?label=ReadTheDocs)](https://kpp.readthedocs.io/en/latest/)

# README for KPP-for-GEOS-Chem
### Description

The GC_updates branch of the KPP repository tracks modifications to KPP that are specific to GEOS-Chem.

### Documentation

Complete documentation for KPP-for-GEOS-Chem may now be found https://kpp.readthedocs.io.
## Versions
  KPP - symbolic chemistry Kinetics PreProcessor, [Version 2.2.3](http://www.cs.vt.edu/~asandu/Software/KPP)

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
