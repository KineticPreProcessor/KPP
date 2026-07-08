.. _Makefile:

============
The Makefile
============

KPP produces a Makefile that allows for an easy compilation of all
KPP-generated source files. The file name is :file:`Makefile_ROOT`. The
Makefile assumes that the selected driver contains the main program.
However, if no driver was selected (i.e. :command:`#DRIVER none`), it is
necessary to add the name of the main program file manually to the
Makefile.
