.. _input-for-kpp:

##############
Input overview
##############

KPP basically handles two types of input files: **Kinetic description
files** and **auxiliary files**. Kinetic description files are in KPP
syntax and described in the following sections. Auxiliary files are
described in the section entitled
:ref:`auxiliary-files-and-the-substitution-preprocessor`.

KPP kinetic description files specify the chemical equations, the
initial values of each of the species involved, the integration
parameters, and many other options. The KPP preprocessor parses the
kinetic description files and generates several output files. Files
that are written in KPP syntax have one of the suffixes :file:`.kpp`,
:file:`.spc`, :file:`.eqn`, or :file:`def`.

The following general rules define the structure of a kinetic
description file:

-  A KPP program is composed of :ref:`kpp-sections`,
   :ref:`kpp-commands`, and :ref:`inlined-code`. Their syntax is
   presented in :ref:`bnf-description`.

-  Comments are either enclosed between the curly braces ":code:`{`"
   and ":code:`}`", or written in a line starting with two slashes and
   a space "// ".

-  Any name given by the user to denote an atom or a species is
   restricted to be less than 32 character in length and can only
   contain letters, numbers, or the underscore character. The first
   character cannot be a number. All names are case insensitive.

The kinetic description files contain a detailed specification of the
chemical model, information about the integration method and the desired
type of results. KPP accepts only one of these files as input, but using
the :ref:`include-cmd` command, code from separate files can be
combined. The include files can be nested up to 10 levels. KPP will
parse these files as if they were a single big file. By carefully
splitting the chemical description, KPP can be configured for a broad
range of users. In this way the users can have direct access to that
part of the model that they are interested in, and all the other details
can be hidden inside several include files. Often, the atom definitions
(:file:`atoms.kpp`) are included first, then species definitions
(:file:`*.spc`), and finally the equations of the chemical mechanism
(:file:`*.eqn`).
