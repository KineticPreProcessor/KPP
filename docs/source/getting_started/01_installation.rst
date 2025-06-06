.. _install:

############
Installation
############

This section can be skipped if KPP is already installed on your system.

========================
Download KPP from Github
========================

Clone the KPP source code from the `KPP Github repository
<https://github.com/KineticPreProcessor/KPP>`_:

.. code-block:: console

   $ cd $HOME
   $ git clone https://github.com/KineticPreProcessor/KPP.git

This will create a directory named KPP in your home directory.

========================================
Define the KPP_HOME environment variable
========================================

Define the :envvar:`$KPP_HOME` environment variable to point to the
complete path where KPP is installed.  Also, add the path of the KPP
executable to the :envvar:`$PATH` environment variable.

If you are using the Unix C-shell (aka :program:`csh`), add
add these statements to your :file:`$HOME/.cshrc` file:

.. code-block:: csh

   setenv KPP_HOME $HOME/kpp
   setenv PATH ${PATH}:$KPP_HOME/bin

and then apply the settings with:

.. code-block:: console

   $ source $HOME/.cshrc

If, on the other hand, you are using the Unix :program:`bash` shell,
add these statements to your :file:`$HOME/.bashrc` file:

.. code-block:: bash

   export KPP_HOME=$HOME/kpp
   export PATH=$PATH:$KPP_HOME/bin

and then apply the settings with:

.. code-block:: console

   $ source $HOME/.bashrc

Now if you type:

.. code-block:: console

   $ which kpp

the path to the executable file (:file:`kpp`) will be displayed. This
path should match the path specified by :file:`$KPP_HOME/bin/kpp`.

.. _test-for-dependencies:

=====================================================
Test if KPP dependencies are installed on your system
=====================================================

KPP depends on several other Unix packages.  Before using KPP for the
first time, test if these are installed on your system.  If any of
these packages are missing, you can install them with your
system's package manager (e.g. :program:`apt` for Ubuntu,
:program:`yum` for Fedora, :program:`homebrew` for MacOS, etc.), or
with `Spack <https://spack.readthedocs.io>`_.

gcc
---

.. important::

   You might have to follow some :ref:`additional configuration
   and installation steps <additional-steps-macosx>` regarding
   :program:`gcc` on MacOS X systems.

KPP uses the `GNU Compiler Collection <https://gcc.gnu.org/>`_ (aka
:program:`gcc`) by default. A version of :program:`gcc` comes
pre-installed with most Linux or MacOS systems. To test if
:program:`gcc` is installed on your system, type:

.. code-block :: console

   $ gcc --version

This will display the version information, such as:

.. code-block:: console

   gcc (GCC) 11.2.0
   Copyright (C) 2021 Free Software Foundation, Inc.
   This is free software; see the source for copying conditions.  There is NO
   warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

sed
---
The :program:`sed` utility is used to search for and replace text
in files.  To test if :program:`sed` has been installed, type:

.. code-block:: console

   $ which sed

This will print the path to :program:`sed` on your system.

bison
-----

The :program:`bison` utility parses the chemical mechanism file into a
computer-readable syntax.  To test :program:`bison` is installed, type:

.. code-block:: console

   $ which bison

This will print the path to :program:`bison` on your system.

.. _flex-dep:

flex
----

.. important::

   You might have to follow some :ref:`additional configuration
   and installation steps <additional-steps-macosx>` regarding
   :program:`flex` on MacOS X systems.

The :program:`flex` (the Fast Lexical Analyzer) creates a scanner that
can recognize the syntax generated by :program:`bison`.  To test if
:program:`flex` is installed, type:

.. code-block:: console

   $ which flex

This will print the path to :program:`flex` on your system.

You will also need to specify the path to the :program:`flex` library
files (:file:`libfl.so` or :file:`libfl.a`) in order to :ref:`build
the KPP executable <build-kpp-exec>`.  This can be done with the
:program:`find` command:

.. code-block:: console

   $ find /usr/ -name "*libfl*" -print

This will generate a list of file paths such as shown below.  Look for
the text :file:`libfl.`:

.. code-block:: console

   /usr/include/libflashrom.h
   /usr/lib/gthumb/extensions/libflicker.so
   /usr/lib/gthumb/extensions/libflicker_utils.so
   /usr/lib/libflashrom.so.1.0.0
   /usr/lib/libfl.so                # <---- This is the flex library file
   # ... etc ...

Once you have located the directory where flex library file
resides (which in this example is :file:`/usr/lib`), use it to define
the :envvar:`KPP_FLEX_LIB_DIR`  environment variable in your
:file:`.bashrc` (or :file:`.bash_aliases` file if you have one):

.. code-block:: bash

   export KPP_FLEX_LIB_DIR=/usr/lib
   export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${KPP_FLEX_LIB_DIR}:"

Then apply the changes with:

.. code-block:: console

   . ~/.bashrc

KPP will use the path specified by :envvar:`KPP_FLEX_LIB_DIR` during
the compilation sequence (described in the next section).

.. _build-kpp-exec:

========================
Build the KPP executable
========================

Change to the KPP/src directory:

.. code-block:: console

   $ cd $KPP_HOME/src

To clean a previously-built KPP installation, delete the KPP object
files and all the examples with:

.. code-block:: console

   $ make clean

To delete a previoulsy-built KPP executable as well, type:

.. code-block:: console

   $ make distclean

KPP will use :program:`gcc` as the default compiler.  If you would
like to use a different compiler (e.g. :program:`icc`), then edit
:file:`src/Makefile.defs` to add your compiler name.

Create the KPP executable with:

.. code-block:: console

   $ make

You should see output similar to:

.. code-block:: console

   gcc -g -Wall -Wno-unused-function -I/usr/include -c code.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c code_c.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c code_f77.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c code_f90.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c code_matlab.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c debug.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c gen.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c kpp.c
   flex -olex.yy.c scan.l
   bison -d -o y.tab.c scan.y
   gcc -g -Wall -Wno-unused-function -I/usr/include -c lex.yy.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c scanner.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c scanutil.c
   gcc -g -Wall -Wno-unused-function -I/usr/include -c y.tab.c
   gcc -g -Wall -Wno-unused-function code.o code_c.o
       code_f77.o code_f90.o code_matlab.o debug.o gen.o kpp.o
       lex.yy.o scanner.o scanutil.o y.tab.o -L/usr/lib -lfl -o kpp

This will create the executable file :file:`$KPP_HOME/bin/kpp`.

.. _additional-steps-macosx:

==============================
Instructions for MacOS X users
==============================

When installing KPP on a MacOS X system, some additional configuration
and installation steps may be necessary.

.. _force-macos-to-recognize-gcc-compiler:

Force MacOS to recognize the gcc compiler
-----------------------------------------

On MacOS X, if you type:

.. code-block:: console

   $ gcc --version

you will probably see output similar to:

.. code-block:: console

   Apple clang version 13.1.6 (clang-1316.0.21.2.5)
   Target: x86_64-apple-darwin21.5.0
   Thread model: posix
   InstalledDir: /Library/Developer/CommandLineTools/usr/bin

This is because MacOS X installs :program:`clang` as :program:`gcc`.
To force MacOS X to recognize the :program:`gcc` compiler, follow
these steps:

#. Use the :program:`homebrew` package manager to install
   :program:`gcc`:

   .. code-block:: console

      $ brew install gcc

#. Type this command:

   .. code-block:: console

      $ ls /usr/local/Cellar/gcc/*/bin/ | grep gcc

   You should see output such as:

   .. code-block:: console

      gcc-11*
      gcc-ar-11*
      gcc-nm-11*
      gcc-ranlib-11*
      # ... etc ...

   This output indicates :program:`gcc` major version 11 has been
   installed, and that the gcc executable is called :code:`gcc-11`.
   (Your version may differ.)

#. Add  the following code block to your :file:`.bashrc` file (or to your
   :file:`.bash_aliases` file if you have one).  This will define
   aliases that will override :program:`clang` with :program:`gcc`.

   .. code-block:: bash

      #============================================================================
      # Compiler settings (MacOS)
      #
      # NOTE: MacOSX installs Clang as /usr/bin/gcc, so we have to manually
      # force reference to gcc-11, g++-11, and gfortran-11, which HomeBrew
      # installs to /usr/local/bin.  (bmy, 10/28/21)
      #============================================================================
      alias gcc=gcc-11
      alias g++=g++-11
      alias gfortran=gfortran-11
      export CC=gcc
      export CXX=g++-11
      export FC=gfortran-11
      export F77=gfortran-11

   Then apply the changes with:

   .. code-block:: console

      $ . ~/.bashrc

#. To check if your shell now recognizes the :program:`gcc` compiler, type:

   .. code-block:: console

      $ gcc --version

   You should see output similar to:

   .. code-block:: console

      gcc-11 (Homebrew GCC 11.3.0_1) 11.3.0
      Copyright (C) 2021 Free Software Foundation, Inc.
      This is free software; see the source for copying conditions.  There is NO
      warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

   This now indicates that your compiler is :program:`gcc` and not
   :program:`clang`.

.. _install-flex-with-homebrew:

Install flex with homebrew
--------------------------

If your MacOS X computer does not have the :program:`flex` library
installed, then you can install it with :program:`homebrew`:

.. code-block:: console

   $ brew install flex

Unlike Linux pacakge managers, which would install the :program:`flex`
library files in the path :file:`/usr/lib/`,
:program:`homebrew` will install it to a path such as
:file:`/usr/local/Cellar/flex/X.Y.Z/lib/`.

To find the version of :program:`flex` that has been installed by
:program:`homebrew`, type:

.. code-block:: console

   $ ls /usr/local/Cellar/flex

and you will get a listing such as:

.. code-block:: console

   2.6.4_2

This indicates that the version of :program:`flex` on your system is
:code:`2.6.4_2` (the :code:`_2` denotes the number of bug-fix updates
since version :code:`2.6.4` was released).

The :program:`flex` library files (:file:`libfl.so` or
:file:`libfl.a`) will be found in :file:`lib/` subfolder.  In this
example, the path will be:

.. code-block:: console

   /usr/local/Cellar/flex/2.6.4_2/lib

Knowing this, you can now define the :envvar:`KPP_FLEX_LIB_DIR`
environment variable :ref:`as described above <flex-dep>`:

.. code-block:: bash

   export FLEX_LIB_DIR=/usr/local/Cellar/flex/2.6.4_2/lib

.. _macosx-limited-stack:

Request maximum stack memory
----------------------------

MacOS X has a hard limit of 65332 bytes for stack memory.  This is
much less memory than what is available on GNU/Linux operating systems
such as Ubuntu, Fedora, etc.

To make sure you are using the maximum amount of stack memory on MacOS
X add this command to your :file:`.bashrc` file:

.. code-block:: bash

   ulimit -s 65532

and then apply the change with:

.. code-block:: console

   $ . ~/.bashrc

This stack memory limit means that KPP will not be able to parse
mechanisms with more than about 2000 equations and 1000 species.
Because of this, we have added an :code:`#ifdef` block to KPP header
file :file:`src/gdata.h` to define the :code:`MAX_EQN` and
:code:`MAX_SPECIES` parameters accordingly:

.. code-block:: C

   #ifdef MACOS
   #define MAX_EQN        2000     // Max number of equations (MacOS only)
   #define MAX_SPECIES    1000     // Max number of species   (MacOS only)
   #else
   #define MAX_EQN       11000     // Max number of equations
   #define MAX_SPECIES    6000     // Max number of species
   #endif

If you find that KPP will not parse your mechanism, you can increase
:code:`MAX_EQN` and decrease :code:`MAX_SPECIES` (or vice-versa) as
needed, and then :ref:`rebuild the KPP executable <build-kpp-exec>`.

.. _macosx-case-insensitive:

Know that MacOS X is case-insenstive
-------------------------------------

If you have two files with identical names except for case
(e.g. :file:`integrator.F90` and :file:`integrator.f90`) then MacOS X
will not be able to tell them apart.  Because of this, you may
encounter an error if you try to commit such files into Git, etc.
