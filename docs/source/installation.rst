.. _install:

############
Installation
############

This section can be skipped if KPP is already installed on your system.
If you work under Linux, you can probably use the precompiled executable
file that is in the directory of the distribution. Then you only have to
define the :envvar:`$KPP_HOME` environment variable.

1. Define the :envvar:`$KPP_HOME` environment variable to point to the
   complete path where KPP is installed. Also, add the path of the KPP
   executable to the :envvar:`$PATH` environment variable. If, for example,
   KPP is installed in :envvar:`$HOME/kpp`, under the C shell you have to edit
   the file :envvar:`$HOME/.cshrc` and add:

.. code-block:: console

   $ setenv KPP_HOME $HOME/kpp
   $ setenv PATH ${PATH}:$KPP_HOME/bin

If you use the bash shell, edit and add:

.. code-block:: bash

   $ export KPP_HOME=$HOME/kpp
   $ export PATH=$PATH:$KPP_HOME/bin

After editing or , start a new shell to make sure these changes are
in effect.

2. Make sure that :program:`sed` is installed on your machine. To test
   this, type:

.. code-block:: console

   $ which sed

3. Make sure that :program:`bison` is installed on your machine. To
   test this, type:

.. code-block:: console

   $ which bison

4. Make sure that the fast lexical analyzer generator is installed on
   your machine. To test this, type:

.. code-block:: console

   $ which flex

   Enter the path where the flex library (:program:`libfl.a` or
   program:`libfl.so`or ) is located into `Makefile.defs`, e.g.

.. code-block:: make

   FLEX_LIB_DIR=/usr/lib

5. Change to the KPP directory:

.. code-block:: console

   $ cd $KPP_HOME

6. To clean the KPP installation, delete the KPP object files and all
   the examples with:

.. code-block:: console

   $ make clean

To delete the KPP executable as well, type:

.. code-block:: console

   $ make distclean

7. If necessary, edit and enter the name of your C compiler. The default
   setting is :program:`gcc`.

#. Create the kpp executable with:

.. code-block:: console

   $ make
