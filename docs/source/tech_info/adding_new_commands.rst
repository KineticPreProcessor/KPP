.. |br| raw:: html

   <br />

.. _adding-new-commands:

#######################
Adding new KPP commands
#######################

To add a new KPP command, the source code has to be edited at several
locations. A short summary is presented here, using :code:`NEWCMD` as an
example:

#. Add the new command to several files in the :file:`src/` directory:

   - :file:`scan.h`: add :code:`void CmdNEWCMD( char *cmd );`
   - :file:`scan.l`: add :code:`{ "NEWCMD", PRM_STATE, NEWCMD },`
   - :file:`scanner.c`: add :code:`void CmdNEWCMD( char *cmd )`
   - :file:`scan.y`:

     - Add :code:`%token NEWCMD`
     - Add :code:`| NEWCMD PARAMETER`
     - Add :code:`{ CmdNEWCMD( $2 ); }` |br| |br|


#. Add :ref:`ci-tests`:

   - Create a new directory :file:`ci-tests/ros_newcmd/ros_newcmd.kpp`
   - Add new :ref:`ci-tests` to the :file:`ci-tests` directory and
     update the scripts in the :file:`.ci-pipelines` directory. |br| |br|

#. Add documentation to user manual :file:`docs/source/*/*.rst`:

   - Add to Table :ref:`table-cmd-defaults`
   - Add a new subsection to :ref:`kpp-commands`
   - Add to the Table :ref:`bnf-description` |br| |br|

#. Add command to :file:`site-lisp/kpp.el` (Emacs configuration file)
