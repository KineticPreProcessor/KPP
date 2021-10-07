.. _creating_fortran_files:

############################################################
Creating Fortran-90 chemical mechanism modules for GEOS-Chem
############################################################

--------------------------------------------------------
Navigate to the KPP folder in your GEOS-Chem source code
--------------------------------------------------------

At this point you can now use :program:`KPP-for-GEOS-Chem` to generate Fortran-90 source code
files that will solve the chemical mechanism in an efficient manner.

Navigate to this folder in your :program:`GEOS-Chem` source code:

-  GEOS-Chem 12.9.3 and prior versions: :file:`KPP`
-  GEOS-Chem 13.0.0 and later versions: :file:`src/GEOS-Chem/KPP`
-  GCHP 13.0.0 and later versions:
   :file:`src/GCHP_GridComp/GEOSChem_GridComp/geos-chem/KPP`

Here you will find two sub-folders: :file:`fullchem` and
:file:`custom`, and a script named :file:`build_mechanism.sh`.

The :file:`custom` folder contains sample chemical mechanism
specification files (:file:`custom.eqn` and :file:`gckpp.kpp`) which
have been copied from fullchem. You can edit these files to define your
own custom mechanism (see subsequent sections for detailed
instructions).

.. note:: The :file:`fullchem` folder contains chemical mechanism
          specification files FlexChem-KPP-generated source code for
	  the default GEOS-Chem mechanism (named fullchem).  You should
          leave these files untouched. This will allow you to revert
	  to the fullchem mechanism if need be.

---------------------------------
Run the build_mechanism.sh script
---------------------------------
	  
Once you are satisfied with your custom mechanism specification you may
now use KPP-for-GEOS-Chem to build the source code files for GEOS-Chem.

Return to the KPP folder containing ``build_mechanism.sh`` and then type:

.. code-block:: console

   $ ./build_mechanism.sh custom

You will see output similar to this:

.. code-block:: none

  This is KPP-2.3.1_gc.
  KPP is parsing the equation file.
  KPP is computing Jacobian sparsity structure.
  KPP is starting the code generation.
  KPP is initializing the code generation.
  KPP is generating the monitor data:
    - gckpp_Monitor
  KPP is generating the utility data:
    - gckpp_Util
  KPP is generating the global declarations:
    - gckpp_Main
  KPP is generating the ODE function:
    - gckpp_Function
  KPP is generating the ODE Jacobian:
    - gckpp_Jacobian
    - gckpp_JacobianSP
  KPP is generating the linear algebra routines:
    - gckpp_LinearAlgebra
  KPP is generating the utility functions:
    - gckpp_Util
  KPP is generating the rate laws:
    - gckpp_Rates
  KPP is generating the parameters:
    - gckpp_Parameters
  KPP is generating the global data:
    - gckpp_Global
  KPP is generating the driver from none.f90:
    - gckpp_Main
  KPP is starting the code post-processing.
  
  KPP has succesfully created the model "gckpp".
  
  Reactivity consists of 172 reactions
  Written to gckpp_Util.F90

If this process is successful, the :file:`custom` folder should now be
populated with several :file:`.F90` source code files:

.. code-block:: none

  CMakeLists.txt*      gckpp_Initialize.F90  gckpp_LinearAlgebra.F90  gckpp_Precision.F90
  custom.eqn           gckpp_Integrator.F90  gckpp.map                gckpp_Rates.F90
  gckpp_Function.F90   gckpp_Jacobian.F90    gckpp_Model.F90          gckpp_Util.F90
  gckpp_Global.F90     gckpp_JacobianSP.F90  gckpp_Monitor.F90        Makefile_gckpp
  gckpp_HetRates.F90@  gckpp.kpp             gckpp_Parameters.F90

These files contain optimized Fortran-90 instructions for solving the chemical
mechanism that you have specified.
