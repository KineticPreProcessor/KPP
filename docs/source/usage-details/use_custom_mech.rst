###########################################
Telling GEOS-Chem to use a custom mechanism
###########################################

:program:`GEOS-Chem` will always use the default mechanism (which is named
:file:`fullchem` or :file:`Standard` depending on which version you
have).  To tell GEOS-Chem to use the :file:`custom` mechanism, follow
these steps.

----------------------------------
For GEOS-Chem versions using CMake
----------------------------------

If your :program:`GEOS-Chem` version uses CMake, thenen navigate to your build directory.

.. note:: GEOS-Chem Classic run directories have a subdirectory named
	  :file:`build` in which you can configure and build
  	  GEOS-Chem.  If you don't have a build directory, you can add one to your run directory.
		
	  For more information about the GEOS-Chem and GCHP
	  configuration process, please see `GEOS-Chem manual
	  <http://wiki.geos-chem.org/Getting_Started_with_GEOS-Chem>`__
	  and `gchp.readthedocs.io <https://gchp.readthedocs.io>`__.

^^^^^^^^^^^^^^^^^^^^^
Configuring GEOS-Chem
^^^^^^^^^^^^^^^^^^^^^

From the build directory, type:

.. code-block:: console

   $ cmake ../CodeDir -DCUSTOMMECH=y -DRUNDIR=..

NOTE: In some later :program:`GEOS-Chem 12` versions, the proper command is:

.. code-block:: console

   $ cmake ../CodeDir -DCHEM=custom

Please consult the `GEOS-Chem manual <https://wiki.geos-chem.org/Getting_Started_with_GEOS-Chem>`__ for more information.

You should see output similar to this written to the screen:

.. code-block:: none

  -- General settings:
     * CUSTOMMECH:  **ON** OFF

This confirms that the custom mechanism has been selected.

^^^^^^^^^^^^^^^^^^^
Compiling GEOS-Chem
^^^^^^^^^^^^^^^^^^^

Once you have configured :program:`GEOS-Chem` to use the :file:`custom` mechanism,
you may build the exectuable.  Type:

.. code-block:: console

   $ make -j
   $ make -j install

The executable file (:file:`gcclassic` or :file:`gchp`, depending on which
mode of GEOS-Chem that you are using) will be placed in the run directory.

-------------------------------------
For GEOS-Chem versions using GNU Make
-------------------------------------

If you are using an older version of :program:`GEOS-Chem` that is only compatible
with GNU Make, then use this command to compile the executable:

.. code-block:: console

   $ make -j CHEM=Custom ...etc other build options...

Please consult the `GEOS-Chem manual <https://wiki.geos-chem.org/Getting_Started_with_GEOS-Chem>`__ for more information.
