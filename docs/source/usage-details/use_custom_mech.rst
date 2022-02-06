###########################################
Telling GEOS-Chem to use a custom mechanism
###########################################

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
you may build the executable.  Type:

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