###########################################
Telling GEOS-Chem to use a custom mechanism
###########################################

---------------------
Configuring GEOS-Chem
---------------------

GEOS-Chem will always use the :file:`fullchem` mechanism by default.  To
configure GEOS-Chem to use the :file:`custom` mechanism instead
of fullchem, navigate to your build directory and type:

.. code-block:: console

   $ cmake ../CodeDir -DCUSTOMMECH=y

.. note:: GEOS-Chem Classic run directories have a subdirectory named
	  :file:`build` in which you can configure and build GEOS-Chem.

	  For more information about the GEOS-Chem and GCHP
	  configuration process, please see `GEOS-Chem manual
	  <http://wiki.geos-chem.org/Getting_Started_with_GEOS-Chem>`__
	  and `gchp.readthedocs.io <https://gchp.readthedocs.io>`__.
  
You should see output such as this written to the screen:

.. code-block:: none

  -- General settings:
     * CUSTOMMECH:  **ON** OFF

This confirms that GEOS-Chem will use the :file:`custom` mechanism.

-------------------
Compiling GEOS-Chem
-------------------

Once you have configured GEOS-Chem to use the :file:`custom` mechanism,
you may build the exectuable.  Type:

.. code-block:: console

   $ make -j
   $ make -j install

The executable file (:file:`gcclassic` or :file:`gchp`, depending on which
mode of GEOS-Chem that you are using) will be placed in the run directory.
