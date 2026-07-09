.. _filter-passive-spc:

==================================================================
Filtering out passive species with an absolute tolerance threshold
==================================================================

In some external models that use KPP (such as GEOS-Chem) it is common
to attach "passive" species to reactions.  This is usually done to
optain production or loss rates for diagnostic purposes.

We have created a simple passive species example C-I test (see the
chemcial mechanism description file
:file:`ci-tests/F90_ros_passivespc/F90_ros_passivespc.kpp`).  In this
example we add 3 passive species (:literal:`PASV1`, :literal:`PASV2`,
:literal:`PASV3`) are to the :literal:`small_strato` mechanism:

.. code-block:: console

   #EQUATIONS                         { SMALL STRATOSPHERIC MECHANISM           }
   <R1>  O2   + hv = 2O               : (2.643E-10) * SUN*SUN*SUN;
   <R2>  O    + O2 = O3               : (8.018E-17);
   <R3>  O3   + hv = O   + O2         : (6.120E-04) * SUN;
   <R4>  O    + O3 = 2O2              : (1.576E-15);
   <R5>  O3   + hv = O1D + O2         : (1.070E-03) * SUN*SUN;
   <R6>  O1D  + M  = O   + M          : (7.110E-11);
   <R7>  O1D  + O3 = 2O2      + PASV1 : (1.200E-10);
   <R8>  NO   + O3 = NO2 + O2 + PASV2 : (6.062E-15);
   <R9>  NO2  + O  = NO  + O2         : (1.069E-11);
   <R10> NO2  + hv = NO  + O  + PASV3 : (1.289E-02) * SUN;

Passive species, however, should not be allowed to influence certain
computations (such as the Rosenbrock error norm), as this may lead to
numerical instability.  In the case of the Rosenbrock error norm, we
have found that results of key species (such as O3 and NO) can vary
depending on the number of passive species that are added to the
mechanism.

In KPP 3.5.0 we have introduced the capability to filter out passive
species.  This feature is currently only available when generating
Fortran90 source code with KPP.  You may do this as follows:

First, give each passive species an extremely high absolute tolerance
(:code:`ATOL`) value:

.. code-block:: Fortran

   USE ROOT_Parameters

   ...

   ATOL(ind_PASV1) = 1.0d25
   ATOL(ind_PASV2) = 1.0d25
   ATOL(ind_PASV3) = 1.0d25

Next, pass the optional argument :code:`PassiveSpc_ATOL_Threshold` to
the :ref:`Initialize <Initialize>` routine (this needs to be done
only once at the model initialization phase):

.. code-block:: Fortran

   CALL Initialize( PassiveSpc_ATOL_Threshold = 1.0d25 )

Any species having a value of :code:`ATOL` greater than or equal to
:code:`PassiveSpc_ATOL_Threshold` will be denoted as a passive
species. Subroutine :ref:`Initialize <Initialize>` will then generate
the number and list of non-passive species in the following variables,
which are stored in the :ref:`Global` module:

.. code-block:: Fortran

   ! NonPassiveSpc_Count - Number of non-passive species in mechanism
     INTEGER :: NonPassiveSpc_Count
   ! NonPassiveSpc_Indices - Indices of non-passive species in mechanism
     INTEGER :: NonPassiveSpc_Indices(NVAR)

You may then use these variables to exclude passive species from
certain computations, such as:

.. code-block:: Fortran

   USE ROOT_Global

   ...

   INTEGER :: NP, I

   ! Loop over the number of non-passive species
   DO NP = 1, NonPassiveSpc_Count

      ! Only use non-passive species in computation
      I   = NonPassiveSpc_Indices(NP)
      Err = ATOL(I) * RTOL(I) * YMAX

      ... etc ...

   ENDDO
