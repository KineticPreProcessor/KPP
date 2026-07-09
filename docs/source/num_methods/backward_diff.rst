.. _back-diff:

################################
Backward differentiation methods
################################

Backward differentiation formulas (BDF) are linear multistep methods
with excellent stability properties for the integration of chemical
systems (cf. :cite:t:`Hairer_and_Wanner_1991`, Section V.1). The
:math:`k`-step BDF method reads

.. math::

   \sum_{i=0}^k \alpha_i y^{n-i} = h_n \beta\; f\left(t^{n},y^{n}\right)
   \label{BDF}

where the coefficients :math:`\alpha_i` and :math:`\beta` are chosen
such that the method has order of consistency :math:`k`.

The KPP library contains two off-the-shelf, highly popular
implementations of BDF methods, described in the following sections:

.. _back-diff-lsode:

=====
LSODE
=====

**Integrator file:** :file:`int/lsode.f90`

LSODE, the Livermore ODE solver
(:cite:t:`Radhakrishnan_and_Hindmarsh_1993`), implements backward
differentiation formula (BDF) methods for stiff problems.  LSODE has
been translated to Fortran90 for the incorporation into the KPP library.

.. attention::

   We have discovered that the current implementation of the LSODE
   integrator is not thread-safe for `OpenMP parallelization
   <https://www.openmp.org/>`_.  When LSODE is called from within
   an OpenMP parallel loop, the integration will fail because key
   internal variables in LSODE will be overwritten by concurrent
   threads.

.. _back-diff-vode:

====
VODE
====

**Integrator file:** :file:`int/dvode.f90`

VODE (:cite:t:`Brown_Byrne_and_Hindmarsh_1989`) uses another formulation
of backward differentiation formulas. The version of VODE present in
the KPP library uses directly the KPP sparse linear algebra routines.

.. _back-diff-beuler:

======
BEULER
======

**Integrator file:** :file:`int/sdirk.f90`

Backward Euler integration method.  To request this method, make sure
you select

.. code-block:: console

   #INTEGRATOR sdirk

in your definition file, and then set :code:`ICNTRL(3) = 6`.
