.. _other-methods:

###############################
Forward differentiation methods
###############################

.. _other-methods-feuler:

FEULER
------

**Integrator file:** :file:`int/feuler.f90`

Forward Euler is an explicit integration method for non-stiff problems.
FEULER computes
:math:`y^{n+1}` as

.. math::

   \begin{aligned}
   y^{n+1} = y^n + hf\left(t^{n},y^{n}\right)
   \end{aligned}
