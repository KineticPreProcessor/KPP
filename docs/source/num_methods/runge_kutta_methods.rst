.. _rk-methods:

############################
Runge-Kutta (aka RK) methods
############################

A general :math:`s`-stage Runge-Kutta method is defined as (see
Section II.1 of :cite:t:`Hairer_Norsett_and_Wanner_1987`)

.. math::

   \begin{aligned}
   y^{n+1} &=& y^n + h \sum_{i=1}^s b_i k_i~,\\
   \nonumber
   T_i &=& t^n + c_i h~, \quad
   Y_i = y^n + h \sum_{j=1}^{s} a_{ij} k_j~,\\
   \nonumber
   k_i &=& f\left( \, T_i, \, Y_i \,\right)~,\end{aligned}

where the coefficients :math:`a_{ij}`, :math:`b_i` and :math:`c_i` are
prescribed for the desired accuracy and stability properties. The stage
derivative values :math:`k_i` are defined implicitly, and require
solving a (set of) nonlinear system(s). Newton-type methods solve
coupled linear systems of dimension (at most) :math:`n \times s`.

The Runge-Kutta methods implemented in KPP are summarized below:

.. _rk-methods-3stage:

===================
3-stage Runge-Kutta
===================

**Integrator file:** :file:`int/runge_kutta.f90`

Fully implicit 3-stage Runge-Kutta methods.  Several variants are available:

- RADAU-2A: order 5
- RADAU-1A: order 5
- Lobatto-3C: order 4
- Gauss: order 6

.. _rk-methods-radau5:

======
RADAU5
======

**Integrator file:** :file:`int/radau5.f90`

This Runge-Kutta method of order 5 based on RADAU-IIA quadrature
is stiffly accurate. The KPP implementation follows the original
implementation of :cite:t:`Hairer_and_Wanner_1991`, Section IV.10. While
RADAU5 is relatively expensive (when compared to the Rosenbrock
methods), it is more robust and is useful to obtain accurate reference
solutions.

.. _rk-methods-sdirk:

=====
SDIRK
=====

**Integrator file:** :file:`int/sdirk.f90`,

SDIRK is an L-stable, singly-diagonally-implicit Runge-Kutta method. The
implementation is based on :cite:t:`Hairer_and_Wanner_1991`. Several
variants are available:

  - Sdirk 2a, 2b: 2 stages, order 2
  - Sdirk 3a: 3 stages, order 2
  - Sdirk 4a, 4b: 5 stages, order 4

.. _rk-methods-sdirk4:

======
SDIRK4
======

**Integrator file:** :file:`int/sdirk4.f90`

SDIRK4 is an L-stable, singly-diagonally-implicit Runge-Kutta method
of order 4. The implementation is based on :cite:t:`Hairer_and_Wanner_1991`.

.. _rk-methods-seulex:

======
SEULEX
======

**Integrator file:** :file:`int/seulex.f90`

SEULEX is a variable  order stiff extrapolation code able to produce
highly accurate solutions. The KPP implementation is based on the
implementation of :cite:t:`Hairer_and_Wanner_1991`.

.. _rk-methods-tlm:

=======================
RK tangent linear model
=======================

The tangent linear method associated with the Runge-Kutta method is

.. math::

   \begin{aligned}
   %y^{n+1} &=& y^n + h \sum_{i=1}^s b_i k_i~,\\
   \delta y^{n+1} &=& \delta y^n + h \sum_{i=1}^s b_i \ell_i~,\\
   \nonumber
   %Y_i &=& y^n + h \sum_{j=1}^{s} a_{ij} k_j~,\\
   \delta Y_i& =& \delta y^n + h \sum_{j=1}^{s} a_{ij} \ell_j~,\\
   \nonumber
   %k_i &=& f\left( \, T_i, \, Y_i \,\right)~,\\
   \ell_i &=& J\left(T_i, \, Y_i \right) \cdot \delta Y_i ~.\end{aligned}

The system is linear and does not require an iterative
procedure. However, even for a SDIRK method (:math:`a_{ij}=0` for
:math:`i>j` and :math:`a_{ii}=\gamma`) each stage requires the LU
factorization of a different matrix.

.. _rk-methods-adj:

=========================
RK discrete adjoint model
=========================

The first order Runge-Kutta adjoint is

.. math::

   \begin{aligned}
   u_i &=& h \, J^T(T_i,Y_i)\cdot
   \left( b_i \lambda^{n+1} + \sum_{j=1}^s a_{ji} u_j \right)\\ %\quad i = 1 \cdots s\\
   \nonumber
   \lambda^{n} &=& \lambda^{n+1} +\sum_{j=1}^s u_j~.\end{aligned}

For :math:`b_i \ne 0` the Runge-Kutta adjoint can be rewritten as
another Runge-Kutta method:

.. math::

   \begin{aligned}
   u_i &=& h \, J^T(T_i,Y_i)\cdot
   \left( \lambda^{n+1} + \sum_{j=1}^s \frac{b_j \,
   a_{ji}}{b_i} u_j \right)\\ %~, \quad i = 1 \cdots s\\
   \nonumber
   \lambda^{n} &=& \lambda^{n+1} +\sum_{j=1}^s b_j \, u_j~.\end{aligned}
