.. _numerical-methods:

#################
Numerical methods
#################

The KPP numerical library contains a set of numerical integrators
selected to be very efficient in the low to medium accuracy regime
(relative errors :math:`\sim 10^{-2} \dots 10^{-5}`). In addition, the
KPP numerical integrators preserve the linear invariants (i.e., mass) of
the chemical system.

KPP implements several Rosenbrock methods: ROS–2
:cite:`1999:Verwer`, ROS–3 :cite:`1997:Sandu_et_al_2`,
RODAS–3 :cite:`1997:Sandu_et_al_2`, ROS–4
:cite:`1991:Hairer_and_Wanner`, and RODAS–4
:cite:`1991:Hairer_and_Wanner`. For each of them KPP implements the
tangent linear model (direct decoupled  ensitivity) and the adjoint
models. The implementations distinguish between sensitivities with
respect to initial values and sensitivities with respect to parameters
for efficiency.

Note that KPP produces the building blocks for the simulation and also
for the sensitivity calculations. It also provides application
programming templates. Some minimal programming may be required from the
users in order to construct their own application from the KPP building
blocks.

In the following sections we introduce the numerical methods implemented
in KPP. The symbols used in the formulas are explained in
:ref:`table-symbols`.

.. _table-symbols:

.. table:: Table 20. Symbols used in numerical methods
   :align: center

   +----------------------------------+----------------------------------+
   | Symbol                           | Description                      |
   +==================================+==================================+
   | :math:`s`                        | Number of stages                 |
   +----------------------------------+----------------------------------+
   | :math:`t^n`                      | Discrete time moment             |
   +----------------------------------+----------------------------------+
   | :math:`h`                        | Time step :math:`h=t^{n+1}-t^n`  |
   +----------------------------------+----------------------------------+
   | :math:`y^n`                      | Numerical solution               |
   |                                  | (concentration) at :math:`t^n`   |
   +----------------------------------+----------------------------------+
   | :math:`\delta y^n`               | tangent linear solution at       |
   |                                  | :math:`t^n`                      |
   +----------------------------------+----------------------------------+
   | :math:`\lambda^n`                | Adjoint numerical solution at    |
   |                                  | :math:`t^n`                      |
   +----------------------------------+----------------------------------+
   | :math:`f(\cdot,\cdot)`           | The ODE derivative function:     |
   |                                  | :math:`y'=f(t,y)`                |
   +----------------------------------+----------------------------------+
   | :math:`f_t(\cdot,\cdot)`         | Partial time derivative          |
   |                                  | :math:`f_t(                      |
   |                                  | t,y)=\partial f(t,y)/\partial t` |
   +----------------------------------+----------------------------------+
   | :math:`J(\cdot,\cdot)`           | The Jacobian                     |
   |                                  | :math:`J(                        |
   |                                  | t,y)=\partial f(t,y)/\partial y` |
   +----------------------------------+----------------------------------+
   | :math:`J_t(\cdot,\cdot)`         | Partial time derivative of       |
   |                                  | Jacobian                         |
   |                                  | :math:`J_t(                      |
   |                                  | t,y)=\partial J(t,y)/\partial t` |
   +----------------------------------+----------------------------------+
   | :math:`A`                        | The system matrix                |
   +----------------------------------+----------------------------------+
   | :math:`H(\cdot,\cdot)`           | The Hessian                      |
   |                                  | :math:`H(t,y)                    |
   |                                  | =\partial^2 f(t,y)/\partial y^2` |
   +----------------------------------+----------------------------------+
   | :math:`T_i`                      | Internal stage time moment for   |
   |                                  | Runge-Kutta and Rosenbrock       |
   |                                  | methods                          |
   +----------------------------------+----------------------------------+
   | :math:`Y_i`                      | Internal stage solution for      |
   |                                  | Runge-Kutta and Rosenbrock       |
   |                                  | methods                          |
   +----------------------------------+----------------------------------+
   | :math:`k_i`, :math:`\ell_i`,     | Internal stage vectors for       |
   | :math:`u_i`, :math:`v_i`         | Runge-Kutta and Rosenbrock       |
   |                                  | methods, their tangent linear    |
   |                                  | and adjoint models               |
   +----------------------------------+----------------------------------+
   | :math:`\alpha_i`,                | Method coefficients              |
   | :math:`\alpha_{ij}`,             |                                  |
   | :math:`a_{ij}`, :math:`b_i`,     |                                  |
   | :math:`c_i`, :math:`c_{ij}`,     |                                  |
   | :math:`e_i`, :math:`m_i`         |                                  |
   +----------------------------------+----------------------------------+

.. _rosenbrock-methods:

==================
Rosenbrock Methods
==================

An :math:`s`-stage Rosenbrock method (cf. Section IV.7 in
:cite:`1991:Hairer_and_Wanner`) computes the next-step solution by the
formulas

.. _alt-rosenbrock:

.. math::

   \begin{aligned}
   y^{n+1} &=& y^n + \sum_{i=1}^s m_i k_i~,
   \quad {\rm Err}^{n+1} = \sum_{i=1}^s e_i k_i\\
   \nonumber
   T_i &=& t^n + \alpha_i h~, \quad
   Y_i =y^n + \sum_{j=1}^{i-1} a_{ij} k_j~,\\
   \nonumber
   A &=& \left[ \frac{1}{h \gamma} - J^T(t^n,y^n) \right]\\
   \nonumber
   A \cdot k_i &=&  f\left( \, T_i,
   \, Y_i \,\right) + \sum_{j=1}^{i-1} \frac{c_{ij}}{h} k_j + h \gamma_i
   f_t\left(t^n,y^n\right)~.
   \end{aligned}

where :math:`s` is the number of stages, :math:`\alpha_i = \sum_j
\alpha_{ij}` and  :math:`\gamma_i = \sum_j \gamma_{ij}`. The formula
coefficients (:math:`a_{ij}` and :math:`\gamma_{ij}`) give the order
of consistency and the stability properties. :math:`A` is the system
matrix (in the linear systems to be solved during implicit
integration, or in the Newton’s method used to solve the nonlinear
systems). It is the scaled identity matrix minus the Jacobian.

The coefficients of the methods implemented in KPP are shown below:

.. _rosenbrock-ros-2:

ROS-2
-----
- Stages (:math:`s`): 2
- Funcion calls: 2
- Order: 2(1)
- Stability properties: L-stable
- Method Coefficients:

.. math::

   \begin{aligned}
   \gamma = 1 + 1/sqrt{2} & \qquad & a_{2,1} = 1/\gamma & \qquad & c_{2,1} = -2/\gamma  &\\
   m_1 = 3/(2\gamma)      & \qquad & m_2 = 1/(2\gamma)  & \qquad & e_1 = 1/(2\gamma)    &\\
   e_2 = 1/(2\gamma)      & \qquad & \alpha_1 = 0       & \qquad & \alpha_2 = 1         &\\
   \gamma_1 = \gamma      & \qquad & \gamma_2 = -\gamma                                 &\\
   \end{aligned}

.. _rosenbrock-ros-3:

ROS-3
-----
- Stages (:math:`s`): 3
- Funcion calls: 2
- Order: 3(2)
- Stability properties: L-stable
- Method Coefficients:

.. math::

   \begin{aligned}
   a_{2,1} = 1       & \qquad & a_{3,1} = 1       & \qquad & a_{3,2} = 0       &\\
   c_{2,1} = -1.015  & \qquad & c_{3,1} = 4.075   & \qquad & c_{3,2} = 9.207   &\\
   m_1 = 1           & \qquad & m_2 = 6.169       & \qquad & m_3 = -0.427      &\\
   e_1 = 0.5         & \qquad & e_2 = -2.908      & \qquad & e_3 = 0.223       &\\
   alpha_1 = 0       & \qquad & \alpha_2 = 0.436  & \qquad & \alpha_3 = 0.436  &\\
   \gamma_1 = 0.436  & \qquad & \gamma_2 = 0.243  & \qquad & \gamma_3 =  2.185 &\\
   \end{aligned}

.. _rosenbrock-ros-4:

ROS-4
-----
- Stages (:math:`s`): 4
- Funcion calls: 3
- Order: 4(3)
- Stability properties: L-stable
- Method Coefficients:

.. math::

   \begin{aligned}
   a_{2,1} = 2        & \qquad & a_{3,1} = 1.868     & \qquad & a_{3,2} = 0.234     &\\
   a_{4,1} = a_{3,1}  & \qquad & a_{4,2} = a_{3,2}   & \qquad & a_{4,3} = 0         &\\
   c_{2,1} = -7.137   & \qquad & c_{3,1} = 2.581     & \qquad & c_{3,2} = 0.652     &\\
   c_{4,1} = -2.137   & \qquad & c_{4,2} = -0.321    & \qquad & c_{4,3} = -0.695    &\\
   m_1 = 2.256        & \qquad & m_2 = 0.287         & \qquad & m_3 = 0.435         &\\
   m_4 = 1.094        & \qquad & e_1 = -0.282        & \qquad & e_2 = -0.073        &\\
   e_3 = -0.108       & \qquad & e_4 = -1.093        & \qquad & \alpha_1 = 0        &\\
   \alpha_2 = 1.146   & \qquad & \alpha_3 = 0.655    & \qquad & \alpha_4 = \alpha_3 &\\
   \gamma_1 = 0.573   & \qquad & \gamma_2 = -1.769   & \qquad & \gamma_3 = 0.759    &\\
   \gamma_4 = -0.104
   \end{aligned}

.. _rosenbrock-rodas-3:

RODAS-3
-------
- Stages (:math:`s`): 4
- Funcion calls: 3
- Order: 3(2)
- Stability properties: Stiffly-accurate
- Method Coefficients:

.. math::

   \begin{aligned}
   a_{2,1} = 0    & \qquad & a_{3,1} = 2     & \qquad & a_{3,2} = 0    &\\
   a_{4,1} = 2    & \qquad & a_{4,2} = 0     & \qquad & a_{4,3} = 1    &\\
   c_{2,1} = 4    & \qquad & c_{3,1} = 1     & \qquad & c_{3,2} = -1   &\\
   c_{4,1} = 1    & \qquad & c_{4,2} = -1    & \qquad & c_{4,3} = -8/3 &\\
   m_1 = 2        & \qquad & m_2 = 0         & \qquad & m_3 = 1        &\\
   m_4 = 1        & \qquad & e_1 = 0         & \qquad & e_2 = 0        &\\
   e_3 = 0        & \qquad & e_4 = 1         & \qquad & \alpha_1 = 0   &\\
   \alpha_2 = 0   & \qquad & \alpha_3 = 1    & \qquad & \alpha_4 = 1   &\\
   \gamma_1 = 0.5 & \qquad & \gamma_2 = 1.5  & \qquad & \gamma_3 = 0   &\\
   \gamma_4 = 0
   \end{aligned}

.. _rosenbrock-rodas-4:

RODAS-4
-------
- Stages (:math:`s`): 6
- Funcion calls: 5
- Order: 4(3)
- Stability properties: Stiffly-accurate
- Method Coefficients:

.. math::

  \begin{aligned}
  \alpha_1 = 0      & \qquad & \alpha_2 = 0.386  & \qquad & \alpha_3 = 0.210  &\\
  \alpha_4 = 0.630  & \qquad & \alpha_5 = 1      & \qquad & \alpha_6 = 1      &\\
  \gamma_1 = 0.25   & \qquad & \gamma_2 = -0.104 & \qquad & \gamma_3 = 0.104  &\\
  \gamma_4 = -0.036 & \qquad & \gamma_5 = 0      & \qquad & \gamma_6 = 0      &\\
  a_{2,1} = 1.544   & \qquad & a_{3,1} = 0.946   & \qquad & a_{3,2} = 0.255   &\\
  a_{4,1} = 3.314   & \qquad & a_{4,2} = 2.896   & \qquad & a_{4,3} = 0.998   &\\
  a_{5,1} = 1.221   & \qquad & a_{5,2} = 6.019   & \qquad & a_{5,3} = 12.537  &\\
  a_{5,4} = -0.687  & \qquad & a_{6,1} = a_{5,1} & \qquad & a_{6,2} = a_{5,2} &\\
  a_{6,3} = a_{5,3} & \qquad & a_{6,4} = a_{5,4} & \qquad & a_{6,5} = 1       &\\
  c_{2,1} = -5.668  & \qquad & c_{3,1} = -2.430  & \qquad & c_{3,2} = -0.206  &\\
  c_{4,1} = -0.107  & \qquad & c_{4,2} = -9.594  & \qquad & c_{4,3} = -20.47  &\\
  c_{5,1} = 7.496   & \qquad & c_{5,2} = -0.124  & \qquad & c_{5,3} = -34     &\\
  c_{5,4} = 11.708  & \qquad & c_{6,1} = 8.083   & \qquad & c_{6,2} = -7.981  &\\
  c_{6,3} = -31.521 & \qquad & c_{6,4} = 16.319  & \qquad & c_{6,5} = -6.058  &\\
  m_1 = a_{5,1}     & \qquad &  m_2 = a_{5,2}    & \qquad & m_3 = a_{5,3}     &\\
  m_4 = a_{5,4}     & \qquad & m_5 = 1           & \qquad & m_6 = 1           &\\
  e_1 = 0           & \qquad & e_2 = 0           & \qquad & e_3 = 0           &\\
  e_4 = 0           & \qquad & e_5 = 0            & \qquad & e_6 = 1          &\\
  \end{aligned}

.. _rosenbrock-tlm:

Rosenbrock Tangent Linear Model
--------------------------------

The Tangent Linear method is combined with the sensitivity
equations. One step of the method reads:

.. math::

   \begin{aligned}
   %y^{n+1} &=& y^n + \sum_{i=1}^s m_i k_i, \qquad
   \delta y^{n+1} &=& \delta y^n + \sum_{i=1}^s m_i \ell_i\\
   \nonumber
   T_i &=& t^n + \alpha_i h~, %\quad Y_i =y^n + \sum_{j=1}^{i-1} a_{ij} k_j~,
   \quad \delta Y_i = \delta y^n + \sum_{j=1}^{i-1} a_{ij} \ell_j\\
   %A &=& \left[ \frac{1}{h \gamma} - J^T(t^n,y^n) \right]\\
   %\nonumber
   %A \cdot k_i &=&
   %           f\left( \, T_i,\, Y_i \,\right)
   %           + \sum_{j=1}^{i-1} \frac{c_{ij}}{h} k_j
   %          + h \gamma_i f_t\left(t^n,y^n\right)~,\\
   \nonumber
   A \cdot \ell_i &=&
           J\left( \, T_i,\, Y_i \,\right)
                 \cdot \delta Y_i
                 + \sum_{j=1}^{i-1} \frac{c_{ij}}{h} \ell_j\\
   \nonumber
   && +
   \left( H( t^n, y^n )\times  k_i \right) \cdot \delta y^n
      + h \gamma_i J_t\left(t^n,y^n\right) \cdot \delta y^n\end{aligned}

The method requires a single `n \times n` LU decomposition per
step to obtain both the concentrations and the sensitivities.

KPP contains tangent linear models (for direct decoupled sensitivity
analysis) for each of the Rosenbrock methods (ROS–2, ROS–3, ROS–4,
RODAS–3, and RODAS–4). The implementations distinguish between
sensitivities with respect to initial values and sensitivities with
respect to parameters for efficiency.

.. _rosenbrock-adjoint:

Rosenbrock Discrete Adjoint Model
---------------------------------

To obtain the adjoint we first differentiate the method with respect to
:math:`y_n`. Here :math:`J` denotes the Jacobian and :math:`H` the
Hessian of the derivative function :math:`f`. The discrete adjoint of
the (non-autonomous) Rosenbrock method is

.. math::

   \begin{aligned}
   \label{Ros_disc_adj}
   %A &=& \left[ \frac{1}{h \gamma} - J^T(t^n,y^n) \right]\\
   %\nonumber
   A \cdot u_i
   &=& m_i \lambda^{n+1} + \sum_{j=i+1}^s \left( a_{ji} v_j + \frac{c_{ji}}{h}
   u_j \right)~,\\
   \nonumber
   v_i &=& J^T(T_i,Y_i)\cdot u_i~, \quad i = s,s-1,\cdots,1~,\\
   \nonumber
   \lambda^n &=& \lambda^{n+1} + \sum_{i=1}^s \left( H(t^n,y^n) \times
   k_i\right)^T
   \cdot u_i\\
   \nonumber
   && + h J^T_t(t^n,y^n) \cdot \sum_{i=1}^s \gamma_i u_i+  \sum_{i=1}^s v_i\end{aligned}

KPP contains adjoint models (for direct decoupled sensitivity analysis)
for each of the Rosenbrock methods (:ref:`rosenbrock-ros-2`,
:ref:`rosenbrock-ros-3`, :ref:`rosenbrock-ros-4`,
:ref:`rosenbrock-rodas3`, :ref:`rosenbrock-rodas-4`).

.. _rk-methods:

============================
Runge-Kutta (aka RK) methods
============================

A general :math:`s`-stage Runge-Kutta method is defined as (see
Section II.1 of :cite:`1993:Hairer_Norsett_and_Wanner`)

.. math::

   \begin{aligned}
   \label{eqn:RungeKutta}
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

.. _rk-method-comparison:

3-stage Runge-Kutta
-------------------

**Integrator file:** :file:`int/runge_kutta.f90`

Fully implicit 3-stage Runge-Kutta methods.  Several variants are available:

- RADAU-2A: order 5
- RADAU-1A: order 5
- Lobatto-3C: order 4
- Gauss: order 6

RADAU5
------
**Integrator files:** :file:`int/kpp_radau5.f90`

This Runge-Kutta method of order 5 based on RADAU-IIA quadrature
is stiffly accurate. The KPP implementation follows the original
implementation of :cite:`1991:Hairer_and_Wanner`, Section IV.10. While
RADAU5 is relatively expensive (when compared to the Rosenbrock
methods), it is more robust and is useful to obtain accurate reference
solutions.

SDIRK
-----
**Integrator file:** :file:`int/sdirk.f90`,

SDIRK is an L-stable, singly-diagonally-implicit Runge-Kutta method. The
implementation is based on :cite:`1991:Hairer_and_Wanner`. Several
variants are available:

  - Sdirk 2a, 2b: 2 stages, order 2
  - Sdirk 3a: 3 stages, order 2
  - Sdirk 4a, 4b: 5 stages, order 4

SDIRK4
------
**Integrator file:** :file:`int/kpp_sdirk4.f90`

SDIRK4 is an L-stable, singly-diagonally-implicit Runge-Kutta method
of order 4. The implementation is based on :cite:`1991:Hairer_and_Wanner`.

SEULEX
------
**Integrator file:** :file:`int/kpp_seulex.f90`

SEULEX is a variable  order stiff extrapolation code able to produce
highly accurate solutions. The KPP implementation is based on the
implementation of :cite:`1991:Hairer_and_Wanner`.

.. _rk-tlm:

RK Tangent Linear Model
------------------------

The tangent linear method associated with the Runge-Kutta method is

.. math::

   \begin{aligned}
   \label{eqn:RK-TLM}
   %y^{n+1} &=& y^n + h \sum_{i=1}^s b_i k_i~,\\
   \delta y^{n+1} &=& \delta y^n + h \sum_{i=1}^s b_i \ell_i~,\\
   \nonumber
   %Y_i &=& y^n + h \sum_{j=1}^{s} a_{ij} k_j~,\\
   \delta Y_i& =& \delta y^n + h \sum_{j=1}^{s} a_{ij} \ell_j~,\\
   \nonumber
   %k_i &=& f\left( \, T_i, \, Y_i \,\right)~,\\
   \ell_i &=& J\left(T_i, \, Y_i \right) \cdot \delta Y_i ~.\end{aligned}

The system (`[eqn:RK-TLM] <#eqn:RK-TLM>`__) is linear and does not
require an iterative procedure. However, even for a SDIRK method
(:math:`a_{ij}=0` for :math:`i>j` and :math:`a_{ii}=\gamma`) each stage
requires the LU factorization of a different matrix.

.. _rk-adj:

RK Discrete Adjoint Model
-------------------------

The first order Runge-Kutta adjoint is

.. math::

   \begin{aligned}
   \label{RK-adj}
   u_i &=& h \, J^T(T_i,Y_i)\cdot
   \left( b_i \lambda^{n+1} + \sum_{j=1}^s a_{ji} u_j \right)\\ %\quad i = 1 \cdots s\\
   \nonumber
   \lambda^{n} &=& \lambda^{n+1} +\sum_{j=1}^s u_j~.\end{aligned}

For :math:`b_i \ne 0` the Runge-Kutta adjoint can be rewritten as
another Runge-Kutta method:

.. math::

   \begin{aligned}
   \label{RK-adj-2}
   u_i &=& h \, J^T(T_i,Y_i)\cdot
   \left( \lambda^{n+1} + \sum_{j=1}^s \frac{b_j \,
   a_{ji}}{b_i} u_j \right)\\ %~, \quad i = 1 \cdots s\\
   \nonumber
   \lambda^{n} &=& \lambda^{n+1} +\sum_{j=1}^s b_j \, u_j~.\end{aligned}

.. _back-diff:

=================================
Backward Differentiation Formulas
=================================

Backward differentiation formulas (BDF) are linear multistep methods
with excellent stability properties for the integration of chemical
systems (cf. :cite:`1991:Hairer_and_Wanner`, Section V.1). The
:math:`k`-step BDF method reads

.. math::

   \sum_{i=0}^k \alpha_i y^{n-i} = h_n \beta\; f\left(t^{n},y^{n}\right)
   \label{BDF}

where the coefficients :math:`\alpha_i` and :math:`\beta` are chosen
such that the method has order of consistency :math:`k`.

The KPP library contains two off-the-shelf, highly popular
implementations of BDF methods, described in the following sections:

LSODE
-----
**Integrator file:** :file:`int/kpp_lsode.f90`

LSODE, the Livermore ODE solver :cite:`1993:LSODE`, implements backward
differentiation formula (BDF) methods for stiff problems.  LSODE has
been translated to Fortran90 for the incorporation into the KPP library.

VODE
-----

**Integrator file:** :file:`int/kpp_dvode.f90`

VODE :cite:`1989:VODE` uses another formulation of backward
differentiation formulas. The version of VODE present in the KPP
library uses directly the KPP sparse linear algebra routines.

.. _integrator-input-output:

=============================
Integrator inputs and outputs
=============================

In order to offer more control over the integrator, the KPP-generated
subroutine provides the :ref:`optional-int-inputs`. Each of them is an
array of 20 elements that allow the fine-tuning of the integrator.

Similarly, to obtain more information about the integration, the
subroutine provides the :ref:`optional-int-outputs`, which are also
also arrays of 20 elements.

.. _optional-int-inputs:

Optional integrator input parameters
------------------------------------

Optional integer (:code:`ICNTRL_U`) and real
(:code:`RCNTRL_U`) input parameters subroutine :code:`INTEGRATE`.
Setting any elements to zero will activate their default values. Array
elements not listed here are either not used or integrator-specific options.
Details can be found in the comment lines of the individual integrator files
:code:`$KPP_HOME/int/*.f90`.

.. option:: ICNTRL_U(1)

   :code:`= 1`: :math:`F = F(y)`, i.e. independent of t (autonomous)

   :code:`= 0`: :math:`F = F(t,y)`, i.e. depends on t (non-autonomous)

   This option is only available for some of the integrators.

.. option:: ICNTRL_U(2)

   The absolute (:code:`ATOL`) and relative (:code:`RTOL`) tolerances
   can be expressed by either a scalar or individually for each
   species in a vector:

   :code:`= 0` : :code:`NVAR` -dimensional vector

   :code:`= 1` : scalar

.. option:: ICNTRL_U(3)

   Selection of a specific method (only available for some of the
   integrators).

.. option:: ICNTRL_U(4)

   Maximum number of integration steps.

.. option:: ICNTRL_U(5)

   Maximum number of Newton iterations (only available for some of the
   integrators).

.. option:: ICNTRL_U(6)

   Starting values of Newton iterations (only avaialble for some of
   the integrators).

   :code:`= 0` : Interpolated

   :code:`= 1` : Zero

.. option:: ICNTRL_U(15)

   This determines which :code:`Update_*` subroutines are called
   within the integrator.

   :code:`= -1` : Do not call any :code:`Update_*` subroutines

   :code:`=  0` :  Use the integrator-specific default values

   :code:`>  1` : A number between 1 and 7, derived by adding up bits
   with values 4, 2, and 1.  The first digit (4) activates
   :code:`Update_SUN`.  The second digit (2) activates
   :code:`Update_PHOTO`.  The third digit (1) activates
   :code:`Update_RCONST`.    |

   For example :code:`ICNTRL(15)=6)` (4+2) will activate the calls to
   :code:`Update_SUN` and :code:`Update_PHOTO`, but not to
   :code:`Update_RCONST`.

.. option:: RCNTRL_U(1)

   :code:`Hmin`, the lower bound of the integration step size. It is
   not recommended to change the default value of zero.

.. option:: RCNTRL_U(2)

   :code:`Hmax`, the upper bound of the integration step size.

.. option:: RCNTRL_U(3)

   :code:`Hstart`, the starting value of the integration step size.

.. option:: RCNTRL_U(4)

   :code:`FacMin`, lower bound on step decrease factor.

.. option:: RCNTRL_U(5)

   :code:`FacMax`, upper bound on step increase factor.

.. option:: RCNTRL_U(6)

   :code:`FacRej`, step decrease factor after multiple rejections.

.. option:: RCNTRL_U(7)

   :code:`FacSafe`, the factor by which the new step is slightly
   smaller than the predicted value.

.. option:: RCNTRL_U(8)

   :code:`ThetaMin`. If the Newton convergence rate is smaller than
   ThetaMin, the Jacobian is not recomputed (only available for some
   of the integrators).

.. option:: RCNTRL_U(9)

   :code:`NewtonTol`, the stopping criterion for Newton’s method (only
   available for some of the integrators).

.. option:: RCNTRL_U(10)

   :code:`Qmin` (only available for some of the integrators).

.. option:: RCNTRL_U(11)

   :code:`Qmax`. If :code:`Qmin < Hnew/Hold < Qmax`, then the step
   size is kept constant and the LU factorization is reused (only
   available for some of the integrators).

.. _optional-int-outputs:

Optional integrator output parameters
-------------------------------------

Optional integer (:code:`ISTATUS_U`) and real (:code:`RSTATUS_U`)
output parameters of subroutine :code:`INTEGRATE`.  Array elements not
listed here are either not used or are integrator-specific options.
Details can be found in the comment lines of the individual integrator files
:code:`$KPP_HOME/int/*.f90`.

.. option:: ISTATUS_U(1)

    Number of function calls.

.. option:: ISTATUS_U(2)`

   Number of Jacobian calls.

.. option:: ISTATUS_U(3)

   Number of steps.

.. option:: ISTATUS_U(4)

   Number of accepted steps.

.. option:: ISTATUS_U(5)

   Number of rejected steps (except at very beginning).

.. option:: ISTATUS_U(6)

   Number of LU decompositions.

.. option:: ISTATUS_U(7)

   Number of forward/backward substitutions.

.. option:: ISTATUS_U(8)

   Number of singular matrix decompositions.

.. option:: RSTATUS_U(1)

   :code:`Texit`, the time corresponding to the computed :math:`Y`
   upon return.

.. option:: RSTATUS_U(2)

  :code:`Hexit`: the last accepted step before exit.

.. option:: RSTATUS_U(3)

   :code:`Hnew`: The last predicted step (not yet taken.  For multiple
   restarts, use :code:`Hnew` as :code:`Hstart` in the subsequent run.
