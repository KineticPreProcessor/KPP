.. _rosenbrock-methods:

##################
Rosenbrock methods
##################

**Integrator file:** :file:`int/rosenbrock.f90`

An :math:`s`-stage Rosenbrock method (cf. Section IV.7 in
:cite:t:`Hairer_and_Wanner_1991`) computes the next-step solution by the
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

=====
ROS-2
=====
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

=====
ROS-3
=====

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

=====
ROS-4
=====

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

=======
RODAS-3
=======

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

.. _rosenbrock-rodas-3-1:

=========
RODAS-3.1
=========

- Reference: :cite:t:`Long_et_al._2026`
- Stages (:math:`s`): 4
- Funcion calls: 3
- Order: 3(2)
- Stability properties: Stiffly-accurate
- Method Coefficients:

.. math::

   \begin{aligned}
   a_{2,1} = 0                  & \qquad & a_{3,1} = 0.646601929740551  &\\
   a_{3,2} = 0.409567801987914  & \qquad & a_{4,1} = 0.646601929740551  &\\
   a_{4,2} = 0.409567801987914  & \qquad & a_{4,3} = 1                  &\\
   c_{2,1} = 4.198495621784201  & \qquad & c_{3,1} = 3.711590161613010  &\\
   c_{3,2} = -1.787771994729384 & \qquad & c_{4,1} = 4.458898153216104  &\\
   c_{4,2} = -2.024095448516552 & \qquad & c_{4,3} = -2.626700600119396 &\\
   m_1 = 0.646601929740551      & \qquad & m_2 = 0.409567801987914      &\\
   m_3 = 1                      & \qquad & m_4 = 1                      &\\
   e_1 = 0                      & \qquad & e_2 = 0                      &\\
   e_3 = 0                      & \qquad & e_4 = 1                      &\\
   \alpha_1 = 0                 & \qquad & \alpha_2 = 0                 &\\
   \alpha_3 = 1                 & \qquad & \alpha_4 = 1                 &\\
   \gamma_1 = 0.515             & \qquad & \gamma_2 = 1.628546001287715 &\\
   \gamma_3 = 0                 & \qquad & \gamma_4 = 0
   \end{aligned}

.. _rosenbrock-rodas-4:

=======
RODAS-4
=======

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

===============================
Rosenbrock tangent linear model
===============================

**Integrator file:** :file:`int/rosenbrock_tlm.f90`

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
analysis) for each of the Rosenbrock methods (:ref:`rosenbrock-ros-2`,
:ref:`rosenbrock-ros-3`, :ref:`rosenbrock-ros-4`,
:ref:`rosenbrock-rodas-3`, and :ref:`rosenbrock-rodas-4`). The
implementations distinguish between sensitivities with respect to
initial values and sensitivities with respect to parameters for
efficiency.

.. _rosenbrock-adjoint:

=================================
Rosenbrock discrete adjoint model
=================================

**Integrator file:** :file:`int/rosenbrock_adj.f90`

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
:ref:`rosenbrock-rodas-3`, :ref:`rosenbrock-rodas-4`).

.. _rosenbrock-autoreduce:

========================================
Rosenbrock with mechanism auto-reduction
========================================

**Integrator file:** :file:`int/rosenbrock_autoreduce.f90`

Mechanism auto-reduction (described in :cite:t:`Lin_et_al._2022`) expands
previous work by :cite:t:`Santillana_et_al._2010` and
:cite:t:`Shen_et_al._2020` to a computationally efficient implementation
in KPP, avoiding memory re-allocation, re-compile of the code, and
on-the-fly mechanism reduction based on dynamically determined
production and loss rate thresholds.

We define a threshold :math:`\delta` which can be fixed (as in
:cite:t:`Santillana_et_al._2010`) or determined by the production and
loss rates of a "target species" scaled by a factor

.. math::

   \begin{aligned}
   \delta = max(P_{target}, L_{target}) * \alpha_{target}.
   \end{aligned}

For each species :math:`i`, the species is partitioned as "slow" iff.

.. math::

   \begin{aligned}
   max(P_i, L_i) < \delta
   \end{aligned}

if the species is partitioned as "slow", it is solved explicitly
(decoupled from the rest of the mechanism) using a first-order
approximation. Otherwise, "fast" species are retained in the implicit
Rosenbrock solver.

.. _rosenbrock-h211b-qssa:

===================================
Rosenbrock with H211b time stepping
===================================

**Integrator file:** :file:`int/rosenbrock_h211b_qssa.f90`

H211b time stepping according to :cite:t:`Soederlind_2003`, as
implemented by :cite:t:`Dreger_2025`.
