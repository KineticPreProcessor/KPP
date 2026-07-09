.. _icntrl-rcntrl:

######################################################
Integrator options (:code:`ICNTRL` and :code:`RCNTRL`)
######################################################

In order to offer more control over the integrator, KPP provides the
arrays :code:`ICNTRL` (integer) and :code:`RCNTRL` (real). Each of them
is an array of 20 elements that allow the fine-tuning of the integrator.
All integrators (except for :code:`tau_leap` and :code:`gillespie`) use
:code:`ICNTRL` and :code:`RCNTRL`. Details can be found in the comment
lines of the individual integrator files in :code:`$KPP_HOME/int/`.

======
ICNTRL
======

.. _table-icntrl:

.. table:: Summary of ICNTRL usage in the f90 integrators.
           Here, Y = used, and s = solver-specific usage.
   :align: center

   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | ICNTRL                 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
   +========================+===+===+===+===+===+===+===+===+===+====+====+====+====+====+====+====+====+====+
   | beuler                 |   | Y | Y | Y | Y | Y | s |   |   |    |    |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | dvode                  |   |   |   |   |   |   |   |   |   |    |    |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | exponential            |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | feuler                 |   |   |   |   |   |   |   |   |   |    |    |    |    |    | Y  | Y  | Y  |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | gillespie              |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | lsode                  |   | Y |   | Y |   |   |   |   |   | s  |    |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | radau5                 |   | Y |   | Y | Y | Y |   |   |   |    | Y  |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | rosenbrock_adj         | Y | Y | Y | Y |   | s | s | s |   |    |    |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | rosenbrock             | Y | Y | Y | Y |   |   |   |   |   |    |    |    |    |    | Y  | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | rosenbrock_tlm         | Y | Y | Y | Y |   |   |   |   |   |    |    | s  |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | rosenbrock_autoreduce  | Y | Y | Y | Y |   |   |   |   |   |    |    | s  | s  | s  | Y  | Y  |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | rosenbrock_h211b_qssa  | Y | Y | Y | Y |   |   |   |   |   |    |    |    |    |    | Y  | Y  |    | s  |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | runge_kutta_adj        |   | Y | Y | Y | Y | s | s | s | s | s  | Y  |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | runge_kutta            |   | Y | Y | Y | Y | Y |   |   |   | s  | Y  |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | runge_kutta_tlm        |   | Y | Y |   | Y | Y | s |   | s | s  | Y  | s  |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | sdirk4                 |   | Y |   | Y |   |   |   |   |   |    |    |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | sdirk_adj              |   | Y | Y | Y | Y | Y | s | s |   |    |    |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | sdirk                  |   | Y | Y | Y | Y | Y |   |   |   |    |    |    |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | sdirk_tlm              |   | Y | Y | Y | Y | Y | s |   | s |    |    | s  |    |    | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | seulex                 | Y | Y |   | Y |   |   |   |   |   | s  | s  | s  | s  | s  | Y  |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+
   | tau_leap               |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+

ICNTRL(1)
---------

Specifies the time-dependence of the function :math:`F` to be integrated:

.. table::
   :align: left

   +------------+----------------------------------------+
   | ICNTRL(1)  | Time dependence of :math:`F`           |
   +============+========================================+
   | 0          | :math:`F = F(t,y)` aka non-autonomous  |
   +------------+----------------------------------------+
   | 1          | :math:`F = F(y)` aka autonomous        |
   +------------+----------------------------------------+

ICNTRL(2)
---------

Specifies the dimensionality of the absolute (:code:`ATOL`) and
relative (:code:`RTOL`) tolerances.  These can be expressed by
either a scalar or individually for each species in a vector.

.. table::
   :align: left

   +------------+-------------------------------------------+
   | ICNTRL(2)  | Dimensionality of ``ATOL`` and ``RTOL``   |
   +============+===========================================+
   | 0          | :code:`NVAR` -dimensional vector          |
   +------------+-------------------------------------------+
   | 1          | scalar                                    |
   +------------+-------------------------------------------+

ICNTRL(3)
---------

Selects a specific integration method.

.. table::
   :align: left

   +-----------------------+------------+-----------------------+
   | Integrator            | ICNTRL(3)  | Method selected       |
   +=======================+============+=======================+
   | rosenbrock            | 0 or 4     | Rodas3 (default)      |
   |                       +------------+-----------------------+
   | rosenbrock_adj        | 1          | Ros2                  |
   |                       +------------+-----------------------+
   | rosenbrock_autoreduce | 2          | Ros3                  |
   |                       +------------+-----------------------+
   | rosenbrock_tlm        | 3          | Ros4                  |
   |                       +------------+-----------------------+
   |                       | 5          | Rodas4                |
   |                       +------------+-----------------------+
   |                       | 6          | Rang                  |
   |                       +------------+-----------------------+
   |                       | 7          | Rodas3.1              |
   +-----------------------+------------+-----------------------+
   | runge_kutta           | 0 or 1     | Radau-2A (default)    |
   |                       +------------+-----------------------+
   | runge_kutta_adj       | 2          | Lobatto-3C            |
   |                       +------------+-----------------------+
   | runge_kutta_tlm       | 3          | Gauss                 |
   |                       +------------+-----------------------+
   |                       | 4          | Radau-1A              |
   |                       +------------+-----------------------+
   |                       | 5          | Lobatto-3A            |
   +-----------------------+------------+-----------------------+
   | sdirk                 | 0 or 1     | Sdirk-2A (default)    |
   |                       +------------+-----------------------+
   | sdirk_adj             | 2          | Sdirk-2B              |
   |                       +------------+-----------------------+
   | sdirk_tlm             | 3          | Sdirk-3A              |
   |                       +------------+-----------------------+
   |                       | 4          | Sdirk-4A              |
   |                       +------------+-----------------------+
   |                       | 5          | Sdirk-4B              |
   |                       +------------+-----------------------+
   |                       | 6          | Backward Euler        |
   +-----------------------+------------+-----------------------+


ICNTRL(4)
---------

Specifies the maximum number of integration steps.

ICNTRL(5)
---------

Specifies the maximum number of Newton iterations.

ICNTRL(6)
---------

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(6)  | Option selected                   |
   +=================+============+===================================+
   | rosenbrock_adj  | 0 thru 6   | Selection of a particular         |
   |                 |            | Rosenbrock method for the         |
   |                 |            | continuous adjoint integration    |
   |                 |            | (see ``ICNTRL(3)``)               |
   +-----------------+------------+-----------------------------------+
   | radau5          | 0          | Starting values for Newton        |
   |                 |            | iterations are interpolated       |
   |                 |            | (default)                         |
   | runge_kutta_adj +------------+-----------------------------------+
   |                 | 1          | Starting values for Newton        |
   | runge_kutta_tlm |            | iterations are zero               |
   |                 |            |                                   |
   | sdirk           |            |                                   |
   |                 |            |                                   |
   | sdirk_adj       |            |                                   |
   |                 |            |                                   |
   | sdirk_tlm       |            |                                   |
   +-----------------+------------+-----------------------------------+

ICNTRL(7)
---------

Selects options for adjoint integrators.

.. table:: Selection of adjoint algorithm
   :align: left

   +-----------------+------------+----------------------------------+
   | Integrator      | ICNTRL(7)  | Adjoint algorithm selected       |
   +=================+============+==================================+
   | rosenbrock_adj  | 0 or 2     | Discrete adjoint with method     |
   |                 |            | ``ICNTRL(3)`` (default)          |
   |                 +------------+----------------------------------+
   |                 | 1          | No adjoint                       |
   |                 +------------+----------------------------------+
   |                 | 3          | Fully adaptive continous adjoint |
   |                 |            | with method ``ICNTRL(6)``        |
   |                 +------------+----------------------------------+
   |                 | 4          | Simplified continuous adjoint    |
   |                 |            | with method ``ICNTRL(6)``        |
   +-----------------+------------+----------------------------------+

.. table:: Method to solve the linear Adj equations
   :align: left

   +-----------------+------------+----------------------------------+
   | Integrator      | ICNTRL(7)  | Method selected                  |
   +=================+============+==================================+
   | runge_kutta_adj | 0 or 1     | Modified Newton re-using LU      |
   |                 |            | (default)                        |
   |                 +------------+----------------------------------+
   |                 | 2          | Direct solution (additional one  |
   |                 |            | LU factorizsation of 3Nx3N       |
   |                 |            | matrix per step; good for        |
   |                 |            | debugging                        |
   |                 +------------+----------------------------------+
   |                 | 3          | Adaptive solution (if Newton     |
   |                 |            | does not converge, switch to     |
   |                 |            | direct)                          |
   +-----------------+------------+----------------------------------+
   | sdirk_adj       | 0          | Modified Newton re-using LU      |
   |                 |            | (default)                        |
   | sdirk_tlm       +------------+----------------------------------+
   |                 | 1          | Direct solution (additional one  |
   |                 |            | LU factorizsation per stage)     |
   +-----------------+------------+----------------------------------+

ICNTRL(8)
---------

Determines if LU factorization will be checkpointed at each step
(for adjoint integrators only).

.. table::
   :align: left

   +-----------------+------------+------------------------------------+
   | Integrator      | ICNTRL(8)  | Option selected                    |
   +=================+============+====================================+
   | rosenbrock_adj  | 0          | Do not save LU factorization at    |
   |                 |            | each step (default)                |
   | runge_kutta_adj +------------+------------------------------------+
   |                 | 1          | Save LU factorization at each step |
   | sdirk_adj       |            |                                    |
   +-----------------+------------+------------------------------------+

ICNTRL(9)
---------

Selects options for adjoint and tangent linear method (TLM)
integrators.

.. table:: Selection of adjoint algorithm
   :align: left

   +-----------------+------------+----------------------------------+
   | Integrator      | ICNTRL(9)  | Adjoint algorithm selected       |
   +=================+============+==================================+
   | runge_kutta_adj | 0 or 2     | Discrete adjoint with method     |
   |                 |            | ``ICNTRL(3)`` (default)          |
   |                 +------------+----------------------------------+
   |                 | 1          | No adjoint                       |
   |                 +------------+----------------------------------+
   |                 | 3          | Fully adaptive continous adjoint |
   |                 |            | with method ``ICNTRL(6)``        |
   |                 +------------+----------------------------------+
   |                 | 4          | Simplified continuous adjoint    |
   |                 |            | with method ``ICNTRL(6)``        |
   +-----------------+------------+----------------------------------+

.. table:: Selection of tangent linear method (TLM) error
	      estimation strategy
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(9)  | Strategy selected                 |
   +=================+============+===================================+
   | runge_kutta_tlm | 0          | Base number of iterations as      |
   |                 |            | forward solution                  |
   | sdirk_tlm       +------------+-----------------------------------+
   |                 | 1          | USE ``ATOL_tlm`` and ``RTOL_tlm`` |
   |                 |            | to calculate error estimation     |
   |                 |            | for TLM at Newton stages          |
   +-----------------+------------+-----------------------------------+

ICNTRL(10)
----------

Selects integrator-specific options.

.. table::
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(10) | Option selected                   |
   +=================+============+===================================+
   | lsode           | user       | Maximum order of the integration  |
   |                 | supplied   | formula allowed                   |
   +-----------------+------------+-----------------------------------+
   | runge_kutta     | 0          | Error estimation strategy:        |
   |                 |            | one additional stage at ``c=0``   |
   | runge_kutta_adj |            | (default)                         |
   |                 +------------+-----------------------------------+
   | runge_kutta_tlm | 1          | Error estimation strategy:        |
   |                 |            | Two additional stages at ``c=0``  |
   |                 |            | and SDIRK at ``c=1``, stiffly     |
   |                 |            | accurate                          |
   |                 +------------+-----------------------------------+
   |                 |            |                                   |
   +-----------------+------------+-----------------------------------+
   | seulex          | 0          | No verbose output                 |
   |                 +------------+-----------------------------------+
   |                 | 1          | Dense verbose output              |
   +-----------------+------------+-----------------------------------+

ICNTRL(11)
----------

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------+------------+-----------------------------------+
   | Integrator      | ICNTRL(11) | Option selected                   |
   +=================+============+===================================+
   | seulex          | user       | Maximum number of columns in the  |
   |                 | supplied   | extrapolation. Default is 12.     |
   +-----------------+------------+-----------------------------------+
   | radau5          | 0          | Gustaffson step size control      |
   |                 +------------+-----------------------------------+
   |                 | 1          | Classical step size control       |
   | runge_kutta     |            |                                   |
   |                 |            |                                   |
   | runge_kutta_adj |            |                                   |
   |                 |            |                                   |
   | runge_kutta_tlm |            |                                   |
   |                 |            |                                   |
   |                 |            |                                   |
   +-----------------+------------+-----------------------------------+

ICNTRL(12)
----------

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------------+------------+-----------------------------------+
   | Integrator            | ICNTRL(12) | Option selected                   |
   +=======================+============+===================================+
   | rosenbrock_autoreduce | 0          | Disable mechanism auto-reduction  |
   |                       |            | (i.e. acts in the same way as the |
   |                       |            | ``rosenbrock`` integrator)        |
   |                       +------------+-----------------------------------+
   |                       | 1          | Enables mechanism auto-reduction, |
   |                       |            | set threshold in ``RCNTRL(`12)``  |
   |                       +------------+-----------------------------------+
   |                       |            |                                   |
   +-----------------------+------------+-----------------------------------+
   | rosenbrock_tlm        | 0          | TLM truncation error is not used  |
   |                       +------------+-----------------------------------+
   | runge_kutta_tlm       | 1          | TLM truncation error is computed  |
   |                       |            | and used                          |
   | sdirk_tlm             |            |                                   |
   |                       +------------+-----------------------------------+
   |                       |            |                                   |
   +-----------------------+------------+-----------------------------------+
   | seulex                | 0          | Nsequence = 2 (default)           |
   |                       +------------+-----------------------------------+
   |                       | 1          | Nsequence =                       |
   |                       |            | 1,2,3,4,6,8,12,16,24,32,48,...    |
   |                       +------------+-----------------------------------+
   |                       | 2          | Nsequence =                       |
   |                       |            | 2,3,4,6,8,12,16,24,32,48,64,...   |
   |                       +------------+-----------------------------------+
   |                       | 3          | Nsequence =                       |
   |                       |            | 1,2,3,4,5,6,7,8,9,10,...          |
   |                       +------------+-----------------------------------+
   |                       | 4          | Nsequence =                       |
   |                       |            | 2,3,4,5,6,7,8,9,10,11,...         |
   +-----------------------+------------+-----------------------------------+

ICNTRL(13)
----------

Selects integrator-specific settings.

.. table::
   :align: left

   +-----------------------+------------+-----------------------------------+
   | Integrator            | ICNTRL(13) | Option selected                   |
   +=======================+============+===================================+
   | rosenbrock_autoreduce | 0          | In auto-reduction, disables       |
   |                       |            | scanning species P and L rates    |
   |                       |            | throughout the internal timesteps |
   |                       |            | of the integrator.                |
   |                       +------------+-----------------------------------+
   |                       | 1          | In auto-reduction, ensables       |
   |                       |            | scanning species P and L rates    |
   |                       |            | throughout the internal timesteps |
   |                       |            | of the integrator, for            |
   |                       |            | repartitioning.                   |
   +-----------------------+------------+-----------------------------------+
   | seulex                | 0, 1       | Sets the ``lambda`` parameter     |
   |                       |            | for verbose output.               |
   +-----------------------+------------+-----------------------------------+

ICNTRL(14)
----------

(Solver-specific for :code:`rosenbrock_autoreduce`) If set to
:code:`> 0`, then the threshold is calculated based on the max of
production and loss rate of the species ID specified in
:code:`ICNTRL(14)` multiplied by :code:`RCNTRL(14)`.

ICNTRL(15)
----------

Determines which :code:`Update_*` subroutines are called within the
integrator.

.. list-table::
   :align: left
   :header-rows: 1

   * - ICNTRL(15)
     - Option selcted
   * - -1
     - Do not call any :code:`Update_*` subroutines
   * - 0
     - Use the integrator-specific default values
   * - 1
     - Call  :code:`Update_RCONST` from within the integrator.
   * - 2
     - Call :code:`Update_PHOTO` from within the integrator
   * - 3
     - Call :code:`Update_RCONST` and :code:`Update_PHOTO` from within
       the integrator.
   * - 4
     - Call :code:`Update_SUN` from within the integrator
   * - 5
     - Call :code:`Update_SUN` and :code:`Update_RCONST` from within
       the integrator
   * - 6
     - Call :code:`Update_SUN` and :code:`Update_PHOTO` from within
       the integrator.
   * - 7
     - Call :code:`Update_SUN`, :code:`Update_PHOTO`, and
       :code:`Update_RCONST` from within the integrator.

Calling :code:`Update_RCONST` may be necessary when reaction rate
coefficients depend on the concentration of a specific species, e.
g.:

.. code-block:: console

   HSO3m + HSO5m + Hp = 2 HSO4m + Hp : k_aqueous( C(ind_Hp) );

This ensures that the concentration :code:`C(ind_Hp)` at the specific
integration time is used when the reaction rate coefficient is
updated within the integrator.

ICNTRL(16)
----------

Specifies how negative values should be handled.

.. list-table::
   :align: left
   :header-rows: 1

   * - ICNTRL(16)
     - Option selcted
   * - 0
     - Leave negative values unchanged
   * - 1
     - Set negative values to zero
   * - 2
     - Print warning and continue
   * - 3
     - Print error message and stop

ICNTRL(17)
----------

Selects the amount of verbose output that will be generated.

.. list-table::
   :align: left
   :header-rows: 1

   * - ICNTRL(17)
     - Option selcted
   * - 0
     - Only return error number
   * - 1
     - Full verbose error output

ICNTRL(18)
----------

Currently not used.

ICNTRL(19)
----------

Currently not used.

ICNTRL(20)
----------

Currently not used.

======
RCNTRL
======

.. _table-rcntrl:

.. table:: Summary of RCNTRL usage in the f90 integrators.
           Here, Y = used, and s = solver-specific usage.
   :align: center

   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | RCNTRL                 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 |
   +========================+===+===+===+===+===+===+===+===+===+====+====+====+====+====+====+====+====+====+====+
   | beuler                 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | dvode                  |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | exponential            |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | feuler                 |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | gillespie              |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | lsode                  | Y | Y | Y |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | radau5                 |   | Y |   | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock_adj         | Y | Y | Y | Y | Y | Y | Y |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock             | Y | Y | Y | Y | Y | Y | Y |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock_tlm         | Y | Y | Y | Y | Y | Y | Y |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock_autoreduce  | Y | Y | Y | Y | Y | Y | Y |   |   |    |    | s  |    | s  |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | rosenbrock_h211b_qssa  | Y | Y | Y | Y | Y | Y | Y |   |   |    |    |    |    |    | s  | s  | s  | s  |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | runge_kutta_adj        | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | runge_kutta            | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | runge_kutta_tlm        | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk4                 | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk_adj              | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk                  | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | sdirk_tlm              | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y  | Y  |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | seulex                 | Y | Y | Y | Y | Y | Y | Y | Y |   | s  | s  | s  | s  | s  | s  | s  | s  | s  | s  |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+
   | tau_leap               |   |   |   |   |   |   |   |   |   |    |    |    |    |    |    |    |    |    |    |
   +------------------------+---+---+---+---+---+---+---+---+---+----+----+----+----+----+----+----+----+----+----+

RCNTRL(1)
---------

:code:`Hmin`, the lower bound of the integration step size. It is
not recommended to change the default value of zero.

RCNTRL(2)
---------

:code:`Hmax`, the upper bound of the integration step size.

RCNTRL(3)
---------

:code:`Hstart`, the starting value of the integration step size.

RCNTRL(4)
---------

:code:`FacMin`, lower bound on step decrease factor.

RCNTRL(5)
---------

:code:`FacMax`, upper bound on step increase factor.

RCNTRL(6)
---------

:code:`FacRej`, step decrease factor after multiple rejections.

RCNTRL(7)
---------

:code:`FacSafe`, the factor by which the new step is slightly
smaller than the predicted value.

RCNTRL(8)
---------

:code:`ThetaMin`. If the Newton convergence rate is smaller than
ThetaMin, the Jacobian is not recomputed.

RCNTRL(9)
---------

:code:`NewtonTol`, the stopping criterion for Newton’s method.

RCNTRL(10)
----------

:code:`Qmin`

Different, solver-specific settings for :code:`seulex`.

RCNTRL(11)
----------

:code:`Qmax`. If :code:`Qmin < Hnew/Hold < Qmax`, then the step
size is kept constant and the LU factorization is reused.

Different, solver-specific settings for :code:`seulex`.

RCNTRL(12)
----------

(Solver-specific for :code:`rosenbrock_autoreduce`) Used to specify
the threshold for auto-reduction partitioning, if :code:`ICNTRL(12) = 1`,
and :code:`ICNTRL(14) = 0`. Will be ignored if :code:`ICNTRL(14) > 0`.

Different, solver-specific settings for :code:`seulex`.

RCNTRL(13)
----------

Solver-specific settings for :code:`seulex`.

RCNTRL(14)
----------

(Solver-specific for :code:`rosenbrock_autoreduce`) Used to specify
the multiplier for threshold for auto-reduction partitioning, if
:code:`ICNTRL(12) = 1`, and :code:`ICNTRL(14) > 0`, :code:`RCNTRL(14)`
is multiplied against max  of production and loss rates of species
:code:`ICNTRL(14)` to produce the partitioning threshold, ignoring
:code:`RCNTRL(12)`.

Different, solver-specific settings for :code:`seulex`.

RCNTRL(15) - RCNTRL(18)
-----------------------

Solver-specific settings for :code:`rosenbrock_h211b_qssa`.

Different, solver-specific settings for :code:`seulex`.

RCNTRL(19)
----------

Solver-specific settings for :code:`seulex`.

RCNTRL(20)
----------

Currently not used.
