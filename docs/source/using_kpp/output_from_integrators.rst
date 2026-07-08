.. _output-from-int:

#################################################################
Output from the Integrators (:code:`ISTATUS` and :code:`RSTATUS`)
#################################################################

In order to obtain more information about the integration, KPP provides
the arrays :code:`ISTATUS` (integer) and :code:`RSTATUS` (real). Each of
them is an array of 20 elements. Array elements not listed here are
currently not used. Details can be found in the comment lines of the
individual integrator files in :code:`$KPP_HOME/int/`.

.. _output-from-int-istatus:

=======
ISTATUS
=======

.. _table-istatus:

.. table:: Summary of ISTATUS usage in the f90 integrators.
           Here, Y = used.
   :align: center

   +----------------------------+---+---+---+---+---+---+---+---+---+
   | ISTATUS                    | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
   +============================+===+===+===+===+===+===+===+===+===+
   | beuler                     | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | dvode                      |   |   |   |   |   |   |   |   |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | exponential                |   |   |   |   |   |   |   |   |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | feuler                     |   |   |   |   |   |   |   |   |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | gillespie                  |   |   |   |   |   |   |   |   |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | lsode                      | Y | Y | Y |   |   |   |   |   |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | radau5                     | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | rosenbrock                 | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | rosenbrock_adj             | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | rosenbrock_autoreduce      | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | rosenbrock_h211b_qssa      | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | rosenbrock_tlm             | Y | Y | Y | Y | Y | Y | Y | Y | Y |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | runge_kutta                | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | runge_kutta_adj            | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | runge_kutta_tlm            | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | sdirk                      | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | sdirk_adj                  | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | sdirk_tlm                  | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | sdirk4                     | Y | Y | Y | Y | Y | Y | Y | Y |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | seulex                     | Y | Y | Y | Y | Y | Y | Y |   |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+
   | tau_leap                   |   |   |   |   |   |   |   |   |   |
   +----------------------------+---+---+---+---+---+---+---+---+---+

ISTATUS(1)
----------

Number of times that the :ref:`function containing the chemical ODE system
<Function>` was called.

ISTATUS(2)
----------

Number times :ref:`the Jacobian matrix <Jacobian-and-JacobianSP>` was constructed.

ISTATUS(3)
----------

Total number of integration timesteps.

ISTATUS(4)
----------

Number of accepted timesteps.

ISTATUS(5)
----------

Number of rejected timesteps (except at very beginning).

ISTATUS(6)
----------

Number of LU decompositions that were performed.

ISTATUS(7)
----------

Number of forward/backward substitutions that were performed.

ISTATUS(8)
----------

Number of singular matrix decompositions that were performed.

ISTATUS(9)
----------

Number of times :ref:`the Hessian matrix <Hessian-and-HessianSP>` was evaluated.

ISTATUS(10) .. ISTATUS(20)
--------------------------

Currently not used.

=======
RSTATUS
=======

.. _table-rstatus:

.. table:: Summary of RSTATUS usage in the f90 integrators.
           Here, Y = used, s = solver specific usage.
   :align: center

   +----------------------------+---+---+---+---+
   | RSTATUS                    | 1 | 2 | 3 | 4 |
   +============================+===+===+===+===+
   | beuler                     | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | dvode                      |   |   |   |   |
   +----------------------------+---+---+---+---+
   | exponential                |   |   |   |   |
   +----------------------------+---+---+---+---+
   | feuler                     | Y |   |   |   |
   +----------------------------+---+---+---+---+
   | gillespie                  |   |   |   |   |
   +----------------------------+---+---+---+---+
   | lsode                      | Y | Y |   |   |
   +----------------------------+---+---+---+---+
   | radau5                     |   |   |   |   |
   +----------------------------+---+---+---+---+
   | rosenbrock                 | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | rosenbrock_adj             | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | rosenbrock_autoreduce      | Y | Y | Y | s |
   +----------------------------+---+---+---+---+
   | rosenbrock_h211b_qssa      | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | rosenbrock_tlm             | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | runge_kutta                | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | runge_kutta_adj            | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | runge_kutta_tlm            | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | sdirk                      | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | sdirk_adj                  | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | sdirk_tlm                  | Y | Y | Y |   |
   +----------------------------+---+---+---+---+
   | sdirk4                     | Y | Y |   |   |
   +----------------------------+---+---+---+---+
   | seulex                     |   |   |   |   |
   +----------------------------+---+---+---+---+
   | tau_leap                   |   |   |   |   |
   +----------------------------+---+---+---+---+

RSTATUS(1)
----------

:code:`Texit`, the time corresponding to the computed :math:`Y`
upon return.

RSTATUS(2)
----------

:code:`Hexit`: the last accepted step before exit.

RSTATUS(3)
----------

:code:`Hnew`: The last predicted step (not yet taken.  For multiple
restarts, use :code:`Hnew` as :code:`Hstart` in the subsequent run.

RSTATUS(4)
----------

(Solver-specific for :code:`rosenbrock_autoreduce`) :code:`AR_thr`:
used to output the calculated (used) auto-reduction threshold for
the integration. Useful when :code:`ICNTRL(10) > 0` where the
threshold is dynamically determined based on a given species.

RSTATUS(5) .. RSTATUS(20)
-------------------------

Currently not used.
