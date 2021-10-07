#################################
Changing the numerical integrator
#################################

Several global options for :program:`KPP` are listed at the top of the
:file:`gckpp.kpp` file: 

.. code-block:: none

   #INTEGRATOR rosenbrock
   #LANGUAGE Fortran90
   #DRIVER none
   #HESSIAN off
   #MEX off
   #STOICMAT off

The :code:`#INTEGRATOR` tag specifies the choice of numerical
integrator that you wish to use with your chemical mechanism.  The
Rosenbrock solver is used by default.

.. important::  We do not recommend changing the value of
		:code:`#INTEGRATOR`.

However, if you wish to use a different integrator for research
purposes, you may select from one of the following options:

   - exponential
   - gillespie
   - kpp_dvode
   - kpp_lsode
   - kpp_radau5
   - kpp_sdirk4
   - kpp_seulex
   - none
   - rosenbrock
   - rosenbrock_adj
   - rosenbrock_split
   - rosenbrock_tlm
   - runge_kutta
   - runge_kutta_adj
   - runge_kutta_tlm
   - sdirk
   - sdirk_adj
   - sdirk_tlm
   - tau_leap
     
The :code:`#LANGUAGE` setting should be set to :code:`Fortran90`.

The other options should be left as they are, as they are not relevant
to :program:`GEOS-Chem`.

For more information about :program:`KPP` settings, please see the
`KPP 2.1 user manual <https://github.com/geoschem/KPP/blob/main/kpp-2.2.3_01/doc/kpp_UserManual.pdf>`__. 
