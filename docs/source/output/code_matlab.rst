.. _matlab-code:

###############
The Matlab code
###############

.. important::

   Some run-time options for Matlab-language integrators (specified in
   the :ref:`ICNTRL and RCNTRL arrays <icntrl-rcntrl>`) do not exactly
   correspond to the Fortran90 run-time options.  We will standardize
   run-time integrator options across all target languages in a future
   KPP release.

`Matlab <http://www.mathworks.com/products/matlab/>`_ provides a
high-level programming environment that allows algorithm development,
numerical computations, and data analysis and visualization. The
KPP-generated Matlab code allows for a rapid prototyping of chemical
kinetic schemes, and for a convenient analysis and visualization of the
results. Differences between different kinetic mechanisms can be easily
understood. The Matlab code can be used to derive reference numerical
solutions, which are then compared against the results obtained with
user-supplied numerical techniques. KPP/Matlab can also be used to teach
students fundamentals of chemical kinetics and chemical numerical
simulations.

Each Matlab function has to reside in a separate m-file. Function calls
use the m-function-file names to reference the function. Consequently,
KPP generates one m-function-file for each of the functions discussed in
the sections entitled :ref:`Function` ,
:ref:`Jacobian-and-JacobianSP`, :ref:`Hessian-and-HessianSP`,
:ref:`Stoichiom-and-StoichiomSP`, :ref:`Util`.  The names of the
m-function-files are the same as the names of the functions (prefixed
by the model name :code:`ROOT`.

The variables of :ref:`table-par` are defined as Matlab :code:`global`
variables and initialized in the file
:file:`ROOT_parameter_defs.m`. The variables of :ref:`table-glob` are
declared as Matlab :code:`global` variables in the file
:file:`ROOT_global_defs.m`. They can be accessed from within each
Matlab function  by using declarations of the variables of interest.

The sparse data structures for the Jacobian (cf. :ref:`table-jac`), the Hessian
(cf. :ref:`table-hess`), the stoichiometric matrix (cf. :ref:`table-sto`),
and the Jacobian of reaction (see :ref:`table-jvrp`) are declared as
Matlab :code:`global` variables in the file
:file:`ROOT_Sparse_defs.m`.  They are initialized in separate m-files,
namely :file:`ROOT_JacobianSP.m`, :file:`ROOT_HessianSP.m`, and
:file:`ROOT_StoichiomSP.m` respectively.

Two wrappers (:file:`ROOT_Fun_Chem.m` and :file:`ROOT_Jac_SP_Chem.m`) are
provided for interfacing the ODE function and the sparse ODE Jacobian
with Matlab’s suite of ODE integrators. Specifically, the syntax of
the wrapper calls matches the syntax required by Matlab’s integrators
like ode15s. Moreover, the Jacobian wrapper converts the sparse KPP
format into a Matlab sparse matrix.

.. _table-matlab:

.. list-table:: List of Matlab model files
   :align: center
   :widths: 35 65
   :header-rows: 1

   * - Function
     - Description
   * - :file:`ROOT.m`
     - Driver
   * - :file:`ROOT_parameter_defs.m`
     - Global parameters
   * - :file:`ROOT_global_defs.m`
     - Global variables
   * - :file:`ROOT_sparse_defs.m`
     - Global sparsity data
   * - :file:`ROOT_Fun_Chem.m`
     - Template for ODE function
   * - :file:`ROOT_Fun.m`
     - ODE function
   * - :file:`ROOT_Jac_Chem.m`
     - Template for ODE Jacobian
   * - :file:`ROOT_Jac_SP.m`
     - Jacobian in sparse format
   * - :file:`ROOT_JacobianSP.m`
     - Sparsity data structures
   * - :file:`ROOT_Hessian.m`
     - ODE Hessian in sparse format
   * - :file:`ROOT_HessianSP.m`
     - Sparsity data structures
   * - :file:`ROOT_Hess_Vec.m`
     - Hessian action on vectors
   * - :file:`ROOT_HessTR_Vec.m`
     - Transposed Hessian action on vectors
   * - :file:`ROOT_stoichiom.m`
     - Derivatives of Fun and Jac w/r/t rate coefficients
   * - :file:`ROOT_stoichiomSP.m`
     - Sparse data
   * - :file:`ROOT_ReactantProd.m`
     - Reactant products
   * - :file:`ROOT_JacReactantProd.m`
     - Jacobian of reactant products
   * - :file:`ROOT_Rates.m`
     - User-defined rate reaction laws
   * - :file:`ROOT_Update_PHOTO.m`
     - Update photolysis rate coefficients
   * - :file:`ROOT_Update_RCONST.m`
     - Update all rate coefficients
   * - :file:`ROOT_Update_SUN.m`
     - Update sola intensity
   * - :file:`ROOT_GetMass.m`
     - Check mass balance for selected atoms
   * - :file:`ROOT_Initialize.m`
     - Set initial values
   * - :file:`ROOT_Shuffle_kpp2user.m`
     - Shuffle concentration vector
   * - :file:`ROOT_Shuffle_user2kpp.m`
     - Shuffle concentration vector
