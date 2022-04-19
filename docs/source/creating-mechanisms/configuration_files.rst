#################################
User-editable configuration files
#################################

In the :file:`KPP/custom` folder within the GEOS-Chem source code, you
will find two files that define the chemical mechanism:

----------
custom.eqn
----------

The :file:`custom.eqn` configuration file contains:

- List of active species
- List of inactive species
- Gas-phase reactions
- Heterogeneous reactions
- Photolysis reactions

---------
gckpp.kpp
---------
  
The :file:`gckpp.kpp` configuration file contains:

- Solver options
- Production and loss family definitions
- Functions to compute reaction rates
- Global definitions
