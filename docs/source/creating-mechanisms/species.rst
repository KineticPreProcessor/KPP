##############
Adding species
##############

-------------------------
Chemically-active species
-------------------------

List chemically-active (aka variable) species in the :code:`#DEFVAR`
section of :file:`custom.eqn`, as shown below:

.. code-block:: none

  #DEFVAR
  A3O2       = IGNORE; {CH3CH2CH2OO; Primary RO2 from C3H8}
  ACET       = IGNORE; {CH3C(O)CH3; Acetone}
  ACTA       = IGNORE; {CH3C(O)OH; Acetic acid}
  ...etc ...

---------------------------
Chemically-inactive species
---------------------------

List species whose concentrations do not change in the :code:`#DEFFIX`
of :file:`custom.eqn`, as shown below:

.. code-block:: none

  #DEFFIX
  H2         = IGNORE; {H2; Molecular hydrogen}
  N2         = IGNORE; {N2; Molecular nitrogen}
  O2         = IGNORE; {O2; Molecular oxygen}
  ... etc ...

Species may be listed in any order, but we have found it convenient to
list them alphabetically.
