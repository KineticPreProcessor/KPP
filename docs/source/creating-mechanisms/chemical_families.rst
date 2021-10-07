###################################
Adding production and loss families
###################################

Certain common families (e.g. POx, LOx) have been pre-defined for you.
You will find the family definitions near the top of the
:file:`gckpp.kpp` file:

.. code-block:: none

  #FAMILIES
  POx : O3 + NO2 + 2NO3 + PAN + PPN + MPAN + HNO4 + 3N2O5 + HNO3 + BrO + HOBr + BrNO2 + 2BrNO3 + MPN + ETHLN + MVKN + MCRHN + MCRHNB + PROPNN + R4N2 + PRN1 + PRPN + R4N1 + HONIT + MONITS + MONITU + OLND + OLNN + IHN1 + IHN2 + IHN3 + IHN4 + INPB + INPD + ICN + 2IDN + ITCN + ITHN + ISOPNOO1 + ISOPNOO2 + INO2B + INO2D + INA + IDHNBOO + IDHNDOO1 + IDHNDOO2 + IHPNBOO + IHPNDOO + ICNOO + 2IDNOO + MACRNO2 + ClO + HOCl + ClNO2 + 2ClNO3 + 2Cl2O2 + 2OClO + O + O1D + IO + HOI + IONO + 2IONO2 + 2OIO + 2I2O2 + 3I2O3 + 4I2O4;
  LOx : O3 + NO2 + 2NO3 + PAN + PPN + MPAN + HNO4 + 3N2O5 + HNO3 + BrO + HOBr + BrNO2 + 2BrNO3 + MPN + ETHLN + MVKN + MCRHN + MCRHNB + PROPNN + R4N2 + PRN1 + PRPN + R4N1 + HONIT + MONITS + MONITU + OLND + OLNN + IHN1 + IHN2 + IHN3 + IHN4 + INPB + INPD + ICN + 2IDN + ITCN + ITHN + ISOPNOO1 + ISOPNOO2 + INO2B + INO2D + INA + IDHNBOO + IDHNDOO1 + IDHNDOO2 + IHPNBOO + IHPNDOO + ICNOO + 2IDNOO + MACRNO2 + ClO + HOCl + ClNO2 + 2ClNO3 + 2Cl2O2 + 2OClO + O + O1D + IO + HOI + IONO + 2IONO2 + 2OIO + 2I2O2 + 3I2O3 + 4I2O4;
  PCO : CO;
  LCO : CO;
  PSO4 : SO4;
  LCH4 : CH4;
  PH2O2 : H2O2;

.. note:: The POx, LOx, PCO, and LCO families are used for computing budgets in
          the GEOS-Chem benchmark simulations.  PSO4 is required for simulations using `TOMAS aerosol microphysics <TOMAS_aerosol_microphysics>`__.

To add a new prod/loss family, add a new line to the :code:`#FAMILIES`
section with the format

.. code-block:: none

  FAM_NAME : MEMBER_1 + MEMBER_2 + ... + MEMBER_N;

The family name must start with :code:`P` or :code:`L` to indicate
whether KPP should calculate a production or a loss rate.

The maximum number of families allowed by KPP is currently set to 300.
Depending on how many prod/loss families you add, you may need to
increase that to a larger number to avoid errors in KPP. You can change
the number for :code:`MAX_FAMILIES` in :file:`KPP/kpp-code/src/gdata.h` and then `rebuild
KPP <FlexChem#KPP_source_code>`__.

.. code-block:: C

   #define MAX_EQN        1500    /* KPP 2.3.0_gc, Bob Yantosca (11 Feb 2021)  */
   #define MAX_SPECIES    1000    /* KPP 2.3.0_gc, Bob Yantosca (11 Feb 2021)  */
   #define MAX_SPNAME       30
   #define MAX_IVAL         40
   #define MAX_EQNTAG       12    /* Max length of equation ID in eqn file     */
   #define MAX_K           150    /* Max length of rate expression in eqn file */
   #define MAX_ATOMS        10
   #define MAX_ATNAME       10
   #define MAX_ATNR        250
   #define MAX_PATH        120
   #define MAX_FILES        20
   #define MAX_FAMILIES    300
   #define MAX_MEMBERS     150
   #define MAX_EQNLEN      200

.. important:: When adding a prod/loss family or changing any of the
	       other settings in :file:`gckpp.kpp`, you must re-run KPP to
	       produce new Fortran-90 files for GEOS-Chem (:ref:`as described in a
	       previous chapter <creating_fortran_files>`).

Production and loss families are archived via the HISTORY diagnostics.
For more information, please see the `Guide to GEOS_Chem History
diagnostics <http://wiki.geos-chem.org/Guide_to_GEOS_Chem_History_diagnostics>`__
on the GEOS-Chem wiki.
