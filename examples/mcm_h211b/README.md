# MCM minimal example (rosenbrock_h211b_qssa)

By Rolf Sander

The Master Chemical Mechanism (MCM) offers the export of a selected
mechanism in KPP format. The downloaded files can be used out-of-the-box
as input for the example model in this directory. Here, the isoprene
degradation mechanism is used as an example.

Before running the example, ensure that KPP is [installed](https://kpp.readthedocs.io/en/stable/getting_started/01_installation.html) properly. The environment variable `$KPP_HOME` must be defined, and the KPP executable `$KPP_HOME/bin/kpp` should be included in the`$PATH`.

## Downloading a mechanism from the MCM web page

1. Browse through the [MCM](https://mcm.york.ac.uk/MCM/browse) and
   select a Subset of the mechanism.

2. Go to the [export page](https://mcm.york.ac.uk/MCM/export) and choose
   KPP as the output format.

3. Click on the `Download` button and rename the downloaded file to
   `mcm_isoprene.eqn`.

4. Download the auxiliary file
   [constants_mcm.f90](https://mcm.york.ac.uk/MCM/export/kpp_constants).


## Executing the MCM minimal example (with the provided script)

Run one of the example scripts provided:
 
   ```console
   $ ./run_example.csh   # If you are using csh or tcsh
   $ ./run_example.sh    # If you are using bash
   ```
## Executing the MCM minimal example (manually)

If you wish to execute the steps manually, follow these instructions:

1. Run KPP to generate Fortran90 solver files for the MCM minimal example:

   ```console
   $ kpp mcm.kpp
   ```
  
2. Compile the KPP-generated source code to an executable:

   ```console
   $ gmake -f Makefile_mcm distclean
   $ gmake -j -f Makefile_mcm EXTERNAL_RATES_F90=constants_mcm.f90
   ```
   NOTE: On some systems, `gmake` may be installed as `make`.

   ALSO NOTE: At present, only a single external F90 module file (with
   rate constants and parameters for the MCM mechanism) can be specified
   with the `EXTERNAL_RATES_F90` environment variable. To specify more
   than one external file you will have to modify the `util/Makfile_f90`
   and/or `util/Makfile_upper_F90` to add additional rules.

3. Run the MCM minimal example executable:

   ```console
   $ ./mcm.exe
   ```

4. If Python 3 and matplotlib are available (e.g., installed in a
   conda/mamba or a pip environment), plot the results with:
   
   ```console
   $ python3 plot_data.py
   ```
   
5. Remove all KPP-generated source code files and output files from the MCM minimal example once you no longer need them:

   ```console
   $ gmake -f Makefile_mcm distclean
   ```
