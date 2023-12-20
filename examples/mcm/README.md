# MCM minimal example 

By Rolf Sander

The Master Chemical Mechanism (MCM) offers the export of a selected
mechanism in KPP format. The downloaded files can be used out-of-the-box
as input for the example model in this directory. Here, the isoprene
degradation mechanism is shown as an example.

## Downloading a mechanism from the MCM web page

1. Browse through the [mcm](https://mcm.york.ac.uk/MCM/browse) and
   select a Subset of the mechanism.

2. Go  to the [export page](https://mcm.york.ac.uk/MCM/export) and
   choose KPP as the output format.

3. Click on the Download button and rename the file to
   `mcm_isoprene.eqn`.

4. Download the auxiliary file
   [constants_mcm.f90](https://mcm.york.ac.uk/MCM/export/kpp_constants).

## Executing the MCM minimal example

1. Ensure that KPP is
   [installed](https://kpp.readthedocs.io/en/stable/getting_started/01_installation.html)
   properly. The environment variable `$KPP_HOME` must be defined, and
   the KPP executable `$KPP_HOME/bin/kpp` should be included in the
   `$PATH`.

2. Run KPP to generate Fortran90 solver files for the MCM minimal example:

   ```console
   $ kpp mcm.kpp
   ```
  
3. Compile the KPP-generated source code to an executable:

   ```console
   $ gmake -f Makefile_mcm clean
   $ gmake -f Makefile_mcm EXTERNAL_RATES_F90=constants_mcm.f90
   ```
   NOTE: On some systems, `gmake` may be installed as `make`.

   ALSO NOTE: At present, only a single external F90 module file (with rate constants and parameters for the MCM mechanism) can be specified with the `EXTERNAL_RATES_F90` environment variable.  To specify more than one external file you will have to modify the `util/Makfile_f90` and/or `util/Makfile_upper_F90` to add additional rules.

4. Run the MCM minimal example executable:

   ```console
   $ ./mcm.exe
   ```

5. If Python 3 and matplotlib are available (e.g., installed in a conda/mamba or a pip environment), plot the results with:
   
   ```console
   $ python3 plot_data.py
   ```
   
6. Remove all KPP-generated source code files and output files from the MCM minimal example once you no longer need them:

   ```console
   $ gmake distclean
   ```
