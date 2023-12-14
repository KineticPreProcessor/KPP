# MCM minimal example 

By Rolf Sander (@RolfSander)

Use these commands to run the MCM minimal example:

1. Run KPP to generate Fortran90 solver files for the MCM minimal example:

   ```console
   $ kpp mcm.kpp
   ```
  
2. Compile the KPP-generated source code to an executable:

   ```console
   $ gmake -f Makefile_mcm clean
   $ gmake -f Makefile_mcm EXTERNAL_RATES_F90=constants_mcm.f90
   ```
   NOTE: On some systems, `gmake` may be installed as `make`.

   ALSO NOTE: At present, only a single external F90 module file (with rate constants and parameters for the MCM mechanism) can be specified with the `EXTERNAL_RATES_F90` environment variable.  To specify more than one external file you will have to modify the `util/Makfile_f90` and/or `util/Makfile_upper_F90` to add additional rules.


3. Run the MCM minimal example executable:

   ```console
   $ ./mcm.exe
   ```

4. Plot results with matplotlib (make sure you have matplotlib installed in a conda/mamba environment or in a pip environment):
   
   ```console
   $ python3 plot_data.py
   ```
   
5. Remove all KPP-generated source code files and output files from the MCM minimal example once you no longer need them:
   ```console
   $ gmake distclean
   ```
