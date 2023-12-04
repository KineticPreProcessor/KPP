# MCM minimal example 

By Rolf Sander (@RolfSander)

Use these commands to run the MCM minimal example:

1. Run KPP to generate Fortran90 solver files for the MCM minimal example:

   ```console
   $ kpp mcm.kpp
   ```
  
2. Compile the KPP-generated source code to an executable:

   ```console
   $ gmake clean
   $ gmake -f Makefile_mcm EXTERNAL_RATES_F90=constants_mcm.f90
   ```
   NOTE: On some systems, `gmake` may be installed as `make`.


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
