#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: Set here the F90 compiler and options
#       Pedefined compilers: INTEL, PGF, HPUX, LAHEY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pick a default compiler if not supplied
ifndef COMPILER
#COMPILER = LAHEY
#COMPILER = INTEL
#COMPILER = PGF
#COMPILER = HPUX
COMPILER = GFORTRAN
endif

FC_LAHEY   = lf95
# More aggressive for production runs:
#FOPT_LAHEY = -Cpp --pca -O
# More checking for debugging:
FOPT_LAHEY = -Cpp --chk a,e,s,u --pca --ap -O0 -g --trap --trace --chkglobal

FC_INTEL   = ifort
# More aggressive for production runs:
#FOPT_INTEL = -cpp -O -fp-model precise -pc80 -prec_div
# More checking for debugging:
FOPT_INTEL = -cpp -O0 -fp-model strict -implicitnone -ftrapuv \
              -debug all -check all -warn all -g

FC_PGF     = pgf90
# More aggressive for production runs:
FOPT_PGF   = -Mpreprocess -O -fast -pc 80 -Kieee
# More checking for debugging:
#FOPT_PGF   = -Mpreprocess -O0 -Mbounds -Mchkfpstk -Mchkptr -Mchkstk \
#             -Ktrap=fp -pc 80 -Kieee

FC_HPUX    = f90
FOPT_HPUX  = -O -u +Oall +check=on

FC_GFORTRAN     = gfortran
FOPT_GFORTRAN   = -cpp -O -g

# define FULL_ALGEBRA for non-sparse integration
FC   = $(FC_$(COMPILER))
FOPT = $(FOPT_$(COMPILER)) # -DFULL_ALGEBRA

LIBS =
#LIBS = -llapack -lblas

#############################################################################
# Test if various files are found on disk or not
# so that we can include/exclude them from the Make process
#############################################################################

# JacobianSP
HAS_JAC_SP =0
ifeq ($(shell test -f mcm_JacobianSP.f90 && echo y),y)
  HAS_JAC_SP =1
endif

# Jacobian
HAS_JAC = 0
ifeq ($(shell test -f mcm_Jacobian.f90 && echo y),y)
  HAS_JAC = 1
endif

# HessianSP
HAS_HES_SP = 0
ifeq ($(shell test -f mcm_HessianSP.f90 && echo y),y)
  HAS_HES_SP = 1
endif

# Hessian
HAS_HES=0
ifeq ($(shell test -f mcm_Hessian.f90 && echo y),y)
  HAS_HES = 1
endif

# StoichiomSP
HAS_STM_SP = 0
ifeq ($(shell test -f  mcm_StoichiomSP.f90 && echo y),y)
  HAS_STM_SP=1
endif

# Stoichiom
HAS_STM = 0
ifeq ($(shell test -f  mcm_Stoichiom.f90 && echo y),y)
  HAS_STM=1
endif


# Command to create Matlab mex gateway routines
# Note: use $(FC) as the mex Fortran compiler
MEX  = mex

GENSRC = mcm_Precision.f90  mcm_Parameters.f90  mcm_Global.f90 constants_mcm.f90
GENOBJ = mcm_Precision.o    mcm_Parameters.o    mcm_Global.o   constants_mcm.o

FUNSRC = mcm_Function.f90
FUNOBJ = mcm_Function.o

JACSRC   =
JACOBJ   =
JACSPOBJ =
ifeq ($(HAS_JAC_SP),1)
  JACSRC   += mcm_JacobianSP.f90
  JACOBJ   += mcm_JacobianSP.o
  JACSPOBJ += mcm_JacobianSP.o
endif
ifeq ($(HAS_JAC),1)
  JACSRC   += mcm_Jacobian.f90
  JACOBJ   += mcm_Jacobian.o
endif

HESSRC   =
HESOBJ   =
HESSPOBJ =
ifeq ($(HAS_HES_SP),1)
  HESSRC   += mcm_HessianSP.f90
  HESOBJ   += mcm_HessianSP.o
  HESSPOBJ += mcm_HessianSP.o
endif
ifeq ($(HAS_HES),1)
  HESSRC   += mcm_Hessian.f90
  HESOBJ   += mcm_Hessian.o
endif

STMSRC   =
STMOBJ   =
STMSPOBJ =
ifeq ($(HAS_STM_SP),1)
  STMSRC   += mcm_StoichiomSP.f90
  STMOBJ   += mcm_StoichiomSP.o
  STMSPOBJ += mcm_StoichiomSP.o
endif
ifeq ($(HAS_HES),1)
  HESSRC   += mcm_Stoichiom.f90
  HESOBJ   += mcm_Stoichiom.o
endif

UTLSRC = mcm_Rates.f90 mcm_Util.f90 mcm_Monitor.f90
UTLOBJ = mcm_Rates.o   mcm_Util.o   mcm_Monitor.o

LASRC  = mcm_LinearAlgebra.f90
LAOBJ  = mcm_LinearAlgebra.o

STOCHSRC = mcm_Stochastic.f90
STOCHOBJ = mcm_Stochastic.o

MAINSRC = mcm_Main.f90        mcm_Initialize.f90 \
          mcm_Integrator.f90  mcm_Model.f90
MAINOBJ = mcm_Main.o          mcm_Initialize.o   \
          mcm_Integrator.o    mcm_Model.o

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: modify the line below to include only the
#       objects needed by your application
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ALLOBJ = $(GENOBJ) $(FUNOBJ) $(JACOBJ) $(HESOBJ) $(STMOBJ) \
	 $(UTLOBJ) $(LAOBJ)

#############################################################################
# Makefile targets and rules
#############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# User: modify the line below to include only the
#       executables needed by your application
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all:    exe

exe:	$(ALLOBJ) $(MAINOBJ)
	$(FC) $(FOPT) $(ALLOBJ) $(MAINOBJ) $(LIBS) -o mcm.exe

stochastic:$(ALLOBJ) $(STOCHOBJ) $(MAINOBJ)
	$(FC) $(FOPT) $(ALLOBJ) $(STOCHOBJ) $(MAINOBJ) $(LIBS) \
	-o mcm_stochastic.exe

mex:    $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O mcm_mex_Fun.f90     $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O mcm_mex_Jac_SP.f90  $(ALLOBJ)
	$(MEX) FC#$(FC) -fortran -O mcm_mex_Hessian.f90 $(ALLOBJ)

clean:
	rm -f *.o *.mod *~ \
	mcm*.dat mcm.exe mcm*.mexglx \
	mcm.log

distclean:
	rm -f *.o *.mod *~ \
	mcm*.dat mcm.exe mcm.log \
	mcm_*.f90 mcm_*.f90 mcm_*.mexglx \
	results*.m Makefile_mcm

mcm_Precision.o: mcm_Precision.f90
	$(FC) $(FOPT) -c $<

constants_mcm.o: constants_mcm.f90 mcm_Precision.o mcm_Global.f90 mcm_Parameters.o
	$(FC) $(FOPT) -c $<

mcm_Parameters.o: mcm_Parameters.f90 mcm_Precision.o
	$(FC) $(FOPT) -c $<

mcm_Monitor.o: mcm_Monitor.f90 mcm_Precision.o
	$(FC) $(FOPT) -c $<

mcm_Global.o: mcm_Global.f90 mcm_Parameters.o mcm_Precision.o
	$(FC) $(FOPT) -c $<

mcm_Initialize.o: mcm_Initialize.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

mcm_Function.o: mcm_Function.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

mcm_Stochastic.o: mcm_Stochastic.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

mcm_JacobianSP.o: mcm_JacobianSP.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

mcm_Jacobian.o: mcm_Jacobian.f90 $(GENOBJ) $(JACSPOBJ)
	$(FC) $(FOPT) -c $<

mcm_LinearAlgebra.o: mcm_LinearAlgebra.f90 $(GENOBJ) $(JACSPOBJ)
	$(FC) $(FOPT) -c $<

mcm_Rates.o: mcm_Rates.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

mcm_HessianSP.o: mcm_HessianSP.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

mcm_Hessian.o:  mcm_Hessian.f90 $(GENOBJ) $(HESSPOBJ)
	$(FC) $(FOPT) -c $<

mcm_StoichiomSP.o: mcm_StoichiomSP.f90 $(GENOBJ)
	$(FC) $(FOPT) -c $<

mcm_Stoichiom.o: mcm_Stoichiom.f90 $(GENOBJ) $(STMSPOBJ)
	$(FC) $(FOPT) -c $<

mcm_Util.o: mcm_Util.f90 $(GENOBJ) mcm_Monitor.o
	$(FC) $(FOPT) -c $<

mcm_Main.o: mcm_Main.f90  $(ALLOBJ) \
                 mcm_Initialize.o mcm_Model.o mcm_Integrator.o
	$(FC) $(FOPT) -c $<

mcm_Model.o: mcm_Model.f90 $(ALLOBJ) mcm_Integrator.o
	$(FC) $(FOPT) -c $<

mcm_Integrator.o: mcm_Integrator.f90 $(ALLOBJ)
	$(FC) $(FOPT) -c $<
