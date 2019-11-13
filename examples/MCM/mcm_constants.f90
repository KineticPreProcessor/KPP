MODULE constants

  USE model_Precision, ONLY: dp
  IMPLICIT NONE
  PRIVATE dp

  INTEGER, PARAMETER :: mnsp=250, mre=2000
  INTEGER i
  ! variables for zenith routine which calculates zenith angle
  REAL(dp) theta, secx, cosx
  ! generic reaction rate variables
  REAL(dp) kro2no, kro2ho2, kapho2, kapno, kro2no3, kno3al, kdec, &
    krosec, kalkoxy, kalkpxy, kroprim
  ! variables for calculation of kfpan and kbpan
  REAL(dp) kfpan, kbpan
  REAL(dp) kc0, kci, krc, fcc, nc, fc
  REAL(dp) kd0, kdi, krd, fcd, ncd, fd
  ! variables for calculation of kmt01
  REAL(dp) kmt01
  REAL(dp) k10, k1i, kr1, fc1, nc1, f1
  ! variables for calculation of kmt02
  REAL(dp) kmt02
  REAL(dp) k20, k2i, kr2, fc2, nc2, f2
  ! variables for calculation of kmt03
  REAL(dp) kmt03
  REAL(dp) k30, k3i, kr3, fc3, nc3, f3
  ! variables for calculation of kmt04
  REAL(dp) kmt04
  REAL(dp) k40, k4i, kr4, fc4, nc4, f4
  ! variables for calculation of kmt05
  REAL(dp) kmt05
  ! variables for calculation of kmt06
  REAL(dp) kmt06
  ! variables for calculation of kmt07
  REAL(dp) kmt07
  REAL(dp) k70, k7i, kr7, fc7, nc7, f7
  ! variables for calculation of kmt08
  REAL(dp) kmt08
  REAL(dp) k80, k8i, kr8, fc8, nc8, f8
  ! variables for calculation of kmt09
  REAL(dp) kmt09
  REAL(dp) k90, k9i, kr9, fc9, nc9, f9
  ! variables for calculation of kmt10
  REAL(dp) kmt10
  REAL(dp) k100, k10i, kr10, fc10, nc10, f10
  ! variables for calculation of kmt11
  REAL(dp) kmt11
  REAL(dp) k1,k2,k3,k4
  ! variables for calculation of kmt12
  REAL(dp) kmt12
  REAL(dp) k0, ki, x, ssign,f
  REAL(dp) k120, k12i, kr12, fc12, nc12, f12
  ! variables for calculation of kmt13
  REAL(dp) kmt13
  REAL(dp) k130, k13i, kr13, nc13, fc13, f13
  ! variables for calculation of kmt14
  REAL(dp) kmt14
  REAL(dp) k140, k14i, kr14, nc14, fc14, f14
  ! variables for calculation of kmt15
  REAL(dp) kmt15
  REAL(dp) k150, k15i, kr15, fc15, nc15, f15
  ! variables for calculation of kmt16
  REAL(dp) kmt16
  REAL(dp) k160, k16i, kr16, nc16, fc16, f16
  ! variables for calculation of kmt17
  REAL(dp) kmt17
  REAL(dp) k170, k17i, kr17, fc17, nc17, f17
  ! variables for calculation of kmt17
  REAL(dp) kmt18

  REAL(dp) kch3o2, k298ch3o2
  
  REAL(dp) kp_rooh, kp_ho2, kp_n2o5, kp_no3, kp_all
  
  INTEGER k

CONTAINS

  !***************************************************************************

  SUBROUTINE mcm_constants(time, temp, M, N2, O2, RO2, H2O)
    ! calculates rate constants from arrhenius informtion
    USE model_Global,  ONLY:LAT, LON, JDAY, SZAS, SVJ_TJ, bs,cs,ds, jfactno2, jfacto1d, &
         DEP_ROOH, DEP_HO2, DEP_N2O5, DEP_NO3, DEP_ALL
    REAL(dp) time, temp, M, N2, O2, RO2, H2O, THETA, TIME2, LAT2
    REAL*8 y,dy,x,tmp(19), tmp2(19),b(19),c(19),d(19)
    integer i,n,jl
    INTEGER LK

 
    ! ************************************************************************
    ! define generic reaction rates.
    ! ************************************************************************

    ! constants used in calculation of reaction rates
    !M  = 2.55E+19
    N2 = 0.79*M
    O2 = 0.2095*M

    kp_rooh = DEP_ROOH
    kp_ho2  = DEP_HO2
    kp_n2o5 = DEP_N2O5
    kp_no3  = DEP_NO3
    kp_all  = DEP_ALL


    ! kro2no : ro2      + no      = ro      + no2
    !        : ro2      + no      = rono2
    ! mcm v3.2
    !kro2no    = 2.7d-12*EXP(360.0/temp)

    ! kro2ho2: ro2      + ho2     = rooh    + o2
    ! mcm protocol v3.2
    !kro2ho2   = 2.91d-13*EXP(1300.0/temp)

    ! kapho2 : rcoo2    + ho2     = products
    ! mcm protocol v3.2
    !kapho2    = 5.20d-13*EXP(980.0/temp)

    ! kapno  : rcoo2    + no      = products
    ! mcm v3.2
    !kapno = 7.50d-12*EXP(290.0/temp)

    ! kro2no3: ro2      + no3     = products
    ! mcm protocol v3.2
    !kro2no3   = 2.30d-12

    ! kno3al : no3      + rcho    = rcoo2   + hno3
    ! mcm protocol v3.2
    !kno3al    = 1.4d-12*EXP(-1860.0/temp)

    ! kdec   : ro                 = products
    ! mcm protocol v3.2
    !kdec      = 1.00d+06
    
    kalkoxy=6.00d-14*EXP(-550.0/temp)*o2
    kalkpxy=1.50d-14*EXP(-200.0/temp)*o2

    !kch3o2 = 1.03d-13*EXP(365.0/TEMP)
    !k298ch3o2 = 3.5d-13

    ! -------------------------------------------------------------------
    ! complex reactions
    ! -------------------------------------------------------------------

    ! kfpan kbpan
    ! formation and decomposition of pan
    ! iupac 2001 (mcmv3.2)
    kc0     = 2.70d-28*m*(temp/300.0)**(-7.1)
    kci     = 1.21d-11*(temp/300.0)**(-0.9)
    krc     = kc0/kci
    fcc     = 0.30
    nc      = 0.75-(1.27*dlog10(fcc))
    fc      = 10**(dlog10(fcc)/(1+((dlog10(krc))/nc)**2))
    kfpan   = (kc0*kci)*fc/(kc0+kci)

    kd0     = 4.90d-03*m*EXP(-12100.0/temp)
    kdi     = 5.40d+16*EXP(-13830.0/temp)
    krd     = kd0/kdi
    fcd     = 0.30
    ncd     = 0.75-(1.27*dlog10(fcd))
    fd      = 10**(dlog10(fcd)/(1+((dlog10(krd))/ncd)**2))
    kbpan   = (kd0*kdi)*fd/(kd0+kdi)

    ! kmt01  : o        + no      = no2
    ! iupac 2001 (mcmv3.2)
    k10     = 1.00d-31*m*(temp/300.0)**(-1.6)

    k1i     = 3.00d-11*(temp/300.0)**(0.3)
    kr1     = k10/k1i
    fc1     = 0.85
    nc1     = 0.75-(1.27*dlog10(fc1))
    f1      = 10**(dlog10(fc1)/(1+((dlog10(kr1)/nc1))**2))
    kmt01   = (k10*k1i)*f1/(k10+k1i)

    ! kmt02  : o        + no2     = no3
    ! iupac 2001 (mcmv3.2)
    k20     = 1.30d-31*m*(temp/300.0)**(-1.5)
    k2i     = 2.30d-11*(temp/300.0)**(0.24)
    kr2     = k20/k2i
    fc2     = 0.6
    nc2     = 0.75-(1.27*dlog10(fc2))
    f2      = 10**(dlog10(fc2)/(1+((dlog10(kr2)/nc2))**2))
    kmt02   = (k20*k2i)*f2/(k20+k2i)

    ! kmt03  : no2      + no3     = n2o5
    ! iupac 2006, mcmv3.2
    k30     = 3.60d-30*m*(temp/300.0)**(-4.1)
    k3i     = 1.90d-12*(temp/300.0)**(0.2)
    kr3     = k30/k3i
    fc3     = 0.35
    nc3     = 0.75-(1.27*dlog10(fc3))
    f3      = 10**(dlog10(fc3)/(1+((dlog10(kr3)/nc3))**2))
    kmt03   = (k30*k3i)*f3/(k30+k3i)

    ! kmt04  : n2o5               = no2     + no3
    ! iupac 2006, mcmv3.2
    k40     = 1.30d-03*m*(temp/300.0)**(-3.5)*EXP(-11000.0/temp)
    k4i     = 9.70d+14*(temp/300.0)**(0.1)*EXP(-11080.0/temp)
    kr4     = k40/k4i
    fc4     = 0.35
    nc4     = 0.75-(1.27*dlog10(fc4))
    f4      = 10**(dlog10(fc4)/(1+((dlog10(kr4)/nc4))**2))
    kmt04   = (k40*k4i)*f4/(k40+k4i)

    ! kmt05  : oh       + co(+o2) = ho2     + co2
    ! iupac 2006
    kmt05  = (1 + (m/4.2d19))

    ! kmt06  : ho2      + ho2     = h2o2    + o2
    ! water enhancement factor
    ! iupac 1992

    kmt06  = 1 + (1.40d-21*EXP(2200.0/temp)*h2o)

    ! kmt06  = 1 + (2.00d-25*EXP(4670.0/temp)*h2o)
    ! S+R 2005 values

    ! kmt07  : oh       + no      = hono

    ! iupac 2006, mcmv3.2
    k70     = 7.40d-31*m*(temp/300.0)**(-2.4)
    k7i     = 3.30d-11*(temp/300.0)**(-0.3)
    kr7     = k70/k7i
    fc7     = EXP(-temp/1420.0)
    nc7     = 0.75-(1.27*dlog10(fc7))
    f7      = 10**(dlog10(fc7)/(1+((dlog10(kr7)/nc7))**2))
    kmt07   = (k70*k7i)*f7/(k70+k7i)

    ! kmt08  : oh       + no2     = hno3

    ! iupac 2006, mcmv3.2
    k80     = 3.30d-30*m*(temp/300.0)**(-3.0)
    k8i     = 4.10d-11
    kr8     = k80/k8i
    fc8     = 0.4
    nc8     = 0.75-(1.27*dlog10(fc8))
    f8      = 10**(dlog10(fc8)/(1+((dlog10(kr8)/nc8))**2))
    kmt08   = (k80*k8i)*f8/(k80+k8i)

    ! kmt09  : ho2      + no2     = ho2no2
    ! iupac 1997, mcmv3.2

    k90     = 1.80d-31*m*(temp/300.0)**(-3.2)
    k9i     = 4.70d-12
    kr9     = k90/k9i
    fc9     = 0.6
    nc9     = 0.75-(1.27*dlog10(fc9))
    f9      = 10**(dlog10(fc9)/(1+((dlog10(kr9)/nc9))**2))
    kmt09   = (k90*k9i)*f9/(k90+k9i)

    ! kmt10  : ho2no2             = ho2     + no2
    ! iupac 1997, mcmv3.2

    k100     = 4.10d-05*m*EXP(-10650.0/temp)
    k10i     = 4.80d+15*EXP(-11170.0/temp)
    kr10     = k100/k10i
    fc10     = 0.6
    nc10     = 0.75-(1.27*dlog10(fc10))
    f10      = 10**(dlog10(fc10)/(1+((dlog10(kr10)/nc10))**2))
    kmt10    = (k100*k10i)*f10/(k100+k10i)

    ! kmt11  : oh       + hno3    = h2o     + no3
    ! iupac 2006, mcmv3.2

    k1     = 2.40d-14*EXP(460.0/temp)
    k3     = 6.50d-34*EXP(1335.0/temp)
    k4     = 2.70d-17*EXP(2199.0/temp)
    k2     = (k3*m)/(1+(k3*m/k4))
    kmt11  = k1 + k2

    ! kmt12 iupac 2006, mcmv3.2

    k120 = 4.50d-31*((temp/300.0)**(-3.9))*m
    k12i = 1.30d-12*((temp/300.0)**(-0.7))
    kr12 = k120/k12i
    fc12 = 0.525
    nc12 = 0.75-(1.27*dlog10(fc12))
    f12  = 10**(dlog10(fc12)/(1.+((dlog10(kr12)/nc12))**2))
    kmt12    = (k120*k12i)*f12/(k120+k12i)

    ! kmt13  : ch3o2    + no2     = ch3o2no2
    ! iupac 2006

    k130     = 2.50d-30*((temp/300.0)**(-5.5))*m
    k13i     = 1.80d-11
    kr13     = k130/k13i
    fc13     = 0.36
    nc13     = 0.75-(1.27*dlog10(fc13))
    f13      = 10**(dlog10(fc13)/(1+((dlog10(kr13)/nc13))**2))
    kmt13    = (k130*k13i)*f13/(k130+k13i)

    ! kmt14  : ch3o2no2           = ch3o2   + no2
    ! iupac 2006, mcmv3.2

    k140     = 9.00d-05*EXP(-9690.0/temp)*m
    k14i     = 1.10d+16*EXP(-10560.0/temp)
    kr14     = k140/k14i
    fc14     = 0.4
    nc14     = 0.75-(1.27*dlog10(fc14))
    f14      = 10**(dlog10(fc14)/(1+((dlog10(kr14)/nc14))**2))
    kmt14    = (k140*k14i)*f14/(k140+k14i)

    ! kmt15 iupac 2006, mcmv3.2

    k150 = 8.60d-29*((temp/300.0)**(-3.1))*m
    k15i = 9.00d-12*((temp/300.0)**(-0.85))
    kr15 = k150/k15i
    fc15 = 0.48
    nc15 = 0.75-(1.27*dlog10(fc15))
    f15  = 10**(dlog10(fc15)/(1+((dlog10(kr15)/nc15))**2))
    kmt15 = (k150*k15i)*f15/(k150+k15i)

    ! kmt16  :  oh  +  c3h6
    ! iupac 2006

    k160     = 8.00d-27*((temp/300.0)**(-3.5))*m
    k16i     = 3.00d-11*((temp/300.0)**(-1.0))
    kr16     = k160/k16i
    fc16     = 0.5
    nc16     = 0.75-(1.27*dlog10(fc16))
    f16      = 10**(dlog10(fc16)/(1+((dlog10(kr16)/nc16))**2))
    kmt16    = (k160*k16i)*f16/(k160+k16i)

    ! kmt17 iupac 2006

    k170 = 5.00d-30*((temp/300.0)**(-1.5))*m
    k17i = 1.00d-12
    kr17 = k170/k17i
    fc17 = (0.17*EXP(-51./temp))+EXP(-temp/204.)
    nc17 = 0.75-(1.27*dlog10(fc17))
    f17  = 10**(dlog10(fc17)/(1+((dlog10(kr17)/nc17))**2))
    kmt17 = (k170*k17i)*f17/(k170+k17i)

    ! kmt18 2011 oh + dms

    !kmt18=(9.5d-39*O2*EXP(5270./TEMP))/(1+7.5d-29*O2*exp(5610./TEMP))

    !       mcm v3.2

    !kroprim  = 2.50d-14*EXP(-300.0/temp)
    !krosec   = 2.50d-14*EXP(-300.0/temp)

  
  END SUBROUTINE mcm_constants

  !***************************************************************************

      
END MODULE constants
