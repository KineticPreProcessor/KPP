! rate constants and functions for the MCM-generated KPP equation file

MODULE constants_mcm

  USE mcm_Precision, ONLY: dp
  USE mcm_Parameters ! ind_*
  USE mcm_Global, ONLY: C, TEMP
  IMPLICIT NONE

  INTEGER, PARAMETER :: J_O3_O1D       =  1 ! O3       -> O(1D) +O2
  INTEGER, PARAMETER :: J_O3_O3P       =  2 ! O3       -> O(3P) +O2
  INTEGER, PARAMETER :: J_H2O2         =  3 ! H2O2     -> OH + OH
  INTEGER, PARAMETER :: J_NO2          =  4 ! NO2      -> NO + O(3P)
  INTEGER, PARAMETER :: J_NO3_NO       =  5 ! NO3      -> NO + O2
  INTEGER, PARAMETER :: J_NO3_NO2      =  6 ! NO3      -> NO2 + O(3P)
  INTEGER, PARAMETER :: J_HONO         =  7 ! HONO     -> NO + OH
  INTEGER, PARAMETER :: J_HNO3         =  8 ! HNO3     -> NO2 + OH
  INTEGER, PARAMETER :: J_HCHO_H       = 11 ! HCHO     -> H + HCO
  INTEGER, PARAMETER :: J_HCHO_H2      = 12 ! HCHO     -> H2 + CO
  INTEGER, PARAMETER :: J_CH3CHO       = 13 ! CH3CHO   -> CH3 + HCO
  INTEGER, PARAMETER :: J_C2H5CHO      = 14 ! C2H5CHO  -> C2H5 + HCO
  INTEGER, PARAMETER :: J_C3H7CHO_HCO  = 15 ! C3H7CHO  -> n-C3H7 + HCO
  INTEGER, PARAMETER :: J_C3H7CHO_C2H4 = 16 ! C3H7CHO  -> C2H4 + CH3CHO
  INTEGER, PARAMETER :: J_IPRCHO       = 17 ! IPRCHO   -> n-C4H9 + HCO
  INTEGER, PARAMETER :: J_MACR_HCO     = 18 ! MACR     -> CH2=CCH3 + HCO
  INTEGER, PARAMETER :: J_MACR_H       = 19 ! MACR     -> CH2=C(CH3)CO + H
  INTEGER, PARAMETER :: J_C5HPALD1     = 20 ! C5HPALD1 -> CH3C(CHO)=CHCH2O + OH
  INTEGER, PARAMETER :: J_CH3COCH3     = 21 ! CH3COCH3 -> CH3CO + CH3
  INTEGER, PARAMETER :: J_MEK          = 22 ! MEK      -> CH3CO + C2H5
  INTEGER, PARAMETER :: J_MVK_CO       = 23 ! MVK      -> CH3CH=CH2 + CO
  INTEGER, PARAMETER :: J_MVK_C2H3     = 24 ! MVK      -> CH3CO + CH2=CH
  INTEGER, PARAMETER :: J_GLYOX_H2     = 31 ! GLYOX    -> CO + CO + H2
  INTEGER, PARAMETER :: J_GLYOX_HCHO   = 32 ! GLYOX    -> HCHO + CO
  INTEGER, PARAMETER :: J_GLYOX_HCO    = 33 ! GLYOX    -> HCO + HCO
  INTEGER, PARAMETER :: J_MGLYOX       = 34 ! MGLYOX   -> CH3CO + HCO
  INTEGER, PARAMETER :: J_BIACET       = 35 ! BIACET   -> CH3CO + CH3CO
  INTEGER, PARAMETER :: J_CH3OOH       = 41 ! CH3OOH   -> CH3O + OH
  INTEGER, PARAMETER :: J_CH3NO3       = 51 ! CH3NO3   -> CH3O + NO2
  INTEGER, PARAMETER :: J_C2H5NO3      = 52 ! C2H5NO3  -> C2H5O + NO2
  INTEGER, PARAMETER :: J_NC3H7NO3     = 53 ! NC3H7NO3 -> n-C3H7O + NO2
  INTEGER, PARAMETER :: J_IC3H7NO3     = 54 ! IC3H7NO3 -> CH3C(O.)CH3 + NO2
  INTEGER, PARAMETER :: J_TC4H9NO3     = 55 ! TC4H9NO3 -> t-C4H9O + NO2
  INTEGER, PARAMETER :: J_NOA          = 56 ! NOA      -> CH3C(O)CH2(O.) + NO2 or CH3CO + HCHO + NO2

  REAL(dp) :: K14ISOM1, K298CH3O2, KAPHO2, KAPNO, KCH3O2, KDEC, KMT05, KMT06, KNO3AL, &
    KRO2HO2, KRO2NO, KRO2NO3, KROPRIM, KROSEC, FCD, KD0, KDI, KRD, NCD, FD, &
    KBPAN, FCC, KC0, KCI, KRC, NC, FC, KFPAN, FC1, K10, K1I, KR1, NC1, F1, &
    KMT01, FC2, K20, K2I, KR2, NC2, F2, KMT02, FC3, K30, K3I, KR3, NC3, F3, &
    KMT03, FC4, K40, K4I, KR4, NC4, F4, KMT04, FC7, K70, K7I, KR7, NC7, F7, &
    KMT07, FC8, K80, K8I, KR8, NC8, F8, KMT08, FC9, K90, K9I, KR9, NC9, F9, &
    KMT09, FC10, K100, K10I, KR10, NC10, F10, KMT10, K3, K4, K1, K2, KMT11, &
    FC12, K120, K12I, KR12, NC12, F12, KMT12, FC13, K130, K13I, KR13, NC13, &
    F13, KMT13, FC14, K140, K14I, KR14, NC14, F14, KMT14, FC15, K150, K15I, &
    KR15, NC15, F15, KMT15, FC16, K160, K16I, KR16, NC16, F16, KMT16
  REAL(dp), DIMENSION(56) :: J
  REAL(dp) :: M, N2, O2, RO2, H2O, zenith

  PUBLIC

CONTAINS

  SUBROUTINE define_constants_mcm()

    RO2 = C(ind_CH3O2) + C(ind_C51O2) + C(ind_CH3COCH2O2) + C(ind_HOCH2CO3) + &
      C(ind_CO2C3CO3) + C(ind_IPRHOCO3) + C(ind_BIACETO2) + C(ind_CH3CO3) + &
      C(ind_HOCH2CH2O2) + C(ind_ETHENO3O2) + C(ind_NO3CH2CO3) + &
      C(ind_HYPROPO2) + C(ind_IPROPOLO2) + C(ind_PRONO3AO2) + &
      C(ind_PRONO3BO2) + C(ind_PRNO3CO3) + C(ind_H13CO2CO3) + &
      C(ind_HCOCH2O2) + C(ind_ACO3) + C(ind_ISOPAO2) + C(ind_HC4CCO3) + &
      C(ind_C57O2) + C(ind_HCOCO3) + C(ind_ISOPBO2) + C(ind_HMVKAO2) + &
      C(ind_HMVKBO2) + C(ind_CO2H3CO3) + C(ind_ISOPCO2) + C(ind_HC4ACO3) + &
      C(ind_C58O2) + C(ind_ISOPDO2) + C(ind_C59O2) + C(ind_MACO3) + &
      C(ind_CH3C2H2O2) + C(ind_MACRO2) + C(ind_NISOPO2) + C(ind_NC4CO3) + &
      C(ind_C510O2) + C(ind_MVKO2) + C(ind_CH3CHOHCO3) + C(ind_C3DIOLO2) + &
      C(ind_MCOCOMOXO2) + C(ind_MVKOHAO2) + C(ind_MVKOHBO2) + &
      C(ind_HCOCOHCO3) + C(ind_MMALANHYO2) + C(ind_C3MCODBCO3) + &
      C(ind_MC3CODBCO3) + C(ind_C4M2ALOHO2) + C(ind_CH3COCO3) + &
      C(ind_CHOMOHCO3) + C(ind_C3MDIALO2) + C(ind_C4CO2O2) + C(ind_PXYFUO2) + &
      C(ind_C23O3CCO3) + C(ind_NPXYFUO2) + C(ind_A2PANOO) + C(ind_OCCOHCO2) + &
      C(ind_ACRO2) + C(ind_MACROHO2) + C(ind_IEACO3) + C(ind_IECCO3) + &
      C(ind_INAO2) + C(ind_INB1O2) + C(ind_INB2O2) + C(ind_INCO2) + &
      C(ind_INDO2) + C(ind_C58NO3CO3) + C(ind_INANCO3) + C(ind_INAHPCO3) + &
      C(ind_INAHCO3) + C(ind_INANCOCO3) + C(ind_MACRNCO3) + C(ind_CONM2CO3) + &
      C(ind_C57NO3CO3) + C(ind_MMALNACO3) + C(ind_MMALNBCO3) + &
      C(ind_INCNCO3) + C(ind_MACRNBCO3) + C(ind_COHM2CO3) + C(ind_INDHPCO3) + &
      C(ind_INDHCO3) + C(ind_CO2N3CO3) + C(ind_INB1NACO3) + &
      C(ind_INB1NBCO3) + C(ind_INB1HPCO3) + C(ind_C524O2) + C(ind_C525O2) + &
      C(ind_NC524O2) + C(ind_HIEB1O2) + C(ind_HIEB2O2) + C(ind_HMACO3) + &
      C(ind_HMACRO2) + C(ind_ISOP34O2) + C(ind_CISOPAO2) + C(ind_CISOPCO2) + &
      C(ind_C526O2) + C(ind_C527O2) + C(ind_C534O2) + C(ind_HPC52O2) + &
      C(ind_C535O2) + C(ind_HPC52CO3) + C(ind_CH2CHCH2O2) + C(ind_NC51O2) + &
      C(ind_C45O2) + C(ind_M3BU3ECO3) + C(ind_C530O2) + C(ind_NC526O2) + &
      C(ind_HC4CO3) + C(ind_C531O2) + C(ind_C31CO3) + C(ind_C533O2) + &
      C(ind_C58AO2) + C(ind_C47CO3) + C(ind_C57AO2) + C(ind_C536O2) + &
      C(ind_C537O2)
    
    K14ISOM1 = 3.00E7*EXP(-5300./TEMP)
    K298CH3O2 = 3.5E-13
    KAPHO2 = 5.2E-13*EXP(980./TEMP)
    KAPNO = 7.5E-12*EXP(290./TEMP)
    KCH3O2 = 1.03E-13*EXP(365./TEMP)
    KDEC = 1.00E+06
    KMT05 = 1.44E-13*(1.+(M/4.2E+19))
    KMT06 = 1. + (1.40E-21*EXP(2200./TEMP)*H2O)
    KNO3AL = 1.44E-12*EXP(-1862./TEMP)
    KRO2HO2 = 2.91E-13*EXP(1300./TEMP)
    KRO2NO = 2.7E-12*EXP(360./TEMP)
    KRO2NO3 = 2.3E-12
    KROPRIM = 2.50E-14*EXP(-300./TEMP)
    KROSEC = 2.50E-14*EXP(-300./TEMP)
    FCD = 0.30
    KD0 = 1.10E-05*M*EXP(-10100./TEMP)
    KDI = 1.90E17*EXP(-14100./TEMP)
    KRD = KD0/KDI
    NCD = 0.75-1.27*(LOG10(FCD))
    FD = 10.**(LOG10(FCD)/(1.+(LOG10(KRD)/NCD)**(2.)))
    KBPAN = (KD0*KDI)*FD/(KD0+KDI)
    FCC = 0.30
    KC0 = 3.28E-28*M*(TEMP/300.)**(-6.87)
    KCI = 1.125E-11*(TEMP/300.)**(-1.105)
    KRC = KC0/KCI
    NC = 0.75-1.27*(LOG10(FCC))
    FC = 10.**(LOG10(FCC)/(1.+(LOG10(KRC)/NC)**(2.)))
    KFPAN = (KC0*KCI)*FC/(KC0+KCI)
    FC1 = 0.85
    K10 = 1.0E-31*M*(TEMP/300.)**(-1.6)
    K1I = 5.0E-11*(TEMP/300.)**(-0.3)
    KR1 = K10/K1I
    NC1 = 0.75-1.27*(LOG10(FC1))
    F1 = 10.**(LOG10(FC1)/(1.+(LOG10(KR1)/NC1)**(2.)))
    KMT01 = (K10*K1I)*F1/(K10+K1I)
    FC2 = 0.6
    K20 = 1.3E-31*M*(TEMP/300.)**(-1.5)
    K2I = 2.3E-11*(TEMP/300.)**(0.24)
    KR2 = K20/K2I
    NC2 = 0.75-1.27*(LOG10(FC2))
    F2 = 10.**(LOG10(FC2)/(1.+(LOG10(KR2)/NC2)**(2.)))
    KMT02 = (K20*K2I)*F2/(K20+K2I)
    FC3 = 0.35
    K30 = 3.6E-30*M*(TEMP/300.)**(-4.1)
    K3I = 1.9E-12*(TEMP/300.)**(0.2)
    KR3 = K30/K3I
    NC3 = 0.75-1.27*(LOG10(FC3))
    F3 = 10.**(LOG10(FC3)/(1.+(LOG10(KR3)/NC3)**(2.)))
    KMT03 = (K30*K3I)*F3/(K30+K3I)
    FC4 = 0.35
    K40 = 1.3E-3*M*(TEMP/300.)**(-3.5)*EXP(-11000./TEMP)
    K4I = 9.7E+14*(TEMP/300.)**(0.1)*EXP(-11080./TEMP)
    KR4 = K40/K4I
    NC4 = 0.75-1.27*(LOG10(FC4))
    F4 = 10.**(LOG10(FC4)/(1.+(LOG10(KR4)/NC4)**(2.)))
    KMT04 = (K40*K4I)*F4/(K40+K4I)
    FC7 = 0.81
    K70 = 7.4E-31*M*(TEMP/300.)**(-2.4)
    K7I = 3.3E-11*(TEMP/300.)**(-0.3)
    KR7 = K70/K7I
    NC7 = 0.75-1.27*(LOG10(FC7))
    F7 = 10.**(LOG10(FC7)/(1.+(LOG10(KR7)/NC7)**(2.)))
    KMT07 = (K70*K7I)*F7/(K70+K7I)
    FC8 = 0.41
    K80 = 3.2E-30*M*(TEMP/300.)**(-4.5)
    K8I = 3.0E-11
    KR8 = K80/K8I
    NC8 = 0.75-1.27*(LOG10(FC8))
    F8 = 10.**(LOG10(FC8)/(1.+(LOG10(KR8)/NC8)**(2.)))
    KMT08 = (K80*K8I)*F8/(K80+K8I)
    FC9 = 0.4
    K90 = 1.4E-31*M*(TEMP/300.)**(-3.1)
    K9I = 4.0E-12
    KR9 = K90/K9I
    NC9 = 0.75-1.27*(LOG10(FC9))
    F9 = 10.**(LOG10(FC9)/(1.+(LOG10(KR9)/NC9)**(2.)))
    KMT09 = (K90*K9I)*F9/(K90+K9I)
    FC10 = 0.4
    K100 = 4.10E-05*M*EXP(-10650./TEMP)
    K10I = 6.0E+15*EXP(-11170./TEMP)
    KR10 = K100/K10I
    NC10 = 0.75-1.27*(LOG10(FC10))
    F10 = 10.**(LOG10(FC10)/(1.+(LOG10(KR10)/NC10)**(2.)))
    KMT10 = (K100*K10I)*F10/(K100+K10I)
    K3 = 6.50E-34*EXP(1335./TEMP)
    K4 = 2.70E-17*EXP(2199./TEMP)
    K1 = 2.40E-14*EXP(460./TEMP)
    K2 = (K3*M)/(1.+(K3*M/K4))
    KMT11 = K1 + K2
    FC12 = 0.53
    K120 = 2.5E-31*M*(TEMP/300.)**(-2.6)
    K12I = 2.0E-12
    KR12 = K120/K12I
    NC12 = 0.75-1.27*(LOG10(FC12))
    F12 = 10.**(LOG10(FC12)/(1.0+(LOG10(KR12)/NC12)**(2.)))
    KMT12 = (K120*K12I*F12)/(K120+K12I)
    FC13 = 0.36
    K130 = 2.5E-30*M*(TEMP/300.)**(-5.5)
    K13I = 1.8E-11
    KR13 = K130/K13I
    NC13 = 0.75-1.27*(LOG10(FC13))
    F13 = 10.**(LOG10(FC13)/(1.+(LOG10(KR13)/NC13)**(2.)))
    KMT13 = (K130*K13I)*F13/(K130+K13I)
    FC14 = 0.36
    K140 = 9.0E-5*EXP(-9690./TEMP)*M
    K14I = 1.1E+16*EXP(-10560./TEMP)
    KR14 = K140/K14I
    NC14 = 0.75-1.27*(LOG10(FC14))
    F14 = 10.**(LOG10(FC14)/(1.+(LOG10(KR14)/NC14)**(2.)))
    KMT14 = (K140*K14I)*F14/(K140+K14I)
    FC15 = 0.48
    K150 = 8.6E-29*M*(TEMP/300.)**(-3.1)
    K15I = 9.0E-12*(TEMP/300.)**(-0.85)
    KR15 = K150/K15I
    NC15 = 0.75-1.27*(LOG10(FC15))
    F15 = 10.**(LOG10(FC15)/(1.+(LOG10(KR15)/NC15)**(2.)))
    KMT15 = (K150*K15I)*F15/(K150+K15I)
    FC16 = 0.5
    K160 = 8.E-27*M*(TEMP/300.)**(-3.5)
    K16I = 3.0E-11*(TEMP/300.)**(-1.)
    KR16 = K160/K16I
    NC16 = 0.75-1.27*(LOG10(FC16))
    F16 = 10.**(LOG10(FC16)/(1.+(LOG10(KR16)/NC16)**(2.)))
    KMT16 = (K160*K16I)*F16/(K160+K16I)
    
    J(1) = 6.073E-05*(cos(zenith)**1.743)*exp(-0.474*(1./cos(zenith)))
    J(2) = 4.775E-04*(cos(zenith)**0.298)*exp(-0.08*(1./cos(zenith)))
    J(3) = 1.041E-05*(cos(zenith)**0.723)*exp(-0.279*(1./cos(zenith)))
    J(4) = 1.165E-02*(cos(zenith)**0.244)*exp(-0.267*(1./cos(zenith)))
    J(5) = 2.485E-02*(cos(zenith)**0.168)*exp(-0.108*(1./cos(zenith)))
    J(6) = 1.747E-01*(cos(zenith)**0.155)*exp(-0.125*(1./cos(zenith)))
    J(7) = 2.644E-03*(cos(zenith)**0.261)*exp(-0.288*(1./cos(zenith)))
    J(8) = 9.312E-07*(cos(zenith)**1.23)*exp(-0.307*(1./cos(zenith)))
    J(11) = 4.642E-05*(cos(zenith)**0.762)*exp(-0.353*(1./cos(zenith)))
    J(12) = 6.853E-05*(cos(zenith)**0.477)*exp(-0.323*(1./cos(zenith)))
    J(13) = 7.344E-06*(cos(zenith)**1.202)*exp(-0.417*(1./cos(zenith)))
    J(15) = 2.792E-05*(cos(zenith)**0.805)*exp(-0.338*(1./cos(zenith)))
    J(17) = 7.914E-05*(cos(zenith)**0.764)*exp(-0.364*(1./cos(zenith)))
    J(18) = 1.482E-06*(cos(zenith)**0.396)*exp(-0.298*(1./cos(zenith)))
    J(19) = 1.482E-06*(cos(zenith)**0.396)*exp(-0.298*(1./cos(zenith)))
    J(20) = 7.600E-04*(cos(zenith)**0.396)*exp(-0.298*(1./cos(zenith)))
    J(21) = 7.992E-07*(cos(zenith)**1.578)*exp(-0.271*(1./cos(zenith)))
    J(22) = 5.804E-06*(cos(zenith)**1.092)*exp(-0.377*(1./cos(zenith)))
    J(23) = 2.4246E-06*(cos(zenith)**0.395)*exp(-0.296*(1./cos(zenith)))
    J(24) = 2.424E-06*(cos(zenith)**0.395)*exp(-0.296*(1./cos(zenith)))
    J(31) = 6.845E-05*(cos(zenith)**0.13)*exp(-0.201*(1./cos(zenith)))
    J(32) = 1.032E-05*(cos(zenith)**0.13)*exp(-0.201*(1./cos(zenith)))
    J(33) = 3.802E-05*(cos(zenith)**0.644)*exp(-0.312*(1./cos(zenith)))
    J(34) = 1.537E-04*(cos(zenith)**0.17)*exp(-0.208*(1./cos(zenith)))
    J(35) = 3.326E-04*(cos(zenith)**0.148)*exp(-0.215*(1./cos(zenith)))
    J(41) = 7.649E-06*(cos(zenith)**0.682)*exp(-0.279*(1./cos(zenith)))
    J(51) = 1.588E-06*(cos(zenith)**1.154)*exp(-0.318*(1./cos(zenith)))
    J(53) = 2.485E-06*(cos(zenith)**1.196)*exp(-0.328*(1./cos(zenith)))
    J(54) = 4.095E-06*(cos(zenith)**1.111)*exp(-0.316*(1./cos(zenith)))
    J(55) = 1.135E-05*(cos(zenith)**0.974)*exp(-0.309*(1./cos(zenith)))
    J(56) = 4.365E-05*(cos(zenith)**1.089)*exp(-0.323*(1./cos(zenith)))

  END SUBROUTINE define_constants_mcm

END MODULE constants_mcm
