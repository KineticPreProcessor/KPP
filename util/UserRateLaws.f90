!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  User-defined Rate Law functions
!
!  NOTE: For computational efficiency, we have created duplicate rate law
!  routines here that take either all single precision or all double precision
!  arguments.  Explicit casts to DBLE are skipped in the functions that take
!  all double precision arguments (as this removes unneeded computations).
!
!  These functions are overloaded by INTERFACE statements, which are located
!  in file UserRateLawsInterfaces.f90.  The UserRateLawsInterfaces.f90 file
!  will be in-lined into the top of the KPP_ROOT_Rates module.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  FUNCTION ARR_abc_dp( a0, b0, c0 ) RESULT( k )
    ! Arrhenius function, for nonzero a0, b0, c0 (dp args)
    REAL(dp), INTENT(IN) :: a0, b0, c0
    KPP_REAL :: k
    k = a0 * EXP(-b0/TEMP) * (TEMP/300.0_dp)**C0
  END FUNCTION ARR_abc_dp

  FUNCTION ARR_abc_sp( a0, b0, c0 ) RESULT( k )
    ! Arrhenius function, for nonzero a0, b0, c0 (sp args)
    REAL(sp), INTENT(IN) :: a0, b0, c0
    KPP_REAL :: k
    k = DBLE(a0) * EXP(-DBLE(b0)/TEMP) * (TEMP/300.0_dp)**DBLE(c0)
  END FUNCTION ARR_abc_sp

  !---------------------------------------------------------------------------

  FUNCTION ARR_ab_dp( a0, b0 ) RESULT( k )
    ! Arrhenius function, for nonzero a0, b0 (dp args)
    REAL(dp), INTENT(IN) :: a0, b0
    KPP_REAL :: k
    k = a0 * EXP(-b0/TEMP)
  END FUNCTION ARR_ab_dp

  FUNCTION ARR_ab_sp( a0, b0 ) RESULT( k )
    ! Arrhenius function, for nonzero a0, b0 (sp args)
    REAL(sp), INTENT(IN) :: a0, b0
    KPP_REAL :: k
    k = DBLE(a0) * EXP(-DBLE(b0)/TEMP)
  END FUNCTION ARR_ab_sp

  !---------------------------------------------------------------------------

  FUNCTION ARR_ac_dp( a0, c0 ) RESULT( k )
    ! Arrhenius function, for nonzero a0, c0 (dp args)
    REAL(dp), INTENT(IN) :: a0, c0
    KPP_REAL :: k
    k = a0 * (TEMP/300.0_dp)**C0
  END FUNCTION ARR_ac_dp

  FUNCTION ARR_ac_sp( a0, c0 ) RESULT( k )
    ! Arrhenius function, for nonzero a0, c0 (sp args)
    REAL(sp), INTENT(IN) :: a0, c0
    KPP_REAL :: k
    k = DBLE(a0) * (TEMP/300.0_dp)**DBLE(c0)
  END FUNCTION ARR_ac_sp

  !---------------------------------------------------------------------------

  FUNCTION EP2_dp( a0, c0, a2, c2, a3, c3 ) RESULT( k )
    ! EP2 function, for saprc99 and saprcnov (dp args)
    REAL(dp), INTENT(IN) :: a0, c0, a2, c2, a3, c3
    REAL(dp) :: k0, k2, k3
    KPP_REAL :: k
    k0 = a0 * EXP(-c0/temp)
    k2 = a2 * EXP(-c2/temp)
    k3 = a3 * EXP(-c3/temp)
    k3 = k3 * CFACTOR * 1.0E6_dp
    k  = k0 + k3/(1.0_dp + k3/k2)
  END FUNCTION EP2_dp

  FUNCTION EP2_sp( a0, c0, a2, c2, a3, c3 ) RESULT( k )
    ! EP2 function, for saprc99 and saprcnov (sp args)
    REAL(sp), INTENT(IN) :: a0, c0, a2, c2, a3, c3
    REAL(dp) :: k0, k2, k3
    KPP_REAL :: k
    k0 = DBLE(a0) * EXP(-DBLE(c0)/TEMP)
    k2 = DBLE(a2) * EXP(-DBLE(c2)/TEMP)
    k3 = DBLE(a3) * EXP(-DBLE(c3)/TEMP)
    k3 = k3 * CFACTOR * 1.0E6_dp
    k  = k0 + K3/(1.0_dp + k3/k2)
  END FUNCTION EP2_sp

  !---------------------------------------------------------------------------

  FUNCTION EP3_dp( a1, c1, a2, c2) RESULT( k )
    ! EP3 function, for saprc99 and saprcnov (dp args)
    REAL(dp), INTENT(IN) :: a1, c1, a2, c2
    REAL(dp) :: k1, k2
    KPP_REAL :: k
    k1 = a1 * EXP(-c1/TEMP)
    k2 = a2 * EXP(-c2/TEMP)
    k  = k1 + k2*(1.0E6_dp * CFACTOR)
  END FUNCTION EP3_dp

  FUNCTION EP3_sp( a1, c1, a2, c2 ) RESULT( k )
    ! EP3 function, for saprc99 and saprcnov (sp args)
    REAL(sp), INTENT(IN) :: a1, c1, a2, c2
    REAL(dp) :: k1, k2
    KPP_REAL :: k
    k1 = DBLE(a1) * EXP(-DBLE(c1)/TEMP)
    k2 = DBLE(a2) * EXP(-DBLE(c2)/TEMP)
    k  = k1 + k2*(1.0E6_dp * CFACTOR)
  END FUNCTION EP3_sp

  !---------------------------------------------------------------------------

  FUNCTION FALL_dp( a0, b0, c0, a1, b1, c1, cf ) RESULT( k )
    ! FALL function, for saprc99 and saprcnov (dp args)
    REAL(dp), INTENT(IN) :: a0, b0, c0, a1, b1, c1, cf
    REAL(dp) :: k0, k1
    KPP_REAL :: k
    k0 = a0 * EXP(-b0/TEMP) * (TEMP/300.0_dp)**c0
    k1 = a1 * EXP(-b1/TEMP) * (TEMP/300.0_dp)**c1
    k0 = k0 * CFACTOR * 1.0E6_dp
    k1 = k0 / k1
    k  = (k0/(1.0_dp+k1)) * cf**(1.0_dp/(1.0_dp+(LOG10(k1))**2))
  END FUNCTION FALL_dp

  FUNCTION FALL_sp( a0, b0, c0, a1, b1, c1, cf ) RESULT( k )
    ! FALL function, for saprc99 and saprcnov (sp args)
    REAL(sp), INTENT(IN) :: a0, b0, c0, a1, b1, c1, cf
    REAL(dp) :: k0, k1
    KPP_REAL :: k
    k0 = DBLE(A0) * EXP(-DBLE(B0)/TEMP)* (TEMP/300.0_dp)**DBLE(C0)
    k1 = DBLE(A1) * EXP(-DBLE(B1)/TEMP)* (TEMP/300.0_dp)**DBLE(C1)
    k0 = k0 * CFACTOR * 1.0E6_dp
    k1 = k0 / k1
    k  = (k0/(1.0_dp+k1)) * DBLE(cf)**(1.0_dp/(1.0_dp+(LOG10(k1))**2))
  END FUNCTION FALL_sp

  !---------------------------------------------------------------------------

  ! JPL (jpldataeval.jpl.nasa.gov) three-body reaction formula:
  
  REAL(dp) FUNCTION k3rd_jpl_dp(cair,k0_300K,n,kinf_300K,m,fc) ! dp args
    INTRINSIC LOG10
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(dp), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(dp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(dp), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(dp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(dp), INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
    REAL(dp) :: zt_help, k0_T, kinf_T, k_ratio
    zt_help  = 300._dp/temp
    k0_T     = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T   = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio  = k0_T/kinf_T
    k3rd_jpl_dp = k0_T/(1._dp+k_ratio)*fc**(1._dp/(1._dp+LOG10(k_ratio)**2))
  END FUNCTION k3rd_jpl_dp

  KPP_REAL FUNCTION k3rd_jpl_sp(cair,k0_300K,n,kinf_300K,m,fc) ! sp args
    INTRINSIC LOG10
    REAL(sp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(sp), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(sp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(sp), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(sp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(sp), INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
    REAL(sp) :: zt_help, k0_T, kinf_T, k_ratio
    zt_help  = 300._dp/temp
    k0_T     = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T   = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio  = k0_T/kinf_T
    k3rd_jpl_sp = k0_T/(1._dp+k_ratio)*fc**(1._dp/(1._dp+LOG10(k_ratio)**2))
  END FUNCTION k3rd_jpl_sp

  ! --------------------------------------------------------------------------

  ! JPL (jpldataeval.jpl.nasa.gov) termolecular chemical activation reaction:

  FUNCTION k3rd_jpl_activation_dp(cair,k0_298K,n,kinf_298K,m,A,B) ! dp args
    INTRINSIC :: LOG10
    REAL(dp), DIMENSION(2) :: k3rd_jpl_activation_dp
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(dp), INTENT(IN) :: k0_298K   ! low pressure limit at 300 K
    REAL(dp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(dp), INTENT(IN) :: kinf_298K ! high pressure limit at 300 K
    REAL(dp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(dp), INTENT(IN) :: A         ! for k_int
    REAL(dp), INTENT(IN) :: B         ! for k_int
    REAL(dp)             :: zt_help, k0_TM, kinf_T, k_ratio, k_f, k_int, k_fCA
    zt_help = 298./temp
    k0_TM   = k0_298K   * zt_help**n * cair ! k_0   at current T * M
    kinf_T  = kinf_298K * zt_help**m        ! k_inf at current T
    k_ratio = k0_TM/kinf_T
    k_f     = k0_TM/(1.+k_ratio)*0.6**(1./(1.+LOG10(k_ratio)**2))
    k_int   = A * exp(-B/temp)
    k_fCA   = k_int * (1. - k_f/kinf_T)
    k3rd_jpl_activation_dp(ASSOC)  = k_f
    k3rd_jpl_activation_dp(DISSOC) = k_fCA
  END FUNCTION k3rd_jpl_activation_dp

  FUNCTION k3rd_jpl_activation_sp(cair,k0_298K,n,kinf_298K,m,A,B) ! sp args
    INTRINSIC :: LOG10
    KPP_REAL, DIMENSION(2) :: k3rd_jpl_activation_sp
    REAL(sp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(sp), INTENT(IN) :: k0_298K   ! low pressure limit at 300 K
    REAL(sp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(sp), INTENT(IN) :: kinf_298K ! high pressure limit at 300 K
    REAL(sp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(sp), INTENT(IN) :: A         ! for k_int
    REAL(sp), INTENT(IN) :: B         ! for k_int
    REAL(sp)             :: zt_help, k0_TM, kinf_T, k_ratio, k_f, k_int, k_fCA
    zt_help = 298./temp
    k0_TM   = k0_298K   * zt_help**n * cair ! k_0   at current T * M
    kinf_T  = kinf_298K * zt_help**m        ! k_inf at current T
    k_ratio = k0_TM/kinf_T
    k_f     = k0_TM/(1.+k_ratio)*0.6**(1./(1.+LOG10(k_ratio)**2))
    k_int   = A * exp(-B/temp)
    k_fCA   = k_int * (1. - k_f/kinf_T)
    k3rd_jpl_activation_sp(ASSOC)  = k_f
    k3rd_jpl_activation_sp(DISSOC) = k_fCA
  END FUNCTION k3rd_jpl_activation_sp

  ! --------------------------------------------------------------------------

  ! IUPAC (iupac.pole-ether.fr) three-body reaction formula:

  REAL(dp) FUNCTION k3rd_iupac_dp(cair,k0_300K,n,kinf_300K,m,fc) ! dp args
    INTRINSIC :: LOG10
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(dp), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(dp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(dp), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(dp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(dp), INTENT(IN) :: fc        ! broadening factor (e.g. 0.45 or 0.6...)
    REAL(dp)             :: nu        ! N
    REAL(dp)             :: zt_help, k0_T, kinf_T, k_ratio
    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    nu      = 0.75-1.27*LOG10(fc)
    k3rd_iupac_dp = k0_T/(1._dp+k_ratio)* &
      fc**(1._dp/(1._dp+(LOG10(k_ratio)/nu)**2))
  END FUNCTION k3rd_iupac_dp

  KPP_REAL FUNCTION k3rd_iupac_sp(cair,k0_300K,n,kinf_300K,m,fc) ! sp args
    INTRINSIC :: LOG10
    REAL(sp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL(sp), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL(sp), INTENT(IN) :: n         ! exponent for low pressure limit
    REAL(sp), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL(sp), INTENT(IN) :: m         ! exponent for high pressure limit
    REAL(sp), INTENT(IN) :: fc        ! broadening factor (e.g. 0.45 or 0.6...)
    REAL(sp)             :: nu        ! N
    REAL(sp)             :: zt_help, k0_T, kinf_T, k_ratio
    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    nu      = 0.75-1.27*LOG10(fc)
    k3rd_iupac_sp = k0_T/(1._dp+k_ratio)* &
      fc**(1._dp/(1._dp+(LOG10(k_ratio)/nu)**2))
  END FUNCTION k3rd_iupac_sp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  End of User-defined Rate Law functions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
