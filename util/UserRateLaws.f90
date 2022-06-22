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
    ! EP3 function, for saprc99 and saprcnov (sp args)
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  End of User-defined Rate Law functions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
