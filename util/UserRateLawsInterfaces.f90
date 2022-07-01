  PRIVATE :: ARR_abc_dp, ARR_abc_sp
  INTERFACE ARR_abc
     MODULE PROCEDURE ARR_abc_dp
     MODULE PROCEDURE ARR_abc_sp
  END INTERFACE ARR_abc

  PRIVATE :: ARR_ab_dp, ARR_ab_sp
  INTERFACE ARR_ab
     MODULE PROCEDURE ARR_ab_dp
     MODULE PROCEDURE ARR_ab_sp
  END INTERFACE ARR_ab

  PRIVATE :: ARR_ac_dp, ARR_ac_sp
  INTERFACE ARR_ac
     MODULE PROCEDURE ARR_ac_dp
     MODULE PROCEDURE ARR_ac_sp
  END INTERFACE ARR_ac

  PRIVATE :: EP2_dp, EP2_sp
  INTERFACE EP2
     MODULE PROCEDURE EP2_dp
     MODULE PROCEDURE EP2_sp
  END INTERFACE EP2

  PRIVATE :: EP3_dp, EP3_sp
  INTERFACE EP3
     MODULE PROCEDURE EP3_dp
     MODULE PROCEDURE EP3_sp
  END INTERFACE EP3

  PRIVATE :: FALL_dp, FALL_sp
  INTERFACE FALL
     MODULE PROCEDURE FALL_dp
     MODULE PROCEDURE FALL_sp
  END INTERFACE FALL

  PRIVATE :: k3rd_jpl_dp, k3rd_jpl_sp
  INTERFACE k3rd_jpl
     MODULE PROCEDURE k3rd_jpl_dp
     MODULE PROCEDURE k3rd_jpl_sp
  END INTERFACE k3rd_jpl

  PRIVATE :: k3rd_jpl_activation_dp, k3rd_jpl_activation_sp
  INTERFACE  k3rd_jpl_activation
     MODULE PROCEDURE  k3rd_jpl_activation_dp
     MODULE PROCEDURE  k3rd_jpl_activation_sp
  END INTERFACE  k3rd_jpl_activation

  PRIVATE :: k3rd_iupac_dp, k3rd_iupac_sp
  INTERFACE k3rd_iupac
     MODULE PROCEDURE k3rd_iupac_dp
     MODULE PROCEDURE k3rd_iupac_sp
  END INTERFACE k3rd_iupac
