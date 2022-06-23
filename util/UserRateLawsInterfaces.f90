  PRIVATE :: ARR_dp, ARR_sp
  INTERFACE ARR
     MODULE PROCEDURE ARR_dp
     MODULE PROCEDURE ARR_sp
  END INTERFACE ARR

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

  PRIVATE :: k_3rd_dp, k_3rd_sp
  INTERFACE k_3rd
     MODULE PROCEDURE k_3rd_dp
     MODULE PROCEDURE k_3rd_sp
  END INTERFACE k_3rd

  PRIVATE :: k_3rd_jpl_activation_dp, k_3rd_jpl_activation_sp
  INTERFACE  k_3rd_jpl_activation
     MODULE PROCEDURE  k_3rd_jpl_activation_dp
     MODULE PROCEDURE  k_3rd_jpl_activation_sp
  END INTERFACE  k_3rd_jpl_activation

  PRIVATE :: k_3rd_iupac_dp, k_3rd_iupac_sp
  INTERFACE k_3rd_iupac
     MODULE PROCEDURE k_3rd_iupac_dp
     MODULE PROCEDURE k_3rd_iupac_sp
  END INTERFACE k_3rd_iupac

  PRIVATE :: k_arr_dp, k_arr_sp
  INTERFACE k_arr
     MODULE PROCEDURE k_arr_dp
     MODULE PROCEDURE k_arr_sp
  END INTERFACE k_arr
