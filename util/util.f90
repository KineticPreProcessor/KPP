! ****************************************************************
!
! InitSaveData - Opens the data file for writing
!   Parameters :
!
! ****************************************************************

      SUBROUTINE InitSaveData ()

      USE KPP_ROOT_Parameters

      open(10, file='KPP_ROOT.dat')

      END SUBROUTINE InitSaveData

! End of InitSaveData function
! ****************************************************************

! ****************************************************************
!                            
! SaveData - Write LOOKAT species in the data file 
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE SaveData ()

      USE KPP_ROOT_Global
      USE KPP_ROOT_Monitor

      INTEGER i

      WRITE(10,999) (TIME-TSTART)/3600.D0,  &
      (C(LOOKAT(i))/CFACTOR, i=1,NLOOKAT)
 999  FORMAT(E24.16,100(1X,E24.16))

      END SUBROUTINE SaveData

! End of SaveData function
! ****************************************************************

! ****************************************************************
!
! CloseSaveData - Close the data file
!   Parameters :
!
! ****************************************************************

      SUBROUTINE CloseSaveData ()

      USE KPP_ROOT_Parameters

      CLOSE(10)

      END SUBROUTINE CloseSaveData

! End of CloseSaveData function
! ****************************************************************

! ****************************************************************
!
! GenerateMatlab - Generates MATLAB file to load the data file
!   Parameters :
!                It will have a character string to prefix each
!                species name with.
!
! ****************************************************************

      SUBROUTINE GenerateMatlab ( PREFIX )

      USE KPP_ROOT_Parameters
      USE KPP_ROOT_Global
      USE KPP_ROOT_Monitor


      CHARACTER(LEN=8) PREFIX
      INTEGER i

      open(20, file='KPP_ROOT.m')
      write(20,*) 'load KPP_ROOT.dat;'
      write(20,990) PREFIX
990   FORMAT(A1,'c = KPP_ROOT;')
      write(20,*) 'clear KPP_ROOT;'
      write(20,991) PREFIX, PREFIX
991   FORMAT(A1,'t=',A1,'c(:,1);')
      write(20,992) PREFIX
992   FORMAT(A1,'c(:,1)=[];')

      do i=1,NLOOKAT
        write(20,993) PREFIX, SPC_NAMES(LOOKAT(i)), PREFIX, i
993     FORMAT(A1,A6,' = ',A1,'c(:,',I2,');')
      end do

      CLOSE(20)

      END SUBROUTINE GenerateMatlab

! End of GenerateMatlab function
! ****************************************************************


! ****************************************************************
!
! Integrator_Update_Options - determine whether to call Update_RCONST,
!    Update_PHOTO, and Update_SUN from within the integrator
!
!   Parameters:
!    option (input)
!        = -1 :  Do not call Update_* functions within the integrator
!        =  0 :  Status quo: Call whichever functions are normally called
!        =  1 :  Call Update_RCONST from within the integrator
!        =  2 :  Call Update_PHOTO from within the integrator
!        =  3 :  Call Update_RCONST and Update_PHOTO from within the int.
!        =  4 :  Call Update_SUN from within the integrator
!        =  5 :  Call Update_SUN and Update_RCONST from within the int.
!        =  6 :  not implemented
!        =  7 :  not implemented
!
!    Do_Update_RCONST (output):
!        =T : Calls Update_RCONST from within the integrator
!        =F : Does not call UPDATE_RCONST from w/in the int.
!
!    Do_Update_PHOTO (output):
!        =T : Calls Update_PHOTO from within the integrator
!        =F : Does not call UPDATE_PHOTO from w/in the int.
!
!    Do_Update_SUN (output):
!        =T : Calls Update_SUN from within the integrator
!        =F : Does not call UPDATE_SUN from w/in the int.
!
! ****************************************************************

      SUBROUTINE Integrator_Update_Options( option,            &
                                            Do_Update_RConst,  &
                                            Do_Update_Photo,   &
                                            Do_Update_Sun     )

      !~~~> Input variables
      INTEGER, INTENT(IN)  :: option

      !~~~> Output variables
      LOGICAL, INTENT(OUT) :: Do_Update_RCONST
      LOGICAL, INTENT(OUT) :: Do_Update_PHOTO
      LOGICAL, INTENT(OUT) :: Do_Update_SUN

      ! Option -1: turn off all Update_* calls within the integrator
      IF ( option == -1 ) THEN
         Do_Update_RCONST = .FALSE.
         Do_Update_PHOTO  = .FALSE.
         Do_Update_SUN    = .FALSE.
         RETURN
      ENDIF

      ! Otherwise determine from the value passed
      Do_Update_RCONST = ( IAND( option, 1 ) > 0 )
      Do_Update_PHOTO  = ( IAND( option, 2 ) > 0 )
      Do_Update_SUN    = ( IAND( option, 4 ) > 0 )

      END SUBROUTINE Integrator_Update_Options

! End of Integrator_Update_Options function
! ****************************************************************
