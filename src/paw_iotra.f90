!.......................................................................
MODULE TRAJECTORYIO_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: TRAJECTORYIO                                               **
!**                                                                   **
!**  PURPOSE: STORES INFORMATION PRODUCED STEP BY STEP IN AN INTERNAL **
!**    ARRAY, AND WRITES IT IN CHUNKS TO FILE                         **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    TRAJECTORYIO$INITIALIZE(NSTEP)                                 **
!**    TRAJECTORYIO$ON(ID,TON)                                        **
!**    TRAJECTORYIO$ADD(ID,NFI,TIME,NSIZE,ARRAY)                      **
!**    TRAJECTORYIO$FLUSH                                             **
!**                                                                   **
!**  TYPICAL USE:                                                     **
!**    1) ADD INFORMATION FOR EACH TIME STEP AND EACH TRAJECTORY      **
!**       USING TRAJECTORY$ADD (NO INITIALIZATION REQUIRED)           **
!**    2) FLUSH ALL TRAJECTORIES BEFORE TERMINATING THE PROCESS       **
!**       (INTERMEDIATE FLUSHES ARE ATOMATICALLY EXECUTED)            **
!**    3) THE ID MUST BE IDENTICAL TO A FILE-ID OF A FILE             **
!**       IN THE FILEHANDLER OBJECT                                   **
!**                                                                   **
!**  REMARKS:                                                         **
!**    1) TRAJECTORYIO$INITIALIZE IS OPTIONAL AND CAN BE CALLED       **
!**       ANY TIME. (IT ONLY RESETS THE SIZE OF WORK ARRAYS)          **
!**                                                                   **
!**    2) TRAJECTORYIO$ON ALLOWS TO SWITCH OF SOME TRAJECTORIES       **
!**                                                                   **
!**  DEPENDENCIES:                                                    **
!**    ERROR                                                          **
!**    FILEHANDLER                                                    **
!**    MPI                                                            **
!**                                                                   **
!***********************************************************************
LOGICAL(4)               :: TINI=.FALSE.
INTEGER(4)               :: NSTEP=100 ! %(TIME STEPS) STORED BEFORE WRITING TO FILE
TYPE XXX_TYPE
  CHARACTER(32)          :: ID        ! IDENTIFIER; MUST BE A FILE-ID OF FILEHANDLER 
  TYPE(XXX_TYPE),POINTER :: NEXT      ! POINTS TO NEXT XXX-TYPE ELEMENT
  LOGICAL(4)             :: TON       ! ON-OFF SWITCH
  INTEGER(4)             :: NSIZE     ! LENGTH OF EACH TRAJECTORY ENTRY
  INTEGER(4)             :: I         ! INDEX OF NEXT FREE ENTRY IN THE POINTER ARRAY
                                      ! I=0 INDICATES THAT ARRAYS ARE NOT ALLOCATED 
  INTEGER(4)    ,POINTER :: ISTEP(:)  ! (NSTEP)         TIME STEP NUMBER
  REAL(8)       ,POINTER :: TIME(:)   ! (NSTEP)         ABSOLUTE TIME IN A.U.
  REAL(8)       ,POINTER :: ARRAY(:,:)! (NSIZE,NSTEP)   ARRAY HOLDING TRAJECTORY INFORMATION 
END TYPE XXX_TYPE
TYPE(XXX_TYPE),POINTER :: XXXSTART
CONTAINS
!     ..................................................................
      SUBROUTINE TRAJECTORYIO_NEW(ID,XXX)
!     ******************************************************************
!     **  CREATES A NEW ID FOR A TRAJECTORY TYPE                      **
!     **  SIZE AND ARRAY ARE NOT YET SET                              **
!     ******************************************************************
      CHARACTER(*)  ,INTENT(IN) :: ID
      TYPE(XXX_TYPE),POINTER    :: XXX
!     ******************************************************************
      ALLOCATE(XXX)
      XXX%ID   =TRIM(ID)
      XXX%TON  =.FALSE.
      XXX%NSIZE=0
      XXX%I    =0         ! I=0 IMPLIES THAT THE ARRAYS ARE NOT ALLOCATED 
      NULLIFY(XXX%NEXT)
      NULLIFY(XXX%ISTEP)
      NULLIFY(XXX%TIME)
      NULLIFY(XXX%ARRAY)
      RETURN
      END SUBROUTINE TRAJECTORYIO_NEW
!
!     ..................................................................
      SUBROUTINE TRAJECTORYIO_FIND(ID,XXX)
!     *******************************************************************
!     *******************************************************************
      CHARACTER(*),INTENT(IN)  :: ID
      TYPE(XXX_TYPE),POINTER   :: XXX           
!     *******************************************************************
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        CALL TRAJECTORYIO_NEW(ID,XXXSTART)
        XXX=>XXXSTART
        RETURN
      END IF
!
      XXX=>XXXSTART
      DO WHILE (TRIM(ID).NE.TRIM(XXX%ID)) 
        IF(ASSOCIATED(XXX%NEXT)) THEN
          XXX=>XXX%NEXT
        ELSE
          CALL TRAJECTORYIO_NEW(ID,XXX%NEXT)
          XXX=>XXX%NEXT
          EXIT
        END IF
      ENDDO
      RETURN
      END SUBROUTINE TRAJECTORYIO_FIND
!
!     ..................................................................
      SUBROUTINE TRAJECTORYIO_FLUSH(XXX)
!     ******************************************************************
!     ** FLUSH                                                        **   
!     **                                                              **   
!     ** WRITES A TRAJECTROY ARRAYS TO FILE                           **   
!     ******************************************************************
      IMPLICIT NONE
      TYPE(XXX_TYPE),INTENT(INOUT):: XXX
      INTEGER(4)                  :: I
      INTEGER(4)                  :: NFIL
!     ******************************************************************
      IF(XXX%I.EQ.0) RETURN
!
!     ==================================================================
!     == WRITE DATA TO FILE                                           ==
!     ==================================================================
      CALL FILEHANDLER$UNIT(TRIM(XXX%ID),NFIL)
      DO I=1,XXX%I-1
        WRITE(NFIL)XXX%ISTEP(I),XXX%TIME(I),XXX%NSIZE,XXX%ARRAY(:,I)
      ENDDO
      CALL FILEHANDLER$CLOSE(TRIM(XXX%ID))
!
!     ==================================================================
!     == RESET ARRAYS                                                 ==
!     ==================================================================
      XXX%I=1
      XXX%ISTEP(:)=0
      XXX%TIME(:)=0
      XXX%ARRAY(:,:)=0.D0
      RETURN
      END SUBROUTINE TRAJECTORYIO_FLUSH
END MODULE TRAJECTORYIO_MODULE
!
!     ..................................................................
      SUBROUTINE TRAJECTORYIO$INITIALIZE(NSTEP_)
!     ******************************************************************
!     **  INITIALIZE                                                  **   
!     **  SETS THE MAXIMUM NUMBER OF STEPS IN THE WORK ARRAYS         **   
!     **  CAN BE CALLED ANY TIME: ARRAYS ARE RESIZED AUTOMATICALLY    **   
!     ******************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NSTEP_   ! #(TIME STEPS) ON WORK ARRAYS
      INTEGER(4)             :: SIZE
      TYPE(XXX_TYPE),POINTER :: XXX           
      REAL(8)   ,ALLOCATABLE :: WORK(:,:)
      INTEGER(4),ALLOCATABLE :: IWORK(:)
      INTEGER(4)             :: ISTEP
!     ******************************************************************
      IF(NSTEP_.EQ.NSTEP) RETURN
!
!     ==================================================================
!     ==  FIRST REALLOCATE ARRAYS IF NECCESARY                        ==
!     ==================================================================
      IF(TINI) THEN
        XXX=>XXXSTART
        DO WHILE (ASSOCIATED(XXX))
          IF(XXX%I.NE.0) THEN   ! XXX%I=0 : ARRAYS NOT ALLOCATED
!           == FLUSH TRAJECTORIES ======================================
            CALL TRAJECTORYIO_FLUSH(XXX)
!           == REALLOCATE ARRAYS =======================================
            SIZE=XXX%NSIZE
            DEALLOCATE(XXX%ARRAY)
            ALLOCATE(XXX%ARRAY(SIZE,NSTEP_))
            DEALLOCATE(XXX%TIME)
            ALLOCATE(XXX%TIME(NSTEP_))
            DEALLOCATE(XXX%ISTEP)
            ALLOCATE(XXX%ISTEP(NSTEP_))
          END IF
          XXX=>XXX%NEXT
        ENDDO
      END IF
!
!     ==================================================================
!     ==  NOW CHANGE THE SIZE VALUE OF THE BUFFER ARRAYS              ==
!     ==================================================================
      NSTEP=NSTEP_
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TRAJECTORYIO$ON(ID_,TON_)
!     ******************************************************************
!     **  SWITCHES A TRAJECTORY ON AND OFF                            **   
!     ******************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*)  ,INTENT(IN) :: ID_
      LOGICAL(4)    ,INTENT(IN) :: TON_
      TYPE(XXX_TYPE),POINTER    :: XXX
!     ******************************************************************
      CALL TRAJECTORYIO_FIND(ID_,XXX)
      XXX%TON=TON_
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TRAJECTORYIO$ADD(ID_,NFI_,TIME_,NSIZE_,ARRAY_)
!     ******************************************************************
!     **  ADD                                                         **   
!     ******************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: ID_
      INTEGER(4)   ,INTENT(IN) :: NFI_
      REAL(8)      ,INTENT(IN) :: TIME_
      INTEGER(4)   ,INTENT(IN) :: NSIZE_
      REAL(8)      ,INTENT(IN) :: ARRAY_(NSIZE_)
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)               :: ISTEP
      TYPE(XXX_TYPE),POINTER   :: XXX
!     ****************************************************************** 
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.GT.1)RETURN
!
!     ==================================================================
!     == FIND TEH CORRECT TRAJECTORY                                  ==
!     ==================================================================
      CALL TRAJECTORYIO_FIND(ID_,XXX)
      IF(.NOT.XXX%TON) RETURN
!
!     ==================================================================
!     == INITIALIZE ARRAY IF NEW                                      ==
!     ==================================================================
      IF(XXX%I.EQ.0) THEN
        XXX%NSIZE=NSIZE_
        XXX%I=1   ! XXX%I POINTS TO THE FIRST FREE ENTRY
        ALLOCATE(XXX%ARRAY(NSIZE_,NSTEP))
        ALLOCATE(XXX%TIME(NSTEP))
        ALLOCATE(XXX%ISTEP(NSTEP))
      END IF
!
!     ==================================================================
!     == PROTECT AGAINST OVERFLOW                                     ==
!     ==================================================================
      IF(NSIZE_.NE.XXX%NSIZE)THEN
        CALL ERROR$MSG('SIZE MISMATCH ERROR')
        CALL ERROR$CHVAL('ID',TRIM(XXX%ID))
        CALL ERROR$I4VAL('NSIZE_',NSIZE_)
        CALL ERROR$I4VAL('XXX%NSIZE',XXX%NSIZE)
        CALL ERROR$STOP('TRAJECTORYIO$ADD')
      ENDIF
      IF(XXX%I.GT.NSTEP) THEN
        CALL ERROR$MSG('STEP NUMBER OUT OF BOUNDS')
        CALL ERROR$CHVAL('ID',TRIM(XXX%ID))
        CALL ERROR$I4VAL('XXX%I',XXX%I)
        CALL ERROR$I4VAL('NSTEP',NSTEP)
        CALL ERROR$STOP('TRAJECTORYIO$ADD')
      END IF
!
!     ==================================================================
!     == MAP THE ARRAY ON TO INTERNAL ARRAY                           ==
!     ==================================================================
      ISTEP=XXX%I
      XXX%ISTEP(ISTEP)=NFI_
      XXX%TIME(ISTEP)=TIME_
      XXX%ARRAY(:,ISTEP)=ARRAY_(:)
      XXX%I=XXX%I+1
!
!     ==================================================================
!     == FLUSH TRAJECTORY IF WORKSPACE IS FULL                        ==
!     ==================================================================
      IF(XXX%I.GT.NSTEP) THEN
        CALL TRAJECTORYIO_FLUSH(XXX)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TRAJECTORYIO$FLUSH
!     ******************************************************************
!     ** FLUSH                                                        **   
!     **                                                              **   
!     ** WRITES ALL TRAJECTROY ARRAYS TO FILE                         **   
!     ******************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      TYPE(XXX_TYPE),POINTER   :: XXX
      INTEGER(4)               :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.GT.1)RETURN
      XXX=>XXXSTART
      DO WHILE (ASSOCIATED(XXX))
        CALL TRAJECTORYIO_FLUSH(XXX)
        XXX=>XXX%NEXT
      ENDDO
      RETURN
      END
