!........1.........2.........3.........4.........5.........6.........7.........8
MODULE TRAJECTORYIO_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: TRAJECTORYIO                                                       **
!**                                                                           **
!**  PURPOSE: STORES INFORMATION PRODUCED STEP BY STEP IN AN INTERNAL         **
!**    ARRAY, AND WRITES IT IN CHUNKS TO FILE                                 **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    TRAJECTORYIO$NEW(ID)                                                   **
!**    TRAJECTORYIO$SELECT(ID)                                                **
!**    TRAJECTORYIO$SETL4(ID,VAL)                                             **
!**    TRAJECTORYIO$SETI4(ID,VAL)                                             **
!**    TRAJECTORYIO$SETCH(ID,VAL)                                             **
!**    TRAJECTORYIO$ADD(ISTEP,TIME,NSIZE,ARRAY)                               **
!**    TRAJECTORYIO$FLUSHALL                                                  **
!**                                                                           **
!**  TYPICAL USE:                                                             **
!**    1) SET UP A NEW INSTANCE:                                              **
!**      NEW(ID): MAKE NEW INSTANCE                                           **
!**      SELECT(ID)                                                           **
!**      SETL4('ON',.TRUE.)   ! DEFAULT IS ONE                                **
!**      SETI4('NRECX',100)   !OPTIONAL                                       **
!**      SETI4('SKIP',0)                                                      **
!**      SETCH('FILE',FILENAME)                                               **
!**      SELECT('NONE')                                                       **
!**                                                                           **
!**    2) ADD INFORMATION TO BUFFER:                                          **
!**      SELECT(ID)                                                           **
!**      ADD(ISTEP,TIME,NSIZE,ARRAY) ADDS A TIME SLICE TO THE BUFFER          **
!**         WRITING TO FILE. THE BUFFER WRITES TO FILE AUTOMATICALLY WHEN FULL**
!**      SELECT('NONE')                                                       **
!**                                                                           **
!**    3) FLUSH EVERYTHING TO FILE (TO BE DONE BEFORE CLOSING DOWN)           **
!**      SELECT(ID)                                                           **
!**      FLUSHALL WRITE ALL BUFFERS TO FILE                                   **
!**      SELECT('NONE')                                                       **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    DEFAULTS:                                                              **
!**       1) THE TRAJECTORY IS SWITCHED OFF BY DEFAULT                        **
!**       2) SKIP IS SET TO ZERO                                              **
!**                                                                           **
!**    A FILE NAME PRECEEDED BY A * REPRESENTS AN EXTENSION OF THE            **
!**       ROOT NAME OF THE FILEHANDLER OBJECT                                 **
!**                                                                           **
!**    THE CRITERION TO SKIP THE FILE IS GOUVERNED BY ISTEP. IF A TIME SLICE  **
!**       IS ADDED THAT SHOULD BE SKIPPED, THE FUNCTION ADD DOES NOT DO       **
!**       ANYTHING. THE NEXT TIME SLICE TO BE ADDED CAN BE REQUESTED BY       **
!**       SETI4('NEXT',ISTEP)                                                 **
!**                                                                           **
!**    THE OBJECT IGNORES ALL REQUESTS FROM TASKS OTHER THAN THE FIRST NODE   **
!**       OF THE MPIGROUP "MONOMER"                                           **
!**                                                                           **
!**                                                                           **
!**  DEPENDENCIES:                                                            **
!**    ERROR                                                                  **
!**    FILEHANDLER                                                            **
!**    MPI                                                                    **
!**                                                                           **
!*******************************************************************************
TYPE THIS_TYPE
  CHARACTER(32)          :: ID   ! IDENTIFIER; MUST BE A FILE-ID OF FILEHANDLER 
  TYPE(THIS_TYPE),POINTER:: NEXT ! POINTS TO NEXT THIS-TYPE ELEMENT
  CHARACTER(512)         :: FILE ! FILE NAME TO WHICH THE TRAJECTORY IS WRITTEN
  LOGICAL(4)             :: TON  ! ON-OFF SWITCH
  INTEGER(4)             :: NSIZE ! LENGTH OF EACH TRAJECTORY ENTRY
  INTEGER(4)             :: NRECX ! #(TIME SLICES) TO BE STORED BEFORE FLUSH
  REAL(8)                :: SKIP=0
  INTEGER(4)             :: NREC ! INDEX OF LAST WRITTEN ENTRY IN  POINTER ARRAY
                                 ! I=0 INDICATES THAT ARRAYS ARE NOT ALLOCATED 
  INTEGER(4)    ,POINTER :: ISTEP(:)  ! (NSTEP)        TIME STEP NUMBER
  REAL(8)       ,POINTER :: TIME(:)   ! (NSTEP)        ABSOLUTE TIME IN A.U.
  REAL(8)       ,POINTER :: ARRAY(:,:)! (NSIZE,NSTEP)  ARRAY HOLDING 
                                      !                TRAJECTORY INFORMATION 
END TYPE THIS_TYPE
LOGICAL(4)               :: TINI=.FALSE.
INTEGER(4)               :: NSTEP=100 ! #(TIME STEPS) STORED BEFORE WRITING TO FILE
TYPE(THIS_TYPE)  ,POINTER :: FIRST_THIS
TYPE(THIS_TYPE)  ,POINTER :: THIS
END MODULE TRAJECTORYIO_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO_SETFILE(FILE)
!     **************************************************************************
!     **  DEFINE THE RELEVANT FILE                                            **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: FILE
      CHARACTER(32)            :: ID
      LOGICAL(4)               :: EXT
      CHARACTER(512)           :: FILE1
!     **************************************************************************
      ID=THIS%ID
      EXT=.FALSE.
      IF(FILE(1:1).EQ.'*') THEN
        EXT=.TRUE.
        FILE1=FILE(2:)
      ELSE
        FILE1=FILE
      END IF
      CALL FILEHANDLER$SETFILE(ID,EXT,FILE1)
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO_FIRSTTASK(TCHK)
!     **************************************************************************
!     **  CHECKS IF ON THE FIRST NODE OF THE RELEVANT TASK GROUP              **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4)   ,INTENT(OUT):: TCHK
      INTEGER(4)               :: NTASKS,THISTASK
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      TCHK=(THISTASK.EQ.1)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$NEW(ID)
!     **************************************************************************
!     **  CREATE A NEW TRAJECTORY INSTANCE                                    **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)              :: TCHK
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        ALLOCATE(FIRST_THIS)
        THIS=>FIRST_THIS
      ELSE
        THIS=>FIRST_THIS
        DO 
          IF(THIS%ID.EQ.ID) THEN
            CALL ERROR$MSG('CANNOT CREATE TRAJECTORY WITH THE SAME NAME')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('TRAJECTORYIO$NEW')
          END IF
          IF(ASSOCIATED(THIS%NEXT)) THEN
            THIS=>THIS%NEXT
          ELSE
            ALLOCATE(THIS%NEXT)
            THIS=>THIS%NEXT
            EXIT
          END IF
        ENDDO
      END IF 
      THIS%ID   =TRIM(ID)
      THIS%TON  =.FALSE.
      THIS%NRECX=100
      THIS%NSIZE=0
      THIS%NREC =0         ! I=0 IMPLIES THAT THE ARRAYS ARE NOT ALLOCATED 
      THIS%SKIP =0         ! AS DEFAULT, NO TIME STEPS ARE SKIPPED
      NULLIFY(THIS%NEXT)
      NULLIFY(THIS%ISTEP)
      NULLIFY(THIS%TIME)
      NULLIFY(THIS%ARRAY)
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$SELECT(ID)
!     **************************************************************************
!     ** SELECT A CERTAIN TRAJECTORYIO INSTANCE OR UNSELECT WITH 'NONE'       **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)              :: TCHK
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('CREATE TRAJECTORYIO BEFORE SELECTING')
        CALL ERROR$STOP('TRAJECTORYIO$SELECT')
      END IF
!
!     == UNSELECT ==============================================================
      IF(ID.EQ.'NONE') THEN
        NULLIFY(THIS)
        RETURN
      END IF
!
!     == FIND TRAJECTORYIO AND SELECT
      THIS=>FIRST_THIS
      DO WHILE (THIS%ID.NE.ID)
        IF(.NOT.ASSOCIATED(THIS%NEXT)) THEN
          CALL ERROR$MSG('TRAJECTORY DOES NOT EXIST')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('TRAJECTORYIO$SELECT')
        END IF
        THIS=>THIS%NEXT
      ENDDO
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$SETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
      LOGICAL(4)              :: TCHK
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO TRAJECTORY SELECTED; USE TRAJECTORYIO$SELECT FIRST')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$SETL4')
      END IF
!
      IF(ID.EQ.'ON') THEN
        THIS%TON=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$SETL4')
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$GETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
      LOGICAL(4)              :: TCHK
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO TRAJECTORY SELECTED; USE TRAJECTORYIO$SELECT FIRST')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$GETL4')
      END IF
!
      IF(ID.EQ.'ON') THEN
        VAL=THIS%TON
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$GETL4')
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$SETI4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
      LOGICAL(4)              :: TCHK
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO TRAJECTORY SELECTED; USE TRAJECTORYIO$SELECT FIRST')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$SETI4')
      END IF
!
      IF(ID.EQ.'NRECX') THEN
        IF(ASSOCIATED(THIS%ARRAY)) THEN
          CALL ERROR$MSG('ARRAY ALREADY ALLOCATED')
          CALL ERROR$MSG('NO CHANGE OF ARRAY-SIZE POSSIBLE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('TRAJECTORYIO$SETI4')
        END IF
        THIS%NRECX=VAL
      ELSE IF (ID.EQ.'SKIP') THEN
        THIS%SKIP=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$SETCH(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
      LOGICAL(4)              :: TCHK
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO TRAJECTORY SELECTED; USE TRAJECTORYIO$SELECT FIRST')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$SETCH')
      END IF
!
      IF(ID.EQ.'FILE') THEN
!       == A PRECEDING STAR AT THE FILE NAME INDICATES THAT IT IS AN EXTENSION =
!       == TO THE ROOT NAME
        THIS%FILE=VAL
        CALL TRAJECTORYIO_SETFILE(THIS%FILE)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAJECTORYIO$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$ADD(ISTEP,TIME,NSIZE,ARRAY)
!     **************************************************************************
!     **  ADD                                                                 **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: ISTEP
      REAL(8)      ,INTENT(IN) :: TIME
      INTEGER(4)   ,INTENT(IN) :: NSIZE
      REAL(8)      ,INTENT(IN) :: ARRAY(NSIZE)
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: I
      INTEGER(4)               :: NREC
      INTEGER(4)               :: NRECX
      INTEGER(4)               :: NFIL
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO TRAJECTORY SELECTED; USE TRAJECTORYIO$SELECT FIRST')
        CALL ERROR$CHVAL('ID',THIS%ID)
        CALL ERROR$STOP('TRAJECTORYIO$ADD')
      END IF
      IF(.NOT.THIS%TON) RETURN
                            CALL TRACE$PUSH('TRAJECTORYIO$ADD')
!
!     ==========================================================================
!     == ALLOCATE ARRAY                                                       ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(THIS%ARRAY)) THEN
        IF(THIS%NRECX.LE.0) THEN
          CALL ERROR$MSG('NUMBER OF RECORDS NOT SPECIFIED')
          CALL ERROR$CHVAL('ID',THIS%ID)
          CALL ERROR$STOP('TRAJECTORYIO$ADD')
        END IF
        THIS%NSIZE=NSIZE
        NRECX=THIS%NRECX
        ALLOCATE(THIS%ARRAY(NSIZE,NRECX))
        ALLOCATE(THIS%ISTEP(NRECX))
        ALLOCATE(THIS%TIME(NRECX))
      END IF
!
!     ==========================================================================
!     == MAP THE ARRAY ON TO INTERNAL ARRAY                                   ==
!     ==========================================================================
      NREC=THIS%NREC
      TCHK=(NREC.EQ.0)
      IF(.NOT.TCHK) TCHK=(ISTEP.GT.THIS%ISTEP(NREC)+THIS%SKIP)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP
        RETURN  ! DO NOT WRITE THIS TIME SLICE
      END IF
      THIS%NREC=THIS%NREC+1
      NREC=THIS%NREC
      THIS%ISTEP(NREC)=ISTEP
      THIS%TIME(NREC)=TIME
      THIS%ARRAY(:,NREC)=ARRAY(:)
!
!     ==========================================================================
!     == FLUSH DATA IF FULL                                                   ==
!     ==========================================================================
      IF(THIS%NREC.EQ.THIS%NRECX) THEN
        CALL FILEHANDLER$UNIT(TRIM(THIS%ID),NFIL)
        DO I=1,THIS%NRECX
          WRITE(NFIL)THIS%ISTEP(I),THIS%TIME(I),THIS%NSIZE,THIS%ARRAY(:,I)
        ENDDO
        CALL FILEHANDLER$CLOSE(TRIM(THIS%ID))
        THIS%NREC=0
        THIS%ISTEP(:)=0
        THIS%TIME(:)=0.D0
        THIS%ARRAY(:,:)=0.D0
      END IF
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORYIO$FLUSHALL
!     **************************************************************************
!     **  CLEAN OUT STORAGE AREA OF ALL TRAJECTORIES                          **
!     **************************************************************************
      USE TRAJECTORYIO_MODULE
      IMPLICIT NONE
      TYPE(THIS_TYPE),POINTER  :: THIS1
      INTEGER(4)               :: I
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: NFIL
!     **************************************************************************
      CALL TRAJECTORYIO_FIRSTTASK(TCHK)
      IF(.NOT.TCHK) RETURN   ! RETURN IF NOT ON TASK 1
                            CALL TRACE$PUSH('TRAJECTORYIO$FLUSHALL')
      THIS1=>FIRST_THIS
      DO WHILE (ASSOCIATED(THIS1))
        IF(.NOT.THIS1%TON) THEN
          THIS1=>THIS1%NEXT
          CYCLE
        END IF
        IF(THIS1%NREC.EQ.0) THEN  ! NEEDED TO AVOID PROBLEM BELOW, 
          THIS1=>THIS1%NEXT       ! WHEN THIS1%ARRAY IS NOT ALLOCATED
          CYCLE
        END IF
        CALL FILEHANDLER$UNIT(TRIM(THIS1%ID),NFIL)
        DO I=1,THIS1%NREC
          WRITE(NFIL)THIS1%ISTEP(I),THIS1%TIME(I),THIS1%NSIZE,THIS1%ARRAY(:,I)
        ENDDO
        CALL FILEHANDLER$CLOSE(TRIM(THIS1%ID))
        THIS1%NREC=0
        THIS1%ISTEP=0
        THIS1%TIME=0.D0
        THIS1%ARRAY(:,:)=0.D0
        THIS1=>THIS1%NEXT
      ENDDO
                            CALL TRACE$POP
      RETURN
      END
