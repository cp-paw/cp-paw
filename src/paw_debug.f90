!**   works for each monomer separately
MODULE DEBUG_MODULE
CHARACTER(32)        :: FILENAME="DEBUGOUT"
INTEGER(8)           :: DUNIT=4243         
LOGICAL(4)           :: ALLTASKS=.TRUE.
END MODULE DEBUG_MODULE
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DEBUG$NAN_R8(LEN,OUTPUT,COMMENT)
       USE MPE_MODULE
       USE DEBUG_MODULE
       IMPLICIT NONE
       INTEGER(4),   INTENT(IN)           ::  LEN
       REAL(8),      INTENT(IN)           ::  OUTPUT(LEN)
       CHARACTER(*), INTENT(IN)           ::  COMMENT
       INTEGER(4)                         ::  NTASKS,THISTASK
       INTEGER(4)                         :: I
!      *************************************************************************
       CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
       IF(.NOT.ALLTASKS) THEN
         OPEN(DUNIT,FILE=TRIM(ADJUSTL(FILENAME)),STATUS="UNKNOWN" &
      &         ,ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE" &
      &         ,POSITION="APPEND")
         WRITE(DUNIT,FMT="(3A)",ADVANCE="YES" ) &
      &              "==============",COMMENT,"==============" 
         IF(THISTASK.EQ.1) THEN
           DO I=1,LEN
             IF(OUTPUT(I).NE.OUTPUT(I)) THEN
                WRITE(DUNIT,FMT="(A,I10,F25.14)") &
      &              "THE FOLLOWING ELEMENT IS DUBIOUS: ",I,OUTPUT(I) 
             END IF
           END DO
           CLOSE(DUNIT)
         END IF
       ELSE
         OPEN (DUNIT,FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN" &
     &                 ,ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE" &
     &                 ,POSITION="APPEND")
       DO I=1,LEN
         IF(OUTPUT(I).NE.OUTPUT(I)) THEN
           WRITE(DUNIT,FMT="(2A,I10,A,I4,F25.14)") &
      &          COMMENT," DUBIOUS ELEMENT: ",I," TASK: ",THISTASK,OUTPUT(I) 
         END IF
       END DO
       CLOSE(DUNIT)
       END IF
       END SUBROUTINE DEBUG$NAN_R8
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DEBUG$NAN_C8(LEN,OUTPUT,COMMENT)
       USE MPE_MODULE
       USE DEBUG_MODULE
       IMPLICIT NONE
       INTEGER(4),   INTENT(IN)           ::  LEN
       COMPLEX(8),   INTENT(IN)           ::  OUTPUT(LEN)
       CHARACTER(*), INTENT(IN)           ::  COMMENT
       INTEGER(4)                         ::  NTASKS,THISTASK
       INTEGER(4)                         :: I
!      *************************************************************************
       CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
       IF(.NOT.ALLTASKS) THEN
         OPEN(DUNIT,FILE=TRIM(ADJUSTL(FILENAME)),STATUS="UNKNOWN" &
      &      ,ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE" &
      &      , POSITION="APPEND")
       WRITE (DUNIT, FMT="(3A)", ADVANCE="YES" ) &
      &      "==============",COMMENT,"==============" 

  IF(THISTASK.EQ.1) THEN
     DO I=1,LEN
        IF(OUTPUT(I).NE.OUTPUT(I)) THEN
           WRITE (DUNIT,FMT="(A,I10,2F25.14)")"THE FOLLOWING ELEMENT IS DUBIOUS: ",I,OUTPUT(I) 
        END IF
     END DO
     CLOSE(DUNIT)
  END IF
ELSE
  OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
     DO I=1,LEN
        IF(OUTPUT(I).NE.OUTPUT(I)) THEN
           WRITE (DUNIT,FMT="(2A,I10,A,I4,2F25.14)")COMMENT," DUBIOUS ELEMENT: ",I," TASK: ",THISTASK,OUTPUT(I) 
        END IF
     END DO
     CLOSE(DUNIT)

END IF
END SUBROUTINE DEBUG$NAN_C8


SUBROUTINE DEBUG$NAN_I4(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  INTEGER(4),   INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  INTEGER(4)                         :: I
 !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
IF(.NOT.ALLTASKS) THEN
  OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
  WRITE (DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 

  IF(THISTASK.EQ.1) THEN
     DO I=1,LEN
        IF(OUTPUT(I).NE.OUTPUT(I)) THEN
           WRITE (DUNIT,FMT="(A,I10,I10)")"THE FOLLOWING ELEMENT IS DUBIOUS: ",I,OUTPUT(I) 
        END IF
     END DO
     CLOSE(DUNIT)
  END IF
ELSE
  OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
       &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")

     DO I=1,LEN
        IF(OUTPUT(I).NE.OUTPUT(I)) THEN
           WRITE (DUNIT,FMT="(2A,I10,A,I4,I10)")COMMENT," DUBIOUS ELEMENT: ",I," TASK: ",THISTASK,OUTPUT(I) 
        END IF
     END DO
     CLOSE(DUNIT)

END IF
END SUBROUTINE DEBUG$NAN_I4






SUBROUTINE DEBUG$WRITE_R8(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  REAL(8),      INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
     WRITE (DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 
!     WRITE (DUNIT,*)OUTPUT 
     WRITE (DUNIT,FMT="(F25.14)")OUTPUT 
     CLOSE(DUNIT)
  END IF
  
END SUBROUTINE DEBUG$WRITE_R8





SUBROUTINE DEBUG$WRITE_I4(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  INTEGER(4),   INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 
     WRITE (DUNIT,*)OUTPUT 
     CLOSE(DUNIT)
  END IF
  
END SUBROUTINE DEBUG$WRITE_I4
!
!  
SUBROUTINE DEBUG$WRITE_L4(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  LOGICAL(4),   INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 
     WRITE (DUNIT,*)OUTPUT 
     CLOSE(DUNIT)
  END IF
  
END SUBROUTINE DEBUG$WRITE_L4


SUBROUTINE DEBUG$WRITE_C8(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  COMPLEX(8),   INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 
!     WRITE (DUNIT,*)OUTPUT 
     WRITE (DUNIT,FMT="(2F25.14)")OUTPUT 
     CLOSE(DUNIT)
  END IF
  
END SUBROUTINE DEBUG$WRITE_C8







!THE FOLLOWING CODE IS HERE FOR HISTORICAL REASONS AND WILL BE DELETED!


SUBROUTINE PAW$DEBUG_R8(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  REAL(8),      INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", FORM="FORMATTED",ACTION="WRITE", POSITION="APPEND")
     WRITE (DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 
!     WRITE (DUNIT,*)OUTPUT 
     WRITE (DUNIT,FMT="(F25.14)")OUTPUT 
     CLOSE(DUNIT)
  END IF
  
END SUBROUTINE PAW$DEBUG_R8





SUBROUTINE PAW$DEBUG_I4(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  INTEGER(4),   INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 
     WRITE (DUNIT,*)OUTPUT 
     CLOSE(DUNIT)
  END IF
  
END SUBROUTINE PAW$DEBUG_I4

SUBROUTINE PAW$DEBUG_C8(LEN,OUTPUT,COMMENT)
  USE MPE_MODULE
  USE DEBUG_MODULE
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)           ::  LEN
  COMPLEX(8),   INTENT(IN)           ::  OUTPUT(LEN)
  CHARACTER(*), INTENT(IN)           ::  COMMENT
  INTEGER(4)                         ::  NTASKS,THISTASK
  !...................................................

  CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=DUNIT, FMT="(3A)", ADVANCE="YES" )"==============",COMMENT,"==============" 
!     WRITE (DUNIT,*)OUTPUT 
     WRITE (DUNIT,FMT="(2F25.14)")OUTPUT 
     CLOSE(DUNIT)
  END IF
END SUBROUTINE PAW$DEBUG_C8
!
!     ..........................................................................
      SUBROUTINE DEBUG$ROUND_R8(LEN,VAL,TOL)
      USE MPE_MODULE
      USE DEBUG_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    ::  LEN
      REAL(8)   ,INTENT(INOUT) ::  VAL(LEN)
      REAL(8)   ,INTENT(IN)    ::  TOL
!     **************************************************************************
      VAL(:)=NINT(VAL(:)/TOL)*TOL
      END SUBROUTINE DEBUG$ROUND_R8
!
!     ..........................................................................
      SUBROUTINE DEBUG$WAVEFUNCTIONS()
      USE WAVES_MODULE
      USE DEBUG_MODULE
      IMPLICIT NONE
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSI1(:,:)
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: IKPT,ISPIN,IB,IDIM,IWAVE
      INTEGER(4)              :: NG1,NGG,NGL,NBH,NB
      REAL(8)                 :: XK(3)
      REAL(8)     ,ALLOCATABLE:: GVECL(:,:)
      REAL(8)     ,ALLOCATABLE:: GVECG(:,:)
      REAL(8)                 :: GBAS(3,3)
      INTEGER(4)              :: NREC1,ISVAR
      CHARACTER(8)            :: KEY
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: A,B
      INTEGER(4)  ,ALLOCATABLE:: IGVECG(:,:)
      INTEGER(4)  ,ALLOCATABLE:: IGVECL(:,:)
      INTEGER(4)              :: IG
      CHARACTER(64)           :: COMMENT='WAVE FUNCTION PRINTOUT'
      REAL(8)      ,PARAMETER :: TOL=1.D-12
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
PRINT*,'ANFANG DEBUG$WAVEFUNCTIONS',THISTASK
  IF(THISTASK.EQ.1) THEN
     OPEN (DUNIT, FILE=TRIM(ADJUSTL(FILENAME)), STATUS="UNKNOWN", &
          &ACCESS="SEQUENTIAL", ACTION="WRITE", POSITION="APPEND")
     WRITE (UNIT=DUNIT, FMT="(3A)", ADVANCE="YES" ) &
 &                 "==============",COMMENT,"==============" 
  END IF
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        NBH=THIS%NBH
        NB=THIS%NB
        CALL PLANEWAVE$GETI4('NGG',NGG)
        ALLOCATE(IGVECL(3,NGL))
        CALL PLANEWAVE$GETI4A('IGVEC',3*NGL,IGVECL)
        ALLOCATE(IGVECG(3,NGG))
        CALL PLANEWAVE$COLLECTI4(3,NGL,IGVECL,NGG,IGVECG)
        DEALLOCATE(IGVECL)
        ALLOCATE(PSIG(NGG,NDIM))
        ALLOCATE(PSI1(NGL,NDIM))
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          DO IB=1,NBH
            DO IDIM=1,NDIM
              PSI1(:,IDIM)=THIS%HPSI(:,IDIM,IB)  !<<<<==========================
              DO IG=1,NGL
                A=REAL(PSI1(IG,IDIM))
                B=AIMAG(PSI1(IG,IDIM))
                A=TOL*REAL(NINT(A/TOL),KIND=8)
                B=TOL*REAL(NINT(B/TOL),KIND=8)
                PSI1(IG,IDIM)=CMPLX(A,B,KIND=8)
              ENDDO
            ENDDO
            DO IDIM=1,NDIM
              CALL PLANEWAVE$COLLECTC8(1,NGL,PSI1(1,IDIM),NGG,PSIG(1,IDIM))
            ENDDO
            IF(THISTASK.EQ.1) THEN
              DO IG=1,NGG
                WRITE(DUNIT,*)IB,IGVECG(:,IG),PSIG(IG,:)
              ENDDO
            END IF
          ENDDO
          DEALLOCATE(PSI1)
          DEALLOCATE(PSIG)
        ENDDO 
        DEALLOCATE(IGVECG)
      ENDDO
      IF(THISTASK.EQ.1) THEN
        CLOSE(DUNIT)
      END IF
PRINT*,'ENDE DEBUG$WAVEFUNCTIONS',THISTASK
      CALL MPE$SYNC('MONOMER')
      CALL ERROR$NORMALSTOP
      END SUBROUTINE DEBUG$WAVEFUNCTIONS
!
!     ..........................................................................
      SUBROUTINE DEBUG$CHKPARACONS_R8A(ID,LEN,VAL)
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LEN
      REAL(8)   ,INTENT(IN)    :: VAL(LEN)
      CHARACTER(*),INTENT(IN) ::  ID
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)   ,ALLOCATABLE   :: VALARR(:,:)
      INTEGER(4)               :: ITASK
      INTEGER(4)               :: ISVARARR(1),ISVAR
!     ==========================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ALLOCATE(VALARR(LEN,NTASKS))
      VALARR(:,1)=VAL(:)
      DO ITASK=2,NTASKS
        IF(ITASK.EQ.THISTASK) THEN
          CALL MPE$SEND('MONOMER',1,ITASK,VAL)
        END IF
        IF(THISTASK.EQ.1) THEN
          CALL MPE$RECEIVE('MONOMER',ITASK,ITASK,VALARR(:,ITASK))
        END IF
      ENDDO
!     =================================================================
!     == ANALYSE                                                     ==
!     =================================================================
      IF(THISTASK.EQ.1) THEN
        DO ITASK=2,NTASKS
          ISVARARR=MAXLOC(ABS(VALARR(:,ITASK)-VALARR(:,1)))
          ISVAR=ISVARARR(1)
          PRINT*,TRIM(ID),ITASK,ISVAR,VALARR(ISVAR,ITASK)-VALARR(ISVAR,1),VALARR(ISVAR,ITASK),VALARR(ISVAR,1)
        ENDDO
      END IF
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE DEBUG$MPEIDENT_R8A(LEN,VAL,IERR)
!     **                                                           **
!     ** THIS ROUTINE CHECKS FOR IF VAL IS IDENTICAL ON ALL NODES  **
!     ** WITH IERR IT PROVIDES THE NUMBER OF TASK, FOR WHICH VAL   **
!     ** DIFFERS FROM THE FIRST TASK. IF IERR=0 ALL ARE IDENTICAL  **
!     **                                                           **
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LEN
      REAL(8)   ,INTENT(IN)    :: VAL(LEN)
      INTEGER(4),INTENT(OUT)   :: IERR
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)   ,ALLOCATABLE   :: VALARR(:,:)
      INTEGER(4)               :: ITASK
      INTEGER(4)               :: ISVARARR(1),ISVAR
!     ******************************************************************
      IERR=0
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ALLOCATE(VALARR(LEN,NTASKS))
      VALARR(:,1)=VAL(:)
      DO ITASK=2,NTASKS
        IF(ITASK.EQ.THISTASK) THEN
          CALL MPE$SEND('MONOMER',1,ITASK,VAL)
        END IF
        IF(THISTASK.EQ.1) THEN
          CALL MPE$RECEIVE('MONOMER',ITASK,ITASK,VALARR(:,ITASK))
        END IF
      ENDDO
!     =================================================================
!     == ANALYSE                                                     ==
!     =================================================================
      IF(THISTASK.EQ.1) THEN
        IERR=0
        DO ITASK=2,NTASKS
          ISVARARR=MAXLOC(ABS(VALARR(:,ITASK)-VALARR(:,1)))
          ISVAR=ISVARARR(1)
          IF(VALARR(ISVAR,ITASK)-VALARR(ISVAR,1).NE.0.D0) IERR=IERR+1
        ENDDO
      END IF
      RETURN
      END

