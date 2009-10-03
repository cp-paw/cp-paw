!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)               :: LL_STP
      CHARACTER(256)              :: SETUPREPORT      
      CHARACTER(16)               :: Selection
      INTEGER(4)                  :: NARGS
      INTEGER(4)                  :: NFIL
!     **************************************************************************
      CALL LIB$NARGS(NARGS)
      IF(NARGS.NE.2) THEN
        CALL ERROR$MSG('INCORRECT NUMBER OF COMMAND LINE ARGUMENTS')
        CALL ERROR$MSG('CALLING SEQUENCE:')
        CALL ERROR$MSG('PAW_STPA.X ARG SETUPREPORTFILENAME')
        CALL ERROR$STOP('MAIN')
      END IF
      CALL LIB$GETARG(1,selection)
      selection=+selection
      CALL LIB$GETARG(2,SETUPREPORT)
!
!     ==========================================================================
!     == READ INPUT FILE                                                      ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('STP_REPORT',.FALSE.,SETUPREPORT)
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('DAT',.FALSE.,-'TMP.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('DAT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('DAT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('DAT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('DAT','FORM','FORMATTED')
!
!     ==========================================================================
!     == READ INPUT FILE                                                      ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_STP)
      CALL FILEHANDLER$UNIT('STP_REPORT',NFIL)
      CALL LINKEDLIST$READ(LL_STP,NFIL,'~')
!
!     ==========================================================================
!     == TAKE CARE OF SCATTERING PROPERTIES                                   ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('DAT',NFIL)
      IF(SELECTION.EQ.'SCATTERING') THEN
        CALL SCATTERING(LL_STP,NFIL)
!
!     ==========================================================================
!     == CONSTRUCT FILE FOR ATOMIC WAVE FUNCTIONS                             ==
!     ==========================================================================
      ELSE IF(SELECTION.EQ.'AEPSI') THEN
        CALL WAVEFUNCTIONS(LL_STP,NFIL,'AEPSI')
      ELSE IF(SELECTION.EQ.'UPSI') THEN
        CALL WAVEFUNCTIONS(LL_STP,NFIL,'UPSI')!
      ELSE IF(SELECTION.EQ.'UPSISM') THEN
        CALL WAVEFUNCTIONS(LL_STP,NFIL,'UPSISM')
!
!     ==========================================================================
!     == CONSTRUCT FILE FOR AUGMENTATION                                      ==
!     ==========================================================================
      ELSE IF(SELECTION.EQ.'AEPHI') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'AEPHI')
      ELSE IF(SELECTION.EQ.'PSPHI') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'PSPHI')
      ELSE IF(SELECTION.EQ.'NLPHI') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'NLPHI')
      ELSE IF(SELECTION.EQ.'QPHI') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'QPHI')
      ELSE IF(SELECTION.EQ.'PRO') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'PRO')
      ELSE IF(SELECTION.EQ.'AEPHIDOT') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'AEPHIDOT')
      ELSE IF(SELECTION.EQ.'PSPHIDOT') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'PSPHIDOT')
      ELSE IF(SELECTION.EQ.'QPHIDOT') THEN
        CALL AUGMENTATION(LL_STP,NFIL,'QPHIDOT')
      ELSE 
        CALL ERROR$MSG('SELECTION NOT RECOGNIZED')
        CALL ERROR$CHVAL('SELECTION',SELECTION)
        CALL ERROR$STOP('MAIN')
      END IF
      CALL FILEHANDLER$CLOSE('DAT')
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE augmentation(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      INTEGER(4)                  :: NR
      REAL(8)      ,ALLOCATABLE   :: R(:)
      REAL(8)      ,ALLOCATABLE   :: PHI(:,:)
      INTEGER(4)                  :: lnx
      INTEGER(4)                  :: ln,IR
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',1)
      CALL LINKEDLIST$GET(LL_STP,'NR',1,NR)
      ALLOCATE(R(NR))
      CALL LINKEDLIST$GET(LL_STP,'R',1,R)
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION',1)
      CALL LINKEDLIST$GET(LL_STP,'LNX',1,LNX)
      ALLOCATE(PHI(NR,LNX))
      DO LN=1,LNX
        IF(ID.EQ.'AEPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'PSPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'PSPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'NLPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'NLPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'QPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'QPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'PRO') THEN
          CALL LINKEDLIST$GET(LL_STP,'PRO',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'AEPHIDOT') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPHIDOT',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'PSPHIDOT') THEN
          CALL LINKEDLIST$GET(LL_STP,'PSPHIDOT',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'NLPHIDOT') THEN
          CALL LINKEDLIST$GET(LL_STP,'NLPHIDOT',LN,PHI(:,LN))
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('augmentation')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PhI(IR,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVEFUNCTIONS(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      INTEGER(4)                  :: NR
      REAL(8)      ,ALLOCATABLE   :: R(:)
      REAL(8)      ,ALLOCATABLE   :: PSI(:,:)
      INTEGER(4)                  :: NB
      INTEGER(4)                  :: IB,IR
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',1)
      CALL LINKEDLIST$GET(LL_STP,'NR',1,NR)
      ALLOCATE(R(NR))
      CALL LINKEDLIST$GET(LL_STP,'R',1,R)
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'ATOM',1)
      CALL LINKEDLIST$GET(LL_STP,'NB',1,NB)
      ALLOCATE(PSI(NR,NB))
      DO IB=1,NB
        IF(ID.EQ.'AEPSI') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPSI',IB,PSI(:,IB))
        ELSE IF(ID.EQ.'UPSI') THEN
          CALL LINKEDLIST$GET(LL_STP,'UPSI',IB,PSI(:,IB))
!          psi(:,ib)=psi(:,ib)/maxval(abs(psi(:,ib)))
        ELSE IF(ID.EQ.'UPSISM') THEN
          CALL LINKEDLIST$GET(LL_STP,'UPSI_SMALL',IB,PSI(:,IB))
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVEFUNCTIONS')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PSI(IR,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCATTERING(LL_STP,NFIL)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      REAL(8)                     :: EMIN,EMAX
      INTEGER(4)                  :: NE
      INTEGER(4)                  :: LX
      REAL(8)     ,ALLOCATABLE    :: AEPHASE(:,:)
      REAL(8)     ,ALLOCATABLE    :: PSPHASE(:,:)
      INTEGER(4)                  :: IE,L
      REAL(8)                     :: E
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SCATTERING',1)
      CALL LINKEDLIST$GET(LL_STP,'EMIN',1,EMIN)
      CALL LINKEDLIST$GET(LL_STP,'EMAX',1,EMAX)
      CALL LINKEDLIST$GET(LL_STP,'NE',1,NE)
      CALL LINKEDLIST$GET(LL_STP,'LX',1,LX)
      ALLOCATE(AEPHASE(NE,LX+1))
      ALLOCATE(PSPHASE(NE,LX+1))
      DO L=0,LX
        CALL LINKEDLIST$GET(LL_STP,'AEPHASE',L+1,AEPHASE(:,L+1))
        CALL LINKEDLIST$GET(LL_STP,'PAWPHASE',L+1,PSPHASE(:,L+1))
      ENDDO
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)
        WRITE(NFIL,*)E,AEPHASE(IE,:),PSPHASE(IE,:)
      ENDDO
      DEALLOCATE(AEPHASE)
      DEALLOCATE(PSPHASE)
      RETURN
      END
