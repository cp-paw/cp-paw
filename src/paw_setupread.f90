!
!............................................................................
MODULE SETUPREAD_MODULE
USE LINKEDLIST_MODULE
INTEGER(4)              :: NFIL=0
TYPE(LL_TYPE)           :: LL_STP
INTEGER(4)              :: GID=0
INTEGER(4)              :: LNX=0
INTEGER(4)              :: NR=0
INTEGER(4)              :: NB    ! #(STATES IN AE ATOMIC CALCULATION)
INTEGER(4)              :: NC    ! #(CORE STATES)
END MODULE SETUPREAD_MODULE
!
!     ......................................................................
      SUBROUTINE SETUPREAD$NEW(NFIL_,TCHK)
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE LINKEDLIST$READ(LL_,NFIL,CID_)
        USE LINKEDLIST_MODULE, ONLY: LL_TYPE 
        TYPE(LL_TYPE),INTENT(IN) :: LL_
        INTEGER(4)   ,INTENT(IN) :: NFIL
        CHARACTER(*) ,INTENT(IN),OPTIONAL :: CID_ ! RELEVANT PROCESSOR GROUP (SEE MPE OBECT)
        END SUBROUTINE LINKEDLIST$READ
      END INTERFACE
      INTEGER(4)  ,INTENT(IN) :: NFIL_
      LOGICAL(4)  ,INTENT(OUT):: TCHK
      CHARACTER(6)            :: STRING
      CHARACTER(16)           :: TYPE
      REAL(8)                 :: R1
      REAL(8)                 :: DEX
!     **********************************************************************
!     == FIRST REMOVE OLD LIST IF NEEDED
      IF(NFIL.NE.0) THEN
        CALL LINKEDLIST$SELECT(LL_STP,'~')
!       CALL LINKEDLIST$DELETE(LL_STP)
        CALL LINKEDLIST$RMLIST(LL_STP,'SETUP')
        GID=0
        LNX=0
        NR=0
        NFIL=0
      END IF
      IF(NFIL_.EQ.0) RETURN
!
!     ======================================================================  
!     == CHECK IF THE FILE IS A LINKEDLIST                                ==  
!     ======================================================================  
!PRINT*,'CHECK'
      TCHK=.TRUE.
      REWIND NFIL_
      READ(NFIL_,*)STRING
      REWIND NFIL_
      IF(STRING.NE.'!SETUP') THEN
        TCHK=.FALSE.  
        RETURN
      END IF
      NFIL=NFIL_
!
!     ======================================================================  
!     == CREATE NEW LINKED LIST AND READ                                  ==  
!     ======================================================================  
!PRINT*,'CREATE LINKEDLIST',STRING
      CALL LINKEDLIST$NEW(LL_STP)
!PRINT*,'READ',NFIL
      CALL LINKEDLIST$READ(LL_STP,NFIL,'MONOMER')
      TCHK=.TRUE.
!
!     ======================================================================  
!     == READ SOME DATA                                                   ==  
!     ======================================================================  
!     == LNX
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
      CALL LINKEDLIST$GET(LL_STP,'LNX',0,LNX)
!     == GRID
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'GRID')
      CALL LINKEDLIST$GET(LL_STP,'TYPE',0,TYPE)
      CALL LINKEDLIST$GET(LL_STP,'R1',0,R1)
      CALL LINKEDLIST$GET(LL_STP,'DEX',0,DEX)
      CALL LINKEDLIST$GET(LL_STP,'NR',0,NR)
      CALL RADIAL$NEW(TRIM(TYPE),GID)
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
!     == NB,NC
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
      CALL LINKEDLIST$GET(LL_STP,'NB',0,NB)
      CALL LINKEDLIST$GET(LL_STP,'NC',0,NC)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE SETUPREAD$GETI4(ID,VAL)
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **********************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      IF(ID.EQ.'LNX') THEN
        VAL=LNX
      ELSE IF(ID.EQ.'NR') THEN
        VAL=NR
      ELSE IF(ID.EQ.'GID') THEN
        VAL=GID
      ELSE IF(ID.EQ.'NC') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        CALL LINKEDLIST$GET(LL_STP,'NC',0,VAL)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUPREAD$GETI4')
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE SETUPREAD$GETR8(ID,VAL)
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     **********************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      IF(ID.EQ.'AEZ') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'GENERIC')
        CALL LINKEDLIST$GET(LL_STP,'AEZ',0,VAL)
      ELSE IF(ID.EQ.'PSZ') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'PSZ',0,VAL)
      ELSE IF(ID.EQ.'RCSM') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'RCSM',0,VAL)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUPREAD$GETR8')
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE SETUPREAD$GETI4A(ID,LEN,VAL)
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)  ,ALLOCATABLE:: IWORK(:)
!     **********************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
!
!     ==  ANGULAR MOMENTA OF THE PARTIAL WAVES  ==========================
      IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.LNX)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETI4A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'LOX',0,VAL)
!
!     ==  ANGULAR MOMENTUM OF THE CORE STATES =============================
      ELSE IF(ID.EQ.'LOFC') THEN
        IF(LEN.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$STOP('SETUPREAD$GETI4A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        ALLOCATE(IWORK(NB))
        CALL LINKEDLIST$GET(LL_STP,'L',0,IWORK)
        VAL(:)=IWORK(1:NC)
        DEALLOCATE(IWORK)
!
!     == #(NODES) OF THE CORE STATES =====================================
      ELSE IF(ID.EQ.'NNOFC') THEN
        IF(LEN.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$STOP('SETUPREAD$GETI4A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        ALLOCATE(IWORK(NB))
        CALL LINKEDLIST$GET(LL_STP,'NN',0,IWORK)
        VAL(:)=IWORK(1:NC)
        DEALLOCATE(IWORK)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUPREAD$GETI4A')
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE SETUPREAD$GETR8A(ID,LEN,VAL)
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      REAL(8)     ,ALLOCATABLE:: WORK(:)
      INTEGER(4)              :: I,I1,I2
!     **********************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
!
      IF(ID.EQ.'PRO') THEN
        IF(LEN.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        DO I=1,LNX
          I1=NR*(I-1)+1
          I2=NR*I
          CALL LINKEDLIST$GET(LL_STP,'PRO',I,VAL(I1:I2))
        ENDDO
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LEN.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        DO I=1,LNX
          I1=NR*(I-1)+1
          I2=NR*I
          CALL LINKEDLIST$GET(LL_STP,'AEPHI',I,VAL(I1:I2))
        ENDDO
      ELSE IF(ID.EQ.'PSPHI') THEN
        IF(LEN.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        DO I=1,LNX
          I1=NR*(I-1)+1
          I2=NR*I
          CALL LINKEDLIST$GET(LL_STP,'PSPHI',I,VAL(I1:I2))
        ENDDO
      ELSE IF(ID.EQ.'AECORE') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'AECORE',0,VAL)
      ELSE IF(ID.EQ.'PSCORE') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'PSCORE',0,VAL)
!
!     =================================================================
      ELSE IF(ID.EQ.'VADD') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'VADD',0,VAL)
!
!     == ONE-CENTER DIFFERENCE KINETIC ENERGY MATRIX ==================
      ELSE IF(ID.EQ.'DT') THEN
        IF(LEN.NE.LNX*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.LNX**2)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'DT',0,VAL)
!
!     == ONE-CENTER DIFFERENCE OVERLAP MATRIX ========================
      ELSE IF(ID.EQ.'DO') THEN
        IF(LEN.NE.LNX*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.LNX**2)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'DO',0,VAL)
!
!     == ONE-CENTER DIFFERENCE HAMILTONIAN OF THE ATOM ==================
      ELSE IF(ID.EQ.'DH') THEN
        IF(LEN.NE.LNX*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.LNX**2)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'DH',0,VAL)
      ELSE IF(ID.EQ.'AEPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        CALL LINKEDLIST$GET(LL_STP,'AEPOT',0,VAL)
!
!     == ENERGIES OF THE CORE STATES ====================================
      ELSE IF(ID.EQ.'EOFC') THEN
        IF(LEN.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        ALLOCATE(WORK(NB))
        CALL LINKEDLIST$GET(LL_STP,'E',0,WORK)
        VAL(:)=WORK(1:NC)
        DEALLOCATE(WORK)
!
!     == OCCUPATIONS OF THE CORE STATES =================================
      ELSE IF(ID.EQ.'FOFC') THEN
        IF(LEN.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        ALLOCATE(WORK(NB))
        CALL LINKEDLIST$GET(LL_STP,'OCC',0,WORK)
        VAL(:)=WORK(1:NC)
        DEALLOCATE(WORK)
!
!     == nodeless partial waves ===================================================
      ELSE IF(ID.EQ.'NDLSPHI') THEN
        IF(LEN.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        DO I=1,LNX
          I1=NR*(I-1)+1
          I2=NR*I
          CALL LINKEDLIST$GET(LL_STP,'NDLSPHI',I,VAL(I1:I2))
        ENDDO
!
!     == KINETIC ENERGY OF NODELESS PARTIAL WAVES =============================
      ELSE IF(ID.EQ.'NDLSTPHI') THEN
        IF(LEN.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LEN.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        DO I=1,LNX
          I1=NR*(I-1)+1
          I2=NR*I
          CALL LINKEDLIST$GET(LL_STP,'NDLSTPHI',I,VAL(I1:I2))
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUPREAD$GETR8A')
      END IF
      RETURN
      END
