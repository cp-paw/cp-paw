!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!
!........1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUPREAD$NEW(NFIL_,TCHK)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL_
      LOGICAL(4)  ,INTENT(OUT):: TCHK
      CHARACTER(6)            :: STRING
      CHARACTER(16)           :: TYPE
      REAL(8)                 :: R1
      REAL(8)                 :: DEX
!     **********************************************************************
      TCHK=.FALSE.
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUPREAD$GETI4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
!
!     ==========================================================================
!     == #(PARTIAL WAVES)                                                     ==
!     ==========================================================================
      IF(ID.EQ.'LNX') THEN
        VAL=LNX
!
!     ==========================================================================
!     == #(RADIAL GRID POINTS)                                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NR') THEN
        VAL=NR
!
!     ==========================================================================
!     == GRID ID                                                              ==
!     ==========================================================================
      ELSE IF(ID.EQ.'GID') THEN
        VAL=GID
!
!     ==========================================================================
!     == NUMBER OF CORE STATES                                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NC') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        CALL LINKEDLIST$GET(LL_STP,'NC',0,VAL)
!
!     ==========================================================================
!     == NUMBER OF ATOMIC STATES AVAILABLE                                    ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NB') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        CALL LINKEDLIST$GET(LL_STP,'NB',0,VAL)
!
!     ==========================================================================
!     ==  WRONG ID                                                            ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUPREAD$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUPREAD$GETR8(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUPREAD$GETI4A(ID,LENG,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      INTEGER(4)  ,INTENT(OUT):: VAL(LENG)
      INTEGER(4)  ,ALLOCATABLE:: IWORK(:)
!     **********************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
!
!     ==  ANGULAR MOMENTA OF THE PARTIAL WAVES  ==========================
      IF(ID.EQ.'LOX') THEN
        IF(LENG.NE.LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETI4A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'LOX',0,VAL)
!
!     ==  ANGULAR MOMENTUM OF THE CORE STATES =============================
      ELSE IF(ID.EQ.'LOFC') THEN
        IF(LENG.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
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
        IF(LENG.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUPREAD$GETR8A(ID,LENG,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE SETUPREAD_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      REAL(8)     ,INTENT(OUT):: VAL(LENG)
      REAL(8)     ,ALLOCATABLE:: WORK(:)
      INTEGER(4)              :: I,I1,I2
      INTEGER(4)  ,ALLOCATABLE:: LOFi(:)
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
!
!     ==========================================================================
!     == PROJECTOR FUNCTIONS                                                  ==
!     ==========================================================================
      IF(ID.EQ.'PRO') THEN
        IF(LENG.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
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
!
!     ==========================================================================
!     == ALL-ELECTRON PARTIAL WAVES                                           ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LENG.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
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
!
!     ==========================================================================
!     == PSEUDO PARTIAL WAVES                                                 ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSPHI') THEN
        IF(LENG.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
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
!
!     ==========================================================================
!     == ALL-ELECTRON CORE DENSITY                                            ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AECORE') THEN
        IF(LENG.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'AECORE',0,VAL)
!
!     ==========================================================================
!     == PSEUDO CORE DENSITY                                                  ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSCORE') THEN
        IF(LENG.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'PSCORE',0,VAL)
!
!     ==========================================================================
!     == POTENTIAL VADD                                                       ==
!     ==========================================================================
      ELSE IF(ID.EQ.'VADD') THEN
        IF(LENG.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'VADD',0,VAL)
!
!     ==========================================================================
!     == ONE-CENTER DIFFERENCE KINETIC ENERGY MATRIX                          ==
!     ==========================================================================
      ELSE IF(ID.EQ.'DT') THEN
        IF(LENG.NE.LNX*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX**2)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'DT',0,VAL)
!
!     ==========================================================================
!     == ONE-CENTER DIFFERENCE OVERLAP MATRIX                                 ==
!     ==========================================================================
      ELSE IF(ID.EQ.'DO') THEN
        IF(LENG.NE.LNX*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX**2)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'DO',0,VAL)
!
!     ==========================================================================
!     == ONE-CENTER DIFFERENCE HAMILTONIAN OF THE ATOM                        ==
!     ==========================================================================
      ELSE IF(ID.EQ.'DH') THEN
        IF(LENG.NE.LNX*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.LNX**2)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('LNX',LNX)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
        CALL LINKEDLIST$GET(LL_STP,'DH',0,VAL)
!
!     ==========================================================================
!     == ALL-ELECTRON ATOMIC POTENTIAL                                        ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPOT') THEN
        IF(LENG.NE.NR) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        CALL LINKEDLIST$GET(LL_STP,'AEPOT',0,VAL)
!
!     ==========================================================================
!     == ENERGIES OF THE CORE STATES                                          ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EOFC') THEN
        IF(LENG.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        ALLOCATE(WORK(NB))
        CALL LINKEDLIST$GET(LL_STP,'E',0,WORK)
        VAL(:)=WORK(1:NC)
        DEALLOCATE(WORK)
!
!     ==========================================================================
!     == OCCUPATIONS OF THE CORE STATES                                       ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FOFC') THEN
        IF(LENG.NE.NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        ALLOCATE(WORK(NB))
        CALL LINKEDLIST$GET(LL_STP,'OCC',0,WORK)
        VAL(:)=WORK(1:NC)
        DEALLOCATE(WORK)
!
!     ==========================================================================
!     == NODELESS PARTIAL WAVES                                               ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NDLSPHI') THEN
        IF(LENG.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
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
!     ==========================================================================
!     == KINETIC ENERGY OF NODELESS PARTIAL WAVES                             ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NDLSTPHI') THEN
        IF(LENG.NE.NR*LNX) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*LNX(')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
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
!
!     ==========================================================================
!     == READ SET OF NODELESS ATOMIC CORE STATES                              ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPSICORE') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
        IF(LENG.NE.NR*NC) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT (LENG.NE.NR*NC)')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$I4VAL('NR',NR)
          CALL ERROR$I4VAL('Nc',Nc)
          CALL ERROR$STOP('SETUPREAD$GETR8A')
        END IF
        DO I=1,NC
          I1=1+(I-1)*NR
          I2=I1-1+NR
          CALL LINKEDLIST$GET(LL_STP,'PHI-NODELESS',I,VAL(I1:I2))
        ENDDO
        ALLOCATE(LOFi(NC))
        CALL SETUPREAD$GETI4A('LOFC',NC,LOFi)
        CALL SETUPREAD_NDLSTONODAL(GID,NR,NC,LOFI,VAL)
        DEALLOCATE(LOFi)
!
!     ==========================================================================
!     == WRONG ID                                                             ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUPREAD$GETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUPREAD_NDLSTONODAL(GID,NR,NC,LOFI,PHI)
!     **************************************************************************
!     ** TRANSFORM NODELESS CORES STATES INTO CORRECT ALL-ELECTRON CORE STATES**
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)    :: GID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: NC
      INTEGER(4)  ,INTENT(IN)    :: LOFI(NC)
      REAL(8)     ,INTENT(INOUT) :: PHI(NR,NC)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: VAL
      INTEGER(4)                 :: IB1,IB2
!     **************************************************************************
      call radial$r(gid,nr,r)
      DO IB1=1,NC
!
!       ========================================================================
!       == NORMALIZE                                                          ==
!       ========================================================================
        AUX(:)=R(:)**2*PHI(:,IB1)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
        PHI(:,IB1)=PHI(:,IB1)/SQRT(VAL)
!
!       ========================================================================
!       == ORTHOGONALIZE HIGHER STATES                                        ==
!       ========================================================================
        DO IB2=IB1+1,NC
          IF(LOFI(IB2).NE.LOFI(IB1)) CYCLE
          AUX(:)=R(:)**2*PHI(:,IB1)*PHI(:,IB2)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
          PHI(:,IB2)=PHI(:,IB2)-PHI(:,IB1)*VAL
        ENDDO
      ENDDO
      RETURN
      END


