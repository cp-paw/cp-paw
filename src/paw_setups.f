!=======================================================================
!=======================================================================
!=======================================================================
!====      FORMFACTORS                                              ====
!=======================================================================
!=======================================================================
!
!.......................................................................
MODULE SETUP_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: SETUP                                                      **
!**                                                                   **
!**  PURPOSE:                                                         **
!**                                                                   **
!**  READS PSEUDOPOTENTIALS AND ATOMIC PSEUDO-WAVEFUNCTIONS           **
!**  GIVEN ON A LINEAR, RADIAL MESH                                   **
!**              AND                                                  **
!**  CALCULATES THE FORMFACTORS FOR LOCAL AND NONLOCAL                **
!**  CONTRIBUTIONS TO THE PSEUDOPOTENTIAL                             **
!**  IN PLANE WAVE REPRESENTATION                                     **
!**                                                                   **
!**  THE LOCAL COMPONENT CONTAINS ALSO THE POTENTIAL OF A             **
!**  (GAUSSIAN SHAPED) CHARGEDENSITY, WHICH COMPENSATES               **
!**  THE IONIC CHARGE.                                                **
!**                                                                   **
!**                                                                   **
!**  RCBG (R) = -ZATOM/(SQRT(PI)*RCRHO)**3 * EXP(-(R/RCRHO)**2)       **
!**                                                                   **
!**  OUTPUT:                                                          **
!**    RCRHO     RADIUS OF GAUSSIAN IN RHOPS                          **
!**    RHOPS     COMPENSATION CHARGE DENSITY                          **
!**    VLOC      LOCAL CONTRIBUTION TO PSEUDOPOTENTIAL                **
!**    WNL       NONLOCAL CONTRIBUTION (PROJECTOR)                    **
!**                                                                   **
!**          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)      **
!***********************************************************************
TYPE THIS_TYPE
INTEGER(4)             :: I
CHARACTER(32)          :: ID
REAL(8)                :: AEZ
REAL(8)                :: PSZ
REAL(8)                :: RCBG
REAL(8)                :: RCSM
INTEGER(4)             :: LNX
INTEGER(4)             :: LMNX
INTEGER(4)             :: LMRX
INTEGER(4),POINTER     :: LOX(:)       !(LNXX)
REAL(8)   ,POINTER     :: VADD(:)      !(NRX)
REAL(8)   ,POINTER     :: AECORE(:)    !(NRX)
REAL(8)   ,POINTER     :: PSCORE(:)    !(NRX)
REAL(8)   ,POINTER     :: PRO(:,:)     !(NRX,LNXX)
REAL(8)   ,POINTER     :: AEPHI(:,:)   !(NRX,LNXX)
REAL(8)   ,POINTER     :: PSPHI(:,:)   !(NRX,LNXX)
REAL(8)   ,POINTER     :: DTKIN(:,:)   !(LNXX,LNXX)
REAL(8)   ,POINTER     :: DOVER(:,:)   !(LNXX,LNXX)
REAL(8)   ,POINTER     :: VADDOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: PSCOREOFG(:) !(NGX)
REAL(8)   ,POINTER     :: VHATOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: NHATPRIMEOFG(:)  !(NGX)
REAL(8)   ,POINTER     :: PROOFG(:,:)  !(NRX,LNXX)
REAL(8)                :: M
REAL(8)                :: ZV
REAL(8)                :: PSG2
REAL(8)                :: PSG4
CHARACTER(16)          :: FILEID
TYPE(THIS_TYPE),POINTER:: NEXT
END TYPE THIS_TYPE
!
REAL(8)                :: R1 =1.056D-4
REAL(8)                :: DEX=5.D-2
INTEGER(4)             :: NR =250
INTEGER(4)             :: NRX=250
REAL(8)   ,PARAMETER   :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
INTEGER(4),PARAMETER   :: NG=256
REAL(8)                :: G1          ! FIRST POINT ON THE RADIAL G-GRID
INTEGER(4)             :: NSP=0
INTEGER(4)             :: LMRXX=0
INTEGER(4)             :: LMNXX=0
INTEGER(4)             :: LNXX=0
TYPE(THIS_TYPE),POINTER :: FIRST
TYPE(THIS_TYPE),POINTER :: THIS
TYPE FASTACCESS_TYPE  
  TYPE(THIS_TYPE),POINTER :: THIS
END TYPE FASTACCESS_TYPE
TYPE(FASTACCESS_TYPE),ALLOCATABLE :: FASTACCESS(:)
END MODULE SETUP_MODULE
!
!     ..................................................................
      SUBROUTINE SETUP$ISELECT(I)
!     ******************************************************************
!     **  SELECTS A SETUP PER INTEGER INDEX                           **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: I
      INTEGER(4)            :: J
!     ******************************************************************
      IF(.NOT.ALLOCATED(FASTACCESS)) THEN
        IF(.NOT.ASSOCIATED(FIRST)) THEN
          CALL ERROR$MSG('NO SETUPS DEFINED')
          CALL ERROR$STOP('SETUP$ISELECT')
        END IF
        THIS=>FIRST
        NSP=1
        DO WHILE(ASSOCIATED(THIS%NEXT))
          NSP=NSP+1
          THIS=>THIS%NEXT
        ENDDO
        ALLOCATE(FASTACCESS(NSP))
        FASTACCESS(1)%THIS=>FIRST
        DO J=2,NSP
          FASTACCESS(J)%THIS=>FASTACCESS(J-1)%THIS%NEXT 
        ENDDO
      END IF
!
      IF(I.GT.NSP) THEN
        CALL ERROR$MSG('INDEX I OUT OF RANGE')
        CALL ERROR$I4VAL('I',I)
        CALL ERROR$I4VAL('NSP',NSP)
        CALL ERROR$STOP('SETUP$ISELECT')
      END IF
!
      THIS=>FASTACCESS(I)%THIS
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$SELECT(ID)
!     ******************************************************************
!     **  SELECTS A SETUP PER ID                                      **
!     **  AND CREATES A NEW, IF IT DOES NOT EXIST                     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
!     ******************************************************************
!
!     == CHECK IF ALREADY SELECTED
      IF(ASSOCIATED(THIS)) THEN
        IF(THIS%ID.EQ.ID) RETURN
      END IF
!
!     == CHECK IF PRESENT ===
      IF(ASSOCIATED(FIRST)) THEN
        THIS=>FIRST
        DO 
          IF(THIS%ID.EQ.ID) RETURN
          IF(.NOT.ASSOCIATED(THIS%NEXT))EXIT
          THIS=>THIS%NEXT
        ENDDO
        ALLOCATE(THIS%NEXT)
        THIS%NEXT%I=THIS%I+1
        THIS=>THIS%NEXT
      ELSE
        ALLOCATE(FIRST)
        THIS=>FIRST
        THIS%I=1
      END IF
!
!     == CREATE NEW
      IF(ALLOCATED(FASTACCESS)) DEALLOCATE(FASTACCESS)
      THIS%ID    =ID
      THIS%AEZ   =0.D0
      THIS%PSZ   =0.D0
      THIS%RCBG  =0.D0
      THIS%RCSM  =0.D0
      THIS%LNX   =0
      THIS%LMNX  =0
      THIS%LMRX  =0
      THIS%PSG2  =0.D0
      THIS%PSG4  =0.D0
      NULLIFY(THIS%LOX)     !(LNX)
      NULLIFY(THIS%VADD)    !(NRX)
      NULLIFY(THIS%AECORE)  !(NRX)
      NULLIFY(THIS%PSCORE)  !(NRX)
      NULLIFY(THIS%PRO)     !(NRX,LNX)
      NULLIFY(THIS%AEPHI)   !(NRX,LNX)
      NULLIFY(THIS%PSPHI)   !(NRX,LNX)
      NULLIFY(THIS%DTKIN)   !(LNXX,LNX)
      NULLIFY(THIS%DOVER)   !(LNXX,LNX)
      NULLIFY(THIS%VADDOFG) !(NGX)
      NULLIFY(THIS%PSCOREOFG) !(NGX)
      NULLIFY(THIS%VHATOFG) !(NGX)
      NULLIFY(THIS%NHATPRIMEOFG) !(NGX)
      NULLIFY(THIS%PROOFG)  !(NGX,LNX)
      NULLIFY(THIS%NEXT)
      WRITE(THIS%FILEID,*)THIS%I
      THIS%FILEID='ATOM'//ADJUSTL(THIS%FILEID)
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE SETUP$GETCH(ID,VAL)
!     ******************************************************************
!     **  COLLECTS INTERNAL DATA                                      **
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      CHARACTER(*),INTENT(OUT) :: VAL
!     ******************************************************************
      IF(ID.EQ.'ID') THEN
        VAL=THIS%ID
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SETUP$GETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$GETI4(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(OUT) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NR') THEN
        VAL=NR
      ELSE IF(ID.EQ.'NRX') THEN
        VAL=NRX
      ELSE IF(ID.EQ.'LNX') THEN
        VAL=THIS%LNX
      ELSE IF(ID.EQ.'LMNX') THEN
        VAL=THIS%LMNX
      ELSE IF(ID.EQ.'LMRX') THEN
        VAL=THIS%LMRX
      ELSE IF(ID.EQ.'LMRXX') THEN
        VAL=LMRXX
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$GETI4A(ID,LEN,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(OUT) :: VAL(len)
!     ******************************************************************
      IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%LOX
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$GETR8(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VAL
!     ******************************************************************
      IF(ID.EQ.'AEZ') THEN
        VAL=THIS%AEZ
      ELSE IF(ID.EQ.'RCSM') THEN
        VAL=THIS%RCSM
      ELSE IF(ID.EQ.'RCBG') THEN
        VAL=THIS%RCBG
      ELSE IF(ID.EQ.'ZV') THEN
        VAL=THIS%ZV
      ELSE IF(ID.EQ.'M') THEN
        VAL=THIS%M
      ELSE IF(ID.EQ.'<G2>') THEN
        VAL=THIS%PSG2
      ELSE IF(ID.EQ.'<G4>') THEN
        VAL=THIS%PSG4
      ELSE IF(ID.EQ.'R1') THEN
        VAL=R1
      ELSE IF(ID.EQ.'DEX') THEN
        VAL=DEX
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$GETR8A(ID,LEN,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      REAL(8)     ,INTENT(OUT) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'PRO') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PRO,(/LEN/))
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPHI,(/LEN/))
      ELSE IF(ID.EQ.'PSPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PSPHI,(/LEN/))
      ELSE IF(ID.EQ.'AECORE') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%AECORE
      ELSE IF(ID.EQ.'PSCORE') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%PSCORE
      ELSE IF(ID.EQ.'VADD') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%VADD
      ELSE IF(ID.EQ.'DEKIN') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DTKIN,(/LEN/))
      ELSE IF(ID.EQ.'DO') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DOVER,(/LEN/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8A')
      END IF
      RETURN
      END  

!
!     ..................................................................
      SUBROUTINE SETUP$GETFOFG(ID,TDER,IND,NG_,G2,CELLVOL,F)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID      !identifier
      LOGICAL(4)  ,INTENT(IN)  :: TDER    ! calculate radial dervative
      INTEGER(4)  ,INTENT(IN)  :: IND     ! selector (used only for id=pro)
      INTEGER(4)  ,INTENT(IN)  :: NG_     ! #(plane waves)
      REAL(8)     ,INTENT(IN)  :: G2(NG_) ! G**2
      REAL(8)     ,INTENT(OUT) :: F(NG_)  
      REAL(8)                  :: FOFG(NG)
      REAL(8)                  :: PI
      INTEGER(4)               :: IG
      INTEGER(4)               :: NGAMMA
      REAL(8)                  :: CELLVOL
      REAL(8)                  :: G
!     ******************************************************************
      if(ng_.eq.0) return
      IF(ID.EQ.'PRO') THEN
        IF(IND.LT.1.OR.IND.GT.THIS%LNX) THEN
          CALL ERROR$MSG('LN OUT OF RANGE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('IND',IND)
          CALL ERROR$I4VAL('LNX',THIS%LNX)
          CALL ERROR$STOP('SETUP$GETFOFG')
        END IF
        FOFG(:)=THIS%PROOFG(:,IND)
      ELSE IF(ID.EQ.'PSCORE') THEN
        FOFG(:)=THIS%PSCOREOFG(:)
      ELSE IF(ID.EQ.'VADD') THEN
        FOFG(:)=THIS%VADDOFG(:)
      ELSE IF(ID.EQ.'V0') THEN
        FOFG(:)=THIS%VHATOFG(:)
      ELSE IF(ID.EQ.'G0') THEN
        FOFG(:)=THIS%NHATPRIMEOFG(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETFOFG')
      END IF
!
!     ==================================================================
!     == INTERPOLATE VALUES FROM RADIAL GRID
!     ==================================================================
      NGAMMA=0
      IF(TDER) THEN
        G=DSQRT(G2(1))
        CALL RADIAL$DERIVATIVE(G1,DEX,NG,FOFG,G,F(1))
        F(1)=G*F(1)
        DO IG=2,NG_
          IF(DABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
            F(IG) =F(IG-1)
          ELSE
            G=DSQRT(G2(IG))
            IF(G.LT.1.D-6) NGAMMA=IG
            CALL RADIAL$DERIVATIVE(G1,DEX,NG,FOFG,G,F(IG))
            F(IG)=G*F(IG)
          END IF
        ENDDO
      ELSE
        G=DSQRT(G2(1))
        CALL RADIAL$VALUE(G1,DEX,NG,FOFG,G,F(1))
        DO IG=2,NG_
          IF(DABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
            F(IG) =F(IG-1)
          ELSE
            G=DSQRT(G2(IG))
            IF(G.LT.1.D-6) NGAMMA=IG
            CALL RADIAL$VALUE(G1,DEX,NG,FOFG,G,F(IG))
          END IF
        ENDDO
      END IF
!
!     ==================================================================
!     == CORRECT EXTRAPOLATION TO THE GAMMA POINT                     ==
!     ==================================================================
      IF(NGAMMA.NE.0) THEN
        PI=4.D0*DATAN(1.D0)
        IF(TDER) THEN
          NGAMMA=0.D0
        ELSE
          IF(ID.EQ.'g0') THEN 
            F(NGAMMA)=4.D0*PI
          ELSE IF(ID.EQ.'V0') THEN
            F(NGAMMA)=PI*(THIS%RCBG**2-THIS%RCSM**2)*4.D0*PI
          END IF
        END IF
      END IF
!
!     ==================================================================
!     == DIVIDE BY CELLVOL                                            ==
!     ==================================================================
      F=F/CELLVOL
      RETURN
      END

!
!     ..................................................................
      SUBROUTINE SETUP_READ
!     ******************************************************************
!     **  READ SELECTED SETUP                                         **
!     **  REQUIRES INFORMATION FROM ATOMTYPELIST                      **
!     **    NAME; LRHOX                                               **
!     **  REQUIRES THE FILEHANDLER TO KNOW THE SETUP FILE             **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      REAL(8)   ,PARAMETER  :: TOL=1.D-6
      INTEGER(4)            :: LMRXCUT
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: ISP
      INTEGER(4)            :: NFIL
      REAL(8)               :: RI
      INTEGER(4)            :: IR
      INTEGER(4)            :: LN
      INTEGER(4)            :: IRCCOR
      REAL(8)               :: XEXP
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: IRMAX
      CHARACTER(16)         :: NAME
      INTEGER(4)            :: L,LX,ISVAR,LNOLD,LNX
      INTEGER(4)            :: ln1,ln2,ln1a,ln2a
      INTEGER(4),ALLOCATABLE:: NPRO(:)
      INTEGER(4),ALLOCATABLE:: IWORK(:)
      REAL(8)   ,ALLOCATABLE:: DWORK(:,:,:) 
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: ERROR       ! ERROR ESTIMATE
      REAL(8)   ,ALLOCATABLE:: FOFR1(:)   
      REAL(8)   ,ALLOCATABLE:: FOFG1(:)
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP$READ')
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('SETUP$READ')
      END IF
!
      CALL ATOMTYPELIST$NAME(THIS%I,THIS%ID)
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETR8('M',THIS%M)
      CALL ATOMTYPELIST$GETR8('ZV',THIS%ZV)
      CALL ATOMTYPELIST$GETR8('PS<G2>',THIS%PSG2)
      CALL ATOMTYPELIST$GETR8('PS<G4>',THIS%PSG4)
!
      CALL FILEHANDLER$UNIT(THIS%FILEID,NFIL)
!
      CALL INPOT$LNX(NFIL,LNX)
      THIS%LNX=LNX
      ALLOCATE(THIS%LOX(LNX))
      ALLOCATE(THIS%VADD(NRX))
      ALLOCATE(THIS%AECORE(NRX))
      ALLOCATE(THIS%PSCORE(NRX))
      ALLOCATE(THIS%PRO(NRX,LNX))
      ALLOCATE(THIS%AEPHI(NRX,LNX))
      ALLOCATE(THIS%PSPHI(NRX,LNX))
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      THIS%VADD=0.D0
      THIS%AECORE=0.D0
      THIS%PSCORE=0.D0
      THIS%PRO=0.D0
      THIS%AEPHI=0.D0
      THIS%PSPHI=0.D0
!     
!     ==================================================================
!     ==  READ PSEUDOPOTENTIALS AND PSEUDO WAVE FUNCTIONS             ==
!     ==================================================================
                            CALL TRACE$PASS('READ SETUP FILES')
      THIS%RCBG=1.D0/DSQRT(0.218D0)
      
      CALL INPOT$READALL(NFIL,NRX,R1,DEX,NR,THIS%LNX,THIS%LOX &
     &         ,THIS%AEZ,THIS%PSZ,THIS%PSPHI,THIS%AEPHI &
     &         ,THIS%VADD,THIS%RCSM,THIS%DTKIN,THIS%DOVER &
     &         ,IRCCOR,THIS%AECORE,THIS%PSCORE,THIS%PRO)
      CALL FILEHANDLER$CLOSE(THIS%FILEID)
!     
!     ==================================================================
!     == LIMIT NUMBER OF PROJECTORS FOR EACH L                        ==
!     ==================================================================
      LX=0
      DO LN=1,THIS%LNX
        LX=MAX(LX,THIS%LOX(LN))
      ENDDO
      ALLOCATE(NPRO(LX+1)) 
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETI4A('NPRO',LX+1,NPRO)
      DO L=0,LX
        ISVAR=0
        DO LN=1,THIS%LNX
          IF(THIS%LOX(LN).NE.L)CYCLE
          ISVAR=ISVAR+1
          IF(ISVAR.GT.NPRO(L+1)) THEN
            THIS%LOX(LN)=-1     ! MARK PROJECTORS TO BE DELETED BY LOX=-1
            ISVAR=ISVAR-1
          END IF
        ENDDO
        NPRO(L+1)=ISVAR
      ENDDO
      CALL ATOMTYPELIST$SETI4A('NPRO',LX+1,NPRO)
      DEALLOCATE(NPRO)
!
      LNOLD=THIS%LNX
      LNX=0
      DO LN=1,LNOLD
        IF(THIS%LOX(LN).NE.-1) LNX=LNX+1
      ENDDO
      THIS%LNX=LNX
!
!     == FOLD DOWN ARRAYS FOR PROJECTORS AND PARTIALWAVES, LOX =========
      ALLOCATE(DWORK(NRX,LNOLD,3))
      ALLOCATE(IWORK(LNOLD))
      DWORK(:,:,1)=THIS%PRO(:,:)
      DWORK(:,:,2)=THIS%AEPHI(:,:)
      DWORK(:,:,3)=THIS%PSPHI(:,:)
      IWORK(:)=THIS%LOX(:)
      DEALLOCATE(THIS%PRO)
      DEALLOCATE(THIS%AEPHI)
      DEALLOCATE(THIS%PSPHI)
      DEALLOCATE(THIS%LOX)
      ALLOCATE(THIS%PRO(NRX,THIS%LNX))
      ALLOCATE(THIS%AEPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%PSPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%LOX(THIS%LNX))
      ISVAR=0
      DO LN=1,LNOLD
        IF(IWORK(LN).EQ.-1) CYCLE
        ISVAR=ISVAR+1
        THIS%PRO(:,ISVAR)=DWORK(:,LN,1)
        THIS%AEPHI(:,ISVAR)=DWORK(:,LN,2)
        THIS%PSPHI(:,ISVAR)=DWORK(:,LN,3)
        THIS%LOX(ISVAR)=IWORK(LN)
      ENDDO
      DEALLOCATE(DWORK)
!
!     == FOLD DOWN ARRAYS FOR DTKIN AND DOVER ==========================
      ALLOCATE(DWORK(LNOLD,LNOLD,2))
      DWORK(:,:,1)=THIS%DTKIN(:,:)
      DWORK(:,:,2)=THIS%DOVER(:,:)
      DEALLOCATE(THIS%DTKIN)
      DEALLOCATE(THIS%dOVER)
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%dOVER(LNX,LNX))
      LN1A=0
      DO LN1=1,LNOLD
        IF(IWORK(LN1).EQ.-1) CYCLE
        LN1A=LN1A+1
        LN2A=0
        DO LN2=1,LNOLD
          IF(IWORK(LN2).EQ.-1) CYCLE
          LN2A=LN2A+1
          THIS%DTKIN(LN1A,LN2A)=DWORK(LN1,LN2,1)
          THIS%dOVER(LN1A,LN2A)=DWORK(LN1,LN2,2)
        ENDDO
      ENDDO
      DEALLOCATE(DWORK)
!
      DEALLOCATE(IWORK)
!     
!     ==================================================================
!     == SET VALUES BEYOND A CERTAIN RADIUS EXACTLY TO ZERO           ==
!     ==================================================================
                            CALL TRACE$PASS('CHECK MAX. RADIUS')
      IRMAX=0
      DO IR=1,NR
        TCHK=(DABS(THIS%VADD(IR)).LT.TOL)
        TCHK=TCHK.AND.(DABS(THIS%PSCORE(IR)-THIS%AECORE(IR)).LT.TOL)
        DO LN=1,THIS%LNX
          TCHK=TCHK.AND. &
     &           (DABS(THIS%AEPHI(IR,LN)-THIS%PSPHI(IR,LN)).LT.TOL)
        ENDDO
        IF(.NOT.TCHK) IRMAX=IR
      ENDDO
      DO IR=IRMAX+1,NR
        THIS%VADD(IR)=0.D0
        DO LN=1,THIS%LNX
          THIS%AEPHI(IR,LN)=0.D0
          THIS%PSPHI(IR,LN)=0.D0
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  DEFINE ARRAYS                                             ==
!     ================================================================
                            CALL TRACE$PASS('DEFINE ARRAYS')
!
!     == SELECT NATURAL VALUES =======================================
      THIS%LMNX=0
      THIS%LMRX=0
      DO LN=1,THIS%LNX
        THIS%LMNX=THIS%LMNX+2*THIS%LOX(LN)+1
        THIS%LMRX=MAX(THIS%LMRX,(2*THIS%LOX(LN)+1)**2)
      ENDDO
!
!     == LIMIT MAX ANGULAR MOMENTUM FOR THE DENSITY TO MAX VALUE =======
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETI4('LRHOX',LMRXCUT)
      LMRXCUT=(LMRXCUT+1)**2
      THIS%LMRX=MIN(THIS%LMRX,LMRXCUT)
      CALL ATOMTYPELIST$UNSELECT
!     
!     ==================================================================
!     ==  UPDATE GLOBAL VARIABLES                                     ==
!     ==================================================================
      LMNXX=MAX(LMNXX,THIS%LMNX)
      LMRXX=MAX(LMRXX,THIS%LMRX)
      LNXX=MAX(LNXX,THIS%LNX)
!     
!     ==================================================================
!     ==  PERFORM BESSELTRANSFORMS                                    ==
!     ==================================================================
      PI=4.D0*DATAN(1.D0)
      FOURPI=4.D0*PI
      G1=GMAX*DEXP(-DEX*DBLE(NR-1))
      IF(NG.LT.NR) THEN
        CALL ERROR$STOP('SETUP_RADTOG')
      END IF
      ALLOCATE(FOFR1(NG))
      ALLOCATE(FOFG1(NG))
      FOFR1(NR+1:NG)=0.D0
!     == VADD (VBAR) ===================================================
      FOFR1(1:NR)=FOURPI*THIS%VADD(1:NR)
      CALL BESSELTRANSFORM(0,NG,R1,G1,DEX,FOFR1,FOFG1,ERROR)
      ALLOCATE(THIS%VADDOFG(NG))
      THIS%VADDOFG(:)=FOFG1(:)
!     == PSCORE (VBAR) =================================================
      FOFR1(1:NR)=FOURPI*THIS%PSCORE(1:NR)
      CALL BESSELTRANSFORM(0,NG,R1,G1,DEX,FOFR1,FOFG1,ERROR)
      ALLOCATE(THIS%PSCOREOFG(NG))
      THIS%PSCOREOFG(:)=FOFG1(:)
!     == PROJECTORS ====================================================
      ALLOCATE(THIS%PROOFG(NG,LNX))
      DO LN=1,LNX
        L=THIS%LOX(LN)
        FOFR1(1:NR)=FOURPI*THIS%PRO(1:NR,LN)
        CALL BESSELTRANSFORM(L,NG,R1,G1,DEX,FOFR1,FOFG1,ERROR)
        THIS%PROOFG(:,LN)=FOFG1(:)
      ENDDO
      DEALLOCATE(FOFR1)
      DEALLOCATE(FOFG1)
!     == COMPENSATION GAUSSIAN =========================================
      ALLOCATE(THIS%NHATPRIMEOFG(NG))
      ALLOCATE(THIS%VHATOFG(NG))
      CALL SETUP_COMPOFG(THIS%RCBG,THIS%RCSM,G1,DEX,NG &
     &                  ,THIS%NHATPRIMEOFG,THIS%VHATOFG)
!      
                            CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$REPORT(NFIL)
!     ******************************************************************
!     **  REPORT SETUP INFORMATION                                    **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NFIL
      INTEGER(4)                 :: L,NPRO,LN,NPROSUM
      TYPE(THIS_TYPE),POINTER    :: THIS1
      CHARACTER(32)              ::STRING
      REAL(8)                    :: U
      INTEGER(4)                 :: THISTASK,NTASKS
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL CONSTANTS$GET('U',U)
      THIS1=>FIRST
      DO 
        WRITE(NFIL,*)
        CALL REPORT$TITLE(NFIL,'ATOMIC SETUP: '//TRIM(THIS1%ID))
        CALL REPORT$R8VAL(NFIL,'ATOMIC NUMBER',THIS1%AEZ,' ')
        CALL REPORT$R8VAL(NFIL,'ATOMIC MASS  ',THIS1%M/U,'U')
        CALL REPORT$R8VAL(NFIL,'#(VALENCE ELECTRONS)',THIS1%ZV,' ')
        L=0
        NPROSUM=0
        DO WHILE (NPROSUM.LT.THIS1%LNX)
          NPRO=0
          DO LN=1,THIS1%LNX
            IF(THIS1%LOX(LN).EQ.L) NPRO=NPRO+1
          ENDDO
          IF(NPRO.NE.0) THEN
             WRITE(STRING,*)L
             STRING='#(PROJECTORS FOR L='//ADJUSTL(STRING)
             STRING=TRIM(STRING)//')'
             CALL REPORT$I4VAL(NFIL,TRIM(STRING),NPRO,' ')
          END IF
          L=L+1
          NPROSUM=NPROSUM+NPRO
        ENDDO        

        CALL REPORT$I4VAL(NFIL,'MAX #(ANGULAR MOMENTA (L,M) FOR 1C-DENSITY)' &
     &                        ,THIS1%LMRX,' ')
        CALL REPORT$R8VAL(NFIL,'GAUSSIAN DECAY FOR COMPENSATION DENSITY ' &
     &                        ,THIS1%RCSM,'ABOHR ')
        CALL REPORT$R8VAL(NFIL,'GAUSSIAN DECAY FOR EXTENDED COMPENSATION DENSITY' &
     &                        ,THIS1%RCBG,'ABOHR ')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G2> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG2,'H')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G4> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG4,'H')
        IF(.NOT.ASSOCIATED(THIS1%NEXT)) EXIT
        THIS1=>THIS1%NEXT
      ENDDO
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$READ
!     ******************************************************************
!     **  READ SETUP                                                  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)           :: NFILO
      INTEGER(4)           :: ISP,NSP1
      CHARACTER(32)        :: NAME
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP$READ')
!
!     ==================================================================
!     ==  CREATE SETUPS                                               ==
!     ==================================================================
      CALL ATOMTYPELIST$LENGTH(NSP1)
      DO ISP=1,NSP1
        CALL ATOMTYPELIST$NAME(ISP,NAME)
        CALL ATOMTYPELIST$SELECT(NAME)
        CALL ATOMTYPELIST$UNSELECT
        CALL SETUP$SELECT(NAME)
      ENDDO
!
!     ==================================================================
!     ==  READ SETUP FILES                                            ==
!     ==================================================================
      DO ISP=1,NSP1
        CALL SETUP$ISELECT(ISP)
        CALL SETUP_READ
      ENDDO
                            CALL TRACE$POP
      RETURN
      END

!
!     ..................................................................
      SUBROUTINE SETUP$RADIALGPRO(ISP_,LNX_,CELLVOL_,NG_,G2_,PROG_,DPROG_)
!     ******************************************************************
!     ** GAUSSIANS AND POTENTIAL FOR THE EWALD TRICK WITH THE         **
!     ** COMPENSATION CHARGE DENSIT TO THE DENSITY G-GRID             **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNX_
      REAL(8)   ,INTENT(IN) :: CELLVOL_
      INTEGER(4),INTENT(IN) :: NG_
      REAL(8)   ,INTENT(IN) :: G2_(NG_)
      REAL(8)   ,INTENT(OUT):: PROG_(NG_,LNX_)
      REAL(8)   ,INTENT(OUT):: DPROG_(NG_,LNX_)
      REAL(8)   ,PARAMETER  :: RMAX=5.D0
      INTEGER(4)            :: IG,LN,L
!     ******************************************************************
      CALL ERROR$MSG('MARKED FOR DELTION: DO NOT USE!')
      CALL ERROR$STOP('SETUP$RADIALGPRO')
                            CALL TRACE$PUSH('SETUP$GPROJECTORS')
!
!     ==================================================================
!     ==  SOME INITIAL CHECKS                                         ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP_)
      IF(LNX_.NE.THIS%LNX) THEN
        CALL ERROR$OVERFLOW('LNX_',LNX_,THIS%LNX)
        CALL ERROR$STOP('SETUP$GPROJECTORS')
      END IF

!     ==================================================================
!     ==  CALCULATE PROG USING A BESSELTRANSFORM                      ==
!     ==================================================================
      DO LN=1,LNX_
        DO IG=1,NG_
          PROG_(IG,LN)=0.D0
        ENDDO
        L=THIS%LOX(LN)
        CALL SETUP_RADTOG(L,R1,DEX,NRX,RMAX,THIS%PRO(1,LN) &
     &                 ,NG_,CELLVOL_,G2_,PROG_(1,LN),DPROG_(1,LN))
      ENDDO
                            CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$COMPENSATION(ISP_,NG_,G2_,CELLVOL_ &
     &                             ,G0_,V0_,DG0_,DV0_)
!     ******************************************************************
!     ** GAUSSIANS AND POTENTIAL FOR THE EWALD TRICK WITH THE         **
!     ** COMPENSATION CHARGE DENSIT TO THE DENSITY G-GRID             **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(IN) :: CELLVOL_
      INTEGER(4),INTENT(IN) :: NG_
      REAL(8)   ,INTENT(IN) :: G2_(NG_)
      REAL(8)   ,INTENT(OUT):: G0_(NG_)
      REAL(8)   ,INTENT(OUT):: V0_(NG_)
      REAL(8)   ,INTENT(OUT):: DG0_(NG_)
      REAL(8)   ,INTENT(OUT):: DV0_(NG_)
!     ******************************************************************
      CALL ERROR$MSG('MARKED FOR DELTION: DO NOT USE!')
      CALL ERROR$STOP('SETUP$COMPENSATION')
                              CALL TRACE$PUSH('SETUP$COMPENSATION')
!
!     ==================================================================
!     ==  SOME INITIAL CHECKS                                         ==
!     ==================================================================
      IF(ISP_.GT.NSP) THEN
        CALL ERROR$OVERFLOW('ISP_',ISP_,NSP)
        CALL ERROR$STOP('SETUP$COMPENSATION')
      END IF
!
!     ==================================================================
!     ==  CALCULATE G0,V0                                             ==
!     ==================================================================
               CALL TIMING$CLOCKON('FORMF-BESSOV COMPENSATION')
      CALL SETUP_COMPENSATIONGAUSSIANS(THIS%RCBG,THIS%RCSM &
     &          ,NG_,G2_,CELLVOL_,G0_,V0_,DG0_,DV0_)
             CALL TIMING$CLOCKOFF('FORMF-BESSOV COMPENSATION')
                            CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$PSCOREG(ISP_,NG_,G2_,CELLVOL_,PSCOREG_,DPSCOREG_)
!     ******************************************************************
!     **  RETURN PSCOREDENSITY                                        **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NG_
      REAL(8)   ,INTENT(IN) :: G2_(NG_)
      REAL(8)   ,INTENT(IN) :: CELLVOL_
      REAL(8)   ,INTENT(OUT):: PSCOREG_(NG_)
      REAL(8)   ,INTENT(OUT):: DPSCOREG_(NG_)
      REAL(8)   ,PARAMETER  :: RMAX=5.D0
!     ******************************************************************
      CALL ERROR$MSG('MARKED FOR DELTION: DO NOT USE!')
      CALL ERROR$STOP('SETUP$pscoreg')
                              CALL TRACE$PUSH('SETUP$PSCOREG')
      CALL SETUP$ISELECT(ISP_)
!
!     ==================================================================
!     ==  CALCULATE PSEUDOCOREDENSITY IN RECIPROCAL SPACE             ==
!     ==================================================================
                              CALL TIMING$CLOCKON('FORMF-BESSOV PSCORE')
      CALL SETUP_RADTOG(0,R1,DEX,NRX,RMAX,THIS%PSCORE,NG_,CELLVOL_,G2_ &
     &               ,PSCOREG_,DPSCOREG_)
                              CALL TIMING$CLOCKOFF('FORMF-BESSOV PSCORE')
                              CALL TRACE$POP
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$VBARG(ISP_,NG_,G2_,CELLVOL_,VBARG_,DVBARG_)
!     ******************************************************************
!     **  CALCULATE FORMFACTOR OF LOCAL PART :<VLOC|G>                **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NG_
      REAL(8)   ,INTENT(IN) :: G2_(NG_)
      REAL(8)   ,INTENT(IN) :: CELLVOL_
      REAL(8)   ,INTENT(OUT):: VBARG_(NG_)
      REAL(8)   ,INTENT(OUT):: DVBARG_(NG_)
      REAL(8)   ,PARAMETER  :: RMAX=5.D0
!     ******************************************************************
      CALL ERROR$MSG('MARKED FOR DELTION: DO NOT USE!')
      CALL ERROR$STOP('SETUP$vbarg')
                            CALL TRACE$PUSH('SETUP$VBARG')
             CALL TIMING$CLOCKON('FORMF-BESSOV VBAR')
      CALL SETUP$ISELECT(ISP_)
      CALL SETUP_RADTOG(0,R1,DEX,NRX,RMAX,THIS%VADD,NG_,CELLVOL_,G2_ &
     &               ,VBARG_,DVBARG_)
             CALL TIMING$CLOCKOFF('FORMF-BESSOV VBAR')
                            CALL TRACE$POP
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$AEPARTIALWAVES(ISP_,NRX_,LNX_,AEPHI_)
!     ******************************************************************
!     **  RETURN AE PARTIAL WAVES ON THE RADIAL GRID                  **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      INTEGER(4),INTENT(IN) :: LNX_
      REAL(8)   ,INTENT(OUT):: AEPHI_(NRX_,LNX_)
      INTEGER(4)            :: LN,IR
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN=1,LNX_
        DO IR=1,NR
          AEPHI_(IR,LN)=THIS%AEPHI(IR,LN)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$PSPARTIALWAVES(ISP_,NRX_,LNX_,PSPHI_)
!     ******************************************************************
!     **  RETURN PS PARTIAL WAVE ON A RADIAL GRID                     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      INTEGER(4),INTENT(IN) :: LNX_
      REAL(8)   ,INTENT(OUT):: PSPHI_(NRX_,LNX_)
      INTEGER(4)            :: IR,LN
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN=1,LNX_
        DO IR=1,NR
          PSPHI_(IR,LN)=THIS%PSPHI(IR,LN)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$AECORE(ISP_,NRX_,AECORE_)
!     ******************************************************************
!     **  RETURN AE CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: AECORE_(NRX_)
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      AECORE_(:)=THIS%AECORE(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$PSCORE(ISP_,NRX_,PSCORE_)
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: PSCORE_(NRX_)
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO IR=1,NR
        PSCORE_(IR)=THIS%PSCORE(IR)
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$VBAR(ISP_,NRX_,VBAR_)
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: VBAR_(NRX_)
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO IR=1,NR
        VBAR_(IR)=THIS%VADD(IR)
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$1COVERLAP(ISP_,LNXX_,DOVER_)
!     ******************************************************************
!     **  RETURN 1-C- OVERLAP OF THE PARTIAL WAVES                    **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNXX_
      REAL(8)   ,INTENT(OUT):: DOVER_(LNXX_,LNXX_)
      INTEGER(4)            :: LN1,LN2
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN1=1,THIS%LNX
        DO LN2=1,THIS%LNX
          DOVER_(LN1,LN2)=THIS%DOVER(LN1,LN2)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$1CKINETIC(ISP_,LNXX_,DTKIN_)
!     ******************************************************************
!     **  RETURN 1-C- KINETIC ENERGY OVERLAP OF THE PARTIAL WAVES     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNXX_
      REAL(8)   ,INTENT(OUT):: DTKIN_(LNXX_,LNXX_)
      INTEGER(4)            :: LN1,LN2
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN1=1,THIS%LNX
        DO LN2=1,THIS%LNX
          DTKIN_(LN1,LN2)=THIS%DTKIN(LN1,LN2)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$RCSM(ISP_,RCSM_)
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: RCSM_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      RCSM_=THIS%RCSM
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$RCBG(ISP_,RCBG_)
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: RCBG_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      RCBG_=THIS%RCBG
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LNX(ISP_,LNX_)
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LNX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LNX_=THIS%LNX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMNX(ISP_,LMNX_)
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LMNX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LMNX_=THIS%LMNX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMNXX(LMNXX_)
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT):: LMNXX_
!     ******************************************************************
      LMNXX_=LMNXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LOFLN(ISP_,LNX_,LOX_)
!     ******************************************************************
!     **  RETURN NUMBER MAIN ANGULAR MOMENTUM OF PARTIAL WAVES        **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNX_
      INTEGER(4),INTENT(OUT):: LOX_(LNX_)
      INTEGER(4)            :: LN
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      IF(LNX_.GT.THIS%LNX) THEN
        CALL ERROR$MSG('LNX ON INPUT TOO SMALL')
        CALL ERROR$OVERFLOW('LNX_',LNX_,THIS%LNX)
        CALL ERROR$STOP('SETUP$LOFLN')
      END IF
      DO LN=1,THIS%LNX
        LOX_(LN)=THIS%LOX(LN)
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$NSPECIES(NSP_)
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NSP_
!     ******************************************************************
      NSP_=NSP
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMRXX(LMRXX_)
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ****************************************************************** 
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: LMRXX_
!     ******************************************************************
      LMRXX_=LMRXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LNXX(LNXX_)
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: LNXX_
!     ******************************************************************
      LNXX_=LNXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMRX(ISP_,LMRX_)
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LMRX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LMRX_=THIS%LMRX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$AEZ(ISP_,AEZ_)
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: AEZ_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      AEZ_=THIS%AEZ
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$RADGRID(ISP_,R1_,DEX_,NR_)
!     ******************************************************************
!     **  RETURN THE  PARAMETERS DETERMINING THE RADIAL GRID          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: R1_
      REAL(8)   ,INTENT(OUT):: DEX_
      INTEGER(4),INTENT(OUT):: NR_
!     ******************************************************************
      R1_=R1
      DEX_=DEX
      NR_=NR
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP_COMPENSATIONGAUSSIANS(RCBG,RCSM &
     &                ,NG,G2,CELLVOL,G0,V0,DG0,DV0)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RCBG
      REAL(8)   ,INTENT(IN) :: RCSM
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(IN) :: G2(NG)
      REAL(8)   ,INTENT(IN) :: CELLVOL
      REAL(8)   ,INTENT(OUT):: G0(NG)
      REAL(8)   ,INTENT(OUT):: V0(NG)
      REAL(8)   ,INTENT(OUT):: DG0(NG)  ! |G|DG0/DG
      REAL(8)   ,INTENT(OUT):: DV0(NG)  ! |G|DV0/DG
      LOGICAL(4)            :: TGAMMA
      REAL(8)   ,PARAMETER  :: EPSILONGAMMA=1.D-7
      REAL(8)               :: PI
      REAL(8)               :: SVAR1,SVAR2,SVAR3,SVAR4
      REAL(8)               :: BGGAUSS,SMGAUSS
      INTEGER(4)            :: IG,NSTART
!     ******************************************************************
      CALL ERROR$MSG('MARKED FOR DELETION: DO NOT USE!')
      CALL ERROR$STOP('SETUP_COMPENSATIONGAUSSIANS')
      PI=4.D0*DATAN(1.D0)
!
!     ==================================================================
!     == CALCULATE GAUSSIANS ETC FOR COMPENSATION CHARGE DENSITY      ==
!     ==================================================================
      SVAR1=-0.25D0*RCBG**2
      SVAR2=-0.25D0*RCSM**2
      SVAR3=4.D0*PI/CELLVOL
      SVAR4=-4.D0*PI*SVAR3
      DO IG=1,NG
        IF(G2(IG).GT.1.D-12) THEN
          BGGAUSS=DEXP(SVAR1*G2(IG))
          SMGAUSS=DEXP(SVAR2*G2(IG))
          G0(IG)=SVAR3*BGGAUSS
          V0(IG)=SVAR4*(BGGAUSS-SMGAUSS)/G2(IG)
          DG0(IG)=2.D0*SVAR1*G2(IG)*G0(IG)
          DV0(IG)=2.D0*SVAR4*(SVAR1*BGGAUSS-SVAR2*SMGAUSS)-2.D0*V0(IG)
        ELSE
          G0(IG)=SVAR3
          V0(IG)=PI*(RCBG**2-RCSM**2)*SVAR3
          V0(IG)=SVAR4*(SVAR1-SVAR2)
          DG0(IG)=0.D0
          DV0(IG)=0.D0
        END IF
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP_COMPOFG(RCBG,RCSM,G1,DEX,NG,G0,V0)
!     ******************************************************************
!     **                                                              **
!     **  COMPENSATION DENSITY AND POTENTIAL ON A RADIAL GRID         **
!     **  IN G-SPACE                                                  **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RCBG
      REAL(8)   ,INTENT(IN) :: RCSM
      REAL(8)   ,INTENT(IN) :: G1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(OUT):: G0(NG)
      REAL(8)   ,INTENT(OUT):: V0(NG)
      LOGICAL(4)            :: TGAMMA
      REAL(8)   ,PARAMETER  :: EPSILONGAMMA=1.D-7
      REAL(8)               :: PI
      REAL(8)               :: SVAR1,SVAR2,SVAR3,SVAR4
      REAL(8)               :: BGGAUSS,SMGAUSS
      INTEGER(4)            :: IG,NSTART
      REAL(8)               :: GI2,xexp2
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      SVAR1=-0.25D0*RCBG**2
      SVAR2=-0.25D0*RCSM**2
      SVAR3=4.D0*PI
      SVAR4=-4.D0*PI*SVAR3
      XEXP2=EXP(DEX)**2
      GI2=G1**2/XEXP2
      DO IG=1,NG
        GI2=GI2*XEXP2
        BGGAUSS=DEXP(SVAR1*GI2)
        SMGAUSS=DEXP(SVAR2*GI2)
        G0(IG)=SVAR3*BGGAUSS
        V0(IG)=SVAR4*(BGGAUSS-SMGAUSS)/GI2
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP_CHECKGAUSS(CELLVOL,RC,TOL,GMAX,RMAX,TCHKR,TCHKG)
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: CELLVOL
      REAL(8)    ,INTENT(IN) :: RC
      REAL(8)    ,INTENT(IN) :: TOL
      REAL(8)    ,INTENT(IN) :: GMAX
      REAL(8)    ,INTENT(IN) :: RMAX
      LOGICAL(4) ,INTENT(OUT):: TCHKR
      LOGICAL(4) ,INTENT(OUT):: TCHKG
      REAL(8)                :: CHECK
!     ******************************************************************
      TCHKG=.TRUE.
      TCHKR=.TRUE.
      CHECK=-1.D0/CELLVOL*DEXP(-0.25D0*(RC*GMAX)**2)
      IF(DABS(CHECK).GT.TOL) TCHKG=.FALSE. 
      CHECK=-DEXP(-(RMAX/RC)**2)
      IF(DABS(CHECK).GT.TOL) TCHKR=.FALSE. 
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP_RADTOG(L,R1,DEX,NR,RMAX,FOFR,NG,CELLVOL,G2,FOFG,DFOFG)
!     **                                                              **
!     **  TRANSFORMS A FUNCTION FOFR GIVEN ON A LOGARITHMIC GRID      **
!     **  INTO G-SPACE                                                **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: R1
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(IN) :: RMAX
      REAL(8)   ,INTENT(IN) :: CELLVOL
      REAL(8)   ,INTENT(IN) :: FOFR(NR)
      REAL(8)   ,INTENT(IN) :: G2(NG)
      REAL(8)   ,INTENT(OUT):: FOFG(NG)
      REAL(8)   ,INTENT(OUT):: DFOFG(NG)   !DF/DG
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: FAC
      REAL(8)               :: G
      REAL(8)               :: RES
      INTEGER(4)            :: IG
!     == NEW VERSION=========
      REAL(8)   ,PARAMETER  :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
      INTEGER(4),PARAMETER  :: NR1=256
      REAL(8)               :: G1          ! FIRST POINT ON THE RADIAL G-GRID
      REAL(8)               :: ERROR       ! ERROR ESTIMATE
      REAL(8)               :: FOFR1(NR1)   
      REAL(8)               :: FOFG1(NR1)
!     ******************************************************************
      CALL ERROR$MSG('MARKED FOR DELETION: DO NOT USE!')
      CALL ERROR$STOP('SETUP_RADTOG')
      PI=4.D0*DATAN(1.D0)
      FOURPI=4.D0*PI
      FAC=FOURPI/CELLVOL
      G1=GMAX*DEXP(-DEX*DBLE(NR-1))
      IF(NR1.LT.NR) THEN
        CALL ERROR$STOP('SETUP_RADTOG')
      END IF
      FOFR1(1:NR)=FOFR(1:NR)
      FOFR1(NR+1:NR1)=0.D0
      CALL BESSELTRANSFORM(L,NR1,R1,G1,DEX,FOFR1,FOFG1,ERROR)
      FOFG1(:)=FAC*FOFG1(:)
!     == REMARK: THE EXTRAPOLATION TO THE GAMMA POINT MAY BE PROBLEMATIC 
!     == BECAUSE OF WRIGGLES IN THE FOURIER TRANSFORMED FUNCTION
!     == AN EARLIER VERSION USED THE VALUE OF THE FIRST RADIAL GRID POINT
      G=DSQRT(G2(1))
      CALL RADIAL$VALUE(G1,DEX,NR1,FOFG1,G,FOFG(1))
      CALL RADIAL$DERIVATIVE(G1,DEX,NR1,FOFG1,G,DFOFG(1))
      DFOFG(1)=G*DFOFG(1)
      DO IG=2,NG
        IF(DABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
          FOFG(IG) =FOFG(IG-1)
          DFOFG(IG)=DFOFG(IG-1)
        ELSE
          G=DSQRT(G2(IG))
          CALL RADIAL$VALUE(G1,DEX,NR1,FOFG1,G,FOFG(IG))
          CALL RADIAL$DERIVATIVE(G1,DEX,NR1,FOFG1,G,DFOFG(IG))
          DFOFG(IG)=G*DFOFG(IG)
        END IF
      ENDDO
      RETURN
      END      
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$LNX(NFIL,LNX)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(OUT) :: LNX
      REAL(8)                :: R1,DEX,NR
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$LNX')
      REWIND NFIL
      READ(NFIL,FMT='(F15.10,F10.5,2I4)')R1,DEX,NR,LNX
!     READ(NFIL,FMT='(F15.10,F10.5,2I4,2F5.2,F20.15,I5)')R1,DEX,NR,LNX
                              CALL TRACE$POP
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$READALL(NFIL,NRX,R1,DEX,NR,LNX,LOX,AEZ,PSZ &
     &         ,PSPHI,AEPHI,VADD,RCSM &
     &         ,DTKIN,DOVER,IRCCOR,RHOCOR,PSCORR,PRO)
!     ******************************************************************
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)  :: NFIL
      INTEGER(4) ,INTENT(IN)  :: NRX
      REAL(8)    ,INTENT(IN)  :: R1
      REAL(8)    ,INTENT(IN)  :: DEX
      INTEGER(4) ,INTENT(IN)  :: NR
      INTEGER(4) ,INTENT(IN)  :: LNX
      REAL(8)    ,INTENT(OUT) :: AEZ
      REAL(8)    ,INTENT(OUT) :: PSZ
      REAL(8)    ,INTENT(OUT) :: RCSM
      REAL(8)    ,INTENT(OUT) :: VADD(NRX)
      REAL(8)    ,INTENT(OUT) :: PRO(NRX,LNX)
      INTEGER(4) ,INTENT(OUT) :: LOX(LNX)
      INTEGER(4) ,INTENT(OUT) :: IRCCOR
      REAL(8)    ,INTENT(OUT) :: DTKIN(LNX,LNX)
      REAL(8)    ,INTENT(OUT) :: DOVER(LNX,LNX)
      REAL(8)    ,INTENT(OUT) :: AEPHI(NRX,LNX)
      REAL(8)    ,INTENT(OUT) :: PSPHI(NRX,LNX)
      REAL(8)    ,INTENT(OUT) :: RHOCOR(NRX)
      REAL(8)    ,INTENT(OUT) :: PSCORR(NRX)
      REAL(8)                 :: R11,DEX1
      INTEGER(4)              :: NR1,LNX1,I,IR,LN1,LN2,LN
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$READALL')
      REWIND NFIL
      READ(NFIL,6000)R11,DEX1,NR1,LNX1,PSZ,AEZ,RCSM,IRCCOR
6000  FORMAT(F15.10,F10.5,2I4,2F5.2,F20.15,I5)
      IF(R11.NE.R1.OR.DEX1.NE.DEX.OR.NR1.NE.NR) THEN
        CALL ERROR$MSG('ONLY ONE TYPE OF RADIAL GRID ALLOWED')
        CALL ERROR$STOP('INPOT')
      END IF
      IF(IRCCOR.LE.0.OR.IRCCOR.GT.NR) THEN
!      PRINT*,'WARNING! NO MT-RADIUS SPECIFIED FOR ATOM WITH Z=',AEZ
       IRCCOR=NR-2
      END IF
      IF(LNX.NE.LNX1) THEN
        CALL ERROR$MSG('LNX OUT OF RANGE')
        CALL ERROR$I4VAL('LNX ON INPUT',LNX)
        CALL ERROR$I4VAL('LNX ON FILE',LNX1)
        CALL ERROR$STOP('INPOT')
      END IF
                              CALL TRACE$PASS('BEFORE LOX')
      READ(NFIL,6020)(LOX(I),I=1,LNX)
6020  FORMAT(14I5)
                              CALL TRACE$PASS('BEFORE VADD')
      READ(NFIL,6100)(VADD(IR),IR=1,NR)
!     ====  RHOCOR = CORE CHARGE DENSITY  ==============================
                              CALL TRACE$PASS('BEFORE RHOCOR')
      READ(NFIL,6100)(RHOCOR(IR),IR=1,NR)
!     ====  PSCORR = PSEUDO CORE CHARGE DENSITY ========================
                              CALL TRACE$PASS('BEFORE PSCORR')
      READ(NFIL,6100)(PSCORR(IR),IR=1,NR)
!     ====  DTKIN = <AEPHI|-DELTA/2|AEPHI> - <PSPHI|-DELTA/2|PSPHI> ====
                              CALL TRACE$PASS('BEFORE DTKIN')
      READ(NFIL,6100)((DTKIN(LN1,LN2),LN1=1,LNX),LN2=1,LNX)
!     PRINT*,'DTKIN ',(DTKIN(LN,LN),LN=1,LNX)
!     ====  DOVER = <AEPHI|AEPHI> - <PSPHI|PSPHI> ======================
                              CALL TRACE$PASS('BEFORE DOVER')
      READ(NFIL,6100)((DOVER(LN1,LN2),LN1=1,LNX),LN2=1,LNX)
                              CALL TRACE$PASS('BEFORE PROJECTORS')
      DO 100 LN=1,LNX
      READ(NFIL,6100)(PRO(IR,LN),IR=1,NR)
      READ(NFIL,6100)(AEPHI(IR,LN),IR=1,NR)
      READ(NFIL,6100)(PSPHI(IR,LN),IR=1,NR)
6100  FORMAT(SP,5E14.8)
100   CONTINUE
                              CALL TRACE$POP
      RETURN
      END


