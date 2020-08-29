!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  OBJECTNAME: AUGMENTATION                                         **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    MOMENTS                                                        **
!**    SPHERE                                                         **
!**    SYNC                                                           **
!**                                                                   **
!**  MODULES USED: MPE                                                **
!**  USES: MPE; SETUP; HYPERFINE, ENERGYLIST, ATOMLIST                **
!**    EXTERNAL1CPOT
!**    CLEBSCH, RADIAL, DFT, GAUSSN                                   **
!**    SELFTEST, TRACE, FILEHANDLER, ERROR                            **
!**                                                                   **
!***********************************************************************
!***********************************************************************
MODULE AUGMENTATION_MODULE
INTEGER(4)  ,PARAMETER :: NE=10
CHARACTER(32)          :: ID(NE)
REAL(8)                :: VAL(NE)
LOGICAL(4)             :: TINI=.FALSE.
END MODULE AUGMENTATION_MODULE
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_INI
!     ******************************************************************
!     **  RESETS AUGMENTATION_MODULE                                  **
!     **  (CALLED BY AUGMENTATION_ADD)                                **
!     **                                                              **
!     ******************************************************************
      USE AUGMENTATION_MODULE
      IMPLICIT NONE
      INTEGER(4) :: I
!     ******************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!
      ID(1)='AE1 EXCHANGE-CORRELATION'
      ID(2)='PS1 EXCHANGE-CORRELATION'
      ID(3)='AE1 ELECTROSTATIC'       
      ID(4)='PS1 ELECTROSTATIC'       
      ID(5)='AE1 BACKGROUND'          
      ID(6)='PS1 BACKGROUND'          
      ID(7)='AE1-PS1 KINETIC'         
      ID(8)='LDA+U EXCHANGE'          
      ID(9)='CORE RELAXATION'          
      ID(10)='EXTERNAL 1CENTER POTENTIAL'          
      DO I=1,NE
        VAL(NE)=0.D0
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADD(ID_,VAL_)
!     ******************************************************************
!     **  ADDS DATA INTO THE TEMPORARY STACK FOR SYNCHRONIZATION      **
!     **  WITH THE ENERGYLIST                                         **
!     ******************************************************************
      USE AUGMENTATION_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_
      INTEGER(4)              :: I
!     ******************************************************************
      IF(.NOT.TINI) CALL AUGMENTATION_INI
      DO I=1,NE
        IF(ID_.EQ.ID(I)) THEN
          VAL(I)=VAL(I)+VAL_
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ID NOT RECOGNIZED')
      CALL ERROR$CHVAL('ID ON INPUT   ',ID_)
      DO I=1,NE
        CALL ERROR$CHVAL('POSSIBLE ID ',ID(I))
      ENDDO
      CALL ERROR$STOP('AUGMENTATION_ADD')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION$SYNC
!     ******************************************************************
!     **  SYNCHRONIZES THE TEMPORARY STACK WITH THE ENERGYLIST        **
!     **  AND RESETS THE TEMPORARY STACK                              **
!     **                                                              **
!     **  MUST BE CALLED AFTER AUGMENTATION$SPHERE HAS BEEN           **
!     **  PROCESSED FOR ALL ATOMS                                     **
!     ******************************************************************
      USE AUGMENTATION_MODULE, ONLY : TINI &
     &                               ,NE &
     &                               ,ID &
     &                               ,VAL
      USE MPE_MODULE, ONLY : MPE$COMBINE
      IMPLICIT NONE
      INTEGER(4) :: I
      REAL(8)    :: SVAR
!     ******************************************************************
      TINI=.FALSE.
      CALL MPE$COMBINE('MONOMER','+',VAL)
      DO I=1,NE
        CALL ENERGYLIST$SET(ID(I),VAL(I))
      ENDDO
!     ==  ID(1)='AE1 EXCHANGE-CORRELATION' =============================
!     ==  ID(2)='PS1 EXCHANGE-CORRELATION' =============================
!     ==  ID(3)='AE1 ELECTROSTATIC' ====================================
!     ==  ID(4)='PS1 ELECTROSTATIC' ====================================
!     ==  ID(5)='AE1 BACKGROUND' =======================================
!     ==  ID(6)='PS1 BACKGROUND' =======================================
!     ==  ID(7)='AE1-PS1 KINETIC' ======================================
!     ==  ID(8)='LDA+U EXCHANGE' =======================================
!     ==  ID(9)='CORE RELAXATION' ======================================
!     ==  ID(10)='EXTERNAL 1CENTER POTENTIAL' ===========================
      SVAR=VAL(7)+VAL(1)-VAL(2)+VAL(3)-VAL(4)+VAL(5)-VAL(6)+VAL(8)+VAL(9) &
     &    +VAL(10)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',SVAR)
      CALL ENERGYLIST$ADD('AE  KINETIC',VAL(7))
      CALL ENERGYLIST$ADD('AE  EXCHANGE-CORRELATION',VAL(1)-VAL(2))
      CALL ENERGYLIST$ADD('    LOCAL CORRELATIONS',VAL(8))
      CALL ENERGYLIST$ADD('AE  ELECTROSTATIC',VAL(3)-VAL(4))
      CALL ENERGYLIST$ADD('BACKGROUND',VAL(5)-VAL(6))
      VAL(:)=0.D0
      RETURN
      END
!
!     .....................................................MOMNTS.......
      SUBROUTINE AUGMENTATION$MOMENTS(ISP,LMNX,DENMAT,LMRX,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ELECTROSTATIC MULTIPOLE MOMENTS FOR THE      **
!     **  COMPENSTATION CHARGE DENSITY                                **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      INTEGER(4),INTENT(IN)  :: ISP
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: LMRX
      COMPLEX(8),INTENT(IN)  :: DENMAT(LMNX,LMNX)
      REAL(8)   ,INTENT(OUT) :: QLM(LMRX)
      REAL(8)   ,ALLOCATABLE :: AERHO(:,:)   !(NR,LMRX)
      REAL(8)   ,ALLOCATABLE :: PSRHO(:,:)   !(NR,LMRX)
      REAL(8)   ,ALLOCATABLE :: AECORE(:)    !(NR)
      REAL(8)   ,ALLOCATABLE :: PSCORE(:)    !(NR)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)   !(NR,LNX)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)   !(NR,LNX)
      INTEGER(4),ALLOCATABLE :: LOX(:)       !(LN)
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)                :: AEZ
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMN1,LMN2
!     ******************************************************************
                          CALL TRACE$PUSH('AUGMENTATION$MOMENTS')
!
!     == PRINT DENSITY MATRIX FOR TEST ========== ======================
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(*,FMT='("TOTAL DENSITY MATRIX")')
        DO LMN1=1,LMNX
          WRITE(*,FMT='(9F10.6)')(DENMAT(LMN1,LMN2),LMN2=1,LMNX)
        ENDDO
      END IF
!     
!     ==================================================================
!     ==   RECEIVE PARTIAL WAVES FROM SETUP OBJECT                    ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
!
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI)
!
!     ==================================================================
!     ==  RECEIVE CORE DENSITY FROM SETUPS OBJECT                     ==
!     ==  RECEIVE ATOMIC NUMBER FROM ATOMLIST                         ==
!     ==================================================================
      ALLOCATE(PSCORE(NR))
      ALLOCATE(AECORE(NR))
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
      CALL SETUP$GETR8A('PSCORE',NR,PSCORE)
      CALL SETUP$GETR8('AEZ',AEZ)
!     
!     ==================================================================
!     ==   CALCULATE ONE CENTER CHARGE DENSITY                        ==
!     ==================================================================
      ALLOCATE(AERHO(NR,LMRX))
      ALLOCATE(PSRHO(NR,LMRX))
      CALL AUGMENTATION_RHO(NR,LNX,LOX,AEPHI,LMNX,DENMAT,LMRX,AERHO)
      CALL AUGMENTATION_RHO(NR,LNX,LOX,PSPHI,LMNX,DENMAT,LMRX,PSRHO)
      DEALLOCATE(AEPHI)
      DEALLOCATE(PSPHI)
      DEALLOCATE(LOX)
!     
!     ==================================================================
!     ==   CALCULATE ONE MULTIPOLE MOMENTS                            ==
!     ==================================================================
      CALL AUGMENTATION_QLM(GID,NR,LMRX,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
      DEALLOCATE(AERHO)
      DEALLOCATE(PSRHO)
      DEALLOCATE(PSCORE)
      DEALLOCATE(AECORE)
      IF(TPR) THEN
        WRITE(*,FMT='(A63)')'MULTIPOLE MOMENTS'
        WRITE(*,FMT='(9E12.5)')QLM(:)
      END IF
      CALL SETUP$UNSELECT()
                   CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION$SPHERE(ISP,IAT,LMNX,NDIMD,DENMAT,EDENMAT &
     &                              ,LMRX,VQLM,RHOB,POTB,DATH,DO)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ** THREE OPTIONS:                                               **
!     **   1 - NON-SPIN POLARIZED (NT)                                **
!     **   2 - SPIN-POLARIZED (NT,NS)                                 **
!     **   4 - NONCOLLINEAR SPIN  (NT,NSX,NSY,NSZ)                    **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: ISP
      INTEGER(4),INTENT(IN)   :: IAT
      INTEGER(4),INTENT(IN)   :: LMNX
      INTEGER(4),INTENT(IN)   :: NDIMD
      COMPLEX(8),INTENT(INOUT):: DENMAT(LMNX,LMNX,NDIMD)
      COMPLEX(8),INTENT(INOUT):: EDENMAT(LMNX,LMNX,NDIMD)
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(INOUT):: VQLM(LMRX)
      REAL(8)   ,INTENT(IN)   :: RHOB ! NEUTRALIZING BACKGROUND DENSITY
!                               RHOB MAY BE SET TO ZERO BY ISOLATE OBJECT
      REAL(8)   ,INTENT(OUT)  :: POTB ! NEG. AV. EL. AUGM. POT.
      COMPLEX(8),INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT)  :: DO(LMNX,LMNX,NDIMD)
      REAL(8)    ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)              :: GID   ! GRID ID
      REAL(8)   ,ALLOCATABLE  :: R(:)  ! RADIAL GRID
      REAL(8)                 :: AEZ
      INTEGER(4)              :: NR
      INTEGER(4)              :: LNX
      INTEGER(4),ALLOCATABLE  :: LOX(:)
      REAL(8)   ,ALLOCATABLE  :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: AECORE(:)
      REAL(8)   ,ALLOCATABLE  :: PSCORE(:)
      REAL(8)                 :: RCSM
      REAL(8)   ,ALLOCATABLE  :: VADD(:)
      REAL(8)   ,ALLOCATABLE  :: DOVER(:,:)
      REAL(8)   ,ALLOCATABLE  :: DTKIN(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DTKIN1(:,:)
      REAL(8)                 :: DELTAH(LMNX,LMNX,NDIMD)  ! SOFTCORE CORRECTION TO DTKIN
      REAL(8)                 :: DELTAO(LMNX,LMNX,NDIMD)  ! SOFTCORE CORRECTION TO DO
      REAL(8)   ,ALLOCATABLE  :: DELTARHO(:,:,:)    ! SOFCORE CORRECTION TO DAERHO
      REAL(8)                 :: DECORE
      REAL(8)   ,ALLOCATABLE  :: AERHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DATP(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DWORK1(:)
      INTEGER(4)              :: IDIM,IR
      INTEGER(4)              :: LM,LMN1
      INTEGER(4)              :: NFILO
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=1
      LOGICAL(4)              :: TBACK
      REAL(8)                 :: DETOT,PSEHARTREE,AEEHARTREE,COREEXC
      REAL(8)                 :: EKINNL,ENL,AEEXC,PSEXC
      CHARACTER(32)           :: ATOM
      REAL(8)                 :: VQLM1(LMRX)
      REAL(8)                 :: QLM(LMRX)
      REAL(8)   ,ALLOCATABLE  :: AEHPOT(:,:)
      REAL(8)   ,ALLOCATABLE  :: AEXCPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSHPOT(:,:)
      REAL(8)   ,ALLOCATABLE  :: PSXCPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: AETOTPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSTOTPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: RHO(:,:,:)
      CHARACTER(32)           :: SOFTCORETYPE
!     ******************************************************************
                            CALL TRACE$PUSH('AUGMENTATION$SPHERE')
!
!     ==================================================================
!     ==  COLLECT ATOM-TYPE SPECIFIC INFORMATION FROM SETUP OBJECT    ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      CALL SETUP$GETR8('AEZ',AEZ)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI)
      ALLOCATE(AECORE(NR))
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
      ALLOCATE(PSCORE(NR))
      CALL SETUP$GETR8A('PSCORE',NR,PSCORE)
      CALL SETUP$GETR8('RCSM',RCSM)
      ALLOCATE(VADD(NR))
      CALL SETUP$GETR8A('VADD',NR,VADD)
!
!     ==================================================================
!     ==  EXPAND OVERLAP AND KINETIC ENERGY MATRIX ELEMENTS           ==
!     ==================================================================
      ALLOCATE(DOVER(LNX,LNX))
      CALL SETUP$GETR8A('DO',LNX*LNX,DOVER)
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DOVER,DO)
      DEALLOCATE(DOVER)
!
      ALLOCATE(DTKIN(LMNX,LMNX,NDIMD))
      ALLOCATE(DTKIN1(LNX,LNX))
      CALL SETUP$GETR8A('DEKIN',LNX*LNX,DTKIN1)
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DTKIN1,DTKIN)
      DEALLOCATE(DTKIN1)
      CALL SETUP$UNSELECT()
!
!     ==================================================================
!     ==  SELF TEST                                                   ==
!     ==================================================================
 1000 CONTINUE
      IF(TTEST) THEN
        IF(ITEST.EQ.1) THEN 
          CALL SELFTEST$START('SPHERE',LMNX*LMNX*NDIMD,DENMAT,1.D-3)
          VQLM(:)=0.D0
          DO IDIM=1,NDIMD
            DENMAT(:,:,IDIM)=0.5D0*(DENMAT(:,:,IDIM)+TRANSPOSE(CONJG(DENMAT(:,:,IDIM))))
          ENDDO
        END IF
      ENDIF
!     == PRINT DENSITY MATRIX FOR TEST
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,FMT='("DENMAT FOR IDIM= ",I2)') IDIM
          DO LMN1=1,LMNX
            WRITE(NFILO,FMT='(10F15.6)')DENMAT(LMN1,:,IDIM)
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     ==  CALCULATE 1-CENTER KINETIC ENERGY                           ==
!     ==================================================================
!     == DENSITY MATRIX IS HERMITEAN FOR EACH SPIN DIRECTION
      EKINNL=REAL(SUM(CONJG(DENMAT)*DTKIN),KIND=8)
      CALL AUGMENTATION_ADD('AE1-PS1 KINETIC',EKINNL)
!
!     ==================================================================
!     ==  CALCULATE 1-CENTER CHARGE DENSITY                           ==
!     ==================================================================
      ALLOCATE(AERHO(NR,LMRX,NDIMD))
      ALLOCATE(PSRHO(NR,LMRX,NDIMD))
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,AEPHI &
     &                  ,LMNX,DENMAT(1,1,IDIM),LMRX,AERHO(1,1,IDIM))
        CALL AUGMENTATION_RHO(NR,LNX,LOX,PSPHI &
     &                  ,LMNX,DENMAT(1,1,IDIM),LMRX,PSRHO(1,1,IDIM))
      ENDDO
!     
!     ================================================================
!     ==  EVALUATE MULTIPOLE MOMENTS                                ==
!     ================================================================
      CALL AUGMENTATION_QLM(GID,NR,LMRX,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
      ALLOCATE(RHO(NR,LMRX,1))
      RHO(:,:,1)=PSRHO(:,:,1)
      RHO(:,1,1)=RHO(:,1,1)+PSCORE(:)   
      CALL AUGMENTATION_ADJUSTVQLM(GID,NR,LMRX,RHO,RCSM,QLM,VQLM1)
      VQLM1(:)=VQLM1(:)+VQLM(:)
      DEALLOCATE(RHO)
!     
!     ================================================================
!     ==  NEW SOFT CORE                                             ==
!     ================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETCH('SOFTCORETYPE',SOFTCORETYPE)
      CALL SETUP$UNSELECT()
      IF(SOFTCORETYPE.EQ.'NONE') THEN
        DECORE=0.D0
      ELSE 
        ALLOCATE(DELTARHO(NR,LMRX,NDIMD))
        DELTAH=0.D0
        DELTAO=0.D0
!        CALL AUGMENTATION_NEWSOFTCORE(SOFTCORETYPE,GID,NR,LMRX,NDIMD,AEZ,LMNX &
!     &          ,DENMAT,EDENMAT,VQLM1,RHOB,DELTAH,DELTAO,DELTARHO,DECORE)
        DTKIN=DTKIN+DELTAH
        DO=DO+DELTAO
        DEALLOCATE(DELTARHO)
      END IF
!
!     =================================================================
!     == HARTREE ENERGY AND POTENTIAL                                ==
!     =================================================================
      ALLOCATE(PSHPOT(NR,LMRX))
      ALLOCATE(AEHPOT(NR,LMRX))
      CALL AUGMENTATION_PSHARTREE(GID,NR,LMRX,PSCORE,PSRHO(:,:,1) &
     &                 ,VADD,RCSM,QLM,VQLM1,RHOB,PSHPOT,PSEHARTREE)
      CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AEZ,AECORE,AERHO(:,:,1) &
     &                               ,VQLM1,RHOB,AEHPOT,AEEHARTREE)
!
!     ================================================================
!     ==   ADD EXCHANGE AND CORRELATION POTENTIAL                   ==
!     ================================================================
      ALLOCATE(AEXCPOT(NR,LMRX,NDIMD))
      ALLOCATE(PSXCPOT(NR,LMRX,NDIMD))
      ALLOCATE(RHO(NR,LMRX,NDIMD))
!     == CORE ONLY EXCHANGE ENERGY ===================================
      CALL AUGMENTATION_XC(GID,NR,1,1,AECORE,COREEXC,AEXCPOT(:,1,1))
      AEXCPOT(:,:,:)=0.D0
!     == AE-EXCHANGE ENERGY AND POTENTIAL ============================
      RHO(:,:,:)=AERHO(:,:,:)
      RHO(:,1,1)=RHO(:,1,1)+AECORE(:)
      CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,AEEXC,AEXCPOT)
      AEEXC=AEEXC-COREEXC
!     == PS-EXCHANGE ENERGY AND POTENTIAL ============================
      RHO(:,:,:)=PSRHO(:,:,:)
      RHO(:,1,1)=RHO(:,1,1)+PSCORE(:)
      CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,PSEXC,PSXCPOT)
      DEALLOCATE(RHO)
!
!     ================================================================
!     ==   ADD UP POTENTIALS                                        ==
!     ================================================================
      ALLOCATE(AETOTPOT(NR,LMRX,NDIMD))
      ALLOCATE(PSTOTPOT(NR,LMRX,NDIMD))
      AETOTPOT=AEXCPOT
      PSTOTPOT=PSXCPOT
      AETOTPOT(:,:,1)=AETOTPOT(:,:,1)+AEHPOT(:,:)
      PSTOTPOT(:,:,1)=PSTOTPOT(:,:,1)+PSHPOT(:,:)
!
!     ================================================================
!     ==   REPORT ENERGIES                                          ==
!     ================================================================
      CALL AUGMENTATION_ADD('AE1 BACKGROUND',0.D0)
      CALL AUGMENTATION_ADD('PS1 BACKGROUND',0.D0)
      CALL AUGMENTATION_ADD('AE1 EXCHANGE-CORRELATION',AEEXC)
      CALL AUGMENTATION_ADD('PS1 EXCHANGE-CORRELATION',PSEXC)
      CALL AUGMENTATION_ADD('AE1 ELECTROSTATIC',AEEHARTREE)
      CALL AUGMENTATION_ADD('PS1 ELECTROSTATIC',PSEHARTREE)
      CALL AUGMENTATION_ADD('CORE RELAXATION',DECORE)
!
!     =================================================================
!     == AVERAGE ELECTROSTATIC ONE-CENTER POTENTIAL                  ==
!     =================================================================
      ALLOCATE(DWORK1(NR))
      DWORK1(:)=AEHPOT(:,1)-(PSHPOT(:,1)-VADD(:))
      DWORK1(:)=-DWORK1(:)*R(:)**2*SQRT(4.D0*PI)  !SQRT(4*PI)=4*PI*Y_0
      CALL RADIAL$INTEGRAL(GID,NR,DWORK1,POTB)
      DEALLOCATE(DWORK1)
!
!     ==========================================================================
!     ==  SEND DENSITIES TO HYPERFINE-PARAMETER OBJECT                        ==
!     ==========================================================================
      CALL AUGMENTATION_FEEDHYPERFINE(GID,NR,LMRX,NDIMD &
     &         ,IAT,AERHO,AECORE,PSRHO,PSCORE,AEHPOT,PSHPOT)
!     
!     ==========================================================================
!     ==  SEND DENSITIES AND POTENTIALS TO GRAPHICS OBJECT                    ==
!     ==========================================================================
      CALL GRAPHICS$SET1CPOT('HARTREE','AE',IAT,GID,NR,NR,LMRX,AEHPOT)
      CALL GRAPHICS$SET1CPOT('HARTREE','PS',IAT,GID,NR,NR,LMRX,PSHPOT)
      CALL GRAPHICS$SET1CPOT('TOT','AE',IAT,GID,NR,NR,LMRX,AETOTPOT(:,:,1))
      CALL GRAPHICS$SET1CPOT('TOT','PS',IAT,GID,NR,NR,LMRX,PSTOTPOT(:,:,1))
!     
!     ==========================================================================
!     ==  EVALUATE CORE SHIFTS                                                ==
!     ==========================================================================
      CALL CORE_CORESHIFTS(IAT,ISP,GID,NR,LMRX,AETOTPOT)
!
      DEALLOCATE(AEHPOT)
      DEALLOCATE(PSHPOT)
      DEALLOCATE(AEXCPOT)
      DEALLOCATE(PSXCPOT)
!
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,*)'AE POTENTIAL FOR SPIN ',IDIM
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(AETOTPOT(IR,LM,IDIM),LM=1,LMRX)
          ENDDO
        ENDDO
        DO IDIM=1,NDIMD
          WRITE(NFILO,*)'PS POTENTIAL FOR ATOM ',IAT,IDIM
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(PSTOTPOT(IR,LM,IDIM),LM=1,LMRX)
          ENDDO
        ENDDO
      ENDIF
!     
!     ==================================================================
!     ==  SELF TEST FOR HAMILTONIAN AND OVERLAP                       ==
!     ==================================================================
 1001 CONTINUE
      IF(TTEST) THEN
        IF(ITEST.EQ.2) THEN
          CALL SELFTEST$START('SPHERE-EKI',LMNX*LMNX*NDIMD,DENMAT,1.D-3)
        END IF
      END IF
!     
!     ==================================================================
!     ==  EVALUATE KINETIC HAMILTONIAN AND OVERLAP                    ==
!     ==================================================================
      DATH(:,:,:)=CMPLX(DTKIN(:,:,:),0.D0,KIND=8)
!
      IF(TTEST) THEN
        IF(ITEST.EQ.2) THEN
          CALL SELFTEST$END('SPHERE-EKI',LMNX*LMNX*NDIMD,DATH,EKINNL,TBACK)
          IF(TBACK) GOTO 1001
          CALL ERROR$MSG('STOP AFTER SELFTEST')
          CALL ERROR$STOP('SPHERE')
        END IF
      END IF
!     
!     
!     ================================================================
!     ==   ADD POTENTIAL ENERGY TO THE ONE-CENTER HAMILTONIAN       ==
!     ==   DATH =DATH + <AEPHI|AEPOT|AEPHI>-<PSPHI|PSPOT|PSPHI>     ==
!     ================================================================
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL AUGMENTATION_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                        ,AETOTPOT,PSTOTPOT,AEPHI,PSPHI,DATP)
      DATH(:,:,:)=DATH(:,:,:)+CMPLX(DATP(:,:,:),0.D0,KIND=8)
      DEALLOCATE(DATP)
!     WRITE(TESTSTRING,FMT='("R8DATH",I2,12(" "))')IAT
!     CALL STOREIT(TESTSTRING,8*LMNX*LMNX*NSPIN,DATH)
!
!
!     ================================================================
!     ==  CI                                                        ==
!     ================================================================
!      CALL PAW_CI(ISP,LMNX,NDIMD,DENMAT)
!
!     ================================================================
!     ==  APPLY EXTERNAL POTENTIAL                                  ==
!     ================================================================
      CALL ATOMLIST$GETCH('NAME',IAT,ATOM)
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATP,DETOT)
      DATH(:,:,:)=DATH(:,:,:)+CMPLX(DATP(:,:,:),0.D0,KIND=8)
      DEALLOCATE(DATP)
      CALL AUGMENTATION_ADD('EXTERNAL 1CENTER POTENTIAL',DETOT)
!     
!     ================================================================
!     ==  SELF-TEST                                                 ==
!     ================================================================
      IF(TTEST) THEN
        IF(ITEST.EQ.1) THEN
          ENL=AEEXC-PSEXC+AEEHARTREE-PSEHARTREE+EKINNL+DETOT
          CALL SELFTEST$END('SPHERE',LMNX*LMNX*NDIMD,DATH,ENL,TBACK)
          IF(TBACK) GOTO 1000
          CALL ERROR$MSG('STOP AFTER SELFTEST')
          CALL ERROR$STOP('SPHERE')
        ENDIF
      END IF
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,FMT='("DATH FOR IDIM= ",I2)') IDIM
          DO LMN1=1,LMNX
            WRITE(NFILO,FMT='(10F15.6)')REAL(DATH(LMN1,:,IDIM))
          ENDDO
        ENDDO
STOP
      END IF
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(AETOTPOT)
      DEALLOCATE(PSTOTPOT)
      DEALLOCATE(AERHO)
      DEALLOCATE(PSRHO)
      DEALLOCATE(LOX)
      DEALLOCATE(AEPHI)
      DEALLOCATE(PSPHI)
      DEALLOCATE(AECORE)
      DEALLOCATE(PSCORE)
      DEALLOCATE(VADD)
      DEALLOCATE(DTKIN)
      DEALLOCATE(R)
                               CALL TRACE$POP()
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_RHO(NR,LNX,LOX,PHI,LMNXX,DENMAT,LMRX,RHOL)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: LMRX
      COMPLEX(8),INTENT(IN) :: DENMAT(LMNXX,LMNXX)
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: RHOL(NR,LMRX)
      INTEGER(4)             :: LMN1,LN1,L1,IM1,LMN2,LN2,L2,IM2
      INTEGER(4)             :: LM1,LM2,LM3
      REAL(8)                :: CG     !CLEBSCH GORDAN COEFFICIENT
      REAL(8)                :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     ==   CALCULATE ONE CENTER CHARGE DENSITY                        ==
!     ==================================================================
!
!     ==================================================================
!     ==  ADD VALENCE CHARGE DENSITY                                  ==
!     ==================================================================
      RHOL(:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LM1=L1**2+IM1
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
!                 == ASSUMES REAL PARTIAL WAVES ==================
                  SVAR=CG*REAL(DENMAT(LMN1,LMN2))
                  RHOL(:,LM3)=RHOL(:,LM3)+SVAR*PHI(:,LN1)*PHI(:,LN2)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_QLM(GID,NR,LMRX &
     &                           ,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE MULTIPOLE MOMENTS OF THE DEVIATION OF             **
!     **  ALL-ELECTRON AND THE PSEUDO DENSITY FROM ONE CENTER         **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)   ,INTENT(IN) :: PSRHO(NR,LMRX)
      REAL(8)   ,INTENT(OUT):: QLM(LMRX)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: DWORK(NR)
      REAL(8)               :: RES
      INTEGER(4)            :: LM,L
      REAL(8)               :: R(NR)
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      QLM(:)=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        DWORK(:)=(AERHO(:,LM)-PSRHO(:,LM))*R(:)**(L+2)
        CALL RADIAL$INTEGRAL(GID,NR,DWORK,QLM(LM))
      ENDDO
!
      DWORK(:)=(AECORE(:)-PSCORE(:))*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,DWORK,RES)
      QLM(1)=QLM(1)+RES-AEZ*Y0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE                          **
!     **  SPHERICAL CONTRIBUTION OF THE DENSITY UP TO QUADRATIC               **
!     **  ORDER IN THE NON-SPHERICAL CONTRIBUTIONS                            **
!     **                                                                      **
!     **  EXC = EXC(XVAL(L=0)*Y0)                                             **
!     **      + 0.5 * D2[EXC]/D[XVAL(L=0)*Y0]**2 * XVAL(L>0)**2               **
!     **                                                                      **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)              **
!     **  IS AN SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.             **
!     **                                                                      **
!     **  DEPENDECIES:                                                        **
!     **    DFT                                                               **
!     **    TIMING                                                            **
!     **    TRACE                                                             **
!     **                                                                      **
!     **  REMARKS: THE GRADIENTS ARE CORRECT ONLY IF DFT SUPPORTS             **
!     **    THIRD DERIVATIVES OF THE XC-ENERGY                                **
!     **   - WHEN USING SELFTEST ON THIS ROUTINE, THEN                        **
!     **     D(EXC)/DRHO(I)=POT(I)*DEX*R(I)**3                                **
!     **     AND THE VALUES AT LARGE RADII MUST BE SURPRESSED                 **
!     **                                                                      **
!     **  REMARK: FOR A COLLINEAR DENSITY THE ROUTINE GIVES DIFFERENT RESULTS **
!     **          WITH NDIMD=2 AND NDIMD=4 DUE TO THE TAYLOR EXPANSION IN     **
!     **          ANGULAR MOMENTUM EXPANSIONS                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1996 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD      ! CAN BE 1,2,4
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      REAL(8)               :: FOURPI
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM
      INTEGER(4)            :: IMAX
      REAL(8)               :: FAC
      REAL(8)               :: CG0LL
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
      REAL(8)   ,PARAMETER  :: XX=1.D0  !BUGFIX 160509 XX=0.5->1.0
!     **************************************************************************
      CALL TRACE$PUSH('AUGMENTATION_XC')
      EXC=0.D0
      VXC(:,:,:)=0.D0
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
      FOURPI=4.D0*PI
      CG0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  OBTAIN SPIN DENSITY                                                 ==
!     ==========================================================================
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      ALLOCATE(RHO(NR,LMRX,NSPIN))
      ALLOCATE(GRHO(NR,LMRX,NSPIN))
      ALLOCATE(VRHO(NR,LMRX,NSPIN))
      RHO(:,:,1)=RHOIN(:,:,1)
      IF(NDIMD.EQ.2) THEN
        RHO(:,:,2)=RHOIN(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == HERE WE NEED TO CALCULATE THE ABSOLUTE VALUE OF THE SPIN DENSITY ====
!       == IN AN ANGULAR MOMENTUM EXPANSION. THIS IS ONLY POSSIBLE APPROXIMATELY
!       == USING A TAYLOR EXPANSION ABOUT THE SPHERICAL PART OF THE SQUARE =====
!       == OF THE SPIN DENSITY. ================================================
        VRHO(:,:,:)=0.D0
        CALL AUGMENTATION_NCOLLTRANS(GID,'RHO',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF
!
!     == IMAX ALLOWS TO RESTRICT SOME LOOPS (1:5) TO (1:IMAX)
      IF(TGRA) THEN
        IF(NSPIN.EQ.2) THEN; IMAX=5; ELSE; IMAX=3; END IF
      ELSE 
        IF(NSPIN.EQ.2) THEN; IMAX=2; ELSE; IMAX=1; END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE RADIAL GRADIENT OF THE DENSITY                            ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE GRADIENTS')
      IF(TGRA) THEN
        GRHO(:,:,:)=0.D0
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
          ENDDO
        ENDDO
      ELSE
        GRHO(:,:,:)=0.D0
      END IF
!
!     ==========================================================================
!     ==  DEFINE VECTOR (RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS)             ==
!     ==========================================================================
      XVAL(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
          XVAL(:,ISPIN,LM)=RHO(:,LM,ISPIN)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN          ! THIS LOOP PUTS T,T->3; S,S->4 ;T,S->5
          DO ISPIN2=ISPIN1,1,-1    ! AND ASSURES CONSISTENCY WITH NSPIN
            II=II+1
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              XVAL(:,II,1)=XVAL(:,II,1) &
        &         +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
        &         +XX*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                      +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY                 ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE DFT')
      WORK1(:)=0.D0
      XDER(:,:,:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ================================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ================================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)
!       == NOW CALCULATE ENERGY DENSITY AND DERIAVTIVES ========================
        WORK1(IR)=FOURPI*EXC1
        XDER(IR,:,1)  =FOURPI*VXC5(:)*Y0
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              WORK1(IR)=WORK1(IR) &
       &              +0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &              +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,WORK1(:)*R(:)**2,EXC)
!
!     ==========================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                             ==
!     ==========================================================================
      ALLOCATE(VGRHO(NR,LMRX,NSPIN))
      VRHO(:,:,:)=0.D0
      VGRHO(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          VRHO(:,LM,ISPIN)=XDER(:,ISPIN,LM)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN
          DO ISPIN2=ISPIN1,1,-1
            II=II+1
!           == FIRST RESOLVE XVAL(:,II,1) ======================================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
!             == THE FOLLOWING LINES DIVIDE BY ZERO IF R=0. ====================
!             == THE VALUES ARE OVERWRITTEN AT THE END OF THE ROUTINE... =======
              VRHO(:,LM,ISPIN1)  =VRHO(:,LM,ISPIN1) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN2)
              VRHO(:,LM,ISPIN2)  =VRHO(:,LM,ISPIN2) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN1) =VGRHO(:,LM,ISPIN1) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN2)
              VGRHO(:,LM,ISPIN2) =VGRHO(:,LM,ISPIN2) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN1)
            ENDDO
!           == NOW RESOLVE XVAL(:,II,LM) =======================================
            DO LM=2,LMRX
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              VGRHO(:,1,ISPIN1) =VGRHO(:,1,ISPIN1) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2) =VGRHO(:,1,ISPIN2) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN2)=VGRHO(:,LM,ISPIN2) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
              VGRHO(:,LM,ISPIN1)=VGRHO(:,LM,ISPIN1) &
      &                 +XX*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
            ENDDO
          ENDDO
        ENDDO               
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==  V = V -1/R**2 D/DR [ R**2 VGRHO ]                                   ==
!     ==  V = V -[2/R VGRHO+ D/DR VGRHO ]                                     ==
!     ==========================================================================
      IF(TGRA) THEN
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
!           == FIRST ALTERNATIVE
!           CALL RADIAL$DERIVE(GID,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(GID,NR,WORK2,WORK1)
            WORK1(2:)=WORK1(2:)/R(2:)**2
            WORK(1)=WORK(2)  ! AVOID DIVIDE BY ZERO
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(VGRHO)
!
!     ==========================================================================
!     ==   CORRECT FOR DIVERGENCE AT THE ORIGIN:                              ==
!     ==   IF A SHIFTED LOGARITHMIC GRID IS USED THE FIRST GRID POINT         ==
!     ==   IS MESSED UP BECAUSE OF FACTORS 1/R                                ==
!     ==========================================================================
      IF(R(1).LT.1.D-5) THEN
        VRHO(1,:,:)=VRHO(2,:,:)
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==========================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=VRHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        CALL AUGMENTATION_NCOLLTRANS(GID,'POT',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF     
      DEALLOCATE(RHO)
      DEALLOCATE(GRHO)
      DEALLOCATE(VRHO)
                      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_XC_PREV(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE                          **
!     **  SPHERICAL CONTRIBUTION OF THE DENSITY UP TO QUADRATIC               **
!     **  ORDER IN THE NON-SPHERICAL CONTRIBUTIONS                            **
!     **                                                                      **
!     **  EXC = EXC(XVAL(L=0)*Y0)                                             **
!     **      + 0.5 * D2[EXC]/D[XVAL(L=0)*Y0]**2 * XVAL(L>0)**2               **
!     **                                                                      **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)              **
!     **  IS AN SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.             **
!     **                                                                      **
!     **  DEPENDECIES:                                                        **
!     **    DFT                                                               **
!     **    TIMING                                                            **
!     **    TRACE                                                             **
!     **                                                                      **
!     **  REMARKS: THE GRADIENTS ARE CORRECT ONLY IF DFT SUPPORTS             **
!     **    THIRD DERIVATIVES OF THE XC-ENERGY                                **
!     **   - WHEN USING SELFTEST ON THIS ROUTINE, THEN                        **
!     **     D(EXC)/DRHO(I)=POT(I)*DEX*R(I)**3                                **
!     **     AND THE VALUES AT LARGE RADII MUST BE SURPRESSED                 **
!     **                                                                      **
!     **  REMARK: FOR A COLLINEAR DENSITY THE ROUTINE GIVES DIFFERENT RESULTS **
!     **          WITH NDIMD=2 AND NDIMD=4 DUE TO THE TAYLOR EXPANSION IN     **
!     **          ANGULAR MOMENTUM EXPANSIONS                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1996 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD      ! CAN BE 1,2,4
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: B(:,:)    
      REAL(8)               :: B0INV(NR) 
      REAL(8)   ,ALLOCATABLE:: C(:,:)    
      REAL(8)   ,ALLOCATABLE:: VB(:,:)    
      REAL(8)   ,ALLOCATABLE:: VC(:,:)    
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      REAL(8)               :: FOURPI
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM,LM1,LM2,LM3
      INTEGER(4)            :: IMAX
      REAL(8)               :: FAC
      REAL(8)               :: CG0LL
      REAL(8)               :: CG
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
      REAL(8)               :: WORK3(NR)
      REAL(8)   ,PARAMETER  :: R8SMALL=1.D-20
!     **************************************************************************
      CALL TRACE$PUSH('AUGMENTATION_XC')
      EXC=0.D0
      VXC(:,:,:)=0.D0
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
      FOURPI=4.D0*PI
      CG0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  OBTAIN SPIN DENSITY                                                 ==
!     ==========================================================================
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      ALLOCATE(RHO(NR,LMRX,NSPIN))
      ALLOCATE(GRHO(NR,LMRX,NSPIN))
      ALLOCATE(VRHO(NR,LMRX,NSPIN))
      RHO(:,:,1)=RHOIN(:,:,1)
      IF(NDIMD.EQ.2) THEN
        RHO(:,:,2)=RHOIN(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == HERE WE NEED TO CALCULATE THE ABSOLUTE VALUE OF THE SPIN DENSITY ====
!       == IN AN ANGULAR MOMENTUM EXPANSION. THIS IS ONLY POSSIBLE APPROXIMATELY
!       == USING A TAYLOR EXPANSION ABOUT THE SPHERICAL PART OF THE SQUARE =====
!       == OF THE SPIN DENSITY. ================================================
        ALLOCATE(B(NR,LMRX))
!       == EVALUATE B, THE SQUARE OF THE SPIN DENSITY ==========================
        B(:,:)=0.D0
        DO LM2=1,LMRX
          DO LM3=LM2,LMRX
            WORK(:)=RHOIN(:,LM2,2)*RHOIN(:,LM3,2) &
     &             +RHOIN(:,LM2,3)*RHOIN(:,LM3,3) &
     &             +RHOIN(:,LM2,4)*RHOIN(:,LM3,4)
            IF(LM2.NE.LM3)WORK(:)=2.D0*WORK(:)
            DO LM1=1,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              B(:,LM1)=B(:,LM1)+CG*WORK(:)
            ENDDO
          ENDDO
        ENDDO
        B0INV=1.D0/(B(:,1)+R8SMALL)
!       ==  TRANSFORM B INTO C=B/(2*B0*Y0) =====================================
!       ==  THE SPHERICAL PART OF C REAMAINS ZERO ==============================
        ALLOCATE(C(NR,LMRX))
        C(:,:)=0.D0
        WORK(:)=0.5D0*B0INV(:)/Y0
        DO LM1=2,LMRX
          C(:,LM1)=B(:,LM1)*WORK(:)
        ENDDO
!       == CALCULATE SPIN DENSITY ==============================================
        RHO(:,:,2)=0.D0
        DO LM2=2,LMRX
          DO LM3=LM2,LMRX
            WORK(:)=C(:,LM2)*C(:,LM3)
            IF(LM2.NE.LM3) WORK(:)=2.D0*WORK(:)
            DO LM1=1,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              RHO(:,LM1,2)=RHO(:,LM1,2)+CG*WORK(:)
            ENDDO
          ENDDO
        ENDDO
        RHO(:,:,2)=-0.5D0*RHO(:,:,2)
        RHO(:,1,2)=RHO(:,1,2)+1.D0/Y0
        DO LM1=2,LMRX
          RHO(:,LM1,2)=RHO(:,LM1,2)+C(:,LM1)
        ENDDO
        WORK(:)=SQRT(B(:,1)*Y0)
        DO LM1=1,LMRX
          RHO(:,LM1,2)=WORK(:)*RHO(:,LM1,2)
        ENDDO
      END IF
!
!     == IMAX ALLOWS TO RESTRICT SOME LOOPS (1:5) TO (1:IMAX)
      IF(TGRA) THEN
        IF(NSPIN.EQ.2) THEN; IMAX=5; ELSE; IMAX=3; END IF
      ELSE 
        IF(NSPIN.EQ.2) THEN; IMAX=2; ELSE; IMAX=1; END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE RADIAL GRADIENT OF THE DENSITY                            ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE GRADIENTS')
      IF(TGRA) THEN
        GRHO(:,:,:)=0.D0
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
          ENDDO
        ENDDO
      ELSE
        GRHO(:,:,:)=0.D0
      END IF
!
!     ==========================================================================
!     ==  DEFINE VECTOR (RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS)             ==
!     ==========================================================================
      XVAL(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
          XVAL(:,ISPIN,LM)=RHO(:,LM,ISPIN)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN          ! THIS LOOP PUTS T,T->3; S,S->4 ;T,S->5
          DO ISPIN2=ISPIN1,1,-1    ! AND ASSURES CONSISTENCY WITH NSPIN
            II=II+1
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              XVAL(:,II,1)=XVAL(:,II,1) &
        &         +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
        &         +0.5D0*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                      +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY                 ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE DFT')
      WORK1(:)=0.D0
      XDER(:,:,:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ================================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ================================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)
!       == NOW CALCULATE ENERGY DENSITY AND DERIAVTIVES ========================
        WORK1(IR)=FOURPI*EXC1
        XDER(IR,:,1)  =FOURPI*VXC5(:)*Y0
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              WORK1(IR)=WORK1(IR) &
       &              +0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &              +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,WORK1(:)*R(:)**2,EXC)
!
!     ==========================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                             ==
!     ==========================================================================
      ALLOCATE(VGRHO(NR,LMRX,NSPIN))
      VRHO(:,:,:)=0.D0
      VGRHO(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          VRHO(:,LM,ISPIN)=XDER(:,ISPIN,LM)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN
          DO ISPIN2=ISPIN1,1,-1
            II=II+1
!           == FIRST RESOLVE XVAL(:,II,1) ======================================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              VRHO(:,LM,ISPIN1)  =VRHO(:,LM,ISPIN1) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN2)
              VRHO(:,LM,ISPIN2)  =VRHO(:,LM,ISPIN2) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN1) =VGRHO(:,LM,ISPIN1) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN2)
              VGRHO(:,LM,ISPIN2) =VGRHO(:,LM,ISPIN2) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN1)
            ENDDO
!           == NOW RESOLVE XVAL(:,II,LM) =======================================
            DO LM=2,LMRX
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              VGRHO(:,1,ISPIN1) =VGRHO(:,1,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2) =VGRHO(:,1,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN2)=VGRHO(:,LM,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
              VGRHO(:,LM,ISPIN1)=VGRHO(:,LM,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
            ENDDO
          ENDDO
        ENDDO               
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==  V = V -1/R**2 D/DR [ R**2 VGRHO ]                                   ==
!     ==  V = V -[2/R VGRHO+ D/DR VGRHO ]                                     ==
!     ==========================================================================
      IF(TGRA) THEN
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
!           == FIRST ALTERNATIVE
!           CALL RADIAL$DERIVE(GID,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(GID,NR,WORK2,WORK1)
            WORK1(:)=WORK1(:)/R(:)**2
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(VGRHO)
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==========================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=VRHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == CALCULATE VC=DE/DC===================================================
        ALLOCATE(VC(NR,LMRX))
        VC(:,:)=0.D0
!       ---- L.NEQ.0 -----------------------------------------------------------
        DO LM2=1,LMRX
          DO LM3=2,LMRX
            WORK(:)=VRHO(:,LM2,2)*C(:,LM3)
            DO LM1=2,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              VC(:,LM1)=VC(:,LM1)+CG*WORK(:)
            ENDDO
          ENDDO
        ENDDO
        WORK(:)=SQRT(B(:,1)*Y0)
        DO LM1=2,LMRX
          VC(:,LM1)=WORK(:)*(VRHO(:,LM1,2)-VC(:,LM1))
        ENDDO
!       ----- L=0 --------------------------------------------------------------
        WORK(:)=0.5D0*B0INV
        DO LM1=1,LMRX
          VC(:,1)=VC(:,1)+WORK(:)*VRHO(:,LM1,2)*RHO(:,LM1,2)
        ENDDO 
!       == CALCULATE VB=DE/DB ==================================================
        ALLOCATE(VB(NR,LMRX))
        VB(:,:)=0.D0
!       ----- L=0 --------------------------------------------------------------
        DO LM1=2,LMRX
          VB(:,1)=VB(:,1)-VC(:,LM1)*C(:,LM1)
        ENDDO
        VB(:,1)=VC(:,1)-B0INV(:)*VB(:,1)
!       ---- L.NEQ.0 -----------------------------------------------------------
        WORK(:)=0.5D0*B0INV(:)/Y0
        DO LM1=2,LMRX
          VB(:,LM1)=VC(:,LM1)*WORK(:)
        ENDDO      
!       == TRANSFORM VB=DE/DB INTO POTENTIAL FOR SPIN DENSITY ==================
        DO LM1=1,LMRX
          VB(:,LM1)=2.D0*VB(:,LM1)
          VXC(:,LM1,2)=0.D0
          VXC(:,LM1,3)=0.D0
          VXC(:,LM1,4)=0.D0
        ENDDO
        DO LM2=1,LMRX
          DO LM3=1,LMRX
            WORK1(:)=VB(:,LM2)*RHOIN(:,LM3,2)
            WORK2(:)=VB(:,LM2)*RHOIN(:,LM3,3)
            WORK3(:)=VB(:,LM2)*RHOIN(:,LM3,4)
            DO LM1=1,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              VXC(:,LM1,2)=VXC(:,LM1,2)+CG*WORK1(:)
              VXC(:,LM1,3)=VXC(:,LM1,3)+CG*WORK2(:)
              VXC(:,LM1,4)=VXC(:,LM1,4)+CG*WORK3(:)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(B)
        DEALLOCATE(VB)
        DEALLOCATE(C)
        DEALLOCATE(VC)
      END IF     
      DEALLOCATE(RHO)
      DEALLOCATE(GRHO)
      DEALLOCATE(VRHO)
!
!     ==========================================================================
!     ==   CORRECT FOR DIVERGENCE AT THE ORIGIN:                              ==
!     ==   IF A SHIFTED LOGARITHMIC GRID IS USED THE FIRST GRID POINT         ==
!     ==   IS MESSED UP BECAUSE OF FACTORS 1/R                                ==
!     ==========================================================================
      IF(R(1).LT.1.D-5) THEN
        VXC(1,:,:)=VXC(2,:,:)
      END IF
                      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS(GID,ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO SECOND ORDER. THUS THE RESULTS      **
!     **  FOR COLLINEAR DENSITIES DO NOT AGREE WITH THOSE OF COLLINEAR        **
!     **  TREATED AS NONCOLLINEAR ONES                                        **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)    :: GID
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)     ,PARAMETER     :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER     :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                    :: Q(NR)
      REAL(8)                    :: VQ(NR)
      REAL(8)                    :: A(NR,LMRX,3)
      REAL(8)                    :: VA(NR,LMRX,3)
      REAL(8)                    :: P(NR,LMRX)
      REAL(8)                    :: VP(NR,LMRX)
      REAL(8)                    :: S(NR,LMRX)
      REAL(8)                    :: VS(NR,LMRX)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: CG   ! GAUNT COEFFICIENT
      REAL(8)     ,PARAMETER     :: SMALL=(1.D-2)**2
      INTEGER(4)                 :: ISIG,LM,LM1,LM2,LM3
      LOGICAL(4)  ,PARAMETER     :: TOTHER=.TRUE.
!     **************************************************************************
      IF(TOTHER) THEN
        CALL AUGMENTATION_NCOLLTRANS_OTHER3(GID,ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z_OTHER.DAT',GID,NR,LMRX,RHO4(:,:,4))
!!$CALL AUGMENTATION_WRITEPHI('RHO2_Z_OTHER.DAT',GID,NR,LMRX,RHO2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT2_Z_OTHER.DAT',GID,NR,LMRX,POT2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT4_Z_OTHER.DAT',GID,NR,LMRX,POT4(:,:,4))
        RETURN
      END IF
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT AUXILIARY VARIABLES                                         ==
!     ==========================================================================
      Q(:)=SQRT(SMALL+(RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)*Y0**2)
!
!     == CONSTRUCT A=SIGMA/Q ===================================================
      DO ISIG=1,3
        DO LM=1,LMRX
          A(:,LM,ISIG)=RHO4(:,LM,ISIG+1)/Q(:)
        ENDDO
      ENDDO
!
!     ==  CONSTRUCT P=A^2 ======================================================
      P(:,:)=0.D0
      DO LM1=1,LMRX
        DO LM2=LM1,LMRX
          AUX(:)=A(:,LM1,1)*A(:,LM2,1) &
     &          +A(:,LM1,2)*A(:,LM2,2) &
     &          +A(:,LM1,3)*A(:,LM2,3)
          IF(LM1.NE.LM2)AUX(:)=2.D0*AUX(:)
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            P(:,LM3)=P(:,LM3)+CG*AUX(:)
          ENDDO
        ENDDO
      ENDDO
!
!     == CONSTRUCT S ===========================================================
      S(:,:)=0.D0
      DO LM1=1,LMRX
        S(:,LM1)=S(:,LM1)+0.75D0*P(:,LM1)
        DO LM2=LM1,LMRX
          AUX(:)=P(:,LM1)*P(:,LM2)
          IF(LM1.NE.LM2)AUX(:)=2.D0*AUX(:)
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            S(:,LM3)=S(:,LM3)-0.125D0*CG*AUX(:)
          ENDDO
        ENDDO
      ENDDO
      S(:,1)=S(:,1)+0.375D0/Y0
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
      RHO2(:,:,1)=RHO4(:,:,1)
      RHO2(:,1,2)=0.D0
      DO LM=1,LMRX
        RHO2(:,LM,2)=RHO2(:,LM,2)+Q(:)*S(:,LM)
      ENDDO
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
!
!     == CONSTRUCT VQ ==========================================================
      VQ(:)=0.D0
      DO LM=1,LMRX
        VQ(:)=VQ(:)+POT2(:,LM,2)*S(:,LM)
        VS(:,LM)=POT2(:,LM,2)*Q(:)
      ENDDO
!
!     == CONSTRUCT VP ==========================================================
      DO LM1=1,LMRX
        VP(:,LM1)=0.75D0*VS(:,LM1)
        DO LM2=1,LMRX
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            VP(:,LM1)=VP(:,LM1)-0.25D0*CG*VS(:,LM2)*P(:,LM3)
          ENDDO
        ENDDO
      ENDDO
!
!     == CONSTRUCT VA ==========================================================
      VA(:,:,:)=0.D0
      DO ISIG=1,3
        DO LM1=1,LMRX
          DO LM2=1,LMRX
            DO LM3=1,LMRX
              CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
              VA(:,LM1,ISIG)=VA(:,LM1,ISIG)+2.D0*CG*VP(:,LM2)*A(:,LM3,ISIG)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     == CONSTRUCT POTENTIAL FOR NON-COLLINEAR SPIN DENSITY ====================
      AUX(:)=0.D0
      DO ISIG=1,3
        DO LM=1,LMRX
          AUX(:)=AUX(:)+VA(:,LM,ISIG)*A(:,LM,ISIG)
        ENDDO
      ENDDO
      AUX(:)=(VQ(:)-AUX(:)/Q(:))*Y0**2
!
      POT4(:,:,:)=0.D0
      DO ISIG=1,3
        POT4(:,1,ISIG+1)=POT4(:,1,ISIG)+AUX(:)*A(:,1,ISIG)
        DO LM1=1,LMRX
          POT4(:,LM1,ISIG+1)=POT4(:,LM1,ISIG+1)+VA(:,LM1,ISIG)/Q(:)
        ENDDO
      ENDDO
!
!     == ADD TOTAL POTENTIAL ===================================================
      POT4(:,:,1)=POT2(:,:,1)
!!$CALL AUGMENTATION_WRITEPHI('Q.DAT',GID,NR,1,Q)
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z.DAT',GID,NR,LMRX,RHO4(:,:,4))
!!$CALL AUGMENTATION_WRITEPHI('RHO2_Z.DAT',GID,NR,LMRX,RHO2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT2_Z.DAT',GID,NR,LMRX,POT2(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('POT4_Z.DAT',GID,NR,LMRX,POT4(:,:,4))
!!$STOP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS_OTHER1(ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO FIRST ORDER ONLY. IT IS CONSTRUCTED **
!     **  SUCH THAT THE RESULTS FOR COLLINEAR DENSITIES DO AGREE WITH         **
!     **  THOSE OF COLLINEAR DENSITIES TREATED AS NONCOLLINEAR ONES           **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)     ,PARAMETER     :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER     :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                    :: A(NR,LMRX,3)
      REAL(8)                    :: VA(NR,LMRX,3)
      REAL(8)                    :: Q(NR)
      REAL(8)                    :: VQ(NR)
      REAL(8)     ,PARAMETER     :: SMALL=1.D-4 ! DO NOT CHOOSE TOO SMALL!
      INTEGER(4)                 :: ISIG,LM
!     **************************************************************************
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT AUXILIARY VARIABLES                                         ==
!     ==========================================================================
      Q(:)=SQRT(SMALL+RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)
      DO ISIG=1,3
        DO LM=1,LMRX
          A(:,LM,ISIG)=RHO4(:,LM,ISIG+1)/Q(:)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
      RHO2(:,:,1)=RHO4(:,:,1)
      RHO2(:,1,2)=Q(:)
      DO LM=2,LMRX
        RHO2(:,LM,2)=0.D0
        DO ISIG=1,3
          RHO2(:,LM,2)=RHO2(:,LM,2)+A(:,1,ISIG)*A(:,LM,ISIG)
        ENDDO
        RHO2(:,LM,2)=RHO2(:,LM,2)*Q(:)
      ENDDO
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
!     -- CALCULATE VQ ----------------------------------------------------------
      VQ(:)=0.D0
      DO LM=1,LMRX
        VQ(:)=VQ(:)+POT2(:,LM,2)*RHO2(:,LM,2)/Q(:)
      ENDDO
!     -- CALCULATE VA ----------------------------------------------------------
      DO ISIG=1,3
        VA(:,1,ISIG)=0.D0
        DO LM=2,LMRX
          VA(:,1,ISIG)=VA(:,1,ISIG)+POT2(:,LM,2)*A(:,LM,ISIG)
          VA(:,LM,ISIG)=POT2(:,LM,2)*A(:,1,ISIG)
        ENDDO
      ENDDO
!     -- NEW VQ ----------------------------------------------------------------
      DO ISIG=1,3
        DO LM=1,LMRX
          VQ(:)=VQ(:)-VA(:,LM,ISIG)*A(:,LM,ISIG)
        ENDDO
      ENDDO
!     --------------------------------------------------------------------------
      POT4(:,:,1)=POT2(:,:,1)
      DO ISIG=1,3
        DO LM=1,LMRX
          POT4(:,LM,ISIG+1)=VA(:,LM,ISIG)
        ENDDO
        POT4(:,1,ISIG+1)=POT4(:,1,ISIG+1)+VQ(:)*A(:,1,ISIG)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS_OTHER2(ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO FIRST ORDER ONLY. IT IS CONSTRUCTED **
!     **  SUCH THAT THE RESULTS FOR COLLINEAR DENSITIES DO AGREE WITH         **
!     **  THOSE OF COLLINEAR DENSITIES TREATED AS NONCOLLINEAR ONES           **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)     ,PARAMETER     :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER     :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)     ,PARAMETER     :: SMALL=1.D-4 ! DO NOT CHOOSE TOO SMALL!
      INTEGER(4)                 :: ISIG,LM
!     **************************************************************************
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
      RHO2(:,:,1)=RHO4(:,:,1)
      RHO2(:,1,2)=SQRT(RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)
      DO LM=2,LMRX
        RHO2(:,LM,2)=Y0*(RHO4(:,1,2)*RHO4(:,LM,2) &
     &                  +RHO4(:,1,3)*RHO4(:,LM,3) &
     &                  +RHO4(:,1,4)*RHO4(:,LM,4)) 
      ENDDO
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
      POT4(:,:,1)=POT2(:,:,1)
      DO ISIG=2,4
        POT4(:,1,ISIG)=POT2(:,1,2)*RHO4(:,1,ISIG)/(SMALL+RHO2(:,1,2))
      ENDDO
      DO LM=2,LMRX
        DO ISIG=2,4
           POT4(:,LM,ISIG)=POT2(:,LM,2)*RHO4(:,1,ISIG)*Y0
           POT4(:,1,ISIG)=POT4(:,1,ISIG)+POT2(:,LM,2)*RHO4(:,LM,ISIG)*Y0
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_NCOLLTRANS_OTHER3(GID,ID,NR,LMRX,RHO4,RHO2,POT2,POT4)
!     **************************************************************************
!     **  CONSTRUCTS A COLLINEAR SPIN DENSITY FROM A NONCOLLINEAR ONE         **
!     **  AND FOR  ID='POT' INSTEAD OF ID='RHO' IT ALSO                       **
!     **  CONSTRUCTS A NONCOLLINEAR POTENTIAL FROM A COLLINEAR ONE AND THE    **
!     **  NONCOLLINEAR DENSITY                                                **
!     **                                                                      **
!     **  THE TRANSFORMATION IS APPROXIMATE BECAUSE IT INVOLVES A TAYLOR      **
!     **  EXPANSION OF THE SQUARE ROOT TO FIRST ORDER ONLY. IT IS CONSTRUCTED **
!     **  SUCH THAT THE RESULTS FOR COLLINEAR DENSITIES DO AGREE WITH         **
!     **  THOSE OF COLLINEAR DENSITIES TREATED AS NONCOLLINEAR ONES           **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: GID
      INTEGER(4)  ,INTENT(IN)    :: NR
      INTEGER(4)  ,INTENT(IN)    :: LMRX
      REAL(8)     ,INTENT(IN)    :: RHO4(NR,LMRX,4)
      REAL(8)     ,INTENT(OUT)   :: RHO2(NR,LMRX,2)
      REAL(8)     ,INTENT(IN)    :: POT2(NR,LMRX,2)
      REAL(8)     ,INTENT(OUT)   :: POT4(NR,LMRX,4)
      REAL(8)     ,PARAMETER     :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER     :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                    :: ARRAY(NR)
      REAL(8)                    :: A(NR,LMRX,3)
      REAL(8)                    :: B(NR,LMRX,3)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: C0LL,CG
      REAL(8)     ,PARAMETER     :: SMALL=1.D-6 ! DO NOT CHOOSE TOO SMALL!
      INTEGER(4)                 :: ISIG,LM,LM1,LM2,LM3
!     **************************************************************************
      C0LL=Y0
      RHO2(:,:,:)=0.D0
      POT4(:,:,:)=0.D0
!
!     ==========================================================================
!     == WORK OUT COLLINEAR DENSITY                                           ==
!     ==========================================================================
!     == TOTAL DENSITY =========================================================
      RHO2(:,:,1)=RHO4(:,:,1)
!
!     == SQUARE OF THE SPIN DENSITY=============================================
      DO LM1=1,LMRX
        AUX(:)=0.D0
        DO LM2=1,LMRX
          DO LM3=LM2,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            IF(LM2.NE.LM3)CG=2.D0*CG
            AUX(:)=AUX(:)+CG*(RHO4(:,LM2,2)*RHO4(:,LM3,2) &
     &                       +RHO4(:,LM2,3)*RHO4(:,LM3,3) &
     &                       +RHO4(:,LM2,4)*RHO4(:,LM3,4)) 
          ENDDO
        ENDDO
        RHO2(:,LM1,2)=AUX(:)
      ENDDO
!
!     == ABSOLUTE VALUE OF THE SPIN DENSITY=====================================
      ARRAY(:)=1.D0/SQRT(SMALL+RHO2(:,1,2)/Y0)
      AUX(:)=ARRAY(:)/(2.D0*Y0)
      DO LM=1,LMRX
        RHO2(:,LM,2)=RHO2(:,LM,2)*AUX(:)
      ENDDO
      RHO2(:,1,2)=2.D0*RHO2(:,1,2)
!!$CALL AUGMENTATION_WRITEPHI('RHO4_X',GID,NR,LMRX,RHO4(:,:,2))
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Y',GID,NR,LMRX,RHO4(:,:,3))
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z',GID,NR,LMRX,RHO4(:,:,4))
!!$CALL AUGMENTATION_WRITEPHI('RHO2',GID,NR,LMRX,RHO2(:,:,2))
!!$AUX(:)=SQRT(RHO4(:,1,2)**2+RHO4(:,1,3)**2+RHO4(:,1,4)**2)
!!$CALL AUGMENTATION_WRITEPHI('AUX',GID,NR,1,AUX)
!
!     ==========================================================================
!     == RETURN IF ONLY DENSITY IS REQUIRED                                   ==
!     ==========================================================================
      IF(ID.EQ.'RHO') RETURN
!
!     ==========================================================================
!     == WORK OUT POTENTIAL                                                   ==
!     ==========================================================================
      IF(ID.NE.'POT') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATION_NCTRANS1')
      END IF
!     == TOTAL POTENTIAL =======================================================
      POT4(:,:,1)=POT2(:,:,1)
!
      DO ISIG=2,4
        DO LM=1,LMRX
          A(:,LM,ISIG-1)=RHO4(:,LM,ISIG)*ARRAY(:)
        ENDDO
      ENDDO
!
      B(:,:,:)=0.D0
      DO LM1=1,LMRX
        DO LM2=1,LMRX
          DO LM3=1,LMRX
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            IF(LM2.EQ.1)CG=2.D0*CG
            DO ISIG=2,4
              B(:,LM1,ISIG-1)=B(:,LM1,ISIG-1)+CG*POT2(:,LM2,2)*A(:,LM3,ISIG-1)        
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      AUX(:)=0.D0
      DO ISIG=2,4
        DO LM=1,LMRX
          AUX(:)=AUX(:)+B(:,LM,ISIG-1)*A(:,LM,ISIG-1)
        ENDDO
      ENDDO
!
      DO ISIG=2,4
        DO LM=1,LMRX
          POT4(:,LM,ISIG)=-0.5D0*AUX(:)*A(:,LM,ISIG-1)+B(:,LM,ISIG-1)
        ENDDO
      ENDDO
      POT4(:,:,2:)=POT4(:,:,2:)/Y0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_PSHARTREE(GID,NR,LMRX,PSRHOC,PSRHO &
     &            ,VADD,RCSM,QLM,VQLM,RHOB,PSPOT,PSEH)
!     **************************************************************************
!     **                                                                      **
!     ** CALCULATES HARTREE ENERGY                                            **
!     ** ADDS COMPENSATION DENSITY TO THE PSEUDO DENSITY                      **
!     **                                                                      **
!     **  1) AEE=AERHO*AEZ/R                                                  **
!     **     PSE=PSRHO*VADD                                                   **
!     **  2) PSRHO=PSRHO+GAUSS(R)*QLM                                         **
!     **  3) AEE=AEE + 0.5*AERHO(R)*AERHO(R')/|R-R'|                          **
!     **     PSE=PSE + 0.5*PSRHO(R)*PSRHO(R')/|R-R'|                          **
!     **     VQLM=VQLM - GAUSS(R)*PSRHO(R')/|R-R'|                            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: PSRHOC(NR)
      REAL(8)    ,INTENT(IN) :: PSRHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: VADD(NR)
      REAL(8)    ,INTENT(IN) :: RCSM
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)    ,INTENT(IN) :: RHOB
      REAL(8)    ,INTENT(OUT):: PSPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: PSEH
      REAL(8)    ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)                :: R(NR)
      REAL(8)                :: RHO1(NR)
      REAL(8)                :: POT(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: G(NR)
      REAL(8)                :: RHOHAT(NR,LMRX)
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      INTEGER(4)             :: LM,L,M,LX
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  INITIALIZE ENERGIES AND POTENTIALS TO  ZERO                 ==
!     ==================================================================
      PSE(:)=0.D0
      PSPOT(:,:)=0.D0
!
!     ==================================================================
!     ==  CONSTRUCT COMPENSATION DENSITY                              ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      LX=INT(SQRT(REAL(LMRX))-1.D0)   ! LMX=(LX+1)**2
      IF((LX+1)**2.NE.LMRX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMRX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMRX',LMRX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('AUGMENTATION_PSHARTREE')
      END IF
      LM=0
      DO L=0,LX
        CALL GAUSSN(L,ALPHA,CL)
        G(:)=CL*EXP(-ALPHA*R(:)**2)*R(:)**L
        DO M=1,2*L+1
          LM=LM+1
          RHOHAT(:,LM)=QLM(LM)*G(:)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==   ADD UNSCREENING POTENTIAL VADD                             ==
!     ==================================================================
      PSE(:)    =VADD(:)*(PSRHO(:,1)+PSRHOC(:))
      PSPOT(:,1)=PSPOT(:,1)+VADD(:)
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL GAUSSN(L,ALPHA,CL)
        SVAR=QLM(LM)*CL
        RHO1(:)=PSRHO(:,LM)+RHOHAT(:,LM)
        IF(L.EQ.0)RHO1(:)=RHO1(:)+PSRHOC(:)
        CALL RADIAL$POISSON(GID,NR,L,RHO1,POT)
        PSE(:)=PSE(:)+0.5D0*POT(:)*RHO1(:)
        PSPOT(:,LM)=PSPOT(:,LM)+POT(:)
      ENDDO
!
!     ==================================================================
!     ==  ADD EXTERNAL POTENTIAL                                      ==
!     ==================================================================
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        POT(:)=VQLM(LM)*R(:)**L
        PSPOT(:,LM)=PSPOT(:,LM)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
        RHO1(:)=PSRHO(:,LM)+RHOHAT(:,LM)
        IF(L.EQ.0) RHO1(:)=RHO1(:)+PSRHOC(:)
        PSE(:)=PSE(:)+RHO1(:)*POT(:)
      ENDDO
!
!     ==================================================================
!     ==  ADD POTENTIAL FROM THE BACKGROUND                           ==
!     ==================================================================
      SVAR=-2.D0*PI*RHOB/3.D0*SQRT(4.D0*PI)
      POT(:)=SVAR*R(:)**2
      PSPOT(:,1)=PSPOT(:,1)+POT(:) ! POTENTIAL OF THE BACKGROUND
      PSE(:)=PSE(:)+(PSRHO(:,1)+RHOHAT(:,1)+PSRHOC(:))*POT(:)
!
!     ==================================================================
!     ==  CALCULATE TOTAL ENERGY                                      ==
!     ==================================================================
      PSE(:)=PSE(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,PSE,PSEH)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_AEHARTREE(GID,NR,LMRX,AEZ,RHOC,AERHO &
     &                                 ,VQLM,RHOB,AEPOT,AEEH)
!     ******************************************************************
!     **  ELECTROSTATIC ENERGY OF THE ALL-ELECTRON ONE-CENTER DENSITY **
!     **  INCLUDING THE EXTERNAL POTENTIAL AND THE POTENTIAL OF THE   **
!     **  COMPENSATING CHARGE BACKGROUND FOR CHARGED SYSTEMS          **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: AEZ        ! ATOMIC NUMBER
      REAL(8)    ,INTENT(IN) :: RHOC(NR)   ! CORE DENSITY
      REAL(8)    ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)    ,INTENT(IN) :: RHOB       ! COMPENSATING BACKGROUND
      REAL(8)    ,INTENT(OUT):: AEPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: AEEH       ! ENERGY
      REAL(8)    ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)    ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX(NR)
      REAL(8)                :: POT(NR)
      REAL(8)                :: AEE(NR)  ! ENERGY DENSITY
      REAL(8)                :: SVAR
      REAL(8)                :: EEXT
      INTEGER(4)             :: LM
      INTEGER(4)             :: L
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  INITIALIZE ENERGIES AND POTENTIALS TO  ZERO                 ==
!     ==================================================================
      AEE(:)=0.D0
      AEPOT(:,:)=0.D0
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC VALENCE-CORE INTERACTION           ==
!     ==   AND VALENCE-NUCLEUS INTERACTION                            ==
!     ==================================================================
      CALL RADIAL$POISSON(GID,NR,0,RHOC,POT)
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX)
      POT(:)  = POT(:)+AUX(:)
      AEE(:)    = POT(:)* AERHO(:,1)
      AEPOT(:,1)= AEPOT(:,1)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
!      AEE(:)=AE(:)+RHOC(:)*VQLM(1)
!AUX(:)=AEE(:)*R(:)**2
!CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!PRINT*,'EL C-V INT', SVAR
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
!AUX(:)=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL RADIAL$POISSON(GID,NR,L,AERHO(:,LM),POT)
        AEE(:)=AEE(:)+0.5D0*POT(:)*AERHO(:,LM)
        AEPOT(:,LM)=AEPOT(:,LM)+POT(:)
!AUX(:)=AUX(:)+0.5D0*POT(:)*AERHO(:,LM)*R(:)**2
      ENDDO
!CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!PRINT*,'EL V-V INT',SVAR
!
!AUX(:)=AEE(:)*R(:)**2
!CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!PRINT*,'EL HARTREE INT ',SVAR
!
!     ==================================================================
!     ==  ADD EXTERNAL POTENTIAL                                      ==
!     ==================================================================
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        POT(:)=VQLM(LM)*R(:)**L
        AEPOT(:,LM)=AEPOT(:,LM)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
        AEE(:)=AEE(:)+AERHO(:,LM)*POT(:)
        IF(L.EQ.0)AEE(:)=AEE(:)+RHOC(:)*POT(:)
      ENDDO
      EEXT=-VQLM(1)*Y0*AEZ  ! WILL BE ADDED TO THE TOTAL ENERGY
!
!     ==================================================================
!     ==  ADD POTENTIAL FROM THE BACKGROUND                           ==
!     ==================================================================
      SVAR=-2.D0*PI*RHOB/3.D0*SQRT(4.D0*PI)
      POT(:)=SVAR*R(:)**2
      AEPOT(:,1)=AEPOT(:,1)+POT(:) ! POTENTIAL OF THE BACKGROUND
      AEE(:)=AEE(:)+(AERHO(:,1)+RHOC(:))*POT(:)
!CALL RADIAL$INTEGRAL(GID,NR,(AERHO(:,1)+RHOC(:))*POT(:)*R**2,SVAR)
!PRINT*,'EL BACKGROUND',SVAR
!
!     ==================================================================
!     ==  CALCULATE ENERGY                                            ==
!     ==================================================================
      AEE(:)=AEE(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AEE,AEEH)
      AEEH=AEEH+EEXT
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADJUSTVQLM(GID,NR,LMRX,PSRHO,RCSM,QLM,VQLM)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      REAL(8)    ,INTENT(IN) :: RCSM
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: PSRHO(NR,LMRX) ! INCLUDES PSEUDO CORE DENSITY
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(OUT):: VQLM(LMRX)
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX(NR)
      REAL(8)                :: RHO(NR)
      REAL(8)                :: G(NR)
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      INTEGER(4)             :: L,M,LM
      INTEGER(4)             :: LX
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      LX=INT(SQRT(REAL(LMRX))-1.D0) ! LMX=(LX+1)**2
      IF((LX+1)**2.NE.LMRX) THEN   
        CALL ERROR$MSG('LMRX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMRX',LMRX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('AUGMENTATION_ADJUSTVQLM')
      END IF
      LM=0
      DO L=0,LX
        CALL GAUSSN(L,ALPHA,CL)
        G(:)=CL*EXP(-ALPHA*R(:)**2)*R(:)**L
        DO M=1,2*L+1
          LM=LM+1
          RHO=PSRHO(:,LM)+QLM(LM)*G(:)
          CALL RADIAL$POISSON(GID,NR,L,RHO,AUX)
          AUX(:)=AUX(:)*G(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          VQLM(LM)=-SVAR
        ENDDO
      ENDDO
      RETURN
      END
!     ................................................................
      SUBROUTINE AUGMENTATION_FEEDHYPERFINE(GID,NR,LMRX,NDIMD &
     &         ,IAT,AERHO,AECORE,PSRHO,PSCORE,AEHPOT,PSHPOT)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: PSRHO(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(IN) :: AEHPOT(NR,LMRX)
      REAL(8)   ,INTENT(IN) :: PSHPOT(NR,LMRX)
      REAL(8)               :: RHO(NR,LMRX)
!     ****************************************************************
!
!     ================================================================
!     ==  TOTAL DENSITY                                             ==
!     ================================================================
      CALL HYPERFINE$SET1CRHO('PS','TOT',IAT,GID,NR,NR,LMRX,PSRHO)
      RHO(:,:)=AERHO(:,:,1)
      RHO(:,1)=RHO(:,1)+AECORE(:)
      CALL HYPERFINE$SET1CRHO('AE','TOT',IAT,GID,NR,NR,LMRX,RHO)
!
!     ================================================================
!     ==  SPIN DENSITY                                              ==
!     ================================================================
      IF(NDIMD.GT.1) THEN
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          RHO(:,:)=PSRHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          RHO(:,:)=SQRT(PSRHO(:,:,2)**2+PSRHO(:,:,3)**2+PSRHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('PS','SPIN',IAT,GID,NR,NR,LMRX,RHO)
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          RHO(:,:)=AERHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          RHO(:,:)=SQRT(AERHO(:,:,2)**2+AERHO(:,:,3)**2+AERHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('AE','SPIN',IAT,GID,NR,NR,LMRX,RHO)
      END IF
!
!     ================================================================
!     ==  TOTAL POTENTIAL FOR ELECTRIC FIELD GRADIENTS              ==
!     ================================================================
      CALL HYPERFINE$SET1CPOT('AE',IAT,GID,NR,NR,LMRX,AEHPOT)
      CALL HYPERFINE$SET1CPOT('PS',IAT,GID,NR,NR,LMRX,PSHPOT)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &           ,AEPOT,PSPOT,AEPHI,PSPHI,DATH)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE EXPECTATION VALUE OF                         **
!     **  THE ONE-CENTER POTENTIALS WITH THE PARTIAL WAVES            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: PSPOT(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: DATH(LMNX,LMNX,NDIMD)
      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LM1,LM2,LM3
      INTEGER(4)            :: L1,L2
      INTEGER(4)            :: IM1,IM2
      INTEGER(4)            :: ISPIN
      REAL(8)               :: AEDMU(NR,NDIMD)
      REAL(8)               :: PSDMU(NR,NDIMD)
      REAL(8)               :: DWORK1(NR)
      REAL(8)               :: CG
      REAL(8)               :: SVAR
      REAL(8)               :: R(NR)
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      DATH(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          LM1=L1**2+IM1
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
!     
!             ==========================================================
!             ==  SUM ALL POTENTIALS THAT ACT ON THE GIVEN PAIR       ==
!             ==  OF PARTIAL WAVES                                    ==
!             ==========================================================
              AEDMU(:,:)=0.D0
              PSDMU(:,:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  DO ISPIN=1,NDIMD
                    AEDMU(:,ISPIN)=AEDMU(:,ISPIN)+CG*AEPOT(:,LM3,ISPIN)
                    PSDMU(:,ISPIN)=PSDMU(:,ISPIN)+CG*PSPOT(:,LM3,ISPIN)
                  ENDDO
                END IF
              ENDDO
!     
!             ==========================================================
!             ==  PERFORM NOW THE INTEGRATION                         ==
!             ==========================================================
              DO ISPIN=1,NDIMD
                DWORK1(:)= &
     &               (AEDMU(:,ISPIN)*AEPHI(:,LN1)*AEPHI(:,LN2) &
     &               -PSDMU(:,ISPIN)*PSPHI(:,LN1)*PSPHI(:,LN2))*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,SVAR)
                DATH(LMN1,LMN2,ISPIN)=DATH(LMN1,LMN2,ISPIN)+SVAR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE AUGMENTATION_CORELEVELS(NDIMD,R1,DEX,NR,LMRX,AEPOT)
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: LMRX
!!$      INTEGER(4),INTENT(IN) :: NDIMD
!!$      REAL(8)   ,INTENT(IN) :: R1
!!$      REAL(8)   ,INTENT(IN) :: DEX
!!$      INTEGER(4),INTENT(IN) :: NR
!!$      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
!!$      INTEGER(4)            :: NC         ! NUMBER OF CORE STATES
!!$      INTEGER(4),ALLOCATABLE:: LOFC(:)    !(NC)
!!$      REAL(8)   ,ALLOCATABLE:: EOFC(:)    !(NC)
!!$      REAL(8)   ,ALLOCATABLE:: VATOM(NR)  !(NC)
!!$      REAL(8)   ,ALLOCATABLE:: PHIC(:,:)  !(NR,NC)
!!$      REAL(8)   ,ALLOCATABLE:: R2(:)   ! R**2
!!$      REAL(8)   ,ALLOCATABLE:: H(:,:)  ! HAMILTONIAN
!!$      COMPLEX(8),ALLOCATABLE:: UC(:,:) ! CORE EIGENVECTORS
!!$      REAL(8)   ,ALLOCATABLE:: EIGC(:) ! CORE EIGENVALUES
!!$      INTEGER(4)            :: NE ! DIMENSIONA OF THE HAMILTONIAN
!!$      INTEGER(4)            :: I1,I2,M1,M2,L1,L2,LM1,LM2,LM3
!!$!     ******************************************************************
!!$      IF(NDIMD.NE.1) THEN
!!$        CALL ERROR$STOP('AUGMENTATION_CORELEVELS')
!!$      END IF
!!$      PI=4.D0*ATAN(1.D0)
!!$      C000=1/SQRT(4.D0*PI)
!!$!
!!$!     ==================================================================
!!$!     ==  READ CORE LEVELS                                            ==
!!$!     ==================================================================
!!$!===  READ NC 
!!$      ALLOCATE(LOFC(NC))
!!$      ALLOCATE(PHIC(NR,NC))
!!$!
!!$!     ==================================================================
!!$!     ==                                                              ==
!!$!     ==================================================================
!!$      NE=0      
!!$      DO I=1,NC
!!$        NE=NE+2*LC(I)+1
!!$      ENDDO
!!$      ALLOCATE(H(NC,NC))
!!$      ALLOCATE(EIGC(NC))
!!$      ALLOCATE(UCC(NC,NC))
!!$!
!!$!     ==================================================================
!!$!     ==  RADIAL GRID R2=R**2                                         ==
!!$!     ==================================================================
!!$      XEXP=EXP(DEX)
!!$      RI=R1/XEXP
!!$      DO IR=1,NR
!!$        RI=RI*XEXP
!!$        R2(IR)=RI**2
!!$      ENDDO
!!$!
!!$!     ==================================================================
!!$!     ==  CALCULATE HAMILTONIAN                                       ==
!!$!     ==================================================================
!!$      DO I1=1,NC
!!$        L1=LOFC(I1)
!!$        DO I2=I1,NC
!!$          L2=LOFC(I2)
!!$          AUX(:)=PHI(:,I1)*PHI(:,I2)*R2(:)
!!$          DO M1=1,2*L1+1
!!$            LM1=L1**2+M1
!!$            DO M2=1,2*L2+1
!!$              LM2=L2**2+M2
!!$              SUM(:)=0.D0
!!$              DO LM3=1,LMRX
!!$                CALL CLEBSCH(LM1,LM2,LM3,CG)
!!$                IF(CG.EQ.0.D0) CYCLE
!!$                SUM(:)=SUM(:)+CG*AEPOT(:,LM3,1)
!!$              ENDDO
!!$              IF(LM1.EQ.LM2) SUM(:)=SUM(:)-VC(:)*C000
!!$              SUM(:)=SUM(:)*AUX(:)
!!$              CALL OLDRADIAL$INTEGRAL(R1,DEX,NR,SUM,H(IE1,IE2))
!!$              H(IE2,IE1)=H(IE1,IE2)
!!$! ADD EIGENVALUES
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==================================================================
!!$!     ==  DIAGONALIZE HAMILTONIAN                                     ==
!!$!     ==================================================================
!!$      CALL LIB$DIAGR8(NE,H,EIGC,UC)
!!$!
!!$!     ==================================================================
!!$!     ==  PRINT EIGENVALUES                                           ==
!!$!     ==================================================================
!!$      CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$      WRITE(NFILO,FMT='("CORE EIGENVALUES"//"================")')
!!$      WRITE(NFILO,FMT='(10F10.3)')EIGC
!!$      RETURN
!!$      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE EXPERTNAL1CPOT_MODULE
!*******************************************************************************
!**                                                                           **
!**  APPLIES AN EXTERNAL POTENTIAL ACTING ON THE ELECTRONS                    **
!**  WITHIN THE AUGMENTATION REGION                                           **
!**                                                                           **
!**    IDIMD EXTERNALLY THIS INDEX REFERS TO                                  **
!**             0: TOTAL  FOR NDIMD=1,2,4                                     **
!**             1: SZ FOR NDIMD=2                                             **
!**             1: SX FOR NDIMD=4                                             **
!**             2: SY FPR NDIMD=4                                             **
!**             3: SZ FPR NDIMD=4                                             **
!**          INTERNALLY THE VALUE IS RAISED BY ONE SO THAT IT                 **
!**          CORRESPONDS TO THE INDEXING                                      **
!**                                                                           **
!*******************************************************************************
TYPE EXTPOT
  REAL(8)          :: VALUE
  CHARACTER(LEN=32):: ATOM
  REAL(8)          :: RC  
  REAL(8)          :: PWR
  CHARACTER(LEN=32):: TYPE   ! CAN BE 'S','P','D','F','ALL'
  INTEGER(4)       :: IDIMD  ! IDIMD=0 REFERES TO ALL SPIN DIRECTIONS
END TYPE EXTPOT
INTEGER(4)             :: NPOT=0
INTEGER(4)             :: NPOTX=0
INTEGER(4),PARAMETER   :: DNPOT=5
TYPE (EXTPOT), ALLOCATABLE :: POT(:)
!*******************************************************************************
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CREATE
      TYPE (EXTPOT),ALLOCATABLE :: TMPPOT(:)
!     **************************************************************************
      IF(NPOT.LT.NPOTX) RETURN
      IF(ALLOCATED(POT)) THEN
        ALLOCATE(TMPPOT(NPOTX))
        TMPPOT=POT
        DEALLOCATE(POT)
        NPOTX=NPOTX+DNPOT
        ALLOCATE(POT(NPOTX))
        POT(1:NPOT)=TMPPOT(1:NPOT)
        DEALLOCATE(TMPPOT)
      ELSE
        NPOTX=DNPOT
        ALLOCATE(POT(NPOTX))
      END IF          
      END SUBROUTINE CREATE
END MODULE EXPERTNAL1CPOT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTERNAL1CPOT$SETPOT(ATOM,TYPE,IDIMD,VALUE,RC,PWR)
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      REAL(8)          ,INTENT(IN) :: VALUE
      CHARACTER(LEN=32),INTENT(IN) :: ATOM
      CHARACTER(LEN=32),INTENT(IN) :: TYPE
      INTEGER(4)       ,INTENT(IN) :: IDIMD !SPIN DIRECTION OR 0
      REAL(8)          ,INTENT(IN) :: RC
      REAL(8)          ,INTENT(IN) :: PWR
!     **************************************************************************
      CALL CREATE
      NPOT=NPOT+1
      POT(NPOT)%ATOM =ATOM
      POT(NPOT)%VALUE=VALUE
      POT(NPOT)%RC=RC
      POT(NPOT)%PWR=PWR
      POT(NPOT)%IDIMD=IDIMD+1
      POT(NPOT)%TYPE =TYPE
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTERNAL1CPOT$REPORT(NFIL)
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IPOT
!     **************************************************************************
      IF(NPOT.EQ.0) RETURN
!     CALL REPORT$TITLE(NFIL,"EXTERNAL POTENTIALS ON ORBITALS")
      WRITE(NFIL,*)
      WRITE(NFIL,FMT='("EXTERNAL POTENTIALS ON ORBITALS"/30("="))')
      DO IPOT=1,NPOT
!       CALL REPORT$I4VAL(NFIL,"POTENTIAL NR:",IPOT,' ')
!       CALL REPORT$CHVAL(NFIL,'ATOM',POT(IPOT)%ATOM)
!       CALL REPORT$CHVAL(NFIL,'TYPE',POT(IPOT)%TYPE)
!       CALL REPORT$R8VAL(NFIL,'VALUE',POT(IPOT)%VALUE,'H')
!       CALL REPORT$I4VAL(NFIL,'SPIN',POT(IPOT)%NSPIN,' ')
        WRITE(NFIL,FMT='("POTENTIAL NR: ",I5)')IPOT
        CALL WRITECH(NFIL,'ATOM',POT(IPOT)%ATOM)
        CALL WRITECH(NFIL,'TYPE',POT(IPOT)%TYPE)
        CALL WRITER8(NFIL,'VALUE',POT(IPOT)%VALUE,'H')
        CALL WRITER8(NFIL,'PWR',POT(IPOT)%PWR,' ')
        CALL WRITER8(NFIL,'RC',POT(IPOT)%RC,'ABOHR')
        CALL WRITEI4(NFIL,'IDIMD ([0=NT],[0=NT,1=NS],[0=NT,1=NX,2=NY,3=NZ])' &
     &                   ,POT(IPOT)%IDIMD-1,' ')
      ENDDO
      RETURN
      CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE WRITECH(NFIL,NAME,VALUE)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        CHARACTER(LEN=*),INTENT(IN) :: VALUE
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,A)')NAME,TRIM(VALUE)
        RETURN
        END SUBROUTINE WRITECH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE WRITER8(NFIL,NAME,VALUE,UNIT)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        REAL(8)          ,INTENT(IN) :: VALUE
        CHARACTER(LEN=*),INTENT(IN) :: UNIT
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,F10.5," ",A)')NAME,VALUE,UNIT
        RETURN
        END SUBROUTINE WRITER8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE WRITEI4(NFIL,NAME,VALUE,UNIT)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        INTEGER(4)       ,INTENT(IN) :: VALUE
        CHARACTER(LEN=*),INTENT(IN) :: UNIT
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,I10," ",A)')NAME,VALUE,UNIT
        RETURN
        END SUBROUTINE WRITEI4
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATH,ETOT)
!     **************************************************************************
!     **  TOTAL ENERGY AND ONE-CENTER HAMILTONIAN CONTRIBUTIONS FROM          **
!     **  EXTERNAL POTENTIAL                                                  **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: ATOM
      INTEGER(4)  ,INTENT(IN)   :: NDIMD
      INTEGER(4)  ,INTENT(IN)   :: LMNX
      COMPLEX(8)  ,INTENT(IN)   :: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(OUT)  :: ETOT
      TYPE(EXTPOT)              :: POT1
      INTEGER(4)                :: LANG
      INTEGER(4)                :: MANG
      INTEGER(4)                :: IAT    ! ATOM INDEX
      INTEGER(4)                :: ISP    ! ATOM TYPE INDEX
      INTEGER(4)                :: NR     ! #(GRID RADIAL POINTS)
      INTEGER(4)                :: LNX    ! #(PARTIAL WAVES (L,N))
      INTEGER(4)   ,ALLOCATABLE :: LOX(:) !(LNX) #(GRID RADIAL POINTS)
      REAL(8)      ,ALLOCATABLE :: AEPHI(:,:) !(NR,LNX) AE PARTIAL WAVES
      REAL(8)      ,ALLOCATABLE :: UONE(:,:)  !(LNX,LNX)
      REAL(8)      ,ALLOCATABLE :: RDEP(:)    !(NR)
      REAL(8)      ,ALLOCATABLE :: AUX(:)    !(NR)
      CHARACTER(32)             :: SPECIES
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: LN1,LN2,LMN1,LMN2,IPOT,IDIMD,L1,L2,I,LM
      INTEGER(4)                :: GID     ! GRID ID
      REAL(8)      ,ALLOCATABLE :: R(:)    ! RADIAL GRID
      LOGICAL(4)                :: TACTIVE
      INTEGER(4)                :: NFILTRACE
!     **************************************************************************
      DATH(:,:,:)=0.D0
      ETOT=0.D0
!
!     ==========================================================================
!     ==  SELECT POTENTIAL                                                    ==
!     ==========================================================================
!     == DETERMINE ID OF THE ATOM TYPE (SPECIES)
      CALL ATOMLIST$INDEX(ATOM,IAT)
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
!
      TACTIVE=.FALSE.
      DO IPOT=1,NPOT
        POT1=POT(IPOT)
        TCHK=(ATOM.EQ.POT1%ATOM.OR.SPECIES.EQ.POT1%ATOM)
        IF(.NOT.TCHK) CYCLE
        TACTIVE=.TRUE.
!
!       ========================================================================
!       ==  UNSCRAMBLE TYPE                                                   ==
!       ========================================================================
!       == LANG=MAIN ANGULAR MOMENTUM  LANG=-1 ALL ANGULAR MOMENTA
!       == MANG="MAGNETIC" QUANTUM NUMBER  MANG=-1 ALL "MAGNETIC" QUANTUM NUMBRS
!       == CAUTION: MANG REFERS TO REAL SPHERICAL HARMONICS
!
        IF(TRIM(POT1%TYPE).EQ.'ALL') THEN
          LANG=-1
!
!       == SELECT COMPLETE ANGULAR MOMENTUM SHELLS =============================
        ELSE IF(TRIM(POT1%TYPE).EQ.'S') THEN
          LANG=0
          MANG=-1
        ELSE IF(TRIM(POT1%TYPE).EQ.'P') THEN
          LANG=1
          MANG=-1
        ELSE IF(TRIM(POT1%TYPE).EQ.'D') THEN
          LANG=2
          MANG=-1
        ELSE IF(TRIM(POT1%TYPE).EQ.'F') THEN
          LANG=3
          MANG=-1
        ELSE   
!         == SELECT INDIVIDUAL REAL SPHERICAL HARMONICS ========================
!         == SPHERICAL$LMBYNAME THROWS AN ERROR IF TYPE IS NOT RECOGNIZED ======
          CALL SPHERICAL$LMBYNAME(TRIM(POT1%TYPE),LM)
          LANG=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
          MANG=LM-LANG**2        
        END IF
!
!       ========================================================================
!       == COLLECT PARTIAL WAVES                                              ==
!       ========================================================================
        CALL ATOMLIST$INDEX(ATOM,IAT)
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETCH('ID',SPECIES)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        ALLOCATE(AEPHI(NR,LNX))
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
        CALL SETUP$UNSELECT()
!
!       ========================================================================
!       == SPECIFY SHAPE OF THE EXTERNAL POTENTIAL                            ==
!       ========================================================================
        ALLOCATE(RDEP(NR))
        ALLOCATE(AUX(NR))
        RDEP(:)=POT1%VALUE*EXP(-(R(:)/POT1%RC)**POT1%PWR)
!
!       ========================================================================
!       == DETERMINE SCALAR PRODUCTS WITH PARTIAL WAVES                       ==
!       == UONE(LN1,LN1)=<AEPHI(LN1)|U(R)|AEPHI(LN2)>                         ==
!       ========================================================================
        ALLOCATE(UONE(LNX,LNX))
        UONE(:,:)=0.D0
        LMN1=0
        DO LN1=1,LNX
          IF(.NOT.(LANG.EQ.LOX(LN1).OR.LANG.EQ.-1)) CYCLE
          DO LN2=1,LNX
            IF(.NOT.(LANG.EQ.LOX(LN2).OR.LANG.EQ.-1)) CYCLE
            IF(LOX(LN2).NE.LOX(LN1)) CYCLE
            IF(LN2.LT.LN1) THEN
              UONE(LN1,LN2)=UONE(LN2,LN1)
              CYCLE
            END IF
            AUX(:)=RDEP(:)*AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,UONE(LN1,LN2))
          ENDDO
        ENDDO
        DEALLOCATE(R)
        DEALLOCATE(RDEP)
        DEALLOCATE(AUX)
        DEALLOCATE(AEPHI)
!
!       ========================================================================
!       ==  CALCULATE ONE-CENTER HAMILTONIAN                                  ==
!       ========================================================================
        IF(POT1%IDIMD.LT.1.OR.POT1%IDIMD.GT.NDIMD) THEN
          CALL ERROR$MSG('VALUE OF POT1%IDIMD NOT RECOGNIZED')
          CALL ERROR$I4VAL('POT1%IDIMD',POT1%IDIMD)
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$STOP('EXTERNAL1CPOT$APPLY')
        END IF
        DO IDIMD=1,NDIMD
          IF(IDIMD.EQ.POT1%IDIMD) THEN
            LMN1=0
            DO LN1=1,LNX
              L1=LOX(LN1)
              LMN2=0
              DO LN2=1,LNX
                L2=LOX(LN2)
                IF(L1.EQ.L2) THEN
                  DO I=1,2*L1+1
                    IF(.NOT.(MANG.EQ.I.OR.MANG.EQ.-1)) CYCLE
                    DATH(LMN1+I,LMN2+I,IDIMD)=DATH(LMN1+I,LMN2+I,IDIMD) &
       &                                     +UONE(LN1,LN2)
                  ENDDO
                END IF
                LMN2=LMN2+2*L2+1
              ENDDO
              LMN1=LMN1+2*L1+1
            ENDDO
          END IF
        ENDDO
        DEALLOCATE(UONE)
        DEALLOCATE(LOX)
      ENDDO
!
!     ==========================================================================
!     ==  RETURN IF NO POTENTIAL FOR THIS ATOM FOUND                          ==
!     ==========================================================================
      IF(.NOT.TACTIVE) RETURN

CALL TRACE$GETL4('TRACEFILE',TCHK)
IF(TCHK) THEN
  CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
  WRITE(NFILTRACE,*)'DENSITY MATRIX IN EXTERNAL1CPOT. ATOM ',POT1%ATOM 
  WRITE(NFILTRACE,*)TRIM(ATOM),TRIM(SPECIES),LNX,LMNX
  DO IDIMD=1,NDIMD
    WRITE(NFILTRACE,*)'IDIMD',IDIMD
    DO LMN1=1,LMNX
      WRITE(NFILTRACE,FMT='(20F10.3)')REAL(DENMAT(LMN1,:,IDIMD))  
    ENDDO
  ENDDO
  WRITE(NFILTRACE,*)'DELTA HAMILTONIAN FROM EXTERNAL1CPOT'
  DO IDIMD=1,NDIMD
    WRITE(NFILTRACE,*)'IDIMD',IDIMD
    DO LMN1=1,LMNX
      WRITE(NFILTRACE,FMT='(20F10.3)')DATH(LMN1,:,IDIMD) 
    ENDDO
  ENDDO
END IF
!
!     ==========================================================================
!     ==  CALCULATE TOTAL ENERGY                                              ==
!     ==========================================================================
      ETOT=0.D0
      DO IDIMD=1,NDIMD
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            ETOT=ETOT+REAL(DENMAT(LMN1,LMN2,IDIMD)*DATH(LMN2,LMN1,IDIMD),KIND=8)
          ENDDO
        ENDDO
      ENDDO
!     PRINT*,'ETOT ',ETOT
      RETURN 
      END
!
!     ..................................................................
      SUBROUTINE GAUSSPOT(RC,L,R,V)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE HARTREE POTENTIAL OF A GENERALIZED           **
!     **  GAUSSIAN  R**L*EXP(-(R/RC)**2)*YLM                          **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RC
      REAL(8)   ,INTENT(IN) :: R
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(OUT):: V
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: FAC,X,SVAR1
!     ******************************************************************
      FAC=4.D0*PI/DBLE(2*L+1)*RC**DBLE(L+2)
      X=R/RC
      CALL DG2NDX(L+1,X,SVAR1)
      V=FAC*(SVAR1*X**DBLE(-L-1)+0.5D0*X**DBLE(L)*EXP(-X**2))
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE DG2NDX(N,X,V)
!     ******************************************************************
!     **                                                              **
!     **  V = INT{DR(0,R)| X**(2*N) *  EXP(-(X)**2) }                 **
!     **                                                              **
!     **  START FROM D/DX{X**(2*N-1)*EXP(-X**2)}                      **
!     **  BUILD RECURSION FOR THE INTEGRAL                            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8),   INTENT(IN) :: X
      REAL(8),   INTENT(OUT):: V
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: X2,GAUSS,FAC1,FAC2
      REAL(8)               :: ERFX
      INTEGER(4)            :: I
!     ******************************************************************
      X2=X*X
      CALL LIB$ERFR8(X,ERFX)
      V=0.5D0*SQRT(PI)*ERFX
      GAUSS=EXP(-X2)
      FAC1=0.5D0
      FAC2=0.5D0*X
      DO I=1,N
        V=FAC1*V - FAC2*GAUSS 
        FAC1=FAC1+1.D0
        FAC2=FAC2*X2
      ENDDO
      RETURN
      END
!!$!
!!$!     .....................................................................
!!$      SUBROUTINE AUGMENTATION_SOFTCORE(ISP,R1,DEX,NR,AERHO1,AEPOT1)
!!$      USE SCHRGL_INTERFACE_MODULE, ONLY : OUTBOUNDARY
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: ISP
!!$      REAL(8)   ,INTENT(IN) :: R1
!!$      REAL(8)   ,INTENT(IN) :: DEX
!!$      INTEGER(4),INTENT(IN) :: NR
!!$      INTEGER(4),INTENT(IN) :: LMRX
!!$      REAL(8)   ,INTENT(IN) :: AERHO1(NR)
!!$      REAL(8)   ,INTENT(IN) :: AEPOT1(NR)
!!$      REAL(8)               :: AECORE0(NR)
!!$      REAL(8)               :: PHI(NR,3)
!!$      REAL(8)   ,PARAMETER  :: TOLB=1.D-6
!!$      INTEGER(4),PARAMETER  :: IBIX=100
!!$      INTEGER(4)            :: NB
!!$      INTEGER(4),ALLOCATABLE:: LOFB(:)  ! MAIN ANGULAR MOMENTUM
!!$      INTEGER(4),ALLOCATABLE:: NNOFB(:) ! #(NODES)
!!$      REAL(8)   ,ALLOCATABLE:: FOFB(:)  ! OCCUPATION
!!$      REAL(8)   ,ALLOCATABLE:: EOFB(:)  ! ENERGY
!!$      REAL(8)   ,ALLOCATABLE:: RCLOFB(:) 
!!$      REAL(8)               :: AEZ
!!$      INTEGER(4)            :: NITER
!!$      INTEGER(4)            :: ITER,IB
!!$!     **********************************************************************
!!$      CALL SETUP$AECORE(ISP,NR,AECORE0)
!!$      CALL SETUP$AEZ(ISP,AEZ)
!!$!     ======================================================================
!!$!     ==  HARDWIRE PARAMETERS FOR IRON                                    ==
!!$!     ======================================================================
!!$      ALLOCATE(LOFB(NB))
!!$      ALLOCATE(FOFB(NB))
!!$      ALLOCATE(EOFB(NB))
!!$      ALLOCATE(RCLOFB(NB))
!!$      NB=5
!!$      LOFB=(/0,0,1,0,1/)
!!$      NNOFB=(/0,1,0,2,1/)
!!$      DO IB=1,NB
!!$        FOFB(IB)=REAL(2*L+1)
!!$        ZEFF=AEZ-SUM(FOFB(1:IB-1))
!!$        RCLOFB(IB)=1.D0/ZEFF/REAL(LOFB(IB)+NNOFB(IB)+1)
!!$        IR=1+NINT(LOG(RCLOFB(IB)/R1)/DEX)
!!$        EOFC(IB)=AEPOT(IR)
!!$      ENDDO
!!$
!!$!     ======================================================================
!!$!     ==  START LOOP                                                      ==
!!$!     ======================================================================
!!$      NITER=2
!!$      DO ITER=1,NITER
!!$!
!!$!       ==============================================================
!!$!       ==  FIND CORE STATES                                        ==
!!$!       ==============================================================
!!$        AECORE(:)=0.D0
!!$        DO IB=1,NLB
!!$          X0=EOFB(IB)
!!$          ISTART=1
!!$          TCONV=.FALSE.
!!$          CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
!!$          DO IBI=1,IBIX
!!$            TFORWARD=.TRUE.
!!$            IR=1+NINT(LOG(RCLOFB(IB)/R1)/DEX)
!!$            CALL SCHROEDER(TFORWARD,.FALSE.,R1,DEX,NR,LOFB(IB),EOFB(IB) &
!!$     &                    ,AEZ,IR,4.D0 &
!!$     &                     ,AEPOT,PHI,GINH,DLG)
!!$            DIN=RCLOFB(IB)*PHI(IR,2)/PHI(IR,1)
!!$            DIN=NNOFB(IB)+0.5D0+ATAN(DIN)/PI            
!!$            TFORWARD=.FALSE.
!!$            CALL SCHROEDER(TFORWARD,.FALSE.,R1,DEX,NR,LOFB(IB),EOFB(IB) &
!!$     &                    ,AEZ,IR,4.D0 &
!!$     &                     ,AEPOT,PHI,GINH,DLG)
!!$            DOUT=RCLOFB(IB)*PHI(IR,2)/PHI(IR,1)
!!$            DIN=0.5D0+ATAN(DOUT)/PI            
!!$            Y0=DIN-DOUT
!!$            TCONV=ABS(DX).LT.TOLB
!!$            IF(TCONV) EXIT
!!$            CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
!!$            EOFB(IB)=X0
!!$          ENDDO
!!$          IF(.NOT.TCONV) THEN
!!$            CALL ERROR$MSG('LOOP NOT CONVERGED')
!!$            CALL ERROR$STOP('AUGMENTATION_SOFTCORE')
!!$          END IF
!!$          AECORE(:)=AECORE(:)+FB(IB)*PHI(:)**2
!!$        ENDDO
!!$!
!!$!       ==============================================================
!!$!       ==  CONSTRUCT DENSITY                                       ==
!!$!       ==============================================================
!!$        RHO(:)=RHO0(:)+AECORE(:)-AECORE0(:)
!!$!
!!$!       ==============================================================
!!$!       ==  CONSTRUCT POTENTIAL                                     ==
!!$!       ==============================================================
!!$
!!$
!!$
!!$      ENDDO
!!$      RETURN
!!$      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADDVQLM(GID,NR,LMRX,VQLM,AEPOT,PSPOT)
!     ******************************************************************
!     **                                                             **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: GID
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(IN)   :: VQLM(LMRX)
      REAL(8)   ,INTENT(INOUT):: AEPOT(NR,LMRX)
      REAL(8)   ,INTENT(INOUT):: PSPOT(NR,LMRX)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: SVAR
      INTEGER(4)              :: LM,L
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        SVAR=VQLM(LM)
        AEPOT(:,LM)=AEPOT(:,LM)+SVAR*R(:)**L
        PSPOT(:,LM)=PSPOT(:,LM)+SVAR*R(:)**L
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADDBACKGROUND(GID,NR,RHOB &
     &                        ,RHO,EB,POT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ONE-CENTER CONTRIBUTION OF THE COMPENSATING  **
!     **  CHARGE BACKGROUND TO TOTAL ENERGY AND POTENTIAL             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: GID
      INTEGER(4),INTENT(IN)   :: NR
      REAL(8)   ,INTENT(IN)   :: RHOB
      REAL(8)   ,INTENT(IN)   :: RHO(NR)
      REAL(8)   ,INTENT(OUT)  :: EB
      REAL(8)   ,INTENT(INOUT):: POT(NR)
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)                 :: SVAR
      REAL(8)                 :: WORK(NR)
      REAL(8)                 :: R(NR)
!     ******************************************************************
      IF(RHOB.EQ.0.D0) THEN
        EB=0.D0
        RETURN
      END IF
      CALL RADIAL$R(GID,NR,R)
      SVAR=-2.D0*PI*RHOB/3.D0*SQRT(4.D0*PI)
      WORK(:)=SVAR*RHO(:)*R(:)**4
      POT(:)=POT(:)+SVAR*R(:)**2 ! POTENTIAL OF THE BACKGROUND
      CALL RADIAL$INTEGRAL(GID,NR,WORK,EB)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_HARTREE(IAT,GID,NR,AEZ,RHOCOR &
     &            ,VADD,RCSM,LMRX,AERHO,PSRHO,QLM,AEPOT,PSPOT,VQLM &
     &            ,AEEHARTREE,PSEHARTREE)
!     ******************************************************************
!     **                                                              **
!     ** CALCULATES HARTREE ENERGY                                    **
!     ** ADDS COMPENSATION DENSITY TO THE PSEUDO DENSITY              **
!     **                                                              **
!     **  1) AEE=AERHO*AEZ/R                                          **
!     **     PSE=PSRHO*VADD                                           **
!     **  2) PSRHO=PSRHO+GAUSS(R)*QLM                                 **
!     **  3) AEE=AEE + 0.5*AERHO(R)*AERHO(R')/|R-R'|                  **
!     **     PSE=PSE + 0.5*PSRHO(R)*PSRHO(R')/|R-R'|                  **
!     **     VQLM=VQLM - GAUSS(R)*PSRHO(R')/|R-R'|                    **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: IAT
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      REAL(8)    ,INTENT(IN) :: AEZ
      REAL(8)    ,INTENT(IN) :: RHOCOR(NR)
      REAL(8)    ,INTENT(IN) :: VADD(NR)
      REAL(8)    ,INTENT(IN) :: RCSM
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)    ,INTENT(INOUT) :: PSRHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(OUT):: VQLM(LMRX)
      REAL(8)    ,INTENT(OUT):: AEPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: PSPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: AEEHARTREE
      REAL(8)    ,INTENT(OUT):: PSEHARTREE
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX1(NR)
      REAL(8)                :: AEDMU(NR)
      REAL(8)                :: PSDMU(NR)
      REAL(8)                :: AEE(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      REAL(8)                :: GAUSSIAN(NR)
      REAL(8)                :: AEH,PSH
      INTEGER(4)             :: LM,IR
      INTEGER(4)             :: L
!     == ARRAYS NEEDED DETAILED REPORT
      LOGICAL(4) ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4)             :: NFILO
!     ******************************************************************
      CALL ERROR$MSG('THIS ROUTINE IS APPARENTLY NOT USED')
      CALL ERROR$MSG('IT IS MARKED FOR DELETION')
      CALL ERROR$STOP('AUGMENTATION_HARTREE')

      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  INITIALIZE ENERGIES AND POTENTIALS TO  ZERO                 ==
!     ==================================================================
      AEEHARTREE=0.D0
      PSEHARTREE=0.D0
      VQLM(:)=0.D0
      AEPOT(:,:)=0.D0
      PSPOT(:,:)=0.D0
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC VALENCE-CORE INTERACTION           ==
!     ==   AND VALENCE-NUCLEUS INTERACTION                            ==
!     ==================================================================
      AEH=0.D0
      PSH=0.D0
      AEE(:)=0.D0
      PSE(:)=0.D0
      AEDMU(:)=0.D0
      CALL RADIAL$POISSON(GID,NR,0,RHOCOR,AEDMU)
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX1)
      AEDMU(:)  = AEDMU(:)+AUX1(:)
      AEE(:)    = AEDMU(:)* AERHO(:,1) *R(:)**2
      AEPOT(:,1)= AEPOT(:,1)+AEDMU(:)
      PSE(:)    = VADD(:) * PSRHO(:,1) *R(:)**2
      PSPOT(:,1)= PSPOT(:,1)+VADD(:)
      CALL RADIAL$INTEGRAL(GID,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(GID,NR,PSE,PSH)
!PRINT*,'AEH C-V INT',AEH
      AEEHARTREE=AEEHARTREE+AEH
      PSEHARTREE=PSEHARTREE+PSH
!
!     ==================================================================
!     ==   ADD COMPENSATION CHARGE DENSITY                            ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL GAUSSN(L,ALPHA,CL)
        SVAR=QLM(LM)*CL
        PSRHO(:,LM)=PSRHO(:,LM)+SVAR*EXP(-ALPHA*R(:)**2)*R(:)**L
      ENDDO
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      AEE(:)=0.D0
      PSE(:)=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL RADIAL$POISSON(GID,NR,L,AERHO(1,LM),AEDMU)
        CALL RADIAL$POISSON(GID,NR,L,PSRHO(1,LM),PSDMU)
        ALPHA=1.D0/RCSM**2
        CALL GAUSSN(L,ALPHA,CL)
!       == ADD COMPENSATION DENSITY AND ITS POTENTIAL ================
        GAUSSIAN(:)=CL*EXP(-ALPHA*R(:)**2)*R(:)**L
!
!       == CALCULATE TOTAL ENERGY AND POTENTIAL ======================
        AEE(:)=AEE(:)+0.5D0*AEDMU(:)*AERHO(:,LM)*R(:)**2
        PSE(:)=PSE(:)+0.5D0*PSDMU(:)*PSRHO(:,LM)*R(:)**2
        AEPOT(:,LM)=AEPOT(:,LM)+AEDMU(:)
        PSPOT(:,LM)=PSPOT(:,LM)+PSDMU(:)
        AUX1(:)=PSDMU(:)*GAUSSIAN*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX1,SVAR)
        VQLM(LM)=VQLM(LM)-SVAR
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(GID,NR,PSE,PSH)
      AEEHARTREE=AEEHARTREE+AEH
      PSEHARTREE=PSEHARTREE+PSH
!PRINT*,'AEH V-V INT',AEH
!
!     ==================================================================
!     ==  PRINTOUT FOR TESTING                                        ==
!     ==================================================================
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'AE POTENTIAL FOR ATOM ',IAT
        DO IR=1,NR,49
          WRITE(NFILO,FMT='(9F10.5)')(AEPOT(IR,LM),LM=1,LMRX)
        ENDDO
        WRITE(NFILO,*)'PS POTENTIAL FOR ATOM ',IAT
        DO IR=1,NR,49
          WRITE(NFILO,FMT='(9F10.5)')(PSPOT(IR,LM),LM=1,LMRX)
        ENDDO
        WRITE(NFILO,*)'PS DENSITY FOR ATOM ',IAT
        DO IR=1,NR,49
          WRITE(NFILO,FMT='(9F10.5)')(PSRHO(IR,LM),LM=1,LMRX)
        ENDDO
      ENDIF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DA1,DA)
!     **                                                              **
!     **  CONVERTS THE MATRIX ELEMENTS OF A SPHERICALLY SYMMETRIC     **
!     **  OPERATOR INTO THE FORM OF THE DENSITY MATRIX                **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LMNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: DA1(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: DA(LMNX,LMNX,NDIMD)
      INTEGER(4)            :: LMN1,LMN2,LN1,LN2,L1,L2,IM
!     ******************************************************************      
      DA(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            DO IM=1,2*L1+1
              DA(LMN1+IM,LMN2+IM,1)=DA1(LN1,LN2)
            ENDDO
          END IF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION_WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI(IR,:))).GT.1.D+3) EXIT
        WRITE(100,FMT='(F15.10,2X,20(E25.10,2X))')R(IR),PHI(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END



