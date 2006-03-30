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
!**    EXTERNAL1CPOT, LDAPLUSU                                        **
!**    CLEBSCH, RADIAL, DFT, GAUSSN                                   **
!**    SELFTEST, TRACE, FILEHANDLER, ERROR                            **
!**                                                                   **
!***********************************************************************
!***********************************************************************
MODULE AUGMENTATION_MODULE
INTEGER(4)  ,PARAMETER :: NE=8
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
      REAL(8)                 :: SVAR
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
      USE AUGMENTATION_MODULE
      USE MPE_MODULE
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
      SVAR=VAL(7)+VAL(1)-VAL(2)+VAL(3)-VAL(4)+VAL(5)-VAL(6)+VAL(8)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',SVAR)
      CALL ENERGYLIST$ADD('AE  KINETIC',VAL(7))
      CALL ENERGYLIST$ADD('AE  EXCHANGE-CORRELATION',VAL(1)-VAL(2)+VAL(8))
      CALL ENERGYLIST$ADD('AE  ELECTROSTATIC',VAL(3)-VAL(4))
      CALL ENERGYLIST$ADD('BACKGROUND',VAL(5)-VAL(6))
      VAL(:)=0.D0
      RETURN
      END
!
!     .....................................................MOMNTS.......
      SUBROUTINE AUGMENTATION$MOMENTS(ISP,LMNX,DENMAT,LMRX,QLM)
!     ******************************************************************
!     **                                                             **
!     **  CALCULATES THE ELECTROSTATIC MULTIPOLE MOMENTS FOR THE     **
!     **  COMPENSTATION CHARGE DENSITY                               **
!     **                                                             **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)**
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      INTEGER(4),INTENT(IN)  :: ISP
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: LMRX
      REAL(8)   ,INTENT(IN)  :: DENMAT(LMNX,LMNX)
      REAL(8)   ,INTENT(OUT) :: QLM(LMRX)
      REAL(8)   ,ALLOCATABLE :: AERHO(:,:)   !(NR,LMRX)
      REAL(8)   ,ALLOCATABLE :: PSRHO(:,:)   !(NR,LMRX)
      REAL(8)   ,ALLOCATABLE :: AECORE(:)    !(NR)
      REAL(8)   ,ALLOCATABLE :: PSCORE(:)    !(NR)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)   !(NR,LNX)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)   !(NR,LNX)
      INTEGER(4),ALLOCATABLE :: LOX(:)       !(LN)
      integer(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)                :: AEZ
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMN1,LMN2,IDIM
!     ******************************************************************
                          CALL TRACE$PUSH('AUGMENTATION$MOMENTS')
!
!     == PRINT DENSITY MATRIX FOR TEST ========== ======================
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
!       WRITE(*,FMT='("TOTAL DENSITY MATRIX")')
        DO LMN1=1,LMNX
!         WRITE(*,FMT='(9F10.6)')(DENMAT(LMN1,LMN2),LMN2=1,LMNX)
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
      CALL SETUP$LNX(ISP,LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$LOFLN(ISP,LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
      CALL SETUP$PSPARTIALWAVES(ISP,NR,LNX,PSPHI)
!
!     ==================================================================
!     ==  RECEIVE CORE DENSITY FROM SETUPS OBJECT                     ==
!     ==  RECEIVE ATOMIC NUMBER FROM ATOMLIST                         ==
!     ==================================================================
      ALLOCATE(PSCORE(NR))
      ALLOCATE(AECORE(NR))
      CALL SETUP$PSCORE(ISP,NR,PSCORE)
      CALL SETUP$AECORE(ISP,NR,AECORE)
      CALL SETUP$AEZ(ISP,AEZ)
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
      CALL AUGMENTATION_QLM(gid,NR,LMRX,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
      DEALLOCATE(AERHO)
      DEALLOCATE(PSRHO)
      DEALLOCATE(PSCORE)
      DEALLOCATE(AECORE)
      IF(TPR) THEN
        WRITE(*,FMT='(A63)')'MULTIPOLE MOMENTS'
        WRITE(*,FMT='(9E12.5)')QLM(:)
      END IF
                   CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION$SPHERE(ISP,IAT,LMNX,NDIMD,DENMAT,DENMATI &
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
      REAL(8)   ,INTENT(INOUT):: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(INOUT):: DENMATI(LMNX,LMNX,NDIMD)
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(INOUT):: VQLM(LMRX)
      REAL(8)   ,INTENT(IN)   :: RHOB ! neutralizing background density
!                               rhob may be set to zero by isolate object
      REAL(8)   ,INTENT(OUT)  :: POTB ! NEG. AV. EL. AUGM. POT.
      REAL(8)   ,INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT)  :: DO(LMNX,LMNX,NDIMD)
      integer(4)              :: gid   ! grid id
      REAL(8)   ,ALLOCATABLE  :: r(:)  ! radial grid
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
      REAL(8)   ,ALLOCATABLE  :: DTKIN(:,:)
      REAL(8)   ,ALLOCATABLE  :: AERHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DATP(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DWORK1(:)
      INTEGER(4)              :: IDIM,LMR,IR
      INTEGER(4)              :: LM,LMN1,LMN2,LMN11,LMN21,LN1,LN2,L1,L2,IM,ISPIN
      INTEGER(4)              :: NSPIN
      INTEGER(4)              :: NFILO
      LOGICAL(4),PARAMETER    :: TPR=.false.
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=1
      LOGICAL(4)              :: TBACK,TSPIN
      REAL(8)                 :: DETOT,PSEHARTREE,AEEHARTREE,COREEXC
      REAL(8)                 :: EKINNL,ENL,AEEXC,PSEXC,HAMUP,HAMDWN
      CHARACTER(32)           :: ATOM
      REAL(8)                 :: VQLM1(LMRX)
      REAL(8)                 :: QLM(LMRX)
      REAL(8)                 :: PI
      LOGICAL(4),PARAMETER    :: TSOFTCORE=.false.
      real(8)   ,allocatable  :: aehpot(:,:)
      real(8)   ,allocatable  :: aexcpot(:,:,:)
      real(8)   ,allocatable  :: pshpot(:,:)
      real(8)   ,allocatable  :: psxcpot(:,:,:)
      real(8)   ,allocatable  :: aetotpot(:,:,:)
      real(8)   ,allocatable  :: pstotpot(:,:,:)
      real(8)   ,allocatable  :: rho(:,:,:)
!     ******************************************************************
                            CALL TRACE$PUSH('AUGMENTATION$SPHERE')
      PI=4.D0*DATAN(1.D0)
print*,'new sphere'
!
!     ==================================================================
!     ==  COLLECT ATOM-TYPE SPECIFIC INFORMATION FROM SETUP OBJECT    ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      allocate(r(nr))
      call radial$r(gid,nr,r)
      CALL SETUP$AEZ(ISP,AEZ)
      CALL SETUP$LNX(ISP,LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$LOFLN(ISP,LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUP$PSPARTIALWAVES(ISP,NR,LNX,PSPHI)
      ALLOCATE(AECORE(NR))
      CALL SETUP$AECORE(ISP,NR,AECORE)
      ALLOCATE(PSCORE(NR))
      CALL SETUP$PSCORE(ISP,NR,PSCORE)
      CALL SETUP$RCSM(ISP,RCSM)
      ALLOCATE(VADD(NR))
      CALL SETUP$VBAR(ISP,NR,VADD)
      ALLOCATE(DOVER(LNX,LNX))
      CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
      ALLOCATE(DTKIN(LNX,LNX))
      CALL SETUP$1CKINETIC(ISP,LNX,DTKIN)
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
            DENMAT(:,:,IDIM)=0.5D0*(DENMAT(:,:,IDIM)+TRANSPOSE(DENMAT(:,:,IDIM)))
          ENDDO
        END IF
      ENDIF
!     == PRINT DENSITY MATRIX FOR TEST
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,FMT='("DENMAT FOR IDIM= ",I2)') IDIM
          DO LMN1=1,LMNX
            WRITE(NFILO,FMT='(9F15.6)')DENMAT(LMN1,:,IDIM)
          ENDDO
        ENDDO
      END IF
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
      rho(:,:,1)=psrho(:,:,1)
      RHO(:,1,1)=RHO(:,1,1)+PSCORE(:)   
      CALL AUGMENTATION_ADJUSTVQLM(GID,NR,LMRX,RHO,RCSM,QLM,VQLM1)
      vqlm1(:)=vqlm1(:)+vqlm(:)
      DEALLOCATE(RHO)
!     
!     ================================================================
!     ==  NEW SOFT CORE                                             ==
!     ================================================================
      IF(TSOFTCORE) THEN
        CALL AUGMENTATION_NEWSOFTCORE(GID,NR,AEZ,AERHO(:,1,1),AECORE &
     &                               ,lmrx,lmnx,denmat(:,:,1),vqlm1,rhob)
      END IF
!
!     =================================================================
!     == HARTREE ENERGY AND POTENTIAL                                ==
!     =================================================================
      ALLOCATE(PSHPOT(NR,LMRX))
      ALLOCATE(AEHPOT(NR,LMRX))
      CALL AUGMENTATION_PSHARTREE(GID,NR,LMRX,PSCORE,PSRHO(:,:,1) &
     &                 ,VADD,RCSM,QLM,VQLM1,RHOB,PSHPOT,PSEHARTREE)
      CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AECORE,AERHO(:,:,1) &
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
!     ==   report energies                                          ==
!     ================================================================
      CALL AUGMENTATION_ADD('AE1 BACKGROUND',0.D0)
      CALL AUGMENTATION_ADD('PS1 BACKGROUND',0.D0)
      CALL AUGMENTATION_ADD('AE1 EXCHANGE-CORRELATION',AEEXC)
      CALL AUGMENTATION_ADD('PS1 EXCHANGE-CORRELATION',PSEXC)
      CALL AUGMENTATION_ADD('AE1 ELECTROSTATIC',AEEHARTREE)
      CALL AUGMENTATION_ADD('PS1 ELECTROSTATIC',PSEHARTREE)
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
!     ================================================================
!     ==  SEND DENSITIES TO HYPERFINE-PARAMETER OBJECT              ==
!     ================================================================
      CALL AUGMENTATION_FEEDHYPERFINE(GID,NR,LMRX,NDIMD &
     &         ,IAT,AERHO,AECORE,PSRHO,PSCORE,AEHPOT,PSHPOT)
!     
!     ================================================================
!     ==  SEND DENSITIES AND POTENTIALS TO GRAPHICS OBJECT          ==
!     ================================================================
      CALL GRAPHICS$SET1CPOT('AE',IAT,gid,NR,NR,LMRX,AEhPOT)
      CALL GRAPHICS$SET1CPOT('PS',IAT,gid,NR,NR,LMRX,PShpOT)
!     
!     ================================================================
!     ==  evaluate core shifts                                      ==
!     ================================================================
      CALL CORE_CORESHIFTS(IAT,ISP,gid,NR,LMRX,AEtotPOT)
!
      DEALLOCATE(AEHPOT)
      DEALLOCATE(PSHPOT)
      DEALLOCATE(AEXCPOT)
      DEALLOCATE(PSXCPOT)
!
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO ISPIN=1,NSPIN
          WRITE(NFILO,*)'AE POTENTIAL FOR SPIN ',ISPIN
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(AEtotPOT(IR,LM,ISPIN),LM=1,LMRX)
          ENDDO
        ENDDO
        DO ISPIN=1,NSPIN
          WRITE(NFILO,*)'PS POTENTIAL FOR ATOM ',IAT,ISPIN
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(PStotPOT(IR,LM,ISPIN),LM=1,LMRX)
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
      EKINNL=0.D0
      DATH(:,:,:)=0.D0
      DO(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            DO IM=1,2*L1+1
              LMN11=LMN1+IM
              LMN21=LMN2+IM
              EKINNL=EKINNL+DENMAT(LMN11,LMN21,1)*DTKIN(LN1,LN2)
              DATH(LMN11,LMN21,1)=DTKIN(LN1,LN2)
              DO(LMN11,LMN21,1)  =DOVER(LN1,LN2)
            ENDDO
          ENDIF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO
      IF(TTEST) THEN
        IF(ITEST.EQ.2) THEN
          CALL SELFTEST$END('SPHERE-EKI',LMNX*LMNX*NDIMD,DATH,EKINNL,TBACK)
          IF(TBACK) GOTO 1001
          CALL ERROR$MSG('STOP AFTER SELFTEST')
          CALL ERROR$STOP('SPHERE')
        END IF
      END IF
!     
      CALL AUGMENTATION_ADD('AE1-PS1 KINETIC',EKINNL)
!     
!     ================================================================
!     ==   ADD POTENTIAL ENERGY TO THE ONE-CENTER HAMILTONIAN       ==
!     ==   DATH =DATH + <AEPHI|AEPOT|AEPHI>-<PSPHI|PSPOT|PSPHI>     ==
!     ================================================================
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL AUGMENTATION_EXPECT(gid,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                        ,AEtotpot,PStotPOT,AEPHI,PSPHI,DATP)
      DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
      DEALLOCATE(DATP)
!     WRITE(TESTSTRING,FMT='("R8DATH",I2,12(" "))')IAT
!     CALL STOREIT(TESTSTRING,8*LMNX*LMNX*NSPIN,DATH)
!
!     ================================================================
!     ==  LDA + U                                                   ==
!     ================================================================
!     ALLOCATE(DATP(LMNX,LMNX,NDIMD))
!     DETOT=0.D0
!     CALL LDAPLUSU(NRX,LNX,LMNX,NSPIN,LOX,LNX &
!     &             ,GID,NR,AEZ,AEPHI,DENMAT,DETOT,DATP)
!     DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
!     DEALLOCATE(DATP)
!     CALL AUGMENTATION_ADD('LDA+U EXCHANGE',DETOT)
!
!     ================================================================
!     ==  APPLY EXTERNAL POTENTIAL                                  ==
!     ================================================================
      CALL ATOMLIST$GETCH('NAME',IAT,ATOM)
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATP,DETOT)
      DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
      DEALLOCATE(DATP)
      CALL AUGMENTATION_ADD('LDA+U EXCHANGE',DETOT)
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
            WRITE(NFILO,FMT='(9F15.6)')DATH(LMN1,:,IDIM)
          ENDDO
        ENDDO
STOP
      END IF
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(AEtotPOT)
      DEALLOCATE(PStotPOT)
      DEALLOCATE(AERHO)
      DEALLOCATE(PSRHO)
      DEALLOCATE(LOX)
      DEALLOCATE(AEPHI)
      DEALLOCATE(PSPHI)
      DEALLOCATE(AECORE)
      DEALLOCATE(PSCORE)
      DEALLOCATE(VADD)
      DEALLOCATE(DOVER)
      DEALLOCATE(DTKIN)
      DEALLOCATE(r)
                               CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_NEWSOFTCORE(GID,NR,AEZ,AERHO,AECORE &
     &                                   ,LMRX,LMNX,DENMAT,VQLM,RHOB)        
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AERHO(NR)
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: DENMAT(LMNX,LMNX)
      REAL(8)   ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)   ,INTENT(IN) :: RHOB
      REAL(8)               :: AEHPOT(NR,LMRX)
      REAL(8)               :: AEXCPOT(NR,LMRX)
      REAL(8)               :: POTNS(NR,LMRX)
      REAL(8)               :: RHONS(NR,LMRX)
      REAL(8)               :: PI,Y0,C0LL
      REAL(8)               :: R(NR)
      LOGICAL(4)            :: CONVG
      REAL(8)               :: RHO(NR)
      REAL(8)               :: POT(NR)
      REAL(8)               :: AEPOT(NR)
      REAL(8)               :: POTIN(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: EH,EXC,DECORE,AEEHARTREE
      REAL(8)   ,ALLOCATABLE:: FOFI(:)      
      REAL(8)   ,ALLOCATABLE:: EOFI(:)      
      INTEGER(4),ALLOCATABLE:: LOFI(:)      
      INTEGER(4),ALLOCATABLE:: SOFI(:)      
      INTEGER(4),ALLOCATABLE:: NNOFI(:)      
      INTEGER(4)            :: LNX
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)   ,ALLOCATABLE:: UPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: TUPHI(:,:)
      REAL(8)               :: XMAX
      REAL(8)               :: SVAR
      INTEGER(4)            :: ITER,I,J,IB
      INTEGER(4),PARAMETER  :: NITER=2000
      REAL(8)   ,PARAMETER  :: TOL=1.D-4
      INTEGER(4)            :: NDIMD=1
      INTEGER(4)            :: NC
      REAL(8)               :: ECORE,ECOREAT
      REAL(8)   ,ALLOCATABLE:: PHIC(:,:)
      REAL(8)   ,ALLOCATABLE:: TPHIC(:,:)
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: TAEPHI(:,:)
      REAL(8)               :: RHOOLD
      INTEGER(4)            :: LN,IC
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ================================================================
!     == COLLECT DATA FROM SETUP OBJECT                             ==
!     ================================================================
      CALL SETUP$GETI4('NC',NC)
      ALLOCATE(EOFI(NC))
      ALLOCATE(LOFI(NC))
      ALLOCATE(SOFI(NC))
      ALLOCATE(FOFI(NC))
      ALLOCATE(NNOFI(NC))
      ALLOCATE(PHIC(NR,NC))
      ALLOCATE(TPHIC(NR,NC))
      CALL SETUP$GETR8A('FOFC',NC,FOFI)
      CALL SETUP$GETR8A('EOFC',NC,EOFI)
      CALL SETUP$GETI4A('LOFC',NC,LOFI)
      SOFI(:)=0.D0
      DO I=1,NC
        NNOFI(I)=0
        DO J=1,I-1
          IF(LOFI(J).EQ.LOFI(I))NNOFI(I)=NNOFI(I)+1
        ENDDO
      ENDDO
      CALL SETUP$GETR8A('ATOMICAEPOT',NR,AEPOT)
!
!     ================================================================
!     == CALCULATE ATOMIC CORE ENERGY                               ==
!     ================================================================
      POT=AEPOT
      CALL SF_AERHO(GID,NR,NC,LOFI,SOFI,FOFI,NNOFI,POT,RHO,EOFI,PHIC,TPHIC)
      CALL SF_SPHERICALCOREETOT(GID,NR,NC,EOFI,FOFI,POT,RHO,ECOREAT,AUX)
PRINT*,'===============ATOMIC CORE========================='
DO IB=1,NC
  PRINT*,'L,E ',LOFI(IB),EOFI(IB)
ENDDO
!
!     =====================================================================
!     == SCF LOOP FOR SPHERICAL POTENTIAL                                ==
!     =====================================================================
      POT=AEPOT ! STARTING POTENTIAL IS THE ATOMIC POTENTIAL 
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR,10,1.D0)
      DO ITER=1,NITER
        CALL SF_AERHO(GID,NR,NC,LOFI,SOFI,FOFI,NNOFI,POT,RHO,EOFI,PHIC,TPHIC)
!
!       ================================================================
!       ==  EXIT IF CONVERGED                                         ==
!       ================================================================
        IF(CONVG) THEN
          CALL BROYDEN$CLEAR
          EXIT
        END IF

!       ====================================================================
!       == CALCULATE OUTPUT POTENTIAL                                     ==
!       ====================================================================
        POTIN=POT
        CALL SF_SPHERICALCOREETOT(GID,NR,NC,EOFI,FOFI,AUX,RHO+AERHO,SVAR,POT)
!       == DO NOT FORGET THE EXTERNAL POTENTIAL AND THE BACKGROUND !!
        SVAR=AEPOT(NR)-POT(NR)
        POT=POT+SVAR
!
!       ================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD     ==
!       ================================================================
        XMAX=MAXVAL(ABS(POT-POTIN)) 
        CALL BROYDEN$STEP(NR,POTIN,(POT-POTIN))
        POT=POTIN
        CONVG=(XMAX.LT.TOL)
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$STOP('AUGMETATION_NEWSOFTCORE')
      END IF
!
!     ================================================================
!     ==  NOW CALCULATE TOTAL ENERGY IN SPHERICAL POTENTIAL         ==
!     ================================================================
      CALL SF_SPHERICALCOREETOT(GID,NR,NC,EOFI,FOFI,POT,RHO,ECORE,AUX)
      DECORE=ECORE-ECOREAT   ! CORE RELAXATION ENERGY
!
PRINT*,'===============SCF CORE========================='
DO IB=1,NC
  PRINT*,'L,E ',LOFI(IB),EOFI(IB)
ENDDO
PRINT*,'SCF CORE ENERGY ',ECORE
PRINT*,'CORE RELAXATION ENERGY ENERGY ',DECORE
!
!     ===============================================================
!     == ADJUST PARTIAL WAVES AND CALCULATE DENSITY                ==
!     ===============================================================
!     == COLLECT NODELESS PARTIAL WAVES ============================
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(UPHI(NR,LNX))
      ALLOCATE(TUPHI(NR,LNX))
      CALL SETUP$GETR8A('NDLSPHI',NR*LNX,UPHI)
      CALL SETUP$GETR8A('NDLSTPHI',NR*LNX,TUPHI)
!
!     ===============================================================
!     == ORTHOGONALIZE NODELESS WAVE FUNCTIONS TO THE CORE STATES  ==
!     ===============================================================
!     ==  CORE STATES ARE NOT ORTHONORMAL BECAUSE OF RELATIVISTIC 
!     ==  EFFECTS. THEY ARE NOT TREATED PROPERLY....
      ALLOCATE(AEPHI(NR,LNX))
      ALLOCATE(TAEPHI(NR,LNX))
      AEPHI=UPHI
      TAEPHI=TUPHI
      DO LN=1,LNX
        DO IC=1,NC
          IF(LOX(LN).NE.LOFI(IC)) CYCLE
          AUX(:)=UPHI(:,LN)*PHIC(:,IC)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          UPHI(:,LN)=UPHI(:,LN)-PHIC(:,IC)*SVAR
          TUPHI(:,LN)=TUPHI(:,LN)-TPHIC(:,IC)*SVAR
        ENDDO
      ENDDO  
      CALL AUGMENTATION_RHO(NR,LNX,LOX,AEPHI,LMNX,DENMAT,LMRX,RHONS)
      RHONS(:,1)=RHONS(:,1)+RHO(:)
      DEALLOCATE(AEPHI)
      DEALLOCATE(TAEPHI)
      AUX(:)=0.D0  ! CORE IS SET TO ZERO, DENSITY CONTAINS CORE
      CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AUX,RHONS &
     &                            ,VQLM,RHOB,AEHPOT,AEEHARTREE)
!     == SPIN SWITCHED OFF      
      CALL AUGMENTATION_XC(GID,NR,LMRX,1,AUX,RHONS,AEXCPOT)
      POTNS=AEHPOT+AEXCPOT
!
!     ===============================================================
!     == SELF-CONSISTENT NONSPHERICAL CORE/VALENCE                 ==
!     ===============================================================
      DO ITER=1,100
        CALL SF_DO(GID,NR,LNX,LOX,UPHI,LMNX,DENMAT,LMRX,POTNS &
     &                ,NC,LOFI,PHIC,TPHIC,2.D0,RHONS)
        AUX(:)=0.D0  ! CORE IS SET TO ZERO, DENSITY CONTAINS CORE
        CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AUX,RHONS &
     &                              ,VQLM,RHOB,AEHPOT,AEEHARTREE)
!       == SPIN SWITCHED OFF      
        CALL AUGMENTATION_XC(GID,NR,LMRX,1,AUX,RHONS,AEXCPOT)
        POTNS=AEHPOT+AEXCPOT
      ENDDO
PRINT*,'MARKE AFTER SF_DO'
      DEALLOCATE(LOX)
      DEALLOCATE(UPHI)
      DEALLOCATE(TUPHI)
      DEALLOCATE(EOFI)
      DEALLOCATE(LOFI)
      DEALLOCATE(SOFI)
      DEALLOCATE(FOFI)
      DEALLOCATE(NNOFI)
      DEALLOCATE(PHIC)
      DEALLOCATE(TPHIC)
PRINT*,'DONE'
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE SF_SPHERICALCOREETOT(GID,NR,NB,E,F,POTIN,RHO,ETOT,POTOUT)
!     **                                                                **
!     ** CALCULATES ENERGY OF A SPHERICAL-NON-SPIN-POLARIZED ATOM       **
!     ** FOR GIVEN ENERGY EIGENVALUES, INPUT POTENTIAL AND DENSITY      **
!     **                                                                **
!     ** remark: uses setup object, which must be properly selected     **
!     **                                                                **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NB
      REAL(8)   ,INTENT(IN) :: E(NB)
      REAL(8)   ,INTENT(IN) :: F(NB)
      REAL(8)   ,INTENT(IN) :: POTIN(NR)
      REAL(8)   ,INTENT(IN) :: RHO(NR)
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: POTOUT(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: potnuc(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: Exc,svar
!     *****************************************************************
      CALL AUGMENTATION_XC(GID,NR,1,1,RHO,EXC,POTOUT)
      potout(1)=potout(2)
      CALL SETUP$GETR8A('NUCPOT',NR,POTNUC)
      POTOUT(:)  = POTOUT(:)+POTNUC(:)
      CALL RADIAL$POISSON(GID,NR,0,RHO,AUX)
      POTOUT(:)  = POTOUT(:)+AUX(:)
      CALL RADIAL$R(GID,NR,R)
      AUX(:)=(-POTIN+POTNUC(:)+0.5D0*AUX(:))*RHO(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      ETOT=SUM(E*F)+EXC+SVAR
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE sf_AErho(GID,NR,NB,LOFI,SO,F,NN,POT,RHO,E,phi,tphi)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: F(NB)     !OCCUPATION
      REAL(8)    ,INTENT(IN)     :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: RHO(NR)   !DENSITY
      REAL(8)    ,INTENT(INOUT)  :: E(NB)     !ONE-PARTICLE ENERGIES
      REAL(8)    ,intent(out)    :: PHI(NR,nb)
      REAL(8)    ,intent(out)    :: tPHI(NR,nb)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)                    :: AUX(NR)
      INTEGER(4),parameter       :: niter=100
      real(8)   ,parameter       :: tol=1.d-4
      INTEGER(4)                 :: IB,IR,i
      REAL(8)                    :: SVAR
      REAL(8)                    :: C0LL,PI,Y0
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RHO(:)=0.D0
      DO IB=1,NB
         do i=1,niter
           svar=e(ib)
           CALL RELATIVISTICCORRECTION(GID,NR,POT,E(IB),DREL)
           aux(:)=0.d0
           CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),DREL,aux,NN(IB),POT &
     &                  ,E(IB),PHI(:,ib))
           if(abs(e(ib)-svar).lt.tol) exit
         enddo
         AUX(:)=(R(:)*PHI(:,ib))**2
         CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
         PHI(:,ib)=PHI(:,ib)/SQRT(SVAR)
         tphi(:,ib)=(e(ib)-pot(:)*y0)*phi(:,ib)
         RHO(:)=RHO(:)+F(IB)*C0LL*PHI(:,ib)**2
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE RELATIVISTICCORRECTION(GID,NR,POT,E,DREL)
!     **                                                                  **
!     **  DREL IS A MEASURE OF THE RELATIVISTIC CORRECTIONS               **
!     **     D:=1/MREL-1/M0                                               **
!     **  WHERE MREL IS THE RELATIVISTIC MASS  MREL=M0+(E-POT)/(2C**2)    **
!     **  AND  M0 IS THE REST MASS                                        **
!     **                                                                  **
!     **  REMARKS:                                                        **
!     **  -  RELATIVISTIC CORRECTIONS FOR EKIN<0 ARE SWITCHED OFF!        **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID         ! GRID ID
      INTEGER(4),INTENT(IN) :: NR          ! #(RADIAL GRID POINTS
      REAL(8)   ,INTENT(IN) :: POT(NR)     ! POTENTIAL (MULTIPLY WITH Y0!)
      REAL(8)   ,INTENT(IN) :: E           ! ONE-PARTICLE ENERGY
      REAL(8)   ,INTENT(OUT):: DREL(NR)    ! RELATIVISTIC CORRECTION 
      INTEGER(4)            :: IR
      REAL(8)               :: C           ! SPEED OF LIGHT
      REAL(8)               :: PI,Y0      
      REAL(8)               :: DRELDOT(NR)
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS$GET('C',C)
      DREL(:)=-1.D0/(1.D0+2.D0*C**2/MAX(E-POT(:)*Y0,0.D0))
      DRELDOT(:)=DREL(:)**2*2.D0*C**2/MAX(E-POT(:)*Y0,0.D0)**2
      DO IR=1,NR
        IF(E.LT.POT(IR)*Y0) THEN 
          DRELDOT(IR)=0.D0
        END IF
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BOUNDSTATE(GID,NR,L,SO,DREL,G,NN,POT,E,PHI)
!     **                                                                  **
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND     **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                 **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G            **
!     **                                                                  **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE             **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.     **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.     **
!     **                                                                  **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND       **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT       **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                     **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,Z0,DX,XM,ZM
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: VAL,DER,DERO,RTEST
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      INTEGER(4)                 :: NN1
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+100 ! MAXIMUM FACTOR IN THE WAVE FUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      integer(4)                 :: nn1m,nn10,nn1p
      real(8)                    :: phip(nr),phim(nr)
      real(8)                    :: rcl
      real(8)                    :: swkb(nr)  ! phi=e^S
      real(8)                    :: aux(nr),aux1(nr)
      integer(4)                 :: irout,ircl
!     *********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RTEST=R(NR)
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!        E=MIN(X0,EMAX)
!
!       =======================================================================
!       ==  CUT OFF THE POTENTIAL                                            ==
!       =======================================================================
        POT1(:)=POT(:)
!       == FIND CLASSICAL TURNING POINT
        RCL=R(NR)
        IF(E.LT.POT(NR)*Y0) THEN
          DO IR=NR-1,1,-1
            IF(E.GT.POT(IR)*Y0) THEN
              IRCL=IR
              RCL=R(IR)-(POT(IR)-E/Y0)/(POT(IR+1)-POT(IR))*(R(IR+1)-R(IR))
! RCL=MIN(RTEST,R(NR))
              EXIT
            END IF
          ENDDO
          RTEST=RCL
        END IF
!
!       == USE WKB SOLUTION FOR THE SCHR.GL. FOR A CONSTANT POTENTIAL AND L=0
!       == TO ESTIMATE FACTOR FROM RTEST TO OUTERMOST POINT
        irout=nr
        IF(RTEST.LT.R(NR)) THEN
          AUX(:IRCL)=0.D0
          AUX(IRCL+1:)=SQRT(REAL(L*(L+1),kind=8)/R(ircl+1:)**2+2.D0*(POT(IRCL+1:)*Y0-E))
          CALL RADIAL$DERIVe(GID,NR,AUX,AUX1)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,Swkb)
!          Swkb(:)=Swkb(:)-0.5D0*LOG(AUX1(:))-LOG(R(:))
          Swkb(:)=Swkb(:)-LOG(R(:))
!         == DETERMINE IROUT WHERE THE WAVE FUNCTION CAN GROW BY A FACTOR 
!         == OF XMAX FROM THE CLASSICAL TURNING POINT
          SVAR=LOG(XMAX)
          DO IR=1,NR
            IF(Swkb(IR).GT.SVAR) THEN
              IROUT=IR-1
              EXIT
            END IF
          ENDDO
        END IF
        svar=pot(irout)
        pot1(:)=pot(:)
        pot1(irout:)=pot(irout)
!print*,'r(irout)',rcl,r(irout),pot1(nr)
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHRODINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(irout).GT.0.OR.PHI(irout).Le.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
!
!       =======================================================================
!       == CALCULATE PHASE SHIFT =========================================
!       =======================================================================
        NN1=0
        DO IR=3,irout
          IF(PHI(IR)*PHI(IR-1).LT.0.D0) NN1=NN1+1
        ENDDO
        Z0=-1.D0
        IF(NN1.Gt.NN) Z0=1.D0
!write(*,fmt='("iter",4i5,4e20.10)')I,L,NN1,NN,E,2.d0*DX,phi(irout),z0
        IF(ABS(2.d0*DX).LE.TOL) EXIT
!       =====================================================================
!       ==  BISECTION                                                      ==
!       =====================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('BOUNDSTATE')
      END IF
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E+2.d0*dx,IDIR,PHIp)
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E-2.d0*dx,IDIR,PHIm)
      NN1p=-nn
      NN1m=-nn
      NN10=-nn
      DO IR=3,irout
        IF(PHI(IR)*PHI(IR-1).LT.0.D0)   NN10=NN10+1
        IF(PHIm(IR)*PHIm(IR-1).LT.0.D0) NN1m=NN1m+1
        IF(PHIp(IR)*PHIp(IR-1).LT.0.D0) NN1p=NN1p+1
      ENDDO
!Print*,'nn1-nn ',nn1m,nn10,nn1p
!Print*,'e ',e-2.d0*dx,e,e+2.d0*dx
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      DO IR=1,NR
         IRMATCH=IR
         IF(POT(IR)-E/y0.GT.0.D0.OR.R(IR).GT.5.D0) EXIT
      ENDDO
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IDIR=-1
      GHOM(:)=0.D0
      if(irout.lt.nr) GHOM(IROUT)=1.D-5
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
      IF(.NOT.THOM) THEN     
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
      ELSE
        PHIINHOM(:)=0.D0
      END IF
!
!     =======================================================================
!     ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!     =======================================================================
      SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
      PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
      CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
      CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
      SVAR=(DERO-DER)/PHI(IRMATCH)
      PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
do ir=1,nr
  if(.not.(phi(ir).gt.0.d0.or.phi(ir).le.0.d0)) then
    print*,'error'
    print*,'phiin',phi(:irmatch-1)
    print*,'phiout',phi(irmatch:)
    print*,'svar ',svar
    call error$stop('boundstate')
  end if
enddo
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE old2BOUNDSTATE(GID,NR,L,SO,DREL,G,NN,POT,E,PHI)
!     **                                                                  **
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND     **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                 **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G            **
!     **                                                                  **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE             **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.     **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.     **
!     **                                                                  **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND       **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT       **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                     **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,Z0,DX,XM,ZM
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: VAL,DER,DERO,RTEST
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      INTEGER(4)                 :: NN1
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+100 ! MAXIMUM FACTOR IN THE WAVE FUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      integer(4)                 :: nn1m,nn10,nn1p
      real(8)                    ::phip(nr),phim(nr)
!     *********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RTEST=R(NR)
      ISTART=1
      X0=E
      DX=1.D0
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!        E=MIN(X0,EMAX)
!
!       =======================================================================
!       ==  CUT OFF THE POTENTIAL                                            ==
!       =======================================================================
!       == FIND CLASSICAL TURNING POINT
        RTEST=R(NR)
        IF(E.lt.POT(NR)*y0) THEN
          DO IR=NR-1,1,-1
            IF(E.gt.POT(IR)*y0) THEN
              RTEST=R(IR)-POT(IR)/(POT(IR+1)-POT(IR))*(R(IR+1)-R(IR))
              RTEST=MIN(RTEST,R(NR))
              EXIT
            END IF
          ENDDO
        END IF
!
!       == USE WKB SOLUTION FOR THE SCHR.GL. FOR A CONSTANT POTENTIAL AND L=0
!       == TO ESTIMATE FACTOR FROM RTEST TO OUTERMOST POINT
        IF(RTEST.LT.R(NR)) THEN
          SVAR=0.5D0*(LOG(XMAX*R(NR)/RTEST)/(R(NR)-RTEST))**2
          POT1(:)=MIN(POT(:),(E+SVAR)/y0)
!print*,'pot changed',svar,e
        ELSE
          POT1(:)=POT(:)
        END IF
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHRODINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(NR).GT.0.OR.PHI(NR).LT.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
!
!       =======================================================================
!       == CALCULATE PHASE SHIFT =========================================
!       =======================================================================
        NN1=0
        DO IR=3,NR
          IF(PHI(IR)*PHI(IR-1).LT.0.D0) NN1=NN1+1
        ENDDO
        Z0=-1.D0
        IF(NN1.Gt.NN) Z0=1.D0
!write(*,fmt='("iter",4i5,4e20.10)')I,L,NN1,NN,E,2.d0*DX,x0-xm,z0
        IF(ABS(2.d0*DX).LE.TOL) EXIT
!       =====================================================================
!       ==  BISECTION                                                      ==
!       =====================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('BOUNDSTATE')
      END IF
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E+2.d0*dx,IDIR,PHIp)
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E-2.d0*dx,IDIR,PHIm)
      NN1p=-nn
      NN1m=-nn
      NN10=-nn
      DO IR=3,NR-1
        IF(PHI(IR)*PHI(IR-1).LT.0.D0)   NN10=NN10+1
        IF(PHIm(IR)*PHIm(IR-1).LT.0.D0) NN1m=NN1m+1
        IF(PHIp(IR)*PHIp(IR-1).LT.0.D0) NN1p=NN1p+1
      ENDDO
!Print*,'nn1-nn ',nn1m,nn10,nn1p
!Print*,'e ',e-2.d0*dx,e,e+2.d0*dx
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      DO IR=1,NR
         IRMATCH=IR
         IF(POT(IR)-E/y0.GT.0.D0.OR.R(IR).GT.5.D0) EXIT
      ENDDO
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IDIR=-1
      GHOM(:)=0.D0
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
      IF(.NOT.THOM) THEN     
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
      ELSE
        PHIINHOM(:)=0.D0
      END IF
!
!     =======================================================================
!     ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!     =======================================================================
      SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
      PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
      CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
      CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
      SVAR=(DERO-DER)/PHI(IRMATCH)
!PRINT*,'SVAR',SVAR,R(IRMATCH),DER/PHI(IRMATCH),DERO/PHI(IRMATCH)
!do ir=1,nr
!  print*,r(ir),phi(ir),phiinhom(ir)
!enddo
!!$      IF(ABS(SVAR).GT.1.D-5) THEN
!!$        CALL ERROR$MSG('DERIVATIVES DO NOT MATCH')
!!$        CALL ERROR$R8VAL('STEP IN LOG. DERIVATIVE',SVAR)
!!$        CALL ERROR$R8VAL('RC',R(IRMATCH))
!!$        CALL ERROR$STOP('BOUNDSTATE')
!!$      END IF
!
      PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE oldBOUNDSTATE(GID,NR,L,SO,DREL,NN,POT,E,PHI)
!     **                                                                  **
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND     **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                 **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G            **
!     **                                                                  **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE             **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.     **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.     **
!     **                                                                  **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND       **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT       **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                     **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      REAL(8)                    :: G(NR)   !INHOMOGENITY
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,Y0,DX,XM,YM
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: PI
      REAL(8)                    :: VAL,DER,RTEST
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      INTEGER(4)                 :: NN1
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      LOGICAL(4)                 :: THOM
!     *********************************************************************
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      g(:)=0.d0
      RTEST=R(NR)
      ISTART=1
      X0=E
      DX=1.D0 
      CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
      DO I=1,NITER
        E=X0
        IDIR=1
        CALL RADIAL$SCHRODINGER(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHI)
        CALL RADIAL$VALUE(GID,NR,PHI,RTEST,VAL)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI,RTEST,DER)
        NN1=0
        DO IR=2,NR-1
          IF(R(IR+1).GT.RTEST) THEN
            IF(VAL*PHI(IR).LT.0.D0) NN1=NN1+1
            EXIT
          END IF
          IF(PHI(IR+1)*PHI(IR).LT.0.D0) NN1=NN1+1
        ENDDO
        Y0=-.5D0-ATAN(DER/VAL)/PI+REAL(NN1-NN,KIND=8)
!       =====================================================================
!       ==  BISECTION                                                      ==
!       =====================================================================
        CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
        IF(ABS(DX).LE.TOL) EXIT
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('BOUNDSTATE')
      END IF
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      DO IR=1,NR
         IRMATCH=IR
         IF(POT(IR)-E.GT.0.D0.OR.R(IR).GT.5.D0) EXIT
      ENDDO
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IDIR=-1
      GHOM(:)=0.D0
      CALL RADIAL$SCHRODINGER(GID,NR,POT,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
      IF(.NOT.THOM) THEN     
        CALL RADIAL$SCHRODINGER(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHIINHOM)
      ELSE
        PHIINHOM(:)=0.D0
      END IF
!
!     =======================================================================
!     ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!     =======================================================================
      SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
      PHI(IRMATCH:)=PHIINHOM(IRMATCH:)+SVAR*PHIHOM(IRMATCH:)
      RETURN
      END
!
!        ....................................................................
         SUBROUTINE SF_DO(GID,NR,LNX,LOX,UPHI,NPHI,THETA,LMRX,POT &
        &                ,LNXCORE,LOXCORE,PHIC,TPHIC,F,RHO)
         IMPLICIT NONE
         INTEGER(4),INTENT(IN)  :: GID
         INTEGER(4),INTENT(IN)  :: NR
         INTEGER(4),INTENT(IN)  :: LNX
         INTEGER(4),INTENT(IN)  :: LOX(LNX)
         REAL(8)   ,INTENT(IN)  :: UPHI(NR,LNX)
         INTEGER(4),INTENT(IN)  :: NPHI   !LMNX
         REAL(8)   ,INTENT(IN)  :: THETA(NPHI,NPHI)
         INTEGER(4),INTENT(IN)  :: LMRX
         REAL(8)   ,INTENT(IN)  :: POT(NR,LMRX)
         INTEGER(4),INTENT(IN)  :: LNXCORE
         INTEGER(4),INTENT(IN)  :: LOXCORE(LNXCORE)
         REAL(8)   ,INTENT(IN)  :: PHIC(NR,LNXCORE)
         REAL(8)   ,INTENT(IN)  :: TPHIC(NR,LNXCORE)
         REAL(8)   ,INTENT(IN)  :: F     ! CORE OCCUPATION (1 OR 2)
         REAL(8)   ,INTENT(OUT) :: RHO(NR,LMRX)
         INTEGER(4)             :: NC
         INTEGER(4)             :: NG
         INTEGER(4)             :: LMOFPHI(NPHI)
         REAL(8)                :: U(NR,NPHI)
         INTEGER(4),ALLOCATABLE :: LMOFG(:)
         REAL(8)   ,ALLOCATABLE :: G(:,:)
         REAL(8)   ,ALLOCATABLE :: TG(:,:)
         REAL(8)   ,ALLOCATABLE :: GG(:,:)
         REAL(8)   ,ALLOCATABLE :: GHG(:,:)
         REAL(8)   ,ALLOCATABLE :: UG(:,:)
         REAL(8)   ,ALLOCATABLE :: UHG(:,:)
         REAL(8)   ,ALLOCATABLE :: Q(:,:)
         integer(4)             :: i,j
real(8) :: aux(nr),svar,r(nr)
real(8),allocatable :: mat(:,:)
!        ***********************************************************************
         CALL SF_EXPAND(NR,LNX,LOX,NPHI,UPHI,LMOFPHI,U)
!
!        =======================================================================
!        ==  use core states as basis functions g                             ==
!        =======================================================================
         NG=SUM(2.D0*LOXCORE(:)+1.D0)
         NC=NG
         ALLOCATE(LMOFG(NG))
         ALLOCATE(G(NR,Ng))
         CALL SF_EXPAND(NR,LNXcore,LOXCORE,NG,PHIC,LMOFG,G)
         ALLOCATE(TG(NR,Ng))
         CALL SF_EXPAND(NR,LNXcore,LOXCORE,NG,TPHIC,LMOFG,TG)
!!$print*,'lmofg',lmofg
!!$call radial$r(gid,nr,r)
!!$do i=1,ng
!!$  do j=1,i
!!$    if(lmofg(i).ne.lmofg(j)) cycle
!!$    aux(:)=g(:,i)*g(:,j)*r(:)**2
!!$    call radial$integral(gid,nr,aux,svar)
!!$    if(i.eq.j) svar=svar-1.d0
!!$    if(abs(svar).gt.1.d-10) then
!!$      print*,'core overlap 1',i,j,lmofg(i),svar
!!$    end if
!!$  enddo
!!$enddo
!
!        =======================================================================
!        ==  calculate matrix elements                                        ==
!        =======================================================================
         ALLOCATE(GG(NG,NG))
         ALLOCATE(GHG(NG,NG))
         ALLOCATE(UG(NPHI,NG))
         ALLOCATE(UHG(NPHI,NG))
         CALL SF_SETUP(GID,NR,NG,NPHI,LMRX,LMOFG,LMOFPHI,G,TG,U,POT &
        &                ,GG,GHG,UG,UHG)
!print*,'ug  ',nphi,ng,ug
!print*,'uhg ',ng,nphi,uhg
!print*,'gg  ',ng,gg
!print*,'ghg ',ng,ghg

!
!        =======================================================================
!        ==  determine expansion coefficients for core states                 ==
!        =======================================================================
         ALLOCATE(q(NG,NC))
         q(:,:)=0.d0
         do i=1,ng
           q(i,i)=1.d0
         enddo
!!$print*,'marke 2',ng,nc
!!$allocate(mat(ng,ng))
!!$mat=matmul(transpose(q),matmul(gg,q))
!!$do i=1,ng
!!$  do j=1,i
!!$    svar=mat(i,j)
!!$    if(i.eq.j) svar=svar-1.d0
!!$    if(abs(svar).gt.1.d-10) then
!!$      print*,'core overlap 2',i,j,svar
!!$    end if
!!$  enddo
!!$enddo

         CALL SF_ITERATE(NPHI,NG,NC,GG,GHG,UG,UHG,THETA,F,Q)
!
!        =======================================================================
!        ==  calculate new one-center density                                 ==
!        =======================================================================
         CALL SF_MAKERHO(NR,NPHI,NG,NC,LMRX,LMOFG,LMOFPHI,U,G,UG,THETA,Q,RHO)
         deALLOCATE(q)
         deALLOCATE(GG)
         deALLOCATE(GHG)
         deALLOCATE(UG)
         deALLOCATE(UHG)
         deALLOCATE(TG)
         deALLOCATE(G)
         deALLOCATE(lmofg)
         RETURN
         END
!
!        ........................................................................
         SUBROUTINE SF_MAKERHO(NR,NPHI,NG,NC,LMRX,LMOFG,LMOFPHI,U,G,UG,THETA,Q,RHO)
         IMPLICIT NONE
         INTEGER(4),INTENT(IN)  :: NR
         INTEGER(4),INTENT(IN)  :: NG
         INTEGER(4),INTENT(IN)  :: NPHI
         INTEGER(4),INTENT(IN)  :: NC
         INTEGER(4),INTENT(IN)  :: LMRX
         INTEGER(4),INTENT(IN)  :: LMOFPHI(NPHI)
         INTEGER(4),INTENT(IN)  :: LMOFG(NG)
         REAL(8)   ,INTENT(IN)  :: U(NR,NPHI)
         REAL(8)   ,INTENT(IN)  :: G(NR,NG)
         REAL(8)   ,INTENT(IN)  :: UG(NPHI,NG)
         REAL(8)   ,INTENT(IN)  :: THETA(NPHI,NPHI)
         REAL(8)   ,INTENT(IN)  :: Q(NG,NC)
         REAL(8)   ,INTENT(OUT) :: RHO(NR,LMRX)
         INTEGER(4)             :: I,J,LM
         REAL(8)                :: CG
         REAL(8)                :: MAT1(NG,NPHI),MAT2(NG,NPHI)
         REAL(8)                :: GUFAC(NG,NPHI)
         REAL(8)                :: GGFAC(NG,NG)
!        ***********************************************************************
         RHO(:,:)=0.D0
!
!        =======================================================================
!        ==   <R|U>D<U|R>                                                     ==
!        =======================================================================
         DO I=1,NPHI
           DO J=1,NPHI
             DO LM=1,LMRX
               CALL CLEBSCH(LMOFPHI(I),LMOFPHI(J),LM,CG)
               IF(CG.EQ.0.D0) CYCLE
               RHO(:,LM)=RHO(:,LM)+THETA(I,J)*CG*U(:,I)*U(:,J)
             ENDDO
           ENDDO
         ENDDO
!
!        =======================================================================
!        ==   <R|U>D<U|R>                                                     ==
!        =======================================================================
!        ==  MAT1=Q^DAGGER Q <G|U> =============================================
         MAT1=MATMUL(TRANSPOSE(Q),MATMUL(Q,TRANSPOSE(UG)))
         MAT2=MATMUL(MAT1,THETA)
         GUFAC=-2.D0*MAT2
         DO I=1,NG
           DO J=1,NPHI
             DO LM=1,LMRX
               CALL CLEBSCH(LMOFG(I),LMOFPHI(J),LM,CG)
               IF(CG.EQ.0.D0) CYCLE
               RHO(:,LM)=RHO(:,LM)+GUFAC(I,J)*CG*G(:,I)*U(:,J)
             ENDDO
           ENDDO
         ENDDO
!
!        =======================================================================
!        ==   <R|U>D<U|R>                                                     ==
!        =======================================================================
         GGFAC=MATMUL(MAT2,TRANSPOSE(MAT1))+MATMUL(TRANSPOSE(Q),Q)
         DO I=1,NG
           DO J=1,NG
             DO LM=1,LMRX
               CALL CLEBSCH(LMOFG(I),LMOFG(J),LM,CG)
               IF(CG.EQ.0.D0) CYCLE
               RHO(:,LM)=RHO(:,LM)+GGFAC(I,J)*CG*G(:,I)*G(:,J)
             ENDDO
           ENDDO
         ENDDO
         RETURN
         END
!
!        ........................................................................
         SUBROUTINE SF_EXPAND(NR,LNX,LOX,LMNX,UPHI,LMOFPHI,U)
!        ** TAKES UPHI STORED WITH ARGUMENT LN AND MAKES THE CORRESPONDING
!        ** FUNCTION WITH ARGUMENT LMN
         IMPLICIT NONE
         INTEGER(4),INTENT(IN)  :: NR
         INTEGER(4),INTENT(IN)  :: LNX
         INTEGER(4),INTENT(IN)  :: LOX(LNX)
         INTEGER(4),INTENT(IN)  :: LMNX
         REAL(8)   ,INTENT(IN)  :: UPHI(NR,LNX)
         INTEGER(4),INTENT(OUT) :: LMOFPHI(LMNX)
         REAL(8)   ,INTENT(OUT) :: U(NR,LMNX)
         INTEGER(4)             :: LMN,L,M,LN
!        ***********************************************************************
         u(:,:)=0.d0
         LMN=0
         DO LN=1,LNX
           L=LOX(LN)
           DO M=1,2*L+1
             LMN=LMN+1
             U(:,LMN)=UPHI(:,LN)
             LMOFPHI(LMN)=L**2+M
           ENDDO
         ENDDO
         RETURN
         END                             
!
!        ........................................................................
         SUBROUTINE SF_SETUP(GID,NR,NG,NPHI,LMX,LMOFG,LMOFPHI,G,TG,U,POT &
        &                ,GG,GHG,UG,UHG)
         IMPLICIT NONE
         INTEGER(4),INTENT(IN)  :: GID
         INTEGER(4),INTENT(IN)  :: NR
         INTEGER(4),INTENT(IN)  :: NG
         INTEGER(4),INTENT(IN)  :: LMOFG(NG)
         INTEGER(4),INTENT(IN)  :: LMX
         REAL(8)   ,INTENT(IN)  :: G(NR,NG)
         REAL(8)   ,INTENT(IN)  :: TG(NR,NG)
         REAL(8)   ,INTENT(IN)  :: POT(NR,LMX)
         INTEGER(4),INTENT(IN)  :: NPHI
         INTEGER(4),INTENT(IN)  :: LMOFPHI(NPHI)
         REAL(8)   ,INTENT(IN)  :: U(NR,NPHI)
         REAL(8)   ,INTENT(OUT) :: GG(NG,NG)
         REAL(8)   ,INTENT(OUT) :: GHG(NG,NG)
         REAL(8)   ,INTENT(OUT) :: UG(NPHI,NG)
         REAL(8)   ,INTENT(OUT) :: UHG(NPHI,NG)
         REAL(8)                :: R(NR)
         REAL(8)                :: AUX(NR)
         REAL(8)                :: VG(NR)
         REAL(8)                :: SVAR
         REAL(8)                :: CG
         INTEGER(4)             :: I,J
         INTEGER(4)             :: LM1,LM2,LM3,LMX1
!        ******************************************************************
         CALL RADIAL$R(GID,NR,R)
         gg(:,:)=0.d0
         ghg(:,:)=0.d0
         ug(:,:)=0.d0
         uhg(:,:)=0.d0
!
!        ==================================================================
!        == <G|G> AND <G|T|G> =============================================
!        ==================================================================
         DO I=1,NG
           LM1=LMOFG(I)
           DO J=1,I
             IF(LM1.eq.LMOFG(J)) THEN
               AUX(:)=G(:,I)*G(:,J)*R(:)**2
               CALL RADIAL$INTEGRAL(GID,NR,AUX,GG(I,J))
               GG(J,I)=GG(i,j)
               AUX(:)=G(:,I)*TG(:,J)*R(:)**2
               CALL RADIAL$INTEGRAL(GID,NR,AUX,GHG(I,J))
               GHG(J,I)=GHG(i,j)
             END IF
           ENDDO
         ENDDO
!
!        ==================================================================
!        == <U|G> AND <U|T|G>
!        ==================================================================
         DO I=1,NPHI
           LM1=LMOFPHI(I)
           DO J=1,NG
             IF(LM1.eq.LMOFG(J)) THEN
               AUX(:)=U(:,I)*G(:,J)*R(:)**2
               CALL RADIAL$INTEGRAL(GID,NR,AUX,UG(I,J))
               AUX(:)=U(:,I)*TG(:,J)*R(:)**2
               CALL RADIAL$INTEGRAL(GID,NR,AUX,UHG(I,J))
             END IF
           ENDDO
         ENDDO
!
!        ==================================================================
!        == ADD <G|V|G> AND <U|V|G>                                      ==
!        ==================================================================
         DO I=1,NG
           LM1=LMOFG(I)
           LMX1=MAX(MAXVAL(LMOFG(1:I)),MAXVAL(LMOFPHI))
           DO LM2=1,LMX1
             VG(:)=0.D0
             DO LM3=1,LMX
               CALL CLEBSCH(LM1,LM2,LM3,CG)
               VG(:)=VG(:)+POT(:,LM2)*G(:,I) 
             ENDDO
             VG(:)=VG(:)*R(:)**2
             DO J=1,I
               IF(LMOFG(J).NE.LM2) CYCLE
               AUX(:)=G(:,J)*VG(:)
               CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
               GHG(I,J)=GHG(I,J)+SVAR
               IF(I.NE.J) GHG(J,I)=GHG(J,I)+SVAR
             ENDDO
             DO J=1,NPHI
               IF(LMOFPHI(J).NE.LM2) CYCLE
               AUX(:)=U(:,J)*VG(:)
               CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
               UHG(J,I)=UHG(J,I)+SVAR
             ENDDO
           ENDDO
         ENDDO
         RETURN
         END
!
!        ..................................................................
         SUBROUTINE SF_ITERATE(NPHI,NG,NC,GG,GHG,UG,UHG,THETA,F,Q)
         IMPLICIT NONE
         INTEGER(4),INTENT(IN)   :: NPHI
         INTEGER(4),INTENT(IN)   :: NG
         INTEGER(4),INTENT(IN)   :: NC
         REAL(8)   ,INTENT(IN)   :: GG(NG,NG)
         REAL(8)   ,INTENT(IN)   :: GHG(NG,NG)
         REAL(8)   ,INTENT(IN)   :: UG(NPHI,NG)
         REAL(8)   ,INTENT(IN)   :: UHG(NPHI,NG)
         REAL(8)   ,INTENT(IN)   :: THETA(NPHI,NPHI)
         REAL(8)   ,INTENT(IN)   :: F    ! CORE OCCUPATION (1.OR 2.)
         REAL(8)   ,INTENT(INOUT):: Q(NG,NC)
         INTEGER(4),PARAMETER    :: NITER=100
         REAL(8)   ,PARAMETER    :: ALPHAMIX=1.D-2
         REAL(8)   ,PARAMETER    :: TOL=1.D-6
         INTEGER(4)              :: ITER
         REAL(8)                 :: E1,E2
         REAL(8)                 :: DEDQ1(NG,NC)
         REAL(8)                 :: DEDQ2(NG,NC)
         REAL(8)                 :: Q2(NG,NC)
         REAL(8)                 :: ALPHANEU
         REAL(8)                 :: SVAR1,SVAR2
REAL(8)                 :: DEDQsave(NG,NC)
integer(4)              :: i
!        ******************************************************************
!rint*,'marke iterate'
         DO ITER=1,NITER
           CALL SF_ORTHO(NG,NC,GG,Q)
           CALL SF_FORCE(NPHI,NG,NC,GG,GHG,UG,UHG,THETA,F &
        &                      ,Q,E1,DEDQ1)
if(iter.eq.1) then
!  dedqsave=dedq1
else
!  dedq1=dedqsave
end if
           SVAR1=SUM(DEDQ1*DEDQ1)
           IF(SQRT(SVAR1/REAL(NC)).le.TOL) then
             print*,'sf iterate converged after ',iter,' iterations'
             return
           end if
           Q2=Q-ALPHAMIX*DEDQ1
           CALL SF_ORTHO(NG,NC,GG,Q2)
           CALL SF_FORCE(NPHI,NG,NC,GG,GHG,UG,UHG,THETA,F &
        &                      ,Q2,E2,DEDQ2)
           SVAR2=SUM(DEDQ1*DEDQ2)
!print*,'e ',e1,e1-alphamix*svar1,e2,e2+alphamix*svar2
           ALPHANEU=ALPHAMIX*SVAR1/(SVAR2-SVAR1)
           Q=Q-ALPHANEU*DEDQ1
         END DO
         CALL ERROR$MSG('LOOP NOT NOT CONVERGED')
         CALL ERROR$STOP('SF_ITERATE')
         RETURN
         END
!
!        ..................................................................
         SUBROUTINE SF_ORTHO(NG,NC,GG,Q)
!        **                                                              **
!        ** ORTHOGONALIZES THE CORE STATES                               **
!        **                                                              **
!        **   QDAGGER*GG*Q=1                                             **
!        **                                                              **
!        ******************************************************************
         IMPLICIT NONE
         INTEGER(4),INTENT(IN)   :: NG
         INTEGER(4),INTENT(IN)   :: NC
         REAL(8)   ,INTENT(IN)   :: GG(NG,NG)      ! 'B'
         REAL(8)   ,INTENT(INOUT):: Q(NG,NC)       ! CORE EXPANSION COEFFICIENTS
         INTEGER(4),PARAMETER    :: NITER=1000
         REAL(8)   ,PARAMETER    :: TOL=1.D-6
         REAL(8)                 :: BQ(NG,NC)
         REAL(8)                 :: MAT(NC,NC)
         REAL(8)                 :: MATINV(NC,NC)
         INTEGER(4)              :: I,ITER,j
         real(8)                 :: svar
!        *********************************************************************
         BQ=MATMUL(GG,Q)
         MAT=MATMUL(TRANSPOSE(BQ),BQ)
         CALL LIB$INVERTR8(NC,MAT,MATINV)
         BQ=0.5D0*MATMUL(BQ,MATINV)
         DO ITER=1,NITER
           MAT=MATMUL(TRANSPOSE(Q),MATMUL(GG,Q))
           DO I=1,NC
             MAT(I,I)=MAT(I,I)-1.D0
           ENDDO
           svar=MAXVAL(ABS(MAT))
           IF(svar.LT.TOL) RETURN
           Q=Q-MATMUL(BQ,MAT)
         ENDDO
!
!        == loop not converged; provide details ==============================
do i=1,nc
  do j=1,i
    if(abs(mat(i,j)).gt.tol.or.abs(mat(j,i)).gt.tol) then
      print*,'i,j ',i,j,mat(i,j),mat(j,i)
    end if
  enddo
enddo
         call error$msg('LOOP NOT CONVERGED')
         call error$stop('sf_ortho')
         RETURN
         END
!
!        ..................................................................
         SUBROUTINE SF_FORCE(NPHI,NG,NC,GG,GHG,UG,UHG,THETA,F,Q,E,DEDQ)
!        **                                                              **
!        ** DETERMINES THE ENERGY AND THE FORCE AND THE CONSTRAINT FORCE **
!        ** ACTING ON THE EXPANSION COEFFICIENTS OF THE CORE STATES      **
!        **                                                              **
!        **   E=CONST + TR[-2*THETA*C*Q*QDAGGER*DDAGGER                  **
!        **                +QDAGGER*A*Q*(F+(C*Q)DAGGER*THETA*C*Q)]       **
!        **                                                              **
!        ******************************************************************
         IMPLICIT NONE
         INTEGER(4),INTENT(IN) :: NPHI
         INTEGER(4),INTENT(IN) :: NG
         INTEGER(4),INTENT(IN) :: NC
         REAL(8)   ,INTENT(IN) :: GG(NG,NG)      ! 'B'
         REAL(8)   ,INTENT(IN) :: GHG(NG,NG)     ! 'A' 
         REAL(8)   ,INTENT(IN) :: UG(NPHI,NG)    ! 'C'
         REAL(8)   ,INTENT(IN) :: UHG(NPHI,NG)   ! 'D'
         REAL(8)   ,INTENT(IN) :: THETA(NPHI,NPHI) ! DENSITY MATRIX
         REAL(8)   ,INTENT(IN) :: F              ! CORE OCCUPATIONS (1 OR 2)
         REAL(8)   ,INTENT(IN) :: Q(NG,NC)       ! CORE EXPANSION COEFFICIENTS
         REAL(8)   ,INTENT(OUT):: E              ! NON-SCF ENERGY
         REAL(8)   ,INTENT(OUT):: DEDQ(NG,NC)    ! ENERGY DERIVATIVE
         REAL(8)               :: FC(NG,NC)      ! FORCE OF CONSTRAINS
         REAL(8)               :: THETAC(NPHI,NG)
         REAL(8)               :: AQQCMINUSD(NG,NPHI)
         REAL(8)               :: MAT1(NG,NG)
         REAL(8)               :: MAT(NC,NC)
         REAL(8)               :: MATINV(NC,NC)
         INTEGER(4)            :: I
!        *********************************************************************
         THETAC=MATMUL(THETA,UG)
         AQQCMINUSD=-TRANSPOSE(UHG) &
        &          +MATMUL(GHG,MATMUL(TRANSPOSE(Q),MATMUL(Q,TRANSPOSE(UG))))
         MAT1=MATMUL(AQQCMINUSD,THETAC)
         MAT1=MAT1+TRANSPOSE(MAT1)+F*GHG
         DEDQ=2.D0*MATMUL(MAT1,Q)
         AQQCMINUSD=-TRANSPOSE(UHG)+AQQCMINUSD
         MAT1=MATMUL(AQQCMINUSD,THETAC)+F*GHG
         E=SUM(Q*MATMUL(MAT1,Q))
!
!        == ORTHOGONALIZE GRADIENT TO EXISTING CORE STATES
         FC=MATMUL(GG,Q)
         MAT=MATMUL(TRANSPOSE(FC),FC)
         CALL LIB$INVERTR8(NC,MAT,MATINV)
         DEDQ=DEDQ-MATMUL(MATMUL(FC,MATMUL(MATINV,TRANSPOSE(FC))),DEDQ)
         RETURN
         END
!
!..................................................................................
MODULE BROYDEN_MODULE
LOGICAL(4)         :: TON=.FALSE.
INTEGER(4)         :: NSTEPX=0
INTEGER(4)         :: NSTEP=0
INTEGER(4)         :: Nx=0
REAL(8)            :: ALPHA
REAL(8),ALLOCATaBLE :: XPREV(:,:)
REAL(8),ALLOCATaBLE :: YPREV(:,:)
END MODULE BROYDEN_MODULE
!      .............................................................................
       subroutine BROYDEN$NEW(NX_,NSTEPX_,ALPHA_)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NX_
       INTEGER(4),INTENT(IN)    :: NSTEPX_
       REAL(8)   ,INTENT(IN)    :: ALPHA_
!      *****************************************************************************
       IF(TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT ALREADY IN USE')
         CALL ERROR$STOP('BROYDEN$NEW')
       END IF
       TON=.TRUE.
       NSTEP=0
       NX=NX_
       NSTEPX=NSTEPX_
       ALPHA=ALPHA_
       ALLOCATE(XPREV(NX,NSTEPX))
       ALLOCATE(YPREV(NX,NSTEPX))
       RETURN
       END
!      .............................................................................
       subroutine BROYDEN$CLEAR
       USE BROYDEN_MODULE
       IMPLICIT NONE
!      *****************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$CLEAR')
       END IF
       TON=.FALSE.
       NSTEPX=0
       NSTEP=0
       NX=0
       ALPHA=0.D0
       DEALLOCATE(XPREV)
       DEALLOCATE(YPREV)
       RETURN
       END
!      .............................................................................
       subroutine BROYDEN$STEP(NX_,X,Y)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NX_
       REAL(8)   ,INTENT(INOUT) :: X(NX_)
       REAL(8)   ,INTENT(IN)    :: Y(NX_)
       REAL(8)   ,ALLOCATABLE   :: DX(:,:)
       REAL(8)   ,ALLOCATABLE   :: DY(:,:)
       REAL(8)   ,ALLOCATABLE   :: B(:,:)
       REAL(8)   ,ALLOCATABLE   :: BINV(:,:)
 REAL(8)   ,ALLOCATABLE   :: w(:,:)
       integer(4)               :: i
!      *****************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
       if(nx_.ne.nx) then
         CALL ERROR$MSG('size inconsistent')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
!print*,'nstep',nstep
!
!      =================================================================
!      == simple mixing in the first step                             ==
!      =================================================================
       if(nstep.eq.0) then
         if(nstepx.gt.0)then
           nstep=1
           XPREV(:,1)=X(:)     
           YPREV(:,1)=Y(:)     
         END IF
         X=X+ALPHA*Y
         return
       end if
!
!      =================================================================
!      == determine inverse hessian alpha+dx otimes dy                ==
!      =================================================================
       ALLOCATE(DX(NX,NSTEP))
       ALLOCATE(DY(NX,NSTEP))
       DO I=1,NSTEP
         DY(:,I)=YPREv(:,I)-Y(:)  
         DX(:,I)=XPREv(:,I)-X(:)+ALPHA*DY(:,I)
       ENDDO
       ALLOCATE(B(NSTEP,NSTEP))
       ALLOCATE(BINV(NSTEP,NSTEP))
       B=MATMUL(TRANSPOSE(DY),DY)   !OVERLAP MATRIX OF DY
!print*,'b',b
       CALL LIB$INVERTR8(NSTEP,B,BINV)           
allocate(w(nx,nstep))
w=MATMUL(DY,BINV)            !NEW DY IS BIORTHOnormal TO OLD DY
!print*,'w ',matmul(transpose(w),dy)
deallocate(w)
       DY=MATMUL(DY,BINV)            !NEW DY IS BIORTHOnormal TO OLD DY
       DEALLOCATE(B)
       DEALLOCATE(BINV)
!
!      =================================================================
!      == store history                                               ==
!      =================================================================
       IF(NSTEP.LT.NSTEPX)NSTEP=NSTEP+1
       DO I=NSTEP,2,-1
         YPREV(:,I)=YPREV(:,I-1)
         XPREV(:,I)=XPREV(:,I-1)
       ENDDO
       XPREV(:,1)=X(:)     
       YPREV(:,1)=Y(:)     
!
!      =================================================================
!      == predict new vector                                          ==
!      =================================================================
       X=X+ALPHA*Y-MATMUL(DX,MATMUL(TRANSPOSE(DY),Y))
       DEALLOCATE(DX)
       DEALLOCATE(DY)
       RETURN
       END

!
!     ..................................................................
      SUBROUTINE AUGMENTATION_RHO(NR,LNX,LOX,PHI,LMNXX &
     &                           ,DENMAT,LMRX,RHOL)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NR
      INTEGER(4), INTENT(IN) :: LNX
      INTEGER(4), INTENT(IN) :: LOX(LNX)
      INTEGER(4), INTENT(IN) :: LMNXX
      INTEGER(4), INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: DENMAT(LMNXX,LMNXX)
      REAL(8)    ,INTENT(IN) :: PHI(NR,LNX)
      REAL(8)    ,INTENT(OUT):: RHOL(NR,LMRX)
      INTEGER(4)             :: LMR
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
                  SVAR=CG*DENMAT(LMN1,LMN2)
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
      SUBROUTINE AUGMENTATION_QLM(gid,NR,LMRX &
     &                           ,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE MULTIPOLE MOMENTS OF THE DEVIATION OF             **
!     **  ALL-ELECTRON AND THE PSEUDO DENSITY FROM ONE CENTER         **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: gid
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)   ,INTENT(IN) :: PSRHO(NR,LMRX)
      REAL(8)   ,INTENT(OUT):: QLM(LMRX)
      REAL(8)               :: DWORK(NR)
      REAL(8)               :: PI,Y0
      REAL(8)               :: RES
      INTEGER(4)            :: LM,L,IR
      real(8)               :: r(nr)
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      call radial$r(gid,nr,r)
      QLM(:)=0.D0
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        DWORK(:)=(AERHO(:,LM)-PSRHO(:,LM))*R(:)**(L+2)
        CALL RADIAL$INTEGRAL(gid,NR,DWORK,QLM(LM))
      ENDDO
!
      DO IR=1,NR
        DWORK(:)=(AECORE(:)-PSCORE(:))*R(:)**2
      ENDDO
      CALL RADIAL$INTEGRAL(gid,NR,DWORK,RES)
      QLM(1)=QLM(1)+RES-AEZ*Y0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_XC(gid,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY              **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID            **
!     **  TIMES REAL SPHERICAL HARMONICS                              **
!     **                                                              **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE                  **
!     **  SPHERICAL CONTRIBUTION OF THE DENSITY UP TO QUADRATIC       **
!     **  ORDER IN THE NON-SPHERICAL CONTRIBUTIONS                    **
!     **                                                              **
!     **  EXC = EXC(XVAL(L=0)*Y0)                                     **
!     **      + 0.5 * D2[EXC]/D[XVAL(L=0)*Y0]**2 * XVAL(L>0)**2       **
!     **                                                              **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)      **
!     **  IS AN SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.     **
!     **                                                              **
!     **  DEPENDECIES:                                                **
!     **    DFT                                                       **
!     **    TIMING                                                    **
!     **    TRACE                                                     **
!     **                                                              **
!     **  REMARKS: THE GRADIENTS ARE CORRECT ONLY IF DFT SUPPORTS     **
!     **    THIRD DERIVATIVES OF THE XC-ENERGY                        **
!     **   - WHEN USING SELFTEST ON THIS ROUTINE, THEN                **
!     **     D(EXC)/DRHO(I)=POT(I)*DEX*R(I)**3                        **
!     **     AND THE VALUES AT LARGE RADII MUST BE SURPRESSED         **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     *********** P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1996) **
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: gid
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: B(:,:)    
      REAL(8)   ,ALLOCATABLE:: VB(:,:)    
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: Y0
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM,LM1,LM2,LM3
      INTEGER(4)            :: IMAX
      REAL(8)               :: RI,FAC
      REAL(8)               :: CG0LL
      REAL(8)               :: CG
      REAL(8)               :: SVAR
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
      REAL(8)               :: WORK3(NR)
      REAL(8)   ,PARAMETER  :: TINY=1.D-300
!     ******************************************************************
      CALL TRACE$PUSH('AUGMENTATION_XC')
      EXC=0.D0
!
!     ==================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                      ==
!     ==================================================================
      CALL DFT$GETL4('GC',TGRA)
      PI=4.D0*DATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/DSQRT(FOURPI)
      CG0LL=Y0
      call radial$r(gid,nr,r)
!
!     ==================================================================
!     ==  OBTAIN SPIN DENSITY                                         ==
!     ==================================================================
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      ALLOCATE(RHO(NR,LMRX,NSPIN))
      ALLOCATE(GRHO(NR,LMRX,NSPIN))
      ALLOCATE(VRHO(NR,LMRX,NSPIN))
      RHO(:,:,1)=RHOIN(:,:,1)
      IF(NDIMD.EQ.2) THEN
        RHO(:,:,2)=RHOIN(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        ALLOCATE(B(NR,LMRX))
!       == EVALUATE SQUARE OF THE SPIN DENSITY =========================
        DO LM1=1,LMRX
          B(:,LM1)=0.D0
        ENDDO
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
!       ==  TRANSFORM B INTO C=B/(2*B0*Y0) =============================
        WORK(:)=0.5D0/(B(:,1)*Y0+TINY)
        DO LM1=2,LMRX
          B(:,LM1)=B(:,LM1)*WORK(:)
        ENDDO
!       == CALCULATE SPIN DENSITY ======================================
        RHO(:,:,2)=0.D0
        DO LM2=2,LMRX
          DO LM3=LM3,LMRX
            WORK(:)=B(:,LM2)*B(:,LM3)
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
          RHO(:,LM1,2)=RHO(:,LM1,2)+B(:,LM1)
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
!     ==================================================================
!     ==  CALCULATE RADIAL GRADIENT OF THE DENSITY                    ==
!     ==================================================================
      CALL TRACE$PASS('BEFORE GRADIENTS')
      IF(TGRA) THEN
        GRHO(:,:,:)=0.D0
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            CALL RADIAL$DERIVE(gid,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
          ENDDO
        ENDDO
      ELSE
        GRHO(:,:,:)=0.D0
      END IF
!
!     ==================================================================
!     ==  DEFINE VECTOR (RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS)     ==
!     ==================================================================
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
              L=INT(SQRT(REAL(LM-1))+1.D-5)
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
!     ==================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY         ==
!     ==================================================================
      CALL TRACE$PASS('BEFORE DFT')
      WORK1(:)=0.D0
      XDER(:,:,:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ========================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ========================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)
!       == NOW CALCULATE ENERGY DENSITY AND DERIAVTIVES =================
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
      CALL RADIAL$INTEGRAL(gid,NR,WORK1(:)*R(:)**2,EXC)
!
!     ==================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                     ==
!     ==================================================================
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
!           == FIRST RESOLVE XVAL(:,II,1) =============================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1))+1.D-5)
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
!           == NOW RESOLVE XVAL(:,II,LM) ==============================
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
!     ==================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS             ==
!     ==  V = V -1/R**2 D/DR [ R**2 VGRHO ]                           ==
!     ==  V = V -[2/R VGRHO+ D/DR VGRHO ]                             ==
!     ==================================================================
      IF(TGRA) THEN
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
!           == FIRST ALTERNATIVE
!           CALL RADIAL$DERIVE(gid,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(gid,NR,WORK2,WORK1)
            WORK1(:)=WORK1(:)/R(:)**2
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(VGRHO)
!
!     ==================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS             ==
!     ==================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=VRHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        ALLOCATE(VB(NR,LMRX))
        VB(:,:)=0.D0
        DO LM2=1,LMRX
          DO LM3=2,LMRX
            WORK(:)=VRHO(:,LM2,2)*B(:,LM3)
            DO LM1=2,LMRX
              CALL CLEBSCH(LM1,LM2,LM3,CG)
              VB(:,LM1)=VB(:,LM1)+CG*WORK(:)
            ENDDO
          ENDDO
        ENDDO
        DO LM1=2,LMRX
          VB(:,LM1)=VRHO(:,LM1,2)-VB(:,LM1)
        ENDDO
        DO LM1=2,LMRX
          VB(:,1)=VB(:,1)+VB(:,LM1)*B(:,LM1)
        ENDDO
        WORK(:)=0.5D0/SQRT(B(:,1)*Y0+TINY)
        DO LM1=1,LMRX
          VB(:,LM1)=VB(:,LM1)*WORK(:)
        ENDDO
        VB(:,1)=VB(:,1)*2.D0*Y0
        WORK(:)=0.D0
        DO LM1=1,LMRX
          WORK(:)=WORK(:)+VRHO(:,LM1,2)*RHO(:,LM1,2)
        ENDDO
        VB(:,1)=VB(:,1)+0.5D0*WORK(:)/(B(:,1)+TINY)
!       == TRANSFORM VB=DE/DB INTO POTENTIAL FOR SPIN DENSITY ==========
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
      END IF     
      DEALLOCATE(RHO)
      DEALLOCATE(GRHO)
      DEALLOCATE(VRHO)
                      CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_PSHARTREE(GID,NR,LMRX,PSRHOC,PSRHO &
     &            ,VADD,RCSM,QLM,VQLM,RHOB,PSPOT,PSEH)
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
      INTEGER(4) ,INTENT(IN) :: gid
      INTEGER(4) ,INTENT(IN) :: NR
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: psRHOC(NR)
      REAL(8)    ,INTENT(IN) :: PSRHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: VADD(NR)
      REAL(8)    ,INTENT(IN) :: RCSM
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(in) :: VQLM(LMRX)
      REAL(8)    ,INTENT(in) :: rhob
      REAL(8)    ,INTENT(OUT):: PSPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: PSEH
      REAL(8)                :: r(NR)
      REAL(8)                :: AUX1(NR)
      REAL(8)                :: rho1(NR)
      REAL(8)                :: pot(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: g(NR)
      REAL(8)                :: rhohat(NR,lmrx)
      REAL(8)                :: PI
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      INTEGER(4)             :: LM,l,m,IR,lx
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$r(gid,NR,r)
!
!     ==================================================================
!     ==  INITIALIZE ENERGIES AND POTENTIALS TO  ZERO                 ==
!     ==================================================================
      pse(:)=0.d0
      PSPOT(:,:)=0.D0
!
!     ==================================================================
!     ==  CONSTRUCT COMPENSATION DENSITY                              ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      LX=INT(SQRT(REAL(LMRX-1)+.01D0))
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
!     ==   add unscreening potential vadd                             ==
!     ==================================================================
      PSE(:)    =VADD(:)*(PSRHO(:,1)+psrhoc(:))
      PSPOT(:,1)=PSPOT(:,1)+VADD(:)
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        CALL GAUSSN(L,ALPHA,CL)
        SVAR=QLM(LM)*CL
        RHO1(:)=PSRHO(:,LM)+RHOHAT(:,LM)
        IF(L.EQ.0)RHO1(:)=RHO1(:)+PSRHOC(:)
        CALL RADIAL$POISSON(GID,NR,L,rho1,pot)
        PSE(:)=PSE(:)+0.5D0*Pot(:)*RHO1(:)
        PSPOT(:,LM)=PSPOT(:,LM)+Pot(:)
      ENDDO
!
!     ==================================================================
!     ==  ADD EXTERNAL POTENTIAL                                      ==
!     ==================================================================
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1)+.01D0))
        POT(:)=VQLM(LM)*R(:)**L
        PSPOT(:,LM)=PSPOT(:,LM)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
!       RHO1(:)=PSRHO(:,LM)+RHOHAT(:,LM)
!       IF(L.EQ.0) RHO1(:)=RHO1(:)+PSRHOC(:)
!       PSE(:)=PSE(:)+PSRHO(:,LM)*POT(:)
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
      SUBROUTINE AUGMENTATION_AEhartree(GID,NR,LMRX,RHOC,AERHO &
     &                                 ,VQLM,RHOB,aepot,AEEH)
!     ******************************************************************
!     **  electrostatic energy of the all-electron one-center density **
!     **  including the external potential and the potential of the   **
!     **  compensating charge background for charged systems          **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID
      INTEGER(4) ,INTENT(IN) :: NR
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: RHOC(NR)   ! core density
      REAL(8)    ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)    ,INTENT(IN) :: RHOB       ! compensating background
      REAL(8)    ,INTENT(OUT):: AEPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: AEEH       ! energy
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX(NR)
      REAL(8)                :: POT(NR)
      REAL(8)                :: AEE(NR)  ! energy density
      REAL(8)                :: svar
      REAL(8)                :: pi
      INTEGER(4)             :: LM
      INTEGER(4)             :: L
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
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
      CALL SETUP$GETR8A('NUCPOT',NR,AUX)
      POT(:)  = POT(:)+AUX(:)
      AEE(:)    = POT(:)* AERHO(:,1)
      AEPOT(:,1)= AEPOT(:,1)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
!      AEE(:)=AE(:)+RHOC(:)*VQLM(1)
aux(:)=aee(:)*r(:)**2
CALL RADIAL$INTEGRAL(GID,NR,aux,svar)
print*,'el c-v int', svar
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
aux(:)=0.d0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1)+.01D0))
        CALL RADIAL$POISSON(GID,NR,L,AERHO(:,LM),POT)
        AEE(:)=AEE(:)+0.5D0*POT(:)*AERHO(:,LM)
        AEPOT(:,LM)=AEPOT(:,LM)+POT(:)
aux(:)=aux(:)+0.5D0*POT(:)*AERHO(:,LM)*r(:)**2
      ENDDO
CALL RADIAL$INTEGRAL(GID,NR,aux,svar)
print*,'el v-v int',svar
!
aux(:)=aee(:)*r(:)**2
CALL RADIAL$INTEGRAL(GID,NR,aux,svar)
print*,'el hartree int ',svar
!
!     ==================================================================
!     ==  ADD EXTERNAL POTENTIAL                                      ==
!     ==================================================================
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1)+.01D0))
        pot(:)=VQLM(LM)*R(:)**L
        AEPOT(:,LM)=AEPOT(:,LM)+pot(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
!       AEE(:)=AE(:)+AERHO(:,LM)*pot(:)
      ENDDO
!
!     ==================================================================
!     ==  ADD POTENTIAL FROM THE BACKGROUND                           ==
!     ==================================================================
      SVAR=-2.D0*PI*RHOB/3.D0*SQRT(4.D0*PI)
      POT(:)=SVAR*R(:)**2
      AEPOT(:,1)=AEPOT(:,1)+POT(:) ! POTENTIAL OF THE BACKGROUND
      AEE(:)=AEE(:)+(AERHO(:,1)+RHOC(:))*POT(:)
CALL RADIAL$INTEGRAL(GID,NR,(AERHO(:,1)+RHOC(:))*POT(:)*r**2,svar)
print*,'el background',svar
!
!     ==================================================================
!     ==  CALCULATE ENERGY                                            ==
!     ==================================================================
      AEE(:)=AEE(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AEE,AEeH)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_adjustvqlm(gid,NR,LMRX,psrho,RCSM,QLM,VQLM)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: gid
      INTEGER(4) ,INTENT(IN) :: NR
      REAL(8)    ,INTENT(IN) :: RCSM
      INTEGER(4) ,INTENT(IN) :: LMRX
      REAL(8)    ,INTENT(IN) :: PSRHO(NR,LMRX) ! INCLUDES PSEUDO CORE DENSITY
      REAL(8)    ,INTENT(IN) :: QLM(LMRX)
      REAL(8)    ,INTENT(OUT):: VQLM(LMRX)
      REAL(8)                :: r(NR)
      REAL(8)                :: AUX(NR)
      REAL(8)                :: RHO(NR)
      REAL(8)                :: g(NR)
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      INTEGER(4)             :: L,M,LM,IR
      INTEGER(4)             :: LX
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      LX=INT(SQRT(REAL(LMRX-1)+.01D0))
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
      subroutine augmentation_feedhyperfine(gid,nr,lmrx,ndimd &
     &         ,iat,aerho,aecore,psrho,pscore,aehpot,pshpot)
      implicit none
      integer(4),intent(in) :: gid
      integer(4),intent(in) :: nr
      integer(4),intent(in) :: lmrx
      integer(4),intent(in) :: ndimd
      integer(4),intent(in) :: iat
      real(8)   ,intent(in) :: aerho(nr,lmrx,ndimd)
      real(8)   ,intent(in) :: psrho(nr,lmrx,ndimd)
      real(8)   ,intent(in) :: aecore(nr)
      real(8)   ,intent(in) :: pscore(nr)
      real(8)   ,intent(in) :: aehpot(nr,lmrx)
      real(8)   ,intent(in) :: pshpot(nr,lmrx)
      real(8)               :: rho(nr,lmrx)
!     ****************************************************************
!
!     ================================================================
!     ==  total density                                             ==
!     ================================================================
      CALL HYPERFINE$SET1CRHO('PS','TOT',IAT,GID,NR,NR,LMRX,PSRHO)
      rho(:,:)=AERHO(:,:,1)
      rho(:,1)=rho(:,1)+AECORE(:)
      CALL HYPERFINE$SET1CRHO('AE','TOT',IAT,GID,NR,NR,LMRX,rho)
!
!     ================================================================
!     ==  spin density                                              ==
!     ================================================================
      IF(NDIMD.GT.1) THEN
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          rho(:,:)=psRHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          rho(:,:)=SQRT(PSRHO(:,:,2)**2+PSRHO(:,:,3)**2+PSRHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('PS','SPIN',IAT,GID,NR,NR,LMRX,rho)
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          rho(:,:)=AERHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          rho(:,:)=SQRT(AERHO(:,:,2)**2+AERHO(:,:,3)**2+AERHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('AE','SPIN',IAT,GID,NR,NR,LMRX,rho)
      END IF
!
!     ================================================================
!     ==  total potential for electric field gradients              ==
!     ================================================================
      CALL HYPERFINE$SET1CPOT('AE',IAT,gid,NR,NR,LMRX,AEhPOT)
      CALL HYPERFINE$SET1CPOT('PS',IAT,gid,NR,NR,LMRX,PShPOT)
      return
      end
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_EXPECT(gid,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &           ,AEPOT,PSPOT,AEPHI,PSPHI,DATH)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE EXPECTATION VALUE OF                         **
!     **  THE ONE-CENTER POTENTIALS WITH THE PARTIAL WAVES            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: gid
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
      REAL(8)               :: r(nr)
!     ******************************************************************
      call radial$r(gid,nr,r)
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
                CALL RADIAL$INTEGRAL(gid,NR,DWORK1,SVAR)
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
!!$      PI=4.D0*DATAN(1.D0)
!!$      C000=1/DSQRT(4.D0*PI)
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
!!$      XEXP=DEXP(DEX)
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
!.......................................................................
MODULE EXPERTNAL1CPOT_MODULE
!**                                                                   **
!**  APPLIES AN EXTERNAL POTENTIAL ACTING ON THE ELECTRONS            **
!**  WITHIN THE AUGMENTATION REGION                                   **
!**                                                                   **
!**    IDIMD EXTERNALLY THIS INDEX REFERS TO                          **
!**             0: TOTAL  FOR NDIMD=1,2,4                             **
!**             1: SZ FOR NDIMD=2                                     **
!**             1: SX FOR NDIMD=4                                     **
!**             2: SY FPR NDIMD=4                                     **
!**             3: SZ FPR NDIMD=4                                     **
!**          INTERNALLY THE VALUE IS RAISED BY ONE SO THAT IT         **
!**          CORRESPONDS TO THE INDEXING                              **
!**                                                                   **
TYPE EXTPOT
  REAL(8)          :: VALUE
  CHARACTER(LEN=32):: ATOM
  REAL(8)          :: RC  
  REAL(8)          :: PWR
  CHARACTER(LEN=32):: TYPE   ! CAN BE 'S','P','D','ALL'
  INTEGER(4)       :: IDIMD  ! IDIMD=0 REFERES TO ALL SPIN DIRECTIONS
END TYPE EXTPOT
INTEGER(4)             :: NPOT=0
INTEGER(4)             :: NPOTX=0
INTEGER(4),PARAMETER   :: DNPOT=5
TYPE (EXTPOT), ALLOCATABLE :: POT(:)
CONTAINS
!  .....................................................................
   SUBROUTINE CREATE
   TYPE (EXTPOT),ALLOCATABLE :: TMPPOT(:)
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
!     ..................................................................
      SUBROUTINE EXTERNAL1CPOT$SETPOT(ATOM,TYPE,IDIMD,VALUE,RC,PWR)
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      REAL(8)          ,INTENT(IN) :: VALUE
      CHARACTER(LEN=32),INTENT(IN) :: ATOM
      CHARACTER(LEN=32),INTENT(IN) :: TYPE
      INTEGER(4)       ,INTENT(IN) :: IDIMD
      REAL(8)          ,INTENT(IN) :: RC
      REAL(8)          ,INTENT(IN) :: PWR
!     ******************************************************************
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
!     ..................................................................
      SUBROUTINE EXTERNAL1CPOT$REPORT(NFIL)
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IPOT
!     ******************************************************************
      IF(NPOT.EQ.0) RETURN
!     CALL REPORT$TITLE(NFIL,"EXTERNAL POTENTIALS ON ORBITALS")
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
!       ................................................................
        SUBROUTINE WRITECH(NFIL,NAME,VALUE)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        CHARACTER(LEN=*),INTENT(IN) :: VALUE
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,A)')NAME,TRIM(VALUE)
        RETURN
        END SUBROUTINE WRITECH
!       ................................................................
        SUBROUTINE WRITER8(NFIL,NAME,VALUE,UNIT)
        INTEGER(4)       ,INTENT(IN) :: NFIL
        CHARACTER(LEN=*),INTENT(IN) :: NAME
        REAL(8)          ,INTENT(IN) :: VALUE
        CHARACTER(LEN=*),INTENT(IN) :: UNIT
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,F10.5," ",A)')NAME,VALUE,UNIT
        RETURN
        END SUBROUTINE WRITER8
!       ................................................................
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
!     ..................................................................
      SUBROUTINE EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATH,ETOT)
!     ******************************************************************
!     ******************************************************************
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: ATOM
      INTEGER(4)  ,INTENT(IN)   :: NDIMD
      INTEGER(4)  ,INTENT(IN)   :: LMNX
      REAL(8)     ,INTENT(IN)   :: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(OUT)  :: ETOT
      TYPE(EXTPOT)              :: POT1
      INTEGER(4)                :: LANG
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
      INTEGER(4)                :: LN1,LN2,LMN1,LMN2,IR,IPOT,IDIMD,L1,L2,I
      INTEGER(4)                :: gid     ! grid id
      real(8)      ,allocatable :: r(:)    ! radial grid
!     ******************************************************************
      DATH(:,:,:)=0.D0
      ETOT=0.D0
!
!     ==================================================================
!     ==  SELECT POTENTIAL                                            ==  
!     ==================================================================
!     == DETERMINE ID OF THE ATOM TYPE (SPECIES)
      CALL ATOMLIST$INDEX(ATOM,IAT)
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
!
      DO IPOT=1,NPOT
        POT1=POT(IPOT)
        TCHK=(ATOM.EQ.POT1%ATOM.OR.SPECIES.EQ.POT1%ATOM)
        IF(.NOT.TCHK) CYCLE
!
!       ==================================================================
!       ==  UNSCRAMBLE TYPE                                             ==  
!       ==================================================================
        IF(TRIM(POT1%TYPE).EQ.'S') THEN
          LANG=0
        ELSE IF(TRIM(POT1%TYPE).EQ.'P') THEN
          LANG=1
        ELSE IF(TRIM(POT1%TYPE).EQ.'D') THEN
          LANG=2
        ELSE IF(TRIM(POT1%TYPE).EQ.'ALL') THEN
          LANG=-1
        ELSE
          CALL ERROR$MSG('VALUE OF TYPE NOT RECOGNIZED')
          CALL ERROR$STOP('EXTERNAL1CPOT$APPLY')
        END IF
!
!       ==================================================================
!       == COLLECT PARTIAL WAVES                                        ==
!       ==================================================================
        CALL ATOMLIST$INDEX(ATOM,IAT)
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETCH('ID',SPECIES)
        CALL RADIAL$GETI4(GID,'NR',NR)
        allocate(r(nr))
        call radial$r(gid,nr,r)
        CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$LOFLN(ISP,LNX,LOX)
        ALLOCATE(AEPHI(NR,LNX))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
!
        ALLOCATE(RDEP(NR))
        ALLOCATE(AUX(NR))
        RDEP(:)=POT1%VALUE*EXP(-(R(:)/POT1%RC)**POT1%PWR)
!
        ALLOCATE(UONE(LNX,LNX))
        UONE(:,:)=0.D0
!
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
            AUX(:)=RDEP(:)*AEPHI(:,LN1)*AEPHI(:,LN2)*r(:)**2
            CALL RADIAL$INTEGRAL(gid,NR,AUX,UONE(LN1,LN2))
          ENDDO
        ENDDO
        DEALLOCATE(R)
        DEALLOCATE(RDEP)
        DEALLOCATE(AUX)
!OPEN(41,FILE='DUMP',FORM='FORMATTED')
!DO IR=1,NR
!  WRITE(41,*)R(ir),(AEPHI(IR,LN1),LN1=1,LNX)
!ENDDO
!PRINT*,'LNX ',LNX,' LOX ',LOX
!DO LN1=1,LNX
!  PRINT*,'UONE ',UONE(LN1,:)
!ENDDO
!STOP
        DEALLOCATE(AEPHI)
!
!       ==============================================================
!       ==  CALCULATE ONE-CENTER HAMILTONIAN                        ==
!       ==============================================================
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
                DO I=1,2*L1+1
                  DATH(LMN1+I,LMN2+I,IDIMD)=DATH(LMN1+I,LMN2+I,IDIMD) &
       &                                   +UONE(LN1,LN2)
                ENDDO
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
!     ==============================================================
!     ==  CALCULATE TOTAL ENERGY                                  ==
!     ==============================================================
      ETOT=0.D0
      DO IDIMD=1,NDIMD
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            ETOT=ETOT+DENMAT(LMN1,LMN2,IDIMD)*DATH(LMN1,LMN2,IDIMD)
          ENDDO
        ENDDO
      ENDDO
!     PRINT*,'ETOT ',ETOT
      RETURN 
      END
!=======================================================================
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
      REAL(8),   INTENT(IN) :: RC
      REAL(8),   INTENT(IN) :: R
      INTEGER(4),INTENT(IN) :: L
      REAL(8),   INTENT(OUT):: V
      REAL(8)               :: PI
      REAL(8)               :: FAC,X,SVAR1
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
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
      REAL(8)               :: PI
      REAL(8)               :: X2,GAUSS,FAC1,FAC2
      REAL(8)               :: ERFX
      INTEGER(4)            :: I
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      X2=X*X
      CALL LIB$ERFR8(X,ERFX)
      V=0.5D0*DSQRT(PI)*ERFX
      GAUSS=DEXP(-X2)
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
!!$            DIN=NNOFB(IB)+0.5D0+DATAN(DIN)/PI            
!!$            TFORWARD=.FALSE.
!!$            CALL SCHROEDER(TFORWARD,.FALSE.,R1,DEX,NR,LOFB(IB),EOFB(IB) &
!!$     &                    ,AEZ,IR,4.D0 &
!!$     &                     ,AEPOT,PHI,GINH,DLG)
!!$            DOUT=RCLOFB(IB)*PHI(IR,2)/PHI(IR,1)
!!$            DIN=0.5D0+DATAN(DOUT)/PI            
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
      SUBROUTINE AUGMENTATION$SPHEREold(ISP,IAT,LMNX,NDIMD,DENMAT,DENMATI &
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
      REAL(8)   ,INTENT(INOUT):: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(INOUT):: DENMATI(LMNX,LMNX,NDIMD)
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(INOUT):: VQLM(LMRX)
      REAL(8)   ,INTENT(IN)   :: RHOB
      REAL(8)   ,INTENT(OUT)  :: POTB ! INTEGRATED ELECTROSTATIC AUGMENTATION POTENTIAL
      REAL(8)   ,INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT)  :: DO(LMNX,LMNX,NDIMD)
      integer(4)              :: gid
      REAL(8)   ,ALLOCATABLE  :: r(:)
      REAL(8)                 :: R1,DEX,XEXP,RI
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
      REAL(8)   ,ALLOCATABLE  :: DTKIN(:,:)
      REAL(8)   ,ALLOCATABLE  :: AERHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: AEPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: PSPOT(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: AEPOT1(:,:)
      REAL(8)   ,ALLOCATABLE  :: PSPOT1(:,:)
      REAL(8)   ,ALLOCATABLE  :: DATP(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DWORK1(:)
      INTEGER(4)              :: IDIM,LMR,IR
      INTEGER(4)              :: LM,LMN1,LMN2,LMN11,LMN21,LN1,LN2,L1,L2,IM,ISPIN
      INTEGER(4)              :: NSPIN
      INTEGER(4)              :: NFILO
      LOGICAL(4),PARAMETER    :: TPR=.false.
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=1
      LOGICAL(4)              :: TBACK,TSPIN
      REAL(8)                 :: DETOT,PSEHARTREE,AEEHARTREE,COREEXC,EKINNL,ENL,AEEXC,PSEXC,HAMUP,HAMDWN
      REAL(8)                 :: AEBACKGROUND,PSBACKGROUND
      CHARACTER(32)           :: ATOM
      REAL(8)                 :: VQLM1(LMRX)
      REAL(8)                 :: QLM(LMRX)
      REAL(8)                 :: PI
      logical(4),parameter    :: tsoftcore=.true.
      real(8)                 :: aetest,pstest
!     ******************************************************************
                            CALL TRACE$PUSH('AUGMENTATION$SPHERE')
      PI=4.D0*DATAN(1.D0)
print*,'old sphere'
!
!     ==================================================================
!     ==  COLLECT ATOM-TYPE SPECIFIC INFORMATION FROM SETUP OBJECT    ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      allocate(r(nr))
      call radial$r(gid,nr,r)
!CALL SETUP$RADGRID(ISP,R1,DEX,NR)
!XEXP=DEXP(DEX)
!
      CALL SETUP$AEZ(ISP,AEZ)
      CALL SETUP$LNX(ISP,LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$LOFLN(ISP,LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUP$PSPARTIALWAVES(ISP,NR,LNX,PSPHI)
      ALLOCATE(AECORE(NR))
      CALL SETUP$AECORE(ISP,NR,AECORE)
      ALLOCATE(PSCORE(NR))
      CALL SETUP$PSCORE(ISP,NR,PSCORE)
      CALL SETUP$RCSM(ISP,RCSM)
      ALLOCATE(VADD(NR))
      CALL SETUP$VBAR(ISP,NR,VADD)
      ALLOCATE(DOVER(LNX,LNX))
      CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
      ALLOCATE(DTKIN(LNX,LNX))
      CALL SETUP$1CKINETIC(ISP,LNX,DTKIN)
!
!     ==================================================================
!     ==  SELF TEST                                                   ==
!     ==================================================================
 1000 CONTINUE
      IF(TTEST) THEN
        IF(ITEST.EQ.1.OR.ITEST.EQ.3) THEN 
          CALL SELFTEST$START('SPHERE',LMNX*LMNX*NDIMD,DENMAT,1.D-3)
          VQLM(:)=0.D0
          DO IDIM=1,NDIMD
            DENMAT(:,:,IDIM)=0.5D0*(DENMAT(:,:,IDIM)+TRANSPOSE(DENMAT(:,:,IDIM)))
          ENDDO
        END IF
      ENDIF
!     == PRINT DENSITY MATRIX FOR TEST
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,FMT='("DENMAT FOR IDIM= ",I2)') IDIM
          DO LMN1=1,LMNX
            WRITE(NFILO,FMT='(9F15.6)')DENMAT(LMN1,:,IDIM)
          ENDDO
        ENDDO
      END IF

!FOR SOFT-CORE ITERATE FROM HERE ...
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
      CALL AUGMENTATION_QLM(gid,NR,LMRX,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
      psrho(:,1,1)=psrho(:,1,1)+pscore(:)   
      CALL AUGMENTATION_adjustvqlm(gid,NR,LMRX,PSRHO(:,:,1),RCSM,QLM,VQLM1)
      psrho(:,1,1)=psrho(:,1,1)-pscore(:)   
print*,'vqlm1 ',vqlm+vqlm1
!     
!     ================================================================
!     ==  new soft core                                             ==
!     ================================================================
      if(tsoftcore) then
!        call augmentation_newsoftcore(gid,nr,aez,aerho(:,1,1),aecore)
      end if

!!$!test
!!$      allocate(pspot1(nr,lmrx))
!!$      allocate(aepot1(nr,lmrx))
!!$      call AUGMENTATION_PSHARTREE(GID,NR,LMRX,PScore,PSRHO(:,:,1) &
!!$     &                 ,VADD,RCSM,QLM,VQLM,RHOB,PSPOT1,PSEHartree)
!!$      call AUGMENTATION_AEhartree(GID,NR,LMRX,aecore,AERHO(:,:,1) &
!!$     &                 ,VQLM,RHOB,aepot1,AEEHartree)
!!$      deallocate(aepot1)
!!$      deallocate(pspot1)
!!$print*,'ehartree ',aeehartree,psehartree
!!$!
!!$!     =================================================================
!!$!     == AVERAGE electrostatic ONE-CENTER POTENTIAL                  ==
!!$!     =================================================================
!!$      ALLOCATE(DWORK1(NR))
!!$      DWORK1(:)=AEPOT(:,1,1)-(PSPOT(:,1,1)-VADD(:))
!!$      DWORK1(:)=DWORK1(:)*R(:)**2*SQRT(4.D0*PI)  !SQRT(4*PI)=4*PI*Y_0
!!$      CALL RADIAL$INTEGRAL(gid,NR,DWORK1,POTB)
!!$      POTB=-POTB
!!$!PEB03 POTB=0.D0
!!$      DEALLOCATE(DWORK1)
!     
!     ================================================================
!     ==  SEND DENSITIES TO HYPERFINE-PARAMETER OBJECT              ==
!     ================================================================
!MOVE THIS PART TO THE END OF THE ROUTINE
!TAKE CARE HOWEVER THAT THE FIRST CALL TAKES THE PSEUDODENSITY
! WITHOUT PSEUDO CORE!
      ALLOCATE(AEPOT(NR,LMRX,1))   ! USED AS AUXILIARY VARIABLE
      CALL HYPERFINE$SET1CRHO('PS','TOT',IAT,gid,NR,NR,LMRX,PSRHO)
      AEPOT(:,:,1)=AERHO(:,:,1)
      AEPOT(:,1,1)=AEPOT(:,1,1)+AECORE(:)
      CALL HYPERFINE$SET1CRHO('AE','TOT',IAT,gid,NR,NR,LMRX,AEPOT)
      IF(NDIMD.GT.1) THEN
        ALLOCATE(PSPOT(NR,LMRX,1))
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          AEPOT(:,:,1)=AERHO(:,:,2)
          PSPOT(:,:,1)=PSRHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          AEPOT(:,:,1)=SQRT(AERHO(:,:,2)**2+AERHO(:,:,3)**2+AERHO(:,:,4)**2)
          PSPOT(:,:,1)=SQRT(PSRHO(:,:,2)**2+PSRHO(:,:,3)**2+PSRHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('PS','SPIN',IAT,gid,NR,NR,LMRX,PSPOT)
        CALL HYPERFINE$SET1CRHO('AE','SPIN',IAT,gid,NR,NR,LMRX,AEPOT)
        DEALLOCATE(PSPOT)
      END IF
      DEALLOCATE(AEPOT)
!     
!     ================================================================
!     ==  ADD CORE CHARGE DENSITY                                   ==
!     ================================================================
      DO IR=1,NR
        PSRHO(IR,1,1)=PSRHO(IR,1,1)+PSCORE(IR)
      ENDDO
!     
!     ================================================================
!     ================================================================
!     ==   CALCULATE 1-CENTER POTENTIAL                             ==
!     ================================================================
!     ================================================================
      ALLOCATE(AEPOT(NR,LMRX,NDIMD))
      ALLOCATE(PSPOT(NR,LMRX,NDIMD))
      AEPOT(:,:,:)=0.D0
      PSPOT(:,:,:)=0.D0
!     
!     ================================================================
!     ==   ADD EXCHANGE AND CORRELATION POTENTIAL                   ==
!     ================================================================
      AEEXC=0.D0
      PSEXC=0.D0
!     == AE-EXCHANGE ENERGY AND POTENTIAL ============================
      AERHO(:,1,1)=AERHO(:,1,1)+AECORE(:)
      CALL AUGMENTATION_XC(gid,NR,LMRX,NDIMD,AERHO,AEEXC,AEPOT)
      AERHO(:,1,1)=AERHO(:,1,1)-AECORE(:)
!     == PS-EXCHANGE ENERGY AND POTENTIAL ============================
      CALL AUGMENTATION_XC(gid,NR,LMRX,NDIMD,PSRHO,PSEXC,PSPOT)
!     == CORE ONLY EXCHANGE ENERGY ===================================
      COREEXC=0.D0
      ALLOCATE(DWORK1(NR))
      CALL AUGMENTATION_XC(gid,NR,1,1,AECORE,COREEXC,DWORK1)
      DEALLOCATE(DWORK1)
      AEEXC=AEEXC-COREEXC
!     
      CALL AUGMENTATION_ADD('AE1 EXCHANGE-CORRELATION',AEEXC)
      CALL AUGMENTATION_ADD('PS1 EXCHANGE-CORRELATION',PSEXC)
!     
!     ================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL                        ==
!     ================================================================
      AEEHARTREE=0.D0
      PSEHARTREE=0.D0
      ALLOCATE(AEPOT1(NR,LMRX))
      ALLOCATE(PSPOT1(NR,LMRX))
      AEPOT1(:,:)=0.D0
      PSPOT1(:,:)=0.D0
      VQLM1(:)=0.D0
!     == NOTE: COMPENSATION DENSITY IS ADDED TO PSRHO ================
      CALL AUGMENTATION_HARTREE(IAT,gid,NR,AEZ,AECORE &
     &          ,VADD,RCSM,LMRX,AERHO,PSRHO,QLM,AEPOT1,PSPOT1,VQLM1 &
     &          ,AEEHARTREE,PSEHARTREE)
      CALL AUGMENTATION_ADD('AE1 ELECTROSTATIC',AEEHARTREE)
      CALL AUGMENTATION_ADD('PS1 ELECTROSTATIC',PSEHARTREE)
print*,'vqlm1 ',vqlm+vqlm1
      VQLM1(:)=VQLM1(:)+VQLM(:)
aetest=aeehartree
pstest=psehartree
print*,'aeehartree without background etc. ',aeehartree
!
!     =================================================================
!     == AVERAGE ONE-CENTER POTENTIAL                                ==
!     =================================================================
      ALLOCATE(DWORK1(NR))
      DWORK1(:)=AEPOT1(:,1)-(PSPOT1(:,1)-VADD(:))
      DWORK1(:)=DWORK1(:)*R(:)**2*SQRT(4.D0*PI)  !SQRT(4*PI)=4*PI*Y_0
      CALL RADIAL$INTEGRAL(gid,NR,DWORK1,POTB)
      POTB=-POTB
!PEB03 POTB=0.D0
      DEALLOCATE(DWORK1)
!     
!     == ADD POTENTIAL FROM CHARGES "OUTSIDE" THE SPHERE =============      
      CALL AUGMENTATION_ADDVQLM(gid,NR,LMRX,VQLM1,AEPOT1,PSPOT1)
!     
!     == ADD POTENTIAL FROM THE NEUTRALIZING BACKGROUND ==============
!     == THE BACKGROUND DOES NOT ADD TO SUM OF ONE-CENTER ENERGIES  ==
!     == AND NEED NOT BE CONSIDERED. IT SHOULD BE INCLUDED, IF THE  ==
!     == POTENTIALS SHOULD BE THE CORRECT ONE-CENTER POTENTIALS     ==
!     == AT THIS POIINT THE DERIVATIVES ARE NOT CORRECT YET         ==
      AEBACKGROUND=0.D0
      PSBACKGROUND=0.D0
print*,'background density ',rhob
      CALL AUGMENTATION_ADDBACKGROUND(gid,NR,RHOB &
     &                      ,AERHO(:,1,1)+AECORE(:),AEBACKGROUND,AEPOT1)
print*,'background energy ',AEBACKGROUND
      CALL AUGMENTATION_ADDBACKGROUND(gid,NR,RHOB &
     &                      ,PSRHO,PSBACKGROUND,PSPOT1)
      CALL AUGMENTATION_ADD('AE1 BACKGROUND',AEBACKGROUND)
      CALL AUGMENTATION_ADD('PS1 BACKGROUND',PSBACKGROUND)
aetest=aetest+aebackground
pstest=pstest+psbackground
print*,'ehartree 2 ',aetest,pstest
!
!     == SOFT CORE ===================================================
!      CALL AUGMENTATION_SOFTCORE(IAT,R1,DEX,NR,LMRX,AEPOT1)
!     
!     == ANALYSIS: ELECTRIC FIELD GRADIENTS ==========================
      CALL HYPERFINE$SET1CPOT('AE',IAT,gid,NR,NR,LMRX,AEPOT1)
      CALL HYPERFINE$SET1CPOT('PS',IAT,gid,NR,NR,LMRX,PSPOT1)
!     
!     == ANALYSIS: POTENTIAL PLOT           ==========================
      CALL GRAPHICS$SET1CPOT('AE',IAT,gid,NR,NR,LMRX,AEPOT1)
      CALL GRAPHICS$SET1CPOT('PS',IAT,gid,NR,NR,LMRX,PSPOT1)
!     
!     == ADD ELECTROSTATIC POTENTIAL TO TOTAL POTENTIAL ==============
      DO LM=1,LMRX
        DO IR=1,NR
          AEPOT(IR,LM,1)=AEPOT(IR,LM,1)+AEPOT1(IR,LM)
          PSPOT(IR,LM,1)=PSPOT(IR,LM,1)+PSPOT1(IR,LM)
        ENDDO
      ENDDO
      DEALLOCATE(AEPOT1)
      DEALLOCATE(PSPOT1)
!
!     ==== THIS IS FOR MARCELLO SANTOS TO CALCULATE CORE LEVEL SHIFTS IN FROZEN CORE
      CALL CORE_CORESHIFTS(IAT,ISP,gid,NR,LMRX,AEPOT)
!
!CALCULATE NEW SOFT-CORE DENSITY HERE
!     
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO ISPIN=1,NSPIN
          WRITE(NFILO,*)'AE POTENTIAL FOR SPIN ',ISPIN
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(AEPOT(IR,LM,ISPIN),LM=1,LMRX)
          ENDDO
        ENDDO
        DO ISPIN=1,NSPIN
          WRITE(NFILO,*)'PS POTENTIAL FOR ATOM ',IAT,ISPIN
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(PSPOT(IR,LM,ISPIN),LM=1,LMRX)
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
        ELSE IF(ITEST.EQ.3) THEN
          call radial$r(gid,'dex',dex,r)
          ENL=0.D0
          DO LM=1,LMRX
            enl=enl+dex*sum(AERHO(:,LM,1)**2*r(:)**3)
            AEPOT(:,LM,1)=2.D0*AERHO(:,LM,1)
            PSPOT(:,LM,1)=0.D0
          ENDDO
          AEEHARTREE=ENL
          PSEHARTREE=0.D0
          AEEXC=0.D0
          PSEXC=0.D0
        END IF
      END IF
!     
!     ==================================================================
!     ==  EVALUATE KINETIC HAMILTONIAN AND OVERLAP                    ==
!     ==================================================================
      EKINNL=0.D0
      DATH(:,:,:)=0.D0
      DO(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            DO IM=1,2*L1+1
              LMN11=LMN1+IM
              LMN21=LMN2+IM
              EKINNL=EKINNL+DENMAT(LMN11,LMN21,1)*DTKIN(LN1,LN2)
              DATH(LMN11,LMN21,1)=DTKIN(LN1,LN2)
              DO(LMN11,LMN21,1)  =DOVER(LN1,LN2)
            ENDDO
          ENDIF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO
      IF(TTEST) THEN
        IF(ITEST.EQ.2) THEN
          CALL SELFTEST$END('SPHERE-EKI',LMNX*LMNX*NDIMD,DATH,EKINNL,TBACK)
          IF(TBACK) GOTO 1001
          CALL ERROR$MSG('STOP AFTER SELFTEST')
          CALL ERROR$STOP('SPHERE')
        END IF
      END IF
!     
      CALL AUGMENTATION_ADD('AE1-PS1 KINETIC',EKINNL)
!     
!     ================================================================
!     ==   ADD POTENTIAL ENERGY TO THE ONE-CENTER HAMILTONIAN       ==
!     ==   DATH =DATH + <AEPHI|AEPOT|AEPHI>-<PSPHI|PSPOT|PSPHI>     ==
!     ================================================================
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL AUGMENTATION_EXPECT(gid,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                        ,AEPOT,PSPOT,AEPHI,PSPHI,DATP)
      DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
      DEALLOCATE(DATP)
!     WRITE(TESTSTRING,FMT='("R8DATH",I2,12(" "))')IAT
!     CALL STOREIT(TESTSTRING,8*LMNX*LMNX*NSPIN,DATH)
!
!     ================================================================
!     ==  LDA + U                                                   ==
!     ================================================================
!     ALLOCATE(DATP(LMNX,LMNX,NDIMD))
!     DETOT=0.D0
!     CALL LDAPLUSU(NRX,LNX,LMNX,NSPIN,LOX,LNX &
!     &             ,GID,NR,AEZ,AEPHI,DENMAT,DETOT,DATP)
!     DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
!     DEALLOCATE(DATP)
!     CALL AUGMENTATION_ADD('LDA+U EXCHANGE',DETOT)
!
!     ================================================================
!     ==  APPLY EXTERNAL POTENTIAL                                  ==
!     ================================================================
      CALL ATOMLIST$GETCH('NAME',IAT,ATOM)
      ALLOCATE(DATP(LMNX,LMNX,NDIMD))
      CALL EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATP,DETOT)
      DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
      DEALLOCATE(DATP)
      CALL AUGMENTATION_ADD('LDA+U EXCHANGE',DETOT)
!     
!     ================================================================
!     ==  SELF-TEST                                                 ==
!     ================================================================
      IF(TTEST) THEN
        IF(ITEST.EQ.1.OR.ITEST.EQ.3) THEN
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
            WRITE(NFILO,FMT='(9F15.6)')DATH(LMN1,:,IDIM)
          ENDDO
        ENDDO
STOP
      END IF
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(AEPOT)
      DEALLOCATE(PSPOT)
      DEALLOCATE(AERHO)
      DEALLOCATE(PSRHO)
      DEALLOCATE(LOX)
      DEALLOCATE(AEPHI)
      DEALLOCATE(PSPHI)
      DEALLOCATE(AECORE)
      DEALLOCATE(PSCORE)
      DEALLOCATE(VADD)
      DEALLOCATE(DOVER)
      DEALLOCATE(DTKIN)
      DEALLOCATE(r)
                               CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADDVQLM(gid,NR,LMRX,VQLM,AEPOT,PSPOT)
!     ******************************************************************
!     **                                                             **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: gid
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
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        SVAR=VQLM(LM)
        AEPOT(:,LM)=AEPOT(:,LM)+SVAR*R(:)**L
        PSPOT(:,LM)=PSPOT(:,LM)+SVAR*R(:)**L
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADDBACKGROUND(gid,NR,RHOB &
     &                        ,RHO,EB,POT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ONE-CENTER CONTRIBUTION OF THE COMPENSATING  **
!     **  CHARGE BACKGROUND TO TOTAL ENERGY AND POTENTIAL             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: gid
      INTEGER(4),INTENT(IN)   :: NR
      REAL(8)   ,INTENT(IN)   :: RHOB
      REAL(8)   ,INTENT(IN)   :: RHO(NR)
      REAL(8)   ,INTENT(OUT)  :: EB
      REAL(8)   ,INTENT(INOUT):: POT(NR)
      REAL(8)                 :: PI
      REAL(8)                 :: SVAR
      REAL(8)                 :: WORK(NR)
      real(8)                 :: r(nr)
!     ******************************************************************
      IF(RHOB.EQ.0.D0) THEN
        EB=0.D0
        RETURN
      END IF
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      SVAR=-2.D0*PI*RHOB/3.D0*DSQRT(4.D0*PI)
      WORK(:)=SVAR*RHO(:)*r(:)**4
      POT(:)=POT(:)+SVAR*r(:)**2 ! POTENTIAL OF THE BACKGROUND
      CALL RADIAL$INTEGRAL(gid,NR,WORK,EB)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_HARTREE(IAT,gid,NR,AEZ,RHOCOR &
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
      INTEGER(4) ,INTENT(IN) :: gid
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
      real(8)                :: r(nr)
      REAL(8)                :: AUX1(NR)
      REAL(8)                :: AEDMU(NR)
      REAL(8)                :: PSDMU(NR)
      REAL(8)                :: AEE(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: PI
      REAL(8)                :: Y0
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      REAL(8)                :: GAUSSIAN(nr)
      REAL(8)                :: AEH,PSH
      INTEGER(4)             :: LM,IR
      INTEGER(4)             :: L
!     == ARRAYS NEEDED DETAILED REPORT
      LOGICAL(4) ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4)             :: NFILO
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      call radial$r(gid,nr,r)
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
      CALL RADIAL$POISSON(gid,NR,0,RHOCOR,AEDMU)
      CALL SETUP$GETR8A('NUCPOT',NR,AUX1)
      AEDMU(:)  = AEDMU(:)+aux1(:)
!      AEDMU(:)  = AEDMU(:)-AEZ/R(:)/Y0
      AEE(:)    = AEDMU(:)* AERHO(:,1) *R(:)**2
      AEPOT(:,1)= AEPOT(:,1)+AEDMU(:)
      PSE(:)    = VADD(:) * PSRHO(:,1) *R(:)**2
      PSPOT(:,1)= PSPOT(:,1)+VADD(:)
      CALL RADIAL$INTEGRAL(gid,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(gid,NR,PSE,PSH)
print*,'aeh c-v int',aeh
      AEEHARTREE=AEEHARTREE+AEH
      PSEHARTREE=PSEHARTREE+PSH
!
!     ==================================================================
!     ==   ADD COMPENSATION CHARGE DENSITY                            ==
!     ==================================================================
      ALPHA=1.D0/RCSM**2
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        CALL GAUSSN(L,ALPHA,CL)
        SVAR=QLM(LM)*CL
        PSRHO(:,LM)=PSRHO(:,LM)+SVAR*DEXP(-ALPHA*r(:)**2)*R(:)**L
      ENDDO
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      AEE(:)=0.D0
      PSE(:)=0.D0
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        CALL RADIAL$POISSON(gid,NR,L,AERHO(1,LM),AEDMU)
        CALL RADIAL$POISSON(gid,NR,L,PSRHO(1,LM),PSDMU)
        ALPHA=1.D0/RCSM**2
        CALL GAUSSN(L,ALPHA,CL)
!       == ADD COMPENSATION DENSITY AND ITS POTENTIAL ================
        GAUSSIAN(:)=CL*DEXP(-ALPHA*R(:)**2)*R(:)**L
!
!       == CALCULATE TOTAL ENERGY AND POTENTIAL ======================
        AEE(:)=AEE(:)+0.5D0*AEdmu(:)*AERHO(:,lm)*R(:)**2
        PSE(:)=PSE(:)+0.5D0*PSdmu(:)*PSRHO(:,lm)*R(:)**2
        AEPOT(:,LM)=AEPOT(:,LM)+AEdmu(:)
        PSPOT(:,LM)=PSPOT(:,LM)+PSdmu(:)
        AUX1(:)=PSdmu(:)*GAUSSIAN*R(:)**2
        CALL RADIAL$INTEGRAL(gid,NR,AUX1,SVAR)
        VQLM(LM)=VQLM(LM)-SVAR
      ENDDO
      CALL RADIAL$INTEGRAL(gid,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(gid,NR,PSE,PSH)
      AEEHARTREE=AEEHARTREE+AEH
      PSEHARTREE=PSEHARTREE+PSH
print*,'aeh v-v int',aeh
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
