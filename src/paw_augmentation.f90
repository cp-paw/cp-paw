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
INTEGER(4)  ,PARAMETER :: NE=9
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
!     ==  ID(9)='CORE RELAXATION' ======================================
      SVAR=VAL(7)+VAL(1)-VAL(2)+VAL(3)-VAL(4)+VAL(5)-VAL(6)+VAL(8)+VAL(9)
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
      complex(8),INTENT(IN)  :: DENMAT(LMNX,LMNX)
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
      CALL AUGMENTATION_QLM(GID,NR,LMRX,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
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
! NEW SPHERE
      SUBROUTINE AUGMENTATION$SPHERE(ISP,IAT,LMNX,NDIMD,DENMAT,eDENMAT &
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
      complex(8),INTENT(INOUT):: DENMAT(LMNX,LMNX,NDIMD)
      complex(8),INTENT(INOUT):: eDENMAT(LMNX,LMNX,NDIMD)
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(INOUT):: VQLM(LMRX)
      REAL(8)   ,INTENT(IN)   :: RHOB ! NEUTRALIZING BACKGROUND DENSITY
!                               RHOB MAY BE SET TO ZERO BY ISOLATE OBJECT
      REAL(8)   ,INTENT(OUT)  :: POTB ! NEG. AV. EL. AUGM. POT.
      REAL(8)   ,INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT)  :: DO(LMNX,LMNX,NDIMD)
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
      INTEGER(4)              :: IDIM,LMR,IR
      INTEGER(4)              :: LM,LMN1,LMN2,LMN11,LMN21,LN1,LN2,L1,L2,IM,ISPIN
      INTEGER(4)              :: NSPIN
      INTEGER(4)              :: NFILO
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=1
      LOGICAL(4)              :: TBACK,TSPIN
      REAL(8)                 :: DETOT,PSEHARTREE,AEEHARTREE,COREEXC
      REAL(8)                 :: EKINNL,ENL,AEEXC,PSEXC,HAMUP,HAMDWN
      CHARACTER(32)           :: ATOM
      REAL(8)                 :: VQLM1(LMRX)
      REAL(8)                 :: QLM(LMRX)
      REAL(8)                 :: PI
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
      PI=4.D0*DATAN(1.D0)
!
!     ==================================================================
!     ==  COLLECT ATOM-TYPE SPECIFIC INFORMATION FROM SETUP OBJECT    ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
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
!
!     ==================================================================
!     ==  EXPAND OVERLAP AND KINETIC ENERGY MATRIX ELEMENTS           ==
!     ==================================================================
      ALLOCATE(DOVER(LNX,LNX))
      CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DOVER,DO)
      DEALLOCATE(DOVER)
!
      ALLOCATE(DTKIN(LMNX,LMNX,NDIMD))
      ALLOCATE(DTKIN1(LNX,LNX))
      CALL SETUP$1CKINETIC(ISP,LNX,DTKIN1)
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,DTKIN1,DTKIN)
      DEALLOCATE(DTKIN1)
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
            DENMAT(:,:,IDIM)=0.5D0*(DENMAT(:,:,IDIM)+TRANSPOSE(conjg(DENMAT(:,:,IDIM))))
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
!     ==  CALCULATE 1-CENTER kinetic energy                           ==
!     ==================================================================
!     == density matrix is hermitean for each spin direction
      EKINNL=real(sum(conjg(denmat)*dtkin),kind=8)
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
      CALL SETUP$GETCH('SOFTCORETYPE',SOFTCORETYPE)
      IF(SOFTCORETYPE.EQ.'NONE') THEN
        DECORE=0.D0
      ELSE 
        ALLOCATE(DELTARHO(NR,LMRX,NDIMD))
        CALL AUGMENTATION_NEWSOFTCORE(SOFTCORETYPE,GID,NR,LMRX,NDIMD,AEZ,LMNX &
     &          ,DENMAT,edenmat,VQLM1,RHOB,DELTAH,DELTAO,DELTARHO,DECORE)
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
!     ================================================================
!     ==  SEND DENSITIES TO HYPERFINE-PARAMETER OBJECT              ==
!     ================================================================
      CALL AUGMENTATION_FEEDHYPERFINE(GID,NR,LMRX,NDIMD &
     &         ,IAT,AERHO,AECORE,PSRHO,PSCORE,AEHPOT,PSHPOT)
!     
!     ================================================================
!     ==  SEND DENSITIES AND POTENTIALS TO GRAPHICS OBJECT          ==
!     ================================================================
      CALL GRAPHICS$SET1CPOT('AE',IAT,GID,NR,NR,LMRX,AEHPOT)
      CALL GRAPHICS$SET1CPOT('PS',IAT,GID,NR,NR,LMRX,PSHPOT)
!     
!     ================================================================
!     ==  EVALUATE CORE SHIFTS                                      ==
!     ================================================================
      CALL CORE_CORESHIFTS(IAT,ISP,GID,NR,LMRX,AETOTPOT)
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
            WRITE(NFILO,FMT='(9F10.5)')(AETOTPOT(IR,LM,ISPIN),LM=1,LMRX)
          ENDDO
        ENDDO
        DO ISPIN=1,NSPIN
          WRITE(NFILO,*)'PS POTENTIAL FOR ATOM ',IAT,ISPIN
          DO IR=1,NR,50
            WRITE(NFILO,FMT='(9F10.5)')(PSTOTPOT(IR,LM,ISPIN),LM=1,LMRX)
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
      DATH(:,:,:)=DTKIN(:,:,:)
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
!      DEALLOCATE(DOVER)
      DEALLOCATE(DTKIN)
      DEALLOCATE(R)
                               CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_NEWSOFTCORE(SOFTCORETYPE,GID,NR,LMRX,NDIMD,AEZ &
     &             ,LMNX,DENMAT,edenmat,VQLM,RHOB,DELTAH,DELTAO,DELTARHO,DECORE)        
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: SOFTCORETYPE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: AEZ
      INTEGER(4),INTENT(IN) :: LMNX
      complex(8),INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD)
      complex(8),INTENT(IN) :: eDENMAT(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)   ,INTENT(IN) :: RHOB
      REAL(8)   ,INTENT(OUT):: DELTAH(LMNX,LMNX,NDIMD) 
      REAL(8)   ,INTENT(OUT):: DELTAO(LMNX,LMNX,NDIMD)   
      REAL(8)   ,INTENT(OUT):: DELTARHO(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: DECORE
      LOGICAL(4)            :: TSPHERICAL=.FALSE.
      LOGICAL(4)            :: TSOFT=.TRUE.
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:)     ! AE PARTIAL WAVES FROM FILE
      REAL(8)   ,ALLOCATABLE:: UPHI(:,:)      ! NODELESS PARTIAL WAVES
      REAL(8)   ,ALLOCATABLE:: TUPHI(:,:)

      REAL(8)               :: AEPOT_AT(NR,LMRX,NDIMD)
      REAL(8)   ,ALLOCATABLE:: AEPHI_AT(:,:)  ! AE PHI FROM UPHI AND CORE FROM FILE
      REAL(8)   ,ALLOCATABLE:: TAEPHI_AT(:,:)
      REAL(8)   ,ALLOCATABLE:: PHIC_AT(:,:)
      REAL(8)   ,ALLOCATABLE:: TPHIC_AT(:,:)
!
      REAL(8)               :: RHOC_AT(NR)
      REAL(8)               :: RHOV_AT(NR,LMRX,NDIMD)
      REAL(8)               :: TKIN_AT(LMNX,LMNX,NDIMD)
      REAL(8)               :: H_AT(LMNX,LMNX,NDIMD)
      REAL(8)               :: O_AT(LMNX,LMNX,NDIMD)
      REAL(8)               :: ECORE_AT,EVALENCE_AT
!
      REAL(8)               :: AEPOT_SPH(NR,LMRX,NDIMD)
      REAL(8)               :: POTIN(NR)
      REAL(8)   ,ALLOCATABLE:: AEPHI_SPH(:,:) ! AE PART. WAVES FROM UPHI AND SPH. SOFTCORE
      REAL(8)   ,ALLOCATABLE:: TAEPHI_SPH(:,:)
      REAL(8)   ,ALLOCATABLE:: PHIC_SPH(:,:)
      REAL(8)   ,ALLOCATABLE:: TPHIC_SPH(:,:)
      REAL(8)               :: RHOC_SPH(NR)
      REAL(8)               :: RHOV_SPH(NR,LMRX,NDIMD)
      REAL(8)               :: TKIN_SPH(LMNX,LMNX,NDIMD)
      REAL(8)               :: H_SPH(LMNX,LMNX,NDIMD)
      REAL(8)               :: O_SPH(LMNX,LMNX,NDIMD)
      REAL(8)               :: ECORE_SPH,EVALENCE_SPH
!
      REAL(8)               :: AEPOT_NSPH(NR,LMRX,NDIMD)
      REAL(8)               :: POTIN_NSPH(NR,LMRX,NDIMD)
      COMPLEX(8),ALLOCATABLE:: AEPHI_NSPH(:,:,:,:) ! AE PART. WAVES FROM UPHI AND NSPH. SOFTCORE
      COMPLEX(8),ALLOCATABLE:: TAEPHI_NSPH(:,:,:,:) ! AE PART. WAVES FROM UPHI AND NSPH. SOFTCORE
      COMPLEX(8),ALLOCATABLE:: sAEPHI_NSPH(:,:,:,:) ! AE PART. WAVES FROM UPHI AND NSPH. SOFTCORE
      COMPLEX(8),ALLOCATABLE:: TsAEPHI_NSPH(:,:,:,:) ! AE PART. WAVES FROM UPHI AND NSPH. SOFTCORE
      INTEGER(4)            :: NC_NSPH
      REAL(8)   ,ALLOCATABLE:: EC_NSPH(:)
      COMPLEX(8),ALLOCATABLE:: PHIC_NSPH(:,:,:,:)
      COMPLEX(8),ALLOCATABLE:: TPHIC_NSPH(:,:,:,:)
      COMPLEX(8),ALLOCATABLE:: SPHIC_NSPH(:,:,:,:)
      COMPLEX(8),ALLOCATABLE:: TSPHIC_NSPH(:,:,:,:)
      REAL(8)               :: RHOC_NSPH(NR,LMRX,NDIMD)
      REAL(8)               :: RHOV_NSPH(NR,LMRX,NDIMD)
      REAL(8)               :: TKIN_NSPH(LMNX,LMNX,NDIMD)
      REAL(8)               :: H_NSPH(LMNX,LMNX,NDIMD)
      REAL(8)               :: O_NSPH(LMNX,LMNX,NDIMD)
      REAL(8)               :: ECORE_NSPH,EVALENCE_NSPH
   
      REAL(8)               :: POTNS(NR,LMRX,NDIMD)
      REAL(8)               :: POTNSOUT(NR,LMRX,NDIMD)
      REAL(8)               :: RHONS(NR,LMRX,NDIMD)

      REAL(8)               :: BROYDENWEIGHT(NR,LMRX,NDIMD)
      REAL(8)               :: AEHPOT(NR,LMRX)
      REAL(8)               :: AEXCPOT(NR,LMRX,NDIMD)
      REAL(8)               :: PI,Y0,C0LL
      REAL(8)               :: R(NR)
      LOGICAL(4)            :: CONVG
      REAL(8)               :: RHO(NR)
      REAL(8)               :: POT(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: EH,EXC,AEEHARTREE
      REAL(8)   ,ALLOCATABLE:: FOFI(:)      
      REAL(8)   ,ALLOCATABLE:: EOFI(:)      
      INTEGER(4),ALLOCATABLE:: LOFI(:)      
      INTEGER(4),ALLOCATABLE:: SOFI(:)      
      INTEGER(4),ALLOCATABLE:: NNOFI(:)      
      INTEGER(4)            :: LNX
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)               :: XMAX
      REAL(8)               :: SVAR,SVAR1,SVAR2
      INTEGER(4)            :: ITER,I,J,IB
      INTEGER(4),PARAMETER  :: NITER=2000
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      INTEGER(4)            :: NC
      REAL(8)               :: EKINC,EKINCAT,EKINV
      INTEGER(4)            :: LN,IC,LN1,LN2,LMN1,LMN2,IM,L1,L2,IBG,L
      INTEGER(4)            :: IDIM,LMN,LM
      REAL(8)               :: DREL(NR)
      REAL(8)               :: G(NR,LMRX)
      REAL(8)               :: E(NR,LMRX)
      INTEGER(4)            :: LMX
      logical(4),allocatable:: tsph(:)
      real(8)               :: pot_hyperfine(nr,lmrx,ndimd)
      real(8)               :: rho_hyperfine(nr,lmrx,ndimd)
!     ******************************************************************
      IF(SOFTCORETYPE.EQ.'SPHERICAL'.OR. &
     &   SOFTCORETYPE.EQ.'NONSPHERICAL'.OR. &
     &   SOFTCORETYPE.EQ.'SPHERICALFIXV'.OR. &
     &   SOFTCORETYPE.EQ.'NONSPHERICALFIXV'.OR. &
     &   SOFTCORETYPE.EQ.'FROZEN') THEN
      ELSE
        CALL ERROR$MSG('ILLEGAL VALUE OF SOFTCORETYPE')
        CALL ERROR$CHVAL('SOFTCORETYPE',SOFTCORETYPE)
        CALL ERROR$STOP('AUGMENTATION_NEWSOFTCORE')
      END IF
!
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ================================================================
!     == reset output data                                          ==
!     ================================================================
      decore=0.d0
      deltao=0.d0
      deltah=0.d0
      deltarho=0.d0
!
!     ================================================================
!     == COLLECT DATA FROM SETUP OBJECT                             ==
!     ================================================================
      CALL SETUP$GETI4('NC',NC)
!     ==  RETURN IF THERE IS NO CORE (I.E. HYDROGEN, HELIUM
      IF(NC.EQ.0) THEN
        PRINT*,'RETURN FROM SOFTCORE BECAUSE NC=0 ',NC,AEZ      
        RETURN
      END IF
!
!     == COLLECT DIMENSIONS FOR PARTIAL WAVES ============================
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
!
!     == COLLECT NODELESS PARTIAL WAVES AND AE PARTIAL WAVES FOR COMPARISON
      ALLOCATE(UPHI(NR,LNX))
      ALLOCATE(TUPHI(NR,LNX))
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('NDLSPHI',NR*LNX,UPHI)
      CALL SETUP$GETR8A('NDLSTPHI',NR*LNX,TUPHI)
!
!     == COLLECT CORE STATES ============================================
      ALLOCATE(EOFI(NC))
      ALLOCATE(LOFI(NC))
      ALLOCATE(SOFI(NC))
      ALLOCATE(FOFI(NC))
      ALLOCATE(NNOFI(NC))
      CALL SETUP$GETR8A('FOFC',NC,FOFI)
      CALL SETUP$GETR8A('EOFC',NC,EOFI)
      CALL SETUP$GETI4A('LOFC',NC,LOFI)
      SOFI(:)=0
      DO I=1,NC
        NNOFI(I)=0
        DO J=1,I-1
          IF(LOFI(J).EQ.LOFI(I))NNOFI(I)=NNOFI(I)+1
        ENDDO
      ENDDO
!
!     == COLLECT CORE STATES ============================================
      allocate(tsph(nc))
      tsph(:)=.false.
print*,'tsph',tsph

!     ==  COLLECT CORE STATES
!      ALLOCATE(PHIC(NR,NC))
!      CALL SETUP$GETR8A('AECOREPSI',NR*NC,PHIC)
!
!     ================================================================
!     == CALCULATE ATOMIC CORE ENERGY                               ==
!     ================================================================
      ALLOCATE(PHIC_AT(NR,NC))
      ALLOCATE(TPHIC_AT(NR,NC))
      ALLOCATE(AEPHI_AT(NR,LNX))
      ALLOCATE(TAEPHI_AT(NR,LNX))
!
      CALL SETUP$GETR8A('ATOMICAEPOT',NR,AEPOT_AT(:,1,1))
      POT(:)=AEPOT_AT(:,1,1)
!
!     == CORE WAVE FUNCTIONS AND CORE DENSITY ============================
      CALL SF_CORESTATES_SPH(GID,NR,NC,LOFI,SOFI,FOFI,NNOFI,POT,RHOC_AT,EOFI,PHIC_AT,TPHIC_AT)
!
!     == CORE ENERGY =====================================================
      CALL SF_ECORE_SPH(GID,NR,NC,EOFI,FOFI,POT,RHOC_AT,ECORE_AT,EKINCAT,AUX)
!
!     == NEW PARTIAL WAVES ===============================================
      CALL SF_PARTWAVES_SPH(GID,NR,LNX,LOX,UPHI,TUPHI &
     &                      ,NC,LOFI,PHIC_AT,TPHIC_AT,AEPHI_AT,TAEPHI_AT)
!
!     == VALENCE DENSITY AND KINETIC ENERGY MATRIX ELEMENTS ==============
      CALL SF_RHOVALENCE_SPH(GID,NR,LNX,LOX,LMNX,NDIMD,DENMAT &
     &                            ,AEPHI_AT,TAEPHI_AT,LMRX,TKIN_AT,RHOV_AT)
!
!     == VALENCE ENERGY AND FULL POTENTIAL
      CALL SF_EVALENCE_SPH(GID,NR,LMNX,LMRX,NDIMD,RHOC_AT,RHOV_AT,DENMAT &
     &                    ,RHOB,VQLM,TKIN_AT,EVALENCE_AT,AEPOT_AT)
!
      CALL SF_POTMATRIX_SPH(GID,NR,LMRX,NDIMD,LNX,LOX,LMNX &
     &                     ,AEPOT_AT,AEPHI_AT,H_AT,O_AT)
!
!
!     =====================================================================
!     ==  hyperfine                                                      ==
!     =====================================================================
       rho_hyperfine(:,:,:)=RHOV_AT(:,:,:)
       rho_hyperfine(:,1,1)=rho_hyperfine(:,1,1)+RHOC_AT(:)
       CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AUX,rho_hyperfine &
      &                           ,VQLM,RHOB,pot_hyperfine,svar)
      CALL SF_HYPERFINE('AT',GID,NR,LMRX,NDIMD,AEZ,rho_hyperfine,pot_hyperfine)
!PRINT*,'O_AT ',O_AT
      H_AT=H_AT+TKIN_AT
!
!     =====================================================================
!     ==  PRINT                                                          ==
!     =====================================================================
PRINT*,'===============ATOMIC CORE========================='
PRINT*,'CORE ENERGY    ',ECORE_AT
PRINT*,'VALENCE ENERGY ',EVALENCE_AT
PRINT*,'TOTAL ENERGY   ',ECORE_AT+EVALENCE_AT
PRINT*,'DECORE         ',0.D0
PRINT*,' '
DO IB=1,NC
  PRINT*,'L,E ',LOFI(IB),EOFI(IB)
ENDDO
PRINT*,' '
!
!     =====================================================================
!     ==  WRAP UP                                                        ==
!     =====================================================================
      IF(SOFTCORETYPE.EQ.'FROZEN') THEN
        DELTARHO(:,:,:)=0.D0
        DELTAH(:,:,:)=0.D0
        DELTAO(:,:,:)=0.D0
        DECORE=0.D0
DO IDIM=1,NDIMD
  DO LMN=1,LMNX
    WRITE(*,FMT='("H_AT",I3,40F10.5)')LMN,H_AT(LMN,:,IDIM)
  ENDDO
ENDDO
DO IDIM=1,NDIMD
  DO LMN=1,LMNX
    WRITE(*,FMT='("O_AT",I3,40F10.5)')LMN,O_AT(LMN,:,IDIM)
  ENDDO
ENDDO
        RETURN
      END IF
!
!     =====================================================================
!     =====================================================================
!     == SCF LOOP FOR SPHERICAL POTENTIAL                                ==
!     =====================================================================
!     =====================================================================
!     =====================================================================
!     == DETERMINE WEIGHT FOR BROYDENS METHOD                             ==
!     =====================================================================
      CALL RADIAL$GETR8(I,'DEX',SVAR1)
      CALL RADIAL$GETR8(I,'R1',SVAR2)
      AUX(:)=svar1*(RHOV_AT(:,1,1)+RHOC_AT(:))*R(:)**3
      BROYDENWEIGHT(:,:,:)=0.D0
      BROYDENWEIGHT(:,1,1)=AUX(:)
      BROYDENWEIGHT(:,:,:)=BROYDENWEIGHT(:,:,:)+1.D0
!
!     =====================================================================
!     == SCF LOOP FOR SPHERICAL POTENTIAL                                ==
!     =====================================================================
      ALLOCATE(AEPHI_SPH(NR,LNX))
      ALLOCATE(TAEPHI_SPH(NR,LNX))
      ALLOCATE(PHIC_SPH(NR,NC))
      ALLOCATE(TPHIC_SPH(NR,NC))
!     == DO NOT FORGET THE EXTERNAL POTENTIAL AND THE BACKGROUND !!
      AEPOT_SPH=AEPOT_AT
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR,10,0.1D0)
      DO ITER=1,NITER
        POTIN=AEPOT_SPH(:,1,1)
!
!       ==================================================================
!       == CORE WAVE FUNCTIONS AND CORE DENSITY =========================
!       ==================================================================
        CALL SF_CORESTATES_SPH(GID,NR,NC,LOFI,SOFI,FOFI,NNOFI,POTIN &
     &      ,RHOC_SPH,EOFI,PHIC_SPH,TPHIC_SPH)
!
!       ==================================================================
!       == CORE ENERGY ===================================================
!       ==================================================================
        CALL SF_ECORE_SPH(GID,NR,NC,EOFI,FOFI,POTIN,RHOC_SPH,ECORE_SPH,EKINCAT,AUX)
!
!       ==================================================================
!       == NEW PARTIAL WAVES                                            ==
!       ==================================================================
        if(SOFTCORETYPE.EQ.'SPHERICALFIXV') THEN
!         == keep original partial waves
          AEPHI_SPH=AEPHI_AT
          TAEPHI_SPH=TAEPHI_AT
        else
!         == recalculate partial waves from nodeless partial waves and new core
          CALL SF_PARTWAVES_SPH(GID,NR,LNX,LOX,UPHI,TUPHI &
     &                      ,NC,LOFI,PHIC_SPH,TPHIC_SPH,AEPHI_SPH,TAEPHI_SPH)
        END IF
!
!       ==================================================================
!       == VALENCE DENSITY AND KINETIC ENERGY MATRIX ELEMENTS ============
!       ==================================================================
        CALL SF_RHOVALENCE_SPH(GID,NR,LNX,LOX,LMNX,NDIMD,DENMAT &
     &                            ,AEPHI_SPH,TAEPHI_SPH,LMRX,TKIN_SPH,RHOV_SPH)
!
!       ==================================================================
!       == VALENCE ENERGY AND FULL POTENTIAL =============================
!       ==================================================================
        CALL SF_EVALENCE_SPH(GID,NR,LMNX,LMRX,NDIMD,RHOC_SPH,RHOV_SPH,DENMAT &
     &                      ,RHOB,VQLM,TKIN_SPH,EVALENCE_SPH,AEPOT_SPH)
!
!       == DO NOT FORGET THE EXTERNAL POTENTIAL AND THE BACKGROUND !!
        SVAR=AEPOT_AT(NR,1,1)-AEPOT_SPH(NR,1,1)
        AEPOT_SPH(:,1,1)=AEPOT_SPH(:,1,1)+SVAR
!
!       ================================================================
!       ==  EXIT IF CONVERGED                                         ==
!       ================================================================
        IF(CONVG) THEN
          CALL BROYDEN$CLEAR
          EXIT
        END IF
!
!       =================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD     ==
!       ================================================================
        XMAX=MAXVAL(ABS(AEPOT_SPH(:,1,1)-POTIN)) 
        CALL BROYDEN$STEP(NR,POTIN,AEPOT_SPH(:,1,1)-POTIN)
        AEPOT_SPH(:,1,1)=POTIN
        CONVG=(XMAX.LT.TOL)
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$STOP('AUGMETATION_NEWSOFTCORE')
      END IF
!
!     ==================================================================
!     ==  HAMILTON MATRIX ELEMENTS                                    ==
!     ==================================================================
      CALL SF_POTMATRIX_SPH(GID,NR,LMRX,NDIMD,LNX,LOX,LMNX &
     &                     ,AEPOT_SPH,AEPHI_SPH,H_SPH,O_SPH)
      H_SPH=H_SPH+TKIN_SPH
!
!     =====================================================================
!     ==  hyperfine                                                      ==
!     =====================================================================
       rho_hyperfine(:,:,:)=RHOV_sph(:,:,:)
       rho_hyperfine(:,1,1)=rho_hyperfine(:,1,1)+RHOC_sph(:)
       CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AUX,rho_hyperfine &
      &                           ,VQLM,RHOB,pot_hyperfine,svar)
      CALL SF_HYPERFINE('sph',GID,NR,LMRX,NDIMD,AEZ,rho_hyperfine,pot_hyperfine)
!
!     =====================================================================
!     ==  PRINT                                                          ==
!     =====================================================================
PRINT*,'===============SCF CORE========================='
PRINT*,'CORE ENERGY    ',ECORE_SPH,ECORE_SPH-ECORE_AT
PRINT*,'VALENCE ENERGY ',EVALENCE_SPH,EVALENCE_SPH-EVALENCE_AT
PRINT*,'TOTAL ENERGY   ',ECORE_SPH+EVALENCE_SPH,ECORE_SPH+EVALENCE_SPH-ECORE_AT-EVALENCE_AT
PRINT*,'DECORE         ',ECORE_SPH+EVALENCE_SPH-ECORE_AT-EVALENCE_AT
PRINT*,' '
DO IB=1,NC
  PRINT*,'L,E ',LOFI(IB),EOFI(IB)
ENDDO
!
!     =====================================================================
!     ==  WRAP UP                                                        ==
!     =====================================================================
      IF(SOFTCORETYPE(1:9).EQ.'SPHERICAL') THEN
        DELTARHO(:,:,:)=RHOV_SPH(:,:,:)-RHOV_AT(:,:,:)
        DELTARHO(:,1,1)=DELTARHO(:,1,1)+RHOC_SPH(:)-RHOC_AT(:)
        DELTAH(:,:,:)=H_SPH-H_AT
        DELTAO(:,:,:)=O_SPH-O_AT
        DECORE=ECORE_SPH+EVALENCE_SPH-ECORE_AT-EVALENCE_AT
        decore=decore-real(sum(conjg(edenmat(:,:,:))*deltao))
DO IDIM=1,NDIMD
  DO LMN=1,LMNX
    WRITE(*,FMT='("DH",I3,40F10.5)')LMN,DELTAH(LMN,:,IDIM)
  ENDDO
ENDDO
DO IDIM=1,NDIMD
  DO LMN=1,LMNX
    WRITE(*,FMT='("DO",I3,40F10.5)')LMN,DELTAO(LMN,:,IDIM)
  ENDDO
ENDDO
        RETURN
      END IF
!
!     =====================================================================
!     =====================================================================
!     ==  NOW NONSPHERICAL                                               ==
!     =====================================================================
!     =====================================================================
      LMX=MAX(LMRX,(MAXVAL(LOFI)+1)**2)
      NC_NSPH=SUM(2*(2*LOFI(:)+1))
      ALLOCATE(EC_NSPH(NC_NSPH))
      ALLOCATE(PHIC_NSPH(NR,LMX,2,NC_NSPH))
      ALLOCATE(TPHIC_NSPH(NR,LMX,2,NC_NSPH))
      ALLOCATE(SPHIC_NSPH(NR,LMX,2,NC_NSPH))
      ALLOCATE(TSPHIC_NSPH(NR,LMX,2,NC_NSPH))
      ALLOCATE(AEPHI_NSPH(NR,LMX,2,2*LMNX))
      ALLOCATE(TAEPHI_NSPH(NR,LMX,2,2*LMNX))
      ALLOCATE(sAEPHI_NSPH(NR,LMX,2,2*LMNX))
      ALLOCATE(TsAEPHI_NSPH(NR,LMX,2,2*LMNX))
      AEPOT_NSPH=AEPOT_SPH ! STARTING POTENTIAL IS RELAXED SPHERICAL POTENTIAL
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR*LMRX*NDIMD,0,0.5D0)
!      CALL BROYDEN$SETWEIGHT(NR*LMRX*NDIMD,BROYDENWEIGHT)
      DO ITER=1,NITER
        POTIN_NSPH=AEPOT_NSPH
!    
!       == CORE WAVE FUNCTIONS AND CORE DENSITY ========================
print*,'timing: nsph calculate core states'
        CALL SF_CORESTATES_NSPH(GID,NR,NDIMD,LMRX,NC,LOFI,SOFI,FOFI,NNOFI,tsph &
     &                  ,POTIN_NSPH,EOFI,LMX,RHOC_NSPH,NC_NSPH,EC_NSPH &
     &                  ,PHIC_NSPH,TPHIC_NSPH,SPHIC_NSPH,TSPHIC_NSPH)
!
PRINT*,'===== CORE STATES DONE===================='
AUX=RHOC_NSPH(:,1,1)*R(:)**2
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'CORE CHARGE ',4.D0*PI*SVAR*Y0,AEZ

IBG=0
DO I=1,NC
  L=LOFI(I)
  DO J=1,2*(2*L+1)
    IBG=IBG+1
    PRINT*,'EC_NSPH ',IBG,EC_NSPH(IBG),EOFI(I)
  ENDDO
ENDDO
!
!       ==================================================================
!       == CORE ENERGY ===================================================
!       ==================================================================
print*,'timing: nsph calculate core energy'
        CALL SF_ECORE_NSPH(GID,NR,NC_NSPH,EC_NSPH,LMRX,NDIMD,POTIN_NSPH &
      &                   ,RHOC_NSPH,ECORE_NSPH)
!
!       ==================================================================
!       == NEW PARTIAL WAVES =============================================
!       ==================================================================
print*,'timing: nsph calculate partial waves'
        IF(SOFTCORETYPE.EQ.'NONSPHERICALFIXV') THEN
          AEPHI_SPH=AEPHI_AT
          TAEPHI_SPH=TAEPHI_AT
        else
          CALL SF_PARTWAVES_NSPH(GID,NR,LNX,LOX,UPHI,TUPHI,LMNX &
      &                ,NC_NSPH,LMX,PHIC_NSPH,TPHIC_NSPH,SPHIC_NSPH,TSPHIC_NSPH &
      &                ,AEPHI_NSPH,TAEPHI_NSPH,sAEPHI_NSPH,TsAEPHI_NSPH)
!!$        call SF_UPDATEPARTWAVES_NSPH(GID,NR,lmx,LMNX,NDIMD,DENMAT &
!!$     &       ,aephi_nsph,taephi_nsph,saephi_nsph,tsaephi_nsph &
!!$     &        ,lmrx,potin_nsph,nc_nsph,phic_nsph,sphic_nsph)
        END IF
!
!       ==================================================================
!       == VALENCE DENSITY AND KINETIC ENERGY MATRIX ELEMENTS ============
!       ==================================================================
print*,'timing: nsph calculate valence density'
        CALL SF_RHOVALENCE_NSPH(GID,NR,LMNX,NDIMD,DENMAT,LMX &
      &                       ,AEPHI_NSPH,TAEPHI_NSPH,sAEPHI_NSPH,TsAEPHI_NSPH &
      &                       ,LMRX,TKIN_NSPH,RHOV_NSPH)
AUX=RHOV_NSPH(:,1,1)*Y0*R(:)**2
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'VALENCE CHARGE ',4.D0*PI*SVAR
AUX=(RHOC_NSPH(:,1,1)+RHOV_NSPH(:,1,1))*R(:)**2
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'TOTAL CHARGE ',4.D0*PI*SVAR*Y0-AEZ
!
!       ==================================================================
!       == VALENCE ENERGY AND FULL POTENTIAL =============================
!       ==================================================================
print*,'timing: nsph calculate valence energy'
        CALL SF_EVALENCE_NSPH(GID,NR,LMNX,LMRX,NDIMD,RHOC_NSPH,RHOV_NSPH,DENMAT &
      &                      ,RHOB,VQLM,TKIN_NSPH,EVALENCE_NSPH,AEPOT_NSPH)
!
!       ================================================================
!       ==  EXIT IF CONVERGED                                         ==
!       ================================================================
        IF(CONVG) THEN
          CALL BROYDEN$CLEAR
          EXIT
        END IF
!       == DO NOT FORGET THE EXTERNAL POTENTIAL AND THE BACKGROUND !!
!        SVAR=AEPOT_AT(NR,1,1)-AEPOT_NSPH(NR,1,1)
!        AEPOT_NSPH(:,1,1)=AEPOT_NSPH(:,1,1)+SVAR
!
!       ================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD     ==
!       ================================================================
print*,'timing: nsph mixing'
        XMAX=MAXVAL(ABS(AEPOT_NSPH-POTIN_NSPH)) 
        CALL BROYDEN$STEP(NR*LMRX*NDIMD,POTIN_NSPH,AEPOT_NSPH-POTIN_NSPH)
        AEPOT_NSPH=POTIN_NSPH
        CONVG=(XMAX.LT.TOL)
PRINT*,'XMAX ',XMAX,TOL,CONVG
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$STOP('AUGMETATION_NEWSOFTCORE')
      END IF
!
!     ==================================================================
!     ==  HAMILTON MATRIX ELEMENTS                                    ==
!     ==================================================================
      CALL SF_POTMATRIX_NSPH(GID,NR,LMRX,NDIMD,LMX,LMNX &
     &                      ,AEPOT_NSPH,AEPHI_NSPH,H_NSPH,O_NSPH)
      H_NSPH=H_NSPH+TKIN_NSPH

!
!     =====================================================================
!     ==  hyperfine                                                      ==
!     =====================================================================
       rho_hyperfine(:,:,:)=RHOV_NSPH(:,:,:)+RHOC_NSPH(:,:,:)
       CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AUX,rho_hyperfine &
      &                           ,VQLM,RHOB,pot_hyperfine,svar)
      CALL SF_HYPERFINE('nsph',GID,NR,LMRX,NDIMD,AEZ,rho_hyperfine,pot_hyperfine)
!
!     =====================================================================
!     ==  WRAP UP                                                        ==
!     =====================================================================
CALL SF_WRITE('RHOV_NSPH.DAT',GID,NR,LMRX*NDIMD,RHOV_NSPH)
DELTARHO=RHOV_NSPH+RHOC_NSPH
CALL SF_WRITE('RHOT_NSPH.DAT',GID,NR,LMRX*NDIMD,DELTARHO)
DELTARHO(:,:,:)=RHOV_NSPH(:,:,:)-RHOV_AT(:,:,:)+RHOC_NSPH(:,:,:)
DELTARHO(:,1,1)=DELTARHO(:,1,1)-RHOC_AT(:)
CALL SF_WRITE('DELTARHO_NSPH.DAT',GID,NR,LMRX*NDIMD,DELTARHO)
!
      DELTARHO(:,:,:)=RHOV_NSPH(:,:,:)-RHOV_AT(:,:,:)+RHOC_NSPH(:,:,:)
      DELTARHO(:,1,1)=DELTARHO(:,1,1)-RHOC_AT(:)
      DELTAH(:,:,:)=H_NSPH-H_AT
      DELTAO(:,:,:)=O_NSPH-O_AT
      DECORE=ECORE_NSPH+EVALENCE_NSPH-ECORE_AT-EVALENCE_AT
      decore=decore-real(sum(conjg(edenmat(:,:,:))*deltao))
print*,'decore before returning from newsoftcore ',decore
DO IDIM=1,NDIMD
  DO LMN=1,LMNX
    WRITE(*,FMT='("DH",I3,40F10.5)')LMN,DELTAH(LMN,:,IDIM)
  ENDDO
ENDDO
DO IDIM=1,NDIMD
  DO LMN=1,LMNX
    WRITE(*,FMT='("DO",I3,40F10.5)')LMN,DELTAO(LMN,:,IDIM)
  ENDDO
ENDDO
      RETURN
!
      DEALLOCATE(LOX)
      DEALLOCATE(UPHI)
      DEALLOCATE(TUPHI)
      DEALLOCATE(AEPHI)
      DEALLOCATE(EOFI)
      DEALLOCATE(LOFI)
      DEALLOCATE(SOFI)
      DEALLOCATE(FOFI)
      DEALLOCATE(NNOFI)
      DEALLOCATE(PHIC_AT)
      DEALLOCATE(TPHIC_AT)
      DEALLOCATE(AEPHI_AT)
      DEALLOCATE(TAEPHI_AT)
!
      DEALLOCATE(AEPHI_SPH)
      DEALLOCATE(TAEPHI_SPH)
      DEALLOCATE(PHIC_SPH)
      DEALLOCATE(TPHIC_SPH)
!
      DEALLOCATE(EC_NSPH)
      DEALLOCATE(PHIC_NSPH)
      DEALLOCATE(TPHIC_NSPH)
      DEALLOCATE(SPHIC_NSPH)
      DEALLOCATE(TSPHIC_NSPH)
      DEALLOCATE(AEPHI_NSPH)
      DEALLOCATE(TAEPHI_NSPH)
!
!     =====================================================================
!     ==WRAP UP                                                          ==
!     =====================================================================
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE SF_WRITE(FILE,GID,NR,NF,F)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NF
      REAL(8)   ,INTENT(IN) :: F(NR,NF)
      REAL(8)               :: R(NR)
      INTEGER(4)            :: I
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE,FORM='FORMATTED')
      DO I=1,NR
        WRITE(100,FMT='(100F20.10)')R(I),F(I,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE SF_HYPERFINE(ID,GID,NR,LMRX,NDIMD,AEZ,RHO,POT)
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: AEZ    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: RHO(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: POT(NR,LMRX,NDIMD)
      REAL(8)               :: R(NR)
      REAL(8)               :: RNUC  !NUCLEAR RADIUS
      REAL(8)               :: PI,Y0
      REAL(8)               :: ISOMER  !ISOMER SHIFT MUST BE MULTIPLIED WITH DRNUC/RNUC
      REAL(8)               :: FERMI   !FERMI CONTACT
      REAL(8)               :: V2(5),EFG(3,3),ANIS(3,3)
      REAL(8)               :: AUX1(NR),AUX2(NR),SVAR1,SVAR2
      REAL(8)               :: EIG(3),U(3,3)
      REAL(8)               :: UNIT
      INTEGER(4)            :: LM
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL PERIODICTABLE$GET(NINT(AEZ),'RNUC',RNUC)
!
!     ======================================================================
!     ==  ISOMER SHIFT                                                    ==
!     ======================================================================
      AUX1(:)=RHO(:,1,1)*Y0*((R(:)/RNUC)**4-(R(:)/RNUC)**2)
      CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
      CALL RADIAL$VALUE(GID,NR,AUX2,RNUC,SVAR1)
      AUX1(:)=((R(:)/RNUC)**4-(R(:)/RNUC)**2)
      CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
      CALL RADIAL$VALUE(GID,NR,AUX2,RNUC,SVAR2)
      ISOMER=svar1/svar2
PRINT*,'ISOMERSHIFT (rho(0)',isomer
!
!     ======================================================================
!     ==  ELECTRIC FIELD GRADIENT                                         ==
!     ======================================================================
      DO LM=5,LMRX
        AUX1(:)=POT(:,LM,1)
        CALL RADIAL$VALUE(GID,NR,AUX1,RNUC,V2(LM-4))
      ENDDO
      V2(:)=V2(:)/RNUC**2

!     ==  FROM SPHERICAL HARMONICS TO CARTESIAN COORDINATES
!         __VLM=V(R)/R**2___VXYZ=D2V/(DI*DJ)____________________________
      CALL DTOXYZ(V2,EFG)  ! SEE EFG
!     == DETERMINE EIGENVALUES AND EIGENVECTORS
      CALL LIB$DIAGR8(3,EFG,EIG,U)
!     == UNITS
      UNIT=1.D-21
      CALL CONSTANTS('VOLT',SVAR1) ; UNIT=UNIT/SVAR1 
      CALL CONSTANTS('METER',SVAR1) ; UNIT=UNIT*SVAR1**2 
      UNIT=-UNIT   ! THE ELEMENTARY CHARGE IS -1 ELECTRON CHARGE
PRINT*,UNIT,'SHOULD BE EQUAL TO ',-9.7175D0
      PRINT*,'EFG EIGENVALUES IN [10**21*V/METER**2] ',EIG*UNIT
!
!     ======================================================================
!     ==  FERMI CONTACT                                                   ==
!     ======================================================================
      IF(NDIMD.EQ.1) THEN
        AUX1(:)=0.D0
      ELSE IF(NDIMD.EQ.2) THEN
        AUX1(:)=4.D0*PI*RHO(:,1,2)*Y0*R(:)**2
      ELSE IF(NDIMD.EQ.4) THEN       
        AUX1(:)=4.D0*PI*SQRT(RHO(:,1,2)**2+RHO(:,1,3)**2+RHO(:,1,4)**2)*Y0*R(:)**2
      END IF  
      CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
      CALL RADIAL$VALUE(GID,NR,AUX2,RNUC,FERMI)
      UNIT=2.D0/3.D0
      CALL CONSTANTS('MU0',SVAR1)          ; UNIT=UNIT*SVAR1 
!     == MULTIPLY WITH MAGNETIC MOMENT OF THE ELECTRON
      CALL CONSTANTS('BOHRMAGNETON',SVAR1) ; UNIT=UNIT*SVAR1
      CALL CONSTANTS('GE',SVAR1)           ; UNIT=UNIT*SVAR1
      CALL CONSTANTS('HBAR',SVAR1)         ; UNIT=UNIT*0.5D0*SVAR1
!     ==   
      CALL CONSTANTS('TESLA',SVAR1)        ; UNIT=UNIT/SVAR1
      PRINT*,'THIS NUMBER SHOULD BE 104.98 : ',UNIT
PRINT*,'FERMI CONTACT TERM [TESLA] ',FERMI*UNIT
!      CALL REPORT$R8VAL(NFILO,'SPIN DENSITY AT THE NUCLEUS OF ATOM ' &
!      &                            //TRIM(ATOMNAME),FERMI,'1/ABOHR^3')
!      CALL REPORT$R8VAL(NFILO,'FERMI CONTACT HYPERFINE FIELD FOR ATOM ' &
!      &                            //TRIM(ATOMNAME),FERMI*UNIT,'TESLA')
!
!     ======================================================================
!     ==  ANISOTROPIC HYPERFINE PARAMETERS                                ==
!     ======================================================================
      IF(LMRX.GE.9) THEN
        DO LM=5,9
          CALL RADIAL$POISSON(GID,NR,2,RHO(1,LM,2),AUX1)
          CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
          AUX2(:)=AUX2(:)*R(:)**2
          CALL RADIAL$VALUE(GID,NR,AUX2,RNUC,V2(LM-4))
        ENDDO
        CALL DTOXYZ(V2,ANIS)  ! SEE EFG OBJECT
      END IF
!
      UNIT=1.D0/(4.D0*PI)
      CALL CONSTANTS('MU0',SVAR1)          ; UNIT=UNIT*SVAR1 
!     ============================================================
      CALL CONSTANTS('BOHRMAGNETON',SVAR1) ; UNIT=UNIT*SVAR1
      CALL CONSTANTS('GE',SVAR1)           ; UNIT=UNIT*SVAR1
      CALL CONSTANTS('HBAR',SVAR1)         ; UNIT=UNIT*0.5D0*SVAR1
!     ============================================================
      CALL CONSTANTS('TESLA',SVAR1)        ; UNIT=UNIT/SVAR1
      PRINT*,'THIS NUMBER SHOULD BE 12.531 : ',UNIT
!
!!$      CALL REPORT$STRING(NFILO,'ANISOTROPIC HYPERFINE FIELD FOR ATOM ' &
!!$     &                //TRIM(ATOMNAME)//' IN UNITS OF TESLA')
!!$      CALL REPORT$STRING(NFILO,'MULTIPLY MATRIX WITH A UNITY VECTOR ' &
!!$     &                //'INTO THE DIRECTION OF THE MAGNETIC FIELD')
      CALL LIB$DIAGR8(3,ANIS,EIG,U)              
!      __ FACTOR 0.5 IS THE ELECTRON SPIN____________
!!$      WRITE(NFILO,FMT='(T5,"VALUE",T20,"DIRECTION")')
!!$      DO I=1,3
!!$        WRITE(NFILO,FMT='(T5,F10.5,T20,3F10.5)')UNIT*EIG(I),(U(J,I),J=1,3)
!!$      ENDDO
      PRINT*,'ANISOTROPIC HYPERFINE[TESLA]',EIG(:)*UNIT
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE SF_CORESTATES_SPH(GID,NR,NB,LOFI,SO,F,NN,POT &
     &                            ,RHO,E,PHI,TPHI)
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
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NB)
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR,NB)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)                    :: AUX(NR)
      INTEGER(4),PARAMETER       :: NITER=100
      REAL(8)   ,PARAMETER       :: TOL=1.D-4
      INTEGER(4)                 :: IB,IR,I
      REAL(8)                    :: SVAR
      REAL(8)                    :: C0LL,PI,Y0
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RHO(:)=0.D0
      DO IB=1,NB
        DO I=1,NITER
          SVAR=E(IB)
          CALL RELATIVISTICCORRECTION(GID,NR,POT,E(IB),DREL)
          AUX(:)=0.D0
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),DREL,AUX,NN(IB),POT &
     &                  ,E(IB),PHI(:,IB))
          IF(ABS(E(IB)-SVAR).LT.TOL) EXIT
        ENDDO
        AUX(:)=(R(:)*PHI(:,IB))**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
        TPHI(:,IB)=(E(IB)-POT(:)*Y0)*PHI(:,IB)
        RHO(:)=RHO(:)+F(IB)*C0LL*PHI(:,IB)**2
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE SF_CORESTATES_NSPH(GID,NR,NDIMD,LMRX,NB,LOFI,SO,F,NN,tsph &
     &                          ,POT,E,LMX,RHOC,NBG,EBG,PHI,TPHI,SPHI,TSPHI)
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)   :: GID
      INTEGER(4) ,INTENT(IN)   :: NR
      INTEGER(4) ,INTENT(IN)   :: NDIMD
      INTEGER(4) ,INTENT(IN)   :: NB        ! #(SHELLS)
      INTEGER(4) ,INTENT(IN)   :: LMRX      ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)   :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)   :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)   :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)   :: F(NB)     !OCCUPATION
      logical(4) ,intent(in)   :: tsph(nb)  ! use spherical potential 
      REAL(8)    ,INTENT(IN)   :: POT(NR,LMRX,NDIMD)   !POTENTIAL
      REAL(8)    ,INTENT(INOUT):: E(NB)     !ONE-PARTICLE ENERGIES
      INTEGER(4) ,INTENT(IN)   :: LMX       ! #(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)   :: NBG       
      REAL(8)    ,INTENT(OUT)  :: EBG(NBG)
      COMPLEX(8) ,INTENT(OUT)  :: PHI(NR,LMX,2,NBG)
      COMPLEX(8) ,INTENT(OUT)  :: TPHI(NR,LMX,2,NBG)
      COMPLEX(8) ,INTENT(OUT)  :: SPHI(NR,LMX,2,NBG)
      COMPLEX(8) ,INTENT(OUT)  :: TSPHI(NR,LMX,2,NBG)
      REAL(8)    ,INTENT(OUT)  :: RHOC(NR,LMRX,NDIMD)
      REAL(8)                  :: GINH(NR,LMX,2)
      LOGICAL(4)               :: TSELECT(2*LMX)
      REAL(8)                  :: DE(2*LMX)
      REAL(8)                  :: R(NR)
      REAL(8)                  :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)                  :: AUX(NR),AUX1(NR)
      COMPLEX(8)               :: CAUX(NR,2,2)
      INTEGER(4),PARAMETER     :: NITER=100
      REAL(8)   ,PARAMETER     :: TOL=1.D-4
      INTEGER(4)               :: IB,IR,I,LM1,LM2,LM3,LMR,IBG,IS1,IS2
      INTEGER(4)               :: I1,I2
      REAL(8)                  :: SVAR
      REAL(8)                  :: C0LL,PI,Y0
      REAL(8)                  :: CG
      LOGICAL(4)               :: TCHK
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      IF(NBG.NE.SUM(2*(2*LOFI(:)+1))) THEN
        CALL ERROR$MSG('NBG INCONSISTENT WITH LOFI AND NB')
        CALL ERROR$STOP('SF_NONSPHAERHO')
      END IF
!
!     ========================================================================
!     == MAP POTENTIAL INTO LOCAL SPIN-POLARIZED ARRAY                      ==
!     ========================================================================
      IBG=0
      DO IB=1,NB
!
!       ======================================================================
!       == DETERMINE SOLUTION IN THE SPHERICAL UNPOLARIZED POTENTIAL        ==
!       ======================================================================
        DO I=1,NITER
          SVAR=E(IB)
          CALL RELATIVISTICCORRECTION(GID,NR,POT(:,1,1),E(IB),DREL)
          AUX(:)=0.D0
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),DREL,AUX,NN(IB),POT(:,1,1) &
     &                 ,E(IB),AUX1)
          IF(ABS(E(IB)-SVAR).LT.TOL) EXIT
        ENDDO
        I1=IBG+1
        I2=IBG+2*(2*LOFI(IB)+1)
!
!       ======================================================================
!       == DETERMINE NONSPHERICAL SOLUTIONS IN THIS SHELL                   ==
!       ======================================================================
        if(tsph(ib)) then
          call sf_soboundstates(gid,nr,lofi(ib),nn(ib),lmx,pot(:,:,1),e(ib) &
    &               ,i2-i1+1,ebg(i1:i2),phi(:,:,:,i1:i2),tphi(:,:,:,i1:i2) &
    &               ,sphi(:,:,:,i1:i2),tsphi(:,:,:,i1:i2))
        else
          GINH=0.D0
          CALL RADIAL$NONSPHBOUND(GID,NR,NDIMD,LMX,LMRX,POT,DREL,GINH,E(IB) &
     &                         ,2*(2*LOFI(IB)+1),EBG(I1:I2) &
     &                         ,PHI(:,:,:,I1:I2),TPHI(:,:,:,I1:I2) &
     &                         ,SPHI(:,:,:,I1:I2),TSPHI(:,:,:,I1:I2),TCHK)
        end if
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('RADIAL$NONSPHBOUND FINISHED WITH ERROR')
          CALL ERROR$STOP('SF_CORESTATES_NSPH')
        END IF
        IBG=I2
PRINT*,'SPHERICAL LOOP ',IB,E(IB),SUM(EBG(I1:I2))/REAL(I2-I1+1),F(IB),I2-I1+1
PRINT*,'ENERGIES       ',EBG(I1:I2)
PRINT*,'-----------------------------------------------------'
      ENDDO
!
!     ======================================================================
!     == CORE DENSITY                                                     ==
!     ======================================================================
PRINT*,'NOW CALCULATE CORE CHARGE'
      RHOC(:,:,:)=0.D0
      DO LMR=1,LMRX
        CAUX(:,:,:)=(0.D0,0.D0)
        DO LM1=1,LMX
          DO LM2=1,LMX
            CALL CLEBSCH(LMR,LM1,LM2,CG)
            IF(CG.EQ.0.D0) CYCLE
            DO IBG=1,NBG
              DO IS1=1,2
                DO IS2=1,2
                  CAUX(:,IS1,IS2)=CAUX(:,IS1,IS2) &
      &                  +CG*CONJG(PHI(:,LM1,IS1,IBG))*PHI(:,LM2,IS2,IBG) &
      &                  +CG*CONJG(sPHI(:,LM1,IS1,IBG))*sPHI(:,LM2,IS2,IBG)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
        RHOC(:,LMR,1)=RHOC(:,LMR,1)+REAL(CAUX(:,1,1)+CAUX(:,2,2))
        IF(NDIMD.EQ.2) THEN
          RHOC(:,LMR,2)=RHOC(:,LMR,2)+REAL(CAUX(:,1,1)-CAUX(:,2,2))
        ELSE IF(NDIMD.EQ.4) THEN
          RHOC(:,LMR,2)=RHOC(:,LMR,2)+REAL(CAUX(:,1,2)+CAUX(:,2,1))
          RHOC(:,LMR,3)=RHOC(:,LMR,3)+AIMAG(CAUX(:,1,2)-CAUX(:,2,1))
          RHOC(:,LMR,4)=RHOC(:,LMR,4)+REAL(CAUX(:,1,1)-CAUX(:,2,2))
        END IF
      ENDDO
AUX=RHOC(:,1,1)*Y0*R(:)**2
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'CORE CHARGE IN CORESTATES',4.D0*PI*SVAR
      RETURN
      END
!
!     ......................................................................
      subroutine sf_soboundstates(gid,nr,l,nn,lmx,pot,enu,nphi,e,phi,tphi,sphi,tsphi)
      implicit none
      integer(4),intent(in) :: gid
      integer(4),intent(in) :: nr
      integer(4),intent(in) :: l
      integer(4),intent(in) :: nn
      integer(4),intent(in) :: lmx
      real(8)   ,intent(in) :: pot(nr)
      real(8)   ,intent(in) :: enu
      integer(4),intent(in) :: nphi
      real(8)   ,intent(out):: e(nphi)
      complex(8),intent(out):: phi(nr,lmx,2,nphi)   !large component
      complex(8),intent(out):: tphi(nr,lmx,2,nphi)  ! kinetic energy * large comp.
      complex(8),intent(out):: sphi(nr,lmx,2,nphi)  ! small component
      complex(8),intent(out):: tsphi(nr,lmx,2,nphi) ! kinetic energy * small comp.
      real(8)               :: phi1(nr,2)
      real(8)               :: tphi1(nr,2)
      real(8)               :: e1(2)
      real(8)               :: drel(nr)
      real(8)               :: aux(nr),svar
      real(8)               :: eprev
      real(8)   ,parameter  :: tol=1.d-6
      integer(4),parameter  :: niter=100
      integer(4)            :: iso,i,iphi,m,lmp,lmm,is,lm
      integer(4)            :: so
      real(8)               :: clm,clmplus1
      complex(8)            :: cp,cm
      real(8)               :: pi,y0
      real(8)               :: r(nr)
!     **********************************************************************
      if(nphi.ne.2*(2*l+1)) then
        call error$msg('number of states requested inconsistent with angular momentum')
        call error$i4val('nphi',nphi)
        call error$i4val('l',l)
        call error$stop('sf_soboundstates')
      end if
      pi=4.d0*datan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
      phi1=0.d0
      e1=0.d0
      CALL RELATIVISTICCORRECTION(GID,NR,POT,Enu,DREL)
!
!     ============================================================================
!     == energies and wave functions for parallel and antiparallel spin-orbit   ==
!     ============================================================================
      do iso=1,2
!       iso=1 -> antiparallel ; iso=2 -> parallel
        e1(iso)=enu
        so=-3+2*iso
        if(l.eq.0) so=0.d0
        if(l.eq.0.and.iso.eq.1) cycle
        DO I=1,NITER
          eprev=E1(iso)
          AUX(:)=0.D0
          CALL BOUNDSTATE(GID,NR,L,SO,DREL,AUX,NN,POT,E1(iso),phi1(:,iso))
          IF(ABS(E1(iso)-eprev).LT.TOL) EXIT
        ENDDO
        tphi1(:,iso)=(e1(iso)-pot(:)*y0)*phi1(:,iso)
      enddo
!
!     ============================================================================
!     ==  large component                                                       ==
!     ============================================================================
      phi(:,:,:,:)=(0.d0,0.d0)
      tphi(:,:,:,:)=(0.d0,0.d0)
      iphi=0
!     == antiparallel spin and orbit (no contribution for l=0)
      do m=-l,l-1
        iphi=iphi+1
        e(iphi)=e1(1)
        clm=sqrt(real(l-m,kind=8)/real(2*l+1,kind=8))
        call sf_torealsphericalharmonics(l,m,lmp,lmm,cp,cm)
        phi(:,lmp,1,iphi)=cp*clm*phi1(:,1)
        tphi(:,lmp,1,iphi)=cp*clm*tphi1(:,1)
        if(m.ne.0) then
          phi(:,lmm,1,iphi)=cm*clm*phi1(:,1)
          tphi(:,lmm,1,iphi)=cm*clm*tphi1(:,1)
        end if
        clmplus1=-sqrt(real(l+m+1,kind=8)/real(2*l+1,kind=8))
        call sf_torealsphericalharmonics(l,m+1,lmp,lmm,cp,cm)
        phi(:,lmp,2,iphi)=cp*clmplus1*phi1(:,1)
        tphi(:,lmp,2,iphi)=cp*clmplus1*tphi1(:,1)
        if(m+1.ne.0) then
          phi(:,lmm,2,iphi)=cm*clmplus1*phi1(:,1)
          tphi(:,lmm,2,iphi)=cm*clmplus1*tphi1(:,1)
        end if
      enddo
!     == parallel spin and orbit (here also l=0)
      do m=-l-1,l
        iphi=iphi+1
        e(iphi)=e1(2)
        if(m.ge.-l) then
          clm=sqrt(real(l+m+1,kind=8)/real(2*l+1,kind=8))
          call sf_torealsphericalharmonics(l,m,lmp,lmm,cp,cm)
          phi(:,lmp,1,iphi)=cp*clm*phi1(:,2)
          tphi(:,lmp,1,iphi)=cp*clm*tphi1(:,2)
          if(m.ne.0) then
            phi(:,lmm,1,iphi)=cm*clm*phi1(:,2)
            tphi(:,lmm,1,iphi)=cm*clm*tphi1(:,2)
          end if
        end if
        if(m+1.le.l) then
          clmplus1=sqrt(real(l-m,kind=8)/real(2*l+1,kind=8))
          call sf_torealsphericalharmonics(l,m+1,lmp,lmm,cp,cm)
          phi(:,lmp,2,iphi)=cp*clmplus1*phi1(:,2)
          Tphi(:,lmp,2,iphi)=cp*clmplus1*Tphi1(:,2)
          if(m+1.ne.0) then
            phi(:,lmm,2,iphi)=cm*clmplus1*phi1(:,2)
            tphi(:,lmm,2,iphi)=cm*clmplus1*tphi1(:,2)
          end if
        end if
      enddo
sphi=(0.d0,0.d0)
tsphi=(0.d0,0.d0)
print*,'warning! small component switched off in sf_soboundstate'
return
!
!     ============================================================================
!     ==  small component                                                       ==
!     ============================================================================
      call radial_smallcomponent(gid,nr,lmx,nphi,drel,phi,sphi)
      do iphi=1,nphi
        aux(:)=e(iphi)-pot(:)*y0
        do is=1,2
          do lm=1,lmx
            tsphi(:,lm,is,iphi)=aux(:)*sphi(:,lm,is,iphi)
          enddo
        enddo
      enddo
!
!     ============================================================================
!     ==  renormalize                                                           ==
!     ============================================================================
      do iphi=1,nphi
        aux(:)=0.d0
        do is=1,2
          do lm=1,lmx
            aux(:)=aux(:)+abs(phi(:,lm,is,iphi))**2+abs(sphi(:,lm,is,iphi))**2
          enddo
        enddo
        aux(:)=aux(:)*r(:)**2
        call radial$integral(gid,nr,aux,svar)
        svar=1.d0/sqrt(svar)
        phi(:,:,:,iphi)=phi(:,:,:,iphi)*svar
        tphi(:,:,:,iphi)=tphi(:,:,:,iphi)*svar
        sphi(:,:,:,iphi)=sphi(:,:,:,iphi)*svar
        tsphi(:,:,:,iphi)=tsphi(:,:,:,iphi)*svar
      enddo
      return
      end
!
!     ....................................................................
      subroutine sf_torealsphericalharmonics(l,m,lmp,lmm,cp,cm)
!     **                                                                **
!     **  an spherical harmonics |l,m>, which is an angular momentum    **
!     **  eigenstate with quantum numbers l and m is represented        **
!     **  in terms of real spherical harmonics as                       **
!     **     |l,m>=cp*ylm(lmp)+cm*ylm(lmm)                              **
!     **  where cp,cm are complex coefficients and lmp and lmm are      **
!     **  an combined index identifying real spherical harmonics        **
!     **                                                                **
!     **  caution: There are different choices of spherical harmonics   **
!     **                                                                **
!     **                                                                **
      implicit none
      integer(4),intent(in) :: l
      integer(4),intent(in) :: m
      integer(4),intent(out):: lmp
      integer(4),intent(out):: lmm
      complex(8),intent(out):: cp
      complex(8),intent(out):: cm
      complex(8),parameter  :: ci=(0.d0,1.d0)
      complex(8)            :: sqr2inv
!     ******************************************************************
      sqr2inv=(1.d0,0.d0)
      sqr2inv=sqr2inv/sqrt(2.d0)
      if(abs(m).gt.l) then
        call error$stop('sf_torealsphericalharmonics')
      end if
      lmp=l**2+l+1+abs(m)
      lmm=l**2+l+1-abs(m)
      if(m.gt.0) then
        cp=sqr2inv
        cm=sqr2inv*ci
      else if(m.lt.0) then
        cp=(-1.d0)**m*sqr2inv
        cm=-(-1.d0)**m*sqr2inv*ci
      else
        cp=(1.d0,0.d0)
        cm=(0.d0,0.d0)
      end if
      if(abs(abs(cp)**2+abs(cm)**2-1.d0).gt.1.d-8) then
        call error$msg('coefficients are not nomalized')
        call error$r8val('norm**2',abs(cp)**2+abs(cm)**2)
        call error$stop('sf_torealsphericalharmonics')
      end if

      return
      end
!
!     ....................................................................
      SUBROUTINE SF_ECORE_SPH(GID,NR,NB,E,F,POTIN,RHO,ETOT,EKIN,POTOUT)
!     **                                                               **
!     ** CALCULATES ENERGY OF A SPHERICAL-NON-SPIN-POLARIZED ATOM       **
!     ** FOR GIVEN ENERGY EIGENVALUES, INPUT POTENTIAL AND DENSITY      **
!     **                                                                **
!     ** REMARK: USES SETUP OBJECT, WHICH MUST BE PROPERLY SELECTED     **
!     **                                                                **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NB
      REAL(8)   ,INTENT(IN) :: E(NB)
      REAL(8)   ,INTENT(IN) :: F(NB)      ! occupations
      REAL(8)   ,INTENT(IN) :: POTIN(NR)  !input potential
      REAL(8)   ,INTENT(IN) :: RHO(NR)    !core density
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: EKIN
      REAL(8)   ,INTENT(OUT):: POTOUT(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: POTNUC(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: EXC,EH,SVAR
!     *****************************************************************
      CALL RADIAL$R(GID,NR,R)
      CALL AUGMENTATION_XC(GID,NR,1,1,RHO,EXC,POTOUT)
      POTOUT(1)=POTOUT(2)
      CALL SETUP$GETR8A('NUCPOT',NR,POTNUC)
      POTOUT(:)  = POTOUT(:)+POTNUC(:)
!
      CALL RADIAL$POISSON(GID,NR,0,RHO,AUX)
      POTOUT(:)  = POTOUT(:)+AUX(:)
      AUX(:)=(POTNUC(:)+0.5D0*AUX(:))*RHO(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,EH)
!
      AUX(:)=POTIN*RHO(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      EKIN=SUM(E*F)-SVAR
      ETOT=EKIN+EXC+EH
PRINT*,'CORE EBAND   ',SUM(E*F)
PRINT*,'CORE EKIN    ',EKIN
PRINT*,'CORE EXC     ',EXC
PRINT*,'CORE EHARTREE',EH
PRINT*,'CORE ECORE   ',ETOT
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE SF_ECORE_NSPH(GID,NR,NB,E,LMRX,NDIMD,POTIN,RHO,ETOT)
!     **                                                                **
!     ** CALCULATES ENERGY OF A SPHERICAL-NON-SPIN-POLARIZED ATOM       **
!     ** FOR GIVEN ENERGY EIGENVALUES, INPUT POTENTIAL AND DENSITY      **
!     **                                                                **
!     ** REMARK: USES SETUP OBJECT, WHICH MUST BE PROPERLY SELECTED     **
!     **                                                                **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: E(NB)
      REAL(8)   ,INTENT(IN) :: POTIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: RHO(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: POTNUC(NR)
      REAL(8)               :: POTOUT(NR,LMRX,NDIMD)
      REAL(8)               :: R(NR)
      REAL(8)               :: EXC,EHARTREE,EKIN
      INTEGER(4)            :: L
      INTEGER(4)            :: IDIMD,LM
!     *****************************************************************
      CALL RADIAL$R(GID,NR,R)
!     == XC POTENTIAL
      CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,EXC,POTOUT)
      POTOUT(1,:,:)=POTOUT(2,:,:)
!
!     == HARTREE ENERGY 
      CALL SETUP$GETR8A('NUCPOT',NR,POTNUC)
      POTOUT(:,1,1)  = POTOUT(:,1,1)+POTNUC(:)
      AUX(:)=POTNUC(:)*RHO(:,1,1)
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1)+1.D-3))
        CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM,1),AUX1)
        POTOUT(:,LM,1)  = POTOUT(:,LM,1)+AUX1(:)
        AUX(:)=AUX(:)+0.5D0*AUX1*RHO(:,LM,1)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,EHARTREE)
!
!     == KINETIC ENERGY
      AUX(:)=0.D0
      DO IDIMD=1,NDIMD
        DO LM=1,LMRX
          AUX(:)=AUX(:)+POTIN(:,LM,IDIMD)*RHO(:,LM,IDIMD)
        ENDDO
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,EKIN)
      EKIN=SUM(E)-EKIN
!
      ETOT=EKIN+EXC+EHARTREE
PRINT*,'CORE EBAND   ',SUM(E)
PRINT*,'CORE EKIN    ',EKIN
PRINT*,'CORE EXC     ',EXC
PRINT*,'CORE EHARTREE',EHARTREE
PRINT*,'CORE ECORE   ',ETOT

      RETURN
      END
!
!     ...................................................................
      SUBROUTINE SF_PARTWAVES_SPH(GID,NR,LNX,LOX,UPHI,TUPHI &
     &                ,NC,LOFC,PHIC,TPHIC,AEPHI,TAEPHI)
!     **                                                               **
!     **  CONSTRUCTS NEW PARTIAL WAVES BY ORTHOGONALIZING WITH         **
!     **  RESPECT TO CORE STATES                                       **
!     **                                                               **
!     **  CORE STATES ARE ASSUMED TO BE ORTHOGONAL TO EACH OTHER       **
!     **                                                               **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: UPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: TUPHI(NR,LNX)
      INTEGER(4),INTENT(IN) :: NC     ! #(CORE STATES)
      INTEGER(4),INTENT(IN) :: LOFC(NC)
      REAL(8)   ,INTENT(IN) :: PHIC(NR,NC)
      REAL(8)   ,INTENT(IN) :: TPHIC(NR,NC)
      REAL(8)   ,INTENT(OUT):: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: TAEPHI(NR,LNX)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: SVAR
      REAL(8)               :: S(NC,LNX)
      INTEGER(4)            :: IC,LN
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
      AEPHI=UPHI
      TAEPHI=TUPHI
      S(:,:)=0.D0
      DO IC=1,NC
        AUX1(:)=PHIC(:,IC)*R(:)**2
        DO LN=1,LNX
          IF(LOX(LN).NE.LOFC(IC)) CYCLE
          AUX(:)=AEPHI(:,LN)*AUX1(:)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          S(IC,LN)=SVAR
        ENDDO
      ENDDO  
      DO IC=1,NC
        DO LN=1,LNX
          IF(LOX(LN).NE.LOFC(IC)) CYCLE
          AEPHI(:,LN) =AEPHI(:,LN) -PHIC(:,IC)*S(IC,LN)
          TAEPHI(:,LN)=TAEPHI(:,LN)-TPHIC(:,IC)*S(IC,LN)
        ENDDO
      ENDDO  
!!$      DO IC=1,NC
!!$        AUX1(:)=PHIC(:,IC)*R(:)**2
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.LOFC(IC)) CYCLE
!!$          AUX(:)=AEPHI(:,LN)*AUX1(:)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$          AEPHI(:,LN) =AEPHI(:,LN) -PHIC(:,IC)*SVAR
!!$          TAEPHI(:,LN)=TAEPHI(:,LN)-TPHIC(:,IC)*SVAR
!!$        ENDDO
!!$      ENDDO  
      RETURN
      END
!
!      ...................................................................
       SUBROUTINE SF_PARTWAVES_NSPH(GID,NR,LNX,LOX,UPHI,TUPHI,LMNX &
      &                ,NC,LMX,PHIC,TPHIC,SPHIC,TSPHIC &
      &                ,AEPHI_NSPH,TAEPHI_NSPH,sAEPHI_NSPH,TsAEPHI_NSPH)
!      **                                                               **
!      **  CALCULATES THE TOTAL DENSITY ASSUMING A NONCOLLINEAR MODEL   **
!      **  WITH COMPLEX WAVE FUNCTIONS AND A COMPLEX DENSITY MATRIX     **
!      **                                                               **
!      **  THE NODE-LESS PARTIAL WAVES ARE STORED REAL                  **
!      **                                                               **
!      **  ASSUMES THAT THE CORE STATES ARE ORTHOGONAL                  **
!      **                                                               **
USE PERIODICTABLE_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: GID
       INTEGER(4),INTENT(IN) :: NR
       INTEGER(4),INTENT(IN) :: LMNX
       INTEGER(4),INTENT(IN) :: LNX
       INTEGER(4),INTENT(IN) :: LOX(LNX)
       REAL(8)   ,INTENT(IN) :: UPHI(NR,LNX)
       REAL(8)   ,INTENT(IN) :: TUPHI(NR,LNX)
       INTEGER(4),INTENT(IN) :: NC    !#(CORE STATES)
       INTEGER(4),INTENT(IN) :: LMX   !#(ANGULAR MOMENTA FOR WAVE FUNCTIONS)
       COMPLEX(8),INTENT(IN) :: PHIC(NR,LMX,2,NC)
       COMPLEX(8),INTENT(IN) :: TPHIC(NR,LMX,2,NC)
       COMPLEX(8),INTENT(IN) :: SPHIC(NR,LMX,2,NC)
       COMPLEX(8),INTENT(IN) :: TSPHIC(NR,LMX,2,NC)
       COMPLEX(8),INTENT(OUT):: AEPHI_NSPH(NR,LMX,2,LMNX*2) !(NR,LMX,NSPIN,NV)
       COMPLEX(8),INTENT(OUT):: TAEPHI_NSPH(NR,LMX,2,LMNX*2) !(NR,LMX,NSPIN,NV)
       COMPLEX(8),INTENT(OUT):: sAEPHI_NSPH(NR,LMX,2,LMNX*2) !(NR,LMX,NSPIN,NV)
       COMPLEX(8),INTENT(OUT):: TsAEPHI_NSPH(NR,LMX,2,LMNX*2) !(NR,LMX,NSPIN,NV)
       COMPLEX(8)            :: S(NC,2*LMNX) ! <PHI_C|U>
       REAL(8)               :: PI,Y0
       COMPLEX(8)            :: CAUX(NR)
       REAL(8)               :: AUX(NR),SVAR1,SVAR2
       INTEGER(4)            :: NV ! 2*LMNX=#(PARTIAL WAVES) EXPANDED
       COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
       INTEGER(4)            :: LMN,LN,L,M,LM,IS,IC,IV,IV1,IV2,LMR,LM1,LM2
       REAL(8)               :: R(NR)
       COMPLEX(8)            :: OV(NC,NC)
INTEGER(4) :: IC1,IC2
!      ***********************************************************************
       PI=4.D0*DATAN(1.D0)
       Y0=1.D0/SQRT(4.D0*PI)
       CALL RADIAL$R(GID,NR,R)
       NV=2*LMNX
!
!      =======================================================================
!      == EVALUATE OVERLAPS WITH UPHI                                       ==
!      =======================================================================
       S(:,:)=(0.D0,0.D0)
       LMN=0
       DO LN=1,LNX
         L=LOX(LN)
         LM=L**2
         DO M=1,2*L+1
           LMN=LMN+1
           LM=LM+1
           IF(LM.GT.LMX) CYCLE
           DO IS=1,2
             DO IC=1,NC
               CAUX(:)=CONJG(PHIC(:,LM,IS,IC))*UPHI(:,LN)*R(:)**2
               AUX(:)=REAL(CAUX(:))
               CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
               AUX(:)=AIMAG(CAUX)
               CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
               S(IC,(IS-1)*LMNX+LMN)=CMPLX(SVAR1,SVAR2)
             ENDDO
           ENDDO
         ENDDO
       ENDDO

!!$       DO IC1=1,NC
!!$         DO IC2=IC1,NC
!!$           CAUX(:)=(0.D0,0.D0)
!!$           DO IS=1,2
!!$             DO LM=1,LMX
!!$               CAUX(:)=CAUX(:)+CONJG(PHIC(:,LM,IS,IC1))*PHIC(:,LM,IS,IC2)
!!$             ENDDO
!!$           ENDDO
!!$           CAUX=CAUX*R**2
!!$           AUX(:)=REAL(CAUX(:))
!!$           CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$           AUX(:)=AIMAG(CAUX)
!!$           CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$           OV(IC1,IC2)=CMPLX(SVAR1,SVAR2)
!!$           OV(IC2,IC1)=CMPLX(SVAR1,-SVAR2)
!!$         ENDDO
!!$PRINT*,'OV ',IC1,OV(IC1,:)
!!$       ENDDO
!!$       CALL LIB$MATRIXSOLVENEWC8(NC,NC,NV,OV,S,S) 
!
!      =======================================================================
!      == ORTHOGONALIZE UPHI TO THE CORE STATES                             ==
!      =======================================================================
       AEPHI_NSPH(:,:,:,:)=(0.D0,0d0)
       TAEPHI_NSPH(:,:,:,:)=(0.D0,0.d0)
       sAEPHI_NSPH(:,:,:,:)=(0.D0,0.d0)
       TsAEPHI_NSPH(:,:,:,:)=(0.D0,0.d0)
       LMN=0
       DO LN=1,LNX
         L=LOX(LN)
         LM=L**2
         DO M=1,2*L+1
           LMN=LMN+1
           LM=LM+1
           DO IS=1,2
             IV=(IS-1)*LMNX+LMN
             AEPHI_NSPH(:,LM,IS,IV)=UPHI(:,LN)
             TAEPHI_NSPH(:,LM,IS,IV)=TUPHI(:,LN)
           ENDDO
         ENDDO
       ENDDO

       DO IV=1,NV
         DO IC=1,NC
           AEPHI_NSPH(:,:,:,IV)  =  AEPHI_NSPH(:,:,:,IV)-  PHIC(:,:,:,IC)*S(IC,IV)
           TAEPHI_NSPH(:,:,:,IV) = TAEPHI_NSPH(:,:,:,IV)- TPHIC(:,:,:,IC)*S(IC,IV)
           sAEPHI_NSPH(:,:,:,IV) = sAEPHI_NSPH(:,:,:,IV)- sPHIC(:,:,:,IC)*S(IC,IV)
           tsAEPHI_NSPH(:,:,:,IV)=tsAEPHI_NSPH(:,:,:,IV)-tsPHIC(:,:,:,IC)*S(IC,IV)
         ENDDO
       ENDDO
      
       RETURN
       END
!
!     ...................................................................
      SUBROUTINE SF_UPDATEPARTWAVES_NSPH(GID,NR,lmx,LMNX,NDIMD,DENMAT &
     &            ,aephi,taephi,saephi,tsaephi,lmrx,pot,nc,phic,sphic)
!     **                                                               **
!     **  VALENCE TOTAL ENERGY FROM GIVEN PARTIAL WAVES                **
!     **                                                               **
!     **  REMARK: ONLY A SPHERICAL POTENTIAL IS CONSTRUCTED            **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LmX
      INTEGER(4),INTENT(IN) :: LmrX
      INTEGER(4),INTENT(IN) :: LmNX
      INTEGER(4),INTENT(IN) :: ndimd
      complex(8),intent(in) :: denmat(lmnx,lmnx,ndimd)
      complex(8),intent(in) :: aephi(nr,lmx,2,2*lmnx)
      complex(8),intent(in) :: taephi(nr,lmx,2,2*lmnx)
      complex(8),intent(in) :: saephi(nr,lmx,2,2*lmnx)
      complex(8),intent(in) :: tsaephi(nr,lmx,2,2*lmnx)
      real(8)   ,intent(in) :: pot(nr,lmrx,ndimd)
      integer(4),intent(in) :: nc
      complex(8),intent(in) :: phic(nr,lmx,2,nc)
      complex(8),intent(in) :: sphic(nr,lmx,2,nc)
      complex(8)            :: theta(2*lmnx,2*lmnx)
      complex(8)            :: pot1(nr,lmrx,2,2)
      complex(8)            :: daephi(nr,lmx,2,2*lmnx)
      complex(8)            :: dsaephi(nr,lmx,2,2*lmnx)
      complex(8)            :: faephi(nr,lmx,2,2*lmnx)
      complex(8)            :: fsaephi(nr,lmx,2,2*lmnx)
      integer(4)            :: lmn,lm1,lm2,lmr,is1,is2,lmn1,lmn2,ir,ic,is,lm
      real(8)               :: cg
      complex(8),parameter  :: ci=(0.d0,1.d0)
      character(64)         :: filename
      complex(8)            :: s(nc)
      complex(8)            :: caux(nr),csvar1,csvar2
      real(8)               :: aux(nr),svar1,svar2
      real(8)               :: r(nr)
      integer(4)            :: irc
!     *******************************************************************
      call radial$r(gid,nr,r)
!
!      =======================================================================
!      == TRANSFORM DENSITY MATRIX IN COMPLEX SPINOR FORM                   ==
!      == use A=0.5*sum_{i=0}^3 sigma_i*Tr[sigma_i*A]                       ==
!      == where sigma_0=1 and for i>0 sigma_i are the pauli matrices        ==
!      =======================================================================
       THETA(:,:)=(0.D0,0.D0)
       IF(NDIMD.EQ.1) THEN
         THETA(:LMNX,:LMNX)    =0.5D0*DENMAT(:,:,1)
         THETA(LMNX+1:,LMNX+1:)=0.5D0*DENMAT(:,:,1)
       ELSE IF(NDIMD.EQ.2) THEN
         THETA(:LMNX,:LMNX)    =0.5D0*(DENMAT(:,:,1)+DENMAT(:,:,2))
         THETA(LMNX+1:,LMNX+1:)=0.5D0*(DENMAT(:,:,1)-DENMAT(:,:,2))
       ELSE IF(NDIMD.EQ.4) THEN
         THETA(:LMNX,:LMNX)    =0.5D0*(DENMAT(:,:,1)+DENMAT(:,:,4))
         THETA(LMNX+1:,LMNX+1:)=0.5D0*(DENMAT(:,:,1)-DENMAT(:,:,4))
         THETA(:LMNX,LMNX+1:)  =0.5D0*(DENMAT(:,:,2)-CI*DENMAT(:,:,3))
         THETA(LMNX+1:,:LMNX)  =0.5D0*(DENMAT(:,:,2)+CI*DENMAT(:,:,3))
       END IF
!
!      =======================================================================
!      == expand potential                                                  ==
!      =======================================================================
       pot1(:,:,:,:)=(0.d0,0.d0)
       pot1(:,:,1,1)=pot1(:,:,1,1)+pot(:,:,1)
       pot1(:,:,2,2)=pot1(:,:,2,2)+pot(:,:,1)
       if(ndimd.eq.2) then
         pot1(:,:,1,1)=pot1(:,:,1,1)+pot(:,:,2)
         pot1(:,:,2,2)=pot1(:,:,2,2)-pot(:,:,2)
       else if(ndimd.eq.4) then
         pot1(:,:,1,2)=pot1(:,:,1,2)+pot(:,:,2)-ci*pot(:,:,3)
         pot1(:,:,2,1)=pot1(:,:,2,1)+pot(:,:,2)+ci*pot(:,:,3)
         pot1(:,:,1,1)=pot1(:,:,1,1)+pot(:,:,4)
         pot1(:,:,2,2)=pot1(:,:,2,2)-pot(:,:,4)
       end if
!
!      =======================================================================
!      == TRANSFORM gradient of partial waves                               ==
!      =======================================================================
      daephi(:,:,:,:)=taephi(:,:,:,:)
      dsaephi(:,:,:,:)=tsaephi(:,:,:,:)
      do lmn=1,2*lmnx
        do is1=1,2
          do lm1=1,lmx
            do is2=1,2
              do lm2=1,lmx
                do lmr=1,lmrx
                  call clebsch(lm1,lm2,lmr,cg)
                  if(cg.eq.0.d0) cycle
                  daephi(:,lm1,is1,lmn)=daephi(:,lm1,is1,lmn) &
     &                           +cg*pot1(:,lmr,is1,is2)*aephi(:,lm2,is2,lmn)
                  dsaephi(:,lm1,is1,lmn)=dsaephi(:,lm1,is1,lmn) &
     &                           +cg*pot1(:,lmr,is1,is2)*saephi(:,lm2,is2,lmn)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
!      =======================================================================
!      == multiply with theta                                              ==
!      =======================================================================
      faephi(:,:,:,:)=0.d0
      fsaephi(:,:,:,:)=0.d0
      do lmn1=1,2*lmnx
        do lmn2=1,2*lmnx
          faephi(:,:,:,lmn1)=faephi(:,:,:,lmn1)+daephi(:,:,:,lmn2)*theta(lmn2,lmn1)
          fsaephi(:,:,:,lmn1)=fsaephi(:,:,:,lmn1)+dsaephi(:,:,:,lmn2)*theta(lmn2,lmn1)
        enddo
      enddo
      daephi=faephi
      dsaephi=fsaephi
!
!      =======================================================================
!      == cancel tail by mixing aephi                                       ==
!      =======================================================================
do ir=1,nr
  if(r(ir).gt.2.5d0) exit
  irc=ir
enddo
       do lmn=1,lmnx
         csvar1=sum(conjg(faephi(irc,:,:,lmn))*aephi(irc,:,:,lmn)) &
     &         +sum(conjg(fsaephi(irc,:,:,lmn))*saephi(irc,:,:,lmn))         
         csvar2=sum(conjg(aephi(irc,:,:,lmn))*aephi(irc,:,:,lmn)) &
     &         +sum(conjg(saephi(irc,:,:,lmn))*saephi(irc,:,:,lmn))
         csvar1=csvar1/csvar2
         faephi(:,:,:,lmn)=faephi(:,:,:,lmn)-aephi(:,:,:,lmn)*csvar1
       enddo
!
!      =======================================================================
!      == orthogonalize to core                                             ==
!      =======================================================================
      do lmn=1,lmnx
        do ic=1,nc
          caux(:)=0.d0
          do lm=1,lmx
            do is=1,2
              caux(:)=caux(:)+conjg(phic(:,lm,is,ic))*faephi(:,lm,is,lmn) &
     &                      +conjg(sphic(:,lm,is,ic))*fsaephi(:,lm,is,lmn)
            enddo
          enddo
          caux(:)=caux(:)*r(:)**2
          aux(:)=real(caux,kind=8)
          call radial$integral(gid,nr,aux,svar1)    
          aux(:)=aimag(caux)
          call radial$integral(gid,nr,aux,svar2)    
          s(ic)=cmplx(svar1,svar2,kind=8)
        enddo
        do ic=1,nc
          faephi(:,:,:,lmn)=faephi(:,:,:,lmn)-phic(:,:,:,ic)*s(ic)
          fsaephi(:,:,:,lmn)=fsaephi(:,:,:,lmn)-sphic(:,:,:,ic)*s(ic)
        enddo
      enddo
!
!      =======================================================================
!      == multiply with theta                                              ==
!      =======================================================================
      do lmn1=1,2*lmnx
        write(filename,*)lmn1
        filename=adjustl(filename)
        filename=trim('daephi'//trim(filename))//".dat"
        open(unit=8,file=filename,form='formatted')
        do ir=1,nr
          write(8,fmt='(100f20.5)')r(ir),real(daephi(ir,:,:,lmn1)),aimag(daephi(ir,:,:,lmn1))
        enddo
        close(8)
      enddo
      do lmn1=1,2*lmnx
        write(filename,*)lmn1
        filename=adjustl(filename)
        filename=trim('faephi'//trim(filename))//".dat"
        open(unit=8,file=filename,form='formatted')
        do ir=1,nr
          write(8,fmt='(100f20.5)')r(ir),real(faephi(ir,:,:,lmn1)),aimag(faephi(ir,:,:,lmn1))
        enddo
        close(8)
      enddo
      do lmn1=1,2*lmnx
        write(filename,*)lmn1
        filename=adjustl(filename)
        filename=trim('fsaephi'//trim(filename))//".dat"
        open(unit=8,file=filename,form='formatted')
        do ir=1,nr
          write(8,fmt='(100f20.5)')r(ir),real(fsaephi(ir,:,:,lmn1)),aimag(fsaephi(ir,:,:,lmn1))
        enddo
        close(8)
      enddo
      do lmn1=1,2*lmnx
        write(filename,*)lmn1
        filename=adjustl(filename)
        filename=trim('aephi'//trim(filename))//".dat"
        open(unit=8,file=filename,form='formatted')
        do ir=1,nr
          write(8,fmt='(100f20.5)')r(ir),real(aephi(ir,:,:,lmn1)),aimag(aephi(ir,:,:,lmn1))
        enddo
        close(8)
      enddo
      do lmn1=1,2*lmnx
        write(filename,*)lmn1
        filename=adjustl(filename)
        filename=trim('taephi'//trim(filename))//".dat"
        open(unit=8,file=filename,form='formatted')
        do ir=1,nr
          write(8,fmt='(100f20.5)')r(ir),real(taephi(ir,:,:,lmn1)),aimag(taephi(ir,:,:,lmn1))
        enddo
        close(8)
      enddo
CALL SF_WRITE('pot.dat',GID,NR,LMRX*NDIMD,pot)
CALL SF_WRITE('pot1r.dat',GID,NR,LMRX*4,real(pot1))
CALL SF_WRITE('pot1i.dat',GID,NR,LMRX*4,aimag(pot1))
stop 'done'
      return
      end
!
!     ...................................................................
      SUBROUTINE SF_RHOVALENCE_SPH(GID,NR,LNX,LOX,LMNX,NDIMD,DENMAT &
     &                            ,AEPHI,TAEPHI,LMRX,TKIN,RHOV)
!     **                                                               **
!     **  VALENCE TOTAL ENERGY FROM GIVEN PARTIAL WAVES                **
!     **                                                               **
!     **  REMARK: ONLY A SPHERICAL POTENTIAL IS CONSTRUCTED            **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNX
      INTEGER(4),INTENT(IN) :: NDIMD
      complex(8),INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: TAEPHI(NR,LNX)
      INTEGER(4),INTENT(IN) :: LMRX
      REAL(8)   ,INTENT(OUT):: TKIN(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT):: RHOV(NR,LMRX,NDIMD)
      REAL(8)               :: TKIN1(LNX,LNX)
      REAL(8)               :: AUX(NR),SVAR
      REAL(8)               :: R(NR)
      INTEGER(4)            :: LN1,LN2,L1,L2,IM1,IM2,LMN1,LMN2,IDIM
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     =================================================================
!     == KINETIC ENERGY MATRIX ELEMENTS                              ==
!     =================================================================
      TKIN1(:,:)=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=AEPHI(:,LN1)*TAEPHI(:,LN2)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          TKIN1(LN1,LN2)=SVAR
        ENDDO
      ENDDO
!          
!     == EXPAND
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,TKIN1,TKIN)
!
!     =================================================================
!     == HARTREE ENERGY AND POTENTIAL                                ==
!     =================================================================
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,AEPHI &
     &                  ,LMNX,DENMAT(:,:,IDIM),LMRX,RHOV(:,:,IDIM))
      ENDDO
      RETURN
      END
!
!      ...................................................................
       SUBROUTINE SF_RHOVALENCE_NSPH(GID,NR,LMNX,NDIMD,DENMAT,LMX &
      &                             ,AEPHI,TAEPHI,saephi,tsaephi,LMRX,TKIN,RHO)
!      **                                                               **
!      **  CALCULATES THE TOTAL DENSITY ASSUMING A NONCOLLINEAR MODEL   **
!      **  WITH COMPLEX WAVE FUNCTIONS AND A COMPLEX DENSITY MATRIX     **
!      **                                                               **
!      **  THE NODE-LESS PARTIAL WAVES ARE STORED REAL                  **
!      **                                                               **
!      **  ASSUMES THAT THE CORE STATES ARE ORTHOGONAL                  **
!      **                                                               **
USE PERIODICTABLE_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: GID
       INTEGER(4),INTENT(IN) :: NR
       INTEGER(4),INTENT(IN) :: LMNX
       INTEGER(4),INTENT(IN) :: NDIMD
       complex(8),INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD)
       INTEGER(4),INTENT(IN) :: LMX   !#(ANGULAR MOMENTA FOR WAVE FUNCTIONS)
       COMPLEX(8),INTENT(IN) :: AEPHI(NR,LMX,2,2*LMNX)
       COMPLEX(8),INTENT(IN) :: TAEPHI(NR,LMX,2,2*LMNX)
       COMPLEX(8),INTENT(IN) :: sAEPHI(NR,LMX,2,2*LMNX)
       COMPLEX(8),INTENT(IN) :: TsAEPHI(NR,LMX,2,2*LMNX)
       INTEGER(4),INTENT(IN) :: LMRX
       REAL(8)   ,INTENT(OUT):: TKIN(LMNX,LMNX,NDIMD)
       REAL(8)   ,INTENT(OUT):: RHO(NR,LMRX,NDIMD)
       COMPLEX(8)            :: PHIVTHETA(NR,LMX,2,2*LMNX)
       COMPLEX(8)            :: sPHIVTHETA(NR,LMX,2,2*LMNX)
       COMPLEX(8)            :: THETA(2*LMNX,2*LMNX)
       COMPLEX(8)            :: CRHO(NR,LMRX,2,2)
       INTEGER(4)            :: NV ! 2*LMNX=#(PARTIAL WAVES) EXPANDED
       COMPLEX(8)            :: CAUX(NR)
       REAL(8)               :: AUX(NR),SVAR,SVAR1,SVAR2
       INTEGER(4)            :: IV1,IV2,IV,LMR,LM1,LM2,IR,IS,LM,LMN1,LMN2
       INTEGER(4)            :: IV1UP,IV1DOWN,IV2UP,IV2DOWN,IS1,IS2
       REAL(8)               :: R(NR)      ! RADIAL GRID
       COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  ! SQRT(-1)
       REAL(8)               :: CG ! GAUNT COEFFICIENT
       COMPLEX(8)            :: CTKIN(2*LMNX,2*LMNX)
REAL(8)       :: PI,Y0
!      ***********************************************************************
PI=4.D0*DATAN(1.D0)
Y0=1.D0/SQRT(4.D0*PI)
       CALL RADIAL$R(GID,NR,R)
       NV=2*LMNX
!
!      =======================================================================
!      == TRANSFORM DENSITY MATRIX IN COMPLEX SPINOR FORM                   ==
!      == use A=0.5*sum_{i=0}^3 sigma_i*Tr[sigma_i*A]                       ==
!      == where sigma_0=1 and for i>0 sigma_i are the pauli matrices        ==
!      =======================================================================
       THETA(:,:)=(0.D0,0.D0)
       IF(NDIMD.EQ.1) THEN
         THETA(:LMNX,:LMNX)    =0.5D0*DENMAT(:,:,1)
         THETA(LMNX+1:,LMNX+1:)=0.5D0*DENMAT(:,:,1)
       ELSE IF(NDIMD.EQ.2) THEN
         THETA(:LMNX,:LMNX)    =0.5D0*(DENMAT(:,:,1)+DENMAT(:,:,2))
         THETA(LMNX+1:,LMNX+1:)=0.5D0*(DENMAT(:,:,1)-DENMAT(:,:,2))
       ELSE IF(NDIMD.EQ.4) THEN
         THETA(:LMNX,:LMNX)    =0.5D0*(DENMAT(:,:,1)+DENMAT(:,:,4))
         THETA(LMNX+1:,LMNX+1:)=0.5D0*(DENMAT(:,:,1)-DENMAT(:,:,4))
         THETA(:LMNX,LMNX+1:)  =0.5D0*(DENMAT(:,:,2)-CI*DENMAT(:,:,3))
         THETA(LMNX+1:,:LMNX)  =0.5D0*(DENMAT(:,:,2)+CI*DENMAT(:,:,3))
       END IF
!
!      =======================================================================
!      == KINETIC ENERGY                                                    ==
!      =======================================================================
       DO IV1=1,NV
         DO IV2=1,NV
           CAUX=0.D0
           DO IS=1,2
             DO LM=1,LMX
               CAUX=CAUX+CONJG(AEPHI(:,LM,IS,IV1))*TAEPHI(:,LM,IS,IV2) &
     &                  +CONJG(SAEPHI(:,LM,IS,IV1))*TSAEPHI(:,LM,IS,IV2)
             ENDDO
           ENDDO
           CAUX=CAUX*R**2
           AUX=REAL(CAUX)
           CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
           AUX=AIMAG(CAUX)
           CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
           CTKIN(IV1,IV2)=CMPLX(SVAR1,SVAR2)
         ENDDO
       ENDDO
       DO LMN1=1,LMNX
         DO LMN2=1,LMNX
           IV1UP=LMN1
           IV1DOWN=LMNX+LMN1
           IV2UP=LMN2
           IV2DOWN=LMNX+LMN2
           TKIN(LMN1,LMN2,1)=0.5D0*REAL(CTKIN(IV1UP,IV2UP)+CTKIN(IV1DOWN,IV2DOWN))
           IF(NDIMD.EQ.2) THEN
             TKIN(LMN1,LMN2,2)=0.5D0*REAL(CTKIN(IV1UP,IV2UP)-CTKIN(IV1DOWN,IV2DOWN))
           ELSE IF(NDIMD.EQ.4) THEN
             TKIN(LMN1,LMN2,2)=0.5D0*REAL(CTKIN(IV1UP,IV2DOWN)+CTKIN(IV1DOWN,IV2UP))
             TKIN(LMN1,LMN2,3)=0.5D0*AIMAG(CTKIN(IV1UP,IV2DOWN)-CTKIN(IV1DOWN,IV2UP))
             TKIN(LMN1,LMN2,4)=0.5D0*REAL(CTKIN(IV1UP,IV2UP)-CTKIN(IV1DOWN,IV2DOWN))
           END IF             
         ENDDO
       ENDDO
!
!      =======================================================================
!      == CALCULATE VALENCE DENSITY                                         ==
!      =======================================================================
       PHIVTHETA(:,:,:,:)=0.D0
       sPHIVTHETA(:,:,:,:)=0.D0
       DO IV1=1,NV
         DO IV2=1,NV
           PHIVTHETA(:,:,:,IV1)=PHIVTHETA(:,:,:,IV1) &
      &                         +AEPHI(:,:,:,IV2)*THETA(IV2,IV1)
           sPHIVTHETA(:,:,:,IV1)=sPHIVTHETA(:,:,:,IV1) &
      &                         +sAEPHI(:,:,:,IV2)*THETA(IV2,IV1)
         ENDDO
       ENDDO
!
       CRHO=(0.D0,0.D0)
       DO IV=1,NV
         DO LMR=1,LMRX
           DO LM1=1,LMX
             DO LM2=1,LMX
               CALL CLEBSCH(LMR,LM1,LM2,CG)
               IF(CG.EQ.0.D0) CYCLE
               DO IS1=1,2
                 DO IS2=1,2
                   CRHO(:,LMR,IS1,IS2)=CRHO(:,LMR,IS1,IS2) &
      &                   +CG*PHIVTHETA(:,LM1,IS1,IV)*CONJG(AEPHI(:,LM2,IS2,IV)) &
      &                   +CG*sPHIVTHETA(:,LM1,IS1,IV)*CONJG(sAEPHI(:,LM2,IS2,IV))
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDDO               
!
!      =======================================================================
!      == MAP DENSITY BACK                                                  ==
!      =======================================================================
       RHO(:,:,1)=REAL(CRHO(:,:,1,1)+CRHO(:,:,2,2))
       IF(NDIMD.EQ.2) THEN
         RHO(:,:,2)=REAL(CRHO(:,:,1,1)-CRHO(:,:,2,2))
       ELSE IF(NDIMD.EQ.4) THEN
         RHO(:,:,2)=REAL(CRHO(:,:,1,2)+CRHO(:,:,2,1))
         RHO(:,:,3)=AIMAG(CRHO(:,:,1,2)-CRHO(:,:,2,1))
         RHO(:,:,4)=REAL(CRHO(:,:,1,1)-CRHO(:,:,2,2))
       END IF
!
!      =======================================================================
!      == WRITE WAVE FUNCTIONS                                              ==
!      =======================================================================
AUX(:)=RHO(:,1,1)*R(:)**2
CALL PERIODICTABLE$GET(80,'R(ASA)',SVAR)
PRINT*,'RASA',SVAR
DO IR=1,NR
  IF(R(IR).GT.SVAR)AUX(IR)=0.D0
ENDDO
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'VALENCE CHARGE ',4.D0*PI*SVAR*Y0
       RETURN
       END
!
!     ...................................................................
      SUBROUTINE SF_EVALENCE_SPH(GID,NR,LMNX,LMRX,NDIMD,RHOC,RHOV,DENMAT &
     &                          ,RHOB,VQLM,TKIN,ETOT,POT)
!     **                                                               **
!     **  VALENCE TOTAL ENERGY FROM GIVEN PARTIAL WAVES                **
!     **                                                               **
!     **  REMARK: ONLY A SPHERICAL POTENTIAL IS CONSTRUCTED            **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMNX
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: RHOC(NR)
      REAL(8)   ,INTENT(IN) :: RHOV(NR,LMRX,NDIMD)
      complex(8),INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(IN) :: RHOB
      REAL(8)   ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)   ,INTENT(IN) :: TKIN(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: POT(NR,LMRX,NDIMD)
      REAL(8)               :: POTH(NR,LMRX) ! HARTREE POTENTIAL
      REAL(8)               :: RHOT(NR,LMRX,NDIMD)
      REAL(8)               :: AUX(NR),SVAR
      REAL(8)               :: R(NR)
      INTEGER(4)            :: LN1,LN2,L1,L2,IM1,IM2,LMN1,LMN2,IDIMD
      REAL(8)               :: EKIN,EH,COREEXC,EXC
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     =================================================================
!     == KINETIC ENERGY                                              ==
!     =================================================================
      EKIN=0.D0
      DO IDIMD=1,NDIMD
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            EKIN=EKIN+real(DENMAT(LMN1,LMN2,IDIMD)*TKIN(LMN2,LMN1,IDIMD),kind=8)
          ENDDO
        ENDDO
      ENDDO
!
!     ================================================================
!     ==   ADD EXCHANGE AND CORRELATION POTENTIAL                   ==
!     ================================================================
!     == CORE ONLY EXCHANGE ENERGY ===================================
      CALL AUGMENTATION_XC(GID,NR,1,1,RHOC,COREEXC,POT(:,1,1))
      POT=0.D0
!     == AE-EXCHANGE ENERGY AND POTENTIAL ============================
      RHOT=RHOV
      RHOT(:,1,1)=RHOT(:,1,1)+RHOC(:)
      CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOT,EXC,POT)
      POT(1,:,:)=POT(2,:,:)
      EXC=EXC-COREEXC
!
!     =================================================================
!     == HARTREE ENERGY AND POTENTIAL                                ==
!     =================================================================
      CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,RHOC,RHOV(:,:,1) &
     &                               ,VQLM,RHOB,POTH,EH)
      POT(:,:,1)=POT(:,:,1)+POTH(:,:)
!
!     ================================================================
!     ==   ADD UP POTENTIALS                                        ==
!     ================================================================
      ETOT=EKIN+EH+EXC
PRINT*,'EKIN    ',EKIN
PRINT*,'EXC     ',EXC,COREEXC
PRINT*,'EHARTREE',EH
PRINT*,'EVALENCE',ETOT
      RETURN
      END
!
!      ...................................................................
       SUBROUTINE SF_EVALENCE_NSPH(GID,NR,LMNX,LMRX,NDIMD,RHOC,RHOV,DENMAT &
      &                           ,RHOB,VQLM,TKIN,ETOT,AEPOT)
!      **                                                               **
!      **  CALCULATES THE TOTAL DENSITY ASSUMING A NONCOLLINEAR MODEL   **
!      **  WITH COMPLEX WAVE FUNCTIONS AND A COMPLEX DENSITY MATRIX     **
!      **                                                               **
!      **  THE NODE-LESS PARTIAL WAVES ARE STORED REAL                  **
!      **                                                               **
!      **  ASSUMES THAT THE CORE STATES ARE ORTHOGONAL                  **
!      **                                                               **
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: GID
       INTEGER(4),INTENT(IN) :: NR
       INTEGER(4),INTENT(IN) :: LMNX
       INTEGER(4),INTENT(IN) :: LMRX
       INTEGER(4),INTENT(IN) :: NDIMD
       REAL(8)   ,INTENT(IN) :: RHOC(NR,LMRX,NDIMD)
       REAL(8)   ,INTENT(IN) :: RHOV(NR,LMRX,NDIMD)
       complex(8),INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD)
       REAL(8)   ,INTENT(IN) :: RHOB
       REAL(8)   ,INTENT(IN) :: VQLM(LMRX)
       REAL(8)   ,INTENT(IN) :: TKIN(LMNX,LMNX,NDIMD)
       REAL(8)   ,INTENT(OUT):: ETOT
       REAL(8)   ,INTENT(OUT):: AEPOT(NR,LMRX,NDIMD)
       REAL(8)               :: R(NR)
       REAL(8)               :: AUX(NR)
       INTEGER(4)            :: IDIMD,LMR,LMN1,LMN2
       REAL(8)               :: EKIN,EXC,EXCCORE,EH,EHCORE
       REAL(8)               :: HPOT(NR,LMRX)
       real(8)               :: VQLM0(LMRX),RHOB0
!      ***********************************************************************
       CALL RADIAL$R(GID,NR,R)
!
!      =======================================================================
!      == KINETIC ENERGY                                                    ==
!      =======================================================================
       EKIN=0.D0
       DO IDIMD=1,NDIMD
         DO LMN1=1,LMNX
           DO LMN2=1,LMNX
             EKIN=EKIN+real(DENMAT(LMN1,LMN2,IDIMD)*TKIN(LMN2,LMN1,IDIMD),kind=8)
           ENDDO
         ENDDO
       ENDDO
!
!      =======================================================================
!      == XC ENERGY                                                         ==
!      =======================================================================
       CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOC,EXCCORE,AEPOT)
       AEPOT=0.D0
       CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOV+RHOC,EXC,AEPOT)
       EXC=EXC-EXCCORE
!
!      =======================================================================
!      == HARTREE ENERGY                                                    ==
!      =======================================================================
       AUX(:)=0.D0  ! CORE IS SET TO ZERO, DENSITY CONTAINS CORE
       RHOB0=0.D0
       VQLM0(:)=0.D0
       CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AUX,RHOC(:,:,1) &
      &                          ,VQLM0,RHOB0,HPOT,EHCORE)
       HPOT=0.D0
       CALL AUGMENTATION_AEHARTREE(GID,NR,LMRX,AUX,RHOV(:,:,1)+RHOC(:,:,1) &
      &                          ,VQLM,RHOB,HPOT,EH)
       EH=EH-EHCORE
       AEPOT(:,:,1)=AEPOT(:,:,1)+HPOT
!
!      =======================================================================
!      == ADD UP ENERGIES AND POTENTIALS                                    ==
!      =======================================================================
       ETOT=EKIN+EH+EXC
PRINT*,'EKIN    ',EKIN
PRINT*,'EXC     ',EXC,EXCCORE
PRINT*,'EHARTREE',EH,EHCORE
PRINT*,'EVALENCE',ETOT
       RETURN
       END
!
!     ...................................................................
      SUBROUTINE SF_POTMATRIX_SPH(GID,NR,LMRX,NDIMD,LNX,LOX,LMNX,AEPOT,AEPHI,HAM,OV)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: HAM(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT):: OV(LMNX,LMNX,NDIMD)
      REAL(8)               :: R(NR)
      INTEGER(4)            :: LN1,LN2,LMN1,LMN2,L1,L2,LM1,LM2,IM1,IM2,LMR,IDIM
      REAL(8)               :: CG ! GAUNT COEFFICIENT
      COMPLEX(8)            :: VMAT(2*LMNX,2*LMNX)
      REAL(8)               :: AUX(NR),SVAR
      REAL(8)               :: OV1(LNX,LNX)
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ====================================================================
!     ==  POTENTIAL TIMES PARTIAL WAVE                                  ==
!     ====================================================================
      HAM=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LM1=L1**2
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LM1=LM1+1
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            LM2=L2**2
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=LM2+1
              DO IDIM=1,NDIMD
                AUX(:)=0.D0
                DO LMR=1,LMRX
                  CALL CLEBSCH(LM1,LM2,LMR,CG)
                  IF(CG.EQ.0.D0) CYCLE
                  AUX(:)=AUX(:)+CG*AEPOT(:,LMR,IDIM)*AEPHI(:,LN1)*AEPHI(:,LN2)
                ENDDO
                AUX(:)=AUX(:)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
                HAM(LMN1,LMN2,IDIM)=SVAR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ====================================================================
!     ==  OVERLAP MATRIX                                                ==
!     ====================================================================
      OV1=0.D0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.NE.L2) CYCLE
          AUX=AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          OV1(LN1,LN2)=SVAR
        ENDDO
      ENDDO
      CALL AUGMENTATION_EXPANDDA(LNX,LMNX,LOX,NDIMD,OV1,OV)
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE SF_POTMATRIX_NSPH(GID,NR,LMRX,NDIMD,LMX,LMNX,AEPOT,AEPHI,HAM,OV)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
      INTEGER(4),INTENT(IN) :: LMX
      INTEGER(4),INTENT(IN) :: LMNX
      COMPLEX(8),INTENT(IN) :: AEPHI(NR,LMX,2,2*LMNX)
      REAL(8)   ,INTENT(OUT):: HAM(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT):: OV(LMNX,LMNX,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: R(NR)
      COMPLEX(8)            :: POT(NR,LMRX,2,2)
      REAL(8)               :: VAEPHI(NR,LMX,2,2*LMNX)
      INTEGER(4)            :: LMN1,LMN2,LM,LM1,LM2,LMR,IS,IS1,IS2
      INTEGER(4)            :: LMN1UP,LMN1DN,LMN2UP,LMN2DN
      REAL(8)               :: CG ! GAUNT COEFFICIENT
      COMPLEX(8)            :: CAUX(NR)
      COMPLEX(8)            :: VMAT(2*LMNX,2*LMNX),OMAT(2*LMNX,2*LMNX)
      REAL(8)               :: AUX(NR),SVAR1,SVAR2
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ====================================================================
!     ==  EXPAND POTENTIAL                                              ==
!     ====================================================================
      POT(:,:,:,:)=(0.D0,0.D0)
      POT(:,:,1,1)=0.5D0*AEPOT(:,:,1)
      POT(:,:,2,2)=0.5D0*AEPOT(:,:,1)
      IF(NDIMD.EQ.2) THEN
        POT(:,:,1,1)=POT(:,:,1,1)+0.5D0*AEPOT(:,:,2)
        POT(:,:,2,2)=POT(:,:,2,2)-0.5D0*AEPOT(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        POT(:,:,1,2)=POT(:,:,1,2)+0.5D0*(AEPOT(:,:,2)-CI*AEPOT(:,:,3))
        POT(:,:,2,1)=POT(:,:,2,1)+0.5D0*(AEPOT(:,:,2)+CI*AEPOT(:,:,3))
        POT(:,:,1,1)=POT(:,:,1,1)+0.5D0*AEPOT(:,:,4)
        POT(:,:,2,2)=POT(:,:,2,2)-0.5D0*AEPOT(:,:,4)
      END IF
      POT=POT*2.D0
!
!     ====================================================================
!     ==  POTENTIAL TIMES PARTIAL WAVE                                  ==
!     ====================================================================
      VAEPHI(:,:,:,:)=(0.D0,0.D0)
      DO LMN1=1,2*LMNX
        DO LM1=1,LMX
          DO LM2=1,LMX
            DO LMR=1,LMRX
              CALL CLEBSCH(LM1,LM2,LMR,CG)
              IF(CG.EQ.0.D0) CYCLE
              DO IS1=1,2
                DO IS2=1,2
                  VAEPHI(:,LM1,IS1,LMN1)=VAEPHI(:,LM1,IS1,LMN1) &
     &                        +CG*POT(:,LMR,IS1,IS2)*AEPHI(:,LM2,IS2,LMN1) 
                ENDDO
              ENDDO
            ENDDO
         ENDDO
        ENDDO
      ENDDO
!
!     ====================================================================
!     ==  MATRIX ELEMENTS                                               ==
!     ====================================================================
      VMAT(:,:)=0.D0
      OMAT(:,:)=0.D0
      DO LMN1=1,2*LMNX
        DO LMN2=1,2*LMNX
          CAUX(:)=(0.D0,0.D0)
          DO LM=1,LMX
            DO IS=1,2
              CAUX(:)=CAUX(:)+CONJG(AEPHI(:,LM,IS,LMN1))*VAEPHI(:,LM,IS,LMN2)
            ENDDO
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          VMAT(LMN1,LMN2)=CMPLX(SVAR1,SVAR2)
!         == NOW OVERLAP
          CAUX(:)=(0.D0,0.D0)
          DO LM=1,LMX
            DO IS=1,2
              CAUX(:)=CAUX(:)+CONJG(AEPHI(:,LM,IS,LMN1))*AEPHI(:,LM,IS,LMN2)
            ENDDO
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          OMAT(LMN1,LMN2)=CMPLX(SVAR1,SVAR2)
        ENDDO
      ENDDO
!
!     ====================================================================
!     ==  TRANSFORM BACK TO SPIN REPRESENTATION                         ==
!     ====================================================================
      DO LMN1=1,LMNX
        LMN1UP=LMN1
        LMN1DN=LMNX+LMN1
        DO LMN2=1,LMNX
          LMN2UP=LMN2
          LMN2DN=LMNX+LMN2
          HAM(LMN1,LMN2,1)=REAL(VMAT(LMN1UP,LMN2UP)+VMAT(LMN1DN,LMN2DN))
          OV(LMN1,LMN2,1)  =REAL(OMAT(LMN1UP,LMN2UP)+OMAT(LMN1DN,LMN2DN))
          IF(NDIMD.EQ.2) THEN
            HAM(LMN1,LMN2,2)=REAL(VMAT(LMN1UP,LMN2UP)-VMAT(LMN1DN,LMN2DN))
            OV(LMN1,LMN2,2)  =REAL(OMAT(LMN1UP,LMN2UP)-OMAT(LMN1DN,LMN2DN))
          ELSE IF(NDIMD.EQ.4) THEN
            HAM(LMN1,LMN2,2)=REAL(VMAT(LMN1UP,LMN2DN)+VMAT(LMN1DN,LMN2UP))
            HAM(LMN1,LMN2,3)=AIMAG(VMAT(LMN1UP,LMN2DN)-VMAT(LMN1DN,LMN2UP))
            HAM(LMN1,LMN2,4)=REAL(VMAT(LMN1UP,LMN2UP)-VMAT(LMN1DN,LMN2DN))
            OV(LMN1,LMN2,2)  =REAL(OMAT(LMN1UP,LMN2DN)+OMAT(LMN1DN,LMN2UP))
            OV(LMN1,LMN2,3)  =AIMAG(OMAT(LMN1UP,LMN2DN)-OMAT(LMN1DN,LMN2UP))
            OV(LMN1,LMN2,4)  =REAL(OMAT(LMN1UP,LMN2UP)-OMAT(LMN1DN,LMN2DN))
          END IF
        ENDDO
      ENDDO
      HAM=0.5D0*HAM
      OV=0.5D0*OV
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
      INTEGER(4)                 :: NN1M,NN10,NN1P
      REAL(8)                    :: PHIP(NR),PHIM(NR)
      REAL(8)                    :: RCL
      REAL(8)                    :: SWKB(NR)  ! PHI=E^S
      REAL(8)                    :: AUX(NR),AUX1(NR)
      INTEGER(4)                 :: IROUT,IRCL
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
        IROUT=NR
        IF(RTEST.LT.R(NR)) THEN
          AUX(:IRCL)=0.D0
          AUX(IRCL+1:)=SQRT(REAL(L*(L+1),KIND=8)/R(IRCL+1:)**2+2.D0*(POT(IRCL+1:)*Y0-E))
          CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,SWKB)
!          SWKB(:)=SWKB(:)-0.5D0*LOG(AUX1(:))-LOG(R(:))
          SWKB(:)=SWKB(:)-LOG(R(:))
!         == DETERMINE IROUT WHERE THE WAVE FUNCTION CAN GROW BY A FACTOR 
!         == OF XMAX FROM THE CLASSICAL TURNING POINT
          SVAR=LOG(XMAX)
          DO IR=1,NR
            IF(SWKB(IR).GT.SVAR) THEN
              IROUT=IR-1
              EXIT
            END IF
          ENDDO
        END IF
        SVAR=POT(IROUT)
        POT1(:)=POT(:)
        POT1(IROUT:)=POT(IROUT)
!PRINT*,'R(IROUT)',RCL,R(IROUT),POT1(NR)
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHRODINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IROUT).GT.0.OR.PHI(IROUT).LE.0)) THEN
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
        DO IR=3,IROUT
          IF(PHI(IR)*PHI(IR-1).LT.0.D0) NN1=NN1+1
        ENDDO
        Z0=-1.D0
        IF(NN1.GT.NN) Z0=1.D0
!WRITE(*,FMT='("ITER",4I5,4E20.10)')I,L,NN1,NN,E,2.D0*DX,PHI(IROUT),Z0
        IF(ABS(2.D0*DX).LE.TOL) EXIT
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
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E+2.D0*DX,IDIR,PHIP)
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E-2.D0*DX,IDIR,PHIM)
      NN1P=-NN
      NN1M=-NN
      NN10=-NN
      DO IR=3,IROUT
        IF(PHI(IR)*PHI(IR-1).LT.0.D0)   NN10=NN10+1
        IF(PHIM(IR)*PHIM(IR-1).LT.0.D0) NN1M=NN1M+1
        IF(PHIP(IR)*PHIP(IR-1).LT.0.D0) NN1P=NN1P+1
      ENDDO
!PRINT*,'NN1-NN ',NN1M,NN10,NN1P
!PRINT*,'E ',E-2.D0*DX,E,E+2.D0*DX
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      DO IR=1,NR
         IRMATCH=IR
         IF(POT(IR)-E/Y0.GT.0.D0.OR.R(IR).GT.5.D0) EXIT
      ENDDO
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IDIR=-1
      GHOM(:)=0.D0
      IF(IROUT.LT.NR) GHOM(IROUT)=1.D-5
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
DO IR=1,NR
  IF(.NOT.(PHI(IR).GT.0.D0.OR.PHI(IR).LE.0.D0)) THEN
    PRINT*,'ERROR'
    PRINT*,'PHIIN',PHI(:IRMATCH-1)
    PRINT*,'PHIOUT',PHI(IRMATCH:)
    PRINT*,'SVAR ',SVAR
    CALL ERROR$STOP('BOUNDSTATE')
  END IF
ENDDO
!
!     =======================================================================
!     ==  normalize solution                                               ==
!     =======================================================================
      call radial$integral(gid,nr,(phi(:)*r(:))**2,svar)
      phi(:)=phi(:)/sqrt(svar)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE OLD2BOUNDSTATE(GID,NR,L,SO,DREL,G,NN,POT,E,PHI)
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
      INTEGER(4)                 :: NN1M,NN10,NN1P
      REAL(8)                    ::PHIP(NR),PHIM(NR)
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
        IF(E.LT.POT(NR)*Y0) THEN
          DO IR=NR-1,1,-1
            IF(E.GT.POT(IR)*Y0) THEN
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
          POT1(:)=MIN(POT(:),(E+SVAR)/Y0)
!PRINT*,'POT CHANGED',SVAR,E
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
        IF(NN1.GT.NN) Z0=1.D0
!WRITE(*,FMT='("ITER",4I5,4E20.10)')I,L,NN1,NN,E,2.D0*DX,X0-XM,Z0
        IF(ABS(2.D0*DX).LE.TOL) EXIT
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
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E+2.D0*DX,IDIR,PHIP)
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E-2.D0*DX,IDIR,PHIM)
      NN1P=-NN
      NN1M=-NN
      NN10=-NN
      DO IR=3,NR-1
        IF(PHI(IR)*PHI(IR-1).LT.0.D0)   NN10=NN10+1
        IF(PHIM(IR)*PHIM(IR-1).LT.0.D0) NN1M=NN1M+1
        IF(PHIP(IR)*PHIP(IR-1).LT.0.D0) NN1P=NN1P+1
      ENDDO
!PRINT*,'NN1-NN ',NN1M,NN10,NN1P
!PRINT*,'E ',E-2.D0*DX,E,E+2.D0*DX
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      DO IR=1,NR
         IRMATCH=IR
         IF(POT(IR)-E/Y0.GT.0.D0.OR.R(IR).GT.5.D0) EXIT
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
!DO IR=1,NR
!  PRINT*,R(IR),PHI(IR),PHIINHOM(IR)
!ENDDO
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
      SUBROUTINE OLDBOUNDSTATE(GID,NR,L,SO,DREL,NN,POT,E,PHI)
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
      G(:)=0.D0
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
!..................................................................................
MODULE BROYDEN_MODULE
LOGICAL(4)         :: TON=.FALSE.
INTEGER(4)         :: NSTEPX=0
INTEGER(4)         :: NSTEP=0
INTEGER(4)         :: NX=0
REAL(8)            :: ALPHA
REAL(8),ALLOCATABLE :: XPREV(:,:)
REAL(8),ALLOCATABLE :: YPREV(:,:)
REAL(8),ALLOCATABLE :: WEIGHT(:)
END MODULE BROYDEN_MODULE
!      .............................................................................
       SUBROUTINE BROYDEN$NEW(NX_,NSTEPX_,ALPHA_)
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
       ALLOCATE(WEIGHT(NX))
       WEIGHT(:)=1.D0
       XPREV(:,:)=0.D0
       YPREV(:,:)=0.D0
       RETURN
       END
!      .............................................................................
       SUBROUTINE BROYDEN$CLEAR
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
       DEALLOCATE(WEIGHT)
       RETURN
       END
!      .............................................................................
       SUBROUTINE BROYDEN$SETWEIGHT(NX_,WEIGHT_)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NX_
       REAL(8)   ,INTENT(IN)  :: WEIGHT_(NX_)
!      *****************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$SETWEIGHT')
       END IF
       IF(NX_.NE.NX) THEN
         CALL ERROR$MSG('SIZE INCONSISTENT')
         CALL ERROR$STOP('BROYDEN$SETWEIGHT')
       END IF
       WEIGHT(:)=WEIGHT_(:)
       RETURN
       END
!      .............................................................................
       SUBROUTINE BROYDEN$STEP(NX_,X,Y)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NX_
       REAL(8)   ,INTENT(INOUT) :: X(NX_)
       REAL(8)   ,INTENT(IN)    :: Y(NX_)
       REAL(8)   ,ALLOCATABLE   :: DX(:,:)
       REAL(8)   ,ALLOCATABLE   :: DY(:,:)
       REAL(8)   ,ALLOCATABLE   :: B(:,:)
       REAL(8)   ,ALLOCATABLE   :: BINV(:,:)
       REAL(8)                  :: WY(NX_)
       INTEGER(4)               :: I
!      *****************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
       IF(NX_.NE.NX) THEN
         CALL ERROR$MSG('SIZE INCONSISTENT')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
!PRINT*,'NSTEP',NSTEP
!
!      =================================================================
!      == APPLY WEIGHTING                                             ==
!      =================================================================
       WY(:)=WEIGHT(:)*Y(:)
!
!      =================================================================
!      == SIMPLE MIXING IN THE FIRST STEP                             ==
!      =================================================================
       IF(NSTEP.EQ.0) THEN
         IF(NSTEPX.GT.0)THEN
           NSTEP=1
           XPREV(:,1)=X(:)     
           YPREV(:,1)=WY(:)     
         END IF
         X=X+ALPHA*WY
         RETURN
       END IF
!
!      =================================================================
!      == DETERMINE INVERSE HESSIAN ALPHA+DX OTIMES DY                ==
!      =================================================================
       ALLOCATE(DX(NX,NSTEP))
       ALLOCATE(DY(NX,NSTEP))
       DO I=1,NSTEP
         DY(:,I)=YPREV(:,I)-WY(:)  
         DX(:,I)=XPREV(:,I)-X(:)+ALPHA*DY(:,I)
       ENDDO
       ALLOCATE(B(NSTEP,NSTEP))
       ALLOCATE(BINV(NSTEP,NSTEP))
       B=MATMUL(TRANSPOSE(DY),DY)   !OVERLAP MATRIX OF DY
!PRINT*,'B',B
       CALL LIB$INVERTR8(NSTEP,B,BINV)           
!ALLOCATE(W(NX,NSTEP))
!W=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
!PRINT*,'W ',MATMUL(TRANSPOSE(W),DY)
!DEALLOCATE(W)
       DY=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
       DEALLOCATE(B)
       DEALLOCATE(BINV)
!
!      =================================================================
!      == STORE HISTORY                                               ==
!      =================================================================
       IF(NSTEP.LT.NSTEPX)NSTEP=NSTEP+1
       DO I=NSTEP,2,-1
         YPREV(:,I)=YPREV(:,I-1)
         XPREV(:,I)=XPREV(:,I-1)
       ENDDO
       XPREV(:,1)=X(:)     
       YPREV(:,1)=WY(:)     
!
!      =================================================================
!      == PREDICT NEW VECTOR                                          ==
!      =================================================================
       X=X+ALPHA*WY-MATMUL(DX,MATMUL(TRANSPOSE(DY),WY))
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
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: LMRX
      complex(8),INTENT(IN) :: DENMAT(LMNXX,LMNXX)
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: RHOL(NR,LMRX)
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
!                 == assumes real partial waves ==================
                  SVAR=CG*real(DENMAT(LMN1,LMN2))
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
      REAL(8)               :: DWORK(NR)
      REAL(8)               :: PI,Y0
      REAL(8)               :: RES
      INTEGER(4)            :: LM,L,IR
      REAL(8)               :: R(NR)
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      QLM(:)=0.D0
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
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
!     ..................................................................
      SUBROUTINE AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
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
!     ****************************************** P.E. BLOECHL, 1996 ****
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
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
      CALL RADIAL$R(GID,NR,R)
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
            CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
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
      CALL RADIAL$INTEGRAL(GID,NR,WORK1(:)*R(:)**2,EXC)
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
!
!     ==================================================================
!     ==   CORRECT FOR DIVERGENCE AT THE ORIGIN:                      ==
!     ==   IF A SHIFTED LOGARITHMIC GRID IS USED THE FIRST GRID POINT ==
!     ==   IS MESSED UP BECAUSE OF FACTORS 1/R                        ==
!     ==================================================================
      IF(R(1).LT.1.D-5) THEN
        VXC(1,:,:)=VXC(2,:,:)
      END IF
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
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX1(NR)
      REAL(8)                :: RHO1(NR)
      REAL(8)                :: POT(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: G(NR)
      REAL(8)                :: RHOHAT(NR,LMRX)
      REAL(8)                :: PI
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      INTEGER(4)             :: LM,L,M,IR,LX
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
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
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
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
        L=INT(SQRT(REAL(LM-1)+.01D0))
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
      SUBROUTINE AUGMENTATION_AEHARTREE(GID,NR,LMRX,RHOC,AERHO &
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
      REAL(8)    ,INTENT(IN) :: RHOC(NR)   ! CORE DENSITY
      REAL(8)    ,INTENT(IN) :: AERHO(NR,LMRX)
      REAL(8)    ,INTENT(IN) :: VQLM(LMRX)
      REAL(8)    ,INTENT(IN) :: RHOB       ! COMPENSATING BACKGROUND
      REAL(8)    ,INTENT(OUT):: AEPOT(NR,LMRX)
      REAL(8)    ,INTENT(OUT):: AEEH       ! ENERGY
      REAL(8)                :: R(NR)
      REAL(8)                :: AUX(NR)
      REAL(8)                :: POT(NR)
      REAL(8)                :: AEE(NR)  ! ENERGY DENSITY
      REAL(8)                :: SVAR
      REAL(8)                :: EEXT
      REAL(8)                :: AEZ
      REAL(8)                :: PI,Y0
      INTEGER(4)             :: LM
      INTEGER(4)             :: L
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
!AUX(:)=AEE(:)*R(:)**2
!CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!PRINT*,'EL C-V INT', SVAR
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
!AUX(:)=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1)+.01D0))
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
        L=INT(SQRT(REAL(LM-1)+.01D0))
        POT(:)=VQLM(LM)*R(:)**L
        AEPOT(:,LM)=AEPOT(:,LM)+POT(:)
!       == THE FOLLOWING ENERGY DENSITY CANCELS WITH THE PSEUDO TERM
!       == CAUTION HOWEVER FOR THE SOFT CORE
        AEE(:)=AEE(:)+AERHO(:,LM)*POT(:)
        IF(L.EQ.0)AEE(:)=AEE(:)+RHOC(:)*POT(:)
      ENDDO
      CALL SETUP$GETR8('AEZ',AEZ)
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
  CHARACTER(LEN=32):: TYPE   ! CAN BE 'S','P','D',F','ALL'
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
      complex(8)  ,INTENT(IN)   :: DENMAT(LMNX,LMNX,NDIMD)
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
      INTEGER(4)                :: GID     ! GRID ID
      REAL(8)      ,ALLOCATABLE :: R(:)    ! RADIAL GRID
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
        ELSE IF(TRIM(POT1%TYPE).EQ.'F') THEN
          LANG=3
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
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
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
            AUX(:)=RDEP(:)*AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,UONE(LN1,LN2))
          ENDDO
        ENDDO
        DEALLOCATE(R)
        DEALLOCATE(RDEP)
        DEALLOCATE(AUX)
!OPEN(41,FILE='DUMP',FORM='FORMATTED')
!DO IR=1,NR
!  WRITE(41,*)R(IR),(AEPHI(IR,LN1),LN1=1,LNX)
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
            ETOT=ETOT+real(DENMAT(LMN1,LMN2,IDIMD)*DATH(LMN2,LMN1,IDIMD),kind=8)
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
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
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
      REAL(8)                 :: PI
      REAL(8)                 :: SVAR
      REAL(8)                 :: WORK(NR)
      REAL(8)                 :: R(NR)
!     ******************************************************************
      IF(RHOB.EQ.0.D0) THEN
        EB=0.D0
        RETURN
      END IF
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      SVAR=-2.D0*PI*RHOB/3.D0*DSQRT(4.D0*PI)
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
      REAL(8)                :: PI
      REAL(8)                :: Y0
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
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
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
      CALL SETUP$GETR8A('NUCPOT',NR,AUX1)
      AEDMU(:)  = AEDMU(:)+AUX1(:)
!      AEDMU(:)  = AEDMU(:)-AEZ/R(:)/Y0
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
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        CALL GAUSSN(L,ALPHA,CL)
        SVAR=QLM(LM)*CL
        PSRHO(:,LM)=PSRHO(:,LM)+SVAR*DEXP(-ALPHA*R(:)**2)*R(:)**L
      ENDDO
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      AEE(:)=0.D0
      PSE(:)=0.D0
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        CALL RADIAL$POISSON(GID,NR,L,AERHO(1,LM),AEDMU)
        CALL RADIAL$POISSON(GID,NR,L,PSRHO(1,LM),PSDMU)
        ALPHA=1.D0/RCSM**2
        CALL GAUSSN(L,ALPHA,CL)
!       == ADD COMPENSATION DENSITY AND ITS POTENTIAL ================
        GAUSSIAN(:)=CL*DEXP(-ALPHA*R(:)**2)*R(:)**L
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



