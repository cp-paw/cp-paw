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
REAL(8)                :: val(NE)
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
        val(NE)=0.D0
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
      CALL MPE$COMBINE('+',VAL)
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
      LOGICAL(4),PARAMETER   :: TPR=.false.
      INTEGER(4),INTENT(IN)  :: ISP
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: LMRX
      REAL(8)   ,INTENT(IN)  :: DENMAT(LMNX,LMNX)
      REAL(8)   ,INTENT(OUT) :: QLM(LMRX)
      REAL(8)   ,ALLOCATABLE :: AERHO(:,:)   !(NRX,LMRX)
      REAL(8)   ,ALLOCATABLE :: PSRHO(:,:)   !(NRX,LMRX)
      REAL(8)   ,ALLOCATABLE :: AECORE(:)    !(NRX)
      REAL(8)   ,ALLOCATABLE :: PSCORE(:)    !(NRX)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)   !(NRX,LNX)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)   !(NRX,LNX)
      INTEGER(4),ALLOCATABLE :: LOX(:)       !(LNX)
      REAL(8)                :: AEZ
      REAL(8)                :: R1,DEX,XEXP
      INTEGER(4)             :: NR
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
!         WRITE(*,FMT='(9F10.6)')(DENMAT(LMN1,lmn2),lmn2=1,lmnx)
        ENDDO
      END IF
!     
!     ==================================================================
!     ==   RECEIVE PARTIAL WAVES FROM SETUP OBJECT                    ==
!     ==================================================================
      CALL SETUP$RADGRID(ISP,R1,DEX,NR)
      XEXP=DEXP(DEX)
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
      CALL AUGMENTATION_QLM(R1,DEX,NR,NR,LMRX &
     &                     ,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
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
      SUBROUTINE AUGMENTATION$SPHERE(ISP,IAT,LMNX,NDIMD,DENMAT &
     &                              ,LMRX,VQLM,RHOB,DATH,DO)
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
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(INOUT):: VQLM(LMRX)
      REAL(8)   ,INTENT(IN)   :: RHOB
      REAL(8)   ,INTENT(OUT)  :: DATH(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT)  :: DO(LMNX,LMNX,ndimd)
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
      REAL(8)   ,ALLOCATABLE  :: HAMUP1(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: HAMDWN1(:,:,:)
      REAL(8)   ,ALLOCATABLE  :: DWORK1(:)
      INTEGER(4)              :: IDIM,LMR,IR
      INTEGER(4)              :: LM,LMN1,LMN2,LMN11,LMN21,LN1,LN2,L1,L2,IM,ISPIN
      INTEGER(4)              :: NSPIN
      INTEGER(4)              :: NFILO
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=1
      LOGICAL(4)              :: TBACK,TSPIN
      REAL(8)                 :: DETOT,PSEHARTREE,AEEHARTREE,COREEXC,EKINNL,ENL,AEEXC,PSEXC,HAMUP,HAMDWN
      REAL(8)                 :: AEBACKGROUND,PSBACKGROUND,ELDA_U 
      CHARACTER(32)           :: ATOM
      REAL(8)                 :: VQLM1(LMRX)
      REAL(8)                 :: QLM(LMRX)
!     ******************************************************************
                            CALL TRACE$PUSH('AUGMENTATION$SPHERE')
!
!     ==================================================================
!     ==  COLLECT ATOM-TYPE SPECIFIC INFORMATION FROM SETUP OBJECT    ==
!     ==================================================================
      CALL SETUP$AEZ(ISP,AEZ)
      CALL SETUP$RADGRID(ISP,R1,DEX,NR)
      XEXP=DEXP(DEX)
      CALL SETUP$LNX(ISP,LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$LOFLN(ISP,LNX,LOX)

! BENGONE MODIFY START 
! 08-04-2002 
! IT SEEMS THAT THE INDICES ARE REVERSED 
! EVERY SUBROUTINES BELOW USE AEPHI(NR,LNX) 
!-> SETUP$AEPARTIALWAVES
!-> AUGMENTATION_RHO
!-> AUGMENTATION_EXPECT 
! DON'T KNOW HOW IT CAN WORK !!!!!!
! SAME FOR PSPHI 

      ALLOCATE(AEPHI(NR,LNX))
!     ALLOCATE(AEPHI(LNX,NR))
! BENGONE MODIFY END  

      CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)

! BENGONE MODIFY START 
! 08-04-2002 
      ALLOCATE(PSPHI(NR,LNX))
!     ALLOCATE(PSPHI(LNX,NR))
! BENGONE MODIFY END  

      CALL SETUP$PSPARTIALWAVES(ISP,NR,LNX,PSPHI)
      ALLOCATE(AECORE(NR))
      CALL SETUP$AECORE(ISP,NR,AECORE)
      ALLOCATE(PSCORE(NR))
      CALL SETUP$PSCORE(ISP,NR,PSCORE)
      CALL SETUP$RCSM(ISP,RCSM)
      ALLOCATE(VADD(NR))
      CALL SETUP$VBAR(ISP,nr,VADD)
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
            WRITE(NFILO,FMT='(9F10.6)')DENMAT(LMN1,:,IDIM)
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
!     ==  ADD CORE CHARGE DENSITY                                   ==
!     ================================================================
      CALL AUGMENTATION_QLM(R1,DEX,NR,NR,LMRX &
     &                   ,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
!     
!     ================================================================
!     ==  SEND DENSITIES TO HYPERFINE-PARAMETER OBJECT              ==
!     ================================================================
      ALLOCATE(AEPOT(NR,LMRX,1))   ! USED AS AUXILIARY VARIABLE
      CALL HYPERFINE$SET1CRHO('PS','TOT',IAT,R1,DEX,NR,NR,LMRX,PSRHO)
      AEPOT(:,:,1)=AERHO(:,:,1)
      AEPOT(:,1,1)=AEPOT(:,1,1)+AECORE(:)
      CALL HYPERFINE$SET1CRHO('AE','TOT',IAT,R1,DEX,NR,NR,LMRX,AEPOT)
      IF(NDIMD.GT.1) THEN
        ALLOCATE(PSPOT(NR,LMRX,1))
        IF(NDIMD.EQ.2) THEN        ! COLLINEAR SPIN POLARIZED 
          AEPOT(:,:,1)=AERHO(:,:,2)
          PSPOT(:,:,1)=PSRHO(:,:,2)
        ELSE IF(NDIMD.EQ.4) THEN   ! NON-COLLINEAR SPIN POLARIZED 
          AEPOT(:,:,1)=SQRT(AERHO(:,:,2)**2+AERHO(:,:,3)**2+AERHO(:,:,4)**2)
          PSPOT(:,:,1)=SQRT(PSRHO(:,:,2)**2+PSRHO(:,:,3)**2+PSRHO(:,:,4)**2)
        END IF
        CALL HYPERFINE$SET1CRHO('PS','SPIN' &
     &                         ,IAT,R1,DEX,NR,NR,LMRX,PSPOT)
        CALL HYPERFINE$SET1CRHO('AE','SPIN' &
     &                         ,IAT,R1,DEX,NR,NR,LMRX,AEPOT)
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
      CALL AUGMENTATION_XC(R1,DEX,NR,LMRX,NDIMD,AERHO,AEEXC,AEPOT)
      AERHO(:,1,1)=AERHO(:,1,1)-AECORE(:)
!     == PS-EXCHANGE ENERGY AND POTENTIAL ============================
      CALL AUGMENTATION_XC(R1,DEX,NR,LMRX,NDIMD,PSRHO,PSEXC,PSPOT)
!     == CORE ONLY EXCHANGE ENERGY
      COREEXC=0.D0
      ALLOCATE(DWORK1(NR))
      CALL AUGMENTATION_XC(R1,DEX,NR,1,1,AECORE,COREEXC,DWORK1)
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
      CALL AUGMENTATION_HARTREE(IAT,R1,DEX,NR,AEZ,AECORE &
     &          ,VADD,RCSM,LMRX,AERHO,PSRHO,QLM,AEPOT1,PSPOT1,VQLM1 &
     &          ,AEEHARTREE,PSEHARTREE)
      CALL AUGMENTATION_ADD('AE1 ELECTROSTATIC',AEEHARTREE)
      CALL AUGMENTATION_ADD('PS1 ELECTROSTATIC',PSEHARTREE)
      VQLM1(:)=VQLM1(:)+VQLM(:)
!     
!     == ADD POTENTIAL FROM CHARGES "OUTSIDE" THE SPHERE =============      
      CALL AUGMENTATION_ADDVQLM(R1,DEX,NR,LMRX,VQLM1,AEPOT1,PSPOT1)
!     
!     == ADD POTENTIAL FROM THE NEUTRALIZING BACKGROUND ==============
!     == THE BACKGROUND DOES NOT ADD TO SUM OF ONE-CENTER ENERGIES  ==
!     == AND NEED NOT BE CONSIDERED. IT SHOULD BE INCLUDED, IF THE  ==
!     == POTENTIALS SHOULD BE THE CORRECT ONE-CENTER POTENTIALS     ==
!     == AT THIS POIINT THE DERIVATIVES ARE NOT CORRECT YET         ==
      AEBACKGROUND=0.D0
      PSBACKGROUND=0.D0
      CALL AUGMENTATION_ADDBACKGROUND(R1,DEX,NR,RHOB &
     &                      ,AERHO,AEBACKGROUND,AEPOT1)
      CALL AUGMENTATION_ADDBACKGROUND(R1,DEX,NR,RHOB &
     &                      ,PSRHO,PSBACKGROUND,PSPOT1)
      CALL AUGMENTATION_ADD('AE1 BACKGROUND',AEBACKGROUND)
      CALL AUGMENTATION_ADD('PS1 BACKGROUND',PSBACKGROUND)
!     
!     == ANALYSIS: ELECTRIC FIELD GRADIENTS ==========================
      CALL HYPERFINE$SET1CPOT('AE',IAT,R1,DEX,NR,NR,LMRX,AEPOT1)
      CALL HYPERFINE$SET1CPOT('PS',IAT,R1,DEX,NR,NR,LMRX,PSPOT1)
!     
!     == ANALYSIS: POTENTIAL PLOT           ==========================
      CALL GRAPHICS$SET1CPOT('AE',IAT,R1,DEX,NR,NR,LMRX,AEPOT1)
      CALL GRAPHICS$SET1CPOT('PS',IAT,R1,DEX,NR,NR,LMRX,PSPOT1)
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
          ENL=0.D0
          DO LM=1,LMRX
            RI=R1/XEXP
            DO IR=1,NR
              RI=RI*XEXP
              ENL=ENL+AERHO(IR,LM,1)**2*DEX*RI**3
              AEPOT(IR,LM,1)=2.D0*AERHO(IR,LM,1)
              PSPOT(IR,LM,1)=0.D0
            ENDDO
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
      CALL AUGMENTATION_EXPECT(R1,DEX,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                        ,AEPOT,PSPOT,AEPHI,PSPHI,DATP)
      DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
      DEALLOCATE(DATP)
!     WRITE(TESTSTRING,FMT='("R8DATH",I2,12(" "))')IAT
!     CALL STOREIT(TESTSTRING,8*LMNX*LMNX*NSPIN,DATH)
!
!     ================================================================
!     ==  LDA + U                                                   ==
!     ================================================================

! BENGONE MODIFY START 
! 08-04-2002 
!     WRITE(*,*) 
!     WRITE(*,*) 'JUST BEFORE LDA+U'    
!     WRITE(*,*) 'DIMENSION DE AEPHI IAT=',IAT  
!     WRITE(*,*) 'DIM 1 = ', SIZE(AEPHI,1)
!     WRITE(*,*) 'DIM 2 = ', SIZE(AEPHI,2)
!     WRITE(*,*) 


      DETOT=0.D0
      CALL LDAPLUSU(IAT,NR,LNX,LMNX,NDIMD,LOX  &
      &             ,R1,DEX,AEZ,AEPHI,DENMAT,ELDA_U,DATH)
!     CALL AUGMENTATION_ADD('LDA+U EXCHANGE',DETOT)

! BENGONE MODIFY END 

 
!     ================================================================
!     ==  APPLY EXTERNAL POTENTIAL                                  ==
!     ================================================================
      CALL ATOMLIST$GETCH('NAME',IAT,ATOM)
      ALLOCATE(DATP(LMNX,LMNX,Ndimd))
      CALL EXTERNAL1CPOT$APPLY(ATOM,LMNX,NDIMD,DENMAT,DATP,DETOT)
      DATH(:,:,:)=DATH(:,:,:)+DATP(:,:,:)
      DEALLOCATE(DATP)

! BENGONE MODIFY START 
! 08-04-2002 
      DETOT = ELDA_U + DETOT 
! BENGONE MODIFY END 
      CALL AUGMENTATION_ADD('LDA+U EXCHANGE',DETOT)
!     
!     ================================================================
!     ==  SELF-TEST                                                 ==
!     ================================================================
      IF(TTEST) THEN
        IF(ITEST.EQ.1.OR.ITEST.EQ.3) THEN
          ENL=AEEXC-PSEXC+AEEHARTREE-PSEHARTREE+EKINNL+DETOT
          CALL SELFTEST$END('SPHERE',LMNX*LMNX*NDIMD,dath,enl,TBACK)
          IF(TBACK) GOTO 1000
          CALL ERROR$MSG('STOP AFTER SELFTEST')
          CALL ERROR$STOP('SPHERE')
        ENDIF
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
                               CALL TRACE$POP
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
      INTEGER(4)             :: LMR,IR
      INTEGER(4)             :: LMN1,LN1,L1,IM1,LMN2,LN2,L2,IM2
      INTEGER(4)             :: LM1,LM2,LM3
      REAL(8)                :: CG     !CLEBSCH GORDAN COEFFICIENT
      REAL(8)                :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     ==   CALCULATE ONE CENTER CHARGE DENSITY                        ==
!     ==================================================================
      DO LMR=1,LMRX
        DO IR=1,NR
          RHOL(IR,LMR)=0.D0
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ADD VALENCE CHARGE DENSITY                                  ==
!     ==================================================================
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
                  DO IR=1,NR
                    RHOL(IR,LM3)=RHOL(IR,LM3)+SVAR*PHI(IR,LN1)*PHI(IR,LN2)
                  ENDDO
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
      SUBROUTINE AUGMENTATION_QLM(R1,DEX,NR,NRX,LMRX &
     &                           ,AEZ,AECORE,PSCORE,AERHO,PSRHO,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE MULTIPOLE MOMENTS OF THE DEVIATION OF             **
!     **  ALL-ELECTRON AND THE PSEUDO DENSITY FROM ONE CENTER         **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NRX
      INTEGER(4),INTENT(IN) :: LMRX
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AECORE(NRX)
      REAL(8)   ,INTENT(IN) :: PSCORE(NRX)
      REAL(8)   ,INTENT(IN) :: AERHO(NRX,LMRX)
      REAL(8)   ,INTENT(IN) :: PSRHO(NRX,LMRX)
      REAL(8)   ,INTENT(OUT):: QLM(LMRX)
      REAL(8)               :: DWORK(NR)
      REAL(8)               :: PI,Y0
      REAL(8)               :: XEXP,RI
      REAL(8)               :: RES
      INTEGER(4)            :: LM,L,IR
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      XEXP=DEXP(DEX)
      QLM(:)=0.D0
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        RI=R1/XEXP
        DO IR=1,NR
          RI=RI*XEXP
          DWORK(IR)=(AERHO(IR,LM)-PSRHO(IR,LM))*RI**(L+2)
        ENDDO
        CALL RADIAL$INTEGRAL(R1,DEX,NR,DWORK,QLM(LM))
      ENDDO
!
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        DWORK(IR)=(AECORE(IR)-PSCORE(IR))*RI**2
      ENDDO
      CALL RADIAL$INTEGRAL(R1,DEX,NR,DWORK,RES)
      QLM(1)=QLM(1)+RES-AEZ*Y0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_XC(R1,DEX,NR,LMRX,NDIMD,RHOIN,EXC,VXC)
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
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
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
      REAL(8)               :: XEXP
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
      real(8)   ,parameter  :: tiny=1.d-300
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
      XEXP=DEXP(DEX)
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        R(IR)=RI
      ENDDO
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
        DO IR=1,NR
          WORK(IR)=0.5D0/(B(IR,1)*Y0+tiny)
        ENDDO
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
            CALL RADIAL$DERIVE(R1,DEX,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
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
      CALL RADIAL$INTEGRAL(R1,DEX,NR,WORK1(:)*R(:)**2,EXC)
!
!     ==================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                     ==
!     ==================================================================
      allocate(vgrho(nr,lmrx,nspin))
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
!           CALL RADIAL$DERIVE(R1,DEX,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(R1,DEX,NR,WORK2,WORK1)
            WORK1(:)=WORK1(:)/R(:)**2
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      deallocate(vgrho)
!
!     ==================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS             ==
!     ==================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=Vrho(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        allocate(vb(nr,lmrx))
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
        WORK(:)=0.5D0/SQRT(B(:,1)*Y0+tiny)
        DO LM1=1,LMRX
          VB(:,LM1)=VB(:,LM1)*WORK(:)
        ENDDO
        VB(:,1)=VB(:,1)*2.d0*Y0
        WORK(:)=0.D0
        DO LM1=1,LMRX
          WORK(:)=WORK(:)+VRHO(:,LM1,2)*RHO(:,LM1,2)
        ENDDO
        VB(:,1)=VB(:,1)+0.5D0*WORK(:)/(B(:,1)+tiny)
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
      SUBROUTINE AUGMENTATION_HARTREE(IAT,R1,DEX,NR,AEZ,RHOCOR &
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
      REAL(8)    ,INTENT(IN) :: R1
      REAL(8)    ,INTENT(IN) :: DEX
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
      REAL(8)                :: AUX1(NR)
      REAL(8)                :: AEDMU(NR)
      REAL(8)                :: PSDMU(NR)
      REAL(8)                :: AEE(NR)
      REAL(8)                :: PSE(NR)
      REAL(8)                :: RI,RI2,RIL
      REAL(8)                :: PI
      REAL(8)                :: Y0
      REAL(8)                :: XEXP
      REAL(8)                :: ALPHA   !1/RCSM
      REAL(8)                :: CL
      REAL(8)                :: SVAR
      REAL(8)                :: GAUSSIAN
      REAL(8)                :: AEPOT1,PSPOT1,AERHO1,PSRHO1,QLM1
      REAL(8)                :: AEH,PSH
      INTEGER(4)             :: LM,IR
      INTEGER(4)             :: L
!     == ARRAYS NEEDED DETAILED REPORT
      LOGICAL(4) ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4)             :: NFILO
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      XEXP=DEXP(DEX)
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
      CALL RADIAL$POISSON(R1,DEX,NR,0,RHOCOR,AEDMU)
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        RI2=RI**2
        AEDMU(IR)  = AEDMU(IR)-AEZ/RI/Y0
        AEE(IR)    = AEDMU(IR)* AERHO(IR,1) *RI2
        AEPOT(IR,1)= AEPOT(IR,1)+AEDMU(IR)
        PSE(IR)    = VADD(IR) * PSRHO(IR,1) *RI2
        PSPOT(IR,1)= PSPOT(IR,1)+VADD(IR)
      ENDDO 
      CALL RADIAL$INTEGRAL(R1,DEX,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(R1,DEX,NR,PSE,PSH)
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
        RI=R1/XEXP
        DO IR=1,NR
          RI=RI*XEXP
          RI2=RI**2
          RIL=RI**L
          PSRHO(IR,LM)=PSRHO(IR,LM)+SVAR*DEXP(-ALPHA*RI2)*RIL
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==   CALCULATE ELECTROSTATIC POTENTIAL AND ENERGY               ==
!     ==================================================================
      AEE(:)=0.D0
      PSE(:)=0.D0
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        CALL RADIAL$POISSON(R1,DEX,NR,L,AERHO(1,LM),AEDMU)
        CALL RADIAL$POISSON(R1,DEX,NR,L,PSRHO(1,LM),PSDMU)
        ALPHA=1.D0/RCSM**2
        CALL GAUSSN(L,ALPHA,CL)
        QLM1=QLM(LM)
        RI=R1/XEXP
        DO IR=1,NR
          RI=RI*XEXP
          RI2=RI**2
          AEPOT1=AEDMU(IR)
          PSPOT1=PSDMU(IR)
          AERHO1=AERHO(IR,LM)
          PSRHO1=PSRHO(IR,LM)
!
!         == ADD COMPENSATION DENSITY AND ITS POTENTIAL ================
          GAUSSIAN=CL*DEXP(-ALPHA*RI2)*RI**L
!
!         == CALCULATE TOTAL ENERGY AND POTENTIAL ======================
          AEE(IR)=AEE(IR)+0.5D0*AEPOT1*AERHO1*RI2
          PSE(IR)=PSE(IR)+0.5D0*PSPOT1*PSRHO1*RI2
          AEPOT(IR,LM)=AEPOT(IR,LM)+AEPOT1
          PSPOT(IR,LM)=PSPOT(IR,LM)+PSPOT1
          AUX1(IR)=PSPOT1*GAUSSIAN*RI2
        ENDDO
        CALL RADIAL$INTEGRAL(R1,DEX,NR,AUX1,SVAR)
        VQLM(LM)=VQLM(LM)-SVAR
      ENDDO
      CALL RADIAL$INTEGRAL(R1,DEX,NR,AEE,AEH)
      CALL RADIAL$INTEGRAL(R1,DEX,NR,PSE,PSH)
      AEEHARTREE=AEEHARTREE+AEH
      PSEHARTREE=PSEHARTREE+PSH
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
      SUBROUTINE AUGMENTATION_ADDVQLM(R1,DEX,NR,LMRX,VQLM,AEPOT,PSPOT)
!     ******************************************************************
!     **                                                             **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: R1
      REAL(8)   ,INTENT(IN)   :: DEX
      INTEGER(4),INTENT(IN)   :: NR
      INTEGER(4),INTENT(IN)   :: LMRX
      REAL(8)   ,INTENT(IN)   :: VQLM(LMRX)
      REAL(8)   ,INTENT(INOUT):: AEPOT(NR,LMRX)
      REAL(8)   ,INTENT(INOUT):: PSPOT(NR,LMRX)
      REAL(8)                 :: XEXP,RI,RIL
      REAL(8)                 :: SVAR
      INTEGER(4)              :: LM,L,IR
!     ******************************************************************
      XEXP=DEXP(DEX)
      DO LM=1,LMRX
        L=INT(DSQRT(DBLE(LM-1)+.01D0))
        SVAR=VQLM(LM)
        RI=R1/XEXP
        DO IR=1,NR
          RI=RI*XEXP
          RIL=RI**L
          AEPOT(IR,LM)=AEPOT(IR,LM)+SVAR*RIL
          PSPOT(IR,LM)=PSPOT(IR,LM)+SVAR*RIL
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_ADDBACKGROUND(R1,DEX,NR,RHOB &
     &                        ,RHO,EB,POT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ONE-CENTER CONTRIBUTION OF THE COMPENSATING  **
!     **  CHARGE BACKGROUND TO TOTAL ENERGY AND POTENTIAL             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: R1
      REAL(8)   ,INTENT(IN)   :: DEX
      INTEGER(4),INTENT(IN)   :: NR
      REAL(8)   ,INTENT(IN)   :: RHOB
      REAL(8)   ,INTENT(IN)   :: RHO(NR)
      REAL(8)   ,INTENT(OUT)  :: EB
      REAL(8)   ,INTENT(INOUT):: POT(NR)
      REAL(8)                 :: PI
      REAL(8)                 :: XEXP
      REAL(8)                 :: SVAR,RI,RI2,SVAR1
      REAL(8)                 :: WORK(NR)
      INTEGER(4)              :: IR
!     ******************************************************************
      IF(RHOB.EQ.0.D0) THEN
        EB=0.D0
        RETURN
      END IF
      PI=4.D0*DATAN(1.D0)
      XEXP=DEXP(DEX)
      SVAR=-2.D0*PI*RHOB/3.D0*DSQRT(4.D0*PI)
      RI=R1/XEXP
print*,'severe warning from augmentation_addbackground'
      DO IR=1,NR
        RI=RI*XEXP
        RI2=RI*RI
        SVAR1=SVAR*RI2
WORK(IR)=SVAR1*RHO(IR)
!       WORK(IR)=SVAR1*RHO(IR)*ri2
        POT(IR)=POT(IR)+SVAR1
      ENDDO
      CALL RADIAL$INTEGRAL(R1,DEX,NR,WORK,EB)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE AUGMENTATION_EXPECT(R1,DEX,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &           ,AEPOT,PSPOT,AEPHI,PSPHI,DATH)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE EXPECTATION VALUE OF                         **
!     **  THE ONE-CENTER POTENTIALS WITH THE PARTIAL WAVES            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
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
      INTEGER(4)            :: IR
      REAL(8)               :: AEDMU(NR,NDIMD)
      REAL(8)               :: PSDMU(NR,NDIMD)
      REAL(8)               :: DWORK1(NR)
      REAL(8)               :: CG
      REAL(8)               :: SVAR
      REAL(8)               :: XEXP,RI
!     ******************************************************************
      XEXP=DEXP(DEX)
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
                    DO IR=1,NR
                      AEDMU(IR,ISPIN)=AEDMU(IR,ISPIN) &
     &                               +CG*AEPOT(IR,LM3,ISPIN)
                      PSDMU(IR,ISPIN)=PSDMU(IR,ISPIN) &
     &                               +CG*PSPOT(IR,LM3,ISPIN)
                    ENDDO
                  ENDDO
                END IF
              ENDDO
!     
!             ==========================================================
!             ==  PERFORM NOW THE INTEGRATION                         ==
!             ==========================================================
              DO ISPIN=1,NDIMD
                RI=R1/XEXP
                DO IR=1,NR
                  RI=RI*XEXP
                  DWORK1(IR)= &
     &               (AEDMU(IR,ISPIN)*AEPHI(IR,LN1)*AEPHI(IR,LN2) &
     &               -PSDMU(IR,ISPIN)*PSPHI(IR,LN1)*PSPHI(IR,LN2))*RI**2
                ENDDO
                CALL RADIAL$INTEGRAL(R1,DEX,NR,DWORK1,SVAR)
                DATH(LMN1,LMN2,ISPIN)=DATH(LMN1,LMN2,ISPIN)+SVAR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      MODULE EXPERTNAL1CPOT_MODULE
      TYPE EXTPOT
        REAL(8)          :: VALUE
        CHARACTER(LEN=32):: ATOM
        REAL(8)          :: RC  
        REAL(8)          :: PWR
        CHARACTER(LEN=32):: TYPE   ! CAN BE 'S','P','D','ALL'
        INTEGER(4)       :: idimd  ! idimd=0 REFERES TO ALL SPIN DIRECTIONS
      END TYPE EXTPOT
      INTEGER(4)             :: NPOT=0
      INTEGER(4)             :: NPOTX=0
      INTEGER(4),PARAMETER   :: DNPOT=5
      TYPE (EXTPOT), ALLOCATABLE :: POT(:)
      CONTAINS
!       ...............................................................
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
      INTEGER(4)       ,INTENT(IN) :: iDIMD
      REAL(8)          ,INTENT(IN) :: RC
      REAL(8)          ,INTENT(IN) :: PWR
!     ******************************************************************
      CALL CREATE
      NPOT=NPOT+1
      POT(NPOT)%ATOM =ATOM
      POT(NPOT)%VALUE=VALUE
      POT(NPOT)%RC=RC
      POT(NPOT)%PWR=PWR

! BENGONE MODIFY START 
! 08-04-2002 
! IF SPIN= 1 or 2 MUST BE SAVED 
! 1 OR 2 OTHERWISE PROBLEM IN 
!-> EXTERNAL1CPOT$APPLY 
      POT(NPOT)%idimd=idimd
!     POT(NPOT)%idimd=idimd+1
! BENGONE MODIFY START 

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
        CALL WRITEI4(NFIL,'IDIMD ([NT],[NT,NS],[NT,NX,NY,NZ])',POT(IPOT)%IDIMD,' ')
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
      SUBROUTINE EXTERNAL1CPOT$APPLY(ATOM,LMNX,ndimd,DENMAT,DATH,ETOT)
!     ******************************************************************
!     ******************************************************************
      USE EXPERTNAL1CPOT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: ATOM
      INTEGER(4)  ,INTENT(IN)   :: Ndimd
      INTEGER(4)  ,INTENT(IN)   :: LMNX
      REAL(8)     ,INTENT(IN)   :: DENMAT(LMNX,LMNX,Ndimd)
      REAL(8)     ,INTENT(OUT)  :: DATH(LMNX,LMNX,Ndimd)
      REAL(8)     ,INTENT(OUT)  :: ETOT
      TYPE(EXTPOT)              :: POT1
      INTEGER(4)                :: LANG
      INTEGER(4)                :: IAT    ! ATOM INDEX
      INTEGER(4)                :: ISP    ! ATOM TYPE INDEX
      REAL(8)                   :: R1     ! FIRST POINT ON RADIAL GRID
      REAL(8)                   :: DEX
      REAL(8)                   :: XEXP   ! EXP(DEX)
      REAL(8)                   :: RI     ! RADIAL GRID POINT
      INTEGER(4)                :: NR     ! #(GRID RADIAL POINTS)
      INTEGER(4)                :: LNX    ! #(PARTIAL WAVES (L,N))
      INTEGER(4)   ,ALLOCATABLE :: LOX(:) !(LNX) #(GRID RADIAL POINTS)
      REAL(8)      ,ALLOCATABLE :: AEPHI(:,:) !(NR,LNX) AE PARTIAL WAVES
      REAL(8)      ,ALLOCATABLE :: UONE(:,:)  !(LNX,LNX)
      REAL(8)      ,ALLOCATABLE :: RDEP(:)    !(NR)
      REAL(8)      ,ALLOCATABLE :: AUX(:)    !(NR)
      CHARACTER(32)             :: SPECIES
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: LN1,LN2,LMN1,LMN2,IR,IPOT,Idimd,L1,L2,I
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
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETCH('ID',SPECIES)
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
        CALL SETUP$RADGRID(ISP,R1,DEX,NR)
        XEXP=EXP(DEX)
        CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$LOFLN(ISP,LNX,LOX)
        ALLOCATE(AEPHI(NR,LNX))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
!
        ALLOCATE(RDEP(NR))
        ALLOCATE(AUX(NR))
        RI=R1/XEXP
        DO IR=1,NR
          RI=RI*XEXP
          RDEP(IR)=POT1%VALUE*EXP(-(RI/POT1%RC)**POT1%PWR)
        ENDDO
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
            RI=R1/XEXP
            DO IR=1,NR
              RI=RI*XEXP
              AUX(IR)=RI**2*RDEP(IR)*AEPHI(IR,LN1)*AEPHI(IR,LN2)
            ENDDO
            CALL RADIAL$INTEGRAL(R1,DEX,NR,AUX,UONE(LN1,LN2))
          ENDDO
        ENDDO
        DEALLOCATE(AEPHI)
        DEALLOCATE(RDEP)
        DEALLOCATE(AUX)
PRINT*,'LNX ',LNX,' LOX ',LOX
PRINT*,'UONE ',UONE
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
      DO Idimd=1,ndimd
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            ETOT=ETOT+DENMAT(LMN1,LMN2,idimd)*DATH(LMN1,LMN2,Idimd)
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
      real(8)               :: erfx
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

