!
!     .................................................................
      MODULE STATEANALYSIS_MODULE
      INTEGER(4)             :: NAT=0
      INTEGER(4)             :: NB=0
      INTEGER(4)             :: LMNXX=0
      INTEGER(4)             :: NKPT=0
      INTEGER(4)             :: NSPIN=0
      REAL(8)   ,ALLOCATABLE :: FNL(:,:,:,:,:) !(NAT,NB,LMNXX,NKPT,NSPIN)
      END MODULE STATEANALYSIS_MODULE
!
!     .................................................................. 
      SUBROUTINE STATEANALYSIS
!     ******************************************************************
!     ******************************************************************
!     ******************************************************************
      USE STATEANALYSIS_MODULE
      IMPLICIT NONE
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:,:)    !(NR,LNXX,NSP)
      INTEGER(4),ALLOCATABLE :: LOX(:,:)        !(LNXX,NSP)
      INTEGER(4),ALLOCATABLE :: LNX(:)          !(NSP)
      INTEGER(4),ALLOCATABLE :: LMNX(:)         !(NSP)
      REAL(8)   ,ALLOCATABLE :: AEZ(:)          !(NSP)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:)     !(NAT)
      REAL(8)   ,ALLOCATABLE :: ZR(:,:,:,:)     !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: HAMLTN(:,:,:,:) !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: EIG(:,:,:)      !(NB,NKPT,NSPIN)
      INTEGER(4)             :: NB1,NSPIN1,NKPT1,LMNXX1
      INTEGER(4)             :: IB,IKPT,ISPIN,ISP
      INTEGER(4)             :: NSP,LNXX,NR,NRX
      REAL(8)                :: R1,DEX
!     ******************************************************************
                              CALL TRACE$PUSH('STATEANALYSIS')
call error$msg('routine marked for deletion')
call error$stop('STATEANALYSIS')
!
      IF(.NOT.ALLOCATED(FNL)) THEN
        CALL ERROR$MSG('FNL NOT SET')
        CALL ERROR$STOP('STATEANALYSIS')
      END IF
      CALL DYNOCC$GETI4('NB',NB1)
      CALL DYNOCC$GETI4('NKPT',NKPT1) 
      CALL DYNOCC$GETI4('NSPIN',NSPIN1)
      IF(NB.NE.NB1.OR.NKPT.NE.NKPT1.OR.NSPIN.NE.NSPIN1) THEN
        CALL ERROR$I4VAL('NB    ',NB)
        CALL ERROR$I4VAL('NB1   ',NB1)
        CALL ERROR$I4VAL('NKPT  ',NKPT)
        CALL ERROR$I4VAL('NKPT1 ',NKPT1)
        CALL ERROR$I4VAL('NSPIN ',NSPIN)
        CALL ERROR$I4VAL('NSPIN1',NSPIN1)
        CALL ERROR$MSG('DIMENSIONS INCONSISTENT')
        CALL ERROR$STOP('STATEANALYSIS')
      END IF
!
!     ==================================================================
!     ==  DIAGONALIZE HAMILTONIAN                                     ==
!     ==================================================================
                              CALL TRACE$PASS('DIAGONALIZE HAMILTONIAN')
      ALLOCATE(HAMLTN(NB,NB,NKPT,NSPIN))
      ALLOCATE(EIG(NB,NKPT,NSPIN))
      ALLOCATE(ZR(NB,NB,NKPT,NSPIN))
      DO ISPIN=1,NSPIN
        CALL WAVES$SETI4('ISPIN',ISPIN)
        DO IKPT=1,NKPT
          CALL WAVES$SETI4('IKPT',IKPT)
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB,HAMLTN(1,1,IKPT,ISPIN))
          CALL LIB$DIAGR8(NB,HAMLTN(1,1,IKPT,ISPIN) &
     &             ,EIG(1,IKPT,ISPIN),ZR(1,1,IKPT,ISPIN))
!CALL DIAG(NB,NB,HAMLTN(1,1,IKPT,ISPIN) &
!&,EIG(1,IKPT,ISPIN),ZR(1,1,IKPT,ISPIN))
!         ZR(:,:,IKPT,ISPIN)=0.D0
!         DO IB=1,NB
!           ZR(IB,IB,IKPT,ISPIN)=1.D0
!         ENDDO
        ENDDO
      ENDDO
      CALL OCCUPATION$SET('EIG',NB*NKPT*NSPIN,EIG)
      DEALLOCATE(HAMLTN)
!
!     ==================================================================
!     ==  CALCULATE OVERLAP MATRIX OF AE PARTIAL WAVES                ==
!     ==================================================================
                              CALL TRACE$PASS('CALCULATE OVERLAP')
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(AEZ(NSP))
      ALLOCATE(LNX(NSP))
      ALLOCATE(LMNX(NSP))
      LNXX=0
      LMNXX1=0
      DO ISP=1,NSP
        CALL SETUP$RADGRID(ISP,R1,DEX,NR)
        CALL SETUP$AEZ(ISP,AEZ(ISP))
        CALL SETUP$LNX(ISP,LNX(ISP))        
        CALL SETUP$LMNX(ISP,LMNX(ISP))        
        LMNXX1=MAX(LMNXX1,LMNX(ISP))
        LNXX=MAX(LNXX,LNX(ISP))
      ENDDO
      IF(LMNXX1.NE.LMNXX) THEN
        CALL ERROR$MSG('INCONSISTENCY OF LMNXX')
        CALL ERROR$STOP('STATEANALYSIS')
      END IF
      ALLOCATE(AEPHI(NR,LNXX,NSP))
      ALLOCATE(LOX(LNXX,NSP))
      DO ISP=1,NSP
        CALL SETUP$LOFLN(ISP,LNX(ISP),LOX(1,ISP))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX(ISP),AEPHI(1,1,ISP))
      ENDDO
      NRX=NR
                              CALL TRACE$PASS('CALCULATE OVERLAPG')
      CALL PDOS(NRX,LNXX,LMNXX,NB,NKPT,NB,NKPT,NSPIN &
     &         ,ISPECIES &
     &         ,NSP,AEZ,R1,DEX,NR,LNX,LOX,AEPHI &
     &         ,ZR,NAT,NB,FNL,EIG)
                              CALL TRACE$PASS('CALCULATE OVERLAPH')
      DEALLOCATE(LOX)
      DEALLOCATE(AEZ)
      DEALLOCATE(LNX)
      DEALLOCATE(LMNX)
      DEALLOCATE(AEPHI)
      DEALLOCATE(ISPECIES)
      DEALLOCATE(ZR)
      DEALLOCATE(FNL)
      DEALLOCATE(EIG)
                              CALL TRACE$POP
      RETURN
      END
!
!     .................................................................
      SUBROUTINE STATEANALYSIS$SETPROJECTIONS(NAT_,NB_,LMNXX_,NKPT_,NSPIN_ &
     &                   ,FNL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE STATEANALYSIS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_
      INTEGER(4),INTENT(IN) :: NB_
      INTEGER(4),INTENT(IN) :: LMNXX_
      INTEGER(4),INTENT(IN) :: NKPT_
      INTEGER(4),INTENT(IN) :: NSPIN_
      REAL(8)   ,INTENT(IN) :: FNL_(NAT_,NB_,LMNXX_,NKPT_,NSPIN_)
      INTEGER(4)            :: ISPIN,IKPT,LMN,IB,IAT
!     ******************************************************************
call error$msg('routine marked for deletion')
call error$stop('STATEANALYSIS$SETPROJECTIONS')
      NAT=NAT_
      NB=NB_
      LMNXX=LMNXX_
      NKPT=NKPT_
      NSPIN=NSPIN_
      IF(.NOT.ALLOCATED(FNL)) THEN
        ALLOCATE(FNL(NAT,NB,LMNXX,NKPT,NSPIN))
      END IF
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO LMN=1,LMNXX
            DO IB=1,NB
              DO IAT=1,NAT
                FNL(IAT,IB,LMN,IKPT,ISPIN)=FNL_(IAT,IB,LMN,IKPT,ISPIN)
              ENDDO 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS(NRX,LNXX,LMNXX,NX,NKPTX &
     &               ,NB,NKPT,NSPIN,ISPECIES &
     &               ,NSP,AEZ,R1,DEX,NR,LNX,LOX,AEPHI &
     &               ,ZR,NAT,NBFNLX,FNL,EIG)
!     ******************************************************************
!     **                                                              **
!     **  PDOS CALCULATES                                             **
!     **  1) THE WAVE FUNCTION AMPLITUDE FOR ENERGY DEPENDENT         **
!     **     PARTIAL WAVES WHICH ARE NORMALIZED WITHIN AN ATOMIC      **
!     **     SPHERE, AND WHOSE SIGN IS DEFINED SUCH THAT THE          **
!     **     PARTIAL WAVES FOR THE VALENCE STATES ARE POSITIVE        **
!     **     IN THE BINDING REGION.                                   **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(32)          :: STRING
      LOGICAL(4)             :: TPDOS=.FALSE.
      LOGICAL(4)             :: TAMP=.TRUE.
      LOGICAL(4)             :: TWEIGHT=.TRUE.
      LOGICAL(4)             :: TAMPFUL=.TRUE.
      LOGICAL(4)             :: TPDOS1
      LOGICAL(4)             :: TAMP1
      LOGICAL(4)             :: TWEIGHT1
      LOGICAL(4)             :: TAMPFUL1
      INTEGER(4) ,INTENT(IN) :: NRX
      INTEGER(4) ,INTENT(IN) :: LNXX
      INTEGER(4) ,INTENT(IN) :: LMNXX
      INTEGER(4) ,INTENT(IN) :: NX
      INTEGER(4) ,INTENT(IN) :: NKPTX
      INTEGER(4) ,INTENT(IN) :: NB
      INTEGER(4) ,INTENT(IN) :: NAT
      INTEGER(4) ,INTENT(IN) :: NKPT
      INTEGER(4) ,INTENT(IN) :: NSPIN
      INTEGER(4) ,INTENT(IN) :: ISPECIES(NAT)
      INTEGER(4) ,INTENT(IN) :: NSP
      REAL(8)    ,INTENT(IN) :: AEZ(NSP)
      REAL(8)    ,INTENT(IN) :: R1
      REAL(8)    ,INTENT(IN) :: DEX
      INTEGER(4) ,INTENT(IN) :: NR
      INTEGER(4) ,INTENT(IN) :: LNX(NSP)
      INTEGER(4) ,INTENT(IN) :: LOX(LNXX,NSP)
      REAL(8)    ,INTENT(IN) :: AEPHI(NRX,LNXX,NSP)
      REAL(8)    ,INTENT(IN) :: ZR(NX,NX,NKPT,NSPIN)
      INTEGER(4) ,INTENT(IN) :: NBFNLX
      REAL(8)    ,INTENT(IN) :: FNL(NAT,NBFNLX,LMNXX,NKPT,NSPIN)
      REAL(8)    ,INTENT(IN) :: EIG(NX,NKPTX,NSPIN)
      REAL(8)    ,ALLOCATABLE:: WOFL(:,:,:,:,:) !(LMXPDOS,NAT,NB,NKPT,NSPIN)
      REAL(8)    ,ALLOCATABLE:: WEIGHT(:)       !(LMXPDOS)
      REAL(8)    ,ALLOCATABLE:: WORK(:)         !(LMXPDOS*NAT)
      INTEGER(4) ,ALLOCATABLE:: IWORK(:)        !(LMXPDOS*NAT)
      REAL(8)                :: RPOS(3) 
      INTEGER(4)             :: I,II,IE,IAT,ISP,LN,ISPIN,LM,LMX,IB,IKPT
      INTEGER(4)             :: NORB,LMAX,LXPDOS,LMXPDOS,NFIL,LISTPERCENT
      INTEGER(4)             :: ISVAR
      INTEGER(4)             :: NE
      REAL(8)                :: EV
      REAL(8)                :: SUM,WIDTH,E,EMIN,EMAX,TOTWEIGHT
      INTEGER(4)             :: NTASKS,THISTASK
!     ******************************************************************
                              CALL TRACE$PUSH('PDOS')
call error$msg('routine marked for deletion')
call error$stop('pdos')
!
!     ==================================================================
!     ==  CALCULATE PARTIAL WAVE AMPLITUDES                           ==
!     ==================================================================
!     -- PDOS WILL BE CALCULATED UP TO ANGULAR MOMENTUM OF LXPDOS
      LXPDOS=2
      LMXPDOS=(LXPDOS+1)**2
      ALLOCATE(WOFL(LMXPDOS,NAT,NB,NKPT,NSPIN)) 
!
!     ==================================================================
!     ==  CALCULATE PARTIAL WAVE AMPLITUDES                           ==
!     ==================================================================
      CALL PDOS_WOFL(NRX,LNXX,LMNXX,NX,NKPTX,NB,NKPT,NSPIN &
     &              ,ISPECIES &
     &              ,NSP,AEZ,R1,DEX,NR,LNX,LOX,AEPHI &
     &              ,ZR,NAT,NBFNLX,FNL,LMXPDOS,WOFL)
!
!     ==================================================================
!     ==  WRITE PARTIAL WAVE AMPLITUDES                               ==
!     ==================================================================
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(TAMPFUL.AND.(THISTASK.EQ.1)) THEN
        CALL FILEHANDLER$UNIT('PDOS',NFIL)
        REWIND NFIL
        WRITE(NFIL,FMT='(4I5)')NAT,NB,NKPT,NSPIN
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          CALL ATOMLIST$GETCH('NAME',IAT,STRING)
          CALL ATOMLIST$GETR8A('R(0)',IAT,3,RPOS)
          LMAX=0
          DO LN=1,LNX(ISP)
            LMAX=MAX(LMAX,LOX(LN,ISP))
          ENDDO 
          LMX=(LMAX+1)**2
          WRITE(NFIL,FMT='(A32,I5,4F10.5)')STRING,LMX,AEZ(ISP),RPOS(:)
          DO IB=1,NB
            DO IKPT=1,NKPT
              DO ISPIN=1,NSPIN
                WRITE(NFIL,FMT='(F10.5,5X,9F10.6)')EIG(IB,IKPT,ISPIN) &
     &              ,(WOFL(LM,IAT,IB,IKPT,ISPIN),LM=1,LMX)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!       CALL FLUSH_(NFIL)
        CALL FILEHANDLER$CLOSE('PDOS')
      END IF
!ALL THE REST IS GARBADGE!!!!!!
!
!     ==================================================================
!     ==  WRITE ANGULAR MOMENTUM WEIGHTS                              ==
!     ==================================================================
!                              CALL TRACE$PASS('TWEIGHT')
!      IF(TWEIGHT) THEN
!        ALLOCATE(WEIGHT(LMXPDOS))
!        DO IAT=1,NAT
!          ISP=ISPECIES(NAT)
!          DO LM=1,LMXPDOS
!            WEIGHT(LM)=0.D0
!          ENDDO
!          DO IKPT=1,NKPT
!            DO ISPIN=1,NSPIN
!              DO IB=1,NB
!                DO LM=1,LMXPDOS
!                  WEIGHT(LM)=WEIGHT(LM)+WOFL(LM,IAT,IB,IKPT,ISPIN)**2
!                ENDDO
!              ENDDO
!            ENDDO
!          ENDDO
!          DO LM=1,LMXPDOS
!            WEIGHT(LM)=WEIGHT(LM)/DBLE(NKPT*NSPIN)
!          ENDDO
!        ENDDO
!        DEALLOCATE(WEIGHT)
!      END IF
!!
!!     ==================================================================
!!     ==  WRITE ANGULAR MOMENTUM WEIGHTS                              ==
!!     ==================================================================
!                              CALL TRACE$PASS('TAMP')
!      IF(TAMP) THEN
!        LISTPERCENT=99
!!       WRITE(NFIL,FMT='("THE WEIGHT OF THE LISTED WAVE FUNCTION "
!!    &                  ,"AMPLITUDES SUM UP TO ",I2," PERCENT "
!!    &                  ,"OF THE TOTAL WEIGHT")')LISTPERCENT
!        ALLOCATE(WORK(LMXPDOS*NAT))
!        ALLOCATE(IWORK(LMXPDOS*NAT))
!        DO ISPIN=1,NSPIN
!          DO IKPT=1,NKPT
!            DO IB=1,NB
!              TOTWEIGHT=0.D0
!              DO IAT=1,NAT
!                DO LM=1,LMXPDOS
!                  TOTWEIGHT=TOTWEIGHT+WOFL(LM,IAT,IB,IKPT,ISPIN)**2
!                ENDDO
!              ENDDO
!              DO IAT=1,NAT
!                DO LM=1,LMXPDOS
!                  II=LMXPDOS*(IAT-1)+LM
!                  WORK(II)=WOFL(LM,IAT,IB,IKPT,ISPIN)**2/TOTWEIGHT
!                ENDDO
!              ENDDO
!              CALL DSORTX(WORK,1,LMXPDOS*NAT,IWORK)
!              NORB=NAT*LMXPDOS
!              SUM=0.D0
!              DO I=1,NAT*LMXPDOS
!                II=NAT*LMXPDOS+1-I                
!                SUM=SUM+WORK(II)
!                IF(SUM.GT.DBLE(LISTPERCENT)*0.01D0) THEN
!                  NORB=I
!                  GOTO 100
!                END IF
!              ENDDO
! 100          CONTINUE
!!             WRITE(NFIL,FMT='(72("=")
!!    &               ,T5,2X,"SPIN=",I2,"IK=",I2,"BAND=",I4
!!    &               ,"TOTAL WEIGHT=",F8.5,2X)')
!!    &               ISPIN,IKPT,IB,TOTWEIGHT
!              DO I=1,NORB
!                II=LMXPDOS*NAT+1-I            
!                ISVAR=IWORK(II)
!                IAT=(ISVAR-1)/LMXPDOS+1
!                LM=ISVAR-LMXPDOS*(IAT-1)
!!               WRITE(NFIL,FMT='("ATOM",I4," LM= ",I4
!!    &               ," AMPLITUDE=",F8.3)')
!!    &               IAT,LM,WOFL(LM,IAT,IB,IKPT,ISPIN)   
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO
!        DEALLOCATE(WORK)
!        DEALLOCATE(IWORK)
!      END IF
!!
!!     ==================================================================
!!     ==  WRITE PROJECTED DENSITY OF STATES                           ==
!!     ==================================================================
!                              CALL TRACE$PASS('TPDOS')
!      IF(TPDOS) THEN
!        CALL CONSTANTS('EV',EV)
!        EMIN=-0.5D0
!        EMAX=0.1D0
!        NE=100
!        WIDTH=0.2D0*EV
!        ALLOCATE(WORK(NE))
!        DO IE=1,NE
!          E=EMIN+(IE-1)/DBLE(NE-1)*(EMAX-EMIN)
!       ENDDO
!       DEALLOCATE(WORK)
!     END IF
!
      DEALLOCATE(WOFL)
                              CALL TRACE$POP
      RETURN
!
!     ==================================================================
!     ==   SELECT PRINTOUT                                            ==
!     ==================================================================
      ENTRY PDOS$SELECT(TPDOS1,TAMP1,TWEIGHT1,TAMPFUL1)
        TPDOS=TPDOS1
        TAMP=TAMP1
        TAMPFUL=TAMPFUL1
        TWEIGHT=TWEIGHT1
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS_WOFL(NRX,LNXX,LMNXX,NX,NKPTX,NB,NKPT,NSPIN &
     &               ,ISPECIES &
     &               ,NSP,AEZ,R1,DEX,NR,LNX,LOX,AEPHI &
     &               ,ZR,NAT,NBFNLX,FNL,LMXPDOS,WOFL)
!     ******************************************************************
!     **                                                              **
!     **  PDOS_WOFL CALCULATES                                        **
!     **  THE WAVE FUNCTION AMPLITUDE FOR ENERGY DEPENDENT            **
!     **     PARTIAL WAVES WHICH ARE NORMALIZED WITHIN AN ATOMIC      **
!     **     SPHERE, AND WHOSE SIGN IS DEFINED SUCH THAT THE          **
!     **     PARTIAL WAVES FOR THE VALENCE STATES ARE POSITIVE        **
!     **     IN THE BINDING REGION.                                   **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NSP
      INTEGER(4),INTENT(IN) :: LNXX
      INTEGER(4),INTENT(IN) :: NRX
      INTEGER(4),INTENT(IN) :: NKPT
      INTEGER(4),INTENT(IN) :: NSPIN
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LMXPDOS
      INTEGER(4),INTENT(IN) :: NBFNLX
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: NKPTX
      INTEGER(4),INTENT(IN) :: LNX(NSP)
      INTEGER(4),INTENT(IN) :: LOX(LNXX,NSP)
      REAL(8)   ,INTENT(IN) :: AEZ(NSP)
      REAL(8)   ,INTENT(IN) :: AEPHI(NRX,LNXX,NSP)
      REAL(8)   ,INTENT(IN) :: ZR(NX,NX,NKPT,NSPIN)
      REAL(8)   ,INTENT(IN) :: FNL(NAT,NBFNLX,LMNXX,NKPTX,NSPIN)
      INTEGER(4),INTENT(IN) :: ISPECIES(NAT)
      REAL(8)   ,INTENT(OUT):: WOFL(LMXPDOS,NAT,NB,NKPT,NSPIN)
      REAL(8)               :: OVCOV(LNXX,LNXX,NSP) 
      REAL(8)               :: AEPHIR1(LNXX,NSP) 
      INTEGER(4)            :: ISPIN,IKPT,IB,IAT,ISP
!     ******************************************************************
                              CALL TRACE$PUSH('PDOS_WOFL')
call error$msg('routine marked for deletion')
call error$stop('STATEANALYSIS_wofl')
!
!     ==================================================================
!     ==  CALCULATE OVERLAP MATRIX OF AE PARTIAL WAVES                ==
!     ==================================================================
      CALL PDOS_OVERCOV(NRX,LNXX &
     &                ,NSP,AEZ,R1,DEX,NR,LNX,LOX,AEPHI,OVCOV,AEPHIR1)
!
!     ==================================================================
!     ==  CALCULATE PARTIAL WAVE AMPLITUDES                           ==
!     ==================================================================
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NB
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              CALL PDOS_STATE(LMXPDOS,NAT,NX,LNXX,LMNXX &
     &                 ,NSP,LNX,LOX,NB &
     &                 ,FNL(1,1,1,IKPT,ISPIN),ZR(1,1,IKPT,ISPIN) &
     &                 ,OVCOV,AEPHIR1 &
     &                 ,ISP,IAT,IB,WOFL(1,IAT,IB,IKPT,ISPIN))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                              CALL TRACE$POP
      RETURN
      END
!
!     ...................................................AUTOPI.........
      SUBROUTINE PDOS_STATE(LMXPDOS,NAT,NX,LNXX,LMNXX &
     &                  ,NSP,LNX,LOX,NB,FNL,ZR,OVCOV,AEPHIR1 &
     &                  ,ISP,IAT,IB,WOFL)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE NUMBER OF ANGULAR MOMENTUM WEIGHTS WOFL      **
!     **  OF STATE "IB" ON ATOM "(ISP,IAT)"                           **
!     **                                                              **
!     **  ZR  IS THE UNITARY TRANSFROMATION ON THE EIGENSTATES        **
!     **  FNL ARE THE PROJECTRION ON THE PARTIAL WAVES                **
!     **                                                              **
!     **  OUTPUT:                                                     **
!     **  -------                                                     **
!     **  WOFL IS THE NUMBER OF ELECTRONS(PER SPIN) IN ORBITAL LM     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMXPDOS
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LNXX
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: NSP
      INTEGER(4),INTENT(IN) :: LNX(NSP)
      INTEGER(4),INTENT(IN) :: LOX(LNXX,NSP)
      INTEGER(4),INTENT(IN) :: NB
      REAL(8)   ,INTENT(IN) :: FNL(NAT,NB,LMNXX)
      REAL(8)   ,INTENT(IN) :: ZR(NX,NX)
      REAL(8)   ,INTENT(IN) :: OVCOV(LNXX,LNXX,NSP)
      REAL(8)   ,INTENT(IN) :: AEPHIR1(LNXX,NSP)
      INTEGER(4),INTENT(IN) :: ISP
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: IB
      REAL(8)   ,INTENT(OUT):: WOFL(LMXPDOS)
      REAL(8)               :: CWAVE(LMNXX)
      INTEGER(4)            :: LN,LMN,IB2,LM,LMN1P,LN1,L1,LMN2P
      INTEGER(4)            :: L2,LMP,LMN1,LMN2,L0,LMNP,L,M0,M,LN2
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: LMAX
      REAL(8)               :: SVAR
      REAL(8)               :: VAL
!     ******************************************************************
call error$msg('routine marked for deletion')
call error$stop('pdos_state')
!
!     ==================================================================
!     == TRANSFORM FNL ON EIGENSTATE                                  ==
!     ==================================================================
      LMNX=0
      DO LN=1,LNX(ISP)
        LMNX=LMNX+2*LOX(LN,ISP)+1
      ENDDO
!
      DO LMN=1,LMNX
        SVAR=0.D0
        DO IB2=1,NB
          SVAR=SVAR+FNL(IAT,IB2,LMN)*ZR(IB2,IB)
        ENDDO
        CWAVE(LMN)=SVAR
      ENDDO
!
!     ==================================================================
!     == CALCULATE NORM                                               ==
!     ==================================================================
!     == WOFL(LM) IS THE NUMBER OF ELECTRONS WITHIN THE ATOMIC SPHERE ==
!     == FROM THE ANGULAR MOMENTUM CHANNEL LM                         ==
!     ==================================================================
      DO LM=1,LMXPDOS
        WOFL(LM)=0.D0
      ENDDO
! 
      LMN1P=0
      DO LN1=1,LNX(ISP)
        L1=LOX(LN1,ISP)
        LMN2P=0
        DO LN2=1,LNX(ISP)
          L2=LOX(LN2,ISP)
          IF(L1.EQ.L2) THEN
            LMP=L1**2
            DO M=1,2*L1+1
              LMN1=LMN1P+M
              LMN2=LMN2P+M
              LM=LMP+M
              IF(LM.LE.LMXPDOS) THEN
                WOFL(LM)=WOFL(LM) &
     &                +OVCOV(LN1,LN2,ISP)*CWAVE(LMN1)*CWAVE(LMN2)
              END IF
            ENDDO
          END IF
          LMN2P=LMN2P+2*L2+1
        ENDDO
        LMN1P=LMN1P+2*L1+1
      ENDDO
!
!     ==================================================================
!     == WOFL(LM) IS CONVERTED TO THE SQUARE ROOT OF THE              ==
!     ==  ANGULAR MOMENTUM WEIGHT FOLLOWING THE CONVENTION            ==
!     == THAT THE PARTIAL WAVE OF A IS POSITIVE IN THE BONDONG REGION ==
!     ==================================================================
      LMAX=0
      DO LN=1,LNX(ISP)
        LMAX=MAX(LMAX,LOX(LN,ISP))
      ENDDO 
!
      LM=0
      DO L0=0,LMAX
        DO M0=1,2*L0+1
          LM=LM+1
          VAL=0.D0
          LMNP=0
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            IF(L.EQ.L0) THEN
              LMN=LMNP+M0
              VAL=VAL+CWAVE(LMN)*AEPHIR1(LN,ISP)
            END IF
            LMNP=LMNP+2*L+1
          ENDDO
          WOFL(LM)=SIGN(SQRT(WOFL(LM)),VAL)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...................................................OVERCOV........
      SUBROUTINE PDOS_OVERCOV(NRX,LNXX &
     &                  ,NSP,AEZ,R1,DEX,NR,LNX,LOX,AEPHI,OVCOV,AEPHIR1)
!     ******************************************************************
!     **                                                              **
!     **  OVERCOV EVALUATES THE NORM OF THE PARTIAL WAVES WITHIN      **
!     **  THE COVALENT RADIUS, AND  CHOOSES A SIGN CONVENTION         **
!     **                                                              **
!     ******************************************************************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NRX
      INTEGER(4),INTENT(IN) :: LNXX
      INTEGER(4),INTENT(IN) :: NSP
      REAL(8)   ,INTENT(IN) :: AEZ(NSP)
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX(NSP)
      INTEGER(4),INTENT(IN) :: LOX(LNXX,NSP) 
      REAL(8)   ,INTENT(IN) :: AEPHI(NRX,LNXX,NSP)
      REAL(8)   ,INTENT(OUT):: OVCOV(LNXX,LNXX,NSP)
      REAL(8)   ,INTENT(OUT):: AEPHIR1(LNXX,NSP)
      LOGICAL(4)            :: TPR=.FALSE.
      CHARACTER(10)         :: CONFIG
      CHARACTER(2)          :: SYMBOL
      REAL(8)               :: AUX(NRX)
      REAL(8)               :: XEXP
      INTEGER(4)            :: ISP,LN,L,LN1,LN2,L1,L2,IR
      INTEGER(4)            :: IZ             ! ATOMIC NUMBER
      INTEGER(4)            :: NODE           ! NUMBER OF NODES 
      REAL(8)               :: RCOV,XCOV
      INTEGER(4)            :: IOUT,XOUT
      REAL(8)               :: RI
      REAL(8)               :: VAL1,VAL2
!     ******************************************************************
call error$msg('routine marked for deletion')
call error$stop('STATEANALYSIS_overcov')
      XEXP=DEXP(DEX)
      DO ISP=1,NSP
        IZ=NINT(AEZ(ISP))
!
!       ================================================================
!       == CALCULATE AEPHIR1                                          ==
!       ================================================================
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(L.EQ.0) THEN
!           CALL PERIODICTABLEI4(IZ,'NODE(S)',NODE)
            CALL PERIODICTABLE$GET(IZ,'#NODES(S)',NODE)
          ELSE IF(L.EQ.1) THEN
!           CALL PERIODICTABLEI4(IZ,'NODE(P)',NODE)
            CALL PERIODICTABLE$GET(IZ,'#NODES(P)',NODE)
          ELSE IF(L.EQ.2) THEN
!           CALL PERIODICTABLEI4(IZ,'NODE(D)',NODE)
            CALL PERIODICTABLE$GET(IZ,'#NODES(D)',NODE)
          ELSE IF(L.EQ.3) THEN
!           CALL PERIODICTABLEI4(IZ,'NODE(F)',NODE)
            CALL PERIODICTABLE$GET(IZ,'#NODES(F)',NODE)
          ELSE
            NODE=0
          END IF
          AEPHIR1(LN,ISP)=AEPHI(1,LN,ISP)*(-1.D0)**NODE
        ENDDO
!
!       ================================================================
!       == CALCULATE OVCOV                                            ==
!       ================================================================
!       CALL PERIODICTABLER8(IZ,'R(ASA)',RCOV)
        CALL PERIODICTABLE$GET(IZ,'R(ASA)',RCOV)
        XCOV=1.D0+DLOG(RCOV/R1)/DEX
        IOUT=INT(XCOV)
        XOUT=DBLE(IOUT)
        DO LN1=1,LNXX
          DO LN2=1,LNXX
            OVCOV(LN1,LN2,ISP)=0.D0
          ENDDO
        ENDDO
        DO LN1=1,LNX(ISP)
          L1=LOX(LN1,ISP)
          DO LN2=LN1,LNX(ISP)
            L2=LOX(LN2,ISP)
            IF(L1.EQ.L2) THEN
              RI=R1/XEXP 
              DO IR=1,NR
                RI=RI*XEXP
                AUX(IR)=RI**2*AEPHI(IR,LN1,ISP)*AEPHI(IR,LN2,ISP)
              ENDDO
              CALL RADIAL$INTEGRAL(R1,DEX,IOUT,AUX,VAL1)
              CALL RADIAL$INTEGRAL(R1,DEX,IOUT+1,AUX,VAL2)
              OVCOV(LN1,LN2,ISP)=VAL1+(XCOV-XOUT)*(VAL2-VAL1)
            ELSE
              OVCOV(LN1,LN2,ISP)=0.D0
            END IF
            OVCOV(LN2,LN1,ISP)=OVCOV(LN1,LN2,ISP)
          ENDDO
        ENDDO   
      ENDDO
!
!     ==================================================================    
!     ==  PRINT                                                       ==    
!     ==================================================================    
      IF(TPR) THEN
        DO ISP=1,NSP
          WRITE(*,FMT='("<PHI|PHI> (R<RCOV) FOR SPECIES ",I5)')ISP
          DO LN1=1,LNX(ISP)
            WRITE(*,FMT='(10F10.5)')(OVCOV(LN1,LN2,ISP),LN2=1,LNX(ISP))
          ENDDO
        ENDDO
      END IF
      RETURN
      END
