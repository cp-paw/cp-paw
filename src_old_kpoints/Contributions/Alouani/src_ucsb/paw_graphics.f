!.......................................................................
MODULE GRAPHICS_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: GRAPHICS                                                   **
!**                                                                   **
!**  PURPOSE: COLLECT WAVE FUNCTIONS OR DENSITIES AND PREPARE         **
!**    FILE FOR PLOTTING                                              **
!**                                                                   **
!**  DOES NOT WORK IF RDYN IS OFF!                                    **
!**                                                                   **
!**  OPTIONS IN FUTURE                                                **
!**    SPINDENSITY                                                    **
!**    TOTALDENSITY                                                   **
!**    STATEDENSITY                                                   **
!**    WAVE                                                           **
!**                                                                   **
!**  LIST STATES IN A SERIES OF STATE RANGES B1,B2,K1,K2,S1,S2        **
!**  (IF B1,K1,OR S1 IS ZERO ALL POSSIBILITIES CONTRIBUTE)            **
!**  FOR A SINGLE STATE B2=B1 ETC                                     **
!**                                                                   **
!**       ORIGINAL VERSION: PETER MARGL                               **
!**       MODIFIED VERSION:                                           **
!**           PETER E. BLOECHL, IBM ZURICH RESEARCH LABORATORY (1996) **
!***********************************************************************
USE LINKEDLIST_MODULE
TYPE STATE_TYPE
 INTEGER(4)              :: IB
 INTEGER(4)              :: IK
 INTEGER(4)              :: IS
 INTEGER(4)              :: F
 TYPE(STATE_TYPE),POINTER:: NEXT
END TYPE STATE_TYPE
TYPE IMAGE_TYPE
  CHARACTER(256)           :: TITLE
  CHARACTER(256)           :: FILE
  CHARACTER(8)             :: TYPE   ! CAN BE 'DENSITY'
  TYPE(STATE_TYPE),POINTER :: STATE
  TYPE(IMAGE_TYPE),POINTER :: NEXT
END TYPE IMAGE_TYPE
!
!=======================================================================
!== NOW THE DATA                                                      ==
!=======================================================================
TYPE (LL_TYPE)         :: LL_GRPHCS
LOGICAL(4)             :: TINI=.FALSE.
LOGICAL(4)             :: Twake=.FALSE.
COMPLEX(8),ALLOCATABLE :: PWPOT(:)
REAL(8),ALLOCATABLE    :: AE1CPOT(:,:,:)
REAL(8),ALLOCATABLE    :: PS1CPOT(:,:,:)
INTEGER(4)             :: LMRXX=0    !initially not set
INTEGER(4)             :: Nrx
INTEGER(4)             :: NGL
END MODULE GRAPHICS_MODULE
!
!     ..................................................................
      SUBROUTINE GRAPHICS$GETLIST(LL_GRPHCS_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE GRAPHICS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(OUT) :: LL_GRPHCS_
      LOGICAL(4)                :: TCHK
!     ******************************************************************
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        CALL LINKEDLIST$NEW(LL_GRPHCS)
        CALL LINKEDLIST$SET(LL_GRPHCS,'LINKEDLISTNAME',0,'GRAPHICS')
      END IF
      LL_GRPHCS_=LL_GRPHCS
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS$SETL4(ID_,VAL_)
      USE GRAPHICS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'WAKE') THEN
        TWAKE=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('GRAPHICS$SETL4')
      END IF
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE GRAPHICS$SETPWPOT(ID,NGL_,VHARTREE)
!     ********************************************************************
!     **  GET SECOND DERIVATIVE OF THE RADIAL POTENTIAL AT THE ORIGIN   **
!     ********************************************************************
      USE GRAPHICS_MODULE
      IMPLICIT NONE
      CHARACTER(32),intent(in) :: ID
      INTEGER(4)   ,intent(in) :: NGL_
      COMPLEX(8)   ,intent(in) :: VHARTREE(NGL_)
      CHARACTER(32)            :: TYPE
      INTEGER(4)               :: IPIC
      INTEGER(4)               :: NPICS
      logical(4)               :: tpot
!     ********************************************************************
      IF (.NOT.Tini) RETURN
      print*,'clemens: INTO SETPWPOT',TWAKE,NGL_
!check if active
      IF (.NOT.TWAKE) RETURN
!set potential
      NGL=NGL_
print*,'marke 1'
      CALL LINKEDLIST$SELECT(LL_GRPHCS,'~')
print*,'marke 2'
      CALL LINKEDLIST$NLISTS(LL_GRPHCS,'IMAGE',NPICS)
print*,'marke 3'
      tpot=.false.
      DO IPIC=1,NPICS
        CALL LINKEDLIST$SELECT(LL_GRPHCS,'~')
        CALL LINKEDLIST$SELECT(LL_GRPHCS,'IMAGE',IPIC)
        CALL LINKEDLIST$GET(LL_GRPHCS,'TYPE',1,TYPE)
        if(type.eq.'POT') then
          if(tpot) cycle     ! ensure that potential will be set only once
          tpot=.true.
          IF(.NOT.ALLOCATED(PWPOT)) THEN
            ALLOCATE(PWPOT(NGL))
          END IF
          PWPOT(:)=VHARTREE(:)
          print*,'clemens: PS contribution set',ID
        end if
      enddo
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE GRAPHICS$SET1CPOT(IDENT_,IAT_,R1,DEX,NR,NRX_,LMRX,POT)
!     ********************************************************************
!     **  USE 1-CENTER POTENTIAL FOR ELECTRIC FIELD GRADIENTS          **
!     ********************************************************************
      USE GRAPHICS_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: IDENT_  ! CAN BE 'AE' OR 'PS' 
      INTEGER(4)   ,INTENT(IN) :: IAT_    ! ATOM INDEX (SEE ATOMLIST)
      REAL(8)      ,INTENT(IN) :: R1
      REAL(8)      ,INTENT(IN) :: DEX
      INTEGER(4)   ,INTENT(IN) :: NR,NRX_
      INTEGER(4)   ,INTENT(IN) :: LMRX
      REAL(8)      ,INTENT(IN) :: POT(NRX_,LMRX)
      REAL(8)                  :: XEXP
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: RI
      INTEGER(4)               :: LM,IR
      CHARACTER(32)            :: TYPE
      INTEGER(4)               :: NPICS
      INTEGER(4)               :: IPIC
      INTEGER(4)               :: NAT
      logical(4)               :: tpot
!     ********************************************************************
      IF(.NOT.Tini) RETURN
!     ===================================================================
!     == CHECK WHETHER ACTION IS REQUIRED                              ==
!     ===================================================================
!check if active
      print*,'clemens: INTO SET1CPOT'
      IF(.NOT.TWAKE) RETURN
!set potential
      CALL LINKEDLIST$SELECT(LL_GRPHCS,'~')
      CALL LINKEDLIST$NLISTS(LL_GRPHCS,'IMAGE',NPICS)
      tpot=.false.
      DO IPIC=1,NPICS
        CALL LINKEDLIST$SELECT(LL_GRPHCS,'~')
        CALL LINKEDLIST$SELECT(LL_GRPHCS,'IMAGE',IPIC)
        CALL LINKEDLIST$GET(LL_GRPHCS,'TYPE',1,TYPE)
        IF(TYPE.EQ.'POT') THEN
          if(tpot) cycle ! ensure that potential will be set only once
          tpot=.true.
          IF(.NOT.ALLOCATED(AE1CPOT)) THEN
            nrx=nrx_
            CALL SETUP$GETI4('LMRXX',LMRXX)
            CALL ATOMLIST$NATOM(NAT)
            ALLOCATE(AE1CPOT(NRX,LMRXX,NAT))
            ALLOCATE(PS1CPOT(NRX,LMRXX,NAT))
            ae1cpot=0.d0
            ps1cpot=0.d0
          END IF
          IF(IDENT_.EQ.'AE') THEN
            AE1CPOT(:,:,IAT_)=0.D0
            AE1CPOT(:,1:LMRX,IAT_)=POT
          ELSE IF(IDENT_.EQ.'PS') THEN
            PS1CPOT(:,:,IAT_)=0.D0
            PS1CPOT(:,1:LMRX,IAT_)=POT
          ELSE
            CALL ERROR$MSG('ID MUST BE WITHER "AE" OR "PS"')
            CALL ERROR$STOP('GRAPHICS$SET1CPOT')
          END IF
        end if
      enddo
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS$PLOT
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      USE GRAPHICS_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NAT
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: NR1,NR1L,NR2,NR3
      INTEGER(4)            :: NNR,NNRL
      INTEGER(4)            :: NR1START
      INTEGER(4)            :: NB
      INTEGER(4)            :: NKPT
      INTEGER(4)            :: NSPIN
      INTEGER(4)            :: LMNXX
      REAL(8)   ,ALLOCATABLE:: OCC(:,:,:)  !(NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE:: WAVE(:,:,:) !(NR1L,NR2,NR3) later refined
      INTEGER(4)            :: NTASKNUM,NTASKID
      REAL(8)   ,ALLOCATABLE:: WORK(:)           
      CHARACTER(32)         :: STRING
      INTEGER(4)            :: NSTATE
      INTEGER(4)            :: IB,IKPT,ISPIN
      INTEGER(4)            :: IB1,IKPT1,ISPIN1,IB2,IKPT2,ISPIN2
      INTEGER(4)            :: IR1,IR2,IR3
      INTEGER(4)            :: ISTATE
      LOGICAL(4)            :: TOCC
      LOGICAL(4)            :: TSPIN
      REAL(8)               :: FAC
      INTEGER(4)            :: IAT,ISP,LN,LMN1,LMN2
      REAL(8)               :: BAS(3)
      REAL(8)               :: SVAR
      INTEGER(4)            :: LNX
      REAL(8)               :: R1
      REAL(8)               :: DEX
      INTEGER(4)            :: NR
      INTEGER(4)            :: LX,LMXX
      INTEGER(4),ALLOCATABLE:: LOX(:)    !LNX
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:) !NR,LNX
      REAL(8)   ,ALLOCATABLE:: PSPHI(:,:) !NR,LNX
      REAL(8)   ,ALLOCATABLE:: PROJ(:)    !LMNXX
      REAL(8)   ,ALLOCATABLE:: DRHOL(:,:) !(NR,LMXX)
      REAL(8)   ,ALLOCATABLE:: DENMAT(:,:) !(LMNXX,LMNXX)
      CHARACTER(32),ALLOCATABLE :: ATOMNAME(:)   !NAT
      REAL(8)   ,ALLOCATABLE:: POS(:,:)   !(3,NAT)
      REAL(8)   ,ALLOCATABLE:: Z(:)
      REAL(8)   ,ALLOCATABLE:: Q(:)
      LOGICAL(4)            :: SAVETRAWSTATES !USED TO RESTORE ORIGINAL
                                              ! STATE OF WAVES OBJECT
      INTEGER(4)            :: NPICS
      INTEGER(4)            :: IPIC
      CHARACTER(256)        :: TITLE
      CHARACTER(256)        :: FILE
      CHARACTER(8)          :: TYPE
      CHARACTER(8)          :: FORMAT
      CHARACTER(8)          :: WEIGHTING
      LOGICAL(4)            :: TIM
      LOGICAL(4)            :: TCHK
      REAL(8)   ,ALLOCATABLE:: WAVEBIG(:,:,:)
      REAL(8)   ,ALLOCATABLE:: WORK1(:,:,:)
      INTEGER(4),parameter  :: FACT=2
!     ******************************************************************
      IF(.NOT.TINI) RETURN
                              CALL TRACE$PUSH('GRAPHICS$PLOT')
!
!     =================================================================
!     ==  GET GENERIC FROM ATOM OBJECT                               ==
!     =================================================================
      CALL ATOMLIST$NATOM(NAT)
      CALL CELL$GETR8A('T(0)',9,RBAS)
      ALLOCATE(POS(3,NAT))
      ALLOCATE(ATOMNAME(NAT))
      ALLOCATE(Z(NAT))
      ALLOCATE(Q(NAT))
      DO IAT=1,NAT
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME(IAT))
        CALL ATOMLIST$GETR8A('R(0)',IAT,3,POS(:,IAT))
        CALL ATOMLIST$GETR8('Z',IAT,Z(IAT))
        CALL ATOMLIST$GETR8('Q',IAT,Q(IAT))
      ENDDO
!
!     =================================================================
!     ==  GET GENERIC INFORMATION ABOUT NUMBER AND SIZE OF THE       ==
!     ==  PSEUDO WAVE FUNCTIONS                                      ==
!     =================================================================
      CALL WAVES$GETI4('NR1',NR1)
      CALL WAVES$GETI4('NR1L',NR1L)
      CALL WAVES$GETI4('NR2',NR2)
      CALL WAVES$GETI4('NR3',NR3)
      NNR=NR1*NR2*NR3
      NNRL=NR1L*NR2*NR3
      CALL WAVES$GETI4('NR1START',NR1START)
      CALL WAVES$GETI4('NB',NB)
      CALL WAVES$GETI4('NKPT',NKPT)
      CALL WAVES$GETI4('NSPIN',NSPIN)
      CALL WAVES$GETL4('RAWSTATES',SAVETRAWSTATES)      
!
!     =================================================================
!     ==  GET GENERIC FROM SETUPS OBJECT                             ==
!     =================================================================
      CALL SETUP$LMNXX(LMNXX)
      ALLOCATE(OCC(NB,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      ALLOCATE(WAVE(NR1L,NR2,NR3))
      ALLOCATE(WAVEBIG(FACT*NR1,FACT*NR2,FACT*NR3))
!
!     =================================================================
!     ==  LOOP OVER ALL IMAGES                                       ==
!     =================================================================
                              CALL TRACE$PASS('MAKE PICTURES')
      CALL LINKEDLIST$SELECT(LL_GRPHCS,'~')
      CALL LINKEDLIST$NLISTS(LL_GRPHCS,'IMAGE',NPICS)
      DO IPIC=1,NPICS
                              CALL TRACE$PASS('NEW PICTURE')
        CALL LINKEDLIST$SELECT(LL_GRPHCS,'~')
        CALL LINKEDLIST$SELECT(LL_GRPHCS,'IMAGE',IPIC)
        CALL LINKEDLIST$GET(LL_GRPHCS,'TITLE',1,TITLE)
        CALL LINKEDLIST$GET(LL_GRPHCS,'FILE',1,FILE)
        CALL LINKEDLIST$GET(LL_GRPHCS,'TYPE',1,TYPE)
        CALL LINKEDLIST$GET(LL_GRPHCS,'FORMAT',1,FORMAT)
        CALL LINKEDLIST$GET(LL_GRPHCS,'WEIGHTING',1,WEIGHTING)
!
!       =================================================================
!       ==  check if it is potential                                   ==
!       =================================================================
        IF(TYPE(1:3).EQ.'POT') THEN
           CALL GRAPHICS_CREATEPOT(IPIC,FACT)
           CYCLE
        END IF
!
!       =================================================================
!       ==  SWITCH WAVES OBJECT TO RAW OR EIGEN STATES                 ==
!       =================================================================
        IF(TRIM(WEIGHTING).EQ.'TOTAL'.OR.TRIM(WEIGHTING).EQ.'SPIN') THEN 
          CALL WAVES$SETL4('RAWSTATES',.TRUE.)
PRINT*,'TITLE ',TRIM(TITLE),.TRUE.
        ELSE
          CALL WAVES$SETL4('RAWSTATES',.FALSE.)
PRINT*,'TITLE ',TRIM(TITLE),.FALSE.
        END IF
!
!       ================================================================
!       ==  PSEUDO WAVE FUNCTIONS/DENSITIES                           ==
!       ================================================================
                              CALL TRACE$PASS('PSEUDO WAVE FUNCTIONS')
        IF(TRIM(TYPE).EQ.'WAVE') THEN
          CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE')
          CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB)
          CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT)
          CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN)
          CALL LINKEDLIST$EXISTD(LL_GRPHCS,'TIM',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_GRPHCS,'TIM',1,TIM)
          ELSE
            TIM=.FALSE.
          END IF
          CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
          CALL WAVES$SETI4('IB',IB)
          CALL WAVES$SETI4('IKPT',IKPT)
          CALL WAVES$SETI4('ISPIN',ISPIN)
          CALL WAVES$SETL4('TIM',TIM)
          CALL WAVES$GETR8A('PSPSI',NNRL,WAVE)
        ELSE IF(TYPE(1:7).EQ.'DENSITY') THEN
          DO IR3=1,NR3            
            DO IR2=1,NR2            
              DO IR1=1,NR1L
                WAVE(IR1,IR2,IR3)=0.D0
              ENDDO
            ENDDO
          ENDDO
          IF(TRIM(WEIGHTING).EQ.'TOTAL') THEN 
            TOCC=.TRUE.
            TSPIN=.FALSE.
          ELSE IF(TRIM(WEIGHTING).EQ.'SPIN') THEN
            TOCC=.TRUE.
            TSPIN=.TRUE.
          ELSE IF(TRIM(WEIGHTING).EQ.'NONE') THEN
            TOCC=.FALSE.
            TSPIN=.FALSE.
          ELSE
            CALL ERROR$MSG('INVALID OPTION OF WEIGHTING')
            CALL ERROR$CHVAL('WEIGHTING',WEIGHTING)
            CALL ERROR$MSG('ALLOWED VALUES: "TOTAL", "SPIN", "NONE".')
            CALL ERROR$STOP('GRAPHICS$PLOT')
          END IF
          IF(TRIM(FORMAT).EQ.'ALL'.OR.TRIM(FORMAT).EQ.'RANGE') THEN
            IF(TRIM(FORMAT).EQ.'RANGE') THEN
              CALL LINKEDLIST$NLISTS(LL_GRPHCS,'STATE',NSTATE)
              IF(NSTATE.NE.2) THEN
                CALL ERROR$MSG('RANGE MUST BE DEFINED BY AT LEAST TWO STATES')
                CALL ERROR$STOP('GRAPHICS$PLOT')
              END IF
              CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE',1)
              CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB1)
              CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT1)
              CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN1)
              CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
              CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE',2)
              CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB2)
              CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT2)
              CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN2)
              CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
            ELSE
              IB1=1
              IB2=NB
              IKPT1=1
              IKPT2=NKPT
              ISPIN1=1
              ISPIN2=NSPIN
            END IF
            DO ISPIN=ISPIN1,ISPIN2
              DO IKPT=IKPT1,IKPT2
                DO IB=IB1,IB2
                  IF(TOCC) THEN
                    FAC=OCC(IB,IKPT,ISPIN)/DBLE(NKPT)
                    IF(TSPIN.AND.ISPIN.EQ.2) FAC=-FAC
                  ELSE
                    FAC=1.D0
                  END IF
                  CALL GRAPHICS_ADDRHO(IB,IKPT,ISPIN,FAC,NNRL,WAVE)
                ENDDO
              ENDDO
            ENDDO
          ELSE  IF(FORMAT(1:4).EQ.'LIST') THEN
            CALL LINKEDLIST$NLISTS(LL_GRPHCS,'STATE',NSTATE)
            DO ISTATE=1,NSTATE
              CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE',ISTATE)
              CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB)
              CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT)
              CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN)
              CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
              IF(TOCC) THEN
                FAC=OCC(IB,IKPT,ISPIN)/DBLE(NKPT)
                IF(TSPIN.AND.ISPIN.EQ.2) FAC=-FAC
              ELSE
                FAC=1.D0
              END IF
              CALL GRAPHICS_ADDRHO(IB,IKPT,ISPIN,FAC,NNR,WAVE)
            ENDDO
          ELSE
            CALL ERROR$MSG('INVALID OPTION OF FORMAT')
            CALL ERROR$CHVAL('FORMAT',FORMAT)
            CALL ERROR$MSG('ALLOWED VALUES: "ALL", "RANGE", "LIST".')
            CALL ERROR$STOP('GRAPHICS$PLOT')
          END IF  
        ELSE
          CALL ERROR$MSG('INVALID OPTION OF TYPE')
          CALL ERROR$CHVAL('TYPE',TYPE)
          CALL ERROR$MSG('ALLOWED VALUES: "DENSITY", "WAVE".')
          CALL ERROR$STOP('GRAPHICS$PLOT')
        END IF  
!
!       ================================================================
!       ==  Expand to a finer r-grid                                  ==
!       ================================================================
        ALLOCATE(WORK1(NR1,NR2,NR3))
        CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WORK1)
        CALL GRAPHICS_REFINEGRID(NR1,NR2,NR3,FACT,WORK1,WAVEBIG)
        DEALLOCATE(WORK1)
!
!       ================================================================
!       ==  ONE-CENTER EXPANSIONS                                     ==
!       ================================================================
                              CALL TRACE$PASS('ONE-CENTER EXPANSIONS')
!         DO I=1,NR1L
!           DO J=1,NR2
!             DO K=1,NR3
!               WAVE(I,J,K)=0.D0
!             ENDDO
!           ENDDO
!         ENDDO

        DO IAT=1,NAT
          CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
!         __ GET POSITIONS______________________________________________
          CALL ATOMLIST$GETR8A('R(0)',IAT,3,BAS)
!         __ GET PARTIAL WAVES__________________________________________
          CALL SETUP$RADGRID(ISP,R1,DEX,NR)
          CALL SETUP$LNX(ISP,LNX)
          ALLOCATE(LOX(LNX))
          CALL SETUP$LOFLN(ISP,LNX,LOX)
          ALLOCATE(AEPHI(NR,LNX))
          ALLOCATE(PSPHI(NR,LNX))
          CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
          CALL SETUP$PSPARTIALWAVES(ISP,NR,LNX,PSPHI)
          CALL WAVES$SETI4('IAT',IAT)
!
!         IF(IAT.EQ.2) THEN
!           PRINT*,'WARNING!! GRAPHICS MODIFIED'
!           DO IR=1,NR
!             DO LN=1,LNX
!               PSPHI(IR,LN)=0.D0
!               AEPHI(IR,LN)=0.D0
!             ENDDO
!           ENDDO
!         END IF
!         DO IR=1,NR
!           DO LN=1,LNX
!             AEPHI(IR,LN)=0.D0
!           ENDDO
!         ENDDO
!
!         __ GET PROJECTIONS____________________________________________
          IF(TRIM(TYPE).EQ.'WAVE') THEN
            CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE')
            CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB)
            CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT)
            CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN)
            CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
            CALL WAVES$SETI4('IB',IB)
            CALL WAVES$SETI4('IKPT',IKPT)
            CALL WAVES$SETI4('ISPIN',ISPIN)
            CALL WAVES$SETL4('TIM',TIM)
            ALLOCATE(PROJ(LMNXX))
            CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROJ)
            LX=0
            DO LN=1,LNX
              LX=MAX(LX,LOX(LN))
            ENDDO
            LMXX=(LX+1)**2
            ALLOCATE(DRHOL(NR,LMXX))
            CALL GRAPHICS_1CWAVE(R1,DEX,NR,LNX,LOX,AEPHI,PSPHI,LMNXX &
     &                            ,PROJ,LMXX,DRHOL)
            DEALLOCATE(PROJ)
          ELSE IF(TRIM(TYPE).EQ.'DENSITY') THEN
            ALLOCATE(DENMAT(LMNXX,LMNXX))
            DO LMN1=1,LMNXX
              DO LMN2=1,LMNXX
                DENMAT(LMN1,LMN2)=0.D0
              ENDDO
            ENDDO
            IF(TRIM(FORMAT).EQ.'ALL'.OR.TRIM(FORMAT).EQ.'RANGE') THEN
              IF(TRIM(FORMAT).EQ.'RANGE') THEN
                CALL LINKEDLIST$NLISTS(LL_GRPHCS,'STATE',NSTATE)
                IF(NSTATE.NE.2) THEN
                  CALL ERROR$STOP('RANGE MUST BE DEFINED &
                                 &BY AT LEAST TWO STATES')
                  CALL ERROR$STOP('GRAPHICS$PLOT')
                END IF
                CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE',1)
                CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB1)
                CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT1)
                CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN1)
                CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
                CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE',2)
                CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB2)
                CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT2)
                CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN2)
                CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
              ELSE
                IB1=1
                IB2=NB
                IKPT1=1
                IKPT2=NKPT
                ISPIN1=1
                ISPIN2=NSPIN
              END IF
              DO ISPIN=ISPIN1,ISPIN2
                DO IKPT=IKPT1,IKPT2
                  DO IB=IB1,IB2
                    IF(TOCC) THEN
                      FAC=OCC(IB,IKPT,ISPIN)/DBLE(NKPT)
                      IF(TSPIN.AND.ISPIN.EQ.2) FAC=-FAC
                    ELSE
                      FAC=1.D0
                    END IF
                    CALL GRAPHICS_ADDDENMAT(IB,IKPT,ISPIN,FAC &
     &                           ,LMNXX,DENMAT)
                  ENDDO
                ENDDO
              ENDDO
            ELSE  IF(TRIM(FORMAT).EQ.'LIST') THEN
              CALL LINKEDLIST$NLISTS(LL_GRPHCS,'STATE',NSTATE)
              DO ISTATE=1,NSTATE
                CALL LINKEDLIST$SELECT(LL_GRPHCS,'STATE')
                CALL LINKEDLIST$GET(LL_GRPHCS,'IB',1,IB)
                CALL LINKEDLIST$GET(LL_GRPHCS,'IKPT',1,IKPT)
                CALL LINKEDLIST$GET(LL_GRPHCS,'ISPIN',1,ISPIN)
                CALL LINKEDLIST$SELECT(LL_GRPHCS,'..')
                IF(TOCC) THEN
                  FAC=OCC(IB,IKPT,ISPIN)/DBLE(NKPT)
                  IF(TSPIN.AND.ISPIN.EQ.2) FAC=-FAC
                ELSE
                  FAC=1.D0
                END IF
                CALL GRAPHICS_ADDDENMAT(IB,IKPT,ISPIN,FAC &
     &                           ,LMNXX,DENMAT)
              ENDDO
            END IF
            LMXX=9
            ALLOCATE(DRHOL(NR,LMXX))
            CALL GRAPHICS_1CRHO(R1,DEX,NR,LNX,LOX,AEPHI,PSPHI,LMNXX &
     &                   ,DENMAT(1,1),LMXX,DRHOL)
            DEALLOCATE(DENMAT)
          END IF
call timing$clockon('--1')
          CALL GRAPHICS_RHOLTOR(RBAS,fact*nr1,fact*NR2,fact*NR3 &
     &           ,1,fact*NR1,wavebig,BAS,R1,DEX,NR,LMXX,DRHOL)
call timing$clockoff('--1')
!          CALL GRAPHICS_RHOLTOR(RBAS,nr1,NR2,NR3 &
!    &           ,NR1START,NR1L,WAVE,BAS &
!    &           ,R1,DEX,NR,LMXX,DRHOL)
          DEALLOCATE(DRHOL)
          DEALLOCATE(LOX)
          DEALLOCATE(AEPHI)
          DEALLOCATE(PSPHI)
        ENDDO  
!
!       ================================================================
!       ==  PRINT WAVE                                                ==
!       ================================================================
                              CALL TRACE$PASS('PRINT')
        CALL MPE$QUERY(NTASKNUM,NTASKID)
        PRINT*,'PRINTING OUT WAVE ',TITLE(1:50)
!       IF(NTASKID.EQ.1) THEN
!         ALLOCATE(WORK(NR1*NR2*NR3))
!       ELSE
!         ALLOCATE(WORK(1))
!       END IF
!       CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WORK)
        IF(NTASKID.EQ.1) THEN
          STRING=' '
          WRITE(STRING,FMT=*)IPIC
          STRING='WAVEPLOT'//ADJUSTL(STRING)
          PRINT*,'FILE ',STRING,'::',TRIM(FILE)
          CALL FILEHANDLER$SETFILE(STRING,.FALSE.,TRIM(FILE))
          CALL FILEHANDLER$SETSPECIFICATION(STRING,'FORM','UNFORMATTED')
          CALL FILEHANDLER$UNIT(STRING,NFIL)
!         CALL WRITEWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,Q,ATOMNAME &
!   &                       ,NR1,NR2,NR3,WORK)
          CALL WRITEWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,Q,ATOMNAME &
    &                       ,NR1*fact,NR2*fact,NR3*fact,wavebig)
          CALL FILEHANDLER$CLOSE(STRING)
        END IF
 9000   CONTINUE
      ENDDO
      CALL WAVES$SETL4('RAWSTATES',SAVETRAWSTATES)
      DEALLOCATE(WAVE)
      DEALLOCATE(WAVEBIG)
      DEALLOCATE(OCC)
      DEALLOCATE(ATOMNAME)
      DEALLOCATE(Z)
      DEALLOCATE(Q)
      DEALLOCATE(POS)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WRITEWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,Q,NAME,NR1,NR2,NR3,WAVE)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      CHARACTER(*) ,INTENT(IN) :: TITLE
      INTEGER(4)   ,INTENT(IN) :: NAT
      REAL(8)      ,INTENT(IN) :: RBAS(3,3)
      REAL(8)      ,INTENT(IN) :: POS(3,NAT)
      REAL(8)      ,INTENT(IN) :: Z(NAT)
      REAL(8)      ,INTENT(IN) :: Q(NAT)
      CHARACTER(32),INTENT(IN) :: NAME(NAT)
      INTEGER(4)   ,INTENT(IN) :: NR1
      INTEGER(4)   ,INTENT(IN) :: NR2
      INTEGER(4)   ,INTENT(IN) :: NR3
      REAL(8)      ,INTENT(IN) :: WAVE(NR1,NR2,NR3)
!     ******************************************************************
      REWIND NFIL
      WRITE(NFIL)'WAVEPLOT',LEN(TITLE)
      WRITE(NFIL)TITLE
      WRITE(NFIL)RBAS,NAT
      WRITE(NFIL)NR1,NR2,NR3
      WRITE(NFIL)NAME
      WRITE(NFIL)Z
      WRITE(NFIL)POS
      WRITE(NFIL)Q
      WRITE(NFIL)WAVE
      WRITE(NFIL)'END OF FILE'
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_ADDDENMAT(IB,IKPT,ISPIN,FAC,LMNXX,DENMAT)
!     ******************************************************************
!     **                                                              **
!     **  GET PROJECTIONS FROM WAVE AND ADD TO 1C-DENSITY MATRIX      **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: IB
      INTEGER(4),INTENT(IN)   :: IKPT
      INTEGER(4),INTENT(IN)   :: ISPIN
      REAL(8)   ,INTENT(IN)   :: FAC
      INTEGER(4),INTENT(IN)   :: LMNXX
      REAL(8)   ,INTENT(INOUT):: DENMAT(LMNXX,LMNXX)
      REAL(8)   ,ALLOCATABLE  :: PROJ(:) ! (LMNXX)  POINTER ($PROJ,PROJ)
      INTEGER(4)              :: LMN1,LMN2
!     ******************************************************************
      CALL WAVES$SETI4('IB',IB)
      CALL WAVES$SETI4('IKPT',IKPT)
      CALL WAVES$SETI4('ISPIN',ISPIN)
      ALLOCATE(PROJ(LMNXX))
      CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROJ)
      DO LMN1=1,LMNXX
        DO LMN2=1,LMNXX
          DENMAT(LMN1,LMN2)=DENMAT(LMN1,LMN2) &
     &                     +FAC*PROJ(LMN1)*PROJ(LMN2)
        ENDDO
      ENDDO
      DEALLOCATE(PROJ)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_ADDRHO(IB,IKPT,ISPIN,FAC,NNR,WAVE)
!     ******************************************************************
!     **                                                              **
!     **  ADD DENSITY OF A GIVEN PS-WAVE FUNCTION TO WAVE             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NNR       ! NUMBER OF R-SPACE GRID POINTS FOR WAVE
      INTEGER(4),INTENT(IN)   :: IB        ! BAND INDEX
      INTEGER(4),INTENT(IN)   :: IKPT      ! K-POINT INDEX
      INTEGER(4),INTENT(IN)   :: ISPIN     ! SPIN INDEX
      REAL(8)   ,INTENT(IN)   :: FAC       ! WEIGHT OF THIS STATE
      REAL(8)   ,INTENT(INOUT):: WAVE(NNR) ! DENSITY
      REAL(8)   ,ALLOCATABLE  :: PSI(:)    ! (NNR)
      INTEGER(4)              :: IR
!     ******************************************************************
      CALL WAVES$SETI4('IB',IB)
      CALL WAVES$SETI4('IKPT',IKPT)
      CALL WAVES$SETI4('ISPIN',ISPIN)
      ALLOCATE(PSI(NNR))
      CALL WAVES$GETR8A('PSPSI',NNR,PSI)
      DO IR=1,NNR            
        WAVE(IR)=WAVE(IR)+FAC*PSI(IR)**2            
      ENDDO
      DEALLOCATE(PSI)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_RHOLTOR(RBAS,NR1,NR2,NR3,NR1START,NR1L,RHO,R0 &
     &           ,R1,DEX,NRX,LMX,DRHOL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER     :: LMXX=36
      REAL(8)   ,PARAMETER     :: TOL=1.D-5
      INTEGER(4),INTENT(IN)    :: NR1,NR2,NR3,NR1START,NR1L
      INTEGER(4),INTENT(IN)    :: NRX
      INTEGER(4),INTENT(IN)    :: LMX
      REAL(8)   ,INTENT(INOUT) :: RHO(NR1L,NR2,NR3)
      REAL(8)   ,INTENT(IN)    :: DRHOL(NRX,LMX)
      REAL(8)   ,INTENT(IN)    :: RBAS(3,3)
      REAL(8)   ,INTENT(IN)    :: R0(3)
      REAL(8)   ,INTENT(IN)    :: R1
      REAL(8)   ,INTENT(IN)    :: DEX
      REAL(8)                  :: DR(3,3)
      REAL(8)                  :: YLM(LMXX)
      REAL(8)                  :: RVEC(3)
      REAL(8)                  :: XEXP
      REAL(8)                  :: RMAX,RMAX2,RI,SVAR,svar1
      INTEGER(4)               :: IR,LM
      INTEGER(4)               :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      REAL(8)                  :: T1,T2,T3
      REAL(8)                  :: X1,X2,X3
      INTEGER(4)               :: I1,I2,I3,I11,I21,I31,I
      REAL(8)                  :: DIS,DIS2
      REAL(8)                  :: XIR
      REAL(8)                  :: W1,W2
!     ******************************************************************
      XEXP=DEXP(DEX)
      IF(LMXX.LT.LMX) THEN
        CALL ERROR$MSG('INCREASE DIMENSION LMXX')
        CALL ERROR$STOP('GRAPHICS_RHOLTOR')
      END IF
!
!     ==================================================================
!     ==  DETERMINE RMAX                                              ==
!     ==================================================================
      RMAX=0.D0
      RI=R1/XEXP
      DO IR=1,NRX
        RI=RI*XEXP
        SVAR=0.D0
        DO LM=1,LMX
          SVAR=MAX(DABS(DRHOL(IR,LM)),SVAR)
        ENDDO
        IF(SVAR.GT.TOL)RMAX=RI
      ENDDO
      RMAX2=RMAX**2
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      DO I=1,3
        DR(I,1)=RBAS(I,1)/DBLE(NR1)
        DR(I,2)=RBAS(I,2)/DBLE(NR2)
        DR(I,3)=RBAS(I,3)/DBLE(NR3)
      ENDDO
      CALL BOXSPH(DR,R0(1),R0(2),R0(3),RMAX &
     &                 ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
      DO I1=MIN1,MAX1
        T1=DBLE(I1)
        I11=MOD(MOD(I1,NR1)+NR1,NR1)+1
        I11=I11-NR1START+1
        IF(I11.GE.1.AND.I11.LE.NR1L) THEN
          DO I2=MIN2,MAX2
            T2=DBLE(I2)
            I21=MOD(MOD(I2,NR2)+NR2,NR2)+1
            X1=DR(1,1)*T1+DR(1,2)*T2-R0(1)
            X2=DR(2,1)*T1+DR(2,2)*T2-R0(2)
            X3=DR(3,1)*T1+DR(3,2)*T2-R0(3)
            DO I3=MIN3,MAX3
              I31=MOD(MOD(I3,NR3)+NR3,NR3)+1
              T3=DBLE(I3)
              RVEC(1)=X1+DR(1,3)*T3
              RVEC(2)=X2+DR(2,3)*T3
              RVEC(3)=X3+DR(3,3)*T3
              DIS2=RVEC(1)**2+RVEC(2)**2+RVEC(3)**2
              IF(DIS2.LE.RMAX2) THEN
                DIS=MAX(1.D-8,DSQRT(DIS2))
                CALL GETYLM(LMX,RVEC,YLM)
                SVAR=0.D0
                DO LM=1,LMX
                  CALL RADIAL$VALUE(R1,DEX,NRX,DRHOL(1,LM),DIS,SVAR1)
                  SVAR=SVAR+SVAR1*YLM(LM)
                ENDDO
                RHO(I11,I21,I31)=RHO(I11,I21,I31)+SVAR
              END IF
            ENDDO
          ENDDO
        END IF
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_1CRHO(R1,DEX,NRX,LNX,LOX,AEPHI,PSPHI,LMNX &
     &                   ,DENMAT,LMX,DRHOL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: LMXX=36
      REAL(8)   ,INTENT(IN)  :: R1
      REAL(8)   ,INTENT(IN)  :: DEX
      INTEGER(4),INTENT(IN)  :: NRX
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      REAL(8)   ,INTENT(IN)  :: AEPHI(NRX,LNX)
      REAL(8)   ,INTENT(IN)  :: PSPHI(NRX,LNX)
      INTEGER(4),INTENT(IN)  :: LMNX
      REAL(8)   ,INTENT(IN)  :: DENMAT(LMNX,LMNX)
      INTEGER(4),INTENT(IN)  :: LMX
      REAL(8)   ,INTENT(OUT) :: DRHOL(NRX,LMX)
      INTEGER(4)             :: LM,IR,LMN1,LN1,LM1,L1,M1,LMN2,LN2,LM2,L2,M2
      REAL(8)                :: CG,DENMAT1,SVAR
      REAL(8)   ,ALLOCATABLE :: WORK(:)
!     ******************************************************************
!
      DO LM=1,LMX
        DO IR=1,NRX
          DRHOL(IR,LM)=0.D0
        ENDDO
      ENDDO
!
      ALLOCATE(WORK(NRX))
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
!
          DO IR=1,NRX
            WORK(IR)=AEPHI(IR,LN1)*AEPHI(IR,LN2) &
     &              -PSPHI(IR,LN1)*PSPHI(IR,LN2)
          ENDDO
          DO M1=1,2*L1+1
            LM1=L1**2+M1
            DO M2=1,2*L2+1
              LM2=L2**2+M2
              DENMAT1=DENMAT(LMN1+M1,LMN2+M2)
!
              DO LM=1,LMX
                CALL CLEBSCH(LM,LM1,LM2,CG)
                IF(CG.NE.0.D0) THEN
                  SVAR=CG*DENMAT1
                  DO IR=1,NRX
                    DRHOL(IR,LM)=DRHOL(IR,LM)+WORK(IR)*SVAR
                  ENDDO
                END IF
              ENDDO
!
            ENDDO
          ENDDO
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO
      DEALLOCATE(WORK)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_1CWAVE(R1,DEX,NRX,LNX,LOX,AEPHI,PSPHI,LMNX &
     &                   ,PROJ,LMX,DRHOL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: R1
      REAL(8)   ,INTENT(IN)  :: DEX
      INTEGER(4),INTENT(IN)  :: NRX
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: LMX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      REAL(8)   ,INTENT(IN)  :: AEPHI(NRX,LNX)
      REAL(8)   ,INTENT(IN)  :: PSPHI(NRX,LNX)
      REAL(8)   ,INTENT(IN)  :: PROJ(LMNX)
      REAL(8)   ,INTENT(OUT) :: DRHOL(NRX,LMX)
      INTEGER(4)             :: LM,IR,LMN,LN,M,L
      REAL(8)                :: SVAR
!     ******************************************************************
      DO LM=1,LMX
        DO IR=1,NRX
          DRHOL(IR,LM)=0.D0
        ENDDO
      ENDDO
!
      LMN=0
      DO LN=1,LNX
        L=LOX(LN)
        DO M=1,2*L+1
          LMN=LMN+1
          LM=L**2+M
          IF(LM.LE.LMX) THEN
            SVAR=PROJ(LMN)
!            PRINT*,'1CWAVE',LN,M,LMN,LM,SVAR
            DO IR=1,NRX
              DRHOL(IR,LM)=DRHOL(IR,LM)+SVAR*(AEPHI(IR,LN)-PSPHI(IR,LN))
            ENDDO
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_REFINEGRID(NR1,NR2,NR3,FACT,WAVE,WAVEBIG)
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: NR1
      INTEGER(4),INTENT(IN)     :: NR2
      INTEGER(4),INTENT(IN)     :: NR3
      INTEGER(4),INTENT(IN)     :: FACT
      REAL(8),INTENT(IN)        :: WAVE(NR1,NR2,NR3)
      REAL(8),INTENT(OUT)       :: WAVEBIG(FACT*NR1,FACT*NR2,FACT*NR3)
      COMPLEX(8),ALLOCATABLE    :: WORKC1(:,:,:)
      COMPLEX(8),ALLOCATABLE    :: WORKC2(:,:,:)
      INTEGER(4)                :: I
      INTEGER(4)                :: J      
      INTEGER(4)                :: K
      INTEGER(4)                :: E1 
      INTEGER(4)                :: E2   
!     ******************************************************************
      ALLOCATE(WORKC1(NR1,NR2,NR3))
      WORKC1(:,:,:)=CMPLX(wave(:,:,:),KIND=8)
      ALLOCATE(WORKC2(NR1,NR2,NR3)) 
      CALL LIB$3DFFTC8('RTOG',NR1,NR2,NR3,WORKC1,WORKC2)
      DEALLOCATE(WORKC1)
      ALLOCATE(WORKC1(FACT*NR1,FACT*NR2,FACT*NR3))      
      WORKC1=(0.D0,0.D0)
      I=NR1/2  !MIND NRX ARE GENERALLY EVEN?
      J=NR2/2
      K=NR3/2
      E1=2*FACT-1   
      E2=2*FACT 
      WORKC1(1:I      ,1:J      ,1:K)      =WORKC2(1:I  ,1:J  ,1:K)
      WORKC1(E1*I:E2*I,1:J      ,1:K)      =WORKC2(I:2*I,1:J  ,1:K)
      WORKC1(1:I      ,E1*J:E2*J,1:K)      =WORKC2(1:I  ,J:2*J,1:K)        
      WORKC1(1:I      ,1:J      ,E1*K:E2*K)=WORKC2(1:I  ,1:J  ,K:2*K)
      WORKC1(E1*I:E2*I,E1*J:E2*J,1:K)      =WORKC2(I:2*I,J:2*J,1:K)
      WORKC1(1:I      ,E1*J:E2*J,E1*K:E2*K)=WORKC2(1:I  ,J:2*J,K:2*K)        
      WORKC1(E1*I:E2*I,1:J      ,E1*K:E2*K)=WORKC2(I:2*I,1:J  ,K:2*K)
      WORKC1(E1*I:E2*I,E1*J:E2*J,E1*K:E2*K)=WORKC2(I:2*I,J:2*J,K:2*K)      
      DEALLOCATE(WORKC2)
      ALLOCATE(WORKC2(FACT*NR1,FACT*NR2,FACT*NR3))
      CALL LIB$3DFFTC8('GTOR',FACT*NR1,FACT*NR2,FACT*NR3,WORKC1,WORKC2)
      WAVEBIG(:,:,:)=REAL(WORKC2(:,:,:),KIND=8)
      DEALLOCATE(WORKC1)
      DEALLOCATE(WORKC2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_CREATEPOT(IPIC,FACT)
!     ******************************************************************
      USE GRAPHICS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)                 :: NPIC
      INTEGER(4)                 :: IPIC
      INTEGER(4)                 :: FACT
      REAL(8)   ,ALLOCATABLE     :: VHARTREE(:)
      REAL(8)   ,ALLOCATABLE     :: WORK(:)
      REAL(8)   ,ALLOCATABLE     :: POTENTIAL(:)
      INTEGER(4)                 :: NRL
      INTEGER(4)                 :: NR1L
      INTEGER(4)                 :: NR1,NR2,NR3
      CHARACTER(32)              :: STRING
      CHARACTER(256)             :: FILE
      REAL(8)                    :: RBAS(3,3)
      CHARACTER(256)             :: TITLE
      INTEGER(4)                 :: NFIL
      INTEGER(4)                 :: NAT
      INTEGER(4)                 :: IAT
      CHARACTER(32),ALLOCATABLE  :: ATOMNAME(:)
      REAL(8)   ,aLLOCATABLE     :: Q(:)
      REAL(8)   ,ALLOCATABLE     :: Z(:)
      REAL(8)   ,ALLOCATABLE     :: POS(:,:)
      INTEGER(4)                 :: NR
      INTEGER(4)                 :: LMRX
      REAL(8)                    :: DEX
      REAL(8)                    :: R1
      REAL(8)   ,ALLOCATABLE     :: AEPOT(:,:)
      REAL(8)   ,ALLOCATABLE     :: PSPOT(:,:)
      INTEGER(4)                 :: NFILO
      real(8)   ,allocatable     :: onecpot(:,:,:)
      integer(4)                 :: ntasks,thistask
!     ******************************************************************
!COLLECTING OF INFORMATION
      CALL MPE$QUERY(NTASKS,THISTASK)
      CALL LINKEDLIST$SELECT(LL_GRPHCS,'~')
      CALL LINKEDLIST$SELECT(LL_GRPHCS,'IMAGE',IPIC)
      CALL LINKEDLIST$GET(LL_GRPHCS,'TITLE',1,TITLE)
      CALL LINKEDLIST$GET(LL_GRPHCS,'FILE',1,FILE)
!
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$GETI4('NRL',NRL)
      CALL PLANEWAVE$GETI4('NR1L',NR1L)
      CALL PLANEWAVE$GETI4('NR1',NR1)
      CALL PLANEWAVE$GETI4('NR2',NR2)
      CALL PLANEWAVE$GETI4('NR3',NR3)
!
!     ==================================================================
!     ==  expand grid by fact and unparallelize                       ==
!     ==================================================================
      ALLOCATE(VHARTREE(NRL))
      ALLOCATE(WORK(NR1*NR2*NR3))
      CALL PLANEWAVE$SUPFFT('GTOR',1,NGL,PWPOT,NRL,VHARTREE)
      CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,VHARTREE,NR1*NR2*NR3,WORK)
      DEALLOCATE(VHARTREE)
      ALLOCATE(POTENTIAL(FACT*NR1*FACT*NR2*FACT*NR3))
      CALL GRAPHICS_REFINEGRID(NR1,NR2,NR3,FACT,WORK,POTENTIAL)
      DEALLOCATE(WORK)
!
!     ==================================================================
!     ==  one-center contributions                                    ==
!     ==================================================================
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ATOMNAME(NAT))
      ALLOCATE(Z(NAT))
      ALLOCATE(Q(NAT))
      ALLOCATE(POS(3,NAT))
      IF(.NOT.ALLOCATED(AE1CPOT)) THEN
        CALL SETUP$GETI4('NRX',NRX)
        CALL SETUP$GETI4('LMRXX',LMRXX)
        CALL ATOMLIST$NATOM(NAT)
        ALLOCATE(AE1CPOT(NRX,LMRXX,NAT))
        ALLOCATE(PS1CPOT(NRX,LMRXX,NAT))
        AE1CPOT=0.D0
        PS1CPOT=0.D0
      END IF
      ALLOCATE(ONECPOT(NRX,LMRXX,NAT))
      onecPOT=AE1CPOT-PS1CPOT
      CALL MPE$COMBINE('+',onecPOT)
      DO IAT=1,NAT
        CALL SETUP$RADGRID(IAT,R1,DEX,NR)
        PRINT*,'CLEMENS: LMRX: ',LMRXx
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME(IAT))
        CALL ATOMLIST$GETR8A('R(0)',IAT,3,POS(:,IAT))
        CALL ATOMLIST$GETR8('Z',IAT,Z(IAT))
        CALL ATOMLIST$GETR8('Q',IAT,Q(IAT))
PRINT*,'INCLUDE AE-CONTRIBUTIONS'
call timing$clockon('graphics 1cpotential')
         CALL GRAPHICS_RHOLTOR(RBAS,FACT*NR1,FACT*NR2,FACT*NR3,1,FACT*NR1 &
      &           ,POTENTIAL,POS(:,IAT),R1,DEX,NRx,LMRXx,oneCPOT(:,:,IAT))
call timing$clockoff('graphics 1cpotential')
PRINT*,'INCLUDEd AE-CONTRIBUTIONS'
      ENDDO
      deallocate(onecpot)
!
!     ==================================================================
!     ==  write to file                                               ==
!     ==================================================================
      if(thistask.eq.1) then
        STRING=' '
        WRITE(STRING,FMT=*)IPIC
        STRING='WAVEPLOT'//ADJUSTL(STRING)
        PRINT*,'FILE ',STRING,'::',TRIM(FILE)
        CALL FILEHANDLER$SETFILE(STRING,.FALSE.,TRIM(FILE))
        CALL FILEHANDLER$SETSPECIFICATION(STRING,'FORM','UNFORMATTED')
        CALL FILEHANDLER$UNIT(STRING,NFIL)
        CALL WRITEWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,Q,ATOMNAME &
     &                  ,fACT*NR1,FACT*NR2,FACT*NR3,POTENTIAL)    
        CALL FILEHANDLER$CLOSE(STRING)
      end if
      DEALLOCATE(POTENTIAL)
      RETURN
      END
