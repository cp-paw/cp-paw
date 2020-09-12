!NICHTKOLLINEAR
!KLEINE KOMPONENTE
!GAMMMA PUNKT NICHT MEHR ERSTES ELEMENT AUF DEM G-GITTER
!TEST EXTRAPOLATE
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE HYPERFINE_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: HYPERFINE_MODULE                                                   **
!**                                                                           **
!**                                                                           **
!**  EVALUATES HYPER FINE PARAMETERS SUCH AS                                  **
!**    ELECTRIC FIELD GRADIENTS                                               **
!**    ISOMER SHIFT                                                           **
!**    FERMI CONTACT TERM                                                     **
!**    ANISOTROPIC HYPER FINE PARAMETER                                       **
!**                                                                           **
!**  STRUCTURE OF THE LINKED LIST                                             **
!**  ============================                                             **
!**  ~                                                                        **
!**   ATOMNAME                                                                **
!**     TEFG   SWITCH FOR ELECTRIC FIELD GRADIENT CALCULATION                 **
!**     TIS    SWITCH FOR ISOMER SHIFT CALCULATION                            **
!**     TFC    SWITCH FOR FERMI CONTACT TERM CALCULATION                      **
!**     TANIS  SWITCH FOR ANISOTROPIC HYPERFINE PARAMETER                     **
!**     AE2VLM                                                                **
!**     PS2VLM                                                                **
!**     PSV                                                                   **
!**     PSRHOTPW                                                              **
!**     AERHOTLM                                                              **
!**     PSRHOTLM                                                              **
!**     PSRHOSPW                                                              **
!**     AERHOSLM                                                              **
!**     PSRHOSLM                                                              **
!**     PSRHOS2PW                                                             **
!**     AERHOS2LM                                                             **
!**     PSRHOS2LM                                                             **
!**                                                                           **
!**  REMARKS                                                                  **
!**     ANISOTROPIC HYPERFINE PARAMETER STILL INCORRECT                       **
!**     ISOMER SHIFT MAY REQUIRE CORE POLARIZATION                            **
!**     EFG NOT YET CONSISTENT PW AND PS-1                                    **
!**     EFG SECOND DERIVATIVES OR PREFACTOR OF THE QUADRATIC TERM?            **
!**                                                                           **
!*******************************************************************************
USE LINKEDLIST_MODULE
LOGICAL(4)   :: TON=.FALSE.    ! ON/OFF SWITCH
LOGICAL(4)   :: TWAKE=.TRUE.   ! ON/OFF SWITCH TEMPORARY
TYPE(LL_TYPE):: LL_HPRFN          ! LINKEDLIST POINTER
LOGICAL(4)   :: TINI=.FALSE.
LOGICAL(4)   :: TPW=.TRUE.     ! ON/OFF SWITCH PLANE WAVE PART
LOGICAL(4)   :: TSPIN=.FALSE.  ! ON/OFF SWITCH SPINPOLARIZED
LOGICAL(4)   :: TPWRHOS=.FALSE.! ON/OFF SWITCH SPINPOLARIZED
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE_NEWLIST
      IMPLICIT NONE
      IF(TINI) RETURN
      CALL LINKEDLIST$NEW(LL_HPRFN)
      TINI=.TRUE.
      RETURN
      END SUBROUTINE HYPERFINE_NEWLIST
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTRAPOLATE(NP,RI,FI_,R0,F0)
!     **************************************************************************
!     **  POLYNOMIAL EXTRAPOLATION OF ORDER NP FROM NP POINTS (RI,FI_)        **
!     **  TO THE POINT (R0,F0)                                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP
      REAL(8)   ,INTENT(IN) :: RI(NP)
      REAL(8)   ,INTENT(IN) :: FI_(NP)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: FI(NP)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J,IP
!     **************************************************************************
      FI(:)=FI_(:)
      F0=0.D0
      DO I=1,NP
        SVAR=1.D0
        DO J=1,I-1
          SVAR=SVAR*(R0-RI(J))/(RI(I)-RI(J))
        ENDDO
        F0=F0+SVAR*FI(I)
        DO IP=I+1,NP
          SVAR=1.D0
          DO J=1,I-1
            SVAR=SVAR*(RI(IP)-RI(J))/(RI(I)-RI(J))
          ENDDO
          FI(IP)=FI(IP)-SVAR*FI(I)
        ENDDO
        FI(I)=0.D0
      ENDDO
      RETURN
      END SUBROUTINE EXTRAPOLATE
END MODULE HYPERFINE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$SELECT(ID)
!     ********************************************************************
!     **  SELECT AN ATOM IN ORDER TO INPUT INFORMATION                  **
!     **                                                                **
!     **                                                                **
!     ********************************************************************
      USE HYPERFINE_MODULE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)              :: TCHK
!     ********************************************************************
      CALL HYPERFINE_NEWLIST
      CALL LINKEDLIST$SELECT(LL_HPRFN,ID)
      IF(ID.EQ.'..'.OR.ID.EQ.'~') RETURN
      CALL LINKEDLIST$EXISTD(LL_HPRFN,'TEFG',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_HPRFN,'TEFG',0,.FALSE.)
      CALL LINKEDLIST$EXISTD(LL_HPRFN,'TIS',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_HPRFN,'TIS',0,.FALSE.)
      CALL LINKEDLIST$EXISTD(LL_HPRFN,'TANIS',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_HPRFN,'TANIS',0,.FALSE.)
      CALL LINKEDLIST$EXISTD(LL_HPRFN,'TFC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_HPRFN,'TFC',0,.FALSE.)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$SETL4(ID,VAL)
!     ********************************************************************
!     **                                                                **
!     **  STORE CONTROLE DATA ON THE LIST                               **
!     **                                                                **
!     ********************************************************************
      USE HYPERFINE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ********************************************************************
      CALL HYPERFINE_NEWLIST
      IF(ID.EQ.'WAKE')THEN
        TWAKE=VAL
        TWAKE=TWAKE.AND.TON
        RETURN
      END IF
      IF((ID.NE.'TEFG').AND.(ID.NE.'TANIS') &
     &                 .AND.(ID.NE.'TIS').AND.(ID.NE.'TFC')) THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('HYPERFINE$SET')
      END IF
      CALL LINKEDLIST$SET(LL_HPRFN,ID,0,VAL)
      TON=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$GETFLAG(IDENT1_,TCHK_)
!     ********************************************************************
!     **                                                                **
!     **  RETURNS TCHK_=.TRUE. IF THE SPIN DENSITY IN G-SPACE IS        **
!     **  IS REQUIRED FOR HYPERFINE PARAMETERS                          **
!     **                                                                **
!     ********************************************************************
      USE HYPERFINE_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN)  :: IDENT1_  !CAN BE 'SPIN' OR 'RHO'
      LOGICAL(4)   ,INTENT(OUT) :: TCHK_
      INTEGER(4)                :: NAT
      INTEGER(4)                :: IAT
      CHARACTER(32)             :: ATOMNAME
!     ********************************************************************
      CALL HYPERFINE_NEWLIST
      CALL ATOMLIST$NATOM(NAT)
      TCHK_=.FALSE.
      IF(.NOT.TWAKE) RETURN
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
        CALL LINKEDLIST$EXISTL(LL_HPRFN,ATOMNAME,1,TCHK_)
        IF(TCHK_) THEN
          CALL LINKEDLIST$SELECT(LL_HPRFN,ATOMNAME)
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'TANIS',1,TCHK_)
          IF(TCHK_) THEN
            CALL LINKEDLIST$GET(LL_HPRFN,'TANIS',1,TCHK_)
            IF(TCHK_.AND.IDENT1_.EQ.'SPIN') RETURN
          END IF
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'TFC',1,TCHK_)
          IF(TCHK_) THEN
            CALL LINKEDLIST$GET(LL_HPRFN,'TFC',1,TCHK_)
            IF(TCHK_.AND.IDENT1_.EQ.'SPIN') RETURN
          END IF
          CALL LINKEDLIST$SELECT(LL_HPRFN,'..')
        END IF
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$SET1CPOT(IDENT_,IAT_,GID,NR,NRX,LMRX,POT)
!     ********************************************************************
!     **  USE 1-CENTER POTENTIAL FOR ELECTRIC FIELD GRADIENTS          **
!     ********************************************************************
      USE HYPERFINE_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: IDENT_  ! CAN BE 'AE' OR 'PS' 
      INTEGER(4)   ,INTENT(IN) :: IAT_    ! ATOM INDEX (SEE ATOMLIST)
      INTEGER(4)   ,INTENT(IN) :: GID
      INTEGER(4)   ,INTENT(IN) :: NR,NRX
      INTEGER(4)   ,INTENT(IN) :: LMRX
      REAL(8)      ,INTENT(IN) :: POT(NRX,LMRX)
      LOGICAL(4)               :: TEFG
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: V2(5)
      INTEGER(4)   ,PARAMETER  :: NP=5
      REAL(8)                  :: YARRAY(NP)
      REAL(8)                  :: R(NR)
      INTEGER(4)               :: LM
      CHARACTER(32)            :: ATOMNAME
!     ********************************************************************
! 
!     ===================================================================
!     == CHECK WHETHER ACTION IS REQUIRED                              ==
!     ===================================================================
      IF(.NOT.TWAKE) RETURN
                                CALL TRACE$PUSH('HYPERFINE$SET1CPOT')
      CALL ATOMLIST$GETCH('NAME',IAT_,ATOMNAME)
      CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
      CALL LINKEDLIST$EXISTL(LL_HPRFN,ATOMNAME,1,TCHK)
      IF(.NOT.TCHK) THEN 
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
      CALL LINKEDLIST$SELECT(LL_HPRFN,ATOMNAME)
      CALL LINKEDLIST$GET(LL_HPRFN,'TEFG',1,TEFG)
      IF(.NOT.TEFG) THEN; CALL TRACE$POP; RETURN; END IF 
! 
!     ===================================================================
!     == EVALUATE ELECTRIC FIELD GRADIENTS                             ==
!     == V(R)=R(I)*V(I,J)*R(J)                                         ==
!     == V2LM=LIM(R->0) V_L(R)/R**2                                    ==
!     ===================================================================
      CALL RADIAL$R(GID,NR,R)
      IF(IDENT_.EQ.'PS') WRITE(*,FMT='("PS1X ",5E20.10)')R(1:NP)
      IF(IDENT_.EQ.'PS') WRITE(*,FMT='("PS1V ",1E20.10)')POT(1,1)
      IF(LMRX.GE.9) THEN
        DO LM=5,9
          YARRAY(1:NP)=POT(1:NP,LM)/R(1:NP)**2
          IF(IDENT_.EQ.'PS') WRITE(*,FMT='("PS1V2 ",5E20.10)')YARRAY
!         CALL EXTRAPOLATE(NP,R(1:NP),YARRAY,0.D0,V2(LM-4))
          V2(LM-4)=YARRAY(1)
        ENDDO
      ELSE
        V2(:)=0.D0
      END IF
      IF(IDENT_.EQ.'AE') THEN
        CALL LINKEDLIST$SET(LL_HPRFN,'AEV2LM',0,V2)
      ELSE IF(IDENT_.EQ.'PS') THEN
        CALL LINKEDLIST$SET(LL_HPRFN,'PSV2LM',0,V2)
      ELSE
        CALL ERROR$MSG('IDENT CAN ONLY BE AE OR PS')
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$STOP('HYPERFINE$SET1CPOT') 
      END IF  
      CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$SET1CRHO(IDENT_,IDENT1_,IAT_,GID,NR,NRX,LMRX,RHO)
!     ********************************************************************
!     **  GET SECOND DERIVATIVE OF THE RADIAL POTENTIAL AT THE ORIGIN   **
!     ********************************************************************
      USE HYPERFINE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_  ! CAN BE 'AE' OR 'PS' 
      CHARACTER(*),INTENT(IN) :: IDENT1_ ! CAN BE 'TOT' OR 'SPIN' 
      INTEGER(4)  ,INTENT(IN) :: IAT_    ! ATOM INDEX (SEE ATOMLIST)
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR,NRX
      INTEGER(4)  ,INTENT(IN) :: LMRX
      REAL(8)     ,INTENT(IN) :: RHO(NRX,LMRX)
      LOGICAL(4)              :: TIS     ! SWITCH FOR ISOMER SHIFT
      LOGICAL(4)              :: TFC     ! SWITCH FOR FERMI CONTACT TERM
      LOGICAL(4)              :: TANIS   ! SWITCH FOR ANISOTROPIC
      LOGICAL(4)              :: TEFG    ! SWITCH FOR EFG
      LOGICAL(4)              :: TCHK
      REAL(8)                 :: RHO0
      INTEGER(4)  ,PARAMETER  :: NP=4
      REAL(8)                 :: YARRAY(NP)
      REAL(8)                 :: WORK(NR)
      REAL(8)                 :: ANIS(5)
      REAL(8)                 :: V2(5)
      INTEGER(4)              :: LM
      CHARACTER(32)           :: ATOMNAME
      REAL(8)                 :: Z
      REAL(8)                 :: R(NR)
!     ********************************************************************
! 
!     ===================================================================
!     == CHECK WHETHER ACTION IS REQUIRED                              ==
!     ===================================================================
      IF(.NOT.TWAKE) RETURN
                                CALL TRACE$PUSH('HYPERFINE$SET1CRHO')
      CALL ATOMLIST$GETCH('NAME',IAT_,ATOMNAME)
      CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
      CALL LINKEDLIST$EXISTL(LL_HPRFN,ATOMNAME,1,TCHK)
      IF(.NOT.TCHK) THEN 
        CALL TRACE$POP 
        RETURN 
      END IF
      CALL LINKEDLIST$SELECT(LL_HPRFN,ATOMNAME)
      IF(IDENT1_.NE.'TOT'.AND.IDENT1_.NE.'SPIN') THEN
        CALL ERROR$MSG('IDENT1_ MUST BE EITHER "TOT" OR "SPIN"')
        CALL ERROR$CHVAL('IDENT1_',IDENT1_)
        CALL ERROR$STOP('HYPERFINE$SET1CRHO')
      END IF
      IF(IDENT_.NE.'AE'.AND.IDENT_.NE.'PS') THEN
        CALL ERROR$MSG('IDENT_ MUST BE EITHER "AE" OR "PS"')
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$STOP('HYPERFINE$SET1CRHO')
      END IF

      TIS=.FALSE.
      TFC=.FALSE.
      TANIS=.FALSE.
      TEFG=.FALSE.
      IF(IDENT1_.EQ.'TOT') THEN
        CALL LINKEDLIST$EXISTD(LL_HPRFN,'TIS',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_HPRFN,'TIS',1,TIS)
        END IF
        CALL LINKEDLIST$EXISTD(LL_HPRFN,'TEFG',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_HPRFN,'TEFG',1,TEFG)
        END IF
      ELSE IF(IDENT1_.EQ.'SPIN') THEN
        TSPIN=.TRUE.
        CALL LINKEDLIST$EXISTD(LL_HPRFN,'TFC',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_HPRFN,'TFC',1,TFC)
        END IF
        CALL LINKEDLIST$EXISTD(LL_HPRFN,'TANIS',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_HPRFN,'TANIS',1,TANIS)
        END IF
      END IF
      IF(.NOT.(TIS.OR.TFC.OR.TANIS)) THEN
        CALL TRACE$POP 
        RETURN 
      END IF
      IF(LMRX.LT.1) THEN; CALL TRACE$POP; RETURN; END IF
      CALL RADIAL$R(GID,NR,R)
! 
!     ===================================================================
!     == ISOMER SHIFT                                                  ==
!     ===================================================================
      IF(TIS) THEN
        CALL EXTRAPOLATE(NP,R(1:NP),RHO(1:NP,1),0.D0,RHO0)
        IF(IDENT_.EQ.'AE') THEN
          CALL LINKEDLIST$SET(LL_HPRFN,'AERHOTLM',0,RHO0)
        ELSE IF(IDENT_.EQ.'PS') THEN
          CALL LINKEDLIST$SET(LL_HPRFN,'PSRHOTLM',0,RHO0)
        END IF  
      END IF 
! 
!     ===================================================================
!     == FERMI CONTACT                                                 ==
!     ===================================================================
      IF(TFC) THEN
        IF(IDENT_.EQ.'AE') THEN
          CALL ATOMLIST$GETR8('Z',IAT_,Z)
          CALL EFG_THOMSON(GID,NR,RHO,Z,RHO0)
          CALL LINKEDLIST$SET(LL_HPRFN,'AERHOSLM',0,RHO0)
        ELSE IF(IDENT_.EQ.'PS') THEN
          CALL EXTRAPOLATE(NP,R(1:NP),RHO(1:NP,1),0.D0,RHO0)
          CALL LINKEDLIST$SET(LL_HPRFN,'PSRHOSLM',0,RHO0)
        END IF  
      END IF
! 
!     ===================================================================
!     == ANISOTROPIC HYPERFINE PARAMETER                               ==
!     ===================================================================
      IF(TANIS) THEN
        IF(LMRX.GE.9) THEN
          DO LM=5,9
            CALL RADIAL$POISSON(GID,NR,2,RHO(1,LM),WORK)
            YARRAY(1:NP)=WORK(1:NP)/R(1:NP)**2
            IF(IDENT_.EQ.'PS') WRITE(*,FMT='("PS1VS2 ",5E20.10)')YARRAY
            CALL EXTRAPOLATE(NP,R(1:NP),YARRAY,0.D0,ANIS(LM-4))
          ENDDO
        ELSE
          ANIS(:)=0.D0
        END IF
        IF(IDENT_.EQ.'AE') THEN
          CALL LINKEDLIST$SET(LL_HPRFN,'AERHOS2LM',0,ANIS)
        ELSE IF(IDENT_.EQ.'PS') THEN
          CALL LINKEDLIST$SET(LL_HPRFN,'PSRHOS2LM',0,ANIS)
        END IF  
      END IF
! 
!     ===================================================================
!     == ELECTRIC FIELD GRADIENT                                       ==
!     == EFG CURRENTLY CALCULATED FROM 1CPOTENTIAL
!     ===================================================================
      IF(TEFG) THEN
        IF(LMRX.GE.9) THEN
          DO LM=5,9
            CALL RADIAL$POISSON(GID,NR,2,RHO(1,LM),WORK)
            YARRAY(1:NP)=WORK(1:NP)/R(1:NP)**2
            IF(IDENT_.EQ.'PS') WRITE(*,FMT='("PS1VS2 ",5E20.10)')YARRAY
            CALL EXTRAPOLATE(NP,R(1:NP),YARRAY,0.D0,V2(LM-4))
          ENDDO
        ELSE
          V2(:)=0.D0
        END IF
!CURRENTLY CALCULATED FROM POTENTIAL (SEE SET1CPOT) (CHECK FACTOR 2!)
!        IF(IDENT_.EQ.'AE') THEN
!          CALL LINKEDLIST$SET(LL_HPRFN,'AEV2LM',0,V2)
!        ELSE IF(IDENT_.EQ.'PS') THEN
!          CALL LINKEDLIST$SET(LL_HPRFN,'PSV2LM',0,V2)
!        END IF  
      END IF
!
      CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$SETPWPOT(NG,G,VOFG)
!     ********************************************************************
!     **  GET SECOND DERIVATIVE OF THE RADIAL POTENTIAL AT THE ORIGIN   **
!     ********************************************************************
      USE HYPERFINE_MODULE
      IMPLICIT NONE
      COMPLEX(8) ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4) ,INTENT(IN) :: NG            ! #(G-VECTORS FOR POT)
      REAL(8)    ,INTENT(IN) :: G(3,NG)       ! G-VECTORS FOR POT
      COMPLEX(8) ,INTENT(IN) :: VOFG(NG)    ! POTENTIAL PLANE WAVE COEFFICIENTS
      CHARACTER(32)          :: ATOMNAME      
      INTEGER(4)             :: IAT,I,J,IG    ! DO LOOP INDICES
      INTEGER(4)             :: NAT           ! NUMBER OF ATOMS
      REAL(8)                :: R(3)          ! ATOM POSITION
      REAL(8)                :: V(3,3)        ! D2V/DR2 AT THE NUCLEUS
      REAL(8)                :: TR            ! TRACE OF SECOND DERIVATIVES
      REAL(8)                :: GR            ! G*R
      LOGICAL(4)             :: TEFG          ! SWITCH FOR EFG CALCULATION
      REAL(8)                :: SVAR1,SVAR2,GLEN,GMAX
      LOGICAL(4)             :: TCHK
REAL(8)   ,ALLOCATABLE :: WORK(:,:)
LOGICAL(4),PARAMETER :: TTEST=.FALSE.
!     ********************************************************************
      IF(.NOT.(TWAKE.AND.TPW)) RETURN
                                CALL TRACE$PUSH('HYPERFINE$SETPWPOT')
!     == GET PLANE WAVE CUTOFF AN CONVERT TO MAX. G-VECTOR LENGTH
      CALL POTENTIAL$GETR8('EPWRHO',GMAX)
      GMAX=SQRT(2.D0*GMAX)
!
      CALL ATOMLIST$NATOM(NAT)
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
        CALL LINKEDLIST$EXISTL(LL_HPRFN,ATOMNAME,1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$SELECT(LL_HPRFN,ATOMNAME)
          CALL ATOMLIST$GETR8A('R(0)',IAT,3,R)
          CALL LINKEDLIST$GET(LL_HPRFN,'TEFG',1,TEFG)
          IF(TEFG) THEN
!           ============================================================
!           == EVALUATE ELECTRIC FIELD GRADIENTS                      ==
!           == V(I,J) = D2V(R)/(DI*DJ)                                ==
!           ============================================================
IF(TTEST) THEN 
  ALLOCATE(WORK(6,NG))
  WORK=0.D0
END IF
            DO I=1,3
              DO J=1,3
                V(I,J)=0.D0
              ENDDO
            ENDDO
            DO IG=1,NG
              GR=R(1)*G(1,IG)+R(2)*G(2,IG)+R(3)*G(3,IG)
!             __-1 FROM DERIVATIVE = I**2
!             __ 2 FROM ADDING COMPLEX CONJUGATE
              SVAR1=0.D0
              GLEN=SQRT(G(1,IG)**2+G(2,IG)**2+G(3,IG)**2)
!             == SVAR DAMPS OUT HIGH FREQUENCY OSCILLATIONS OF THE POTENTIAL 
              SVAR1=1.D0/(1.D0+EXP(GLEN+3.D0-GMAX))  !MAY NOT BE A GOOD CHOICE
IF(TTEST) SVAR1=1.D0
              SVAR1=-2.D0*REAL(EXP(+CI*GR)*VOFG(IG),KIND=8)*SVAR1
              SVAR2=G(1,IG)*SVAR1
              V(1,1)=V(1,1)+G(1,IG)*SVAR2
              V(2,1)=V(2,1)+G(2,IG)*SVAR2
              V(3,1)=V(3,1)+G(3,IG)*SVAR2
IF(TTEST) THEN
  WORK(1,IG)=G(1,IG)*SVAR2
  WORK(2,IG)=G(2,IG)*SVAR2
  WORK(3,IG)=G(3,IG)*SVAR2
END IF
              SVAR2=G(2,IG)*SVAR1
              V(2,2)=V(2,2)+G(2,IG)*SVAR2        
              V(3,2)=V(3,2)+G(3,IG)*SVAR2        
IF(TTEST) THEN
  WORK(4,IG)=G(2,IG)*SVAR2
  WORK(5,IG)=G(3,IG)*SVAR2
END IF
              SVAR2=G(3,IG)*SVAR1
              V(3,3)=V(3,3)+G(3,IG)*SVAR2        

IF(TTEST) THEN
  WORK(6,IG)=G(3,IG)*SVAR2
END IF
            ENDDO
IF(TTEST) THEN
!THIS IS TO TEST THE PLANE WAVE CONVERGENCE OF THE EFG AND THE FILTER
 DO IG=1,NG
    TR=WORK(1,IG)+WORK(4,IG)+WORK(6,IG)
    WORK(1,IG)=WORK(1,IG)-TR/3.D0
    WORK(4,IG)=WORK(4,IG)-TR/3.D0
    WORK(6,IG)=WORK(6,IG)-TR/3.D0
  ENDDO
  DO IG=2,NG
    WORK(:,IG)=WORK(:,IG)+WORK(:,IG-1)
  ENDDO
  SVAR2=0.D0
  PRINT*,'HERE COMES THE PLANE WAVE PART OF THE EFG'
  DO IG=1,NG-1
    SVAR1=G(1,IG)**2+G(2,IG)**2+G(3,IG)**2
    SVAR1=SQRT(SVAR1)
    IF(SVAR1.GT.SVAR2) WRITE(6,FMT='(7F20.8)')SVAR1,WORK(:,IG)
    IF(SVAR1.LT.SVAR2) PRINT*,'WARNING'
    SVAR2=SVAR1
  ENDDO
  DEALLOCATE(WORK)
  CALL ERROR$STOP('EFG$SETPWPOT: FORCE STOP')
END IF
!           == SYMMETRIZE ================================================
            V(1,2)=V(2,1)
            V(1,3)=V(3,1)
            V(2,3)=V(3,2)
!           __SUBTRACT TRACE______________________________________________
            TR=V(1,1)+V(2,2)+V(3,3)
            V(1,1)=V(1,1)-TR/3.D0
            V(2,2)=V(2,2)-TR/3.D0
            V(3,3)=V(3,3)-TR/3.D0
            CALL LINKEDLIST$SET(LL_HPRFN,'PSV',0,V)
          END IF
        END IF
      ENDDO
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$SETPWRHO(IDENT_,NG,G,RHOG)
!     ********************************************************************
!     **  GET SECOND DERIVATIVE OF THE RADIAL POTENTIAL AT THE ORIGIN   **
!     ********************************************************************
      USE HYPERFINE_MODULE
      IMPLICIT NONE
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      CHARACTER(*),INTENT(IN) :: IDENT_  ! CAN BE 'TOT' OR 'SPIN'
      INTEGER(4)  ,INTENT(IN) :: NG
      REAL(8)     ,INTENT(IN) :: G(3,NG)
      COMPLEX(8)  ,INTENT(IN) :: RHOG(NG)
      CHARACTER(32)           :: ATOMNAME
      INTEGER(4)              :: IAT,I,J,IG
      REAL(8)                 :: R(3)
      REAL(8)                 :: RHO0
      REAL(8)                 :: GR 
      REAL(8)                 :: TR
      REAL(8)                 :: G2
      REAL(8)                 :: PI
      REAL(8)                 :: FAC
      REAL(8)                 :: V(3,3)
      LOGICAL(4)              :: TIS   ! ON/OFF SWITCH ISOMER SHIFT 
      LOGICAL(4)              :: TFC   ! ON/OFF SWITCH FERMI CONTACT
      LOGICAL(4)              :: TANIS ! ON/OFF SWITCH ANISTROPIC HYPERF.
      LOGICAL(4)              :: TCHK
      REAL(8)                 :: SVAR1,SVAR2
      INTEGER(4)              :: NFILO
      INTEGER(4)              :: NAT
!     ********************************************************************
      IF(.NOT.(TWAKE.AND. TPW)) RETURN
                                CALL TRACE$PUSH('HYPERFINE$SETPWRHO')
      PI=4.D0*ATAN(1.D0)
      IF(IDENT_.NE.'TOT'.AND.IDENT_.NE.'SPIN') THEN
        CALL ERROR$MSG('IDENT_ MUST BE EITHER TOT OR SPIN')
        CALL ERROR$CHVAL('IDENT_',IDENT_)
        CALL ERROR$STOP('HYPERFINE$SETPWRHO')
      END IF
      IF(IDENT_.EQ.'SPIN') TSPIN=.TRUE.
      CALL ATOMLIST$NATOM(NAT)
      CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
!     CALL FILEHANDLER$UNIT('PROT',NFILO)
!     CALL LINKEDLIST$REPORT(LL_HPRFN,NFILO)
      DO IAT=1,NAT
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
        CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
        CALL LINKEDLIST$EXISTL(LL_HPRFN,ATOMNAME,1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$SELECT(LL_HPRFN,ATOMNAME)
          CALL LINKEDLIST$GET(LL_HPRFN,'TIS',1,TIS)
          CALL LINKEDLIST$GET(LL_HPRFN,'TFC',1,TFC)
          IF((TIS.AND.IDENT_.EQ.'TOT').OR.(TFC.AND.IDENT_.EQ.'SPIN')) THEN
!           ============================================================
!           == EVALUATE VALUE OF THE PS-DENSITY AT THE NCULEAR SITE   ==
!           == FOR  ISOMER SHIFT OR FERMI CONTACT INTERACTION         ==
!           ============================================================
            CALL ATOMLIST$GETR8A('R(0)',IAT,3,R)
            RHO0=0.D0
            DO IG=1,NG
              GR=R(1)*G(1,IG)+R(2)*G(2,IG)+R(3)*G(3,IG)
!             __ 2 FROM ADDING COMPLEX CONJUGATE
              RHO0=RHO0+2.D0*REAL(EXP(+CI*GR)*RHOG(IG),KIND=8)
            ENDDO
!           __ SUBTRACT DOUBLE COUNTING OF THE GAMMA POINT____
!ERROR! ASSUMPTION THAT GAMMA POINT RESIDES ON FIRST ELEMENT INCORRECT
!REMARK: USE PARAMETER NGAMMA IF POSSIBLE
            IF(G(1,1)**2+G(2,1)**2+G(3,1)**2.LT.1.D-6) THEN
              RHO0=RHO0-REAL(RHOG(1),KIND=8)
            END IF
            IF(IDENT_.EQ.'SPIN') THEN             
              CALL LINKEDLIST$SET(LL_HPRFN,'PSRHOSPW',0,RHO0)
            ELSE 
              CALL LINKEDLIST$SET(LL_HPRFN,'PSRHOTPW',0,RHO0)
            END IF
          END IF
          CALL LINKEDLIST$GET(LL_HPRFN,'TANIS',1,TANIS)
          IF(TANIS.AND.(IDENT_.EQ.'SPIN')) THEN
            CALL ATOMLIST$GETR8A('R(0)',IAT,3,R)
!           ============================================================
!           == EVALUATE THE MAGNETIC FIELD GRADIENT                   ==
!           == V(R)=R(I)*V(I,J)*R(J)                                  ==
!           == WHERE V(R)=RHOS(RPRIME)/|R-RPRIME|                     ==
!           ============================================================
            DO I=1,3
              DO J=1,3
                V(I,J)=0.D0
              ENDDO
            ENDDO
            FAC=-2.D0*4.D0*PI
!ERROR! ASSUMES GAMMA APOINT RESIDES ON FIRST ELEMENT
            DO IG=2,NG
              GR=R(1)*G(1,IG)+R(2)*G(2,IG)+R(3)*G(3,IG)
              G2=G(1,IG)**2+G(2,IG)**2+G(3,IG)**2
!             __ 2 FROM ADDING COMPLEX CONJUGATE
              SVAR1=FAC*REAL(EXP(+CI*GR)*RHOG(IG),KIND=8)/G2
              SVAR2=G(1,IG)*SVAR1
              V(1,1)=V(1,1)+G(1,IG)*SVAR2
              V(2,1)=V(2,1)+G(2,IG)*SVAR2
              V(3,1)=V(3,1)+G(3,IG)*SVAR2
              SVAR2=G(2,IG)*SVAR1
              V(2,2)=V(2,2)+G(2,IG)*SVAR2        
              V(3,2)=V(3,2)+G(3,IG)*SVAR2        
              SVAR2=G(3,IG)*SVAR1
              V(3,3)=V(3,3)+G(3,IG)*SVAR2        
            ENDDO
!           == SYMMETRIZE ================================================
            V(1,2)=V(2,1)
            V(1,3)=V(3,1)
            V(2,3)=V(3,2)
!           __SUBTRACT TRACE______________________________________________
            TR=V(1,1)+V(2,2)+V(3,3)
            V(1,1)=V(1,1)-TR/3.D0
            V(2,2)=V(2,2)-TR/3.D0
            V(3,3)=V(3,3)-TR/3.D0
            CALL LINKEDLIST$SET(LL_HPRFN,'PSRHOS2PW',0,V)
            CALL TRACE$PASS('HYPERFINE$SETPWRHO PW ANISOTROPIC HYPERFINE PARM')
          END IF
        END IF
      ENDDO
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HYPERFINE$PRINT
!     **************************************************************************
!     **  CALCULATE ELECTRIC FIELD GRADIENTS AND PRINT RESULTS ON FILE        **
!     **************************************************************************
      USE HYPERFINE_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4)           :: TPLANEWAVE
      CHARACTER(32)        :: ATOMNAME
      INTEGER(4)           :: NAT
      INTEGER(4)           :: NTASKID,NTASKNUM
      INTEGER(4)           :: I,J,IAT
      INTEGER(4)           :: NFILO
      LOGICAL(4)           :: TCHK
      LOGICAL(4)           :: TEFG
      LOGICAL(4)           :: TIS
      LOGICAL(4)           :: TFC
      LOGICAL(4)           :: TANIS
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
      REAL(8)              :: EFG(3,3),PSEFG(3,3),AEEFG1(3,3),PSEFG1(3,3)
      REAL(8)              :: ANIS(3,3)
      REAL(8)              :: RHOT0
      REAL(8)              :: RHOS0,PSRHOS0,AERHO1S0,PSRHO1S0
      REAL(8)              :: SVAR,SVAR1
      REAL(8)              :: VXYZ(3,3)
      REAL(8)              :: V2(5)
      REAL(8)              :: EIG(3)
      REAL(8)              :: U(3,3)
      REAL(8)              :: Y0
      REAL(8)              :: PI
!     **************************************************************************
      IF(.NOT.TWAKE) RETURN
                                CALL TRACE$PUSH('HYPERFINE$PRINT')
!
!     ==========================================================================
!     ==  PREPARE SOME CONSTANTS NEEDED LATER                                 ==
!     ==========================================================================
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL MPE$QUERY('MONOMER',NTASKNUM,NTASKID)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==========================================================================
!     ==  LOOP THROUGH ALL ATOMS                                              ==
!     ==========================================================================
                         CALL TRACE$PASS('BEFORE LOOP')
      CALL REPORT$TITLE(NFILO,'HYPERFINE PARAMETERS') 
      CALL ATOMLIST$NATOM(NAT)
      DO IAT=1,NAT
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
        CALL LINKEDLIST$SELECT(LL_HPRFN,'~')
        CALL LINKEDLIST$EXISTL(LL_HPRFN,ATOMNAME,1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$SELECT(LL_HPRFN,ATOMNAME)
          CALL LINKEDLIST$GET(LL_HPRFN,'TEFG',1,TEFG)
          CALL LINKEDLIST$GET(LL_HPRFN,'TIS',1,TIS)
          CALL LINKEDLIST$GET(LL_HPRFN,'TFC',1,TFC)
          CALL LINKEDLIST$GET(LL_HPRFN,'TANIS',1,TANIS)
        ELSE
          TEFG=.FALSE.
          TFC=.FALSE.
          TIS=.FALSE.
          TANIS=.FALSE.
        END IF
!
!       ========================================================================
!       ==  ISOMER SHIFT                                                      ==
!       ========================================================================
                         CALL TRACE$PASS('BEFORE IS:'//TRIM(ATOMNAME))
        IF(TIS) THEN
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'AERHOTLM',1,TCHK)
          IF(.NOT.TCHK) THEN
            IF(NTASKNUM.EQ.1) THEN
              CALL ERROR$MSG('AERHOTLM NOT PRESENT')
              CALL ERROR$STOP('HYPERFINE$PRINT')            
            END IF
            SVAR=0.D0
          ELSE
            CALL LINKEDLIST$GET(LL_HPRFN,'AERHOTLM',1,SVAR)
          END IF
!
          RHOT0=SVAR*Y0
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSRHOTPW',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_HPRFN,'PSRHOTPW',1,SVAR)
            RHOT0=RHOT0+SVAR
            CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSRHOTLM',1,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_HPRFN,'PSRHOTLM',1,SVAR)
              RHOT0=RHOT0-SVAR*Y0
            END IF
          END IF
          CALL MPE$COMBINE('MONOMER','+',RHOT0)
                         CALL TRACE$PASS('BEFORE PRINT IS')
          IF(NTASKID.EQ.1) THEN
            CALL REPORT$R8VAL(NFILO &
      &                      ,'ELECTRON DENSITY AT THE NUCLEUS OF ATOM ' &
      &                       //TRIM(ATOMNAME),RHOT0,'1/ABOHR^3')
            CALL REPORT$R8VAL(NFILO,'ISOMERSHIFT FOR ATOM '//TRIM(ATOMNAME) &
       &                            ,RHOT0,'E/A_0^3')
          END IF
        END IF
!
!       ========================================================================
!       ==  FERMI CONTACT                                                     ==
!       ========================================================================
                         CALL TRACE$PASS('BEFORE FC')
        IF(TFC.AND.TSPIN) THEN
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'AERHOSLM',1,TCHK)
          IF(.NOT.TCHK) THEN
            IF(NTASKNUM.EQ.1) THEN
              CALL ERROR$MSG('AERHOSLM NOT PRESENT')
              CALL ERROR$STOP('HYPERFINE$PRINT')            
            END IF
            SVAR=0.D0
          ELSE
            CALL LINKEDLIST$GET(LL_HPRFN,'AERHOSLM',1,SVAR)
          END IF
          AERHO1S0=SVAR*Y0
          RHOS0=SVAR*Y0
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSRHOSPW',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_HPRFN,'PSRHOSPW',1,SVAR)
            PSRHOS0=SVAR
            RHOS0=RHOS0+SVAR
            CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSRHOSLM',1,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_HPRFN,'PSRHOSLM',1,SVAR)
               PSRHO1S0=SVAR*Y0
               RHOS0=RHOS0-SVAR*Y0
            END IF
            CALL MPE$COMBINE('MONOMER','+',AERHO1S0)
            CALL MPE$COMBINE('MONOMER','+',PSRHO1S0)
            CALL MPE$COMBINE('MONOMER','+',PSRHOS0)
            CALL MPE$COMBINE('MONOMER','+',RHOS0)
          END IF
          IF(NTASKID.EQ.1) THEN
            SVAR=2.D0/3.D0
            CALL CONSTANTS('MU0',SVAR1)          ; SVAR=SVAR*SVAR1 
!           == MULTIPLY WITH MAGNETIC MOMENT OF THE ELECTRON
            CALL CONSTANTS('BOHRMAGNETON',SVAR1) ; SVAR=SVAR*SVAR1
            CALL CONSTANTS('GE',SVAR1)           ; SVAR=SVAR*SVAR1
            CALL CONSTANTS('HBAR',SVAR1)         ; SVAR=SVAR*0.5D0*SVAR1
!           ==   
            CALL CONSTANTS('TESLA',SVAR1)        ; SVAR=SVAR/SVAR1
            PRINT*,'THIS NUMBER SHOULD BE 104.98 : ',SVAR
            CALL REPORT$R8VAL(NFILO,'SPIN DENSITY AT THE NUCLEUS OF ATOM ' &
      &                            //TRIM(ATOMNAME),RHOS0,'1/ABOHR^3')
            CALL REPORT$R8VAL(NFILO,'FERMI CONTACT HYPERFINE FIELD FOR ATOM ' &
      &                            //TRIM(ATOMNAME),SVAR*RHOS0,'TESLA')
            CALL REPORT$R8VAL(NFILO,'PSEUDO PART OF FERMI CONTACT ' &
      &                            //'HYPERFINE FIELD FOR ATOM ' &
      &                            //TRIM(ATOMNAME),SVAR*PSRHOS0,'TESLA')
            CALL REPORT$R8VAL(NFILO,'ONE-CENTER PART OF FERMI CONTACT ' &
      &                            //'HYPERFINE FIELD FOR ATOM ' &
      &                            //TRIM(ATOMNAME),SVAR*AERHO1S0,'TESLA')
            CALL REPORT$R8VAL(NFILO,'ONE-CENTER-PSEUDO PART OF FERMI CONTACT ' &
      &                            //'HYPERFINE FIELD FOR ATOM ' &
      &                            //TRIM(ATOMNAME),SVAR*PSRHO1S0,'TESLA')
          END IF
        END IF
        IF(TFC.AND.(.NOT.TSPIN).AND.NTASKID.EQ.1) THEN
          CALL REPORT$STRING(NFILO &
      &                     ,'NO FERMI CONTACT HYPERFINE FIELD FOR ATOM ' &
      &                     //TRIM(ATOMNAME)//' BECAUSE NO SPIN POLARIZATION')
        END IF
!
!       ========================================================================
!       ==  ELECTRIC FIELD GRADIENTS                                          ==
!       ========================================================================
                         CALL TRACE$PASS('BEFORE EFG')
        IF(TEFG) THEN
!         == AE PART OF THE EFG ================================================
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'AEV2LM',1,TCHK)
          IF(.NOT.TCHK) THEN
            IF(NTASKNUM.EQ.1) THEN
              CALL ERROR$MSG('AE2VLM NOT PRESENT')
              CALL ERROR$STOP('HYPERFINE$PRINT')            
            END IF
            V2(:)=0.D0
          ELSE
            CALL LINKEDLIST$GET(LL_HPRFN,'AEV2LM',1,V2)
          END IF
!         __TRANSFORM LM REPRESENTATION INTO CARTESIAN COORDINATES______________
!         __VLM=V(R)/R**2___VXYZ=D2V/(DI*DJ)____________________________________
          CALL DTOXYZ(V2,VXYZ)
          AEEFG1(:,:)=VXYZ(:,:)
          EFG(:,:)=VXYZ(:,:)
!         == PS-PS1 PART  OF THE EFG ===========================================
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSV',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_HPRFN,'PSV',1,VXYZ)
            PSEFG(:,:)=VXYZ(:,:)
            EFG(:,:)=EFG(:,:)+VXYZ(:,:)
            CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSV2LM',1,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_HPRFN,'PSV2LM',1,V2)
!             __TRANSFORM LM REPRESENTATION INTO CARTESIAN COORDINATES__________
              CALL DTOXYZ(V2,VXYZ)
              PSEFG1(:,:)=VXYZ(:,:)
              EFG(:,:)=EFG(:,:)-VXYZ(:,:)
            END IF
          END IF
!         ==  PARALLELIZATION: SUM OVER ALL NODES ==============================
                         CALL TRACE$PASS('IN EFG 4')
          CALL MPE$COMBINE('MONOMER','+',EFG)
          CALL MPE$COMBINE('MONOMER','+',PSEFG)
          CALL MPE$COMBINE('MONOMER','+',PSEFG1)
          CALL MPE$COMBINE('MONOMER','+',AEEFG1)
                         CALL TRACE$PASS('IN EFG 5')
!         ==  CONVERSION AND PRINTOUT FOR ELECTRIC FIELD GRADIENTS =============
          IF(NTASKID.EQ.1) THEN
            SVAR=1.D-21
            CALL CONSTANTS('VOLT',SVAR1) ; SVAR=SVAR/SVAR1 
            CALL CONSTANTS('METER',SVAR1) ; SVAR=SVAR*SVAR1**2 
            SVAR=-SVAR   ! THE ELEMENTARY CHARGE IS -1 ELECTRON CHARGE
            PRINT*,SVAR,'SHOULD BE EQUAL TO ',-9.7175D0
!
!           == TOTAL ELECTIC FIELD GRADIENT ====================================
            CALL REPORT$STRING(NFILO,'ELECTRIC FIELD GRADIENT FOR ATOM ' &
      &                      //TRIM(ATOMNAME)//' IN [10**21*V/M**2]')
            CALL LIB$DIAGR8(3,EFG,EIG,U)
            WRITE(NFILO,FMT='(T5,"VALUE",T20,"DIRECTION")')
            DO I=1,3
              WRITE(NFILO,FMT='(T5,F10.5,T20,3F10.5)') &
     &                                        SVAR*EIG(I),(U(J,I),J=1,3)
            ENDDO
!
!           == PRINT PLANE WAVE CONTRIBUTION TO THE EFG ========================
            IF(TPR) THEN
              CALL REPORT$STRING(NFILO &
      &                         ,'PSEUDO ELECTRIC FIELD GRADIENT FOR ATOM ' &
      &                         //TRIM(ATOMNAME)//' IN [10**21*V/M**2]')
              CALL LIB$DIAGR8(3,PSEFG,EIG,U)
              WRITE(NFILO,FMT='(T5,"VALUE",T20,"DIRECTION")')
              DO I=1,3
                WRITE(NFILO,FMT='(T5,F10.5,T20,3F10.5)') &
     &                                                SVAR*EIG(I),(U(J,I),J=1,3)
              ENDDO
            END IF
!
!           == PRINT ONE-CENTER CONTRIBUTION TO THE EFG ========================
            IF(TPR) THEN 
              CALL REPORT$STRING(NFILO &
      &                         ,'ONE-CENTER ELECTIC FIELD GRADIENT FOR ATOM ' &
      &                         //TRIM(ATOMNAME)//' IN [10**21*V/M**2]')
              CALL LIB$DIAGR8(3,AEEFG1,EIG,U)
              WRITE(NFILO,FMT='(T5,"VALUE",T20,"DIRECTION")')
              DO I=1,3
                WRITE(NFILO,FMT='(T5,F10.5,T20,3F10.5)') &
                                SVAR*EIG(I),(U(J,I),J=1,3)
              ENDDO
            END IF
                         CALL TRACE$PASS('IN EFG 8')
          END IF
        END IF
!
!       ========================================================================
!       ==  ANISOTROPIC HYPERFINE PARAMETER                                   ==
!       ========================================================================
                         CALL TRACE$PASS('BEFORE ANIS')
        IF(TANIS.AND.TSPIN) THEN
!         == AE PART OF THE EFG ================================================
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'AERHOS2LM',1,TCHK)
          IF(.NOT.TCHK) THEN
            IF(NTASKNUM.EQ.1) THEN
              CALL ERROR$MSG('AERHOS2LM NOT PRESENT')
              CALL ERROR$STOP('HYPERFINE$PRINT')            
            END IF
            V2(:)=0.D0
          ELSE
            CALL LINKEDLIST$GET(LL_HPRFN,'AERHOS2LM',1,V2)
          END IF
!         __TRANSFORM LM REPRESENTATION INTO CARTESIAN COORDINATES__
          CALL DTOXYZ(V2,VXYZ)
          ANIS(:,:)=VXYZ(:,:)
!         == PS-PS1 PART  OF THE EFG ===========================================
          CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSRHOS2PW',1,TPLANEWAVE)
          IF(TPLANEWAVE) THEN
            CALL LINKEDLIST$GET(LL_HPRFN,'PSRHOS2PW',1,VXYZ)
            ANIS(:,:)=ANIS(:,:)+VXYZ(:,:)
            CALL LINKEDLIST$EXISTD(LL_HPRFN,'PSRHOS2LM',1,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_HPRFN,'PSRHOS2LM',1,V2)
!             __TRANSFORM LM REPRESENTATION INTO CARTESIAN COORDINATES__________
              CALL DTOXYZ(V2,VXYZ)
              ANIS(:,:)=ANIS(:,:)-VXYZ(:,:)
            END IF
          END IF
!         == PARALLELIZATION : SUM OVER ALL NODES ==============================
          CALL MPE$COMBINE('MONOMER','+',ANIS)
!         == CONVERSION AND PRINTOUT FOR ANISOTROPIC HYPERFINE PARAMETERS
          IF(NTASKID.EQ.1) THEN
            SVAR=1.D0/(4.D0*PI)
            CALL CONSTANTS('MU0',SVAR1)          ; SVAR=SVAR*SVAR1 
!           ====================================================================
            CALL CONSTANTS('BOHRMAGNETON',SVAR1) ; SVAR=SVAR*SVAR1
            CALL CONSTANTS('GE',SVAR1)           ; SVAR=SVAR*SVAR1
            CALL CONSTANTS('HBAR',SVAR1)         ; SVAR=SVAR*0.5D0*SVAR1
!           ====================================================================
            CALL CONSTANTS('TESLA',SVAR1)        ; SVAR=SVAR/SVAR1
            PRINT*,'THIS NUMBER SHOULD BE 12.531 : ',SVAR
            CALL REPORT$STRING(NFILO,'ANISOTROPIC HYPERFINE FIELD FOR ATOM ' &
     &                      //TRIM(ATOMNAME)//' IN UNITS OF TESLA')
            CALL REPORT$STRING(NFILO,'MULTIPLY MATRIX WITH A UNITY VECTOR ' &
     &                      //'INTO THE DIRECTION OF THE MAGNETIC FIELD')
            CALL LIB$DIAGR8(3,ANIS,EIG,U)              
!            __ FACTOR 0.5 IS THE ELECTRON SPIN____________
            WRITE(NFILO,FMT='(T5,"VALUE",T20,"DIRECTION")')
            DO I=1,3
              WRITE(NFILO,FMT='(T5,F10.5,T20,3F10.5)')SVAR*EIG(I),(U(J,I),J=1,3)
            ENDDO
          END IF
        END IF
        IF(TANIS.AND. (.NOT.TSPIN).AND.NTASKID.EQ.1)  THEN
          WRITE(NFILO,FMT='("ANISOTROPIC HYPERFINE FIELD FOR ATOM : ",A)')&
     &                    ATOMNAME
          WRITE(NFILO, &
     &         FMT='("NOT OBTAINED BECAUSE CALCULATION IS NON SPIN POLARIZED")')
        END IF
      ENDDO
                                CALL TRACE$POP
      RETURN
      END
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DTOXYZ(VLM_,V2)
!     ****************************************************************
!     **                                                            **
!     **  GIVEN THE VALUE OF F(R)/R**2 FOR THE 5 D-TYPE RADIAL      **
!     **  FUNCTIONS, CALCULATED THE MATRIX OF THE SECOND            **
!     **  DERIVATIVES AT THE ORIGIN                                 **
!     **                                                            **
!     **  THE FIRST NINE REAL SPHERICAL HARMONICS ARE:              **
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2    **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2    **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2    **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2    **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2    **
!     ****************************************************************
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: VLM_(5)
      REAL(8), INTENT(OUT):: V2(3,3)
      REAL(8)             :: VLM(5)
      REAL(8)             :: PI
      REAL(8)             :: SQ15
      REAL(8)             :: SQ5
      REAL(8)             :: SQ60
      REAL(8)             :: SQ16PI
!     ****************************************************************
      PI=4.D0*ATAN(1.D0)
      SQ15=SQRT(15.D0)
      SQ5=SQRT(5.D0)
      SQ60=SQRT(60.D0)
      SQ16PI=SQRT(16.D0*PI)
      VLM(1)=VLM_(1)*SQ15/SQ16PI
      VLM(2)=VLM_(2)*SQ60/SQ16PI
      VLM(3)=VLM_(3)*SQ5/SQ16PI
      VLM(4)=VLM_(4)*SQ60/SQ16PI
      VLM(5)=VLM_(5)*SQ60/SQ16PI
      V2(1,1)=+2.D0*VLM(1)-2.D0*VLM(3)
      V2(1,2)=+VLM(5)
      V2(1,3)=+VLM(2)
      V2(2,2)=-2.D0*VLM(1)-2.D0*VLM(3)
      V2(2,3)=+VLM(4)
      V2(3,3)=+4.D0*VLM(3)
      V2(2,1)=V2(1,2)
      V2(3,1)=V2(1,3)
      V2(3,2)=V2(2,3)
      RETURN
      END SUBROUTINE DTOXYZ
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EFG_THOMSON(GID,NR,RHO,Z,RHO0)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE HYPERFINE_MODULE, ONLY: EXTRAPOLATE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID          ! GRID ID
      INTEGER(4),INTENT(IN) :: NR           ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: RHO(NR)      ! DENSITY
      REAL(8)   ,INTENT(IN) :: Z            ! ATOMIC NUMBER
      REAL(8)   ,INTENT(OUT):: RHO0         ! WEIGHTED INTEGRAL
      REAL(8)               :: ALPHA        ! FINE STRUCTURE CONSTANT
      REAL(8)               :: RT           ! THOMSON RADIUS
      REAL(8)               :: RTH          ! THOMSON RADIUS/2
      REAL(8)               :: LAMBDA       ! PSI(R<<1)=C*R**(LAMBDA-1)
      REAL(8)               :: DLAM         ! 2*(LAMBDA-1)
      REAL(8)               :: R(NR)        !RADIAL GRID
      REAL(8)               :: X1,CI(4)     ! POLYNOMIAL
      REAL(8)               :: RA,RB        ! INTEGRATION BOUNDS
      INTEGER(4)            :: I,IR,J
      REAL(8)               :: Y(4)
      REAL(8)               :: SVAR,FAC,R2,SUM
!     ******************************************************************
      CALL CONSTANTS$GET('ALPHA',ALPHA) ! FINE STRUCTURE CONSTANT
      RT=Z*ALPHA**2    !THOMSON RADIUS
      RTH=0.5D0*RT
      LAMBDA=SQRT(1.D0-(ALPHA*Z)**2)
      DLAM=2.D0*(LAMBDA-1.D0)   ! ASSUMES S-LIKE BEHAVIOR OF WAVE FUNCTIONS
      IF(Z.LT.1) DLAM=0.D0     
!     == PREPARE RADIAL GRID ===========================================
      CALL RADIAL$R(GID,NR,R)
!
!     ================================================================== 
!     == INTEGRATE TO THE SECOND K-POINT                              == 
!     ================================================================== 
!     ==  OBTAIN INTERPOLATING POLYNOMIAL 
      DO IR=1,4
        Y(IR)=RHO(IR)*(R(IR)/RTH)**(-DLAM)
      ENDDO      
      CALL POLYNOM$COEFF(4,X1,CI,R,Y)
      CALL POLYNOM$SHIFTORIGIN(4,X1,CI,0.D0)
!     ==  CLASSICAL EXPRESSION USED FOR Z<1
      IF(Z.LT.1.D0) THEN
        RHO0=CI(1)
        RETURN
      END IF
!     == NOW INTEGRATE FROM THE ORIGIN TO TEH SECOND GRID POINT
      DO I=1,4
        CI(I)=CI(I)*RTH**(I-1)
      ENDDO
      R2=R(2)
      RHO0=0.D0
      FAC=1.D0/(1.D0+RTH/R2)
      DO I=1,4
        SVAR=DLAM+REAL(I,KIND=8)
        CI(I)=CI(I)*(R2/RTH)**SVAR/(SVAR*(1.D0+R2/RTH)**2)
        SUM=CI(I)
        DO J=1,100
          CI(I)=CI(I)*FAC/(1.D0+SVAR/REAL(J,KIND=8))
          SUM=SUM+REAL(J+1,KIND=8)*CI(I)
        ENDDO
        RHO0=RHO0+SUM
      ENDDO
!
!     ================================================================== 
!     == INTEGRATE FROM THE SECOND R-POINT                            == 
!     ================================================================== 
      DO IR=2,NR-3
!       == OBTAIN COEFFICENTS FOR Y=SUM:AI*(X-X2)**I
        CALL POLYNOM$COEFF(4,X1,CI,R(IR-1),RHO(IR-1))
        CALL POLYNOM$SHIFTORIGIN(4,X1,CI,-RTH)
        RA=R(IR)  +RTH
        RB=R(IR+1)+RTH
        RHO0=RHO0+RTH*(CI(1)*(1.D0/RA-1.D0/RB) &
    &                 +CI(2)*(LOG(RB)-LOG(RA)) &
    &                 +CI(3)*(RB-RA) &
    &                 +CI(4)*0.5D0*(RB**2-RA**2))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EFG_THOMSON_OLD(R1,DEX,NR,RHO,Z,RHO0)
!     ****************************************************************
!     **                                                            **
!     ****************************************************************
      USE HYPERFINE_MODULE, ONLY: EXTRAPOLATE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1           ! FIRST GRID POINT
      REAL(8)   ,INTENT(IN) :: DEX          ! LOG-SPACING OF GRID POINTS
      INTEGER(4),INTENT(IN) :: NR           ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: RHO(NR)      ! DENSITY
      REAL(8)   ,INTENT(IN) :: Z            ! ATOMIC NUMBER
      REAL(8)   ,INTENT(OUT):: RHO0         ! WEIGHTED INTEGRAL
      REAL(8)               :: ALPHA        ! FINE STRUCTURE CONSTANT
      REAL(8)               :: RT           ! THOMSON RADIUS
      REAL(8)               :: RTH          ! THOMSON RADIUS/2
      REAL(8)               :: LAMBDA       ! PSI(R<<1)=C*R**(LAMBDA-1)
      INTEGER(4)            :: I,IR
      REAL(8)               :: Y(4),X(4)
      REAL(8)               :: R(NR)        !RADIAL GRID
      REAL(8)               :: RI,XEXP
      REAL(8)               :: A0,A1,A2,A3  ! POWER EXPANSION COEFFICIENTS 
      REAL(8)               :: S0,S1,S2,S3  ! INTEGRATION OF POWERS
      REAL(8)               :: X12,X13,X32,X42,X43,Y13,Y43
      REAL(8)               :: SVAR,FAC
      REAL(8)               :: RA,RB ! INTEGRATION BOUNDS
      LOGICAL(4), PARAMETER :: TTEST=.FALSE.
!     ****************************************************************
      CALL CONSTANTS$GET('ALPHA',ALPHA) ! FINE STRUCTURE CONSTANT
      RT=Z*ALPHA**2    !THOMSON RADIUS
      RTH=0.5D0*RT
      LAMBDA=SQRT(1.D0-(ALPHA*Z)**2)
      XEXP=EXP(DEX)
      IF(TTEST) THEN
        PRINT*,'Z      ',Z
        PRINT*,'RT     ',RT
        PRINT*,'LAMBDA ',LAMBDA
        PRINT*,'ALPHA  ',ALPHA
        PRINT*,'R1/RTH ',R1/RTH
      END IF
      IF(R1/RTH.GT.1.D0) THEN
        RI=R1/XEXP
        DO IR=1,4
          RI=RI*XEXP
          X(IR)=RI
          Y(IR)=RHO(IR)
        ENDDO      
        CALL EXTRAPOLATE(4,X,Y,0.D0,RHO0)
        RETURN
      END IF
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        R(IR)=RI
      ENDDO
!
!     ================================================================== 
!     == INTEGRATE TO THE FIRST K-POINT                               == 
!     ================================================================== 
      RHO0=RHO(1)*R1/(R1+RTH)
      FAC=-RHO(1)
      DO I=0,100
        FAC=FAC*(-R1/RTH)
        SVAR=FAC*(1.D0/(1.D0+2*(LAMBDA-1.D0)/DBLE(I+1))-1.D0)
        RHO0=RHO0+SVAR
        IF(ABS(SVAR).LT.1.D-10) EXIT
      ENDDO
!
!     ================================================================== 
!     == INTEGRATE FROM THE FIRST K-POINT                             == 
!     ================================================================== 
      DO IR=2,NR-3
!       == MAP ARRAY
        DO I=1,4
          Y(I)=RHO(IR-2+I)
          X(I)=R(IR-2+I)
        ENDDO   
!       == OBTAIN COEFFICENTS FOR Y=SUM:AI*(X-X2)**I
        A0=Y(2)
        DO I=1,4
          IF(I.EQ.2) CYCLE   ! SECOND GRID POINT IS NOT USED ANY MORE
          Y(I)=(Y(I)-A0)/(X(I)-X(2))
        ENDDO
        X12=X(1)-X(2)
        X13=X(1)-X(3)
        X32=X(3)-X(2)
        X42=X(4)-X(2)
        X43=X(4)-X(3)
        Y13=Y(1)-Y(3)
        Y43=Y(4)-Y(3)
        A3=(Y43/X43-Y13/X13)/((X42**2-X32**2)/X43-(X12**2-X32**2)/X13)
        A2=(Y43-A3*(X42**2-X32**2))/X43
        A1=Y(3)-A2*X32-A3*X32**2
        IF(TTEST) THEN
          PRINT*,'AI ',A0,A1,A2,A3
          PRINT*,'XI2 ',X12,X32,X42
          PRINT*,'TEST1 ',RHO(IR-1),A0+A1*X12+A2*X12**2+A3*X12**3
          PRINT*,'TEST2 ',RHO(IR),A0
          PRINT*,'TEST3 ',RHO(IR+1),A0+A1*X32+A2*X32**2+A3*X32**3
          PRINT*,'TEST4 ',RHO(IR+2),A0+A1*X42+A2*X42**2+A3*X42**3
        END IF
!
!       == CALCULATE BARE INTEGRALS OF SI=INT:S(R)*(R+RTH)**I ==========
        RA=X(2)
        RB=X(3)
        IF(IR.EQ.2) THEN ! THE FIRST INTERVAL EXTENDS FROM R(1) TO R(3)
          RA=X(1)
          A2=A2+3.D0*A3*X12
          A1=A1+2.D0*A2*X12+3.D0*A3*X12**2
          A0=A0+     A1*X12+     A2*X12**2+A3*X12**3
        END IF
        S0=RTH*(1.D0/(RA+RTH)-1.D0/(RB+RTH))
        S1=RTH*LOG((RB+RTH)/(RA+RTH))
        S2=RTH*(RB-RA)
        S3=RTH*0.5D0*((RB+RTH)**2-(RA+RTH)**2)
!
!       == TRANSFORM TO SI=INT: S(R)*(X-RA)**I =======================
        SVAR=-(RA+RTH)
        S3=S3+3.D0*SVAR*S2+3.D0*SVAR**2*S1+SVAR**3*S0
        S2=S2+2.D0*SVAR*S1+SVAR**2*S0
        S1=S1+SVAR*S0
!
!       ===
        RHO0=RHO0+A0*S0+A1*S1+A2*S2+A3*S3
      ENDDO
      RETURN
      END



