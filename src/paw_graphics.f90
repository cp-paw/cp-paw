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
!.......................................................................
MODULE GRAPHICS_MODULE
TYPE WAVEPLOT_TYPE
 CHARACTER(512)          :: FILE  ! FILE NAME
 CHARACTER(128)          :: TITLE ! IMAGE TITLE 
 REAL(8)                 :: DR    ! STEP SIZE OF THE GRID
 INTEGER(4)              :: IB    ! BAND
 INTEGER(4)              :: IKPT  ! K-POINT
 INTEGER(4)              :: ISPIN ! SPIN INDEX
 LOGICAL(4)              :: TIMAG  ! IMAGINARY OR REAL PART
END TYPE WAVEPLOT_TYPE
TYPE DENSITYPLOT_TYPE
 CHARACTER(512)          :: FILE  ! FILE NAME
 CHARACTER(128)          :: TITLE ! IMAGE TITLE 
 REAL(8)                 :: DR    ! STEP SIZE OF THE GRID
 CHARACTER(8)            :: TYPE  ! CAN BE 'TOTAL','SPIN','UP','DOWN'
 LOGICAL(4)              :: TDIAG
 LOGICAL(4)              :: TOCC
 LOGICAL(4)              :: TCORE
 CHARACTER(16)           :: SELECTOR
 REAL(8)                 :: EMIN
 REAL(8)                 :: EMAX
 INTEGER(4)              :: IBMIN
 INTEGER(4)              :: IBMAX
END TYPE DENSITYPLOT_TYPE
TYPE POTPLOT_TYPE
 CHARACTER(512)          :: FILE  ! FILE NAME
 CHARACTER(128)          :: TITLE ! IMAGE TITLE 
 REAL(8)                 :: DR    ! STEP SIZE OF THE GRID
END TYPE POTPLOT_TYPE
type onecrho_type
  real(8),pointer        :: ae1cpot(:,:)
  real(8),pointer        :: ps1cpot(:,:)
end type onecrho_type
TYPE(DENSITYPLOT_TYPE),ALLOCATABLE :: DENSITYPLOT(:)
TYPE(WAVEPLOT_TYPE)   ,ALLOCATABLE :: WAVEPLOT(:)
TYPE(POTPLOT_TYPE)                 :: POTPLOT
LOGICAL(4)                :: TINI=.FALSE.
LOGICAL(4)                :: TWAKE=.FALSE.
INTEGER(4)                :: IWAVEPTR=0
INTEGER(4)                :: IDENSITYPTR=0
INTEGER(4)                :: IPOTPTR=0
COMPLEX(8),ALLOCATABLE    :: PWPOT(:)
type(onecrho_type),allocatable :: onecpotarray(:)  ! shALL REPLACE AE1CPOT AND PS1CPOT
REAL(8)   ,ALLOCATABLE    :: AE1CPOT(:,:,:)  !(NRX,LMRX,NAT)
REAL(8)   ,ALLOCATABLE    :: PS1CPOT(:,:,:)  !(NRX,LMRX,NAT)
INTEGER(4)                :: LMRXX=0    !INITIALLY NOT SET
INTEGER(4)                :: NRX=0
INTEGER(4)                :: NGL=0
INTEGER(4)                :: NWAVE=0
INTEGER(4)                :: NDENSITY=0
LOGICAL(4)                :: TPOT=.FALSE.
END MODULE GRAPHICS_MODULE
!
!     ..................................................................
      SUBROUTINE GRAPHICS$SETI4(ID,VAL)
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      USE GRAPHICS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NWAVE') THEN 
        IF(ALLOCATED(WAVEPLOT)) THEN
          CALL ERROR$MSG('WAVEPLOT ALREADY ALLOCATED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETI4')
        END IF
        NWAVE=VAL
        ALLOCATE(WAVEPLOT(NWAVE))
!
      ELSE IF(ID.EQ.'NDENSITY') THEN 
        IF(ALLOCATED(DENSITYPLOT)) THEN
          CALL ERROR$MSG('DENSITYPLOT ALREADY ALLOCATED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETI4')
        END IF
        NDENSITY=VAL
        ALLOCATE(DENSITYPLOT(NDENSITY))
!
      ELSE IF(ID.EQ.'IWAVE') THEN 
        IDENSITYPTR=0
        IPOTPTR=0
        IWAVEPTR=VAL
        TINI=.TRUE.
!
      ELSE IF(ID.EQ.'IDENSITY') THEN 
        IWAVEPTR=0
        IPOTPTR=0
        IDENSITYPTR=VAL
        TINI=.TRUE.
!
      ELSE IF(ID.EQ.'IB') THEN 
        IF(IWAVEPTR.EQ.0) THEN
          CALL ERROR$MSG('WAVE PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETI4')
        END IF
        WAVEPLOT(IWAVEPTR)%IB=VAL
!
      ELSE IF(ID.EQ.'IKPT') THEN 
        IF(IWAVEPTR.EQ.0) THEN
          CALL ERROR$MSG('WAVE PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETI4')
        END IF
        WAVEPLOT(IWAVEPTR)%IKPT=VAL
!
      ELSE IF(ID.EQ.'ISPIN') THEN 
        IF(IWAVEPTR.EQ.0) THEN
          CALL ERROR$MSG('WAVE PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETI4')
        END IF
        WAVEPLOT(IWAVEPTR)%ISPIN=VAL
!
      ELSE IF(ID.EQ.'IBMIN') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETI4')
        END IF
        DENSITYPLOT(IDENSITYPTR)%IBMIN=VAL
!
      ELSE IF(ID.EQ.'IBMAX') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETI4')
        END IF
        DENSITYPLOT(IDENSITYPTR)%IBMAX=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GRAPHICS$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS$SETL4(ID,VAL)
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      USE GRAPHICS_MODULE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'WAKE') THEN 
        TWAKE=VAL
!
      ELSE IF(ID.EQ.'IMAG') THEN 
        IF(IWAVEPTR.EQ.0) THEN
          CALL ERROR$MSG('WAVE PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETL4')
        END IF
        WAVEPLOT(IWAVEPTR)%TIMAG=VAL
!
      ELSE IF(ID.EQ.'TDIAG') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETL4')
        END IF
        DENSITYPLOT(IDENSITYPTR)%TDIAG=VAL
!
      ELSE IF(ID.EQ.'TOCC') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETL4')
        END IF
        DENSITYPLOT(IDENSITYPTR)%TOCC=VAL
!
      ELSE IF(ID.EQ.'TCORE') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETL4')
        END IF
        DENSITYPLOT(IDENSITYPTR)%TCORE=VAL
!
      ELSE IF(ID.EQ.'TPOT') THEN 
        TPOT=VAL
        IDENSITYPTR=0
        IWAVEPTR=0
        IPOTPTR=0
        IF(TPOT)IPOTPTR=1
        TINI=.TRUE.
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GRAPHICS$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS$SETCH(ID,VAL)
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      USE GRAPHICS_MODULE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'FILE') THEN 
        IF(IWAVEPTR.NE.0) THEN
          WAVEPLOT(IWAVEPTR)%FILE=VAL
        ELSE IF(IDENSITYPTR.NE.0) THEN
          DENSITYPLOT(IDENSITYPTR)%FILE=VAL
        ELSE IF(IPOTPTR.NE.0) THEN
          POTPLOT%FILE=VAL
        ELSE
          CALL ERROR$MSG('NEITHER WAVE OR DENSITY OR POT PLOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETCH')
        END IF
!
      ELSE IF(ID.EQ.'TITLE') THEN 
        IF(IWAVEPTR.NE.0) THEN
          WAVEPLOT(IWAVEPTR)%TITLE=VAL
        ELSE IF(IDENSITYPTR.NE.0) THEN
          DENSITYPLOT(IDENSITYPTR)%TITLE=VAL
        ELSE IF(IPOTPTR.NE.0) THEN
          POTPLOT%TITLE=VAL
        ELSE
          CALL ERROR$MSG('NEITHER WAVE OR DENSITY OR POT PLOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETCH')
        END IF
!
      ELSE IF(ID.EQ.'TYPE') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETCH')
        END IF
        DENSITYPLOT(IDENSITYPTR)%TYPE=VAL
!
      ELSE IF(ID.EQ.'SELECT') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETCH')
        END IF
        DENSITYPLOT(IDENSITYPTR)%SELECTOR=VAL
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GRAPHICS$SETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS$SETR8(ID,VAL)
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      USE GRAPHICS_MODULE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
      REAL(8)                 :: PI,Y0
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      IF(ID.EQ.'DR') THEN 
        IF(IWAVEPTR.NE.0) THEN
          WAVEPLOT(IWAVEPTR)%DR=VAL
        ELSE IF(IDENSITYPTR.NE.0) THEN
          DENSITYPLOT(IDENSITYPTR)%DR=VAL
        ELSE IF (IPOTPTR.NE.0) THEN
          POTPLOT%DR=VAL
        ELSE
          CALL ERROR$MSG('NEITHER WAVE OR DENSITY PLOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETR8')
        END IF
!
      ELSE IF(ID.EQ.'EMIN') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETR8')
        END IF
        DENSITYPLOT(IDENSITYPTR)%EMIN=VAL
!
      ELSE IF(ID.EQ.'EMAX') THEN 
        IF(IDENSITYPTR.EQ.0) THEN
          CALL ERROR$MSG('DENSITY PLOT NOT SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GRAPHICS$SETR8')
        END IF
        DENSITYPLOT(IDENSITYPTR)%EMAX=VAL
!
      ELSE IF(ID.EQ.'POTSHIFT') THEN 
        IF(IPOTPTR.EQ.0) RETURN
        IF(ALLOCATED(PWPOT)) THEN
          PWPOT(:)=PWPOT(:)+VAL
        END IF
        IF(ALLOCATED(AE1CPOT).AND.ALLOCATED(PS1CPOT)) THEN
          AE1CPOT(:,1,:)=AE1CPOT(:,1,:)+VAL/Y0
          PS1CPOT(:,1,:)=PS1CPOT(:,1,:)+VAL/Y0
        END IF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GRAPHICS$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS$PLOT
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      USE GRAPHICS_MODULE
      CHARACTER(512) ::  FILE
      CHARACTER(128) :: TITLE
      REAL(8)        :: DR
      LOGICAL(4)     :: TIMAG
      LOGICAL(4)     :: TDIAG
      LOGICAL(4)     :: TOCC
      LOGICAL(4)     :: TCORE
      CHARACTER(16)  :: SELECTOR
      CHARACTER(16)  :: TYPE
      REAL(8)        :: E,ELAST
      INTEGER(4)     :: I,IB,IKPT,ISPIN
      INTEGER(4),ALLOCATABLE :: IBMIN(:,:),IBMAX(:,:)
      REAL(8)   ,ALLOCATABLE :: EB(:,:,:) !(NB,NKPT,NSPIN) E-EXPECTATIONV.
      REAL(8)   ,ALLOCATABLE :: EIGVAL(:,:,:) !(NB,NKPT,NSPIN) E-EIGENVALUES
      REAL(8)                :: EMIN,EMAX
!     ******************************************************************
      IF(.NOT.TINI) RETURN
                              CALL TRACE$PUSH('GRAPHICS$PLOT')
!
!     ==================================================================
!     == PLOT POTENTIAL                                               ==
!     ==================================================================
      IF(TPOT) THEN
        FILE=POTPLOT%FILE
        TITLE=POTPLOT%TITLE
        DR=POTPLOT%DR
        CALL GRAPHICS_CREATEPOT(FILE,TITLE,DR)
      END IF
!
!     ==================================================================
!     == PLOT WAVE FUNCTIONS                                          ==
!     ==================================================================
      DO I=1,NWAVE
        FILE=WAVEPLOT(I)%FILE
        TITLE=WAVEPLOT(I)%TITLE
        DR=WAVEPLOT(I)%DR
        IB=WAVEPLOT(I)%IB
        IKPT=WAVEPLOT(I)%IKPT
        ISPIN=WAVEPLOT(I)%ISPIN
        TIMAG=WAVEPLOT(I)%TIMAG
        CALL GRAPHICS_WAVEPLOT(FILE,TITLE,DR,IB,IKPT,ISPIN,TIMAG)
      ENDDO
!
!     ==================================================================
!     == PLOT DENSITIES                                               ==
!     ==================================================================
      IF(.NOT.ALLOCATED(IBMIN)) THEN
        CALL WAVES$GETI4('NB',NB)
        CALL WAVES$GETI4('NKPT',NKPT)
        CALL WAVES$GETI4('NSPIN',NSPIN)
        ALLOCATE(IBMIN(NKPT,NSPIN))
        ALLOCATE(IBMAX(NKPT,NSPIN))
      END IF
      ALLOCATE(EB(NB,NKPT,NSPIN))
      ALLOCATE(EIGVAL(NB,NKPT,NSPIN))
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB,EB(:,IKPT,ISPIN))
          CALL WAVES$GETR8A('EIGVAL',NB,EIGVAL(:,IKPT,ISPIN))
        ENDDO
      ENDDO
!
      DO I=1,NDENSITY
        FILE=DENSITYPLOT(I)%FILE
        TITLE=DENSITYPLOT(I)%TITLE
        DR=DENSITYPLOT(I)%DR
        TYPE=DENSITYPLOT(I)%TYPE
        TDIAG=DENSITYPLOT(I)%TDIAG
        TOCC=DENSITYPLOT(I)%TOCC
        TCORE=DENSITYPLOT(I)%TCORE
        SELECTOR=DENSITYPLOT(I)%SELECTOR
        IF(SELECTOR.EQ.'EMINMAX') THEN
          EMIN=DENSITYPLOT(I)%EMIN
          EMAX=DENSITYPLOT(I)%EMAX
          IBMIN(:,:)=NB
          IBMAX(:,:)=1
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              ELAST=-1.D+10
              DO IB=1,NB
                IF(TDIAG) THEN
                  E=EIGVAL(IB,IKPT,ISPIN)
                ELSE
                  E=EB(IB,IKPT,ISPIN)
                END IF
                IF(E.GE.EMIN.AND.E.LE.EMAX) THEN
                  IBMIN(IKPT,ISPIN)=MIN(IB,IBMIN(IKPT,ISPIN))
                  IBMAX(IKPT,ISPIN)=MAX(IB,IBMIN(IKPT,ISPIN))
                END IF
                IF(E.LT.ELAST) THEN
                  CALL ERROR$MSG('BAND ENERGIES NOT IN INCREASING ORDER ')
                  CALL ERROR$STOP('GRAPHICS$PLOT')
                END IF
                ELAST=E
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(SELECTOR.EQ.'BMINMAX') THEN
          IBMIN(:,:)=DENSITYPLOT(I)%IBMIN
          IBMAX(:,:)=DENSITYPLOT(I)%IBMAX
!
        ELSE IF(SELECTOR.EQ.'NONE') THEN
          IBMIN(:,:)=1
          IBMAX(:,:)=NB
        ELSE
          CALL ERROR$MSG('SELECTOR NOT RECOGNIZED')
          CALL ERROR$CHVAL('SELECTOR',SELECTOR)
          CALL ERROR$STOP('GRAPHICS$PLOT')
        END IF
        CALL GRAPHICS_DENSITYPLOT(FILE,TITLE,DR,TYPE,TDIAG,TOCC,TCORE &
     &                               ,NKPT,NSPIN,IBMIN,IBMAX)
      ENDDO
      DEALLOCATE(IBMIN)
      DEALLOCATE(IBMAX)
      DEALLOCATE(EB)
      DEALLOCATE(EIGVAL)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_WAVEPLOT(FILE,TITLE,DR,IB,IKPT,ISPIN,TIMAG)
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(512),INTENT(IN) :: FILE
      CHARACTER(128),INTENT(IN) :: TITLE
      REAL(8)       ,INTENT(IN) :: DR    ! STEP SIZE OF THE GRID
      INTEGER(4)    ,INTENT(IN) :: IB    ! BAND
      INTEGER(4)    ,INTENT(IN) :: IKPT  ! K-POINT
      INTEGER(4)    ,INTENT(IN) :: ISPIN ! SPIN INDEX
      LOGICAL(4)    ,INTENT(IN) :: TIMAG  ! IMAGINARY OR REAL PART
      INTEGER(4)                :: NTASKS,THISTASK
      INTEGER(4)                :: NR1,NR1L,NR2,NR3
      INTEGER(4)                :: NNR,NNRL
      INTEGER(4)                :: NR1START
      INTEGER(4)                :: NB
      INTEGER(4)                :: NKPT
      INTEGER(4)                :: NSPIN
      LOGICAL(4)                :: SAVETRAWSTATES !USED TO RESTORE ORIGINAL
                                              ! STATE OF WAVES OBJECT
      REAL(8)                   :: GBAS(3,3),DET ! DUMMY
      INTEGER(4)                :: FACT
      REAL(8)       ,ALLOCATABLE:: WAVE(:,:,:) !(NR1L,NR2,NR3) LATER REFINED
      REAL(8)       ,ALLOCATABLE:: WAVEBIG(:,:,:)
      REAL(8)       ,ALLOCATABLE:: WORK1(:,:,:)
      INTEGER(4)                :: NAT
      REAL(8)                   :: RBAS(3,3)
      CHARACTER(32) ,ALLOCATABLE:: ATOMNAME(:)   !NAT
      REAL(8)       ,ALLOCATABLE:: POS(:,:)   !(3,NAT)
      REAL(8)       ,ALLOCATABLE:: Z(:)
      REAL(8)       ,ALLOCATABLE:: Q(:)
      INTEGER(4)                :: LMNXX
      REAL(8)                   :: R1,DEX
      INTEGER(4)                :: NR
      INTEGER(4)                :: LNX
      INTEGER(4)    ,ALLOCATABLE:: LOX(:)    !LNX
      REAL(8)       ,ALLOCATABLE:: AEPHI(:,:) !NR,LNX
      REAL(8)       ,ALLOCATABLE:: PSPHI(:,:) !NR,LNX
      REAL(8)       ,ALLOCATABLE:: PROJ(:)    !LMNXX
      INTEGER(4)                :: LX,LMXX
      REAL(8)       ,ALLOCATABLE:: DRHOL(:,:) !(NR,LMXX)
      REAL(8)       ,ALLOCATABLE:: DENMAT(:,:) !(LMNXX,LMNXX)
      INTEGER(4)                :: IAT,ISP,LN
      INTEGER(4)                :: NFIL
      INTEGER(4)                :: NR1B,NR2B,NR3B
      INTEGER(4)                :: gid    ! grid id
!     ******************************************************************
                              CALL TRACE$PUSH('GRAPHICS_WAVEPLOT')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     =================================================================
!     =================================================================
!     == PLANE WAVE PART         
!     =================================================================
!     =================================================================
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
      CALL WAVES$GETL4('RAWSTATES',SAVETRAWSTATES)  ! USED TO RESET WAVE OBJECT INTO ORIGINAL STATE
!
      CALL WAVES$SETL4('RAWSTATES',.FALSE.)
      CALL WAVES$SETI4('IB',IB)
      CALL WAVES$SETI4('IKPT',IKPT)
      CALL WAVES$SETI4('ISPIN',ISPIN)
      CALL WAVES$SETL4('TIM',TIMAG)
!
!     =================================================================
!     ==  FACTOR FOR GRID REFINEMENT                                 ==
!     =================================================================
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL GBASS(RBAS,GBAS,DET)
      FACT=NINT(((DET/DR**3)/REAL(NR1*NR2*NR3,KIND=8))**(1./3.))
      FACT=MAX(1,FACT)
!
!     =================================================================
!     ==  GET GENERIC FROM SETUPS OBJECT                             ==
!     =================================================================
      ALLOCATE(WAVE(NR1L,NR2,NR3))
      CALL WAVES$GETR8A('PSPSI',NNRL,WAVE)
!     
!     ================================================================
!     ==  EXPAND TO A FINER R-GRID                                  ==
!     ================================================================
      NR1B=FACT*NR1
      NR2B=FACT*NR2
      NR3B=FACT*NR3
!     == TAKE VALUES COMPATIBLE WITH FFT ROUTINES ====================
      CALL LIB$FFTADJUSTGRD(NR1B)
      CALL LIB$FFTADJUSTGRD(NR2B)
      CALL LIB$FFTADJUSTGRD(NR3B)
      ALLOCATE(WAVEBIG(NR1B,NR2B,NR3B))
      IF(FACT.EQ.1) THEN
        CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WAVEBIG)
      ELSE
        ALLOCATE(WORK1(NR1,NR2,NR3))
        CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WORK1)
        CALL GRAPHICS_REFINEGRID(NR1,NR2,NR3,NR1B,NR2B,NR3B,WORK1,WAVEBIG)
        DEALLOCATE(WORK1)
      END IF
      DEALLOCATE(WAVE)
!
!     =================================================================
!     =================================================================
!     == ONE-CENTER CONTRIBUTION
!     =================================================================
!     =================================================================
!
!     =================================================================
!     ==  GET GENERIC FROM ATOM OBJECT                               ==
!     =================================================================
      CALL ATOMLIST$NATOM(NAT)
      CALL SETUP$LMNXX(LMNXX)
      ALLOCATE(POS(3,NAT))
      ALLOCATE(ATOMNAME(NAT))
      ALLOCATE(Z(NAT))
      ALLOCATE(Q(NAT))
      DO IAT=1,NAT
CALL TRACE$PASS('MARKE 3')
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME(IAT))
        CALL ATOMLIST$GETR8A('R(0)',IAT,3,POS(:,IAT))
        CALL ATOMLIST$GETR8('Z',IAT,Z(IAT))
        CALL ATOMLIST$GETR8('Q',IAT,Q(IAT))
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
!       __ GET PARTIAL WAVES__________________________________________
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
        ALLOCATE(PROJ(LMNXX))
        CALL WAVES$SETI4('IAT',IAT)
        CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROJ)
        LX=0
        DO LN=1,LNX
          LX=MAX(LX,LOX(LN))
        ENDDO
        LMXX=(LX+1)**2
        ALLOCATE(DRHOL(NR,LMXX))
        CALL GRAPHICS_1CWAVE(NR,LNX,LOX,AEPHI,PSPHI,LMNXX &
     &                        ,PROJ,LMXX,DRHOL)
        DEALLOCATE(PROJ)
        CALL GRAPHICS_RHOLTOR(RBAS,NR1B,NR2B,NR3B &
     &         ,1,NR1B,WAVEBIG,POS(:,IAT),gid,NR,LMXX,DRHOL)
        DEALLOCATE(DRHOL)
        DEALLOCATE(LOX)
        DEALLOCATE(AEPHI)
        DEALLOCATE(PSPHI)
      ENDDO  
!     
!     ================================================================
!     ==  PRINT WAVE                                                ==
!     ================================================================
                            CALL TRACE$PASS('PRINT')
      PRINT*,'PRINTING OUT WAVE ',TITLE(1:50)
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$SETFILE('WAVEPLOT',.FALSE.,TRIM(FILE))
        CALL FILEHANDLER$SETSPECIFICATION('WAVEPLOT','FORM','UNFORMATTED')
        CALL FILEHANDLER$UNIT('WAVEPLOT',NFIL)
        CALL WRITEWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,Q,ATOMNAME &
    &                     ,NR1B,NR2B,NR3B,WAVEBIG)
        CALL FILEHANDLER$CLOSE('WAVEPLOT')
      END IF
!
!     ==================================================================
!     == CLOSE DOWN                                                   ==
!     ==================================================================
      CALL WAVES$SETL4('RAWSTATES',SAVETRAWSTATES)
      DEALLOCATE(WAVEBIG)
      DEALLOCATE(ATOMNAME)
      DEALLOCATE(Z)
      DEALLOCATE(Q)
      DEALLOCATE(POS)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_DENSITYPLOT(FILE,TITLE,DR,TYPE,TDIAG,TOCC,TCORE &
     &                               ,NKPT,NSPIN,IBMIN,IBMAX)
!     ******************************************************************
!     **  PLOT                                                        **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(512),INTENT(IN) :: FILE
      CHARACTER(128),INTENT(IN) :: TITLE
      REAL(8)       ,INTENT(IN) :: DR    ! STEP SIZE OF THE GRID
      LOGICAL(4)    ,INTENT(IN) :: TCORE ! INCLUDE CORE STATES
      LOGICAL(4)    ,INTENT(IN) :: TDIAG ! USE EIGENSTATES INSTEAD OF RAW STATES
      LOGICAL(4)    ,INTENT(IN) :: TOCC  ! MULTIPLY WITH OCCUPATIONS
      CHARACTER(4)  ,INTENT(IN) :: TYPE  ! CAN BE 'TOTAL','SPIN',UP,DOWN
      INTEGER(4)    ,INTENT(IN) :: NKPT
      INTEGER(4)    ,INTENT(IN) :: NSPIN
      INTEGER(4)    ,INTENT(IN) :: IBMIN(NKPT,NSPIN)
      INTEGER(4)    ,INTENT(IN) :: IBMAX(NKPT,NSPIN)
      INTEGER(4)                :: NTASKS,THISTASK
      INTEGER(4)                :: NR1,NR1L,NR2,NR3
      INTEGER(4)                :: NNR,NNRL
      INTEGER(4)                :: NR1START
      INTEGER(4)                :: NB
      INTEGER(4)                :: NSPIN_,NKPT_ ! DUMMY
      REAL(8)       ,ALLOCATABLE:: OCC(:,:,:)  !(NB,NKPT,NSPIN)
      LOGICAL(4)                :: SAVETRAWSTATES !USED TO RESTORE ORIGINAL
                                              ! STATE OF WAVES OBJECT
      REAL(8)                   :: RBAS(3,3)
      REAL(8)                   :: GBAS(3,3),DET ! DUMMY
      INTEGER(4)                :: FACT
      REAL(8)       ,ALLOCATABLE:: WAVE(:,:,:) !(NR1L,NR2,NR3) LATER REFINED
      REAL(8)       ,ALLOCATABLE:: WAVEBIG(:,:,:)
      REAL(8)       ,ALLOCATABLE:: WORK1(:,:,:)
      INTEGER(4)                :: NAT
      CHARACTER(32) ,ALLOCATABLE:: ATOMNAME(:)   !NAT
      REAL(8)       ,ALLOCATABLE:: POS(:,:)   !(3,NAT)
      REAL(8)       ,ALLOCATABLE:: Z(:)
      REAL(8)       ,ALLOCATABLE:: Q(:)
      INTEGER(4)                :: LMNXX
      INTEGER(4)                :: NR
      INTEGER(4)                :: LNX
      INTEGER(4)    ,ALLOCATABLE:: LOX(:)    !LNX
      REAL(8)       ,ALLOCATABLE:: AEPHI(:,:) !NR,LNX
      REAL(8)       ,ALLOCATABLE:: PSPHI(:,:) !NR,LNX
      REAL(8)       ,ALLOCATABLE:: PROJ(:)    !LMNXX
      INTEGER(4)                :: LX,LMXX
      REAL(8)       ,ALLOCATABLE:: DRHOL(:,:) !(NR,LMXX)
      REAL(8)       ,ALLOCATABLE:: DENMAT(:,:) !(LMNXX,LMNXX)
      INTEGER(4)                :: IAT,ISP,LN,ISPIN,IKPT,IB
      REAL(8)                   :: FAC,SPINFAC
      INTEGER(4)                :: NFIL
      REAL(8)      ,ALLOCATABLE :: AECORE(:)
      REAL(8)                   :: SVAR,XEXP,RI
      INTEGER(4)                :: I
      INTEGER(4)                :: NR1B,NR2B,NR3B
      INTEGER(4)                :: gid
!     ******************************************************************
                              CALL TRACE$PUSH('GRAPHICS_DENSITYPLOT')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
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
      CALL WAVES$GETI4('NKPT',NKPT_)
      CALL WAVES$GETI4('NSPIN',NSPIN_)
      IF(MAXVAL(IBMAX).GT.NB.OR.MINVAL(IBMIN).LT.1) THEN
        CALL ERROR$MSG('BANDS OUT OF RANGE')
        CALL ERROR$STOP('GRAPHICS_DENSITYPLOT')
      END IF
      IF(NKPT_.NE.NKPT.OR.NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('NKPT OR NSPIN INCONSISTENT')
        CALL ERROR$STOP('GRAPHICS_DENSITYPLOT')
      END IF
      CALL WAVES$GETL4('RAWSTATES',SAVETRAWSTATES)      
      ALLOCATE(OCC(NB,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
!
!     =================================================================
!     ==  SWITCH WAVES OBJECT TO RAW OR EIGEN STATES                 ==
!     =================================================================
      CALL WAVES$SETL4('RAWSTATES',.NOT.TDIAG)
!
!     =================================================================
!     ==  FACTOR FOR GRID REFINEMENT                                 ==
!     =================================================================
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL GBASS(RBAS,GBAS,DET)
      FACT=NINT(((DET/DR**3)/REAL(NR1*NR2*NR3,KIND=8))**(1./3.))
      FACT=MAX(1,FACT)
!     
!     ================================================================
!     ==  PSEUDO WAVE FUNCTIONS/DENSITIES                           ==
!     ================================================================
                            CALL TRACE$PASS('PSEUDO WAVE FUNCTIONS')
      ALLOCATE(WAVE(NR1L,NR2,NR3))
      WAVE(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        IF(TYPE.EQ.'UP'.AND.ISPIN.EQ.2) CYCLE 
        IF(TYPE.EQ.'DOWN'.AND.ISPIN.EQ.1) CYCLE 
        SPINFAC=1.D0
        IF(TYPE.EQ.'SPIN'.AND.ISPIN.EQ.2) SPINFAC=-1.D0
        DO IKPT=1,NKPT
          DO IB=IBMIN(IKPT,ISPIN),IBMAX(IKPT,ISPIN)
            FAC=SPINFAC
            IF(TOCC)FAC=SPINFAC*OCC(IB,IKPT,ISPIN)
            CALL GRAPHICS_ADDRHO(IB,IKPT,ISPIN,FAC,NNRL,WAVE)
          ENDDO
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  EXPAND TO A FINER R-GRID                                  ==
!     ================================================================
      NR1B=FACT*NR1
      NR2B=FACT*NR2
      NR3B=FACT*NR3
!     == TAKE VALUES COMPATIBLE WITH FFT ROUTINES ====================
      CALL LIB$FFTADJUSTGRD(NR1B)
      CALL LIB$FFTADJUSTGRD(NR2B)
      CALL LIB$FFTADJUSTGRD(NR3B)
      ALLOCATE(WAVEBIG(NR1B,NR2B,NR3B))
      IF(FACT.EQ.1) THEN
        CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WAVEBIG)
      ELSE
        ALLOCATE(WORK1(NR1,NR2,NR3))
        CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WORK1)
        CALL GRAPHICS_REFINEGRID(NR1,NR2,NR3,NR1B,NR2B,NR3B,WORK1,WAVEBIG)
        DEALLOCATE(WORK1)
      END IF
      DEALLOCATE(WAVE)
!     
!     ================================================================
!     ================================================================
!     ==  ONE-CENTER EXPANSIONS                                     ==
!     ================================================================
!     ================================================================
!
!     =================================================================
!     ==  GET GENERIC FROM ATOM OBJECT                               ==
!     =================================================================
      CALL SETUP$LMNXX(LMNXX)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(POS(3,NAT))
      ALLOCATE(ATOMNAME(NAT))
      ALLOCATE(Z(NAT))
      ALLOCATE(Q(NAT))
      DO IAT=1,NAT
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME(IAT))
        CALL ATOMLIST$GETR8A('R(0)',IAT,3,POS(:,IAT))
        CALL ATOMLIST$GETR8('Z',IAT,Z(IAT))
        CALL ATOMLIST$GETR8('Q',IAT,Q(IAT))
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$LOFLN(ISP,LNX,LOX)
        ALLOCATE(AEPHI(NR,LNX))
        ALLOCATE(PSPHI(NR,LNX))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
        CALL SETUP$PSPARTIALWAVES(ISP,NR,LNX,PSPHI)
!     
!       __ GET PROJECTIONS____________________________________________
        CALL WAVES$SETI4('IAT',IAT)
        ALLOCATE(DENMAT(LMNXX,LMNXX))
        DENMAT(:,:)=0.D0
        DO ISPIN=1,NSPIN
          IF(TYPE.EQ.'UP'.AND.ISPIN.EQ.2) CYCLE 
          IF(TYPE.EQ.'DOWN'.AND.ISPIN.EQ.1) CYCLE 
          SPINFAC=1.D0
          IF(TYPE.EQ.'SPIN'.AND.ISPIN.EQ.2) SPINFAC=-1.D0
          DO IKPT=1,NKPT
            DO IB=IBMIN(IKPT,ISPIN),IBMAX(IKPT,ISPIN)
              FAC=SPINFAC
              IF(TOCC)FAC=SPINFAC*OCC(IB,IKPT,ISPIN)
              CALL GRAPHICS_ADDDENMAT(IB,IKPT,ISPIN,FAC,LMNXX,DENMAT)
            ENDDO
          ENDDO
        ENDDO
        LMXX=9
        ALLOCATE(DRHOL(NR,LMXX))
        CALL GRAPHICS_1CRHO(NR,LNX,LOX,AEPHI,PSPHI,LMNXX &
     &               ,DENMAT(1,1),LMXX,DRHOL)
        DEALLOCATE(DENMAT)
        CALL GRAPHICS_RHOLTOR(RBAS,NR1B,NR2B,NR3B &
    &         ,1,NR1B,WAVEBIG,POS(:,IAT),gid,NR,LMXX,DRHOL)
!     
!     ================================================================
!     ==  ADD CORE CHARGE DENSITY                                   ==
!     ================================================================
        IF (TCORE.AND.(.NOT.(TYPE.EQ.'SPIN'))) THEN
          ALLOCATE(AECORE(NR))
          CALL SETUP$AECORE(ISP,NR,AECORE)
          IF(TYPE.EQ.'UP'.OR.TYPE.EQ.'DOWN') AECORE(:)=0.5D0*AECORE(:)
          CALL GRAPHICS_RHOLTOR(RBAS,NR1B,NR2B,NR3B &
    &           ,1,NR1B,WAVEBIG,POS(:,IAT),gid,NR,1,AECORE)
          DEALLOCATE(AECORE)
        END IF
        DEALLOCATE(DRHOL)
        DEALLOCATE(LOX)
        DEALLOCATE(AEPHI)
        DEALLOCATE(PSPHI)
      ENDDO  
!     
!     ================================================================
!     ==  PRINT WAVE                                                ==
!     ================================================================
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$SETFILE('WAVEPLOT',.FALSE.,TRIM(FILE))
        CALL FILEHANDLER$SETSPECIFICATION('WAVEPLOT','FORM','UNFORMATTED')
        CALL FILEHANDLER$UNIT('WAVEPLOT',NFIL)
        CALL WRITEWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,Q,ATOMNAME &
    &                     ,NR1B,NR2B,NR3B,WAVEBIG)
        CALL FILEHANDLER$CLOSE('WAVEPLOT')
      END IF
!     
!     ================================================================
!     ==  CLOSE DOWN                                                ==
!     ================================================================
      CALL WAVES$SETL4('RAWSTATES',SAVETRAWSTATES)
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
      CALL WAVES$SETL4('TIM',.TRUE.)
      CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROJ)
      DO LMN1=1,LMNXX
        DO LMN2=1,LMNXX
          DENMAT(LMN1,LMN2)=DENMAT(LMN1,LMN2) &
     &                     +FAC*PROJ(LMN1)*PROJ(LMN2)
        ENDDO
      ENDDO
      CALL WAVES$SETL4('TIM',.FALSE.)
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
      CALL WAVES$SETL4('TIM',.TRUE.)
      CALL WAVES$GETR8A('PSPSI',NNR,PSI)
      DO IR=1,NNR            
        WAVE(IR)=WAVE(IR)+FAC*PSI(IR)**2            
      ENDDO
      CALL WAVES$SETL4('TIM',.FALSE.)
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
     &           ,gid,NR,LMX,DRHOL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER     :: LMXX=36
      REAL(8)   ,PARAMETER     :: TOL=1.D-5
      INTEGER(4),INTENT(IN)    :: NR1,NR2,NR3,NR1START,NR1L
      INTEGER(4),INTENT(IN)    :: NR
      INTEGER(4),INTENT(IN)    :: LMX
      REAL(8)   ,INTENT(INOUT) :: RHO(NR1L,NR2,NR3)
      REAL(8)   ,INTENT(IN)    :: DRHOL(NR,LMX)
      REAL(8)   ,INTENT(IN)    :: RBAS(3,3)
      REAL(8)   ,INTENT(IN)    :: R0(3)
      integer(4),INTENT(IN)    :: gid
      REAL(8)                  :: DR(3,3)
      REAL(8)                  :: YLM(LMXX)
      REAL(8)                  :: RVEC(3)
      REAL(8)                  :: RMAX,RMAX2,RI,SVAR,SVAR1
      INTEGER(4)               :: IR,LM
      INTEGER(4)               :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      REAL(8)                  :: T1,T2,T3
      REAL(8)                  :: X1,X2,X3
      INTEGER(4)               :: I1,I2,I3,I11,I21,I31,I
      REAL(8)                  :: DIS,DIS2
      REAL(8)                  :: XIR
      REAL(8)                  :: W1,W2
      REAL(8)                  :: r(nr)
!     ******************************************************************
      call radial$r(gid,nr,r)
      IF(LMXX.LT.LMX) THEN
        CALL ERROR$MSG('INCREASE DIMENSION LMXX')
        CALL ERROR$STOP('GRAPHICS_RHOLTOR')
      END IF
!
!     ==================================================================
!     ==  DETERMINE RMAX                                              ==
!     ==================================================================
      RMAX=0.D0
      DO IR=1,NR
        SVAR=0.D0
        DO LM=1,LMX
          SVAR=MAX(DABS(DRHOL(IR,LM)),SVAR)
        ENDDO
        IF(SVAR.GT.TOL)RMAX=R(ir)
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
                  CALL RADIAL$VALUE(gid,NR,DRHOL(1,LM),DIS,SVAR1)
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
      SUBROUTINE GRAPHICS_1CRHO(NR,LNX,LOX,AEPHI,PSPHI,LMNX &
     &                   ,DENMAT,LMX,DRHOL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: LMXX=36
      INTEGER(4),INTENT(IN)  :: NR
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      REAL(8)   ,INTENT(IN)  :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN)  :: PSPHI(NR,LNX)
      INTEGER(4),INTENT(IN)  :: LMNX
      REAL(8)   ,INTENT(IN)  :: DENMAT(LMNX,LMNX)
      INTEGER(4),INTENT(IN)  :: LMX
      REAL(8)   ,INTENT(OUT) :: DRHOL(NR,LMX)
      INTEGER(4)             :: LM,LMN1,LN1,LM1,L1,M1,LMN2,LN2,LM2,L2,M2
      REAL(8)                :: CG,DENMAT1,SVAR
      REAL(8)   ,ALLOCATABLE :: WORK(:)
!     ******************************************************************
      drhol(:,:)=0.d0
      ALLOCATE(WORK(NR))
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
!
          WORK(:)=AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2)
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
                  DRHOL(:,LM)=DRHOL(:,LM)+WORK(:)*SVAR
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
      SUBROUTINE GRAPHICS_1CWAVE(NR,LNX,LOX,AEPHI,PSPHI,LMNX &
     &                   ,PROJ,LMX,DRHOL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NR
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: LMX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      REAL(8)   ,INTENT(IN)  :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN)  :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(IN)  :: PROJ(LMNX)
      REAL(8)   ,INTENT(OUT) :: DRHOL(NR,LMX)
      INTEGER(4)             :: LM,IR,LMN,LN,M,L
      REAL(8)                :: SVAR
!     ******************************************************************
      drhol(:,:)=0.d0
      LMN=0
      DO LN=1,LNX
        L=LOX(LN)
        DO M=1,2*L+1
          LMN=LMN+1
          LM=L**2+M
          IF(LM.LE.LMX) THEN
            SVAR=PROJ(LMN)
!            PRINT*,'1CWAVE',LN,M,LMN,LM,SVAR
            DRHOL(:,LM)=DRHOL(:,LM)+SVAR*(AEPHI(:,LN)-PSPHI(:,LN))
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_REFINEGRID(NR1,NR2,NR3,NR1B,NR2B,NR3B,WAVE,WAVEBIG)
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: NR1
      INTEGER(4),INTENT(IN)     :: NR2
      INTEGER(4),INTENT(IN)     :: NR3
      INTEGER(4),INTENT(IN)     :: NR1B
      INTEGER(4),INTENT(IN)     :: NR2B
      INTEGER(4),INTENT(IN)     :: NR3B
      REAL(8),INTENT(IN)        :: WAVE(NR1,NR2,NR3)
      REAL(8),INTENT(OUT)       :: WAVEBIG(NR1B,NR2B,NR3B)
      COMPLEX(8),ALLOCATABLE    :: WORKC1(:,:,:)
      COMPLEX(8),ALLOCATABLE    :: WORKC2(:,:,:)
      INTEGER(4)                :: I
      INTEGER(4)                :: J      
      INTEGER(4)                :: K
      INTEGER(4)                :: E1 
      INTEGER(4)                :: E2   
!     ******************************************************************
      ALLOCATE(WORKC1(NR1,NR2,NR3))
      WORKC1(:,:,:)=CMPLX(WAVE(:,:,:),KIND=8)
      ALLOCATE(WORKC2(NR1,NR2,NR3)) 
      CALL LIB$3DFFTC8('RTOG',NR1,NR2,NR3,WORKC1,WORKC2)
      DEALLOCATE(WORKC1)
      ALLOCATE(WORKC1(NR1B,NR2B,NR3B))      
      WORKC1=(0.D0,0.D0)
      I=NR1/2  !MIND: REQUIRES THAT ONLY EVEN NUMBERS ARE USED (-> ASSURED BY LIB$FFTADJUSTGRD)
      J=NR2/2
      K=NR3/2
!
      WORKC1(1:I          ,1:J          ,1:K)          =WORKC2(1:I    ,1:J    ,1:K)
      WORKC1(NR1B-I+1:NR1B,1:J          ,1:K)          =WORKC2(I+1:2*I,1:J    ,1:K)
      WORKC1(1:I          ,NR2B-J+1:NR2B,1:K)          =WORKC2(1:I    ,J+1:2*J,1:K)        
      WORKC1(1:I          ,1:J          ,NR3B-K+1:NR3B)=WORKC2(1:I    ,1:J    ,K+1:2*K)
      WORKC1(NR1B-I+1:NR1B,NR2B-J+1:NR2B,1:K)          =WORKC2(I+1:2*I,J+1:2*J,1:K)
      WORKC1(1:I          ,NR2B-J+1:NR2B,NR3B-K+1:NR3B)=WORKC2(1:I    ,J+1:2*J,K+1:2*K)        
      WORKC1(NR1B-I+1:NR1B,1:J          ,NR3B-K+1:NR3B)=WORKC2(I+1:2*I,1:J    ,K+1:2*K)
      WORKC1(NR1B-I+1:NR1B,NR2B-J+1:NR2B,NR3B-K+1:NR3B)=WORKC2(I+1:2*I,J+1:2*J,K+1:2*K)
!
      DEALLOCATE(WORKC2)
      ALLOCATE(WORKC2(NR1B,NR2B,NR3B))
      CALL LIB$3DFFTC8('GTOR',NR1B,NR2B,NR3B,WORKC1,WORKC2)
      WAVEBIG(:,:,:)=REAL(WORKC2(:,:,:),KIND=8)
      DEALLOCATE(WORKC1)
      DEALLOCATE(WORKC2)
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
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: NGL_
      COMPLEX(8)  ,INTENT(IN) :: VHARTREE(NGL_)
!     ********************************************************************
      IF(.NOT.TINI) RETURN
      IF(.NOT.TWAKE) RETURN
      IF(.NOT.TPOT) RETURN
      IF(NGL.NE.0.AND.NGL.NE.NGL_) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$I4VAL('NGL',NGL)
        CALL ERROR$I4VAL('NGL_',NGL_)
        CALL ERROR$STOP('GRAPHICS$SETPWPOT')
      END IF
      NGL=NGL_
      IF(.NOT.ALLOCATED(PWPOT)) THEN
        ALLOCATE(PWPOT(NGL))
      END IF
      PWPOT(:)=VHARTREE(:)
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE GRAPHICS$SET1CPOT(IDENT_,IAT_,gid,NR,NRX_,LMRX,POT)
!     ********************************************************************
!     **  USE 1-CENTER POTENTIAL FOR ELECTRIC FIELD GRADIENTS          **
!     ********************************************************************
      USE GRAPHICS_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: IDENT_  ! CAN BE 'AE' OR 'PS' 
      INTEGER(4)   ,INTENT(IN) :: IAT_    ! ATOM INDEX (SEE ATOMLIST)
      INTEGER(4)   ,INTENT(IN) :: gid
      INTEGER(4)   ,INTENT(IN) :: NR,NRX_
      INTEGER(4)   ,INTENT(IN) :: LMRX
      REAL(8)      ,INTENT(IN) :: POT(NRX_,LMRX)
      INTEGER(4)               :: NAT
!     ********************************************************************
      IF(.NOT.TINI) RETURN
!     ===================================================================
!     == CHECK WHETHER ACTION IS REQUIRED                              ==
!     ===================================================================
      IF(.NOT.TWAKE) RETURN
      IF(.NOT.TPOT) RETURN
      IF(.NOT.ALLOCATED(AE1CPOT)) THEN
        NRX=NRX_
        CALL SETUP$GETI4('LMRXX',LMRXX)
        CALL ATOMLIST$NATOM(NAT)
        ALLOCATE(AE1CPOT(NRX,LMRXX,NAT))
        ALLOCATE(PS1CPOT(NRX,LMRXX,NAT))
        AE1CPOT=0.D0
        PS1CPOT=0.D0
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
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRAPHICS_CREATEPOT(FILE,TITLE,DR)
!     ******************************************************************
      USE GRAPHICS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: TITLE
      CHARACTER(*),INTENT(IN)    :: FILE
      REAL(8)     ,INTENT(IN)    :: DR
      INTEGER(4)                 :: FACT
      REAL(8)                    :: RBAS(3,3),GBAS(3,3),DET
      REAL(8)   ,ALLOCATABLE     :: VHARTREE(:)
      REAL(8)   ,ALLOCATABLE     :: WORK(:)
      REAL(8)   ,ALLOCATABLE     :: POTENTIAL(:,:,:)
      INTEGER(4)                 :: NRL
      INTEGER(4)                 :: NR1L
      INTEGER(4)                 :: NR1,NR2,NR3
      CHARACTER(32)              :: STRING
      INTEGER(4)                 :: NFIL
      INTEGER(4)                 :: NAT
      INTEGER(4)                 :: IAT
      CHARACTER(32),ALLOCATABLE  :: ATOMNAME(:)
      REAL(8)   ,ALLOCATABLE     :: Q(:)
      REAL(8)   ,ALLOCATABLE     :: Z(:)
      REAL(8)   ,ALLOCATABLE     :: POS(:,:)
      INTEGER(4)                 :: NR
      INTEGER(4)                 :: LMRX
      REAL(8)                    :: DEX
      REAL(8)                    :: R1
      REAL(8)   ,ALLOCATABLE     :: AEPOT(:,:)
      REAL(8)   ,ALLOCATABLE     :: PSPOT(:,:)
      INTEGER(4)                 :: NFILO
      REAL(8)   ,ALLOCATABLE     :: ONECPOT(:,:,:)
      INTEGER(4)                 :: NTASKS,THISTASK
      INTEGER(4)                 :: NR1B,NR2B,NR3B
      INTEGER(4)                 :: isp
      INTEGER(4)                 :: gid
!     ******************************************************************
!COLLECTING OF INFORMATION
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$GETI4('NRL',NRL)
      CALL PLANEWAVE$GETI4('NR1L',NR1L)
      CALL PLANEWAVE$GETI4('NR1',NR1)
      CALL PLANEWAVE$GETI4('NR2',NR2)
      CALL PLANEWAVE$GETI4('NR3',NR3)
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL GBASS(RBAS,GBAS,DET)
      FACT=NINT(((DET/DR**3)/REAL(NR1*NR2*NR3,KIND=8))**(1./3.))
      FACT=MAX(1,FACT)
!
!     ==================================================================
!     ==  EXPAND GRID BY FACT AND UNPARALLELIZE                       ==
!     ==================================================================
      NR1B=FACT*NR1
      NR2B=FACT*NR2
      NR3B=FACT*NR3
!     == TAKE VALUES COMPATIBLE WITH FFT ROUTINES ====================
      CALL LIB$FFTADJUSTGRD(NR1B) 
      CALL LIB$FFTADJUSTGRD(NR2B) 
      CALL LIB$FFTADJUSTGRD(NR3B) 
      ALLOCATE(POTENTIAL(NR1B,NR2B,NR3B))
      ALLOCATE(VHARTREE(NRL))
      CALL PLANEWAVE$SUPFFT('GTOR',1,NGL,PWPOT,NRL,VHARTREE)
      IF(FACT.EQ.1) THEN
        CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,VHARTREE,NR1*NR2*NR3,POTENTIAL)
      ELSE
        ALLOCATE(WORK(NR1*NR2*NR3))
        CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,VHARTREE,NR1*NR2*NR3,WORK)
        CALL GRAPHICS_REFINEGRID(NR1,NR2,NR3,NR1B,NR2B,NR3B,WORK,POTENTIAL)
        DEALLOCATE(WORK)
      END IF
      DEALLOCATE(VHARTREE)
!
!     ==================================================================
!     ==  ONE-CENTER CONTRIBUTIONS                                    ==
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
      ONECPOT=AE1CPOT-PS1CPOT
      CALL MPE$COMBINE('MONOMER','+',ONECPOT)
      DO IAT=1,NAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        call setup$geti4('gid',gid)                
        call radial$geti4(gid,'nr',nr)
        IF(NRX.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT GRID SIZE')
          CALL ERROR$MSG('ERROR ENTERED WHILE ALLOWING ATOM SPECIFIC RADIAL GRIDS')
          CALL ERROR$STOP('GRAPHICS_CREATEPOT')
        END IF
!        CALL SETUP$RADGRID(IAT,R1,DEX,NR)
        CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME(IAT))
        CALL ATOMLIST$GETR8A('R(0)',IAT,3,POS(:,IAT))
        CALL ATOMLIST$GETR8('Z',IAT,Z(IAT))
        CALL ATOMLIST$GETR8('Q',IAT,Q(IAT))
PRINT*,'INCLUDE AE-CONTRIBUTIONS'
CALL TIMING$CLOCKON('GRAPHICS 1CPOTENTIAL')
         CALL GRAPHICS_RHOLTOR(RBAS,NR1B,NR2B,NR3B,1,NR1B &
      &           ,POTENTIAL,POS(:,IAT),gid,NRX,LMRXX,ONECPOT(:,:,IAT))
CALL TIMING$CLOCKOFF('GRAPHICS 1CPOTENTIAL')
PRINT*,'INCLUDED AE-CONTRIBUTIONS'
      ENDDO
      DEALLOCATE(ONECPOT)
!
!     ==================================================================
!     ==  WRITE TO FILE                                               ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$SETFILE('WAVEPLOT',.FALSE.,TRIM(FILE))
        CALL FILEHANDLER$SETSPECIFICATION('WAVEPLOT','FORM','UNFORMATTED')
        CALL FILEHANDLER$UNIT('WAVEPLOT',NFIL)
        CALL WRITEWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,Q,ATOMNAME &
     &                  ,NR1B,NR2B,NR3B,POTENTIAL)    
        CALL FILEHANDLER$CLOSE('WAVEPLOT')
      END IF
      DEALLOCATE(POTENTIAL)
      RETURN
      END
