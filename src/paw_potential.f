      MODULE POTENTIAL_MODULE
      LOGICAL(4) :: TINI=.FALSE.  ! PLANE WAVE GRID DEFINED OR NOT
      LOGICAL(4) :: TSET=.FALSE.  ! MODULE DATA DEFINED OR NOT
      REAL(8)    :: EPWRHO     ! PLANE WAVE CUTOFF FOR THE DENSITY
      INTEGER(4) :: NR1GLOB    ! #(R-POINTS IN FIRST DIRECTION; ALL TASK)
      INTEGER(4) :: NR1L       ! #(R-POINTS IN FIRST DIRECTION; THIS TASK)
      INTEGER(4) :: NR1START   ! 
      INTEGER(4) :: NR2,NR3    ! #(R-POINTS IN SECOND AND THIRD DIRECTION)
      INTEGER(4) :: NGL        ! #(G-POINTS; FOR THIS TASK)
      INTEGER(4) :: NSP        ! #(ELEMENTS)
      INTEGER(4) :: LMRXX      ! MAX #(ANGULAR MOMENTA IN THE 1-CENTER MOMENTS)
      REAL(8)    ,ALLOCATABLE :: V0(:,:)      !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: DV0(:,:)     !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: G0(:,:)      !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: DG0(:,:)     !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: PSCOREG(:,:) !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: DPSCOREG(:,:)!(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: VBARG(:,:)   !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: DVBARG(:,:)  !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: YLMOFG(:,:)  !(LMRXX,NGL)
      INTEGER(4) ,ALLOCATABLE :: LMRX(:)      !(NSP)
      LOGICAL(4)              :: TSTRESS=.FALSE.
      LOGICAL(4)              :: TFORCE=.FALSE.
      END MODULE POTENTIAL_MODULE
!
!     ..................................................................
      SUBROUTINE POTENTIAL$GETR8(ID,VALUE)
!     ******************************************************************
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VALUE
!     ******************************************************************
      IF(ID.EQ.'EPWRHO') THEN
        VALUE=EPWRHO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('POTENTIAL$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL$SETL4(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(in) :: VAL
!     ******************************************************************
      IF(ID.EQ.'STRESS') THEN
        TSTRESS=VAL
      ELSE IF(ID.EQ.'FORCE') THEN
        TFORCE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('POTENTIAL$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL$SETR8A(ID,LEN,VAL)
!     ******************************************************************
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
      INTEGER(4)              :: I,J,IND
!     ******************************************************************
      IF(ID.EQ.'RBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('SIZE INCONSISTENCY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('POTENTIAL$SETR8A')
        END IF
        CALL POTENTIAL_UPDATE(VAL)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('POTENTIAL$SETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL$REPORT(NFIL)
!     ******************************************************************
!     **                                                              **
!     **  REPORT ON THE POTENTIAL OBJECT                              **
!     **                                                              **
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      INTEGER(4)    :: NFIL
      REAL(8)       :: MBYTE
      REAL(8)       :: RY
      REAL(8)       :: MEMORY
      INTEGER(4)    :: NGG
!     ******************************************************************
      CALL CONSTANTS('RY',RY)
      CALL CONSTANTS('MBYTE',MBYTE)
      CALL REPORT$TITLE(NFIL,'POTENTIAL')
      IF(.NOT.TINI) THEN
        CALL REPORT$STRING(NFIL,'POTENTIAL OBJECT NOT YET INITIALIZED')
        RETURN
      END IF
      CALL REPORT$R8VAL(NFIL,'DENSITY PLANE WAVE CUTOFF',EPWRHO/RY,'RY')
      CALL PLANEWAVE$SELECT('DENSITY')      
      CALL PLANEWAVE$GETI4('NGG',NGG)      
      CALL REPORT$I4VAL(NFIL,'#(G-VECTORS FOR DENSITY)',NGG,' ') 
      MEMORY=REAL(8*(NGL*NSP)*4   & ! V0,G0,PSCOREG,VBARG
     &           +8*(NGL*LMRXX)*1 & ! YLMOFG
     &           +4*(NSP)*1       & ! LMRX
     &           ,KIND=8)/MBYTE
      CALL REPORT$R8VAL(NFIL,'MEMORY REQUIRED BIG PERMANENT ARRAYS' &
     &           ,MEMORY/MBYTE,'MBYTE')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL$CLEAR
!     ******************************************************************
!     **  CLEAR PERMANENT MEMORY OF THE POTENTIAL OBJECT              **
!     **  (WILL AUTOMATICALLY BE RECREATED WHEN REQUIRED)             **
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
!     ******************************************************************
                                CALL TRACE$PUSH('POTENTIAL$CLEAR')
      IF(ALLOCATED(YLMOFG)) DEALLOCATE(YLMOFG)
      IF(ALLOCATED(PSCOREG))DEALLOCATE(PSCOREG)
      IF(ALLOCATED(VBARG))  DEALLOCATE(VBARG)
      IF(ALLOCATED(G0))     DEALLOCATE(G0)
      IF(ALLOCATED(V0))     DEALLOCATE(V0)
      IF(ALLOCATED(LMRX))   DEALLOCATE(LMRX)
                                CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL$INITIALIZE(EPWRHO_,NR1START_,NR1L_,NR2_,NR3_)
!     ******************************************************************
!     **  DEFINE GRID FOR FOURIER TRANSFORM                           **
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: EPWRHO_
      INTEGER(4),INTENT(IN) :: NR1START_
      INTEGER(4),INTENT(IN) :: NR1L_
      INTEGER(4),INTENT(IN) :: NR2_
      INTEGER(4),INTENT(IN) :: NR3_
      REAL(8)               :: RBAS(3,3)     ! LATTICE VECTORS
      REAL(8)               :: CELLVOL       ! UNIT CELL VOLUME
      REAL(8)   ,PARAMETER  :: K(3)=(/0.D0,0.D0,0.D0/)
      INTEGER(4)            :: ISP
!     ******************************************************************
                                CALL TRACE$PUSH('POTENTIAL$INITIALIZE')
      IF(TINI) THEN
        CALL ERROR$MSG('POTENTIAL OBJECT IS ALREADY INITIALIZED')
        CALL ERROR$STOP('POTENTIAL$INITIALIZE')
      END IF
      TINI=.TRUE.
      TSET=.FALSE.
!     
!     ================================================================
!     ==  MAP INPUT DATA ONTO MODULE                                ==
!     ================================================================
      EPWRHO=EPWRHO_
      NR1START=NR1START_
      NR1L=NR1L_
      NR2=NR2_
      NR3=NR3_
!     
!     ================================================================
!     ==  INITIALIZE PLANEWAVE OBJECT                               ==
!     ================================================================
      CALL CELL$GETR8A('TREF',9,RBAS)
      CALL PLANEWAVE$INITIALIZE('DENSITY',RBAS,K,.TRUE.,EPWRHO &
     &                         ,NR1START,NR1L,NR2,NR3)
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$SETL4('SUPER',.TRUE.)
      CALL PLANEWAVE$GETI4('NGL',NGL)
      CALL SETUP$LMRXX(LMRXX)
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(LMRX(NSP))
      DO ISP=1,NSP
        CALL SETUP$LMRX(ISP,LMRX(ISP))
      ENDDO
      CALL POTENTIAL_UPDATE(RBAS)
                                CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_RESET
!     ******************************************************************
!     **   RESET PERMANENT DATA OF POTENTIAL_MODULE                   **
!     **   USED TO UPDATE LATTICE VECTORS                             **
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
!     ******************************************************************
      TSET=.FALSE.
      DEALLOCATE(YLMOFG)
      DEALLOCATE(PSCOREG)
      DEALLOCATE(VBARG)
      DEALLOCATE(G0)
      DEALLOCATE(V0)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_UPDATE(RBAS)
!     ******************************************************************
!     **  COLLECTS INTERNAL ARRAYS                                    **
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,ALLOCATABLE:: GVEC(:,:)
      REAL(8)   ,ALLOCATABLE:: G2(:)
      INTEGER(4)            :: ISP
      REAL(8)               :: GBAS(3,3)
      REAL(8)               :: CELLVOL
      CHARACTER(64)         :: PLANEWAVESELECTION
!     ******************************************************************
                                CALL TRACE$PUSH('POTENTIAL_UPDATE')
      CALL PLANEWAVE$GETCH('ID',PLANEWAVESELECTION)
!
!     ================================================================
!     == ALLOCATE ARRAYS                                            ==
!     ================================================================
      IF(.NOT.ALLOCATED(YLMOFG))  ALLOCATE(YLMOFG(LMRXX,NGL))           
      IF(.NOT.ALLOCATED(PSCOREG)) ALLOCATE(PSCOREG(NGL,NSP))            
      IF(.NOT.ALLOCATED(DPSCOREG))ALLOCATE(DPSCOREG(NGL,NSP))           
      IF(.NOT.ALLOCATED(VBARG))   ALLOCATE(VBARG(NGL,NSP))              
      IF(.NOT.ALLOCATED(DVBARG))  ALLOCATE(DVBARG(NGL,NSP))             
      IF(.NOT.ALLOCATED(G0))      ALLOCATE(G0(NGL,NSP)) ;G0(:,:)=0.D0   
      IF(.NOT.ALLOCATED(DG0))     ALLOCATE(DG0(NGL,NSP));DG0(:,:)=0.D0 
      IF(.NOT.ALLOCATED(V0))      ALLOCATE(V0(NGL,NSP))                 
      IF(.NOT.ALLOCATED(DV0))     ALLOCATE(DV0(NGL,NSP))                
!
!     ================================================================
!     == RESET LATTICE CONSTANT OF PLANEWAVE MODULE                 ==
!     ================================================================
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
!
!     ================================================================
!     == CALCULATE YLM(G)*G**L FOR VOFRHO                           ==
!     ================================================================
      CALL GBASS(RBAS,GBAS,CELLVOL)
      ALLOCATE(GVEC(3,NGL))
      ALLOCATE(G2(NGL))
      CALL PLANEWAVE$GETR8A('G2',NGL,G2)
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
      CALL POTENTIAL_YLMOFG(LMRXX,LMRXX,NGL,GVEC,YLMOFG)
!
!     ================================================================
!     == OBTAIN ATOM SPECIFIC ARRAYS FROM SETUP                     ==
!     ================================================================
      DO ISP=1,NSP
          CALL SETUP$ISELECT(ISP)
          CALL SETUP$GETFOFG('G0',.FALSE.,1,NGL,G2,CELLVOL,G0(1,ISP))
          CALL SETUP$GETFOFG('G0',.TRUE. ,1,NGL,G2,CELLVOL,DG0(1,ISP))
          CALL SETUP$GETFOFG('V0',.FALSE.,1,NGL,G2,CELLVOL,V0(1,ISP))
          CALL SETUP$GETFOFG('V0',.TRUE. ,1,NGL,G2,CELLVOL,DV0(1,ISP))
          CALL SETUP$GETFOFG('PSCORE',.FALSE.,1,NGL,G2,CELLVOL,PSCOREG(1,ISP))
          CALL SETUP$GETFOFG('PSCORE',.TRUE. ,1,NGL,G2,CELLVOL,DPSCOREG(1,ISP))
          CALL SETUP$GETFOFG('VADD',.FALSE.,1,NGL,G2,CELLVOL,VBARG(1,ISP))
          CALL SETUP$GETFOFG('VADD',.TRUE. ,1,NGL,G2,CELLVOL,DVBARG(1,ISP))
!
!        CALL SETUP$COMPENSATION(ISP,NGL,G2,CELLVOL &
!     &                      ,G0(1,ISP),V0(1,ISP),DG0(1,ISP),DV0(1,ISP))
!        CALL SETUP$PSCOREG(ISP,NGL,G2,CELLVOL,PSCOREG(1,ISP),DPSCOREG(1,ISP))
!        CALL SETUP$VBARG(ISP,NGL,G2,CELLVOL,VBARG(1,ISP),DVBARG(1,ISP))
      ENDDO
      DEALLOCATE(G2)
      DEALLOCATE(GVEC)
      CALL PLANEWAVE$SELECT(PLANEWAVESELECTION)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL$VOFRHO(NRL,NDIMD,RHO,LMRXX_,NAT_,QLM,VQLM &
     &                           ,R0,FORCE,RBAS,STRESS,RHOB)
!     ******************************************************************
!     **  CALCULATE POTENTIAL                                         **
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NRL
      INTEGER(4),INTENT(IN)    :: NDIMD
      INTEGER(4),INTENT(IN)    :: NAT_
      INTEGER(4),INTENT(IN)    :: LMRXX_
      REAL(8)   ,INTENT(IN)    :: QLM(LMRXX_,NAT_)
      REAL(8)   ,INTENT(OUT)   :: VQLM(LMRXX_,NAT_)
      REAL(8)   ,INTENT(IN)    :: R0(3,NAT_)
      REAL(8)   ,INTENT(OUT)   :: FORCE(3,NAT_)
      REAL(8)   ,INTENT(IN)    :: RBAS(3,3)
      REAL(8)   ,INTENT(OUT)   :: STRESS(3,3)
      REAL(8)   ,INTENT(OUT)   :: RHOB
      REAL(8)   ,INTENT(INOUT) :: RHO(NRL,NDIMD)
      INTEGER(4),ALLOCATABLE   :: ISPECIES(:) !(NAT) 
      REAL(8)   ,ALLOCATABLE   :: G2(:)       !(NGL)
      REAL(8)   ,ALLOCATABLE   :: GVEC(:,:)     !(3,NGL)
      INTEGER(4)               :: IAT,IR
      INTEGER(4)               :: NAT
      INTEGER(4)               :: NSPIN
      INTEGER(4)               :: ISVAR
      LOGICAL(4)               :: TGRA
      REAL(8)   ,ALLOCATABLE   :: RHOTEMP(:,:)
      REAL(8)                  :: SVAR
!     ******************************************************************
                                CALL TRACE$PUSH('POTENTIAL$VOFRHO')
                                CALL TIMING$CLOCKON('POTENTIAL')
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('POTENTIAL OBJECT IS NOT INITIALIZED')
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
      IF(LMRXX_.NE.LMRXX) THEN
        CALL ERROR$MSG('INPUT DATA INCONSISTENT WITH MODULE')
        CALL ERROR$I4VAL('LMRXX_',LMRXX_)
        CALL ERROR$I4VAL('LMRXX',LMRXX)
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
!     == COLLECT FROM ATOMLIST =======================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      IF(NAT_.NE.NAT) THEN
        CALL ERROR$MSG('INPUT DATA INCONSISTENT WITH MODULE')
        CALL ERROR$I4VAL('NAT_',NAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
  
!     
!     == PLOT DENSITY ================================================
!     IF(TDENSITYTODX(1).AND.TDENSITYTODX(2)) THEN
!       CALL FILEHANDLER$UNIT('DENSITY.DX',NFIL)
!       CALL  WRITEDX_DENSITY
!    &       (NFIL,NR1GLOB,NR2,NR3,NR1GLOB,NR2,NR3,RBAS,RHOE)
!       CALL FILEHANDLER$CLOSE('DENSITY.DX')
!       TDENSITYTODX(2)=.FALSE.
!     END IF
!     
!     == CONVERT DENSITY INTO POTENTIAL ==============================
      CALL SETUP$NSPECIES(NSP)

!     == COLLECT FROM PLANEWAVE OBJECT =================================
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$GETI4('NGL',NGL)
      ALLOCATE(G2(NGL))
      ALLOCATE(GVEC(3,NGL))
      CALL PLANEWAVE$GETR8A('G2',NGL,G2)
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
      CALL PLANEWAVE$GETI4('NRL',ISVAR)
      IF(ISVAR.NE.NRL) THEN
        CALL ERROR$MSG('#(GRID POINTS INCONSISTENT)')
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
      CALL PLANEWAVE$GETI4('NR1',NR1GLOB)
!     == COLLECT FROM DFT OBJECT =======================================
      CALL DFT$GETL4('GC',TGRA)
!
!     ==================================================================
!     == POTENTIAL; BRANCH INTO COLLINEAR AND NON-COLLINEAR CASE      ==
!     ==================================================================
      IF(NDIMD.EQ.4) THEN
        NSPIN=2
        ALLOCATE(RHOTEMP(NRL,2))
        DO IR=1,NRL
          RHOTEMP(IR,1)=RHO(IR,1)
          RHOTEMP(IR,2)=SQRT(RHO(IR,2)**2+RHO(IR,3)**2+RHO(IR,4)**2)
        ENDDO
        CALL POTENTIAL_VOFRHO(LMRXX,NRL,NSP,NAT,ISPECIES,R0,FORCE &
     &                        ,NR1GLOB,NR1L,NR2,NR3,RHOTEMP,NSPIN,RBAS &
     &       ,PSCOREG,DPSCOREG,VBARG,DVBARG,YLMOFG,G0,V0,QLM,VQLM,LMRX &
     &                        ,NGL,GVEC,G2,RHOB,TSTRESS,DG0,DV0,STRESS)
        DO IR=1,NRL
          RHO(IR,1)=RHOTEMP(IR,1)
          SVAR=SQRT(RHO(IR,2)**2+RHO(IR,3)**2+RHO(IR,4)**2)
          SVAR=RHOTEMP(IR,2)/(SVAR+1.D-300)
          RHO(IR,2)=SVAR*RHO(IR,2)
          RHO(IR,3)=SVAR*RHO(IR,3)
          RHO(IR,4)=SVAR*RHO(IR,4)
        ENDDO
        DEALLOCATE(RHOTEMP)
      ELSE
        NSPIN=NDIMD
        CALL POTENTIAL_VOFRHO(LMRXX,NRL,NSP,NAT,ISPECIES,R0,FORCE &
     &                            ,NR1GLOB,NR1L,NR2,NR3,RHO,NSPIN,RBAS &
     &       ,PSCOREG,DPSCOREG,VBARG,DVBARG,YLMOFG,G0,V0,QLM,VQLM,LMRX &
     &                        ,NGL,GVEC,G2,RHOB,TSTRESS,DG0,DV0,STRESS)
      END IF
!
IF(TSTRESS) THEN
  WRITE(*,FMT='("C+XC STRESS ",3F15.7)')STRESS(1,:)
  WRITE(*,FMT='("C+XC STRESS ",3F15.7)')STRESS(2,:)
  WRITE(*,FMT='("C+XC STRESS ",3F15.7)')STRESS(3,:)
END IF
!
!     ==================================================================
!     == CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(ISPECIES)
      DEALLOCATE(G2)
      DEALLOCATE(GVEC)
                                CALL TIMING$CLOCKOFF('POTENTIAL')
                                CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_VOFRHO(LMRXX,NRL &
     &                    ,NSP,NAT,ISPECIES,TAU0,FION &
     &                    ,NR1GLOB,NR1,NR2,NR3,RHOE,NDIMD,RBAS &
     &         ,PSCORG,DPSCORG,VBARG,DVBARG,YLMOFG,G0,V0,QLM,VQLM,LMRX &
     &                    ,NGL,GVEC,G2,RHOB,TSTRESS,DG0,DV0,STRESS)
!     ******************************************************************
!     **                                                              **
!     **  THE PLANE WAVE CONTRIBUTION TO TOTAL ENERGY AND FORCES      **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    RHOE           CHARGE DENSITY IN REAL SPACE               **
!     **                   (SPIN UP AND SPIN DOWN                     **
!     **                                                              **
!     **  OUTPUT:                                                     **
!     **    RHOE           PLANE WAVE POTENTIAL IN REAL SPACE         **
!     **                   (SPIN UP AND SPIN DOWN)                    **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     ******************************************************************
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4)              :: TSPIN
      LOGICAL(4)              :: TCHK
!     == DUMMY   ARRAYS ================================================
      INTEGER(4),INTENT(IN)   :: LMRXX
      INTEGER(4),INTENT(IN)   :: NRL
      INTEGER(4),INTENT(IN)   :: NSP
      INTEGER(4),INTENT(IN)   :: NAT
      INTEGER(4),INTENT(IN)   :: NGL
      INTEGER(4),INTENT(IN)   :: NDIMD
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT)
      REAL(8)   ,INTENT(INOUT):: TAU0(3,NAT)     !<-
      REAL(8)   ,INTENT(OUT)  :: FION(3,NAT)
      INTEGER(4),INTENT(IN)   :: NR1GLOB
      INTEGER(4),INTENT(IN)   :: NR1
      INTEGER(4),INTENT(IN)   :: NR2
      INTEGER(4),INTENT(IN)   :: NR3
      REAL(8)   ,INTENT(INOUT):: RHOE(NRL,NDIMD)
      REAL(8)   ,INTENT(IN)   :: RBAS(3,3)
      REAL(8)   ,INTENT(IN)   :: PSCORG(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: DPSCORG(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: VBARG(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: DVBARG(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: YLMOFG(LMRXX,NGL)
      REAL(8)   ,INTENT(IN)   :: G0(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: DG0(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: V0(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: DV0(NGL,NSP)
      REAL(8)   ,INTENT(INOUT):: QLM(LMRXX,NAT)   !<-
      REAL(8)   ,INTENT(OUT)  :: VQLM(LMRXX,NAT)
      INTEGER(4),INTENT(IN)   :: LMRX(NSP)
      LOGICAL(4),INTENT(IN)   :: TSTRESS
      REAL(8)   ,INTENT(OUT)  :: STRESS(3,3)
      REAL(8)   ,INTENT(OUT)  :: RHOB
      LOGICAL(4)              :: TGRA
      INTEGER(4)              :: NSPIN
      REAL(8)                 :: GVEC(3,NGL)
      REAL(8)                 :: G2(NGL)
      REAL(8)   ,ALLOCATABLE  :: GRHO(:,:,:) !(NNR,3,NSPIN)
      COMPLEX(8),ALLOCATABLE  :: VHARTREE(:) !(NGL)
      COMPLEX(8),ALLOCATABLE  :: RHOG(:,:)   !(NGL,NSPIN)
      COMPLEX(8),ALLOCATABLE  :: CWORK(:,:)  !(NGL,NSPIN)
      REAL(8)                 :: RCSM(NSP)   ! COMPENSATION GAUSSIAN DECAY
      REAL(8)                 :: RCBG(NSP)   ! LARGE GAUSSIAN DECAY
      REAL(8)                 :: GBAS(3,3)   ! RECIPROCAL LATTICE VECTORS
      REAL(8)                 :: FIONT(3,NAT)      
      REAL(8)                 :: VQLMT(LMRXX,NAT)  
      REAL(8)                 :: FORCE1(3,NAT)      
      REAL(8)                 :: VQLM1(LMRXX,NAT)  
      REAL(8)                 :: EHARTREE
      REAL(8)                 :: EPAIR
      REAL(8)                 :: EISOLATE
      REAL(8)                 :: EXC
      REAL(8)                 :: RHOBT
      COMPLEX(8)              :: RHOTOT,RHOSPN,RHOUP,RHODWN
      REAL(8)                 :: CELLVOL
      REAL(8)                 :: SVAR
      INTEGER(4)              :: NNR,IR,LRXX,IAT,ISPIN,I,LM
      INTEGER(4)              :: IC,IG,ISP
      REAL(8)                 :: STRESST(3,3),STRESS1(3,3)
!     ==  VARIABLES FOR SELF TEST ======================================
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=2
      COMPLEX(8),ALLOCATABLE  :: RHO_SELFTEST(:,:)
      REAL(8)   ,ALLOCATABLE  :: RHELP(:)
      LOGICAL(4)              :: TBACK
      INTEGER(4)              :: NGAMMA
!     ******************************************************************
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL GBASS(RBAS,GBAS,CELLVOL)
      CALL DFT$GETL4('GC',TGRA)
      NNR=NR1*NR2*NR3
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      TSPIN=(NSPIN.EQ.2)
      LRXX=INT(SQRT(REAL(LMRXX-1)+1.E-5))
!
!     ==================================================================
!     == TRANSFORM PSEUDO DENSITY INTO FOURIER SPACE                  ==
!     ==================================================================
!CALL PLANEWAVE$GETI4('NGAMMA',NGAMMA)
!svar=0.d0
!do ir=1,nrl
!  svar=svar+rhoe(ir,1)
!enddo
!svar=svar*cellvol/real(nrl,kind=8)
!print*,'charge before transform ',svar
!print*,'nspin,ngl,nrl,cellvol,ngamma ',nspin,ngl,nrl,cellvol,ngamma

      ALLOCATE(RHOG(NGL,NSPIN))
      CALL PLANEWAVE$SUPFFT('RTOG',NSPIN,NGL,RHOG,NRL,RHOE)

!print*,'rhog(1) a ',rhog(ngamma,1)*cellvol
!stop
!
!     ==================================================================
!     ==  ADD PSEUDO CORE DENSITY                                     ==
!     ==================================================================
      CALL POTENTIAL_ADDCORE(NGL,NSP,NAT,ISPECIES,TAU0,PSCORG,RHOG(1,1))
!
!     ==================================================================
!     == SELFTEST                                                     ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(RHO_SELFTEST(NGL,NSPIN))
!       == USE THE FOLLOWING LINE TO AVOID PROBLEMS WITH TEH LOW-DENSITY
!       == CUTOFF OF THE XC ENERGY
        RHOG(1,:)=RHOG(1,:)+(1.D-1,0.D0) 
        RHO_SELFTEST(:,:)=RHOG(:,:)
      END IF
 1000 CONTINUE
      IF(TTEST) THEN
        IF(ITEST.EQ.1) THEN
          ALLOCATE(RHELP(2*NGL*NSPIN))
          RHELP=TRANSFER(RHOG,RHELP)
          CALL SELFTEST$START('POTENTIAL-RHO',2*NGL*NSPIN,RHELP,1.D-5)
          RHOG(:,:)=RESHAPE(TRANSFER(RHELP,RHOG),(/NGL,NSPIN/))
          DEALLOCATE(RHELP)
          EHARTREE=0.D0
          EXC=0.D0
          EISOLATE=0.D0
        ELSE IF(ITEST.EQ.2) THEN
!         == REDEFINE QLM AS INTENT(INOUT)
          RHOG(:,:)=RHO_SELFTEST(:,:)
          CALL SELFTEST$START('POTENTIAL-QLM',LMRXX*NAT,QLM,1.D-1)
          VQLM(:,:)=0.D0
        ELSE IF(ITEST.EQ.3) THEN
!         == REDEFINE TAU0 AS INTENT(INOUT)
          RHOG(:,:)=RHO_SELFTEST(:,:)
          CALL SELFTEST$START('POTENTIAL-FORCE',3*NAT,TAU0,1.D-4)
          FION(:,:)=0.D0
        END IF
      END IF
!
!     ==================================================================
!     ==  RESET ARRAYS FOR PARALLEL PROCESSING                        ==
!     ==================================================================
      FION(:,:)=0.D0
      VQLM(:,:)=0.D0
      STRESS(:,:)=0.D0
      FIONT(:,:)=0.D0
      VQLMT(:,:)=0.D0
      STRESST(:,:)=0.D0
!
!     ==================================================================
!     ==  EVALUATE HARTREE ENERGY IN G-SPACE                          ==
!     ==================================================================
                                CALL TIMING$CLOCKON('VOFRHO: HARTREE')
      EHARTREE=0.D0
      FORCE1(:,:)=0.D0
      VQLM1(:,:)=0.D0
      ALLOCATE(VHARTREE(NGL))
      VHARTREE(:)=(0.D0,0.D0)
      CALL POTENTIAL_HARTREE(LMRXX,LRXX,NGL,NSP,NAT,ISPECIES,LMRX &
     &            ,RHOG(1,1),VHARTREE,QLM,VQLM1,TAU0,FORCE1,EHARTREE &
     &            ,CELLVOL,G2,GVEC,YLMOFG,VBARG,DVBARG,G0,V0,RHOBT &
     &            ,TSTRESS,DG0,DV0,STRESS1)
      CALL MPE$COMBINE('+',RHOBT)
      FIONT(:,:)=FIONT(:,:)+FORCE1(:,:)
      VQLMT(:,:)=VQLMT(:,:)+VQLM1(:,:)
      STRESST(:,:)=STRESST(:,:)+STRESS1(:,:)
                                CALL TIMING$CLOCKOFF('VOFRHO: HARTREE')
!
!     ==================================================================
!     ==  CALCULATE PAIR-POTENTIAL                                    ==
!     ==================================================================
                                CALL TIMING$CLOCKON('VOFRHO: PAIRP')
      DO ISP=1,NSP
        CALL SETUP$RCSM(ISP,RCSM(ISP))
        CALL SETUP$RCBG(ISP,RCBG(ISP))
      ENDDO
      EPAIR=0.D0
      FORCE1(:,:)=0.D0
      VQLM1(:,:)=0.D0
      CALL POTENTIAL_PAIRP(RBAS,NSP,NAT,ISPECIES,LMRXX,LMRX,RCSM,RCBG &
     &                          ,TAU0,FORCE1,QLM,VQLM1,EPAIR,STRESS1)
      FION(:,:)=FION(:,:)+FORCE1(:,:)
      VQLM(:,:)=VQLM(:,:)+VQLM1(:,:)
      STRESS(:,:)=STRESS(:,:)+STRESS1(:,:)
                                CALL TIMING$CLOCKOFF('VOFRHO: PAIRP')
PRINT*,'PAIRP. ENERGY ',EPAIR
WRITE(*,FMT='("PAIRP STRESS ",3F15.7)')STRESS1(1,:)
WRITE(*,FMT='("PAIRP STRESS ",3F15.7)')STRESS1(2,:)
WRITE(*,FMT='("PAIRP STRESS ",3F15.7)')STRESS1(3,:)
!
!     ==================================================================
!     ==  ISOLATE CHARGE DENSITY FROM PERIODIC IMAGES                 ==
!     ==================================================================
!     CALL DIPOLE(RBAS,NGL,NGL,GVEC,RHOG(1,1),NAX,NSP,NA,TAU0,LMRXX,LMRX,QLM)
                                CALL TIMING$CLOCKON('VOFRHO: ISOLATE')
      EISOLATE=0.D0
      FORCE1(:,:)=0.D0
      VQLM1(:,:)=0.D0
      CALL ISOLATE(NSP,NAT,ISPECIES,TAU0,RBAS,FORCE1,EISOLATE &
     &            ,LMRXX,LMRX,QLM,VQLM1,RHOBT,NGL,RHOG(1,1),VHARTREE)
      FION(:,:)=FION(:,:)+FORCE1(:,:)
      VQLM(:,:)=VQLM(:,:)+VQLM1(:,:)
                                CALL TIMING$CLOCKOFF('VOFRHO: ISOLATE')
!
!     ==================================================================
!     ==  HYPERFINE PARAMETERS                                        ==
!     ==================================================================
!     == CHECK IF NTILDE+NTILDEC (LIKE HERE) OR ONLY NTILDE SHALL BE USED      
      CALL HYPERFINE$SETPWRHO('TOT',NGL,GVEC,RHOG(1,1))
      CALL HYPERFINE$GETFLAG('SPIN',TCHK)
      IF(TSPIN.AND.TCHK) THEN
        CALL HYPERFINE$SETPWRHO('SPIN',NGL,GVEC,RHOG(1,NSPIN))
      END IF
      CALL HYPERFINE$SETPWPOT(NGL,GVEC,VHARTREE)
!
!     ==================================================================
!     ==  GRAPHICS  PARAMETERS                                        ==
!     ==================================================================
      CALL GRAPHICS$SETPWPOT('ELECTROSTATIC',NGL,VHARTREE)
!
!     ==================================================================
!     ==  TRANSFORM  TOTAL CHARGE DENSITY AND POTENTIAL TO REAL SPACE ==
!     ==================================================================
                            CALL TIMING$CLOCKON('VOFRHO: FFT POTENTIAL')
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$SUPFFT('GTOR',NSPIN,NGL,RHOG,NRL,RHOE)
                          CALL TIMING$CLOCKOFF('VOFRHO: FFT POTENTIAL')
!
!     ==================================================================
!     ==  CALCULATE |NABLA|RHO IN G-SPACE AND TRANSFORM TO REAL SPACE ==
!     ==================================================================
      IF(TGRA) THEN
                        CALL TIMING$CLOCKON('VOFRHO: NABLA-RHO AND FFT')
        ALLOCATE(GRHO(NNR,3,NSPIN))
        CALL POTENTIAL_GRADIENT('RHO',NDIMD,NGL,GVEC,RHOG,NRL,GRHO) 
                       CALL TIMING$CLOCKOFF('VOFRHO: NABLA-RHO AND FFT')
      ELSE
        ALLOCATE(GRHO(1,1,1))
      END IF
!
!     ==================================================================
!     == CALCULATE XC-ENERGY AND POTENTIAL                            ==
!     ==================================================================
                           CALL TIMING$CLOCKON('VOFRHO: XC-POTENTIAL')
      EXC=0.D0
      CALL POTENTIAL_XC(TGRA,NSPIN,NRL,NRL,NR1GLOB*NR2*NR3,CELLVOL &
     &                 ,RHOE,GRHO,EXC,TSTRESS,STRESS1)
      STRESST(:,:)=STRESST(:,:)+STRESS1(:,:)
WRITE(*,FMT='("XC STRESS ",3F15.7)')STRESS1(1,:)
WRITE(*,FMT='("XC STRESS ",3F15.7)')STRESS1(2,:)
WRITE(*,FMT='("XC STRESS ",3F15.7)')STRESS1(3,:)
                           CALL TIMING$CLOCKOFF('VOFRHO: XC-POTENTIAL')
!
!     ==================================================================
!     == FOURIER TRANSFORM XC POTENTIALS                              ==
!     ================================================================== 
      IF(TGRA) THEN
                             CALL TIMING$CLOCKON('VOFRHO: GRAD-COR-POT')
        CALL POTENTIAL_GRADIENT('POT',NSPIN,NGL,GVEC,RHOG,NRL,GRHO)
                           CALL TIMING$CLOCKOFF('VOFRHO: GRAD-COR-POT')
      ELSE
        RHOG(:,:)=(0.D0,0.D0)
      END IF
      DEALLOCATE(GRHO)
      CALL PLANEWAVE$SELECT('DENSITY')
      ALLOCATE(CWORK(NGL,NSPIN))
      CALL PLANEWAVE$SUPFFT('RTOG',NSPIN,NGL,CWORK,NRL,RHOE)
      DO ISPIN=1,NSPIN
        DO IG=1,NGL
          RHOG(IG,ISPIN)=RHOG(IG,ISPIN)+CWORK(IG,ISPIN)
        ENDDO
      ENDDO        
      DEALLOCATE(CWORK)
!
!     ==================================================================
!     == ADD HARTREE POTENTIAL TO EXCHANGE POTENTIAL                  ==
!     ==================================================================
      DO IG=1,NGL
         RHOG(IG,1)=RHOG(IG,1)+VHARTREE(IG)
      ENDDO
      DEALLOCATE(VHARTREE)
!
!     ==================================================================
!     == CALCULATE FORCE ON THE PSEUDO CORE                           ==
!     ==================================================================
      CALL TIMING$CLOCKON('VOFRHO: FORCE ON PS-CORE')
      CALL POTENTIAL_FPSCORE(NSP,NAT,ISPECIES,RBAS,TAU0,FORCE1,STRESS1 &
     &              ,NGL,G2,GVEC,RHOG(1,1),PSCORG,DPSCORG)
      FIONT=FIONT+FORCE1
      STRESST=STRESST+STRESS1
WRITE(*,FMT='("CORE STRESS ",3F15.7)')STRESS1(1,:)
WRITE(*,FMT='("CORE STRESS ",3F15.7)')STRESS1(2,:)
WRITE(*,FMT='("CORE STRESS ",3F15.7)')STRESS1(3,:)
      CALL TIMING$CLOCKOFF('VOFRHO: FORCE ON PS-CORE')
!
!     ==================================================================
!     == SEND POTENTIAL TO OPTICS CODE                                ==
!     ==================================================================
!     CALL OPTIC3$VOFG(NR1,NR2,NR3,NGL,RHOE,NSPIN,INB1,INB2,INB3)
!
!     ==================================================================
!     ==  COMBINE ARRAYS FOR PARALELL PROCESSING                      ==
!     ==================================================================
!     ++ PARALLEL START ++++++++++++++++++++++++++++++++++++++++++++++++
      CALL MPE$COMBINE('+',EHARTREE)
      CALL MPE$COMBINE('+',FIONT)
      CALL MPE$COMBINE('+',VQLMT)
      CALL MPE$COMBINE('+',STRESST)
      CALL MPE$COMBINE('+',EXC)
!     ++ PARALLEL END ++++++++++++++++++++++++++++++++++++++++++++++++++
      STRESS(:,:)=STRESS(:,:)+STRESST(:,:)
      DO IAT=1,NAT
        DO I=1,3
          FION(I,IAT)=FION(I,IAT)+FIONT(I,IAT)
        ENDDO
        DO LM=1,LMRXX
          VQLM(LM,IAT)=VQLM(LM,IAT)+VQLMT(LM,IAT)
        ENDDO
      ENDDO
      RHOB=RHOBT
!
!     ==================================================================
!     ==  COMMUNICATE ENERGIES TO ENERGYLIST                          ==
!     ==================================================================
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EHARTREE)
      CALL ENERGYLIST$ADD('AE  ELECTROSTATIC',EHARTREE)
      CALL ENERGYLIST$ADD('PS  ELECTROSTATIC',EHARTREE)
!
      CALL ENERGYLIST$SET('PAIRPOTENTIAL',EPAIR)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EPAIR)
      CALL ENERGYLIST$ADD('AE  ELECTROSTATIC',EPAIR)
      CALL ENERGYLIST$ADD('PS  ELECTROSTATIC',EPAIR)
!   
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EISOLATE)
      CALL ENERGYLIST$ADD('AE  ELECTROSTATIC',EISOLATE)
      CALL ENERGYLIST$ADD('PS  ELECTROSTATIC',EISOLATE)
      CALL ENERGYLIST$SET('ISOLATE ENERGY',EISOLATE)
!
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EXC)
      CALL ENERGYLIST$ADD('AE  EXCHANGE-CORRELATION',EXC)
      CALL ENERGYLIST$ADD('PS  EXCHANGE-CORRELATION',EXC)
!
!     ==================================================================
!     ==  SELFTEST                                                    ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR=EHARTREE+EPAIR+EISOLATE+EXC
        IF(ITEST.EQ.1) THEN
          ALLOCATE(RHELP(2*NGL*NSPIN))
          RHELP=TRANSFER(RHOG,RHELP)*CELLVOL*2.D0
          CALL PLANEWAVE$GETI4('NGAMMA',NGAMMA)
          IF(NGAMMA.NE.0) THEN
            RHELP(NGAMMA:NGAMMA+1)=0.5D0*RHELP(NGAMMA:NGAMMA+1)
          END IF
          CALL SELFTEST$END('POTENTIAL-RHO',2*NGL*NSPIN,RHELP,SVAR,TBACK)
          DEALLOCATE(RHELP)
        ELSE IF(ITEST.EQ.2) THEN
          CALL SELFTEST$END('POTENTIAL-QLM',LMRXX*NAT,VQLM,SVAR,TBACK)
        ELSE IF(ITEST.EQ.3) THEN
          CALL SELFTEST$END('POTENTIAL-FORCE',3*NAT,FION,-SVAR,TBACK)
        END IF
        IF(TBACK) GOTO 1000
        CALL ERROR$MSG('NORMAL STOP AFTER SELFTEST')
        CALL ERROR$STOP('POTENTIAL_VOFRHO')
        DEALLOCATE(RHO_SELFTEST) 
      END IF
!
!     ==================================================================
!     ==  FOURIER TRANSFORM                                           ==
!     ==================================================================
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$SUPFFT('GTOR',NSPIN,NGL,RHOG,NRL,RHOE)
      DEALLOCATE(RHOG)
      RETURN
      END
!
!     ..................................................FPSCORE.........
      SUBROUTINE POTENTIAL_FPSCORE(NSP,NAT,ISPECIES,RBAS,TAU0,FORCE,STRESS &
     &                            ,NGL,G2,GVEC,VTEMP,PSCORG,DPSCORG)
!     ******************************************************************
!     **                                                              **
!     **  FORCE ON THE PSEUDO CORE                                    **
!     **                                                              **
!     **  F=F-V*NABLA*PSCORE                                          **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NSP
      INTEGER(4),INTENT(IN)   :: NAT
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT)
      REAL(8)   ,INTENT(IN)   :: RBAS(3,3)
      REAL(8)   ,INTENT(IN)   :: TAU0(3,NAT)
      REAL(8)   ,INTENT(OUT)  :: FORCE(3,NAT)
      REAL(8)   ,INTENT(OUT)  :: STRESS(3,3) 
      INTEGER(4),INTENT(IN)   :: NGL
      REAL(8)   ,INTENT(IN)   :: G2(NGL)
      REAL(8)   ,INTENT(IN)   :: GVEC(3,NGL)
      REAL(8)   ,INTENT(IN)   :: PSCORG(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: DPSCORG(NGL,NSP)
      COMPLEX(8),INTENT(IN)   :: VTEMP(NGL)
      REAL(8)                 :: GBAS(3,3)
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)    !(NGL)
      INTEGER(4),PARAMETER    :: NBLOCK=100
      REAL(8)   ,PARAMETER    :: TINY=1.D-300
      REAL(8)                 :: FAC,SVAR,SVAR1,SVAR2
      REAL(8)                 :: CELLVOL
      REAL(8)                 :: PI,Y0
      INTEGER(4)              :: IG,ISP,IAT,I,J
      REAL(8)                 :: FORCE1(3)
      REAL(8)                 :: STRESS1(3,3)
      COMPLEX(8)              :: CSVAR
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      FORCE(:,:)=0.D0
      STRESS(:,:)=0.D0
      FORCE1(:)=0.D0
      STRESS1(:,:)=0.D0      
      ALLOCATE(EIGR(NGL))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL PLANEWAVE$STRUCTUREFACTOR(TAU0(1,IAT),NGL,EIGR)
        DO IG=1,NGL
          CSVAR=CONJG(VTEMP(IG))*PSCORG(IG,ISP)*EIGR(IG)
!         == FORCES ====================================================
          SVAR=-AIMAG(CSVAR)
          DO I=1,3
            FORCE1(I)=FORCE1(I)+SVAR*GVEC(I,IG)
          ENDDO
!         == STRESSES =================================================
          SVAR1=-REAL(CONJG(VTEMP(IG))*DPSCORG(IG,ISP)*EIGR(IG))/(G2(IG)+TINY)
          DO I=1,3
            DO J=1,3
              STRESS1(I,J)=STRESS1(I,J)+GVEC(I,IG)*GVEC(J,IG)*SVAR1
            ENDDO
          ENDDO
!         ==  NOW TAKE CARE OF BLOCKING ================================ 
          IF(MOD(IG,NBLOCK).EQ.0) THEN
            FORCE(:,IAT)=FORCE(:,IAT)+FORCE1
            STRESS=STRESS+STRESS1
            FORCE1=0.D0
            STRESS1=0.D0
          END IF
        ENDDO
        FORCE(:,IAT)=FORCE(:,IAT)+FORCE1(:)
        STRESS=STRESS+STRESS1
        FORCE1=0.D0
        STRESS1=0.D0
      ENDDO
      DEALLOCATE(EIGR)
      FAC=Y0*2.D0*CELLVOL
      FORCE=FORCE*FAC
      STRESS=STRESS*FAC
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_GRADIENT(IDENT,NDIMD,NGL,GVEC,RHOOFG,NRL,GRHOOFR)
!     ******************************************************************
!     **  ID='RHO': COMPUTES THE GRADIENT OF A DENSITY PROVIDED       **
!     **            IN RECIPROCAL SPACE AND TRANSFORMS THE RESULT TO  **
!     **            REAL SPACE                                        **
!     **  ID='POT': TRANSFORMS A VECTOR ARRAY IN TO RECIPROCAL SPACE  **
!     **            AND COMPUTES THE DIVERGENCE                       **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      COMPLEX(8)  ,PARAMETER    :: CI=(0.D0,1.D0)
      CHARACTER(3),INTENT(IN)   :: IDENT  ! CAN BE 'RHO' OR 'POT'
      INTEGER(4)  ,INTENT(IN)   :: NDIMD
      INTEGER(4)  ,INTENT(IN)   :: NRL
      INTEGER(4)  ,INTENT(IN)   :: NGL
      REAL(8)     ,INTENT(IN)   :: GVEC(3,NGL)
      COMPLEX(8)  ,INTENT(INOUT):: RHOOFG(NGL,NDIMD)
      REAL(8)     ,INTENT(INOUT):: GRHOOFR(NRL,3,NDIMD)
      COMPLEX(8)  ,ALLOCATABLE  :: CWORK(:,:,:)    !(NGL,3,NGL)
      COMPLEX(8)                :: CSVAR
      INTEGER(4)                :: I,IG,IDIM
!     ******************************************************************
      CALL PLANEWAVE$SELECT('DENSITY')
!
!     ==================================================================
!     ==  CALCULATE GRADIENT OF RHO,GRHO AND OBTAIN RESULT IN R-SPACE ==
!     ==================================================================
      ALLOCATE(CWORK(NGL,3,NDIMD))
      IF(IDENT.EQ.'RHO') THEN
        DO IDIM=1,NDIMD
          DO I=1,3
            DO IG=1,NGL
              CWORK(IG,I,IDIM)=RHOOFG(IG,IDIM)*CI*GVEC(I,IG)
            ENDDO
          ENDDO
        ENDDO
        CALL PLANEWAVE$SUPFFT('GTOR',3*NDIMD,NGL,CWORK,NRL,GRHOOFR)
      ELSE IF(IDENT.EQ.'POT') THEN
!       ================================================================
!       == CALCULATE DIVERGENCE OF THE GRADIENTPOTENTIAL IN G-SPACE   ==
!       ================================================================
        CALL PLANEWAVE$SUPFFT('RTOG',3*NDIMD,NGL,CWORK,NRL,GRHOOFR)
        RHOOFG(:,:)=(0.D0,0.D0)
        DO IDIM=1,NDIMD
          DO I=1,3
            DO IG=1,NGL
              CSVAR=-CI*GVEC(I,IG)
              RHOOFG(IG,IDIM)=RHOOFG(IG,IDIM)+CWORK(IG,I,IDIM)*CSVAR
            ENDDO
          ENDDO
        ENDDO      
      ELSE
        CALL ERROR$MSG('INCORRECT IDENTIFIER')
        CALL ERROR$CHVAL('IDENT',IDENT)
        CALL ERROR$STOP ('POTENTIAL_GRADIENT')
      END IF
      DEALLOCATE(CWORK)         
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_ADDCORE(NGL,NSP,NAT,ISPECIES,RAT,PSCORG,RHO)
!     ******************************************************************
!     ** ADD THE PSEUDO CORE DENSITY TO THE PSEUDO DENSITY RHO        **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NGL            ! #(G-VECTORS)
      INTEGER(4),INTENT(IN)   :: NAT            ! #(ATOMS)
      INTEGER(4),INTENT(IN)   :: NSP            ! #(ELEMENTS)
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT)  !
      REAL(8)   ,INTENT(IN)   :: RAT(3,NAT)     ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(IN)   :: PSCORG(NGL,NSP)! PS-CORE
      COMPLEX(8),INTENT(INOUT):: RHO(NGL)       ! IN PS-DENSITY; OUT: +PS-CORE
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)        !(NGL) STRUCTURE FACTOR
      REAL(8)                 :: PI,Y0
      INTEGER(4)              :: IAT,ISP,IG
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      ALLOCATE(EIGR(NGL))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL PLANEWAVE$STRUCTUREFACTOR(RAT(1,IAT),NGL,EIGR)

        DO IG=1,NGL
          RHO(IG)=RHO(IG)+PSCORG(IG,ISP)*Y0*EIGR(IG)
        ENDDO
      ENDDO
      DEALLOCATE(EIGR)
      RETURN 
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_HARTREE(LMRXX,LRXX &
     &            ,NGL,NSP,NAT,ISPECIES,LMRX &
     &            ,RHO,POT,QLM,VQLM,TAU0,FORCE,EHARTREE &
     &            ,GWEIGHT,G2,GVEC,YLMOFG,VBARG,DVBARG,G0,V0,RHOB &
     &            ,TSTRESS,DG0,DV0,STRESS)
!     ******************************************************************
!     **                                                              **
!     ** CALCULATES THE HARTREE ENERGY IN RECIPROCAL SPACE            **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    RHO       DENSITY OF THE PS WAVE FUNCTIONS (G-SPACE)     **
!     **                                                              **
!     **  OUTPUT:                                                     **
!     **    POT      HARTREE POTENTIAL (G-SPACE)                    **
!     **    RHO       DENSITY OF PS WVAE FUNCTIONS AND PS CORE       **
!     **    FORCE       HARTREE FORCE ON THE ATOMS (ADDED TO INPUT)    **
!     **    EHARTREE   PS HARTREE ENERGY (WITHOUT PAIRPOTENTIAL)      **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
      INTEGER(4),INTENT(IN)   :: LMRXX
      INTEGER(4),INTENT(IN)   :: LRXX 
      INTEGER(4),INTENT(IN)   :: NGL
      INTEGER(4),INTENT(IN)   :: NAT             ! #(ATOMS)
      INTEGER(4),INTENT(IN)   :: NSP             ! #(ELEMENTS)
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT)
      INTEGER(4),INTENT(IN)   :: LMRX(NSP)
      COMPLEX(8),INTENT(IN)   :: RHO(NGL)        ! IN PS-DENSITY; OUT: +PS-CORE
      COMPLEX(8),INTENT(OUT)  :: POT(NGL)        ! HARTREE-POTENTIAL
      REAL(8)   ,INTENT(IN)   :: QLM(LMRXX,NAT)  ! MOMENTS OF 1C-DENSITY
      REAL(8)   ,INTENT(OUT)  :: VQLM(LMRXX,NAT) ! DE/DQLM
      REAL(8)   ,INTENT(IN)   :: TAU0(3,NAT)     ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(OUT)  :: FORCE(3,NAT)    ! FORCE ON ATOMS (ADDED)
      REAL(8)   ,INTENT(OUT)  :: EHARTREE        ! PS-HARTREE ENERGY (WITHOUT PAIRPOT)
      REAL(8)   ,INTENT(IN)   :: GWEIGHT         ! UNIT CELL VOLUME
      REAL(8)   ,INTENT(IN)   :: G2(NGL)         ! G**2
      REAL(8)   ,INTENT(IN)   :: GVEC(3,NGL)     ! G-VECTORS
      REAL(8)   ,INTENT(IN)   :: YLMOFG(LMRXX,NGL)! SPHERICAL HARMONICS
      REAL(8)   ,INTENT(IN)   :: G0(NGL,NSP)     ! WIDE GAUSSIANS OF COMPENSATION DENSITY
      REAL(8)   ,INTENT(IN)   :: V0(NGL,NSP)     ! POTENTIAL OF KOMPENSATION DENSITIES
      REAL(8)   ,INTENT(IN)   :: VBARG(NGL,NSP)  ! V-HAT
      REAL(8)   ,INTENT(IN)   :: DVBARG(NGL,NSP) ! |G|D(V-HAT)/DG
      REAL(8)   ,INTENT(OUT)  :: RHOB            ! DENSITY OF THE COMPENSATING BACKGROUND
      LOGICAL(4),INTENT(IN)   :: TSTRESS
      REAL(8)   ,INTENT(IN)   :: DG0(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: DV0(NGL,NSP)
      REAL(8)   ,INTENT(OUT)  :: STRESS(3,3)
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)         !(NGL)
      COMPLEX(8),ALLOCATABLE  :: RHO2(:)         !(NGL) psrho+compensation
      COMPLEX(8),ALLOCATABLE  :: VG(:)           !(NGL)
      COMPLEX(8)              :: CFAC(LMRXX)
      INTEGER(4)              :: NGAMMA  
      COMPLEX(8)              :: CSVAR,CSVAR0,AVAL,BVAL,CVAL,DVAL,CSUM
      COMPLEX(8)              :: EIGR1
      REAL(8)                 :: VB               !VHAT(GAMMA)
      INTEGER(4)              :: M,L,LM,ISP,IAT,IG,I,J,LM1,LM2
      REAL(8)                 :: SVAR,SVAR1,SVAR2,FPIBG2
      REAL(8)                 :: PI,Y0,FPI
      REAL(8)                 :: CC(3,3)
      COMPLEX(8)              :: P0(3,3,LMRXX),PM(3,3,LMRXX),PT(3,3)
      REAL(8)   ,PARAMETER    :: TINY=1.D-300
      REAL(8)                 :: STRESSA(3,3),STRESSB(3,3),STRESSC(3,3)
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FPI=4.D0*PI
      Y0=1.D0/DSQRT(FPI)
      CALL PLANEWAVE$GETI4('NGAMMA',NGAMMA)
!
      ALLOCATE(EIGR(NGL))
      ALLOCATE(RHO2(NGL))
      ALLOCATE(VG(NGL))
!
!     ==================================================================
!     ==  ADD COMPENSATION CHARGE DENSITY                             ==
!     ==  RHO = PS-RHO + PS-CORE                                      ==
!     ==  RHO2 = PS-RHO + RHO-HAT + PS-CORE + RHO-BACKGROUND          ==
!     ==  POT = VBAR + VHAT                                           ==
!     ==  VB   = VHAT(GAMMA)                                          ==
!     ==================================================================
      VB=0.D0
      POT(:)=(0.D0,0.D0)
      RHO2(:)=RHO(:)
!
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LM=0
        CSVAR=CI       
        DO L=0,LRXX
          CSVAR=-CSVAR*CI/DBLE(2*L+1)   
          DO M=-L,L
            LM=LM+1
            CFAC(LM)=CSVAR*QLM(LM,IAT)
          ENDDO
        ENDDO
        CALL PLANEWAVE$STRUCTUREFACTOR(TAU0(1,IAT),NGL,EIGR)
        IF(NGAMMA.NE.0) THEN
          VB=VB+QLM(1,IAT)*Y0*V0(NGAMMA,ISP)
        END IF
        DO IG=1,NGL
          CSVAR=(0.D0,0.D0)
          DO LM=1,LMRX(ISP)
            CSVAR=CSVAR+CFAC(LM)*YLMOFG(LM,IG)
          ENDDO
          CSVAR=CSVAR*EIGR(IG)
          RHO2(IG)    =RHO2(IG)    +CSVAR*G0(IG,ISP)
          POT(IG)=POT(IG)+CSVAR*V0(IG,ISP)+VBARG(IG,ISP)*Y0*EIGR(IG)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ADD BACKGROUND                                              ==
!     ==================================================================
      IF(NGAMMA.NE.0) THEN
        RHOB=-REAL(RHO2(NGAMMA))
PRINT*,'RHOB  FROM VOFRHO_HARTREE',RHOB*GWEIGHT
PRINT*,'Q(RHO) FROM VOFRHO_HARTREE',RHO(NGAMMA)*GWEIGHT
PRINT*,'Q(RHO2) FROM VOFRHO_HARTREE',RHO2(NGAMMA)*GWEIGHT
        RHO2(NGAMMA)=(0.D0,0.D0)   !
      ELSE 
        RHOB=0.D0
      END IF
!
!     ==================================================================
!     ==   CALCULATE ENERGY                                           ==
!     ==================================================================
      CSUM=(0.D0,0.D0)
      CSVAR=0.D0
      DO IG=1,NGL
        IF(IG.NE.NGAMMA) THEN
          VG(IG)=RHO2(IG)*FPI/G2(IG)
          CSUM=CSUM+0.5D0*CONJG(RHO2(IG))*VG(IG)+CONJG(RHO(IG))*POT(IG)
CSVAR=CSVAR+0.5D0*CONJG(RHO2(IG))*VG(IG)
          POT(IG)=POT(IG)+VG(IG)
        ELSE 
          VG(IG)=-VB
          CSUM=CSUM+0.5D0*( CONJG(RHO(IG))*POT(IG)+RHOB*VB )
          POT(IG)=POT(IG)+VG(IG)
        END IF
      ENDDO
      EHARTREE=REAL(CSUM)
PRINT*,'PURE HARTREE ENERGY ',REAL(CSVAR)*2.D0*GWEIGHT
PRINT*,'RHO*POT ',REAL(CSUM-CSVAR)*2.D0*GWEIGHT
PRINT*,'STRESS PREDICTED ',-REAL(CSVAR)*2.D0*GWEIGHT/3.D0-REAL(CSUM-CSVAR)*2.D0*GWEIGHT
!
!     == NOW TAKE CARE OF STRESSES =====================================
!     == INTEGRATION WEIGHT IS DONE LATER ==============================
      STRESS(:,:)=0.D0
      IF(TSTRESS) THEN
        DO IG=1,NGL
          IF(IG.NE.NGAMMA) THEN         
            SVAR=VG(IG)*CONJG(RHO2(IG))
            SVAR1=SVAR/G2(IG)
            DO I=1,3
              DO J=1,3
                STRESS(I,J)=STRESS(I,J)+GVEC(I,IG)*GVEC(J,IG)*SVAR1
              ENDDO
            ENDDO
          END IF
        ENDDO
        DO I=1,3
          STRESS(I,I)=STRESS(I,I)-EHARTREE
        ENDDO
      ENDIF

WRITE(*,FMT='("MAXWELL STRESS ",3F15.7)')STRESS(1,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL STRESS ",3F15.7)')STRESS(2,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL STRESS ",3F15.7)')STRESS(3,:)*2.D0*GWEIGHT
STRESSA=0.D0
STRESSB=0.D0
STRESSC=0.D0
!
!     ==================================================================
!     ==   CALCULATE FORCES AND TOTAL ENERGY                          ==
!     ==      VBARG  = V_ADD                                          ==
!     ==      V0 = POTENTIAL FROM THE DIFFERENCE                      ==
!     ==           BETWEEN THE COMPENSATION CHARGES                   ==
!     ==================================================================
!     == N_A=RHOTILDE+N_B ==============================================
      DO IG=1,NGL
        RHO2(IG)=RHO(IG)
      ENDDO
      IF(NGAMMA.NE.0) RHO2(NGAMMA)=RHO2(NGAMMA)+RHOB
!     == (-I)**L/(2L+1)!! ==============================================
      LM=0
      CSVAR=CI
      DO L=0,LRXX
        CSVAR=-CSVAR*CI/DBLE(2*L+1)
        DO M=-L,L
          LM=LM+1
          CFAC(LM)=CSVAR
        ENDDO
      ENDDO
      VQLM(:,:)=0.D0
      FORCE(:,:)=0.D0
      PT(:,:)=0.D0
      DO IAT=1,NAT
        IF(TSTRESS) THEN
          P0(:,:,:)=0.D0
          PM(:,:,:)=0.D0
          DO LM1=1,LMRX(ISP)
            DO I=1,10
              CALL SPHERICAL$CCMAT0(LM1,I,LM2,CC)
              IF(LM2.EQ.0) EXIT
              P0(:,:,LM2)=P0(:,:,LM2)+CFAC(LM1)*QLM(LM1,IAT)*CC(:,:)
            ENDDO
            DO I=1,10
              CALL SPHERICAL$CCMATM(LM1,I,LM2,CC)
              IF(LM2.EQ.0) EXIT
              PM(:,:,LM2)=PM(:,:,LM2)+CFAC(LM1)*QLM(LM1,IAT)*CC(:,:)
            ENDDO
          ENDDO
        END IF

        ISP=ISPECIES(IAT)
        CALL PLANEWAVE$STRUCTUREFACTOR(TAU0(1,IAT),NGL,EIGR)
        DO IG=1,NGL
          IF(IG.EQ.NGAMMA) THEN
            CSVAR0=(G0(IG,ISP)*CONJG(VG(IG))+V0(IG,ISP)*CONJG(RHO2(IG)))
            CSVAR=CSVAR0*CFAC(1)*YLMOFG(1,IG)*EIGR(IG)
            VQLM(1,IAT)=VQLM(1,IAT)-0.5D0*REAL(CSVAR)
          END IF
!         == MULTIPOLE POTENTIALS ======================================
          AVAL=( G0(IG,ISP)*CONJG(VG(IG)) + V0(IG,ISP)*CONJG(RHO2(IG)) )
          AVAL=AVAL*EIGR(IG)
          CSUM=0.D0
          DO LM=1,LMRX(ISP)
            CSVAR=CFAC(LM)*YLMOFG(LM,IG)
            VQLM(LM,IAT)=VQLM(LM,IAT)+REAL(CSVAR*AVAL)
            CSUM=CSUM   + QLM(LM,IAT)*CSVAR
          ENDDO
!         == FORCES ====================================================
          CVAL=VBARG(IG,ISP)*Y0*EIGR(IG)*CONJG(RHO(IG))
          SVAR=AIMAG(CVAL+CSUM*AVAL)
          FORCE(1,IAT)=FORCE(1,IAT)-GVEC(1,IG)*SVAR
          FORCE(2,IAT)=FORCE(2,IAT)-GVEC(2,IG)*SVAR
          FORCE(3,IAT)=FORCE(3,IAT)-GVEC(3,IG)*SVAR
!         == NOW STRESSES ==============================================
          IF(TSTRESS) THEN
            DO LM=1,LMRX(ISP)
              PT(:,:)=PT(:,:)+(P0(:,:,LM)+PM(:,:,LM)*G2(IG))*YLMOFG(LM,IG) 
            ENDDO
            BVAL=DG0(IG,ISP)*CONJG(VG(IG))+DV0(IG,ISP)*CONJG(RHO2(IG)) 
            BVAL=BVAL*EIGR(IG)
            DVAL=DVBARG(IG,ISP)*Y0*EIGR(IG)*CONJG(RHO2(IG))
            SVAR=REAL(DVAL+BVAL*CSUM)/(G2(IG)+TINY)
            DO I=1,3
              DO J=1,3
                STRESS(I,J)=STRESS(I,J)- REAL(PT(I,J)*AVAL) &
     &                                 - GVEC(I,IG)*GVEC(J,IG)*SVAR
!
                STRESSA(I,J)=STRESSA(I,J)- REAL(PT(I,J)*AVAL) 
                STRESSB(I,J)=STRESSB(I,J)- GVEC(I,IG)*GVEC(J,IG)*REAL(DVAL)/(G2(IG)+TINY)
                STRESSC(I,J)=STRESSC(I,J)- GVEC(I,IG)*GVEC(J,IG)*REAL(BVAL*CSUM)/(G2(IG)+TINY)
              ENDDO
            ENDDO
            PT(:,:)=0.D0  ! RESET FOR NEXT G-VECTOR
          END IF
        ENDDO
      ENDDO
WRITE(*,FMT='("MAXWELL++  STRESS ",3F15.7)')STRESSA(1,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL++  STRESS ",3F15.7)')STRESSA(2,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL++  STRESS ",3F15.7)')STRESSA(3,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL+++ STRESS ",3F15.7)')STRESSB(1,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL+++ STRESS ",3F15.7)')STRESSB(2,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL+++ STRESS ",3F15.7)')STRESSB(3,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL+++ STRESS ",3F15.7)')STRESSC(1,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL+++ STRESS ",3F15.7)')STRESSC(2,:)*2.D0*GWEIGHT
WRITE(*,FMT='("MAXWELL+++ STRESS ",3F15.7)')STRESSC(3,:)*2.D0*GWEIGHT
!
!     ==================================================================
!     ==  MULTIPLY WITH INTEGRATION WEIGHT                            ==
!     ==  FACTOR TWO SINCE ONLY HALF OF THE G-VECTORS ARE TREATED     ==
!     ==================================================================
      EHARTREE=EHARTREE*2.D0*GWEIGHT
      VQLM=VQLM*2.D0*GWEIGHT
      FORCE=FORCE*2.D0*GWEIGHT
      STRESS=STRESS*2.D0*GWEIGHT

      DEALLOCATE(VG)
      DEALLOCATE(EIGR)
      DEALLOCATE(RHO2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_XC &
     &     (TGRA,NSPIN,NNRX,NNR,NNRSCAL,CELLVOL,RHO,GRHO,EXC,TSTRESS,STRESS)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES EXCHANGE AND CORRELATION ENERGY OF THE           **
!     **  PSEUDO DENSITY                                              **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **  1)RHO AND GRHO ARE OVERWRITTEN WITH THE RESPECTIVE          **
!     **    POTENTIALS                                                **
!     **  2) GRHO MUST NOT BE USED IF TGRA=.FALSE.                    **
!     **  3) THE SUM OVER R-POINTS IS BLOCKED FOR BETTER ACCURACY     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)   :: TGRA         ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4),INTENT(IN)   :: NSPIN        ! NUMBER OF SPIN DIRECTIONS (1 OR 2)
      INTEGER(4),INTENT(IN)   :: NNRX         ! DIMENSION FOR RHO AND GRHO
      INTEGER(4),INTENT(IN)   :: NNR          ! LOCAL #(R-POINTS) (THIS PROCESSOR)
      INTEGER(4),INTENT(IN)   :: NNRSCAL      ! GLOBAL #(R-POINTS) (ALL PROCESSORS)
      REAL(8)   ,INTENT(IN)   :: CELLVOL      ! VOLUME OF THE UNIT CELL
      REAL(8)   ,INTENT(OUT)  :: EXC          ! EXCHANGE AND CORRELATION ENERGY
      REAL(8)   ,INTENT(INOUT):: RHO(NNRX,NSPIN)
      REAL(8)   ,INTENT(INOUT):: GRHO(NNR,3,NSPIN) ! GRHO MUST NOT BE USED IF TGRA=.TRUE.
      LOGICAL(4),INTENT(IN)   :: TSTRESS
      REAL(8)   ,INTENT(OUT)  :: STRESS(3,3)
      INTEGER(4),PARAMETER    :: NBLOCK=100    ! BLOCKS THE SUMMATION OVER GRID POINTS
      LOGICAL(4)              :: TSPIN
      REAL(8)                 :: SXC,SSXC
      INTEGER(4)              :: IR,I,J
      REAL(8)                 :: EXC1
      REAL(8)                 :: RHOT                   ! TOTAL DENSITY (RHOUP+RHODOWN)
      REAL(8)                 :: RHOS                   ! SPIN DENSITY  (RHOUP-RHODOWN)
      REAL(8)                 :: GRHOTX,GRHOTY,GRHOTZ   ! GRAD(RHOT)
      REAL(8)                 :: GRHOSX,GRHOSY,GRHOSZ   ! GRAD(RHOS)
      REAL(8)                 :: GRHOT2                 ! |GRAD(RHOT)|**2
      REAL(8)                 :: GRHOS2                 ! |GRAD(RHOS)|**2
      REAL(8)                 :: GRHOST                 ! GRAD(RHOS)*GRAD(RHOT)
      REAL(8)                 :: VT                     ! DE/DRHOT
      REAL(8)                 :: VS                     ! DE/DRHOS
      REAL(8)                 :: GVT2                   ! DE/D(GRHOT2)
      REAL(8)                 :: GVS2                   ! DE/D(GRHOS2)
      REAL(8)                 :: GVST                   ! DE/D(GRHOST)
      REAL(8)                 :: SVAR,SVAR1,SVAR2
      REAL(8)                 :: STRESS1(3,3),DIAGSTRESS,DIAGSTRESS1
!     ******************************************************************
      TSPIN=(NSPIN.EQ.2) 
!
      SXC=0.D0
      RHOS=0.D0
      GRHOT2=0.D0
      GRHOS2=0.D0
      GRHOST=0.D0
!
!     ==================================================================
!     ==  BLOCKED LOOP OVER R-POINTS                                  ==
!     ==================================================================
      SSXC=0.D0    ! FINAL SUM OF EXCHANGE ENERGIES
      SXC=0.D0     ! HOLDS INTERMEDIATE SUM
      STRESS=0.D0
      STRESS1=0.D0
      DIAGSTRESS=0.D0
      DIAGSTRESS1=0.D0
      DO IR=1,NNR
!
!       ================================================================
!       == EVALUATE DENSITIES AND GRADIENTS                           ==
!       ================================================================
        RHOT=RHO(IR,1)
        IF(TSPIN) THEN
          RHOS=RHO(IR,NSPIN)
        ELSE 
          RHOS=0.D0
        END IF
        IF(TGRA) THEN
          GRHOTX=GRHO(IR,1,1)
          GRHOTY=GRHO(IR,2,1)
          GRHOTZ=GRHO(IR,3,1)
          GRHOT2=GRHOTX**2+GRHOTY**2+GRHOTZ**2 ! |GRAD(RHOT)|**2
          IF(TSPIN) THEN
            GRHOSX=GRHO(IR,1,NSPIN)
            GRHOSY=GRHO(IR,2,NSPIN)
            GRHOSZ=GRHO(IR,3,NSPIN)
            GRHOST=GRHOSX*GRHOTX+GRHOSY*GRHOTY+GRHOSZ*GRHOTZ ! GRAD(RHOS)*GRAD(RHOT)
            GRHOS2=GRHOSX**2+GRHOSY**2+GRHOSZ**2             ! |GRAD(RHOS)|**2
          ELSE
            GRHOS2=0.D0
            GRHOST=0.D0
          END IF
        ELSE
          GRHOT2=0.D0
          GRHOS2=0.D0
          GRHOST=0.D0
        END IF
!
!       ================================================================
!       == EVALUATE DENSITY FUNCTIONAL                                ==
!       ================================================================
        CALL DFT(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST,EXC1,VT,VS,GVT2,GVS2,GVST)
!
!       ================================================================
!       == ADD TO TOTAL ENERGY AND POTENTIAL                          ==
!       ================================================================
        SXC=SXC+EXC1
        RHO(IR,1)=VT
        IF(TSPIN) THEN
          RHO(IR,NSPIN)=VS
        END IF
        IF(TGRA) THEN
          SVAR=2.D0*GVT2
          GRHO(IR,1,1)=GRHOTX*SVAR       
          GRHO(IR,2,1)=GRHOTY*SVAR       
          GRHO(IR,3,1)=GRHOTZ*SVAR       
          IF(TSPIN) THEN
            SVAR=2.D0*GVS2
            GRHO(IR,1,NSPIN)=GRHOSX*SVAR+GRHOTX*GVST
            GRHO(IR,2,NSPIN)=GRHOSY*SVAR+GRHOTY*GVST
            GRHO(IR,3,NSPIN)=GRHOSZ*SVAR+GRHOTZ*GVST
            GRHO(IR,1,1)=GRHO(IR,1,1)+GRHOSX*GVST
            GRHO(IR,2,1)=GRHO(IR,2,1)+GRHOSY*GVST
            GRHO(IR,3,1)=GRHO(IR,3,1)+GRHOSZ*GVST
          END IF
        END IF
!
!       ================================================================
!       == NOW ADD TO STRESSES                                        ==
!       ================================================================
        IF(TSTRESS) THEN
          DIAGSTRESS1=DIAGSTRESS1+EXC1-RHOT*VT
          IF(TSPIN) THEN
            DIAGSTRESS1=DIAGSTRESS1-RHOS*VS
          END IF
          IF(TGRA) THEN
            SVAR=2.D0*GVT2
            DIAGSTRESS1=DIAGSTRESS1-GRHOT2*SVAR
            STRESS1(1,1)=STRESS1(1,1)-GRHOTX*GRHOTX*SVAR
            STRESS1(1,2)=STRESS1(1,2)-GRHOTX*GRHOTY*SVAR
            STRESS1(1,3)=STRESS1(1,3)-GRHOTX*GRHOTZ*SVAR
            STRESS1(2,2)=STRESS1(2,2)-GRHOTY*GRHOTY*SVAR
            STRESS1(2,3)=STRESS1(2,3)-GRHOTY*GRHOTZ*SVAR
            STRESS1(3,3)=STRESS1(3,3)-GRHOTZ*GRHOTZ*SVAR
            IF(TSPIN) THEN
              SVAR1=2.D0*GVS2
              SVAR2=GVST
              DIAGSTRESS=DIAGSTRESS-GRHOS2*SVAR1-2.D0*GRHOST*SVAR2
              STRESS1(1,1)=STRESS1(1,1)-GRHOSX*GRHOSX*SVAR1 &
     &                            -2.D0*GRHOSX*GRHOTX*SVAR2
              STRESS1(2,2)=STRESS1(2,2)-GRHOSY*GRHOSY*SVAR1 &
     &                            -2.D0*GRHOSY*GRHOTY*SVAR2
              STRESS1(3,3)=STRESS1(3,3)-GRHOSZ*GRHOSZ*SVAR1 &
     &                            -2.D0*GRHOSZ*GRHOTZ*SVAR2
              STRESS1(1,2)=STRESS1(1,2)-GRHOSX*GRHOSY*SVAR1 &
     &                 -(GRHOTX*GRHOSY+GRHOSX*GRHOTY)*SVAR2
              STRESS1(1,3)=STRESS1(1,3)-GRHOSX*GRHOSZ*SVAR1 &
     &                 -(GRHOTX*GRHOSZ+GRHOSX*GRHOTZ)*SVAR2
              STRESS1(2,3)=STRESS1(2,3)-GRHOSY*GRHOSZ*SVAR1 &
     &                 -(GRHOTY*GRHOSZ+GRHOSY*GRHOTZ)*SVAR2
            END IF
          END IF
        END IF
!
!       ================================================================
!       == TAKE CARE OF BLOCKING                                      ==
!       ================================================================
        IF(MOD(IR,NBLOCK).EQ.0) THEN
          SSXC=SSXC+SXC               ! SWAP INTERMEDIATE SUM 
          DIAGSTRESS=DIAGSTRESS+DIAGSTRESS1
          STRESS=STRESS+STRESS1
          SXC=0.D0                    ! INTO FINAL SUM
          STRESS1=0.D0
          DIAGSTRESS1=0.D0
        END IF
      ENDDO
!     == CLEAN OUT INTERMEDIATE SUM ====================================
      SSXC=SSXC+SXC  ! 
      DIAGSTRESS=DIAGSTRESS+DIAGSTRESS1
      STRESS=STRESS+STRESS1
!     == SET TO ZERO FOR AESTETICAL REASONS ============================
      SXC=0.D0                                      
      STRESS1=0.D0
      DIAGSTRESS1=0.D0
!
!     ==================================================================
!     == SCALE EXCHANGE ENERGY WITH THE INTEGRATION WEIGHT            ==
!     ==================================================================
      EXC=SSXC*CELLVOL/DBLE(NNRSCAL)
      IF(TSTRESS) THEN
        SVAR=CELLVOL/DBLE(NNRSCAL)
        DO I=1,3
          STRESS(I,I)=(STRESS(I,I)+DIAGSTRESS)*SVAR
          DO J=I+1,3
            STRESS(I,J)=STRESS(I,J)*SVAR
            STRESS(J,I)=STRESS(I,J)
          ENDDO
        ENDDO
      END IF
!
      IF(EXC.NE.EXC) THEN
        CALL ERROR$MSG('POTENTIAL_XC')
        CALL ERROR$R8VAL('EXC',EXC)
        CALL ERROR$STOP('POTENTIAL_XC')
      END IF
      RETURN 
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_YLMOFG(LMXX,LMX,NG,GVEC,YLMOFG)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES Y_LM(G)*G**L                                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: LMXXX=36
      INTEGER(4),INTENT(IN) :: LMXX
      INTEGER(4),INTENT(IN) :: LMX
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(IN) :: GVEC(3,NG)
      REAL(8)   ,INTENT(OUT):: YLMOFG(LMXX,NG)
      REAL(8)               :: YLM(LMXXX)
      INTEGER(4)            :: LX
      INTEGER(4)            :: LM,IG,L,M
      INTEGER(4)            :: ISTART
      REAL(8)               :: GTOL,GLENG
      REAL(8)               :: Y0
!     ******************************************************************
      Y0=1.D0/SQRT(16.D0*DATAN(1.D0))  ! =1/SQRT(4*PI)
      LX=INT(SQRT(REAL(LMX-1))+1.E-8)
      IF(LMX.NE.(LX+1)**2) THEN
        CALL ERROR$MSG('LMX MUST CONTAIN ALL MAGNETIC QUANTUM NUMBER')
        CALL ERROR$MSG('FOR EACH MAIN ANGULAR MOMENTUM')
        CALL ERROR$I4VAL('LMX ',LMX)
        CALL ERROR$STOP('POTENTIAL_YLMOFG')
      END IF
      IF(LMXX.GT.LMXXX) THEN
        CALL ERROR$STOP('POTENTIAL_YLMOFG')
      END IF
!
!     ==================================================================
!     ==  CALCULATE REAL SPHERICAL HARMONICS TIME G**L                ==
!     ==================================================================
      DO LM=1,LMXX
        DO IG=1,NG
          YLMOFG(LM,IG)=0.D0
        ENDDO
      ENDDO
      DO IG=1,NG
        GLENG=DSQRT(GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2)
        IF(GLENG.LT.1.D-7) THEN
          YLMOFG(:,IG)=0.D0
          YLMOFG(1,IG)=Y0
          CYCLE
        END IF
        CALL GETYLM(LMX,GVEC(1,IG),YLM)
        GTOL=1.D0
        LM=0
        DO L=0,LX
          DO M=1,2*L+1
            LM=LM+1
            YLMOFG(LM,IG)=YLM(LM)*GTOL
          ENDDO
          GTOL=GTOL*GLENG
        ENDDO
      ENDDO
      RETURN
      END
