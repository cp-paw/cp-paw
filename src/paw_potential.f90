!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE POTENTIAL_MODULE
      LOGICAL(4) :: TINI=.FALSE.  ! PLANE WAVE GRID DEFINED OR NOT
      LOGICAL(4) :: TSET=.FALSE.  ! MODULE DATA DEFINED OR NOT
      LOGICAL(4) :: TCONFINE=.FALSE. ! CONFINING POTENTIAL OR NOT
      REAL(8)    :: VCONFINE=0.D0
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$SETR8(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN):: VAL
!     **************************************************************************
      IF(ID.EQ.'VCONFINE') THEN
        VCONFINE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('POTENTIAL$SETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$GETR8(ID,VALUE)
!     **************************************************************************
!     **************************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VALUE
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$SETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'STRESS') THEN
        TSTRESS=VAL
      ELSE IF(ID.EQ.'FORCE') THEN
        TFORCE=VAL
      ELSE IF(ID.EQ.'CONFINE') THEN
!       == APPLIES A REPULSIVE POTENTIAL AROUND THE MOLECULE DESCRIBING THE ====
!       == PAULI REPULSION OF A SOLVENT ENVIRONMENT ============================
        TCONFINE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('POTENTIAL$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$SETR8A(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$REPORT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **  REPORT ON THE POTENTIAL OBJECT                                      **
!     **                                                                      **
!     **************************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      INTEGER(4)    :: NFIL
      REAL(8)       :: MBYTE
      REAL(8)       :: RY
      REAL(8)       :: MEMORY
      INTEGER(4)    :: NGG
!     **************************************************************************
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
      IF(TCONFINE) THEN
        CALL REPORT$STRING(NFIL,'SOLVENT PAULI-REPULSION INCLUDED')
        CALL REPORT$R8VAL(NFIL,'HEIGHT OF CONFINING POTENTIAL',VCONFINE,'H')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$CLEAR
!     **************************************************************************
!     **  CLEAR PERMANENT MEMORY OF THE POTENTIAL OBJECT                      **
!     **  (WILL AUTOMATICALLY BE RECREATED WHEN REQUIRED)                     **
!     **************************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$INITIALIZE(EPWRHO_,NR1START_,NR1L_,NR2_,NR3_)
!     **************************************************************************
!     **  DEFINE GRID FOR FOURIER TRANSFORM                                   **
!     **************************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: EPWRHO_
      INTEGER(4),INTENT(IN) :: NR1START_
      INTEGER(4),INTENT(IN) :: NR1L_
      INTEGER(4),INTENT(IN) :: NR2_
      INTEGER(4),INTENT(IN) :: NR3_
      REAL(8)               :: RBAS(3,3)     ! LATTICE VECTORS
      REAL(8)   ,PARAMETER  :: K(3)=(/0.D0,0.D0,0.D0/)
      INTEGER(4)            :: ISP
!     **************************************************************************
                                CALL TRACE$PUSH('POTENTIAL$INITIALIZE')
      IF(TINI) THEN
        CALL ERROR$MSG('POTENTIAL OBJECT IS ALREADY INITIALIZED')
        CALL ERROR$STOP('POTENTIAL$INITIALIZE')
      END IF
      TINI=.TRUE.
      TSET=.FALSE.
!     
!     ==========================================================================
!     ==  MAP INPUT DATA ONTO MODULE                                          ==
!     ==========================================================================
      EPWRHO=EPWRHO_
      NR1START=NR1START_
      NR1L=NR1L_
      NR2=NR2_
      NR3=NR3_
!     
!     ==========================================================================
!     ==  INITIALIZE PLANEWAVE OBJECT                                         ==
!     ==========================================================================
      CALL CELL$GETR8A('TREF',9,RBAS)
      CALL PLANEWAVE$INITIALIZE('DENSITY','MONOMER',RBAS,K,.TRUE.,EPWRHO &
     &                         ,NR1START,NR1L,NR2,NR3)
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$SETL4('SUPER',.TRUE.)
      CALL PLANEWAVE$GETI4('NGL',NGL)
      CALL SETUP$GETI4('LMRXX',LMRXX) !FORMER CALL SETUP$LMRXX(LMRXX)
      CALL SETUP$GETI4('NSP',NSP)     !FORMER CALL SETUP$NSPECIES(NSP)
      ALLOCATE(LMRX(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LMRX',LMRX(ISP)) !FORMERCALL SETUP$LMRX(ISP,LMRX(ISP))
        CALL SETUP$UNSELECT()
      ENDDO
      CALL POTENTIAL_UPDATE(RBAS)
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
          CALL SETUP$UNSELECT()
      ENDDO
      DEALLOCATE(G2)
      DEALLOCATE(GVEC)
      CALL PLANEWAVE$SELECT(PLANEWAVESELECTION)
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL$VOFRHO(NRL,NDIMD,RHO,LMRXX_,NAT_,QLM,VQLM &
     &                           ,R0,FORCE,RBAS,STRESS,RHOB)
!     **************************************************************************
!     **  MAIN INTERFACE FOR POTENTIAL OBJECT                                 **
!     **  CALCULATE POTENTIAL FROM DENSITY                                    **
!     **************************************************************************
      USE POTENTIAL_MODULE
      USE MPE_MODULE
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
      INTEGER(4)               :: IR
      INTEGER(4)               :: NAT
      INTEGER(4)               :: NSPIN
      INTEGER(4)               :: ISVAR
      LOGICAL(4)               :: TGRA
      REAL(8)   ,ALLOCATABLE   :: RHOTEMP(:,:)
      REAL(8)                  :: SVAR
      REAL(8)   ,ALLOCATABLE   :: VEXT(:)
      REAL(8)   ,ALLOCATABLE   :: FORCEEXT (:,:), STRESSEXT(:,:)
      REAL(8)                  :: EEXT
!     **************************************************************************
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
!     == COLLECT FROM ATOMLIST =================================================
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
!     == PLOT DENSITY ==========================================================
!     IF(TDENSITYTODX(1).AND.TDENSITYTODX(2)) THEN
!       CALL FILEHANDLER$UNIT('DENSITY.DX',NFIL)
!       CALL  WRITEDX_DENSITY
!    &       (NFIL,NR1GLOB,NR2,NR3,NR1GLOB,NR2,NR3,RBAS,RHOE)
!       CALL FILEHANDLER$CLOSE('DENSITY.DX')
!       TDENSITYTODX(2)=.FALSE.
!     END IF
!     
!     == CONVERT DENSITY INTO POTENTIAL ========================================
      CALL SETUP$GETI4('NSP',NSP) !FORMER CALL SETUP$NSPECIES(NSP)

!     == COLLECT FROM PLANEWAVE OBJECT =========================================
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
!     == COLLECT FROM DFT OBJECT ===============================================
      CALL DFT$GETL4('GC',TGRA)
!
!     ==========================================================================
!     == CALCULATE CONFINING POTENTIAL                                        ==
!     ==========================================================================
      IF(TCONFINE) THEN
         ALLOCATE(VEXT(NRL))
         ALLOCATE(FORCEEXT(3,NAT))
         ALLOCATE(STRESSEXT(3,3))
         CALL POTENTIAL_CONFINE(NAT, NR1GLOB, NR1START, NR1L, NR2, NR3, &
     &                          VCONFINE, R0, RBAS, RHO, &
     &                          EEXT, FORCEEXT, STRESSEXT, VEXT)
         CALL MPE$COMBINE('MONOMER','+',EEXT)
         CALL ENERGYLIST$SET('SOLVENT PAULI REPULSION', EEXT)
         CALL ENERGYLIST$ADD('TOTAL ENERGY', EEXT)
      END IF 
!
!     ==========================================================================
!     == POTENTIAL; BRANCH INTO COLLINEAR AND NON-COLLINEAR CASE              ==
!     ==========================================================================
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
!IF(TSTRESS) THEN
!  WRITE(*,FMT='("C+XC STRESS ",3F15.7)')STRESS(1,:)
!  WRITE(*,FMT='("C+XC STRESS ",3F15.7)')STRESS(2,:)
!  WRITE(*,FMT='("C+XC STRESS ",3F15.7)')STRESS(3,:)
!END IF
!
!     ==========================================================================
!     == ADD EXTERNAL POTENTIAL                                               ==
!     ==========================================================================
      IF(TCONFINE) THEN
         RHO(:,1) = RHO(:,1) + VEXT(:)
         FORCE  = FORCE  + FORCEEXT
         STRESS = STRESS + STRESSEXT
         DEALLOCATE(VEXT)
         DEALLOCATE(FORCEEXT)
         DEALLOCATE(STRESSEXT)
      END IF 
!
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(ISPECIES)
      DEALLOCATE(G2)
      DEALLOCATE(GVEC)
                                CALL TIMING$CLOCKOFF('POTENTIAL')
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL_VOFRHO(LMRXX,NRL &
     &                    ,NSP,NAT,ISPECIES,TAU0,FION &
     &                    ,NR1GLOB,NR1,NR2,NR3,RHOE,NDIMD,RBAS &
     &         ,PSCORG,DPSCORG,VBARG,DVBARG,YLMOFG,G0,V0,QLM,VQLM,LMRX &
     &                    ,NGL,GVEC,G2,RHOB,TSTRESS,DG0,DV0,STRESS)
!     **************************************************************************
!     **                                                                      **
!     **  THE PLANE WAVE CONTRIBUTION TO TOTAL ENERGY AND FORCES              **
!     **                                                                      **
!     **  INPUT:                                                              **
!     **    RHOE           CHARGE DENSITY IN REAL SPACE                       **
!     **                   (SPIN UP AND SPIN DOWN                             **
!     **                                                                      **
!     **  OUTPUT:                                                             **
!     **    RHOE           PLANE WAVE POTENTIAL IN REAL SPACE                 **
!     **                   (SPIN UP AND SPIN DOWN)                            **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
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
      REAL(8)   ,INTENT(IN)   :: TAU0(3,NAT)     !<-
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
      REAL(8)   ,INTENT(IN)   :: QLM(LMRXX,NAT)   !<-
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
      REAL(8)                 :: CELLVOL
      REAL(8)                 :: SVAR
      INTEGER(4)              :: NNR,LRXX,IAT,ISPIN,I,LM
      INTEGER(4)              :: IG,ISP
      REAL(8)                 :: STRESST(3,3),STRESS1(3,3)
!     ==  VARIABLES FOR SELF TEST ======================================
      LOGICAL(4),PARAMETER    :: TTEST=.FALSE.
      INTEGER(4),PARAMETER    :: ITEST=2
      COMPLEX(8),ALLOCATABLE  :: RHO_SELFTEST(:,:)
      REAL(8)   ,ALLOCATABLE  :: RHELP(:)
      LOGICAL(4)              :: TBACK
      INTEGER(4)              :: NGAMMA
      LOGICAL(4)              :: TOPTIC
!     **************************************************************************
!      CALL OPTICS$GETL4('ON',TOPTIC)
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
!     == HAND PSEUDO DENSITY OF VALENCE ELECTRONS TO OPTICS CODE      ==
!     ==================================================================
!      IF(TOPTIC)CALL OPTICS$VXCG_VALENCE(TGRA,NSPIN,NRL,NRL,NR1GLOB*NR2*NR3 &
!     &                        ,NR1,NR2,NR3,NGL,CELLVOL,RHOE)
!
!     ==================================================================
!     == TRANSFORM PSEUDO DENSITY INTO FOURIER SPACE                  ==
!     ==================================================================
!CALL PLANEWAVE$GETI4('NGAMMA',NGAMMA)
!PRINT*,'NSPIN,NGL,NRL,CELLVOL,NGAMMA ',NSPIN,NGL,NRL,CELLVOL,NGAMMA
!SVAR=0.D0
!DO IR=1,NRL
!  SVAR=SVAR+RHOE(IR,1)
!ENDDO
!SVAR=SVAR !*CELLVOL
!CALL MPE$COMBINE('MONOMER','+',SVAR)
!NGAMMA=NRL
!CALL MPE$COMBINE('MONOMER','+',NGAMMA)
!SVAR=SVAR/REAL(NGAMMA,KIND=8)
!PRINT*,'CHARGE BEFORE TRANSFORM ',SVAR

      ALLOCATE(RHOG(NGL,NSPIN))
      CALL PLANEWAVE$SUPFFT('RTOG',NSPIN,NGL,RHOG,NRL,RHOE)

!CALL PLANEWAVE$GETI4('NGAMMA',NGAMMA)
!IF(NGAMMA.NE.0)PRINT*,'RHOG(1) A ',RHOG(NGAMMA,1) !*CELLVOL
!STOP
!
!     ==================================================================
!     ==  ADD PSEUDO CORE DENSITY                                     ==
!     ==================================================================
      CALL POTENTIAL_ADDCORE(NGL,NSP,NAT,ISPECIES,TAU0,PSCORG,RHOG(1,1))
!
!     ==================================================================
!     == HAND TOTAL PSEUDO DENSITY TO OPTICS CODE                     ==
!     ==================================================================
!      IF(TOPTIC)CALL OPTICS$VXCG_CORE(TGRA,NSPIN,NRL,NRL,NR1GLOB*NR2*NR3 &
!     &                     ,NR1,NR2,NR3,NGL,CELLVOL,RHOG)
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
      CALL MPE$COMBINE('MONOMER','+',RHOBT)
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
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('RCSM',RCSM(ISP))
        CALL SETUP$GETR8('RCBG',RCBG(ISP))
        CALL SETUP$ISELECT(0)
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
!
!     ==================================================================
!     ==  ISOLATE CHARGE DENSITY FROM PERIODIC IMAGES                 ==
!     ==================================================================
!     CALL DIPOLE(RBAS,NGL,NGL,GVEC,RHOG(1,1),NAX,NSP,NA,TAU0,LMRXX,LMRX,QLM)
                                CALL TIMING$CLOCKON('VOFRHO: ISOLATE')
      EISOLATE=0.D0
      FORCE1(:,:)=0.D0
      VQLM1(:,:)=0.D0
!     __ ISOLATE MAY SET RHOBT TO ZERO _________________________________
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
      CALL GRAPHICS$SETPWPOT('HARTREE',NGL,VHARTREE)
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
!     == SEND POTENTIAL TO OPTICS CODE                                ==
!     ==================================================================
!     IF(TOPTIC)CALL OPTICS$VXCG(NR1,NR2,NR3,NGL,RHOG,NSPIN)
!
!     ==================================================================
!     == ADD HARTREE POTENTIAL TO EXCHANGE POTENTIAL                  ==
!     ==================================================================
      DO IG=1,NGL
         RHOG(IG,1)=RHOG(IG,1)+VHARTREE(IG)
      ENDDO
      DEALLOCATE(VHARTREE)
!
!     == PASS COMPLETE EFFECTIVE POTENTIAL TO THE GRAPHICS MODULE=======
      CALL GRAPHICS$SETPWPOT('TOT',NGL,RHOG(1,1))
!
!     ==================================================================
!     == CALCULATE FORCE ON THE PSEUDO CORE                           ==
!     ==================================================================
      CALL TIMING$CLOCKON('VOFRHO: FORCE ON PS-CORE')
      CALL POTENTIAL_FPSCORE(NSP,NAT,ISPECIES,RBAS,TAU0,FORCE1,STRESS1 &
     &              ,NGL,G2,GVEC,RHOG(1,1),PSCORG,DPSCORG)
      FIONT=FIONT+FORCE1
      STRESST=STRESST+STRESS1
      CALL TIMING$CLOCKOFF('VOFRHO: FORCE ON PS-CORE')
!
!     ==================================================================
!     == SEND POTENTIAL TO OPTICS CODE                                ==
!     ==================================================================
!     IF(TOPTIC)CALL OPTICS$VOFG(NR1,NR2,NR3,NGL,RHOG,NSPIN)
!
!     ==================================================================
!     ==  COMBINE ARRAYS FOR PARALELL PROCESSING                      ==
!     ==================================================================
!     ++ PARALLEL START ++++++++++++++++++++++++++++++++++++++++++++++++
      CALL MPE$COMBINE('MONOMER','+',EHARTREE)
      CALL MPE$COMBINE('MONOMER','+',FIONT)
      CALL MPE$COMBINE('MONOMER','+',VQLMT)
      CALL MPE$COMBINE('MONOMER','+',STRESST)
      CALL MPE$COMBINE('MONOMER','+',EXC)
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
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EXC)
      CALL ENERGYLIST$ADD('AE  EXCHANGE-CORRELATION',EXC)
      CALL ENERGYLIST$ADD('PS  EXCHANGE-CORRELATION',EXC)
!   
!     == EISOLATE CONTAINS THE ENERGY CONTRIBUTIONS FROM ISOLATE,      ==
!     == COSMO, QMMM ETC. THE CORRESPONDING ENERGYLIST CONTRIBUTIONS   ==
!     == ARE ALREADY TAKEN ARE OFF IN ISOLATE_INTERFACE                ==
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
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER    :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                 :: GBAS(3,3)
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)    !(NGL)
      INTEGER(4),PARAMETER    :: NBLOCK=100
      REAL(8)   ,PARAMETER    :: R8SMALL=1.D-20
      REAL(8)                 :: FAC,SVAR,SVAR1
      REAL(8)                 :: CELLVOL
      INTEGER(4)              :: IG,ISP,IAT,I,J
      REAL(8)                 :: FORCE1(3)
      REAL(8)                 :: STRESS1(3,3)
      COMPLEX(8)              :: CSVAR
!     ******************************************************************
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
          SVAR1=-REAL(CONJG(VTEMP(IG))*DPSCORG(IG,ISP)*EIGR(IG),KIND=8)/(G2(IG)+R8SMALL)
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
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER    :: Y0=1.D0/SQRT(4.D0*PI)
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)        !(NGL) STRUCTURE FACTOR
      INTEGER(4)              :: IAT,ISP,IG
!     ******************************************************************
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
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER    :: FPI=4.D0*PI
      REAL(8)   ,PARAMETER    :: Y0=1.D0/SQRT(4.D0*PI)
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)         !(NGL)
      COMPLEX(8),ALLOCATABLE  :: RHO2(:)         !(NGL) PSRHO+COMPENSATION
      COMPLEX(8),ALLOCATABLE  :: VG(:)           !(NGL)
      COMPLEX(8)              :: CFAC(LMRXX)
      INTEGER(4)              :: NGAMMA  
      COMPLEX(8)              :: CSVAR,CSVAR0,AVAL,BVAL,CVAL,DVAL,CSUM
      REAL(8)                 :: VB               !VHAT(GAMMA)
      INTEGER(4)              :: M,L,LM,ISP,IAT,IG,I,J,LM1,LM2
      REAL(8)                 :: SVAR,SVAR1
      REAL(8)                 :: CC(3,3)
      COMPLEX(8)              :: P0(3,3,LMRXX),PM(3,3,LMRXX),PT(3,3)
      REAL(8)   ,PARAMETER    :: R8SMALL=1.D-20
!      REAL(8)                 :: STRESSA(3,3),STRESSB(3,3),STRESSC(3,3)
!     ******************************************************************
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
        RHOB=-REAL(RHO2(NGAMMA),KIND=8)
!!$PRINT*,'GWEIGHT',GWEIGHT,'RHOB',RHOB
!!$PRINT*,'RHOB  FROM VOFRHO_HARTREE',RHOB*GWEIGHT
!!$PRINT*,'Q(RHO) FROM VOFRHO_HARTREE',RHO(NGAMMA)*GWEIGHT
!!$PRINT*,'Q(RHO2) FROM VOFRHO_HARTREE',RHO2(NGAMMA)*GWEIGHT
        RHO2(NGAMMA)=(0.D0,0.D0)   !
      ELSE 
        RHOB=0.D0
      END IF
!
!     ==================================================================
!     ==   CALCULATE ENERGY                                           ==
!     ==================================================================
      CSUM=(0.D0,0.D0)
!!$CSVAR=0.D0
      DO IG=1,NGL
        IF(IG.NE.NGAMMA) THEN
          VG(IG)=RHO2(IG)*FPI/G2(IG)
          CSUM=CSUM+0.5D0*CONJG(RHO2(IG))*VG(IG)+CONJG(RHO(IG))*POT(IG)
!!$CSVAR=CSVAR+0.5D0*CONJG(RHO2(IG))*VG(IG)
          POT(IG)=POT(IG)+VG(IG)
        ELSE 
          VG(IG)=-VB
          CSUM=CSUM+0.5D0*( CONJG(RHO(IG))*POT(IG)+RHOB*VB )
          POT(IG)=POT(IG)+VG(IG)
        END IF
      ENDDO
      EHARTREE=REAL(CSUM,KIND=8)
!!$PRINT*,'PURE HARTREE ENERGY ',REAL(CSVAR,KIND=8)*2.D0*GWEIGHT
!!$PRINT*,'RHO*POT ',REAL(CSUM-CSVAR,KIND=8)*2.D0*GWEIGHT
!!$PRINT*,'STRESS PREDICTED ',-REAL(CSVAR,KIND=8)*2.D0*GWEIGHT/3.D0-REAL(CSUM-CSVAR,KIND=8)*2.D0*GWEIGHT
!
!     == NOW TAKE CARE OF STRESSES =====================================
!     == INTEGRATION WEIGHT IS DONE LATER ==============================
      STRESS(:,:)=0.D0
      IF(TSTRESS) THEN
        DO IG=1,NGL
          IF(IG.NE.NGAMMA) THEN         
            SVAR=REAL(VG(IG)*CONJG(RHO2(IG)))
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
            VQLM(1,IAT)=VQLM(1,IAT)-0.5D0*REAL(CSVAR,KIND=8)
          END IF
!         == MULTIPOLE POTENTIALS ======================================
          AVAL=( G0(IG,ISP)*CONJG(VG(IG)) + V0(IG,ISP)*CONJG(RHO2(IG)) )
          AVAL=AVAL*EIGR(IG)
          CSUM=0.D0
          DO LM=1,LMRX(ISP)
            CSVAR=CFAC(LM)*YLMOFG(LM,IG)
            VQLM(LM,IAT)=VQLM(LM,IAT)+REAL(CSVAR*AVAL,KIND=8)
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
! CAUTION PT IS NOT EXACTLY SYMMETRIC
            DO LM=1,LMRX(ISP)
              PT(:,:)=PT(:,:)+(P0(:,:,LM)+PM(:,:,LM)*G2(IG))*YLMOFG(LM,IG) 
            ENDDO
            BVAL=DG0(IG,ISP)*CONJG(VG(IG))+DV0(IG,ISP)*CONJG(RHO2(IG)) 
            BVAL=BVAL*EIGR(IG)
            DVAL=DVBARG(IG,ISP)*Y0*EIGR(IG)*CONJG(RHO2(IG))
            SVAR=REAL(DVAL+BVAL*CSUM,KIND=8)/(G2(IG)+R8SMALL)
            DO I=1,3
              DO J=1,3
                STRESS(I,J)=STRESS(I,J)- REAL(PT(I,J)*AVAL,KIND=8) &
     &                                 - GVEC(I,IG)*GVEC(J,IG)*SVAR
!
!                STRESSA(I,J)=STRESSA(I,J)- REAL(PT(I,J)*AVAL,KIND=8) 
!                STRESSB(I,J)=STRESSB(I,J)- GVEC(I,IG)*GVEC(J,IG)*REAL(DVAL,KIND=8)/(G2(IG)+R8SMALL)
!                STRESSC(I,J)=STRESSC(I,J)- GVEC(I,IG)*GVEC(J,IG)*REAL(BVAL*CSUM,KIND=8)/(G2(IG)+R8SMALL)
              ENDDO
            ENDDO
            PT(:,:)=0.D0  ! RESET FOR NEXT G-VECTOR
          END IF
        ENDDO
      ENDDO
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
!       __ SAFEGUARD AGAINST NEGATIVE DENSITIES DUE TO FOURIER TRANSFORM________
        RHOT=MAX(RHOT,1.D-8)

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
      INTEGER(4),PARAMETER  :: LMXXX=81
      INTEGER(4),INTENT(IN) :: LMXX
      INTEGER(4),INTENT(IN) :: LMX
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(IN) :: GVEC(3,NG)
      REAL(8)   ,INTENT(OUT):: YLMOFG(LMXX,NG)
      REAL(8)               :: YLM(LMXXX)
      INTEGER(4)            :: LX
      INTEGER(4)            :: LM,IG,L,M
      REAL(8)               :: GTOL,GLENG
      REAL(8)               :: Y0
!     ******************************************************************
      Y0=1.D0/SQRT(16.D0*ATAN(1.D0))  ! =1/SQRT(4*PI)
      LX=INT(SQRT(REAL(LMX))-1.D0)
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$MSG('LMX MUST CONTAIN ALL MAGNETIC QUANTUM NUMBER')
        CALL ERROR$MSG('FOR EACH MAIN ANGULAR MOMENTUM')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('POTENTIAL_YLMOFG')
      END IF
      IF(LMXX.GT.LMXXX) THEN
        CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED')
        CALL ERROR$MSG('LMXX DIFFERS FROM LMXXX')
        CALL ERROR$I4VAL('LMXXX',LMXXX)
        CALL ERROR$I4VAL('LMXX',LMXX)
        CALL ERROR$MSG('INCREASE HARDWIRED LIMIT IN PAW_POTENTIAL.F90')
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
        GLENG=SQRT(GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2)
        IF(GLENG.LT.1.D-7) THEN
          YLMOFG(:,IG)=0.D0
!         == WITH NSP=0, YLMOFG(1,NG) DOES NOT EXIST ===========================
          IF(SIZE(YLMOFG).NE.0)YLMOFG(1,IG)=Y0
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIAL_CONFINE(NAT,NR1G,NR1START,NR1L,NR2,NR3,V0, &
     &                             RAT,RBAS,RHO,ETOT,FORCE,STRESS,POT)
!     **************************************************************************
!     **                                                                      **
!     ** SUBROUTINE FOR CALCULATION OF AN EMPIRICAL POTENTIAL FOR COSMO       **
!     ** CALCULATIONS THAT DESCRIBES PAULI REPULSION BETWEEN THE ELECTRONS    **
!     ** OF SOLVENT AND SOLUTE.                                               **
!     **                                                                      **
!     **************************************************************************
      USE PERIODICTABLE_MODULE, ONLY: PERIODICTABLE$GET
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)  :: NAT                 ! #(ATOMS)
      INTEGER(4), INTENT(IN)  :: NR1G                ! GLOBAL #(GRID PLANES)
      INTEGER(4), INTENT(IN)  :: NR1START            ! FIRST LOCAL GRID PLANE
      INTEGER(4), INTENT(IN)  :: NR1L                ! LOCAL #(GRID PLANES)
      INTEGER(4), INTENT(IN)  :: NR2   ! #(DIVISIONS ALONG 2ND LATTICE VECTOR)
      INTEGER(4), INTENT(IN)  :: NR3   ! #(DIVISIONS ALONG 3RD LATTICE VECTOR)
      REAL(8),    INTENT(IN)  :: V0                  ! MAXIMUM POTENTIAL
      REAL(8),    INTENT(IN)  :: RAT(3,NAT)          ! ATOMIC POSITIONS
      REAL(8),    INTENT(IN)  :: RBAS(3, 3)             ! LATTICE VECTORS
      REAL(8),    INTENT(IN)  :: RHO(NR1L, NR2, NR3) ! ELECTRON DENSITY

      REAL(8),    INTENT(OUT) :: ETOT                ! TOTAL ENERGY
      REAL(8),    INTENT(OUT) :: FORCE(3, NAT)           ! FORCE
      REAL(8),    INTENT(OUT) :: STRESS(3, 3)        ! STRESS
      REAL(8),    INTENT(OUT) :: POT(NR1L, NR2, NR3) ! POTENTIAL
      CHARACTER(2)            :: ATOM_NAME           ! TO GET VDW RADIUS
      REAL(8)                 :: R1(NAT), R2(NAT)    ! VDW RADIUS AND 2*R1
      REAL(8)                 :: GRID_POINT_VOLUME   ! VOLUME OF ONE GRID POINT
      REAL(8)                 :: GRIDBAS(3,3)        ! GRID LATTICE VECTORS
      INTEGER(4)              :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      REAL(8)                 :: XP(3)               ! RELATIVE GRID POSITION      
      REAL(8)                 :: RP(3)               ! ABSOLUTE GRID POSITION      
      REAL(8)                 :: DIS                 ! DISTANCE FROM CENTER      
      REAL(8)                 :: ARG                 ! ARGUMENT FOR CUTOFF FUNC.
      REAL(8)                 :: SVAR                ! AUXILARY VARIABLE
      REAL(8)                 :: VEC(3)              !                      
      INTEGER(4)              :: IAT,I1,I2,I3,J1,J2,J3 ! ITERATION VARIABLES
!     **************************************************************************
!     ==========================================================================
!     == DEFAULT VALUES FOR RADII                                             ==
!     ==========================================================================
      DO IAT=1,NAT
        CALL ATOMLIST$GETCH('NAME',IAT,ATOM_NAME)
        CALL PERIODICTABLE$GET(ATOM_NAME,'R(VDW)',R1(IAT))
        R2(IAT)=2.D0*R1(IAT)
      END DO 
      GRIDBAS(:,1)=RBAS(:,1)/REAL(NR1G,KIND=8)
      GRIDBAS(:,2)=RBAS(:,2)/REAL(NR2,KIND=8)
      GRIDBAS(:,3)=RBAS(:,3)/REAL(NR3,KIND=8)
      GRID_POINT_VOLUME=(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) * RBAS(1,1) &
     &                 +(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) * RBAS(2,1) &
     &                 +(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)) * RBAS(3,1)  
      GRID_POINT_VOLUME=ABS(GRID_POINT_VOLUME)/REAL(NR1G*NR2*NR3,KIND=8)

!     =========================================================================
!     ==  CALCULATE POTENTIAL AND TOTAL ENERGY                               ==
!     =========================================================================
      POT(:,:,:)=V0
      DO IAT=1,NAT
        CALL BOXSPH(GRIDBAS,RAT(1,IAT),RAT(2,IAT),RAT(3,IAT),R2(IAT) &
     &             ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!PRINT*,'MINMAX ',MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
        DO I1=MIN1,MAX1
          J1=MODULO(I1,NR1G)+1
          J1=J1-(NR1START-1)
          IF(J1.LT.1) CYCLE
          IF(J1.GT.NR1L) CYCLE
          XP(1)=REAL(I1,KIND=8)
          DO I2=MIN2,MAX2
            J2=MODULO(I2,NR2)+1
            XP(2)=REAL(I2,KIND=8)
            DO I3=MIN3,MAX3
              J3=MODULO(I3,NR3)+1
              IF(POT(J1,J2,J3).EQ.0.D0) CYCLE
              XP(3)=REAL(I3,KIND=8)          
              RP=MATMUL(GRIDBAS,XP)          ! GRID POINT IN ABSOLUTE COORDINATES
              DIS=SQRT(SUM((RP(:)-RAT(:,IAT))**2))
              IF(DIS.GT.R2(IAT)) CYCLE
              IF(DIS.LT.R1(IAT)) THEN
                POT(J1,J2,J3)=0.D0
                CYCLE
              END IF
              ARG=(DIS-R1(IAT))/(R2(IAT)-R1(IAT))
              SVAR=(3.D0-2.D0*ARG)*ARG**2
              POT(J1,J2,J3)=POT(J1,J2,J3)*SVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
      ETOT=SUM(POT*RHO)*GRID_POINT_VOLUME
!
!     =========================================================================
!     == CALCULATE FORCES                                                    ==
!     =========================================================================
      FORCE(:,:)=0.D0
      DO IAT=1,NAT
        CALL BOXSPH(GRIDBAS,RAT(1,IAT),RAT(2,IAT),RAT(3,IAT),R2(IAT) &
     &             ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
        DO I1=MIN1,MAX1
          J1=MODULO(I1,NR1G)
          J1=J1+1-NR1START+1
          IF(J1.LT.1) CYCLE
          IF(J1.GT.NR1L) CYCLE
          XP(1)=REAL(I1,KIND=8)
          DO I2=MIN2,MAX2
            J2=MODULO(I2,NR2)+1
            XP(2)=REAL(I2,KIND=8)
            DO I3=MIN3,MAX3
              J3=MODULO(I3,NR3)+1
              IF(POT(J1,J2,J3).EQ.0.D0) CYCLE
              IF(POT(J1,J2,J3).EQ.V0) CYCLE
              XP(3)=REAL(I3,KIND=8)
              RP=MATMUL(GRIDBAS,XP)              ! GRID POINT IN ABSOLUTE COORDS
              DIS=SQRT(SUM((RP(:)-RAT(:,IAT))**2))
              IF(DIS.GE.R2(IAT)) CYCLE
              IF(DIS.LE.R1(IAT)) CYCLE
              ARG=(DIS-R1(IAT))/(R2(IAT)-R1(IAT))
              IF(ARG.LT.1.D-8) CYCLE ! AVOID OVERFLOW
              SVAR=6.D0*(1.D0-ARG)/(3.D0-2.D0*ARG)/ARG  !DF/F
              SVAR=RHO(J1,J2,J3)*POT(J1,J2,J3)*SVAR/(R2(IAT)-R1(IAT))
              VEC(:)=(RP(:)-RAT(:,IAT))/DIS
              FORCE(:,IAT)=FORCE(:,IAT)-SVAR*VEC(:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
      FORCE=FORCE*GRID_POINT_VOLUME

!     =========================================================================
!     == CALCULATE STRESS                                                    ==
!     =========================================================================
      STRESS = 0D0
      RETURN
      END SUBROUTINE POTENTIAL_CONFINE
!
!     .........................................................................
      SUBROUTINE POTENTIAL_CONFINE_GRIEGER(NAT,NR1G,NR1START,NR1L,NR2,NR3,V0, &
     &                             R,T,RHO,ETOT,F,STRESS,POT)
!     *************************************************************************
!     **                                                                     **
!     ** SUBROUTINE FOR CALCULATION OF AN EMPIRICAL POTENTIAL FOR COSMO      **
!     ** CALCULATIONS THAT DESCRIBES PAULI REPULSION BETWEEN THE ELECTRONS   **
!     ** OF SOLVENT AND SOLUTE.                                              **
!     **                                                                     **
!     *************************************************************************
      USE PERIODICTABLE_MODULE, ONLY: PERIODICTABLE$GET
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)  :: NAT                 ! #(ATOMS)
      INTEGER(4), INTENT(IN)  :: NR1G                ! GLOBAL #(GRID PLANES)
      INTEGER(4), INTENT(IN)  :: NR1START            ! FIRST LOCAL GRID PLANE
      INTEGER(4), INTENT(IN)  :: NR1L                ! LOCAL #(GRID PLANES)
      INTEGER(4), INTENT(IN)  :: NR2   ! #(DIVISIONS ALONG 2ND LATTICE VECTOR)
      INTEGER(4), INTENT(IN)  :: NR3   ! #(DIVISIONS ALONG 3RD LATTICE VECTOR)
      REAL(8),    INTENT(IN)  :: V0                  ! MAXIMUM POTENTIAL
      REAL(8),    INTENT(IN)  :: R(3, NAT)           ! ATOMIC POSITIONS
      REAL(8),    INTENT(IN)  :: T(3, 3)             ! LATTICE VECTORS
      REAL(8),    INTENT(IN)  :: RHO(NR1L, NR2, NR3) ! ELECTRON DENSITY

      REAL(8),    INTENT(OUT) :: ETOT                ! TOTAL ENERGY
      REAL(8),    INTENT(OUT) :: F(3, NAT)           ! FORCE
      REAL(8),    INTENT(OUT) :: STRESS(3, 3)        ! STRESS
      REAL(8),    INTENT(OUT) :: POT(NR1L, NR2, NR3) ! POTENTIAL

!     INTEGER(4)              :: IPIV(3)             ! PIVOT INDICES
!     INTEGER(4)              :: INFO                ! INFO VARIABLE
      INTEGER(4)              :: I, ATOM_NR          ! ITERATION VARIABLES
      CHARACTER(2)            :: ATOM_NAME           ! TO GET VDW RADIUS

!     REAL(8)                 :: T_LU(3, 3)          ! LU FACT. OF T
      REAL(8)                 :: R1(NAT), R2(NAT)    ! VDW RADIUS AND 2*R1
      REAL(8)                 :: R_GRID(3, NAT)      ! ATOMIC GRID POSITIONS
      REAL(8)                 :: LENGTH_T(3)         ! LENGTH OF LATTICE VECS

      REAL(8) :: GRID_POINT_VOLUME          ! VOLUME OF ONE GRID POINT
!     *************************************************************************
!     =========================================================================
!     == DEFAULT VALUES FOR RADII                                            ==
!     =========================================================================
      DO I=1,NAT
        CALL ATOMLIST$GETCH('NAME',I,ATOM_NAME)
        CALL PERIODICTABLE$GET(ATOM_NAME,'R(VDW)',R1(I))
        R2(I)=2.D0*R1(I)
      END DO 

!     =========================================================================
!     == CALCULATE ATOMIC GRID POSITIONS                                     ==
!     =========================================================================
      CALL LIB$MATRIXSOLVER8(3,3,3,T,R_GRID,R)
      R_GRID(1,:) = R_GRID(1,:) * NR1G + 1
      R_GRID(2,:) = R_GRID(2,:) * NR2  + 1
      R_GRID(3,:) = R_GRID(3,:) * NR3  + 1

!     =========================================================================
!     == LENGTH OF LATTICE VECTORS                                           ==
!     =========================================================================
      LENGTH_T(:) = SQRT(T(1, :)**2 + T(2, :)**2 + T(3, :)**2)

!     =========================================================================
!     == VOLUME OF ONE GRID POINT                                            ==
!     =========================================================================
      GRID_POINT_VOLUME = ABS ( &
     &  ((T(2, 2)*T(3, 3) - T(3, 2)*T(2, 3)) * T(1, 1)  + &
     &   (T(3, 2)*T(1, 3) - T(1, 2)*T(3, 3)) * T(2, 1)  + &
     &   (T(1, 2)*T(2, 3) - T(2, 2)*T(1, 3)) * T(3, 1)) / &
     &  (NR1G * NR2 * NR3))
!
!     =========================================================================
!     == CALCULATE FORCES                                                    ==
!     =========================================================================
      DO ATOM_NR=1,NAT
        DO I =1,3
          CALL POTENTIAL_CONFINE_CALCULATE(NAT, NR1G, NR1START, NR1L, NR2, &
     &                   NR3, V0, R1, R2, R_GRID, T, LENGTH_T, ATOM_NR, I, POT)
          F(I, ATOM_NR)=-SUM(POT*RHO)*GRID_POINT_VOLUME
        END DO
      END DO
!
!     =========================================================================
!     == DO ACTUAL POTENTIAL CALCULATION                                     ==
!     =========================================================================
      CALL POTENTIAL_CONFINE_CALCULATE(NAT, NR1G, NR1START, NR1L, NR2, NR3, &
     &                              V0, R1, R2, R_GRID, T, LENGTH_T, 0, 0, POT)
!
!     =========================================================================
!     == CALCULATE ENERGY                                                    ==
!     =========================================================================
      ETOT=SUM(POT*RHO)*GRID_POINT_VOLUME

!     =========================================================================
!     == CALCULATE STRESS                                                    ==
!     =========================================================================
      STRESS = 0D0
      RETURN
      END SUBROUTINE POTENTIAL_CONFINE_GRIEGER

!     .........................................................................
      SUBROUTINE POTENTIAL_CONFINE_CALCULATE(NAT, NR1G, NR1START, NR1L, &
     &   NR2, NR3, V0, R1, R2, R_GRID, T, LENGTH_T, &
     &   DERIVE_ATOM, DERIVE_DIRECTION, POT)
!     *************************************************************************
!     **                                                                     **
!     ** CALCULATE THE POTENTIAL OBTAINED BY SUBROUTINE CONFINE_POTENTIAL,   **
!     ** AS WELL AS ITS SPATIAL DERIVATIVES WITH RESPECT TO ATOM DERIVE_ATOM **
!     ** AND SPATIAL DIRECTION DERIVE_DIRECTION.                             **
!     ** USE DERIVE_ATOM == 0 IF YOU DO NOT WANT DERIVATIVES!                **
!     **                                                                     **
!     *************************************************************************
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)  :: NAT                 ! #(ATOMS)
      INTEGER(4), INTENT(IN)  :: NR1G                ! GLOBAL #(GRID PLANES)
      INTEGER(4), INTENT(IN)  :: NR1START            ! FIRST LOCAL GRID PLANE
      INTEGER(4), INTENT(IN)  :: NR1L                ! LOCAL #(GRID PLANES)
      INTEGER(4), INTENT(IN)  :: NR2   ! #(DIVISIONS ALONG 2ND LATTICE VECTOR)
      INTEGER(4), INTENT(IN)  :: NR3   ! #(DIVISIONS ALONG 3RD LATTICE VECTOR)

      REAL(8),    INTENT(IN)  :: V0                  ! MAXIMUM POTENTIAL
      REAL(8),    INTENT(IN)  :: R1 (NAT)            ! VAN-DER-WAALS RADIUS
      REAL(8),    INTENT(IN)  :: R2 (NAT)            ! R2
      REAL(8)                 :: R12(NAT)            ! R2 - R1

      REAL(8),    INTENT(IN)  :: R_GRID(3, NAT)      ! ATOMIC GRID POSITIONS
      REAL(8),    INTENT(IN)  :: T(3, 3)             ! LATTICE VECTORS
      REAL(8),    INTENT(IN)  :: LENGTH_T(3)         ! LENGTH OF LATTICE VECS

      INTEGER(4), INTENT(IN)  :: DERIVE_ATOM         ! ATOM NR.  TO DERIVE
      INTEGER(4), INTENT(IN)  :: DERIVE_DIRECTION    ! DIRECTION TO DERIVE

      REAL(8),    INTENT(OUT) :: POT(NR1L, NR2, NR3) ! POTENTIAL

      INTEGER(4)              :: I, J, K, L, ATOM_NR ! ITERATION VARIABLES
      INTEGER(4)              :: II, JJ, KK, ATOM_NR2! DITTO
      INTEGER(4)              :: XMIN, XMAX          ! BARRIERS
      INTEGER(4)              :: YMIN, YMAX          ! DITTO
      INTEGER(4)              :: ZMIN(3), ZMAX(3)    ! DITTO
      REAL(8)                 :: C, D, D1, D2, Y     ! DISTANCES TO ATOMIC POS.
      REAL(8)                 :: P1, P2, P3          ! AXIS INTERCEPTS
      REAL(8)                 :: T_UNIT(3, 3)        ! T AS UNIT VECTORS

      REAL(8) :: COS12, COS13, COS23    ! ANGLES BETWEEN LATTICE VECTORS
      REAL(8) :: DX, DY, DZ             ! GRID INCREMENT: LENGTH_T(X) / NRX
!     *************************************************************************
!     =========================================================================
!     == INITIALIZATIONS                                                     ==
!     =========================================================================
      R12 = R2 - R1
      IF (DERIVE_ATOM .EQ. 0) THEN
         POT=V0
      ELSE
         POT=0.D0
      END IF
      COS12 = DOT_PRODUCT (T(:,1),T(:,2)) / (LENGTH_T(1) * LENGTH_T(2))
      COS13 = DOT_PRODUCT (T(:,1),T(:,3)) / (LENGTH_T(1) * LENGTH_T(3))
      COS23 = DOT_PRODUCT (T(:,2),T(:,3)) / (LENGTH_T(2) * LENGTH_T(3))
      DX = LENGTH_T(1) / REAL(NR1G,KIND=8)
      DY = LENGTH_T(2) / REAL(NR2,KIND=8)
      DZ = LENGTH_T(3) / REAL(NR3,KIND=8)
      DO I = 1, 3
        T_UNIT(:, I) = T(:, I) / LENGTH_T(I)
      END DO

!     =========================================================================
!     == ITERATION OVER GRID POINTS FOR D < R2                               ==
!     =========================================================================
!     == ATOM_NR2 .EQ. 0 IS ADDED FOR THE ATOM TO BE DERIVED ==
ATOMS: DO ATOM_NR2 = 0, NAT
         IF (ATOM_NR2.EQ.DERIVE_ATOM) THEN
            CYCLE ATOMS
         ELSE IF (ATOM_NR2.EQ.0) THEN
            ATOM_NR = DERIVE_ATOM
         ELSE
            ATOM_NR = ATOM_NR2
         END IF
         XMIN = CEILING(R_GRID(1,ATOM_NR) - R2(ATOM_NR)/DX )
         XMAX =   FLOOR(R_GRID(1,ATOM_NR) + R2(ATOM_NR)/DX )
         II = MODULO (XMIN, NR1G) - NR1START
         IF (II.EQ.-NR1START) II = NR1G - NR1START

ILOOP:   DO I = XMIN, XMAX
            IF (II .GT. NR1G - NR1START) II = 1 - NR1START
            II = II + 1
            IF (II .LE. 0 .OR. II .GT. NR1L) CYCLE ILOOP
            P1 = DBLE(I - R_GRID(1, ATOM_NR)) * DX
            D2 = SQRT(R2(ATOM_NR)**2 + P1**2 * (COS12**2 - 1)) / DY
            C  = R_GRID(2, ATOM_NR) - P1 * COS12 / DY
            YMIN = CEILING(C - D2)
            YMAX = FLOOR(C + D2)

            JJ = MODULO (YMIN, NR2)
            IF (JJ .EQ. 0) JJ = NR2

JLOOP:   DO J = YMIN, YMAX
            P2 = DBLE(J - R_GRID(2, ATOM_NR)) * DY

            D2 = SQRT(R2(ATOM_NR)**2 + P1**2 * (COS13**2 - 1) &
     &         + P2**2 * (COS23**2 - 1) + 2*P1*P2 * (COS13*COS23 - COS12)) / DZ
            C = R_GRID(3, ATOM_NR) - (P1*COS13 + P2*COS23) / DZ
            ZMIN(1) = CEILING(C - D2)

            D1 = R1(ATOM_NR)**2 + P1**2 * (COS13**2 - 1) &
     &          + P2**2 * (COS23**2 - 1) + 2*P1*P2 * (COS13*COS23 - COS12)
            IF (D1 .LT. 0) THEN
               ZMAX(1) = FLOOR(C + D2)
               ZMIN(3) = 1
               ZMAX(3) = 0
               ZMIN(2) = 1
               ZMAX(2) = 0
            ELSE
               D1 = SQRT(D1) / DZ
               ZMAX(1) = FLOOR(C - D1)
               ZMIN(2) = CEILING(C + D1)
               ZMAX(2) = FLOOR(C + D2)
               ZMIN(3) = ZMAX(1) + 1
               ZMAX(3) = ZMIN(2) - 1
            END IF

KPRELOOP:DO L = 1, 2
            KK = MODULO (ZMIN(L), NR3)
            IF (KK .EQ. 0) KK = NR3

KLOOP:   DO K = ZMIN(L), ZMAX(L)
            P3 = DBLE(K - R_GRID(3, ATOM_NR)) * DZ

            D = SQRT (P1**2 + P2**2 + P3**2 &
     &              + 2*P1*P2*COS12 + 2*P1*P3*COS13 + 2*P2*P3*COS23)
            Y = (D - R1(ATOM_NR)) / R12(ATOM_NR)

            IF (ATOM_NR2 .NE. 0) THEN ! 0 MEANS: DERIVE
               POT(II, JJ, KK) = &
     &         POT(II, JJ, KK) * (3D0 * Y**2 - 2D0 * Y**3)
            ELSE
               POT(II, JJ, KK) = - (Y - Y**2) * (6D0 / (R12(ATOM_NR) * D)) * &
     &            (P1 * T_UNIT(DERIVE_DIRECTION, 1) + &
     &             P2 * T_UNIT(DERIVE_DIRECTION, 2) + &
     &             P3 * T_UNIT(DERIVE_DIRECTION, 3)) * V0
            END IF

         KK = KK + 1
         IF (KK .GT. NR3) KK = 1
         END DO KLOOP
         END DO KPRELOOP

         IF (ATOM_NR2 .NE. 0) THEN ! IF ZERO, VALUES ARE INITIALIZED TO ZERO
            KK = MODULO (ZMIN(3), NR3)
            IF (KK .EQ. 0) KK = NR3

K0LOOP:     DO K = ZMIN(3), ZMAX(3)
               POT(II, JJ, KK) = 0D0
               KK = KK + 1
               IF (KK .GT. NR3) KK = 1
            END DO K0LOOP
         END IF

         JJ = JJ + 1
         IF (JJ .GT. NR2) JJ = 1
         END DO JLOOP

         END DO ILOOP
      END DO ATOMS
      RETURN
      END SUBROUTINE POTENTIAL_CONFINE_CALCULATE

!!$!
!!$!     ..................................................................
!!$      SUBROUTINE POTENTIAL_DENSITYINTERFACE(NR1START,NR1L,NR1G,NR2,NR3 &
!!$     &                         ,NRL,NDIMD,RHO,POT,R0,FORCE,RBAS,STRESS)
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN)    :: NR1START
!!$      INTEGER(4),INTENT(IN)    :: NR1L
!!$      INTEGER(4),INTENT(IN)    :: NR1G
!!$      INTEGER(4),INTENT(IN)    :: NR2
!!$      INTEGER(4),INTENT(IN)    :: NR3
!!$      INTEGER(4),INTENT(IN)    :: NRL
!!$      INTEGER(4),INTENT(IN)    :: NDIMD
!!$      INTEGER(4),INTENT(IN)    :: NAT
!!$      REAL(8)   ,INTENT(IN)    :: RHO(NRL,NDIMD)
!!$      REAL(8)   ,INTENT(OUT)   :: POT(NRL,NDIMD)
!!$      REAL(8)   ,INTENT(IN)    :: R0(3,NAT)
!!$      REAL(8)   ,INTENT(OUT)   :: FORCE(3,NAT)
!!$      REAL(8)   ,INTENT(IN)    :: RBAS(3,3)
!!$      REAL(8)   ,INTENT(OUT)   :: STRESS(3,3)    
!!$      REAL(8)                  :: DT1(3),DT2(3),DT3(3)
!!$      REAL(8)                  :: RGRID(3)
!!$      REAL(8)                  :: VAL
!!$      REAL(8)                  :: DR(3)
!!$      REAL(8)                  :: F1(3,NAT)
!!$      INTEGER(4)               :: I1,I2,I3
!!$      LOGICAL(4)               :: TZERO,TONE
!!$      REAL(8)                  :: RVDW(NAT)
!!$      REAL(8)                  :: RVDWSOLV
!!$      REAL(8)                  :: VDW(NAT)
!!$      REAL(8)                  :: RMAX2(NAT)
!!$!     ******************************************************************
!!$      DT1(:)=RBAS(:,1)/REAL(NR1G-1)
!!$      DT2(:)=RBAS(:,2)/REAL(NR2-1)
!!$      DT3(:)=RBAS(:,3)/REAL(NR3-1)
!!$      RMAX2(:)=(RVDW(:)+RVDWSOLV)**2
!!$      POT(:,1)=1.D0
!!$      DO I1=IR1START,IR1START_1+NR1L
!!$        DO I2=1,NR2
!!$          DO I3=1,NR3
!!$            RGRID(:)=DT1(:)*REAL(I1-1)+DT2(:)*REAL(I2-1)+DT3(:)*REAL(I3-1)
!!$            IR=IR1-IR1START+1+NR2*(IR2-1+NR3*(IR3-1))
!!$            VAL=1.D0
!!$            F1(:,:)=1.D0
!!$            TZERO=.FALSE.
!!$            TONE=.TRUE.
!!$            DO IAT=1,NAT
!!$              DR(:)=RGRID(:)-R0(:,IAT)
!!$              DIS2=SUM((DR(:))*2)
!!$              IF(DIS2.GT.RMAX(IAT)) CYCLE
!!$              DIS=SQRT(DIS2)
!!$              IF(DIS.LE.RVDW(IAT)) THEN
!!$                VAL=0.D0
!!$                F1(:,:)=0.D0
!!$                TZERO=.TRUE.
!!$                EXIT
!!$              ELSE 
!!$                X=(DIS-RMIN(IAT))/RVDWSOLV
!!$                SVAR=(3.D0-2.D0*X)*X**2
!!$                SVAR2=6.D0*(1.D0-X)*X/RVDWSOLV/DIS
!!$                VAL=VAL*SVAR
!!$                F1(:,:IAT-1)=F1(:,:IAT-1)*SVAR
!!$                F1(:,IAT)=F1(:,IAT)*DR(:)*SVAR2
!!$                F1(:,IAT+1:)=F1(:,IAT+1:)*SVAR
!!$                TONE=.FALSE.
!!$              END IF
!!$            ENDDO
!!$            IF(TZERO) CYCLE
!!$            IF(TONE) THEN
!!$              POT(IR,1)=1.D0
!!$            ELSE
!!$              POT(IR,1)=VAL
!!$              FORCE(:,:)=FORCE(:,:)+F1(:,:)*RHO(IR,1)
!!$            END IF
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      RETURN
!!$      END


