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
      REAL(8)    ,ALLOCATABLE :: V0(:,:)     !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: G0(:,:)     !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: PSCOREG(:,:)!(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: VBARG(:,:)  !(NGL,NSP)
      REAL(8)    ,ALLOCATABLE :: YLMOFG(:,:) !(LMRXX,NGL)
      INTEGER(4) ,ALLOCATABLE :: LMRX(:)     !(NSP)
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
        CALL ERROR$STOP('POTENTIAL$GVECTORS')
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
      SUBROUTINE POTENTIAL_SET
!     ******************************************************************
!     **  DEFINE GRID FOR FOURIER TRANSFORM                           **
!     ******************************************************************
      USE POTENTIAL_MODULE
      IMPLICIT NONE
      REAL(8)   ,ALLOCATABLE:: GVEC(:,:)
      REAL(8)   ,ALLOCATABLE:: G2(:)
      INTEGER(4)            :: ISP
      REAL(8)               :: RBAS(3,3)
      REAL(8)               :: GBAS(3,3)
      REAL(8)               :: CELLVOL
!     ******************************************************************
                                CALL TRACE$PUSH('POTENTIAL_SET')
      IF(TSET) RETURN
      TSET=.TRUE.
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      ALLOCATE(GVEC(3,NGL))
      ALLOCATE(G2(NGL))
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
      CALL PLANEWAVE$GETR8A('G2',NGL,G2)
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
!
!     ================================================================
!     == ALLOCATE ARRAYS                                            ==
!     ================================================================
      ALLOCATE(YLMOFG(LMRXX,NGL))
      ALLOCATE(PSCOREG(NGL,NSP))
      ALLOCATE(VBARG(NGL,NSP))
      ALLOCATE(G0(NGL,NSP)) ;G0(:,:)=0.D0
      ALLOCATE(V0(NGL,NSP))
!
!     ================================================================
!     == CALCULATE YLM(G)*G**L FOR VOFRHO                           ==
!     ================================================================
      CALL POTENTIAL_YLMOFG(LMRXX,LMRXX,NGL,GVEC,YLMOFG)
!     
!     == CALCULATE PS CORE DENSITY IN G-SPACE=========================
      DO ISP=1,NSP
        CALL SETUP$COMPENSATION(ISP,NGL,G2,CELLVOL,G0(1,ISP),V0(1,ISP))
        CALL SETUP$PSCOREG(ISP,NGL,G2,CELLVOL,PSCOREG(1,ISP))
        CALL SETUP$VBARG(ISP,NGL,G2,CELLVOL,VBARG(1,ISP))
      ENDDO
      DEALLOCATE(G2)
      DEALLOCATE(GVEC)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL$VOFRHO(NRL,NDIMD,RHO,LMRXX_,NAT_,QLM,VQLM,RHOB)
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
      REAL(8)   ,INTENT(OUT)   :: RHOB
      REAL(8)   ,INTENT(INOUT) :: RHO(NRL,NDIMD)
      REAL(8)                  :: RBAS(3,3)
      INTEGER(4),ALLOCATABLE   :: ISPECIES(:) !(NAT) 
      REAL(8)   ,ALLOCATABLE   :: R0(:,:)     !(3,NAT)
      REAL(8)   ,ALLOCATABLE   :: FORCE(:,:)  !(3,NAT)
      REAL(8)   ,ALLOCATABLE   :: FORCE1(:,:)  !(3,NAT)
      REAL(8)   ,ALLOCATABLE   :: G2(:)       !(NGL)
      REAL(8)   ,ALLOCATABLE   :: GVEC(:,:)     !(3,NGL)
      INTEGER(4)               :: IAT
      INTEGER(4)               :: NAT
      INTEGER(4)               :: NSPIN
      INTEGER(4)               :: ISVAR
      LOGICAL(4)               :: TGRA
!     ******************************************************************
                                CALL TRACE$PUSH('POTENTIAL$VOFRHO')
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('POTENTIAL OBJECT IS NOT INITIALIZED')
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
      CALL POTENTIAL_SET
      IF(LMRXX_.NE.LMRXX) THEN
        CALL ERROR$MSG('INPUT DATA INCONSISTENT WITH MODULE')
        CALL ERROR$I4VAL('LMRXX_',LMRXX_)
        CALL ERROR$I4VAL('LMRXX',LMRXX)
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
      CALL ATOMLIST$NATOM(NAT)
      IF(NAT_.NE.NAT) THEN
        CALL ERROR$MSG('INPUT DATA INCONSISTENT WITH MODULE')
        CALL ERROR$I4VAL('NAT_',NAT_)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
  
      CALL ATOMLIST$LATTICE(RBAS)
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
!     == COLLECT FROM ATOMLIST =======================================
      ALLOCATE(ISPECIES(NAT))
      ALLOCATE(R0(3,NAT))    
      ALLOCATE(FORCE(3,NAT)) 
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCE)
      IF(NDIMD.EQ.1) THEN
        NSPIN=1
      ELSE IF(NDIMD.EQ.2) THEN
        NSPIN=2
      ELSE IF(NDIMD.EQ.4) THEN
        CALL ERROR$MSG('NONCOLINEAR DESCRIPTION NOT COMPLETED')
        CALL ERROR$STOP('POTENTIAL$VOFRHO')
      END IF
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
!     == CONVERT DENSITY INTO POTENTIAL                               ==
!     ==================================================================
      CALL POTENTIAL_VOFRHO(LMRXX,NRL,NSP,NAT,ISPECIES,R0,FORCE &
     &                  ,NR1GLOB,NR1L,NR2,NR3,RHO,NSPIN,RBAS &
     &                  ,PSCOREG,VBARG,YLMOFG,G0,V0,QLM,VQLM,LMRX &
     &                  ,NGL,GVEC,G2,RHOB)
!
!     ==================================================================
!     ==  ADD FORCES AND SEND FORCES TO ATOMLIST                      ==
!     ==================================================================
      ALLOCATE(FORCE1(3,nat))
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCE1)
      FORCE(:,:)=FORCE(:,:)+FORCE1(:,:)
      DEALLOCATE(FORCE1)
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
!
!     ==================================================================
!     == CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(ISPECIES)
      DEALLOCATE(R0)      
      DEALLOCATE(FORCE)   
      DEALLOCATE(G2)
      DEALLOCATE(GVEC)
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_VOFRHO(LMRXX,NRL &
     &                    ,NSP,NAT,ISPECIES,TAU0,FION &
     &                    ,NR1GLOB,NR1,NR2,NR3,RHOE,NDIMD,RBAS &
     &                    ,PSCORG,VBARG,YLMOFG,G0,V0,QLM,VQLM,LMRX &
     &                    ,NGL,GVEC,G2,RHOB)
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
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT)
      REAL(8)   ,INTENT(INOUT):: TAU0(3,NAT)     !<-
      REAL(8)   ,INTENT(OUT)  :: FION(3,NAT)
      INTEGER(4),INTENT(IN)   :: NR1GLOB
      INTEGER(4),INTENT(IN)   :: NR1
      INTEGER(4),INTENT(IN)   :: NR2
      INTEGER(4),INTENT(IN)   :: NR3
      REAL(8)   ,INTENT(INOUT):: RHOE(NRL,NDIMD)
      INTEGER(4),INTENT(IN)   :: NDIMD
      REAL(8)   ,INTENT(IN)   :: RBAS(3,3)
      REAL(8)   ,INTENT(IN)   :: PSCORG(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: VBARG(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: YLMOFG(LMRXX,NGL)
      REAL(8)   ,INTENT(IN)   :: G0(NGL,NSP)
      REAL(8)   ,INTENT(IN)   :: V0(NGL,NSP)
      REAL(8)   ,INTENT(INOUT):: QLM(LMRXX,NAT)   !<-
      REAL(8)   ,INTENT(OUT)  :: VQLM(LMRXX,NAT)
      INTEGER(4),INTENT(IN)   :: LMRX(NSP)
      INTEGER(4),INTENT(IN)   :: NGL
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
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      ALLOCATE(RHOG(NGL,NSPIN))
      CALL PLANEWAVE$SUPFFT('RTOG',NSPIN,NGL,RHOG,NRL,RHOE)
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
      FIONT(:,:)=0.D0
      VQLMT(:,:)=0.D0
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
     &            ,CELLVOL,G2,GVEC,YLMOFG,VBARG,G0,V0,RHOBT)
      CALL MPE$BROADCAST(1,RHOBT)
      FIONT(:,:)=FIONT(:,:)+FORCE1(:,:)
      VQLMT(:,:)=VQLMT(:,:)+VQLM1(:,:)
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
     &                          ,TAU0,FORCE1,QLM,VQLM1,EPAIR)
      FION(:,:)=FION(:,:)+FORCE1(:,:)
      VQLM(:,:)=VQLM(:,:)+VQLM1(:,:)
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

!PRINT*,'WARNING IN VOFRHO XC POTENTIAL SWITCHED OFF!!!!!'
!RHOE=0.D0
!GRHO=0.D0
      EXC=0.D0
      CALL POTENTIAL_XC &
    &     (TGRA,NSPIN,NRL,NRL,NR1GLOB*NR2*NR3,CELLVOL,RHOE,GRHO,EXC)
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
      CALL POTENTIAL_FPSCORE(NSP,NAT,ISPECIES,RBAS,TAU0,FIONT &
     &              ,NGL,GVEC,RHOG(1,1),PSCORG)
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
      CALL MPE$COMBINE('+',EXC)
!     ++ PARALLEL END ++++++++++++++++++++++++++++++++++++++++++++++++++
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
      SUBROUTINE POTENTIAL_FPSCORE(NSP,NAT,ISPECIES,RBAS,TAU0,FION &
     &                            ,NGL,GVEC,VTEMP,PSCORG)
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
      REAL(8)   ,INTENT(INOUT):: FION(3,NAT)
      INTEGER(4),INTENT(IN)   :: NGL
      REAL(8)   ,INTENT(IN)   :: GVEC(3,NGL)
      REAL(8)   ,INTENT(IN)   :: PSCORG(NGL,NSP)
      COMPLEX(8),INTENT(IN)   :: VTEMP(NGL)
      REAL(8)                 :: GBAS(3,3)
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)    !(NGL)
      REAL(8)                 :: SUM1,SUM2,SUM3
      INTEGER(4),PARAMETER    :: NBLOCK=100
      REAL(8)                 :: SSUM1,SSUM2,SSUM3
      REAL(8)                 :: FAC,SVAR
      REAL(8)                 :: CELLVOL
      REAL(8)                 :: PI,Y0
      INTEGER(4)              :: IG,ISP,IAT
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      ALLOCATE(EIGR(NGL))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL PLANEWAVE$STRUCTUREFACTOR(TAU0(1,IAT),NGL,EIGR)
        SSUM1=0.D0
        SSUM2=0.D0
        SSUM3=0.D0
        SUM1=0.D0
        SUM2=0.D0
        SUM3=0.D0
        DO IG=1,NGL
          SVAR=-AIMAG(VTEMP(IG)*CONJG(EIGR(IG)*PSCORG(IG,ISP)))
          SUM1=SUM1+SVAR*GVEC(1,IG)
          SUM2=SUM2+SVAR*GVEC(2,IG)
          SUM3=SUM3+SVAR*GVEC(3,IG)
          IF(MOD(IG,NBLOCK).EQ.0) THEN
            SSUM1=SSUM1+SUM1
            SSUM2=SSUM2+SUM2
            SSUM3=SSUM3+SUM3
            SUM1=0.D0
            SUM2=0.D0
            SUM3=0.D0
          END IF
        ENDDO
        SUM1=SUM1+SSUM1
        SUM2=SUM2+SSUM2
        SUM3=SUM3+SSUM3
        FAC=Y0*2.D0*CELLVOL
        FION(1,IAT)=FION(1,IAT)-SUM1*FAC
        FION(2,IAT)=FION(2,IAT)-SUM2*FAC
        FION(3,IAT)=FION(3,IAT)-SUM3*FAC
      ENDDO
      DEALLOCATE(EIGR)
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
          RHO(IG)=RHO(IG)+PSCORG(IG,ISP)*y0*EIGR(IG)
        ENDDO
      ENDDO
      DEALLOCATE(EIGR)
      RETURN 
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_HARTREE(LMRXX,LRXX &
     &            ,NGL,NSP,NAT,ISPECIES,LMRX &
     &            ,RHO1,VHARTREE,QLM,VQLM,TAU0,FION,EHARTREE &
     &            ,CELLVOL,G2,GVEC,YLMOFG,VBARG,G0,V0,RHOB)
!     ******************************************************************
!     **                                                              **
!     ** CALCULATES THE HARTREE ENERGY IN RECIPROCAL SPACE            **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    RHO1       DENSITY OF THE PS WAVE FUNCTIONS (G-SPACE)     **
!     **                                                              **
!     **  OUTPUT:                                                     **
!     **    VHARTREE      HARTREE POTENTIAL (G-SPACE)                    **
!     **    RHO1       DENSITY OF PS WVAE FUNCTIONS AND PS CORE       **
!     **    FION       HARTREE FORCE ON THE ATOMS (ADDED TO INPUT)    **
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
      COMPLEX(8),INTENT(IN)   :: RHO1(NGL)       ! IN PS-DENSITY; OUT: +PS-CORE
      COMPLEX(8),INTENT(OUT)  :: VHARTREE(NGL)   ! HARTREE-POTENTIAL
      REAL(8)   ,INTENT(IN)   :: QLM(LMRXX,NAT)  ! MOMENTS OF 1C-DENSITY
      REAL(8)   ,INTENT(OUT)  :: VQLM(LMRXX,NAT) ! DE/DQLM
      REAL(8)   ,INTENT(IN)   :: TAU0(3,NAT)     ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(OUT)  :: FION(3,NAT)     ! FORCE ON ATOMS (ADDED)
      REAL(8)   ,INTENT(OUT)  :: EHARTREE        ! PS-HARTREE ENERGY (WITHOUT PAIRPOT)
      REAL(8)   ,INTENT(IN)   :: CELLVOL         ! UNIT CELL VOLUME
      REAL(8)   ,INTENT(IN)   :: G2(NGL)         ! G**2
      REAL(8)   ,INTENT(IN)   :: GVEC(3,NGL)       ! G-VECTORS
      REAL(8)   ,INTENT(IN)   :: YLMOFG(LMRXX,NGL)! SPHERICAL HARMONICS
      REAL(8)   ,INTENT(IN)   :: G0(NGL,NSP)     ! WIDE GAUSSIANS OF COMPENSATION DENSITY
      REAL(8)   ,INTENT(IN)   :: V0(NGL,NSP)     ! POTENTIAL OF KOMPENSATION DENSITIES
      REAL(8)   ,INTENT(IN)   :: VBARG(NGL,NSP)  ! V-HAT
      REAL(8)   ,INTENT(OUT)  :: RHOB           ! DENSITY OF THE COMPENSATING BACKGROUND
      COMPLEX(8),ALLOCATABLE  :: EIGR(:)         !(NGL)
      COMPLEX(8),ALLOCATABLE  :: RHO2(:)         !(NGL)
      COMPLEX(8)              :: CFAC(LMRXX)
      INTEGER(4)              :: NGAMMA  
      COMPLEX(8)              :: CSVAR
      COMPLEX(8)              :: CSVAR1
      COMPLEX(8)              :: EIGR1
      real(8)                 :: Vb
      INTEGER(4)              :: M,L,NSTART,LM,ISP,IAT,IG
      REAL(8)                 :: SVAR,FPIBG2
      REAL(8)                 :: PI,Y0,FPI
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FPI=4.D0*PI
      Y0=1.D0/DSQRT(FPI)
!
      ALLOCATE(EIGR(NGL))
      ALLOCATE(RHO2(NGL))
      CALL PLANEWAVE$GETI4('NGAMMA',NGAMMA)
!
!     ==================================================================
!     ==  ADD COMPENSATION CHARGE DENSITY                             ==
!     ==  RHO1 = PS-RHO + PS-CORE                                     ==
!     ==  RHO2 = PS-RHO + RHO-HAT + PS-CORE + RHO-BACKGROUND          ==
!     ==  VHARTREE =                                                     ==
!     ==================================================================
      RHOB=0.D0
      VB=0.D0
      DO IG=1,NGL
        RHO2(IG)=RHO1(IG)
        VHARTREE(IG)=(0.D0,0.D0)
      ENDDO
      IF(NGAMMA.NE.0) THEN
        RHOB=-REAL(RHO1(NGAMMA))
      ENDIF
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
          RHOB=RHOB-real(CFAC(1))*Y0*G0(NGAMMA,ISP)
          VB  =VB  +real(CFAC(1))*Y0*V0(NGAMMA,ISP)
        END if
        DO IG=1,NGL
          EIGR1=EIGR(IG)
          CSVAR=(0.D0,0.D0)
          DO LM=1,LMRX(ISP)
            CSVAR=CSVAR+CFAC(LM)*YLMOFG(LM,IG)
          ENDDO
          CSVAR=CSVAR*EIGR1
          RHO2(IG)    =RHO2(IG)    +CSVAR*G0(IG,ISP)
          VHARTREE(IG)=VHARTREE(IG)+CSVAR*V0(IG,ISP) &
     &                             +VBARG(IG,ISP)*Y0*EIGR1
        ENDDO
      ENDDO
      IF(NGAMMA.NE.0) THEN
!PRINT*,'rhob  FROM VOFRHO_HARTREE',RHOB*CELLVOL
!PRINT*,'Q(rho1) FROM VOFRHO_HARTREE',RHO1(NGAMMA)*CELLVOL
!PRINT*,'Q(rho2) FROM VOFRHO_HARTREE',RHO2(NGAMMA)*CELLVOL
        RHOB=-REAL(RHO2(NGAMMA))
        RHO2(NGAMMA)=(0.D0,0.D0)   !
      END IF
!
!     ==================================================================
!     ==   CALCULATE ENERGY                                           ==
!     ==================================================================
      CSVAR=(0.D0,0.D0)
      DO IG=1,NGL
        IF(IG.NE.NGAMMA) THEN
          FPIBG2=FPI/G2(IG)
          CSVAR1=FPIBG2*RHO2(IG)
          CSVAR=CSVAR+2.D0*CONJG(RHO1(IG))*VHARTREE(IG) &
     &               +CONJG(RHO2(IG))*CSVAR1 
          VHARTREE(IG)=VHARTREE(IG)+CSVAR1
        ELSE 
          RHO2(IG)=(0.D0,0.D0)
          CSVAR=CSVAR+VHARTREE(IG)*CONJG(RHO1(IG))+RHOB*VB
        END IF
      ENDDO
      EHARTREE=REAL(CSVAR)*CELLVOL
!
!     ==================================================================
!     ==   CALCULATE FORCES AND TOTAL ENERGY                          ==
!     ==      VBARG  = V_ADD                                          ==
!     ==      V0 = POTENTIAL FROM THE DIFFERENCE                      ==
!     ==           BETWEEN THE COMPENSATION CHARGES                   ==
!     ==================================================================
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL PLANEWAVE$STRUCTUREFACTOR(TAU0(1,IAT),NGL,EIGR)
!
        LM=0
        CSVAR=CI
        DO L=0,LRXX
          CSVAR=-CSVAR*CI/DBLE(2*L+1)
          DO M=-L,L
            LM=LM+1
            CFAC(LM)=CSVAR
          ENDDO
        ENDDO
        DO IG=1,NGL
          IF(IG.NE.NGAMMA) THEN
            FPIBG2=FPI/G2(IG)
            EIGR1=EIGR(IG)*2.D0*CELLVOL
            CSVAR1 = ( G0(IG,ISP)*FPIBG2*CONJG(RHO2(IG)) &
     &               + V0(IG,ISP)*CONJG(RHO1(IG)))*EIGR1
            CSVAR=(0.D0,0.D0)
            DO LM=1,LMRX(ISP)
              VQLM(LM,IAT)=VQLM(LM,IAT) &
     &                       +CSVAR1*CFAC(LM)*YLMOFG(LM,IG)
              CSVAR  = CSVAR+CFAC(LM)*QLM(LM,IAT)*YLMOFG(LM,IG)
            ENDDO
            SVAR = AIMAG(VBARG(IG,ISP)*Y0*EIGR1*CONJG(RHO1(IG)) &
     &                                         + CSVAR*CSVAR1 )
            FION(1,IAT)=FION(1,IAT)-SVAR*GVEC(1,IG)
            FION(2,IAT)=FION(2,IAT)-SVAR*GVEC(2,IG)
            FION(3,IAT)=FION(3,IAT)-SVAR*GVEC(3,IG)
          ELSE
            CSVAR=EIGR(IG)*CFAC(1)*Y0*CELLVOL
            VQLM(1,IAT)=VQLM(1,IAT) &
     &               +CSVAR*(-VB*G0(IG,ISP)+(RHO1(IG)+RHOB)*V0(IG,ISP))
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(EIGR)
      DEALLOCATE(RHO2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE POTENTIAL_XC &
     &     (TGRA,NSPIN,NNRX,NNR,NNRSCAL,CELLVOL,RHO,GRHO,EXC)
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
      INTEGER(4),PARAMETER    :: NBLOCK=100    ! BLOCKS THE SUMMATION OVER GRID POINTS
      LOGICAL(4)              :: TSPIN
      REAL(8)                 :: SXC,SSXC
      INTEGER(4)              :: IR
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
      REAL(8)                 :: SVAR
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
!       == TAKE CARE OF BLOCKING                                      ==
!       ================================================================
        IF(MOD(IR,NBLOCK).EQ.0) THEN
          SSXC=SSXC+SXC               ! SWAP INTERMEDIATE SUM 
          SXC=0.D0                    ! INTO FINAL SUM
        END IF
      ENDDO
      SSXC=SSXC+SXC  ! CLEAN OUT INTERMEDIATE SUM
      SXC=0.D0       ! THIS LAST STATEMENT IS JUST TO AVOID CONFUSION
!
!     ==================================================================
!     == SCALE EXCHANGE ENERGY WITH THE INTEGRATION WEIGHT            ==
!     ==================================================================
      EXC=SSXC*CELLVOL/DBLE(NNRSCAL)
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
      real(8)               :: y0
!     ******************************************************************
      Y0=1.D0/SQRT(16.D0*DATAN(1.D0))  ! =1/sqrt(4*pi)
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
