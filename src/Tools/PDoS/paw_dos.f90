!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**  NAME: PDOS                                                               **
!**                                                                           **
!**  PURPOSE: ANALYSIS TOOL FOR DENSITY OF STATES                             **
!**                                                                           **
!**  PDOS_MODULE (SEE OBJECT PAW_PDOS.F90 OF THE SIMULATION CODE)             **
!**                                                                           **
!**                                                                           **
!**                                                                           **
!**  NEWORBITAL_MODULE: DESCRIBES A GENERAL ORBITAL.                          **
!**    EACH ENTRY OF AN ORBITAL IS (ATOM-INDEX IAT, LATTICE TRANSLATION IT,   **
!**    ARRAY ORB OF ORBITAL COEFFICIENTS, ARRAY IPRO OF ORBITAL INDICES       **
!**    THAT POINT TO A SPECIFIC POSITION ON THE STATE VECTOR                  **
!**                                                                           **
!**  NEWSET_MODULE: ENCODES A PROJECTION OR COOP                              **
!**    A SET SPECIFIES EXACTLY ONE OF THE FOLLOWING                           **
!**    1) COOP 2)WGHT 3) ORBITALS, 4) SPECIALS: (TOTAL,EMPTY,ALL)             **
!**    THE NEWSET REFERS TO THE NE NEWORBITALS                                **
!**                                                                           **
!**                                                                           **
!**                                                                           **
!*******************************************************************************
!*******************************************************************************
!**
!*******************************************************************************
!**  COPY OF THE PDOS_MODULE FROM THE PAW_PDOS.F90 OBJECT OF THE MAIN CODE    **
!*******************************************************************************
!**  MODULE PDOS_MODULE
!**    TYPE STATE_TYPE
!**      INTEGER(4)                        :: NB
!**      REAL(8)   ,POINTER                :: EIG(:)
!**      REAL(8)   ,POINTER                :: OCC(:)
!**      COMPLEX(8),POINTER                :: VEC(:,:,:)
!**    END TYPE STATE_TYPE 
!**    CHARACTER(6)                        :: FLAG
!**    INTEGER(4)                          :: NAT
!**    INTEGER(4)                          :: NSP
!**    INTEGER(4)                          :: NKPT
!**    INTEGER(4)                          :: NSPIN
!**    INTEGER(4)                          :: NDIM
!**    INTEGER(4)                          :: NPRO
!**    INTEGER(4)                          :: NKDIV(3)
!**    INTEGER(4)                          :: ISHIFT(3)
!**    REAL(8)                             :: RNTOT
!**    REAL(8)                             :: NEL
!**    LOGICAL(4)                          :: TINV
!**    INTEGER(4)                          :: LNXX
!**    REAL(8)                             :: RBAS(3,3)
!**    REAL(8)   ,ALLOCATABLE              :: R(:,:)
!**    INTEGER(4),ALLOCATABLE              :: LNX(:)
!**    INTEGER(4),ALLOCATABLE              :: LOX(:,:)
!**    INTEGER(4),ALLOCATABLE              :: ISPECIES(:)
!**    CHARACTER(16),ALLOCATABLE           :: ATOMID(:)
!**    INTEGER(4),ALLOCATABLE              :: IZ(:)
!**    REAL(8)   ,ALLOCATABLE              :: RAD(:)
!**    REAL(8)   ,ALLOCATABLE              :: PHIOFR(:,:)
!**    REAL(8)   ,ALLOCATABLE              :: DPHIDR(:,:)
!**    REAL(8)   ,ALLOCATABLE              :: OV(:,:,:)
!**    REAL(8)   ,ALLOCATABLE              :: XK(:,:)
!**    REAL(8)   ,ALLOCATABLE              :: WKPT(:)
!**    TYPE(STATE_TYPE),ALLOCATABLE,TARGET :: STATEARR(:,:)
!**    TYPE(STATE_TYPE),POINTER            :: STATE
!**    INTEGER(4)                          :: SPACEGROUP
!**    LOGICAL(4)                          :: TSHIFT
!**  END MODULE PDOS_MODULE
!*******************************************************************************
!*******************************************************************************
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SPINDIR_MODULE
REAL(8)      ,ALLOCATABLE :: SPINDIR(:,:)
SAVE
END MODULE SPINDIR_MODULE
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE DOS_WGHT_MODULE
USE BRILLOUIN_MODULE, ONLY: EWGHT_TYPE
REAL(8),ALLOCATABLE          :: WGHT(:,:)
TYPE(EWGHT_TYPE),ALLOCATABLE :: EWGHT(:,:)
REAL(8)                      :: EF
INTEGER(4)                   :: SPACEGROUP
SAVE
END MODULE DOS_WGHT_MODULE
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE READCNTL_MODULE
USE LINKEDLIST_MODULE,ONLY : LL_TYPE
TYPE(LL_TYPE)   :: LL_CNTL
SAVE
END MODULE READCNTL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM PDOS
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **  - PDOS$READ                                                         **
!     **  - READCNTL$ORBITAL                                                  **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE,ONLY : LINKEDLIST$REPORT_UNUSED &
     &                            ,LINKEDLIST$SELECT
      USE READCNTL_MODULE  ,ONLY : LL_CNTL
      USE SPINDIR_MODULE   ,ONLY : SPINDIR !(IS ONLY ALLOCATED)
      IMPLICIT NONE
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NAT
      INTEGER(4)                :: NB
      INTEGER(4)                :: NKPT
      INTEGER(4)                :: NSPIN
      INTEGER(4)                :: NDIM !=2 FOR SPINOR WF; OTHERWISE =1
      INTEGER(4)                :: LENG
      INTEGER(4)                :: NSET
      REAL(8)      ,ALLOCATABLE :: RPOS(:,:)
      REAL(8)                   :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)      ,ALLOCATABLE :: EIG(:,:)
      REAL(8)      ,ALLOCATABLE :: OCC(:,:)
      REAL(8)      ,ALLOCATABLE :: SET(:,:,:,:)
      CHARACTER(16),ALLOCATABLE :: ATOMID(:)
      REAL(8)                   :: EMIN
      REAL(8)                   :: EMAX
      INTEGER(4)                :: NE
      REAL(8)                   :: EBROAD
      REAL(8)                   :: ELSCALE
      INTEGER(4)                :: NFILIN
      INTEGER(4)                :: NPRO
      INTEGER(4)                :: NBB             ! #(SPIN STATES PER K-POINTS)
      INTEGER(4)                :: IKPT,ISPIN,IB
      INTEGER(4)   ,ALLOCATABLE :: NBARR(:,:)
      CHARACTER(6)              :: FLAG
      CHARACTER(32)             :: MODE
      CHARACTER(256)            :: PREFIX !OPTIONAL PREFIX FOR DOS AND NOS FILES
      CHARACTER(128)            :: FORMAT
!     **************************************************************************
      CALL TRACE$PUSH('MAIN')
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!     ==========================================================================
!     ==  RESOLVE ARGUMENTLIST AND INITIALIZE FILE HANDLER                    ==
!     ==========================================================================
      CALL INITIALIZEFILEHANDLER
!
!     ==========================================================================
!     ==  ANALYZE CONTROL FILE                                                ==
!     ==========================================================================
      CALL READCNTL
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==========================================================================
!     ==  WRITE HEADER                                                        ==
!     ==========================================================================
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(80("*"),T15,A50)') &
     &            '              PDOS ANALYSIS TOOL                '
      WRITE(NFILO,FMT='(80("*"),T15,A50)') &
     &             '    FOR THE PROJECTOR-AUGMENTED WAVE METHOD   ' 
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(T10 &
     &           ,"P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY")')
      WRITE(NFILO,FMT='(T10 &
     &            ,"DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3")')
      WRITE(NFILO,*)
!
!     ==========================================================================
!     ==  READ PDOSFILE                                                       ==
!     == ( DONE BEFORE READING FROM DCNTL TO SUGGEST RANGE FOR ENERGY GRID)   ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PDOS',NFILIN)
      REWIND(NFILIN)
      CALL PDOS$READ(NFILIN)
      CALL PDOS$GETI4('NAT',NAT)
      CALL PDOS$GETI4('NKPT',NKPT)
      CALL PDOS$GETI4('NSPIN',NSPIN)
      CALL PDOS$GETI4('NDIM',NDIM)
      CALL PDOS$GETI4('NPRO',NPRO)
      ALLOCATE(NBARR(NKPT,NSPIN))
      CALL PDOS$GETI4A('NB',NKPT*NSPIN,NBARR)
      NB=MAXVAL(NBARR)
      DEALLOCATE(NBARR)
      LENG=NPRO
      CALL PDOS$GETR8A('RBAS',3*3,RBAS)
      ALLOCATE(RPOS(3,NAT))
      CALL PDOS$GETR8A('R',3*NAT,RPOS)
      ALLOCATE(ATOMID(NAT))
      CALL PDOS$GETCHA('ATOMID',NAT,ATOMID)
!
!     ==========================================================================
!     == CHECK IF PDOS FILE CONTAINS DATA FOR THE TETRAHEDRON METHOD ===========
!     ==========================================================================
      CALL PDOS$GETCH('FLAG',FLAG)
      CALL REPORT$CHVAL(NFILO,"FLAG OF PDOS FILE=",FLAG)
      IF(MODE.EQ.'TETRA'.AND.FLAG.NE.'181213')THEN 
        CALL ERROR$MSG('THE PDOS-FILE IS TOO OLD FOR MODE=TETRA')
        CALL ERROR$CHVAL('FLAG ',FLAG)
        CALL ERROR$STOP('PDOS MAIN')
      ENDIF
                            CALL TRACE$PASS('AFTER READPDOS')
!
!     ==========================================================================
!     ==  TRANSFORM STATES IN PDOS OBJECT ONTO AN ORTHORMAL BASISSET          ==
!     ==  AFTER THAT THE VARIABLE OV EQUALS THE UNIT MATRIX.                  ==
!     ==  CHANGES STATE%VEC AND OV IN PDOS_MODULE.
!     ==========================================================================
      CALL ORTHONORMALIZESTATES()
!
!     ==========================================================================
!     ==  CALCULATE ANGULAR MOMENTUM WEIGHTS AND SPINS                        ==
!     ==  AND WRITE RESULT TO THE PROTOCOL FILE.                              ==
!     ==  THE VARIABLE SPIN DIR GIVES THE LOCAL SPIN AXIS AND IS KEPT FOR LATER=
!     ==========================================================================
      ALLOCATE(SPINDIR(3,NAT))
      CALL REPORT(NFILO)
                            CALL TRACE$PASS('AFTER REPORT')
!
!     ==========================================================================
!     ==  READ PREDEFINED ORBITALS  (DATA -> NEWORBITAL_MODULE)               ==
!     ==========================================================================
      CALL READCNTL$ORBITAL(LENG,NAT,RBAS,RPOS,ATOMID)
      CALL NEWORBITAL$REPORT(NFILO)
                            CALL TRACE$PASS('AFTER READCNTL$ORBITAL')
!
!     ==========================================================================
!     ==  SELECT MATRIXELEMENTS   (->NEWSET_MODULE)                           ==
!     ==========================================================================
      CALL READCNTL$SETNUMBER(NSET)
      CALL READCNTL$SETS_NEW(NKPT,NSET,NAT,RBAS,ATOMID,RPOS)
      CALL NEWSET$REPORT(NFILO)
                            CALL TRACE$PASS('AFTER READCNTL$SETS_NEW')
!
!     ==========================================================================
!     ==  DETERMINE SELECTED MATRIX ELEMENTS FOR EACH STATE                   ==
!     ==  DATA WILL BE ENCODED BY (EIG,OCC,SET)
!     ==========================================================================
!     == NBB=#(SPIN-STATES PER K-POINT) 
      NBB=NB                         
      IF(NDIM.EQ.1)NBB=2*NB          ! COLLINEAR AND NON-SPIN-POLARIZED
      ALLOCATE(SET(NBB,NKPT,2,NSET)) ! SET WORKS ALWAYS WITH TWO SPINS
      ALLOCATE(EIG(NBB,NKPT))
      ALLOCATE(OCC(NBB,NKPT))
!     __ FILL IN EIGENVALUES AND OCCUPATIONS____________________________________
      CALL SET$ENOCC(NBB,NKPT,EIG,OCC)
      CALL NEWSET$PROCESS(NBB,NKPT,NSET,SET)
                            CALL TRACE$PASS('AFTER NEWSET$PROCESS')
!
!     ==========================================================================
!     ==  READ GENERAL INFORMATION FROM CONTROL FILE
!     ==========================================================================
      CALL READCNTL$GENERIC(MODE,PREFIX)
!     == DEFAULT VALUES FOR RANGE OF ENERGY GRID ===============================
      CALL READCNTL$GRID(EMIN,EMAX,NE,EBROAD)
      CALL READCNTL$REPORT1(MODE,PREFIX,EMIN,EMAX,NE,EBROAD)
!
!     ==========================================================================
!     ==  MAKE PLOTS                                                          ==
!     ==========================================================================
!     ==  CALCULAT WEIGHTS FOR DOS USING THE TETRAHEDRON METHOD               ==
      IF(MODE.EQ.'TETRA')THEN
        ELSCALE=1.D0
        IF(NSPIN.EQ.1.AND.NDIM.EQ.1) ELSCALE=2.D0
        CALL GENERATE_TETRA_WGHT(NFILO,NBB,NKPT,EMAX,EMIN,NE,RBAS,EIG,ELSCALE)
      ENDIF
!     == WRITE FILES ===========================================================
      CALL READCNTL$OUTPUT(EMIN,EMAX,NE,EBROAD,PREFIX &
     &                    ,NBB,NKPT,EIG,OCC,NSET,SET,MODE)
                            CALL TRACE$PASS('AFTER READCNTL$OUTPUT')
!
!     ==========================================================================
!     ==  CLOSING                                                             ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$REPORT_UNUSED(LL_CNTL,NFILO)
!
      CALL FILEHANDLER$REPORT(NFILO,'USED')
      WRITE(NFILO,FMT='(80("="))')
      WRITE(NFILO,FMT='(80("="),T20,"  PAW_DOS TOOL FINISHED  ")')
      WRITE(NFILO,FMT='(80("="))')
                            CALL TRACE$PASS('AFTER CLOSING')
!
!     ==========================================================================
!     ==  CLOSE FILES                                                         ==
!     ==========================================================================
      CALL FILEHANDLER$CLOSEALL
                            CALL TRACE$PASS('AFTER FILEHANDLER$CLOSEALL')
      CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      STOP
      END PROGRAM PDOS
!      
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE REPORT(NFILO) 
!     **************************************************************************
!     **  WRITES PROJECTED CHARGES AND SPINS FOR EACH ATOM TO                 **
!     **  DPROT FILE AND CALCULATES THE SPIN DIRECTIONS                       **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NSPIN &
     &                          ,NKPT &
     &                          ,NDIM &
     &                          ,STATEARR,STATE &
     &                          ,NAT &
     &                          ,ISPECIES &
     &                          ,LNX,LOX &
     &                          ,ATOMID &
     &                          ,R &
     &                          ,OV
      USE SPINDIR_MODULE, ONLY : SPINDIR
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFILO
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: ISPIN,IKPT,IB,IDIM
      INTEGER(4)              :: IPRO0,IPRO1,IPRO2
      INTEGER(4)              :: IAT,ISP,IAT1,IAT2,ITEN
      INTEGER(4)              :: IDIR,L1,L2,M,LN,LN1,LN2,NATSPINANGLE
      INTEGER(4) ,ALLOCATABLE :: IATSPINANGLE(:)
      REAL(8)                 :: SUM_(3),SPIN(3,NAT),TOTALSPIN(3)
      REAL(8)    ,ALLOCATABLE :: ANGWGHT(:,:,:) ! (LOX,2,NAT)
      REAL(8)                 :: SUML,ANGLE(NAT),PI
      REAL(8)                 :: SIGMA
      REAL(8)                 :: SVAR
!     **************************************************************************
                                   CALL TRACE$PUSH('REPORT')
      ALLOCATE(ANGWGHT(MAXVAL(LOX)+1,2,NAT))
      ANGWGHT(:,:,:)=0.D0
      SPIN(:,:)=0.D0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          SIGMA=REAL(3-2*ISPIN,KIND=8)
          IPRO0=0
          DO IAT=1,NAT
            ISP=ISPECIES(IAT)
            IPRO1=IPRO0
            DO LN1=1,LNX(ISP)
              L1=LOX(LN1,ISP)
              IPRO2=IPRO0
              DO LN2=1,LNX(ISP)
                L2=LOX(LN2,ISP)
                IF(L1.NE.L2) THEN
                  IPRO2=IPRO2+2*L2+1
                  CYCLE
                END IF
                DO IB=1,STATE%NB
                  SUML=0.D0
                  SUM_(:)=0.D0
                  DO M=1,2*L1+1
                      DO IDIM=1,NDIM
                        SUML=SUML+REAL(CONJG(STATE%VEC(IDIM,IPRO1+M,IB)) &
     &                                      *STATE%VEC(IDIM,IPRO2+M,IB))
                      END DO
                      IF(NSPIN.EQ.2.AND.NDIM.EQ.1) THEN
                        SUM_(3)=SUM_(3)+SIGMA &
     &                                  *REAL(CONJG(STATE%VEC(1,IPRO1+M,IB)) &
     &                                             *STATE%VEC(1,IPRO2+M,IB))
                      ELSE IF(NDIM.EQ.2) THEN
                        SUM_(1)=SUM_(1)+2.D0 &
     &                                 *REAL(CONJG(STATE%VEC(1,IPRO1+M,IB)) &
     &                                            *STATE%VEC(2,IPRO2+M,IB)) 
                        SUM_(2)=SUM_(2)+2.D0 &
     &                                 *AIMAG(CONJG(STATE%VEC(1,IPRO1+M,IB)) &
     &                                             *STATE%VEC(2,IPRO2+M,IB))
                        SUM_(3)=SUM_(3)+REAL(CONJG(STATE%VEC(1,IPRO1+M,IB)) &
     &                                            *STATE%VEC(1,IPRO2+M,IB) &
     &                                      -CONJG(STATE%VEC(2,IPRO1+M,IB)) &
     &                                             *STATE%VEC(2,IPRO2+M,IB))
                      END IF
                  ENDDO
                  DO IDIR=1,3
                    SPIN(IDIR,IAT)=SPIN(IDIR,IAT) &
     &                            +SUM_(IDIR)*STATE%OCC(IB)*OV(LN1,LN2,ISP)
                  END DO
                  ANGWGHT(L1+1,ISPIN,IAT)=ANGWGHT(L1+1,ISPIN,IAT) &
     &                                   +SUML*STATE%OCC(IB)*OV(LN1,LN2,ISP)
                ENDDO  !END OF LOOP OVER BANDS
                IPRO2=IPRO2+2*L2+1
              ENDDO            
              IPRO1=IPRO1+2*L1+1
            ENDDO
            DO LN=1,LNX(ISP)
              IPRO0=IPRO0+2*LOX(LN,ISP)+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CALCULATE SPIN DIRECTION FOR EACH ATOM                               ==
!     ==========================================================================
      IF(NDIM.EQ.2.OR.NSPIN.EQ.2) THEN
        DO IAT=1,NAT
          SUML=SPIN(1,IAT)**2+SPIN(2,IAT)**2+SPIN(3,IAT)**2
          DO IDIR=1,3
            SPINDIR(IDIR,IAT)=SPIN(IDIR,IAT)/SQRT(SUML)
          END DO
        END DO
      END IF
!
!     ==========================================================================
!     == REPORT IN CHARGE DISTRIBUTION                                        ==
!     ==========================================================================
      CALL REPORT$TITLE(NFILO,'PROJECTED CHARGE ANALYSIS')
      WRITE(NFILO,FMT='(T11,"|",T28,"CHARGE[-E] PROJECTED ON")')
      WRITE(NFILO &
     &     ,FMT='(T3,"ATOM",T11,"|",T18,"ALL",T28,"S",T38,"P",T48,"D",T58,"F")')
      WRITE(NFILO,FMT='(65("-"))')
      DO IAT=1,NAT
        L1=1+MAXVAL(LOX(:,ISPECIES(IAT)))
        WRITE(NFILO,FMT='(A10,"|",10F10.3)')  &
     &       ATOMID(IAT),SUM(ANGWGHT(:L1,:,IAT)) &
     &                  ,(SUM(ANGWGHT(IDIR,:,IAT)),IDIR=1,L1)
      END DO
!
!     ==========================================================================
!     == SPIN REPORT FOR COLLINEAR CALCULATION                                ==
!     ==========================================================================
      IF(NSPIN.EQ.2) THEN
        CALL REPORT$TITLE(NFILO,'PROJECTED SPIN ANALYSIS')
        WRITE(NFILO,FMT='(T11,"|",T28,"SPIN[HBAR/2] PROJECTED ON")')
        WRITE(NFILO &
     &     ,FMT='(T3,"ATOM",T11,"|",T18,"ALL",T28,"S",T38,"P",T48,"D",T58,"F")')
        WRITE(NFILO,FMT='(65("-"))')
        DO IAT=1,NAT
          L1=1+MAXVAL(LOX(:,ISPECIES(IAT)))
          WRITE(NFILO,FMT='(A10,"|",10F10.3)') &
     &         ATOMID(IAT),SUM(ANGWGHT(:,1,IAT))-SUM(ANGWGHT(:,2,IAT)) &
     &                    ,(ANGWGHT(IDIR,1,IAT)-ANGWGHT(IDIR,2,IAT),IDIR=1,L1)
        END DO
      END IF
!
!     ==========================================================================
!     == SPIN REPORT FOR NON-COLLINEAR CALCULATION                            ==
!     ==========================================================================
      IF(NDIM.EQ.2) THEN
        TOTALSPIN(:)=0.D0
        WRITE(NFILO,FMT='(82("="))')
        WRITE(NFILO,FMT='(82("="),T30," SPIN ANALYSIS  ")')
        WRITE(NFILO,FMT='(82("="))')
        WRITE(NFILO,FMT='(A,T40,"X",T50,"Y",T60,"Z",T69,"TOTAL")') &
     &                   'TOTAL SPIN PROJECTED ON'
!
!       ==  SPIN DIRECTIONS ====================================================
        DO IAT=1,NAT
          WRITE(NFILO,FMT='("SPIN[HBAR/2] ON ATOM",T23,A10,":",4F10.3)')  &
     &          ATOMID(IAT),SPIN(1,IAT),SPIN(2,IAT),SPIN(3,IAT) &
     &                     ,SQRT(SPIN(1,IAT)**2+SPIN(2,IAT)**2+SPIN(3,IAT)**2)
          DO IDIR=1,3
            TOTALSPIN(IDIR)=TOTALSPIN(IDIR)+SPIN(IDIR,IAT)
          END DO
        END DO
        WRITE(NFILO,FMT='(" TOTAL PROJECTED SPIN:      ",4F10.3)') &
     &         TOTALSPIN(1),TOTALSPIN(2),TOTALSPIN(3),SQRT(TOTALSPIN(1)**2+&
     &         TOTALSPIN(2)**2+TOTALSPIN(3)**2)
!
!       ==  PRINT ANGLES BETWEEN THE SPINS =====================================
        CALL CONSTANTS('PI',PI)

!       CHOSE ATOMS WITH SPIN GREATER THAN 0.1 HBAR/2
        NATSPINANGLE=0
        ALLOCATE(IATSPINANGLE(NAT))
        DO IAT1=1,NAT
          IF(SQRT(SUM(SPIN(:,IAT1)**2)).LE.0.1D0) CYCLE
          NATSPINANGLE=NATSPINANGLE+1
          IATSPINANGLE(NATSPINANGLE)=IAT1
        END DO
        IF(NAT.GE.2) THEN
          WRITE(NFILO,FMT='(82("="),T20,A)') &
    &              " ANGLES [DEG] BETWEEN THE SPINS > 0.1 ON THE ATOMS "
          WRITE(NFILO,FMT='(T9,10(1X,A7))') &
     &                       (ATOMID(IATSPINANGLE(IAT2)),IAT2=NATSPINANGLE,2,-1)
          DO IAT1=1,NATSPINANGLE !SENKRECHT
            DO IAT2=NATSPINANGLE,IAT1+1,-1   !WAAGRECHT
              ANGLE(IAT2)=180.D0/PI*ACOS(SUM(SPINDIR(:,IATSPINANGLE(IAT1)) &
     &                                      *SPINDIR(:,IATSPINANGLE(IAT2))))
            END DO
            ITEN=NATSPINANGLE
            DO WHILE (IAT1+1.LE.ITEN)
              WRITE(NFILO,FMT='(A6,10F8.1)')ATOMID(IATSPINANGLE(IAT1)) &
     &                         ,(ANGLE(IAT2),IAT2=ITEN,MAX(ITEN-9,IAT1+1),-1)
              ITEN=ITEN-10
            ENDDO
          END DO
        END IF
        DEALLOCATE(IATSPINANGLE)
!
!       ========================================================================
!       ==  THIS BLOCK IS INTENDED AS INPUTFILE FOR MOLDEN                    ==
!       ==  MOLDEN WILL PLOT THE SPIN DISTRIBUTION                            ==
!       ========================================================================
        CALL FILEHANDLER$UNIT('MOL',NFIL)
        WRITE(NFIL,*)'[MOLDEN FORMAT]'
        WRITE(NFIL,*)'[GEOMETRIES] XYZ'
        WRITE(NFIL,*)'    ',NAT
        WRITE(NFIL,*)' '
        DO IAT=1,NAT
          WRITE(NFIL,FMT='(A2,F10.5,F10.5,F10.5)')ATOMID(IAT),R(:,IAT)
        END DO
        WRITE(NFIL,*)' '
        WRITE(NFIL,*)'[FREQ]'
        WRITE(NFIL,*)'4.'
        WRITE(NFIL,*)'[FR-COORD]'
        DO IAT=1,NAT
          WRITE(NFIL,FMT='(A2,F10.5,F10.5,F10.5)')ATOMID(IAT),R(:,IAT)
        END DO
        WRITE(NFIL,*)'[FR-NORM-COORD]'
        WRITE(NFIL,*)'VIBRATION 1'
        DO IAT=1,NAT
          WRITE(NFIL,FMT='(F10.5,F10.5,F10.5)')SPIN(:,IAT)
        END DO
      END IF
                           CALL TRACE$POP
      RETURN
     END SUBROUTINE REPORT
!      
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
!     **************************************************************************
!     **************************************************************************
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: PDOSINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
!     **************************************************************************
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE PDOS TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,PDOSINNAME)
      ISVAR=INDEX(PDOSINNAME,-'.DCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=PDOSINNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('DCNTL',.FALSE.,PDOSINNAME)
      RETURN
      END SUBROUTINE INITIALIZEFILEHANDLER
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STANDARDFILES
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: ID
!     **************************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!  
!     ==========================================================================
!     == SET STANDARD FILENAMES                                               ==
!     ==========================================================================
!
!     ==  ERROR FILE ===========================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOL FILE ========================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =====================================================
      ID=+'DCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  STRUCTURE FILE   =====================================================
      ID=+'PDOS'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PDOS')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  DENSITY OF STATES FILE PRODUCES AS OUTPUT ============================
!     ==  WILL BE ATTACHED TO DIFFERENT FILES DURING EXECUTION =================
      ID=+'PDOSOUT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PDOSOUT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  NUMBER OF STATES FILE PRODUCES AS OUTPUT =============================
!     ==  WILL BE ATTACHED TO DIFFERENT FILES DURING EXECUTION =================
      ID=+'PNOSOUT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PNOSOUT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  SPIN GRAPHICS FILE   =================================================
      ID=+'MOL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.MOL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ORTHONORMALIZESTATES()
!     **************************************************************************
!     ** TRANSFORMS THE LOCAL BASIS SET SO THAT IT IS ORTHONORMAL.            **
!     ** A CHOLESKY DECOMPOSITION OF THE LOCAL ORBITALS TRANSFORMS THE        **
!     ** OVERLAP INTO A UNIT MATRIX                                           **
!     **                                                                      **
!     ** CHANGES STATE%VEC AND OV ON PDOS_MODULE                              **
!     **                                                                      **
!     ** THE FIRST ORBITAL FOR EACH ANGULAR MOMENTUM IS ONLY NORMALIZED       **
!     ** THE SECOND ORBITAL FOR EACH ANGULAR MOMENTUM IS ORTHONORMALIZED      **
!     **     TO THE FIRST ORBITAL AND THEN NORMALIZED                         **
!     **************************************************************************
      USE PDOS_MODULE, ONLY : NSP &
     &                       ,LNX &      !(NSP)
     &                       ,LOX &      !(LNXX,NSP)
     &                       ,NAT &
     &                       ,NDIM &
     &                       ,ISPECIES & !(NAT)
     &                       ,OV  &      !(LNXX,LNXX,NSP)
     &                       ,NKPT  & 
     &                       ,NSPIN  & 
     &                       ,STATEARR & !(IKPT,ISPIN)
     &                       ,STATE 
      IMPLICIT NONE
      INTEGER(4)             :: LX
      INTEGER(4)             :: LMNX(NSP)
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: TRANSFORM(:,:,:)
      INTEGER(4)             :: ISP,LN1,LN2,L1,L2,LMN1,LMN2,IM,ISVAR,IDIM,LN
      INTEGER(4)             :: IPRO,IKPT,ISPIN,IAT
!     **************************************************************************
      LX=MAXVAL(LOX)
      DO ISP=1,NSP
        LMNX(ISP)=SUM(2*LOX(:LNX(ISP),ISP)+1)
      ENDDO
      ISVAR=MAXVAL(LMNX)
      ALLOCATE(TRANSFORM(ISVAR,ISVAR,NSP))
      TRANSFORM(:,:,:)=(0.D0,0.D0)
!
!     ==========================================================================
!     == DETERMINE TRANSFORMATION MATRIX USING CHOLESKY DECOMPOSITION         ==
!     ==========================================================================
      DO ISP=1,NSP
!!$        WRITE(*,FMT='(80("="),T10," OV FOR ISP=",I2,"  ")')ISP
!!$        DO LN=1,LNX(ISP)
!!$          WRITE(*,FMT='("OV=",10F10.5)')OV(LN,:LNX(ISP),ISP)
!!$        ENDDO 
        ALLOCATE(AMAT(LNX(ISP),LNX(ISP)))
        ALLOCATE(BMAT(LNX(ISP),LNX(ISP)))
        AMAT=OV(:LNX(ISP),:LNX(ISP),ISP)
!       ==  AMAT=MATMUL(BMAT,TRANSPOSE(BMAT)) WITH BMAT(I,J)=0 FOR J>I =========
        CALL CHOLESKY(LNX(ISP),AMAT,BMAT)
        BMAT=TRANSPOSE(BMAT)
        DEALLOCATE(AMAT)
!
!!$        WRITE(*,FMT='(80("="),T10," G FOR ISP=",I2,"  ")')ISP
!!$        DO LN1=1,LNX(ISP)
!!$          WRITE(*,FMT='(10F10.5)')BMAT(:,LN1)
!!$        ENDDO 
!
!       == UPFOLD BMAT AND PLACE INTO TRANSFORM ================================
        LMN1=0
        DO LN1=1,LNX(ISP)
          L1=LOX(LN1,ISP)
          LMN2=0
          DO LN2=1,LNX(ISP)
            L2=LOX(LN2,ISP)
            IF(L2.EQ.L1) THEN
              DO IM=1,2*L1+1
                TRANSFORM(LMN1+IM,LMN2+IM,ISP)=CMPLX(BMAT(LN1,LN2),KIND=8)
              ENDDO
            ENDIF
            LMN2=LMN2+2*L2+1
          ENDDO
          LMN1=LMN1+2*L1+1
        ENDDO
        DEALLOCATE(BMAT)
      ENDDO !ISP
!
!     ==========================================================================
!     == TRANSSFORM STATE ARRAY                                               ==
!     ==========================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          IPRO=0
          DO IAT=1,NAT
            ISP=ISPECIES(IAT)
            DO IDIM=1,NDIM
              STATE%VEC(IDIM,IPRO+1:IPRO+LMNX(ISP),:) &
   &                            =MATMUL(TRANSFORM(:LMNX(ISP),:LMNX(ISP),ISP) &
   &                                   ,STATE%VEC(IDIM,IPRO+1:IPRO+LMNX(ISP),:))
            ENDDO
            IPRO=IPRO+LMNX(ISP)
          ENDDO
        ENDDO
      ENDDO        
      DEALLOCATE(TRANSFORM)
!
!     ==========================================================================
!     == SET OV TO UNITY                                                      ==
!     ==========================================================================
      OV(:,:,:)=0.D0
      DO ISP=1,NSP
        DO LN=1,LNX(ISP)
          OV(LN,LN,ISP)=1.D0
        ENDDO
      ENDDO
!
!!$      DO ISP=1,NSP
!!$        WRITE(*,FMT='(80("="),T10," OV(SET TO 1) FOR ISP=",I2,"  ")')ISP
!!$        DO LN=1,LNX(ISP)
!!$          WRITE(*,FMT='(10F10.5)')OV(LN,:LNX(ISP),ISP)
!!$        ENDDO 
!!$      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CHOLESKY(N,A,G)
!     **************************************************************************
!     ** CHOLESKY DECOMPOSITION OF THE POSITIVE-DEFINITE, SYMMETRIC MATRIX A  **
!     **   A=MATMUL(G,TRANSPOSE(G))                                           **
!     **   G HAS VALUES ONLY ON LOWER TRIANGULAR, I.E.  G(I,J)=0 FOR J>I      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: A(N,N)
      REAL(8)   ,INTENT(OUT):: G(N,N)
      INTEGER(4)           :: I,J,K
      REAL(8)               :: S
!     **************************************************************************
      G=A
      DO I=1,N
        DO J=1,I
          S=G(I,J)
          DO K=1,J-1
            S=S-G(I,K)*G(J,K)
          ENDDO
          IF(I.GT.J) THEN
            G(I,J)=S/G(J,J)
          ELSE
            IF(S.LE.0.D0) THEN
               STOP 'MATRIX NOT POSITIVE DEFINITE'
            END IF
            G(I,I)=SQRT(S)
          END IF
        ENDDO
      ENDDO
      DO I=1,N
        DO J=I+1,N
          G(I,J)=0.D0
        ENDDO
      ENDDO
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE NEWORBITAL_MODULE
!*******************************************************************************
!** CONTAINER FOR ORBITALS USED TO PROJECT THE DENSITY OF STATES              **
!**                                                                           **
!** EACH ORBITAL CONSISTS OF ENTRIES THAT HOLD THE CONTRIBUTION OF A GIVEN    **
!** SITE IDENTIFIED BY AN ATOM INDEX IAT AND A TRANSLATION IT(3).             **
!**                                                                           **
!** 1) CREATE ORBITAL WITH NEWORBITAL(NAME)                                   **
!** 2) SELECT AND UNSELECT ORBITAL WITH SELECT(NAME) OR ISELECT(IORB).        **
!**    UNSELECT WITH SELECT('NONE') OR ISELECT(0)                             **
!**                                                                           **
!*******************************************************************************
IMPLICIT NONE
TYPE ORBITALENTRY_TYPE
  ! EACH ENTRY DESCRIBES THE ORBITALS CENTERED ON A GIVEN SITE
  INTEGER(4)             :: IAT
  INTEGER(4)             :: IT(3)
  INTEGER(4)             :: LMNX
  COMPLEX(8),ALLOCATABLE :: ORB(:)  !(LMNX)
  INTEGER(4),ALLOCATABLE :: IPRO(:) !(LMNX)
END TYPE ORBITALENTRY_TYPE
!
TYPE ORBITAL_TYPE
  CHARACTER(32)           :: NAME
  INTEGER(4)              :: NENTRY=0
  TYPE(ORBITALENTRY_TYPE) :: ENTRY(1000)
END TYPE ORBITAL_TYPE
!
INTEGER(4)        ,PARAMETER   :: NORBX=1000
INTEGER(4)                     :: NORB=0
INTEGER(4)                     :: IORB=0
TYPE(ORBITAL_TYPE)             :: NEWORBITAL(NORBX) !(NORBX)
END MODULE NEWORBITAL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$GETI4(ID,IVAL)
!     **************************************************************************
!     **  COLLECT INFORMATION (INTEGER(4))                                    **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : NORB,IORB,NEWORBITAL
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: IVAL
!     **************************************************************************
      IF(ID.EQ.'NORB') THEN
        IVAL=NORB
!
      ELSE IF(ID.EQ.'IORB') THEN
        IVAL=IORB
!
      ELSE IF(ID.EQ.'NENTRY') THEN
        IF(IORB.EQ.0) THEN
          CALL ERROR$MSG('NO ORBITAL SELECTED')
          CALL ERROR$MSG('ITEM "NENTRY" IS NOT AVAILABLE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('NEWORBITAL$GETI4')
        END IF
        IVAL=NEWORBITAL(IORB)%NENTRY

      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('NEWORBITAL$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$IORB(NAME,IORB1)
!     **************************************************************************
!     **  IDENTIFY THE INDEX OF THE ORBITAL WITH SPECIFIED NAME               **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : NORB,NEWORBITAL
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(OUT):: IORB1
      INTEGER(4)              :: I
!     **************************************************************************
      DO I=1,NORB
        IF(NAME.EQ.NEWORBITAL(I)%NAME) THEN
          IORB1=I
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ORBITAL INDEX NOT FOUND')
      CALL ERROR$CHVAL('NAME',NAME)
      CALL ERROR$STOP('NEWORBITAL$IORB')
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$REPORT(NFIL)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : NORB,NEWORBITAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NENTRY
      INTEGER(4)            :: I,J
!     **************************************************************************
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='(80("="),T20,"  NEWORBITAL REPORT  ")')
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='(40("."),":",T1,A,T42,I3)')'NUMBER OF ORBITALS',NORB
      DO I=1,NORB
        WRITE(NFIL,FMT='("ORBITAL ",A)')NEWORBITAL(I)%NAME
        NENTRY=NEWORBITAL(I)%NENTRY
        DO J=1,NENTRY
          WRITE(NFIL,FMT='("IAT=",I4," IT=",3I2," ORB=",99("(",2F10.5,")"))') &
    &                NEWORBITAL(I)%ENTRY(J)%IAT,NEWORBITAL(I)%ENTRY(J)%IT &
    &               ,NEWORBITAL(I)%ENTRY(J)%ORB
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$NEWORBITAL(NAME)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **  THE NEW ORBITAL HAS A NAME BUT NO ENTRIES YET                       **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : NORBX,NORB,NEWORBITAL
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)              :: I
!     **************************************************************************
!
!     ==========================================================================
!     == CHECK IF NAME HAS BEEN USED BEFORE                                   ==
!     ==========================================================================
      DO I=1,NORB
        IF(NAME.EQ.NEWORBITAL(I)%NAME) THEN
          CALL ERROR$MSG('ORBITAL WITH THE NAME ALREADY EXISTS')
          CALL ERROR$CHVAL('NAME',NAME)
          CALL ERROR$STOP('NEWORBITAL$NEWORBITAL')
        END IF
      ENDDO
!
!     ==========================================================================
!     == CHECK IF ARRAY SIZE IS SUFFICIENT                                    ==
!     ==========================================================================
      IF(NORB+1.GT.NORBX) THEN
        CALL ERROR$MSG('MAX. NR. OF ORBITALS EXCEEDED')
        CALL ERROR$MSG('INCREASE PARAMETER "NORBX" IN NEWORBITAL_MODULE')
        CALL ERROR$I4VAL('NORBX',NORBX)
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('NEWORBITAL$NEWORBITAL')
      END IF
!
!     ==========================================================================
!     == CHECK IF ARRAY SIZE IS SUFFICIENT                                    ==
!     ==========================================================================
      NORB=NORB+1
      NEWORBITAL(NORB)%NAME=NAME
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$SELECT(NAME)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : IORB,NORB,NEWORBITAL
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)              :: I
!     **************************************************************************
      IF(NAME.EQ.'NONE') THEN
        IORB=0
        RETURN
      END IF
!
!     ==========================================================================
!     == CHECK IF NAME HAS BEEN USED BEFORE                                   ==
!     ==========================================================================
      IF(IORB.NE.0) THEN
         CALL ERROR$MSG('ATTEMPT TO SELECT ORBITAL WHILE ANOTHER IS ACTIVE')
         CALL ERROR$MSG('UNSELECT FIRST')
         CALL ERROR$CHVAL('NAME',NAME)
         CALL ERROR$STOP('NEWORBITAL$SELECT')
      END IF
     
      DO I=1,NORB
        IF(NAME.EQ.NEWORBITAL(I)%NAME) THEN
          IORB=I 
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ORBITAL NOT FOUND')
      CALL ERROR$CHVAL('NAME',NAME)
      CALL ERROR$STOP('NEWORBITAL$SELECT')
      STOP
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$ISELECT(IORB1)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : IORB,NORB,NEWORBITAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IORB1
!     **************************************************************************
      IF(IORB1.EQ.0) THEN
        IORB=0
        RETURN
      END IF
!
!     ==========================================================================
!     == CHECK IF ORBITAL POINTER IS IN ALLOWED RANGE                         ==
!     ==========================================================================
      IF(IORB1.LT.1.OR.IORB1.GT.NORB) THEN
        CALL ERROR$MSG('SPECIFIED ORBITAL INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('IORB1',IORB1)
        CALL ERROR$I4VAL('NORB',NORB)
        CALL ERROR$STOP('NEWORBITAL$ISELECT')
      END IF
!
!     ==========================================================================
!     == CHECK IF NAME HAS BEEN USED BEFORE                                   ==
!     ==========================================================================
      IF(IORB.NE.0) THEN
        CALL ERROR$MSG('ATTEMPT TO SELECT ORBITAL WHILE ANOTHER IS ACTIVE')
        CALL ERROR$MSG('UNSELECT FIRST')
        CALL ERROR$CHVAL('NAME ACTIVE',NEWORBITAL(IORB)%NAME)
        CALL ERROR$CHVAL('NAME REQUESTED',NEWORBITAL(IORB1)%NAME)
        CALL ERROR$STOP('NEWORBITAL$ISELECT')
      END IF
!
!     ==========================================================================
!     == SET ORBITAL POINTER                                                  ==
!     ==========================================================================
      IORB=IORB1
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$ADDENTRY(IAT,IT,LMNX,ORB)
!     **************************************************************************
!     **  ADD AN ENTRY TO AN ORBITAL, OR, IF ONE ALREADY EXISTS FOR THAT SITE,**
!     **  ADD COEFFICIENTS TO THOSE IN THAT ENTRY                             **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : IORB,NEWORBITAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: IT(3)
      INTEGER(4),INTENT(IN) :: LMNX
      COMPLEX(8),INTENT(IN) :: ORB(LMNX)
      INTEGER(4)            :: IENTRY
      INTEGER(4)            :: NENTRY
      INTEGER(4)            :: I
!     **************************************************************************
      IF(IORB.EQ.0) THEN
        CALL ERROR$MSG('ORBITAL NOT SELECTED')
        CALL ERROR$I4VAL('IORB',IORB)
        CALL ERROR$STOP('NEWORBITAL$ADDENTRY')
      END IF
!
!     ==========================================================================
!     == FIND INDEX OF EXISTING ENTRY                                         ==
!     ==========================================================================
      NENTRY=NEWORBITAL(IORB)%NENTRY
      IENTRY=0
      DO I=1,NENTRY
        IF(IAT.NE.NEWORBITAL(IORB)%ENTRY(I)%IAT) CYCLE
        IF(IT(1).NE.NEWORBITAL(IORB)%ENTRY(I)%IT(1)) CYCLE
        IF(IT(2).NE.NEWORBITAL(IORB)%ENTRY(I)%IT(2)) CYCLE
        IF(IT(3).NE.NEWORBITAL(IORB)%ENTRY(I)%IT(3)) CYCLE
        IENTRY=I
        IF(LMNX.NE.NEWORBITAL(IORB)%ENTRY(IENTRY)%LMNX) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE OF ORBITAL ARRAY FOR SITE')
          CALL ERROR$I4VAL('ORBITAL ARRAY SIZE OFFERED',LMNX)
          CALL ERROR$I4VAL('ORBITAL ARRAY SIZE ACCEPTED' &
    &                     ,NEWORBITAL(IORB)%ENTRY(IENTRY)%LMNX)
          CALL ERROR$I4VAL('IENTRY',IENTRY)
          CALL ERROR$STOP('NEWORBITAL$ADDENTRY')
        END IF
        NEWORBITAL(IORB)%ENTRY(IENTRY)%ORB=NEWORBITAL(IORB)%ENTRY(IENTRY)%ORB &
     &                                    +ORB
        EXIT
      ENDDO
!
!     ==========================================================================
!     == CREATE NEW ENTRY IF IT DOES NOT EXIST
!     ==========================================================================
      IF(IENTRY.EQ.0) THEN
        NENTRY=NENTRY+1
        IF(NENTRY.GT.SIZE(NEWORBITAL(IORB)%ENTRY)) THEN
          CALL ERROR$MSG('NUMBER OF ENTRIES FOR ORBITAL EXCEEDED')
          CALL ERROR$CHVAL('ORBITAL NAME',NEWORBITAL(IORB)%NAME)
          CALL ERROR$I4VAL('NENTRY',NENTRY)
          CALL ERROR$I4VAL('SIZE OF ENTRY',SIZE(NEWORBITAL(IORB)%ENTRY))
          CALL ERROR$STOP('NEWORBITAL$ADDENTRY')
        END IF
        NEWORBITAL(IORB)%NENTRY=NENTRY
        NEWORBITAL(IORB)%ENTRY(NENTRY)%IAT=IAT
        NEWORBITAL(IORB)%ENTRY(NENTRY)%IT=IT
        NEWORBITAL(IORB)%ENTRY(NENTRY)%LMNX=LMNX
        ALLOCATE(NEWORBITAL(IORB)%ENTRY(NENTRY)%ORB(LMNX))
        NEWORBITAL(IORB)%ENTRY(NENTRY)%ORB=ORB
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$GETENTRY(IENTRY,IAT,IT,LMNXX,LMNX,ORB)
!     **************************************************************************
!     **  RETURN THE ORBITAL DESCRIPTORS AND COEFFICIENTS OF AN ENTRY         **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : IORB,NEWORBITAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IENTRY
      INTEGER(4),INTENT(OUT):: IAT
      INTEGER(4),INTENT(OUT):: IT(3)
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(OUT):: LMNX
      COMPLEX(8),INTENT(OUT):: ORB(LMNXX)
      INTEGER(4)            :: IND
!     **************************************************************************
      IF(IORB.EQ.0) THEN
        CALL ERROR$MSG('ORBITAL NOT SELECTED')
        CALL ERROR$STOP('NEWORBITAL$GETENTRY')
      END IF
      IF(IENTRY.GT.NEWORBITAL(IORB)%NENTRY) THEN
        CALL ERROR$MSG('ENTRY INDEX OUT OF RANGE')
        CALL ERROR$STOP('NEWORBITAL$GETENTRY')
      END IF
      IF(LMNXX.LT.NEWORBITAL(IORB)%ENTRY(IENTRY)%LMNX) THEN
        CALL ERROR$MSG('INCONSISTENT SIZE OF ORBITAL ARRAY')
        CALL ERROR$STOP('NEWORBITAL$GETENTRY')
      END IF
      IAT=NEWORBITAL(IORB)%ENTRY(IENTRY)%IAT
      IT=NEWORBITAL(IORB)%ENTRY(IENTRY)%IT
      LMNX=NEWORBITAL(IORB)%ENTRY(IENTRY)%LMNX
      IND=MIN(LMNX,LMNXX)
      ORB(:IND)=NEWORBITAL(IORB)%ENTRY(IENTRY)%ORB(:IND)
      ORB(IND+1:)=(0.D0,0.D0)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$SCALEORBITAL(FACTOR)
!     **************************************************************************
!     **  SCALE SELECTED ORBITAL (IORB) BY A FACTOR                           **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : IORB,NEWORBITAL
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN) :: FACTOR
      INTEGER(4)            :: NENTRY
      INTEGER(4)            :: I
!     **************************************************************************
      NENTRY=NEWORBITAL(IORB)%NENTRY
      DO I=1,NENTRY
        NEWORBITAL(IORB)%ENTRY(I)%ORB=NEWORBITAL(IORB)%ENTRY(I)%ORB*FACTOR
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$SHIFTORBITAL(IT)
!     **************************************************************************
!     **  ADD AN ORBITAL TO THE SELECTED ONE. THE ORBITAL CAN BE DISPLACED    **
!     **  WITH LATICE TRANSLATIONS GIVEN BY IT AND MULTIPLIED WITH A FACTOR   **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : IORB,NEWORBITAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IT(3)
      INTEGER(4)            :: NENTRY
      INTEGER(4)            :: I
!     **************************************************************************
      NENTRY=NEWORBITAL(IORB)%NENTRY
      DO I=1,NENTRY
        NEWORBITAL(IORB)%ENTRY(I)%IT=NEWORBITAL(IORB)%ENTRY(I)%IT+IT
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$ADDORBITAL(IORB2,IT2,FACTOR2)
!     **************************************************************************
!     **  ADD AN ORBITAL TO THE SELECTED ONE. THE ORBITAL CAN BE DISPLACED    **
!     **  WITH LATICE TRANSLATIONS GIVEN BY IT AND MULTIPLIED WITH A FACTOR   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IORB2
      INTEGER(4),INTENT(IN) :: IT2(3)
      COMPLEX(8),INTENT(IN) :: FACTOR2
      INTEGER(4),PARAMETER  :: LMNXX=16
      COMPLEX(8)            :: ORB(LMNXX)
      INTEGER(4)            :: NENTRY2
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: IAT
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: IORB1
      INTEGER(4)            :: I
!     **************************************************************************
      CALL NEWORBITAL$GETI4('IORB',IORB1)
      CALL NEWORBITAL$ISELECT(0)
      CALL NEWORBITAL$ISELECT(IORB2)
      CALL NEWORBITAL$GETI4('NENTRY',NENTRY2)
      DO I=1,NENTRY2
        CALL NEWORBITAL$ISELECT(0)
        CALL NEWORBITAL$ISELECT(IORB2)
        CALL NEWORBITAL$GETENTRY(I,IAT,IT,LMNXX,LMNX,ORB)
        ORB(:LMNX)=ORB(:LMNX)*FACTOR2
        CALL NEWORBITAL$ISELECT(0)
        CALL NEWORBITAL$ISELECT(IORB1)
        CALL NEWORBITAL$ADDENTRY(IAT,IT+IT2,LMNX,ORB(:LMNX))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL$ORBPRO(TINV,NBB,NKPT,NORB,ORBPRO)
!     **************************************************************************
!     ** CALCULATE THE PROJECTION OF THE KOHN-SHAM WAVE FUNCTIONS ONTO THE    **
!     ** DEFINED ORBITALS.                                                    **
!     **                                                                      **
!     ** WITH TINV=T, THE WAVE FUNCTIONS ARE TRANSFORMED THROUGH              **
!     ** TIME-INVERSION SYMMETRY.                                             **
!     **   FOR NDIM=1, PSI(S,-K)=CONJG(PSI(S,K)).                             **
!     **   FOR NDIM=2, PSI(UP,-K)=-CONJG(PSI(DN,K)) AND                       **
!     **               PSI(DN,-K)=+CONJG(PSI(UP,K))                           **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NDIM &
     &                          ,NSPIN &
     &                          ,NKPT1=>NKPT &
     &                          ,XK &          !(3,NKPT)
     &                          ,STATEARR &    !(IKPT,ISPIN)
     &                          ,STATE 
      USE NEWORBITAL_MODULE, ONLY : NORB1=>NORB &
     &                             ,NEWORBITAL
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TINV
      INTEGER(4),INTENT(IN) :: NBB    !=2*NB FOR COLLINEAR AND NON-SPIN CALC.
      INTEGER(4),INTENT(IN) :: NKPT
      INTEGER(4),INTENT(IN) :: NORB
      COMPLEX(8),INTENT(OUT):: ORBPRO(2,NBB,NKPT,NORB)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: NENTRY
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: NB
      INTEGER(4)            :: IKPT,ISPIN,IORB,IENTRY,LMN,IPRO
      INTEGER(4)            :: IAT
      INTEGER(4)            :: IT(3)
      COMPLEX(8)            :: EIK1,EIK2,EIK3,EIKT
      COMPLEX(8)            :: CFAC
!     **************************************************************************
      CALL NEWORBITAL_COMPLETEIPRO()
!
      ORBPRO=(0.D0,0.D0)
      DO IKPT=1,NKPT
!       == E^(IKT)=EXP(I* XK * GBAS^T * RBAS * IT )
!       ==        = EIK1^IT(1) * EIK2^IT(2) * EIK3^IT(3)
        EIK1=EXP(-CI*2.D0*PI*XK(1,IKPT))
        EIK2=EXP(-CI*2.D0*PI*XK(2,IKPT))
        EIK3=EXP(-CI*2.D0*PI*XK(3,IKPT))
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          NB=STATE%NB
          DO IORB=1,NORB
            NENTRY=NEWORBITAL(IORB)%NENTRY
            DO IENTRY=1,NENTRY
              IAT=NEWORBITAL(IORB)%ENTRY(IENTRY)%IAT
              IT=NEWORBITAL(IORB)%ENTRY(IENTRY)%IT
              EIKT=(EIK1**IT(1))*(EIK2**IT(2))*(EIK3**IT(3))
              LMNX=NEWORBITAL(IORB)%ENTRY(IENTRY)%LMNX
              DO LMN=1,LMNX
                IPRO=NEWORBITAL(IORB)%ENTRY(IENTRY)%IPRO(LMN)
                IF(IPRO.EQ.0) CYCLE
                IF(.NOT.TINV) THEN
                  CFAC=CONJG(NEWORBITAL(IORB)%ENTRY(IENTRY)%ORB(LMN)*EIKT)
                  IF(NDIM.EQ.2) THEN
                    ORBPRO(:,:,IKPT,IORB)=ORBPRO(:,:,IKPT,IORB) &
     &                                 +CFAC*STATE%VEC(:,IPRO,:)
                  ELSE IF(NDIM.EQ.1) THEN
                    IF(ISPIN.EQ.1) THEN ! FIRST SPIN DIRECTION
                      ORBPRO(1,:NB,IKPT,IORB)=ORBPRO(1,:NB,IKPT,IORB) &
     &                                     +CFAC*STATE%VEC(1,IPRO,:)
                    END IF
                    IF(ISPIN.EQ.2.OR.NSPIN.EQ.1) THEN ! SECOND SPIN DIRECTION
                    ORBPRO(2,NB+1:2*NB,IKPT,IORB)=ORBPRO(2,NB+1:2*NB,IKPT,IORB)&
     &                                       +CFAC*STATE%VEC(1,IPRO,:)
                    END IF
                  END IF ! NDIM=1 OR 2
                ELSE  !TIME INVERTED WAVE FUNCTION (-K)
                  CFAC=CONJG( NEWORBITAL(IORB)%ENTRY(IENTRY)%ORB(LMN) &
     &                       *CONJG(EIKT))
                  IF(NDIM.EQ.2) THEN
                     ORBPRO(1,:,IKPT,IORB)=ORBPRO(1,:,IKPT,IORB) &
     &                                 -CFAC*CONJG(STATE%VEC(2,IPRO,:))
                     ORBPRO(2,:,IKPT,IORB)=ORBPRO(2,:,IKPT,IORB) &
     &                                 +CFAC*CONJG(STATE%VEC(1,IPRO,:))
                   ELSE IF(NDIM.EQ.1) THEN
                     IF(ISPIN.EQ.1) THEN ! FIRST SPIN DIRECTION
                       ORBPRO(1,:NB,IKPT,IORB)=ORBPRO(1,:NB,IKPT,IORB) &
     &                                        +CFAC*CONJG(STATE%VEC(1,IPRO,:))
                    END IF
                    IF(ISPIN.EQ.2.OR.NSPIN.EQ.1) THEN ! SECOND SPIN DIRECTION
                      ORBPRO(2,NB+1:2*NB,IKPT,IORB) &
     &                                       =ORBPRO(2,NB+1:2*NB,IKPT,IORB) &
     &                                       +CFAC*CONJG(STATE%VEC(1,IPRO,:))
                    END IF
                  END IF ! NDIM=1 OR 2
                END IF  ! TINV=T OR F
              ENDDO  ! LMN
            ENDDO ! IENTRY
          ENDDO ! IORB
        ENDDO ! ISPIN
      ENDDO ! IKPT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWORBITAL_COMPLETEIPRO()
!     **************************************************************************
!     ** CREATES FOR EACH ORBITAL AN ARRAY IPRO POINTING TO THE POSITION      **
!     ** ON THE STATE ARRAY                                                   **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NAT &         ! #(ATOMS)
     &                          ,ISPECIES &    !(NAT)  SPECIES INDEX
     &                          ,LNX &         !(NSP)
     &                          ,LOX           !(LNXX,NSP) MAIN ANGULAR MOMENTUM
      USE NEWORBITAL_MODULE, ONLY : NORB,NEWORBITAL
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NX=3
      INTEGER(4),PARAMETER :: LX=3
      INTEGER(4)           :: IPOINT(NX,LX+1,NAT)
      INTEGER(4)           :: IORB,IENTRY,IAT,IT(3),LMNX,L,M,LMN,N
!     **************************************************************************
!
!     ==========================================================================
!     ==  FILL POINT ARRAY IPOINT(N,L+1,NAT)                                  ==
!     ==  IPOINT CONTAINS THE FIRST VALUE IPRO FOR A A GIVEN (N,L,IAT)        ==
!     ==========================================================================
      CALL MAKEIPOINT(NX,LX,NAT,IPOINT)
!
!     ==========================================================================
!     == COMPLETE IPRO ARRAY                                                  ==
!     ==========================================================================
      DO IORB=1,NORB
        DO IENTRY=1,NEWORBITAL(IORB)%NENTRY
          IAT=NEWORBITAL(IORB)%ENTRY(IENTRY)%IAT
          LMNX=NEWORBITAL(IORB)%ENTRY(IENTRY)%LMNX
          IF(ALLOCATED(NEWORBITAL(IORB)%ENTRY(IENTRY)%IPRO)) CYCLE
          ALLOCATE(NEWORBITAL(IORB)%ENTRY(IENTRY)%IPRO(LMNX))
          NEWORBITAL(IORB)%ENTRY(IENTRY)%IPRO(:)=0
          LMN=0
          DO L=0,LX+1
            DO M=1,2*L+1
              LMN=LMN+1
              IF(LMN.GT.LMNX) EXIT
              N=1  ! MORE THAN ONE PROJECTOR NOT YET SUPPORTED!!!
              IF(ABS(NEWORBITAL(IORB)%ENTRY(IENTRY)%ORB(LMN)).LT.1.D-7) CYCLE
              NEWORBITAL(IORB)%ENTRY(IENTRY)%IPRO(LMN)=IPOINT(N,L+1,IAT)-1+M
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE NEWSET_MODULE
TYPE SET_TYPE
  CHARACTER(32)             :: ID
  INTEGER(4)                :: NCOOP=0
  CHARACTER(32),ALLOCATABLE :: COOP(:,:)    !(2,NCOOP)
  INTEGER(4)   ,ALLOCATABLE :: ICOOP(:,:)   !(2,NCOOP)
  INTEGER(4)                :: NORB=0
  CHARACTER(32),ALLOCATABLE :: ORB(:)       !(NORB)
  INTEGER(4)   ,ALLOCATABLE :: IORB(:)      !(NORB)
  INTEGER(4)                :: NWGHT=0
  INTEGER(4)   ,ALLOCATABLE :: IAT(:)       !(NWGHT)
  INTEGER(4)   ,ALLOCATABLE :: L(:)         !(NWGHT)
  CHARACTER(32)             :: SPECIAL=' '  ! SPECIAL TYPE LIKE TOTAL AND EMPTY
  CHARACTER(32)             :: LEGEND=' '   ! TEXT FOR ANNOTATION IN THE FIGURE
  CHARACTER(32)             :: SPINID='+Z'  ! ID FOR SPIN PROJECTION
END TYPE SET_TYPE
INTEGER(4)     ,PARAMETER   :: NSETX=1000
INTEGER(4)                  :: NSET=0
INTEGER(4)                  :: ISET=0
TYPE(SET_TYPE)              :: NEWSET(NSETX) 
END MODULE NEWSET_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$GETI4(ID,IVAL)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : NSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: IVAL
!     **************************************************************************
      IF(ID.EQ.'NSET') THEN
        IVAL=NSET
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('NEWSET$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$SETCH(ID,VAL)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : NSET,ISET,NEWSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'SPECIAL') THEN
        IF(ISET.EQ.0) THEN
          CALL ERROR$MSG('NO SET SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('NEWSET$SETCH')
        END IF
        NEWSET(ISET)%SPECIAL=VAL
!
      ELSE IF(ID.EQ.'SPINID') THEN
        IF(ISET.EQ.0) THEN
          CALL ERROR$MSG('NO SET SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('NEWSET$SETCH')
        END IF
        NEWSET(ISET)%SPINID=VAL
!
      ELSE IF(ID.EQ.'LEGEND') THEN
        IF(ISET.EQ.0) THEN
          CALL ERROR$MSG('NO SET SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('NEWSET$SETCH')
        END IF
        NEWSET(ISET)%LEGEND=VAL
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('NEWSET$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$REPORT(NFIL)
!     **************************************************************************
!     **  REPORT ALL DEFINED SETS                                             **
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : NSET,NEWSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NCOOP
      INTEGER(4)            :: NWGHT
      INTEGER(4)            :: NORB
      INTEGER(4)            :: I,J,K,J1,J2
      CHARACTER(16)         :: TYPEID
      CHARACTER(218)        :: FMT
!     **************************************************************************
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='(80("="),T20,"  NEWSET REPORT  ")')
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='(40("."),":",T1,A,T42,I3)')'NUMBER OF SETS',NSET
      DO I=1,NSET
!
!       ========================================================================
!       == HEADER LINE AND GENERIC INFORMATION                                ==
!       ========================================================================
        FMT='(80("="))'
        WRITE(NFIL,FMT=FMT) 
        FMT='(20("."),":",T1,A,T25,A)'
        WRITE(NFIL,FMT=FMT)'ID',TRIM(NEWSET(I)%ID)
        WRITE(NFIL,FMT=FMT)'LEGEND',TRIM(NEWSET(I)%LEGEND)
        WRITE(NFIL,FMT=FMT)'SPINID',NEWSET(I)%SPINID
!
!       ========================================================================
!       == REPORT FOR SPECIFIC TYPES                                          ==
!       ========================================================================
        NCOOP=NEWSET(I)%NCOOP
        NORB=NEWSET(I)%NORB
        NWGHT=NEWSET(I)%NWGHT
        IF(LEN_TRIM(NEWSET(I)%SPECIAL).NE.0) THEN
          WRITE(NFIL,FMT=FMT)'TYPE','SPECIAL:'//NEWSET(I)%SPECIAL
!
        ELSE IF(NCOOP.GT.0) THEN 
          WRITE(NFIL,FMT=FMT)'TYPE','COOP'
          DO J=1,NCOOP
            WRITE(NFIL,FMT='("ORBITAL1=",A32," ORBITAL2=",A21)') &
    &                      NEWSET(I)%COOP(:,J)
          ENDDO
!
        ELSE IF(NWGHT.GT.0) THEN
          WRITE(NFIL,FMT=FMT)'TYPE','ANGULAR-MOMENTUM WEIGHTS'
!         == 5 ENTRIES PER LINE=================================================
          FMT='(5(" [IAT=",I5," L=",I1,"] "))'
          DO J=1,INT(REAL(NWGHT)/5.)
            J1=1+5*(J-1)
            J2=5*(J-1)
            WRITE(NFIL,FMT=FMT)(NEWSET(I)%IAT(K),NEWSET(I)%L(K),K=J1,J2)
          ENDDO
          J1=1+5*INT(REAL(NWGHT)/5.)
          FMT='(T1'
          DO J=J1,NWGHT
            FMT=TRIM(ADJUSTL(FMT))//'," [IAT=",I5," L=",I1,"] "'
          ENDDO
          FMT=TRIM(ADJUSTL(FMT))//')'
          WRITE(NFIL,FMT=FMT)(NEWSET(I)%IAT(J),NEWSET(I)%L(J),J=J1,NWGHT)
!
        ELSE IF(NORB.GT.0) THEN
          WRITE(NFIL,FMT=FMT)'TYPE','ORBITAL-WEIGHTS'
          WRITE(NFIL,FMT='("ORB=",A,T40," ORB=",A)') &
    &                                          (TRIM(NEWSET(I)%ORB(J)),J=1,NORB)
!
        ELSE
          CALL ERROR$MSG('TYPE ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('TYPEID',TYPEID)
          CALL ERROR$STOP('NEWSET$REPORT')
        ENDIF
      ENDDO
!
      WRITE(NFIL,*)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$NEWSET(ID)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **  THE NEW ORBITAL HAS A NAME BUT NO ENTRIES YET                       **
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : NSETX,NSET,NEWSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)              :: I
!     **************************************************************************
!
!     ==========================================================================
!     == CHECK IF NAME HAS BEEN USED BEFORE                                   ==
!     ==========================================================================
      DO I=1,NSET
        IF(ID.EQ.NEWSET(I)%ID) THEN
          CALL ERROR$MSG('SET WITH THE NAME ALREADY EXISTS')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('NEWSET$NEWSET')
        END IF
      ENDDO
!
!     ==========================================================================
!     == CHECK IF ARRAY SIZE IS SUFFICIENT                                    ==
!     ==========================================================================
      IF(NSET+1.GT.NSETX) THEN
        CALL ERROR$MSG('MAX. NR. OF SETS EXCEEDED')
        CALL ERROR$MSG('INCREASE PARAMETER "NSETX" IN NEWSET_MODULE')
        CALL ERROR$I4VAL('NSETX',NSETX)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('NEWSET$NEWSET')
      END IF
!
!     ==========================================================================
!     == CHECK IF ARRAY SIZE IS SUFFICIENT                                    ==
!     ==========================================================================
      NSET=NSET+1
      NEWSET(NSET)%ID=ID
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$SELECT(NAME)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : ISET,NSET,NEWSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)              :: I
!     **************************************************************************
      IF(NAME.EQ.'NONE') THEN
        ISET=0
        RETURN
      END IF
!
!     ==========================================================================
!     == CHECK IF NAME HAS BEEN USED BEFORE                                   ==
!     ==========================================================================
      IF(ISET.NE.0) THEN
         CALL ERROR$MSG('ATTEMPT TO SELECT SER WHILE ANOTHER IS ACTIVE')
         CALL ERROR$MSG('UNSELECT FIRST')
         CALL ERROR$CHVAL('NAME',NAME)
         CALL ERROR$STOP('NEWSET$SELECT')
      END IF
     
      DO I=1,NSET
        IF(NAME.EQ.NEWSET(I)%ID) THEN
          ISET=I 
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('SET NOT FOUND')
      CALL ERROR$CHVAL('NAME',NAME)
      CALL ERROR$STOP('NEWSET$SELECT')
      STOP
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$ISELECT(ISET1)
!     **************************************************************************
!     **  CREATE NEW ORBITAL                                                  **
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : ISET,NSET,NEWSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISET1
!     **************************************************************************
      IF(ISET1.EQ.0) THEN
        ISET=0
        RETURN
      END IF
!
!     ==========================================================================
!     == CHECK IF ORBITAL POINTER IS IN ALLOWED RANGE                         ==
!     ==========================================================================
      IF(ISET1.LT.1.OR.ISET1.GT.NSET) THEN
        CALL ERROR$MSG('SPECIFIED SET INDEX OUT OF RANGE')
        CALL ERROR$I4VAL('ISET1',ISET1)
        CALL ERROR$I4VAL('NSET',NSET)
        CALL ERROR$STOP('NEWSET$ISELECT')
      END IF
!
!     ==========================================================================
!     == CHECK IF NAME HAS BEEN USED BEFORE                                   ==
!     ==========================================================================
      IF(ISET.NE.0) THEN
        CALL ERROR$MSG('ATTEMPT TO SELECT SET WHILE ANOTHER IS ACTIVE')
        CALL ERROR$MSG('UNSELECT FIRST')
        CALL ERROR$CHVAL('NAME ACTIVE',NEWSET(ISET)%ID)
        CALL ERROR$CHVAL('NAME REQUESTED',NEWSET(ISET1)%ID)
        CALL ERROR$STOP('NEWSET$ISELECT')
      END IF
!
!     ==========================================================================
!     == SET ORBITAL POINTER                                                  ==
!     ==========================================================================
      ISET=ISET1
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$ADDLWEIGHT(IAT,L)
!     **************************************************************************
!     ** ADD AN ANGULAR MOMENTUM WEIGHT TO THE SET
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : ISET,NEWSET
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: IAT
      INTEGER(4)  ,INTENT(IN) :: L
      INTEGER(4)  ,PARAMETER  :: NWGHTX=10
      INTEGER(4)              :: NWGHT
      INTEGER(4) ,ALLOCATABLE :: TMP(:)
!     **************************************************************************
!
!     ==========================================================================
!     == CHECK IF SET IS SELECTED                                             ==
!     ==========================================================================
      IF(ISET.EQ.0) THEN
        CALL ERROR$MSG('NO SET SELECTED')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$STOP('NEWSET$ADDLWGHT')
      END IF
!
!     ==========================================================================
!     == RESIZE ARRAY                                                         ==
!     ==========================================================================
      NWGHT=NEWSET(ISET)%NWGHT
      IF(NWGHT.EQ.0) THEN
        ALLOCATE(NEWSET(ISET)%IAT(NWGHTX))
        ALLOCATE(NEWSET(ISET)%L(NWGHTX))
      ELSE
        IF(NWGHT.GE.SIZE(NEWSET(ISET)%IAT)) THEN
          ALLOCATE(TMP(NWGHT))
!
          TMP=NEWSET(ISET)%IAT(:NWGHT)
          DEALLOCATE(NEWSET(ISET)%IAT)
          ALLOCATE(NEWSET(ISET)%IAT(NWGHT+NWGHTX))
          NEWSET(ISET)%IAT(:NWGHT)=TMP
          NEWSET(ISET)%IAT(NWGHT+1:)=-99
!
          TMP=NEWSET(ISET)%L(:NWGHT)
          DEALLOCATE(NEWSET(ISET)%L)
          ALLOCATE(NEWSET(ISET)%L(NWGHT+NWGHTX))
          NEWSET(ISET)%L(:NWGHT)=TMP
          NEWSET(ISET)%L(NWGHT+1:)=-99
! 
          DEALLOCATE(TMP)
        END IF
      END IF
!
!     ==========================================================================
!     == SET DATA                                                             ==
!     ==========================================================================
      NWGHT=NWGHT+1
      NEWSET(ISET)%NWGHT=NWGHT
      NEWSET(ISET)%IAT(NWGHT)=IAT
      NEWSET(ISET)%L(NWGHT)=L
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$ADDORBWGHT(ORBITAL)
!     **************************************************************************
!     ** ADD AN ORBITAL TO THE SELECTED SET
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : ISET,NEWSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ORBITAL
      INTEGER(4)  ,PARAMETER  :: NORBX=10
      INTEGER(4)              :: NORB
      CHARACTER(32),ALLOCATABLE :: TMP(:)
!     **************************************************************************
!
!     ==========================================================================
!     == CHECK IF SET IS SELECTED                                             ==
!     ==========================================================================
      IF(ISET.EQ.0) THEN
        CALL ERROR$MSG('NO SET SELECTED')
        CALL ERROR$CHVAL('ORBITAL',ORBITAL)
        CALL ERROR$STOP('NEWSET$ADDORBWGHT')
      END IF
!
!     ==========================================================================
!     == RESIZE ARRAY                                                         ==
!     ==========================================================================
      NORB=NEWSET(ISET)%NORB
      IF(NORB.EQ.0) THEN
        ALLOCATE(NEWSET(ISET)%ORB(NORBX))
      ELSE
        IF(NORB.GE.SIZE(NEWSET(ISET)%ORB)) THEN
          ALLOCATE(TMP(NORB))
          TMP=NEWSET(ISET)%ORB(:NORB)
          DEALLOCATE(NEWSET(ISET)%ORB)
          ALLOCATE(NEWSET(ISET)%ORB(NORB+NORBX))
          NEWSET(ISET)%ORB(:NORB)=TMP
          NEWSET(ISET)%ORB(NORB+1:)=' '
          DEALLOCATE(TMP)
        END IF
      END IF
!
!     ==========================================================================
!     == SET DATA                                                             ==
!     ==========================================================================
      NORB=NORB+1
      NEWSET(ISET)%NORB=NORB
      NEWSET(ISET)%ORB(NORB)=ORBITAL
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$ADDCOOP(ORBITAL1,ORBITAL2)
!     **************************************************************************
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : ISET,NEWSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ORBITAL1
      CHARACTER(*),INTENT(IN) :: ORBITAL2
      INTEGER(4)  ,PARAMETER  :: NCOOPX=10
      INTEGER(4)              :: NCOOP
      CHARACTER(32),ALLOCATABLE :: TMP(:,:)
!     **************************************************************************
!
!     ==========================================================================
!     == CHECK IF SET IS SELECTED                                             ==
!     ==========================================================================
      IF(ISET.EQ.0) THEN
        CALL ERROR$MSG('NO SET SELECTED')
        CALL ERROR$CHVAL('ORBITAL1',ORBITAL1)
        CALL ERROR$CHVAL('ORBITAL2',ORBITAL2)
        CALL ERROR$STOP('NEWSET$ADDCOOP')
      END IF
!
!     ==========================================================================
!     == RESIZE ARRAY                                                         ==
!     ==========================================================================
      NCOOP=NEWSET(ISET)%NCOOP
      IF(NCOOP.EQ.0) THEN
        ALLOCATE(NEWSET(ISET)%COOP(2,NCOOPX))
      ELSE
        IF(NCOOP.GE.SIZE(NEWSET(ISET)%COOP,2)) THEN
          ALLOCATE(TMP(2,NCOOP))
          TMP=NEWSET(ISET)%COOP(:,:NCOOP)
          DEALLOCATE(NEWSET(ISET)%COOP)
          ALLOCATE(NEWSET(ISET)%COOP(2,NCOOP+NCOOPX))
          NEWSET(ISET)%COOP(:,:NCOOP)=TMP
          NEWSET(ISET)%COOP(:,NCOOP+1:)=' '
          DEALLOCATE(TMP)
        END IF
      END IF
!
!     ==========================================================================
!     == SET DATA                                                             ==
!     ==========================================================================
      NCOOP=NCOOP+1
      NEWSET(ISET)%NCOOP=NCOOP
      NEWSET(ISET)%COOP(1,NCOOP)=ORBITAL1
      NEWSET(ISET)%COOP(2,NCOOP)=ORBITAL2
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$CLEANUP()
!     **************************************************************************
!     **  TRANSLATES ORBITAL NAMES INTO ORBITAL INDICES                       **
!     **************************************************************************
      USE NEWORBITAL_MODULE, ONLY : NORB &
     &                             ,NEWORBITAL
      USE NEWSET_MODULE, ONLY : NSET &
     &                         ,NEWSET
      USE SPINDIR_MODULE, ONLY : SPINDIR ! (3,NAT)
      IMPLICIT NONE
      INTEGER(4)       :: ISET
      INTEGER(4)       :: NCOOP
      INTEGER(4)       :: NSETORB
      INTEGER(4)       :: I,J,IORB
!     **************************************************************************
      DO ISET=1,NSET
!
!       ========================================================================
!       == TRANSLATE ORBITAL NAMES %COOP INTO ORBITAL INDICES %ICOOP          ==
!       ========================================================================
        NCOOP=NEWSET(ISET)%NCOOP
        IF(NCOOP.GT.0) THEN
          ALLOCATE(NEWSET(ISET)%ICOOP(2,NCOOP))
          NEWSET(ISET)%ICOOP(:,:)=0
          DO I=1,NCOOP
            DO J=1,2
              DO IORB=1,NORB
                IF(NEWSET(ISET)%COOP(J,I).NE.NEWORBITAL(IORB)%NAME) CYCLE
                NEWSET(ISET)%ICOOP(J,I)=IORB
                EXIT
              ENDDO  !IORB
              IF(NEWSET(ISET)%ICOOP(J,I).EQ.0) THEN
                CALL ERROR$MSG('ORBITAL NOT FOUND')
                CALL ERROR$STOP('NEWSET$CLEANUP')
              END IF
            ENDDO
          ENDDO
        END IF
!
!       ========================================================================
!       == TRANSLATE ORBITAL NAMES %ORB INTO ORBITAL INDICES %IORB            ==
!       ========================================================================
        NSETORB=NEWSET(ISET)%NORB
        IF(NSETORB.GT.0) THEN
          ALLOCATE(NEWSET(ISET)%IORB(NSETORB))
          NEWSET(ISET)%IORB(:)=0
          DO I=1,NSETORB
            DO IORB=1,NORB
              IF(NEWSET(ISET)%ORB(I).NE.NEWORBITAL(IORB)%NAME) CYCLE
              NEWSET(ISET)%IORB(I)=IORB
              EXIT
            ENDDO  !IORB
            IF(NEWSET(ISET)%IORB(I).EQ.0) THEN
              CALL ERROR$MSG('ORBITAL NOT FOUND')
              CALL ERROR$STOP('NEWSET$CLEANUP')
            END IF
          ENDDO !I
        END IF
      ENDDO !ISET        
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET$PROCESS(NBB,NKPT,NSET,SET)
!     **************************************************************************
!     ** CONSTRUCTS THE NEW SETS, THAT IS THE CONTRIBUTION FROM EACH STATE    **
!     ** TO THE PROJECTED DENSITY OF STATES. THERE IS ONE SET FOR EACH        **
!     ** PROJECTED DENSITY OF STATES                                          **
!     **                                                                      **
!     ** THIS ROUTINE FIRST CONSTRUCTS ALL SPIN COMPONENTS FROM A SET OF      **
!     ** ORBITALS AND THEN PROJECTS ONTO SPIN CONTRIBUTIONS                   **
!     **                                                                      **
!     ** THE INFORMATION ON THE WAVE FNCTIONS IS ENCODED IN                   **
!     **   STATE(IKPT,ISPIN)%NB
!     **                    %VEC(IDIM,IPRO,IB)
!     **                    %OCC(IB)
!     **                    %EIG(IB)
!     **                                                                      **
!     **                                                                      **
!     **  TREATMENT OF SPINS:                                                 **
!     **    NON SPIN-POLARIZED:            NDIM=1,NSPIN=1                     **
!     **    COLLINEAR, SPIN POLARIZED:     NDIM=1,NSPIN=2                     **
!     **    NON-COLLINEAR, SPIN POLARIZED: NDIM=1,NSPIN=2                     **
!     **                                                                      **
!     **    THE VARIABLE ORBPRO AND ORBPROCC ARE FORMULATED IN THE            **
!     **    NON-COLLINEARMODEL: 
!     **    - THE BAND INDEX (NBB) PER K-POINT COUNTS ALWAYS SPIN-UP AND      **
!     **      SPIN-DOWN STATES                                                **
!     **    - EACH STATE HAS SPIN-UP AND SPIN-DOWN COMPONENTS                 **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NDIM &
     &                          ,NSPIN &
     &                          ,NKPT1=>NKPT &
     &                          ,STATEARR !(IKPT,ISPIN)
      USE NEWSET_MODULE, ONLY : NEWSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NBB
      INTEGER(4),INTENT(IN)  :: NKPT
      INTEGER(4),INTENT(IN)  :: NSET
      REAL(8)   ,INTENT(OUT) :: SET(NBB,NKPT,2,NSET)
      INTEGER(4)             :: NBB1
      INTEGER(4)             :: NORB
      INTEGER(4)             :: NSET1
      INTEGER(4)             :: IKPT,ISPIN,I,ISET
      REAL(8)   ,ALLOCATABLE :: MATEL(:,:,:,:)
      COMPLEX(8),ALLOCATABLE :: ORBPRO(:,:,:,:)   !(2,NBB,NKPT,NORB))
      COMPLEX(8),ALLOCATABLE :: ORBPROCC(:,:,:,:) !(2,NBB,NKPT,NORB))
!     **************************************************************************
!
!     ==========================================================================
!     ==  ALLOCATE MATRIX ELEMENTS MATEL(NDIM,NBB,NKPT,NSET)                  ==
!     ==========================================================================
      NBB1=0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          NBB1=MAX(NBB1,STATEARR(1,1)%NB)
        ENDDO
      ENDDO
      IF(NDIM.EQ.1)NBB1=2*NBB1
      IF(NBB1.NE.NBB) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$STOP('NEWSET$PROCESS')
      END IF

      CALL NEWSET$GETI4('NSET',NSET1)
      IF(NSET1.NE.NSET) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$STOP('NEWSET$PROCESS')
      END IF
!
!     ==========================================================================
!     ==  RESOLVE ORBITALS                                                    ==
!     ==========================================================================
      CALL NEWORBITAL$GETI4('NORB',NORB)
      ALLOCATE(ORBPRO(2,NBB,NKPT,NORB))
      ALLOCATE(ORBPROCC(2,NBB,NKPT,NORB))
      CALL NEWORBITAL$ORBPRO(.FALSE.,NBB,NKPT,NORB,ORBPRO)
      CALL NEWORBITAL$ORBPRO(.TRUE.,NBB,NKPT,NORB,ORBPROCC)
!
!     ==========================================================================
!     ==  RESOLVE SETS                                                        ==
!     ==========================================================================
      CALL NEWSET$CLEANUP()
      ALLOCATE(MATEL(4,NBB,NKPT,NSET))
      CALL NEWSET_COLLECTMATEL(NBB,NKPT,NSET,NORB,ORBPRO,ORBPROCC,MATEL)
!
!     ==========================================================================
!     ==  PROJECT ONTO UP AND DOWN SPIN                                       ==
!     ==========================================================================
      CALL NEWSET_PROJECTSPIN(NBB,NKPT,NSET,MATEL)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO I=1,2
        SET(:,:,I,:)=MATEL(I,:,:,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET_PROJECTSPIN(NBB,NKPT,NSET,MATEL)
!     **************************************************************************
!     ** RESOLVE SPIN PROJECTIONS                                             **
!     **                                                                      **
!     ** ON INPUT, MATEL HOLDS THE (RHOT,SPIN_X,SPIN_Y,SPIN_Z) COMPONENTS     **
!     **                                                                      **
!     ** ON OUTPUT:                                                           **
!     ** - THE FIRST TWO ELEMENTS OF MATEL ARE THE SPIN-UP AND SPIN-DOWN      **
!     **   COMPONENTS FOR THE SELECTED SPIN AXIS                              **
!     ** - FOR COOPS, ONLY ONE THE SELECTED SPIN-UP DIRECTION IS MAINTAINED.  **
!     ** - FOR SPINID='TOTAL' THE TWO SPIN DIRECTIONS ARE SUMMED IN THE FIRST **
!     **   AND THE SECOND IS SET TO ZERO.                                     **
!     ** - THE SPIN DIRECTION 'MAIN' IS ONLY ALLOWED IF ANGULAR MOMENTUM      **
!     **   WEIGHTS ON A SINGLE ATOM ARE SELECTED                              **
!     **                                                                      **
!     ** SPINID MAY HAVE THE VALUES:                                          **
!     **   'TOTAL','MAIN','X','+X','-X','Y','+Y','-Y','Z','+Z','-Z'           **
!     **   AND A STRING OF THREE REAL NUMBERS IN FREE FORMAT 'X Y Z'          **
!     **   DEFINING THE SPIN DIRECTION                                        **
!     ************************************P. BLOECHL, GOSLAR, DEC. 30,2015******
      USE NEWSET_MODULE, ONLY : NSET1=>NSET &
     &                         ,NEWSET
      USE SPINDIR_MODULE, ONLY : SPINDIR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NBB
      INTEGER(4),INTENT(IN) :: NKPT
      INTEGER(4),INTENT(IN) :: NSET
      REAL(8)   ,INTENT(OUT):: MATEL(4,NBB,NKPT,NSET)
      INTEGER(4)            :: ISET,I
      INTEGER(4)            :: IATP
      CHARACTER(32)         :: SPINID
      REAL(8)               :: SPINVEC(3)
!     **************************************************************************
      IF(NSET1.NE.NSET) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZES')
        CALL ERROR$STOP('NEWSET_PROJECTSPIN')
      END IF
      DO ISET=1,NSET
        SPINID=NEWSET(ISET)%SPINID
!       == TOTAL SPIN ==========================================================
        IF(SPINID.EQ.'TOTAL') THEN
          MATEL(1,:,:,ISET)=2.D0*MATEL(1,:,:,ISET)
          MATEL(2:,:,:,ISET)=0.D0
!
!       == SPECIFIED SPIN DIRECTIONS ALONG CARTESIAN AXES ======================
        ELSE IF(SPINID.EQ.'Z'.OR.SPINID.EQ.'+Z'.OR.SPINID.EQ.'-Z'.OR. &
     &          SPINID.EQ.'X'.OR.SPINID.EQ.'+X'.OR.SPINID.EQ.'-X'.OR. &
     &          SPINID.EQ.'Y'.OR.SPINID.EQ.'+Y'.OR.SPINID.EQ.'-Y') THEN
          IF(SPINID.EQ.'Z'.OR.SPINID.EQ.'+Z') THEN
          ELSE IF(SPINID.EQ.'-Z') THEN
            MATEL(4,:,:,ISET)=-MATEL(4,:,:,ISET)
          ELSE IF(SPINID.EQ.'X'.OR.SPINID.EQ.'+X') THEN
            MATEL(4,:,:,ISET)=MATEL(2,:,:,ISET)
          ELSE IF(SPINID.EQ.'-X') THEN
            MATEL(4,:,:,ISET)=-MATEL(2,:,:,ISET)
          ELSE IF(SPINID.EQ.'Y'.OR.SPINID.EQ.'+Y') THEN
            MATEL(4,:,:,ISET)=MATEL(3,:,:,ISET)
          ELSE IF(SPINID.EQ.'-Y') THEN
            MATEL(4,:,:,ISET)=-MATEL(3,:,:,ISET)
          END IF
          MATEL(2,:,:,ISET)=MATEL(1,:,:,ISET)-MATEL(4,:,:,ISET)
          MATEL(1,:,:,ISET)=MATEL(1,:,:,ISET)+MATEL(4,:,:,ISET)
!
!       == MAIN SPIN DIRECTION AS SPIN AXIS ====================================
        ELSE IF(SPINID.EQ.'MAIN') THEN
          IF(NEWSET(ISET)%NWGHT.EQ.0) THEN
            CALL ERROR$MSG('SPIN="MAIN" ONLY COMPATIBLE')
            CALL ERROR$MSG('WITH ANGULAR MOMENTUM WEIGHTS')
            CALL ERROR$STOP('NEWSET_PROJECTSPIN')
          END IF
          IATP=NEWSET(ISET)%IAT(1)
          IF(IATP.EQ.-1) THEN
            CALL ERROR$MSG('SPIN="MAIN" ONLY COMPATIBLE')
            CALL ERROR$MSG('WITH WEIGHTS ON A SINGLE ATOM')
            CALL ERROR$STOP('NEWSET_PROJECTSPIN')
          END IF
          DO I=1,NEWSET(ISET)%NWGHT
            IF(NEWSET(ISET)%IAT(I).NE.IATP) THEN
              CALL ERROR$MSG('SPIN="MAIN" ONLY COMPATIBLE')
              CALL ERROR$MSG('WITH WEIGHTS ON A SINGLE ATOM')
              CALL ERROR$STOP('NEWSET_PROJECTSPIN')
            END IF
          ENDDO
          SPINVEC=SPINDIR(:,IATP)
          SPINVEC=SPINVEC/SQRT(SUM(SPINVEC**2))
          MATEL(4,:,:,ISET)=SPINVEC(1)*MATEL(2,:,:,ISET) &
       &                   +SPINVEC(2)*MATEL(3,:,:,ISET) &
       &                   +SPINVEC(3)*MATEL(4,:,:,ISET) 
          MATEL(2,:,:,ISET)=MATEL(1,:,:,ISET)-MATEL(4,:,:,ISET)
          MATEL(1,:,:,ISET)=MATEL(1,:,:,ISET)+MATEL(4,:,:,ISET)
!
!       == SPIN AXIS AS NUMBER TRIPLE ==========================================
        ELSE 
          READ(SPINID,*)SPINVEC
          SPINVEC=SPINVEC/SQRT(SUM(SPINVEC**2))
          MATEL(4,:,:,ISET)=SPINVEC(1)*MATEL(2,:,:,ISET) &
     &                     +SPINVEC(2)*MATEL(3,:,:,ISET) &
     &                     +SPINVEC(3)*MATEL(4,:,:,ISET) 
          MATEL(2,:,:,ISET)=MATEL(1,:,:,ISET)-MATEL(4,:,:,ISET)
          MATEL(1,:,:,ISET)=MATEL(1,:,:,ISET)+MATEL(4,:,:,ISET)
        END IF
        MATEL(3:4,:,:,ISET)=0.D0
!
!       ========================================================================
!       == REMOVE SPIN DOWN DIRECTION FOR COOPS                               ==
!       ========================================================================
        IF(NEWSET(ISET)%NCOOP.GT.0) THEN
          MATEL(2,:,:,ISET)=0.D0
        END IF
      ENDDO
      RETURN
      END
     
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET_COLLECTMATEL(NBB,NKPT,NSET,NORB,ORB,ORBCC,MATEL)
!     **************************************************************************
!     ** THE FIRST INDEX OF MATEL SPECIFIES (RHOT,MAG_X,MAG_Y,MAG_Z)          **
!     ** ORB AND ORBC ARE FORMULATED IN THE NON-COLLINEAR SPIN MODEL          **
!     **************************************************************************
      USE NEWSET_MODULE, ONLY : NSET1=>NSET &
     &                         ,NEWSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NBB
      INTEGER(4),INTENT(IN) :: NKPT
      INTEGER(4),INTENT(IN) :: NORB
      INTEGER(4),INTENT(IN) :: NSET
      COMPLEX(8),INTENT(IN) :: ORB(2,NBB,NKPT,NORB)
      COMPLEX(8),INTENT(IN) :: ORBCC(2,NBB,NKPT,NORB) ! -K
      REAL(8)   ,INTENT(OUT):: MATEL(4,NBB,NKPT,NSET)
      INTEGER(4)            :: IKPT,ISPIN
      INTEGER(4)            :: NCOOP
      INTEGER(4)            :: NSETORB
      INTEGER(4)            :: ISET,ICOOP,I,IORB1,IORB2
!     **************************************************************************
      IF(NSET1.NE.NSET) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZES')
        CALL ERROR$STOP('NEWSET_COLLECTMATEL')
      END IF
!
!     ==========================================================================
!     == CHECK CONSISTENCY                                                    ==
!     ==========================================================================
      DO ISET=1,NSET
        I=0
        IF(NEWSET(ISET)%NWGHT.GT.0) I=I+1
        IF(NEWSET(ISET)%NCOOP.GT.0) I=I+1
        IF(NEWSET(ISET)%NORB.GT.0) I=I+1
        IF(NEWSET(ISET)%SPECIAL.EQ.'TOTAL') I=I+1
        IF(NEWSET(ISET)%SPECIAL.EQ.'EMPTY') I=I+1
        IF(NEWSET(ISET)%SPECIAL.EQ.'ALL') I=I+1
        IF(I.NE.1) THEN
          CALL ERROR$MSG('INCOMPATIBLE SELECTIONS IN NEWSET OBJECT')
          CALL ERROR$I4VAL('ISET',ISET)
          CALL ERROR$CHVAL('NEWSET%ID',NEWSET(ISET)%ID)
          CALL ERROR$I4VAL('NCOOP',NEWSET(ISET)%NCOOP)
          CALL ERROR$I4VAL('NWGHT',NEWSET(ISET)%NWGHT)
          CALL ERROR$I4VAL('NORB',NEWSET(ISET)%NORB)
          CALL ERROR$CHVAL('SPECIAL',NEWSET(ISET)%SPECIAL)
          CALL ERROR$STOP('NEWSET_COLLECTMATEL')
        END IF
      ENDDO
!
!     ==========================================================================
!     == ANGULAR MOMENTUM WEIGHTS AND SPECIAL (TOTAL,ALL,EMPTY)               ==
!     ==========================================================================
      MATEL=0.D0
      CALL NEWSET_ANGMOMWEIGHTS(NBB,NKPT,NSET,MATEL)
!
!     ==========================================================================
!     == COOPS                                                                ==
!     ==========================================================================
      DO ISET=1,NSET
        NCOOP=NEWSET(ISET)%NCOOP
        DO ICOOP=1,NCOOP
          IORB1=NEWSET(ISET)%ICOOP(1,ICOOP)
          IORB2=NEWSET(ISET)%ICOOP(2,ICOOP)
!         == 0.5*TRACE(D*SIGMA)=0.5*SUM_I,J D_JI * SIGMA_IJ ====================
!         == D_IJ=<I|PSI1>... <PSI2|J> =========================================
!         == REAL PART TO ACCOUNT FOR -K (WAVE VECTOR)
!         == ORBCC ACCOUNTS FOR -K (AVERAGE FO +K AND -K IS CALCULATED)
!         == REAL(C) INHERITS THE KIND PARAMETER IF THE ARGUMENT C IS COMPLEX
          MATEL(1,:,:,ISET)=MATEL(1,:,:,ISET) &
     &               +0.25D0*REAL(ORB(1,:,:,IORB1)*CONJG(ORB(1,:,:,IORB2)) &
     &                           +ORB(2,:,:,IORB1)*CONJG(ORB(2,:,:,IORB2)))&
     &               +0.25D0*REAL(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB2)) &
     &                           +ORBCC(2,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB2)))
          MATEL(2,:,:,ISET)=MATEL(2,:,:,ISET) &
     &               +0.25D0*REAL(ORB(1,:,:,IORB1)*CONJG(ORB(2,:,:,IORB2)) &
     &                          +ORB(2,:,:,IORB1)*CONJG(ORB(1,:,:,IORB2))) &
     &               +0.25D0*REAL(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB2)) &
     &                           +ORBCC(2,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB2)))
          MATEL(3,:,:,ISET)=MATEL(3,:,:,ISET) &
     &               +0.25D0*AIMAG(ORB(1,:,:,IORB1)*CONJG(ORB(2,:,:,IORB2)) &
     &                            -ORB(2,:,:,IORB1)*CONJG(ORB(1,:,:,IORB2))) &
     &               +0.25D0*AIMAG(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB2))&
     &                            -ORBCC(2,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB2)))
          MATEL(4,:,:,ISET)=MATEL(4,:,:,ISET) &
     &               +0.25D0*REAL(ORB(1,:,:,IORB1)*CONJG(ORB(1,:,:,IORB2)) &
     &                           -ORB(2,:,:,IORB1)*CONJG(ORB(2,:,:,IORB2))) &
     &               +0.25D0*REAL(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB2)) &
     &                           -ORBCC(2,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB2)))
        ENDDO !ICOOP
      ENDDO
!
!     ==========================================================================
!     == ORBITAL WEIGHTS                                                      ==
!     ==========================================================================
      DO ISET=1,NSET
        NSETORB=NEWSET(ISET)%NORB
        DO I=1,NSETORB
          IORB1=NEWSET(ISET)%IORB(I)
          MATEL(1,:,:,ISET)=MATEL(1,:,:,ISET) &
     &               +0.25D0*REAL(ORB(1,:,:,IORB1)*CONJG(ORB(1,:,:,IORB1)) &
     &                           +ORB(2,:,:,IORB1)*CONJG(ORB(2,:,:,IORB1))) &
     &               +0.25D0*REAL(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB1)) &
     &                           +ORBCC(2,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB1)))
          MATEL(2,:,:,ISET)=MATEL(2,:,:,ISET) &
     &               +0.25D0*REAL(ORB(1,:,:,IORB1)*CONJG(ORB(2,:,:,IORB1)) &
     &                          +ORB(2,:,:,IORB1)*CONJG(ORB(1,:,:,IORB1))) &
     &               +0.25D0*REAL(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB1)) &
     &                           +ORBCC(2,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB1)))
          MATEL(3,:,:,ISET)=MATEL(3,:,:,ISET) &
     &              +0.25D0*AIMAG(ORB(1,:,:,IORB1)*CONJG(ORB(2,:,:,IORB1)) &
     &                           -ORB(2,:,:,IORB1)*CONJG(ORB(1,:,:,IORB1))) &
     &              +0.25D0*AIMAG(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB1)) &
     &                           -ORBCC(2,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB1)))
          MATEL(4,:,:,ISET)=MATEL(4,:,:,ISET) &
     &               +0.25D0*REAL(ORB(1,:,:,IORB1)*CONJG(ORB(1,:,:,IORB1)) &
     &                           -ORB(2,:,:,IORB1)*CONJG(ORB(2,:,:,IORB1))) &
     &               +0.25D0*REAL(ORBCC(1,:,:,IORB1)*CONJG(ORBCC(1,:,:,IORB1)) &
     &                           -ORBCC(2,:,:,IORB1)*CONJG(ORBCC(2,:,:,IORB1)))
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWSET_ANGMOMWEIGHTS(NBB,NKPT,NSET,MATEL)
!     **************************************************************************
!     ** STATE%VEC(NDIM,NPRO,NB)
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NAT &
     &                          ,NDIM &
     &                          ,NSPIN &
     &                          ,NKPT1=>NKPT &
     &                          ,ISPECIES & !(NAT)
     &                          ,LNX & !(NSP)
     &                          ,LOX & !(LNX,NSP)
     &                          ,NPRO &
     &                          ,STATEARR,STATE  !(IKPT,ISPIN)
      USE NEWSET_MODULE, ONLY : NSET1=>NSET &
     &                         ,NEWSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NBB
      INTEGER(4),INTENT(IN) :: NKPT
      INTEGER(4),INTENT(IN) :: NSET
      REAL(8)   ,INTENT(OUT):: MATEL(4,NBB,NKPT,NSET)
      INTEGER(4)            :: NWGHT
      INTEGER(4)            :: NB
      INTEGER(4)            :: IATP,LP,NP
      INTEGER(4)            :: ISET,IWGHT,IPRO,IM,IKPT,ISPIN,LN
      INTEGER(4)            :: ISP,IAT,L,N
!     **************************************************************************
      IF(NKPT1.NE.NKPT) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZES')
        CALL ERROR$STOP('NEWSET$ANGMOMWEIGHTS')
      END IF
      IF(NSET1.NE.NSET) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZES')
        CALL ERROR$STOP('NEWSET$ANGMOMWEIGHTS')
      END IF
      MATEL(:,:,:,:)=0.D0
!
!     ==========================================================================
!     == PROCESS SPECIAL                                                      ==
!     ==========================================================================
      DO ISET=1,NSET
!
        IF(    NEWSET(ISET)%SPECIAL.EQ.'ALL'    &
     &     .OR.NEWSET(ISET)%SPECIAL.EQ.'EMPTY') THEN
!         == CONSTRUCT SUMM OF ALL ANGULAR MOMENTUM WEIGHTS ====================
!         == THE RESULT FOR EMPTY WILL BE ADJUSTED BELOW =======================
          DO IKPT=1,NKPT
            DO ISPIN=1,NSPIN
              STATE=>STATEARR(IKPT,ISPIN)
              NB=STATE%NB
              CALL XXX(NDIM,NPRO,NB,ISPIN,NSPIN,STATE%VEC &
     &                ,NBB,MATEL(:,:,IKPT,ISET))
            ENDDO
          ENDDO
        END IF

        IF(NEWSET(ISET)%SPECIAL.EQ.'EMPTY') THEN
          MATEL(:,:,:,ISET)=-MATEL(:,:,:,ISET)
!          MATEL(1,:,:,ISET)=MATEL(1,:,:,ISET)+1.D0
        END IF

        IF(    NEWSET(ISET)%SPECIAL.EQ.'TOTAL' &
     &     .OR.NEWSET(ISET)%SPECIAL.EQ.'EMPTY') THEN
!         == MATEL IS THE PREFACTOR OF THE UNIT MATRIX
!         == WHICH HAS TRACE TWO. THEREFORE THE FACTOR 0.5
          MATEL(1,:,:,ISET)=MATEL(1,:,:,ISET)+0.5D0
!         == ADD SPIN INFORMATION FOR COLLINEAR SPIN-POLARIZED CALC
!         == NO SPIN INFORMATION IS AVAILABLE FOR TOTAL IN A NON-COLLINEAR CALC.
          IF(NSPIN.EQ.2.AND.NDIM.EQ.1) THEN
            DO IKPT=1,NKPT
              DO ISPIN=1,NSPIN
                STATE=>STATEARR(IKPT,ISPIN)
                NB=STATE%NB
                IF(ISPIN.EQ.1) THEN
                  MATEL(4,1:NB ,IKPT,ISET)=MATEL(4,1:NB     ,IKPT,ISET)+0.5D0
                ELSE IF(ISPIN.EQ.2.OR.NSPIN.EQ.1) THEN
                  MATEL(4,NB+1:,IKPT,ISET)=MATEL(4,NB+1:2*NB,IKPT,ISET)-0.5D0
                END IF
              ENDDO
            ENDDO
          END IF
        END IF
      ENDDO         
!
!     ==========================================================================
!     == PROCESS ANGULAR MOMENTUM WEIGHTS                                     ==
!     ==========================================================================
      DO ISET=1,NSET
        NWGHT=NEWSET(ISET)%NWGHT
        DO IWGHT=1,NWGHT
          IATP=NEWSET(ISET)%IAT(IWGHT)
          LP=NEWSET(ISET)%L(IWGHT)
          NP=-1
          IPRO=0
          DO IAT=1,NAT
            ISP=ISPECIES(IAT)
            DO LN=1,LNX(ISP)
              L=LOX(LN,ISP)
              IF(IAT.EQ.IATP.OR.IATP.EQ.-1) THEN
                IF(L.EQ.LP.OR.LP.EQ.-1) THEN
! N/NP THIS IS FOR LATER WHEN ALSO DIFFERENT ORBITALS PER L SHALL BE ADDRESSED
!                  IF(N.EQ.NP.OR.NP.EQ.-1) THEN 
                    DO IKPT=1,NKPT
                      DO ISPIN=1,NSPIN
                        STATE=>STATEARR(IKPT,ISPIN)
                        NB=STATE%NB
                        CALL XXX(NDIM,2*L+1,NB,ISPIN,NSPIN &
        &                                   ,STATE%VEC(:,IPRO+1:IPRO+2*L+1,:) &
        &                                   ,NBB,MATEL(:,:,IKPT,ISET))
                      ENDDO ! ISPIN
                    ENDDO !IKPT
!                 END IF !N
                END IF ! L
              END IF ! IAT
              IPRO=IPRO+2*L+1
            ENDDO ! LN
          ENDDO ! IAT
        ENDDO ! IWGHT
      ENDDO ! ISET
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE XXX(NDIM,NPRO,NB,ISPIN,NSPIN,C,NBB,MATEL)
!     **************************************************************************
!     ** CALCULATES GENERALIZED ANGULAR MOMENTUM WEIGHTS.                     **
!     ** MATEL ADDS THE SUM OVER THE WEIGHTS FROM THE SET OF ORBITALS SUPPLIED**
!     ** THE RESULT IS RESOLVED INTO THE (TOTAL,SX,SY,SZ) SPIN COMPONENTS     **
!     ** THE CALLING ROUTINE SPECIFIES THE ORBITAL RANGE.                     **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIM
      INTEGER(4),INTENT(IN) :: NPRO
      INTEGER(4),INTENT(IN) :: NSPIN
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: ISPIN
      COMPLEX(8),INTENT(IN) :: C(NDIM,NPRO,NB)
      INTEGER(4),INTENT(IN)    :: NBB
      REAL(8)   ,INTENT(INOUT) :: MATEL(4,NBB)
      INTEGER(4)               :: I
!     **************************************************************************
      DO I=1,NPRO
        IF(NDIM.EQ.2) THEN
          MATEL(1,:NB)=MATEL(1,:NB)+0.5D0* REAL(C(1,I,:)*CONJG(C(1,I,:)) &
     &                                         +C(2,I,:)*CONJG(C(2,I,:)))
          MATEL(2,:NB)=MATEL(2,:NB)+0.5D0* REAL(C(1,I,:)*CONJG(C(2,I,:)) &
     &                                         +C(2,I,:)*CONJG(C(1,I,:)))
          MATEL(3,:NB)=MATEL(3,:NB)+0.5D0*AIMAG(C(1,I,:)*CONJG(C(2,I,:)) &
     &                                         -C(2,I,:)*CONJG(C(1,I,:)))
          MATEL(4,:NB)=MATEL(4,:NB)+0.5D0* REAL(C(1,I,:)*CONJG(C(1,I,:)) &
     &                                         -C(2,I,:)*CONJG(C(2,I,:)))
        ELSE IF(NDIM.EQ.1) THEN
          IF(ISPIN.EQ.1) THEN
            MATEL(1,:NB)=MATEL(1,:NB)+0.5D0*REAL(C(1,I,:)*CONJG(C(1,I,:)))
            MATEL(4,:NB)=MATEL(4,:NB)+0.5D0*REAL(C(1,I,:)*CONJG(C(1,I,:)))
          END IF
          IF(ISPIN.EQ.2.OR.NSPIN.EQ.1) THEN
            MATEL(1,NB+1:2*NB)=MATEL(1,NB+1:2*NB) &
     &                        +0.5D0*REAL(C(1,I,:)*CONJG(C(1,I,:)))
            MATEL(4,NB+1:2*NB)=MATEL(4,NB+1:2*NB) &
     &                        -0.5D0*REAL(C(1,I,:)*CONJG(C(1,I,:)))
          END IF
        END IF
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEIPOINT(NX,LX,NAT,IPOINT)
!     **************************************************************************
!     ** DETERMINES THE IPRO INDEX FOR THE FIRST ORBITAL WITH GIVEN           **
!     ** ANGULAR MOMENTUM AND "MAIN QUANTUM NUMBER" N                         **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NAT1=>NAT &
     &                          ,ISPECIES &  !(NAT)
     &                          ,LNX &       !(LNX)
     &                          ,LOX         !(LNXX,ISP)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LX
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(OUT):: IPOINT(NX,LX+1,NAT)
      INTEGER(4)            :: NSP    !#(SPECIES)
      INTEGER(4)            :: ISP,IAT,L,N,LN,IPRO
!     **************************************************************************
!
!     ==========================================================================
!     ==  CHECK ARRAY SIZE                                                    ==
!     ==========================================================================
      IF(NAT.NE.NAT1) THEN
        CALL ERROR$MSG('NAT DIFFERS BETWEEN ARGUMENT AND MODULE')
        CALL ERROR$STOP('MAKEIPOINT')
      END IF
      NSP=MAXVAL(ISPECIES)
      DO ISP=1,NSP
        IF(MAXVAL(LOX(:LNX(ISP),ISP)).GT.LX) THEN
          CALL ERROR$MSG('LX TOO SMALL')
          CALL ERROR$STOP('MAKEIPOINT')
        END IF
        DO L=0,LX
          N=0
          DO LN=1,LNX(ISP)
            IF(LOX(LN,ISP).EQ.L) N=N+1
          ENDDO
          IF(N.GT.NX) THEN
            CALL ERROR$MSG('NX TOO SMALL')
            CALL ERROR$STOP('MAKEIPOINT')
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ALLOCATE AND FILL POINT ARRAY IPOINT(N,L+1,NAT)                     ==
!     ==  IPOINT CONTAINS THE FIRST VALUE IPRO FOR A A GIVEN (N,L,IAT)        ==
!     ==========================================================================
      IPOINT(:,:,:)=0
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          N=1
          DO WHILE (IPOINT(N,L+1,IAT).NE.0)
            N=N+1
          ENDDO
          IPOINT(N,L+1,IAT)=IPRO+1
          IPRO=IPRO+2*L+1
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READORBENTRY(LL_CNTL,RBAS,NAT,ATOMID,RPOS,IAT,IT,LMXX,ORB)
!     **************************************************************************
!     ** READ ONE !ORB BLOCK IN THE !DCNTL FILE WITH ATOM=' ' SPECIFIED       **
!     ** AND RETURN A THE DATA FOR AN ORBITAL ENTRY                           **
!     **                                                                      **
!     ** 1) ALL ATOMS CAN BE SPECIFIED IN EXTENDED ATOM NOTATION              **
!     ** 2) A SPECIFIED TRANSLATION VECTOR IT SHIFTS THE FINAL ORBITAL        **
!     **    BY AN ADDITIONAL LATTICE TRANSLATION                              **
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL
      REAL(8)      ,INTENT(IN) :: RBAS(3,3)   ! LATTICE VECTORS 
      INTEGER(4)   ,INTENT(IN) :: NAT         ! #(ATOMS)
      CHARACTER(16),INTENT(IN) :: ATOMID(NAT) ! ATOM NAMES
      REAL(8)      ,INTENT(IN) :: RPOS(3,NAT) ! ATOM POSITIONS
      INTEGER(4)   ,INTENT(OUT):: IAT         ! ATOM INDEX
      INTEGER(4)   ,INTENT(OUT):: IT(3)       ! UNIT-CELL INDEX
      INTEGER(4)   ,INTENT(IN) :: LMXX        ! MAX LENGTH OF ORBITAL ARRAY
      COMPLEX(8)   ,INTENT(OUT):: ORB(LMXX)   ! ORBITAL COEFFICIENT VECTOR
      CHARACTER(8)             :: TYPE  ! ORBITAL TYPE IN LOCAL COORDINATES
      INTEGER(4)               :: IAT2
      CHARACTER(32)            :: ATOMEX      ! ATOM NAME IN EXTENDED NOTATION
      CHARACTER(32)            :: ATOM
      CHARACTER(32)            :: ATOM2
      INTEGER(4)               :: IT2(3)
      INTEGER(4)               :: I
      COMPLEX(8)               :: CFAC
      LOGICAL(4)               :: TCHK,TCHK1
      REAL(8)                  :: DRZ(3)
      REAL(8)                  :: DRX(3)
      REAL(8)                  :: ROT(3,3)
      REAL(8)                  :: YLMROT(LMXX,LMXX)
!     **************************************************************************
      IAT=0
      IT(:)=0
      ORB(:)=(0.D0,0.D0)
!
!     ==========================================================================
!     ==  CHECK DATA COMPATIBILITY                                            ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('INCOMPATIBLE DATA IDENTIFIER "NAME" ENCOUNTERED')
        CALL ERROR$STOP('READORBENTRY')
      END IF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK)
      IF(.NOT.TCHK) THEN   
        CALL ERROR$MSG('NO ATOM IDENTIFIER "ATOM" ENCOUNTERED')
        CALL ERROR$STOP('READORBENTRY')
      END IF
!
!     ==========================================================================
!     ==  SPECIFY ATOM SITE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,ATOMEX)
      CALL RESOLVEEXTENDEDNAME(ATOMEX,ATOM,IT)
      IAT=0
      DO I=1,NAT
        IF(ATOM.NE.ATOMID(I)) CYCLE
        IAT=I
        EXIT
      ENDDO
      IF(IAT.EQ.0) THEN
        CALL ERROR$MSG('ATOM NAME NOT FOUND')
        CALL ERROR$CHVAL('ATOM',ATOM)
        CALL ERROR$STOP('READORBENTRY')
      END IF
!     
!     ==========================================================================
!     ==  SPECIFY LOCAL Z-COORDINATE DRZ                                      ==
!     ==========================================================================
!     == DEFAULT: LOCAL COORDINATES EQUAL GLOBAL COORDINATE 
      DRZ(:)=(/0.D0,0.D0,1.D0/)
!
!     == SPECIFY LOCAL AXIS IN GLOBAL COORDINATES ==============================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'Z',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'Z',1,DRZ)
!
!     == SPECIFY LOCAL Z-COORDINATE BY NEIGBOR ATOM ============================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NNZ',1,TCHK)
      IF(TCHK) THEN
!       == CHECK INCOMPATIBILITY OF DATA SPECIFICATIONS ========================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'Z',1,TCHK1)
        IF(TCHK1) THEN
CALL LINKEDLIST$REPORT(LL_CNTL,6)
          CALL ERROR$MSG('DO NOT SPECIFY BOTH "Z" AND "NNZ"')
          CALL ERROR$STOP('READORBENTRY')
        END IF
!
!       == COLLECT NEIGHBOR VECTOR =============================================
        CALL LINKEDLIST$GET(LL_CNTL,'NNZ',1,ATOMEX)
        CALL RESOLVEEXTENDEDNAME(ATOMEX,ATOM,IT2)
        IAT2=0
        DO I=1,NAT
          IF(ATOM.NE.ATOMID(I)) CYCLE
          IAT2=I
          EXIT
        ENDDO
        IF(IAT2.EQ.0) THEN
          CALL ERROR$MSG('NAME OF ATOM DEFINING LOCAL Z-AXIS NOT FOUND')
          CALL ERROR$CHVAL('ATOM',ATOM)
          CALL ERROR$STOP('READORBENTRY')
        END IF
!
        DRZ(:)=RPOS(:,IAT2)-RPOS(:,IAT) &
     &        +RBAS(:,1)*REAL(IT2(1)-IT(1),KIND=8) &
     &        +RBAS(:,2)*REAL(IT2(2)-IT(2),KIND=8) &
     &        +RBAS(:,3)*REAL(IT2(3)-IT(3),KIND=8)
      END IF
!       
!     ==========================================================================
!     ==  SPECIFY LOCAL X-COORDINATE DRX                                      ==
!     ==========================================================================
!     == DEFAULT: LOCAL COORDINATES EQUAL GLOBAL COORDINATE 
      DRX(:)=(/1.D0,0.D0,0.D0/)
!
!     == SPECIFY LOCAL AXIS IN GLOBAL COORDINATES ==============================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'X',1,TCHK)
      IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'X',1,DRX(:))
!
!     == SPECIFY LOCAL Z-COORDINATE BY NEIGBOR ATOM ============================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NNX',1,TCHK)
      IF(TCHK) THEN
!       == CHECK INCOMPATIBILITY OF DATA SPECIFICATIONS ========================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'X',1,TCHK1)
        IF(TCHK1) THEN
          CALL ERROR$MSG('DO NOT SPECIFY BOTH "X" AND "NNX"')
          CALL ERROR$STOP('READORBENTRY')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'NNX',1,ATOMEX)
        CALL RESOLVEEXTENDEDNAME(ATOMEX,ATOM,IT2)
        IAT2=0
        DO I=1,NAT
          IF(ATOM.NE.ATOMID(I)) CYCLE
          IAT2=I
          EXIT
        ENDDO
        IF(IAT2.EQ.0) THEN
          CALL ERROR$MSG('NAME OF ATOM DEFINING LOCAL X-AXIS NOT FOUND')
          CALL ERROR$CHVAL('ATOM',ATOM)
          CALL ERROR$STOP('READORBENTRY')
        END IF
!
        DRX(:)=RPOS(:,IAT2)-RPOS(:,IAT) &
     &        +RBAS(:,1)*REAL(IT2(1)-IT(1),KIND=8) &
     &        +RBAS(:,2)*REAL(IT2(2)-IT(2),KIND=8) &
     &        +RBAS(:,3)*REAL(IT2(3)-IT(3),KIND=8)
      END IF
!       
!     ==========================================================================
!     ==  COLLECT ORBITAL TYPE                                                ==
!     ==========================================================================
      CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
      TYPE=+TYPE
      CALL RESOLVETYPE(LMXX,TYPE,ORB)
!       
!     ==========================================================================
!     ==  ROTATE ORBITAL INTO LOCAL COORDINATE AXIS                           ==
!     ==========================================================================
      CALL RESOLVEROTATION(DRZ,DRX,ROT)
      CALL ROTATEYLM(LMXX,ROT,YLMROT)
      ORB=MATMUL(CMPLX(YLMROT,KIND=8),ORB)
!
!     ==========================================================================
!     ==  MULTYPLY WITH OPTIONAL PREFACTOR                                    ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FAC',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$CONVERT(LL_CNTL,'FAC',1,'C(8)')
        CALL LINKEDLIST$GET(LL_CNTL,'FAC',1,CFAC)
        ORB=ORB*CFAC
      END IF
!
!     ==========================================================================
!     ==  SHIFT BY OPTIONAL LATTICE TRANSLATION                               ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'IT',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'IT',1,IT2)
        IT=IT+IT2
      END IF
      RETURN
      END SUBROUTINE READORBENTRY
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RESOLVETYPE(LMX,TYPE,ORBITAL)
!     **************************************************************************
!     ** CONSTRUCTS THE PREFACTORS OF AN ORBITAL IN AN EXPANSION OF           **
!     ** REAL SPHERICAL HARMONICS                                             **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: LMX
      CHARACTER(*),INTENT(IN)   :: TYPE
      COMPLEX(8)  ,INTENT(OUT)  :: ORBITAL(LMX)
      REAL(8)                   :: ORB(16)
      INTEGER(4)                :: LM
!     **************************************************************************
      ORB(:)=0.D0
!
!     ==========================================================================
!     ==  SET ORBITAL COEFFICIENTS                                            ==
!     ==========================================================================
      IF(TRIM(TYPE).EQ.'SP1') THEN
!       == SP1 HYBRID ORBITAL IN Z DIRECTION ===================================
        ORB(1)=SQRT(1.D0/2.D0)  ! S
        ORB(3)=SQRT(1.D0/2.D0)  ! PZ
!
      ELSE IF(TRIM(TYPE).EQ.'SP2') THEN
!       == SP2 HYBRID ORBITAL IN Z DIRECTION ===================================
        ORB(1)=SQRT(1.D0/3.D0)
        ORB(3)=SQRT(2.D0/3.D0)
!
      ELSE IF(TRIM(TYPE).EQ.'SP3') THEN
!       == SP3 HYBRID ORBITAL IN Z DIRECTION ===================================
        ORB(1)=SQRT(1.D0/4.D0)
        ORB(3)=SQRT(3.D0/4.D0)
!
      ELSE
!       == RESOLVE PURE REAL SPHERICAL HARMONICS ===============================
        CALL SPHERICAL$LMBYNAME(TYPE,LM)
        ORB(LM)=1.D0
      END IF
!
!     ==========================================================================
!     == CHECK WHETHER THE ORBITAL VECTOR FITS INTO THE ARRAY SUPPLIED        ==
!     ==========================================================================
      DO LM=LMX+1,16
        IF(ORB(LM).NE.0) THEN
          CALL ERROR$MSG('LMX TOO SMALL')
          CALL ERROR$CHVAL('TYPE',TRIM(TYPE))
          CALL ERROR$I4VAL('LMX',LMX)
          CALL ERROR$STOP('RESOLVETYPE')
        END IF
      ENDDO
!
!     ==========================================================================
!     == MAP ORB VECTOR ONTO OUTPUT                                           ==
!     ==========================================================================
      IF(LMX.LE.16) THEN
        ORBITAL(:)=CMPLX(ORB(:LMX),KIND=8)
      ELSE
        ORBITAL(1:16)=CMPLX(ORB,KIND=8)
        ORBITAL(17:)=(0.D0,0.D0)
      END IF
      RETURN
      END SUBROUTINE RESOLVETYPE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEORBITAL(ATOM,LMXX,ORB,NPRO_,ORBITAL)
!     **************************************************************************
!     **  ONLY THE FIRST PARTIAL WAVE PER ANGULAR MOMENTUM IS                 **
!     **  CONSIDERED                                                          **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE, ONLY : NAT &
     &                       ,ISPECIES &
     &                       ,ATOMID &
     &                       ,LNX &
     &                       ,LOX 
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ATOM
      INTEGER(4)  ,INTENT(IN) :: LMXX
      COMPLEX(8)  ,INTENT(IN) :: ORB(LMXX)
      INTEGER(4)  ,INTENT(IN) :: NPRO_
      COMPLEX(8)  ,INTENT(OUT):: ORBITAL(NPRO_)
      INTEGER(4)              :: IPRO,IAT,ISP,LN,L,LM,M
      LOGICAL(4)              :: TCHK
      LOGICAL(4)              :: LCHK(10)
!     **************************************************************************
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        TCHK=(ATOM.EQ.ATOMID(IAT))
        LCHK=.TRUE.
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(TCHK.AND.LCHK(L+1)) THEN
            LCHK(L+1)=.FALSE.
!PB:WARNING: PICKS ONLY THE WEIGHT FROM THE FIRST PARTIAL WAVE
            LM=L**2
            DO M=1,2*L+1
              IPRO=IPRO+1
              LM=LM+1
              ORBITAL(IPRO)=ORB(LM)
            ENDDO
          ELSE
            IPRO=IPRO+2*L+1
          END IF
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MAKEORBITAL
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RESOLVEROTATION(DZ_,DX_,ROT)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: DZ_(3)
      REAL(8)   ,INTENT(IN)  :: DX_(3)
      REAL(8)   ,INTENT(OUT) :: ROT(3,3)
      REAL(8)                :: DZ(3)
      REAL(8)                :: DX(3)
      REAL(8)                :: DY(3)
      REAL(8)                :: DZLEN,DXLEN
!     **************************************************************************
      DX=DX_
      DZ=DZ_
!     
!     ==========================================================================
!     == NORMALIZE AND COMPLETE VECTORS                                       ==
!     ==========================================================================
!     == SET DZ ================================================================
      DZLEN=SQRT(DZ(1)**2+DZ(2)**2+DZ(3)**2)
      IF(DZLEN.EQ.0.D0) THEN
        DZ=(/0.D0,0.D0,1.D0/)
      ELSE
        DZ=DZ/DZLEN
      END IF
!     == SET DX ================================================================
      DXLEN=SQRT(DX(1)**2+DX(2)**2+DX(3)**2)
      IF(DXLEN.EQ.0.D0) THEN
        DX=(/1.D0,0.D0,0.D0/)
      ELSE
        DX=DX/DXLEN
      END IF
      DX=DX-DZ*DOT_PRODUCT(DZ,DX)
      DXLEN=SQRT(DX(1)**2+DX(2)**2+DX(3)**2)
      IF(DXLEN.EQ.0.D0) THEN   ! DX=0 OR PARALLEL TO DZ
        DX=(/0.D0,1.D0,0.D0/)  ! CHOOSE Y DIRECTION AS ALTERNATIVE
        DX=DX-DZ*DOT_PRODUCT(DZ,DX) !JO AB DA
        DXLEN=SQRT(DX(1)**2+DX(2)**2+DX(3)**2) 
        IF(DXLEN.EQ.0.D0) THEN
          DX=(/0.D0,0.D0,1.D0/) ! CHOOSE Z DIRECTION AS ALTERNATIVE
          DX=DX-DZ*DOT_PRODUCT(DZ,DX) 
          DXLEN=SQRT(DX(1)**2+DX(2)**2+DX(3)**2)   
        END IF       !JO BIS DA
      END IF
      DX=DX/DXLEN
!     == SET DY ================================================================
      DY(1)=DZ(2)*DX(3)-DZ(3)*DX(2)
      DY(2)=DZ(3)*DX(1)-DZ(1)*DX(3)
      DY(3)=DZ(1)*DX(2)-DZ(2)*DX(1)
!     
!     ==========================================================================
!     == NORMALIZE AND COMPLETE VECTORS                                       ==
!     ==========================================================================
      ROT(:,3)=DZ(:)
      ROT(:,2)=DY(:)
      ROT(:,1)=DX(:)
!     WRITE(*,FMT='("ROT",3F10.5)')ROT(1,:)
!     WRITE(*,FMT='("ROT",3F10.5)')ROT(2,:)
!     WRITE(*,FMT='("ROT",3F10.5/)')ROT(3,:)
      RETURN
      END SUBROUTINE RESOLVEROTATION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NORMALIZEORBITAL(NPRO_,ORBITAL)
!     **************************************************************************
!     **************************************************************************
!     USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)    :: NPRO_
      COMPLEX(8)  ,INTENT(INOUT) :: ORBITAL(NPRO_)
      INTEGER(4)                 :: I
      REAL(8)                    :: SUM
!     **************************************************************************
      SUM=0.D0
      DO I=1,NPRO_
        SUM=SUM+REAL(CONJG(ORBITAL(I))*ORBITAL(I))
      ENDDO
      IF(SUM.LE.0.D0) THEN
        CALL ERROR$MSG('ORBITAL COEFFICIENTS VANISH')
        CALL ERROR$STOP('NORMALIZEORBITAL')
      END IF
      SUM=1.D0/SQRT(SUM)
      DO I=1,NPRO_
        ORBITAL(I)=ORBITAL(I)*SUM
      ENDDO
      RETURN
      END SUBROUTINE NORMALIZEORBITAL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RESOLVEATOM(ATOMEX,IAT,IT)
!     **************************************************************************
!     **  RESOLVES THE EXTENDED ATOM NAME NOTATION, WHICH INCLUDES            **
!     **  A LATTICE TRANSLATION                                               **
!     **                                                                      **
!     **  THE EXTENDED NOTATION INCLUDES AN INTEGER LATTICE TRANSLATIONS      **
!     **  IN THE ATOM NAME FOLLOWING A COLON                                  **
!     **                                                                      **
!     **   'O_23:+1-1+1'  ATOM 'O_23' SHIFTED BY RBAS(:,1)-RBAS(:,2)+RBAS(:,3)**
!     **                                                                      **
!     **   THE '+'SIGNS ARE NOT REQUIRED.                                     **
!     **   ONLY SINGLE-DIGIT TRANSLATIONS ARE PERMITTED                       **
!     **                                                                      **
!     ** ANALOGOUS TO STRCIN_RESOLVEEXTENDEDNAME(XNAME,NAME,IT)
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE, ONLY : NAT,ATOMID
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ATOMEX
      INTEGER(4)  ,INTENT(OUT) :: IAT
      INTEGER(4)  ,INTENT(OUT) :: IT(3) !INTEGER TRANSLATION
      INTEGER(4)               :: I
      CHARACTER(32)            :: ATOM
!     ******************************************************************
      CALL RESOLVEEXTENDEDNAME(ATOMEX,ATOM,IT)
!
      DO I=1,NAT
        IF(ATOM.EQ.ATOMID(I)) THEN
          IAT=I
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ATOM NAME NOT FOUND')
      CALL ERROR$CHVAL('ATOM',ATOM)
      CALL ERROR$STOP('RESOLVEATOM')
      RETURN
      END SUBROUTINE RESOLVEATOM
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RESOLVEEXTENDEDNAME(XNAME,NAME,IT)
!     **************************************************************************
!     **  RESOLVES THE EXTENDED ATOM NAME NOTATION, WHICH INCLUDES            **
!     **  A LATTICE TRANSLATION                                               **
!     **                                                                      **
!     **  THE EXTENDED NOTATION INCLUDES AN INTEGER LATTICE TRANSLATIONS      **
!     **  IN THE ATOM NAME FOLLOWING A COLON                                  **
!     **                                                                      **
!     **   'O_23:+1-1+1'  ATOM 'O_23' SHIFTED BY RBAS(:,1)-RBAS(:,2)+RBAS(:,3)**
!     **                                                                      **
!     **   THE '+'SIGNS ARE NOT REQUIRED.                                     **
!     **   ONLY SINGLE-DIGIT TRANSLATIONS ARE PERMITTED                       **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: XNAME  ! EXTENDED ATOM NAME
      CHARACTER(*),INTENT(OUT):: NAME   ! NON-EXTENDED ATOM NAME
      INTEGER(4)  ,INTENT(OUT):: IT(3)  ! INTEGER LATTICE TRANSLATIONS
      INTEGER(4)              :: ICOLON ! POSITION OF THE COLON IN XNAME
      INTEGER(4)              :: IPOS,IND,SGN
      INTEGER(4)              :: ICH    ! ASCII NUMBER OF THE SELECTED LETTER
!     **************************************************************************
      ICOLON=INDEX(XNAME,':')
!     == RETURN IF NO TRANSLATION VECTOR GIVEN =================================
      IF(ICOLON.EQ.0) THEN
        NAME=XNAME
        IT(:)=0
        RETURN
      END IF
!
!     ==========================================================================
!     == RESOLVE EXTENDED ATOM NAME                                           ==
!     ==========================================================================
      NAME=XNAME(:ICOLON-1)
      IPOS=ICOLON+1
      IND=0
      SGN=+1
      DO WHILE(IND.LT.3) 
        ICH=IACHAR(XNAME(IPOS:IPOS))
!       ==  IACHAR('+')=43; IACHAR('-')=45; IACHAR('0')=48; IACHAR('1')=49;...
        IF(ICH.GE.48.AND.ICH.LE.57) THEN ! IF "0,1,...,9"
          IND=IND+1
          IT(IND)=SGN*(ICH-48)
          SGN=+1
        ELSE IF(ICH.EQ.43) THEN   ! IF "+"
          SGN=+1
        ELSE IF(ICH.EQ.45) THEN   ! IF "-"
          SGN=-1
        ELSE
          CALL ERROR$MSG('ILLEGAL CHARACTER IN EXTENDED ATOM NOTATION')  
          CALL ERROR$CHVAL('EXT. NAME ',XNAME)
          CALL ERROR$CHVAL('ILLEGAL CHARACTER ',XNAME(IPOS:IPOS))
          CALL ERROR$STOP('STRCIN_RESOLVEEXTENDEDNAME')
        END IF
        IPOS=IPOS+1
      ENDDO
      IF(XNAME(IPOS:).NE.' ') THEN
        CALL ERROR$MSG('LETTERS FOUND BEYOND END OF EXTENDED ATOM NOTATION')  
        CALL ERROR$CHVAL('EXT. NAME ',XNAME)
        CALL ERROR$CHVAL('ADDITIONAL LETTERS ',XNAME(IPOS:))
        CALL ERROR$STOP('STRCIN_RESOLVEEXTENDEDNAME')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE, ONLY : LL_CNTL
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
      LOGICAL(4)           :: TCHK
      INTEGER(4)           :: NFIL
      CHARACTER(32)        :: ID 
      CHARACTER(256)       :: FILENAME
      INTEGER(4)           :: ITH
      INTEGER(4)           :: NUM
      INTEGER(4)           :: NFILO
!     **************************************************************************
                          CALL TRACE$PUSH('READCNTL')
!
!     ==========================================================================
!     ==  READ CONTROL FILE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('DCNTL',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)

      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO) 
        CALL LINKEDLIST$REPORT(LL_CNTL,NFILO)
      END IF
!
!     ==========================================================================
!     ==  !PDOSIN!FILES!FILE                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILES',1,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',ITH)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
!       ==  READ ACTUAL VALUES  ======================================
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILENAME)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        CALL FILEHANDLER$SETFILE(ID,TCHK,FILENAME)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$GENERIC(MODE,PREFIX)
!     **************************************************************************
!     ** READ !DCNTL!GENERIC BLOCK FROM THE CONTROL FILE                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE, ONLY: LL_CNTL
      IMPLICIT NONE
      CHARACTER(*),INTENT(OUT) :: MODE    ! "SAMPLE" OR "TETRA"
      CHARACTER(*),INTENT(OUT) :: PREFIX  ! PREFIX FOR DOS AND NOS FILES
      LOGICAL(4)               :: TCHK
!     **************************************************************************
!     ==========================================================================
!     == SET DEFAULT VALUES                                                   ==
!     ==========================================================================
      MODE='TETRA'
      PREFIX=' '
!
!     ==========================================================================
!     == READ GENERIC BLOCK                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'GENERIC',1,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MODE',1,TCHK)
      IF(TCHK)THEN
        CALL LINKEDLIST$GET(LL_CNTL,'MODE',1,MODE)
        IF(MODE.NE.'SAMPLE'.AND.MODE.NE.'TETRA')THEN
          CALL ERROR$MSG('MODE UNKNOWN (SHOULD BE "SAMPLE" OR "TETRA")')
          CALL ERROR$CHVAL('MODE ',MODE)
          CALL ERROR$STOP('READCNTL$GENERIC')
        ENDIF
        IF(MODE.EQ.'SAMPLE')THEN
          CALL ERROR$MSG('MODE SAMPLE NOT PROPERLY IMPLEMENTED')
          CALL ERROR$MSG('STOPPING...')
          CALL ERROR$CHVAL('MODE ',MODE)
          CALL ERROR$STOP('READCNTL$GENERIC')
        ENDIF
      ENDIF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'PREFIX',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'PREFIX',1,PREFIX)
!
!     ==========================================================================
!     == OBSOLETE INPUT DATA                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DOS',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('!DCNTL!GENERIC:DOS IS OBSOLETE')
        CALL ERROR$MSG('THE DENSITY OF STATES ID CALCULATED AS DEFAULT')
        CALL ERROR$STOP('READCNTL$GENERIC')
      END IF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NOS',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('!DCNTL!GENERIC:NOS IS OBSOLETE')
        CALL ERROR$MSG('THE DENSITY OF STATES ID CALCULATED AS DEFAULT')
        CALL ERROR$MSG('THE NUMBER OF STATES IS OBTAINED BY INTEGRATION')
        CALL ERROR$MSG('THE ENERGY GRID COVERS ALL CORE STATES')
        CALL ERROR$STOP('READCNTL$GENERIC')
      END IF
      RETURN
      END SUBROUTINE READCNTL$GENERIC
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$GRID(EMIN,EMAX,NE,EBROAD)
!     **************************************************************************
!     ** READ !DCNTL!GRID FROM CONTROL FILE                                   **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NSPIN &
     &                          ,NKPT &
     &                          ,STATEARR,STATE 
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE, ONLY: LL_CNTL
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT)   :: EMIN   ! MINIMUM OF ENERGY GRID
      REAL(8)   ,INTENT(OUT)   :: EMAX   ! MAXIMMUM OF ENERGY GRID
      INTEGER(4),INTENT(OUT)   :: NE     ! NUMBER OF ENERGY GRID POINTS
      REAL(8)   ,INTENT(OUT)   :: EBROAD ! THERMAL ENERGY BROADENING
      REAL(8)                  :: EV     ! ELECTRON VOLT
      REAL(8)                  :: KB     ! BOLTZMANN CONSTANT
      REAL(8)                  :: DE
      LOGICAL(4)               :: TCHK,TCHK1
      INTEGER(4)               :: ISPIN,IKPT
      INTEGER(4)               :: NB
!     **************************************************************************
      EMIN=HUGE(EMIN)
      EMAX=HUGE(EMAX)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
          NB=STATE%NB
          EMIN=MIN(EMIN,MINVAL(STATE%EIG))
          EMAX=MIN(EMAX,STATE%EIG(NB))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == SET DEFAULT VALUES                                                   ==
!     == (DEFAULT FOR EMIN/EMAX IS INPUT)                                     ==
!     ==========================================================================
      CALL CONSTANTS('EV',EV)
      CALL CONSTANTS('KB',KB)
      EBROAD=KB*300.D0
      DE=1.D-2*EV
      EMIN=EMIN-20.D0*EBROAD
      EMAX=EMAX+20.D0*EBROAD
      NE=INT((EMAX-EMIN)/DE)+1
!
!     ==========================================================================
!     == READ GRID BLOCK                                                      ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GRID')
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DE[EV]',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'DE[EV]',1,DE)
        DE=DE*EV
        NE=INT((EMAX-EMIN)/DE)+1
      END IF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'BROADENING[EV]',1,TCHK)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'BROADENING[K]',1,TCHK1)
      IF(TCHK.AND.TCHK1) THEN
        CALL ERROR$MSG('SELECT ONLY ONE OF THE ALTERNATIVE OPTIONS, ')
        CALL ERROR$MSG('EITHER !DCNTL!GRID:BROADENING[EV]')
        CALL ERROR$MSG('    OR !DCNTL!GRID:BROADENING[K]')
        CALL ERROR$STOP('READCNTL$GRID')
      END IF
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'BROADENING[EV]',1,EBROAD)
        EBROAD=EBROAD*EV
      END IF
      IF(TCHK1) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'BROADENING[K]',1,EBROAD)
        EBROAD=EBROAD*KB
      END IF
!
!     ==========================================================================
!     == OBSOLETE INPUT DATA                                                  ==
!     ==========================================================================
!     ==  READ ACTUAL VALUES  ==================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMIN[EV]',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('!DCNTL!GRID:EMIN[EV] IS OBSOLETE')
        CALL ERROR$MSG('THE ENTIRE ENERGY REGION IS USED ')
        CALL ERROR$STOP('READCNTL$GRID')
      END IF
!!$      IF(TCHK) THEN
!!$        CALL LINKEDLIST$GET(LL_CNTL,'EMIN[EV]',1,EMIN)
!!$        EMIN=EMIN*EV
!!$      END IF
!!$!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMAX[EV]',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('!DCNTL!GRID:EMAX[EV] IS OBSOLETE')
        CALL ERROR$MSG('THE ENTIRE ENERGY REGION IS USED ')
        CALL ERROR$STOP('READCNTL$GRID')
      END IF
!!$      IF(TCHK) THEN
!!$        CALL LINKEDLIST$GET(LL_CNTL,'EMAX[EV]',1,EMAX)
!!$        EMAX=EMAX*EV
!!$      END IF

      CALL LINKEDLIST$EXISTD(LL_CNTL,'SCALEY',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('!DCNTL!GRID:SCALEY IS OBSOLETE')
        CALL ERROR$MSG('DOS AND NOS FILES ARE PRODUCED AS SEPARATE FILES')
        CALL ERROR$MSG('DEPENDING ON THE SETTING OF "DOS" AND "NOS"')
        CALL ERROR$STOP('READCNTL$GRID')
      END IF
!
      RETURN
      END SUBROUTINE READCNTL$GRID
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$REPORT1(MODE,PREFIX,EMIN,EMAX,NE,EBROAD)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MODE    ! "SAMPLE" OR "TETRA"
      CHARACTER(*),INTENT(IN) :: PREFIX  ! PREFIX FOR DOS AND NOS FILES
      REAL(8)     ,INTENT(IN) :: EMIN   ! MINIMUM OF ENERGY GRID
      REAL(8)     ,INTENT(IN) :: EMAX   ! MAXIMMUM OF ENERGY GRID
      INTEGER(4)  ,INTENT(IN) :: NE     ! NUMBER OF ENERGY GRID POINTS
      REAL(8)     ,INTENT(IN) :: EBROAD ! THERMAL ENERGY BROADENING
      INTEGER(4)              :: NFILO
      REAL(8)                 :: EV
!     **************************************************************************
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL CONSTANTS('EV',EV)
      CALL REPORT$CHVAL(NFILO,'INTEGRATION METHOD',MODE)
      CALL REPORT$CHVAL(NFILO,'PREFIX FOR DOS AND NOS FILES',PREFIX)
      CALL REPORT$R8VAL(NFILO,'ENERGY GRID STARTS AT',EMIN/EV,'EV')
      CALL REPORT$R8VAL(NFILO,'ENERGY GRID ENDS AT',EMAX/EV,'EV')
      CALL REPORT$R8VAL(NFILO,'SPACING OF THE ENERGY GRID' &
     &                       ,(EMAX-EMIN)/REAL(NE-1)/EV,'EV')
      WRITE(NFILO,*)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$ORBITAL(NPRO,NAT,RBAS,RPOS,ATOMID)
!     **************************************************************************
!     ** READ AND CONSTRUCT PREDEFINED ORBITALS                               **
!     ** 
!     ** 
!     ** 
!     ** REMARK: IN THE NEW VERSION THE ORBITAL IS NO MORE NORMALIZED.        **
!     ** 
!     **************************************************************************
      USE LINKEDLIST_MODULE 
      USE READCNTL_MODULE   ,ONLY : LL_CNTL
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NPRO
      INTEGER(4)   ,INTENT(IN) :: NAT
      REAL(8)      ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)      ,INTENT(IN) :: RPOS(3,NAT)  !ATOMIC POSITIONS
      CHARACTER(16),INTENT(IN) :: ATOMID(NAT)
      INTEGER(4)   ,PARAMETER  :: LMXX=25 
      CHARACTER(32)            :: ORBITALNAME
      CHARACTER(32)            :: ORBITALNAME1
      COMPLEX(8)               :: ORBITAL(NPRO)
      COMPLEX(8)               :: ORBITAL1(NPRO)
      INTEGER(4)               :: IORB,ITH,IORB2
      INTEGER(4)               :: IAT,IT(3)
      COMPLEX(8)               :: ORB(LMXX)
      INTEGER(4)               :: NUM
      INTEGER(4)               :: NORB
      COMPLEX(8)               :: CFAC
      LOGICAL(4)               :: TCHK
!     **************************************************************************
                          CALL TRACE$PUSH('READCNTL$ORBITAL')
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'ORBITAL',NORB)
      DO IORB=1,NORB
        CALL LINKEDLIST$SELECT(LL_CNTL,'ORBITAL',IORB)
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,ORBITALNAME)
        CALL NEWORBITAL$NEWORBITAL(ORBITALNAME)
        CALL NEWORBITAL$SELECT(ORBITALNAME)
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB',NUM)
        DO ITH=1,NUM
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB',ITH)
          CALL READONENEWORB(LL_CNTL,RBAS,NAT,ATOMID,RPOS)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL NEWORBITAL$SELECT('NONE')
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL$ORBITAL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READONENEWORB(LL_CNTL,RBAS,NAT,ATOMID,RPOS)
!     **************************************************************************
!     **  READS CONTENTS OF BLOCK !ORB FROM LINKEDLIST AND ADDS IT TO THE     **
!     **  CURRENTLY SELECTED ORBITAL OF THE NEWORBITAL OBJECT                 **
!     **                                                                      **
!     **  THE !ORB BLOCK MAY SELECT A PREDEFINED ORBITAL WITH KEYWORD "NAME=" **
!     **  OR IT SELECTS AN NEW ENTRY WITH THE KEYWORD "ATOM="                 **
!     **                                                                      **
!     **  PRECONDITIONS:                                                      **
!     **   - THE LINKED LIST LL_CNTL MUST BE POSITIONED INSIDE THE !ORB BLOCK **
!     **   - ONE ORBITAL IS SELECTED IN THE NEWORBITAL OBJECT                 **
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE, ONLY : LL_TYPE &
     &                             ,LINKEDLIST$CONVERT &
     &                             ,LINKEDLIST$EXISTD &
     &                             ,LINKEDLIST$GET
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL      ! LINKED LIST WITH INPUT DATA
      REAL(8)      ,INTENT(IN) :: RBAS(3,3)    ! LATTICE VECTORS
      INTEGER(4)   ,INTENT(IN) :: NAT          ! #(ATOMS)
      CHARACTER(16),INTENT(IN) :: ATOMID(NAT)  ! ATOM NAMES
      REAL(8)      ,INTENT(IN) :: RPOS(3,NAT)  ! ATOMIC POSITIONS
      INTEGER(4)   ,PARAMETER  :: LMXX=16
      LOGICAL(4)               :: TCHK,TCHK1
      CHARACTER(32)            :: ORBITALNAME1
      INTEGER(4)               :: IORB2
      INTEGER(4)               :: IAT
      INTEGER(4)               :: IT(3)
      COMPLEX(8)               :: CFAC
      COMPLEX(8)               :: ORB(LMXX)
!     **************************************************************************
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK1)
      IF(TCHK1.EQV.TCHK) THEN
        CALL ERROR$MSG('!ORB:NAME AND !ORB:ATOM ARE MUTUALLY EXCLUSIVE')
        CALL ERROR$MSG('ONE OF THE TWO NEEDS TO BE SPECIFIED')
        CALL ERROR$L4VAL('NAME PRESENT',TCHK)
        CALL ERROR$L4VAL('ATOM PRESENT',TCHK1)
        CALL ERROR$STOP('READONENEWORB')
      END IF
!
!     ==========================================================================
!     == RESOLVE ORBITAL (!ORB CONTAINS NAME= )                               ==
!     ==========================================================================
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,ORBITALNAME1)

!
!       ========================================================================
!       ==  GET PREFACTOR                                                     ==
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FAC',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$CONVERT(LL_CNTL,'FAC',1,'C(8)')
          CALL LINKEDLIST$GET(LL_CNTL,'FAC',1,CFAC)
        ELSE
          CFAC=(1.D0,0.D0)
        END IF
!
!       ========================================================================
!       ==  GET LATTICE TRANSLATION                                           ==
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'IT',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'IT',1,IT)
        ELSE
          IT=(/0,0,0/)
        END IF
!
!       ========================================================================
!       ==  ADD ORBITAL TO CURRENT ORBITAL                                    ==
!       ========================================================================
        CALL NEWORBITAL$IORB(ORBITALNAME1,IORB2)
        CALL NEWORBITAL$ADDORBITAL(IORB2,IT(:),CFAC)
!
!     ==========================================================================
!     ==  ADD NEW ENTRY (!ORB CONTAINS KEYWORD "ATOM=" )                      ==
!     ==========================================================================
      ELSE
        CALL READORBENTRY(LL_CNTL,RBAS,NAT,ATOMID,RPOS,IAT,IT,LMXX,ORB)
        CALL NEWORBITAL$ADDENTRY(IAT,IT,LMXX,ORB)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$SETNUMBER(NSET)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE, ONLY : LL_CNTL
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NSET
      INTEGER(4)             :: I
!     **************************************************************************
                          CALL TRACE$PUSH('READCNTL$SETNUMBER')
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'COOP',I)
      NSET=I
      CALL LINKEDLIST$NLISTS(LL_CNTL,'WEIGHT',I)
      NSET=NSET+I
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL$SETNUMBER
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SET$ENOCC(NBB,NKPT_,EIG,OCC)
!     **************************************************************************
!     **  ENERGIES AND OCCUPATIONS ARE TREATED IN THE FOLLOWING DATA MODEL:   **
!     **                                                                      **
!     **  THERE ARE INDEPENDENT ARRAYS FOR THE TWO SPIN DIRECTIONS, EACH      **
!     **  WITH ITS OWN ENERGIES AN OCCUPATIONS.                               **
!     **  1) FOR NDIM=1,NSPIN=1: EACH STATE IS DUPLICATED SO THAT BOTH SPIN   **
!     **     CHANNELS ARE FILLED                                              **
!     **  2) FOR NDIM=1,NSPIN=2: EIGENVALUES AND OCCUPATIONS DIFFER FOR BOTH  **
!     **     SPIN DIRECTIONS.                                                 **
!     **  3) FOR NDIM=2,NSPIN=1: ENERGIES AND OCCUPATIONS ARE DUPLICATED.     **
!     **     BUT THE ENTRIES IN SET ARE SPIN DEPENDENT                        **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NDIM &
     &                          ,NSPIN &
     &                          ,NKPT &
     &                          ,STATEARR,STATE 
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NBB  ! X#(SPIN STATES PER KPOINT)
      INTEGER(4),INTENT(IN) :: NKPT_
      REAL(8)   ,INTENT(OUT):: EIG(NBB,NKPT_)
      REAL(8)   ,INTENT(OUT):: OCC(NBB,NKPT_)
      INTEGER(4)            :: ISPIN,IKPT
      REAL(8)               :: EMAX
      INTEGER(4)            :: NB
!     **************************************************************************
!
!     ==========================================================================
!     == EXTRACT UPPER BOUND OF THE SPECTRUM AND THE NUMBER OF BANDS FOR TEST ==
!     ==========================================================================
      EMAX=-HUGE(EMAX)
      NB=0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
          NB=MAX(NB,STATE%NB)
          EMAX=MAX(EMAX,MAXVAL(STATE%EIG))
        ENDDO
      ENDDO
      IF(NB*2/NDIM.NE.NBB.OR.NKPT_.NE.NKPT) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NBB',NBB)
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$I4VAL('NKPT_',NKPT_)
        CALL ERROR$I4VAL('NKPT',NKPT)
        CALL ERROR$STOP('SET$ENOCC')
      END IF
!      
!     ==========================================================================
!     == FILL IN EIGENVALUES AND OCCUPATIONS                                  ==
!     ==========================================================================
      EIG=EMAX+0.1D0
      OCC=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
          NB=STATE%NB
          EIG(1+NB*(ISPIN-1):NB*ISPIN,IKPT)=STATE%EIG(:)
          OCC(1+NB*(ISPIN-1):NB*ISPIN,IKPT)=STATE%OCC(:)
          IF(NDIM*NSPIN.EQ.1) THEN ! FOR NON-SPIN-POLARIZED CASE
            EIG(NB+1:2*NB,IKPT)=EIG(1:NB,IKPT)
            OCC(1:NB,IKPT)=0.5D0*OCC(1:NB,IKPT)
            OCC(NB+1:2*NB,IKPT)=OCC(1:NB,IKPT)
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$SETS_NEW(NKPT,NSET,NAT,RBAS,ATOMID,RPOS)
!     **************************************************************************
!     ** REMARK: THE STRING VARIABLE SPIN DESCRIBES THE SPIN PROJECTION FOR   **
!     ** NON-COLLINEAR CALCULATIONS. IN THIS CASE, NSPIN=1.                   **
!     **    SPIN='TOTAL' CALCULATES TOTAL ELECTRON DENSITY                    **
!     **    SPIN='MAIN'  CALCULATES THE SPIN DENSITY PROJECTED ONTO THE       **
!     **                 MAIN SPIN AXIS                                       **
!     **    SPIN='X','Y','Z' PRODUCES ONE OF THE CARTESIAN COMPONENTS OF THE  **
!     **                 SPIN DENSITY                                         **
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)  :: NKPT
      INTEGER(4)   ,INTENT(IN)  :: NSET
      INTEGER(4)   ,INTENT(IN)  :: NAT          ! #(ATOMS)
      REAL(8)      ,INTENT(IN)  :: RBAS(3,3)    ! LATTICE VECTORS
      CHARACTER(16),INTENT(IN)  :: ATOMID(NAT)  ! ATOM NAMES
      REAL(8)      ,INTENT(IN)  :: RPOS(3,NAT)  ! ATOMIC POSITIONS
      CHARACTER(32)             :: SETID
      CHARACTER(32)             :: LEGEND
      INTEGER(4)                :: ISET   ! SET COUNTER
      INTEGER(4)                :: IAT
      INTEGER(4)                :: IAT0
      CHARACTER(8)              :: TYPE
      CHARACTER(32)             :: NAME,SPIN
      LOGICAL(4)                :: TCHK,TCHK1
      INTEGER(4)                :: ITH, NUM
      INTEGER(4)                :: JTH, NUMJTH
      INTEGER(4)                :: IORB1,IORB2,IORB
      INTEGER(4)                :: NORB1,NORB2,NORB
      CHARACTER(32)             :: ORBITAL1NAME
      CHARACTER(32)             :: ORBITAL2NAME
!     **************************************************************************
                          CALL TRACE$PUSH('READCNTL$SETS_NEW')
      ISET=0 !COUNTER USE TO SET DEFAULT NAMES FOR SETS
!
!     ==========================================================================
!     ==========================================================================
!     ==  SCAN COOPS                                                          ==
!     ==========================================================================
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'COOP',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNTL,'COOP',ITH)
        ISET=ISET+1
        WRITE(SETID,FMT='("COOP:SET",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ID',0,SETID)
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,SETID)
        CALL NEWSET$NEWSET(SETID)
        CALL NEWSET$SELECT(SETID)
!
!       ========================================================================
!       == READ LEGEND                                                        ==
!       ========================================================================
        WRITE(LEGEND,FMT='("COOP:SET",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'LEGEND',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'LEGEND',1,LEGEND)
        CALL NEWSET$SETCH('LEGEND',LEGEND)
!
!       ========================================================================
!       == DEFINE SPIN AXIS                                                   ==
!       == CAN BE '+Z', '-Z'
!       ========================================================================
        SPIN='+Z'
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',1,TCHK1)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,SPIN)
          SPIN=+SPIN
        END IF
        CALL NEWSET$SETCH('SPINID',SPIN)
!
!       ========================================================================
!       == THE SUM OF COOPS BETWEEN ALL ORB1 AND ALL ORB2 IS EQUAL TO THE     ==
!       == COOP OF THE SUM OF ORBITALS ORB1 WITH THE SUM OF ALL ORBITALS ORB2 ==
!       ========================================================================
        WRITE(ORBITAL1NAME,*)ISET
        ORBITAL2NAME='COOP'//TRIM(ADJUSTL(ORBITAL1NAME))//'ORB2'
        ORBITAL1NAME='COOP'//TRIM(ADJUSTL(ORBITAL1NAME))//'ORB1'
        CALL NEWORBITAL$NEWORBITAL(ORBITAL1NAME)
        CALL NEWORBITAL$NEWORBITAL(ORBITAL2NAME)
        CALL NEWSET$ADDCOOP(ORBITAL1NAME,ORBITAL2NAME)
!
        CALL NEWORBITAL$SELECT(ORBITAL1NAME)
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB1',NORB1)
        DO IORB1=1,NORB1
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB1',IORB1) 
          CALL READONENEWORB(LL_CNTL,RBAS,NAT,ATOMID,RPOS)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL NEWORBITAL$SELECT('NONE')
!
        CALL NEWORBITAL$SELECT(ORBITAL2NAME)
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB2',NORB2)
        DO IORB2=1,NORB2
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB2',IORB2) 
          CALL READONENEWORB(LL_CNTL,RBAS,NAT,ATOMID,RPOS)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL NEWORBITAL$SELECT('NONE')
        CALL NEWSET$SELECT('NONE')
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$PASS('COOPS FINISHED')
!
!     ==========================================================================
!     ==========================================================================
!     ==  SCAN WEIGHTS                                                        ==
!     ==========================================================================
!     ==========================================================================
      CALL LINKEDLIST$NLISTS(LL_CNTL,'WEIGHT',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNTL,'WEIGHT',ITH)
        ISET=ISET+1
        WRITE(SETID,FMT='("WEIGHT",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ID',0,SETID)
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,SETID)
        CALL NEWSET$NEWSET(SETID)
        CALL NEWSET$SELECT(SETID)
!
!       ========================================================================
!       == READ LEGEND                                                        ==
!       ========================================================================
        WRITE(LEGEND,FMT='("WEIGHT",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'LEGEND',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'LEGEND',1,LEGEND)
        CALL NEWSET$SETCH('LEGEND',LEGEND)
!
!       ========================================================================
!       ==  SELECT SPIN AXIS                                                  ==
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',1,TCHK1)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,SPIN)
          SPIN=+SPIN
          CALL NEWSET$SETCH('SPINID',SPIN)
        ELSE
!         ==  SPIN=TOTAL MEANS "SUM OVER ALL SPINOR COMPONENTS". =============
!         ==  FOR COLLINEAR CALCULATIONS, SPIN='Z' PRODUCES THE TWO ==========
!         ==  SPIN DENSITIES. ================================================
          SPIN='Z' 
        END IF
!
!       ========================================================================
!       ==  SELECT TYPE                                                       ==
!       ========================================================================
!
!       == TYPE MAY BE 'TOTAL', 'ALL' OR 'EMPTY' ===============================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TYPE',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
          TYPE=+TYPE
!         == CHECK SYNTAX ======================================================
          CALL LINKEDLIST$NLISTS(LL_CNTL,'ATOM',NUMJTH)
          CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB',NORB)
          IF(NUMJTH.GT.0.OR.NORB.GT.0) THEN
            CALL ERROR$MSG('!DCNTL!WEIGHT:TYPE IS INCOMPATIBLE WITH')
            CALL ERROR$MSG('!DCNTL!WEIGHT!ATOM AND !DCNTL!WEIGHT!ORB')
            CALL ERROR$CHVAL('TYPE',TYPE)
            CALL ERROR$I4VAL('NUMJTH',NUMJTH)
            CALL ERROR$I4VAL('NORB',NORB)
            CALL ERROR$STOP('READCNTL$SETS_NEW')
          END IF
!
!         ======================================================================
!         ==  'TOTAL' = TOTAL DENSITY OF STATES                               ==
!         ======================================================================
          IF(TRIM(TYPE).EQ.'TOTAL') THEN
            CALL NEWSET$SETCH('SPECIAL',TYPE)
!
!         ======================================================================
!         ==  'EMPTY' = VACCUM DENSITY OF STATES                              ==
!         ======================================================================
          ELSE IF(TRIM(TYPE).EQ.'EMPTY') THEN
            CALL NEWSET$SETCH('SPECIAL',TYPE)
!
!         ======================================================================
!         ==  'ALL' = ALL PROJECTED DENSITY OF STATES                         ==
!         ======================================================================
          ELSE IF(TRIM(TYPE).EQ.'ALL') THEN
            CALL NEWSET$ADDLWEIGHT(-1,-1)

!         ======================================================================
!         ==  NOT RECOGNIZED                                                  ==
!         ======================================================================
          ELSE 
            CALL ERROR$MSG('TYPE NOT RECOGNIZED')
            CALL ERROR$MSG('MUST BE EITHER TOTAL,ALL OR EMPTY')
            CALL ERROR$CHVAL('TYPE ',TYPE)
            CALL ERROR$STOP('READCNTL$SETS_NEW')
          END IF
        END IF
!
!       ========================================================================
!       == COLLECT CONTRIBUTIONS FROM INDIVIDUAL ATOMS                        ==
!       ========================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ATOM',NUMJTH)
        DO JTH=1,NUMJTH
          CALL LINKEDLIST$SELECT(LL_CNTL,'ATOM',JTH)
!
!         ==  ATOM SPECIFIER IAT0: IAT0<0 IMPLIES ALL ATOMS              
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,NAME)
          IAT0=-1
          DO IAT=1,NAT
            IF(NAME.EQ.ATOMID(IAT)) THEN
              IAT0=IAT
              EXIT
            END IF
          ENDDO
          IF(IAT0.EQ.-1) THEN
            IF(NAME.NE.'ALL') THEN
              CALL ERROR$MSG('ATOMID NOT RECOGNIZED')
              CALL ERROR$CHVAL('ATOMID_',NAME)
              CALL ERROR$STOP('READCNTL$SETS_NEW')
            END IF
          END IF
!
!         == SELECT TYPE =======================================================
          CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
          TYPE=+TYPE
          IF(TYPE.EQ.'ALL') THEN
            CALL NEWSET$ADDLWEIGHT(IAT,-1)
          ELSE
            IF(SCAN(TYPE,'S').NE.0) CALL NEWSET$ADDLWEIGHT(IAT0,0)
            IF(SCAN(TYPE,'P').NE.0) CALL NEWSET$ADDLWEIGHT(IAT0,1)
            IF(SCAN(TYPE,'D').NE.0) CALL NEWSET$ADDLWEIGHT(IAT0,2)
            IF(SCAN(TYPE,'F').NE.0) CALL NEWSET$ADDLWEIGHT(IAT0,3)

            IF(SCAN(TYPE,'ABCEGHIJKLMNOQRTUVWXYZ0123456789').NE.0) THEN
              CALL ERROR$MSG('ALLOWED VALUES FOR VARIABLE TYPE ARE')
              CALL ERROR$MSG('ONLY "ALL" OR A COMBINATION OF "S", "P", "D","F"')
              CALL ERROR$CHVAL('ATOM',NAME)
              CALL ERROR$CHVAL('TYPE',TYPE)
              CALL ERROR$STOP('READCNTL$SETS_NEW')
            END IF
          END IF
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
!       ========================================================================
!       == ADD CONTRIBUTION FROM ORBITALS                                     ==
!       ========================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB',NORB)
        DO IORB=1,NORB
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB',IORB)
          WRITE(ORBITAL1NAME,*)ISET
          WRITE(ORBITAL2NAME,*)IORB
          ORBITAL1NAME='WGHT'//TRIM(ADJUSTL(ORBITAL1NAME))//'ORB'
          ORBITAL1NAME=TRIM(ADJUSTL(ORBITAL1NAME))//TRIM(ADJUSTL(ORBITAL2NAME))
          CALL NEWSET$ADDORBWGHT(ORBITAL1NAME)
          CALL NEWORBITAL$NEWORBITAL(ORBITAL1NAME)
!
          CALL NEWORBITAL$SELECT(ORBITAL1NAME)
          CALL READONENEWORB(LL_CNTL,RBAS,NAT,ATOMID,RPOS)
          CALL NEWORBITAL$SELECT('NONE')
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        CALL NEWSET$SELECT('NONE')
      ENDDO
                          CALL TRACE$POP
      RETURN
    END SUBROUTINE READCNTL$SETS_NEW
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINBRA(X,Y,Z,CP,CM)
!     **************************************************************************
!     ** BRAS OF SPIN EIGENSTATES FOR A SPIN-AXIS (X,Y,Z)                     **
!     **                                                                      **
!     ** THE TWO CASES EZ>0 AND EZ<0 ARE TREATED DIFFERENTLY TO AVOID A       **
!     ** DIVIDE BY ZERO                                                       **
!     **                                                                      **
!     ** THE EIGENSTATES ARE DEFINED ONLY UP TO A PHASE FACTOR.               **
!     ** IT IS FIXED BY:                                                      **
!     **   (1) FOR AN AXIS IN Z DIRECTION, THE STATES ARE (1,0) AND (0,1)    **
!     **   (2) FOR Z>0 THE CP(1) AND CM(2) IS REAL AND POSITIVE
!     **   (3) FOR Z<0 THE CP(2) AND CM(1) IS REAL AND POSITIVE
!     **   (4) CP(-X,-Y,-Z)=CM(X,Y,Z) AND VICE VERSA                          **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: X,Y,Z
      COMPLEX(8),INTENT(OUT) :: CP(2)
      COMPLEX(8),INTENT(OUT) :: CM(2)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)                :: EX,EY,EZ
      REAL(8)                :: SVAR
!     **************************************************************************
      SVAR=SQRT(X**2+Y**2+Z**2)
      EX=X/SVAR
      EY=Y/SVAR
      EZ=Z/SVAR
!     == EIGENVECTORS OF EZ*SIGMAX+EY*SIGMAY+EZ*SIGMAZ
!     == CP IS EIGENVECTOR FOR EIGENVALUE +1
!     == CM IS EIGENVECTOR FOR EIGENVALUE -1
      IF(EZ.GT.0.D0) THEN
        CP(1)=1.D0+EZ
        CP(2)=+EX+CI*EY
        CM(1)=-EX+CI*EY
        CM(2)=1.D0+EZ
      ELSE
        CP(1)=+EX-CI*EY
        CP(2)=1.D0-EZ
        CM(1)=1.D0-EZ
        CM(2)=-EX-CI*EY
      END IF
!     == NORMALIZE =============================================================
      CP=CP/SQRT(ABS(CP(1))**2+ABS(CP(2))**2)
      CM=CM/SQRT(ABS(CM(1))**2+ABS(CM(2))**2)
!
!     == TURN KETS INTO BRAS =================================================
      CP=CONJG(CP)
      CM=CONJG(CM)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$OUTPUT(EMIN,EMAX,NE,EBROAD,PREFIX &
     &                    ,NBB,NKPT,EIG,OCC,NSET,SET,MODE)
!     **************************************************************************
!     **************************************************************************
      USE STRINGS_MODULE
      USE READCNTL_MODULE
      USE DOS_WGHT_MODULE
      USE NEWSET_MODULE  ,ONLY : NEWSET
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NE   ! #(ENERGY-GRID POINTS)
      REAL(8)      ,INTENT(IN) :: EMIN ! MIN ENERGY ON GRID
      REAL(8)      ,INTENT(IN) :: EMAX
      REAL(8)      ,INTENT(IN) :: EBROAD
      CHARACTER(*) ,INTENT(IN) :: PREFIX ! FOR DOS AND NOS FILES
      INTEGER(4)   ,INTENT(IN) :: NBB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      INTEGER(4)   ,INTENT(IN) :: NSET
      REAL(8)      ,INTENT(IN) :: EIG(NBB,NKPT)
      REAL(8)      ,INTENT(IN) :: OCC(NBB,NKPT)
      REAL(8)      ,INTENT(IN) :: SET(NBB,NKPT,2,NSET)
      CHARACTER(32),INTENT(IN) :: MODE
      LOGICAL(4)   ,PARAMETER  :: TDOS=.TRUE.  ! DENSITY OF STATES PRINTED ?
      LOGICAL(4)   ,PARAMETER  :: TNOS=.FALSE.  ! NUMBER OF STATES PRINTED ?
      CHARACTER(32)        :: LEGEND
      CHARACTER(256)       :: FILE
      INTEGER(4)           :: NFILDOS  ! UNIT FOR DENSITY OF STATES FILE
      INTEGER(4)           :: NFILO
      REAL(8)              :: EV
      REAL(8)              :: E1,E2
      INTEGER(4)           :: IE1,IE2
      LOGICAL(4),ALLOCATABLE :: DEADZONE(:) !(NE,ISPIN)
      INTEGER(4)           :: IB,ISET
!     **************************************************************************
                         CALL TRACE$PUSH('READCNTL$OUTPUT')
      CALL CONSTANTS('EV',EV)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==========================================================================
!     ==  DETERMINE DEAD ZONES. IN A DEADZONE THE DENSITY OF STATES VANISHES  ==
!     == AND NEED NOT BE PRINTED. DEADZONES MUST BE THE SAME FOR ALL SETS     ==
!     == TO ALLOW DATASET OPERATIONS, FOR EXAMPLE WITH XMGRACE.               ==
!     ==========================================================================
      ALLOCATE(DEADZONE(NE))
      DEADZONE(:)=.TRUE.
      DO IB=1,NBB
        E1=MINVAL(EIG(IB,:))-1.D0*EV-5.D0*EBROAD
        E2=MAXVAL(EIG(IB,:))+1.D0*EV+5.D0*EBROAD
        IE1=MAX(1, 1+NINT((E1-EMIN)/(EMAX-EMIN)*REAL(NE-1)))
        IE2=MIN(NE,1+NINT((E2-EMIN)/(EMAX-EMIN)*REAL(NE-1)))
        DEADZONE(IE1:IE2)=.FALSE.
      ENDDO
!
!     ==========================================================================
!     ==========================================================================
!     == PRINT ALL COOPS, DOS AND NOS                                         ==
!     ==========================================================================
!     ==========================================================================
      DO ISET=1,NSET
!
!       ========================================================================
!       ==  SPECIFY OUTPUT FILE                                               ==
!       ========================================================================
        FILE=TRIM(PREFIX)//TRIM(NEWSET(ISET)%ID)//-'.DOS'
        CALL FILEHANDLER$SETFILE('PDOSOUT',.FALSE.,FILE)
        CALL FILEHANDLER$UNIT('PDOSOUT',NFILDOS)
!
!       ========================================================================
!       ==  EXTRACT LEGEND                                                    ==
!       ========================================================================
        LEGEND=NEWSET(ISET)%LEGEND
!
!       ========================================================================
!       ==  WRITE DOS AND INTEGRATED DOS ON FILE                              ==
!       ========================================================================
        IF(MODE.EQ.'SAMPLE')THEN
          CALL PUTONGRID_SAMPLE(NFILDOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
      &                 ,NBB,NKPT,EIG,OCC,SET(:,:,:,ISET),LEGEND)
        ELSE IF(MODE.EQ.'TETRA')THEN
!!$          IF((MAXVAL(SET(:,:,:,ISET)).NE.1.0D0.OR.&
!!$      &       MINVAL(SET(:,:,:,ISET)).NE.1.0D0).AND.SPACEGROUP.NE.1)THEN
!!$            CALL ERROR$MSG('TETRAEDRON METHOD ONLY IMPLEMENTED FOR')
!!$            CALL ERROR$MSG('SPACEGROUP=1 OR TOTAL DENSITY OF STATES')
!!$            CALL ERROR$I4VAL('ISET',ISET)   
!!$            CALL ERROR$I4VAL('SPACEGROUP',SPACEGROUP)   
!!$            CALL ERROR$STOP('READCNTL$OUTPUT')
!!$          ENDIF
!
          CALL PUTONGRID_TETRA(NFILDOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
      &                 ,NBB,NKPT,EIG,SET(:,:,:,ISET),LEGEND)
        ELSE
          CALL ERROR$MSG('VALUE FOR "MODE" NOT RECOGNIZED')
          CALL ERROR$MSG('MUST BE "SAMPLE" OR "TETRA"')
          CALL ERROR$CHVAL('MDOE',MODE)   
          CALL ERROR$STOP('READCNTL$OUTPUT')
        ENDIF
!
!       == CLOSE DOWN ==========================================================
        IF(NFILDOS.GE.0)CALL FILEHANDLER$CLOSE('PDOSOUT')
      ENDDO  ! END OF LOOP OVER SETS
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL$OUTPUT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PUTONGRID_SAMPLE(NFILDOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
     &                    ,NBB,NKPT,EIG,OCC,SET,LEGEND)
!     **************************************************************************
!     **  MAPS THE CONTRIBUTION FROM EACH STATE ONTO AN ENERGY GRID,          **
!     **  CONSTRUCTS DOS AND NOS AND WRITES THE RESULT ON FILE                **
!     **                                                                      **
!     **  SET IS THE CONTRIBUTION FROM EACH STATE (WITHOUT K-POINT WEIGHT     **
!     **       AND OCCUPATIONS MULTIPLIED TO THEM)                            **
!     **  STATE(IK,IS)%OCC(IB) IS THE OCCUPATION MULTIPLIED WITH THE          **
!     **     WITH THE K-POINT WEIGHT                                          **
!     **  NOS(IE,ISPIN,1) IS MULTIPLIED WITH MAX OCCUPATION OF ALL BANDS      **
!     **  NOS(IE,ISPIN,2) IS MULTIPLIED WITH ACTUAL OCCUPATION OF EACH STATE  **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE, ONLY: STATE,STATEARR,WKPT,NSPIN
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NE
      REAL(8)      ,INTENT(IN) :: EMIN
      REAL(8)      ,INTENT(IN) :: EMAX
      REAL(8)      ,INTENT(IN) :: EBROAD
      INTEGER(4)   ,INTENT(IN) :: NBB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      REAL(8)      ,INTENT(IN) :: EIG(NBB,NKPT)
      REAL(8)      ,INTENT(IN) :: OCC(NBB,NKPT)
      REAL(8)      ,INTENT(IN) :: SET(NBB,NKPT,2)
      INTEGER(4)   ,INTENT(IN) :: NFILDOS  
      CHARACTER(32),INTENT(IN) :: LEGEND
      LOGICAL(4)   ,INTENT(IN) :: DEADZONE(NE)
      REAL(8)      ,PARAMETER  :: TOL=1.D-2
      REAL(8)                  :: DE
      INTEGER(4)               :: IE1,IE2,IDE
      INTEGER(4)               :: ND,IOCC
      REAL(8)                  :: NOS(NE,2,2)   
      REAL(8)                  :: DOS(NE,2,2)
      REAL(8)      ,ALLOCATABLE:: SMEAR(:)
      REAL(8)                  :: EV
      REAL(8)                  :: W1,W2,X
      REAL(8)                  :: WGHTX
      INTEGER(4)               :: IKPT,ISPIN,IE,IB
      REAL(8)                  :: SVAR
      REAL(8)                  :: SIG
      REAL(8)                  :: E
      REAL(8)                  :: OCC1
!     **************************************************************************
                                 CALL TRACE$PUSH('PUTONGRID_SAMPLE')
      CALL CONSTANTS('EV',EV)
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)   !STEP OF THE ENERGY GRID
!
      IF(SUM(WKPT).EQ.0.D0) THEN
        CALL ERROR$MSG('NO K-POINT WEIGHTS AVAILABLE: OLD VERSION OF PDOS FILE')
        CALL ERROR$MSG('RERUN ONE PAW ITERATION WITH NEW CODE')
        CALL ERROR$STOP('PUTONGRID_SAMPLE')
      END IF
!
!     ==========================================================================
!     ==  MAP CONTRIBUTION FROM EACH STATE) ONTO THE ENERGY GRID.             ==
!     ==  (IT IS DIVIDED PROPORTIONALLY TO THE TWO ENCLOSING GRID POINTS)     ==
!     == THE SUM OVER ALL NOS-DATA ADDS TO THE TOTAL NUMBER OF STATES         ==
!     == STATES LYING BELOW EMIN ARE MAPPED INTO NOSSMALL                     ==
!     ==========================================================================
      NOS(:,:,:)=0.D0
      DO ISPIN=1,2
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
!         == CAUTION: HERE I ESTIMATE THE WEIGHT AND SPIN-DEGENERACY FACTOR 
!         == FROM THE MAX OCCUPATION, WHICH MAY BE INCORRECT
          WGHTX=WKPT(IKPT)
          DO IB=1,NBB
            X=(EIG(IB,IKPT)-EMIN)/DE+1.D0
            IE1=INT(X)
            IE2=IE1+1
            W2=(X-REAL(IE1,KIND=8))
            W1=1.D0-W2
            OCC1=OCC(IB,IKPT)
            IF(IE1.LE.NE.AND.IE1.GE.1) THEN 
              NOS(IE1,ISPIN,1)=NOS(IE1,ISPIN,1)+W1*SET(IB,IKPT,ISPIN)*WGHTX
              NOS(IE1,ISPIN,2)=NOS(IE1,ISPIN,2)+W1*SET(IB,IKPT,ISPIN)*OCC1
            END IF
            IF(IE2.LE.NE.AND.IE2.GE.1) THEN
              NOS(IE2,ISPIN,1)=NOS(IE2,ISPIN,1)+W2*SET(IB,IKPT,ISPIN)*WGHTX
              NOS(IE2,ISPIN,2)=NOS(IE2,ISPIN,2)+W2*SET(IB,IKPT,ISPIN)*OCC1
            END IF
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE DOS                                                       ==
!     ==========================================================================
!     == DETERMINE SMEARING FUNCTION (DERIVATIVE OF FERMI FUNCTION) ============
      ND=NINT(EBROAD/DE * LOG(4.D0/TOL))
      ALLOCATE(SMEAR(-ND:ND))
      DO IDE=-ND,ND
        SMEAR(IDE)=EBROAD/( 0.5D0*COSH(0.5D0*REAL(IDE,KIND=8)*DE/EBROAD) )**2
      ENDDO
      SMEAR=SMEAR/SUM(SMEAR)  ! RENORMALIZE TO MAINTAIN SUM RULES
!
!     ==  SMEAR DENSITY OF STATES ==============================================
      DOS(:,:,:)=0.D0
      DO ISPIN=1,2
        DO IOCC=1,2
          DO IDE=-ND,ND
            IE1=MAX(1,1-IDE)
            IE2=MIN(NE,NE-IDE)
            W1=SMEAR(IDE)
            DO IE=IE1,IE2
              DOS(IE,ISPIN,IOCC)=DOS(IE,ISPIN,IOCC)+NOS(IE+IDE,ISPIN,IOCC)*W1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  RESCALE DOS                                                         ==
!     ==========================================================================
      DOS=DOS/DE
!
!     ==========================================================================
!     ==  ROUND INSIGNIFICANT VALUES TO ZERO                                  ==
!     ==========================================================================
      DO ISPIN=1,2
        DO IOCC=1,2
          DO IE=1,NE
            IF(ABS(DOS(IE,ISPIN,IOCC)).LE.1.D-99)DOS(IE,ISPIN,IOCC)=0.D0
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  WRITE RESULT ON PSODOUT                                             ==
!     ==========================================================================
!     == WRITE DENSITY OF STATES ===============================================
      DO ISPIN=1,2
        SIG=REAL(3-2*ISPIN) ! +1 FOR ISPIN=1 / -1 FOR ISPIN=2
        WRITE(NFILDOS,FMT='(F14.8,2F20.8)')EMIN/EV,0.D0,0.D0
        DO IE=1,NE
          IF(DEADZONE(IE)) CYCLE
          E=EMIN+(EMAX-EMIN)*REAL(IE-1)/REAL(NE-1)
          WRITE(NFILDOS,FMT='(F14.8,2F20.8)')E/EV,SIG*DOS(IE,ISPIN,1)*EV &
                                                 ,SIG*DOS(IE,ISPIN,2)*EV
        ENDDO
        WRITE(NFILDOS,FMT='(F14.8,2F20.8)')EMAX/EV,0.D0,0.D0
      ENDDO
      WRITE(NFILDOS,FMT='("# THIS WAS: ",A)')LEGEND
                                 CALL TRACE$POP
      RETURN
      END SUBROUTINE PUTONGRID_SAMPLE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PUTONGRID_TETRA(NFILDOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
     &                          ,NBB,NKPT,EIG,SET,LEGEND)
!     **************************************************************************
!     **  MAPS THE CONTRIBUTION FROM EACH STATE ONTO AN ENERGY GRID,          **
!     **  CONSTRUCTS DOS AND NOS WITH THE TETRAHEDRON METHOD AND WRITES THE   **
!     **  RESULT ON FILE                                                      **
!     **                                                                      **
!     **  NOS(IE,ISPIN,1) IS MULTIPLIED WITH MAX OCCUPATION OF ALL BANDS      **
!     **  NOS(IE,ISPIN,2) IS MULTIPLIED WITH ACTUAL OCCUPATION OF EACH STATE  **
!     **                                                                      **
!     **  EWGHT ARE ENERGY-RESOLVED DENSITY OF STATES FROM A SPECIFIC K-POINT **
!     **  AND BAND INDEX. (ENERGY DERIVATIVE OF INTEGRATION WEIGHTS)          **
!     **  THE STRUCTURE EWGHT_TYPE CONTAINS THE DATA WGHT ON A SUBSECTION     **
!     **  (I1,...,I2) OF THE ENERGY GRID                                      **
!     **************************************************************************
      USE DOS_WGHT_MODULE, ONLY: EF &
     &                          ,EWGHT &
     &                          ,WGHT
      USE PDOS_MODULE    , ONLY: STATEARR
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NE     ! #(ENERGY GRID POINTS)
      REAL(8)      ,INTENT(IN) :: EMIN   ! MIN ENERGY ON GRID
      REAL(8)      ,INTENT(IN) :: EMAX   ! MAX ENERGY ON GRID
      REAL(8)      ,INTENT(IN) :: EBROAD ! KBT OF THERMAL BROADENIUNG
      INTEGER(4)   ,INTENT(IN) :: NBB    ! #(BANDS PER K-POINT)
                                         ! NBB=2*NB FOR NSPIN=1,2
      INTEGER(4)   ,INTENT(IN) :: NKPT   ! #(K-POINTS)
      REAL(8)      ,INTENT(IN) :: EIG(NBB,NKPT)   ! ONE-PARTICLE ENERGIES
      REAL(8)      ,INTENT(IN) :: SET(NBB,NKPT,2) ! MATRIX ELEMENTS
      INTEGER(4)   ,INTENT(IN) :: NFILDOS ! UNIT FOR DOS FILE OR "-1"
      CHARACTER(32),INTENT(IN) :: LEGEND
      LOGICAL(4)   ,INTENT(IN) :: DEADZONE(NE) !INDICATES POINTS TO BE SKIPPED
      REAL(8)      ,PARAMETER  :: TOL=1.D-2
      REAL(8)                  :: NOS(NE,2)
      REAL(8)                  :: DOS(NE,2)
      REAL(8)    ,ALLOCATABLE  :: SMEAR(:)
      REAL(8)                  :: DE
      INTEGER(4)               :: IE1,IE2,IDE,I1,I2
      INTEGER(4)               :: ND
      REAL(8)                  :: EV
      REAL(8)                  :: W1
      INTEGER(4)               :: IKPT,ISPIN,IE,IB
      REAL(8)                  :: SVAR
      REAL(8)                  :: SIG
      REAL(8)                  :: E
!     **************************************************************************
                                 CALL TRACE$PUSH('PUTONGRID_TETRA')
      CALL CONSTANTS('EV',EV)
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)  !STEP OF THE ENERGY GRID
!
!     ==========================================================================
!     ==  CALCULATE DENSITY OF STATES                                         ==
!     ==========================================================================
      DOS(:,:)=0.D0
      DO ISPIN=1,2   ! THE DATA MODEL HAS BEEN EXPANDED TO 2 SPIN COMPONENTS
        DO IKPT=1,NKPT
          DO IB=1,NBB
            SVAR=SET(IB,IKPT,ISPIN)*DE
            I1=EWGHT(IB,IKPT)%I1
            I2=EWGHT(IB,IKPT)%I2
            DO IE=I1,I2
              DOS(IE,ISPIN)=DOS(IE,ISPIN)+EWGHT(IB,IKPT)%WGHT(IE)*SVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  BROADEN RESULT                                                      ==
!     ==========================================================================
!     == DETERMINE SMEARING FUNCTION (DERIVATIVE OF FERMI FUNCTION) ============
!     == BROADENING EXTENDS OVER 2N+1 STEPS
      ND=NINT(EBROAD/DE * LOG(4.D0/TOL))
      ALLOCATE(SMEAR(-ND:ND))
      DO IDE=-ND,ND
        SMEAR(IDE)=EBROAD/( 0.5D0*COSH(0.5D0*REAL(IDE,KIND=8)*DE/EBROAD) )**2
      ENDDO
      SMEAR=SMEAR/SUM(SMEAR)  ! RENORMALIZE TO MAINTAIN SUM RULES
!
!     ==  SMEAR OUT THE DENSITY OF STATES. (NOS IS ONLY A SUPPORT ARRAY.) ======
      NOS=DOS
      DOS(:,:)=0.D0
      DO ISPIN=1,2
        DO IDE=-ND,ND
          IE1=MAX(1,1-IDE)
          IE2=MIN(NE,NE-IDE)
          W1=SMEAR(IDE)
          DO IE=IE1,IE2
            DOS(IE,ISPIN)=DOS(IE,ISPIN)+NOS(IE+IDE,ISPIN)*W1
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  CONVERT DELTA-NOS INTO DOS                                          ==
!     ==========================================================================
      DOS=DOS/DE
!
!     ==========================================================================
!     ==  ROUND INSIGNIFICANT VALUES TO ZERO                                  ==
!     ==========================================================================
      DO ISPIN=1,2
        DO IE=1,NE
          IF(ABS(DOS(IE,ISPIN)).LE.1.D-99) DOS(IE,ISPIN)=0.D0
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  WRITE RESULT ON FILE WITH ID PDOSOUT                                ==
!     ==  1) THE DOS OF THE FIRST SPIN DIRECTION IS WRITTEN WITH POSTIVE SIGN.==
!     ==  2) THE VALUES JUMP TO (EMAX,0) AND (EMIN,0).                        ==
!     ==  3) THE DOS FOR THE SECOND SPIN DIRECTION IS DRAWN WITH NEGATIVE SIGN==
!     ==========================================================================
      DO ISPIN=1,2
        SIG=REAL(3-2*ISPIN) ! +1 FOR ISPIN=1 / -1 FOR ISPIN=2
        WRITE(NFILDOS,FMT='(F14.8,2F20.8)')EMIN/EV,0.D0,0.D0
        DO IE=1,NE
          IF(DEADZONE(IE)) CYCLE
          E=EMIN+(EMAX-EMIN)*REAL(IE-1,KIND=8)/REAL(NE-1)
          SVAR=1.D0
          IF(E.GT.EF) SVAR=0.D0
          WRITE(NFILDOS,FMT='(F14.8,2F20.8)')E/EV,SIG*DOS(IE,ISPIN)*EV &
                                                 ,SIG*DOS(IE,ISPIN)*SVAR*EV
        ENDDO
        WRITE(NFILDOS,FMT='(F14.8,2F20.8)')EMAX/EV,0.D0,0.D0
      ENDDO
      WRITE(NFILDOS,FMT='("# THIS WAS: ",A)')LEGEND
                                 CALL TRACE$POP
      RETURN
      END SUBROUTINE PUTONGRID_TETRA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GENERATE_TETRA_WGHT(NFILO,NBB,NKPT,EMAX,EMIN,NE,RBAS,EIG &
     &                              ,ELSCALE)
!     **************************************************************************
!     **  EWGHT ARE ENERGY-RESOLVED DENSITY OF STATES FROM A SPECIFIC K-POINT **
!     **  AND BAND INDEX. (ENERGY DERIVATIVE OF INTEGRATION WEIGHTS)          **
!     **  THE STRUCTURE EWGHT_TYPE CONTAINS THE DATA WGHT ON A SUBSECTION     **
!     **  (I1,...,I2) OF THE ENERGY GRID                                      **
!     **                                                                      **
!     **************************************************************************
      USE DOS_WGHT_MODULE , ONLY: EF,SPACEGROUP,WGHT,EWGHT
      USE BRILLOUIN_MODULE, ONLY: EWGHT_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)        :: NFILO
      INTEGER(4),INTENT(IN)        :: NBB
      INTEGER(4),INTENT(IN)        :: NKPT
      REAL(8)   ,INTENT(IN)        :: EMIN
      REAL(8)   ,INTENT(IN)        :: EMAX
      INTEGER(4),INTENT(IN)        :: NE
      REAL(8)   ,INTENT(IN)        :: RBAS(3,3)
      REAL(8)   ,INTENT(IN)        :: EIG(NBB,NKPT)
      REAL(8)   ,INTENT(IN)        :: ELSCALE
      INTEGER(4)                   :: ISHIFT(3)
      INTEGER(4)                   :: NKDIV(3)
      INTEGER(4)                   :: NKPT2
      LOGICAL(4)                   :: TINV
      REAL(8)                      :: RNTOT,NEL
      REAL(8)                      :: SUMA(2),SVAR
      REAL(8),ALLOCATABLE          :: A(:,:)
      INTEGER(4)                   :: IB,IKPT,ISPIN
      REAL(8)                      :: A0,B0,C0,ALPHA,BETA,GAMMA
      INTEGER(4)                   :: NSYM,NOP
      INTEGER(4),PARAMETER         :: NOPX=48
      INTEGER(4)                   :: IARB(3)
      CHARACTER(3)                 :: BRAVAIS
      INTEGER(4)                   :: IIO(3,3,NOPX)
      REAL(8)                      :: C(3,NOPX)
      INTEGER(4)                   :: ISYM
      LOGICAL(4)                   :: TSHIFT
      REAL(8)                      :: EV
!     **************************************************************************
                          CALL TRACE$PUSH('GENERATE_TETRA_WGHT')
      CALL CONSTANTS('EV',EV)
!     
!     ==========================================================================
!     == COLLECT TETRAHEDRON METHOD RELATED INFORMATION FROM PDOS FILE        ==
!     ==========================================================================
      CALL PDOS$GETI4A('NKDIV',3,NKDIV)
      CALL PDOS$GETI4A('ISHIFT',3,ISHIFT)
!     == FOR NON-SPIN POLARIZED CALCULATIONS, RNTOT IS ONLY ONE-HALF OF THE ====
!     == NUMBER OF ELECTRONS (SET IN PAW_WAVES2.F90). HERE IT IS BLOWN UP ======
!     == TO THE FULL NUMBER OF ELECTRONS =======================================
      CALL PDOS$GETR8('RNTOT',RNTOT)
      RNTOT=RNTOT*ELSCALE
      CALL PDOS$GETR8('NEL',NEL)
      CALL PDOS$GETL4('TINV',TINV)
      CALL PDOS$GETI4('SPACEGROUP',SPACEGROUP)
      CALL PDOS$GETL4('TSHIFT',TSHIFT)
      CALL SPACEGROUP$SETI4('SPACEGROUP',SPACEGROUP)
      CALL SPACEGROUP$GETCH('BRAVAIS',BRAVAIS)
!     
!     ==========================================================================
!     == REPORT INFORMATION ON THE PDOS FILE RELATED TO THE TETRAHEDRON METHOD==
!     ==========================================================================
      CALL REPORT$TITLE(NFILO,'TETRAHEDRON METHOD')
      CALL REPORT$STRING(NFILO,'WARNING!!')
      CALL REPORT$STRING(NFILO,'PROJECTED DENSITY OF STATES WILL BE IN ERROR')
      CALL REPORT$STRING(NFILO,'UNLESS THEY TRANSFORM LIKE THE IDENTITY')
      CALL REPORT$STRING(NFILO,'UNDER THE POINTGROUP OF THE CRYSTAL EMPLOYED')
      CALL REPORT$STRING(NFILO,'IN THE CALCULATION')
      CALL REPORT$I4VAL(NFILO,"SPACEGROUP",SPACEGROUP,' ')
      CALL REPORT$CHVAL(NFILO,"BRAVAIS LATTICE",BRAVAIS)
      CALL REPORT$L4VAL(NFILO,"TIME REVERSAL SYMMETRY EXPLOITED?",TINV)
      CALL REPORT$L4VAL(NFILO,"TSHIFT",TSHIFT)
      CALL REPORT$I4VAL(NFILO,"NUMBER OF K-GRID POINTS ALONG G1",NKDIV(1),' ')
      CALL REPORT$I4VAL(NFILO,"NUMBER OF K-GRID POINTS ALONG G2",NKDIV(2),' ')
      CALL REPORT$I4VAL(NFILO,"NUMBER OF K-GRID POINTS ALONG G3",NKDIV(3),' ')
      CALL REPORT$I4VAL(NFILO,"SHIFT ALONG G1",ISHIFT(1),' ')
      CALL REPORT$I4VAL(NFILO,"SHIFT ALONG G2",ISHIFT(2),' ')
      CALL REPORT$I4VAL(NFILO,"SHIFT ALONG G3",ISHIFT(3),' ')
!     
!     ==========================================================================
!     == CONSTRUCT SYMMETRY OPERATIONS                                        ==
!     ==========================================================================
      CALL BRILLOUIN$CHECKRBAS(BRAVAIS,A0,B0,C0,ALPHA,BETA,GAMMA,RBAS)
      CALL SPACEGROUP$GENERATORS('RECI',NOPX,NOP,IIO,C)
!     
!     ==========================================================================
!     == CONSTRUCT IRREDUCIBLE K-POINTS                                       ==
!     ==========================================================================
!!$      IF(BRAVAIS.EQ.'GH'.OR.BRAVAIS.EQ.'GQ'.OR.BRAVAIS.EQ.'GOB') THEN
!!$        IARB=(/1,0,0/)
!!$      ELSE IF(BRAVAIS.EQ.'GOF'.OR.BRAVAIS.EQ.'GO'.OR.BRAVAIS.EQ.'GM') THEN
!!$        IARB=(/0,0,0/)
!!$      ELSE IF(BRAVAIS.EQ.'GMB') THEN
!!$        IARB=(/0,1,0/)
!!$      ELSE
!!$        IARB=(/1,1,1/)
!!$      ENDIF 
!!$      NKPT2=(NKDIV(1)+1)*(NKDIV(2)+1)*(NKDIV(3)+1)
!!$      CALL BRILLOUIN$MSH(RBAS,NKPT2,NOP,IIO,IARB,TSHIFT)
!!$      CALL BRILLOUIN$MSHNOSYM(TINV,RBAS,NKDIV,ISHIFT)
      CALL BRILLOUIN$MSHSYM(RBAS,NKDIV,ISHIFT,TINV,NOP,IIO)
      CALL BRILLOUIN$GETI4('NK',NKPT2)
!
      IF(NKPT2.NE.NKPT)THEN
        CALL ERROR$MSG('NUMBER OF KPOINTS INCONSISTENT')
        CALL ERROR$I4VAL('NKPT FROM PDOS',NKPT)
        CALL ERROR$I4VAL('NKPT FROM BRILLOUIN',NKPT2)
        CALL ERROR$I4VAL('NKDIV(1)',NKDIV(1))
        CALL ERROR$I4VAL('NKDIV(2)',NKDIV(2))
        CALL ERROR$I4VAL('NKDIV(3)',NKDIV(3))
        CALL ERROR$STOP('DOS')
      ENDIF
!     
!     ==========================================================================
!     == CONSTRUCT K-INTEGRATION WEIGHTS                                      ==
!     ==========================================================================
      ALLOCATE(WGHT(NBB,NKPT))
      CALL BRILLOUIN$DOS(NBB,NKPT,EIG,WGHT,RNTOT,EF)
!!$!
!!$!     ==========================================================================
!!$!     ==  PERFORM BRILLOUIN ZONE INTEGRATION OF A(K) FOR TESTING              ==
!!$!     ==========================================================================
!!$      !FIXME TOTAL DENSITY FOR TESTING
!!$      ALLOCATE(A(NB*NSPIN,NKPT))
!!$      A=1.D0
!!$      SUMA(:)=0.D0
!!$      DO IB=1,NB
!!$        DO ISPIN=1,NSPIN
!!$          SUMA(ISPIN)=0.0D0
!!$          DO IKPT=1,NKPT
!!$            SUMA(ISPIN)=SUMA(ISPIN) &
!!$     &                 +WGHT(IB+NB*(ISPIN-1),IKPT)*A(IB+NB*(ISPIN-1),IKPT)
!!$          ENDDO
!!$          PRINT*,"IB=",IB," ISPIN=",ISPIN," SUMA=",SUMA(ISPIN)
!!$        ENDDO
!!$      ENDDO
!!$      
!!$      A=1.0D0
!!$      SUMA(:)=0.D0
!!$      DO IB=1,NB
!!$        DO ISPIN=1,NSPIN
!!$          DO IKPT=1,NKPT
!!$            SUMA(ISPIN)=SUMA(ISPIN) &
!!$     &                 +WGHT(IB+NB*(ISPIN-1),IKPT)*A(IB+NB*(ISPIN-1),IKPT)
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      PRINT*,'INTEGRAL OF A=1 :             ',SUM(SUMA(:)),' SHOULD BE ',RNTOT 
!!$      PRINT*,'INTEGRAL OF A=1 (SPIN UP) :   ',SUMA(1) 
!!$      PRINT*,'INTEGRAL OF A=1 (SPIN DOWN) : ',SUMA(2)
!     
!     ==========================================================================
!     == CALCULATE ENERGY DEPENDENT INTEGRATION WEIGHTS                       ==
!     ==========================================================================
      ALLOCATE(EWGHT(NBB,NKPT))
      CALL BRILLOUIN$WGHT(NKPT,NBB,EMIN,EIG,WGHT)
      CALL BRILLOUIN$EWGHT(NKPT,NBB,EIG,EMIN,EMAX,NE,EWGHT)
!     
!     ==========================================================================
!     == REPORT INFORMATION                                                   ==
!     ==========================================================================
      CALL REPORT$R8VAL(NFILO,"FERMI LEVEL",EF/EV,'EV')
                          CALL TRACE$POP()
      RETURN
      END SUBROUTINE GENERATE_TETRA_WGHT

