!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**  NAME: PDOS                                                               **
!**                                                                           **
!**  PURPOSE: ANALYSIS TOOL FOR DENSITY OF STATES                             **
!**                                                                           **
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
USE LINKEDLIST_MODULE
TYPE(LL_TYPE)   :: LL_CNTL
SAVE
END MODULE READCNTL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM PDOS
!     **************************************************************************
!     **************************************************************************
      USE PDOS_MODULE, ONLY: STATE,STATEARR
      USE SPINDIR_MODULE
      USE DOS_WGHT_MODULE
      IMPLICIT NONE
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NAT
      INTEGER(4)                :: NB
      INTEGER(4)                :: NKPT
      INTEGER(4)                :: NSPIN
      INTEGER(4)                :: NDIM !=2 FOR SPINOR WF; OTHERWISE =1
      INTEGER(4)                :: LENG
      INTEGER(4)                :: NSET
      INTEGER(4)   ,ALLOCATABLE :: LMX(:)
      REAL(8)      ,ALLOCATABLE :: RPOS(:,:)
      REAL(8)                   :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)      ,ALLOCATABLE :: EIG(:,:)
      REAL(8)      ,ALLOCATABLE :: OCC(:,:)
      REAL(8)      ,ALLOCATABLE :: SET(:,:,:,:)
      CHARACTER(32),ALLOCATABLE :: LEGEND(:)
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
      LOGICAL(4)                :: TCHK
      CHARACTER(256)            :: PREFIX !OPTIONAL PREFIX FOR DOS AND NOS FILES
!      REAL(8)      ,ALLOCATABLE :: SPINDIR(:,:)
!     **************************************************************************
      CALL TRACE$PUSH('MAIN')
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
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"           PDOS ANALYSIS TOOL                ")')
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"    FOR THE PROJECTOR-AUGMENTED WAVE METHOD  ")')
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(T20 &
     &               ," P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY ")')
      WRITE(NFILO,FMT='(T20 &
     &      ,"(C) CLAUSTHAL UNIVERSITY OF TECHNOLOGY (CUT), GERMANY " &
     &      /T20,"ANY USE REQUIRES WRITTEN LICENSE FROM CUT")')
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
      ALLOCATE(LMX(NAT))
      CALL PDOS$GETR8A('RBAS',3*3,RBAS)
      ALLOCATE(RPOS(3,NAT))
      CALL PDOS$GETR8A('R',3*NAT,RPOS)
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
!     ==  READ PREDEFINED ORBITALS                                            ==
!     ==========================================================================
!
!     __DETERMINE ANGULAR MOMENTUM WEIGHTS AND SPIN CONTRIBUTIONS AND WRITE_____
!     __TO PROTOCOLL FILE.
      ALLOCATE(SPINDIR(3,NAT))
      CALL REPORT(NFILO)
                            CALL TRACE$PASS('AFTER REPORT')
!
!     ==========================================================================
!     ==  READ PREDEFINED ORBITALS                                            ==
!     ==========================================================================
      CALL READCNTL$ORBITAL(LENG,NAT,LMX,RBAS,RPOS)
                            CALL TRACE$PASS('AFTER READCNTL$ORBITAL')
!
!     ==========================================================================
!     ==  SELECT MATRIXELEMENTS                                               ==
!     ==========================================================================
      NBB=NB    ! #(STATES PER K-POINT) NOS SPIN DEGERACY
      IF(NDIM.EQ.1)NBB=2*NB
      CALL READCNTL$SETNUMBER(NSET)
      ALLOCATE(LEGEND(NSET))
      ALLOCATE(SET(NBB,NKPT,2,NSET))  ! SET WORKS ALWAYS WITH TWO SPINS
      ALLOCATE(EIG(NBB,NKPT))
      ALLOCATE(OCC(NBB,NKPT))
!     __ FILL IN EIGENVALUES AND OCCUPATIONS____________________________________
      CALL SET$ENOCC(NBB,NKPT,EIG,OCC)
!     __ FILL IN WEIGHT FOR EACH STATE__________________________________________
      CALL READCNTL$SETS(NBB,NKPT,NSET,NAT,LMX,RBAS,RPOS,LENG,SET,LEGEND)
      DEALLOCATE(LMX)
                            CALL TRACE$PASS('AFTER READCNTL$SETS')
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
     &                    ,NBB,NKPT,EIG,OCC,NSET,SET,LEGEND,MODE)
                            CALL TRACE$PASS('AFTER READCNTL$OUTPUT')
!
!     ==========================================================================
!     ==  CLOSING                                                             ==
!     ==========================================================================
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
        WRITE(NFIL,*)''
        DO IAT=1,NAT
          WRITE(NFIL,FMT='(A2,F10.5,F10.5,F10.5)')ATOMID(IAT),R(:,IAT)
        END DO
        WRITE(NFIL,*)''
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
      CALL LIB$NARGS(NARGS)
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE PDOS TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL LIB$GETARG(1,PDOSINNAME)
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
      LOGICAL(4),PARAMETER :: T=.TRUE.
      LOGICAL(4),PARAMETER :: F=.FALSE.
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
      CALL FILEHANDLER$SETFILE(ID,T,-'.DERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOLL FILE========================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.DPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =====================================================
      ID=+'DCNTL'
      CALL FILEHANDLER$SETFILE(ID,T,-'.DCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  STRUCTURE FILE   =====================================================
      ID=+'PDOS'
      CALL FILEHANDLER$SETFILE(ID,T,-'.PDOS')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  DENSITY OF STATES FILE PRODUCES AS OUTPUT ============================
!     ==  WILL BE ATTACHED TO DIFFERENT FILES DURING EXECUTION =================
      ID=+'PDOSOUT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.PDOSOUT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  NUMBER OF STATES FILE PRODUCES AS OUTPUT =============================
!     ==  WILL BE ATTACHED TO DIFFERENT FILES DURING EXECUTION =================
      ID=+'PNOSOUT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.PNOSOUT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  SPIN GRAPHICS FILE   =================================================
      ID=+'MOL'
      CALL FILEHANDLER$SETFILE(ID,T,-'.MOL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE ORBITALS_MODULE
PRIVATE
PUBLIC ORBITALS$SETORB
PUBLIC ORBITALS$GETORB
LOGICAL(4)                :: TINI=.FALSE.
INTEGER(4)                :: NORB=0      ! #(ORBITALS)
INTEGER(4)                :: NORBX=0     ! MAX#(ORBITALS)
INTEGER(4)    ,PARAMETER  :: NORBSTEP=10 ! STEP IN #(ORBITALS)
INTEGER(4)                :: LENG=0      ! LENGTH OF ORBITALVECTOR
COMPLEX(8)   ,ALLOCATABLE :: ORBITAL(:,:)
CHARACTER(21),ALLOCATABLE :: ORBITALNAME(:)
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ORBITALS$SETORB(NAME_,LENG_,ORBITAL_)
      CHARACTER(*),INTENT(IN) :: NAME_
      INTEGER(4)  ,INTENT(IN) :: LENG_
      COMPLEX(8)  ,INTENT(IN) :: ORBITAL_(LENG_)
!     **************************************************************************
      IF(.NOT.TINI) THEN
        LENG=LENG_
        TINI=.TRUE.
      END IF
      IF(LENG.NE.LENG_) THEN
        CALL ERROR$MSG('LENGTH OF ORBITAL VECTOR INCONSISTENT')
        CALL ERROR$I4VAL('LENG_',LENG_)
        CALL ERROR$I4VAL('LENG',LENG)
        CALL ERROR$STOP('ORBITALS$SETORB')
      END IF
      CALL RESIZE
      NORB=NORB+1
      ORBITAL(:,NORB)=ORBITAL_(:)
      ORBITALNAME(NORB)=NAME_
      RETURN
      END SUBROUTINE ORBITALS$SETORB
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ORBITALS$GETORB(NAME_,LENG_,ORBITAL_)
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: NAME_
      INTEGER(4)   ,INTENT(IN) :: LENG_
      COMPLEX(8)   ,INTENT(OUT):: ORBITAL_(LENG_)
      INTEGER(4)               :: IORB
!     **************************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('ORBITALS MODULE NOT YET INITIALIZED')
        CALL ERROR$STOP('ORBITALS$GETORB')
      END IF
      IF(LENG.NE.LENG_) THEN
        CALL ERROR$MSG('LENGTH OF ORBITAL VECTOR INCONSISTENT')
        CALL ERROR$I4VAL('LENG_',LENG_)
        CALL ERROR$I4VAL('LENG',LENG)
        CALL ERROR$STOP('ORBITALS$GETORB')
      END IF
      DO IORB=1,NORB
        IF(TRIM(ORBITALNAME(IORB)).EQ.TRIM(NAME_)) THEN
          ORBITAL_(:)=ORBITAL(:,IORB)
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ORBITAL NAME NOT RECOGNIZED')
      CALL ERROR$CHVAL('NAME_',NAME_)
      CALL ERROR$STOP('ORBITALS$GETORB')
      RETURN
      END SUBROUTINE ORBITALS$GETORB
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RESIZE
      COMPLEX(8)  ,ALLOCATABLE :: TMPORBITAL(:,:)
      CHARACTER(32),ALLOCATABLE :: TMPNAME(:)
!     **************************************************************************
      IF(NORB+1.LT.NORBX) RETURN
      IF(NORB.GT.0) THEN
        ALLOCATE(TMPORBITAL(LENG,NORB))
        ALLOCATE(TMPNAME(NORB))
        TMPORBITAL(:,:)=ORBITAL(:,1:NORB)
        TMPNAME(:)     =ORBITALNAME(1:NORB)
        DEALLOCATE(ORBITAL)
        DEALLOCATE(ORBITALNAME)
      END IF
      NORBX=NORB+NORBSTEP
      ALLOCATE(ORBITAL(LENG,NORBX))
      ALLOCATE(ORBITALNAME(NORBX))
      IF(NORB.GT.0) THEN
        ORBITAL(:,1:NORB)=TMPORBITAL(:,:)
        ORBITALNAME(1:NORB)=TMPNAME(:)
        DEALLOCATE(TMPORBITAL)
        DEALLOCATE(TMPNAME)
      END IF
      ORBITAL(:,NORB+1:NORBX)=0.D0
      ORBITALNAME(NORB+1:NORBX)=' '
      END SUBROUTINE RESIZE
END MODULE ORBITALS_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READONEORB(LL_CNTL,NAT,LMX,RBAS,RPOS,NPRO,ORBITAL)
!     **************************************************************************
!     ** READ  AN ORBITAL BLOCK FROM LIST LL_CNTL AND RETURN                  **
!     ** A VECTOR DEFINING THAT ORBITAL                                       **
!     **                                                                      **
!     **************************************************************************
      USE ORBITALS_MODULE, ONLY : ORBITALS$GETORB
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NAT
      INTEGER(4)   ,INTENT(IN) :: LMX(NAT)
      REAL(8)      ,INTENT(IN) :: RBAS(3,3)
      REAL(8)      ,INTENT(IN) :: RPOS(3,NAT)
      INTEGER(4)   ,INTENT(IN) :: NPRO
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL
      COMPLEX(8)   ,INTENT(OUT):: ORBITAL(NPRO)
      INTEGER(4)   ,PARAMETER  :: LMXX=16
      REAL(8)                  :: ORB(LMXX)
      INTEGER(4)               :: IAT,IAT2
      CHARACTER(32)            :: ATOM1,ATOMZ,ATOMX
      CHARACTER(32)            :: ORBITALNAME1 
      CHARACTER(8)             :: TYPE
      INTEGER(4)               :: IT3(3),IT3Z(3),IT3X(3)
      COMPLEX(8)               :: CFAC
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: DRZ(3)
      REAL(8)                  :: DRX(3)
      REAL(8)                  :: ROT(3,3)
      REAL(8)      ,ALLOCATABLE:: YLMROT(:,:)
!     **************************************************************************
      ORBITAL(:)=0.D0
!
!     ==========================================================================
!     ==  GET PREFACTOR                                                       ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FAC',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$CONVERT(LL_CNTL,'FAC',1,'C(8)')
        CALL LINKEDLIST$GET(LL_CNTL,'FAC',1,CFAC)
      ELSE
        CFAC=(1.D0,0.D0)
      END IF
!
!     ==========================================================================
!     ==  SEARCH PREDEFINED ORBITALS                                          ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,ORBITALNAME1)
        CALL ORBITALS$GETORB(ORBITALNAME1,NPRO,ORBITAL)
      END IF
!
!     ==========================================================================
!     ==  BUILD NEW ORBITAL COMPONENT FROM BASIC BUILDING BLOCKS              ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK)
      IF(TCHK) THEN   
        CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,ATOM1)
        CALL RESOLVEATOM(ATOM1,IAT,IT3)
        CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
        TYPE=+TYPE
        CALL RESOLVETYPE(LMXX,TYPE,ORB)
!       
!       ========================================================================
!       ==  FIND NEAREST NEIGHBOUR DIRECTIONS                                 ==
!       ========================================================================
        DRZ(:)=0.D0
        DRZ(3)=1.D0
        DRX(:)=0.D0
        DRX(1)=1.D0
        CALL PDOS$GETR8A('RBAS',9,RBAS) !ADDED
        CALL LINKEDLIST$EXISTD(LL_CNTL,'Z',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'Z',1,DRZ(:))
        CALL LINKEDLIST$EXISTD(LL_CNTL,'X',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'X',1,DRX(:))
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NNZ',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NNZ',1,ATOMZ)
          CALL RESOLVEATOM(ATOMZ,IAT2,IT3Z)
          DRZ(:)=RPOS(:,IAT2)-RPOS(:,IAT) &
       &        +RBAS(:,1)*REAL(IT3Z(1)-IT3(1),KIND=8) &
       &        +RBAS(:,2)*REAL(IT3Z(2)-IT3(2),KIND=8) &
       &        +RBAS(:,3)*REAL(IT3Z(3)-IT3(3),KIND=8)
          CALL LINKEDLIST$SET(LL_CNTL,'Z',0,DRZ(:))
        END IF
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NNX',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NNX',1,ATOMX)
          CALL RESOLVEATOM(ATOMX,IAT2,IT3X)
          DRX(:)=RPOS(:,IAT2)-RPOS(:,IAT) &
       &        +RBAS(:,1)*REAL(IT3X(1)-IT3(1),KIND=8) &
       &        +RBAS(:,2)*REAL(IT3X(2)-IT3(2),KIND=8) &
       &        +RBAS(:,3)*REAL(IT3X(3)-IT3(3),KIND=8)
          CALL LINKEDLIST$SET(LL_CNTL,'X',0,DRX(:))
        END IF
!       
!       ========================================================================
!       ==  MAKE ORBITAL                                                      ==
!       ========================================================================
        CALL RESOLVEROTATION(DRZ,DRX,ROT)
        ALLOCATE(YLMROT(LMXX,LMXX))
        CALL ROTATEYLM(LMXX,ROT,YLMROT)
        ORB=MATMUL(YLMROT,ORB)
        DEALLOCATE(YLMROT)
        CALL MAKEORBITAL(ATOM1,LMXX,ORB,NPRO,ORBITAL)
        ORBITAL=ORBITAL*CFAC
      END IF
      RETURN
      CONTAINS
!
!       .1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE RESOLVETYPE(LMX,TYPE,ORBITAL)
!       ************************************************************************
!       ** CONSTRUCTS THE PREFACTORS OF AN ORBITAL IN AN EXPANSION OF         **
!       ** REAL SPHERICAL HARMONICS                                           **
!       ************************************************************************
        IMPLICIT NONE
        INTEGER(4)  ,INTENT(IN)   :: LMX
        CHARACTER(*),INTENT(IN)   :: TYPE
        REAL(8)     ,INTENT(OUT)  :: ORBITAL(LMX)
        REAL(8)                   :: ORB(9)
        INTEGER(4)                :: LM
!       ************************************************************************
        ORB(:)=0.D0
!
!       ========================================================================
!       ==  SET ORBITAL COEFFICIENTS                                          ==
!       ========================================================================
        IF(TRIM(TYPE).EQ.'S') THEN
          ORB(1)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'PX') THEN
          ORB(2)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'PZ') THEN
          ORB(3)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'PY') THEN
          ORB(4)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'SP1') THEN
          ORB(1)=SQRT(1.D0/2.D0)
          ORB(3)=SQRT(1.D0/2.D0)
        ELSE IF(TRIM(TYPE).EQ.'SP2') THEN
          ORB(1)=SQRT(1.D0/3.D0)
          ORB(3)=SQRT(2.D0/3.D0)
        ELSE IF(TRIM(TYPE).EQ.'SP3') THEN
          ORB(1)=SQRT(1.D0/4.D0)
          ORB(3)=SQRT(3.D0/4.D0)
        ELSE IF(TRIM(TYPE).EQ.'DX2-Y2') THEN
          ORB(5)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'DXZ') THEN
          ORB(6)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'D3Z2-R2') THEN
          ORB(7)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'DYZ') THEN
          ORB(8)=1.D0
        ELSE IF(TRIM(TYPE).EQ.'DXY') THEN
          ORB(9)=1.D0
        ELSE
          CALL ERROR$MSG('TYPE NOT IDENTIFIED')
          CALL ERROR$CHVAL('TYPE',TRIM(TYPE))
          CALL ERROR$STOP('RESOLVETYPE')
        END IF
        DO LM=LMX+1,9
          IF(ORB(LM).NE.0) THEN
            CALL ERROR$MSG('LMX TOO SMALL')
            CALL ERROR$CHVAL('TYPE',TRIM(TYPE))
            CALL ERROR$I4VAL('LMX',LMX)
            CALL ERROR$STOP('RESOLVETYPE')
          END IF
        ENDDO
        ORBITAL(1:9)=ORB
        ORBITAL(10:)=0.D0
        RETURN
        END SUBROUTINE RESOLVETYPE
      END SUBROUTINE READONEORB
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEORBITAL(ATOM,LMXX,ORB,NPRO_,ORBITAL)
!     **************************************************************************
!     **  ONLY THE FIRST PARTIAL WAVE PER ANGULAR MOMENTUM IS                 **
!     **  CONSIDERED                                                          **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ATOM
      INTEGER(4)  ,INTENT(IN) :: LMXX
      REAL(8)     ,INTENT(IN) :: ORB(LMXX)
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
      IF(DXLEN.EQ.0.D0) THEN
        DX=(/0.D0,1.D0,0.D0/)
        DX=DX-DZ*DOT_PRODUCT(DZ,DX) !JO AB DA
        DXLEN=SQRT(DX(1)**2+DX(2)**2+DX(3)**2) 
        IF(DXLEN.EQ.0.D0) THEN
          DX=(/1.D0,0.D0,0.D0/)
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
      USE PDOS_MODULE
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
      USE READCNTL_MODULE
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
      PREFIX=''
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
          CALL ERROR$STOP('PDOS MAIN')
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
      SUBROUTINE READCNTL$ORBITAL(NPRO,NAT,LMX,RBAS,RPOS)
!     **************************************************************************
!     **************************************************************************
      USE READCNTL_MODULE
      USE ORBITALS_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NPRO
      INTEGER(4)   ,INTENT(IN) :: NAT
      INTEGER(4)   ,INTENT(IN) :: LMX(NAT)
      REAL(8)      ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)      ,INTENT(IN) :: RPOS(3,NAT)  !ATOMIC POSITIONS
      CHARACTER(32)            :: ORBITALNAME
      COMPLEX(8)               :: ORBITAL(NPRO)
      COMPLEX(8)               :: ORBITAL1(NPRO)
      INTEGER(4)               :: IORB,ITH
      INTEGER(4)               :: NUM
      INTEGER(4)               :: NORB
!     **************************************************************************
                          CALL TRACE$PUSH('READCNTL$ORBITAL')
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'ORBITAL',NORB)
      DO IORB=1,NORB
        CALL LINKEDLIST$SELECT(LL_CNTL,'ORBITAL',IORB)
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,ORBITALNAME)
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB',NUM)
        ORBITAL(:)=0.D0
        DO ITH=1,NUM
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB',ITH)
          CALL READONEORB(LL_CNTL,NAT,LMX,RBAS,RPOS,NPRO,ORBITAL1)
          ORBITAL(:)=ORBITAL(:)+ORBITAL1(:)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL NORMALIZEORBITAL(NPRO,ORBITAL)
        CALL ORBITALS$SETORB(ORBITALNAME,NPRO,ORBITAL)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL$ORBITAL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$SETNUMBER(NSET)
!     **************************************************************************
!     **************************************************************************
      USE READCNTL_MODULE
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
      SUBROUTINE READCNTL$SETS(NBB,NKPT,NSET &
     &                        ,NAT,LMX,RBAS,RPOS,LENG,SET,LEGEND)
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
      USE READCNTL_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)  :: NBB
      INTEGER(4)   ,INTENT(IN)  :: NKPT
      INTEGER(4)   ,INTENT(IN)  :: NSET
      INTEGER(4)   ,INTENT(IN)  :: NAT
      INTEGER(4)   ,INTENT(IN)  :: LMX(NAT)
      INTEGER(4)   ,INTENT(IN)  :: LENG
      REAL(8)      ,INTENT(IN)  :: RBAS(3,3)
      REAL(8)      ,INTENT(IN)  :: RPOS(3,NAT)
      REAL(8)      ,INTENT(OUT) :: SET(NBB,NKPT,2,NSET) ! ALWAYS WITH 2 SPINS 
      CHARACTER(32),INTENT(OUT) :: LEGEND(NSET)
      REAL(8)      ,ALLOCATABLE :: SET1(:,:,:) !(NBB,NKPT,2)
      INTEGER(4)                :: ISET   ! SET COUNTER
      INTEGER(4)                :: ISPIN   ! SET COUNTER
      CHARACTER(8)              :: TYPE
      CHARACTER(32)             :: NAME,SPIN
      LOGICAL(4)                :: TCHK,TCHK1
      INTEGER(4)                :: ITH, NUM
      INTEGER(4)                :: JTH, NUMJTH
      INTEGER(4)                :: IORB1,IORB2,IORB
      INTEGER(4)                :: NORB1,NORB2,NORB
      COMPLEX(8)                :: ORBITAL1(LENG)
      COMPLEX(8)                :: ORBITAL2(LENG)
      COMPLEX(8)                :: ORBITALI(LENG)
!     **************************************************************************
                          CALL TRACE$PUSH('READCNTL$SETS')
      SET(:,:,:,:)=0.D0
      ISET=0
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
!
!       == LOOK UP ORBITALS ====================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB1',NORB1)
        ORBITAL1=0.D0 
        DO IORB1=1,NORB1
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB1',IORB1) 
          CALL READONEORB(LL_CNTL,NAT,LMX,RBAS,RPOS,LENG,ORBITALI)
          ORBITAL1(:)=ORBITAL1(:)+ORBITALI(:)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB2',NORB2)
        ORBITAL2=0.D0 
        DO IORB2=1,NORB2
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB2',IORB2)
          CALL READONEORB(LL_CNTL,NAT,LMX,RBAS,RPOS,LENG,ORBITALI)
          ORBITAL2(:)=ORBITAL2(:)+ORBITALI(:)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',1,TCHK1)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,SPIN)
          SPIN=+SPIN
        ELSE
!         ==  SPIN=TOTAL MEANS "SUM OVER ALL SPINOR COMPONENTS". =============
!         ==  FOR COLLINEAR CALCULATIONS, SPIN='TOTAL' PRODUCES THE TWO ======
!         ==  SPIN DENSITIES. ================================================
          SPIN='TOTAL' 
        END IF
        CALL SET$PROJECT(LENG,NBB,NKPT,ORBITAL1,ORBITAL2,SPIN,SET(:,:,:,ISET))
        WRITE(LEGEND(ISET),FMT='("COOP:SET",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ID',0,LEGEND(ISET))
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,LEGEND(ISET))
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
        TYPE=''
        SET(:,:,:,ISET)=0.D0
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',1,TCHK1)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,SPIN)
          SPIN=+SPIN
        ELSE
!         ==  SPIN=TOTAL MEANS "SUM OVER ALL SPINOR COMPONENTS". =============
!         ==  FOR COLLINEAR CALCULATIONS, SPIN='Z' PRODUCES THE TWO ==========
!         ==  SPIN DENSITIES. ================================================
          SPIN='Z' 
        END IF
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
            CALL ERROR$STOP('READCNTL$SETS')
          END IF
!
!         ======================================================================
!         ==  'TOTAL' = TOTAL DENSITY OF STATES                               ==
!         ======================================================================
          IF(TRIM(TYPE).EQ.'TOTAL') THEN
            CALL SET$WEIGHT('TOTAL','',NBB,NKPT,SPIN,SET(1,1,1,ISET))
!
!         ======================================================================
!         ==  'ALL' = ALL PROJECTED DENSITY OF STATES                         ==
!         ======================================================================
          ELSE IF(TRIM(TYPE).EQ.'ALL') THEN
            CALL SET$WEIGHT('ALL','',NBB,NKPT,SPIN,SET(1,1,1,ISET))
!
!         ======================================================================
!         ==  'EMPTY' = VACCUM DENSITY OF STATES                              ==
!         ======================================================================
          ELSE IF(TRIM(TYPE).EQ.'EMPTY') THEN
            ALLOCATE(SET1(NBB,NKPT,2))
            SET1=0.D0
            CALL SET$WEIGHT('TOTAL','',NBB,NKPT,SPIN,SET1)
            SET(:,:,:,ISET)=SET(:,:,:,ISET)+SET1
            SET1=0.D0
            CALL SET$WEIGHT('ALL','',NBB,NKPT,SPIN,SET1)
            SET(:,:,:,ISET)=SET(:,:,:,ISET)-SET1
            DEALLOCATE(SET1)

!         ======================================================================
!         ==  NOT RECOGNIZED                                                  ==
!         ======================================================================
          ELSE 
            CALL ERROR$MSG('TYPE NOT RECOGNIZED')
            CALL ERROR$MSG('MUST BE EITHER TOTAL,ALL OR EMPTY')
            CALL ERROR$CHVAL('TYPE ',TYPE)
            CALL ERROR$STOP('READCNTL$SETS')
          END IF
        END IF
!
!       ========================================================================
!       == COLLECT CONTRIBUTIONS FROM INDIVIDUAL ATOMS                        ==
!       ========================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ATOM',NUMJTH)
        DO JTH=1,NUMJTH
          CALL LINKEDLIST$SELECT(LL_CNTL,'ATOM',JTH)
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,NAME)
          CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',1,TCHK)
!
!         == SELECT TYPE =======================================================
          CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
          TYPE=+TYPE
          IF(TYPE.EQ.'ALL') THEN
            CALL SET$WEIGHT(NAME,'ALL',NBB,NKPT,SPIN,SET(1,1,1,ISET))
          ELSE
            IF(SCAN(TYPE,'S').NE.0) THEN
              CALL SET$WEIGHT(NAME,'S',NBB,NKPT,SPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'P').NE.0) THEN
              CALL SET$WEIGHT(NAME,'P',NBB,NKPT,SPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'D').NE.0) THEN
              CALL SET$WEIGHT(NAME,'D',NBB,NKPT,SPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'F').NE.0) THEN
              CALL SET$WEIGHT(NAME,'F',NBB,NKPT,SPIN,SET(1,1,1,ISET))
            END IF

            IF(SCAN(TYPE,'ABCEGHIJKLMNOQRTUVWXYZ0123456789').NE.0) THEN
              CALL ERROR$MSG('ALLOWED VALUES FOR VARIABLE TYPE ARE')
              CALL ERROR$MSG('ONLY "ALL" OR A COMBINATION OF "S", "P", "D","F"')
              CALL ERROR$CHVAL('ATOM',NAME)
              CALL ERROR$CHVAL('TYPE',TYPE)
              CALL ERROR$STOP('READCNTL$SETS')
            END IF
          END IF
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
!       ========================================================================
!       == ADD CONTRIBUTION FROM PREDEFINED ORBITALS ===========================
!       ========================================================================
!
!       == LOOK UP ORBITALS ====================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB',NORB)
        DO IORB=1,NORB
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB',IORB)
          CALL READONEORB(LL_CNTL,NAT,LMX,RBAS,RPOS,LENG,ORBITALI)
          CALL SET$PROJECT(LENG,NBB,NKPT,ORBITALI,ORBITALI,SPIN,SET(1,1,1,ISET))
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        WRITE(LEGEND(ISET),FMT='("WEIGHT",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ID',0,LEGEND(ISET))
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,LEGEND(ISET))
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL$SETS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SET$WEIGHT(ATOMID_,ORBITALID,NBB,NKPT_,SPIN,SET)
!     **************************************************************************
!     ** ADDS A CONTRIBUTION TO THE SET, WHICH CONTAINS THE CONTRIBUTION      **
!     ** TO THE WEIGHT FROM ALL THE STATES                                    **
!     **   SPIN: MAY BE 'TOTAL', 'MAIN', 'X', 'Y', 'Z'                        **
!     **         TOTAL IS THE TOTAL DENSITY OF STATES                         **
!     **         MAIN IS THE PROJECTION ON SPINDIR(I,IAT)                     **
!     **   ORBITALID MAY BE:  'TOTAL','ALL','S','P','D','F'                   **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE   , ONLY : NDIM &
     &                          ,NSPIN &
     &                          ,NKPT &
     &                          ,STATEARR,STATE &
     &                          ,NAT &
     &                          ,ATOMID &
     &                          ,ISPECIES &
     &                          ,LNX &
     &                          ,LOX &
     &                          ,OV
      USE SPINDIR_MODULE, ONLY : SPINDIR
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: ATOMID_
      CHARACTER(*),INTENT(IN)   :: ORBITALID !MAY BE 'ALL','S','P','D','F'
      INTEGER(4)  ,INTENT(IN)   :: NBB
      INTEGER(4)  ,INTENT(IN)   :: NKPT_
      CHARACTER(*),INTENT(IN)   :: SPIN   
      REAL(8)     ,INTENT(INOUT):: SET(NBB,NKPT_,2) !ALWAYS TWO SPIN COMPONENTS
      INTEGER(4)                :: ISPIN,IKPT,IB,IDIM
      INTEGER(4)                :: IPRO0,IPRO1,IPRO2,IP1,IP2
      INTEGER(4)                :: IAT,IAT0,ISP
      INTEGER(4)                :: L,L1,L2,M,LN,LN1,LN2
      INTEGER(4)                :: IBB
      INTEGER(4)                :: NB
      REAL(8)                   :: SUM,SUM_(3)
      COMPLEX(8),PARAMETER      :: CI=(0.D0,1.D0)
      COMPLEX(8),PARAMETER      :: CONE=(1.D0,0.D0)
      COMPLEX(8)                :: CSVAR1P,CSVAR1M,CSVAR2P,CSVAR2M 
      REAL(8)                   :: SVARP,SVARM      
      COMPLEX(8)                :: CP(2),CM(2) ! SPIN EIGENSTATES
      REAL(8)                   :: SVAR
!     **************************************************************************
                                 CALL TRACE$PUSH('SET$WEIGHT')
!
!     ==========================================================================
!     ==  DETERMINE SPINOR EIGENSTATES                                        ==
!     ==========================================================================
      IF(TRIM(SPIN).EQ.'X') THEN
        CALL SPINBRA(1.D0,0.D0,0.D0,CP,CM)
      ELSE IF(TRIM(SPIN).EQ.'Y') THEN
        CALL SPINBRA(0.D0,1.D0,0.D0,CP,CM)
      ELSE IF(TRIM(SPIN).EQ.'Z'.OR.TRIM(SPIN).EQ.'TOTAL') THEN
        CALL SPINBRA(0.D0,0.D0,1.D0,CP,CM)
      ELSE IF(TRIM(SPIN).EQ.'MAIN') THEN
        CALL SPINBRA(SPINDIR(1,IAT0),SPINDIR(2,IAT0),SPINDIR(3,IAT0),CP,CM)
      ELSE
        CALL ERROR$MSG('SPIN COMPONENT NOT RECOGNIZED')
        CALL ERROR$MSG('SPIN SHOULD BE TOTAL,X,Y,Z OR MAIN')
        CALL ERROR$CHVAL('SPIN',SPIN)
        CALL ERROR$STOP('SET$WEIGHT')
      END IF

      IF(ATOMID_.EQ.'TOTAL') THEN
!        SVAR=2.D0/REAL(NDIM*NSPIN,KIND=8)
        svar=1.d0
        IF(TRIM(SPIN).EQ.'TOTAL') THEN
!         == WEIGHT OF BOTH SPIN COMPONENTS IS ADDED TO THE FIRST SPIN CHANNEL
          SET(:,:,1)=SET(:,:,1)+SVAR
          RETURN
        ELSE IF(TRIM(SPIN).EQ.'Z'.AND.NSPIN.EQ.2) THEN
          DO IKPT=1,NKPT
            DO ISPIN=1,NSPIN
              STATE=>STATEARR(IKPT,ISPIN)
              NB=STATE%NB
              SET(:NB  ,IKPT,1)=SET(:NB  ,IKPT,1)+0.5D0*SVAR
              SET(NB+1:,IKPT,2)=SET(NB+1:,IKPT,2)+0.5D0*SVAR
            ENDDO
          ENDDO
        ELSE IF(TRIM(SPIN).EQ.'MAIN') THEN
          CALL ERROR$MSG('SPIN COMPONENT NOT RECOGNIZED')
          CALL ERROR$MSG('SPIN=MAIN INCOMPATIBLE WITH ATOMID=TOTAL')
          CALL ERROR$CHVAL('SPIN',SPIN)
          CALL ERROR$CHVAL('ATOMID',ATOMID)
          CALL ERROR$STOP('SET$WEIGHT')
        ELSE
!         == FOR NON-COLLINEAR CALCULATION THERE IS NO PREFERRED SPIN AXIS =====
!         == THEREFORE THE WEIGHT IS EQUALLY DISTRIBUTED TO BOTH SPIN DIRECTIONS
          SET(:,:,:)=SET(:,:,:)+0.5D0*SVAR
        END IF
        RETURN
      END IF
!
!     ==========================================================================
!     ==  ATOM SPECIFIER IAT0: IAT0<0 IMPLIES ALL ATOMS                       ==
!     ==========================================================================
      IF(ATOMID_.EQ.'ALL'.OR.ATOMID_.EQ.'') THEN
        IAT0=-1
      ELSE
        IAT0=-1
        DO IAT=1,NAT
          IF(ATOMID_.EQ.ATOMID(IAT)) THEN
            IAT0=IAT
            EXIT
          END IF
        ENDDO
        IF(IAT0.EQ.-1) THEN
          CALL ERROR$MSG('ATOMID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ATOMID_',ATOMID_)
          CALL ERROR$STOP('SET$WEIGHT')
        END IF
      END IF
!
      IF(TRIM(SPIN).EQ.'MAIN'.AND.IAT0.EQ.-1.AND.NAT.GT.1) THEN
        CALL ERROR$MSG('MAIN SPIN FOR ALL ATOMS NOT POSSIBLE')
        CALL ERROR$CHVAL('SPIN',SPIN)
        CALL ERROR$CHVAL('ATOMID_',ATOMID_)
        CALL ERROR$STOP('SET$WEIGHT')
      END IF
!
!     ==========================================================================
!     ==  ANGULAR MOMENTUM SPECIFIER L: L<0 MEANS ALL ANGULAR MOMENTA         ==
!     ==========================================================================
      IF(ORBITALID.EQ.'S') THEN
        L=0
      ELSE IF(ORBITALID.EQ.'P') THEN
        L=1
      ELSE IF(ORBITALID.EQ.'D') THEN
        L=2
      ELSE IF(ORBITALID.EQ.'F') THEN
        L=3
      ELSE IF(ORBITALID.EQ.'ALL'.OR.ORBITALID.EQ.'') THEN
        L=-1
      ELSE
        CALL ERROR$MSG('ORBITALID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ORBITALID',ORBITALID)
        CALL ERROR$STOP('SET$WEIGHT')
      END IF
!
!     ==========================================================================
!     ==  NOW EVALUATE WEIGHTS                                                ==
!     ==========================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          IPRO0=0
          DO IAT=1,NAT
            ISP=ISPECIES(IAT)
            IF(IAT.EQ.IAT0.OR.IAT0.LE.0) THEN ! ALL ATOMS FOR IAT0<1
              IPRO1=IPRO0
              DO LN1=1,LNX(ISP)
                L1=LOX(LN1,ISP)
                IPRO2=IPRO0
                DO LN2=1,LNX(ISP)
                  L2=LOX(LN2,ISP)
!                 == SELECT ANGULAR MOMENTA (ALL FOR L<0)
                  IF(L1.NE.L2.OR.(L.GE.0.AND.L1.NE.L)) THEN  
                    IPRO2=IPRO2+2*L2+1
                    CYCLE
                  END IF
                  NB=STATE%NB
                  DO IB=1,NB
                    IBB=IB+NB*(ISPIN-1)
                    DO M=1,2*L1+1
                      IP1=IPRO1+M
                      IP2=IPRO2+M
                      IF(NDIM.EQ.1) THEN
                        CSVAR1P=CP(ISPIN)*STATE%VEC(1,IP1,IB)
                        CSVAR2P=CP(ISPIN)*STATE%VEC(1,IP2,IB)
                        CSVAR1M=CM(ISPIN)*STATE%VEC(1,IP1,IB)
                        CSVAR2M=CM(ISPIN)*STATE%VEC(1,IP2,IB)
                        SVARP=REAL(CONJG(CSVAR1P)*CSVAR2P,KIND=8)
                        SVARM=REAL(CONJG(CSVAR1M)*CSVAR2M,KIND=8)
                        IF(SPIN.EQ.'TOTAL') THEN
                          SVARP=SVARP+SVARM
                          SVARM=0.D0
                        END IF
                        SET(IBB,IKPT,1)=SET(IBB,IKPT,1)+SVARP*OV(LN1,LN2,ISP)
                        SET(IBB,IKPT,2)=SET(IBB,IKPT,2)+SVARM*OV(LN1,LN2,ISP)
                        IF(ISPIN.EQ.1.AND.NSPIN.EQ.1) THEN
                          CSVAR1P=CP(2)*STATE%VEC(1,IP1,IB)
                          CSVAR2P=CP(2)*STATE%VEC(1,IP2,IB)
                          CSVAR1M=CM(2)*STATE%VEC(1,IP1,IB)
                          CSVAR2M=CM(2)*STATE%VEC(1,IP2,IB)
                          SVARP=REAL(CONJG(CSVAR1P)*CSVAR2P,KIND=8)
                          SVARM=REAL(CONJG(CSVAR1M)*CSVAR2M,KIND=8)
                          IF(SPIN.EQ.'TOTAL') THEN
                            SVARP=SVARP+SVARM
                            SVARM=0.D0
                          END IF
                          SET(IB+NB,IKPT,1)=SET(IB+NB,IKPT,1) &
      &                                                  +SVARP*OV(LN1,LN2,ISP)
                          SET(IB+NB,IKPT,2)=SET(IB+NB,IKPT,2) &
      &                                                  +SVARM*OV(LN1,LN2,ISP)
                        END IF
                      ELSE IF(NDIM.EQ.2) THEN
                        CSVAR1P=(0.D0,0.D0)
                        CSVAR1M=(0.D0,0.D0)
                        CSVAR2P=(0.D0,0.D0)
                        CSVAR2M=(0.D0,0.D0)
                        DO IDIM=1,NDIM
                          CSVAR1P=CSVAR1P+CP(IDIM)*STATE%VEC(IDIM,IP1,IB)
                          CSVAR1M=CSVAR1M+CM(IDIM)*STATE%VEC(IDIM,IP1,IB)
                          CSVAR2P=CSVAR2P+CP(IDIM)*STATE%VEC(IDIM,IP2,IB)
                          CSVAR2M=CSVAR2M+CM(IDIM)*STATE%VEC(IDIM,IP2,IB)
                        ENDDO
                        SVARP=REAL(CONJG(CSVAR1P)*CSVAR2P,KIND=8)
                        SVARM=REAL(CONJG(CSVAR1M)*CSVAR2M,KIND=8)
                        IF(SPIN.EQ.'TOTAL') THEN
                          SVARP=SVARP+SVARM
                          SVARM=0.D0
                        END IF
                        SET(IB,IKPT,1)=SET(IB,IKPT,1)+SVARP*OV(LN1,LN2,ISP)
                        SET(IB,IKPT,2)=SET(IB,IKPT,2)+SVARM*OV(LN1,LN2,ISP)
                      END IF
                    ENDDO
                  ENDDO  ! END OF LOOP OVER BANDS (IB)
                  IPRO2=IPRO2+2*L2+1
                ENDDO    !END OF LOOP OVER LN2
                IPRO1=IPRO1+2*L1+1
              ENDDO      !END OF LOOP OVER LN1
            END IF
            DO LN=1,LNX(ISP)
              IPRO0=IPRO0+2*LOX(LN,ISP)+1
            ENDDO
          ENDDO  ! END LOOP OVER ATOMS (IAT)
        ENDDO  ! END LOOP OVER SPINS (ISPIN)
      ENDDO  ! END LOOP OVER K-POINTS (IKPT)
                                 CALL TRACE$POP
      RETURN
      END SUBROUTINE SET$WEIGHT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINBRA(X,Y,Z,CP,CM)
!     **************************************************************************
!     ** BRAS OF SPIN EIGENSTATES FOR A SPIN-AXIS (X,Y,Z)                     **
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
      CP(1)=1.D0+EZ
      CP(2)=EX+CI*EY
      CM(1)=EX-CI*EY
      CM(2)=-(1.D0+EZ)
      CP=CP/SQRT(2.D0*(1.D0+EZ))
      CM=CM/SQRT(2.D0*(1.D0+EZ))
!
!     == TURN KETS INTO BRAS =================================================
      CP=CONJG(CP)
      CM=CONJG(CM)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SET$PROJECT(NPRO_,NBB,NKPT_,ORBITAL1,ORBITAL2,SPIN,SET)
!     **************************************************************************
!     **  ADDS PROJECTION ONTO MATRIX ELEMENTS OF PREDEFINED ORBITALS         **
!     **  SPECIFIED BY ORBITAL1 AND ORBITAL2 TO SET                           **
!     **                                                                      **
!     **  ASSUMES THAT ORBITAL1 AND 2 ARE NONZERO ONLY FOR THE FIRST PARTIAL  **
!     **    WAVE PER SITE AND ANGULAR MOMENTUM
!     **                                                                      **
!     **  CAUTION: SET IS AN INOUT VARIABLE!                                  **
!     **************************************************************************
      USE PDOS_MODULE, ONLY : NSPIN &
     &                       ,NKPT &
     &                       ,NDIM &
     &                       ,STATEARR,STATE &
     &                       ,NPRO &
     &                       ,NAT &
     &                       ,ISPECIES &
     &                       ,LNX &
     &                       ,LOX &
     &                       ,OV
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: NPRO_
      INTEGER(4)  ,INTENT(IN)   :: NBB
      INTEGER(4)  ,INTENT(IN)   :: NKPT_
      COMPLEX(8)  ,INTENT(IN)   :: ORBITAL1(NPRO)
      COMPLEX(8)  ,INTENT(IN)   :: ORBITAL2(NPRO)
      CHARACTER(*),INTENT(IN)   :: SPIN   
      REAL(8)     ,INTENT(INOUT):: SET(NBB,NKPT_,2) 
      COMPLEX(8)                :: ORBITAL1P(2,NPRO)
      COMPLEX(8)                :: ORBITAL2P(2,NPRO)
      COMPLEX(8)                :: ORBITAL1M(2,NPRO)
      COMPLEX(8)                :: ORBITAL2M(2,NPRO)
      INTEGER(4)                :: ISPIN,IKPT,IB,IPRO,IDIM,IBB
      INTEGER(4)                :: IAT,ISP,LN,L,M
      INTEGER(4)                :: NB
      COMPLEX(8)                :: CSVAR,CSVAR1,CSVAR2
      COMPLEX(8)                :: CVEC(2),CP(2),CM(2)
      REAL(8)                   :: SVAR
!     **************************************************************************
                                 CALL TRACE$PUSH('SET$PROJECT')
!
!     ==========================================================================
!     ==  DETERMINE SPINOR EIGENSTATES                                        ==
!     ==========================================================================
      IF(TRIM(SPIN).EQ.'X') THEN
        CALL SPINBRA(1.D0,0.D0,0.D0,CP,CM)
      ELSE IF(TRIM(SPIN).EQ.'Y') THEN
        CALL SPINBRA(0.D0,1.D0,0.D0,CP,CM)
      ELSE IF(TRIM(SPIN).EQ.'Z'.OR.TRIM(SPIN).EQ.'TOTAL') THEN
        CALL SPINBRA(0.D0,0.D0,1.D0,CP,CM)
      ELSE
        CALL ERROR$MSG('SPIN COMPONENT NOT RECOGNIZED')
        CALL ERROR$MSG('SPIN SHOULD BE TOTAL,X,Y,Z OR MAIN')
        CALL ERROR$CHVAL('SPIN',SPIN)
        CALL ERROR$STOP('SET$PROJECT')
      END IF
!
!     ==========================================================================
!     ==  DETERMINE SPINOR EIGENSTATES                                        ==
!     ==========================================================================
      ORBITAL1P(1,:)=CP(1)*CONJG(ORBITAL1(:))
      ORBITAL1P(2,:)=CP(2)*CONJG(ORBITAL1(:))
      ORBITAL1M(1,:)=CM(1)*CONJG(ORBITAL1(:))
      ORBITAL1M(2,:)=CM(2)*CONJG(ORBITAL1(:))
      ORBITAL2P(1,:)=CP(1)*CONJG(ORBITAL2(:))
      ORBITAL2P(2,:)=CP(2)*CONJG(ORBITAL2(:))
      ORBITAL2M(1,:)=CM(1)*CONJG(ORBITAL2(:))
      ORBITAL2M(2,:)=CM(2)*CONJG(ORBITAL2(:))
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          DO M=1,2*L+1
            IPRO=IPRO+1
            ORBITAL1P(:,IPRO)=ORBITAL1P(:,IPRO)*OV(LN,LN,ISP)
            ORBITAL1M(:,IPRO)=ORBITAL1M(:,IPRO)*OV(LN,LN,ISP)
          ENDDO
        ENDDO
      ENDDO

      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
          NB=STATE%NB
          DO IB=1,NB
            IBB=IB+NB*(ISPIN-1)
            IF(NDIM.EQ.2) THEN
              CSVAR1=SUM(ORBITAL1P(:,:)*STATE%VEC(:,:,IB))
              CSVAR2=SUM(ORBITAL2P(:,:)*STATE%VEC(:,:,IB))
              SET(IBB,IKPT,1)=SET(IBB,IKPT,1)+REAL(CONJG(CSVAR1)*CSVAR2)
              CSVAR1=SUM(ORBITAL1M(:,:)*STATE%VEC(:,:,IB))
              CSVAR2=SUM(ORBITAL2M(:,:)*STATE%VEC(:,:,IB))
              SET(IBB,IKPT,2)=SET(IBB,IKPT,2)+REAL(CONJG(CSVAR1)*CSVAR2)
            ELSE
              CSVAR1=SUM(ORBITAL1P(ISPIN,:)*STATE%VEC(1,:,IB))
              CSVAR2=SUM(ORBITAL2P(ISPIN,:)*STATE%VEC(1,:,IB))
              SET(IBB,IKPT,1)=SET(IBB,IKPT,1)+REAL(CONJG(CSVAR1)*CSVAR2)
              CSVAR1=SUM(ORBITAL1M(ISPIN,:)*STATE%VEC(1,:,IB))
              CSVAR2=SUM(ORBITAL2M(ISPIN,:)*STATE%VEC(1,:,IB))
              SET(IBB,IKPT,2)=SET(IBB,IKPT,2)+REAL(CONJG(CSVAR1)*CSVAR2)
              IF(NSPIN.EQ.1) THEN
                CSVAR1=SUM(ORBITAL1M(2,:)*STATE%VEC(1,:,IB))
                CSVAR2=SUM(ORBITAL2M(2,:)*STATE%VEC(1,:,IB))
                SET(IBB,IKPT,2)=SET(IBB,IKPT,2)+REAL(CONJG(CSVAR1)*CSVAR2)
              END IF
            END IF
          ENDDO   ! END OF LOOP OVER BANDS IB
        ENDDO  ! END OF LOOP OVER K-POINTS IKPT
      ENDDO  ! END OF LOOP OVER SPINS ISPIN
!
!     ==========================================================================
!     == ADD UP CONTRIBUTIONS FOR SPIN=TOTAL                                  ==
!     ==========================================================================
      IF(SPIN.EQ.'TOTAL') THEN
        SET(:,:,1)=SET(:,:,1)+SET(:,:,2)
        SET(:,:,2)=0.D0
      END IF
                                 CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$OUTPUT(EMIN,EMAX,NE,EBROAD,PREFIX &
     &                    ,NBB,NKPT,EIG,OCC,NSET,SET,LEGEND,MODE)
!     **************************************************************************
!     **************************************************************************
      USE STRINGS_MODULE
      USE READCNTL_MODULE
      USE ORBITALS_MODULE
      USE DOS_WGHT_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NE
      REAL(8)      ,INTENT(IN) :: EMIN
      REAL(8)      ,INTENT(IN) :: EMAX
      REAL(8)      ,INTENT(IN) :: EBROAD
      CHARACTER(*) ,INTENT(IN) :: PREFIX ! FOR DOS AND NOS FILES
      INTEGER(4)   ,INTENT(IN) :: NBB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      REAL(8)      ,INTENT(IN) :: EIG(NBB,NKPT)
      REAL(8)      ,INTENT(IN) :: OCC(NBB,NKPT)
      INTEGER(4)   ,INTENT(IN) :: NSET
      REAL(8)      ,INTENT(IN) :: SET(NBB,NKPT,2,NSET)
      CHARACTER(32),INTENT(IN) :: LEGEND(NSET)
      CHARACTER(32),INTENT(IN) :: MODE
      LOGICAL(4)   ,PARAMETER  :: TDOS=.TRUE.  ! DENSITY OF STATES PRINTED ?
      LOGICAL(4)   ,PARAMETER  :: TNOS=.FALSE.  ! NUMBER OF STATES PRINTED ?
      CHARACTER(32)        :: LEGEND1
      CHARACTER(256)       :: FILE
      LOGICAL(4)           :: TCHK
      INTEGER(4)           :: NOUT
      INTEGER(4)           :: NFIL
      INTEGER(4)           :: NFILDOS  ! UNIT FOR DENSITY OF STATES FILE
      INTEGER(4)           :: NFILO
      REAL(8)              :: EV
      REAL(8)              :: E1,E2
      INTEGER(4)           :: IE1,IE2
      LOGICAL(4),ALLOCATABLE :: DEADZONE(:) !(NE,ISPIN)
      INTEGER(4)           :: IB,IKPT,I,ISET
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
        LEGEND1=LEGEND(ISET)
!
!       ========================================================================
!       ==  SPECIFY OUTPUT FILE ================================================
!       ========================================================================
        FILE=TRIM(PREFIX)//TRIM(LEGEND1)//-'.DOS'
        CALL FILEHANDLER$SETFILE('PDOSOUT',.FALSE.,FILE)
        CALL FILEHANDLER$UNIT('PDOSOUT',NFILDOS)
!
!       ========================================================================
!       ==  WRITE DOS AND INTEGRATED DOS ON FILE                              ==
!       ========================================================================
        IF(MODE.EQ.'SAMPLE')THEN
          CALL PUTONGRID_SAMPLE(NFILDOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
      &                 ,NBB,NKPT,EIG,OCC,SET(:,:,:,ISET),LEGEND(ISET))
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
      &                 ,NBB,NKPT,EIG,SET(:,:,:,ISET),LEGEND(ISET))
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
!     **  SCALEY MAY BE 'DOS' OR 'NOS' OR 'NONE'                              **
!     **     SCALEY='DOS' RESCALES THE DENSITY OF STATES TO FIT WINDOW        **
!     **     SCALEY='NOS' RESCALES THE NUMBER OF STATES TO FIT WINDOW         **
!     **  NOS(IE,ISPIN,1) IS MULTIPLIED WITH MAX OCCUPATION OF ALL BANDS      **
!     **  NOS(IE,ISPIN,2) IS MULTIPLIED WITH ACTUAL OCCUPATION OF EACH STATE  **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE, ONLY: STATE,STATEARR,WKPT
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
      REAL(8)                  :: W1,W2,X,FAC
      REAL(8)                  :: WGHTX
      INTEGER(4)               :: IKPT,ISPIN,IE,IB
      REAL(8)                  :: SVAR
      REAL(8)                  :: SIG
      REAL(8)                  :: E
      REAL(8)                  :: SPINDEG
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
      ND=EBROAD/DE * LOG(4.D0/TOL)
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
     &                    ,NBB,NKPT,EIG,SET,LEGEND)
!     **************************************************************************
!     **  MAPS THE CONTRIBUTION FROM EACH STATE ONTO AN ENERGY GRID,          **
!     **  CONSTRUCTS DOS AND NOS WITH THE TETRAHEDRON METHOD AND WRITES THE   **
!     **  RESULT ON FILE                                                      **
!     **                                                                      **
!     **  SCALEY MAY BE 'DOS' OR 'NOS' OR 'NONE'                              **
!     **     SCALEY='DOS' RESCALES THE DENSITY OF STATES TO FIT WINDOW        **
!     **     SCALEY='NOS' RESCALES THE NUMBER OF STATES TO FIT WINDOW         **
!     **  NOS(IE,ISPIN,1) IS MULTIPLIED WITH MAX OCCUPATION OF ALL BANDS      **
!     **  NOS(IE,ISPIN,2) IS MULTIPLIED WITH ACTUAL OCCUPATION OF EACH STATE  **
!     **                                                                      **
!     **************************************************************************
      USE DOS_WGHT_MODULE, ONLY: EF,EWGHT,WGHT
      USE PDOS_MODULE, ONLY: STATE,STATEARR,WKPT
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NE
      REAL(8)      ,INTENT(IN) :: EMIN
      REAL(8)      ,INTENT(IN) :: EMAX
      REAL(8)      ,INTENT(IN) :: EBROAD
      INTEGER(4)   ,INTENT(IN) :: NBB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      REAL(8)      ,INTENT(IN) :: EIG(NBB,NKPT)
      REAL(8)      ,INTENT(IN) :: SET(NBB,NKPT,2)
      INTEGER(4)   ,INTENT(IN) :: NFILDOS ! UNIT FOR DOS FILE OR "-1"
      CHARACTER(32),INTENT(IN) :: LEGEND
      LOGICAL(4)   ,INTENT(IN) :: DEADZONE(NE)
      REAL(8)      ,PARAMETER  :: TOL=1.D-2
      REAL(8)                  :: NOS(NE,2)
      REAL(8)                  :: NOSMIN(2)
      REAL(8)                  :: DOS(NE,2)
      REAL(8)    ,ALLOCATABLE  :: SMEAR(:)
      REAL(8)                  :: DE
      INTEGER(4)               :: IE1,IE2,IDE,I1,I2,IBIS
      INTEGER(4)               :: ND,IOCC
      REAL(8)                  :: EV
      REAL(8)                  :: W1,W2,X
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
!     ==  CALCULATE NUMBER OF STATES FOR EACH ENERGY INTERVAL                 ==
!     ==========================================================================
      DOS(:,:)=0.D0
      NOSMIN(:)=0.D0
      DO ISPIN=1,2   ! THE DATA MODEL HAS BEEN EXPANDED TO 2 COMPONENTS
        DO IKPT=1,NKPT
          DO IB=1,NBB
            SVAR=SET(IB,IKPT,ISPIN)*DE
            NOSMIN(ISPIN)=NOSMIN(ISPIN)+WGHT(IB,IKPT)*SVAR
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
      ND=EBROAD/DE * LOG(4.D0/TOL)
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
          DO IE=1-IDE,IE1-1
            NOSMIN(ISPIN)=NOSMIN(ISPIN)+NOS(IE+IDE,ISPIN)*W1
          ENDDO
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
!     ==  WRITE RESULT ON PSODOUT                                             ==
!     ==========================================================================
!     == WRITE DENSITY OF STATES ===============================================
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
      SUBROUTINE GENERATE_TETRA_WGHT(NFILO,NBB,NKPT,EMAX,EMIN,NE,RBAS,EIG,ELSCALE)
!     **************************************************************************
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
      CALL REPORT$I4VAL(NFILO,"SPACEGROUP",SPACEGROUP,'')
      CALL REPORT$CHVAL(NFILO,"BRAVAIS LATTICE",BRAVAIS)
      CALL REPORT$L4VAL(NFILO,"TIME REVERSAL SYMMETRY EXPLOITED?",TINV)
      CALL REPORT$L4VAL(NFILO,"TSHIFT",TSHIFT)
      CALL REPORT$I4VAL(NFILO,"NUMBER OF K-GRID POINTS ALONG G1",NKDIV(1),'')
      CALL REPORT$I4VAL(NFILO,"NUMBER OF K-GRID POINTS ALONG G2",NKDIV(2),'')
      CALL REPORT$I4VAL(NFILO,"NUMBER OF K-GRID POINTS ALONG G3",NKDIV(3),'')
      CALL REPORT$I4VAL(NFILO,"SHIFT ALONG G1",ISHIFT(1),'')
      CALL REPORT$I4VAL(NFILO,"SHIFT ALONG G2",ISHIFT(2),'')
      CALL REPORT$I4VAL(NFILO,"SHIFT ALONG G3",ISHIFT(3),'')
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
PRINT*,'EFERMI[EV] ',EF*27.211D0
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

