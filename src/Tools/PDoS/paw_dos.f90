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
      REAL(8)      ,ALLOCATABLE :: EIG(:,:,:)
      REAL(8)      ,ALLOCATABLE :: SET(:,:,:,:)
      CHARACTER(32),ALLOCATABLE :: LEGEND(:)
      REAL(8)                   :: EMIN
      REAL(8)                   :: EMAX
      INTEGER(4)                :: NE
      REAL(8)                   :: EBROAD
      LOGICAL(4)                :: TDOS  ! DENSITY OF STATES FILE PRINTED ?
      LOGICAL(4)                :: TNOS  ! NUMBER  OF STATES FILE PRINTED ?
      INTEGER(4)                :: NFILIN
      INTEGER(4)                :: NPRO
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
      ALLOCATE(EIG(NB,NKPT,NSPIN))
      ALLOCATE(SPINDIR(3,NAT))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          DO IB=1,NB
            EIG(IB,IKPT,ISPIN)=STATE%EIG(IB)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  READ GENERAL INFORMATION FROM CONTROL FILE
!     ==========================================================================
      CALL READCNTL$GENERIC(MODE,TDOS,TNOS,PREFIX)
!     == DEFAULT VALUES FOR RANGE OF ENERGY GRID ===============================
      EMIN=MINVAL(EIG)-1.D-1
      EMAX=MINVAL(EIG(NB,:,:)) 
      CALL READCNTL$GRID(EMIN,EMAX,NE,EBROAD)
      CALL READCNTL$REPORT1(MODE,TDOS,TNOS,PREFIX,EMIN,EMAX,NE,EBROAD)
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
      CALL REPORT(NFILO,EIG)
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
      CALL READCNTL$SETNUMBER(NSET)
      ALLOCATE(SET(NB,NKPT,NSPIN,NSET))
      ALLOCATE(LEGEND(NSET))
      CALL READCNTL$SETS(NB,NKPT,NSPIN,NSET,NAT,LMX,RBAS,RPOS,LENG,SET,LEGEND)
      DEALLOCATE(LMX)
                            CALL TRACE$PASS('AFTER READCNTL$SETS')
!
!     ==========================================================================
!     ==  MAKE PLOTS                                                          ==
!     ==========================================================================
!     ==  CALCULAT WEIGHTS FOR DOS USING THE TETRAHEDRON METHOD               ==
      IF(MODE.EQ.'TETRA')THEN
        CALL GENERATE_TETRA_WGHT(NFILO,NB,NSPIN,NKPT,EMAX,EMIN,NE,RBAS,EIG)
      ENDIF
!     == WRITE FILES ===========================================================
      CALL READCNTL$OUTPUT(EMIN,EMAX,NE,EBROAD,PREFIX,TDOS,TNOS &
     &                    ,NB,NKPT,NSPIN,NDIM,EIG,NSET,SET,LEGEND,MODE)
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
      SUBROUTINE REPORT(NFILO,EIG) 
!     **************************************************************************
!     **  WRITES PROJECTED CHARGES AND SPINS FOR EACH ATOM TO                 **
!     **  DPROT FILE AND CALCULATES THE SPIN DIRECTIONS                       **
!     **************************************************************************
      USE PDOS_MODULE
      USE SPINDIR_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFILO
      REAL(8)     ,INTENT(IN) :: EIG(STATE%NB,NKPT,NSPIN)
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
!PRINT*,'STATE%VEC',IKPT,ISPIN,STATE%VEC(1,:,:)
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
        WRITE(NFILO,FMT='(A,T38,"X",T48,"Y",T58,"Z",T66,"TOTAL")') &
     &                   'TOTAL SPIN PROJECTED ON'
!
!       ==  SPIN DIRECTIONS ====================================================
        DO IAT=1,NAT
          WRITE(NFILO,FMT='("SPIN[HBAR/2] ON ATOM",T20,A10,":",4F10.3)')  &
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
          WRITE(NFILO,FMT='(T8,10(4X,A6))') &
     &                       (ATOMID(IATSPINANGLE(IAT2)),IAT2=NATSPINANGLE,2,-1)
          DO IAT1=1,NATSPINANGLE !SENKRECHT
            DO IAT2=NATSPINANGLE,IAT1+1,-1   !WAAGRECHT
              ANGLE(IAT2)=180.D0/PI*ACOS(SUM(SPINDIR(:,IATSPINANGLE(IAT1)) &
     &                                      *SPINDIR(:,IATSPINANGLE(IAT2))))
            END DO
            ITEN=NATSPINANGLE
            DO WHILE (IAT1+1.LE.ITEN)
              WRITE(NFILO,FMT='(A6,10F8.2)')ATOMID(IATSPINANGLE(IAT1)) &
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
      USE ORBITALS_MODULE
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
      SUBROUTINE READCNTL$GENERIC(MODE,TDOS,TNOS,PREFIX)
!     **************************************************************************
!     ** READ !DCNTL!GENERIC BLOCK FROM THE CONTROL FILE                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE, ONLY: LL_CNTL
      IMPLICIT NONE
      CHARACTER(*),INTENT(OUT) :: MODE    ! "SAMPLE" OR "TETRA"
      CHARACTER(*),INTENT(OUT) :: PREFIX  ! PREFIX FOR DOS AND NOS FILES
      LOGICAL(4)  ,INTENT(OUT) :: TDOS    ! DENSITY OF STATES WILL BE PRINTED
      LOGICAL(4)  ,INTENT(OUT) :: TNOS    ! NUMBER OF STATES WILL BE PRINTED
      LOGICAL(4)               :: TCHK
!     **************************************************************************
!     ==========================================================================
!     == SET DEFAULT VALUES                                                   ==
!     ==========================================================================
      MODE='TETRA'
      TDOS=.TRUE.
      TNOS=.FALSE.
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
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DOS',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'DOS',1,TDOS)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NOS',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'NOS',1,TNOS)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'PREFIX',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'PREFIX',1,PREFIX)
      RETURN
      END SUBROUTINE READCNTL$GENERIC
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$GRID(EMIN,EMAX,NE,EBROAD)
!     **************************************************************************
!     ** READ !DCNTL!GRID FROM CONTROL FILE                                   **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE, ONLY: LL_CNTL
      IMPLICIT NONE
      REAL(8)   ,INTENT(INOUT) :: EMIN   ! MINIMUM OF ENERGY GRID
      REAL(8)   ,INTENT(INOUT) :: EMAX   ! MAXIMMUM OF ENERGY GRID
      INTEGER(4),INTENT(OUT)   :: NE     ! NUMBER OF ENERGY GRID POINTS
      REAL(8)   ,INTENT(OUT)   :: EBROAD ! THERMAL ENERGY BROADENING
      REAL(8)                  :: EV     ! ELECTRON VOLT
      REAL(8)                  :: KB     ! BOLTZMANN CONSTANT
      REAL(8)                  :: DE
      LOGICAL(4)               :: TCHK,TCHK1
!     **************************************************************************
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
!     ==  READ ACTUAL VALUES  ==================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMIN[EV]',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'EMIN[EV]',1,EMIN)
        EMIN=EMIN*EV
      END IF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMAX[EV]',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'EMAX[EV]',1,EMAX)
        EMAX=EMAX*EV
      END IF
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
        CALL ERROR$MSG('EITHER !DCNTL!GENERIC:BROADENING[EV]')
        CALL ERROR$MSG('    OR !DCNTL!GENERIC:BROADENING[K]')
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

      CALL LINKEDLIST$EXISTD(LL_CNTL,'SCALEY',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('!DCNTL!GENERIC:SCALEY IS OBSOLETE')
        CALL ERROR$MSG('DOS AND NOS FILES ARE PRODUCED AS SEPARATE FILES')
        CALL ERROR$MSG('DEPENDING ON THE SETTING OF "DOS" AND "NOS"')
        CALL ERROR$STOP('READCNTL$GRID')
      END IF
!
      RETURN
      END SUBROUTINE READCNTL$GRID
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$REPORT1(MODE,TDOS,TNOS,PREFIX,EMIN,EMAX,NE,EBROAD)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MODE    ! "SAMPLE" OR "TETRA"
      CHARACTER(*),INTENT(IN) :: PREFIX  ! PREFIX FOR DOS AND NOS FILES
      LOGICAL(4)  ,INTENT(IN) :: TDOS    ! DENSITY OF STATES WILL BE PRINTED
      LOGICAL(4)  ,INTENT(IN) :: TNOS    ! NUMBER OF STATES WILL BE PRINTED
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
      CALL REPORT$L4VAL(NFILO,'DENSITY OF STATES ARE WRITTEN?',TDOS)
      CALL REPORT$L4VAL(NFILO,'NUMBER OF STATES ARE WRITTEN?',TNOS)
      CALL REPORT$R8VAL(NFILO,'ENERGY GRID STARTS AT',EMIN/EV,'EV')
      CALL REPORT$R8VAL(NFILO,'ENERGY GRID ENDS AT',EMAX/EV,'EV')
      CALL REPORT$R8VAL(NFILO,'SPACING OF THE ENERGY GRID',(EMAX-EMIN)/REAL(NE-1)/EV,'EV')
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
      SUBROUTINE READCNTL$SETS(NB,NKPT,NSPIN,NSET &
     &                        ,NAT,LMX,RBAS,RPOS,LENG,SET,LEGEND)
!     **************************************************************************
!     **                                                                      **
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
      INTEGER(4)   ,INTENT(IN)  :: NB
      INTEGER(4)   ,INTENT(IN)  :: NKPT
      INTEGER(4)   ,INTENT(IN)  :: NSPIN
      INTEGER(4)   ,INTENT(IN)  :: NSET
      INTEGER(4)   ,INTENT(IN)  :: NAT
      INTEGER(4)   ,INTENT(IN)  :: LMX(NAT)
      INTEGER(4)   ,INTENT(IN)  :: LENG
      REAL(8)      ,INTENT(IN)  :: RBAS(3,3)
      REAL(8)      ,INTENT(IN)  :: RPOS(3,NAT)
      REAL(8)      ,INTENT(OUT) :: SET(NB,NKPT,NSPIN,NSET)
      CHARACTER(32),INTENT(OUT) :: LEGEND(NSET)
      INTEGER(4)                :: ISET   ! SET COUNTER
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
        CALL SET$PROJECT(LENG,NB,NKPT,NSPIN,ORBITAL1,ORBITAL2,SET(:,:,:,ISET))
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
          CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',1,TCHK1)
          IF(TCHK1) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,SPIN)
            SPIN=+SPIN
            IF(NSPIN.EQ.2) THEN
              CALL ERROR$MSG('!DCNTL!WEIGHT:TYPE IS ONLY COMPATIBLE WITH')
              CALL ERROR$MSG('NON-COLLINEAR CALCULATIONS.')
              CALL ERROR$MSG('FOR COLLINEAR CALCULATIONS, (NSPIN=2)')
              CALL ERROR$MSG('BOTH SPIN COMPONENTS ARE PRODUCED SIMULTANEOUSLY')
              CALL ERROR$CHVAL('TYPE',TYPE)
              CALL ERROR$CHVAL('SPIN',SPIN)
              CALL ERROR$I4VAL('NSPIN',NSPIN)
              CALL ERROR$STOP('READCNTL$SETS')
            END IF
          ELSE
!           ==  SPIN=TOTAL MEANS "SUM OVER ALL SPINOR COMPONENTS". =============
!           ==  FOR COLLINEAR CALCULATIONS, SPIN='TOTAL' PRODUCES THE TWO ======
!           ==  SPIN DENSITIES. ================================================
            SPIN='TOTAL' 
          END IF
!
!         ======================================================================
!         ==  'TOTAL' = TOTAL DENSITY OF STATES                               ==
!         ======================================================================
          IF(TRIM(TYPE).EQ.'TOTAL') THEN
            SET(:,:,:,ISET)=1.D0
!
!         ======================================================================
!         ==  'ALL' = ALL PROJECTED DENSITY OF STATES                         ==
!         ======================================================================
          ELSE IF(TRIM(TYPE).EQ.'ALL') THEN
            CALL SET$WEIGHT('ALL','',NB,NKPT,NSPIN,SPIN,SET(1,1,1,ISET))
!
!         ======================================================================
!         ==  'EMPTY' = VACCUM DENSITY OF STATES                              ==
!         ======================================================================
          ELSE IF(TRIM(TYPE).EQ.'EMPTY') THEN
            CALL SET$WEIGHT('ALL','',NB,NKPT,NSPIN,SPIN,SET(1,1,1,ISET))
            SET(:,:,:,ISET)=1.D0-SET(:,:,:,ISET)

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
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,SPIN)
          ELSE
            SPIN='TOTAL'
          END IF
!
!         == SELECT TYPE =======================================================
          CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
          TYPE=+TYPE
          IF(TYPE.EQ.'ALL') THEN
            CALL SET$WEIGHT(NAME,'ALL',NB,NKPT,NSPIN,SPIN,SET(1,1,1,ISET))
          ELSE
            IF(SCAN(TYPE,'S').NE.0) THEN
              CALL SET$WEIGHT(NAME,'S',NB,NKPT,NSPIN,SPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'P').NE.0) THEN
              CALL SET$WEIGHT(NAME,'P',NB,NKPT,NSPIN,SPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'D').NE.0) THEN
              CALL SET$WEIGHT(NAME,'D',NB,NKPT,NSPIN,SPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'F').NE.0) THEN
              CALL SET$WEIGHT(NAME,'F',NB,NKPT,NSPIN,SPIN,SET(1,1,1,ISET))
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
          CALL SET$PROJECT(LENG,NB,NKPT,NSPIN,ORBITALI,ORBITALI,SET(1,1,1,ISET))
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
      SUBROUTINE SET$WEIGHT(ATOMID_,ORBITALID,NB_,NKPT_,NSPIN_,SPIN,SET)
!     **************************************************************************
!     ** ADDS A CONTRIBUTION TO THE SET, WHICH CONTAINS THE CONTRIBUTION      **
!     ** TO THE WEIGHT FROM ALL THE STATES                                    **
!     **   SPIN: MAY BE 'TOTAL', 'MAIN', 'X', 'Y', 'Z'                        **
!     **         TOTAL IS THE TOTAL DENSITY OF STATES                         **
!     **         MAIN IS THE PROJECTION ON SPINDIR(I,IAT)                     **
!     **   ORBITALID MAY BE:  'ALL','S','P','D','F'                           **
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE
      USE SPINDIR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: ATOMID_
      CHARACTER(*),INTENT(IN)   :: ORBITALID !MAY BE 'ALL','S','P','D','F'
      INTEGER(4)  ,INTENT(IN)   :: NB_
      INTEGER(4)  ,INTENT(IN)   :: NKPT_
      INTEGER(4)  ,INTENT(IN)   :: NSPIN_
      CHARACTER(*),INTENT(INOUT):: SPIN   
      REAL(8)     ,INTENT(INOUT):: SET(NB_,NKPT_,NSPIN_)
      INTEGER(4)                :: ISPIN,IKPT,IB,IDIM
      INTEGER(4)                :: IPRO0,IPRO1,IPRO2
      INTEGER(4)                :: IAT,IAT0,ISP
      INTEGER(4)                :: L,L1,L2,M,LN,LN1,LN2
      REAL(8)                   :: SUM,SUM_(3)
!     **************************************************************************
                                 CALL TRACE$PUSH('SET$WEIGHT')
      IF(NDIM.EQ.1) SPIN='TOTAL'   ! OVERWRITE IF ONLY ONE CHOICE POSSIBLE
!
!     ==========================================================================
!     ==  SELECT ATOM                                                         ==
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
      IF(TRIM(SPIN).EQ.'MAIN'.AND.IAT0.EQ.-1.AND.NAT.GT.1) THEN
        CALL ERROR$MSG('MAIN SPIN FOR ALL ATOMS NOT POSSIBLE')
        CALL ERROR$CHVAL('SPIN',SPIN)
        CALL ERROR$CHVAL('ATOMID_',ATOMID_)
        CALL ERROR$STOP('SET$WEIGHT')
      END IF
!
!     ==========================================================================
!     ==  SELECT ORBITAL                                                      ==
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
            IF(IAT.EQ.IAT0.OR.IAT0.EQ.-1) THEN
!
              IPRO1=IPRO0
              DO LN1=1,LNX(ISP)
                L1=LOX(LN1,ISP)
                IPRO2=IPRO0
                DO LN2=1,LNX(ISP)
                  L2=LOX(LN2,ISP)
                  IF(L1.NE.L2.OR.(L.GE.0.AND.L1.NE.L)) THEN
                    IPRO2=IPRO2+2*L2+1
                    CYCLE
                  END IF
                  DO IB=1,STATE%NB
                    SUM=0.D0
                    SUM_(:)=0.D0
                    DO M=1,2*L1+1
                      DO IDIM=1,NDIM
                        SUM=SUM+REAL(CONJG(STATE%VEC(IDIM,IPRO1+M,IB)) &
     &                                    *STATE%VEC(IDIM,IPRO2+M,IB))
                      ENDDO
                      IF(NDIM.EQ.2.AND.TRIM(SPIN).NE.'TOTAL') THEN
                        SUM_(1)=SUM_(1) +2.D0*REAL( &
     &                                          CONJG(STATE%VEC(1,IPRO1+M,IB)) &
     &                                                *STATE%VEC(2,IPRO2+M,IB)) 
                        SUM_(2)=SUM_(2)+2.D0*AIMAG( &
     &                                          CONJG(STATE%VEC(1,IPRO1+M,IB)) &
     &                                               *STATE%VEC(2,IPRO2+M,IB))
                        SUM_(3)=SUM_(3)+REAL( &
     &                                          CONJG(STATE%VEC(1,IPRO1+M,IB)) &
     &                                               *STATE%VEC(1,IPRO2+M,IB)  &
     &                                         -CONJG(STATE%VEC(2,IPRO1+M,IB)) &
     &                                               *STATE%VEC(2,IPRO2+M,IB))
                      END IF
                    ENDDO
                    SUM    =SUM    *OV(LN1,LN2,ISP)
                    SUM_(:)=SUM_(:)*OV(LN1,LN2,ISP)
                    IF(TRIM(SPIN).EQ.'TOTAL') THEN
                      SET(IB,IKPT,ISPIN)=SET(IB,IKPT,ISPIN)+SUM
                    ELSE IF(TRIM(SPIN).EQ.'X') THEN
                      SET(IB,IKPT,ISPIN)=SET(IB,IKPT,ISPIN)+SUM_(1)
                    ELSE IF(TRIM(SPIN).EQ.'Y') THEN
                      SET(IB,IKPT,ISPIN)=SET(IB,IKPT,ISPIN)+SUM_(2)
                    ELSE IF(TRIM(SPIN).EQ.'Z') THEN
                      SET(IB,IKPT,ISPIN)=SET(IB,IKPT,ISPIN)+SUM_(3)
                    ELSE IF(TRIM(SPIN).EQ.'MAIN') THEN
                      SET(IB,IKPT,ISPIN)=SET(IB,IKPT,ISPIN) &
     &                                                 +SUM_(1)*SPINDIR(1,IAT) &
     &                                                 +SUM_(2)*SPINDIR(2,IAT) &
     &                                                 +SUM_(3)*SPINDIR(3,IAT)
                    ELSE
                      CALL ERROR$MSG('SPIN COMPONENT NOT RECOGNIZED')
                      CALL ERROR$MSG('SPIN SHOULD BE TOTAL,X,Y,Z OR MAIN')
                      CALL ERROR$CHVAL('SPIN',SPIN)
                      CALL ERROR$STOP('SET$WEIGHT')
                    END IF
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
      SUBROUTINE SET$PROJECT(NPRO_,NB_,NKPT_,NSPIN_,ORBITAL1,ORBITAL2,SET)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NPRO_
      INTEGER(4),INTENT(IN) :: NB_
      INTEGER(4),INTENT(IN) :: NKPT_
      INTEGER(4),INTENT(IN) :: NSPIN_
      COMPLEX(8),INTENT(IN) :: ORBITAL1(NPRO)
      COMPLEX(8),INTENT(IN) :: ORBITAL2(NPRO)
      REAL(8)   ,INTENT(INOUT):: SET(NB_,NKPT_,NSPIN_)
      INTEGER(4)            :: ISPIN,IKPT,IB,IPRO,IDIM
      COMPLEX(8)            :: CSVAR1,CSVAR2
!     **************************************************************************
                                 CALL TRACE$PUSH('SET$PROJECT')
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
          DO IB=1,NB_
            CSVAR1=0.D0
            CSVAR2=0.D0
            DO IPRO=1,NPRO
              DO IDIM=1,NDIM
                CSVAR1=CSVAR1+CONJG(ORBITAL1(IPRO))*STATE%VEC(IDIM,IPRO,IB)
                CSVAR2=CSVAR2+CONJG(ORBITAL2(IPRO))*STATE%VEC(IDIM,IPRO,IB)
              ENDDO                      
            ENDDO
!JOHANNES:  SET(IB,IKPT,ISPIN)=REAL(CONJG(CSVAR1)*CSVAR2)
            SET(IB,IKPT,ISPIN)=SET(IB,IKPT,ISPIN)+REAL(CONJG(CSVAR1)*CSVAR2)
          ENDDO
        ENDDO
      ENDDO
                                 CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL$OUTPUT(EMIN,EMAX,NE,EBROAD,PREFIX,TDOS,TNOS &
     &                    ,NB,NKPT,NSPIN,NDIM,EIG,NSET,SET,LEGEND,MODE)
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
      LOGICAL(4)   ,INTENT(IN) :: TDOS  ! DENSITY OF STATES PRINTED ?
      LOGICAL(4)   ,INTENT(IN) :: TNOS  ! NUMBER OF STATES PRINTED ?
      INTEGER(4)   ,INTENT(IN) :: NB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      INTEGER(4)   ,INTENT(IN) :: NSPIN
      INTEGER(4)   ,INTENT(IN) :: NDIM  !#(SPINOR COMPONENTS)
      REAL(8)      ,INTENT(IN) :: EIG(NB,NKPT,NSPIN)
      INTEGER(4)   ,INTENT(IN) :: NSET
      REAL(8)      ,INTENT(IN) :: SET(NB,NKPT,NSPIN,NSET)
      CHARACTER(32),INTENT(IN) :: LEGEND(NSET)
      CHARACTER(32),INTENT(IN) :: MODE
      CHARACTER(32)        :: LEGEND1
      CHARACTER(256)       :: FILE
      LOGICAL(4)           :: TIB,TE,TIK,TIS,TCHK
      INTEGER(4)           :: IB,IKPT,ISPIN,I,IOUT
      INTEGER(4)           :: ISET
      INTEGER(4)           :: NOUT
      INTEGER(4)           :: NFIL
      INTEGER(4)           :: NFILDOS  ! UNIT FOR DENSITY OF STATES FILE
      INTEGER(4)           :: NFILNOS  ! UNIT FOR NUMBER OF STATES FILE
      INTEGER(4)           :: NFILO
      INTEGER(4)           :: IB0,IK0,IS0
      REAL(8)              :: EV
      REAL(8)              :: ENERGY
      REAL(8)              :: E1,E2
      INTEGER(4)           :: IE1,IE2
      REAL(8)              :: SUM,SUMS,SVAR 
      LOGICAL(4),ALLOCATABLE :: DEADZONE(:,:) !(NE,ISPIN)
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
      ALLOCATE(DEADZONE(NE,NSPIN))
      DEADZONE(:,:)=.TRUE.
      DO ISPIN=1,NSPIN
        DO IB=1,NB
          E1=MINVAL(EIG(IB,:,ISPIN))-1.D0*EV
          E2=MAXVAL(EIG(IB,:,ISPIN))+1.D0*EV
          IE1=max(1, 1+NINT((E1-EMIN)/(EMAX-EMIN)*REAL(NE-1)))
          IE2=min(ne,1+NINT((E2-EMIN)/(EMAX-EMIN)*REAL(NE-1)))
          DEADZONE(IE1:IE2,ISPIN)=.FALSE.
        ENDDO
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
        FILE=TRIM(PREFIX)//TRIM(LEGEND1)//-'.NOS'
        CALL FILEHANDLER$SETFILE('PNOSOUT',.FALSE.,FILE)
        NFILDOS=-1
        NFILNOS=-1
        IF(TDOS) CALL FILEHANDLER$UNIT('PDOSOUT',NFILDOS)
        IF(TNOS)CALL FILEHANDLER$UNIT('PNOSOUT',NFILNOS)
!
!       ========================================================================
!       ==  WRITE DOS AND INTEGRATED DOS ON FILE                              ==
!       ========================================================================
        IF(MODE.EQ.'SAMPLE')THEN
          CALL PUTONGRID_SAMPLE(NFILDOS,NFILNOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
      &                 ,NB,NKPT,NSPIN,NDIM,EIG,SET(:,:,:,ISET),LEGEND(ISET))
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
          CALL PUTONGRID_TETRA(NFILDOS,NFILNOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
      &                 ,NB,NKPT,NSPIN,NDIM,EIG,SET(:,:,:,ISET),LEGEND(ISET))
        ELSE
          CALL ERROR$MSG('VALUE FOR "MODE" NOT RECOGNIZED')
          CALL ERROR$MSG('MUST BE "SAMPLE" OR "TETRA"')
          CALL ERROR$CHVAL('MDOE',MODE)   
          CALL ERROR$STOP('READCNTL$OUTPUT')
        ENDIF
!
!       == CLOSE DOWN ==========================================================
        IF(NFILDOS.GE.0)CALL FILEHANDLER$CLOSE('PDOSOUT')
        IF(NFILNOS.GE.0)CALL FILEHANDLER$CLOSE('PNOSOUT')
      ENDDO  ! END OF LOOP OVER SETS
!
!     ==========================================================================
!     ==========================================================================
!     == SCAN SPECIFIC OUTPUT BLOCKS                                          ==
!     ==========================================================================
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'OUTPUT',NOUT)
      DO IOUT=1,NOUT
        CALL LINKEDLIST$SELECT(LL_CNTL,'OUTPUT',IOUT)
                          CALL TRACE$PASS('NEXT IOUT')
!
!       ==  SPECIFY SET ========================================================
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,LEGEND1)
        ISET=0
        DO I=1,NSET
          IF(TRIM(LEGEND1).EQ.TRIM(LEGEND(I))) THEN
            ISET=I ; EXIT
          END IF
        ENDDO
        IF(ISET.EQ.0) THEN
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          DO I=1,NSET
            PRINT*,'LEGEND ',I,TRIM(LEGEND(I)),'::'
          ENDDO
          CALL ERROR$CHVAL('LEGEND1',LEGEND1)
          CALL ERROR$STOP('READCNTL$OUTPUT')
        END IF
 !
!       ==  OUTPUT TYPE ========================================================
!       ==  (B,K,S) SPECIFIES A SPECIFIC STATE REPORTED IN THE PROTOCOLL
!       ==  E[EV]    SPECIFIES AN ENERGY
!       ==  OTHERWISE THE INFORMATION ON THE GRID IS WRITTEN TO FILE
        CALL LINKEDLIST$EXISTD(LL_CNTL,'B',1,TIB)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'E[EV]',1,TE)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'K',1,TIK)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'S',1,TIS)
        IF((.NOT.(TIB.OR.TE)).AND.(TIK.OR.TIS)) THEN
          CALL ERROR$MSG('K-POINT AND SPIN MUST NOT BE SELECTED') 
          CALL ERROR$MSG('                FOR DENSITY OF STATES')
          CALL ERROR$STOP('READCNTL$OUTPUT') 
        END IF
        IK0=0
        IS0=0
        IF(TIK)CALL LINKEDLIST$GET(LL_CNTL,'K',1,IK0)
        IF(TIS)CALL LINKEDLIST$GET(LL_CNTL,'S',1,IS0)
        IF(TIB.AND.TE) THEN
          CALL ERROR$MSG('NB AND E ARE INCOMPATIBLE')
          CALL ERROR$STOP('READCNTL$OUTPUT')
        END IF
!
!       ========================================================================
!       ==  WRITE PROJECTION OF A GIVEN STATE ON FILE                         ==
!       ========================================================================
        IF(TIB) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'B',1,IB0)
          DO ISPIN=1,NSPIN
            IF(TIS.AND.(ISPIN.NE.IS0)) CYCLE
            DO IKPT=1,NKPT
              IF(TIK.AND.(IKPT.NE.IK0)) CYCLE
              WRITE(NFILO,FMT='(A32,"; B=",I3,"; K=",I2,"; S=",I1 &
     &         ,"; E[EV]=",F10.5,";TOTAL PRO=",F10.5)') &
     &                LEGEND(ISET),IB0,IKPT,ISPIN &
     &               ,EIG(IB0,IKPT,ISPIN)/EV,SET(IB0,IKPT,ISPIN,ISET)
            ENDDO
          ENDDO
        END IF
!
!       ========================================================================
!       ==  WRITE INTEGRATED DENSITY OF STATES AT A GIVEN ENERGY              ==
!       ========================================================================
        IF(TIB) THEN
          ENERGY=-1.D+10
          SUM=0.D0
          SUMS=0.D0
          DO ISPIN=1,NSPIN
            IF(TIS.AND.(ISPIN.NE.IS0)) CYCLE
            IF(ISPIN.EQ.1) SVAR=1.D0
            IF(ISPIN.EQ.2) SVAR=-1.D0
            DO IKPT=1,NKPT
              IF(TIK.AND.(IKPT.NE.IK0)) CYCLE
              DO IB=1,IB0
                SUM=SUM+SET(IB,IKPT,ISPIN,ISET)
                SUMS=SUMS+SET(IB,IKPT,ISPIN,ISET)*SVAR
                ENERGY=MAX(ENERGY,EIG(IB,IKPT,ISPIN))
              ENDDO
            ENDDO        
          ENDDO
          IF(.NOT.TIS) SUM=SUM*2.D0/DBLE(NSPIN)  ! NO K-POINT WEIGHT ? NO NDIM ?
          WRITE(NFILO,FMT='(A32,"; B=",I3,"; K=",I2,"; S=",I1 &
     &                  ,"; E[EV]=",F10.5,"; TOTAL NOS=",F10.5)') &
     &                  LEGEND(ISET),IB0,IK0,IS0,ENERGY/EV,SUM
!         IF(NSPIN.EQ.2.AND.(.NOT.TIS)) THEN
!           WRITE(NFIL,FMT='(A32,"; B=",I3,"; K=",I2,"; S=",I1 &
!    &                  ,";       ",10X,  "; SOS=",F10.5)') &
!    &                  LEGEND(ISET),IB0,IK0,IS0,SUMS
!         END IF
        END IF
!
!       ========================================================================
!       ==  WRITE INTEGRATED DENSITY OF STATES AT A GIVEN ENERGY              ==
!       ========================================================================
        IF(TE) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'E[EV]',1,ENERGY)
          ENERGY=ENERGY*EV
          SUM=0.D0
          SUMS=0.D0
          DO ISPIN=1,NSPIN
            IF(TIS.AND.(ISPIN.NE.IS0)) CYCLE
            IF(ISPIN.EQ.1) SVAR=1.D0
            IF(ISPIN.EQ.2) SVAR=-1.D0
            DO IKPT=1,NKPT
              IF(TIK.AND.(IKPT.NE.IK0)) CYCLE
              DO IB=1,NB
                IF(EIG(IB,IKPT,ISPIN).LT.ENERGY) THEN
                  SUM=SUM+SET(IB,IKPT,ISPIN,ISET)
                  SUMS=SUMS+SET(IB,IKPT,ISPIN,ISET)*SVAR
                END IF
              ENDDO
            ENDDO        
          ENDDO
          IF(.NOT.TIS) SUM=SUM*2.D0/DBLE(NSPIN)  ! NO K-POINT WEIGHT ? NO NDIM ?
          WRITE(NFILO,FMT='(A32,8X," K=",I2,"; S=",I1 &
     &                  ,"; E[EV]=",F10.5,"; TOTAL NOS=",F10.5)') &
     &                  LEGEND(ISET),IK0,IS0,ENERGY/EV,SUM
          IF(NSPIN.EQ.2.AND.(.NOT.TIS)) THEN
            WRITE(NFILO,FMT='(A32,8X,"; K=",I2,"; S=",I1 &
     &                  ,"; E[EV]=",F10.5,"; TOTAL SOS=",F10.5)') &
     &                  LEGEND(ISET),IK0,IS0,ENERGY/EV,SUMS
          END IF
        END IF
!
!       == CLOSE DOWN ==========================================================
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO  ! END OF LOOP OVER OUTPUT BLOCKS
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READCNTL$OUTPUT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PUTONGRID_SAMPLE(NFILDOS,NFILNOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
     &                    ,NB,NKPT,NSPIN,NDIM,EIG,SET,LEGEND)
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
      INTEGER(4)   ,INTENT(IN) :: NB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      INTEGER(4)   ,INTENT(IN) :: NSPIN
      INTEGER(4)   ,INTENT(IN) :: NDIM
      REAL(8)      ,INTENT(IN) :: EIG(NB,NKPT,NSPIN)
      REAL(8)      ,INTENT(IN) :: SET(NB,NKPT,NSPIN)
      INTEGER(4)   ,INTENT(IN) :: NFILDOS  ! UNIT FOR DOS FILE OR "-1"
      INTEGER(4)   ,INTENT(IN) :: NFILNOS  ! UNIT FOR NOS FILE OR "-1"
      CHARACTER(32),INTENT(IN) :: LEGEND
      LOGICAL(4)   ,INTENT(IN) :: DEADZONE(NE,NSPIN)
      REAL(8)      ,PARAMETER  :: TOL=1.D-2
      REAL(8)              :: DE
      INTEGER(4)           :: IE1,IE2,IDE
      INTEGER(4)           :: ND,IOCC
      REAL(8)              :: NOS(NE,NSPIN,2)
      REAL(8)              :: DOS(NE,NSPIN,2)
      REAL(8),ALLOCATABLE  :: SMEAR(:)
      REAL(8)              :: EV
      REAL(8)              :: W1,W2,X,FAC
      REAL(8)              :: NOSSMALL(NSPIN,2)
      REAL(8)              :: WGHTX
!      REAL(8)              :: WKPT(NKPT)
      INTEGER(4)           :: IKPT,ISPIN,IE,IB
      REAL(8)              :: SVAR
      REAL(8)              :: SIG
      REAL(8)              :: E
      REAL(8)              :: SPINDEG
!     **************************************************************************
                                 CALL TRACE$PUSH('PUTONGRID_SAMPLE')
      CALL CONSTANTS('EV',EV)
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)   !STEP OF THE ENERGY GRID
      SPINDEG=1.D0
      IF(NSPIN.EQ.1.AND.NDIM.EQ.1) SPINDEG=2.D0
!
!     == THIS IS A DIRTY FIX THAT WAS NECESSARY BEFORE THE K-POINT WEIGHT WAS
!     == AVAILABLE ON THE PDOS FILE 
      IF(SUM(WKPT).EQ.0.D0) THEN
        CALL ERROR$MSG('NO K-POINT WEIGHTS AVAILABLE: OLD VERSION OF PDOS FILE')
        CALL ERROR$MSG('RERUN ONE PAW ITERATION WITH NEW CODE')
        CALL ERROR$STOP('PUTONGRID_SAMPLE')
        DO IKPT=1,NKPT
          WKPT(IKPT)=0.D0
          DO ISPIN=1,NSPIN
            STATE=>STATEARR(IKPT,ISPIN)
!           == CAUTION: HERE I ESTIMATE THE WEIGHT AND SPIN-DEGENERACY FACTOR 
!           == FROM THE MAX OCCUPATION, WHICH MAY BE INCORRECT
            WKPT(IKPT)=MAX(WKPT(IKPT),MAXVAL(STATE%OCC(:)))
          ENDDO
          IF(WKPT(IKPT).EQ.0.D0) THEN
            CALL ERROR$MSG('NO ELECTRONS FOR THIS K-POINT. GOT CONFUSED')
            CALL ERROR$STOP('PUTONGRID_SAMPLE')
          END IF
        ENDDO       
      END IF
!
!     ==========================================================================
!     ==  MAP CONTRIBUTION FROM EACH STATE) ONTO THE ENERGY GRID.             ==
!     ==  (IT IS DIVIDED PROPORTIONALLY TO THE TWO ENCLOSING GRID POINTS)     ==
!     == THE SUM OVER ALL NOS-DATA ADDS TO THE TOTAL NUMBER OF STATES         ==
!     == STATES LYING BELOW EMIN ARE MAPPED INTO NOSSMALL                     ==
!     ==========================================================================
      NOS(:,:,:)=0.D0
      NOSSMALL(:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
!         == CAUTION: HERE I ESTIMATE THE WEIGHT AND SPIN-DEGENERACY FACTOR 
!         == FROM THE MAX OCCUPATION, WHICH MAY BE INCORRECT
          WGHTX=WKPT(IKPT)*SPINDEG
          DO IB=1,NB
            X=(EIG(IB,IKPT,ISPIN)-EMIN)/DE+1.D0
            IE1=INT(X)
            IE2=IE1+1
            W2=(X-REAL(IE1,KIND=8))
            W1=1.D0-W2
            IF (IE1.LT.1) THEN
              NOSSMALL(ISPIN,1)=NOSSMALL(ISPIN,1)+SET(IB,IKPT,ISPIN)*WGHTX
              NOSSMALL(ISPIN,2)=NOSSMALL(ISPIN,2) &
     &                                         +SET(IB,IKPT,ISPIN)*STATE%OCC(IB)
            ELSE
              IF(IE1.LE.NE.AND.IE1.GE.1) THEN 
                NOS(IE1,ISPIN,1)=NOS(IE1,ISPIN,1)+W1*SET(IB,IKPT,ISPIN)*WGHTX
                NOS(IE1,ISPIN,2)=NOS(IE1,ISPIN,2) &
     &                                      +W1*SET(IB,IKPT,ISPIN)*STATE%OCC(IB)
              END IF
              IF(IE2.LE.NE.AND.IE2.GE.1) THEN
                NOS(IE2,ISPIN,1)=NOS(IE2,ISPIN,1)+W2*SET(IB,IKPT,ISPIN)*WGHTX
                NOS(IE2,ISPIN,2)=NOS(IE2,ISPIN,2) &
     &                                      +W2*SET(IB,IKPT,ISPIN)*STATE%OCC(IB)
              END IF
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
!     ==  SMEAR DENSITY OF STATES ==============================================
      DOS(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IOCC=1,2
          DO IDE=-ND,ND
            IE1=MAX(1,1-IDE)
            IE2=MIN(NE,NE-IDE)
            W1=SMEAR(IDE)
            DO IE=-IDE+1,IE1-1
              NOSSMALL(ISPIN,IOCC)=NOSSMALL(ISPIN,IOCC) &
    &                             +NOS(IE+IDE,ISPIN,IOCC)*W1
            ENDDO
            DO IE=IE1,IE2
              DOS(IE,ISPIN,IOCC)=DOS(IE,ISPIN,IOCC)+NOS(IE+IDE,ISPIN,IOCC)*W1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  DETERMINE NOS BY SUMMATION                                          ==
!     ==========================================================================
      NOS=DOS
      DO ISPIN=1,NSPIN
        DO IOCC=1,2
          NOS(1,ISPIN,IOCC)=NOS(1,ISPIN,IOCC)+NOSSMALL(ISPIN,IOCC)
          DO IE=2,NE
            NOS(IE,ISPIN,IOCC)=NOS(IE,ISPIN,IOCC)+NOS(IE-1,ISPIN,IOCC)
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
      DO ISPIN=1,NSPIN
        DO IOCC=1,2
          DO IE=1,NE
            IF(ABS(NOS(IE,ISPIN,IOCC)).LE.1.D-99)NOS(IE,ISPIN,IOCC)=0.D0
            IF(ABS(DOS(IE,ISPIN,IOCC)).LE.1.D-99)DOS(IE,ISPIN,IOCC)=0.D0
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  WRITE RESULT ON PSODOUT                                             ==
!     ==========================================================================
!     == WRITE DENSITY OF STATES ===============================================
      IF(NFILDOS.GE.0) THEN
        DO ISPIN=1,NSPIN
          SIG=REAL(3-2*ISPIN) ! +1 FOR ISPIN=1 / -1 FOR ISPIN=2
          WRITE(NFILDOS,FMT='(F14.8,2F14.8)')EMIN/EV,0.D0,0.D0
          DO IE=1,NE
            IF(DEADZONE(IE,ISPIN)) CYCLE
            E=EMIN+(EMAX-EMIN)*REAL(IE-1)/REAL(NE-1)
            WRITE(NFILDOS,FMT='(F14.8,2F14.8)')E/EV,SIG*DOS(IE,ISPIN,1)*EV &
                                                   ,SIG*DOS(IE,ISPIN,2)*EV
          ENDDO
          WRITE(NFILDOS,FMT='(F14.8,2F14.8)')EMAX/EV,0.D0,0.D0
        ENDDO
        WRITE(NFILDOS,FMT='("# THIS WAS: ",A)')LEGEND
      END IF
!
!     == WRITE NUMBER OF STATES ================================================
      IF(NFILNOS.GE.0) THEN
        DO ISPIN=1,NSPIN
          SIG=REAL(3-2*ISPIN) ! +1 FOR ISPIN=1 / -1 FOR ISPIN=2
          WRITE(NFILNOS,FMT='(F14.8,2F14.8)')EMIN/EV,0.D0,0.D0
          DO IE=1,NE
            IF(DEADZONE(IE,ISPIN)) CYCLE
            E=EMIN+(EMAX-EMIN)*REAL(IE-1)/REAL(NE-1)
            WRITE(NFILNOS,FMT='(F14.8,2F14.8)')E/EV,SIG*NOS(IE,ISPIN,1) &
                                                   ,SIG*NOS(IE,ISPIN,2)
          ENDDO
          WRITE(NFILNOS,FMT='(F14.8,2F14.8)')EMAX/EV,0.D0,0.D0
        ENDDO
        WRITE(NFILNOS,FMT='("# THIS WAS: ",A)')LEGEND
      END IF
                                 CALL TRACE$POP
      RETURN
    END SUBROUTINE PUTONGRID_SAMPLE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PUTONGRID_TETRA(NFILDOS,NFILNOS,EMIN,EMAX,NE,EBROAD,DEADZONE &
     &                    ,NB,NKPT,NSPIN,NDIM,EIG,SET,LEGEND)
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
      INTEGER(4)   ,INTENT(IN) :: NB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      INTEGER(4)   ,INTENT(IN) :: NSPIN
      INTEGER(4)   ,INTENT(IN) :: NDIM
      REAL(8)      ,INTENT(IN) :: EIG(NB,NKPT,NSPIN)
      REAL(8)      ,INTENT(IN) :: SET(NB,NKPT,NSPIN)
      INTEGER(4)   ,INTENT(IN) :: NFILDOS ! UNIT FOR DOS FILE OR "-1"
      INTEGER(4)   ,INTENT(IN) :: NFILNOS ! UNIT FOR NOS FILE OR "-1"
      CHARACTER(32),INTENT(IN) :: LEGEND
      LOGICAL(4)   ,INTENT(IN) :: DEADZONE(NE,NSPIN)
      REAL(8)      ,PARAMETER  :: TOL=1.D-2
      REAL(8)    ,ALLOCATABLE  :: NOS(:,:)
      REAL(8)    ,ALLOCATABLE  :: NOSMIN(:)
      REAL(8)    ,ALLOCATABLE  :: DOS(:,:)
      REAL(8)    ,ALLOCATABLE  :: SMEAR(:)
      REAL(8)                  :: DE
      INTEGER(4)               :: IE1,IE2,IDE,I1,I2,IBIS
      INTEGER(4)               :: ND,IOCC
      REAL(8)                  :: EV
      REAL(8)                  :: W1,W2,X,FAC
      REAL(8)                  :: NOSSMALL(NSPIN,2)
      REAL(8)                  :: WGHTX
      INTEGER(4)               :: IKPT,ISPIN,IE,IB
      REAL(8)                  :: SVAR
      REAL(8)                  :: SIG
      REAL(8)                  :: E
      LOGICAL(4)               :: TSKIP
      REAL(8)                  :: SPINDEG
!     **************************************************************************
                                 CALL TRACE$PUSH('PUTONGRID_TETRA')
      CALL CONSTANTS('EV',EV)
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)  !STEP OF THE ENERGY GRID
      SPINDEG=1.D0
      IF(NSPIN.EQ.1.AND.NDIM.EQ.1) SPINDEG=2.D0
!
!     ==========================================================================
!     ==  CALCULATE NUMBER OF STATES FOR EACH ENERGY INTERVAL                 ==
!     ==========================================================================
      ALLOCATE(DOS(NE,NSPIN))
      ALLOCATE(NOSMIN(NSPIN))
      ALLOCATE(NOS(NE,NSPIN))
      DOS(:,:)=0.D0
      NOSMIN(:)=0.D0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            STATE=>STATEARR(IKPT,ISPIN)
            IBIS=IB+NB*(ISPIN-1)
            SVAR=SPINDEG*SET(IB,IKPT,ISPIN)*DE
            NOSMIN(ISPIN)=NOSMIN(ISPIN)+WGHT(IBIS,IKPT)*SVAR
            I1=EWGHT(IB+NB*(ISPIN-1),IKPT)%I1
            I2=EWGHT(IB+NB*(ISPIN-1),IKPT)%I2
            DO IE=I1,I2
              DOS(IE,ISPIN)=DOS(IE,ISPIN)+EWGHT(IBIS,IKPT)%WGHT(IE)*SVAR
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
      DO ISPIN=1,NSPIN
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
!     ==  DETERMINE NOS BY INTEGRATION                                        ==
!     ==========================================================================
      DO ISPIN=1,NSPIN
        NOS(1,ISPIN)=NOSMIN(ISPIN)
        DO IE=2,NE
          NOS(IE,ISPIN)=NOS(IE-1,ISPIN)+DOS(IE,ISPIN)*DE
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ROUND INSIGNIFICANT VALUES TO ZERO                                  ==
!     ==========================================================================
      DO ISPIN=1,NSPIN
        DO IE=1,NE
          IF(ABS(NOS(IE,ISPIN)).LE.1.D-99) NOS(IE,ISPIN)=0.D0
          IF(ABS(DOS(IE,ISPIN)).LE.1.D-99) DOS(IE,ISPIN)=0.D0
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  WRITE RESULT ON PSODOUT                                             ==
!     ==========================================================================
!     == WRITE DENSITY OF STATES ===============================================
      IF(NFILDOS.GE.0) THEN
        DO ISPIN=1,NSPIN
          SIG=REAL(3-2*ISPIN) ! +1 FOR ISPIN=1 / -1 FOR ISPIN=2
          WRITE(NFILDOS,FMT='(F14.8,2F14.8)')EMIN/EV,0.D0,0.D0
          DO IE=1,NE
            IF(DEADZONE(IE,ISPIN)) CYCLE
            E=EMIN+(EMAX-EMIN)*REAL(IE-1,KIND=8)/REAL(NE-1)
            SVAR=1.D0
            IF(E.GT.EF) SVAR=0.D0
            WRITE(NFILDOS,FMT='(F14.8,2F14.8)')E/EV,SIG*DOS(IE,ISPIN)*EV &
                                              ,SIG*DOS(IE,ISPIN)*SVAR*EV
          ENDDO
          WRITE(NFILDOS,FMT='(F14.8,2F14.8)')EMAX/EV,0.D0,0.D0
        ENDDO
        WRITE(NFILDOS,FMT='("# THIS WAS: ",A)')LEGEND
      END IF
!
!     == WRITE NUMBER OF STATES ================================================
      IF(NFILNOS.GE.0) THEN
        DO ISPIN=1,NSPIN
          SIG=REAL(3-2*ISPIN) ! +1 FOR ISPIN=1 / -1 FOR ISPIN=2
          WRITE(NFILNOS,FMT='(F14.8,2F14.8)')EMIN/EV,0.D0,0.D0
          IE1=1
          DO IE=1,NE
            IF(DEADZONE(IE,ISPIN)) CYCLE
            E=EMIN+(EMAX-EMIN)*REAL(IE-1,KIND=8)/REAL(NE-1)
            IF(E.LT.EF)IE1=IE
            WRITE(NFILNOS,FMT='(F14.8,2F14.8)')E/EV,SIG*NOS(IE,ISPIN) &
      &                                            ,SIG*NOS(IE1,ISPIN)
          ENDDO
          WRITE(NFILNOS,FMT='(F14.8,2F14.8)')EMAX/EV,0.D0,0.D0
        ENDDO
        WRITE(NFILNOS,FMT='("# THIS WAS: ",A)')LEGEND
      END IF
                                 CALL TRACE$POP
      RETURN
      END SUBROUTINE PUTONGRID_TETRA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GENERATE_TETRA_WGHT(NFILO,NB,NSPIN,NKPT,EMAX,EMIN,NE,RBAS,EIG)
      USE DOS_WGHT_MODULE, ONLY: EF,SPACEGROUP,WGHT,EWGHT
      USE BRILLOUIN_MODULE, ONLY: EWGHT_TYPE
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)        :: NFILO
      INTEGER(4),INTENT(IN)        :: NB
      INTEGER(4),INTENT(IN)        :: NSPIN
      INTEGER(4),INTENT(IN)        :: NKPT
      REAL(8),INTENT(IN)           :: EMIN
      REAL(8),INTENT(IN)           :: EMAX
      INTEGER(4),INTENT(IN)        :: NE
      REAL(8),INTENT(IN)           :: RBAS(3,3)
      REAL(8),INTENT(IN)           :: EIG(NB,NKPT,NSPIN)
      INTEGER(4)                   :: ISHIFT(3)
      INTEGER(4)                   :: NKDIV(3)
      INTEGER(4)                   :: NKPT2
      LOGICAL(4)                   :: TINV
      REAL(8),ALLOCATABLE          :: EB(:,:)
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
      CALL PDOS$GETR8('RNTOT',RNTOT)
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
      ALLOCATE(EB(NB*NSPIN,NKPT))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          EB(1+NB*(ISPIN-1):NB+NB*(ISPIN-1),IKPT)=EIG(1:NB,IKPT,ISPIN) 
        ENDDO
      ENDDO
      ALLOCATE(WGHT(NB*NSPIN,NKPT))
      CALL BRILLOUIN$DOS(NSPIN*NB,NKPT,EB,WGHT,RNTOT,EF)
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
      ALLOCATE(EWGHT(NB*NSPIN,NKPT))
      CALL BRILLOUIN$WGHT(NKPT,NB*NSPIN,EMIN,EB,WGHT)
      CALL BRILLOUIN$EWGHT(NKPT,NB*NSPIN,EB,EMIN,EMAX,NE,EWGHT)
!     
!     ==========================================================================
!     == REPORT INFORMATION                                                   ==
!     ==========================================================================
      CALL REPORT$R8VAL(NFILO,"FERMI LEVEL",EF/EV,'EV')
                          CALL TRACE$POP()
      RETURN
      END SUBROUTINE GENERATE_TETRA_WGHT

