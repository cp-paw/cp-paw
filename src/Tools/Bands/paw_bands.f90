!........1.........2.........3.........4.........5.........6.........7.........8
!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**  NAME: BANDS                                                              **
!**                                                                           **
!**  PURPOSE: ANALYSIS TOOL FOR BAND STRUCTURE                                **
!**                                                                           **
!*******************************************************************************
!*******************************************P.E.BLOECHL, GOSLAR JUNE 20,2010****
MODULE REFCELL_MODULE
REAL(8)      :: GREF(3,3)
CHARACTER(3) :: BRAVAIS ! ID OF TH THE BRAVAIS LATTICE
                        ! SEE BRADLEY CRACKNELL TABLE 3.3
INTEGER(4)               :: NP=0
REAL(8)     ,ALLOCATABLE :: P(:,:)   !(3,NP)  GVECTOR OF HIGHSYMMETRY POINT
CHARACTER(1),ALLOCATABLE :: PNAME(:) !(NP)
END MODULE REFCELL_MODULE

      PROGRAM MAIN
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)             :: LL_CNTL
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NFIL
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: THISTASK,NTASKS
      CHARACTER(32)             :: MODE
      INTEGER(4)                :: NBANDS,NPDOS
!     ******************************************************************
!
!     ==========================================================================
!     == INITIALIZE MPE ROUTINE FOR PARALLEL PROCESSING                       ==
!     ==========================================================================
      CALL MPE$INIT
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      CALL PAW_VERSION
      CALL TIMING$START
      CALL TRACE$PUSH('MAIN')
!
!     ==========================================================================
!     ==  RESOLVE ARGUMENTLIST AND INITIALIZE FILE HANDLER                    ==
!     ==========================================================================
      CALL INITIALIZEFILEHANDLER
!
!     ==========================================================================
!     ==  READ CONTROL FILE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('BCNTL',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
!
      IF(THISTASK.EQ.1)THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!       ========================================================================
!       ==  WRITE HEADER                                                      ==
!       ========================================================================
        WRITE(NFILO,FMT='(80("*"))')
        WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"             BANDS ANALYSIS TOOL                ")')
        WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"      FOR THE PROJECTOR-AUGMENTED WAVE METHOD   ")')
        WRITE(NFILO,FMT='(80("*"))')
        WRITE(NFILO,FMT='(T10 &
     &           ,"P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY")')
        WRITE(NFILO,FMT='(T10 &
     &            ,"DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3")')
        WRITE(NFILO,*)
      ENDIF  
!
!     ==========================================================================
!     ==  SELECT METHOD                                                       ==
!     ==========================================================================
      !DO BANDSTRUCTURE
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'BANDSTRUCTURE',NBANDS)
      IF(NBANDS.EQ.1)THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
        !READ MODE
        CALL LINKEDLIST$EXISTD(LL_CNTL,'MODE',0,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'MODE',1,MODE)
        ELSE
          MODE='LINEARINTERPOLATION'
        ENDIF
        IF(MODE.EQ.'DIAGONALISATION')THEN
          CALL BANDS_BANDSTRUCTURE_DIAG(LL_CNTL,NFIL,NFILO)
        ELSE IF(MODE.EQ.'LINEARINTERPOLATION')THEN
          IF(THISTASK.EQ.1)THEN
            CALL BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION(LL_CNTL,NFIL,NFILO)
          ENDIF
        ELSE
          CALL ERROR$MSG('BANDS: MODE UNKNOWN')
          CALL ERROR$CHVAL('MODE',MODE)
          CALL ERROR$STOP('PAW_BANDS')
        ENDIF
      ELSE IF(NBANDS.GT.1)THEN
        CALL ERROR$MSG('TOO MANY BANDSTRCTURE-BLOCKS GIVEN')
        CALL ERROR$MSG('JUST ONE BLOCK ALLOWED')
        CALL ERROR$I4VAL('NBANDS',NBANDS)
        CALL ERROR$STOP('PAW_BANDS')
      ENDIF
      
      !DO PDOS
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'PDOS',NPDOS)
      IF(NPDOS.EQ.1)THEN
        CALL BANDS_PDOS(LL_CNTL,NFIL,NFILO)
      ELSE IF(NPDOS.GT.1)THEN
        CALL ERROR$MSG('TOO MANY PDOS-BLOCKS GIVEN, JUST ONE BLOCK ALLOWED')
        CALL ERROR$I4VAL('NPDOS',NPDOS)
        CALL ERROR$STOP('PAW_BANDS')
      ENDIF
!
!     ==========================================================================
!     ==  CLOSING                                                             ==
!     ==========================================================================
      IF(THISTASK.EQ.1)THEN
        CALL FILEHANDLER$REPORT(NFILO,'USED')
        WRITE(NFILO,FMT='(72("="))')
        WRITE(NFILO,FMT='(72("="),T20,"  PAW_BANDS TOOL FINISHED  ")')
        WRITE(NFILO,FMT='(72("="))')
                            CALL TRACE$PASS('AFTER CLOSING')
      ENDIF
!
!     ==========================================================================
!     ==  CLOSE FILES                                                         ==
!     ==========================================================================
      CALL FILEHANDLER$CLOSEALL
                            CALL TRACE$PASS('AFTER FILEHANDLER$CLOSEALL')
      CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      END PROGRAM MAIN
!      
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL_REFCELL(LL_CNTL)
!     **************************************************************************
!     ** READS A REFERENCE UNIT CELL FOR THE INTERPRETATION OF K-POINT        **
!     ** COORDINATES. K-POINTS ARE SPECIFIED IN RELATIVE COORDINATES          **
!     ** AND WILL BE TRANSFORMED INTO CARTESIAN COORDINATES BY MULTIPLICATION **
!     ** WITH GREF                                                            **
!
!     !REFCELL LUNIT= T= G= ID= !END
!     !KPOINT ID='G' GREL= 0. 0. 0.5 !END
!     !LINE POINTS= 'G' 'X' 'N' '/' 'U' 'G' !END
!
!     **************************************************************************
      USE REFCELL_MODULE, ONLY : GREF,BRAVAIS,NP,P,PNAME
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)             :: LL_CNTL
      LOGICAL(4)                :: TCHK,TCHK1,TCHK2
      REAL(8)                   :: TREF(3,3)
      REAL(8)                   :: LUNIT
      REAL(8)                   :: DET
      INTEGER(4)                :: I
      INTEGER(4),PARAMETER      :: NPX=20
      REAL(8)                   :: P1(3,NPX)
      CHARACTER(1)              :: PNAME1(NPX)
!     **************************************************************************
!
!     ==========================================================================
!     ==  SET DEFAULTS                                                        ==
!     ==========================================================================
      TREF=0.D0
      GREF=0.D0
      DO I=1,3
        TREF(I,I)=1.D0
        GREF(I,I)=1.D0
      ENDDO
!
!     ==========================================================================
!     == READ  REFERENCE UNIT CELL GREF                                       ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'REFCELL',0,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'REFCELL')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'LUNIT',0,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'LUNIT',1,LUNIT)

      CALL LINKEDLIST$EXISTD(LL_CNTL,'T',0,TCHK1)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'G',0,TCHK2)
      IF(TCHK1.AND.TCHK2) THEN
        CALL ERROR$MSG('!BCNTL!REFCELL:T AND G ARE MUTUALLY EXCLUSIVE')
        CALL ERROR$STOP('READCNTL_REFCELL')
      END IF
      IF(TCHK1) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'T',1,TREF)
        TREF=TREF*LUNIT
        CALL GBASS(TREF,GREF,DET)
      ELSE IF(TCHK2) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'G',1,GREF)
        GREF=GREF/LUNIT
      END IF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'BRAVAIS',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,BRAVAIS)
        CALL SPACEGROUP$SYMMETRYPOINTS(BRAVAIS,NPX,NP,PNAME1,P1)
        ALLOCATE(PNAME(NP))
        ALLOCATE(P(3,NP))
        PNAME(:)=PNAME1(:NP)
        DO I=1,NP
          P(:,I)=MATMUL(GREF,P1(:,I))
        ENDDO
      ELSE
        NP=0
      END IF
      RETURN
      END     
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPACEGROUP$SYMMETRYPOINTS(BRAVAIS,NPX,NP,NAME,G)
!     **************************************************************************
!     ** NOT FINISHED!!!!                                                     **
!     ** HIGH-SYMMETRY POINTS IN THE BRILLOUIN ZONE IN RELATIVE COORDINATES   **
!     ** FOLLOWING BRADLEY CRACKNELL FIG.3.2-FIG 3.15                         **
!     **                                                                      **
!     ** GAMMA IS DENOTE G                                                    **
!     ** AN ELEMENT '/' IN PATH IDENTIFIES A INTERRUPTION OF THE PATH         **
!     **                                                                      **
!     ** THE PATH IS SELECTED ACCORDING TO                                    **
!     **                HTTP://EN.WIKIPEDIA.ORG/WIKI/BRILLOUIN_ZONE           **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: BRAVAIS
      INTEGER(4)  ,INTENT(IN) :: NPX
      INTEGER(4)  ,INTENT(OUT):: NP
      CHARACTER(*),INTENT(OUT):: NAME(NPX)
      REAL(8)     ,INTENT(OUT):: G(3,NPX)
      CHARACTER(1)            :: PATH(20) ! PATH FROM HIGH SYMMETRY POINTS
!     **************************************************************************
      NP=1;  NAME='G' ; G(:,NP)=(/0.D0,0.D0,0.D0/)
      IF(BRAVAIS.EQ.'GT') THEN
        NP=NP+1; NAME(NP)='B';  G(:,NP)=(/0.5D0,0.D0,0.D0/)
        NP=NP+1; NAME(NP)='F';  G(:,NP)=(/0.D0,0.5D0,0.D0/)
        NP=NP+1; NAME(NP)='G';  G(:,NP)=(/0.D0,0.0D0,0.5D0/)
      ELSE IF(BRAVAIS.EQ.'GM') THEN
        NP=NP+1; NAME(NP)='B';  G(:,NP)=(/-0.5D0,0.D0,0.D0/)
        NP=NP+1; NAME(NP)='Y';  G(:,NP)=(/0.D0,0.5D0,0.D0/)
        NP=NP+1; NAME(NP)='Z';  G(:,NP)=(/0.D0,0.0D0,0.5D0/)
        NP=NP+1; NAME(NP)='C';  G(:,NP)=(/0.D0,0.5D0,0.5D0/)
        NP=NP+1; NAME(NP)='D';  G(:,NP)=(/-0.5D0,0.D0,0.5D0/)
        NP=NP+1; NAME(NP)='A';  G(:,NP)=(/-0.5D0,0.5D0,0.D0/)
        NP=NP+1; NAME(NP)='E';  G(:,NP)=(/-0.5D0,0.5D0,0.5D0/)
      ELSE IF(BRAVAIS.EQ.'GMB') THEN
        NP=NP+1; NAME(NP)='A';  G(:,NP)=(/-0.5D0,0.D0,0.D0/)
        NP=NP+1; NAME(NP)='Z';  G(:,NP)=(/0.D0,-0.5D0,0.5D0/)
        NP=NP+1; NAME(NP)='M';  G(:,NP)=(/-0.5D0,-0.5D0,0.5D0/)
        NP=NP+1; NAME(NP)='L';  G(:,NP)=(/-0.5D0,0.D0,0.5D0/)
        NP=NP+1; NAME(NP)='V';  G(:,NP)=(/0.D0,0.D0,0.5D0/)
      ELSE IF(BRAVAIS.EQ.'GO') THEN
        NP=NP+1; NAME(NP)='Y';  G(:,NP)=(/-0.5D0,0.D0,0.D0/)
        NP=NP+1; NAME(NP)='X';  G(:,NP)=(/0.D0,0.5D0,0.D0/)
        NP=NP+1; NAME(NP)='Z';  G(:,NP)=(/0.D0,0.D0,0.5D0/)
        NP=NP+1; NAME(NP)='U';  G(:,NP)=(/0.D0,0.5D0,0.5D0/)
        NP=NP+1; NAME(NP)='T';  G(:,NP)=(/-0.5D0,0.D0,0.5D0/)
        NP=NP+1; NAME(NP)='S';  G(:,NP)=(/-0.5D0,0.5D0,0.D0/)
        NP=NP+1; NAME(NP)='R';  G(:,NP)=(/-0.5D0,0.5D0,0.5D0/)




      ELSE IF(BRAVAIS.EQ.'GC') THEN
        NP=NP+1; NAME(NP)='X';  G(:,NP)=(/0.D0,0.5D0,0.D0/)
        NP=NP+1; NAME(NP)='M';  G(:,NP)=(/0.5D0,0.5D0,0.D0/)
        NP=NP+1; NAME(NP)='R';  G(:,NP)=(/0.5D0,0.5D0,0.5D0/)
        PATH(1:9)=(/'G','X','M','G','R','X','/','M','R'/)
      ELSE IF(BRAVAIS.EQ.'GCF') THEN
        NP=NP+1; NAME(NP)='X';  G(:,NP)=(/0.5D0,0.D0,0.5D0/)
        NP=NP+1; NAME(NP)='L';  G(:,NP)=(/0.5D0,0.5D0,0.5D0/)
        NP=NP+1; NAME(NP)='W';  G(:,NP)=(/0.5D0,0.25D0,0.75D0/)
        PATH(1:13)=(/'G','X','W','K','G','L','U','W','L','K','/','U','X'/)
      ELSE IF(BRAVAIS.EQ.'GCV') THEN
        NP=NP+1; NAME(NP)='H';  G(:,NP)=(/0.5D0,-0.5D0,0.5D0/)
        NP=NP+1; NAME(NP)='P';  G(:,NP)=(/0.25D0,0.25D0,0.25D0/)
        NP=NP+1; NAME(NP)='N';  G(:,NP)=(/0.D0,0.D0,0.5D0/)
        PATH(1:7)=(/'G','H','N','G','/','P','N'/)
      ELSE
        CALL ERROR$MSG('ID FOR BRAVAIS LATTICE NOT RECOGNIZED')
        CALL ERROR$MSG('SPACEGROUP$SYMMETRYPOINTS')
      END IF
      RETURN
      END     
        
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: PDOSINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
!     ******************************************************************
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE PDOS TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,PDOSINNAME)
      ISVAR=INDEX(PDOSINNAME,-'.BCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=PDOSINNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES_BANDS
      CALL FILEHANDLER$SETFILE('BCNTL',.FALSE.,PDOSINNAME)
      RETURN
      END SUBROUTINE INITIALIZEFILEHANDLER
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STANDARDFILES_BANDS
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: ID
!     **************************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES_BANDS')
!  
!     ==================================================================
!     == SET STANDARD FILENAMES                                       ==
!     ==================================================================
!
!     ==  ERROR FILE ===================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.BERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOL FILE ================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.BPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =============================================
      ID=+'BCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.BCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  STRUCTURE FILE   =============================================
      ID=+'PDOS'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PDOS')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES_BANDS
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION(LL_CNTL,NFIL,NFILO)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE PDOS_MODULE, ONLY: STATE,STATEARR
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT):: LL_CNTL
      INTEGER(4)   ,INTENT(IN)  :: NFIL
      INTEGER(4)   ,INTENT(IN)  :: NFILO
      INTEGER(4)                :: NAT
      INTEGER(4)                :: NB
      INTEGER(4)                :: NKPT
      INTEGER(4)                :: NKPTBIG
      INTEGER(4)                :: NSPIN
      INTEGER(4)                :: NDIM !=2 FOR SPINOR WF; OTHERWISE =1
      REAL(8)                   :: RBAS(3,3) ! REAL-SPACE LATTICE VECTORS
      REAL(8)                   :: GBAS(3,3) ! G-SPACE LATTICE VECTORS
      REAL(8)                   :: GBASINV(3,3) ! INVERTED G-SPACE LATTICE VECTORS
      REAL(8)                   :: GUNIT ! 2*PI/LUNIT
      REAL(8)      ,ALLOCATABLE :: RPOS(:,:)
      REAL(8)      ,ALLOCATABLE :: XK(:,:) ! K IN RELATIVE COORDINATES
      REAL(8)      ,ALLOCATABLE :: EB(:,:,:)
      REAL(8)      ,ALLOCATABLE :: PROPSI(:,:,:)
      REAL(8)      ,ALLOCATABLE :: EBBIG(:,:,:)
      INTEGER(4)   ,ALLOCATABLE :: MAP(:,:,:)
      INTEGER(4)   ,ALLOCATABLE :: MAPBIG(:,:,:)
      REAL(8)      ,ALLOCATABLE :: WORK1(:,:,:)
      REAL(8)      ,ALLOCATABLE :: WORK2(:,:,:)
      INTEGER(4)                :: NFILIN
      INTEGER(4)                :: NPRO
      INTEGER(4)                :: IKPT,ISPIN,IB,I,J,K,J1,J2
      INTEGER(4)                :: N1,N2,N3
      INTEGER(4)                :: N1B,N2B,N3B
      INTEGER(4)                :: IND
      INTEGER(4)   ,ALLOCATABLE :: NBARR(:,:)
      REAL(8)                   :: XK1(3),XK2(3),XQ(3)
      REAL(8)                   :: KVEC1(3),KVEC2(3)
      CHARACTER(512)            :: FILE
      INTEGER(4)                :: NLINE,ILINE
      INTEGER(4)                :: NQ
      INTEGER(4)                :: NP
      REAL(8)                   :: X1,X2
      REAL(8)                   :: SVAR,SVAR1
      LOGICAL(4)                :: TCHK,TCHK1,TCHK2
      LOGICAL(4)                :: TAPPEND
      INTEGER(4)                :: NBMAX   ! LIMITS THE NUMBER OF BANDS PLOTTED
      REAL(8)      ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)                   :: ANGSTROM
      LOGICAL(4),PARAMETER      :: TREFINE=.FALSE.
!     **************************************************************************
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
!
!     ==========================================================================
!     ==  READ PDOSFILE                                                       ==
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
      CALL PDOS$GETR8A('RBAS',3*3,RBAS)
      ALLOCATE(RPOS(3,NAT))
      CALL PDOS$GETR8A('R',3*NAT,RPOS)
      ALLOCATE(EB(NB,NKPT,NSPIN))
      ALLOCATE(PROPSI(NB,NKPT,NSPIN))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          DO IB=1,NB
            EB(IB,IKPT,ISPIN)=STATE%EIG(IB)
            !STATE%VEC(NDIM,NPRO,IB)
            PROPSI(IB,IKPT,ISPIN)=ABS(STATE%VEC(1,1,IB))
          ENDDO
        ENDDO
      ENDDO
      ALLOCATE(XK(3,NKPT))
      CALL PDOS$GETR8A('XK',3*NKPT,XK)
      CALL GBASS(RBAS,GBAS,SVAR)
                            CALL TRACE$PASS('AFTER READPDOS')
!
!!$PRINT*,'RBAS ',RBAS
!!$PRINT*,'GBAS ',GBAS
!!$PRINT*,'NKPT ',NKPT
!
!      =========================================================================
!      ==  EXTRACT DIVISION OF LATTICE VECTORS                                ==
!      =========================================================================
       DO I=1,3
         SVAR=1.D0
         DO J1=1,NKPT
           DO J2=J1+1,NKPT
             SVAR1=ABS(MODULO(XK(I,J1)-XK(I,J2)+0.5D0,1.D0)-0.5D0)
             IF(SVAR1.NE.0.D0) SVAR=MIN(SVAR,SVAR1)
           ENDDO
         ENDDO
         SVAR=1.D0/SVAR
         IF(I.EQ.1) THEN
           N1=NINT(SVAR)
!!$PRINT*,'N1 ',N1,SVAR
         ELSE IF(I.EQ.2) THEN
           N2=NINT(SVAR)
!!$PRINT*,'N2 ',N2,SVAR
         ELSE IF(I.EQ.3) THEN
           N3=NINT(SVAR)
!!$PRINT*,'N3 ',N3,SVAR
         END IF
       ENDDO
!
!      =========================================================================
!      ==  CONSTRUCT MAPPING ONTO K-POINT ARRAY                               ==
!      =========================================================================
       ALLOCATE(MAP(N1,N2,N3))       
       CALL BANDS_GETMAP(GBAS,N1,N2,N3,NKPT,XK,MAP)
!
!      =========================================================================
!      ==  INTERPOLATE TO A FINER MESH                                        ==
!      =========================================================================
       IF(TREFINE) THEN
         N1B=32
         N2B=32
         N3B=32
         CALL LIB$FFTADJUSTGRD(N1B)
         CALL LIB$FFTADJUSTGRD(N2B)
         CALL LIB$FFTADJUSTGRD(N3B)
PRINT*,'N1B ',N1B,N2B,N3B
         ALLOCATE(MAPBIG(N1B,N2B,N3B))
         IND=0
         DO K=1,N3B
           DO J=1,N2B
             DO I=1,N1B
               IND=IND+1
               MAPBIG(I,J,K)=IND
             ENDDO
           ENDDO
         ENDDO    
         NKPTBIG=IND
         ALLOCATE(EBBIG(NB,NKPTBIG,NSPIN))
         ALLOCATE(WORK1(N1,N2,N3))
         ALLOCATE(WORK2(N1B,N2B,N3B))
         DO ISPIN=1,NSPIN
           DO IB=1,NB
             DO I=1,N1
               DO J=1,N2
                 DO K=1,N3
                   WORK1(I,J,K)=EB(IB,MAP(I,J,K),ISPIN)
                 ENDDO
               ENDDO
             ENDDO    
             CALL REFINEGRID(N1,N2,N3,N1B,N2B,N3B,WORK1,WORK2)
             DO I=1,N1B
               DO J=1,N2B
                 DO K=1,N3B
                   EBBIG(IB,MAPBIG(I,J,K),ISPIN)=WORK2(I,J,K)
                 ENDDO
               ENDDO
             ENDDO    
           ENDDO
         ENDDO  
         DEALLOCATE(WORK1)
         DEALLOCATE(WORK2)
       END IF
!
!      =========================================================================
!      ==  CONSTRUCT BAND STRUCTURE                                           ==
!      =========================================================================
!      == DEFAULT OUTPUT FILE ==================================================
       CALL FILEHANDLER$SETFILE('BANDS',.TRUE.,-'.BANDS')
       CALL FILEHANDLER$SETSPECIFICATION('BANDS','STATUS','REPLACE')
       CALL FILEHANDLER$SETSPECIFICATION('BANDS','POSITION','APPEND')
       CALL FILEHANDLER$SETSPECIFICATION('BANDS','ACTION','WRITE')
       CALL FILEHANDLER$SETSPECIFICATION('BANDS','FORM','FORMATTED')
!
       CALL LINKEDLIST$SELECT(LL_CNTL,'~')
       CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
       CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
       CALL LINKEDLIST$NLISTS(LL_CNTL,'LINE',NLINE)
       X2=0.D0
       DO ILINE=1,NLINE
         WRITE(NFILO,FMT='(72("="))')
         WRITE(NFILO,*)'LINE-BLOCK:',ILINE,' OF ',NLINE
         WRITE(NFILO,FMT='(72("="))')
!
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
         CALL LINKEDLIST$SELECT(LL_CNTL,'LINE',ILINE)
!
         CALL LINKEDLIST$EXISTD(LL_CNTL,'TAPPEND',0,TCHK)
         IF(TCHK) THEN
           CALL LINKEDLIST$GET(LL_CNTL,'TAPPEND',1,TAPPEND)
         ELSE
           TAPPEND=.FALSE.
         ENDIF
         CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',0,TCHK)
         IF(TCHK) THEN
           CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FILE)
           CALL FILEHANDLER$SETFILE('BANDS',.FALSE.,FILE)
           IF(TAPPEND)THEN
             WRITE(NFILO,*)'ATTACHING OUTPUT FOR LINE-BLOCK',ILINE &
      &                   ,' TO EXISTING FILE ',TRIM(FILE)
             CALL FILEHANDLER$SETSPECIFICATION('BANDS','STATUS','OLD')
           ELSE
             WRITE(NFILO,*)'WRITING OUTPUT FOR LINE-BLOCK',ILINE &
      &                   ,' TO NEW FILE ',TRIM(FILE)
             CALL FILEHANDLER$SETSPECIFICATION('BANDS','STATUS','REPLACE')
           ENDIF
           CALL FILEHANDLER$SETSPECIFICATION('BANDS','POSITION','APPEND')
           CALL FILEHANDLER$SETSPECIFICATION('BANDS','ACTION','WRITE')
           CALL FILEHANDLER$SETSPECIFICATION('BANDS','FORM','FORMATTED')
           WRITE(NFILO,*)'OUTPUT FILE: ',TRIM(FILE)
         ELSE
           CALL ERROR$MSG('NO OUTPUT FILE GIVEN')
           CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
           CALL ERROR$STOP('BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION')
         END IF
!
!        =======================================================================
!        ==  READ K-VECTORS                                                   ==
!        =======================================================================
!
!        == SCALE FACTOR =======================================================
         CALL LINKEDLIST$EXISTD(LL_CNTL,'LUNIT',0,TCHK1)
         CALL LINKEDLIST$EXISTD(LL_CNTL,'LUNIT[AA]',0,TCHK2)
         IF(TCHK1.AND.TCHK2) THEN
           CALL ERROR$MSG('!BCNTL!BANDSTRUCTURE!LINE:LUNIT AND LUNIT[AA] ')
           CALL ERROR$MSG('MUST NOT BE SPECIFIED SIMULTANEOUSLY')
           CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
           CALL ERROR$STOP('BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION')
         ELSE           
           GUNIT=1.D0
           IF(TCHK1) THEN
             CALL LINKEDLIST$GET(LL_CNTL,'LUNIT',1,SVAR)
             GUNIT=2.D0*PI/SVAR
           ELSE IF(TCHK2) THEN
             CALL LINKEDLIST$GET(LL_CNTL,'LUNIT[AA]',1,SVAR)
             SVAR=SVAR*ANGSTROM
             GUNIT=2.D0*PI/SVAR
           ELSE
!            == THIS OPTION IS OBSOLETE AND IS THERE FOR BACKWARD COMPATIBILITY
             CALL LINKEDLIST$EXISTD(LL_CNTL,'KVECSCALE',0,TCHK)
             IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'KVECSCALE',1,SVAR)
               GUNIT=2.D0*PI/SVAR
             END IF
           END IF
         END IF
!
         CALL LIB$INVERTR8(3,GBAS,GBASINV) 

         CALL LINKEDLIST$EXISTD(LL_CNTL,'XK1',0,TCHK1)
         CALL LINKEDLIST$EXISTD(LL_CNTL,'KVEC1',0,TCHK2)
         IF(TCHK1.AND.TCHK2) THEN
           CALL ERROR$MSG('!BCNTL!BANDSTRUCTURE!LINE:XK1 AND KVEC1')
           CALL ERROR$MSG('MUST NOT BE SPECIFIED SIMULTANEOUSLY')
           CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
           CALL ERROR$STOP('BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION')
         ELSE
           IF(TCHK1)THEN
             CALL LINKEDLIST$GET(LL_CNTL,'XK1',1,XK1)
             WRITE(NFILO,*)'FIRST K-POINT: GIVEN IN RELATIVE COORDINATES'
             KVEC1=MATMUL(GBAS,XK1)
           ELSE IF(TCHK2)THEN
             CALL LINKEDLIST$GET(LL_CNTL,'KVEC1',1,KVEC1)
             WRITE(NFILO,*)'FIRST K-POINT: GIVEN IN ABSOLUTE COORDINATES'
             KVEC1=GUNIT*KVEC1
             XK1=MATMUL(GBASINV,KVEC1)
           ELSE
             CALL ERROR$MSG('NO XK1 OR KVEC1 GIVEN')
             CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
             CALL ERROR$STOP('BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION')
           ENDIF
         END IF

         CALL LINKEDLIST$EXISTD(LL_CNTL,'XK2',0,TCHK1)
         CALL LINKEDLIST$EXISTD(LL_CNTL,'KVEC2',0,TCHK2)
         IF(TCHK1)THEN
           CALL LINKEDLIST$GET(LL_CNTL,'XK2',1,XK2)
           WRITE(NFILO,*)'LAST K-POINT: GIVEN IN RELATIVE COORDINATES'
           KVEC2=MATMUL(GBAS,XK2)
         ELSE IF (TCHK2)THEN
           CALL LINKEDLIST$GET(LL_CNTL,'KVEC2',1,KVEC2)
           WRITE(NFILO,*)'LAST K-POINT: GIVEN IN ABSOLUTE COORDINATES'
           KVEC2=KVEC2*GUNIT
           XK2=MATMUL(GBASINV,KVEC2)
         ELSE
           CALL ERROR$MSG('NO XK2 OR KVEC2 GIVEN')
           CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
           CALL ERROR$STOP('BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION')
         ENDIF

         WRITE(NFILO,*)'FIRST K-POINT IN RELATIVE COORDINATES:',XK1
         WRITE(NFILO,*)'FIRST K-POINT IN ABSOLUTE COORDINATES:',KVEC1 
         WRITE(NFILO,*)'LAST K-POINT IN RELATIVE COORDINATES:',XK2
         WRITE(NFILO,*)'LAST K-POINT IN ABSOLUTE COORDINATES:',KVEC2
         IF(TAPPEND)THEN
           X1=X2
           X2=X2+SQRT(SUM((KVEC2-KVEC1)**2))/GUNIT
         ELSE
           X1=0.0D0
           X2=SQRT(SUM((KVEC2-KVEC1)**2))/GUNIT
         ENDIF
         
         NP=100
         CALL LINKEDLIST$EXISTD(LL_CNTL,'NK',0,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NK',1,NP)
!
!        == NBMAX LIMITS THE NUMBER OF BANDS PRINTED ===========================
         NBMAX=20
         CALL LINKEDLIST$EXISTD(LL_CNTL,'NB',0,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NB',1,NBMAX)

!        == PROJECTION FOR 2D- BANDSTRUCTURE ===================================
         CALL LINKEDLIST$EXISTD(LL_CNTL,'XKPROJECT',0,TCHK)
         IF(TCHK) THEN
           CALL LINKEDLIST$GET(LL_CNTL,'XKPROJECT',1,XQ)
           NQ=11
           CALL LINKEDLIST$EXISTD(LL_CNTL,'NPROJECT',0,TCHK1)
           IF(TCHK1)CALL LINKEDLIST$GET(LL_CNTL,'NPROJECT',1,NQ)
           IF(MODULO(N1,2).NE.1) THEN
!            == FOR 2-D- BAND STRUCTURES ONE GOES FORWARD-BACKWARD-FORWARD =====
!            == FOR EVEN NQ A STRAIGHT LINE IS DRAWN TOWARDS THE NEXT SEGMENT ==
!            == WHEN SEVERAL LINES ARE WRITTEN IN THE SAME FILE=================
             CALL ERROR$MSG('VARIABLE !LINE:NPROJECT MUST BE AN ODD NUMBER')
             CALL ERROR$I4VAL('NPROJECT',NQ)
             CALL ERROR$STOP('MAIN')
           END IF
         ELSE
           XQ(:)=0.D0
           NQ=1
         END IF
!        == SPECIFY SPIN DIRECTION =============================================
         ISPIN=1
         CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',0,TCHK)
         IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,ISPIN)
         IF(ISPIN.GT.NSPIN) THEN
           CALL ERROR$MSG('ISPIN EXCEEDS RANGE')
           CALL ERROR$STOP('MAIN')
         END IF
!
!        == PLOT BANDS =========================================================
         IF(TREFINE) THEN         
!           CALL BANDS_PLOTBANDS(GBAS,N1B,N2B,N3B,MAPBIG &
!      &                        ,NKPTBIG,NB,EBBIG(:,:,ISPIN) &
!      &                        ,NP,X1,XK1,X2,XK2,NQ,XQ)
         ELSE
           CALL BANDS_PLOTBANDS(GBAS,N1,N2,N3,MAP,NKPT &
      &                        ,NBMAX,EB(:NBMAX,:,ISPIN) &
      &                        ,PROPSI(:NBMAX,:,ISPIN) &
      &                        ,NP,X1,XK1,X2,XK2,NQ,XQ)
         END IF
!
         CALL LINKEDLIST$SELECT(LL_CNTL,'..')
       ENDDO
       IF(TREFINE) THEN
         PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         PRINT*,'CAUTION! THE EXPERIMENTAL OPTION TREFINE HAS BEEN USED'
         PRINT*,'THIS IMPLEMENTATION IS NOT READY FOR USE'
         PRINT*,'THE INTERPOLATION TOA FINE GRID DOES NOT WORK'
         PRINT*,'PROBABLY ALSO THE BAND CROSSINGS NEED TO BE ACCOUNTED FOR' 
         PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       END IF
       RETURN
       END SUBROUTINE BANDS_BANDSTRUCTURE_LINEAR_INTERPOLATION
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE BANDS_PLOTBANDS(GBAS,N1,N2,N3,MAP,NK,NB,EB,PROPSI &
      &                          ,NP,X1,XK1,X2,XK2,NQ,XQ)
!      *************************************************************************
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: GBAS(3,3)
       INTEGER(4),INTENT(IN) :: N1,N2,N3
       INTEGER(4),INTENT(IN) :: NK
       INTEGER(4),INTENT(IN) :: NB
       REAL(8)   ,INTENT(IN) :: EB(NB,NK)
       REAL(8)   ,INTENT(IN) :: PROPSI(NB,NK)
       INTEGER(4),INTENT(IN) :: MAP(N1,N2,N3)
       INTEGER(4),INTENT(IN) :: NP
       INTEGER(4),INTENT(IN) :: NQ
       REAL(8)   ,INTENT(IN) :: XK1(3)  !INITIAL K-POINT IN RELATIVE COORDINATES
       REAL(8)   ,INTENT(IN) :: XK2(3)  !FINAL K-POINT IN RELATIVE COORDINATES
       REAL(8)   ,INTENT(IN) :: XQ(3)   !PROJECTION DIRECTION IN RELATIVE COORD.
       REAL(8)   ,INTENT(IN) :: X1
       REAL(8)   ,INTENT(IN) :: X2
       REAL(8)               :: D1,D2,D3
       REAL(8)               :: KI(3)
       REAL(8)               :: EBI(NB)
       REAL(8)               :: PROPSII(NB)
       INTEGER(4)            :: IP,IQ
       INTEGER(4)            :: NFIL=11
       CHARACTER(64)         :: FMTSTRING
!      *************************************************************************
                                            CALL TRACE$PUSH('BANDS_PLOTBANDS')
       PRINT*,'NB ',NB
       PRINT*,'X1X2 ',X1,X2
       WRITE(FMTSTRING,*)NB+1
       FMTSTRING='('//TRIM(ADJUSTL(FMTSTRING))//'F10.5)'
       CALL FILEHANDLER$UNIT('BANDS',NFIL)
       DO IQ=1,NQ
         IF(NQ.EQ.1) THEN
           D3=0.D0
         ELSE
           D3=REAL(IQ-1)/REAL(NQ-1)
         END IF
         DO IP=1,NP
           D2=REAL(IP-1)/REAL(NP-1)
           IF((-1)**IQ.GT.0)D2=1.D0-D2
           D1=1.D0-D2
!          == KI IS THE K-POINT IN CARTESIAN COORDINATES =======================
           KI=MATMUL(GBAS,XK1*D1+XK2*D2+XQ*D3)
           CALL BANDS_GETBAND(GBAS,N1,N2,N3,MAP,NK,NB,EB,PROPSI,KI,EBI,PROPSII)
           WRITE(NFIL,FMT=FMTSTRING)X1+(X2-X1)*D2,EBI 
         ENDDO
       ENDDO
       FLUSH(NFIL)
       CALL FILEHANDLER$CLOSE('BANDS')
                                            CALL TRACE$POP()
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE BANDS_GETMAP(GBAS,N1,N2,N3,NK,XK,MAP)
!      *************************************************************************
!      ** MAPPING FROM THE K-POINT GRID TO THE POSITION IN THE K-POINT ARRAY  **
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: GBAS(3,3)
       INTEGER(4),INTENT(IN) :: N1,N2,N3
       INTEGER(4),INTENT(IN) :: NK
       REAL(8)   ,INTENT(IN) :: XK(3,NK)
       INTEGER(4),INTENT(OUT):: MAP(N1,N2,N3)
       INTEGER(4)            :: IK,I1,I2,I3
!      *************************************************************************
       MAP(:,:,:)=-1
       DO IK=1,NK
         I1=NINT(XK(1,IK)*REAL(N1,KIND=8))
         I2=NINT(XK(2,IK)*REAL(N2,KIND=8))
         I3=NINT(XK(3,IK)*REAL(N3,KIND=8))
         I1=MODULO(I1,N1)
         I2=MODULO(I2,N2)
         I3=MODULO(I3,N3)
         MAP(I1+1,I2+1,I3+1)=IK
         I1=MODULO(-I1,N1)
         I2=MODULO(-I2,N2)
         I3=MODULO(-I3,N3)
         MAP(I1+1,I2+1,I3+1)=IK
       ENDDO
!
!      == TEST IF ALL SLOTS ARE OCCUPIED
       DO I1=1,N1
         DO I2=1,N2
           DO I3=1,N3
             IF(MAP(I1,I2,I3).EQ.-1) THEN
               STOP 'MAP IS INCOMPLETE '
             END IF
           ENDDO
         ENDDO
       ENDDO
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE BANDS_GETBAND(GBAS,N1,N2,N3,MAP,NK,NB,EB,PROPSI,K1,EB1,PROPSI1)
!      *************************************************************************
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: GBAS(3,3)
       INTEGER(4),INTENT(IN) :: N1,N2,N3
       INTEGER(4),INTENT(IN) :: MAP(N1,N2,N3)
       INTEGER(4),INTENT(IN) :: NK
       INTEGER(4),INTENT(IN) :: NB
       REAL(8)   ,INTENT(IN) :: EB(NB,NK)
       REAL(8)   ,INTENT(IN) :: PROPSI(NB,NK)
       REAL(8)   ,INTENT(IN) :: K1(3)
       REAL(8)   ,INTENT(OUT):: EB1(NB)
       REAL(8)   ,INTENT(OUT):: PROPSI1(NB)
       REAL(8)               :: GBASIN(3,3)
       REAL(8)               :: XK(3)
       INTEGER(4)            :: I1,I2,I3
       INTEGER(4)            :: I,J,K,N
       INTEGER(4)            :: I1A,I2A,I3A
       INTEGER(4)            :: IND
       INTEGER(4)            :: CORNER_MAP(4)
       REAL(8)               :: CORNER_XK(3,4)
       INTEGER(4)            :: IND4(4)
       REAL(8)               :: AMAT(3,3),AMATIN(3,3)
       REAL(8)               :: WGHT(4)
       REAL(8)               :: S1,S2,S3
!      *************************************************************************
!       CALL LIB__INVERTR8(3,GBAS,GBASIN)
       GBASIN(1,1)=GBAS(2,2)*GBAS(3,3)-GBAS(2,3)*GBAS(3,2)
       GBASIN(2,1)=GBAS(3,2)*GBAS(1,3)-GBAS(3,3)*GBAS(1,2)
       GBASIN(3,1)=GBAS(1,2)*GBAS(2,3)-GBAS(1,3)*GBAS(2,2)
       GBASIN(1,2)=GBAS(2,3)*GBAS(3,1)-GBAS(2,1)*GBAS(3,3)
       GBASIN(2,2)=GBAS(3,3)*GBAS(1,1)-GBAS(3,1)*GBAS(1,3)
       GBASIN(3,2)=GBAS(1,3)*GBAS(2,1)-GBAS(1,1)*GBAS(2,3)
       GBASIN(1,3)=GBAS(2,1)*GBAS(3,2)-GBAS(2,2)*GBAS(3,1)
       GBASIN(2,3)=GBAS(3,1)*GBAS(1,2)-GBAS(3,2)*GBAS(1,1)
       GBASIN(3,3)=GBAS(1,1)*GBAS(2,2)-GBAS(1,2)*GBAS(2,1)
       GBASIN(:,:)=GBASIN(:,:)/DOT_PRODUCT(GBAS(:,1),GBASIN(:,1))
       GBASIN=TRANSPOSE(GBASIN)
!
       XK=MATMUL(GBASIN,K1) ! RELATIVE COORDINATE IN K-SPACE
!      == MAP INTO FIRST UNIT CELL =============================================
       XK(1)=MODULO(XK(1),1.D0)*REAL(N1,KIND=8) 
       XK(2)=MODULO(XK(2),1.D0)*REAL(N2,KIND=8) 
       XK(3)=MODULO(XK(3),1.D0)*REAL(N3,KIND=8) 
!      == SELECT LEFT-LOWER-FRONT (LLF) CORNER CONTAINING THE K-POINT
       I1=INT(XK(1))
       I2=INT(XK(2))
       I3=INT(XK(3))
!      == DISPLACEMENT RELATIVE TO THE LLF CORNER
       XK(1)=XK(1)-REAL(I1,KIND=8)
       XK(2)=XK(2)-REAL(I2,KIND=8)
       XK(3)=XK(3)-REAL(I3,KIND=8)
!
!      =========================================================================
!      ==  SELECT THE FOUR CLOSEST CORNERS  SPANNING A TETRAHEDRON            ==
!      ==    1=(0,0,0)   2=(1,0,0)  3=(0,1,0)  4=(1,1,0)                      ==
!      ==    5=(0,0,1)   6=(1,0,1)  7=(0,1,1)  8=(1,1,1)                      ==
!      =========================================================================
       S1=XK(2)-XK(3)
       S2=XK(3)-XK(1)       
       S3=XK(1)-XK(2)       
       IND4(1)=1
       IND4(2)=8
       IF(S1.GE.0.AND.S2.LE.0.AND.S3.GE.0) THEN
         IND4(3)=2
         IND4(4)=4
       ELSE IF(S1.GE.0.AND.S2.LE.0.AND.S3.LE.0) THEN
         IND4(3)=4
         IND4(4)=3
       ELSE IF(S1.GE.0.AND.S2.GE.0.AND.S3.LE.0) THEN
         IND4(3)=3
         IND4(4)=7
       ELSE IF(S1.LE.0.AND.S2.GE.0.AND.S3.LE.0) THEN
         IND4(3)=7
         IND4(4)=5
       ELSE IF(S1.LE.0.AND.S2.GE.0.AND.S3.GE.0) THEN
         IND4(3)=5
         IND4(4)=6
       ELSE IF(S1.LE.0.AND.S2.LE.0.AND.S3.GE.0) THEN
         IND4(3)=6
         IND4(4)=2
       ELSE
         STOP 'ERROR!!'
       END IF
!
!      == COLLECT THE FOUR CORNERS =============================================
       IND=0
       DO K=0,1
         DO J=0,1
           DO I=0,1
             IND=IND+1
             DO N=1,4
               IF(IND4(N).NE.IND) CYCLE
               I1A=1+MODULO(I1+I,N1)
               I2A=1+MODULO(I2+J,N2)
               I3A=1+MODULO(I3+K,N3)
               CORNER_MAP(N)=MAP(I1A,I2A,I3A)
               CORNER_XK(1,N)=(REAL(I,KIND=8)-XK(1))/REAL(N1,KIND=8) 
               CORNER_XK(2,N)=(REAL(J,KIND=8)-XK(2))/REAL(N2,KIND=8)
               CORNER_XK(3,N)=(REAL(K,KIND=8)-XK(3))/REAL(N3,KIND=8)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
!
!      =========================================================================
!      ==  DETERMINE INTERPOLATION WEIGHTS                                    ==
!      =========================================================================
       DO J=1,3
         AMAT(:,J)=CORNER_XK(:,J+1)-CORNER_XK(:,1)
       ENDDO
       AMATIN(1,1)=AMAT(2,2)*AMAT(3,3)-AMAT(2,3)*AMAT(3,2)
       AMATIN(2,1)=AMAT(3,2)*AMAT(1,3)-AMAT(3,3)*AMAT(1,2)
       AMATIN(3,1)=AMAT(1,2)*AMAT(2,3)-AMAT(1,3)*AMAT(2,2)
       AMATIN(1,2)=AMAT(2,3)*AMAT(3,1)-AMAT(2,1)*AMAT(3,3)
       AMATIN(2,2)=AMAT(3,3)*AMAT(1,1)-AMAT(3,1)*AMAT(1,3)
       AMATIN(3,2)=AMAT(1,3)*AMAT(2,1)-AMAT(1,1)*AMAT(2,3)
       AMATIN(1,3)=AMAT(2,1)*AMAT(3,2)-AMAT(2,2)*AMAT(3,1)
       AMATIN(2,3)=AMAT(3,1)*AMAT(1,2)-AMAT(3,2)*AMAT(1,1)
       AMATIN(3,3)=AMAT(1,1)*AMAT(2,2)-AMAT(1,2)*AMAT(2,1)
       AMATIN(:,:)=AMATIN(:,:)/DOT_PRODUCT(AMAT(:,1),AMATIN(:,1))
       AMATIN=TRANSPOSE(AMATIN)
       WGHT(2:4)=-MATMUL(AMATIN,CORNER_XK(:,1))
       WGHT(1)=1.D0-SUM(WGHT(2:4))
!
!      =========================================================================
!      ==  INTERPOLATE BANDS                                                  ==
!      =========================================================================
       EB1(:)=0.D0
       PROPSI1(:)=0.D0
       DO I=1,4
         EB1(:)=EB1(:)+WGHT(I)*EB(:,CORNER_MAP(I))*27.211D0
         PROPSI1(:)=PROPSI1(:)+WGHT(I)*PROPSI(:,CORNER_MAP(I))
       ENDDO       
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE REFINEGRID(NR1,NR2,NR3,NR1B,NR2B,NR3B,WAVE,WAVEBIG)
!     **************************************************************************
!     ** COPY OF GRAPHICS_REFINEGRID.
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: NR1
      INTEGER(4),INTENT(IN)     :: NR2
      INTEGER(4),INTENT(IN)     :: NR3
      INTEGER(4),INTENT(IN)     :: NR1B
      INTEGER(4),INTENT(IN)     :: NR2B
      INTEGER(4),INTENT(IN)     :: NR3B
      REAL(8)   ,INTENT(IN)     :: WAVE(NR1,NR2,NR3)
      REAL(8)   ,INTENT(OUT)    :: WAVEBIG(NR1B,NR2B,NR3B)
      COMPLEX(8),ALLOCATABLE    :: WORKC1(:,:,:)
      COMPLEX(8),ALLOCATABLE    :: WORKC2(:,:,:)
      INTEGER(4)                :: I,J,K,IM,JM,KM
!     **************************************************************************
                           CALL TRACE$PUSH('REFINEGRID')
!
!     ==========================================================================
!     ==  FOURIER TRANSFORM INTO G-SPACE                                      ==
!     ==========================================================================
      ALLOCATE(WORKC1(NR1,NR2,NR3))
      WORKC1(:,:,:)=CMPLX(WAVE(:,:,:),KIND=8)
      ALLOCATE(WORKC2(NR1,NR2,NR3)) 
      CALL LIB$3DFFTC8('RTOG',NR1,NR2,NR3,WORKC1,WORKC2)
      DEALLOCATE(WORKC1)
!
!     ==========================================================================
!     ==  EXPAND ARRAY AND PAD HIGH-G WITH ZEROS                              ==
!     ==========================================================================
      ALLOCATE(WORKC1(NR1B,NR2B,NR3B))      
      WORKC1=(0.D0,0.D0)
      I=NR1/2  !MIND: REQUIRES THAT ONLY EVEN NUMBERS ARE USED (-> ASSURED BY LIB$FFTADJUSTGRD)
      J=NR2/2
      K=NR3/2
      IM=NR1B-I+1
      JM=NR2B-J+1
      KM=NR3B-K+1
!
      WORKC1(1:I    ,1:J    ,1:K)    =WORKC2(1:I    ,1:J    ,1:K)
      WORKC1(IM:NR1B,1:J    ,1:K)    =WORKC2(I+1:2*I,1:J    ,1:K)
      WORKC1(1:I    ,JM:NR2B,1:K)    =WORKC2(1:I    ,J+1:2*J,1:K)        
      WORKC1(1:I    ,1:J    ,KM:NR3B)=WORKC2(1:I    ,1:J    ,K+1:2*K)
      WORKC1(IM:NR1B,JM:NR2B,1:K)    =WORKC2(I+1:2*I,J+1:2*J,1:K)
      WORKC1(1:I    ,JM:NR2B,KM:NR3B)=WORKC2(1:I    ,J+1:2*J,K+1:2*K)        
      WORKC1(IM:NR1B,1:J    ,KM:NR3B)=WORKC2(I+1:2*I,1:J    ,K+1:2*K)
      WORKC1(IM:NR1B,JM:NR2B,KM:NR3B)=WORKC2(I+1:2*I,J+1:2*J,K+1:2*K)
      DEALLOCATE(WORKC2)
!
!     ==========================================================================
!     ==  FOURIER BACKTRANSFORM                                               ==
!     ==========================================================================
      ALLOCATE(WORKC2(NR1B,NR2B,NR3B))
      CALL LIB$3DFFTC8('GTOR',NR1B,NR2B,NR3B,WORKC1,WORKC2)
      WAVEBIG(:,:,:)=REAL(WORKC2(:,:,:),KIND=8)
      DEALLOCATE(WORKC1)
      DEALLOCATE(WORKC2)
                                                 CALL TRACE$POP()
      RETURN
      END

      MODULE KPOINTDIAG_MODULE
      IMPLICIT NONE
      CONTAINS

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_KINDEP(NG,GVEC,G2,GWEIGHT,GIDG_PROTO,TI_H,TI_S)
!     **************************************************************************
!     **  CONSTRUCT K-INDEPENDENT PART OF HAMILTONIAN                         **
!     **************************************************************************
      USE STRINGS_MODULE
      USE BANDDATA_MODULE, ONLY : NRG,NDIM,RBAS,GBAS,VOFR,NR1,NR2,NR3,NSPIN &
     &                          ,TYPEID_PROTO,G1_PROTO,DEX_PROTO,NG_PROTO,EPW 
      USE WAVES_MODULE, ONLY : GSET_TYPE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)              :: NG
      REAL(8),ALLOCATABLE,INTENT(INOUT)   :: GVEC(:,:)   !G
      REAL(8),ALLOCATABLE,INTENT(INOUT)   :: G2(:)
      REAL(8),INTENT(OUT)                 :: GWEIGHT
      INTEGER(4),INTENT(OUT)              :: GIDG_PROTO
      COMPLEX(8),ALLOCATABLE,INTENT(INOUT):: TI_H(:,:,:)
      COMPLEX(8),ALLOCATABLE,INTENT(INOUT):: TI_S(:,:,:)
      REAL(8)                             :: KVEC(3)
      TYPE(GSET_TYPE)                     :: GSET
      COMPLEX(8),ALLOCATABLE              :: PSI1(:,:,:)
      COMPLEX(8),ALLOCATABLE              :: HPSI1(:,:,:)
      INTEGER(4)                          :: ISPIN,I,IDIM1,IDIM2
      LOGICAL(4)                          :: TPRINT=.TRUE.
      REAL(8)                             :: GMAX
      REAL(8)                             :: GDIFF(3,3,3)
      INTEGER(4)                          :: I1,I2,I3,J
      REAL(8)                             :: EPW2
      INTEGER(4)                          :: THISTASK,NTASKS,NGL
      INTEGER(4),ALLOCATABLE              :: ICOLOR(:)
      CHARACTER(512)                      :: FFTCID
      CHARACTER(512)                      :: FFTID
      LOGICAL(4),PARAMETER                :: TFIXSYMMETRYBREAK=.TRUE.
      INTEGER(4),SAVE                     :: COUNT=0
!     **************************************************************************
      COUNT=COUNT+1 
!
!     == GENERATE G-GRID =======================================================
      KVEC(1)=0.0D0
      KVEC(2)=0.0D0
      KVEC(3)=0.0D0
      
      IF(TFIXSYMMETRYBREAK)THEN
        IF(TPRINT)PRINT*,"EPW FROM INPUT=",EPW
        GMAX=SQRT(2.0D0*EPW)
        !FIND LONGEST DIAGONAL IN GSPACE GRID
        DO I1=-1,1
          DO I2=-1,1
            DO I3=-1,1
              GDIFF(I1+2,I2+2,I3+2)=SQRT(SUM((REAL(I1,KIND=8)*GBAS(:,1)&
      &              +REAL(I2,KIND=8)*GBAS(:,2)+REAL(I3,KIND=8)*GBAS(:,3))**2))
            ENDDO
          ENDDO
        ENDDO
        !IF(TPRINT)PRINT*,"LONGEST DIAGONAL=",MAXVAL(GDIFF(:,:,:)),MAXLOC(GDIFF(:,:,:))
        EPW2=0.5D0*(GMAX+MAXVAL(GDIFF(:,:,:)))**2
        IF(TPRINT)PRINT*,"EPW MODIFIED FOR PROPER CUTOFF=",EPW2
      ELSE
        EPW2=EPW
      ENDIF

      !CREATE COMMUNICATOR WITH ONE TASK
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      ALLOCATE(ICOLOR(NTASKS))
      ICOLOR(1:NTASKS)=(/ (I,I=1,NTASKS) /)
      WRITE(FFTCID,'(A,I4,A,I4)')"FFT ",COUNT," ",THISTASK
      CALL MPE$NEW('~',FFTCID,NTASKS,ICOLOR)
      CALL MPE$QUERY(FFTCID,NTASKS,THISTASK)
      
      WRITE(FFTID,'(A,I4)')"WAVE GAMMA ",COUNT
       
      CALL PLANEWAVE$INITIALIZE(FFTID,FFTCID,RBAS,KVEC,.FALSE.,EPW2 &
     &                               ,0,NR1,NR2,NR3)
      CALL PLANEWAVE$SELECT(FFTID)     

      CALL PLANEWAVE$GETI4('NGG',NG)
      IF(TPRINT)PRINT*,"EPW2 NG",EPW2,NG
      ALLOCATE(GVEC(3,NG))
      ALLOCATE(G2(NG))
      
      !COLLECT GVECTORS
      NGL=NG
      GSET%NGL=NG

      GVEC(:,:)=0.0D0
      G2(:)=0.0D0

      CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
      CALL PLANEWAVE$GETR8A('G2',NGL,G2)
      
      CALL GBASS(RBAS,GBAS,GWEIGHT)

!     == GENERATE G-GRID FOR PROJECTORS ========================================
      CALL RADIAL$NEW(TYPEID_PROTO,GIDG_PROTO)
      !DEX_PROTO=LOG(GMAX_PROTO/G1_PROTO)/REAL(NG_PROTO-1,KIND=8)
      CALL RADIAL$SETI4(GIDG_PROTO,'NR',NG_PROTO)
      CALL RADIAL$SETR8(GIDG_PROTO,'R1',G1_PROTO)
      CALL RADIAL$SETR8(GIDG_PROTO,'DEX',DEX_PROTO)

!     == ALLOCATE K-INDEPENDET ARRAYS ==========================================
      ALLOCATE(TI_H(NG*NDIM,NG*NDIM,NSPIN))
      ALLOCATE(TI_S(NG*NDIM,NG*NDIM,NSPIN))
      TI_H=0.0D0
      TI_S=0.0D0
!
!     ==========================================================================
!     == COMPUTE K-INDEPENDENT ARRAYS                                         ==
!     ==========================================================================
                            CALL TRACE$PUSH('COMPUTE K-INDEPENDENT ARRAYS')
      CALL TIMING$CLOCKON('PS_POTENTIAL')
      ALLOCATE(HPSI1(NG,NDIM,1))
      ALLOCATE(PSI1(NG,NDIM,1))
!
!     ==  EVALUATE MATRIX-ELEMENTS OF PS-POTENTIAL =============================
!FIXME: INTRODUCE G-G' AS OPTIMIZATION
      DO ISPIN=1,NSPIN
        DO I=1,NG
          DO IDIM1=1,NDIM
            I1=I*(NDIM)+(IDIM1-1)
            PSI1(:,:,:)=0.0D0
            PSI1(I,IDIM1,1)=CMPLX(1.D0,0.D0,KIND=8)
            HPSI1(:,:,:)=0.0D0
            CALL WAVES_VPSI(GSET,NGL,NDIM,1,NRG,PSI1(:,:,:) &
      &                                        ,VOFR(:,ISPIN),HPSI1(:,:,:))
            DO IDIM2=1,NDIM
              TI_H(I1,1+(IDIM2-1)*NG:IDIM2*NG,ISPIN)=HPSI1(:,IDIM2,1)
            ENDDO
            TI_S(I1,I1,ISPIN)=1.0D0
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(HPSI1)
      DEALLOCATE(PSI1)
      !REMOVE KINETIC ENERGY FROM TI_H
      DO I=1,NG
        DO IDIM1=1,NDIM
          I1=I*(NDIM)+(IDIM1-1)
          !EKIN
          TI_H(I1,I1,:)=TI_H(I1,I1,:)-0.5D0*G2(I)
        ENDDO
      ENDDO
      CALL TIMING$CLOCKOFF('PS_POTENTIAL')
                            CALL TRACE$POP
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_KPOINT(NG,NB,ISPIN,METHOD,GIDG_PROTO,TPROJ,&
     &               KVEC,GVEC,TI_H,TI_S,E,PROJ)
!     ******************************************************************
!     **                                                              **
!     ** COMPUTES EIGENVALUES[H] AND PROJECTIONS FOR A GIVEN KVEC     **
!     **                                                              **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NG
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: ISPIN
      INTEGER(4),INTENT(IN)  :: METHOD
      INTEGER(4),INTENT(IN)  :: GIDG_PROTO
      LOGICAL(4),INTENT(IN)  :: TPROJ
      REAL(8)   ,INTENT(IN)  :: KVEC(3)
      REAL(8)   ,INTENT(IN)  :: GVEC(3,NG)
      COMPLEX(8),INTENT(IN)  :: TI_H(NG*NDIM,NG*NDIM,NSPIN)
      COMPLEX(8),INTENT(IN)  :: TI_S(NG*NDIM,NG*NDIM,NSPIN)
      REAL(8)   ,INTENT(OUT),ALLOCATABLE :: E(:)
      COMPLEX(8),INTENT(OUT) :: PROJ(NAT,NB,LMNXX)
      COMPLEX(8),ALLOCATABLE :: TI_HK(:,:)
      COMPLEX(8),ALLOCATABLE :: TI_SK(:,:)
      COMPLEX(8),ALLOCATABLE :: U(:,:)
      REAL(8)                :: GVECPK(3,NG)
      REAL(8)                :: G2(NG)
      COMPLEX(8),ALLOCATABLE :: EIGR(:)
      COMPLEX(8),ALLOCATABLE :: FOFG(:)
      COMPLEX(8),ALLOCATABLE :: PRO(:,:,:)
      REAL(8)   ,ALLOCATABLE :: BAREPRO(:,:)
      REAL(8)   ,ALLOCATABLE :: YLM(:,:)
      REAL(8)                :: YLM_(LMX)
      INTEGER(4)             :: LOX_(LNXX)
      INTEGER(4)             :: ISP,IND,LN,IAT,IDIM1,IDIM2,IBPRO
      INTEGER(4)             :: LMN1,LMN2,IB
      INTEGER(4)             :: LNX_,LMNX_,IPRO
      INTEGER(4)             :: I,I1,J,K,L,M,N
      REAL(8)                :: CELLVOL,GWEIGHT
      REAL(8)                :: SVAR
      COMPLEX(8)             :: CSVAR,CSVAR_H,CSVAR_S
      LOGICAL(4)             :: TTIMING
      INTEGER(4)             :: NG2
      INTEGER(4)             :: NGARR(NG)
      LOGICAL(4),PARAMETER   :: TTEST_POSITIVE_DEFINITE=.FALSE.
      LOGICAL(4),PARAMETER   :: TFIXSYMMETRYBREAK=.TRUE.
!     ******************************************************************
      TTIMING=.FALSE.
      IF(TTIMING)CALL TRACE$PUSH('COMPUTE K_DEPENDENT ARRAY')
      
      CALL GBASS(RBAS,GBAS,GWEIGHT)
      CELLVOL=GWEIGHT

      !COMPUTE G+K AND (G+K)^2
      IF(TTIMING) CALL TIMING$CLOCKON('G+K,G2')
      
      IF(TFIXSYMMETRYBREAK)THEN
        IF(NDIM.NE.1)THEN
          CALL ERROR$MSG('NDIM.NE.1 NOT PROPERLY IMPLEMENTED YET')
          CALL ERROR$STOP('BANDS_KPOINT')
        ENDIF
        NG2=1
        DO I=1,NG
          IF(0.5D0*SUM((GVEC(:,I)+KVEC(:))**2).LE.EPW)THEN
            GVECPK(:,NG2)=GVEC(:,I)+KVEC(:)
            G2(NG2)=SUM(GVECPK(:,NG2)**2)
            NGARR(I)=NG2
            NG2=NG2+1
          ELSE
            NGARR(I)=0
          ENDIF
        ENDDO
        NG2=NG2-1

        ALLOCATE(TI_HK(NG2*NDIM,NG2*NDIM))
        ALLOCATE(TI_SK(NG2*NDIM,NG2*NDIM))
        N=1
        DO I=1,NG
          IF(NGARR(I).NE.0)THEN
            M=1
            DO J=1,NG
              IF(NGARR(J).NE.0)THEN
                TI_HK(N,M)=TI_H(I,J,ISPIN)
                TI_SK(N,M)=TI_S(I,J,ISPIN)
                M=M+1
              ENDIF
            ENDDO
            N=N+1
          ENDIF
        ENDDO
      ELSE
        DO I=1,NG
          GVECPK(:,I)=GVEC(:,I)+KVEC(:)
          G2(I)=SUM(GVECPK(:,I)**2)
        ENDDO
        NG2=NG
        ALLOCATE(TI_HK(NG2*NDIM,NG2*NDIM))
        ALLOCATE(TI_SK(NG2*NDIM,NG2*NDIM))
        TI_HK(1:NG2*NDIM,1:NG2*NDIM)=TI_H(1:NG2*NDIM,1:NG2*NDIM,ISPIN)
        TI_SK(1:NG2*NDIM,1:NG2*NDIM)=TI_S(1:NG2*NDIM,1:NG2*NDIM,ISPIN)
      ENDIF
      ALLOCATE(E(NG2*NDIM))    
      ALLOCATE(EIGR(NG2))      
      ALLOCATE(FOFG(NG2))      
      ALLOCATE(PRO(NAT,NG2,LMNXX))
      ALLOCATE(BAREPRO(NG2,NBAREPRO))
      ALLOCATE(YLM(NG2,LMX)) 

      PRINT*,"DIMENSION OF THE RAW POTENTIAL MATRIX",NG
      PRINT*,"DIMENSION OF EIGENPROBLEM",NG2

      IF(TTIMING)CALL TIMING$CLOCKOFF('G+K,G2')

      IF(TTIMING)CALL TIMING$CLOCKON('PROJECTORS')
      IND=0
      DO ISP=1,NSP
        DO LN=1,LNX(ISP)
          IND=IND+1
          !BANDS_GETFOFG IS A MODIFIED VERSION OF SETUP$GETFOFG
          CALL BANDS_GETFOFG(GIDG_PROTO,LN,ISP,NG2,G2(1:NG2) &
     &                      ,CELLVOL,BAREPRO(1:NG2,IND))
        ENDDO
      ENDDO
      
      DO I=1,NG2
        CALL GETYLM(LMX,GVECPK(:,I),YLM_)
        YLM(I,:)=YLM_(:) 
      ENDDO        
      
      PRO(:,:,:)=0.0D0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LNX_=LNX(ISP)
        LMNX_=LMNX(ISP)
        LOX_=LOX(:,ISP)
        IBPRO=1+SUM(LNX(1:ISP-1))
        DO I=1,NG2
          SVAR=SUM(GVECPK(:,I)*R(:,IAT))
          EIGR(I)=CMPLX(COS(SVAR),-SIN(SVAR),KIND=8)
        ENDDO

!PRINT*,"LOX",IAT,ISP,LNX_,LMNX_,LOX_(1:LNX_),IBPRO,R(:,IAT)
        CALL WAVES_EXPANDPRO(LNX_,LOX_,LMNX_,NG2,GVECPK(1:3,1:NG2) &
     &                      ,BAREPRO(1:NG2,IBPRO:IBPRO+LNX_-1),LMX &
     &                      ,YLM,EIGR(1:NG2),PRO(IAT,1:NG2,1:LMNX_))
      ENDDO
      IF(TTIMING)CALL TIMING$CLOCKOFF('PROJECTORS')

!      !NOTE: DTKIN/DH FROM SETUP_MAKEPARTIALWAVES NOW 
!      !(AFTER COMMIT 85F857388191C58F050ACF5BE749412FCFE54D12) PRODUCE A HERMITIAN HAMILTONIAN
!      DO IAT=1,NAT
!        DH(:,:,ISPIN,IAT)=0.5D0*(DH(:,:,ISPIN,IAT)+&
!  &                    TRANSPOSE(CONJG(DH(:,:,ISPIN,IAT))))  
!      ENDDO

      !ADD (G+K)^2/2 AND AUGMENTATION TO TI_HK
      DO I=1,NG2
        DO IDIM1=1,NDIM
          I1=I*(NDIM)+(IDIM1-1)
          TI_HK(I1,I1)=TI_HK(I1,I1)+0.5D0*G2(I)
        ENDDO
      ENDDO
      IF(TTIMING)CALL TIMING$CLOCKON('AUGMENTATION')
      PRO(:,:,:)=SQRT(GWEIGHT)*PRO(:,:,:)
      !FIXME: OPTIMIZE THIS BLOCK!!!
!PRINT*,"GWEIGHT",GWEIGHT
      DO I=1,NG2
        DO J=1,I
          CSVAR_H=0.0D0
          CSVAR_S=0.0D0
          DO M=1,NAT
            ISP=ISPECIES(M)
            LMNX_=LMNX(ISP)
            DO K=1,LMNX_
              DO L=1,LMNX_
              !FIXME: THIS WILL NOT WORK FOR NON-COLLINEAR CALC.
!                    PROTMP(1:I)=CONJG(PRO(M,I,K))*PRO(M,1:I,L)*GWEIGHT
!                    TI_HK(I,1:I)=TI_HK(I,1:I)+&
!       &               DH(K,L,ISPIN,M)*PROTMP
!                    TI_SK(I,1:I)=TI_SK(I,1:I)+&
!       &               DO(K,L,ISPIN,M)*PROTMP
                CSVAR_H=CSVAR_H+&
 &                 DH(K,L,ISPIN,M)*CONJG(PRO(M,I,K))*PRO(M,J,L)
                CSVAR_S=CSVAR_S+&
 &                 DO(K,L,ISPIN,M)*CONJG(PRO(M,I,K))*PRO(M,J,L)
              ENDDO
            ENDDO
          ENDDO
          TI_HK(I,J)=TI_HK(I,J)+CSVAR_H
          TI_SK(I,J)=TI_SK(I,J)+CSVAR_S
        ENDDO
      ENDDO
     
!      !CHECK IF TI_HK IS HERMITIAN
!      DO I=1,NG2
!        DO J=I,NG2
!          IF(ABS(TI_HK(I,J)-CONJG(TI_HK(J,I))).GT.1D-10)THEN
!            PRINT*,"H",I,J
!            STOP
!          ENDIF
!          IF(ABS(TI_SK(I,J)-CONJG(TI_SK(J,I))).GT.1D-10)THEN
!            PRINT*,"S",I,J
!            STOP
!          ENDIF
!        ENDDO
!      ENDDO

      !COMPLETE TI_HK AND TI_SK
      DO I=1,NG2
        DO J=I+1,NG2
          TI_HK(I,J)=CONJG(TI_HK(J,I))
          TI_SK(I,J)=CONJG(TI_SK(J,I))
        ENDDO
      ENDDO
      
      IF(TTIMING)CALL TIMING$CLOCKOFF('AUGMENTATION')
      IF(TTIMING)CALL TRACE$POP

      ALLOCATE(U(NG2*NDIM,NG2*NDIM))
      IF(TTEST_POSITIVE_DEFINITE)THEN
        CALL LIB$DIAGC8(NG2*NDIM,TI_SK,E,U)
        PRINT*,"EIGENVALUES OF OVERLAP MATRIX"
        DO I=1,NG2
          PRINT*,I,E(I)
          IF(E(I).LT.0.0D0)THEN
            CALL ERROR$MSG('OVERLAP MATRIX NOT POSITIVE DEFINITE')
            CALL ERROR$STOP('BANDS_KPOINT')        
          ENDIF
        ENDDO
      ENDIF
      !SOLVE GENERALIZED EIGENVALUE PROBLEM
      IF(TTIMING)CALL TIMING$CLOCKON('DIAG')
      IF(METHOD.EQ.1)THEN
        CALL LAPACKOPTIONS$SETCH('GENERALEIGENVALUEC8_MODE','ZHEGVD')
!
      ELSE IF(METHOD.EQ.2)THEN
        CALL LAPACKOPTIONS$SETCH('GENERALEIGENVALUEC8_MODE','ZHEGV')

!     ==========================================================================
!     == USE FEAST HTTP://WWW.ECS.UMASS.EDU/~POLIZZI/FEAST/
!     ==========================================================================
      ELSE IF(METHOD.EQ.3)THEN
#IF DEFINED(CPPVAR_FEAST)
        IF(TTIMING)CALL TRACE$PUSH('FEAST_ZHEGV')
        CALL FEAST_ZHEGV(NG*NDIM,NB,-15.0D0,5.0D0,.FALSE.,TI_HK,TI_SK,E,U)
        IF(TTIMING)CALL TRACE$POP
#ELSE
        CALL ERROR$MSG('LIBRARY FEAST FOR MATRIX DIAGONALIZATION MISSING')
        CALL ERROR$MSG('SEE HTTP://WWW.ECS.UMASS.EDU/~POLIZZI/FEAST/')
        CALL ERROR$I4VAL('METHOD',METHOD)
        CALL ERROR$STOP('BANDS_KPOINT')
#ENDIF
!
!     ==========================================================================
!     == USE SLEPC HTTPS://SLEPC.UPV.ES/                                     ==
!     ==========================================================================
      ELSE IF(METHOD.EQ.4)THEN
#IF DEFINED(CPPVAR_SLEPC)
        IF(TTIMING)CALL TRACE$PUSH('SLEPC_ZHEGV')
        CALL SLEPC_ZHEGV(NG,NB,TI_HK,TI_SK,E)           
        IF(TTIMING)CALL TRACE$POP
#ELSE
        CALL ERROR$MSG('LIBRARY SLEPC FOR MATRIX DIAGONALIZATION MISSING')
        CALL ERROR$MSG('SEE HTTPS://SLEPC.UPV.ES/')
        CALL ERROR$I4VAL('METHOD',METHOD)
        CALL ERROR$STOP('BANDS_KPOINT')
#ENDIF
!
!     ==========================================================================
!     == USE JADAMILU (HTTP://HOMEPAGES.ULB.AC.BE/~JADAMILU/)
!     ==========================================================================
      ELSE IF(METHOD.EQ.5)THEN
#IF DEFINED(CPPVAR_JADAMILU)
        IF(TTIMING)CALL TRACE$PUSH('JADAMILU_ZHEGV')
        CALL JADAMILU_ZHEGV(NG,NB,TI_HK,TI_SK,E)            
        IF(TTIMING)CALL TRACE$POP
#ELSE
        CALL ERROR$MSG('LIBRARY JADAMILU FOR MATRIX DIAGONALIZATION MISSING')
        CALL ERROR$MSG('SEE HTTP://HOMEPAGES.ULB.AC.BE/~JADAMILU/')
        CALL ERROR$I4VAL('METHOD',METHOD)
        CALL ERROR$STOP('BANDS_KPOINT')
#ENDIF
!
      ELSE 
        CALL ERROR$MSG('METHOD_DIAG NOT IMPLEMENTED')
        CALL ERROR$I4VAL('METHOD_DIAG',METHOD)
        CALL ERROR$STOP('BANDS_KPOINT')
      ENDIF
      
      IF(TTIMING)CALL TRACE$PUSH('LIB$GENERALEIGENVALUEC8_ZHEGVD')
!      CALL LIB$GENERALEIGENVALUEC8(NG2*NDIM,TI_HK,TI_SK,E,U,'ZHEGVD')
      CALL LIB$GENERALEIGENVALUEC8(NG2*NDIM,TI_HK,TI_SK,E,U)
      IF(TTIMING)CALL TRACE$POP

      IF(TTIMING)CALL TIMING$CLOCKOFF('DIAG')
      IF(TPROJ)THEN
        IF(TTIMING)CALL TIMING$CLOCKON('PROJECTIONS')
        PRO(:,:,:)=0.0D0
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          LNX_=LNX(ISP)
          LMNX_=LMNX(ISP)
          LOX_=LOX(:,ISP)
          IBPRO=1+SUM(LNX(1:ISP-1))
          DO I=1,NG2
            SVAR=SUM(GVECPK(:,I)*R(:,IAT))
            EIGR(I)=CMPLX(COS(SVAR),SIN(SVAR),KIND=8)
          ENDDO

          CALL WAVES_EXPANDPRO(LNX_,LOX_,LMNX_,NG2,GVECPK(1:3,1:NG2) &
       &                      ,BAREPRO(1:NG2,IBPRO:IBPRO+LNX_-1),LMX,YLM &
       &                      ,EIGR(1:NG2),PRO(IAT,1:NG2,1:LMNX_))
        ENDDO
        PRO(:,:,:)=SQRT(GWEIGHT)*PRO(:,:,:)
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          LNX_=LNX(ISP)
          LMNX_=LMNX(ISP)
          DO LN=1,LMNX_
            DO IB=1,NB
              !FIXME: COMPUTE MULTIPLE SCALARPRODUCTS WITH ONE CALL
              CALL LIB$SCALARPRODUCTC8(.FALSE.,NG2*NDIM,1,PRO(IAT,:,LN),1,U(:,IB),CSVAR)
              !CSVAR=SUM(CONJG(PRO(IAT,1:NG2,LN))*U(1:NG2,IB))
              PROJ(IAT,IB,LN)=CSVAR
!PRINT*,"PROJ",KVEC,IAT,ISP,IB,LN,ABS(CSVAR)
            ENDDO
          ENDDO
        ENDDO
        IF(TTIMING)CALL TIMING$CLOCKOFF('PROJECTIONS')
      ENDIF
      DEALLOCATE(TI_HK)
      DEALLOCATE(TI_SK)
      DEALLOCATE(U)
      DEALLOCATE(EIGR)      
      DEALLOCATE(FOFG)      
      DEALLOCATE(PRO)
      DEALLOCATE(BAREPRO)
      DEALLOCATE(YLM) 
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_GETFOFG(GIDG,IND,ISP,NG_,G2,CELLVOL,F)
!     ******************************************************************
!     **                                                              **
!     **  THIS IS A MODIFIED VERSION OF SETUP$GETFOFGD                **
!     **                                                              **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: GIDG
      INTEGER(4)  ,INTENT(IN)  :: IND
      INTEGER(4)  ,INTENT(IN)  :: ISP
      INTEGER(4)  ,INTENT(IN)  :: NG_     ! #(PLANE WAVES)
      REAL(8)     ,INTENT(IN)  :: G2(NG_) ! G**2
      REAL(8)     ,INTENT(OUT) :: F(NG_)  
      REAL(8)     ,ALLOCATABLE :: FOFG(:)  !(NG)
      REAL(8)                  :: PI
      INTEGER(4)               :: IG
      INTEGER(4)               :: NG
      INTEGER(4)               :: NGAMMA
      REAL(8)                  :: CELLVOL
      REAL(8)                  :: G
!     ******************************************************************
      IF(NG_.EQ.0) RETURN
      CALL RADIAL$GETI4(GIDG,'NR',NG)
      ALLOCATE(FOFG(NG))
      FOFG(:)=PROOFG(:,IND,ISP)
!
!     ==================================================================
!     == INTERPOLATE VALUES FROM RADIAL GRID
!     ==================================================================
      NGAMMA=0
      G=SQRT(G2(1))
      IF(G.LT.1.D-6) NGAMMA=1
      CALL RADIAL$VALUE(GIDG,NG,FOFG,G,F(1))
      DO IG=2,NG_
        IF(ABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
          F(IG) =F(IG-1)
        ELSE
          G=SQRT(G2(IG))
          IF(G.LT.1.D-6) NGAMMA=IG
          CALL RADIAL$VALUE(GIDG,NG,FOFG,G,F(IG))
        END IF
      ENDDO
      DEALLOCATE(FOFG)
!
!     ==================================================================
!     == DIVIDE BY CELLVOL                                            ==
!     ==================================================================
      F=F/CELLVOL
      RETURN
      END SUBROUTINE
END MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_BANDSTRUCTURE_DIAG(LL_CNTL,NFIL,NFILO)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE BANDDATA_MODULE
      USE RADIAL_MODULE
      USE WAVES_MODULE, ONLY : GSET_TYPE
      USE KPOINTDIAG_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT):: LL_CNTL
      INTEGER(4)   ,INTENT(IN)  :: NFIL
      INTEGER(4)   ,INTENT(IN)  :: NFILO
      INTEGER(4)                :: NFILN
      LOGICAL(4)                :: TPRINT=.FALSE.
      LOGICAL(4)                :: TCHK,TCHK1,TCHK2
      INTEGER(4)                :: METHOD_DIAG
      REAL(8)                   :: XK1(3)    ! INITIAL K-POINT IN REL. COORD.
      REAL(8)                   :: XK2(3)    ! FINAL K-POINT IN REL. COORD.
      REAL(8)                   :: KVEC1(3)  ! INITIAL K-POINT IN ABS. COORD.
      REAL(8)                   :: KVEC2(3)  ! FINAL K-POINT IN ABS. COORD.
      REAL(8)                   :: GUNIT     ! 2*PI/LUNIT
      REAL(8)                   :: GBASINV(3,3)
      INTEGER(4)                :: NK,NKDIAG
      INTEGER(4)                :: NB
      INTEGER(4)                :: ILINE,NLINE,NFILE
      INTEGER(4)                :: ISPIN
      REAL(8)                   :: X1,X2
      CHARACTER(512)            :: BANDDATAFILE,ID
      CHARACTER(512)            :: FILE
      LOGICAL(4)                :: TPROJ

      INTEGER(4)                :: I,J,I1,J1,IDIM1,IDIM2
      REAL(8)                   :: XK(3)
      REAL(8)                   :: KVEC(3)
      REAL(8)                   :: T1,T2,T3
      REAL(8)                   :: GWEIGHT
      INTEGER(4)                :: NG,IG
      INTEGER(4)                :: IKDIAG
      REAL(8)                   :: CELLVOL
      TYPE(GSET_TYPE)           :: GSET
      REAL(8)   ,ALLOCATABLE    :: GVEC(:,:)   !G
      REAL(8)   ,ALLOCATABLE    :: G2(:)
      INTEGER(4)                :: GIDG_PROTO
      COMPLEX(8),ALLOCATABLE    :: TI_H(:,:,:)
      COMPLEX(8),ALLOCATABLE    :: TI_S(:,:,:)
      COMPLEX(8),ALLOCATABLE    :: PROJ(:,:,:)
      INTEGER(4)                :: IAT,IBPRO,LN,IND,IPRO,ISP,K,L,M
      INTEGER(4)                :: LMNX_,LNX_,LMN1,LMN2
      REAL(8),ALLOCATABLE       :: E(:)
      REAL(8),ALLOCATABLE       :: FATBANDVAL(:,:,:)
      REAL(8),ALLOCATABLE       :: EIGVAL(:,:)
      REAL(8),ALLOCATABLE       :: XKVAL(:,:)
      REAL(8),ALLOCATABLE       :: KVECVAL(:,:)
      INTEGER(4)                :: IB
      COMPLEX(8)                :: CSVAR
      CHARACTER(64)             :: ATOM
      INTEGER(4)                :: NFILBAND
      INTEGER(4)                :: NFILFATBAND
      INTEGER(4)                :: NFATBAND,IFATBAND,ILINESPERBAND
      INTEGER(4),ALLOCATABLE    :: FATBANDIPRO(:),FATBANDIAT(:),LINESPERBAND(:)
      REAL(8)   ,ALLOCATABLE    :: FATBANDMAXWIDTH(:)
      CHARACTER(512),ALLOCATABLE:: FATBANDFILE(:)
      REAL(8)                   :: FATBANDMAX,SVAR,SVAR1,SVAR2
      LOGICAL(4)                :: TAPPEND
      REAL(8)                   :: COORD,COORD0
      INTEGER(4)                :: THISTASK,NTASKS,INDEX
      REAL(8)   ,PARAMETER      :: PI=4.D0*ATAN(1.D0)
      REAL(8)                   :: ANGSTROM
!     **************************************************************************
                            CALL TRACE$PUSH('BANDS_BANDSTRUCTURE_DIAG')
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)

      CALL MPE$QUERY('~',NTASKS,THISTASK)
      CALL TIMING$START
!
!     ==========================================================================
!     ==  SET BANDDATAFILE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'INPUTFILE',NFILE)
      IF(NFILE.EQ.1)THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'INPUTFILE',NFILE)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(TCHK)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,BANDDATAFILE)
          ID=+'BANDDATAIN'
          CALL FILEHANDLER$SETFILE(ID,.FALSE.,BANDDATAFILE)
          CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
        ELSE
          CALL ERROR$MSG('!BCNTL!INPUTFILE!NAME NOT FOUND')
          CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
        ENDIF
      ELSE IF (NFILE.GT.1)THEN
        CALL ERROR$MSG('MULTIPLE INPUT FILES GIVEN (!BCNTL!INPUTFILE)')
        CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
      ELSE
        ID=+'BANDDATAIN'
        CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.BANDDATA')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
      ENDIF
!
!     ==========================================================================
!     ==  READ BANDDATAFILE                                                   ==
!     ==========================================================================
      CALL BANDDATA$READFILE
!
!     ==========================================================================
!     ==  CHECKS                                                              ==
!     ==========================================================================
      IF(NDIM.EQ.2)THEN
        CALL ERROR$MSG('ONLY NDIM=1 AND NSPIN=1 OR 2 IMPLEMENTED AND TESTED')
        CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
      ENDIF
!
!     ==========================================================================
!     ==  READ BCNTL                                                          ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'METHOD_DIAG',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'METHOD_DIAG',1,METHOD_DIAG)
      ELSE
        METHOD_DIAG=1
      ENDIF
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EPWPSI',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'EPWPSI',1,EPW)
        IF(THISTASK.EQ.1)WRITE(NFILO,*)'WARNING: EPWPSI SET TO',EPW,' RY!'
        !INPUT IS IN RYDBERG, WE NEED HARTREE
        EPW=0.5D0*EPW
      ENDIF
!
!     =========================================================================
!     ==  CONSTRUCT K-INDEPENDENT PART OF HAMILTONIAN                        ==
!     =========================================================================
      CALL BANDS_KINDEP(NG,GVEC,G2,GWEIGHT,GIDG_PROTO,TI_H,TI_S)
!
!     =========================================================================
!     ==  CONSTRUCT BAND STRUCTURE                                           ==
!     =========================================================================
!     == DEFAULT OUTPUT FILE ==================================================
      CALL FILEHANDLER$SETFILE('BANDS',.TRUE.,-'.BANDS')
      CALL FILEHANDLER$SETSPECIFICATION('BANDS','STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION('BANDS','POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION('BANDS','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('BANDS','FORM','FORMATTED')
!
!     ==========================================================================
!     ==  LOOP OVER ALL LINES IN RECIPROCAL SPACE                             ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'LINE',NLINE)
      DO ILINE=1,NLINE
        IF(THISTASK.EQ.1)WRITE(NFILO,FMT='(72("="))')
        IF(THISTASK.EQ.1)WRITE(NFILO,*)'LINE-BLOCK:',ILINE,' OF ',NLINE
        IF(THISTASK.EQ.1)WRITE(NFILO,FMT='(72("="))')
!        
        CALL LINKEDLIST$SELECT(LL_CNTL,'~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
        CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
        CALL LINKEDLIST$SELECT(LL_CNTL,'LINE',ILINE)
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TAPPEND',0,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'TAPPEND',1,TAPPEND)
        ELSE
          TAPPEND=.FALSE.
        ENDIF

        IF(THISTASK.EQ.1)THEN
          CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',0,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FILE)
            CALL FILEHANDLER$SETFILE('BANDS',.FALSE.,FILE)
            IF(TAPPEND)THEN
              WRITE(NFILO,*)'ATTACHING OUTPUT FOR LINE-BLOCK',ILINE &
     &                    ,' TO EXISTING FILE ',TRIM(FILE)
              CALL FILEHANDLER$SETSPECIFICATION('BANDS','STATUS','OLD')
            ELSE
              WRITE(NFILO,*)'WRITING OUTPUT FOR LINE-BLOCK',ILINE &
     &                    ,' TO NEW FILE ',TRIM(FILE)
              CALL FILEHANDLER$SETSPECIFICATION('BANDS','STATUS','REPLACE')
            ENDIF
            CALL FILEHANDLER$SETSPECIFICATION('BANDS','POSITION','APPEND')
            CALL FILEHANDLER$SETSPECIFICATION('BANDS','ACTION','WRITE')
            CALL FILEHANDLER$SETSPECIFICATION('BANDS','FORM','FORMATTED')
            WRITE(NFILO,*)'OUTPUT FILE: ',TRIM(FILE)
          ELSE
            CALL ERROR$MSG('NO OUTPUT FILE GIVEN')
            CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
            CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
          END IF
        ENDIF
!
!
!       ========================================================================
!       ==  READ K-VECTORS                                                    ==
!       ========================================================================
!       == SCALE FACTOR =======================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'LUNIT',0,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'LUNIT[AA]',0,TCHK2)
        IF(TCHK1.AND.TCHK2) THEN
          CALL ERROR$MSG('!BCNTL!BANDSTRUCTURE!LINE:LUNIT AND LUNIT[AA] ')
          CALL ERROR$MSG('MUST NOT BE SPECIFIED SIMULTANEOUSLY')
          CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
          CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
        ELSE           
          GUNIT=1.D0
          IF(TCHK1) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'LUNIT',1,SVAR)
            GUNIT=2.D0*PI/SVAR
          ELSE IF(TCHK2) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'LUNIT[AA]',1,SVAR)
            SVAR=SVAR*ANGSTROM
            GUNIT=2.D0*PI/SVAR
           ELSE
!            == THIS OPTION IS OBSOLETE AND IS THERE FOR BACKWARD COMPATIBILITY
             CALL LINKEDLIST$EXISTD(LL_CNTL,'KVECSCALE',0,TCHK)
             IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'KVECSCALE',1,SVAR)
               GUNIT=2.D0*PI/SVAR
             END IF
          END IF
        END IF
!
        CALL LIB$INVERTR8(3,GBAS,GBASINV) 

        CALL LINKEDLIST$EXISTD(LL_CNTL,'XK1',0,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'KVEC1',0,TCHK2)
        IF(TCHK1.AND.TCHK2) THEN
          CALL ERROR$MSG('!BCNTL!BANDSTRUCTURE!LINE:XK1 AND KVEC1')
          CALL ERROR$MSG('MUST NOT BE SPECIFIED SIMULTANEOUSLY')
          CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
          CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
        ELSE
          IF(TCHK1)THEN
            CALL LINKEDLIST$GET(LL_CNTL,'XK1',1,XK1)
            KVEC1=MATMUL(GBAS,XK1)
          ELSE IF(TCHK2)THEN
            CALL LINKEDLIST$GET(LL_CNTL,'KVEC1',1,KVEC1)
            KVEC1=GUNIT*KVEC1
            XK1=MATMUL(GBASINV,KVEC1)
          ELSE
            CALL ERROR$MSG('NEITHER XK1 NOR KVEC1 GIVEN')
            CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
            CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
          ENDIF
        END IF

        CALL LINKEDLIST$EXISTD(LL_CNTL,'XK2',0,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'KVEC2',0,TCHK2)
        IF(TCHK1)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'XK2',1,XK2)
          KVEC2=MATMUL(GBAS,XK2)
        ELSE IF (TCHK2)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'KVEC2',1,KVEC2)
          KVEC2=KVEC2*GUNIT
          XK2=MATMUL(GBASINV,KVEC2)
        ELSE
          CALL ERROR$MSG('NEITHER XK2 NOR KVEC2 GIVEN')
          CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
          CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
        ENDIF

        IF(THISTASK.EQ.1)THEN
          WRITE(NFILO,*)'FIRST K-POINT IN RELATIVE COORDINATES:',XK1
          WRITE(NFILO,*)'FIRST K-POINT IN ABSOLUTE COORDINATES:',KVEC1
          WRITE(NFILO,*)'LAST K-POINT IN RELATIVE COORDINATES:',XK2
          WRITE(NFILO,*)'LAST K-POINT IN ABSOLUTE COORDINATES:',KVEC2
        ENDIF
        
        !NK IS NUMBER OF POINTS FOR OUTPUT
        !IF NKDIAG<NK: DO 1D-FFT INTERPOALTION
        !IF NKDIAG=NK: USE RESULTS FROM DIAGONALISATION DIRECTLY
        !IF NKDIAG>NK: ERROR
        NB=20
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NB',0,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NB',1,NB)
        NK=10
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NK',0,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NK',1,NK)
        !NKDIAG IS NUMBER OF K-POINTS IN DIAGONALISATION
        NKDIAG=NK
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NKDIAG',0,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NKDIAG',1,NKDIAG)
        IF(NKDIAG.GT.NK)THEN
          CALL ERROR$MSG('NKDIAG>NK')
          CALL ERROR$I4VAL('NK',NK)
          CALL ERROR$I4VAL('NKDIAG',NKDIAG)
          CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
        ENDIF
        !FIXME: IMPLEMENT 1D FFT INTERPOLATION
        IF(NKDIAG.LT.NK)THEN
          CALL ERROR$MSG('1D FFT INTERPOLATION NOT YET IMPLEMENTED!!')
          CALL ERROR$I4VAL('NK',NK)
          CALL ERROR$I4VAL('NKDIAG',NKDIAG)
          CALL ERROR$STOP('BANDS_BANDSTRUCTURE_DIAG')
        ENDIF
        IF(TPRINT)PRINT*,XK1,XK2,NB,NK,NKDIAG,X1,X2,TRIM(FILE)
!       == SPECIFY SPIN DIRECTION =============================================
        ISPIN=1
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',0,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,ISPIN)
PRINT*,'+++++++++++++++++ SPIN=',ISPIN,'++++++++++++++++++'
        IF(ISPIN.GT.NSPIN) THEN
          CALL ERROR$MSG('ISPIN EXCEEDS RANGE')
          CALL ERROR$STOP('MAIN')
        END IF
!       == SPECIFY FATBAND ===================================================
        CALL LINKEDLIST$SELECT(LL_CNTL,'~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
        CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
        CALL LINKEDLIST$SELECT(LL_CNTL,'LINE',ILINE)
        CALL LINKEDLIST$NLISTS(LL_CNTL,'FATBAND',NFATBAND)

        IF(NFATBAND.GE.1)THEN
          ALLOCATE(FATBANDIPRO(NFATBAND))
          ALLOCATE(FATBANDIAT(NFATBAND))
          ALLOCATE(LINESPERBAND(NFATBAND))
          ALLOCATE(FATBANDMAXWIDTH(NFATBAND))
          ALLOCATE(FATBANDFILE(NFATBAND))

          DO IFATBAND=1,NFATBAND
            CALL LINKEDLIST$SELECT(LL_CNTL,'~')
            CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
            CALL LINKEDLIST$SELECT(LL_CNTL,'BANDSTRUCTURE')
            CALL LINKEDLIST$SELECT(LL_CNTL,'LINE',ILINE)
            CALL LINKEDLIST$SELECT(LL_CNTL,'FATBAND',IFATBAND)
            !FILE
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',0,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FATBANDFILE(IFATBAND))
            ELSE
              CALL ERROR$I4VAL('FATBAND',IFATBAND)
              CALL ERROR$MSG('FILE NOT GIVEN')
              CALL ERROR$STOP('MAIN')
            END IF
            !ATOM
            CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',0,TCHK)
            IF(TCHK)THEN
              CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,ATOM)
            ELSE
              CALL ERROR$I4VAL('FATBAND',IFATBAND)
              CALL ERROR$MSG('ATOM NOT GIVEN')
              CALL ERROR$STOP('MAIN')
            ENDIF
            !FIND NUMBER OF ATOM IN ATOMID
            FATBANDIAT(IFATBAND)=0
            DO IAT=1,NAT
              IF(ATOMID(IAT).EQ.ATOM)THEN
                FATBANDIAT(IFATBAND)=IAT
                EXIT
              ENDIF
            ENDDO
            IF(FATBANDIAT(IFATBAND).EQ.0)THEN
              CALL ERROR$MSG('ATOM NOT FOUND')
              CALL ERROR$CHVAL('ATOM',ATOM)
              CALL ERROR$I4VAL('FATBAND',IFATBAND)
              CALL ERROR$STOP('MAIN')
            ENDIF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'IPRO',0,TCHK)
            IF(TCHK)THEN
              CALL LINKEDLIST$GET(LL_CNTL,'IPRO',1,FATBANDIPRO(IFATBAND))
            ELSE
              CALL ERROR$MSG('IPRO NOT GIVEN')
              CALL ERROR$I4VAL('FATBAND',IFATBAND)
              CALL ERROR$STOP('MAIN')
            ENDIF

            CALL LINKEDLIST$EXISTD(LL_CNTL,'MAXWIDTH',0,TCHK)
            IF(TCHK)THEN
              CALL LINKEDLIST$GET(LL_CNTL,'MAXWIDTH',1,FATBANDMAXWIDTH(IFATBAND))
            ELSE
              FATBANDMAXWIDTH(IFATBAND)=1.0D0
            ENDIF
            CALL LINKEDLIST$EXISTD(LL_CNTL,'LINESPERBAND',0,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_CNTL,'LINESPERBAND',1,LINESPERBAND(IFATBAND))
            ELSE
              LINESPERBAND(IFATBAND)=10
            ENDIF
          ENDDO
        ENDIF
                            CALL TRACE$PASS('AFTER READ BNCTL')

        ALLOCATE(EIGVAL(NB,NKDIAG))
        ALLOCATE(KVECVAL(3,NKDIAG))
        ALLOCATE(XKVAL(3,NKDIAG))
        ALLOCATE(FATBANDVAL(NFATBAND,NKDIAG,NB))
        TPROJ=NFATBAND.GE.1
        EIGVAL(:,:)=0.0D0
        KVECVAL(:,:)=0.0D0
        XKVAL(:,:)=0.0D0
        FATBANDVAL(:,:,:)=0.0D0

!       == ITERATE K-POINTS =================================================
        CALL MPE$SYNC('~')
                            CALL TRACE$PUSH('ITERATE KPOINTS')
        INDEX=1
        DO IKDIAG=0,NKDIAG-1
          INDEX=INDEX+1
          IF(INDEX.GT.NTASKS)INDEX=1
          IF(INDEX.NE.THISTASK)CYCLE

          ALLOCATE(PROJ(NAT,NB,LMNXX))
          ALLOCATE(E(NG*NDIM))
          IF(TPRINT)PRINT*,"IKDIAG",IKDIAG
          XK=XK1+(XK2-XK1)*REAL(IKDIAG,KIND=8)/REAL(MAX(NKDIAG-1,1),KIND=8)
!HERE THE K-POINT 
          KVEC=MATMUL(GBAS,XK)
          KVECVAL(:,IKDIAG+1)=KVEC
          XKVAL(:,IKDIAG+1)=XK
          
          IF(TPRINT) THEN
            WRITE(NFILO,*)'LINE ',ILINE,' OF ',NLINE,' K-POINT ',IKDIAG &
       &                 ,' OF ',NKDIAG-1,' IN RELATIVE COORDINATES ',XK
            WRITE(NFILO,*)'LINE ',ILINE,' OF ',NLINE,' K-POINT ',IKDIAG &
       &                 ,' OF ',NKDIAG-1,' IN ABSOLUTE COORDINATES ',KVEC
          END IF
!
!         == DIAGONALIZE HAMILTONIAN FOR A GIVEN K-POINT AND SPIN ==============
          CALL BANDS_KPOINT(NG,NB,ISPIN,METHOD_DIAG,GIDG_PROTO,TPROJ,KVEC &
      &                    ,GVEC,TI_H,TI_S,E,PROJ)
          EIGVAL(1:NB,IKDIAG+1)=E(1:NB)*27.21139D0

          IF(TPRINT) THEN
            PRINT*,'DIAGBANDS_EIG',SQRT(SUM(KVEC**2)),E(1:NB)*27.21139D0
          END IF

!
!         ====================================================================
!         ==  FAT BANDS                                                     ==
!         ====================================================================
          IF(NFATBAND.GE.1)THEN
            DO IFATBAND=1,NFATBAND
              IF(TPRINT)PRINT*,"FATBAND IFATBAND",IFATBAND
              IF(TPRINT)PRINT*,"FATBAND IAT",FATBANDIAT(IFATBAND)
              IF(TPRINT)PRINT*,"FATBAND IPRO",FATBANDIPRO(IFATBAND)
              DO IB=1,NB
                !FIXME: INTEGRATE PDOS INTERFACE HERE
                CSVAR=PROJ(FATBANDIAT(IFATBAND),IB,FATBANDIPRO(IFATBAND))
                FATBANDVAL(IFATBAND,IKDIAG+1,IB)=ABS(CSVAR)**2
              ENDDO
            ENDDO
          ELSE
            FATBANDVAL(:,:,:)=0.0D0 
          ENDIF
          DEALLOCATE(PROJ)
          DEALLOCATE(E)
        ENDDO
        CALL MPE$SYNC('~')
        CALL MPE$COMBINE('~','+',EIGVAL)
        CALL MPE$COMBINE('~','+',KVECVAL)
        CALL MPE$COMBINE('~','+',XKVAL)
        CALL MPE$COMBINE('~','+',FATBANDVAL)
                            CALL TRACE$POP()
!
!       ========================================================================
!       ==  WRITE BANDS                                                       ==
!       ========================================================================
        IF(THISTASK.EQ.1)THEN
                              CALL TRACE$PUSH('WRITE_BANDS')
          CALL FILEHANDLER$UNIT('BANDS',NFILBAND)
          IF(NB.GT.100) THEN
            CALL ERROR$MSG('NUMBER OF BANDS EXCEEDS LIMIT OF 100')
            CALL ERROR$I4VAL('NB ',NB)
            CALL ERROR$STOP('BANDS_PLOTBANDS')
          END IF

          !ITERATE K-POINTS
          WRITE(NFILBAND,*)'#DATA FOR LINE BLOCK',ILINE
          WRITE(NFILBAND,*)'#NUMBER OF K-POINTS',NKDIAG
          WRITE(NFILBAND,*)'#FIRST K-POINT IN ABSOLUTE COORDINATES ',KVECVAL(:,1)
          WRITE(NFILBAND,*)'#FIRST K-POINT IN RELATIVE COORDINATES ',XKVAL(:,1)
          WRITE(NFILBAND,*)'#LAST  K-POINT IN ABSOLUTE COORDINATES ',KVECVAL(:,NKDIAG)
          WRITE(NFILBAND,*)'#LAST  K-POINT IN RELATIVE COORDINATES ',XKVAL(:,NKDIAG)
          DO IKDIAG=1,NKDIAG
            !WRITE EIGENVALUES
            IF(IKDIAG.EQ.1.AND..NOT.TAPPEND)COORD0=0.0D0
            COORD=COORD0+SQRT(SUM((KVECVAL(:,IKDIAG)-KVECVAL(:,1))**2))/GUNIT
            WRITE(NFILBAND,FMT='(103F10.5)')COORD,EIGVAL(:,IKDIAG)
          ENDDO
          FLUSH(NFILBAND)
          CALL FILEHANDLER$CLOSE('BANDS')
          COORD0=COORD
          !WRITE FATBANDS
          DO IFATBAND=1,NFATBAND
            CALL FILEHANDLER$SETFILE('FATBANDS',.FALSE.,FATBANDFILE(IFATBAND))
            CALL FILEHANDLER$SETSPECIFICATION('FATBANDS','STATUS','REPLACE')
            CALL FILEHANDLER$SETSPECIFICATION('FATBANDS','POSITION','APPEND')
            CALL FILEHANDLER$SETSPECIFICATION('FATBANDS','ACTION','WRITE')
            CALL FILEHANDLER$SETSPECIFICATION('FATBANDS','FORM','FORMATTED')
            CALL FILEHANDLER$UNIT('FATBANDS',NFILFATBAND)
            !ITERATE K-POINTS
            DO IKDIAG=1,NKDIAG
              WRITE(NFILFATBAND,FMT='(103F10.5)',ADVANCE="NO") &
       &                 SQRT(SUM((KVECVAL(:,IKDIAG)-KVECVAL(:,1))**2))/GUNIT

              FATBANDMAX=MAXVAL(FATBANDVAL(IFATBAND,:,:))
              IF(TPRINT)PRINT*,"FATBAND FATBANDMAX",FATBANDMAX
              DO ILINESPERBAND=0,LINESPERBAND(IFATBAND)-1
                SVAR1=2.0D0/REAL(LINESPERBAND(IFATBAND)-1,KIND=8) &
       &                   *REAL(ILINESPERBAND,KIND=8)-1.0D0
                SVAR2=0.5D0*SVAR1*FATBANDMAXWIDTH(IFATBAND)/FATBANDMAX
                WRITE(NFILFATBAND,FMT='(100F10.5)',ADVANCE="NO") &
       &                    EIGVAL(:,IKDIAG)+SVAR2*FATBANDVAL(IFATBAND,IKDIAG,:)
              ENDDO
              WRITE(NFILFATBAND,FMT='(A)')" "
            ENDDO
          ENDDO
                              CALL TRACE$POP
        ENDIF
        DEALLOCATE(EIGVAL)
        DEALLOCATE(FATBANDVAL)
        DEALLOCATE(KVECVAL)
        DEALLOCATE(XKVAL)
        IF(NFATBAND.GE.1)THEN
          DEALLOCATE(FATBANDIPRO)
          DEALLOCATE(FATBANDIAT)
          DEALLOCATE(LINESPERBAND)
          DEALLOCATE(FATBANDMAXWIDTH)
          DEALLOCATE(FATBANDFILE)
        ENDIF
      ENDDO
                            CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDS_BANDSTRUCTURE_DIAG

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_PDOS(LL_CNTL,NFIL,NFILO)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE BANDDATA_MODULE
      USE RADIAL_MODULE
      USE WAVES_MODULE, ONLY : GSET_TYPE
      USE KPOINTDIAG_MODULE
      USE BRILLOUIN_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT)  :: LL_CNTL
      INTEGER(4)   ,INTENT(IN)     :: NFIL
      INTEGER(4)   ,INTENT(IN)     :: NFILO
      INTEGER(4)                   :: ILINE,NLINE,NFILE
      CHARACTER(512)               :: BANDDATAFILE,ID
      LOGICAL(4)                   :: TCHK
      INTEGER(4)                   :: IB
      INTEGER(4)                   :: NKP
      INTEGER(4)                   :: IKP
      LOGICAL(4)                   :: TINV=.TRUE.
      INTEGER(4)                   :: NKDIV(3)
      INTEGER(4)                   :: ISHIFT(3)
      INTEGER(4)                   :: NB
      INTEGER(4)                   :: I,J,I1,I2,IE,ISPIN

      REAL(8)      ,ALLOCATABLE    :: BK(:,:)
      REAL(8)      ,ALLOCATABLE    :: XK(:,:)
      REAL(8)                      :: RBASINV(3,3)
      REAL(8),ALLOCATABLE          :: EB(:,:)
      REAL(8),ALLOCATABLE          :: EBTMP(:,:)
      REAL(8),ALLOCATABLE          :: WGHT(:,:)
      TYPE(EWGHT_TYPE),ALLOCATABLE :: EWGHT(:,:)
      REAL(8),ALLOCATABLE          :: A(:,:)
      REAL(8)                      :: RNTOT
      REAL(8)                      :: EF
      REAL(8)                      :: SUMA
      REAL(8)                      :: EMIN,EMAX,EMIN2,EMAX2
      INTEGER(4)                   :: NE
      REAL(8),ALLOCATABLE          :: DOS(:,:)
     
      INTEGER(4)                   :: NG
      INTEGER(4)                   :: METHOD_DIAG 
      INTEGER(4)                   :: GIDG_PROTO
      LOGICAL(4)                   :: TPROJ=.TRUE.
      REAL(8)                      :: KVEC(3)
      REAL(8),ALLOCATABLE          :: GVEC(:,:)
      REAL(8),ALLOCATABLE          :: G2(:)
      COMPLEX(8),ALLOCATABLE       :: TI_H(:,:,:)
      COMPLEX(8),ALLOCATABLE       :: TI_S(:,:,:)
      REAL(8),ALLOCATABLE          :: E(:)
      COMPLEX(8),ALLOCATABLE       :: PROJ(:,:,:)
      COMPLEX(8),ALLOCATABLE       :: PROJK(:,:,:,:)
      REAL(8)                      :: GWEIGHT

      INTEGER(4)                   :: SPACEGROUP
      REAL(8)                      :: A0,B0,C0,ALPHA,BETA,GAMMA,EI
      INTEGER(4)                   :: NSYM,NOP
      INTEGER(4),PARAMETER         :: NOPX=48
      INTEGER(4)                   :: IARB(3)
      CHARACTER(3)                 :: BRAVAIS
      INTEGER(4)                   :: IIO(3,3,NOPX)
      REAL(8)                      :: C(3,NOPX)
      INTEGER(4)                   :: ISYM
      LOGICAL(4)                   :: TSHIFT

      INTEGER(4)                   :: NFILOUT,NFILIN
      REAL(8),ALLOCATABLE          :: OCC(:,:)
      REAL(8),ALLOCATABLE          :: WKPT(:)
      COMPLEX(8),ALLOCATABLE       :: PROJTMP(:,:,:)
      INTEGER(4)                   :: LMN,INDEX,IAT

      INTEGER(4)                   :: NPDOS
      CHARACTER(512)               :: FILE
     
      LOGICAL(4)                   :: TPRINT=.FALSE.

      INTEGER(4)                   :: THISTASK,NTASKS
      INTEGER(4),ALLOCATABLE       :: JOBLIST(:,:)

      INTEGER(4)                   :: LMNX_,LNX_,ISP

      REAL(8)                      :: T1,T2
!     **************************************************************************
      IF(NDIM.EQ.2)THEN
        CALL ERROR$MSG('ONLY NDIM=1 AND NSPIN=1 OR 2 IMPLEMENTED AND TESTED')
        CALL ERROR$STOP('BANDS_PDOS')
      ENDIF
                            CALL TRACE$PUSH('BANDS_PDOS')
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      CALL TIMING$START
!
!     ==========================================================================
!     ==  SET BANDDATAFILE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'INPUTFILE',NFILE)
      IF(NFILE.EQ.1)THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'INPUTFILE',NFILE)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(TCHK)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,BANDDATAFILE)
          ID=+'BANDDATAIN'
          CALL FILEHANDLER$SETFILE(ID,.FALSE.,BANDDATAFILE)
          CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
        ELSE
          CALL ERROR$MSG('!BCNTL!INPUTFILE!NAME NOT FOUND')
          CALL ERROR$STOP('BANDS_PDOS')
        ENDIF
      ELSE IF (NFILE.GT.1)THEN
        CALL ERROR$MSG('MULTIPLE INPUT FILES GIVEN (!BCNTL!INPUTFILE)')
        CALL ERROR$STOP('BANDS_PDOS')
      ELSE
        ID=+'BANDDATAIN'
        CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.BANDDATA')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
      ENDIF
!
!     ==========================================================================
!     ==  READ BANDDATAFILE                                                   ==
!     ==========================================================================
      CALL BANDDATA$READFILE
!
!     ==========================================================================
!     ==  CHECKS                                                              ==
!     ==========================================================================
      IF(NDIM.EQ.2)THEN
        CALL ERROR$MSG('ONLY NDIM=1 AND NSPIN=1 OR 2 IMPLEMENTED AND TESTED')
        CALL ERROR$STOP('BANDS_PDOS')
      ENDIF

!
!     ==========================================================================
!     ==  READ BCNTL                                                          ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'PDOS')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'METHOD_DIAG',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'METHOD_DIAG',1,METHOD_DIAG)
      ELSE
        METHOD_DIAG=1
      ENDIF
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EPWPSI',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'EPWPSI',1,EPW)
        IF(THISTASK.EQ.1)WRITE(NFILO,*)'WARNING: EPWPSI SET TO',EPW,' RY!'
        !INPUT IS IN HARTREE
        EPW=0.5D0*EPW
      ENDIF
!
!     ==========================================================================
!     ==  READ BCNTL: PDOS BLOCK                                              ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'PDOS',NPDOS)
      IF(NPDOS.GT.1)THEN
        CALL ERROR$MSG('TOO MANY PDOS-BLOCKS GIVEN, JUST ONE BLOCK ALLOWED')
        CALL ERROR$I4VAL('NPDOS',NPDOS)
        CALL ERROR$STOP('BANDS_PDOS')
      ENDIF
      IF(NPDOS.EQ.0)THEN
        CALL ERROR$MSG('NO PDOS-BLOCK GIVEN')
        CALL ERROR$I4VAL('NPDOS',NPDOS)
        CALL ERROR$STOP('BANDS_PDOS')
      ENDIF
      CALL LINKEDLIST$SELECT(LL_CNTL,'PDOS',NPDOS)

      !PDOSINFILE
      CALL LINKEDLIST$EXISTD(LL_CNTL,'PDOSINFILE',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'PDOSINFILE',1,FILE)
        CALL FILEHANDLER$SETFILE('PDOS',.FALSE.,FILE)
      ELSE
        CALL FILEHANDLER$SETFILE('PDOS',.TRUE.,-'.PDOS')
      ENDIF
      CALL FILEHANDLER$SETSPECIFICATION('PDOS','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('PDOS','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PDOS','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('PDOS','FORM','UNFORMATTED')

      !PDOSOUTFILE
      IF(THISTASK.EQ.1)THEN
        CALL LINKEDLIST$EXISTD(LL_CNTL,'PDOSOUTFILE',0,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'PDOSOUTFILE',1,FILE)
          CALL FILEHANDLER$SETFILE('PDOSOUT',.FALSE.,FILE)
        ELSE
          CALL FILEHANDLER$SETFILE('PDOSOUT',.TRUE.,-'.PDOSOUT')
        ENDIF
        CALL FILEHANDLER$SETSPECIFICATION('PDOSOUT','STATUS','UNKNOWN')
        CALL FILEHANDLER$SETSPECIFICATION('PDOSOUT','POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION('PDOSOUT','ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION('PDOSOUT','FORM','UNFORMATTED')
        CALL FILEHANDLER$FILENAME('PDOSOUT',FILE)
        WRITE(NFILO,*)'WRITING NEW PDOS-FILE TO FILE:',TRIM(FILE)
      ENDIF
      
      !NB
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NB',0,TCHK)
      IF(TCHK)THEN
        CALL LINKEDLIST$GET(LL_CNTL,'NB',1,NB)
      ELSE
        NB=20
      ENDIF
      
      !NKDIV
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NKDIV',1,TCHK)
      IF(TCHK)THEN
        CALL LINKEDLIST$GET(LL_CNTL,'NKDIV',1,NKDIV)
      ELSE
        CALL ERROR$MSG('NKDIV NOT GIVEN')
        CALL ERROR$STOP('BANDS_PDOS')
      ENDIF
      
      !ISHIFT
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ISHIFT',1,TCHK)
      IF(TCHK)THEN
        CALL LINKEDLIST$GET(LL_CNTL,'ISHIFT',1,ISHIFT)
      ELSE
        ISHIFT(1)=0
        ISHIFT(2)=0
        ISHIFT(3)=0
      ENDIF

      !SPACEGROUP
      CALL LINKEDLIST$EXISTD(LL_CNTL,'SPACEGROUP',0,TCHK)
      IF(TCHK)THEN
        CALL LINKEDLIST$GET(LL_CNTL,'SPACEGROUP',1,SPACEGROUP)
      ELSE
        SPACEGROUP=1 !NO SYMMETRY
      ENDIF
      
      !TSHIFT
      CALL LINKEDLIST$EXISTD(LL_CNTL,'TSHIFT',0,TCHK)
      IF(TCHK)THEN
        CALL LINKEDLIST$GET(LL_CNTL,'TSHIFT',1,TSHIFT)
      ELSE
        TSHIFT=.FALSE.
      ENDIF

      IF(NSPIN.EQ.1)THEN
        RNTOT=0.5D0*NEL
      ELSE
        RNTOT=NEL
      ENDIF

      !CHECK IF NB*NSPIN>RNTOT
      IF(NB*NSPIN.LT.RNTOT)THEN
        CALL ERROR$MSG('NB*NPSIN IS LESS THAN THE NUMBER OF BANDS, PLEASE INCREASE NB')
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$R8VAL('NUMBER OF BANDS',RNTOT)
        CALL ERROR$STOP('BANDS_PDOS')
      ENDIF

                            CALL TRACE$PASS('AFTER READ BNCTL')
      IF(.TRUE.)THEN
!     ==========================================================================
!     ==  FIND IRREDUCIBLE K-POINTS AND TETRAHEDRA                            ==
!     ==========================================================================
        NKP=NKDIV(1)*NKDIV(2)*NKDIV(3)
        CALL SPACEGROUP$SETI4('SPACEGROUP',SPACEGROUP)
        CALL SPACEGROUP$GETCH('BRAVAIS',BRAVAIS)
        CALL BRILLOUIN$CHECKRBAS(BRAVAIS,A0,B0,C0,ALPHA,BETA,GAMMA,RBAS)
        IARB=1
        IF(BRAVAIS.EQ.'GH'.OR.BRAVAIS.EQ.'GQ'.OR.BRAVAIS.EQ.'GOB') THEN
         IARB(1)=1
         IARB(2)=0
         IARB(3)=0
        ENDIF 
        IF(BRAVAIS.EQ.'GOF'.OR.BRAVAIS.EQ.'GO'.OR.BRAVAIS.EQ.'GM') THEN
          IARB=0
        ENDIF 
        IF(BRAVAIS.EQ.'GMB') THEN
         IARB(1)=0
         IARB(2)=1
         IARB(3)=0
        ENDIF 
                            CALL TRACE$PUSH('BRILLOUIN$MSH')

        CALL SPACEGROUP$GENERATORS('RECI',NOPX,NOP,IIO,C)
        IF(THISTASK.EQ.1)WRITE(*,FMT='(82("="),T10," GENERATORS OF THE GROUP ")')
        DO I=1,NOP
          IF(THISTASK.EQ.1)WRITE(*,FMT='(I5,T20,3("|",3I5,"|"))')I,IIO(:,:,I)
        ENDDO
        CALL BRILLOUIN$MSH(RBAS,NKP,NOP,IIO,IARB,TSHIFT)
                            CALL TRACE$POP()
      ELSE
                            CALL TRACE$PUSH('BRILLOUIN$MSHNOSYM')
!     ==========================================================================
!     ==  FIND K-POINTS AND TETRAHEDRA                                        ==
!     ==========================================================================
        CALL BRILLOUIN$MSHNOSYM(TINV,RBAS,NKDIV,ISHIFT)
                            CALL TRACE$POP()
      ENDIF
! 
!     ==========================================================================
!     ==  CALCULATE ENERGIES AT THE IRREDUCIBLE K-POINTS                      ==
!     ==========================================================================
      CALL BRILLOUIN$GETI4('NK',NKP)
!     ==========================================================================
!     ==  CALCULATE WEIGHTS                                                   ==
!     ==========================================================================
      IF(TPRINT)WRITE(*,*)' NKP NACH BRILLOUIN$GETI4',NKP  
      ALLOCATE(BK(3,NKP))
      ALLOCATE(XK(3,NKP))
      ALLOCATE(WKPT(NKP))
      CALL BRILLOUIN$GETR8A('K',3*NKP,BK)
      CALL BRILLOUIN$GETR8A('XK',3*NKP,XK)
      CALL BRILLOUIN$GETR8A('WKPT',NKP,WKPT)
!
!     =========================================================================
!     ==  CONSTRUCT K-INDEPENDENT PART OF HAMILTONIAN                        ==
!     =========================================================================
      CALL BANDS_KINDEP(NG,GVEC,G2,GWEIGHT,GIDG_PROTO,TI_H,TI_S)

      ALLOCATE(EB(NB*NSPIN,NKP))        
      ALLOCATE(WGHT(NB*NSPIN,NKP))
      EB(:,:)=0.0D0 
      IF(TPROJ)THEN
        ALLOCATE(PROJK(NAT,NB*NSPIN,LMNXX,NKP))
        PROJK(:,:,:,:)=CMPLX(0.D0,0.D0,KIND=8)
      ENDIF

!     ==========================================================================
!     ==  READ BEGINNING OF PDOS FILE                                         ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PDOS',NFILIN)
      REWIND(NFILIN)
      CALL PDOS$READ(NFILIN)
      CALL PDOS$SETI4('NSP',NSP)
      CALL PDOS$SETI4('NAT',NAT)
      CALL PDOS$SETI4('NKPT',NKP)
      CALL PDOS$SETI4('NSPIN',NSPIN)
      CALL PDOS$SETI4('NDIM',NDIM)
      CALL PDOS$SETI4('NPRO',NPRO)
      CALL PDOS$SETI4A('NKDIV',3,NKDIV)
      CALL PDOS$SETI4A('ISHIFT',3,ISHIFT)
      CALL PDOS$SETR8('RNTOT',RNTOT)
      CALL PDOS$SETR8('NEL',NEL)
      CALL PDOS$SETL4('TINV',TINV)
      CALL PDOS$SETI4('SPACEGROUP',SPACEGROUP)
      CALL PDOS$SETL4('TSHIFT',TSHIFT)
!     ==========================================================================
!     ==  WRITE BEGINNING OF PDOS FILE                                        ==
!     ==========================================================================
      
      IF(THISTASK.EQ.1)THEN
        ALLOCATE(OCC(NKP,NB))
        OCC(:,:)=1.0D0
        CALL FILEHANDLER$UNIT('PDOSOUT',NFILOUT)
        REWIND NFILOUT
        CALL PDOS$WRITE(NFILOUT,'181213')
      ENDIF
!       == ITERATE K-POINTS =================================================
                            CALL TRACE$PUSH('ITERATE KPOINTS')
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      ALLOCATE(JOBLIST(NKP,NSPIN))
      INDEX=1
      DO IKP=1,NKP
        DO ISPIN=1,NSPIN
          JOBLIST(IKP,ISPIN)=INDEX
          INDEX=INDEX+1
          IF(INDEX.GT.NTASKS)INDEX=1
        ENDDO
      ENDDO
      PRINT*,"JOBLIST",JOBLIST(:,:)
      


      ALLOCATE(PROJ(NAT,NB,LMNXX))
      ALLOCATE(E(NG*NDIM))
      DO IKP=1,NKP
        DO ISPIN=1,NSPIN
          IF(THISTASK.NE.JOBLIST(IKP,ISPIN))CYCLE
          KVEC=BK(:,IKP)
          IF(TPRINT)WRITE(NFILO,*)'K-POINT ',IKP,'(',INDEX,') OF ',NKP,&
     &        ' FOR SPIN ',ISPIN,' IN ABSOLUTE COORDINATES ',KVEC
          CALL CPU_TIME(T1)
          CALL BANDS_KPOINT(NG,NB,ISPIN,METHOD_DIAG,GIDG_PROTO,TPROJ,KVEC,&
     &              GVEC,TI_H,TI_S,E,PROJ)
          CALL CPU_TIME(T2)
          PRINT*,'CPU-TIME FOR KPOINT: NODE=',THISTASK,' KPOINT:',IKP,&
     &         ' SPIN=',ISPIN,'TIME=',T2-T1
          EB(1+NB*(ISPIN-1):NB+NB*(ISPIN-1),IKP)=E(1:NB)
          IF(TPROJ)THEN
            DO IAT=1,NAT
              DO LMN=1,LMNXX
                PROJK(IAT,1+NB*(ISPIN-1):NB+NB*(ISPIN-1),LMN,IKP)=PROJ(IAT,1:NB,LMN)
              ENDDO         
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(PROJ)
      DEALLOCATE(E)
      CALL MPE$COMBINER8R2('~','+',EB)
      CALL MPE$COMBINEC8R4('~','+',PROJK)
                            CALL TRACE$POP()
      !
      IF(THISTASK.NE.1)RETURN

      DO IKP=1,NKP
        DO I=1,NB
          PRINT*,'EB(',I,',',IKP,')=',EB(I,IKP)          
        ENDDO
      ENDDO
        
!     ==========================================================================
!     ==  CALCULATE WEIGHTS                                                   ==
!     ==========================================================================
      IF(NSPIN.EQ.2) WRITE(NFILO,*)'WARNING: USING THE SAME FERMI ENERGY FOR&
     &      BOTH SPIN DIRECTIONS.'
      CALL BRILLOUIN$DOS(NSPIN*NB,NKP,EB,WGHT,RNTOT,EF)
                            
                        CALL TRACE$PUSH('PDOS WRITEK LOOP')
!     ==========================================================================
!     ==  WRITE TO PDOS FILE                                                  ==
!     ==========================================================================
      ALLOCATE(PROJTMP(NDIM,NPRO,NB))
      DO IKP=1,NKP
        DO ISPIN=1,NSPIN
          KVEC=BK(:,IKP)
          WRITE(NFILO,*)'K-POINT ',IKP,' OF ',NKP*NSPIN,' FOR SPIN ',ISPIN,&
     &               ' IN ABSOLUTE COORDINATES ',KVEC
          IF(TPROJ)THEN
            !REORDER PROJECTIONS
            INDEX=1
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              LMNX_=LMNX(ISP)
              DO LMN=1,LMNX_
                PROJTMP(1,INDEX,:)=PROJK(IAT,1+NB*(ISPIN-1):NB+NB*(ISPIN-1),&
      &                            LMN,IKP)
                INDEX=INDEX+1 
              ENDDO
            ENDDO
          ELSE
            PROJTMP(:,:,:)=0.0D0
          ENDIF
          DO IB=1,NB
            IF(EB(IB+NB*(ISPIN-1),IKP).LE.EF)THEN
              OCC(IKP,IB)=WKPT(IKP)*2.0D0/REAL(NSPIN,KIND=8)
            ELSE
              OCC(IKP,IB)=0.0D0
            ENDIF
          ENDDO
          CALL PDOS$WRITEK(NFILOUT,XK(:,IKP),NB,NDIM,NPRO,&
      &  WKPT(IKP),EB(1+NB*(ISPIN-1):NB+NB*(ISPIN-1),IKP),OCC(IKP,1:NB),PROJTMP)
        ENDDO
      ENDDO
      DEALLOCATE(PROJTMP)
!     ==========================================================================
!     ==  CLOSE PDOS FILE                                                     ==
!     ==========================================================================
                            CALL TRACE$PUSH('PDOS CLOSE')
      FLUSH(NFILOUT)
      CALL FILEHANDLER$CLOSE('PDOSOUT')
                            CALL TRACE$POP()     
!$!
!$!     ==========================================================================
!$!     ==  CALCULATE AND WRITE TOTAL DENSITY OF STATES                         ==
!$!     ==  THE FOLLOWING CODE WILL BE REMOVED, WHEN SYMMETRIES ARE DONE!       ==
!$!     ==========================================================================
!$      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
!$      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
!$      CALL LINKEDLIST$NLISTS(LL_CNTL,'DOS',NPDOS)
!$      IF(NPDOS.EQ.1)THEN
!$                        CALL TRACE$PUSH('TOTAL DOS')
!$        CALL LINKEDLIST$SELECT(LL_CNTL,'DOS',1)
!$        EMIN=MINVAL(EB(:,:))
!$        EMAX=MAXVAL(EB(:,:))
!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'EMIN',0,TCHK)
!$        IF(TCHK)THEN
!$          CALL LINKEDLIST$GET(LL_CNTL,'EMIN',1,EMIN2)
!$        ELSE
!$          EMIN2=EMIN*27.211D0
!$        ENDIF
!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'EMAX',0,TCHK)
!$        IF(TCHK)THEN
!$          CALL LINKEDLIST$GET(LL_CNTL,'EMAX',1,EMAX2)
!$        ELSE
!$          EMAX2=EMAX*27.211D0
!$        ENDIF
!$        
!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'NE',0,TCHK)
!$        IF(TCHK)THEN
!$          CALL LINKEDLIST$GET(LL_CNTL,'NE',1,NE)
!$        ELSE
!$          NE=1000
!$        ENDIF
!$        
!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',0,TCHK)
!$        IF(TCHK)THEN
!$          CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FILE)
!$        ELSE
!$          CALL ERROR$MSG('FILE IN DOS-BLOCK NOT GIVEN')
!$          CALL ERROR$STOP('BANDS_PDOS')
!$        ENDIF
!$        ALLOCATE(EWGHT(NSPIN*NB,NKP))
!$        CALL BRILLOUIN$EWGHT(NKP,NB*NSPIN,EB,EMIN,EMAX,NE,EWGHT)
!$
!$        ALLOCATE(DOS(NE,NSPIN))
!$        DOS(:,:)=0.0D0
!$        DO IKP=1,NKP
!$          DO IB=1,NB
!$            DO ISPIN=1,NSPIN
!$              I1=EWGHT(IB+NB*(ISPIN-1),IKP)%I1
!$              I2=EWGHT(IB+NB*(ISPIN-1),IKP)%I2
!$              DO IE=I1,I2
!$                DOS(IE,ISPIN)=DOS(IE,ISPIN)+EWGHT(IB+NB*(ISPIN-1),IKP)%WGHT(IE)
!$              ENDDO
!$            ENDDO
!$          ENDDO
!$        ENDDO
!$
!$        WRITE(NFILO,*)'WRITING TOTAL DOS TO FILE:', FILE
!$        OPEN(101,FILE=FILE)
!$        
!$        DO IE=1,NE
!$          EI=27.211D0*(EMIN+REAL(IE-1,KIND=8)*(EMAX-EMIN)/REAL(NE-1,KIND=8))
!$          IF((EI.GE.EMIN2).AND.(EI.LE.EMAX2))THEN
!$            WRITE(101,*)EI,(DOS(IE,ISPIN),ISPIN=1,NSPIN)
!$          ENDIF
!$        ENDDO
!$        CLOSE(101)        
!$                        CALL TRACE$POP()
!$      ENDIF
                            CALL TRACE$POP()
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LAPACK_ZHEGVD(N,JOBZ,H,S,E)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: N
      CHARACTER(1),INTENT(IN)  :: JOBZ
      COMPLEX(8),INTENT(INOUT) :: H(N,N)
      COMPLEX(8),INTENT(INOUT) :: S(N,N)
      REAL(8),INTENT(OUT)      :: E(N)
      INTEGER(4)               :: LWORK,LRWORK,LIWORK,INFO
      COMPLEX(8),ALLOCATABLE   :: WORK(:)
      REAL(8)   ,ALLOCATABLE   :: RWORK(:)
      INTEGER(4),ALLOCATABLE   :: IWORK(:)
!     **************************************************************************

      LWORK=-1
      LRWORK=-1
      LIWORK=-1
      ALLOCATE(WORK(1))!COMPLEX(8)
      ALLOCATE(RWORK(1))!REAL(8)
      ALLOCATE(IWORK(1))!INTEGER(4)
      
      CALL ZHEGVD(1,JOBZ,'U',N,H,N,S,N,E &
     &           ,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
      LWORK=WORK(1)
      LRWORK=RWORK(1)
      LIWORK=IWORK(1)
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(IWORK)
      ALLOCATE(WORK(LWORK)) 
      ALLOCATE(RWORK(LRWORK)) 
      ALLOCATE(IWORK(LIWORK)) 
      CALL ZHEGVD(1,JOBZ,'U',N,H,N,S,N,E &
     &           ,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(IWORK)

      IF(INFO.NE.0)THEN
        CALL ERROR$MSG('ZHEGVD FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LAPACK_ZHEGVD')
      ENDIF
      RETURN
      END
!
#IF DEFINED(CPPVAR_FEAST)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FEAST_ZHEGV(N,NB,EMIN,EMAX,TINITIALGUESS,H,S,E,U)
!     **************************************************************************
!     ** INTERFACE TO THE SPARSE-MATRIX DIAGONALIZATION                       **
!     **                                                                      **
!     **  HTTP://WWW.ECS.UMASS.EDU/~POLIZZI/FEAST/                            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: N
      INTEGER(4)  ,INTENT(IN)  :: NB
      REAL(8),INTENT(IN)       :: EMIN,EMAX
      LOGICAL(4),INTENT(IN)    :: TINITIALGUESS
      COMPLEX(8),INTENT(INOUT) :: H(N,N)
      COMPLEX(8),INTENT(INOUT) :: S(N,N)
      REAL(8),INTENT(OUT)      :: E(N)
      COMPLEX(8),INTENT(INOUT) :: U(N,N)
      INTEGER(4)               :: NMAT
      CHARACTER*1,PARAMETER    :: UPLO='F'
      INTEGER(4)               :: FPM(64),M
      REAL(8)                  :: EPSOUT
      INTEGER(4)               :: INFO,LOOP,M0
      REAL(8)                  :: RES(N)
!     **************************************************************************
      M0=NINT(2.0D0*REAL(NB,KIND=8))
      NMAT=N

      CALL FEASTINIT(FPM)
      FPM(1)=1 !PRINT RUNTIME COMMENTS
      FPM(2)=3 !NUMBER OF CONTOUR POINTS 
      FPM(3)=0 !STOPPING CRITERIA FOR DOUBLE PRECISION EPS=10^(-FPM(3))
!     FPM(4)=0 !MAXIMUM NUMBER OF REFINEMENT LOOPS
      IF(TINITIALGUESS)THEN
        FPM(5)=1 !PROVIDE INITIAL GUESS SUBSPACE
      ELSE
        FPM(5)=0 !PROVIDE INITIAL GUESS SUBSPACE
        U=CMPLX(0.D0,0.D0,KIND=8)
      ENDIF
!     __CHOOSE CONVERGENCE CRITERIA:____________________________________________
!     _________0 RELATIVE ERROR ON TRACE (EPSOUT<EPS)
!     _________1 RELATIVE RESIDUAL (RES<EPS)
      FPM(6)=0

      M=0
      E=0.D0
      RES=0.D0
      EPSOUT=0.D0
      CALL ZFEAST_HEGV(UPLO,NMAT,H,NMAT,S,NMAT,FPM,EPSOUT&
                   &,LOOP,EMIN/27.211D0, EMAX/27.211D0,M0,E,U,M,RES,INFO) 

      IF(INFO.NE.0)THEN
        CALL ERROR$MSG('FEAST_ZHEGV FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('FEAST_ZHEGV')
      ENDIF
      RETURN
      END
#ENDIF
!
#IF DEFINED(CPPVAR_JADAMILU)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE JADAMILU_ZHEGV(N,NB,H,S,E)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: N
      INTEGER(4)  ,INTENT(IN)  :: NB
      COMPLEX(8),INTENT(INOUT) :: H(N,N)
      COMPLEX(8),INTENT(INOUT) :: S(N,N)
      REAL(8)   ,INTENT(OUT)   :: E(N)
      INTEGER(4)               :: INFO
      COMPLEX(8),ALLOCATABLE   :: H_CSR(:)
      INTEGER(4)               :: IH_CSR(N+1)
      INTEGER(4),ALLOCATABLE   :: JH_CSR(:)
      COMPLEX(8),ALLOCATABLE   :: S_CSR(:)
      INTEGER(4)               :: IS_CSR(N+1)
      INTEGER(4),ALLOCATABLE   :: JS_CSR(:)
      INTEGER(4)               :: IERR,IPRINT,NEIG,NINIT,ITER,ICNTL(5)
      REAL(8)                  :: TOL,RES(NB),SIGMA,SHIFT,MEM,DROPTOL,GAP
      INTEGER(4)               :: LX,MAXSP,ISEARCH,MADSPACE
      COMPLEX(8),ALLOCATABLE   :: X(:)
      REAL(8)                  :: EIGS(NB)
      INTEGER(4)               :: I,J,NCSR
!     ******************************************************************
      !CONVERT DENSE MATRIX TO CSR
      NCSR=(N*N+N)/2
      ALLOCATE(H_CSR(NCSR))
      ALLOCATE(S_CSR(NCSR))
      ALLOCATE(JH_CSR(NCSR))
      ALLOCATE(JS_CSR(NCSR))

      CALL DNSCSR (1D-6, N, N, N*N, H, N, H_CSR, JH_CSR, IH_CSR, IERR )
      CALL DNSCSR (1D-6, N, N, N*N, S, N, S_CSR, JS_CSR, IS_CSR, IERR )
      PRINT*, IERR
      !DIAGONAL PREPROCESSING
      JH_CSR(1)=-1
      JS_CSR(1)=-1

      IPRINT=6
      NEIG=NB
      NINIT=0
      ITER=10000
      TOL=1.0D-4
      ICNTL(1)=0
      ICNTL(2)=0
      ICNTL(3)=0
      ICNTL(4)=0
      ICNTL(5)=1
      MAXSP=20
      ISEARCH=1
      SIGMA=0.0D0
      SHIFT=SIGMA
      MEM=20.0D0
      DROPTOL=1.D-3
      MADSPACE=MAXSP

      LX=N*(4*MAXSP+2*NEIG+1)+4*MAXSP*MAXSP
      ALLOCATE(X(LX))

      CALL ZPJD_GEP(N,H_CSR,JH_CSR,IH_CSR,S_CSR,JS_CSR,IS_CSR,&
        &   EIGS,RES,X,LX,NEIG,SIGMA,ISEARCH,NINIT,&
        &   MADSPACE,ITER,TOL,SHIFT,DROPTOL,MEM,ICNTL,IPRINT,INFO,GAP)
      E(1:NB)=EIGS
      IF(INFO.NE.0)THEN
        CALL ERROR$MSG('ZHEGV FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('JADAMILU_ZHEGV')
      ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DNSCSR(TOL,NROW,NCOL,NZMAX,DNS,NDNS,A,JA,IA,IERR )
!     **************************************************************************
!     ** DNSCSR CONVERTS DENSE TO COMPRESSED ROW SPARSE FORMAT.               **
!     **                                                                      **
!     ** DISCUSSION:                                                          **
!     **                                                                      **
!     ** THIS ROUTINE CONVERTS A DENSELY STORED MATRIX INTO A ROW ORIENTIED   **
!     ** COMPACTLY SPARSE MATRIX.  IT IS THE REVERSE OF CSRDNS.               **
!     **                                                                      **
!     ** THIS ROUTINE DOES NOT CHECK WHETHER AN ELEMENT IS SMALL. IT CONSIDERS**
!     ** THAT A(I,J) IS ZERO ONLY IF IT IS EXACTLY EQUAL TO ZERO.             **
!     **                                                                      **
!     ** MODIFIED: 07 JANUARY 2004                                            **
!     ** AUTHOR:   YOUCEF SAAD                                                **
!     ** PART OF JADAMILU                                                    **
!     **                                                                      **
!     ** PARAMETERS:                                                          **
!     **                                                                      **
!     **  INPUT,INTEGER(KIND=4) NROW: THE ROW DIMENSION OF THE MATRIX.        **
!     **                                                                      **
!     **  INPUT,INTEGER(KIND= 4) NCOL: THE COLUMN DIMENSION OF THE MATRIX.    **
!     **                                                                      **
!     **  INPUT,INTEGER(KIND=4) NZMAX: THE MAXIMUM NUMBER OF NONZERO ELEMENTS **
!     **                               ALLOWED. THIS SHOULD BE SET TO BE THE  **
!     **                               LENGTHS OF THE ARRAYS A AND JA.        **
!     **                                                                      **
!     **  INPUT,REAL DNS(NDNS,NCOL): AN NROW BY NCOL DENSE MATRIX.            **
!     **                                                                      **
!     **  INPUT,INTEGER(KIND=4) NDNS: THE FIRST DIMENSION OF DNS, WHICH MUST  **
!     **                              BE AT LEAST NROW.                       **
!     **                                                                      **
!     **  OUTPUT,REAL A(*), INTEGER(KIND=4) JA(*), IA(NROW+1): THE MATRIX IN  **
!     **                                CSR COMPRESSED SPARSE ROW FORMAT.     **
!     **                                                                      **
!     **  OUTPUT, INTEGER ( KIND = 4 ) IERR: ERROR INDICATOR.                 **
!     **            0 MEANS NORMAL RETURN;
!     **            I, MEANS THAT THE THE CODE STOPPED WHILE PROCESSING ROW I,**
!     **               BECAUSE THERE WAS NO SPACE LEFT IN A AND JA,           **
!     **               AS DEFINED BY NZMAX.                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER ( KIND = 4 ),INTENT(IN) :: NCOL
      INTEGER ( KIND = 4 ),INTENT(IN) :: NDNS
      INTEGER ( KIND = 4 ),INTENT(IN) :: NROW
      COMPLEX ( KIND = 8 ),INTENT(IN) :: DNS(NDNS,NCOL)
      COMPLEX ( KIND = 8 ),INTENT(OUT)::  A(*)
      INTEGER ( KIND = 4 ),INTENT(OUT):: IA(NROW+1)
      INTEGER ( KIND = 4 ),INTENT(OUT):: JA(*)
      INTEGER ( KIND = 4 ) :: IERR
      INTEGER ( KIND = 4 ) :: I
      INTEGER ( KIND = 4 ) :: J
      INTEGER ( KIND = 4 ) :: NEXT
      INTEGER ( KIND = 4 ) :: NZMAX
      REAL ( KIND = 8 )    :: TOL
!     **************************************************************************
      IERR = 0
      NEXT = 1
      IA(1) = 1
      DO I = 1, NROW
        DO J = I, NCOL
          IF ( ABS(DNS(I,J)).GT.TOL ) THEN
            IF ( NZMAX < NEXT ) THEN
              IERR = I
              RETURN
            END IF
            JA(NEXT) = J
            A(NEXT) = DNS(I,J)
            NEXT = NEXT + 1
          END IF
        END DO
        IA(I+1) = NEXT
      END DO
      RETURN
      END
#ENDIF
!
#IF DEFINED(CPPVAR_SLEPC)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SLEPC_ZHEGV(DIM,NB,HK,SK,E)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
#INCLUDE <FINCLUDE/SLEPCEPSDEF.H>
      USE SLEPCEPS
      USE ISO_C_BINDING
      IMPLICIT NONE
! FOR USAGE WITHOUT MODULES, UNCOMMENT THE FOLLOWING LINES AND REMOVE 
! THE PREVIOUS LINES BETWEEN 'PROGRAM MAIN' AND 'IMPLICIT NONE'
!
!#INCLUDE <FINCLUDE/PETSC.H>
!#INCLUDE <FINCLUDE/SLEPC.H>
      INTEGER(4)  ,INTENT(IN)  :: DIM
      INTEGER(4)  ,INTENT(IN)  :: NB
      COMPLEX(8),INTENT(INOUT) :: HK(DIM,DIM)
      COMPLEX(8),INTENT(INOUT) :: SK(DIM,DIM)
      REAL(8),INTENT(OUT)      :: E(DIM)
      INTEGER(4)               :: INFO=0
#IF DEFINED(PETSC_USE_FORTRAN_DATATYPES)
      TYPE(MAT)   ::   H2
      TYPE(MAT)   ::   S2
      TYPE(EPS)   ::   SOLVER
#ELSE
      MAT         ::   H2
      MAT         ::   S2
      EPS         ::   SOLVER
#ENDIF
      EPSTYPE     ::   TNAME
      PETSCINT    ::   N, I,J, ISTART, IEND, ONE, TWO, THREE
      PETSCINT    ::   NEV
      PETSCINT    ::   ROW(1)
      PETSCINT    ::   COL(DIM)
      PETSCMPIINT ::   RANK
      PETSCERRORCODE :: IERR
      PETSCBOOL     :: FLG
      PETSCINT       :: NNZ(DIM)
      COMPLEX(8)     :: CSVAR1,CSVAR2
      REAL(8)        :: TOL=1D-3
!     **************************************************************************
      ONE = 1
      TWO = 2
      THREE = 3
      CALL SLEPCINITIALIZE(PETSC_NULL_CHARACTER,IERR)
      CALL MPI_COMM_RANK(PETSC_COMM_WORLD,RANK,IERR)
      N = DIM
      !CALL PETSCOPTIONSGETINT(PETSC_NULL_CHARACTER,'-N',N,FLG,IERR)
      
      CALL MATCREATE(PETSC_COMM_WORLD,H2,IERR)
      CALL MATSETSIZES(H2,PETSC_DECIDE,PETSC_DECIDE,N,N,IERR)
      CALL MATSETFROMOPTIONS(H2,IERR)
      CALL MATSETUP(H2,IERR)
      DO I=1,DIM
        NNZ(I)=0
        DO J=1,DIM
          IF(ABS(HK(I,J)).GT.TOL)NNZ(I)=NNZ(I)+1
        ENDDO
      ENDDO
      PRINT*,"SUM",SUM(NNZ(:)),DIM*DIM &
     &            ,REAL(SUM(NNZ(:)),KIND=8)/REAL(DIM*DIM,KIND=8)
      CALL MATSEQAIJSETPREALLOCATION(H2,DIM,NNZ,IERR)
      
      CALL MATGETOWNERSHIPRANGE(H2,ISTART,IEND,IERR)
      PRINT*, DIM, ISTART,IEND
!      DO I=1,DIM
!        COL(I)=I-1
!      ENDDO
      DO I=ISTART,IEND-1
        !ROW(1)=I
        !CALL MATSETVALUES(H2,ONE,ROW,DIM,COL,HK(:,I+1),INSERT_VALUES,IERR)
        DO J=0,DIM-1
          IF(ABS(HK(I+1,J+1)).GT.TOL)THEN
            CALL MATSETVALUE(H2,I,J,HK(I+1,J+1),INSERT_VALUES,IERR)
          ENDIF
        ENDDO
      ENDDO
      CALL MATASSEMBLYBEGIN(H2,MAT_FINAL_ASSEMBLY,IERR)
      CALL MATASSEMBLYEND(H2,MAT_FINAL_ASSEMBLY,IERR)
      
      
      CALL MATCREATE(PETSC_COMM_WORLD,S2,IERR)
      CALL MATSETSIZES(S2,PETSC_DECIDE,PETSC_DECIDE,N,N,IERR)
      CALL MATSETFROMOPTIONS(S2,IERR)
      CALL MATSETUP(S2,IERR)
      DO I=1,DIM
        NNZ(I)=0
        DO J=1,DIM
          IF(ABS(SK(I,J)).GT.TOL)NNZ(I)=NNZ(I)+1
        ENDDO
      ENDDO
      PRINT*,"SUM",SUM(NNZ(:)),DIM*DIM &
    &             ,REAL(SUM(NNZ(:)),KIND=8)/REAL(DIM*DIM,KIND=8)
      CALL MATSEQAIJSETPREALLOCATION(S2,DIM,NNZ,IERR)
      
      CALL MATGETOWNERSHIPRANGE(S2,ISTART,IEND,IERR)
      PRINT*, DIM, ISTART,IEND
!      DO I=1,DIM
!        COL(I)=I-1
!      ENDDO
      DO I=ISTART,IEND-1
        !ROW(1)=I
        !CALL MATSETVALUES(S2,ONE,ROW,DIM,COL,HK(:,I+1),INSERT_VALUES,IERR)
        DO J=0,DIM-1
          IF(ABS(SK(I+1,J+1)).GT.TOL)THEN
            CALL MATSETVALUE(S2,I,J,SK(I+1,J+1),INSERT_VALUES,IERR)
          ENDIF
        ENDDO
      ENDDO
      CALL MATASSEMBLYBEGIN(S2,MAT_FINAL_ASSEMBLY,IERR)
      CALL MATASSEMBLYEND(S2,MAT_FINAL_ASSEMBLY,IERR)
      


!      CALL MATVIEW(H2,      PETSC_VIEWER_STDOUT_SELF)
!      CALL MATVIEW(S2,      PETSC_VIEWER_STDOUT_SELF)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     CREATE THE EIGENSOLVER AND DISPLAY INFO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     ** CREATE EIGENSOLVER CONTEXT
      CALL EPSCREATE(PETSC_COMM_WORLD,SOLVER,IERR)
      CALL EPSSETTYPE(SOLVER,EPSARPACK,IERR)
!      CALL EPSSETTOLERANCES(SOLVER,1D-3,0,IERR)
      CALL EPSSETWHICHEIGENPAIRS(SOLVER,EPS_SMALLEST_REAL,IERR)
      CALL EPSSETDIMENSIONS(SOLVER,NB,PETSC_DECIDE,PETSC_DECIDE,IERR)
!     ** SET OPERATORS. IN THIS CASE, IT IS A STANDARD EIGENVALUE PROBLEM
      CALL EPSSETOPERATORS(SOLVER,H2,S2,IERR)
      CALL EPSSETPROBLEMTYPE(SOLVER,EPS_GHEP,IERR)
!     ** SET SOLVER PARAMETERS AT RUNTIME
      !CALL EPSSETFROMOPTIONS(SOLVER,IERR)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     SOLVE THE EIGENSYSTEM
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      PRINT*,"SOLVING..."
      CALL EPSSOLVE(SOLVER,IERR) 

!     ** OPTIONAL: GET SOME INFORMATION FROM THE SOLVER AND DISPLAY IT
      CALL EPSGETTYPE(SOLVER,TNAME,IERR)
      IF (RANK .EQ. 0) THEN
        WRITE(*,120) TNAME
      ENDIF
 120  FORMAT (' SOLUTION METHOD: ',A)
      CALL EPSGETDIMENSIONS(SOLVER,NEV,PETSC_NULL_INTEGER &
     &                     ,PETSC_NULL_INTEGER,IERR)
      IF (RANK .EQ. 0) THEN
        WRITE(*,130) NEV
      ENDIF
 130  FORMAT (' NUMBER OF REQUESTED EIGENVALUES:',I4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     DISPLAY SOLUTION AND CLEAN UP
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CALL EPSPRINTSOLUTION(SOLVER,PETSC_NULL_OBJECT,IERR)
      CALL EPSGETCONVERGED(SOLVER,NEV,IERR)
      DO I=0,NEV-1
        CALL EPSGETEIGENVALUE(SOLVER,I,CSVAR1,CSVAR2,IERR)
        E(I+1)=REAL(CSVAR1,KIND=8)
      ENDDO
      CALL EPSDESTROY(SOLVER,IERR)
      CALL MATDESTROY(H2,IERR)
      CALL MATDESTROY(S2,IERR)

      CALL SLEPCFINALIZE(IERR)

      IF(INFO.NE.0)THEN
        CALL ERROR$MSG('ZHEGV FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('SLEPC_ZHEGV')
      ENDIF
      END SUBROUTINE SLEPC_ZHEGV
#ENDIF
