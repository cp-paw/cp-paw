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
      PROGRAM MAIN
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE PDOS_MODULE, ONLY: STATE,STATEARR
      IMPLICIT NONE
      TYPE(LL_TYPE)   :: LL_CNTL
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NAT
      INTEGER(4)                :: NB
      INTEGER(4)                :: NKPT
      INTEGER(4)                :: NKPTBIG
      INTEGER(4)                :: NSPIN
      INTEGER(4)                :: NDIM !=2 FOR SPINOR WF; OTHERWISE =1
      REAL(8)                   :: RBAS(3,3) ! REAL-SPACE LATTICE VECTORS
      REAL(8)                   :: GBAS(3,3) ! REAL-SPACE LATTICE VECTORS
      REAL(8)      ,ALLOCATABLE :: RPOS(:,:)
      REAL(8)      ,ALLOCATABLE :: XK(:,:) ! K IN RELATIVE COORDINATES
      REAL(8)      ,ALLOCATABLE :: EB(:,:,:)
      REAL(8)      ,ALLOCATABLE :: EBBIG(:,:,:)
      INTEGER(4)   ,ALLOCATABLE :: MAP(:,:,:)
      INTEGER(4)   ,ALLOCATABLE :: MAPBIG(:,:,:)
      REAL(8)      ,ALLOCATABLE :: WORK1(:,:,:)
      REAL(8)      ,ALLOCATABLE :: WORK2(:,:,:)
      INTEGER(4)                :: NFILIN
      INTEGER(4)                :: NPRO
      INTEGER(4)                :: IKPT,ISPIN,IB,I,J,K,j1,j2
      INTEGER(4)                :: N1,N2,N3
      INTEGER(4)                :: N1B,N2B,N3B
      INTEGER(4)                :: IND
      INTEGER(4)   ,ALLOCATABLE :: NBARR(:,:)
      REAL(8)                   :: XK1(3),XK2(3),XQ(3)
      CHARACTER(512)            :: FILE
      INTEGER(4)                :: NLINE,ILINE
      INTEGER(4)                :: NQ
      INTEGER(4)                :: NP
      INTEGER(4)                :: NFIL
      REAL(8)                   :: X1,X2
      REAL(8)                   :: SVAR,svar1
      LOGICAL(4)                :: TCHK,TCHK1
      LOGICAL(4),PARAMETER      :: TREFINE=.FALSE.
!     **************************************************************************
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
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==========================================================================
!     ==  WRITE HEADER                                                        ==
!     ==========================================================================
      WRITE(NFILO,FMT='(72("*"))')
      WRITE(NFILO,FMT='(72("*"),T15 &
     &             ,"           BANDS ANALYSIS TOOL                ")')
      WRITE(NFILO,FMT='(72("*"),T15 &
     &             ,"    FOR THE PROJECTOR-AUGMENTED WAVE METHOD  ")')
      WRITE(NFILO,FMT='(72("*"))')
      WRITE(NFILO,FMT='(T28 &
     &           ," P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY ")')
      WRITE(NFILO,FMT='(T10 &
     &      ,"(C) CLAUSTHAL UNIVERSITY OF TECHNOLOGY (CUT), GERMANY " &
     &      ,"* ANY USE REQUIRES WRITTEN LICENSE FROM CUT")')
      WRITE(NFILO,*)
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
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          DO IB=1,NB
            EB(IB,IKPT,ISPIN)=STATE%EIG(IB)
!PRINT*,'STATE ',STATE%VEC(:,:,IB)
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
         svar=1.d0
         do j1=1,nkpt
           do j2=j1+1,nkpt
             svar1=abs(modulo(xk(i,j1)-xk(i,j2)+0.5d0,1.d0)-0.5d0)
             if(svar1.ne.0.d0) svar=min(svar,svar1)
           enddo
         enddo
         svar=1.d0/svar
         IF(I.EQ.1) THEN
           N1=NINT(SVAR)
!!$PRINT*,'N1 ',N1,svar
         ELSE IF(I.EQ.2) THEN
           N2=NINT(SVAR)
!!$PRINT*,'N2 ',N2,svar
         ELSE IF(I.EQ.3) THEN
           N3=NINT(SVAR)
!!$PRINT*,'N3 ',N3,svar
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
       CALL LINKEDLIST$NLISTS(LL_CNTL,'LINE',NLINE)
       X2=0.D0
       DO ILINE=1,NLINE
         CALL LINKEDLIST$SELECT(LL_CNTL,'LINE',ILINE)
!
         CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',0,TCHK)
         IF(TCHK) THEN
           CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FILE)
           CALL FILEHANDLER$SETFILE('BANDS',.FALSE.,FILE)
           CALL FILEHANDLER$SETSPECIFICATION('BANDS','STATUS','REPLACE')
           CALL FILEHANDLER$SETSPECIFICATION('BANDS','POSITION','APPEND')
           CALL FILEHANDLER$SETSPECIFICATION('BANDS','ACTION','WRITE')
           CALL FILEHANDLER$SETSPECIFICATION('BANDS','FORM','FORMATTED')
           X1=0.D0
           X2=1.D0
         ELSE
           X1=X2
           X2=X1+1.D0
         END IF
!
         CALL LINKEDLIST$GET(LL_CNTL,'XK1',1,XK1)
         CALL LINKEDLIST$GET(LL_CNTL,'XK2',1,XK2)
         NP=100
         CALL LINKEDLIST$EXISTD(LL_CNTL,'NK',0,TCHK)
         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NK',1,NP)
!        == PROJECTION FOR 2D- BANDSTRUCTURE ===================================
         CALL LINKEDLIST$EXISTD(LL_CNTL,'XKPROJECT',0,TCHK)
         IF(TCHK) THEN
           CALL LINKEDLIST$GET(LL_CNTL,'XKPROJECT',1,XQ)
           NQ=11
           CALL LINKEDLIST$EXISTD(LL_CNTL,'NPROJECT',0,TCHK1)
           IF(TCHK1)CALL LINKEDLIST$GET(LL_CNTL,'NPROJECT',1,NQ)
           IF(MODULO(N1,2).NE.1) THEN
!            == for 2-d- band structures one goes forward-backward-forward =====
!            == for even nq a straight line is drawn towards the next segment ==
!            == when several lines are written in the same file=================
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
           CALL BANDS_PLOTBANDS(GBAS,N1B,N2B,N3B,MAPBIG &
      &                        ,NKPTBIG,NB,EBBIG(:,:,ISPIN) &
      &                        ,NP,X1,XK1,X2,XK2,NQ,XQ)
         ELSE
           CALL BANDS_PLOTBANDS(GBAS,N1,N2,N3,MAP,NKPT,NB,EB(:,:,ISPIN) &
      &                      ,NP,X1,XK1,X2,XK2,NQ,XQ)
         END IF
!
         CALL LINKEDLIST$SELECT(LL_CNTL,'..')
       ENDDO
!
!     ==========================================================================
!     ==  CLOSING                                                             ==
!     ==========================================================================
      CALL FILEHANDLER$REPORT(NFILO,'USED')
      WRITE(NFILO,FMT='(72("="))')
      WRITE(NFILO,FMT='(72("="),T20,"  PAW_BANDS TOOL FINISHED  ")')
      WRITE(NFILO,FMT='(72("="))')
                            CALL TRACE$PASS('AFTER CLOSING')
!
!     ==========================================================================
!     ==  CLOSE FILES                                                         ==
!     ==========================================================================
      IF(TREFINE) THEN
        PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        PRINT*,'CAUTION! THE EXPERIMENTAL OPTION TREFINE HAS BEEN USED'
        PRINT*,'THIS IMPLEMENTATION IS NOT READY FOR USE'
        PRINT*,'THE INTERPOLATION TOA FINE GRID DOES NOT WORK'
        PRINT*,'PROBABLY ALSO THE BAND CROSSINGS NEED TO BE ACCOUNTED FOR' 
        PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     END IF
      CALL FILEHANDLER$CLOSEALL
                            CALL TRACE$PASS('AFTER FILEHANDLER$CLOSEALL')
      CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      STOP
      END PROGRAM MAIN
!      
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: PDOSINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
!     ******************************************************************
      CALL LIB$NARGS(NARGS)
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE PDOS TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL LIB$GETARG(1,PDOSINNAME)
      ISVAR=INDEX(PDOSINNAME,-'.BCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=PDOSINNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('BCNTL',.FALSE.,PDOSINNAME)
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
!     ==================================================================
!     == SET STANDARD FILENAMES                                       ==
!     ==================================================================
!
!     ==  ERROR FILE ===================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,T,-'.BERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOLL FILE================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.BPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =============================================
      ID=+'BCNTL'
      CALL FILEHANDLER$SETFILE(ID,T,-'.BCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  STRUCTURE FILE   =============================================
      ID=+'PDOS'
      CALL FILEHANDLER$SETFILE(ID,T,-'.PDOS')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE BANDS_PLOTBANDS(GBAS,N1,N2,N3,MAP,NK,NB,EB &
      &                          ,NP,X1,XK1,X2,XK2,NQ,XQ)
!      *************************************************************************
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: GBAS(3,3)
       INTEGER(4),INTENT(IN) :: N1,N2,N3
       INTEGER(4),INTENT(IN) :: NK
       INTEGER(4),INTENT(IN) :: NB
       REAL(8)   ,INTENT(IN) :: EB(NB,NK)
       INTEGER(4),INTENT(IN) :: MAP(N1,N2,N3)
       INTEGER(4),INTENT(IN) :: NP
       INTEGER(4),INTENT(IN) :: NQ
       REAL(8)   ,INTENT(IN) :: XK1(3)  !initial k-point in relative coordinates
       REAL(8)   ,INTENT(IN) :: XK2(3)  !final k-point in relative coordinates
       REAL(8)   ,INTENT(IN) :: XQ(3)   !projection direction in relative coord.
       REAL(8)   ,INTENT(IN) :: X1
       REAL(8)   ,INTENT(IN) :: X2
       REAL(8)               :: D1,D2,D3
       REAL(8)               :: KI(3)
       REAL(8)               :: EBI(NB)
       INTEGER(4)            :: IP,IQ
       INTEGER(4)            :: NFIL=11
!      *************************************************************************
print*,'nb ',nb
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
!          == ki is the k-point in cartesian coordinates =======================
           KI=MATMUL(GBAS,XK1*D1+XK2*D2+XQ*D3)
           CALL BANDS_GETBAND(GBAS,N1,N2,N3,MAP,NK,NB,EB,KI,EBI)
           IF(NB.GT.100) THEN
             CALL ERROR$MSG('NUMBER OF BANDS EXCEEDS LIMIT OF 100')
             CALL ERROR$I4VAL('NB ',NB)
             CALL ERROR$STOP('BANDS_PLOTBANDS')
           END IF
           WRITE(NFIL,FMT='(100F9.5)')X1+D2*(X2-X1),EBI
         ENDDO
       ENDDO
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
       SUBROUTINE BANDS_GETBAND(GBAS,N1,N2,N3,MAP,NK,NB,EB,K1,EB1)
!      *************************************************************************
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: GBAS(3,3)
       INTEGER(4),INTENT(IN) :: N1,N2,N3
       INTEGER(4),INTENT(IN) :: MAP(N1,N2,N3)
       INTEGER(4),INTENT(IN) :: NK
       INTEGER(4),INTENT(IN) :: NB
       REAL(8)   ,INTENT(IN) :: EB(NB,NK)
       REAL(8)   ,INTENT(IN) :: K1(3)
       REAL(8)   ,INTENT(OUT):: EB1(NB)
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
       DO I=1,4
         EB1(:)=EB1(:)+WGHT(I)*EB(:,CORNER_MAP(I))
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
