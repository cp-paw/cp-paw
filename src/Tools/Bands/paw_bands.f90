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
      IMPLICIT NONE
      TYPE(LL_TYPE)             :: LL_CNTL
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NFIL
      INTEGER(4)                :: METHOD
      LOGICAL(4)                :: TCHK
!     ******************************************************************
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
!     ==  SELECT METHOD                                                       ==
!     ==========================================================================
      METHOD=1
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'METHOD',1)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,METHOD)
      ENDIF
      IF(METHOD.eq.0)THEN
        CALL BANDS_LINEAR_INTERPOLATION(LL_CNTL,NFIL,NFILO)
      ELSE IF(METHOD.eq.1)then
        CALL BANDS_DIAG(LL_CNTL,NFIL,NFILO)
      ELSE
        CALL ERROR$MSG('!METHOD!ID NOT IMPLEMENTED')
        CALL ERROR$I4VAL('ID',METHOD)
        CALL ERROR$STOP('MAIN')
      ENDIF
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
      CALL FILEHANDLER$CLOSEALL
                            CALL TRACE$PASS('AFTER FILEHANDLER$CLOSEALL')
      CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      END PROGRAM MAIN
!      
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: PDOSINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: IARGC
!     ******************************************************************
      IF(IARGC().LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE PDOS TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GETARG(1,PDOSINNAME)
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
      LOGICAL(4),PARAMETER :: T=.TRUE.
      LOGICAL(4),PARAMETER :: F=.FALSE.
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
      END SUBROUTINE STANDARDFILES_BANDS
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_LINEAR_INTERPOLATION(LL_CNTL,NFIL,NFILO)
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
      REAL(8)                   :: GBAS(3,3) ! REAL-SPACE LATTICE VECTORS
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
      REAL(8)                   :: X1,X2
      REAL(8)                   :: SVAR,svar1
      LOGICAL(4)                :: TCHK,TCHK1
      LOGICAL(4),PARAMETER      :: TREFINE=.FALSE.
!     **************************************************************************
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
            PROPSI(IB,IKPT,ISPIN)=abs(STATE%VEC(1,1,IB))
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
!           CALL BANDS_PLOTBANDS(GBAS,N1B,N2B,N3B,MAPBIG &
!      &                        ,NKPTBIG,NB,EBBIG(:,:,ISPIN) &
!      &                        ,NP,X1,XK1,X2,XK2,NQ,XQ)
         ELSE
           CALL BANDS_PLOTBANDS(GBAS,N1,N2,N3,MAP,NKPT,NB,EB(:,:,ISPIN),PROPSI(:,:,ISPIN) &
      &                      ,NP,X1,XK1,X2,XK2,NQ,XQ)
         END IF
!
         CALL LINKEDLIST$SELECT(LL_CNTL,'..')
       ENDDO
      IF(TREFINE) THEN
        PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        PRINT*,'CAUTION! THE EXPERIMENTAL OPTION TREFINE HAS BEEN USED'
        PRINT*,'THIS IMPLEMENTATION IS NOT READY FOR USE'
        PRINT*,'THE INTERPOLATION TOA FINE GRID DOES NOT WORK'
        PRINT*,'PROBABLY ALSO THE BAND CROSSINGS NEED TO BE ACCOUNTED FOR' 
        PRINT*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      END IF
      RETURN
      END SUBROUTINE BANDS_LINEAR_INTERPOLATION
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
       REAL(8)   ,INTENT(IN) :: XK1(3)  !initial k-point in relative coordinates
       REAL(8)   ,INTENT(IN) :: XK2(3)  !final k-point in relative coordinates
       REAL(8)   ,INTENT(IN) :: XQ(3)   !projection direction in relative coord.
       REAL(8)   ,INTENT(IN) :: X1
       REAL(8)   ,INTENT(IN) :: X2
       REAL(8)               :: D1,D2,D3
       REAL(8)               :: KI(3)
       REAL(8)               :: EBI(NB)
       REAL(8)               :: PROPSII(NB)
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
           CALL BANDS_GETBAND(GBAS,N1,N2,N3,MAP,NK,NB,EB,PROPSI,KI,EBI,PROPSII)
           IF(NB.GT.100) THEN
             CALL ERROR$MSG('NUMBER OF BANDS EXCEEDS LIMIT OF 100')
             CALL ERROR$I4VAL('NB ',NB)
             CALL ERROR$STOP('BANDS_PLOTBANDS')
           END IF
           WRITE(NFIL,FMT='(100F9.5)')sqrt(sum(KI**2)),X1+D2*(X2-X1),sqrt(sum(KI**2)),EBI,PROPSII
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_DIAG(LL_CNTL,NFIL,NFILO)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE BANDDATA_MODULE
      USE RADIAL_MODULE
      USE WAVES_MODULE, ONLY : GSET_TYPE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT):: LL_CNTL
      INTEGER(4)   ,INTENT(IN)  :: NFIL
      INTEGER(4)   ,INTENT(IN)  :: NFILO
      INTEGER(4)                :: NFILN
      LOGICAL(4)                :: TPRINT=.FALSE.
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: METHOD_DIAG !1=LAPACK,2=LANCZOS
      REAL(8)                   :: XK1(3)  !initial k-point in relative coordinates
      REAL(8)                   :: XK2(3)  !final k-point in relative coordinates
      INTEGER(4)                :: NK,NKDIAG
      INTEGER(4)                :: NB
      INTEGER(4)                :: ILINE,NLINE,NFILE
      INTEGER(4)                :: ISPIN
      REAL(8)                   :: X1,X2
      CHARACTER(512)            :: BANDDATAFILE,ID
      CHARACTER(512)            :: FILE

      INTEGER(4)                :: I,J,I1,J1,IDIM1,IDIM2
      REAL(8)                   :: XK(3)
      REAL(8)                   :: KVEC(3)
      REAL(8)                   :: T1,T2,T3
      REAL(8)                   :: GWEIGHT
      INTEGER(4)                :: NG,IG
      INTEGER(4)                :: IKDIAG
      REAL(8)                   :: CELLVOL
      TYPE(GSET_TYPE)           :: GSET
      REAL(8)   ,allocatable    :: GVEC(:,:)   !G
      REAL(8)   ,allocatable    :: GVECPK(:,:) !G+K
      REAL(8)   ,allocatable    :: G2(:)
      INTEGER(4)                :: GIDG_PROTO
      COMPLEX(8),allocatable    :: PSI1(:,:,:)
      COMPLEX(8),allocatable    :: HPSI1(:,:,:)
      COMPLEX(8),allocatable    :: EIGR(:)
      COMPLEX(8),allocatable    :: TI_H(:,:,:)
      COMPLEX(8),allocatable    :: TI_S(:,:,:)
      COMPLEX(8),allocatable    :: TI_HK(:,:)
      COMPLEX(8),allocatable    :: TI_SK(:,:)
      COMPLEX(8),allocatable    :: FOFG(:)
      COMPLEX(8) ,ALLOCATABLE   :: PRO(:,:,:)   !(NAT,NGG,LMNXX)
      REAL(8),allocatable       :: BAREPRO(:,:)
      REAL(8),allocatable       :: YLM(:,:)
      REAL(8),allocatable       :: YLM_(:)
      INTEGER(4)                :: IAT,IBPRO,LN,IND,IPRO,ISP,K,L,M
      INTEGER(4)                :: LMNX_,LNX_,LMN1,LMN2
      INTEGER(4),allocatable    :: LOX_(:)
      REAL(8)                   :: SVAR
      REAL(8),allocatable       :: E(:)
      REAL(8),allocatable       :: FATBANDVAL(:,:,:)
      REAL(8),allocatable       :: EIGVAL(:,:)
      REAL(8),allocatable       :: XKVAL(:,:)
      REAL(8),allocatable       :: KVECVAL(:,:)
      COMPLEX(8),allocatable    :: U(:,:)
      INTEGER(4)                :: IB
      COMPLEX(8)                :: CSVAR
      CHARACTER(64)             :: ATOM
      INTEGER(4)                :: NFILBAND
      INTEGER(4)                :: NFILFATBAND
      INTEGER(4)                :: NFATBAND,IFATBAND,ILINESPERBAND
      INTEGER(4),ALLOCATABLE    :: FATBANDIPRO(:),FATBANDIAT(:),LINESPERBAND(:)
      REAL(8),ALLOCATABLE       :: FATBANDMAXWIDTH(:)
      CHARACTER(512),ALLOCATABLE:: FATBANDFILE(:)
      REAL(8)                   :: FATBANDMAX,SVAR1,SVAR2
!     **************************************************************************
                            CALL TRACE$PUSH('BANDS_DIAG')
!
!     ==========================================================================
!     ==  SET BANDDATAFILE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'INPUTFILE',NFILE)
      IF(NFILE.eq.1)THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'INPUTFILE',NFILE)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(TCHK)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,BANDDATAFILE)
          ID=+'BANDDATAIN'
          CALL FILEHANDLER$SETFILE(ID,.false.,-BANDDATAFILE)
          CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
        ELSE
          CALL ERROR$MSG('!BCNTL!INPUTFILE!NAME NOT FOUND')
          CALL ERROR$STOP('BANDS_DIAG')
        ENDIF
      ELSE IF (NFILE.gt.1)then
        CALL ERROR$MSG('MULTIPLE INPUT FILES GIVEN (!BCNTL!INPUTFILE)')
        CALL ERROR$STOP('BANDS_DIAG')
      ELSE
        ID=+'BANDDATAIN'
        CALL FILEHANDLER$SETFILE(ID,.true.,-'.BANDDATA')
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
      IF(NDIM.eq.2)THEN
        CALL ERROR$MSG('ONLY NDIM=1 AND NSPIN=1 OR 2 IMPLEMENTED AND TESTED')
        CALL ERROR$STOP('BANDS_DIAG')
      ENDIF
!
!     ==========================================================================
!     ==  READ BCNTL                                                          ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'METHOD_DIAG',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'METHOD_DIAG',1,METHOD_DIAG)
      ELSE
        METHOD_DIAG=1
      ENDIF
      IF(METHOD_DIAG.ne.1)THEN
        CALL ERROR$MSG('ONLY LAPACK IMPLEMENTED (METHOD_DIAG=1)')
        CALL ERROR$STOP('BANDS_DIAG')
      ENDIF
!
!     =========================================================================
!     ==  CONSTRUCT K-INDEPENDENT PART OF HAMILTONIAN                        ==
!     =========================================================================
!     == GENERATE G-GRID ==================================================
      KVEC(1)=0.0D0
      KVEC(2)=0.0D0
      KVEC(3)=0.0D0
!
!     ==========================================================================
!     == INITIALIZE MPE ROUTINE FOR PARALLEL PROCESSING                       ==
!     ==========================================================================
      CALL MPE$INIT

      CALL PLANEWAVE$INITIALIZE('WAVE GAMMA','~',RBAS,KVEC,.false.,EPW &
     &                               ,0,NR1,NR2,NR3)
      CALL PLANEWAVE$SELECT('WAVE GAMMA')      

      CALL PLANEWAVE_COUNTG(GBAS,KVEC,EPW,NG)
      GSET%NGL=NG
      IF(TPRINT)PRINT*,"EPW NG",EPW,NG
      
      ALLOCATE(GVEC(3,NG))
      ALLOCATE(GVECPK(3,NG))
      ALLOCATE(G2(NG))
      ALLOCATE(EIGR(NG))
      ALLOCATE(PRO(NAT,NG,LMNXX))
      CALL PLANEWAVE$GETR8A('GVEC',3*NG,GVEC)
      CALL PLANEWAVE$GETR8A('G2',NG,G2)
      CALL GBASS(RBAS,GBAS,GWEIGHT)
      CELLVOL=GWEIGHT

!     == GENERATE G-GRID FOR PROJECTORS ========================================
      CALL RADIAL$NEW(GRIDTYPE(TYPE_PROTO),GIDG_PROTO)
      !DEX_PROTO=LOG(GMAX_PROTO/G1_PROTO)/REAL(NG_PROTO-1,KIND=8)
      CALL RADIAL$SETI4(GIDG_PROTO,'NR',NG_PROTO)
      CALL RADIAL$SETR8(GIDG_PROTO,'R1',G1_PROTO)
      CALL RADIAL$SETR8(GIDG_PROTO,'DEX',DEX_PROTO)

!     == ALLOCATE K-INDEPENDET ARRAYS ==========================================
      allocate(TI_H(NG*NDIM,NG*NDIM,NSPIN))
      allocate(TI_S(NG*NDIM,NG*NDIM,NSPIN))
      TI_H=0.0D0
      TI_S=0.0D0
       
!     == COMPUTE K_INDEPENDET ARRAYS ===========================================
      allocate(HPSI1(NG,NDIM,1))
      allocate(PSI1(NG,NDIM,1))
      !EVALUATE MATRIX-ELEMENTS OF PS-POTENTIAL
      !FIXME: introduce G-G', this part could use some work
      DO ISPIN=1,NSPIN
        DO I=1,NG
          DO IDIM1=1,NDIM
            I1=I*(NDIM)+(IDIM1-1)
            PSI1(:,:,:)=0.0D0
            PSI1(I,IDIM1,1)=1.0D0
            HPSI1(:,:,:)=0.0D0
            call WAVES_VPSI(GSET,NG,NDIM,1,NRL,PSI1,VOFR(:,ISPIN),HPSI1)
         
            DO IDIM2=1,NDIM
              TI_H(I1,1+(IDIM2-1)*NG:IDIM2*NG,ISPIN)=HPSI1(:,IDIM2,1)
            ENDDO
            TI_S(I1,I1,ISPIN)=1.0D0
          ENDDO
        ENDDO
      ENDDO
      !REMOVE KINETIC ENERGY FROM TI_H
      DO I=1,NG
        DO IDIM1=1,NDIM
          I1=I*(NDIM)+(IDIM1-1)
          !EKIN
          TI_H(I1,I1,:)=TI_H(I1,I1,:)-0.5D0*sum(GVEC(:,I)**2)
        ENDDO
      ENDDO
      allocate(TI_HK(NG*NDIM,NG*NDIM))
      allocate(TI_SK(NG*NDIM,NG*NDIM))
      ALLOCATE(U(NG*NDIM,NG*NDIM))
      ALLOCATE(E(NG*NDIM))
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
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'LINE',NLINE)
      
      X2=0.D0
      DO ILINE=1,NLINE
        CALL LINKEDLIST$SELECT(LL_CNTL,'~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
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
        !NK is number of points for output
        !if NKDIAG<NK: do 1D-fft interpoaltion
        !if NKDIAG=NK: use results from diagonalisation directly
        !if NKDIAG>NK: error
        NB=10
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NB',0,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NB',1,NB)
        NK=10
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NK',0,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NK',1,NK)
        !NKDIAG is number of k-points in diagonalisation
        NKDIAG=NK
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NKDIAG',0,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NKDIAG',1,NKDIAG)
        IF(NKDIAG.gt.NK)THEN
          CALL ERROR$MSG('NKDIAG>NK')
          CALL ERROR$I4VAL('NK',NK)
          CALL ERROR$I4VAL('NKDIAG',NKDIAG)
          CALL ERROR$STOP('BANDS_DIAG')
        ENDIF
        !FIXME: IMPLEMENT 1D FFT INTERPOATION
        IF(NKDIAG.lt.NK)THEN
          CALL ERROR$MSG('1D FFT INTERPOLATION NOT YET IMPLEMENTED!!')
          CALL ERROR$I4VAL('NK',NK)
          CALL ERROR$I4VAL('NKDIAG',NKDIAG)
          CALL ERROR$STOP('BANDS_DIAG')
        ENDIF
        IF(TPRINT)PRINT*,XK1,XK2,NB,NK,NKDIAG,X1,X2,FILE
!       == SPECIFY SPIN DIRECTION =============================================
        ISPIN=1
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SPIN',0,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'SPIN',1,ISPIN)
        IF(ISPIN.GT.NSPIN) THEN
          CALL ERROR$MSG('ISPIN EXCEEDS RANGE')
          CALL ERROR$STOP('MAIN')
        END IF
!       == SPECIFY FATBAND ===================================================
        CALL LINKEDLIST$SELECT(LL_CNTL,'~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
        CALL LINKEDLIST$SELECT(LL_CNTL,'LINE',ILINE)
        CALL LINKEDLIST$NLISTS(LL_CNTL,'FATBAND',NFATBAND)

        IF(NFATBAND.GE.1)then
          ALLOCATE(FATBANDIPRO(NFATBAND))
          ALLOCATE(FATBANDIAT(NFATBAND))
          ALLOCATE(LINESPERBAND(NFATBAND))
          ALLOCATE(FATBANDMAXWIDTH(NFATBAND))
          ALLOCATE(FATBANDFILE(NFATBAND))

          DO IFATBAND=1,NFATBAND
            CALL LINKEDLIST$SELECT(LL_CNTL,'~')
            CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
            CALL LINKEDLIST$SELECT(LL_CNTL,'LINE',ILINE)
            CALL LINKEDLIST$SELECT(LL_CNTL,'FATBAND',IFATBAND)
            !FILE
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',0,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FATBANDFILE(IFATBAND))
            ELSE
              CALL ERROR$I4VAL('FATBAND',IFATBAND)
              CALL ERROR$MSG('FILE not given')
              CALL ERROR$STOP('MAIN')
            END IF
            !ATOM
            CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',0,TCHK)
            IF(TCHK)THEN
              CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,ATOM)
            ELSE
              CALL ERROR$I4VAL('FATBAND',IFATBAND)
              CALL ERROR$MSG('ATOM not given')
              CALL ERROR$STOP('MAIN')
            ENDIF
            !FIND NUMBER OF ATOM IN ATOMID
            FATBANDIAT(IFATBAND)=0
            DO IAT=1,NAT
              IF(ATOMID(IAT).eq.ATOM)THEN
                FATBANDIAT(IFATBAND)=IAT
                exit
              ENDIF
            ENDDO
            IF(FATBANDIAT(IFATBAND).eq.0)THEN
              CALL ERROR$MSG('ATOM not found')
              CALL ERROR$CHVAL('ATOM',ATOM)
              CALL ERROR$I4VAL('FATBAND',IFATBAND)
              CALL ERROR$STOP('MAIN')
            ENDIF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'IPRO',0,TCHK)
            IF(TCHK)THEN
              CALL LINKEDLIST$GET(LL_CNTL,'IPRO',1,FATBANDIPRO(IFATBAND))
            ELSE
              CALL ERROR$MSG('IPRO not given')
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
        ALLOCATE(FATBANDVAL(NFATBAND,NKDIAG,NB))
        ALLOCATE(KVECVAL(3,NKDIAG))
        ALLOCATE(XKVAL(3,NKDIAG))

        IF(METHOD_DIAG.eq.1)THEN  
!         == ITERATE K-POINTS =================================================
          DO IKDIAG=0,NKDIAG-1
            XK=XK1+(XK2-XK1)*REAL(IKDIAG,KIND=8)/REAL(MAX(NKDIAG-1,1),KIND=8)
            KVEC=MATMUL(GBAS,XK)
            KVECVAL(:,IKDIAG+1)=KVEC
            XKVAL(:,IKDIAG+1)=XK
            
            TI_HK=TI_H(:,:,ISPIN)
            TI_SK=TI_S(:,:,ISPIN)

            !COMPUTE G+K and (G+K)^2
            DO I=1,NG
              GVECPK(:,I)=GVEC(:,I)+KVEC(:)
              G2(I)=sum(GVECPK(:,I)**2)
            ENDDO

            IF(.NOT.ALLOCATED(BAREPRO)) ALLOCATE(BAREPRO(NG,NBAREPRO))
            IND=0
            DO ISP=1,NSP
              DO LN=1,LNX(ISP)
                IND=IND+1
                !BANDS_GETFOFG IS A MODIFIED VERSION OF SETUP$GETFOFG
                CALL BANDS_GETFOFG(GIDG_PROTO,LN,ISP,NG,G2,CELLVOL,BAREPRO(:,IND))
              ENDDO
            ENDDO
            
            IF(.NOT.ALLOCATED(YLM)) ALLOCATE(YLM(NG,LMX))
            ALLOCATE(YLM_(LMX))
            DO I=1,NG
              CALL GETYLM(LMX,GVECPK(1,I),YLM_)
              YLM(I,:)=YLM_(:) 
            ENDDO        
            DEALLOCATE(YLM_)
            
            ALLOCATE(LOX_(LNXX))    
            IPRO=1
            PRO(:,:,:)=0.0D0
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              LNX_=LNX(ISP)
              LMNX_=LMNX(ISP)
              LOX_=LOX(ISP,:)
              IBPRO=1+SUM(LNX(1:ISP-1))
              DO I=1,NG
                SVAR=SUM(GVECPK(:,I)*R(:,IAT))
                EIGR(I)=CMPLX(cos(SVAR),-sin(SVAR))
              ENDDO
              CALL WAVES_EXPANDPRO(LNX_,LOX_,LMNX_,NG,GVECPK &
           &         ,BAREPRO(:,IBPRO:IBPRO+LNX_-1),LMX,YLM,EIGR,PRO(IAT,:,:))
            ENDDO
            deallocate(LOX_)

            !FIXME: DH IS NOT HERMITIAN, SEE paw_setups.f90 setup_MAKEPARTIALWAVES
            !here we make it symmetric, so that ZHEGV can be used instead of
            !ZGGEV. 
            DO IAT=1,NAT
              DH(:,:,ISPIN,IAT)=0.5D0*(DH(:,:,ISPIN,IAT)+&
        &                   TRANSPOSE(CONJG(DH(:,:,ISPIN,IAT))))  
            ENDDO

            !add (G+k)^2/2 and augmentation to TI_HK
            DO I=1,NG
              DO IDIM1=1,NDIM
                DO J=1,NG
                  DO IDIM2=1,NDIM
                    I1=I*(NDIM)+(IDIM1-1)
                    J1=J*(NDIM)+(IDIM2-1)
                    !EKIN
                    IF(I1.eq.J1)TI_HK(I1,J1)=TI_HK(I1,J1)+0.5D0*G2(I)
    
                    !Augmentation LMNXX,LMNXX,NDIMD,NAT
                    DO K=1,LMNXX
                      DO L=1,LMNXX
                        DO M=1,NAT
                          TI_HK(I1,J1)=TI_HK(I1,J1)+&
       &                    DH(K,L,ISPIN,M)*CONJG(PRO(M,I,K))*PRO(M,J,L)*GWEIGHT
                          TI_SK(I1,J1)=TI_SK(I1,J1)+&
       &                    DO(K,L,ISPIN,M)*CONJG(PRO(M,I,K))*PRO(M,J,L)*GWEIGHT
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            !SOLVE EIGENVALUE PROBLEM WITH LAPACK ROUTINES
            CALL LIB$GENERALEIGENVALUEC8(NG*NDIM,TI_HK,TI_SK,E,U)
            IF(TPRINT)PRINT*,'DIAGBANDS_EIG',sqrt(sum(KVEC**2)),E(1:NB)*27.21139D0
            EIGVAL(1:NB,IKDIAG+1)=E(1:NB)*27.21139D0
!
!           ====================================================================
!           ==  FAT BANDS                                                     ==
!           ====================================================================
            IF(NFATBAND.GE.1)THEN
              DO IFATBAND=1,NFATBAND
                IF(TPRINT)PRINT*,"FATBAND IFATBAND",IFATBAND
                IF(TPRINT)PRINT*,"FATBAND IAT",FATBANDIAT(IFATBAND)
                IF(TPRINT)PRINT*,"FATBAND IPRO",FATBANDIPRO(IFATBAND)
                DO IB=1,NB
                  CALL PLANEWAVE$SCALARPRODUCT(' ',NG,NDIM,1,PRO(FATBANDIAT(IFATBAND),:,FATBANDIPRO(IFATBAND)),1,U(:,IB),CSVAR)
                  FATBANDVAL(IFATBAND,IKDIAG+1,IB)=ABS(CSVAR)**2
                ENDDO
              ENDDO
            ELSE
              FATBANDVAL(:,:,:)=0.0D0 
            ENDIF
          ENDDO
        ENDIF
!
!       =========================================================================
!       ==  WRITE BANDS                                                        ==
!       =========================================================================
        CALL FILEHANDLER$UNIT('BANDS',NFILBAND)
        IF(NB.GT.100) THEN
          CALL ERROR$MSG('NUMBER OF BANDS EXCEEDS LIMIT OF 100')
          CALL ERROR$I4VAL('NB ',NB)
          CALL ERROR$STOP('BANDS_PLOTBANDS')
        END IF

        !ITERATE K-POINTS
        DO IKDIAG=1,NKDIAG
          !WRITE EIGENVALUES
          WRITE(NFILBAND,FMT='(103F9.5)')sqrt(sum(XKVAL(:,IKDIAG)**2)),sqrt(sum(KVECVAL(:,IKDIAG)**2)),&
    &           REAL(IKDIAG-1,KIND=8)/REAL(MAX(NKDIAG-1,1),KIND=8),EIGVAL(:,IKDIAG)
        ENDDO
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
            WRITE(NFILFATBAND,FMT='(103F9.5)',advance="no")sqrt(sum(XKVAL(:,IKDIAG)**2)),sqrt(sum(KVECVAL(:,IKDIAG)**2)),&
    &           REAL(IKDIAG-1,KIND=8)/REAL(MAX(NKDIAG-1,1),KIND=8)

            FATBANDMAX=MAXVAL(FATBANDVAL(IFATBAND,:,:))
            IF(TPRINT)PRINT*,"FATBAND FATBANDMAX",FATBANDMAX
            DO ILINESPERBAND=0,LINESPERBAND(IFATBAND)-1
              SVAR1=2.0D0/REAL(LINESPERBAND(IFATBAND)-1,KIND=8)*REAL(ILINESPERBAND,KIND=8)-1.0D0
              SVAR2=0.5D0*SVAR1*FATBANDMAXWIDTH(IFATBAND)/FATBANDMAX
              WRITE(NFILFATBAND,FMT='(100F9.5)',advance="no")EIGVAL(:,IKDIAG)+SVAR2*FATBANDVAL(IFATBAND,IKDIAG,:) 
            ENDDO
            WRITE(NFILFATBAND,FMT='(a)')" "
          ENDDO
        ENDDO
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
      END SUBROUTINE BANDS_DIAG
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
      END
