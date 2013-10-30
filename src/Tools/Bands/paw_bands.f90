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
      USE OMP_LIB
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
      ELSE IF(METHOD.eq.2)THEN
        CALL BANDS_DOS(LL_CNTL,NFIL,NFILO)
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

      ID=+'PDOSOUT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.PDOSOUT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
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
      USE BANDDATA_MODULE
      USE RADIAL_MODULE
      USE WAVES_MODULE, ONLY : GSET_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)              :: NG
      REAL(8),allocatable,INTENT(INOUT)   :: GVEC(:,:)   !G
      REAL(8),allocatable,INTENT(INOUT)   :: G2(:)
      REAL(8),INTENT(OUT)                 :: GWEIGHT
      INTEGER(4),INTENT(OUT)              :: GIDG_PROTO
      COMPLEX(8),allocatable,INTENT(INOUT):: TI_H(:,:,:)
      COMPLEX(8),allocatable,INTENT(INOUT):: TI_S(:,:,:)
      REAL(8)                             :: KVEC(3)
      TYPE(GSET_TYPE)                     :: GSET
      COMPLEX(8),allocatable              :: PSI1(:,:,:)
      COMPLEX(8),allocatable              :: HPSI1(:,:,:)
      INTEGER(4)                          :: ISPIN,I,IDIM1,IDIM2,I1
      LOGICAL(4)                          :: TPRINT=.true.
!
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
      ALLOCATE(G2(NG))
      CALL PLANEWAVE$GETR8A('GVEC',3*NG,GVEC)
      CALL PLANEWAVE$GETR8A('G2',NG,G2)
      CALL GBASS(RBAS,GBAS,GWEIGHT)

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
                            CALL TRACE$PUSH('COMPUTE K_INDEPENDET ARRAYS')
      CALL TIMING$CLOCKON('PS_POTENTIAL')
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
      deallocate(HPSI1)
      deallocate(PSI1)
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
      SUBROUTINE BANDS_KPOINT(NG,NB,ISPIN,METHOD,GIDG_PROTO,TPROJ,KVEC,GVEC,TI_H,TI_S,E,PROJ)
!     ******************************************************************
!     **                                                              **
!     ** COMPUTES EIGENVALUES[H] AND PROJECTIONS FOR A GIVEN KVEC     **
!     **                                                              **
!     ******************************************************************
      USE BANDDATA_MODULE
      use omp_lib
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: NG
      INTEGER(4)  ,INTENT(IN)  :: NB
      INTEGER(4)  ,INTENT(IN)  :: ISPIN
      INTEGER(4)  ,INTENT(IN)  :: METHOD
      INTEGER(4)  ,INTENT(IN)  :: GIDG_PROTO
      LOGICAL(4)  ,INTENT(IN)  :: TPROJ
      REAL(8)     ,INTENT(IN)  :: KVEC(3)
      REAL(8)     ,INTENT(IN)  :: GVEC(3,NG)
      COMPLEX(8)  ,INTENT(IN)  :: TI_H(NG*NDIM,NG*NDIM,NSPIN)
      COMPLEX(8)  ,INTENT(IN)  :: TI_S(NG*NDIM,NG*NDIM,NSPIN)
      REAL(8)     ,INTENT(OUT) :: E(NG*NDIM)
      COMPLEX(8)  ,INTENT(OUT) :: PROJ(NAT,NB,LMNXX)
      COMPLEX(8),allocatable   :: TI_HK(:,:)
      COMPLEX(8),allocatable   :: TI_SK(:,:)
      COMPLEX(8)               :: U(NG*NDIM,NG*NDIM)
      REAL(8)                  :: GVECPK(3,NG)
      REAL(8)                  :: G2(NG)
      COMPLEX(8)               :: EIGR(NG)
      COMPLEX(8)               :: FOFG(NG)
      COMPLEX(8)               :: PRO(NAT,NG,LMNXX)
      REAL(8)                  :: BAREPRO(NG,NBAREPRO)
      REAL(8)                  :: YLM(NG,LMX)
      REAL(8)                  :: YLM_(LMX)
      INTEGER(4)               :: LOX_(LNXX)
      INTEGER(4)               :: ISP,IND,LN,IAT,IDIM1,IDIM2,IBPRO
      INTEGER(4)               :: LMN1,LMN2,IB
      INTEGER(4)               :: LNX_,LMNX_,IPRO
      INTEGER(4)               :: I,I1,J,K,L,M,N
      REAL(8)                  :: CELLVOL,GWEIGHT
      REAL(8)                  :: SVAR
      COMPLEX(8)               :: CSVAR,CSVAR_H,CSVAR_S
      LOGICAL(4)               :: TTIMING
!     ******************************************************************
      if(omp_get_num_threads().eq.1)then
        TTIMING=.true.
      else
        TTIMING=.false.
      endif
      IF(TTIMING)CALL TRACE$PUSH('COMPUTE K_DEPENDENT ARRAY')
      allocate(TI_HK(NG*NDIM,NG*NDIM))
      allocate(TI_SK(NG*NDIM,NG*NDIM))
      TI_HK(:,:)=TI_H(:,:,ISPIN)
      TI_SK(:,:)=TI_S(:,:,ISPIN)
      
      CALL GBASS(RBAS,GBAS,GWEIGHT)
      CELLVOL=GWEIGHT

      !COMPUTE G+K and (G+K)^2
      IF(TTIMING) CALL TIMING$CLOCKON('G+K,G2')
      DO I=1,NG
        GVECPK(:,I)=GVEC(:,I)+KVEC(:)
        G2(I)=sum(GVECPK(:,I)**2)
      ENDDO
      IF(TTIMING)CALL TIMING$CLOCKOFF('G+K,G2')

      IF(TTIMING)CALL TIMING$CLOCKON('PROJECTORS')
      IND=0
      DO ISP=1,NSP
        DO LN=1,LNX(ISP)
          IND=IND+1
          !BANDS_GETFOFG IS A MODIFIED VERSION OF SETUP$GETFOFG
          CALL BANDS_GETFOFG(GIDG_PROTO,LN,ISP,NG,G2,CELLVOL,BAREPRO(:,IND))
        ENDDO
      ENDDO
      
      DO I=1,NG
        CALL GETYLM(LMX,GVECPK(1,I),YLM_)
        YLM(I,:)=YLM_(:) 
      ENDDO        
      
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
      IF(TTIMING)CALL TIMING$CLOCKOFF('PROJECTORS')

!      !NOTE: DTKIN/DH FROM SETUP_MAKEPARTIALWAVES NOW 
!      !(AFTER COMMIT 85f857388191c58f050acf5be749412fcfe54d12) PRODUCE A HERMITIAN HAMILTONIAN
!      DO IAT=1,NAT
!        DH(:,:,ISPIN,IAT)=0.5D0*(DH(:,:,ISPIN,IAT)+&
!  &                    TRANSPOSE(CONJG(DH(:,:,ISPIN,IAT))))  
!      ENDDO

      !add (G+k)^2/2 and augmentation to TI_HK
      DO I=1,NG
        DO IDIM1=1,NDIM
          I1=I*(NDIM)+(IDIM1-1)
          TI_HK(I1,I1)=TI_HK(I1,I1)+0.5D0*G2(I)
        ENDDO
      ENDDO
      IF(TTIMING)CALL TIMING$CLOCKON('AUGMENTATION')
      PRO(:,:,:)=sqrt(GWEIGHT)*PRO(:,:,:)
      !FIXME: OPTIMIZE THIS BLOCK!!!

      DO I=1,NG
        DO J=1,I
          CSVAR_H=0.0D0
          CSVAR_S=0.0D0
          DO M=1,NAT
            ISP=ISPECIES(M)
            LMNX_=LMNX(ISP)
            DO K=1,LMNX_
              DO L=1,LMNX_
              !FIXME: this will not work for non-collinear calc.
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
     
!      !check if TI_HK is hermitian
!      DO I=1,NG
!        DO J=I,NG
!          IF(abs(TI_HK(I,J)-CONJG(TI_HK(J,I))).gt.1d-10)THEN
!            PRINT*,"H",I,J
!            STOP
!          ENDIF
!          IF(abs(TI_SK(I,J)-CONJG(TI_SK(J,I))).gt.1d-10)THEN
!            PRINT*,"S",I,J
!            STOP
!          ENDIF
!        enddo
!      ENDDO

      !complete TI_HK and TI_SK
      DO I=1,NG
        DO J=I+1,NG
          TI_HK(I,J)=CONJG(TI_HK(J,I))
          TI_SK(I,J)=CONJG(TI_SK(J,I))
        ENDDO
      ENDDO

      IF(TTIMING)CALL TIMING$CLOCKOFF('AUGMENTATION')
      IF(TTIMING)CALL TRACE$POP

      !SOLVE GENERALIZED EIGENVALUE PROBLEM
      IF(TTIMING)CALL TIMING$CLOCKON('DIAG')
      IF(METHOD.eq.1)then
        IF(TTIMING)CALL TRACE$PUSH('LAPACK_ZHEGVD')
        IF(TPROJ)THEN
          U=TI_HK
          CALL LAPACK_ZHEGVD(NG,'V',U,TI_SK,E)            
        ELSE
          CALL LAPACK_ZHEGVD(NG,'N',TI_HK,TI_SK,E)            
        ENDIF
        IF(TTIMING)CALL TRACE$POP
      ELSE IF(METHOD.eq.2)then
        IF(TTIMING)CALL TRACE$PUSH('LIB$GENERALEIGENVALUEC8')
        CALL LIB$GENERALEIGENVALUEC8(NG*NDIM,TI_HK,TI_SK,E,U)
        IF(TTIMING)CALL TRACE$POP
      ELSE 
        CALL ERROR$MSG('METHOD_DIAG NOT IMPLEMENTED')
        CALL ERROR$I4VAL('METHOD_DIAG',METHOD)
        CALL ERROR$STOP('BANDS_KPOINT')
      ENDiF

      IF(TTIMING)CALL TIMING$CLOCKOFF('DIAG')
!$omp critical
      IF(TPROJ)THEN
        IF(TTIMING)CALL TIMING$CLOCKON('PROJECTIONS')
        DO IAT=1,NAT
          DO IB=1,NB
            DO LN=1,LMNXX
              CALL PLANEWAVE$SCALARPRODUCT(' ',NG,NDIM,1,PRO(IAT,:,LN),1,U(:,IB),CSVAR)
              PROJ(IAT,IB,LN)=CSVAR
            ENDDO
          ENDDO
        ENDDO
        IF(TTIMING)CALL TIMING$CLOCKOFF('PROJECTIONS')
      ENDIF
!$omp end critical
      DEALLOCATE(TI_HK)
      DEALLOCATE(TI_SK)
      RETURN
      end subroutine
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LAPACK_ZHEGVD(N,JOBZ,H,S,E)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: N
      Character(1),INTENT(IN)  :: JOBZ
      COMPLEX(8),INTENT(INOUT) :: H(N,N)
      COMPLEX(8),INTENT(INOUT) :: S(N,N)
      REAL(8),INTENT(OUT)      :: E(N)
      INTEGER(4)               :: LWORK,LRWORK,LIWORK,info
      COMPLEX(8),ALLOCATABLE   :: WORK(:)
      REAL(8)   ,ALLOCATABLE   :: RWORK(:)
      INTEGER(4),ALLOCATABLE   :: IWORK(:)
!     ******************************************************************

      LWORK=-1
      LRWORK=-1
      LIWORK=-1
      Allocate(WORK(1))!complex(8)
      Allocate(RWORK(1))!real(8)
      Allocate(IWORK(1))!integer(4)
      
      CALL ZHEGVD(1,JOBZ,'U',N,H,N,S,N,E,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
      LWORK=WORK(1)
      LRWORK=RWORK(1)
      LIWORK=IWORK(1)
      deallocate(WORK)
      deallocate(RWORK)
      deallocate(IWORK)
      allocate(WORK(LWORK)) 
      allocate(RWORK(LRWORK)) 
      allocate(IWORK(LIWORK)) 
      CALL ZHEGVD(1,JOBZ,'U',N,H,N,S,N,E,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
      deallocate(WORK)
      deallocate(RWORK)
      deallocate(IWORK)

      IF(INFO.ne.0)THEN
        CALL ERROR$MSG('ZHEGVD FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LAPACK_ZHEGVD')
      ENDIF
      return
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LAPACK_ZGGEV(N,JOBZ,H,S,E)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: N
      Character(1),INTENT(IN)  :: JOBZ
      COMPLEX(8),INTENT(INOUT) :: H(N,N)
      COMPLEX(8),INTENT(INOUT) :: S(N,N)
      REAL(8),INTENT(OUT)      :: E(N)
      INTEGER(4)               :: LWORK,LRWORK,LIWORK,info
      COMPLEX(8),ALLOCATABLE   :: WORK(:)
      REAL(8)   ,ALLOCATABLE   :: RWORK(:)
      Character(1)             :: JOBVL,JOBVR
      COMPLEX(8)               :: ALPHA(N),BETA(N)
      COMPLEX(8)               :: VL(N,N),VR(N,N)
!     ******************************************************************

      LWORK=-1
      Allocate(WORK(1))!complex(8)
      Allocate(RWORK(8*N))!complex(8)
      JOBVL='N'
      JOBVR=JOBZ
!      ( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
!         22 *                         VL, LDVL, VR, LDVR, WORK, LWORK, RWORK,
!         INFO )
      CALL ZGGEV(JOBVL,JOBVR,N,H,N,S,N,ALPHA,BETA,VL,N,VR,N,WORK,LWORK,RWORK,INFO)
      LWORK=WORK(1)
      deallocate(WORK)
      allocate(WORK(LWORK)) 
      CALL ZGGEV(JOBVL,JOBVR,N,H,N,S,N,ALPHA,BETA,VL,N,VR,N,WORK,LWORK,RWORK,INFO)
      deallocate(WORK)

!      IF(INFO.ne.0)THEN
!        CALL ERROR$MSG('ZGGEV FAILED')
!        CALL ERROR$I4VAL('INFO',INFO)
!        CALL ERROR$STOP('LAPACK_ZGGEV')
!      ENDIF
      E=REAL(ALPHA/BETA)
      return
      END SUBROUTINE

END MODULE
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
      USE KPOINTDIAG_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT):: LL_CNTL
      INTEGER(4)   ,INTENT(IN)  :: NFIL
      INTEGER(4)   ,INTENT(IN)  :: NFILO
      INTEGER(4)                :: NFILN
      LOGICAL(4)                :: TPRINT=.FALSE.
      LOGICAL(4)                :: TCHK,TCHK1,TCHK2
      INTEGER(4)                :: METHOD_DIAG !1=LAPACK,2=LANCZOS
      REAL(8)                   :: XK1(3)  !initial k-point in relative coordinates
      REAL(8)                   :: XK2(3)  !final k-point in relative coordinates
      REAL(8)                   :: KVEC1(3)  !initial k-point in absolute coordinates
      REAL(8)                   :: KVEC2(3)  !final k-point in abolute coordinates
      REAL(8)                   :: KVECSCALE
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
      REAL(8)   ,allocatable    :: GVEC(:,:)   !G
      REAL(8)   ,allocatable    :: G2(:)
      INTEGER(4)                :: GIDG_PROTO
      COMPLEX(8),allocatable    :: TI_H(:,:,:)
      COMPLEX(8),allocatable    :: TI_S(:,:,:)
      COMPLEX(8),allocatable    :: PROJ(:,:,:)
      INTEGER(4)                :: IAT,IBPRO,LN,IND,IPRO,ISP,K,L,M
      INTEGER(4)                :: LMNX_,LNX_,LMN1,LMN2
      REAL(8),allocatable       :: E(:)
      REAL(8),allocatable       :: FATBANDVAL(:,:,:)
      REAL(8),allocatable       :: EIGVAL(:,:)
      REAL(8),allocatable       :: XKVAL(:,:)
      REAL(8),allocatable       :: KVECVAL(:,:)
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
      CALL TIMING$START
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
      CALL LINKEDLIST$SELECT(LL_CNTL,'METHOD')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'METHOD_DIAG',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'METHOD_DIAG',1,METHOD_DIAG)
      ELSE
        METHOD_DIAG=1
      ENDIF
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EPWPSI',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'EPWPSI',1,EPW)
        WRITE(NFILO,*)'WARNING: EPWPSI set to',EPW,' Ry!'
        !input is in Hartree
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
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'BCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'LINE',NLINE)
      
      DO ILINE=1,NLINE
        WRITE(NFILO,FMT='(72("="))')
        WRITE(NFILO,*)'LINE-BLOCK:',ILINE,' of ',NLINE
        WRITE(NFILO,FMT='(72("="))')
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
          WRITE(NFILO,*)'OUTPUT FILE: ',trim(FILE)
        ELSE
          CALL ERROR$MSG('NO OUTPUT FILE GIVEN')
          CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
          CALL ERROR$STOP('BANDS_DIAG')
        END IF
!
        !INPUT OF K-VECTORS
        CALL LIB$INVERTR8(3,GBAS,GBASINV) 

        CALL LINKEDLIST$EXISTD(LL_CNTL,'XK1',0,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'KVEC1',0,TCHK2)
        IF(TCHK1)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'XK1',1,XK1)
          WRITE(NFILO,*)'FIRST K-POINT: GIVEN IN RELATIVE COORDINATES'
          KVEC1=MATMUL(GBAS,XK1)
        ELSE IF (TCHK2)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'KVEC1',1,KVEC1)
          WRITE(NFILO,*)'FIRST K-POINT: GIVEN IN ABSOLUTE COORDINATES'

          CALL LINKEDLIST$EXISTD(LL_CNTL,'KVECSCALE',0,TCHK)
          IF(TCHK)THEN
            CALL LINKEDLIST$GET(LL_CNTL,'KVECSCALE',1,KVECSCALE)
            WRITE(NFILO,*)'INPUT SCALE FACTOR (KVEC=KVECINPUT*(2*PI)/SCALE): ',KVECSCALE
            KVEC1=2.0D0*(4.0D0*ATAN(1.0D0))/KVECSCALE*KVEC1
          ENDIF
          XK1=MATMUL(GBASINV,KVEC1)
        ELSE
          CALL ERROR$MSG('NO XK1 OR KVEC1 GIVEN')
          CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
          CALL ERROR$STOP('BANDS_DIAG')
        ENDIF

        CALL LINKEDLIST$EXISTD(LL_CNTL,'XK2',0,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'KVEC2',0,TCHK2)
        IF(TCHK1)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'XK2',1,XK2)
          WRITE(NFILO,*)'LAST K-POINT: GIVEN IN RELATIVE COORDINATES'
          KVEC2=MATMUL(GBAS,XK2)
        ELSE IF (TCHK2)THEN
          CALL LINKEDLIST$GET(LL_CNTL,'KVEC2',1,KVEC2)
          WRITE(NFILO,*)'LAST K-POINT: GIVEN IN ABSOLUTE COORDINATES'

          CALL LINKEDLIST$EXISTD(LL_CNTL,'KVECSCALE',0,TCHK)
          IF(TCHK)THEN
            CALL LINKEDLIST$GET(LL_CNTL,'KVECSCALE',1,KVECSCALE)
            WRITE(NFILO,*)'INPUT SCALE FACTOR (KVEC=KVECINPUT*(2*PI)/SCALE): ',KVECSCALE
            KVEC2=2.0D0*(4.0D0*ATAN(1.0D0))/KVECSCALE*KVEC2
          ENDIF
          XK2=MATMUL(GBASINV,KVEC2)
        ELSE
          CALL ERROR$MSG('NO XK2 OR KVEC2 GIVEN')
          CALL ERROR$I4VAL('LINE BLOCK NUMBER',ILINE)
          CALL ERROR$STOP('BANDS_DIAG')
        ENDIF

        WRITE(NFILO,*)'FIRST K-POINT IN RELATIVE COORDINATES:',XK1
        WRITE(NFILO,*)'FIRST K-POINT IN ABSOLUTE COORDINATES:',KVEC1

        WRITE(NFILO,*)'LAST K-POINT IN RELATIVE COORDINATES:',XK2
        WRITE(NFILO,*)'LAST K-POINT IN ABSOLUTE COORDINATES:',KVEC2
        
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
        IF(TPRINT)PRINT*,XK1,XK2,NB,NK,NKDIAG,X1,X2,TRIM(FILE)
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
        ALLOCATE(KVECVAL(3,NKDIAG))
        ALLOCATE(XKVAL(3,NKDIAG))
        ALLOCATE(FATBANDVAL(NFATBAND,NKDIAG,NB))
        TPROJ=NFATBAND.GE.1

!       == ITERATE K-POINTS =================================================
!$omp parallel do private(IKDIAG,XK,KVEC,E,PROJ)
        DO IKDIAG=0,NKDIAG-1
          ALLOCATE(PROJ(NAT,NB,LMNXX))
          ALLOCATE(E(NG*NDIM))
!$omp critical
          IF(TPRINT)PRINT*,"IKDIAG",IKDIAG
          XK=XK1+(XK2-XK1)*REAL(IKDIAG,KIND=8)/REAL(MAX(NKDIAG-1,1),KIND=8)
          KVEC=MATMUL(GBAS,XK)
          KVECVAL(:,IKDIAG+1)=KVEC
          XKVAL(:,IKDIAG+1)=XK
          
          WRITE(NFILO,*)'LINE ',ILINE,' OF ',NLINE,' K-POINT ',IKDIAG,' OF ',&
            &NKDIAG-1,' IN RELATIVE COORDINATES ',XK
          WRITE(NFILO,*)'LINE ',ILINE,' OF ',NLINE,' K-POINT ',IKDIAG,' OF ',&
            &NKDIAG-1,' IN ABSOLUTE COORDINATES ',KVEC
!$omp end critical

          CALL BANDS_KPOINT(NG,NB,ISPIN,METHOD_DIAG,GIDG_PROTO,TPROJ,KVEC,&
      &              GVEC,TI_H,TI_S,E,PROJ)
!$omp critical
          IF(TPRINT)PRINT*,'DIAGBANDS_EIG',sqrt(sum(KVEC**2)),E(1:NB)*27.21139D0

          EIGVAL(1:NB,IKDIAG+1)=E(1:NB)*27.21139D0

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
                !FIXME: integrate pdos interface here
                CSVAR=PROJ(FATBANDIAT(IFATBAND),IB,FATBANDIPRO(IFATBAND))
                FATBANDVAL(IFATBAND,IKDIAG+1,IB)=ABS(CSVAR)**2
              ENDDO
            ENDDO
          ELSE
            FATBANDVAL(:,:,:)=0.0D0 
          ENDIF
!$omp end critical
          DEALLOCATE(PROJ)
          DEALLOCATE(E)
        ENDDO
!$omp end parallel do
!
!       =========================================================================
!       ==  WRITE BANDS                                                        ==
!       =========================================================================
                            CALL TRACE$PUSH('WRITE_BANDS')
        CALL FILEHANDLER$UNIT('BANDS',NFILBAND)
        IF(NB.GT.100) THEN
          CALL ERROR$MSG('NUMBER OF BANDS EXCEEDS LIMIT OF 100')
          CALL ERROR$I4VAL('NB ',NB)
          CALL ERROR$STOP('BANDS_PLOTBANDS')
        END IF

        !ITERATE K-POINTS
        DO IKDIAG=1,NKDIAG
          !WRITE EIGENVALUES
          WRITE(NFILBAND,FMT='(103F10.5)')&
            &  sqrt(sum((KVECVAL(:,IKDIAG)-KVECVAL(:,1))**2)),&
            &  EIGVAL(:,IKDIAG)
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
            WRITE(NFILFATBAND,FMT='(103F10.5)',advance="no")&
            &  sqrt(sum((KVECVAL(:,IKDIAG)-KVECVAL(:,1))**2))

            FATBANDMAX=MAXVAL(FATBANDVAL(IFATBAND,:,:))
            IF(TPRINT)PRINT*,"FATBAND FATBANDMAX",FATBANDMAX
            DO ILINESPERBAND=0,LINESPERBAND(IFATBAND)-1
              SVAR1=2.0D0/REAL(LINESPERBAND(IFATBAND)-1,KIND=8)*REAL(ILINESPERBAND,KIND=8)-1.0D0
              SVAR2=0.5D0*SVAR1*FATBANDMAXWIDTH(IFATBAND)/FATBANDMAX
              WRITE(NFILFATBAND,FMT='(100F10.5)',advance="no")EIGVAL(:,IKDIAG)+SVAR2*FATBANDVAL(IFATBAND,IKDIAG,:) 
            ENDDO
            WRITE(NFILFATBAND,FMT='(a)')" "
          ENDDO
        ENDDO
                            CALL TRACE$POP
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
      CALL TIMING$PRINT('ALL',NFILO)
                            CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDS_DIAG

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BANDS_DOS(LL_CNTL,NFIL,NFILO)
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
      TYPE(EWGHT_TYPE),allocatable :: EWGHT(:,:)
      REAL(8),ALLOCATABLE          :: A(:,:)
      REAL(8)                      :: RNTOT
      REAL(8)                      :: EF
      REAL(8)                      :: SUMA
      REAL(8)                      :: EMIN,EMAX
      INTEGER(4)                   :: NE
      REAL(8),ALLOCATABLE          :: DOS(:,:)
     
      INTEGER(4)                   :: NG
      INTEGER(4)                   :: METHOD_DIAG 
      INTEGER(4)                   :: GIDG_PROTO
      LOGICAL(4)                   :: TPROJ
      REAL(8)                      :: KVEC(3)
      REAL(8),ALLOCATABLE          :: GVEC(:,:)
      REAL(8),ALLOCATABLE          :: G2(:)
      COMPLEX(8),ALLOCATABLE       :: TI_H(:,:,:)
      COMPLEX(8),ALLOCATABLE       :: TI_S(:,:,:)
      REAL(8),ALLOCATABLE          :: E(:)
      COMPLEX(8),ALLOCATABLE       :: PROJ(:,:,:)
      COMPLEX(8),ALLOCATABLE       :: PROJK(:,:,:,:)
      REAL(8)                      :: GWEIGHT

      LOGICAL(4)                   :: TUSESYM
      INTEGER(4)                   :: SPACEGROUP
      REAL(8)                      :: A0,B0,C0,ALPHA,BETA,GAMMA
      INTEGER(4)                   :: NSYM,NOP
      INTEGER(4),parameter         :: NSYMX=48
      INTEGER(4),parameter         :: NOPX=48
      INTEGER(4)                   :: IARB(3)
      CHARACTER(3)                 :: BRAVAIS
      INTEGER(4)                   :: IIO(3,3,NOPX)
      REAL(8)                      :: C(3,NOPX)
      INTEGER(4)                   :: ISYM
      LOGICAL(4)                   :: TSHIFT

      INTEGER(4)                   :: NFILOUT,NFILIN
     
      LOGICAL(4)                   :: TPRINT=.FALSE.
!     **************************************************************************
      IF(NDIM.eq.2)THEN
        CALL ERROR$MSG('ONLY NDIM=1 AND NSPIN=1 OR 2 IMPLEMENTED AND TESTED')
        CALL ERROR$STOP('BANDS_DOS')
      ENDIF
                            CALL TRACE$PUSH('BANDS_DOS')
      CALL TIMING$START
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
          CALL ERROR$STOP('BANDS_DOS')
        ENDIF
      ELSE IF (NFILE.gt.1)then
        CALL ERROR$MSG('MULTIPLE INPUT FILES GIVEN (!BCNTL!INPUTFILE)')
        CALL ERROR$STOP('BANDS_DOS')
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
        CALL ERROR$STOP('BANDS_DOS')
      ENDIF

      !FIXME: TO BE READ FROM BCNTL/DCNTL
      NKDIV(1)=20
      NKDIV(2)=20
      NKDIV(3)=20
      ISHIFT(1)=0
      ISHIFT(2)=0
      ISHIFT(3)=0
      NB=10
      !EPW=0.5D0*15.D0
      TPROJ=.FALSE.
      IF(NSPIN.eq.1)THEN
        RNTOT=0.5D0*NEL
      ELSE
        RNTOT=NEL
      ENDIF

      !check if NB*NSPIN>RNTOT
      IF(NB*NSPIN.lt.RNTOT)THEN
        CALL ERROR$MSG('NB*NPSIN IS LESS THAN THE NUMBER OF BANDS, PLEASE INCREASE NB')
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$R8VAL('NUMBER OF BANDS',RNTOT)
        CALL ERROR$STOP('BANDS_DOS')
      ENDIF

      METHOD_DIAG=2
      NE=1000
      TUSESYM=.false.
      SPACEGROUP=229
      TSHIFT=.FALSE. 


                            CALL TRACE$PASS('AFTER READ BNCTL')
      IF(TUSESYM)THEN
!     ==========================================================================
!     ==  FIND IRREDUCIBLE K-POINTS AND TETRAHEDRA                            ==
!     ==========================================================================
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
        NKP=NKDIV(1)*NKDIV(2)*NKDIV(3)
        CALL SPACEGROUP$SETI4('SPACEGROUP',SPACEGROUP)
        CALL SPACEGROUP$GETCH('BRAVAIS',BRAVAIS)
        CALL BRILLOUIN$CHECKRBAS(BRAVAIS,A0,B0,C0,ALPHA,BETA,GAMMA,RBAS)

!        CALL BRILLOUIN_TESTCOMPLETE()
        CALL SPACEGROUP$GENERATORS('RECI',NOPX,NOP,IIO,C)
        WRITE(*,FMT='(82("="),T10," GENERATORS OF THE GROUP ")')
        DO I=1,NOP
          WRITE(*,FMT='(I5,T20,3("|",3I5,"|"))')I,IIO(:,:,I)
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
      IF(TPRINT)write(*,*)' NKP nach BRILLOUIN$GETI4',NKP  
      ALLOCATE(BK(3,NKP))
      ALLOCATE(XK(3,NKP))
      CALL BRILLOUIN$GETR8A('K',3*NKP,BK)

      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      XK=MATMUL(RBASINV,BK)

!
!     =========================================================================
!     ==  CONSTRUCT K-INDEPENDENT PART OF HAMILTONIAN                        ==
!     =========================================================================
      CALL BANDS_KINDEP(NG,GVEC,G2,GWEIGHT,GIDG_PROTO,TI_H,TI_S)
      ALLOCATE(EB(NB*NSPIN,NKP))        
      ALLOCATE(WGHT(NB*NSPIN,NKP))
      
      IF(TPROJ)THEN
        ALLOCATE(PROJK(NAT,NB*NSPIN,LMNXX,IKP))
      ENDIF
!       == ITERATE K-POINTS =================================================
                            CALL TRACE$PUSH('ITERATE KPOINTS')
!$omp parallel do private(IKP,KVEC,E,PROJ,ISPIN)
      DO IKP=1,NKP
        DO ISPIN=1,NSPIN
          KVEC=BK(:,IKP)
!$omp critical
          WRITE(NFILO,*)'K-POINT ',IKP,' OF ',NKP,' FOR SPIN ',ISPIN,' IN ABSOLUTE COORDINATES ',KVEC
!$omp end critical
          ALLOCATE(PROJ(NAT,NB,LMNXX))
          ALLOCATE(E(NG*NDIM))
          CALL BANDS_KPOINT(NG,NB,ISPIN,METHOD_DIAG,GIDG_PROTO,TPROJ,KVEC,&
      &              GVEC,TI_H,TI_S,E,PROJ)
!$omp critical
          EB(1+NB*(ISPIN-1):NB+NB*(ISPIN-1),IKP)=E(1:NB)*27.21139D0
          IF(TPROJ)THEN
            PROJK(:,1+NB*(ISPIN-1):NB+NB*(ISPIN-1),:,IKP)=PROJ(:,1:NB,:)
          ENDIF
!$omp end critical
          DEALLOCATE(PROJ)
          DEALLOCATE(E)
        ENDDO
      ENDDO
!$omp end parallel do
                            CALL TRACE$POP()
      DO IKP=1,NKP
        DO I=1,NB
          PRINT*,"EB(",I,",",IKP,")=",EB(I,IKP)          
        ENDDO
      ENDDO

!bcc Fe, 555
!ALLOCATE(EBTMP(10,4))        
!EBTMP(1,1)=8.90519471374423D0
!EBTMP(2,1)=15.9434033733550D0
!EBTMP(3,1)=15.9437318985426D0
!EBTMP(4,1)=15.9439782149963D0
!EBTMP(5,1)=17.4331464451587D0
!EBTMP(6,1)=17.4332575780743D0
!EBTMP(7,1)=42.1940822516814D0
!EBTMP(8,1)=42.1948986745993D0
!EBTMP(9,1)=42.1955108417937D0
!EBTMP(10,1)=47.3335853129882D0
!EBTMP(1,2)=12.3673809615346D0
!EBTMP(2,2)=15.1409086127278D0
!EBTMP(3,2)=15.7856215226624D0
!EBTMP(4,2)=17.2826678202824D0
!EBTMP(5,2)=17.5522773569203D0
!EBTMP(6,2)=17.7871577224451D0
!EBTMP(7,2)=29.6606949980102D0
!EBTMP(8,2)=37.2043315559138D0
!EBTMP(9,2)=39.9091911884596D0
!EBTMP(10,2)=41.1678488884605D0
!EBTMP(1,3)=14.5347471996221D0
!EBTMP(2,3)=14.9176142640354D0
!EBTMP(3,3)=14.9185269579434D0
!EBTMP(4,3)=17.5706279784390D0
!EBTMP(5,3)=17.5748786681166D0
!EBTMP(6,3)=24.2776710779338D0
!EBTMP(7,3)=24.2962112814693D0
!EBTMP(8,3)=27.1274063016807D0
!EBTMP(9,3)=36.6869486713712D0
!EBTMP(10,3)=45.8183345525610D0
!EBTMP(1,4)=13.8692342529125D0
!EBTMP(2,4)=14.6195957197494D0
!EBTMP(3,4)=17.0029133195168D0
!EBTMP(4,4)=17.0710170479733D0
!EBTMP(5,4)=17.8472980812975D0
!EBTMP(6,4)=20.5994303387914D0
!EBTMP(7,4)=31.1501089743446D0
!EBTMP(8,4)=31.3120829681898D0
!EBTMP(9,4)=35.6501004998864D0
!EBTMP(10,4)=36.2712404230226D0
!
!EB(:,1)=EBTMP(:,1)
!EB(:,2)=EBTMP(:,2)
!EB(:,3)=EBTMP(:,2)
!EB(:,4)=EBTMP(:,3)
!EB(:,5)=EBTMP(:,2)
!EB(:,6)=EBTMP(:,2)
!EB(:,7)=EBTMP(:,3)
!EB(:,8)=EBTMP(:,2)
!EB(:,9)=EBTMP(:,3)
!EB(:,10)=EBTMP(:,3)
!EB(:,11)=EBTMP(:,4)
!EB(:,12)=EBTMP(:,2)
!EB(:,13)=EBTMP(:,4)
!EB(:,14)=EBTMP(:,4)

!     ==========================================================================
!     ==  WRITE PDOS FILE                                                     ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PDOS',NFILIN)
      REWIND(NFILIN)
      CALL PDOS$READ(NFILIN)
      CALL PDOS$SETI4('NKPT',NKP)
      CALL PDOS$SETI4('NB',NB)
      CALL PDOS$SETR8A('XK',3*NKP,XK)
      
      CALL FILEHANDLER$UNIT('PDOSOUT',NFILOUT)
      REWIND NFILOUT
      CALL PDOS$WRITE(NFILOUT)
      CALL LIB$FLUSHFILE(NFILOUT)
      CALL FILEHANDLER$CLOSE('PDOSOUT')
      

STOP

!     ==========================================================================
!     ==  CALCULATE WEIGHTS                                                   ==
!     ==========================================================================
      IF(NSPIN.eq.2) WRITE(NFILO,*)'WARNING: USING THE SAME FERMI ENERGY FOR&
&      BOTH SPIN DIRECTIONS.'
      CALL BRILLOUIN$DOS(NSPIN*NB,NKP,EB,WGHT,RNTOT,EF)
!
!     ==========================================================================
!     ==  PERFORM BRILLOUIN ZONE INTEGRATION OF A(K)                          ==
!     ==========================================================================
      !FIXME TOTAL DENSITY for testing
      ALLOCATE(A(NB*NSPIN,NKP))
      A=1
      SUMA=0.D0
      DO IB=1,NB
        DO ISPIN=1,NSPIN
          SUMA=0.0D0
          DO IKP=1,NKP
            SUMA=SUMA+WGHT(IB+NB*(ISPIN-1),IKP)*A(IB+NB*(ISPIN-1),IKP)
          ENDDO
          PRINT*,"IB=",IB," ISPIN=",ISPIN," SUMA=",SUMA
        ENDDO
      ENDDO
      
      A=1
      SUMA=0.D0
      DO IB=1,NB
        DO ISPIN=1,NSPIN
          DO IKP=1,NKP
            SUMA=SUMA+WGHT(IB+NB*(ISPIN-1),IKP)*A(IB+NB*(ISPIN-1),IKP)
          ENDDO
        ENDDO
      ENDDO
      PRINT*,'INTEGRAL OF A : ',SUMA,' should be ',RNTOT 
 
      !TOTAL DENSITY OF STATES
      EMIN=minval(EB(:,:))
      EMAX=maxval(EB(:,:))
      allocate(EWGHT(NSPIN*NB,NKP))
      CALL BRILLOUIN$EWGHT(NKP,NB*NSPIN,EB,EMIN,EMAX,NE,EWGHT)

      A(:,:)=1.0d0
      ALLOCATE(DOS(NE,NSPIN))
      DOS(:,:)=0.0D0
      DO IKP=1,NKP
        DO IB=1,NB
          DO ISPIN=1,NSPIN
            I1=EWGHT(IB+NB*(ISPIN-1),IKP)%I1
            I2=EWGHT(IB+NB*(ISPIN-1),IKP)%I2
            DO IE=I1,I2
              DOS(IE,ISPIN)=DOS(IE,ISPIN)+EWGHT(IB+NB*(ISPIN-1),IKP)%WGHT(IE)*A(IB+NB*(ISPIN-1),IKP)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      !FIXME: make output nice
      open(101,file='dos.dat')
      do ie=1,ne
      write(101,*)emin+real(ie-1,kind=8)*(emax-emin)/real(ne-1,kind=8),(DOS(ie,ISPIN),ISPIN=1,NSPIN)
      enddo
      close(101)
                            CALL TRACE$POP()
      RETURN
      END SUBROUTINE



