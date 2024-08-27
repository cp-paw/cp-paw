!************************************************************************
!**                                                                    **
!**  NAME: WAVEPLOT                                                    **
!**                                                                    **
!**  PURPOSE: READS A WAVE FUNCTION FILE AND WRITES INPUT FILES        **
!**    FOR VISUALIZATION PROGRAMS                                      **
!**                                                                    **
!************************************************************************
       PROGRAM MAIN
       USE STRINGS_MODULE
       USE LINKEDLIST_MODULE
       IMPLICIT NONE
       TYPE(LL_TYPE)             :: LL_CNTL
       TYPE(LL_TYPE)             :: LL_STRC
       INTEGER(4)                :: NFIL
       INTEGER(4)                :: NFILO
       CHARACTER(256)            :: FILE
       CHARACTER(32)             :: ID
       LOGICAL(4)                :: TCHK
       INTEGER(4)                :: I
       INTEGER(4)                :: NLISTS
!      ******************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT

       I=COMMAND_ARGUMENT_COUNT()
       IF(I.NE.1) THEN
         WRITE(*,FMT='(A)')'CORRECT USAGE: PAW_WAVE.X [FILE]'
         WRITE(*,FMT='(A)') &
    &                  'WHERE [FILE] IS THE FILE CONTROL FILE NAME OF THIS TOOL'
         WRITE(*,FMT='(A)')'SEE PAW_MANUAL FOR FURTHER DESCRIPTION'
         CALL ERROR$NORMALSTOP
       END IF
       CALL  GET_COMMAND_ARGUMENT(1,FILE)
       IF(FILE.EQ.'?'.OR.FILE.EQ.-'-H') THEN
         WRITE(*,FMT='(A)')'CORRECT USAGE: PAW_WAVE [FILE]'
         WRITE(*,FMT='(A)') &
    &                  'WHERE [FILE] IS THE FILE CONTROL FILE NAME OF THIS TOOL'
         WRITE(*,FMT='(A)')'SEE PAW_MANUAL FOR FURTHER DESCRIPTION'
         CALL ERROR$NORMALSTOP
       END IF
!      ==================================================================
!      ==  INITIALIZE FILEHANDLER                                      ==
!      ==================================================================
      I=INDEX(FILE,'.',BACK=.TRUE.)
      CALL FILEHANDLER$SETROOT(FILE(1:I-1))
      CALL FILEHANDLER$SETFILE('CNTL',.FALSE.,TRIM(FILE))
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC_OUT')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('PROT',.TRUE.,-'.WPROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('WAVE',.TRUE.,-'.WV')
      CALL FILEHANDLER$SETSPECIFICATION('WAVE','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('WAVE','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('WAVE','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('WAVE','FORM','UNFORMATTED')
      CALL FILEHANDLER$SETFILE('WAVEDX',.TRUE.,-'.DX')
      CALL FILEHANDLER$SETSPECIFICATION('WAVEDX','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('WAVEDX','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('WAVEDX','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('WAVEDX','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('CUBE',.TRUE.,-'.CUB')
      CALL FILEHANDLER$SETSPECIFICATION('CUBE','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('CUBE','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('CUBE','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('CUBE','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('VRML',.TRUE.,-'.WRL')
      CALL FILEHANDLER$SETSPECIFICATION('VRML','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('VRML','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('VRML','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('VRML','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('GNUCONTOUR',.TRUE.,-'_C.GNU')
      CALL FILEHANDLER$SETSPECIFICATION('GNUCONTOUR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('GNUCONTOUR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('GNUCONTOUR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('GNUCONTOUR','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('GNURUBBERSHEET',.TRUE.,-'_R.GNU')
      CALL FILEHANDLER$SETSPECIFICATION('GNURUBBERSHEET','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('GNURUBBERSHEET','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('GNURUBBERSHEET','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('GNURUBBERSHEET','FORM','FORMATTED')
!
!     ==========================================================================
!     ==  READ CNTL FILE TO LINKEDLIST                                        ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CNTL',NFIL)
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)
!
!     ==========================================================================
!     ==  RESET FILE NAME FOR STRUCTURE FILE IF REQUESTED                     ==
!     ==========================================================================
      CALL TRACE$PASS('RESET FILENAME FOR STRC FILE')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'WCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NLISTS)
      DO I=1,NLISTS
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',I)
!
!       == GET FILE ID =================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('KEYWORD ID IS MANDATORY')
          CALL ERROR$STOP('MAIN')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
        IF(ID.NE.'STRC'.AND.ID.NE.'WAVE' &
     &                 .AND.ID.NE.'WAVEDX' &
     &                 .AND.ID.NE.'CUBE' &
     &                 .AND.ID.NE.'WRL' &
     &                 .AND.ID.NE.'GNUCONTOUR' &
     &                 .AND.ID.NE.'GNURUBBERSHEET') THEN
         CALL ERROR$MSG('!WCNTL!FILE!FILES:ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('MAIN')
        END IF
!
!       == GET FILE NAME ================================================
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
!
!       == GET EXTENSION FLAG ==========================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        END IF
        CALL FILEHANDLER$SETFILE(ID,TCHK,FILE)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
      CALL TRACE$PASS('FILENAME FOR STRC FILE SET')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
!
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL FILEHANDLER$REPORT(NFILO,'ALL')
!
!     ==================================================================
!     ==  READ STRC FILE TO LINKEDLIST                                ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$NEW(LL_STRC)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$MARK(LL_STRC,1)
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL MWAVE(LL_CNTL,LL_STRC)
      
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'WCNTL')
      CALL LINKEDLIST$REPORT_UNUSED(LL_CNTL,NFILO)

      !!LL_STRC IS ONLY USED IN GETIZ
      !CALL LINKEDLIST$SELECT(LL_STRC,'~')
      !CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      !CALL LINKEDLIST$REPORT_UNUSED(LL_STRC,NFILO)
!
      CALL ERROR$NORMALSTOP()
      STOP
      END
!
!      ..................................................................
       SUBROUTINE MWAVE(LL_CNTL_,LL_STRC_)
       USE STRINGS_MODULE
       USE LINKEDLIST_MODULE
       USE PERIODICTABLE_MODULE
       IMPLICIT NONE
       TYPE(LL_TYPE), INTENT(IN) :: LL_CNTL_
       TYPE(LL_TYPE), INTENT(IN) :: LL_STRC_
       TYPE(LL_TYPE)             :: LL_CNTL
       TYPE(LL_TYPE)             :: LL_STRC
       CHARACTER(256)            :: TITLE
       INTEGER(4)                :: NFIL
       INTEGER(4)                :: NFILO !FORTRAN UNIT OF PROTOCOL FILE
       REAL(8)                   :: RBAS(3,3)
       INTEGER(4)                :: NAT
       REAL(8)      ,ALLOCATABLE :: POS(:,:)   !(3,NAT)
       REAL(8)      ,ALLOCATABLE :: Z(:)       !(NAT)
       CHARACTER(32),ALLOCATABLE :: ATOMNAME(:)
       INTEGER(4)                :: NR1,NR2,NR3
       REAL(8)      ,ALLOCATABLE :: WAVE(:,:,:)
       COMPLEX(8)   ,ALLOCATABLE :: CWAVE(:,:,:)
       REAL(8)                   :: BOXR0(3)
       REAL(8)                   :: BOXVEC(3,3)
       REAL(8)                   :: PLANER0(3)
       REAL(8)                   :: PLANEVEC(3,2)
       REAL(8)      ,ALLOCATABLE :: POSM(:,:)  ! POSITIONS
       INTEGER(4)   ,ALLOCATABLE :: MAP(:)     !MAPPING FROM MODEL 
       REAL(8)      ,ALLOCATABLE :: RAD(:)   
       REAL(8)      ,ALLOCATABLE :: COLOR(:,:)
       INTEGER(4)                :: NBOND
       INTEGER(4)   ,ALLOCATABLE :: BOND(:,:)
       LOGICAL(4)                :: TCHK,TCHK1
       INTEGER(4)                :: NATM
       INTEGER(4)                :: IAT
       INTEGER(4)                :: IVEC(3)
       REAL(8)      ,ALLOCATABLE :: SCALEDRAD(:)
       LOGICAL(4)                :: TPLANE
       LOGICAL(4)                :: TC
       REAL(8)                   :: XK(3)
       REAL(8)                   :: SVAR
       INTEGER(4)                :: I1,I2,I3
!      *************************************************************************
       CALL TRACE$PUSH('MWAVE')
       LL_CNTL=LL_CNTL_
       LL_STRC=LL_STRC_
       CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!      =========================================================================
!      ==  READ INPUT FILE                                                    ==
!      =========================================================================
       CALL TRACE$PASS('READING WAVEPLOT')
       CALL FILEHANDLER$UNIT('WAVE',NFIL)
       CALL READWAVEPLOTDIM(NFIL,TITLE,NAT,NR1,NR2,NR3,TC)
       ALLOCATE(POS(3,NAT))
       ALLOCATE(Z(NAT))
       ALLOCATE(ATOMNAME(NAT))
       IF(TC) THEN
         ALLOCATE(CWAVE(NR1,NR2,NR3))
         CALL READWAVEPLOTC(NFIL,TITLE,RBAS,NAT,POS,Z,ATOMNAME,XK &
      &                    ,NR1,NR2,NR3,CWAVE)
         ALLOCATE(WAVE(NR1,NR2,NR3))
         WAVE=REAL(CWAVE)
       ELSE
         ALLOCATE(WAVE(NR1,NR2,NR3))
         CALL READWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,ATOMNAME,XK &
      &                   ,NR1,NR2,NR3,WAVE)
       END IF
       CALL TRACE$PASS('WAVEPLOT FILE READ')
!
!      =========================================================================
!      ==  SET MINIMUM AND MAXIMUM VALUES                                     ==
!      =========================================================================
       CALL LINKEDLIST$SELECT(LL_CNTL,'~')
       CALL LINKEDLIST$SELECT(LL_CNTL,'WCNTL')
       CALL LINKEDLIST$EXISTL(LL_CNTL,'GENERIC',1,TCHK)
       IF(TCHK) THEN
         CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
         CALL LINKEDLIST$EXISTD(LL_CNTL,'MIN',1,TCHK1)
         IF(TCHK1) THEN
           CALL LINKEDLIST$GET(LL_CNTL,'MIN',1,SVAR)
           DO I1=1,NR1
             DO I2=1,NR2
               DO I3=1,NR3
                 WAVE(I1,I2,I3)=MAX(WAVE(I1,I2,I3),SVAR)
               ENDDO
             ENDDO
           ENDDO
           IF(TC)THEN
             DO I1=1,NR1
               DO I2=1,NR2
                 DO I3=1,NR3
                   CWAVE(I1,I2,I3)=CMPLX(MAX(REAL(CWAVE(I1,I2,I3)),SVAR) &
      &                                ,MAX(AIMAG(CWAVE(I1,I2,I3)),SVAR),KIND=8)
                 ENDDO
               ENDDO
             ENDDO
           ENDIF
         END IF     
         CALL LINKEDLIST$EXISTD(LL_CNTL,'MAX',1,TCHK1)
         IF(TCHK1) THEN
           CALL LINKEDLIST$GET(LL_CNTL,'MAX',1,SVAR)
           DO I1=1,NR1
             DO I2=1,NR2
               DO I3=1,NR3
                 WAVE(I1,I2,I3)=MIN(WAVE(I1,I2,I3),SVAR)
               ENDDO
             ENDDO
           ENDDO
           IF(TC)THEN
             DO I1=1,NR1
               DO I2=1,NR2
                 DO I3=1,NR3
                   CWAVE(I1,I2,I3)=CMPLX(MIN(REAL(CWAVE(I1,I2,I3)),SVAR) &
     &                                 ,MIN(AIMAG(CWAVE(I1,I2,I3)),SVAR),KIND=8)
                 ENDDO
               ENDDO
             ENDDO
           ENDIF
         END IF     
       END IF
!
!      =========================================================================
!      ==  GET VIEWBOX                                                        ==
!      =========================================================================
       CALL LINKEDLIST$SELECT(LL_CNTL,'~')
       CALL LINKEDLIST$SELECT(LL_CNTL,'WCNTL')
       CALL LINKEDLIST$SELECT(LL_CNTL,'VIEWBOX')
       BOXR0(:)=0.D0
       BOXVEC(:,:)=RBAS(:,:)
       CALL LINKEDLIST$EXISTD(LL_CNTL,'T',1,TCHK)
       IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'T',1,BOXVEC)
       CALL LINKEDLIST$EXISTD(LL_CNTL,'O',1,TCHK)
       IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'O',1,BOXR0)
       CALL LINKEDLIST$EXISTD(LL_CNTL,'C',1,TCHK1)
       IF(TCHK1) THEN
         IF(TCHK) THEN
           CALL ERROR$MSG('!WCNTL!VIEWBOX:O AND :C ARE MUTUALLY EXCLUSIVE')
           CALL ERROR$STOP('MWAVE')
         END IF
         CALL LINKEDLIST$GET(LL_CNTL,'C',1,BOXR0)
         BOXR0=BOXR0-0.5D0*MATMUL(BOXVEC,(/1.D0,1.D0,1.D0/))
       END IF
!
!      ==================================================================
!      ==  GET PLANE FOR RUBBERSHEET                                   ==
!      ==================================================================
       CALL LINKEDLIST$SELECT(LL_CNTL,'~')
       CALL LINKEDLIST$SELECT(LL_CNTL,'WCNTL')
       CALL LINKEDLIST$EXISTL(LL_CNTL,'PLANE',1,TPLANE)
       PLANER0(:)=BOXR0
       PLANEVEC(:,:)=BOXVEC(:,:2)
       IF(TPLANE) THEN
         CALL LINKEDLIST$SELECT(LL_CNTL,'PLANE')
         CALL LINKEDLIST$EXISTD(LL_CNTL,'T',1,TCHK)
         IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'T',1,PLANEVEC)
!        == MAKE SECOND PLANE VECTOR PERPENDICULAR TO THE FIRST ONE ============
         PLANEVEC(:,2)=PLANEVEC(:,2)-PLANEVEC(:,1) &
      &                       /DOT_PRODUCT(PLANEVEC(:,1),PLANEVEC(:,1)) &
      &                       *DOT_PRODUCT(PLANEVEC(:,1),PLANEVEC(:,2)) 

         CALL LINKEDLIST$EXISTD(LL_CNTL,'O',1,TCHK)
         IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'O',1,PLANER0)
         CALL LINKEDLIST$EXISTD(LL_CNTL,'C',1,TCHK1)
         IF(TCHK1) THEN
           IF(TCHK) THEN
             CALL ERROR$MSG('!WCNTL!VIEWBOX:O AND :C ARE MUTUALLY EXCLUSIVE')
             CALL ERROR$STOP('MWAVE')
           END IF
           CALL LINKEDLIST$GET(LL_CNTL,'C',1,PLANER0)
           BOXR0=BOXR0-0.5D0*MATMUL(PLANEVEC,(/1.D0,1.D0/))
         END IF
       END IF
!
!      ==================================================================
!      ==  WRITE DENSITY TO CUBE FILE                                  ==
!      ==================================================================
!      == ATTENTION! ONLY REAL PART IS USED
       CALL TRACE$PASS('WRITE WAVE TO CUBE-FILE')
       CALL FILEHANDLER$UNIT('CUBE',NFIL)
       REWIND NFIL
       CALL MAKECUBE(NFIL,NAT,Z,POS,RBAS,NR1,NR2,NR3,WAVE,BOXR0,BOXVEC)
       CALL FILEHANDLER$CLOSE('CUBE')
!
!      =========================================================================
!      ==  WRITE RUBBERSHEET TO GNU FILE                                      ==
!      =========================================================================
       IF(TPLANE) THEN
         CALL TRACE$PASS('MAKE GNU FILE FOR CONTOUR PLOT')
!         == ATTENTION! ONLY REAL PART IS USED =================================
         CALL FILEHANDLER$UNIT('GNUCONTOUR',NFIL)
         REWIND NFIL
         CALL MAKEGNU(NFIL,'CONTOUR',NAT,Z,POS,RBAS,NR1,NR2,NR3,WAVE &
      &              ,PLANER0,PLANEVEC)
         CALL FILEHANDLER$CLOSE('GNUCONTOUR')
         CALL TRACE$PASS('MAKE GNU FILE FOR RUBBERSHEET')
         CALL FILEHANDLER$UNIT('GNURUBBERSHEET',NFIL)
         REWIND NFIL
         CALL MAKEGNU(NFIL,'SURFACE',NAT,Z,POS,RBAS,NR1,NR2,NR3,WAVE &
      &              ,PLANER0,PLANEVEC)
         CALL FILEHANDLER$CLOSE('GNURUBBERSHEET')
       END IF
!
!      ==================================================================
!      ==  WRITE VRML SCENE TO WRL FILE                                ==
!      ==================================================================
       CALL TRACE$PASS('WRITE WAVE TO VRML-FILE')
       CALL FILEHANDLER$UNIT('VRML',NFIL)
       REWIND NFIL
       IF(TC) THEN
         CALL MAKEVRMLC(NFIL,NAT,Z,POS,RBAS,XK,NR1,NR2,NR3,CWAVE &
      &             ,BOXR0,BOXVEC,TPLANE,PLANER0,PLANEVEC)
       ELSE
         CALL MAKEVRML(NFIL,NAT,Z,POS,RBAS,NR1,NR2,NR3,WAVE &
      &             ,BOXR0,BOXVEC,TPLANE,PLANER0,PLANEVEC)
       END IF
       CALL FILEHANDLER$CLOSE('VRML')
!
!      ==================================================================
!      ==  WRITE DENSITY TO DATAEXPLORER FILE                          ==
!      ==================================================================
       CALL TRACE$PASS('WRITE WAVE TO DX-FILE')
       CALL FILEHANDLER$UNIT('WAVEDX',NFIL)
       REWIND NFIL
       CALL DX_DENSITY(NFIL,NR1,NR2,NR3,RBAS,WAVE,BOXR0,BOXVEC)
       DEALLOCATE(WAVE)
       CALL TRACE$PASS('WAVE WRITTEN')
!
!      ==================================================================
!      ==  MAP ATOMS INTO VIEWBOX                                      ==
!      ==================================================================
       CALL MODEL$NATM(BOXR0,BOXVEC,NAT,POS,RBAS,NATM)
       ALLOCATE(POSM(3,NATM))
       ALLOCATE(MAP(NATM))
       CALL MODEL$ATOMS(BOXR0,BOXVEC,NAT,POS,RBAS,NATM,POSM,MAP)
       ALLOCATE(RAD(NATM))
       ALLOCATE(SCALEDRAD(NATM))
       ALLOCATE(COLOR(3,NATM))
       DO IAT=1,NATM
         CALL PERIODICTABLE$GET(NINT(Z(MAP(IAT))),'R(COV)',RAD(IAT))
         CALL ATOMCOLOR(NINT(Z(MAP(IAT))),IVEC)
         COLOR(:,IAT)=REAL(IVEC,KIND=8)/200.D0
       ENDDO
       SCALEDRAD(:)=1.2D0*RAD(:)
       CALL MODEL$NBONDM(NATM,POSM,SCALEDRAD,NBOND)
       ALLOCATE(BOND(2,NBOND))
       CALL MODEL$BONDS(NATM,POSM,SCALEDRAD,NBOND,BOND)
!
!      ==================================================================
!      ==  WRITE BALLSTICK MODEL TO DATAEXPLORER FILE                  ==
!      ==================================================================
       CALL TRACE$PASS('BEFORE BALLSTICK')
       SCALEDRAD(:)=0.5D0*RAD(:)
       CALL DXBALLSTICK(NFIL,NATM,COLOR,SCALEDRAD,POSM,NBOND,BOND,BOXR0,BOXVEC)
       CALL TRACE$PASS('AFTER BALLSTICK')
!
!      ==================================================================
!      ==  FINISH DATAEXPLORER FILE                                    ==
!      ==================================================================
       WRITE(NFIL,*)-'END'
       CALL FILEHANDLER$CLOSE('WAVEDX')
!
!      ==================================================================
!      ==  CLOSE DOWN                                                  ==
!      ==================================================================
       DEALLOCATE(POS)
       DEALLOCATE(Z)
       DEALLOCATE(ATOMNAME)
       CALL TRACE$POP
       RETURN
       END
!
!      ..................................................................
       SUBROUTINE GETIZ(LL_STRC_,NAT,IZ)
       USE LINKEDLIST_MODULE
       USE PERIODICTABLE_MODULE
       TYPE(LL_TYPE),INTENT(IN)  :: LL_STRC_
       TYPE(LL_TYPE)             :: LL_STRC
       INTEGER(4)   ,INTENT(IN)  :: NAT
       INTEGER(4)   ,INTENT(OUT) :: IZ(NAT)
       CHARACTER(16)             :: SP(NAT)
       CHARACTER(16)             :: SPNAME
       INTEGER(4)                :: I,IAT,ISVAR
       INTEGER(4)                :: NSP
!      ******************************************************************
       LL_STRC=LL_STRC_
       CALL LINKEDLIST$SELECT(LL_STRC,'~')
       CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
       CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',ISVAR)
       IF(ISVAR.NE.NAT) THEN
         CALL ERROR$MSG('INCONSISTENT WAVEPLOT FILE AND STRC_OUT FILE')
         CALL ERROR$STOP('MWAVE')
       ENDIF
       DO I=1,NAT
         CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',I)
         CALL LINKEDLIST$GET(LL_STRC,'INDEX',1,IAT)
         CALL LINKEDLIST$GET(LL_STRC,'SP',1,SP(IAT))
         CALL LINKEDLIST$SELECT(LL_STRC,'..')
       ENDDO
       CALL LINKEDLIST$NLISTS(LL_STRC,'SPECIES',NSP)
       DO I=1,NSP
         CALL LINKEDLIST$SELECT(LL_STRC,'SPECIES',I)
         CALL LINKEDLIST$GET(LL_STRC,'NAME',1,SPNAME)
         CALL PERIODICTABLE$GET(SPNAME(1:2),'Z',ISVAR)
         DO IAT=1,NAT
           IF(SPNAME.EQ.SP(IAT)) IZ(IAT)=ISVAR
         ENDDO
         CALL LINKEDLIST$SELECT(LL_STRC,'..')
       ENDDO
       RETURN
       END 
!
!      ..................................................................
       SUBROUTINE MODEL$NATM(BOXR0,BOXVEC,NAT,POS,RBAS,NATM)
!      ******************************************************************
!      **                                                              **
!      ******************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: BOXR0(3)
       REAL(8)   ,INTENT(IN) :: BOXVEC(3,3)
       INTEGER(4),INTENT(IN) :: NAT
       REAL(8)   ,INTENT(IN) :: POS(3,NAT)
       REAL(8)   ,INTENT(IN) :: RBAS(3,3)
       INTEGER(4),INTENT(OUT):: NATM
       INTEGER(4)            :: IAT,I1,I2,I3
       REAL(8)               :: BOXVECIN(3,3)
       REAL(8)               :: DET
       REAL(8)               :: R(3)
       REAL(8)               :: XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3
       INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
       REAL(8)               :: T1,T2,T3
       LOGICAL(4)            :: TCHK
!      ******************************************************************
!
!      ==================================================================
!      ==  INVERSE OF BOXVEC                                           ==
!      ==================================================================
       BOXVECIN(1,1)=BOXVEC(2,2)*BOXVEC(3,3)-BOXVEC(3,2)*BOXVEC(2,3)
       BOXVECIN(2,1)=BOXVEC(3,2)*BOXVEC(1,3)-BOXVEC(1,2)*BOXVEC(3,3)
       BOXVECIN(3,1)=BOXVEC(1,2)*BOXVEC(2,3)-BOXVEC(2,2)*BOXVEC(1,3)
       BOXVECIN(1,2)=BOXVEC(2,3)*BOXVEC(3,1)-BOXVEC(3,3)*BOXVEC(2,1)
       BOXVECIN(2,2)=BOXVEC(3,3)*BOXVEC(1,1)-BOXVEC(1,3)*BOXVEC(3,1)
       BOXVECIN(3,2)=BOXVEC(1,3)*BOXVEC(2,1)-BOXVEC(2,3)*BOXVEC(1,1)
       BOXVECIN(1,3)=BOXVEC(2,1)*BOXVEC(3,2)-BOXVEC(3,1)*BOXVEC(2,2)
       BOXVECIN(2,3)=BOXVEC(3,1)*BOXVEC(1,2)-BOXVEC(1,1)*BOXVEC(3,2)
       BOXVECIN(3,3)=BOXVEC(1,1)*BOXVEC(2,2)-BOXVEC(2,1)*BOXVEC(1,2)
       DET=DOT_PRODUCT(BOXVECIN(:,1),BOXVEC(:,1))
       BOXVECIN(:,:)=BOXVECIN(:,:)/DET
       BOXVECIN=TRANSPOSE(BOXVECIN)
!
!      ==================================================================
!      ==  CALCULATE #(ATOMS) IN THE MODEL                             ==
!      ==================================================================
       NATM=0
       DO IAT=1,NAT
         R(:)=POS(:,IAT)
         CALL BOXBOX(RBAS,BOXR0-R,BOXVEC,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
         MIN1=NINT(XMIN1)-1
         MAX1=NINT(XMAX1)+1
         MIN2=NINT(XMIN2)-1
         MAX2=NINT(XMAX2)+1
         MIN3=NINT(XMIN3)-1
         MAX3=NINT(XMAX3)+1
         DO I1=MIN1,MAX1
           T1=REAL(I1,KIND=8)
           DO I2=MIN2,MAX2
             T2=REAL(I2,KIND=8)
             DO I3=MIN3,MAX3
               T3=REAL(I3,KIND=8)
               R=POS(:,IAT)+RBAS(:,1)*T1+RBAS(:,2)*T2+RBAS(:,3)*T3-BOXR0(:)
               R=MATMUL(BOXVECIN,R)
               TCHK=(R(1).GE.0.D0.AND.R(1).LE.1.D0).AND. &
                    (R(2).GE.0.D0.AND.R(2).LE.1.D0).AND. &
                    (R(3).GE.0.D0.AND.R(3).LE.1.D0)
               IF(TCHK) NATM=NATM+1
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       RETURN
       END
!
!      ..................................................................
       SUBROUTINE MODEL$ATOMS(BOXR0,BOXVEC,NAT,POS,RBAS,NATM,POSM,MAP)
!      ******************************************************************
!      **                                                              **
!      ******************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: BOXR0(3)
       REAL(8)   ,INTENT(IN) :: BOXVEC(3,3)
       INTEGER(4),INTENT(IN) :: NAT
       REAL(8)   ,INTENT(IN) :: POS(3,NAT)
       REAL(8)   ,INTENT(IN) :: RBAS(3,3)
       INTEGER(4),INTENT(IN) :: NATM
       REAL(8)   ,INTENT(OUT):: POSM(3,NATM)
       INTEGER(4),INTENT(OUT):: MAP(NATM)
       INTEGER(4)            :: IAT,IATM,I1,I2,I3
       REAL(8)               :: BOXVECIN(3,3)
       REAL(8)               :: DET
       REAL(8)               :: R(3),XR(3)
       REAL(8)               :: XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3
       INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
       REAL(8)               :: T1,T2,T3
       LOGICAL(4)            :: TCHK
!      ******************************************************************
!
!      ==================================================================
!      ==  INVERSE OF BOXVEC                                           ==
!      ==================================================================
       BOXVECIN(1,1)=BOXVEC(2,2)*BOXVEC(3,3)-BOXVEC(3,2)*BOXVEC(2,3)
       BOXVECIN(2,1)=BOXVEC(3,2)*BOXVEC(1,3)-BOXVEC(1,2)*BOXVEC(3,3)
       BOXVECIN(3,1)=BOXVEC(1,2)*BOXVEC(2,3)-BOXVEC(2,2)*BOXVEC(1,3)
       BOXVECIN(1,2)=BOXVEC(2,3)*BOXVEC(3,1)-BOXVEC(3,3)*BOXVEC(2,1)
       BOXVECIN(2,2)=BOXVEC(3,3)*BOXVEC(1,1)-BOXVEC(1,3)*BOXVEC(3,1)
       BOXVECIN(3,2)=BOXVEC(1,3)*BOXVEC(2,1)-BOXVEC(2,3)*BOXVEC(1,1)
       BOXVECIN(1,3)=BOXVEC(2,1)*BOXVEC(3,2)-BOXVEC(3,1)*BOXVEC(2,2)
       BOXVECIN(2,3)=BOXVEC(3,1)*BOXVEC(1,2)-BOXVEC(1,1)*BOXVEC(3,2)
       BOXVECIN(3,3)=BOXVEC(1,1)*BOXVEC(2,2)-BOXVEC(2,1)*BOXVEC(1,2)
       DET=DOT_PRODUCT(BOXVECIN(:,1),BOXVEC(:,1))
       BOXVECIN(:,:)=BOXVECIN(:,:)/DET
       BOXVECIN=TRANSPOSE(BOXVECIN)
!
!      ==================================================================
!      ==  CALCULATE #(ATOMS) IN THE MODEL                             ==
!      ==================================================================
       IATM=0
       DO IAT=1,NAT
         R(:)=POS(:,IAT)
         CALL BOXBOX(RBAS,BOXR0-R,BOXVEC,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
         MIN1=NINT(XMIN1)-1
         MAX1=NINT(XMAX1)+1
         MIN2=NINT(XMIN2)-1
         MAX2=NINT(XMAX2)+1
         MIN3=NINT(XMIN3)-1
         MAX3=NINT(XMAX3)+1
         DO I1=MIN1,MAX1
           T1=REAL(I1,KIND=8)
           DO I2=MIN2,MAX2
             T2=REAL(I2,KIND=8)
             DO I3=MIN3,MAX3
               T3=REAL(I3,KIND=8)
               R=POS(:,IAT)+RBAS(:,1)*T1+RBAS(:,2)*T2+RBAS(:,3)*T3
               XR=MATMUL(BOXVECIN,R-BOXR0)
               TCHK=(XR(1).GE.0.D0.AND.XR(1).LE.1.D0).AND. &
                    (XR(2).GE.0.D0.AND.XR(2).LE.1.D0).AND. &
                    (XR(3).GE.0.D0.AND.XR(3).LE.1.D0)
               IF(TCHK) THEN
                 IATM=IATM+1
                 POSM(:,IATM)=R(:)
                 MAP(IATM)=IAT
               END IF
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       RETURN
       END
!
!      ..................................................................
       SUBROUTINE MODEL$NBONDM(NAT,POS,RAD,NBOND)
!      ******************************************************************
!      **                                                              **
!      ******************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NAT
       REAL(8)   ,INTENT(IN) :: POS(3,NAT)
       REAL(8)   ,INTENT(IN) :: RAD(NAT)
       INTEGER(4),INTENT(OUT):: NBOND
       INTEGER(4)            :: IAT1,IAT2
       REAL(8)               :: DR(3)
       REAL(8)               :: DIS
!      ******************************************************************
       NBOND=0
       DO IAT1=1,NAT
         DO IAT2=IAT1+1,NAT
           DR(:)=POS(:,IAT1)-POS(:,IAT2)
           DIS=SQRT(DR(1)**2+DR(2)**2+DR(3)**2)
           IF(DIS.LT.RAD(IAT1)+RAD(IAT2)) NBOND=NBOND+1
         ENDDO
       ENDDO
       RETURN
       END
!
!      ..................................................................
       SUBROUTINE MODEL$BONDS(NAT,POS,RAD,NBOND,BOND)
!      ******************************************************************
!      **                                                              **
!      ******************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NAT
       REAL(8)   ,INTENT(IN) :: POS(3,NAT)
       REAL(8)   ,INTENT(IN) :: RAD(NAT)
       INTEGER(4),INTENT(IN) :: NBOND
       INTEGER(4),INTENT(OUT):: BOND(2,NBOND)
       INTEGER(4)            :: IAT1,IAT2,IBOND
       REAL(8)               :: DR(3)
       REAL(8)               :: DIS
!      ******************************************************************
       IBOND=0
       DO IAT1=1,NAT
         DO IAT2=IAT1+1,NAT
           DR(:)=POS(:,IAT1)-POS(:,IAT2)
           DIS=SQRT(DOT_PRODUCT(DR,DR))
           IF(DIS.LT.RAD(IAT1)+RAD(IAT2)) THEN
             IBOND=IBOND+1
             BOND(1,IBOND)=IAT1
             BOND(2,IBOND)=IAT2
           END IF
         ENDDO
       ENDDO
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READWAVEPLOTDIM(NFIL,TITLE,NAT,NR1,NR2,NR3,TC)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)   :: NFIL
      CHARACTER(*) ,INTENT(OUT ) :: TITLE
      INTEGER(4)   ,INTENT(OUT)  :: NAT
      LOGICAL(4)   ,INTENT(OUT)  :: TC
      INTEGER(4)   ,INTENT(OUT)  :: NR1
      INTEGER(4)   ,INTENT(OUT)  :: NR2
      INTEGER(4)   ,INTENT(OUT)  :: NR3
      REAL(8)                    :: RBAS(3,3)
      CHARACTER(8)               :: BEGINID
      INTEGER(4)                 :: LENTITLE
!     ******************************************************************
      REWIND NFIL
      READ(NFIL)BEGINID,LENTITLE
      TC=(BEGINID.EQ.'CWAVEPLO')
      LENTITLE=MIN(LENTITLE,LEN(TITLE))
      READ(NFIL)TITLE(1:LENTITLE)
      READ(NFIL)RBAS,NAT
      READ(NFIL)NR1,NR2,NR3
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,NAME,XK,NR1,NR2,NR3,WAVE)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)   :: NFIL
      CHARACTER(*) ,INTENT(OUT ) :: TITLE
      INTEGER(4)   ,INTENT(IN)   :: NAT
      REAL(8)      ,INTENT(OUT)  :: RBAS(3,3)
      REAL(8)      ,INTENT(OUT)  :: POS(3,NAT)
      REAL(8)      ,INTENT(OUT)  :: Z(NAT)
      REAL(8)      ,INTENT(OUT)  :: XK(3)
      CHARACTER(32),INTENT(OUT)  :: NAME(NAT)
      INTEGER(4)   ,INTENT(IN)   :: NR1
      INTEGER(4)   ,INTENT(IN)   :: NR2
      INTEGER(4)   ,INTENT(IN)   :: NR3
      REAL(8)      ,INTENT(OUT)  :: WAVE(NR1,NR2,NR3)
      REAL(8)                    :: Q(NAT) !POINT CHARGES (NOT YET USED)
      INTEGER(4)                 :: NR1_,NR2_,NR3_,NAT_
      CHARACTER(8)               :: BEGINID
      INTEGER(4)                 :: LENTITLE
      CHARACTER(11)              :: ENDID
!     ******************************************************************
      REWIND NFIL
      READ(NFIL)BEGINID,LENTITLE
      IF(BEGINID.NE.'WAVEPLOT') THEN
        CALL ERROR$MSG('INCORRECT FILE TYPE. MUST BE "WAVEPLOT".')
        CALL ERROR$CHVAL('BEGINID',BEGINID)
        CALL ERROR$STOP('READWAVEPLOT')
      END IF
      LENTITLE=MIN(LENTITLE,LEN(TITLE))
      READ(NFIL)TITLE(1:LENTITLE)
      READ(NFIL)RBAS,NAT_
      READ(NFIL)NR1_,NR2_,NR3_
      IF(NAT_.NE.NAT.OR.NR1_.NE.NR1.OR.NR2_.NE.NR2.OR.NR3_.NE.NR3) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$I4VAL('NAT ON FILE ',NAT_)
        CALL ERROR$I4VAL('NAT ON INPUT',NAT)
        CALL ERROR$I4VAL('NR1 ON FILE ',NR1_)
        CALL ERROR$I4VAL('NR1 ON INPUT',NR1)
        CALL ERROR$I4VAL('NR2 ON FILE ',NR2_)
        CALL ERROR$I4VAL('NR2 ON INPUT',NR2)
        CALL ERROR$I4VAL('NR3 ON FILE ',NR3_)
        CALL ERROR$I4VAL('NR3 ON INPUT',NR3)
        CALL ERROR$STOP('READWAVEPLOT')
      END IF
      IF(NAT.NE.0) THEN
        READ(NFIL)NAME
        READ(NFIL)Z
        READ(NFIL)POS
        READ(NFIL)Q
      END IF
      XK(:)=0.D0
      READ(NFIL)WAVE
      READ(NFIL)ENDID 
      IF(ENDID.NE.'END OF FILE') THEN
        CALL ERROR$MSG('DATA FILE CORRUPTED')
        CALL ERROR$CHVAL('ENDID',ENDID)
        CALL ERROR$STOP('READWAVEPLOT')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READWAVEPLOTC(NFIL,TITLE,RBAS,NAT,POS,Z,NAME,XK,NR1,NR2,NR3,WAVE)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)   :: NFIL
      CHARACTER(*) ,INTENT(OUT ) :: TITLE
      INTEGER(4)   ,INTENT(IN)   :: NAT
      REAL(8)      ,INTENT(OUT)  :: RBAS(3,3)
      REAL(8)      ,INTENT(OUT)  :: POS(3,NAT)
      REAL(8)      ,INTENT(OUT)  :: Z(NAT)
      REAL(8)      ,INTENT(OUT)  :: XK(3)
      CHARACTER(32),INTENT(OUT)  :: NAME(NAT)
      INTEGER(4)   ,INTENT(IN)   :: NR1
      INTEGER(4)   ,INTENT(IN)   :: NR2
      INTEGER(4)   ,INTENT(IN)   :: NR3
      COMPLEX(8)   ,INTENT(OUT)  :: WAVE(NR1,NR2,NR3)
      REAL(8)                    :: Q(NAT) !POINT CHARGES (NOT YET USED)
      INTEGER(4)                 :: NR1_,NR2_,NR3_,NAT_
      CHARACTER(8)               :: BEGINID
      INTEGER(4)                 :: LENTITLE
      CHARACTER(11)              :: ENDID
!     ******************************************************************
      REWIND NFIL
      READ(NFIL)BEGINID,LENTITLE
      IF(BEGINID.NE.'CWAVEPLO') THEN
        CALL ERROR$MSG('INCORRECT FILE TYPE. MUST BE "CWAVEPLO".')
        CALL ERROR$CHVAL('BEGINID',BEGINID)
        CALL ERROR$STOP('READWAVEPLOTC')
      END IF
      LENTITLE=MIN(LENTITLE,LEN(TITLE))
      READ(NFIL)TITLE(1:LENTITLE)
      READ(NFIL)RBAS,NAT_
      READ(NFIL)NR1_,NR2_,NR3_
      IF(NAT_.NE.NAT.OR.NR1_.NE.NR1.OR.NR2_.NE.NR2.OR.NR3_.NE.NR3) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$I4VAL('NAT ON FILE ',NAT_)
        CALL ERROR$I4VAL('NAT ON INPUT',NAT)
        CALL ERROR$I4VAL('NR1 ON FILE ',NR1_)
        CALL ERROR$I4VAL('NR1 ON INPUT',NR1)
        CALL ERROR$I4VAL('NR2 ON FILE ',NR2_)
        CALL ERROR$I4VAL('NR2 ON INPUT',NR2)
        CALL ERROR$I4VAL('NR3 ON FILE ',NR3_)
        CALL ERROR$I4VAL('NR3 ON INPUT',NR3)
        CALL ERROR$STOP('READWAVEPLOTC')
      END IF
      IF(NAT.NE.0) THEN
        READ(NFIL)NAME
        READ(NFIL)Z
        READ(NFIL)POS
        READ(NFIL)Q
      END IF
      READ(NFIL)XK
      READ(NFIL)WAVE
      READ(NFIL)ENDID 
      IF(ENDID.NE.'END OF FILE') THEN
        CALL ERROR$MSG('DATA FILE CORRUPTED')
        CALL ERROR$CHVAL('ENDID',ENDID)
        CALL ERROR$STOP('READWAVEPLOTC')
      END IF
      RETURN
      END
!                                                                       
!     .....................................................BONDS .......
      SUBROUTINE DX_DENSITY(NFIL,NR1,NR2,NR3,RBAS,DENSITYIN,BOXR0,BOXVEC)
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      INTEGER(4),INTENT(IN) :: NR3
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: DENSITYIN(NR1,NR2,NR3)
      REAL(8)   ,INTENT(IN) :: BOXR0(3)
      REAL(8)   ,INTENT(IN) :: BOXVEC(3,3)
      REAL(8)               :: SBAS(3,3)
      REAL(8)               :: XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3
      INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: I1,I2,I3
      REAL(8)               :: T1,T2,T3
      REAL(8)               :: X,Y,Z
      REAL(8)               :: SH1,SH2,SH3
      INTEGER(4)            :: I,J,K
      REAL(8)  ,PARAMETER   :: MAXVALUE=9999.D0
      REAL(8)               :: DENSITY(NR1,NR2,NR3)
      REAL(8)               :: SVAR
!     *******************************************************************
      DENSITY(:,:,:)=DENSITYIN(:,:,:)
      DO K=1,NR3
        DO J=1,NR2
          DO I=1,NR1
            SVAR=MIN(MAXVALUE,DENSITY(I,J,K))
            DENSITY(I,J,K)=MAX(-MAXVALUE,SVAR)
          END DO
        END DO
      END DO
      SBAS(:,1)=RBAS(:,1)/REAL(NR1,KIND=8)
      SBAS(:,2)=RBAS(:,2)/REAL(NR2,KIND=8)
      SBAS(:,3)=RBAS(:,3)/REAL(NR3,KIND=8)
      CALL BOXBOX(SBAS,BOXR0,BOXVEC,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
      MIN1=INT(ABS(XMIN1+1.D0))*NINT(SIGN(1.D0,XMIN1))
      MAX1=INT(ABS(XMAX1+1.D0))*NINT(SIGN(1.D0,XMAX1))
      MIN2=INT(ABS(XMIN2+1.D0))*NINT(SIGN(1.D0,XMIN2))
      MAX2=INT(ABS(XMAX2+1.D0))*NINT(SIGN(1.D0,XMAX2))
      MIN3=INT(ABS(XMIN3+1.D0))*NINT(SIGN(1.D0,XMIN3))
      MAX3=INT(ABS(XMAX3+1.D0))*NINT(SIGN(1.D0,XMAX3))
      N1=MAX1-MIN1+1
      N2=MAX2-MIN2+1
      N3=MAX3-MIN3+1
!
!     ==================================================================
!     ==================================================================
!     ==   WRITE DX FILE                                              ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     ==   DATA ARRAY: DENSITY                                        ==
!     ==================================================================
      WRITE(NFIL,FMT='(A/5X,A,5X,A,5X,A,5X,A,I10/A)') &
     &                -"OBJECT 1" &
     &               ,-"CLASS ARRAY",-"TYPE FLOAT",-"RANK 0",-"ITEMS ",N1*N2*N3 &
     &               ,-"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10F12.5)') &
     &     (((DENSITY(MOD(I+1000*NR1,NR1)+1 &
     &               ,MOD(J+1000*NR2,NR2)+1 &
     &               ,MOD(K+1000*NR3,NR3)+1) &
     &               ,K=MIN3,MAX3) &
     &               ,J=MIN2,MAX2) &
     &               ,I=MIN1,MAX1)
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   POSITIONS ARRAY: GRID POSITIONS                            ==
!     ==================================================================
      WRITE(NFIL,FMT='(A/A,3I10)') &
     &               -"OBJECT 2" &
     &              ,-"CLASS GRIDPOSITIONS COUNTS ",N1,N2,N3
      SH1=SBAS(1,1)*MIN1+SBAS(1,2)*MIN2+SBAS(1,3)*MIN3
      SH2=SBAS(2,1)*MIN1+SBAS(2,2)*MIN2+SBAS(2,3)*MIN3
      SH3=SBAS(3,1)*MIN1+SBAS(3,2)*MIN2+SBAS(3,3)*MIN3
      WRITE(NFIL,FMT='(A,3F10.5)')-"ORIGIN",SH1,SH2,SH3
      WRITE(NFIL,FMT='(A,3F10.5)')-'DELTA',SBAS(:,1)
      WRITE(NFIL,FMT='(A,3F10.5)')-'DELTA',SBAS(:,2)
      WRITE(NFIL,FMT='(A,3F10.5)')-'DELTA',SBAS(:,3)
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   CONNECTIONS ARRAY: CONNECTIONS                             ==
!     ==================================================================
      WRITE(NFIL,FMT='(A/A,3I10)') &
     &              -"OBJECT 3" &
     &             ,-"CLASS GRIDCONNECTIONS COUNTS ",N1,N2,N3
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""ELEMENT TYPE"" STRING ""CUBES"""
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""REF"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!
!     ================================================================
!     ==   CLIPBOX                                                  ==
!     ================================================================
      WRITE(NFIL,FMT='(A,I10/A/A)') &
     &             -"OBJECT ",4 &
     &            ,-"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS 8" &
     &            ,-"DATA FOLLOWS"
      DO I3=0,1 
         T3=DBLE(I3)
        DO I2=0,1
          T2=DBLE(I2)
          DO I1=0,1
            T1=DBLE(I1)
            X=BOXR0(1)+BOXVEC(1,1)*T1+BOXVEC(1,2)*T2+BOXVEC(1,3)*T3
            Y=BOXR0(2)+BOXVEC(2,1)*T1+BOXVEC(2,2)*T2+BOXVEC(2,3)*T3
            Z=BOXR0(3)+BOXVEC(3,1)*T1+BOXVEC(3,2)*T2+BOXVEC(3,3)*T3
            WRITE(NFIL,FMT='(3F10.5)')X,Y,Z
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   OBJECT MOLECULE:                                           ==
!     ==================================================================
      WRITE(NFIL,FMT='(A)')-"OBJECT ""DENSITY"""
      WRITE(NFIL,FMT='(A)')-"CLASS FIELD"
      WRITE(NFIL,FMT='(A)')-"COMPONENT ""DATA"" VALUE 1"
      WRITE(NFIL,FMT='(A)')-"COMPONENT ""POSITIONS"" VALUE 2"
      WRITE(NFIL,FMT='(A)')-"COMPONENT ""CONNECTIONS"" VALUE 3"
      WRITE(NFIL,FMT='(A)')-"COMPONENT ""BOX"" VALUE 4"
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""NAME"" STRING ""DENSITY"""
      WRITE(NFIL,FMT='(A)')-"#"
!
!     ==================================================================
!     ==   END                                                        ==
!     ==================================================================
!     WRITE(NFIL,FMT='("END")')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE DXBALLSTICK(NFIL,NAT,COLOR,RAD,POS,NBOND,IBOND,BOXR0,BOXVEC)
!     **                                                              **
!     **  WRITES A DATAEXPLORER FILE FOR A BALLSTICK MODEL            **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: COLOR(3,NAT)
      REAL(8)   ,INTENT(IN) :: RAD(NAT)
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      INTEGER(4),INTENT(IN) :: NBOND
      INTEGER(4),INTENT(IN) :: IBOND(2,NBOND)
      REAL(8)   ,INTENT(IN) :: BOXR0(3)
      REAL(8)   ,INTENT(IN) :: BOXVEC(3,3)
      INTEGER(4)            :: IOBJECT0=4
      REAL(8)               :: X,Y,Z
      REAL(8)               :: T1,T2,T3
      INTEGER(4)            :: I1,I2,I3
!     ******************************************************************
!     
!     ==================================================================
!     ==   DATA ARRAY: SIZE OF THE SPHERES                            ==
!     ==================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"SPHERE SIZE"/"#")') !COMMENT
      WRITE(NFIL,FMT='(A,I10,/5X,A,5X,A,5X,A,5X,A,I10/A)') &
     &                -"OBJECT ",IOBJECT0+1 &
     &               ,-"CLASS ARRAY",-"TYPE FLOAT",-"RANK 0",-"ITEMS ",NAT &
     &               ,-"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10F10.5)')RAD(:)
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   COLORS                                                   ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"COLORS"/"#")') !COMMENT
      WRITE(NFIL,FMT='(A,I10/A,I10/A)') &
     &             -"OBJECT ",IOBJECT0+2 &
     &            ,-"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",NAT &
     &            ,-"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10F10.5)')COLOR(:,:)
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   POSITIONS ARRAY: ATOMIC POSITIONS                        ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"ATOMIC POSITIONS"/"#")') !COMMENT
      WRITE(NFIL,FMT='(A,I10/A,I10/A)') &
     &            -"OBJECT ",IOBJECT0+3 &
     &           ,-"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",NAT &
     &           ,-"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10F10.5)')POS(:,:)
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   CONNECTIONS ARRAY: BONDS                                 ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"BONDS"/"#")') !COMMENT
      IF(NBOND.EQ.0) THEN
        WRITE(NFIL,FMT='(A,I10/A,I10)') &
     &             -"OBJECT ",IOBJECT0+4 &
     &            ,-"CLASS ARRAY TYPE INT RANK 1 SHAPE 2 ITEMS ",NBOND 
      ELSE
        WRITE(NFIL,FMT='(A,I10/A,I10/A)') &
     &             -"OBJECT ",IOBJECT0+4 &
     &            ,-"CLASS ARRAY TYPE INT RANK 1 SHAPE 2 ITEMS ",NBOND &
     &            ,-"DATA FOLLOWS"
      END IF
      WRITE(NFIL,FMT='(10I5)')IBOND(:,:)-1
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""REF"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""ELEMENT TYPE"" STRING ""LINES"""
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   BOX: LATTICE VECTORS                                     ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"LATTICE VECTORS"/"#")')
      WRITE(NFIL,FMT='(A,I10/A/A)') &
     &            -"OBJECT ",IOBJECT0+5 &
     &           ,-"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS 8" &
     &           ,-"DATA FOLLOWS"
      DO I3=0,1 
         T3=DBLE(I3)
        DO I2=0,1
          T2=DBLE(I2)
          DO I1=0,1
            T1=DBLE(I1)
            X=BOXR0(1)+BOXVEC(1,1)*T1+BOXVEC(1,2)*T2+BOXVEC(1,3)*T3
            Y=BOXR0(2)+BOXVEC(2,1)*T1+BOXVEC(2,2)*T2+BOXVEC(2,3)*T3
            Z=BOXR0(3)+BOXVEC(3,1)*T1+BOXVEC(3,2)*T2+BOXVEC(3,3)*T3
            WRITE(NFIL,FMT='(3F10.5)')X,Y,Z
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   OBJECT MOLECULE:                                         ==
!     ================================================================
      WRITE(NFIL,FMT='(A)')-"OBJECT ""BALLSTICK"""
      WRITE(NFIL,FMT='(A)')-"CLASS FIELD"
      WRITE(NFIL,FMT='(A,I10)')-"COMPONENT ""DATA"" VALUE ",IOBJECT0+1
      WRITE(NFIL,FMT='(A,I10)')-"COMPONENT ""POSITIONS"" VALUE ",IOBJECT0+3
      WRITE(NFIL,FMT='(A,I10)')-"COMPONENT ""CONNECTIONS"" VALUE ",IOBJECT0+4
      WRITE(NFIL,FMT='(A,I10)')-"COMPONENT ""BOX"" VALUE ",IOBJECT0+5
      WRITE(NFIL,FMT='(A,I10)')-"COMPONENT ""COLORS"" VALUE ",IOBJECT0+2
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""NAME"" STRING ""ATOMS"""
      WRITE(NFIL,FMT='("#")')
      RETURN
      END
!        
!     ..................................................................
      SUBROUTINE ATOMCOLOR(IZ,ICOLOR)
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IZ
      INTEGER(4),INTENT(OUT) :: ICOLOR(3)
      INTEGER(4)             :: ICOLORSTANDARD(3,106)
      INTEGER(4)             :: I
!     ******************************************************************
      DATA (ICOLORSTANDARD(I,  1),I=1,3)/135,125,131/
      DATA (ICOLORSTANDARD(I,  2),I=1,3)/255,228,196/!BISQUE
      DATA (ICOLORSTANDARD(I,  3),I=1,3)/240,248,255/!ALICE BLUE
      DATA (ICOLORSTANDARD(I,  4),I=1,3)/200,200,  0/!YELLOW
      DATA (ICOLORSTANDARD(I,  5),I=1,3)/192,102,624/!CORAL
      DATA (ICOLORSTANDARD(I,  6),I=1,3)/ 50, 50, 50/
      DATA (ICOLORSTANDARD(I,  7),I=1,3)/ 10,200, 10/
      DATA (ICOLORSTANDARD(I,  8),I=1,3)/255,  0,  0/!RED
      DATA (ICOLORSTANDARD(I,  9),I=1,3)/ 50,205, 50/!LIME GREEN
      DATA (ICOLORSTANDARD(I, 10),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 11),I=1,3)/ 60,  1,  1/!REDISH BLACK
      DATA (ICOLORSTANDARD(I, 12),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 13),I=1,3)/ 91,189,255/!LIGHT BLUE
      DATA (ICOLORSTANDARD(I, 14),I=1,3)/135,125,131/!DARK BLUE
      DATA (ICOLORSTANDARD(I, 15),I=1,3)/230,171, 17/
      DATA (ICOLORSTANDARD(I, 16),I=1,3)/240,240,  0/
      DATA (ICOLORSTANDARD(I, 17),I=1,3)/ 60,180,  0/!YELLOW GREEN
      DATA (ICOLORSTANDARD(I, 18),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 19),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 20),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 21),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 22),I=1,3)/ 50,255, 50/!LIME GREEN
      DATA (ICOLORSTANDARD(I, 23),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 24),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 25),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 26),I=1,3)/176, 48, 96/
      DATA (ICOLORSTANDARD(I, 27),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 28),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 29),I=1,3)/278,134, 34/!FIREBRICK
      DATA (ICOLORSTANDARD(I, 30),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 31),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 32),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 33),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 34),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 35),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 36),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 37),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 38),I=1,3)/239,255, 18/
      DATA (ICOLORSTANDARD(I, 39),I=1,3)/ 40, 40,140/
      DATA (ICOLORSTANDARD(I, 40),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 41),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 42),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 43),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 44),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 45),I=1,3)/230, 51, 41/
      DATA (ICOLORSTANDARD(I, 46),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 47),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 48),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 49),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 50),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 51),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 52),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 53),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 54),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 55),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 56),I=1,3)/100, 40,  0/
      DATA (ICOLORSTANDARD(I, 57),I=1,3)/255,125,255/ !PURPLE
      DATA (ICOLORSTANDARD(I, 58),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 59),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 60),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 61),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 62),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 63),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 64),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 65),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 66),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 67),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 68),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 69),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 70),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 71),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 72),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 73),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 74),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 75),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 76),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 77),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 78),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 79),I=1,3)/255,215,  0/!GOLD
      DATA (ICOLORSTANDARD(I, 80),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 81),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 82),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 83),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 84),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 85),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 86),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 87),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 88),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 89),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 90),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 91),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 92),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 93),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 94),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 95),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 96),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 97),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 98),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 99),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,100),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,101),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,102),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,103),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,104),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,105),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,106),I=1,3)/  1,  1,  1/   
      DO I=1,3
        ICOLOR(I)=ICOLORSTANDARD(I,IZ)
      ENDDO
      RETURN
      END
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITECMCV(NFIL,TITLE,RBAS,N1,N2,N3,FIELD)
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: TITLE
      REAL(8)     ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4)  ,INTENT(IN) :: N1,N2,N3
      REAL(8)     ,INTENT(IN) :: FIELD(N1,N2,N3)
      INTEGER(4)              :: I,J,K
!     **************************************************************************
      REWIND(NFIL)
      WRITE(NFIL,*)TITLE
      WRITE(NFIL,*)N1,N2,N3
      WRITE(NFIL,*)RBAS
      DO K=1,N3
        DO J=1,N2
          DO I=1,N1
            WRITE(NFIL,*)FIELD(I,J,K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKECUBE(NFIL,NAT,Z,R,RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     ** IT THEN PLACES THE DATA INTO A FILE WITH THE GAUSSIAN CUBE FORMAT.   **
!     ** A VERY CLEAR DESCRIPTION OF THE GAUSSIAN CUBE FORMAT HAS BEEN GIVEN  **
!     **  ON HTTP://LOCAL.WASP.UWA.EDU.AU/~PBOURKE/DATAFORMATS/CUBE/          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) !LATTICE VECTORS OF THE PERIODIC DATA
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3 ! # GRID POINTS ON THE PERIODIC GRID
      REAL(8)   ,INTENT(IN) :: WAVE(NR1,NR2,NR3) ! PERIODIC DENSITY DATA
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! CORNER OF THE NEW GRID
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! VECTORS SPANNING THE NEW GRID
      INTEGER(4),PARAMETER  :: N1=80      ! DISPLACEMENT
      INTEGER(4),PARAMETER  :: N2=80
      INTEGER(4),PARAMETER  :: N3=80
      REAL(8)   ,ALLOCATABLE:: DATA(:,:,:) !(N1,N2,N3)
      REAL(8)               :: RTOX(3,3)
      REAL(8)               :: VOL
      INTEGER(4)            :: NRVIEW,NRVIEWX
      REAL(8)   ,ALLOCATABLE:: RVIEW(:,:)
      REAL(8)   ,ALLOCATABLE:: ZVIEW(:)
      REAL(8)               :: PI
!     **************************************************************************
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE FUNCTION ONTO VIEWING GRID                         ==
!     ==========================================================================
      ALLOCATE(DATA(N1,N2,N3))
      CALL MAPFIELDTOGRID(RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX,N1,N2,N3,DATA)
 !
!     ==========================================================================
!     ==  MAP ATOMS INTO VIEW BOX                                             ==
!     ==========================================================================
      PI=4.D0*ATAN(1.D0)
      CALL GBASS(BOX,RTOX,VOL)      
      NRVIEWX=2*NINT(VOL/(4.D0*PI/3.D0*1.4D0**3)) ! ESTIMATE MAX NUMBER OF ATOMS
                                                ! INSIDE THE VIEWBOX
      ALLOCATE(RVIEW(3,NRVIEWX)) ! POSITIONS
      ALLOCATE(ZVIEW(NRVIEWX))   ! ATOMIC NUMBERS
      CALL MAKEATOMSINTOBOX(NAT,Z,R,RBAS,ORIGIN,BOX,NRVIEWX,NRVIEW,RVIEW,ZVIEW)
!
!     ==========================================================================
!     ==  WRITE CUBE FILE                                                     ==
!     ==========================================================================
      CALL WRITECUBEFILE(NFIL,NRVIEW,ZVIEW,RVIEW,ORIGIN,BOX,N1,N2,N3,DATA)
      DEALLOCATE(DATA)
!
      DEALLOCATE(RVIEW)
      DEALLOCATE(ZVIEW)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEVRML(NFIL,NAT,Z,R,RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX &
     &                   ,TPLANE,PLANER0,PLANEVEC)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     ** IT THEN PLACES THE DATA INTO A FILE WITH THE GAUSSIAN CUBE FORMAT.   **
!     ** A VERY CLEAR DESCRIPTION OF THE GAUSSIAN CUBE FORMAT HAS BEEN GIVEN  **
!     **  ON HTTP://LOCAL.WASP.UWA.EDU.AU/~PBOURKE/DATAFORMATS/CUBE/          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) !LATTICE VECTORS OF THE PERIODIC DATA
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3 ! # GRID POINTS ON THE PERIODIC GRID
      REAL(8)   ,INTENT(IN) :: WAVE(NR1,NR2,NR3) ! PERIODIC DENSITY DATA
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! CORNER OF THE NEW GRID
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! VECTORS SPANNING THE NEW GRID
      LOGICAL(4),INTENT(IN) :: TPLANE
      REAL(8)   ,INTENT(IN) :: PLANER0(3)   ! 
      REAL(8)   ,INTENT(IN) :: PLANEVEC(3,2)  ! 
      INTEGER(4),PARAMETER  :: N1=160      ! DISPLACEMENT
      INTEGER(4),PARAMETER  :: N2=160
      INTEGER(4),PARAMETER  :: N3=1
      REAL(8)   ,ALLOCATABLE:: DATA(:,:,:) !(N1,N2,N3) TOO LARGE FOR STACK
      REAL(8)               :: RTOX(3,3)
      REAL(8)               :: VOL
      INTEGER(4)            :: NRVIEW,NRVIEWX
      REAL(8)   ,ALLOCATABLE:: RVIEW(:,:)
      REAL(8)   ,ALLOCATABLE:: ZVIEW(:)
      REAL(8)               :: PI
      REAL(8)               :: PLANEBOX(3,3)
!     **************************************************************************
                          CALL TRACE$PUSH('MAKEVRML')
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE FUNCTION ONTO VIEWING GRID                         ==
!     ==========================================================================
      IF(TPLANE) THEN
        ALLOCATE(DATA(N1,N2,N3))
        PLANEBOX(:,:)=0.D0
        PLANEBOX(:,:2)=PLANEVEC(:,:)
        CALL MAPFIELDTOGRID(RBAS,NR1,NR2,NR3,WAVE &
     &                                     ,PLANER0,PLANEBOX,N1,N2,N3,DATA)
      END IF
!
!     ==========================================================================
!     ==  MAP ATOMS INTO VIEW BOX                                             ==
!     ==========================================================================
      PI=4.D0*ATAN(1.D0)
      CALL GBASS(BOX,RTOX,VOL)      
      NRVIEWX=2*NINT(VOL/(4.D0*PI/3.D0*1.4D0**3)) ! ESTIMATE MAX NUMBER OF ATOMS
                                                ! INSIDE THE VIEWBOX
      ALLOCATE(RVIEW(3,NRVIEWX)) ! POSITIONS
      ALLOCATE(ZVIEW(NRVIEWX))   ! ATOMIC NUMBERS
      CALL MAKEATOMSINTOBOX(NAT,Z,R,RBAS,ORIGIN,BOX,NRVIEWX,NRVIEW,RVIEW,ZVIEW)
!
!     ==========================================================================
!     ==  WRITE CUBE FILE                                                     ==
!     ==========================================================================
      WRITE(NFIL,FMT='("#VRML V2.0 UTF8")')
      WRITE(NFIL,*)'BACKGROUND {'
      WRITE(NFIL,*)'  SKYCOLOR 0.933 0.863 0.51'
      WRITE(NFIL,*)'  GROUNDCOLOR 0.933 0.863 0.51'
      WRITE(NFIL,*)'}'
      WRITE(NFIL,*)'VIEWPOINT {'
      WRITE(NFIL,*)'  POSITION 0 0 40 '
      WRITE(NFIL,*)'  FIELDOFVIEW 0.8 '
      WRITE(NFIL,*)'} # CLOSE VIEWPOINT'
      CALL VRML$ADDBALLSTICK(NFIL,NRVIEW,ZVIEW,RVIEW)
      IF(TPLANE) THEN
        CALL VRML$ADDRUBBERSHEET(NFIL,PLANER0,PLANEBOX,N1,N2,DATA)
        DEALLOCATE(DATA)
      END IF
!
      DEALLOCATE(RVIEW)
      DEALLOCATE(ZVIEW)

                          CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEVRMLC(NFIL,NAT,Z,R,RBAS,XK,NR1,NR2,NR3,CWAVE,ORIGIN,BOX &
     &                   ,TPLANE,PLANER0,PLANEVEC)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     ** IT THEN PLACES THE DATA INTO A FILE WITH THE GAUSSIAN CUBE FORMAT.   **
!     ** A VERY CLEAR DESCRIPTION OF THE GAUSSIAN CUBE FORMAT HAS BEEN GIVEN  **
!     **  ON HTTP://LOCAL.WASP.UWA.EDU.AU/~PBOURKE/DATAFORMATS/CUBE/          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    ! ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(IN) :: XK(3)     ! K-POINT IN RELATIVE COORDINATES
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS OF THE PERIODIC DATA
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3 ! # GRID POINTS ON THE PERIODIC GRID
      COMPLEX(8),INTENT(IN) :: CWAVE(NR1,NR2,NR3) ! PERIODIC DENSITY DATA
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! CORNER OF THE NEW GRID
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! VECTORS SPANNING THE NEW GRID
      LOGICAL(4),INTENT(IN) :: TPLANE
      REAL(8)   ,INTENT(IN) :: PLANER0(3) ! 
      REAL(8)   ,INTENT(IN) :: PLANEVEC(3,2)  ! 
      INTEGER(4),PARAMETER  :: N1=80      ! DISPLACEMENT
      INTEGER(4),PARAMETER  :: N2=80
      INTEGER(4),PARAMETER  :: N3=1
      COMPLEX(8),ALLOCATABLE:: CDATA(:,:,:) !(N1,N2,N3)
      REAL(8)   ,ALLOCATABLE:: DATA(:,:,:)  !(N1,N2,N3)
      REAL(8)               :: RTOX(3,3)
      REAL(8)               :: VOL
      INTEGER(4)            :: NRVIEW,NRVIEWX
      REAL(8)   ,ALLOCATABLE:: RVIEW(:,:)
      REAL(8)   ,ALLOCATABLE:: ZVIEW(:)
      REAL(8)               :: PI
      REAL(8)               :: PLANEBOX(3,3)
!     **************************************************************************
                          CALL TRACE$PUSH('MAKEVRML')
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE FUNCTION ONTO VIEWING GRID                         ==
!     ==========================================================================
      IF(TPLANE) THEN
        ALLOCATE(CDATA(N1,N2,N3))
        PLANEBOX(:,:)=0.D0
        PLANEBOX(:,:2)=PLANEVEC(:,:)
        CALL MAPFIELDTOGRIDC(RBAS,XK,NR1,NR2,NR3,CWAVE &
     &                                     ,PLANER0,PLANEBOX,N1,N2,N3,CDATA)
      END IF
!
!     ==========================================================================
!     ==  MAP ATOMS INTO VIEW BOX                                             ==
!     ==========================================================================
      PI=4.D0*ATAN(1.D0)
      CALL GBASS(BOX,RTOX,VOL)      
      NRVIEWX=2*NINT(VOL/(4.D0*PI/3.D0*1.4D0**3)) ! ESTIMATE MAX NUMBER OF ATOMS
                                                ! INSIDE THE VIEWBOX
      ALLOCATE(RVIEW(3,NRVIEWX)) ! POSITIONS
      ALLOCATE(ZVIEW(NRVIEWX))   ! ATOMIC NUMBERS
      CALL MAKEATOMSINTOBOX(NAT,Z,R,RBAS,ORIGIN,BOX,NRVIEWX,NRVIEW,RVIEW,ZVIEW)
!
!     ==========================================================================
!     ==  WRITE CUBE FILE                                                     ==
!     ==========================================================================
      WRITE(NFIL,FMT='("#VRML V2.0 UTF8")')
      WRITE(NFIL,*)'BACKGROUND {'
      WRITE(NFIL,*)'  SKYCOLOR 0.933 0.863 0.51'
      WRITE(NFIL,*)'  GROUNDCOLOR 0.933 0.863 0.51'
      WRITE(NFIL,*)'}'
      WRITE(NFIL,*)'VIEWPOINT {'
      WRITE(NFIL,*)'  POSITION 0 0 40 '
      WRITE(NFIL,*)'  FIELDOFVIEW 0.8 '
      WRITE(NFIL,*)'} # CLOSE VIEWPOINT'
      CALL VRML$ADDBALLSTICK(NFIL,NRVIEW,ZVIEW,RVIEW)
      IF(TPLANE) THEN
        ALLOCATE(DATA(N1,N2,N3))
        DATA=REAL(CDATA)
        CALL VRML$ADDRUBBERSHEET(NFIL,PLANER0,PLANEBOX,N1,N2,DATA)
        DEALLOCATE(CDATA)
        DEALLOCATE(DATA)
      END IF
!
      DEALLOCATE(RVIEW)
      DEALLOCATE(ZVIEW)

                          CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEGNU(NFIL,TYPE,NAT,Z,R,RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) !LATTICE VECTORS OF THE PERIODIC DATA
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3 ! # GRID POINTS ON THE PERIODIC GRID
      REAL(8)   ,INTENT(IN) :: WAVE(NR1,NR2,NR3) ! PERIODIC DENSITY DATA
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! CORNER OF THE NEW GRID
      REAL(8)   ,INTENT(IN) :: BOX(3,2)   ! VECTORS SPANNING THE NEW GRID
      INTEGER(4),PARAMETER  :: N1=60      ! DISPLACEMENT
      INTEGER(4),PARAMETER  :: N2=60
      INTEGER(4),PARAMETER  :: N3=1
      REAL(8)               :: BOX3D(3,3)
      REAL(8)               :: DATA(N1,N2,N3)
      REAL(8)               :: XMAX,XMIN,YMAX,YMIN
      INTEGER(4)            :: I,J
!     **************************************************************************
                              CALL TRACE$PUSH('MAKEGNU')
      IF(TYPE.EQ.'SURFACE') THEN
      ELSE IF(TYPE.EQ.'CONTOUR') THEN
      ELSE
        CALL ERROR$MSG('TYPE NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$MSG('PERMITTED VALUES ARE "CONTOUR" AND "SURFACE"')
        CALL ERROR$STOP('MAKEGNU')
      END IF
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE FUNCTION ONTO VIEWING GRID                         ==
!     ==========================================================================
      BOX3D=0.D0  ! ONLY NEEDED BECAUSE MAPFIELDTOGRID EXPECTS A 3D BOX
      BOX3D(:,:2)=BOX  
      CALL MAPFIELDTOGRID(RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX3D,N1,N2,N3,DATA)
!
      XMIN=-0.5D0*SQRT(SUM(BOX(:,1)**2)) 
      XMAX=+0.5D0*SQRT(SUM(BOX(:,1)**2)) 
      YMIN=-0.5D0*SQRT(SUM(BOX(:,2)**2)) 
      YMAX=+0.5D0*SQRT(SUM(BOX(:,2)**2)) 
 !
!     ==========================================================================
!     ==  WRITE GNUPLOT FILE RUBBERSHEET.GNU                                  ==
!     ==========================================================================
!     ==  REMARK "SET DATA STYLE LINES" IS NOT RECOGNIZED ANY MORE
!     ==  REMARK "SET DGRID3D  60,60,1" DOES NOT POPERLY INTERPOLATE
!     ==  Z-RANGE (DATA) REMARK "SET DGRID3D  60,60,1" DOES NOT POPERLY 
!     ==  INTERPOLATE
      WRITE(NFIL,*)'#'                                                     
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DATA SECTION TO BE CHANGED BY THE USER                 =='
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)-'XMIN=',XMIN
      WRITE(NFIL,*)-'XMAX=',XMAX
      WRITE(NFIL,*)-'YMIN=',YMIN
      WRITE(NFIL,*)-'YMAX=',YMAX
      WRITE(NFIL,*)-'ZMIN=',MINVAL(DATA)
      WRITE(NFIL,*)-'ZMAX=',MAXVAL(DATA)
      IF(TYPE.EQ.'SURFACE') THEN
        WRITE(NFIL,*)-'ROT_X=30'
        WRITE(NFIL,*)-'ROT_Z=20'
      ELSE IF(TYPE.EQ.'CONTOUR') THEN
        WRITE(NFIL,*)-'ROT_X=0'
        WRITE(NFIL,*)-'ROT_Z=0'
      END IF
      WRITE(NFIL,*)-'SCALE=1.8'
      WRITE(NFIL,*)-'SCALE_Z=1.'
      WRITE(NFIL,*)'#'                                                     
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DEFINE LINE STYLES TO BE USED WITH LS IN SPLOT         =='
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DEFINE A LINESTYLE FOR SPLOT (USED WITH LS)'         
      WRITE(NFIL,*)-'SET STYLE LINE 1 LT 1 LC RGB "BLACK" LW 1'             
      WRITE(NFIL,*)'# MAP HIGHT VALUES TO COLORS'                          
      WRITE(NFIL,*)-'SET PALETTE RGBFORMULA 22,13,-31'                      
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DATA RELATED STATEMENTS                                =='
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)-'SET XRANGE [XMIN:XMAX]'                                
      WRITE(NFIL,*)-'SET YRANGE [YMIN:YMAX]'                                
      WRITE(NFIL,*)-'SET ZRANGE [ZMIN:ZMAX]'                                
      WRITE(NFIL,*)-'SET DGRID3D  60,60,1','     #SAMPLE DATA ONTO A 60X60 GRID'
      WRITE(NFIL,*)'#'                                                     
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# SURFACE PLOT                                           =='
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)-'SET SURFACE'                                           
      WRITE(NFIL,*)-'#SET HIDDEN3D'                                          
      WRITE(NFIL,*)-'SET XYPLANE AT 0.','      # PLACE Z=0 INTO THE XY PLANE'
      WRITE(NFIL,*)-'UNSET BORDER','           # REMOVE AXES'
      WRITE(NFIL,*)-'UNSET XTICS','            # REMOVE TICS FROM THE AXES'
      WRITE(NFIL,*)-'UNSET YTICS','            # REMOVE TICS FROM THE AXES'
      WRITE(NFIL,*)-'UNSET ZTICS','            # REMOVE TICS FROM THE AXES'
      WRITE(NFIL,*)-'SET KEY OFF','            # REMOVE TITLE'
      IF(TYPE.EQ.'SURFACE') THEN
        WRITE(NFIL,*)-'SET PM3D HIDDEN3D 1','  # LINESTYLE FOR THE SURFACE GRID'
      ELSE IF(TYPE.EQ.'CONTOUR') THEN
        WRITE(NFIL,*)-'#SET PM3D HIDDEN3D 1',' # LINESTYLE FOR THE SURFACE GRID'
      END IF
      WRITE(NFIL,*)'#'                                                     
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DEFINE CONTOUR                                         =='
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)-'SET CONTOUR SURFACE',' # DRAW CONTOUR ONTO THE SURFACE'
      WRITE(NFIL,*)-'SET CNTRPARAM LEVELS INCREMENTAL -0.5,0.05,0.5'
      WRITE(NFIL,*)-'SET CNTRPARAM BSPLINE'
      WRITE(NFIL,*)-'SET CNTRPARAM ORDER 6'
      WRITE(NFIL,*)-'UNSET CLABEL'       ,' # NO AUTOCOLORING OF CONTOURS'
      WRITE(NFIL,*)'#'                                                     
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DEFINE ORIENTATION                                     =='
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)-'SET VIEW EQUAL XY # AVOIDS DISTORTIONS OF XY-PLANE'
      WRITE(NFIL,*)'#  ANGLE, ANGLE, OVERALL SCALE, SCALE Z-AXIS'        
      WRITE(NFIL,*)-'SET VIEW ROT_X,ROT_Z,SCALE,SCALE_Z'                  
      WRITE(NFIL,*)'#'                                                   
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# SET TERMINALS (UNCOMMENT ONE)                          =='
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'#----USE POSTSCRIPT TERMINAL FOR EPS FILES-----------------'
      WRITE(NFIL,*)-'SET TERMINAL POSTSCRIPT EPS ' &
                  ,-'SIZE (XMAX-XMIN)/2,(YMAX-YMIN)/2 ' &
     &            ,-" ENHANCED COLOR CLIP FONT 'HELVETICA,20' LINEWIDTH 1 "
      WRITE(NFIL,*)'#----USE WXT TERMINAL FOR LINUX SCREEN---------------------'
      WRITE(NFIL,*)-'# SET TERMINAL WXT SIZE 350,262 ENHANCED ' &
     &            ,-" FONT 'VERDANA,10' PERSIST "
      WRITE(NFIL,*)'#----USE PDF TERMINAL FOR PDF FILES------------------------'
      WRITE(NFIL,*)-'# SET TERMINAL PDF COLOR ENHANCED '
      WRITE(NFIL,*)'#----USE AQUA TERMINAL FOR OSX SCREEN----------------------'
      WRITE(NFIL,*)-"# SET TERMINAL AQUA ENHANCED SOLID FONT 'HELVETICA,20'" &
     &            ,-"TITLE 'CONTOUR'"
      WRITE(NFIL,*)'#'                                                   
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# PERFORM PLOT                                           =='
      WRITE(NFIL,*)'#=========================================================='
      IF(TYPE.EQ.'SURFACE') THEN
        WRITE(NFIL,*)-'UNSET CONTOUR'
        WRITE(NFIL,*)-'# SET BORDER'
        WRITE(NFIL,*)-'# UNSET SURFACE'
        WRITE(NFIL,*)-'# UNSET COLORBOX'
      ELSE IF(TYPE.EQ.'CONTOUR') THEN
        WRITE(NFIL,*)-'# UNSET CONTOUR'
        WRITE(NFIL,*)-'SET BORDER'
        WRITE(NFIL,*)-'# UNSET SURFACE'
        WRITE(NFIL,*)-'UNSET COLORBOX'
      END IF
      WRITE(NFIL,*)-"SPLOT '-' WITH PM3D LINEWIDTH 3 LINECOLOR RGB '#000000'"
      WRITE(NFIL,*)'#'                                                     
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DATA SECTION                                           =='
      WRITE(NFIL,*)'#=========================================================='
      DO I=1,N1
        DO J=1,N2
          WRITE(NFIL,*)(REAL(I-1)/REAL(N1-1)-0.5D0)*SQRT(SUM(BOX(:,1)**2)) &
     &              ,(REAL(J-1)/REAL(N2-1)-0.5D0)*SQRT(SUM(BOX(:,2)**2)) &
     &              ,DATA(I,J,1)
        ENDDO
      ENDDO
      WRITE(NFIL,*)'#=========================================================='
      WRITE(NFIL,*)'# DATA SECTION FINISHED                                  =='
      WRITE(NFIL,*)'#=========================================================='
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAPFIELDTOGRID(RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX,N1,N2,N3,DATA)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) !LATTICE VECTORS OF THE PERIODIC DATA
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3 ! # GRID POINTS ON THE PERIODIC GRID
      REAL(8)   ,INTENT(IN) :: WAVE(NR1,NR2,NR3) ! PERIODIC DENSITY DATA
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! CORNER OF THE NEW GRID
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! VECTORS SPANNING THE NEW GRID
      INTEGER(4),INTENT(IN) :: N1,N2,N3
      REAL(8)   ,INTENT(OUT):: DATA(N1,N2,N3)
      REAL(8)               :: POS(3)
      REAL(8)               :: DT1(3),DT2(3),DT3(3)
      REAL(8)               :: XTOR(3,3)
      REAL(8)               :: RTOX(3,3)
      REAL(8)               :: SVAR1,SVAR2
      REAL(8)               :: F00,F01,F10,F11,G0,G1
      REAL(8)               :: XDIS(3)
      INTEGER(4)            :: I1M,I1P,I2M,I2P,I3M,I3P
      INTEGER(4)            :: I,J,K
!     **************************************************************************
                                CALL TRACE$PUSH('MAPFIELDTOGRID')
      XTOR(:,1)=RBAS(:,1)/REAL(NR1,KIND=8)
      XTOR(:,2)=RBAS(:,2)/REAL(NR2,KIND=8)
      XTOR(:,3)=RBAS(:,3)/REAL(NR3,KIND=8)
      CALL LIB$INVERTR8(3,XTOR,RTOX)
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE ONTO NEW GRID                                      ==
!     ==========================================================================
      DT1(:)=BOX(:,1)/REAL(N1-1,KIND=8)
      DT2(:)=BOX(:,2)/REAL(N2-1,KIND=8)
      IF(N3.GT.1) THEN
        DT3(:)=BOX(:,3)/REAL(N3-1,KIND=8)
      ELSE                   ! N3=1 IS FOR A 2D GRID
        DT3(:)=0.D0
      END IF
      DO I=1,N1
        DO J=1,N2
          POS(:)=ORIGIN(:)+DT1*REAL(I-1,KIND=8)+DT2*REAL(J-1,KIND=8)-DT3(:)
          DO K=1,N3
            POS(:)=POS(:)+DT3(:)  ! POSITION ON THE NEW GRID
!           == MAP POSITION INTO RELATIVE COORDINATES AND INTO 1ST UNIT CELL
            XDIS=MATMUL(RTOX,POS)
            XDIS(1)=MODULO(XDIS(1),REAL(NR1,KIND=8))
            XDIS(2)=MODULO(XDIS(2),REAL(NR2,KIND=8))
            XDIS(3)=MODULO(XDIS(3),REAL(NR3,KIND=8))
!           == AVOID COMPILER BUG. FOR VERY SMALL NEGATIVE VALUES A RESULT  ==
!           == EQUAL TO THE RIGHT BOUNDARY CAN OCCUR
            IF(XDIS(1).EQ.REAL(NR1,KIND=8))XDIS(1)=XDIS(1)-REAL(NR1,KIND=8)
            IF(XDIS(2).EQ.REAL(NR2,KIND=8))XDIS(2)=XDIS(2)-REAL(NR2,KIND=8)
            IF(XDIS(3).EQ.REAL(NR3,KIND=8))XDIS(3)=XDIS(3)-REAL(NR3,KIND=8)
            I1M=1+INT(XDIS(1))
            I1P=1+MODULO(I1M,NR1)
            I2M=1+INT(XDIS(2))
            I2P=1+MODULO(I2M,NR2)
            I3M=1+INT(XDIS(3))
            I3P=1+MODULO(I3M,NR3)
            XDIS(1)=MODULO(XDIS(1),1.D0)
            XDIS(2)=MODULO(XDIS(2),1.D0)
            XDIS(3)=MODULO(XDIS(3),1.D0)
!           == INTERPOLATE ALONG THIRD DIRECTION ===========================
            SVAR2=XDIS(3)
            SVAR1=1.D0-SVAR2
            F00=WAVE(I1M,I2M,I3M)*SVAR1+WAVE(I1M,I2M,I3P)*SVAR2
            F10=WAVE(I1P,I2M,I3M)*SVAR1+WAVE(I1P,I2M,I3P)*SVAR2
            F01=WAVE(I1M,I2P,I3M)*SVAR1+WAVE(I1M,I2P,I3P)*SVAR2
            F11=WAVE(I1P,I2P,I3M)*SVAR1+WAVE(I1P,I2P,I3P)*SVAR2
!           == INTERPOLATE ALONG SECOND DIRECTION ===========================
            SVAR2=XDIS(2)
            SVAR1=1.D0-SVAR2
            G0=F00*SVAR1+F01*SVAR2
            G1=F10*SVAR1+F11*SVAR2
!           == INTERPOLATE ALONG SECOND DIRECTION ===========================
            SVAR2=XDIS(1)
            SVAR1=1.D0-SVAR2
            DATA(I,J,K)=G0*SVAR1+G1*SVAR2
          ENDDO
        ENDDO
      ENDDO
                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAPFIELDTOGRIDC(RBAS,XK,NR1,NR2,NR3,WAVE &
     &                           ,ORIGIN,BOX,N1,N2,N3,DATA)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) !LATTICE VECTORS OF THE PERIODIC DATA
      REAL(8)   ,INTENT(IN) :: XK(3)     !K-POINT IN RELATIVE COORDINATES
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3 ! # GRID POINTS ON THE PERIODIC GRID
      COMPLEX(8),INTENT(IN) :: WAVE(NR1,NR2,NR3) ! PERIODIC DENSITY DATA
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! CORNER OF THE NEW GRID
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! VECTORS SPANNING THE NEW GRID
      INTEGER(4),INTENT(IN) :: N1,N2,N3
      COMPLEX(8),INTENT(OUT):: DATA(N1,N2,N3)
      REAL(8)               :: POS(3)
      REAL(8)               :: DT1(3),DT2(3),DT3(3)
      REAL(8)               :: XTOR(3,3)
      REAL(8)               :: RTOX(3,3)
      REAL(8)               :: SVAR1,SVAR2
      COMPLEX(8)            :: CSVAR1,CSVAR2
      COMPLEX(8)            :: EIK1M,EIK1P,EIK2M,EIK2P,EIK3M,EIK3P
      COMPLEX(8)            :: F00,F01,F10,F11,G0,G1
      REAL(8)               :: XDIS(3),XDIS0(3)
      INTEGER(4)            :: I1M,I1P,I2M,I2P,I3M,I3P
      INTEGER(4)            :: I,J,K
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: PI
      COMPLEX(8)            :: C2PII
!     **************************************************************************
                                CALL TRACE$PUSH('MAPFIELDTOGRIDC')
      PI=4.D0*ATAN(1.D0)
      C2PII=(0.D0,1.D0)*2.D0*PI
      XTOR(:,1)=RBAS(:,1)/REAL(NR1,KIND=8)
      XTOR(:,2)=RBAS(:,2)/REAL(NR2,KIND=8)
      XTOR(:,3)=RBAS(:,3)/REAL(NR3,KIND=8)
      CALL LIB$INVERTR8(3,XTOR,RTOX)
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE ONTO NEW GRID                                      ==
!     ==========================================================================
      DT1(:)=BOX(:,1)/REAL(N1-1,KIND=8)
      DT2(:)=BOX(:,2)/REAL(N2-1,KIND=8)
      IF(N3.GT.1) THEN
        DT3(:)=BOX(:,3)/REAL(N3-1,KIND=8)
      ELSE                   ! N3=1 IS FOR A 2D GRID
        DT3(:)=0.D0
      END IF
      DO I=1,N1
        DO J=1,N2
          POS(:)=ORIGIN(:)+DT1*REAL(I-1,KIND=8)+DT2*REAL(J-1,KIND=8)-DT3(:)
          DO K=1,N3
            POS(:)=POS(:)+DT3(:)  ! POSITION ON THE NEW GRID
!           == MAP POSITION INTO RELATIVE COORDINATES AND INTO 1ST UNIT CELL ===
            XDIS0=MATMUL(RTOX,POS)
            XDIS(1)=MODULO(XDIS0(1),REAL(NR1,KIND=8))
            XDIS(2)=MODULO(XDIS0(2),REAL(NR2,KIND=8))
            XDIS(3)=MODULO(XDIS0(3),REAL(NR3,KIND=8))
!           == AVOID COMPILER BUG. FOR VERY SMALL NEGATIVE VALUES A RESULT  ====
!           == EQUAL TO THE RIGHT BOUNDARY CAN OCCUR ===========================
            IF(XDIS(1).EQ.REAL(NR1,KIND=8))XDIS(1)=XDIS(1)-REAL(NR1,KIND=8)
            IF(XDIS(2).EQ.REAL(NR2,KIND=8))XDIS(2)=XDIS(2)-REAL(NR2,KIND=8)
            IF(XDIS(3).EQ.REAL(NR3,KIND=8))XDIS(3)=XDIS(3)-REAL(NR3,KIND=8)
            I1M=1+INT(XDIS(1))
            I1P=1+MODULO(I1M,NR1)
            I2M=1+INT(XDIS(2))
            I2P=1+MODULO(I2M,NR2)
            I3M=1+INT(XDIS(3))
            I3P=1+MODULO(I3M,NR3)
            EIK1M=EXP(-C2PII*XK(1)*(XDIS(1)-XDIS0(1))/REAL(NR1,KIND=8))
            EIK2M=EXP(-C2PII*XK(2)*(XDIS(2)-XDIS0(2))/REAL(NR2,KIND=8))
            EIK3M=EXP(-C2PII*XK(3)*(XDIS(3)-XDIS0(3))/REAL(NR3,KIND=8))
            EIK1P=EIK1M*EXP(-C2PII*XK(1)*REAL(I1P-1-I1M)/REAL(NR1,KIND=8))
            EIK2P=EIK2M*EXP(-C2PII*XK(2)*REAL(I2P-1-I2M)/REAL(NR2,KIND=8))
            EIK3P=EIK3M*EXP(-C2PII*XK(3)*REAL(I3P-1-I3M)/REAL(NR3,KIND=8))
!
            XDIS(1)=MODULO(XDIS(1),1.D0)
            XDIS(2)=MODULO(XDIS(2),1.D0)
            XDIS(3)=MODULO(XDIS(3),1.D0)
!           == INTERPOLATE ALONG THIRD DIRECTION ===============================
            SVAR2=XDIS(3)
            SVAR1=1.D0-SVAR2
            CSVAR1=SVAR1*EIK3M
            CSVAR2=SVAR2*EIK3P
            F00=WAVE(I1M,I2M,I3M)*CSVAR1+WAVE(I1M,I2M,I3P)*CSVAR2
            F10=WAVE(I1P,I2M,I3M)*CSVAR1+WAVE(I1P,I2M,I3P)*CSVAR2
            F01=WAVE(I1M,I2P,I3M)*CSVAR1+WAVE(I1M,I2P,I3P)*CSVAR2
            F11=WAVE(I1P,I2P,I3M)*CSVAR1+WAVE(I1P,I2P,I3P)*CSVAR2
!           == INTERPOLATE ALONG SECOND DIRECTION ==============================
            SVAR2=XDIS(2)
            SVAR1=1.D0-SVAR2
            CSVAR1=SVAR1*EIK2M
            CSVAR2=SVAR2*EIK2P
            G0=F00*CSVAR1+F01*CSVAR2
            G1=F10*CSVAR1+F11*CSVAR2
!           == INTERPOLATE ALONG SECOND DIRECTION ==============================
            SVAR2=XDIS(1)
            SVAR1=1.D0-SVAR2
            CSVAR1=SVAR1*EIK1M
            CSVAR2=SVAR2*EIK1P
            DATA(I,J,K)=G0*CSVAR1+G1*CSVAR2
          ENDDO
        ENDDO
      ENDDO
                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEATOMSINTOBOX(NAT,Z,R,RBAS,ORIGIN,BOX &
     &                           ,NRVIEWX,NRVIEW,RVIEW,ZVIEW)
!     **************************************************************************
!     ** MAPS PERIODIC LATTICE INTO VIEWBOX                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    ! ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! ATOMIC POSITIONS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS OF THE PERIODIC DATA
      REAL(8)   ,INTENT(IN) :: ORIGIN(3) ! CORNER OF THE NEW GRID
      REAL(8)   ,INTENT(IN) :: BOX(3,3)  ! VECTORS SPANNING THE NEW GRID
      INTEGER(4),INTENT(IN) :: NRVIEWX   ! X#(ATOMS IN VIEWBOX)
      INTEGER(4),INTENT(OUT):: NRVIEW    ! #(ATOMS IN VIEWBOX)
      REAL(8)   ,INTENT(OUT):: RVIEW(3,NRVIEWX)
      REAL(8)   ,INTENT(OUT):: ZVIEW(NRVIEWX)
      REAL(8)               :: POS(3)
      REAL(8)               :: RTOX(3,3)
      REAL(8)               :: XDIS(3)
      REAL(8)               :: T1(3),T2(3),T3(3)
      INTEGER(4)            :: IT1,IT2,IT3,IAT
      INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
!     **************************************************************************
!
!     ==========================================================================
!     ==  MAP ATOMS INTO THE VIEWBOX                                          ==
!     ==========================================================================
      CALL LIB$INVERTR8(3,BOX,RTOX)
      MIN1=-5
      MAX1=5
      MIN2=-5
      MAX2=5
      MIN3=-5
      MAX3=5
      NRVIEW=0
      DO IT1=MIN1,MAX1
        T1(:)=RBAS(:,1)*REAL(IT1,KIND=8)
        DO IT2=MIN2,MAX2
          T2(:)=RBAS(:,2)*REAL(IT2,KIND=8)
          DO IT3=MIN3,MAX3
            T3(:)=RBAS(:,3)*REAL(IT3,KIND=8)
            DO IAT=1,NAT
              POS(:)=R(:,IAT)+T1+T2+T3
              XDIS(:)=MATMUL(RTOX,POS-ORIGIN(:))
              IF(XDIS(1).LT.-0.05D0.OR.XDIS(1).GT.1.05D0) CYCLE
              IF(XDIS(2).LT.-0.05D0.OR.XDIS(2).GT.1.05D0) CYCLE
              IF(XDIS(3).LT.-0.05D0.OR.XDIS(3).GT.1.05D0) CYCLE
              NRVIEW=NRVIEW+1
              IF(NRVIEW.GT.NRVIEWX) THEN
                CALL ERROR$MSG('NUMBER OF ATOMS IN THE BOX OUT OF RANGE')
                CALL ERROR$I4VAL('NRVIEWX',NRVIEWX)
                CALL ERROR$I4VAL('NRVIEW',NRVIEW)
                CALL ERROR$STOP('MAKECUBE')
              END IF
              RVIEW(:,NRVIEW)=POS(:)
              ZVIEW(NRVIEW)=Z(IAT)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITECUBEFILE(NFIL,NAT,Z,R,ORIGIN,BOX,N1,N2,N3,DATA)
!     **************************************************************************
!     **  WRITE A GAUSSIAN CUBE FILE (EXTENSION .CUB) WITH VOLUMINETRIC DATA  **
!     **                                                                      **
!     ** A VERY CLEAR DESCRIPTION OF THE GAUSSIAN CUBE FORMAT HAS BEEN GIVEN  **
!     **  ON HTTP://LOCAL.WASP.UWA.EDU.AU/~PBOURKE/DATAFORMATS/CUBE/          **
!     **                                                                      **
!     ** REMARK:                                                              **
!     ** UNITS WRITTEN ARE ABOHR, CONSISTENT WITH AVOGADROS IMPLEMENTATION.   **
!     ** THE SPECS REQUIRE N1,N2,N3 TO BE MULTIPLIED BY -1 IF ABOHR ARE USED  **
!     ** AND ANGSTROM IS THE UNIT IF THEY ARE POSITIVE.                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    ! ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)
      REAL(8)   ,INTENT(IN) :: BOX(3,3)
      INTEGER(4),INTENT(IN) :: N1,N2,N3
      REAL(8)   ,INTENT(IN) :: DATA(N1,N2,N3)
      REAL(8)               :: ANGSTROM
      REAL(8)               :: SCALE
      INTEGER(4)            :: SUNIT=1   ! =1 FOR BOHR RADII/ =-1 FOR ANGSTROM
      INTEGER(4)            :: IAT,I,J,K
!     **************************************************************************
      IF(SUNIT.EQ.1) THEN       ! LENGTHS IN BOHR RADII
        SCALE=1.D0
      ELSE IF(SUNIT.EQ.-1) THEN ! LENGTHS IN ANGSTROM
        CALL CONSTANTS('ANGSTROM',ANGSTROM)
        SCALE=1.D0/ANGSTROM
      ELSE
        CALL ERROR$MSG('ILLEGAL VALUE OF SUNIT. MUST BE EITHER +1 OR -1')
        CALL ERROR$I4VAL('SUNIT',SUNIT)
        CALL ERROR$STOP('WRITECUBEFILE')
      END IF
      WRITE(NFIL,FMT='("CP-PAW CUBE FILE")')
      WRITE(NFIL,FMT='("NOCHN KOMMENTAR")')
      WRITE(NFIL,FMT='(I5,3F12.6)')NAT,ORIGIN*SCALE
!     == SHEARED GRIDS ARE PERMITTED BUT OFTEN NOT SUPPORTED ===================
!     == MINUS SIGN INDICATS ANGSTROM-UNIT. OTHERWISE BOHR RADII ===============
      WRITE(NFIL,FMT='(I5,3F12.6)')SUNIT*N1,BOX(:,1)/REAL(N1-1,KIND=8)*SCALE
      WRITE(NFIL,FMT='(I5,3F12.6)')SUNIT*N2,BOX(:,2)/REAL(N2-1,KIND=8)*SCALE
      WRITE(NFIL,FMT='(I5,3F12.6)')SUNIT*N3,BOX(:,3)/REAL(N3-1,KIND=8)*SCALE
!
      DO IAT=1,NAT
        WRITE(NFIL,FMT='(I5,4F12.6)')NINT(Z(IAT)),0.D0,R(:,IAT)*SCALE
      ENDDO  
!     == DATA ARE WRITTEN WITH Z AS INNERMOST AND X AS OUTERMOST LOOP
      WRITE(NFIL,FMT='(6(E12.6," "))')(((DATA(I,J,K),K=1,N3),J=1,N2),I=1,N1)
!!$DO K=1,N3
!!$  SCALE=0.D0
!!$  DO I=1,N1
!!$    DO J=1,N2
!!$      SCALE=SCALE+DATA(I,J,K)**2
!!$    ENDDO
!!$  ENDDO
!!$PRINT*,'SCALE ',ORIGIN(:)+BOX(:,3)/REAL(N3,KIND=8)*REAL(K-1,KIND=8),SCALE
!!$ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VRML$ADDRUBBERSHEET(NFIL,ORIGIN,PLANE,N1,N2,DATA)
!     **************************************************************************
!     **  ADD A BALL-STICK MODEL TO A VRML SCENE                              **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)
      REAL(8)   ,INTENT(IN) :: PLANE(3,2)
      REAL(8)   ,INTENT(IN) :: DATA(N1,N2)
      REAL(8)   ,PARAMETER  :: SCALE=10.D0
      REAL(8)               :: ROT(3,3)
      REAL(8)               :: AXIS(3)
      REAL(8)               :: ANGLE
      REAL(8)               :: D1,D2
      REAL(8)               :: SVAR1,SVAR2,X
      REAL(8)               :: COLOR(3,N1,N2)
      INTEGER(4)            :: I,J
      INTEGER(4)            :: R,G,B
      LOGICAL(4)            :: TINV
!     **************************************************************************
                             CALL TRACE$PUSH('VRML$ADDRUBBERSHEET')
      D1=SQRT(SUM(PLANE(:,1)**2))
      D2=SQRT(SUM(PLANE(:,2)**2))
!     == CONSTRUCT ROTATION MATRIX =============================================
!     == FIRST DIRECTION IS MAPPED ONTO X-AXIS
!     == SECOND DIRECTION IS MAPPED ONTO Z-AXIS
      ROT(:,1)=PLANE(:,1)/D1
      ROT(:,3)=PLANE(:,2)/D2
      ROT(1,2)=ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)
      ROT(2,2)=ROT(3,3)*ROT(1,1)-ROT(1,3)*ROT(3,1)
      ROT(3,2)=ROT(1,3)*ROT(2,1)-ROT(2,3)*ROT(1,1)
!     == EXTRACT ROTATION ANGLE ================================================
!     == SEE HTTP://WWW.EUCLIDEANSPACE.COM/MATHS/GEOMETRY/ROTATIONS           ==
!     ==                                           /CONVERSIONS/MATRIXTOANGLE ==
      CALL ROTATION$MATRIXTOANGLE(ROT,AXIS,TINV)
      ANGLE=SQRT(SUM(AXIS**2))
      IF(ANGLE.NE.0.D0) THEN
         AXIS=AXIS/ANGLE
      ELSE
         AXIS(:)=(/0.D0,0.D0,1.D0/)
      END IF
!
!     == CREATE COLOR ===== ====================================================
      COLOR(:,:,:)=0.D0
      SVAR1=MINVAL(DATA)
      SVAR2=MAXVAL(DATA)
      DO I=1,N1
        DO J=1,N1
          X=(DATA(I,J)-SVAR1)/(SVAR2-SVAR1)
          CALL RGB$RAINBOW(X,R,G,B)
          COLOR(1,I,J)=REAL(R)/REAL(255)
          COLOR(2,I,J)=REAL(G)/REAL(255)
          COLOR(3,I,J)=REAL(B)/REAL(255)
        ENDDO
      ENDDO
!
!     == CREATE WRITE MODEL TO SCENE ===========================================
      WRITE(NFIL,*)'TRANSFORM{TRANSLATION ',ORIGIN,' ROTATION ',AXIS,ANGLE
      WRITE(NFIL,*)'CHILDREN['
!
!     == CREATE RUBBERSHEET ====================================================
      WRITE(NFIL,*)'SHAPE{'
      WRITE(NFIL,*)'  APPEARANCE APPEARANCE {'
      WRITE(NFIL,*)'    MATERIAL MATERIAL {'
      WRITE(NFIL,*)'      DIFFUSECOLOR 0.8 0.4 0.0'
      WRITE(NFIL,*)'      SHININESS 1'
      WRITE(NFIL,*)'    } # CLOSE MATERIAL'
      WRITE(NFIL,*)'  }   # CLOSE APPEARANCE'
      WRITE(NFIL,*)'  GEOMETRY ELEVATIONGRID {'
      WRITE(NFIL,*)'    XDIMENSION ',N1
      WRITE(NFIL,*)'    ZDIMENSION ',N2
      WRITE(NFIL,*)'    XSPACING ',D1/REAL(N1-1)
      WRITE(NFIL,*)'    ZSPACING ',D2/REAL(N2-1)
      WRITE(NFIL,*)'    SOLID FALSE # MAKES SHEET VISIBLE FROM BOTH SIDES'   
      WRITE(NFIL,*)'    CREASEANGLE 3 # SMOOTHENING'
      WRITE(NFIL,*)'    COLOR COLOR {COLOR ['
      WRITE(NFIL,FMT='(5(3F10.5," , "))')COLOR
      WRITE(NFIL,*)'      ]  # END OF COLOR'
      WRITE(NFIL,*)'    }    # END OF COLOR'
      WRITE(NFIL,*)'    HEIGHT ['
      WRITE(NFIL,FMT='(20F10.5)')-DATA*SCALE
      WRITE(NFIL,*)'    ]}}]}'  !END OF HEIGHT,GEOMETRY,SHAPE,CHILDREN,TRANSFORM
                             CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VRML$ADDBALLSTICK(NFIL,NAT,Z,R)
!     **************************************************************************
!     **  ADD A BALL-STICK MODEL TO A VRML SCENE                              **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,PARAMETER  :: SCALEBALL=0.6D0
      REAL(8)   ,PARAMETER  :: BONDDIAMETER=0.2D0
      REAL(8)   ,PARAMETER  :: BONDLENGTHCUTOFF=1.3D0
      REAL(8)   ,PARAMETER  :: BONDCOLOR(3)=(/0.3D0,0.3D0,0.3D0/)
      INTEGER(4)            :: IAT,IAT1,IAT2
      REAL(8)               :: RAD,RAD1,RAD2
      REAL(8)               :: COLOR(3)
      INTEGER(4)            :: IVEC(3)
      REAL(8)               :: DR(3),R1(3),R2(3)
      REAL(8)               :: DIS
      REAL(8)               :: AXIS(3),ANGLE
!     **************************************************************************
!
!      == CREATE BALLS =========================================================
       DO IAT=1,NAT
         CALL PERIODICTABLE$GET(NINT(Z(IAT)),'R(COV)',RAD)
         CALL ATOMCOLOR(NINT(Z(IAT)),IVEC)
         COLOR(:)=REAL(IVEC,KIND=8)/255.D0
         WRITE(NFIL,*)'TRANSFORM{'
         WRITE(NFIL,*)'  TRANSLATION ',R(:,IAT)
         WRITE(NFIL,*)'  CHILDREN['
         WRITE(NFIL,*)'    SHAPE{'
         WRITE(NFIL,*)'      APPEARANCE APPEARANCE {'
         WRITE(NFIL,*)'        MATERIAL MATERIAL {'
         WRITE(NFIL,*)'          DIFFUSECOLOR ',COLOR
         WRITE(NFIL,*)'          SHININESS 1'
         WRITE(NFIL,*)'          TRANSPARENCY 0.5'
         WRITE(NFIL,*)'        } # CLOSE MATERIAL'
         WRITE(NFIL,*)'      } # CLOSE APPEARANCE'
         WRITE(NFIL,*)'      GEOMETRY SPHERE{ RADIUS ',RAD*SCALEBALL,'}'
         WRITE(NFIL,*)'    }   # CLOSE SHAPE'
         WRITE(NFIL,*)'  ]     # CLOSE CHILDREN'
         WRITE(NFIL,*)'}       # CLOSE TRANSFORM'
       ENDDO
!
!      == CREATE BONDS =========================================================
       DO IAT1=1,NAT
         CALL PERIODICTABLE$GET(NINT(Z(IAT1)),'R(COV)',RAD1)
         DO IAT2=IAT1+1,NAT
           CALL PERIODICTABLE$GET(NINT(Z(IAT2)),'R(COV)',RAD2)
           DR(:)=R(:,IAT2)-R(:,IAT1)
           DIS=SQRT(SUM(DR**2))
           IF(DIS/(RAD1+RAD2).GT.BONDLENGTHCUTOFF) CYCLE           
           DR(:)=DR(:)/DIS
           R1(:)=R(:,IAT1)+RAD1*SCALEBALL*1.001D0*DR(:)
           R2(:)=R(:,IAT2)-RAD2*SCALEBALL*1.001D0*DR(:)
           ANGLE=-ACOS(DR(2))   ! ORIGINAL ORIENTATION (0,1,0)
           IF(ABS(DR(2))-1.D0.GT.1.D-8) THEN
             AXIS(1)=-DR(3)
             AXIS(2)=0.D0
             AXIS(3)=DR(1)
             AXIS=AXIS/SQRT(SUM(AXIS**2))
           ELSE  ! EXCEPTION TO AVOID DIVIDE-BY-ZERO
             AXIS(:)=0.D0  
             AXIS(1)=1.D0
           END IF
           WRITE(NFIL,*)'TRANSFORM{'
           WRITE(NFIL,*)'  TRANSLATION ',0.5D0*(R1(:)+R2(:))
           WRITE(NFIL,*)'  ROTATION    ',AXIS,ANGLE
           WRITE(NFIL,*)'  CHILDREN ['
           WRITE(NFIL,*)'    SHAPE{'
           WRITE(NFIL,*)'      APPEARANCE APPEARANCE {'
           WRITE(NFIL,*)'        MATERIAL MATERIAL {'
           WRITE(NFIL,*)'          DIFFUSECOLOR ',BONDCOLOR
           WRITE(NFIL,*)'          SHININESS 1'
           WRITE(NFIL,*)'          TRANSPARENCY 0.7'
           WRITE(NFIL,*)'        } # CLOSE MATERIAL'
           WRITE(NFIL,*)'      }   # CLOSE APPEARANCE'
           WRITE(NFIL,*)'      GEOMETRY CYLINDER{'
           WRITE(NFIL,*)'        RADIUS ',BONDDIAMETER
           WRITE(NFIL,*)'        HEIGHT ',SQRT(SUM((R2-R1)**2))
           WRITE(NFIL,*)'      }   # CLOSE GEOMETRY'
           WRITE(NFIL,*)'    }     # CLOSE SHAPE'
           WRITE(NFIL,*)'  ]       # CLOSE CHILDREN'
           WRITE(NFIL,*)'}         # CLOSE TRANSFORM'
         ENDDO
       ENDDO
      RETURN
      END
