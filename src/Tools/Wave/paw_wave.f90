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
       CHARACTER(256)            :: FILE
       CHARACTER(32)             :: ID
       LOGICAL(4)                :: TCHK
       INTEGER(4)                :: I
       INTEGER(4)                :: NLISTS
       INTEGER(4)                :: IARGC
!      ******************************************************************
       I=IARGC()
       IF(I.NE.1) THEN
         WRITE(*,FMT='(A)')'CORRECT USAGE: paw_wave.x [FILE]'
         WRITE(*,FMT='(A)')'WHERE [FILE] IS THE FILE CONTROL FILE NAME OF THIS TOOL'
         WRITE(*,FMT='(A)')'SEE PAW_MANUAL FOR FURTHER DESCRIPTION'
         CALL ERROR$NORMALSTOP
       END IF
       CALL GETARG(1,FILE)
       IF(FILE.EQ.'?'.OR.FILE.EQ.-'-H') THEN
         WRITE(*,FMT='(A)')'CORRECT USAGE: PAW_WAVE [FILE]'
         WRITE(*,FMT='(A)')'WHERE [FILE] IS THE FILE CONTROL FILE NAME OF THIS TOOL'
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
      CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
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
!
!     ==================================================================
!     ==  READ CNTL FILE TO LINKEDLIST                                ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('CNTL',NFIL)
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
!
!     ==================================================================
!     ==  RESET FILE NAME FOR STRUCTURE FILE IF REQUESTED             ==
!     ================================================================== 
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
!     ==================================================================
!     ==  READ STRC FILE TO LINKEDLIST                                ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$NEW(LL_STRC)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL MWAVE(LL_CNTL,LL_STRC)
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
       REAL(8)                   :: RBAS(3,3)
       INTEGER(4)                :: NAT
       REAL(8)      ,ALLOCATABLE :: POS(:,:)   !(3,NAT)
       REAL(8)      ,ALLOCATABLE :: Z(:)       !(NAT)
       CHARACTER(32),ALLOCATABLE :: ATOMNAME(:)
       INTEGER(4)                :: NR1,NR2,NR3
       REAL(8)      ,ALLOCATABLE :: WAVE(:,:,:)
       complex(8)   ,ALLOCATABLE :: cWAVE(:,:,:)
       REAL(8)                   :: BOXR0(3)
       REAL(8)                   :: BOXVEC(3,3)
       REAL(8)                   :: planer0(3)
       REAL(8)                   :: planeVEC(3,2)
       REAL(8)      ,ALLOCATABLE :: POSM(:,:)  ! POSITIONS
       INTEGER(4)   ,ALLOCATABLE :: MAP(:)     !MAPPING FROM MODEL 
       REAL(8)      ,ALLOCATABLE :: RAD(:)   
       REAL(8)      ,ALLOCATABLE :: COLOR(:,:)
       INTEGER(4)                :: NBOND
       INTEGER(4)   ,ALLOCATABLE :: BOND(:,:)
       LOGICAL(4)                :: TCHK
       INTEGER(4)                :: NATM
       INTEGER(4)                :: IAT
       INTEGER(4)                :: IVEC(3)
       REAL(8)      ,ALLOCATABLE :: SCALEDRAD(:)
       logical(4)                :: tplane
       logical(4)                :: tc
       real(8)                   :: xk(3)
!      ******************************************************************
       CALL TRACE$PUSH('MWAVE')
       LL_CNTL=LL_CNTL_
       LL_STRC=LL_STRC_
!
!      ==================================================================
!      ==  READ INPUT FILE                                             ==
!      ==================================================================
       CALL TRACE$PASS('READING WAVEPLOT')
       CALL FILEHANDLER$UNIT('WAVE',NFIL)
       call READWAVEPLOTdim(NFIL,TITLE,NAT,nr1,nr2,nr3,tc)
       ALLOCATE(POS(3,NAT))
       ALLOCATE(Z(NAT))
       ALLOCATE(ATOMNAME(NAT))
       if(tc) then
         ALLOCATE(CWAVE(NR1,NR2,NR3))
         CALL READWAVEPLOTC(NFIL,TITLE,RBAS,NAT,POS,Z,ATOMNAME,xk,NR1,NR2,NR3,cWAVE)
       else
         ALLOCATE(WAVE(NR1,NR2,NR3))
         CALL READWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,ATOMNAME,xk,NR1,NR2,NR3,WAVE)
       end if
       CALL TRACE$PASS('WAVEPLOT FILE READ')
!
!      ==================================================================
!      ==  GET VIEWBOX                                                 ==
!      ==================================================================
       CALL LINKEDLIST$SELECT(LL_CNTL,'~')
       CALL LINKEDLIST$SELECT(LL_CNTL,'WCNTL')
       CALL LINKEDLIST$SELECT(LL_CNTL,'VIEWBOX')
       BOXR0(:)=0.D0
       BOXVEC(:,:)=RBAS(:,:)
       CALL LINKEDLIST$EXISTD(LL_CNTL,'O',1,TCHK)
       IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'O',1,BOXR0)
       CALL LINKEDLIST$EXISTD(LL_CNTL,'T',1,TCHK)
       IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'T',1,BOXVEC)
!
!      ==================================================================
!      ==  GET plane for rubbersheet                                   ==
!      ==================================================================
       CALL LINKEDLIST$SELECT(LL_CNTL,'~')
       CALL LINKEDLIST$SELECT(LL_CNTL,'WCNTL')
       CALL LINKEDLIST$EXISTL(LL_CNTL,'PLANE',1,TPLANE)
       PLANER0(:)=BOXR0
       PLANEVEC(:,:)=BOXVEC(:,:2)
       IF(TPLANE) THEN
         CALL LINKEDLIST$SELECT(LL_CNTL,'PLANE')
         CALL LINKEDLIST$EXISTD(LL_CNTL,'O',1,TCHK)
         IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'O',1,PLANER0)
         CALL LINKEDLIST$EXISTD(LL_CNTL,'T',1,TCHK)
         IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'T',1,PLANEVEC)
       END IF
       PLANEVEC(:,2)=PLANEVEC(:,2)-PLANEVEC(:,1) &
      &                       /DOT_PRODUCT(PLANEVEC(:,1),PLANEVEC(:,1)) &
      &                       *DOT_PRODUCT(PLANEVEC(:,1),PLANEVEC(:,2)) 
!
!      ==================================================================
!      ==  WRITE DENSITY TO CUBE FILE                                  ==
!      ==================================================================
!       if(.not.tc) then
         CALL TRACE$PASS('WRITE WAVE TO CUBE-FILE')
         CALL FILEHANDLER$UNIT('CUBE',NFIL)
         REWIND NFIL
allocate(wave(nr1,nr2,nr3))
wave=real(cwave)
         CALL MAKECUBE(NFIL,NAT,Z,POS,RBAS,NR1,NR2,NR3,WAVE,BOXR0,BOXVEC)
         CALL FILEHANDLER$CLOSE('CUBE')
!       end if
!
!      ==================================================================
!      ==  WRITE vrml scene TO wrl FILE                                ==
!      ==================================================================
       CALL TRACE$PASS('WRITE WAVE TO vrml-FILE')
       CALL FILEHANDLER$UNIT('VRML',NFIL)
       REWIND NFIL
       if(tplane.and.(.not.tc)) then
         CALL MAKEgnu(NFIL,NAT,Z,POS,RBAS,NR1,NR2,NR3,WAVE,BOXR0,BOXVEC)
       end if
       if(tc) then
print*,'xk ',xk
         CALL MAKEVRMLc(NFIL,NAT,Z,POS,RBAS,xk,NR1,NR2,NR3,cWAVE &
      &             ,BOXR0,BOXVEC,tplane,planer0,planevec)
       else
         CALL MAKEVRML(NFIL,NAT,Z,POS,RBAS,NR1,NR2,NR3,WAVE &
      &             ,BOXR0,BOXVEC,tplane,planer0,planevec)
       end if
       CALL FILEHANDLER$CLOSE('VRML')
stop 'forced'
!
!      ==================================================================
!      ==  WRITE DENSITY TO DATAEXPLORER FILE                          ==
!      ==================================================================
       CALL TRACE$PASS('WRITE WAVE TO dx-FILE')
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
      SUBROUTINE READWAVEPLOTdim(NFIL,TITLE,NAT,nr1,nr2,nr3,tc)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)   :: NFIL
      CHARACTER(*) ,INTENT(OUT ) :: TITLE
      INTEGER(4)   ,INTENT(OUT)  :: NAT
      logical(4)   ,INTENT(OUT)  :: tc
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
      SUBROUTINE READWAVEPLOT(NFIL,TITLE,RBAS,NAT,POS,Z,NAME,xk,NR1,NR2,NR3,WAVE)
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
      REAL(8)      ,INTENT(OUT)  :: xk(3)
      CHARACTER(32),INTENT(OUT)  :: NAME(NAT)
      INTEGER(4)   ,INTENT(IN)   :: NR1
      INTEGER(4)   ,INTENT(IN)   :: NR2
      INTEGER(4)   ,INTENT(IN)   :: NR3
      REAL(8)      ,INTENT(OUT)  :: WAVE(NR1,NR2,NR3)
      REAL(8)                    :: Q(NAT) !POINT CHARGES (NOT YET USED)
      integer(4)                 :: nr1_,nr2_,nr3_,nat_
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
      if(nat_.ne.nat.or.nr1_.ne.nr1.or.nr2_.ne.nr2.or.nr3_.ne.nr3) then
        CALL ERROR$MSG('INCOnsistent array size')
        CALL ERROR$i4VAL('nat on file ',nat_)
        CALL ERROR$i4VAL('nat on input',nat)
        CALL ERROR$i4VAL('nr1 on file ',nr1_)
        CALL ERROR$i4VAL('nr1 on input',nr1)
        CALL ERROR$i4VAL('nr2 on file ',nr2_)
        CALL ERROR$i4VAL('nr2 on input',nr2)
        CALL ERROR$i4VAL('nr3 on file ',nr3_)
        CALL ERROR$i4VAL('nr3 on input',nr3)
        CALL ERROR$STOP('READWAVEPLOT')
      end if
      IF(NAT.ne.0) THEN
        READ(NFIL)NAME
        READ(NFIL)Z
        READ(NFIL)POS
        READ(NFIL)Q
      END IF
      xk(:)=0.d0
      READ(NFIL)WAVE
      READ(NFIL)ENDID 
      IF(ENDID.NE.'END OF FILE') THEN
        call error$msg('data file corrupted')
        CALL ERROR$CHVAL('ENDID',endid)
        call error$stop('READWAVEPLOT')
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
      IF(NAT.ne.0) THEN
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
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! atomic positions
      REAL(8)   ,INTENT(IN) :: Rbas(3,3) !lattice vectors of the periodic data
      INTEGER(4),INTENT(IN) :: Nr1,Nr2,Nr3 ! # grid points on the periodic grid
      REAL(8)   ,INTENT(IN) :: wave(Nr1,Nr2,Nr3) ! periodic density data
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! corner of the new grid
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! vectors spanning the new grid
      integer(4),parameter  :: n1=80      ! displacement
      integer(4),parameter  :: n2=80
      integer(4),parameter  :: n3=80
      real(8)               :: data(n1,n2,n3)
      real(8)               :: rtox(3,3)
      real(8)               :: vol
      integer(4)            :: nrview,nrviewx
      real(8)   ,allocatable:: rview(:,:)
      real(8)   ,allocatable:: zview(:)
      real(8)               :: pi
!     **************************************************************************
!
!     ==========================================================================
!     ==  interpolate wave function onto viewing grid                         ==
!     ==========================================================================
      call mapfieldtogrid(RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX,n1,n2,n3,data)
 !
!     ==========================================================================
!     ==  map atoms into view box                                             ==
!     ==========================================================================
      pi=4.d0*atan(1.d0)
      call gbass(box,rtox,vol)      
      nrviewx=2*nint(vol/(4.d0*pi/3.d0*1.4d0**3)) ! estimate max number of atoms
                                                ! inside the viewbox
      allocate(rview(3,nrviewx)) ! positions
      allocate(zview(nrviewx))   ! atomic numbers
      call MAKEatomsintobox(NAT,Z,R,RBAS,ORIGIN,BOX,nrviewx,nrview,rview,zview)
!
!     ==========================================================================
!     ==  write cube file                                                     ==
!     ==========================================================================
      call WRITECUBEFILE(NFIL,Nrview,Zview,Rview,ORIGIN,BOX,N1,N2,N3,DATA)
!
      deallocate(rview)
      deallocate(zview)
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEvrml(NFIL,NAT,Z,R,RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX &
     &                   ,tplane,planer0,planevec)
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
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! atomic positions
      REAL(8)   ,INTENT(IN) :: Rbas(3,3) !lattice vectors of the periodic data
      INTEGER(4),INTENT(IN) :: Nr1,Nr2,Nr3 ! # grid points on the periodic grid
      REAL(8)   ,INTENT(IN) :: wave(Nr1,Nr2,Nr3) ! periodic density data
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! corner of the new grid
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! vectors spanning the new grid
      logical(4),intent(in) :: tplane
      REAL(8)   ,INTENT(IN) :: planer0(3)   ! 
      REAL(8)   ,INTENT(IN) :: planevec(3,2)  ! 
      integer(4),parameter  :: n1=160      ! displacement
      integer(4),parameter  :: n2=160
      integer(4),parameter  :: n3=1
      real(8)               :: data(n1,n2,n3)
      real(8)               :: rtox(3,3)
      real(8)               :: vol
      integer(4)            :: nrview,nrviewx
      real(8)   ,allocatable:: rview(:,:)
      real(8)   ,allocatable:: zview(:)
      real(8)               :: pi
      real(8)               :: planebox(3,3)
!     **************************************************************************
                          call trace$push('MAKEvrml')
!
!     ==========================================================================
!     ==  interpolate wave function onto viewing grid                         ==
!     ==========================================================================
      if(tplane) then
        planebox(:,:)=0.d0
        planebox(:,:2)=planevec(:,:)
        call mapfieldtogrid(RBAS,NR1,NR2,NR3,WAVE &
     &                                     ,planer0,planeBOX,n1,n2,n3,data)
      end if
!
!     ==========================================================================
!     ==  map atoms into view box                                             ==
!     ==========================================================================
      pi=4.d0*atan(1.d0)
      call gbass(box,rtox,vol)      
      nrviewx=2*nint(vol/(4.d0*pi/3.d0*1.4d0**3)) ! estimate max number of atoms
                                                ! inside the viewbox
      allocate(rview(3,nrviewx)) ! positions
      allocate(zview(nrviewx))   ! atomic numbers
      call MAKEatomsintobox(NAT,Z,R,RBAS,ORIGIN,BOX,nrviewx,nrview,rview,zview)
!
!     ==========================================================================
!     ==  write cube file                                                     ==
!     ==========================================================================
      write(nfil,fmt='("#VRML V2.0 utf8")')
      write(nfil,*)'Background {'
      write(nfil,*)'  skyColor 0.933 0.863 0.51'
      write(nfil,*)'  groundColor 0.933 0.863 0.51'
      write(nfil,*)'}'
      write(nfil,*)'Viewpoint {'
      write(nfil,*)'  position 0 0 40 '
      write(nfil,*)'  fieldOfView 0.8 '
      write(nfil,*)'} # close Viewpoint'
      call vrml$addballstick(nfil,Nrview,Zview,Rview)
      if(tplane) then
        call vrml$addrubbersheet(nfil,planer0,planebox,n1,n2,data)
      end if
!
      deallocate(rview)
      deallocate(zview)

                          call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEvrmlC(NFIL,NAT,Z,R,RBAS,xk,NR1,NR2,NR3,cWAVE,ORIGIN,BOX &
     &                   ,tplane,planer0,planevec)
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
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! atomic positions
      REAL(8)   ,INTENT(IN) :: xk(3)     ! k-point in relative coordinates
      REAL(8)   ,INTENT(IN) :: Rbas(3,3) ! lattice vectors of the periodic data
      INTEGER(4),INTENT(IN) :: Nr1,Nr2,Nr3 ! # grid points on the periodic grid
      complex(8),INTENT(IN) :: cwave(Nr1,Nr2,Nr3) ! periodic density data
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! corner of the new grid
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! vectors spanning the new grid
      logical(4),intent(in) :: tplane
      REAL(8)   ,INTENT(IN) :: planer0(3) ! 
      REAL(8)   ,INTENT(IN) :: planevec(3,2)  ! 
      integer(4),parameter  :: n1=80      ! displacement
      integer(4),parameter  :: n2=80
      integer(4),parameter  :: n3=1
      complex(8)            :: cdata(n1,n2,n3)
      real(8)               :: data(n1,n2,n3)
      real(8)               :: rtox(3,3)
      real(8)               :: vol
      integer(4)            :: nrview,nrviewx
      real(8)   ,allocatable:: rview(:,:)
      real(8)   ,allocatable:: zview(:)
      real(8)               :: pi
      real(8)               :: planebox(3,3)
!     **************************************************************************
                          call trace$push('MAKEvrml')
!
!     ==========================================================================
!     ==  interpolate wave function onto viewing grid                         ==
!     ==========================================================================
      if(tplane) then
        planebox(:,:)=0.d0
        planebox(:,:2)=planevec(:,:)
        call mapfieldtogridc(RBAS,xk,NR1,NR2,NR3,cWAVE &
     &                                     ,planer0,planeBOX,n1,n2,n3,cdata)
      end if
!
!     ==========================================================================
!     ==  map atoms into view box                                             ==
!     ==========================================================================
      pi=4.d0*atan(1.d0)
      call gbass(box,rtox,vol)      
      nrviewx=2*nint(vol/(4.d0*pi/3.d0*1.4d0**3)) ! estimate max number of atoms
                                                ! inside the viewbox
      allocate(rview(3,nrviewx)) ! positions
      allocate(zview(nrviewx))   ! atomic numbers
      call MAKEatomsintobox(NAT,Z,R,RBAS,ORIGIN,BOX,nrviewx,nrview,rview,zview)
!
!     ==========================================================================
!     ==  write cube file                                                     ==
!     ==========================================================================
      write(nfil,fmt='("#VRML V2.0 utf8")')
      write(nfil,*)'Background {'
      write(nfil,*)'  skyColor 0.933 0.863 0.51'
      write(nfil,*)'  groundColor 0.933 0.863 0.51'
      write(nfil,*)'}'
      write(nfil,*)'Viewpoint {'
      write(nfil,*)'  position 0 0 40 '
      write(nfil,*)'  fieldOfView 0.8 '
      write(nfil,*)'} # close Viewpoint'
      call vrml$addballstick(nfil,Nrview,Zview,Rview)
      if(tplane) then
        data=real(cdata)
        call vrml$addrubbersheet(nfil,planer0,planebox,n1,n2,data)
      end if
!
      deallocate(rview)
      deallocate(zview)

                          call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEgnu(NFIL,NAT,Z,R,RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! atomic positions
      REAL(8)   ,INTENT(IN) :: Rbas(3,3) !lattice vectors of the periodic data
      INTEGER(4),INTENT(IN) :: Nr1,Nr2,Nr3 ! # grid points on the periodic grid
      REAL(8)   ,INTENT(IN) :: wave(Nr1,Nr2,Nr3) ! periodic density data
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! corner of the new grid
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! vectors spanning the new grid
      integer(4),parameter  :: n1=30      ! displacement
      integer(4),parameter  :: n2=30
      integer(4),parameter  :: n3=1
      real(8)               :: data(n1,n2,n3)
      integer(4)            :: i,j
!     **************************************************************************
                              call trace$push('makegnu')
!
!     ==========================================================================
!     ==  interpolate wave function onto viewing grid                         ==
!     ==========================================================================
      call mapfieldtogrid(RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX,n1,n2,n3,data)

open(99,file='gnudata')
do i=1,n1
  do j=1,n2
    write(99,*)(real(i-1)/real(n1-1)-0.5d0)*sqrt(sum(box(:,1)**2)) &
   &          ,(real(j-1)/real(n2-1)-0.5d0)*sqrt(sum(box(:,2)**2)) &
   &          ,data(i,j,1)
  enddo
enddo
close(99)
      open(99,file='rubbersheet.gnu')
      write(99,*)'#'                                                     
      write(99,*)'#================================================='      
      write(99,*)'# data section to be changed by the user        =='      
      write(99,*)'#================================================='      
      write(99,*)'xmin=',-0.5d0*sqrt(sum(box(:,1)**2)) 
      write(99,*)'xmax=',+0.5d0*sqrt(sum(box(:,1)**2)) 
      write(99,*)'ymin=',-0.5d0*sqrt(sum(box(:,2)**2)) 
      write(99,*)'ymax=',+0.5d0*sqrt(sum(box(:,2)**2)) 
      write(99,*)'zmin=',minval(data)
      write(99,*)'zmax=',maxval(data)
      write(99,*)'rot_x=30'
      write(99,*)'rot_z=20'
      write(99,*)'scale=1.0'
      write(99,*)'scale_z=1.'
      write(99,*)'infile="gnudata"'                                  
      write(99,*)'#'                                                     
      write(99,*)'#================================================='      
      write(99,*)'# define line styles to be used with ls in splot=='      
      write(99,*)'#================================================='      
      write(99,*)'# define a linestyle for splot (used with ls)'         
      write(99,*)'set style line 1 lt 1 lc rgb "black" lw 1'             
      write(99,*)'# map hight values to colors'                          
      write(99,*)'set palette rgbformula 22,13,-31'                      
      write(99,*)'#'
      write(99,*)'#================================================='
      write(99,*)'# data related statements                       =='       
      write(99,*)'#================================================='
      write(99,*)'set xrange [xmin:xmax]'                                
      write(99,*)'set yrange [ymin:ymax]'                                
      write(99,*)'set zrange [zmin:zmax]'                                
      write(99,*)'# sample input data onto a 60x60 grid'                 
      write(99,*)'set dgrid3d  60,60,1'                                  
      write(99,*)'#'                                                     
      write(99,*)'#================================================='
      write(99,*)'# surface plot                                  =='      
      write(99,*)'#================================================='
      write(99,*)'set surface'                                           
      write(99,*)'set hidden3d'                                          
      write(99,*)'set data style lines'                                  
      write(99,*)'# places the zero of the z-axis into the xy plane'     
      write(99,*)'set xyplane at 0.'                                     
      write(99,*)'# remove axes'                                         
      write(99,*)'unset border'                                          
      write(99,*)'#  remove tics from the axes'                          
      write(99,*)'unset xtics'                                           
      write(99,*)'unset ytics'                                           
      write(99,*)'unset ztics'                                           
      write(99,*)'# no title written'                                    
      write(99,*)'set key off'                                           
      write(99,*)'# place contours onto the surface'                     
      write(99,*)'set pm3d explicit hidden3d 1'                          
      write(99,*)'#'                                                   
      write(99,*)'#  angle, angle, overall scale, scale z-axis'        
      write(99,*)'set view rot_x,rot_z,scale,scale_z'                  
      write(99,*)'#'                                                   
      write(99,*)'splot  infile with pm3d'                             
      close(99)
                              call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE mapfieldtogrid(RBAS,NR1,NR2,NR3,WAVE,ORIGIN,BOX,n1,n2,n3,data)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: Rbas(3,3) !lattice vectors of the periodic data
      INTEGER(4),INTENT(IN) :: Nr1,Nr2,Nr3 ! # grid points on the periodic grid
      REAL(8)   ,INTENT(IN) :: wave(Nr1,Nr2,Nr3) ! periodic density data
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! corner of the new grid
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! vectors spanning the new grid
      integer(4),intent(in) :: n1,n2,n3
      real(8)   ,intent(out):: data(n1,n2,n3)
      real(8)               :: pos(3)
      real(8)               :: dt1(3),dt2(3),dt3(3)
      real(8)               :: xtor(3,3)
      real(8)               :: rtox(3,3)
      real(8)               :: svar1,svar2
      real(8)               :: f00,f01,f10,f11,g0,g1
      real(8)               :: xdis(3)
      INTEGER(4)            :: i1m,i1p,i2m,i2p,i3m,i3p
      INTEGER(4)            :: I,J,K
!     **************************************************************************
                                call trace$push('mapfieldtogrid')
      xtor(:,1)=rbas(:,1)/real(nr1,kind=8)
      xtor(:,2)=rbas(:,2)/real(nr2,kind=8)
      xtor(:,3)=rbas(:,3)/real(nr3,kind=8)
      call lib$invertr8(3,xtor,rtox)
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE ONTO NEW GRID                                      ==
!     ==========================================================================
      DT1(:)=BOX(:,1)/REAL(N1-1,KIND=8)
      DT2(:)=BOX(:,2)/REAL(N2-1,KIND=8)
      if(n3.gt.1) then
        DT3(:)=BOX(:,3)/REAL(N3-1,KIND=8)
      else                   ! n3=1 is for a 2d grid
        DT3(:)=0.d0
      end if
      DO I=1,N1
        DO J=1,N2
          POS(:)=ORIGIN(:)+DT1*REAL(I-1,KIND=8)+DT2*REAL(J-1,KIND=8)-DT3(:)
          DO K=1,N3
            POS(:)=POS(:)+DT3(:)  ! POSITION ON THE NEW GRID
!           == map position into relative coordinates and into 1st unit cell
            XDIS=MATMUL(RTOX,POS)
            xdis(1)=MODULO(XDIS(1),real(NR1,kind=8))
            xdis(2)=MODULO(XDIS(2),real(NR2,kind=8))
            xdis(3)=MODULO(XDIS(3),real(NR3,kind=8))
!           == avoid compiler bug. for very small negative values a result  ==
!           == equal to the right boundary can occur
            if(xdis(1).eq.real(nr1,kind=8))xdis(1)=xdis(1)-real(NR1,kind=8)
            if(xdis(2).eq.real(nr2,kind=8))xdis(2)=xdis(2)-real(NR2,kind=8)
            if(xdis(3).eq.real(nr3,kind=8))xdis(3)=xdis(3)-real(NR3,kind=8)
            I1M=1+int(XDIS(1))
            I1P=1+MODULO(I1M,NR1)
            I2M=1+int(XDIS(2))
            I2P=1+MODULO(I2M,NR2)
            I3M=1+int(XDIS(3))
            I3P=1+MODULO(I3M,NR3)
            XDIS(1)=modulo(XDIS(1),1.d0)
            XDIS(2)=modulo(XDIS(2),1.d0)
            XDIS(3)=modulo(XDIS(3),1.d0)
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
            DATA(I,j,k)=G0*SVAR1+G1*SVAR2
          ENDDO
        ENDDO
      ENDDO
                                call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE mapfieldtogridc(RBAS,xk,NR1,NR2,NR3,WAVE &
     &                           ,ORIGIN,BOX,n1,n2,n3,data)
!     **************************************************************************
!     ** INTERPOLATES PERIODIC DENSITY DATA ON AN ARBITRARY GRID ONTO         **
!     ** AN INDEPENDENT GRID SPECIFIED BY ORIGIN AND BOX.                     **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: Rbas(3,3) !lattice vectors of the periodic data
      REAL(8)   ,INTENT(IN) :: xk(3)     !k-point in relative coordinates
      INTEGER(4),INTENT(IN) :: Nr1,Nr2,Nr3 ! # grid points on the periodic grid
      complex(8),INTENT(IN) :: wave(Nr1,Nr2,Nr3) ! periodic density data
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)  ! corner of the new grid
      REAL(8)   ,INTENT(IN) :: BOX(3,3)   ! vectors spanning the new grid
      integer(4),intent(in) :: n1,n2,n3
      complex(8),intent(out):: data(n1,n2,n3)
      real(8)               :: pos(3)
      real(8)               :: dt1(3),dt2(3),dt3(3)
      real(8)               :: xtor(3,3)
      real(8)               :: rtox(3,3)
      real(8)               :: svar1,svar2
      complex(8)            :: csvar1,csvar2
      complex(8)            :: eik1m,eik1p,eik2m,eik2p,eik3m,eik3p
      complex(8)            :: f00,f01,f10,f11,g0,g1
      real(8)               :: xdis(3),xdis0(3)
      INTEGER(4)            :: i1m,i1p,i2m,i2p,i3m,i3p
      INTEGER(4)            :: I,J,K
      complex(8),parameter  :: ci=(0.d0,1.d0)
      real(8)               :: pi
      complex(8)            :: c2pii
!     **************************************************************************
                                call trace$push('mapfieldtogridc')
      pi=4.d0*atan(1.d0)
      c2pii=(0.d0,1.d0)*2.d0*pi
      xtor(:,1)=rbas(:,1)/real(nr1,kind=8)
      xtor(:,2)=rbas(:,2)/real(nr2,kind=8)
      xtor(:,3)=rbas(:,3)/real(nr3,kind=8)
      call lib$invertr8(3,xtor,rtox)
!
!     ==========================================================================
!     ==  INTERPOLATE WAVE ONTO NEW GRID                                      ==
!     ==========================================================================
      DT1(:)=BOX(:,1)/REAL(N1-1,KIND=8)
      DT2(:)=BOX(:,2)/REAL(N2-1,KIND=8)
      if(n3.gt.1) then
        DT3(:)=BOX(:,3)/REAL(N3-1,KIND=8)
      else                   ! n3=1 is for a 2d grid
        DT3(:)=0.d0
      end if
      DO I=1,N1
        DO J=1,N2
          POS(:)=ORIGIN(:)+DT1*REAL(I-1,KIND=8)+DT2*REAL(J-1,KIND=8)-DT3(:)
          DO K=1,N3
            POS(:)=POS(:)+DT3(:)  ! POSITION ON THE NEW GRID
!           == map position into relative coordinates and into 1st unit cell ===
            XDIS0=MATMUL(RTOX,POS)
            xdis(1)=MODULO(XDIS0(1),real(NR1,kind=8))
            xdis(2)=MODULO(XDIS0(2),real(NR2,kind=8))
            xdis(3)=MODULO(XDIS0(3),real(NR3,kind=8))
!           == avoid compiler bug. for very small negative values a result  ====
!           == equal to the right boundary can occur ===========================
            if(xdis(1).eq.real(nr1,kind=8))xdis(1)=xdis(1)-real(NR1,kind=8)
            if(xdis(2).eq.real(nr2,kind=8))xdis(2)=xdis(2)-real(NR2,kind=8)
            if(xdis(3).eq.real(nr3,kind=8))xdis(3)=xdis(3)-real(NR3,kind=8)
            I1M=1+int(XDIS(1))
            I1P=1+MODULO(I1M,NR1)
            I2M=1+int(XDIS(2))
            I2P=1+MODULO(I2M,NR2)
            I3M=1+int(XDIS(3))
            I3P=1+MODULO(I3M,NR3)
            eik1m=exp(-c2pii*xk(1)*(xdis(1)-xdis0(1))/real(nr1,kind=8))
            eik2m=exp(-c2pii*xk(2)*(xdis(2)-xdis0(2))/real(nr2,kind=8))
            eik3m=exp(-c2pii*xk(3)*(xdis(3)-xdis0(3))/real(nr3,kind=8))
            eik1p=eik1m*exp(-c2pii*xk(1)*real(i1p-1-i1m)/real(nr1,kind=8))
            eik2p=eik2m*exp(-c2pii*xk(2)*real(i2p-1-i2m)/real(nr2,kind=8))
            eik3p=eik3m*exp(-c2pii*xk(3)*real(i3p-1-i3m)/real(nr3,kind=8))
!
            XDIS(1)=modulo(XDIS(1),1.d0)
            XDIS(2)=modulo(XDIS(2),1.d0)
            XDIS(3)=modulo(XDIS(3),1.d0)
!           == INTERPOLATE ALONG THIRD DIRECTION ===============================
            SVAR2=XDIS(3)
            SVAR1=1.D0-SVAR2
            csvar1=svar1*eik3m
            csvar2=svar2*eik3p
            F00=WAVE(I1M,I2M,I3M)*cSVAR1+WAVE(I1M,I2M,I3P)*cSVAR2
            F10=WAVE(I1P,I2M,I3M)*cSVAR1+WAVE(I1P,I2M,I3P)*cSVAR2
            F01=WAVE(I1M,I2P,I3M)*cSVAR1+WAVE(I1M,I2P,I3P)*cSVAR2
            F11=WAVE(I1P,I2P,I3M)*cSVAR1+WAVE(I1P,I2P,I3P)*cSVAR2
!           == INTERPOLATE ALONG SECOND DIRECTION ==============================
            SVAR2=XDIS(2)
            SVAR1=1.D0-SVAR2
            csvar1=svar1*eik2m
            csvar2=svar2*eik2p
            G0=F00*cSVAR1+F01*cSVAR2
            G1=F10*cSVAR1+F11*cSVAR2
!           == INTERPOLATE ALONG SECOND DIRECTION ==============================
            SVAR2=XDIS(1)
            SVAR1=1.D0-SVAR2
            csvar1=svar1*eik1m
            csvar2=svar2*eik1p
            DATA(I,j,k)=G0*cSVAR1+G1*cSVAR2
          ENDDO
        ENDDO
      ENDDO
                                call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MAKEatomsintobox(NAT,Z,R,RBAS,ORIGIN,BOX &
     &                           ,nrviewx,nrview,rview,zview)
!     **************************************************************************
!     ** maps periodic lattice into viewbox                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    ! ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)  ! atomic positions
      REAL(8)   ,INTENT(IN) :: Rbas(3,3) ! lattice vectors of the periodic data
      REAL(8)   ,INTENT(IN) :: ORIGIN(3) ! corner of the new grid
      REAL(8)   ,INTENT(IN) :: BOX(3,3)  ! vectors spanning the new grid
      INTEGER(4),INTENT(IN) :: nrviewx   ! x#(atoms in viewbox)
      INTEGER(4),INTENT(out):: nrview    ! #(atoms in viewbox)
      real(8)   ,intent(out):: rview(3,nrviewx)
      real(8)   ,intent(out):: zview(nrviewx)
      real(8)               :: pos(3)
      real(8)               :: rtox(3,3)
      real(8)               :: xdis(3)
      real(8)               :: t1(3),t2(3),t3(3)
      INTEGER(4)            :: it1,it2,it3,iat
      INTEGER(4)            :: min1,max1,min2,max2,min3,max3
!     **************************************************************************
!
!     ==========================================================================
!     ==  map atoms into the viewbox                                          ==
!     ==========================================================================
      call lib$invertr8(3,box,rtox)
      min1=-5
      max1=5
      min2=-5
      max2=5
      min3=-5
      max3=5
      nrview=0
      do it1=min1,max1
        t1(:)=rbas(:,1)*real(it1,kind=8)
        do it2=min2,max2
          t2(:)=rbas(:,2)*real(it2,kind=8)
          do it3=min3,max3
            t3(:)=rbas(:,3)*real(it3,kind=8)
            do iat=1,nat
              pos(:)=r(:,iat)+t1+t2+t3
              xdis(:)=matmul(rtox,pos-origin(:))
              if(xdis(1).lt.-0.05d0.or.xdis(1).gt.1.05d0) cycle
              if(xdis(2).lt.-0.05d0.or.xdis(2).gt.1.05d0) cycle
              if(xdis(3).lt.-0.05d0.or.xdis(3).gt.1.05d0) cycle
              nrview=nrview+1
              IF(NRVIEW.GT.NRVIEWX) THEN
                CALL ERROR$MSG('NUMBER OF ATOMS IN THE BOX OUT OF RANGE')
                CALL ERROR$I4VAL('NRVIEWX',NRVIEWX)
                CALL ERROR$I4VAL('NRVIEW',NRVIEW)
                CALL ERROR$STOP('MAKECUBE')
              END IF
              rview(:,nrview)=pos(:)
              zview(nrview)=z(iat)
            enddo
          enddo
        enddo
      enddo
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITECUBEFILE(NFIL,NAT,Z,R,ORIGIN,BOX,N1,N2,N3,DATA)
!     **************************************************************************
!     **  write a Gaussian Cube file (extension .cub) with voluminetric data  **
!     **                                                                      **
!     ** remark:                                                              **
!     ** units written are abohr, consistent with avogadro's implementation.  **
!     ** the specs require N1,N2,N3 to be multiplied by -1 if abohr are used  **
!     ** and angstrom is the unit if they are positive.                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)
      REAL(8)   ,INTENT(IN) :: BOX(3,3)
      INTEGER(4),INTENT(IN) :: N1,N2,N3
      REAL(8)   ,INTENT(IN) :: DATA(N1,N2,N3)
      REAL(8)               :: ANGSTROM
      REAL(8)               :: scale
      INTEGER(4)            :: IAT,I,J,K
!     **************************************************************************
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      scale=1.d0
!      SCALE=1.D0/ANGSTROM
      WRITE(NFIL,FMT='("CP-PAW CUBE FILE")')
      WRITE(NFIL,FMT='("NOCHN KOMMENTAR")')
      WRITE(NFIL,FMT='(I5,3F12.6)')NAT,ORIGIN*scale
      WRITE(NFIL,FMT='(I5,3F12.6)')-N1,BOX(:,1)/REAL(N1-1,KIND=8)*scale
      WRITE(NFIL,FMT='(I5,3F12.6)')-N2,BOX(:,2)/REAL(N2-1,KIND=8)*scale
      WRITE(NFIL,FMT='(I5,3F12.6)')-N3,BOX(:,3)/REAL(N3-1,KIND=8)*scale
      DO IAT=1,NAT
        WRITE(NFIL,FMT='(I5,4F12.6)')NINT(Z(IAT)),0.D0,R(:,IAT)*scale
      ENDDO  
      WRITE(NFIL,FMT='(6(E12.6," "))')(((DATA(I,J,K),K=1,N3),J=1,N2),I=1,N1)
!!$do k=1,n3
!!$  scale=0.d0
!!$  do i=1,n1
!!$    do j=1,n2
!!$      scale=scale+data(i,j,k)**2
!!$    enddo
!!$  enddo
!!$print*,'scale ',origin(:)+box(:,3)/real(n3,kind=8)*real(k-1,kind=8),scale
!!$enddo
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE vrml$addrubbersheet(nfil,origin,plane,n1,n2,data)
!     **************************************************************************
!     **  add a ball-stick model to a vrml scene                              **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      integer(4),intent(IN) :: nfil
      integer(4),intent(IN) :: n1
      integer(4),intent(IN) :: n2
      real(8)   ,intent(in) :: origin(3)
      real(8)   ,intent(in) :: plane(3,2)
      real(8)   ,intent(in) :: data(n1,n2)
      real(8)   ,parameter  :: scale=10.d0
      real(8)               :: rot(3,3)
      real(8)               :: axis(3)
      real(8)               :: angle
      real(8)               :: d1,d2
      real(8)               :: svar1,svar2,x
      real(8)               :: color(3,n1,n2)
      integer(4)            :: i,j
      integer(4)            :: r,g,b
      logical(4)            :: tinv
!     **************************************************************************
                             call trace$push('vrml$addrubbersheet')
      d1=sqrt(sum(plane(:,1)**2))
      d2=sqrt(sum(plane(:,2)**2))
!     == construct rotation matrix =============================================
!     == first direction is mapped onto x-axis
!     == second direction is mapped onto z-axis
      rot(:,1)=plane(:,1)/d1
      rot(:,3)=plane(:,2)/d2
      rot(1,2)=rot(2,3)*rot(3,1)-rot(3,3)*rot(2,1)
      rot(2,2)=rot(3,3)*rot(1,1)-rot(1,3)*rot(3,1)
      rot(3,2)=rot(1,3)*rot(2,1)-rot(2,3)*rot(1,1)
!     == extract rotation angle ================================================
!     == see http://www.euclideanspace.com/maths/geometry/rotations           ==
!     ==                                           /conversions/matrixToAngle ==
      call rotation$matrixtoangle(rot,axis,tinv)
      angle=sqrt(sum(axis**2))
      axis=axis/angle
!
!     == create color ===== ====================================================
      color(:,:,:)=0.d0
      svar1=minval(data)
      svar2=maxval(data)
      do i=1,n1
        do j=1,n1
          x=(data(i,j)-svar1)/(svar2-svar1)
          call rgb$rainbow(x,r,g,b)
          color(1,i,j)=real(r)/real(255)
          color(2,i,j)=real(g)/real(255)
          color(3,i,j)=real(b)/real(255)
        enddo
      enddo
!
!     == create write model to scene ===========================================
      write(nfil,*)'Transform{translation ',Origin,' rotation ',axis,angle
      write(nfil,*)'children['
!
!     == create rubbersheet ====================================================
      write(nfil,*)'Shape{'
      write(nfil,*)'  appearance Appearance {'
      write(nfil,*)'    material Material {'
      write(nfil,*)'      diffuseColor 0.8 0.4 0.0'
      write(nfil,*)'      shininess 1'
      write(nfil,*)'    } # close Material'
      write(nfil,*)'  }   # close Appearance'
      write(nfil,*)'  geometry ElevationGrid {'
      write(nfil,*)'    xDimension ',n1
      write(nfil,*)'    zDimension ',n2
      write(nfil,*)'    xSpacing ',d1/real(n1-1)
      write(nfil,*)'    zSpacing ',d2/real(n2-1)
      write(nfil,*)'    solid FALSE # makes sheet visible from both sides'   
      write(nfil,*)'    creaseAngle 3 # smoothening'
      write(nfil,*)'    color Color {color ['
      write(nfil,fmt='(5(3f10.5," , "))')color
      write(nfil,*)'      ]  # end of color'
      write(nfil,*)'    }    # end of Color'
      write(nfil,*)'    height ['
      write(nfil,fmt='(20f10.5)')-data*scale
      write(nfil,*)'    ]}}]}'  !end of height,geometry,Shape,children,Transform
                             call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE vrml$addballstick(nfil,NAT,Z,R)
!     **************************************************************************
!     **  add a ball-stick model to a vrml scene                              **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: nfil
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      real(8)   ,parameter  :: scaleball=0.6d0
      real(8)   ,parameter  :: bonddiameter=0.2d0
      real(8)   ,parameter  :: bondlengthcutoff=1.3d0
      real(8)   ,parameter  :: bondcolor(3)=(/0.3d0,0.3d0,0.3d0/)
      integer(4)            :: iat,iat1,iat2
      real(8)               :: rad,rad1,rad2
      real(8)               :: color(3)
      integer(4)            :: ivec(3)
      real(8)               :: dr(3),r1(3),r2(3)
      real(8)               :: dis
      real(8)               :: axis(3),angle
!     **************************************************************************
!
!      == create balls =========================================================
       do iat=1,nat
         CALL PERIODICTABLE$GET(NINT(Z(IAT)),'R(COV)',RAD)
         CALL ATOMCOLOR(NINT(Z(IAT)),IVEC)
         COLOR(:)=REAL(IVEC,KIND=8)/255.D0
         write(nfil,*)'Transform{'
         write(nfil,*)'  translation ',r(:,iat)
         write(nfil,*)'  children['
         write(nfil,*)'    Shape{'
         write(nfil,*)'      appearance Appearance {'
         write(nfil,*)'        material Material {'
         write(nfil,*)'          diffuseColor ',color
         write(nfil,*)'          shininess 1'
         write(nfil,*)'          transparency 0.5'
         write(nfil,*)'        } # close Material'
         write(nfil,*)'      } # close appearance'
         write(nfil,*)'      geometry Sphere{ radius ',rad*scaleball,'}'
         write(nfil,*)'    }   # close Shape'
         write(nfil,*)'  ]     # close children'
         write(nfil,*)'}       # close Transform'
       enddo
!
!      == create bonds =========================================================
       do iat1=1,nat
         CALL PERIODICTABLE$GET(NINT(Z(IAT1)),'R(COV)',RAD1)
         do iat2=iat1+1,nat
           CALL PERIODICTABLE$GET(NINT(Z(IAT2)),'R(COV)',RAD2)
           dr(:)=r(:,iat2)-r(:,iat1)
           dis=sqrt(sum(dr**2))
           if(dis/(rad1+rad2).gt.bondlengthcutoff) cycle           
           dr(:)=dr(:)/dis
           r1(:)=r(:,iat1)+rad1*scaleball*1.001d0*dr(:)
           r2(:)=r(:,iat2)-rad2*scaleball*1.001d0*dr(:)
           angle=-acos(dr(2))   ! original orientation (0,1,0)
           axis(1)=-dr(3)
           axis(2)=0.d0
           axis(3)=dr(1)
           axis=axis/sqrt(sum(axis**2))
           write(nfil,*)'Transform{'
           write(nfil,*)'  translation ',0.5d0*(r1(:)+r2(:))
           write(nfil,*)'  rotation    ',axis,angle
           write(nfil,*)'  children ['
           write(nfil,*)'    Shape{'
           write(nfil,*)'      appearance Appearance {'
           write(nfil,*)'        material Material {'
           write(nfil,*)'          diffuseColor ',bondcolor
           write(nfil,*)'          shininess 1'
           write(nfil,*)'          transparency 0.7'
           write(nfil,*)'        } # close Material'
           write(nfil,*)'      }   # close Appearance'
           write(nfil,*)'      geometry Cylinder{'
           write(nfil,*)'        radius ',bonddiameter
           write(nfil,*)'        height ',sqrt(sum((r2-r1)**2))
           write(nfil,*)'      }   # close geometry'
           write(nfil,*)'    }     # close Shape'
           write(nfil,*)'  ]       # close children'
           write(nfil,*)'}         # close Transform'
         enddo
       enddo
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE WRITEgnuFILE(NFIL,l1,l2,n1,n2,data)
!!$!     **************************************************************************
!!$!     **  write a gnuplot file                                                **
!!$!     **                                                                      **
!!$!     ** remark:                                                              **
!!$!     ** units written are abohr, consistent with avogadro's implementation.  **
!!$!     ** the specs require N1,N2,N3 to be multiplied by -1 if abohr are used  **
!!$!     ** and angstrom is the unit if they are positive.                       **
!!$!     **************************************************************************
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: NFIL
!!$      INTEGER(4),INTENT(IN) :: N1,N2
!!$      REAL(8)   ,INTENT(IN) :: l1,l2
!!$      REAL(8)   ,INTENT(IN) :: DATA(N1,N2)
!!$      REAL(8)               :: ANGSTROM
!!$      INTEGER(4)            :: I,J,K
!!$!     **************************************************************************
!!$      write(nfil,fmt='("set style line 1 lt 1 lc rgb "black" lw 1")')
!!$      write(nfil,fmt='("set palette rgbformula 22,13,-31")')
!!$      write(nfil,fmt='("set xrange [",f10.5,":",f10.5]")')0.d0,l1
!!$      write(nfil,fmt='("set yrange [",f10.5,":",f10.5"]")')0.d0,l2
!!$      write(nfil,fmt='("set zrange [-0.015:0.3]")')
!!$      write(nfil,fmt='("set dgrid3d  60,60,1")')
!!$      write(nfil,fmt='("set cntrparam levels incremental -0.2,0.01,0.2")')
!!$      write(nfil,fmt='("set surface")')
!!$      write(nfil,fmt='("set hidden3d")')
!!$      write(nfil,fmt='("set data style lines")')
!!$      write(nfil,fmt='("set view 30, 20, 1.0, 2.5")')
!!$      write(nfil,fmt='("set xyplane at 0.")')
!!$      write(nfil,fmt='("unset border")')
!!$      write(nfil,fmt='("unset xtics")')
!!$      write(nfil,fmt='("unset ytics"0')
!!$      write(nfil,fmt='("unset ztics")')
!!$      write(nfil,fmt='("set key off"0')
!!$      write(nfil,fmt='("set pm3d explicit hidden3d 1")')
!!$      write(nfil,fmt='("splot "ntb2d_iat1s.dat" with pm3d")')
!!$      RETURN
!!$      END
