! THESE ROUTINES ALLOW TO WRITE AND READ THE OLD FILE FORMATS OF THE
! RESTART FILE. AFTER READING HOWEVER THE WRITE ROUTINE IS NOT ATTACHED.
!
! THE RESTART FILES ARE CONVERT.RSTRTIN AND CONVERT.RSTRTOUT IN THE
! CURRENT DIRECTORY
!
!     ..................................................................
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER(4)      :: NFIL
      INTEGER(4)      :: NFI,NSP1,NB1,NKPT1,NSPIN1
      CHARACTER(80)   :: TEXT
      CHARACTER(32)   :: STYLE      
      REAL(8)         :: DELTAT=10.D0
      INTEGER(4)      :: IOS
      CHARACTER(256)  :: INFILE,OUTFILE
!     *****************************************************************
!
!     =================================================================
!     == CONNECT OLD AN NEW FILES                                    ==
!     =================================================================
      CALL FILEHANDLER$SETROOT('./CONVERT')
!
!     ==  RESTART_OUT FILE =============================================
      ID=+'RESTART_IN'
      CALL FILEHANDLER$SETFILE(ID,T,-'.RSTRT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  RESTART_OUT FILE =============================================
      ID=+'RESTART_OUT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.RSTRT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     =================================================================
!     == DETERMINE STYLE                                             ==
!     =================================================================
      CALL FILEHANDLER$UNIT('RESTART_IN',NFIL)
      REWIND NFIL
      READ(NFIL,END=9000,IOSTAT=IOS)NFI,NSP1,NB1,NKPT1,NSPIN1,TEXT
      IF(TEXT.EQ.'TESTING OF ROUTINE RDWRTE FOR IO') THEN
        STYLE='STYLE1'
      ELSE
!       TEXT='STYLE(19JAN1995)'
        STYLE='STYLE2'
      END IF
!
!     =================================================================
!     == NOW READ 
!     =================================================================
      IF(STYLE(1:6).EQ.'STYLE1') THEN
        CALL READSTYLE1(DELTAT)
      ELSE IF(STYLE(1:6).EQ.'STYLE2') THEN
        CALL READSTYLE2(DELTAT)
      END IF
                        CALL TRACE$POP
      RETURN
!
 9000 CONTINUE
      CALL IOSTATMESSAGE(IOS,NFIL)
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$STOP('WAVEIO')
!
      END
!
!     ..................................................................
      SUBROUTINE READSTYLE2(DELTAT)
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: DELTAT
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NFI
      INTEGER(4)             :: NSP
      CHARACTER(80)          :: TEXT
      INTEGER(4),ALLOCATABLE :: NA(:) !(NSP)
      INTEGER(4),ALLOCATABLE :: NGWG(:) !(NSP)
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NB
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: NSPIN
      INTEGER(4)             :: NGWGX
      REAL(8)   ,ALLOCATABLE :: RLAM0(:,:,:,:) !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: RLAMM(:,:,:,:) !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: H0(:,:,:,:)    !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: HM(:,:,:,:)    !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)        !(3,NAT)
      REAL(8)   ,ALLOCATABLE :: RM(:,:)        !(3,NAT)
      COMPLEX(8),ALLOCATABLE :: C0(:,:,:,:)    !(NG,NB,NKPT,NSPIN)
      COMPLEX(8),ALLOCATABLE :: CM(:,:,:,:)    !(NG,NB,NKPT,NSPIN)
      INTEGER(4)             :: IOS
      REAL(8)   ,ALLOCATABLE :: EIG(:,:,:)
      REAL(8)                :: X0,XM
      REAL(8)                :: XE0,XEM,XE2M
      REAL(8)                :: DUMMY
      INTEGER(4)             :: ISP,IKPT,IB,IB1,IB2,ISPIN,I,J,K
!         
!     ==================================================================
!     ==================================================================
!     ==  READ FROM FILE                                              ==
!     ==================================================================
!     ==================================================================
      CALL FILEHANDLER$UNIT('RESTART_IN',NFIL)
      REWIND NFIL
      READ(NFIL,END=9000,IOSTAT=IOS)NFI,NSP,NB,NKPT,NSPIN,TEXT
      ALLOCATE(NA(NSP))
      ALLOCATE(NGWG(NKPT))
      READ(NFIL,END=9000,IOSTAT=IOS)(NA(I),I=1,NSP),(NGWG(K),K=1,NKPT)
      NAT=0
      DO ISP=1,NSP
        NAT=NAT+NA(ISP)
      ENDDO
      NGWGX=0
      DO IKPT=1,NKPT
        NGWGX=MAX(NGWGX,NGWG(IKPT))
      ENDDO
!     
!     ===============================================================
!     ==   READ COORDINATES
!     ===============================================================
      ALLOCATE(R0(3,NAT))
      ALLOCATE(RM(3,NAT))
      READ(NFIL,END=9000,IOSTAT=IOS)((R0(I,J),RM(I,J),I=1,3),J=1,NAT)
!     
!     ===============================================================
!     ==   READ WAVE FUNCTIONS                                     ==
!     ===============================================================
      ALLOCATE(C0(NGWGX,NB,NKPT,NSPIN))
      ALLOCATE(CM(NGWGX,NB,NKPT,NSPIN))
      C0=(0.D0,0.D0)
      CM=(0.D0,0.D0)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NB
            READ(NFIL,END=9000,IOSTAT=IOS) &
     &        (C0(I,IB,IKPT,ISPIN),CM(I,IB,IKPT,ISPIN),I=1,NGWG(IKPT))
          ENDDO
        ENDDO
      ENDDO
!     
!     ===============================================================
!     ==   READ NOSE VARIABLE                                      ==
!     ===============================================================
      READ(NFIL,END=9000,IOSTAT=IOS)X0,XM,XE0,XEM,XE2M
!     
!     ===============================================================
!     ==   READ LAGRANGE PARAMETERS                                ==
!     ===============================================================
      ALLOCATE(RLAM0(NB,NB,NKPT,NSPIN))
      ALLOCATE(RLAMM(NB,NB,NKPT,NSPIN))
      RLAM0(:,:,:,:)=0.D0
      RLAMM(:,:,:,:)=0.D0
      READ(NFIL,END=1000)((((RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                      ,RLAMM(IB1,IB2,IKPT,ISPIN) &
     &                      ,DUMMY &
     &         ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
!     
!     ===============================================================
!     ==   READ HAMILTON MATRIX                                    ==
!     ===============================================================
      ALLOCATE(H0(NB,NB,NKPT,NSPIN))
      ALLOCATE(HM(NB,NB,NKPT,NSPIN))
      H0(:,:,:,:)=0.D0
      HM(:,:,:,:)=0.D0
      READ(NFIL,END=1000)((((H0(IB1,IB2,IKPT,ISPIN) &
     &                      ,HM(IB1,IB2,IKPT,ISPIN) &
     &         ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
!     
!     ===============================================================
!     ==   CONTINUE HERE (IF END ENCOUNTERED)                      ==
!     ===============================================================
 1000 CONTINUE
      CALL FILEHANDLER$CLOSE('RESTART_IN')
!     
!     ==================================================================
!     == TRANSFORM X(0),V(-1/2) INTO X(0),V(-1/2)                     ==
!     ==================================================================
      RM(:,:)=R0(:,:)-RM(:,:)*DELTAT
      CM(:,:,:,:)=C0(:,:,:,:)-CM(:,:,:,:)*DELTAT
      HM(:,:,:,:)=H0(:,:,:,:)-HM(:,:,:,:)*DELTAT
      RLAMM(:,:,:,:)=RLAM0(:,:,:,:)-RLAMM(:,:,:,:)*DELTAT
      CALL WAVES$SETR8A('RLAM(0)',NB*NB*NKPT*NSPIN,RLAM0)
      CALL WAVES$SETR8A('RLAM(-)',NB*NB*NKPT*NSPIN,RLAMM)
      ALLOCATE(EIG(NB,NKPT,NSPIN))
      DO IB=1,NB
        EIG(IB,:,:)=H0(IB,IB,:,:)
      ENDDO
                          CALL TRACE$POP
      RETURN
!       

 9000 CONTINUE
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$I4VAL('IOS ',IOS)
      CALL ERROR$STOP('RDWRTE')
      END
!
!     ..................................................................
      SUBROUTINE READSTYLE1(DELTAT)
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: DELTAT
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NFI
      INTEGER(4)             :: NSP,NB,NKPT,NSPIN
      CHARACTER(80)          :: TEXT
      INTEGER(4),ALLOCATABLE :: NA(:) !(NSP)
      INTEGER(4),ALLOCATABLE :: NGWG(:) !(NSP)
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NGWGX
      REAL(8)   ,ALLOCATABLE :: RLAM0(:,:,:,:) !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: RLAMM(:,:,:,:) !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: H0(:,:,:,:)    !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: HM(:,:,:,:)    !(NB,NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)        !(3,NAT)
      REAL(8)   ,ALLOCATABLE :: RM(:,:)        !(3,NAT)
      COMPLEX(8),ALLOCATABLE :: C0(:,:,:,:)    !(NG,NB,NKPT,NSPIN)
      COMPLEX(8),ALLOCATABLE :: CM(:,:,:,:)    !(NG,NB,NKPT,NSPIN)
      INTEGER(4)             :: IOS
      REAL(8)   ,ALLOCATABLE :: EIG(:,:,:)
      REAL(8)                :: X0,XM
      REAL(8)                :: XE0,XEM,XE2M
      REAL(8)                :: DUMMY
      INTEGER(4)             :: ISP,IKPT,IB,IB1,IB2,ISPIN,I,J,K
!         
!     ==================================================================
!     ==  READ FROM FILE                                              ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('RESTART_IN',NFIL)
      REWIND NFIL
      READ(NFIL,END=9000,IOSTAT=IOS)NFI,NSP,NB,NKPT,NSPIN,TEXT
      ALLOCATE(NA(NSP))
      ALLOCATE(NGWG(NKPT))
      READ(NFIL,END=9000,IOSTAT=IOS)(NA(I),I=1,NSP),(NGWG(K),K=1,NKPT)
      NAT=0
      DO ISP=1,NSP
        NAT=NAT+NA(ISP)
      ENDDO
      NGWGX=0
      DO IKPT=1,NKPT
        NGWGX=MAX(NGWGX,NGWG(IKPT))
      ENDDO
      ALLOCATE(R0(3,NAT))
      ALLOCATE(RM(3,NAT))
      ALLOCATE(C0(NGWGX,NB,NKPT,NSPIN))
      ALLOCATE(CM(NGWGX,NB,NKPT,NSPIN))
      ALLOCATE(RLAM0(NB,NB,NKPT,NSPIN))
      ALLOCATE(RLAMM(NB,NB,NKPT,NSPIN))
      ALLOCATE(H0(NB,NB,NKPT,NSPIN))
      ALLOCATE(HM(NB,NB,NKPT,NSPIN))
!     
!     ==   READ COORDINATES
      READ(NFIL,END=9000,IOSTAT=IOS) &
     &          ((R0(I,J),RM(I,J),I=1,3),J=1,NAT) &
     &         ,((((C0(I,J,K,ISPIN),CM(I,J,K,ISPIN) &
     &            ,I=1,NGWG(K)),J=1,NB),K=1,NKPT),ISPIN=1,NSPIN) &
     &         ,X0,XM,XE0,XEM,XE2M  
      READ(NFIL,END=1000)((((RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                      ,RLAMM(IB1,IB2,IKPT,ISPIN) &
     &                      ,DUMMY &
     &         ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
      READ(NFIL,END=1000)((((H0(IB1,IB2,IKPT,ISPIN) &
     &                      ,HM(IB1,IB2,IKPT,ISPIN) &
     &         ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
 1000 CONTINUE
      
      CALL FILEHANDLER$CLOSE('RESTART_IN')
!
!     ==================================================================
!     == TRANSFORM X(0),V(-1/2) INTO X(0),V(-1/2)                     ==
!     ==================================================================
      RM(:,:)=R0(:,:)-RM(:,:)*DELTAT
      CM(:,:,:,:)=C0(:,:,:,:)-CM(:,:,:,:)*DELTAT
      HM(:,:,:,:)=H0(:,:,:,:)-HM(:,:,:,:)*DELTAT
      RLAMM(:,:,:,:)=RLAM0(:,:,:,:)-RLAMM(:,:,:,:)*DELTAT
      ALLOCATE(EIG(NB,NKPT,NSPIN))
      DO IB=1,NB
        EIG(IB,:,:)=H0(IB,IB,:,:)
      ENDDO

 9000 CONTINUE
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$STOP('WAVEIOSTYLE1')
      END
!
!     ..................................................................
      SUBROUTINE WRITESTYLE2(DELTAT,NFI,NAT,NB,NKPT,NSPIN &
     &                      ,NGWX,NGWG,C0,CM,X0,XM,XE0,XEM,XE2M &
     &                      ,R0,RM,RLAM0,RLAMM,H0,HM)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NFI
      INTEGER(4),INTENT(IN)   :: NAT
      REAL(8)   ,INTENT(IN)   :: R0(3,NAT)
      REAL(8)   ,INTENT(INOUT):: RM(3,NAT)
      REAL(8)   ,INTENT(IN)   :: DELTAT
      INTEGER(4),INTENT(IN)   :: NB
      INTEGER(4),INTENT(IN)   :: NKPT
      INTEGER(4),INTENT(IN)   :: NSPIN
      INTEGER(4),INTENT(IN)   :: NGWX
      INTEGER(4),INTENT(IN)   :: NGWG(NKPT)
      COMPLEX(8),INTENT(IN)   :: C0(NGWX,NB,NKPT,NSPIN)
      COMPLEX(8),INTENT(INOUT):: CM(NGWX,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(IN)   :: RLAM0(NB,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(INOUT):: RLAMM(NB,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(IN)   :: H0(NB,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(INOUT):: HM(NB,NB,NKPT,NSPIN)
      INTEGER(4)              :: NFIL
      REAL(8)   ,INTENT(IN)   :: XE0,XEM,XE2M
      REAL(8)   ,INTENT(IN)   :: X0,XM
      INTEGER(4)              :: IB,IKPT,ISPIN,IB1,IB2,I,J,K
      INTEGER(4)              :: IOS
      CHARACTER(80)           :: TEXT
!         
!     ================================================================
!     ==  TRANSFORM TO VELOCITIES
!     ================================================================
      RM(:,:)=(R0(:,:)-RM(:,:))/DELTAT
      RLAMM(:,:,:,:)=(RLAM0(:,:,:,:)-RLAMM(:,:,:,:))/DELTAT
      HM(:,:,:,:)=(H0(:,:,:,:)-HM(:,:,:,:))/DELTAT
      CM(:,:,:,:)=(C0(:,:,:,:)-CM(:,:,:,:))/DELTAT
!     
!     ================================================================
!     ==  WRITE ON FILE                                             ==
!     ================================================================
      CALL FILEHANDLER$UNIT('RESTART_OUT',NFIL)
      REWIND NFIL
      WRITE(NFIL,ERR=9999,IOSTAT=IOS)NFI,1,NB,NKPT,NSPIN,TEXT
      WRITE(NFIL,ERR=9999,IOSTAT=IOS)NAT,NGWG(:)
      WRITE(NFIL,ERR=9999,IOSTAT=IOS)((R0(I,J),RM(I,J),I=1,3),J=1,NAT)
!     
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NB
            WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &         (C0(I,IB,IKPT,ISPIN),CM(I,IB,IKPT,ISPIN),I=1,NGWG(IKPT))
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,ERR=9999,IOSTAT=IOS)X0,XM,XE0,XEM,XE2M
      WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &           ((((RLAM0(IB1,IB2,IKPT,ISPIN) &
     &              ,RLAMM(IB1,IB2,IKPT,ISPIN) &
     &              ,0.D0 &
     &              ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
      WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &           ((((H0(IB1,IB2,IKPT,ISPIN) &
     &              ,HM(IB1,IB2,IKPT,ISPIN) &
     &              ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
      CALL FILEHANDLER$CLOSE('RESTART_OUT')
 9999 CONTINUE
      CALL IOSTATMESSAGE(IOS,NFIL)
      CALL ERROR$MSG('ERROR WHILE WRITING TO FILE')
      CALL ERROR$STOP('WAVEIOSTYLE1')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WRITESTYLE1(DELTAT,NFI,NAT,NB,NKPT,NSPIN &
     &                      ,NGWX,NGWG,C0,CM,X0,XM,XE0,XEM,XE2M &
     &                      ,R0,RM,RLAM0,RLAMM,H0,HM)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NFI
      INTEGER(4),INTENT(IN)   :: NAT
      REAL(8)   ,INTENT(IN)   :: R0(3,NAT)
      REAL(8)   ,INTENT(INOUT):: RM(3,NAT)
      REAL(8)   ,INTENT(IN)   :: DELTAT
      INTEGER(4),INTENT(IN)   :: NB
      INTEGER(4),INTENT(IN)   :: NKPT
      INTEGER(4),INTENT(IN)   :: NSPIN
      INTEGER(4),INTENT(IN)   :: NGWX
      INTEGER(4),INTENT(IN)   :: NGWG(NKPT)
      COMPLEX(8),INTENT(IN)   :: C0(NGWX,NB,NKPT,NSPIN)
      COMPLEX(8),INTENT(INOUT):: CM(NGWX,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(IN)   :: RLAM0(NB,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(INOUT):: RLAMM(NB,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(IN)   :: H0(NB,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(INOUT):: HM(NB,NB,NKPT,NSPIN)
      REAL(8)   ,INTENT(IN)   :: XE0,XEM,XE2M
      REAL(8)   ,INTENT(IN)   :: X0,XM
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: IB,IKPT,ISPIN,IB1,IB2,I,J,K
      INTEGER(4)              :: IOS
      CHARACTER(80)           :: TEXT
!     ****************************************************************
!         
!     ================================================================
!     ==  TRANSFORM TO VELOCITIES
!     ================================================================
      RM(:,:)=(R0(:,:)-RM(:,:))/DELTAT
      RLAMM(:,:,:,:)=(RLAM0(:,:,:,:)-RLAMM(:,:,:,:))/DELTAT
      HM(:,:,:,:)=(H0(:,:,:,:)-HM(:,:,:,:))/DELTAT
      CM(:,:,:,:)=(C0(:,:,:,:)-CM(:,:,:,:))/DELTAT

!         
!     ================================================================
!     ==  WRITE ON FILE                                             ==
!     ================================================================
      CALL FILEHANDLER$UNIT('RESTART_OUT',NFIL)
      REWIND NFIL
      WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &           NFI,1,NB,NKPT,NSPIN,TEXT
      WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &           NAT,(NGWG(K),K=1,NKPT)
      WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &            ((R0(I,J),RM(I,J),I=1,3),J=1,NAT) &
     &           ,((((C0(I,J,K,ISPIN),CM(I,J,K,ISPIN) &
     &                    ,I=1,NGWG(K)),J=1,NB),K=1,NKPT), &
     &     ISPIN=1,NSPIN),X0,XM,XE0,XEM,XE2M
      WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &           ((((RLAM0(IB1,IB2,IKPT,ISPIN) &
     &              ,RLAMM(IB1,IB2,IKPT,ISPIN) &
     &              ,0.D0 &
     &              ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
      WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &           ((((H0(IB1,IB2,IKPT,ISPIN) &
     &              ,HM(IB1,IB2,IKPT,ISPIN) &
     &              ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
      CALL FILEHANDLER$CLOSE('RESTART_OUT')
 9999 CONTINUE
      CALL IOSTATMESSAGE(IOS,NFIL)
      CALL ERROR$MSG('ERROR WHILE WRITING TO FILE')
      CALL ERROR$STOP('WAVEIOSTYLE1')
      RETURN
      END


