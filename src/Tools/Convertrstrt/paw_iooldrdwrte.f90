!
!     ..................................................................
      SUBROUTINE WAVEIO
      CHARACTER(32) STYLE
      CHARACTER*(*) STRING
      CHARACTER*80 TEXT,TEXT1
      LOGICAL(4)  :: TCHK
!     
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      ENTRY WAVEIO$WRITELOCAL(TLOCAL_)
        CALL WAVEIOSTYLE1$WRITELOCAL(TLOCAL_)
        CALL WAVEIOSTYLE2$WRITELOCAL(TLOCAL_)
        RETURN
!     
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      ENTRY WAVEIO$RDWRTE(STRING,TNWSTR,NFI,TEXT,DELTAT)
                          CALL TRACE$PUSH('RDWRTE')
!
!\BEGIN(TEST)
      IF(STRING.EQ.'PUT') THEN
        CALL WRITERESTART
        CALL TRACE$POP
        RETURN
      ELSE IF(STRING.EQ.'GET') THEN
        CALL READRESTART(TCHK)
        IF(TCHK) then
          CALL TRACE$POP
          RETURN
        end if
      END IF
!\END(TEST)

        STYLE='NONE'
        IF(STRING.EQ.'PUT') THEN
          STYLE='STYLE2'
!         STYLE='STYLE1'
        ELSE IF(STRING.EQ.'GET') THEN
          CALL MPE$QUERY(NTASKNUM,NTASKID)
          IF(NTASKID.EQ.1) THEN
            CALL FILEHANDLER$UNIT('RESTART_IN',NFIL)
            REWIND NFIL
            READ(NFIL,END=9000,IOSTAT=IOS)NFI,NSP1,NB1,NKPT1,NSPIN1,TEXT
            CALL MPE$BCAST(TEXT,80,1)
          ELSE 
            CALL MPE$BCAST(TEXT,80,1)
          END IF
          TEXT1='TESTING OF ROUTINE RDWRTE FOR IO'
          IF(TEXT.EQ.TEXT1) THEN
            STYLE='STYLE1'
          ELSE
!           TEXT='STYLE(19JAN1995)'
            STYLE='STYLE2'
          END IF
        ELSE
          CALL ERROR$MSG('STRING MUST BE EITHER "PUT" OR "GET"')
          CALL ERROR$STOP('WAVEIO')
        END IF
!        
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,FMT='("RESTART STYLE=",A)')STYLE
        IF(STYLE(1:6).EQ.'STYLE1') THEN
          CALL WAVEIOSTYLE1$RDWRTE(STRING,TNWSTR,NFI,TEXT &
     &                 ,DELTAT,X0,XM,XE0,XEM,XE2M)
        ELSE IF(STYLE(1:6).EQ.'STYLE2') THEN
          CALL WAVEIOSTYLE2$RDWRTE(STRING,TNWSTR,NFI,TEXT &
     &                 ,DELTAT,X0,XM,XE0,XEM,XE2M)
        END IF
                          CALL TRACE$POP
        RETURN
!
 9000 CONTINUE
      WRITE(NFILO,FMT='("ERROR WHILE READING FROM FILE IN RDWRTE")')
      CALL IOSTATMESSAGE(IOS,NFIL)
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$STOP('WAVEIO')
!
      END
!
!     ..................................................................
      SUBROUTINE WAVEIOSTYLE2
      IMPLICIT INTEGER ($)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(256):: STRINGDUMMY
      CHARACTER(3)  ::STRING
      CHARACTER(80) :: TEXT
      LOGICAL      :: TLOCAL=.FALSE.
      LOGICAL      :: TLOCAL_
      LOGICAL      :: TNWSTR
      INTEGER(4)    :: NGWL(NKPT)              ;POINTER ($NGWL,NGWL)
      INTEGER(4)    :: NGWG(NKPT)              ;POINTER ($NGWG,NGWG)
      INTEGER(4)    :: NGWG1(NKPT1)            ;POINTER ($NGWG1,NGWG1)
      COMPLEX(8)   :: C0(NGWLX,NB,NKPT,NSPIN) ;POINTER ($C0,C0)
      COMPLEX(8)   :: CM(NGWLX,NB,NKPT,NSPIN) ;POINTER ($CM,CM)
      COMPLEX(8)   :: C0_TMP(NGWGX1)          ;POINTER ($C0_TMP,C0_TMP)
      COMPLEX(8)   :: CM_TMP(NGWGX1)          ;POINTER ($CM_TMP,CM_TMP)
      INTEGER(4)    :: NA1(NSP1)               ;POINTER ($NA1,NA1)
      REAL(8)       :: R0(3,NAT)               ;POINTER ($R0,R0)
      REAL(8)       :: RM(3,NAT)               ;POINTER ($RM,RM)
      REAL(8)       :: RV(3,NAT)               ;POINTER ($RV,RV)
      REAL(8)       :: RLAM0(NB,NB,NKPT,NSPIN) ;POINTER ($RLAM0,RLAM0)
      REAL(8)       :: RLAMM(NB,NB,NKPT,NSPIN) ;POINTER ($RLAMM,RLAMM)
      REAL(8)       :: H0(NB,NB,NKPT,NSPIN)    ;POINTER ($H0,H0)
      REAL(8)       :: HM(NB,NB,NKPT,NSPIN)    ;POINTER ($HM,HM)
      REAL(8)       :: RLAM0_TMP(NB1,NB1,NKPT1,NSPIN1) ;POINTER ($RLAM0_TMP,RLAM0_TMP)
      REAL(8)       :: RLAMM_TMP(NB1,NB1,NKPT1,NSPIN1) ;POINTER ($RLAMM_TMP,RLAMM_TMP)
      REAL(8)       :: H0_TMP(NB1,NB1,NKPT1,NSPIN1)    ;POINTER ($H0_TMP,H0_TMP)
      REAL(8)       :: HM_TMP(NB1,NB1,NKPT1,NSPIN1)    ;POINTER ($HM_TMP,HM_TMP)
      REAL(8),ALLOCATABLE :: EIG(:,:,:)
      CALL ERROR$MSG('WAVEIO MUST BE ACCESSED VIA ENTRIES')
      CALL ERROR$STOP('WAVEIO')
      RETURN
!     
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      ENTRY WAVEIOSTYLE2$WRITELOCAL(TLOCAL_)
        TLOCAL=TLOCAL_
        RETURN
!     
!     ******************************************************************
!     **  READ ATOMIC POSITIONS AND WAVE FUNCTIONS ETC. FROM FILE     **
!     ******************************************************************
      ENTRY WAVEIOSTYLE2$RDWRTE(STRING,TNWSTR,NFI,TEXT &
     &                 ,DELTAT,X0,XM,XE0,XEM,XE2M)
                          CALL TRACE$PUSH('WAVEIOSTYLE2$RDWRTE')
!
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      TEXT='STYLE(19JAN1995)'
      IF(STRING.NE.'PUT'.AND.STRING.NE.'GET') THEN
        CALL ERROR$MSG('STRING IS NOT GET OR PUT')
        CALL ERROR$CHVAL('STRING',STRING)
        CALL ERROR$STOP('RDWRTE')
      END IF
!
      CALL DYNOCC$GETI4('NB',NB)
      CALL DYNOCC$GETI4('NKPT',NKPT) 
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
      CALL WAVES$SETWAVE(NGWLX_,NX_,NKPTX_,NSPINX_,$C0,$CM)
      NGWLX=NGWLX_
      IF(NKPTX_.NE.NKPT.OR.NSPINX_.NE.NSPIN.OR.NX_.NE.NB) THEN
        CALL ERROR$MSG('DIMENSIONS FROM WAVES$SETWAVE INCONSISTENT')
        CALL ERROR$MSG('WITH ACTUAL VALUES')
        CALL ERROR$STOP('WAVEIO$RDWRTE')
      END IF
!
      CALL ATOMLIST$NATOM(NAT)
      CALL STACK$ALLOCATE(8,$R0,3*NAT)
      CALL STACK$ALLOCATE(8,$RM,3*NAT)
      CALL STACK$ALLOCATE(8,$RV,3*NAT)
      CALL STACK$ALLOCATE(8,$RLAM0,NB*NB*NKPT*NSPIN)
      CALL STACK$ALLOCATE(8,$RLAMM,NB*NB*NKPT*NSPIN)
      CALL STACK$ALLOCATE(8,$H0,NB*NB*NKPT*NSPIN)
      CALL STACK$ALLOCATE(8,$HM,NB*NB*NKPT*NSPIN)
!     == PARALLEL ======================================================
      CALL MPE$QUERY(NTASKNUM,NTASKID)
      CALL STACK$ALLOCATE(4,$NGWL,NKPT)
      CALL STACK$ALLOCATE (4,$NGWG,NKPT)
      NGWGX=-1
      DO IKPT=1,NKPT
        CALL PLANEWAVE$SELECT('WAVE',IKPT)
        CALL PLANEWAVE$GSIZE(NGWL(IKPT),NGWG(IKPT),NGW_MX_)
        NGWGX=MAX(NGWGX,NGWG(IKPT))
      ENDDO        
!     PRINT*,'NGWGX ',NGWGX,NB,NKPT,NSPIN
!
!     == PARALLEL ======================================================
!
!     ==================================================================
!     ==================================================================
!     == WRITE ON FILE                                                ==
!     ==================================================================
!     ==================================================================
      IF(STRING.EQ.'PUT') THEN 
        CALL PLANEWAVE$SELECT('WAVE',1)
        CALL PLANEWAVE$GSIZE(NGW,NGX,NGW_MX_)
        NB1=NB
        NKPT1=NKPT
        NSPIN1=NSPIN
        IF(NTASKID.EQ.1) THEN
!         ==============================================================
!         == TRANSFORM X(0),X(-) INTO X(0),V(-1/2)                    ==
!         ==============================================================
          CALL ATOMLIST$GETALL('R(0)',8*3,NAT,R0)
          CALL ATOMLIST$GETALL('R(-)',8*3,NAT,RM)
          DO IAT=1,NAT
            DO I=1,3
              RV(I,IAT)=(R0(I,IAT)-RM(I,IAT))/DELTAT
            ENDDO
          ENDDO
        
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,RLAM0)
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,RLAMM)
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,H0)
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,HM)
!         CALL OCCUPATION$GET('RLAM(0)',8*NB*NB*NKPT*NSPIN,RLAM0)
!         CALL OCCUPATION$GET('RLAM(-)',8*NB*NB*NKPT*NSPIN,RLAMM)
!         CALL OCCUPATION$GET('H(0)',8*NB*NB*NKPT*NSPIN,H0)
!         CALL OCCUPATION$GET('H(-)',8*NB*NB*NKPT*NSPIN,HM)
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              DO IB1=1,NB
                DO IB2=1,NB
                  RLAMM(IB1,IB2,IKPT,ISPIN)=(RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                               -RLAMM(IB1,IB2,IKPT,ISPIN))/DELTAT
                  HM(IB1,IB2,IKPT,ISPIN)=(H0(IB1,IB2,IKPT,ISPIN) &
     &                                  -HM(IB1,IB2,IKPT,ISPIN))/DELTAT
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL THERMOSTAT$SELECT('ATOMS')
          CALL THERMOSTAT$GETR8('XNOS0',X0)
          CALL THERMOSTAT$GETR8('XNOSM',XM)
!         
!         ================================================================
!         ==  WRITE ON FILE                                             ==
!         ================================================================
          CALL FILEHANDLER$UNIT('RESTART_OUT',NFIL)
          WRITE(NFILO,*)'WAVEIOSTYLE2 NFIL:',NFIL
          CALL FLUSH_(NFILO)
          INQUIRE(UNIT=NFIL,NAME=STRINGDUMMY)
          WRITE(NFILO,*)'WAVEIOSTYLE2 NAME:',STRINGDUMMY
          CALL FLUSH_(NFILO)
          REWIND NFIL
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               NFI,1,NB,NKPT,NSPIN,TEXT
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               NAT,(NGWG(K),K=1,NKPT)
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &                ((R0(I,J),RV(I,J),I=1,3),J=1,NAT)
          CALL FLUSH_(NFILO)
!
          NGWGX1=NGWGX
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1)
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              CALL PLANEWAVE$SELECT('WAVE',IKPT)
              DO IB=1,NB
                CALL PLANEWAVE$COLLECT(16,NGWL(IKPT) &
     &                          ,C0(1,IB,IKPT,ISPIN),NGWG(IKPT),C0_TMP)
                CALL PLANEWAVE$COLLECT(16,NGWL(IKPT) &
     &                         ,CM(1,IB,IKPT,ISPIN),NGWG(IKPT),CM_TMP)
                DO IG=1,NGWG(IKPT)
                  CM_TMP(IG)=(C0_TMP(IG)-CM_TMP(IG))/DELTAT
                ENDDO
                WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               (C0_TMP(I),CM_TMP(I),I=1,NGWG(IKPT))
              ENDDO
            ENDDO
          ENDDO
          CALL STACK$FREE($C0_TMP)
          CALL STACK$FREE($CM_TMP)
          WRITE(NFIL,ERR=9999,IOSTAT=IOS)X0,XM,XE0,XEM,XE2M
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               ((((RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                  ,RLAMM(IB1,IB2,IKPT,ISPIN) &
     &                  ,0.D0 &
     &                  ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               ((((H0(IB1,IB2,IKPT,ISPIN) &
     &                  ,HM(IB1,IB2,IKPT,ISPIN) &
     &                  ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
!         PRINT *,' DISK WRITTEN',NFIL,' NFI = ',NFI
          CALL FILEHANDLER$CLOSE('RESTART_OUT')
        ELSE
!       == ON ALL OTHER TASKS JUST COLLECT ==
          NGWGX1=NGWGX
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1)
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              CALL PLANEWAVE$SELECT('WAVE',IKPT)
              DO IB=1,NB
                CALL PLANEWAVE$COLLECT(16,NGWL(IKPT) &
     &                        ,C0(1,IB,IKPT,ISPIN),NGWG(IKPT),C0_TMP)
                CALL PLANEWAVE$COLLECT(16,NGWL(IKPT) &
     &                        ,CM(1,IB,IKPT,ISPIN),NGWG(IKPT),CM_TMP)
              ENDDO
            ENDDO
          ENDDO
          CALL STACK$FREE($C0_TMP)
          CALL STACK$FREE($CM_TMP)
        END IF

      END IF
!         
!     ==================================================================
!     ==================================================================
!     ==  READ FROM FILE                                              ==
!     ==================================================================
!     ==================================================================
      IF(STRING.EQ.'GET') THEN
        IF(NTASKID.EQ.1) THEN
          CALL FILEHANDLER$UNIT('RESTART_IN',NFIL)
          REWIND NFIL
          READ(NFIL,END=9000,IOSTAT=IOS)NFI,NSP1,NB1,NKPT1,NSPIN1,TEXT
          CALL STACK$ALLOCATE(4,$NA1,NSP1)
          CALL STACK$ALLOCATE(4,$NGWG1,NKPT1)
          READ(NFIL,END=9000,IOSTAT=IOS) &
     &        (NA1(I),I=1,NSP1),(NGWG1(K),K=1,NKPT1)
          NAT1=0
          DO ISP=1,NSP1
            NAT1=NAT1+NA1(ISP)
          ENDDO
          CALL STACK$FREE($NA1)
          NGWGX1=NGWGX
          DO IKPT=1,NKPT1
            NGWGX1=MAX(NGWGX1,NGWG1(IKPT))
          ENDDO
          CALL STACK$ALLOCATE(16,$RLAM0_TMP,NB1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$RLAMM_TMP,NB1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$H0_TMP,NB1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$HM_TMP,NB1*NB1*NKPT1*NSPIN1)
          DO ISPIN=1,NSPIN1
            DO IKPT=1,NKPT1
              DO IB1=1,NB1
                DO IB2=1,NB1
                  RLAM0_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                  RLAMM_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                  H0_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                  HM_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!               
          IF((NAT1.NE.NAT).AND.(.NOT.TNWSTR)) THEN
            CALL ERROR$MSG('NUMBER OF ATOMS ON FILE IS TOO LARGE')
            CALL ERROR$OVERFLOW('NAT ON FILE',NAT1,NAT)
            CALL ERROR$STOP('WAVEIO$RDWRTE')
          END IF
!
!         ==   READ COORDINATES
          READ(NFIL,END=9000,IOSTAT=IOS) &
     &              ((R0(I,J),RV(I,J),I=1,3),J=1,NAT1)
          CALL MPE$BCAST(NSPIN1,4,1)
          CALL MPE$BCAST(NKPT1 ,4,1)
          CALL MPE$BCAST(NB1   ,4,1)
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1)
          DO ISPIN=1,MAX(NSPIN,NSPIN1)
            DO IKPT=1,MAX(NKPT,NKPT1)
              DO IB=1,MAX(NB,NB1)
                DO IG=1,NGWGX
                  C0_TMP(IG)=(0.D0,0.D0)
                  CM_TMP(IG)=(0.D0,0.D0)
                ENDDO
                IF(ISPIN.LE.NSPIN1.AND.IKPT.LE.NKPT1.AND.IB.LE.NB1) THEN
                  READ(NFIL,END=9000,IOSTAT=IOS) &
     &              (C0_TMP(I),CM_TMP(I),I=1,NGWG1(IKPT))
                ELSE
                  IS=ISPIN
                  IK=IKPT
                  I=IB
                  IF(ISPIN.GT.NSPIN1) IS=1
                  IF(IKPT.GT.NKPT1) IK=1
                  IF(IB.GT.NB1) I=1
                  CALL PLANEWAVE$SELECT('WAVE',IK)
                  CALL PLANEWAVE$COLLECT(16,NGWL(IK) &
     &                        ,C0(1,I,IK,IS),NGWG(IK),C0_TMP)
                  CALL PLANEWAVE$COLLECT(16,NGWL(IK) &
     &                       ,CM(1,I,IK,IS),NGWG(IK),CM_TMP)
                END IF
                IF(ISPIN.LE.NSPIN) THEN
                  IF(IKPT.LE.NKPT) THEN
                    IF(IB.LE.NB) THEN
                      IF(IKPT.LE.NKPT)CALL PLANEWAVE$SELECT('WAVE',IKPT)
                      CALL PLANEWAVE$DISTRIBUTE(16,NGWG(IKPT),C0_TMP &
     &                            ,NGWL(IKPT),C0(1,IB,IKPT,ISPIN))
                      CALL PLANEWAVE$DISTRIBUTE(16,NGWG(IKPT),CM_TMP &
     &                            ,NGWL(IKPT),CM(1,IB,IKPT,ISPIN))
                    END IF
                  END IF
                END IF
              ENDDO
            ENDDO
          ENDDO
          CALL STACK$FREE($C0_TMP)
          CALL STACK$FREE($CM_TMP)
          READ(NFIL,END=9000,IOSTAT=IOS)X0,XM,XE0,XEM,XE2M
          READ(NFIL,END=1000)((((RLAM0_TMP(IB1,IB2,IKPT,ISPIN) &
     &                          ,RLAMM_TMP(IB1,IB2,IKPT,ISPIN) &
     &                          ,DUMMY &
     &             ,IB1=1,NB1),IB2=1,NB1),IKPT=1,NKPT1),ISPIN=1,NSPIN1)
          READ(NFIL,END=1000)((((H0_TMP(IB1,IB2,IKPT,ISPIN) &
     &                          ,HM_TMP(IB1,IB2,IKPT,ISPIN) &
     &             ,IB1=1,NB1),IB2=1,NB1),IKPT=1,NKPT1),ISPIN=1,NSPIN1)
 1000     CONTINUE
!         WRITE(NFILO,*)' DISK READ  ',NFIL,' NFI = ',NFI
          
          CALL FILEHANDLER$CLOSE('RESTART_IN')
!
!         ==  COMPLETE MISSING K-POINTS
          DO ISPIN=1,NSPIN
            ISPIN1=ISPIN
            IF(NSPIN1.LT.NSPIN)ISPIN1=1
            DO IKPT=1,NKPT
              IKPT1=IKPT
              IF(NKPT1.LT.NKPT)IKPT1=1
              DO IB1=1,MIN(NB,NB1)
                DO IB2=1,MIN(NB,NB1)
                  RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                 =RLAM0_TMP(IB1,IB2,IKPT1,ISPIN1)
                  RLAMM(IB1,IB2,IKPT,ISPIN) &
     &                 =RLAMM_TMP(IB1,IB2,IKPT1,ISPIN1)
                  H0(IB1,IB2,IKPT,ISPIN)=H0_TMP(IB1,IB2,IKPT1,ISPIN1)
                  HM(IB1,IB2,IKPT,ISPIN)=HM_TMP(IB1,IB2,IKPT1,ISPIN1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL STACK$FREE($RLAMM_TMP)
          CALL STACK$FREE($RLAM0_TMP)
          CALL STACK$FREE($HM_TMP)
          CALL STACK$FREE($H0_TMP)
          CALL STACK$FREE($NGWG1)
        ELSE
          NGWGX1=1
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1)
          CALL MPE$BCAST(NSPIN1,4,1)
          CALL MPE$BCAST(NKPT1 ,4,1)
          CALL MPE$BCAST(NB1   ,4,1)
          DO ISPIN=1,MAX(NSPIN,NSPIN1)
            DO IKPT=1,MAX(NKPT,NKPT1)
              DO IB=1,MAX(NB,NB1)
                DO IG=1,NGWGX1
                  C0_TMP(IG)=(0.D0,0.D0)
                  CM_TMP(IG)=(0.D0,0.D0)
                ENDDO
                IF(ISPIN.LE.NSPIN1.AND.IKPT.LE.NKPT1.AND.IB.LE.NB1) THEN
                ELSE
                  IS=ISPIN
                  IK=IKPT
                  I=IB
                  IF(ISPIN.GT.NSPIN1) IS=1
                  IF(IKPT.GT.NKPT1) IK=1
                  IF(IB.GT.NB1) I=1
                  CALL PLANEWAVE$SELECT('WAVE',IK)
                  CALL PLANEWAVE$COLLECT(16,NGWL(IK) &
     &                        ,C0(1,I,IK,IS),NGWG(IK),C0_TMP)
                  CALL PLANEWAVE$COLLECT(16,NGWL(IK) &
     &                       ,CM(1,I,IK,IS),NGWG(IK),CM_TMP)
                END IF
                IF(ISPIN.LE.NSPIN) THEN
                  IF(IKPT.LE.NKPT) THEN
                    IF(IB.LE.NB) THEN
                      IF(IKPT.LE.NKPT)CALL PLANEWAVE$SELECT('WAVE',IKPT)
                      CALL PLANEWAVE$DISTRIBUTE(16,NGWG(IKPT),C0_TMP &
     &                            ,NGWL(IKPT),C0(1,IB,IKPT,ISPIN))
                      CALL PLANEWAVE$DISTRIBUTE(16,NGWG(IKPT),CM_TMP &
     &                            ,NGWL(IKPT),CM(1,IB,IKPT,ISPIN))
                    END IF
                  END IF
                END IF
              ENDDO
            ENDDO
          ENDDO
          CALL STACK$FREE($C0_TMP)
          CALL STACK$FREE($CM_TMP)
        END IF
!
!       == BROADCAST NOSE ==============================================
        CALL MPE$BCAST(X0  ,8,1)
        CALL MPE$BCAST(XM  ,8,1)
        CALL MPE$BCAST(XE0 ,8,1)
        CALL MPE$BCAST(XEM ,8,1)
        CALL MPE$BCAST(XE2M,8,1)
        CALL THERMOSTAT$SELECT('ATOMS')
        CALL THERMOSTAT$SETR8('XNOS0',X0)
        CALL THERMOSTAT$SETR8('XNOSM',XM)
!       ================================================================
        CALL MPE$BCAST(R0   ,8*3*NAT,1)
        CALL MPE$BCAST(RV   ,8*3*NAT,1)
        CALL MPE$BCAST(RLAM0,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(RLAMM,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(H0   ,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(HM   ,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(NFI   ,4,1)
!
!       == PARALLEL ======================================================
!
!       ==================================================================
!       == TRANSFORM X(0),V(-1/2) INTO X(0),V(-1/2)                     ==
!       ==================================================================

        DO IAT=1,NAT
          DO I=1,3
            RM(I,IAT)=R0(I,IAT)-RV(I,IAT)*DELTAT
          ENDDO
        ENDDO
        IF(.NOT.TNWSTR) THEN
          CALL ATOMLIST$SETALL('R(0)',8*3,NAT,R0)
          CALL ATOMLIST$SETALL('R(-)',8*3,NAT,RM)
        ELSE
          WRITE(NFILO,*)'WARNING: STRUCTURE NOT READ FROM RESTART FILE'
        END IF
!       
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NB
              DO IG=1,NGWL(IKPT)
                CM(IG,IB,IKPT,ISPIN)=C0(IG,IB,IKPT,ISPIN) &
     &                              -CM(IG,IB,IKPT,ISPIN)*DELTAT
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB1=1,NB
              DO IB2=1,NB
                RLAMM(IB1,IB2,IKPT,ISPIN)=RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                               -RLAMM(IB1,IB2,IKPT,ISPIN)*DELTAT
                HM(IB1,IB2,IKPT,ISPIN)=H0(IB1,IB2,IKPT,ISPIN) &
     &                                -HM(IB1,IB2,IKPT,ISPIN)*DELTAT
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL WAVES$SETR8A('RLAM(0)',NB*NB*NKPT*NSPIN,RLAM0)
        CALL WAVES$SETR8A('RLAM(-)',NB*NB*NKPT*NSPIN,RLAMM)
!       CALL OCCUPATION$SET('RLAM(0)',8*NB*NB*NKPT*NSPIN,RLAM0)
!       CALL OCCUPATION$SET('RLAM(-)',8*NB*NB*NKPT*NSPIN,RLAMM)
!       CALL OCCUPATION$SET('H(0)',8*NB*NB*NKPT*NSPIN,H0)
!       CALL OCCUPATION$SET('H(-)',8*NB*NB*NKPT*NSPIN,HM)
        ALLOCATE(EIG(NB,NKPT,NSPIN))
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NB
              EIG(IB,IKPT,ISPIN)=H0(IB,IB,IKPT,ISPIN)
            ENDDO
          ENDDO
        ENDDO
        CALL DYNOCC$SETR8A('EPSILON',NB*NKPT*NSPIN,EIG)
        DEALLOCATE(EIG)
      ENDIF
!
!     ==================================================================
!     ==  FREE ARRAYS                                                 ==
!     ==================================================================
      CALL STACK$FREE($RLAM0)
      CALL STACK$FREE($RLAMM)
      CALL STACK$FREE($H0)
      CALL STACK$FREE($HM)
      CALL STACK$FREE($R0)
      CALL STACK$FREE($RM)
      CALL STACK$FREE($RV)
      CALL STACK$FREE($NGWG)
      CALL STACK$FREE($NGWL)
                          CALL TRACE$POP
      RETURN
!       

 9000 CONTINUE
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$I4VAL('IOS ',IOS)
      CALL ERROR$STOP('RDWRTE')
 9999 CONTINUE
      WRITE(NFILO,FMT='("ERROR WHILE WRITING TO FILE IN RDWRTE")')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL IOSTATMESSAGE(IOS,NFIL)
      CALL ERROR$MSG('ERROR WHILE WRITING TO FILE')
      CALL ERROR$STOP('RDWRTE')
      END
!
!     ..................................................................
      SUBROUTINE WAVEIOSTYLE1
      IMPLICIT INTEGER ($)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 STRING
      CHARACTER*80 TEXT
!     == DUMMY ARRAYS ++++++++++++++++++++++++++++++++++++++++++++++++++
      LOGICAL TLOCAL_,TNWSTR
!     ==  POINTER ARRAYS
      POINTER   ($NGWL,NGWL),($NGWG,NGWG),($NGWG1,NGWG1)
      DIMENSION NGWL(NKPT),NGWG(NKPT),NGWG1(NKPT1)
      POINTER ($C0,C0),($CM,CM)
      COMPLEX(8) C0(NGWLX,NB,NKPT,NSPIN)
      COMPLEX(8) CM(NGWLX,NB,NKPT,NSPIN)
      POINTER ($C0_TMP,C0_TMP),($CM_TMP,CM_TMP)
      COMPLEX(8) C0_TMP(NGWGX1,NB1,NKPT1,NSPIN1)
      COMPLEX(8) CM_TMP(NGWGX1,NB1,NKPT1,NSPIN1)
!
      POINTER ($NA1,NA1)
      DIMENSION NA1(NSP1)
      POINTER ($R0,R0),($RM,RM),($RV,RV)
      DIMENSION R0(3,NAT),RM(3,NAT),RV(3,NAT)
!
      POINTER  ($RLAM0,RLAM0),($RLAMM,RLAMM)
      POINTER  ($H0,H0),($HM,HM)
      DIMENSION RLAM0(NB,NB,NKPT,NSPIN),RLAMM(NB,NB,NKPT,NSPIN)
      DIMENSION H0(NB,NB,NKPT,NSPIN),HM(NB,NB,NKPT,NSPIN)
!
      POINTER  ($RLAM0_TMP,RLAM0_TMP),($RLAMM_TMP,RLAMM_TMP)
      POINTER  ($H0_TMP,H0_TMP),($HM_TMP,HM_TMP)
      DIMENSION RLAM0_TMP(NB1,NB1,NKPT1,NSPIN1) &
     &         ,RLAMM_TMP(NB1,NB1,NKPT1,NSPIN1)
      DIMENSION H0_TMP(NB1,NB1,NKPT1,NSPIN1) &
     &         ,HM_TMP(NB1,NB1,NKPT1,NSPIN1)
      REAL(8),ALLOCATABLE :: EIG(:,:,:)

      LOGICAL TLOCAL
      DATA TLOCAL /.FALSE./
      CALL ERROR$MSG('WAVEIO MUST BE ACCESSED VIA ENTRIES')
      CALL ERROR$STOP('WAVEIO')
      RETURN
!     
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      ENTRY WAVEIOSTYLE1$WRITELOCAL(TLOCAL_)
        TLOCAL=TLOCAL_
        RETURN
!     
!     ******************************************************************
!     **  READ ATOMIC POSITIONS AND WAVE FUNCTIONS ETC. FROM FILE     **
!     ******************************************************************
      ENTRY WAVEIOSTYLE1$RDWRTE(STRING,TNWSTR,NFI,TEXT &
     &                 ,DELTAT,X0,XM,XE0,XEM,XE2M)
                          CALL TRACE$PUSH('WAVEIOSTYLE1$RDWRTE')
!
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      TEXT='TESTING OF ROUTINE RDWRTE FOR IO'
      IF(STRING.NE.'PUT'.AND.STRING.NE.'GET') THEN
        CALL ERROR$MSG('STRING IS NOT GET OR PUT')
        CALL ERROR$CHVAL('STRING',STRING)
        CALL ERROR$STOP('RDWRTE')
      END IF
!
      CALL DYNOCC$GETI4('NB',NB)
      CALL DYNOCC$GETI4('NKPT',NKPT) 
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
      CALL WAVES$SETWAVE(NGWLX_,NX_,NKPTX_,NSPINX_,$C0,$CM)
      NGWLX=NGWLX_
      IF(NKPTX_.NE.NKPT.OR.NSPINX_.NE.NSPIN.OR.NX_.NE.NB) THEN
        CALL ERROR$MSG('DIMENSIONS FROM WAVES$SETWAVE INCONSISTENT')
        CALL ERROR$MSG('WITH ACTUAL VALUES')
        CALL ERROR$STOP('WAVEIO$RDWRTE')
      END IF
!
      CALL ATOMLIST$NATOM(NAT)
      CALL STACK$ALLOCATE(8,$R0,3*NAT)
      CALL STACK$ALLOCATE(8,$RM,3*NAT)
      CALL STACK$ALLOCATE(8,$RV,3*NAT)
      CALL STACK$ALLOCATE(8,$RLAM0,NB*NB*NKPT*NSPIN)
      CALL STACK$ALLOCATE(8,$RLAMM,NB*NB*NKPT*NSPIN)
      CALL STACK$ALLOCATE(8,$H0,NB*NB*NKPT*NSPIN)
      CALL STACK$ALLOCATE(8,$HM,NB*NB*NKPT*NSPIN)
!     == PARALLEL ======================================================
      CALL MPE$QUERY(NTASKNUM,NTASKID)
      CALL STACK$ALLOCATE(4,$NGWL,NKPT)
      CALL STACK$ALLOCATE (4,$NGWG,NKPT)
      NGWGX=-1
      DO IKPT=1,NKPT
        CALL PLANEWAVE$SELECT('WAVE',IKPT)
        CALL PLANEWAVE$GSIZE(NGWL(IKPT),NGWG(IKPT),NGW_MX_)
        NGWGX=MAX(NGWGX,NGWG(IKPT))
      ENDDO        
!     PRINT*,'NGWGX ',NGWGX,NB,NKPT,NSPIN
!
!     == PARALLEL ======================================================
!
!     ==================================================================
!     ==================================================================
!     == WRITE ON FILE                                                ==
!     ==================================================================
!     ==================================================================
      IF(STRING.EQ.'PUT') THEN 
        CALL PLANEWAVE$SELECT('WAVE',1)
        CALL PLANEWAVE$GSIZE(NGW,NGX,NGW_MX_)
        NB1=NB
        NKPT1=NKPT
        NSPIN1=NSPIN
        IF (NTASKID.EQ.1) THEN
          NGWGX1=NGWGX
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
        ELSE
          NGWGX1=1
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
        ENDIF
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            CALL PLANEWAVE$SELECT('WAVE',IKPT)
            DO IB=1,NB
              CALL PLANEWAVE$COLLECT(16,NGWL(IKPT),C0(1,IB,IKPT,ISPIN) &
     &                          ,NGWG(IKPT),C0_TMP(1,IB,IKPT,ISPIN))
              CALL PLANEWAVE$COLLECT(16,NGWL(IKPT),CM(1,IB,IKPT,ISPIN) &
     &                         ,NGWG(IKPT),CM_TMP(1,IB,IKPT,ISPIN))
            ENDDO
          ENDDO
        ENDDO
      
        IF(NTASKID.EQ.1) THEN
!         ==============================================================
!         == TRANSFORM X(0),X(-) INTO X(0),V(-1/2)                    ==
!         ==============================================================
          CALL ATOMLIST$GETALL('R(0)',8*3,NAT,R0)
          CALL ATOMLIST$GETALL('R(-)',8*3,NAT,RM)
          DO IAT=1,NAT
            DO I=1,3
              RV(I,IAT)=(R0(I,IAT)-RM(I,IAT))/DELTAT
            ENDDO
          ENDDO
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              DO I=1,NB
                DO IG=1,NGWG(IKPT)
                  CM_TMP(IG,I,IKPT,ISPIN)=(C0_TMP(IG,I,IKPT,ISPIN) &
     &                                -CM_TMP(IG,I,IKPT,ISPIN))/DELTAT
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,RLAM0)
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,RLAMM)
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,H0)
          CALL WAVES$GETR8A('<PSI|H|PSI>',NB*NB*NKPT*NSPIN,HM)
!         CALL OCCUPATION$GET('RLAM(0)',8*NB*NB*NKPT*NSPIN,RLAM0)
!         CALL OCCUPATION$GET('RLAM(-)',8*NB*NB*NKPT*NSPIN,RLAMM)
!         CALL OCCUPATION$GET('H(0)',8*NB*NB*NKPT*NSPIN,H0)
!         CALL OCCUPATION$GET('H(-)',8*NB*NB*NKPT*NSPIN,HM)
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              DO IB1=1,NB
                DO IB2=1,NB
                  RLAMM(IB1,IB2,IKPT,ISPIN)=(RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                               -RLAMM(IB1,IB2,IKPT,ISPIN))/DELTAT
                  HM(IB1,IB2,IKPT,ISPIN)=(H0(IB1,IB2,IKPT,ISPIN) &
     &                                  -HM(IB1,IB2,IKPT,ISPIN))/DELTAT
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL THERMOSTAT$SELECT('ATOMS')
          CALL THERMOSTAT$GETR8('XNOS0',X0)
          CALL THERMOSTAT$GETR8('XNOSM',XM)
          CALL THERMOSTAT$GETR8('XNOSMM',XMM)
!         
!         ================================================================
!         ==  WRITE ON FILE                                             ==
!         ================================================================
          CALL FILEHANDLER$UNIT('RESTART_OUT',NFIL)
          REWIND NFIL
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               NFI,1,NB,NKPT,NSPIN,TEXT
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               NAT,(NGWG(K),K=1,NKPT)
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &                ((R0(I,J),RV(I,J),I=1,3),J=1,NAT) &
     &               ,((((C0_TMP(I,J,K,ISPIN),CM_TMP(I,J,K,ISPIN) &
     &                        ,I=1,NGWG(K)),J=1,NB),K=1,NKPT), &
     &         ISPIN=1,NSPIN),X0,XM,XE0,XEM,XE2M
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               ((((RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                  ,RLAMM(IB1,IB2,IKPT,ISPIN) &
     &                  ,0.D0 &
     &                  ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
          WRITE(NFIL,ERR=9999,IOSTAT=IOS) &
     &               ((((H0(IB1,IB2,IKPT,ISPIN) &
     &                  ,HM(IB1,IB2,IKPT,ISPIN) &
     &                  ,IB1=1,NB),IB2=1,NB),IKPT=1,NKPT),ISPIN=1,NSPIN)
!         PRINT *,' DISK WRITTEN',NFIL,' NFI = ',NFI
          
          CALL FILEHANDLER$CLOSE('RESTART_OUT')
        END IF
        CALL STACK$FREE($C0_TMP)
        CALL STACK$FREE($CM_TMP)
      END IF
!         
!     ==================================================================
!     ==  READ FROM FILE                                              ==
!     ==================================================================
      IF(STRING.EQ.'GET') THEN
        IF(NTASKID.EQ.1) THEN
          CALL FILEHANDLER$UNIT('RESTART_IN',NFIL)
          REWIND NFIL
          READ(NFIL,END=9000,IOSTAT=IOS)NFI,NSP1,NB1,NKPT1,NSPIN1,TEXT
          CALL STACK$ALLOCATE(4,$NA1,NSP1)
          CALL STACK$ALLOCATE(4,$NGWG1,NKPT1)
          READ(NFIL,END=9000,IOSTAT=IOS) &
     &        (NA1(I),I=1,NSP1),(NGWG1(K),K=1,NKPT1)
          NAT1=0
          DO ISP=1,NSP1
            NAT1=NAT1+NA1(ISP)
          ENDDO
          CALL STACK$FREE($NA1)
          NGWGX1=NGWGX
          DO IKPT=1,NKPT1
            NGWGX1=MAX(NGWGX1,NGWG1(IKPT))
          ENDDO
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$RLAM0_TMP,NB1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$RLAMM_TMP,NB1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$H0_TMP,NB1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$HM_TMP,NB1*NB1*NKPT1*NSPIN1)
          DO ISPIN=1,NSPIN1
            DO IKPT=1,NKPT1
              DO IB=1,NB1
                DO IG=1,NGWGX1
                  C0_TMP(IG,IB,IKPT,ISPIN)=(0.D0,0.D0)
                  CM_TMP(IG,IB,IKPT,ISPIN)=(0.D0,0.D0)
                ENDDO
              ENDDO
              DO IB1=1,NB1
                DO IB2=1,NB1
                  RLAM0_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                  RLAMM_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                  H0_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                  HM_TMP(IB1,IB2,IKPT,ISPIN)=0.D0
                ENDDO
              ENDDO
            ENDDO
          ENDDO
!               
          IF(NAT1.NE.NAT) THEN
            CALL ERROR$MSG('NUMBER OF ATOMS ON FILE IS TOO LARGE')
            CALL ERROR$OVERFLOW('NAT ON FILE',NAT1,NAT)
            CALL ERROR$STOP('WAVEIO$RDWRTE')
          END IF
!
!         ==   READ COORDINATES
          READ(NFIL,END=9000,IOSTAT=IOS) &
     &              ((R0(I,J),RV(I,J),I=1,3),J=1,NAT) &
     &             ,((((C0_TMP(I,J,K,ISPIN),CM_TMP(I,J,K,ISPIN) &
     &                ,I=1,NGWG1(K)),J=1,NB1),K=1,NKPT1),ISPIN=1,NSPIN1) &
     &             ,X0,XM,XE0,XEM,XE2M  
          READ(NFIL,END=1000)((((RLAM0_TMP(IB1,IB2,IKPT,ISPIN) &
     &                          ,RLAMM_TMP(IB1,IB2,IKPT,ISPIN) &
     &                          ,DUMMY &
     &             ,IB1=1,NB1),IB2=1,NB1),IKPT=1,NKPT1),ISPIN=1,NSPIN1)
          READ(NFIL,END=1000)((((H0_TMP(IB1,IB2,IKPT,ISPIN) &
     &                          ,HM_TMP(IB1,IB2,IKPT,ISPIN) &
     &             ,IB1=1,NB1),IB2=1,NB1),IKPT=1,NKPT1),ISPIN=1,NSPIN1)
 1000     CONTINUE
!         WRITE(NFILO,*)' DISK READ  ',NFIL,' NFI = ',NFI
          
          CALL FILEHANDLER$CLOSE('RESTART_IN')
!
!
!         ==  COMPLETE MISSING K-POINTS
          DO ISPIN=1,NSPIN
            ISPIN1=ISPIN
            IF(NSPIN1.LT.NSPIN)ISPIN1=1
            DO IKPT=1,NKPT
              IKPT1=IKPT
              IF(NKPT1.LT.NKPT)IKPT1=1
              DO IB1=1,MIN(NB,NB1)
                DO IB2=1,MIN(NB,NB1)
                  RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                 =RLAM0_TMP(IB1,IB2,IKPT1,ISPIN1)
                  RLAMM(IB1,IB2,IKPT,ISPIN) &
     &                 =RLAMM_TMP(IB1,IB2,IKPT1,ISPIN1)
                  H0(IB1,IB2,IKPT,ISPIN)=H0_TMP(IB1,IB2,IKPT1,ISPIN1)
                  HM(IB1,IB2,IKPT,ISPIN)=HM_TMP(IB1,IB2,IKPT1,ISPIN1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL STACK$FREE($RLAMM_TMP)
          CALL STACK$FREE($RLAM0_TMP)
          CALL STACK$FREE($HM_TMP)
          CALL STACK$FREE($H0_TMP)
          CALL STACK$FREE($NGWG1)
        END IF
!
!       == BROADCAST NOSE ==============================================
        CALL MPE$BCAST(X0  ,8,1)
        CALL MPE$BCAST(XM  ,8,1)
        CALL MPE$BCAST(XE0 ,8,1)
        CALL MPE$BCAST(XEM ,8,1)
        CALL MPE$BCAST(XE2M,8,1)
        CALL THERMOSTAT$SELECT('ATOMS')
        CALL THERMOSTAT$SETR8('XNOS0',X0)
        CALL THERMOSTAT$SETR8('XNOSM',XM)
!       ================================================================
        CALL MPE$BCAST(R0   ,8*3*NAT,1)
        CALL MPE$BCAST(RV   ,8*3*NAT,1)
        CALL MPE$BCAST(RLAM0,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(RLAMM,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(H0   ,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(HM   ,8*NB*NB*NKPT*NSPIN,1)
        CALL MPE$BCAST(NSPIN1,4,1)
        CALL MPE$BCAST(NKPT1 ,4,1)
        CALL MPE$BCAST(NB1   ,4,1)
        CALL MPE$BCAST(NFI   ,4,1)
        IF(NTASKID.NE.1) THEN
          NGWGX1=1
          CALL STACK$ALLOCATE(16,$C0_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
          CALL STACK$ALLOCATE(16,$CM_TMP,NGWGX1*NB1*NKPT1*NSPIN1)
        ENDIF
        DO ISPIN=1,NSPIN
          ISPIN1=ISPIN
          IF(ISPIN.GT.NSPIN1)ISPIN1=1
          DO IKPT=1,NKPT
            IKPT1=IKPT
            IF(IKPT.GT.NKPT1)IKPT1=1
            CALL PLANEWAVE$SELECT('WAVE',IKPT)
            DO IB=1,MIN(NB,NB1)
              CALL PLANEWAVE$DISTRIBUTE(16 &
     &              ,NGWG(IKPT),C0_TMP(1,IB,IKPT1,ISPIN1) &
     &              ,NGWL(IKPT),C0(1,IB,IKPT,ISPIN))
              CALL PLANEWAVE$DISTRIBUTE(16 &
     &             ,NGWG(IKPT),CM_TMP(1,IB,IKPT1,ISPIN1) &
     &             ,NGWL(IKPT),CM(1,IB,IKPT,ISPIN))
            ENDDO
          ENDDO
        ENDDO
        CALL STACK$FREE($C0_TMP)
        CALL STACK$FREE($CM_TMP)
!       == PARALLEL ======================================================
!
!       ==================================================================
!       == TRANSFORM X(0),V(-1/2) INTO X(0),V(-1/2)                     ==
!       ==================================================================
        DO IAT=1,NAT
          DO I=1,3
            RM(I,IAT)=R0(I,IAT)-RV(I,IAT)*DELTAT
          ENDDO
        ENDDO
        IF(.NOT.TNWSTR) THEN
          CALL ATOMLIST$SETALL('R(0)',8*3,NAT,R0)
          CALL ATOMLIST$SETALL('R(-)',8*3,NAT,RM)
        ELSE
          WRITE(NFILO,*)'WARNING: STRUCTURE NOT READ FROM RESTART FILE'
        END IF
!       
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NB
              DO IG=1,NGWL(IKPT)
                CM(IG,IB,IKPT,ISPIN)=C0(IG,IB,IKPT,ISPIN) &
     &                              -CM(IG,IB,IKPT,ISPIN)*DELTAT
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB1=1,NB
              DO IB2=1,NB
                RLAMM(IB1,IB2,IKPT,ISPIN)=RLAM0(IB1,IB2,IKPT,ISPIN) &
     &                               -RLAMM(IB1,IB2,IKPT,ISPIN)*DELTAT
                HM(IB1,IB2,IKPT,ISPIN)=H0(IB1,IB2,IKPT,ISPIN) &
     &                                -HM(IB1,IB2,IKPT,ISPIN)*DELTAT
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!       CALL OCCUPATION$SET('RLAM(0)',8*NB*NB*NKPT*NSPIN,RLAM0)
!       CALL OCCUPATION$SET('RLAM(-)',8*NB*NB*NKPT*NSPIN,RLAMM)
!       CALL OCCUPATION$SET('H(0)',8*NB*NB*NKPT*NSPIN,H0)
!       CALL OCCUPATION$SET('H(-)',8*NB*NB*NKPT*NSPIN,HM)
        CALL WAVES$SETR8A('RLAM(0)',NB*NB*NKPT*NSPIN,RLAM0)
        CALL WAVES$SETR8A('RLAM(-)',NB*NB*NKPT*NSPIN,RLAMM)
        ALLOCATE(EIG(NB,NKPT,NSPIN))
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NB
              EIG(IB,IKPT,ISPIN)=H0(IB,IB,IKPT,ISPIN)
            ENDDO
          ENDDO
        ENDDO
        CALL DYNOCC$SETR8A('EPSILON',NB*NKPT*NSPIN,EIG)
        DEALLOCATE(EIG)
      ENDIF
!
!     ==================================================================
!     ==  FREE ARRAYS                                                 ==
!     ==================================================================
      CALL STACK$FREE($RLAM0)
      CALL STACK$FREE($RLAMM)
      CALL STACK$FREE($H0)
      CALL STACK$FREE($HM)
      CALL STACK$FREE($R0)
      CALL STACK$FREE($RM)
      CALL STACK$FREE($RV)
      CALL STACK$FREE($NGWG)
      CALL STACK$FREE($NGWL)
                          CALL TRACE$POP
      RETURN
!       

 9000 CONTINUE
      WRITE(NFILO,FMT='("ERROR WHILE WRITING TO FILE IN RDWRTE")')
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$STOP('WAVEIOSTYLE1')
 9999 CONTINUE
      WRITE(NFILO,FMT='("ERROR WHILE WRITING TO FILE IN RDWRTE")')
      CALL IOSTATMESSAGE(IOS,NFIL)
      CALL ERROR$MSG('ERROR WHILE WRITING TO FILE')
      CALL ERROR$STOP('WAVEIOSTYLE1')
      END


