!
!.......................................................................
      MODULE BANDDATA_MODULE
        INTEGER(4)             :: NDIM
        LOGICAL(4)             :: NDIM_SET
        INTEGER(4)             :: NSPIN
        LOGICAL(4)             :: NSPIN_SET
        INTEGER(4)             :: NDIMD
        LOGICAL(4)             :: NDIMD_SET

        REAL(8)                :: EPW !IN HARTREE
        LOGICAL(4)             :: EPW_SET
        REAL(8)                :: RBAS(3,3)
        LOGICAL(4)             :: RBAS_SET
        REAL(8)                :: GBAS(3,3)
        LOGICAL(4)             :: GBAS_SET !PLANEWAVE$INITIALIZE

        INTEGER(4)             :: NRL
        LOGICAL(4)             :: NRL_SET
        INTEGER(4)             :: NR1
        LOGICAL(4)             :: NR1_SET
        INTEGER(4)             :: NR2
        LOGICAL(4)             :: NR2_SET
        INTEGER(4)             :: NR3
        LOGICAL(4)             :: NR3_SET
        REAL(8),allocatable    :: VOFR(:,:)

        INTEGER(4)             :: NAT
        LOGICAL(4)             :: NAT_SET
        INTEGER(4)             :: NSP
        LOGICAL(4)             :: NSP_SET
        INTEGER(4)             :: NPRO
        LOGICAL(4)             :: NPRO_SET
        INTEGER(4)             :: LNXX
        LOGICAL(4)             :: LNXX_SET
        INTEGER(4)             :: LMNXX
        LOGICAL(4)             :: LMNXX_SET
        INTEGER(4)             :: LMX
        LOGICAL(4)             :: LMX_SET
        INTEGER(4)             :: NBAREPRO
        LOGICAL(4)             :: NBAREPRO_SET
        
        !THE FOLLOWING IS INCLUDED IN CASE OF GRID CHANGES
        INTEGER(4)             :: NG_PROTO
        LOGICAL(4)             :: NG_PROTO_SET
        INTEGER(4)             :: TYPE_PROTO
        LOGICAL(4)             :: TYPE_PROTO_SET
        REAL(8)                :: GMAX_PROTO
        LOGICAL(4)             :: GMAX_PROTO_SET
        REAL(8)                :: G1_PROTO
        LOGICAL(4)             :: G1_PROTO_SET
        REAL(8)                :: DEX_PROTO
        LOGICAL(4)             :: DEX_PROTO_SET
        
        INTEGER(4),allocatable :: ISPECIES(:) !NAT
        INTEGER(4),allocatable :: LNX(:) !NSP
        INTEGER(4),allocatable :: LMNX(:) !NSP
        INTEGER(4),allocatable :: LOX(:,:) !NSP,LNXX
        REAL(8)   ,allocatable :: R(:,:) !NAT
        CHARACTER(16),ALLOCATABLE :: ATOMID(:) !NAT
        REAL(8)   ,allocatable :: PROOFG(:,:,:) !
        COMPLEX(8),allocatable :: DH(:,:,:,:)!LMNXX,LMNXX,NDIMD,NAT
        REAL(8)   ,allocatable :: DO(:,:,:,:)!LMNXX,LMNXX,NDIMD,NAT

        REAL(8)                :: NEL
        LOGICAL(4)             :: NEL_SET
        REAL(8)                :: SPIN
        LOGICAL(4)             :: SPIN_SET
      END MODULE BANDDATA_MODULE
!
!     ..................................................................
      SUBROUTINE BANDDATA$GETI4(ID,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'NG_PROTO') THEN
        IF(NG_PROTO_SET)THEN
          VAL=NG_PROTO
        ELSE
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$GETI4')        
        ENDIF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$GETI4')
      END IF
      RETURN
      END
!!
!!     ..................................................................
!      SUBROUTINE BANDDATA$GETR8(ID,VAL)
!      USE BANDDATA_MODULE
!      IMPLICIT NONE
!      CHARACTER(*),INTENT(IN) :: ID
!      REAL(8)     ,INTENT(OUT):: VAL
!!     ******************************************************************
!      IF(ID.EQ.'EPW') THEN
!        IF(EPW_SET)THEN
!          VAL=EPW
!        ELSE
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETI4')        
!        ENDIF
!      ELSE
!        CALL ERROR$MSG('ID NOT RECOGNIZED')
!        CALL ERROR$CHVAL('ID',ID)
!        CALL ERROR$STOP('BANDDATA$GETR8')
!      END IF
!      RETURN
!      END
!!
!!     ..................................................................
!      SUBROUTINE BANDDATA$GETI4A(ID,LEN,VAL)
!      USE BANDDATA_MODULE
!      IMPLICIT NONE
!      CHARACTER(*),INTENT(IN) :: ID
!      INTEGER(4)  ,INTENT(IN) :: LEN
!      INTEGER(4)  ,INTENT(OUT):: VAL(LEN)
!      INTEGER(4)              :: ISPIN,IKPT,I
!!     ******************************************************************
!      IF(ID.EQ.'LNX') THEN
!        IF(LEN.NE.NSP) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$I4VAL('LEN',LEN)
!          CALL ERROR$I4VAL('NSP',NSP)
!          CALL ERROR$STOP('BANDDATA$GETI4')
!        END IF
!        IF(.NOT.ALLOCATED(LNX)) THEN
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETI4')
!        END IF
!        VAL=LNX
!      ELSE IF(ID.EQ.'LOX') THEN
!        IF(LEN.NE.NSP*LNXX) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETI4')
!        END IF
!        IF(.NOT.ALLOCATED(LOX)) THEN 
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETI4')
!        END IF
!        VAL=RESHAPE(LOX,(/LNXX*NSP/))
!      ELSE IF(ID.EQ.'ISPECIES') THEN
!        IF(LEN.NE.NAT) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETI4')
!        END IF
!        IF(.NOT.ALLOCATED(ISPECIES)) THEN
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETI4')
!        END IF
!        VAL=ISPECIES
!      ELSE
!        CALL ERROR$MSG('ID NOT RECOGNIZED')
!        CALL ERROR$CHVAL('ID',ID)
!        CALL ERROR$STOP('BANDDATA$GETI4A')
!      END IF
!      RETURN
!      END
!!
!!     ..................................................................
!      SUBROUTINE BANDDATA$GETR8A(ID,LEN,VAL)
!      USE BANDDATA_MODULE
!      IMPLICIT NONE
!      CHARACTER(*),INTENT(IN) :: ID
!      INTEGER(4)  ,INTENT(IN) :: LEN
!      REAL(8)     ,INTENT(OUT):: VAL(LEN)
!!     ******************************************************************
!      IF(ID.EQ.'RBAS') THEN
!        IF(LEN.NE.9) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        IF(RBAS_SET)THEN
!          VAL=RESHAPE(RBAS,(/9/))
!        ELSE
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR84')        
!        ENDIF
!      ELSE IF(ID.EQ.'GBAS') THEN
!        IF(LEN.NE.9) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        IF(GBAS_SET)THEN
!          VAL=RESHAPE(GBAS,(/9/))
!        ELSE
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR84')        
!        ENDIF
!      ELSE IF(ID.EQ.'VOFR') THEN
!        IF(LEN.NE.NRL*NDIMD) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        IF(.NOT.ALLOCATED(VOFR)) THEN
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        VAL=RESHAPE(VOFR,(/NRL*NDIMD/))
!      ELSE IF(ID.EQ.'R') THEN
!        IF(LEN.NE.3*NAT) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        IF(.NOT.ALLOCATED(R)) THEN
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        VAL=RESHAPE(R,(/3*NAT/))
!      ELSE
!        CALL ERROR$MSG('ID NOT RECOGNIZED')
!        CALL ERROR$CHVAL('ID',ID)
!        CALL ERROR$STOP('BANDDATA$GETR8A')
!      END IF
!      RETURN
!      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$SETI4(ID,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NDIM') THEN
        NDIM=VAL
        NDIM_SET=.true.
      ELSE IF(ID.EQ.'NSPIN') THEN
        NSPIN=VAL
        NSPIN_SET=.true.
      ELSE IF(ID.EQ.'NDIMD') THEN
        NDIMD=VAL
        NDIMD_SET=.true.
      ELSE IF(ID.EQ.'NR1') THEN
        NR1=VAL
        NR1_SET=.true.
      ELSE IF(ID.EQ.'NR2') THEN
        NR2=VAL
        NR2_SET=.true.
      ELSE IF(ID.EQ.'NR3') THEN
        NR3=VAL
        NR3_SET=.true.
      ELSE IF(ID.EQ.'NRL') THEN
        NRL=VAL
        NRL_SET=.true.
      ELSE IF(ID.EQ.'NAT') THEN
        NAT=VAL
        NAT_SET=.true.
      ELSE IF(ID.EQ.'NSP') THEN
        NSP=VAL
        NSP_SET=.true.
      ELSE IF(ID.EQ.'NPRO') THEN
        NPRO=VAL
        NPRO_SET=.true.
      ELSE IF(ID.EQ.'LNXX') THEN
        LNXX=VAL
        LNXX_SET=.true.
      ELSE IF(ID.EQ.'LMNXX') THEN
        LMNXX=VAL
        LMNXX_SET=.true.
      ELSE IF(ID.EQ.'LMX') THEN
        LMX=VAL
        LMX_SET=.true.
      ELSE IF(ID.EQ.'NBAREPRO') THEN
        NBAREPRO=VAL
        NBAREPRO_SET=.true.
      ELSE IF(ID.EQ.'NG_PROTO') THEN
        NG_PROTO=VAL
        NG_PROTO_SET=.true.
      ELSE IF(ID.EQ.'TYPE_PROTO') THEN
        TYPE_PROTO=VAL
        TYPE_PROTO_SET=.true.
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$SETI4A(ID,LEN,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'ISPECIES') THEN
        IF(.not.NAT_SET)then
          CALL ERROR$MSG('NAT NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        endif
        IF(LEN.NE.NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        END IF
        IF(.NOT.ALLOCATED(ISPECIES))ALLOCATE(ISPECIES(NAT))
        ISPECIES=VAL
      ELSE IF(ID.EQ.'LNX') THEN
        IF(.not.NSP_SET)then
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        endif
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        END IF
        IF(.NOT.ALLOCATED(LNX))ALLOCATE(LNX(NSP))
        LNX=VAL
      ELSE IF(ID.EQ.'LMNX') THEN
        IF(.not.NSP_SET)then
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        endif
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        END IF
        IF(.NOT.ALLOCATED(LMNX))ALLOCATE(LMNX(NSP))
        LMNX=VAL
      ELSE IF(ID.EQ.'LOX') THEN
        IF(.not.NSP_SET)then
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4AA')
        endif
        IF(.not.LNXX_SET)then
          CALL ERROR$MSG('LNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        endif
        IF(LEN.NE.NSP*LNXX) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$I4VAL('LNXX',LNXX)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        END IF
        IF(.NOT.ALLOCATED(LOX))ALLOCATE(LOX(LNXX,NSP))
        LOX=RESHAPE(VAL,(/LNXX,NSP/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$SETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$SETR8(ID,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'EPW') THEN
        EPW=VAL
        EPW_SET=.true.
      ELSE IF(ID.EQ.'GMAX_PROTO') THEN
        GMAX_PROTO=VAL
        GMAX_PROTO_SET=.true.
      ELSE IF(ID.EQ.'G1_PROTO') THEN
        G1_PROTO=VAL
        G1_PROTO_SET=.true.
      ELSE IF(ID.EQ.'DEX_PROTO') THEN
        DEX_PROTO=VAL
        DEX_PROTO_SET=.true.
      ELSE IF(ID.EQ.'NEL') THEN
        NEL=VAL
        NEL_SET=.true.
      ELSE IF(ID.EQ.'SPIN') THEN
        SPIN=VAL
        SPIN_SET=.true.
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$SETR8A(ID,LEN,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'RBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        RBAS(:,1)=VAL(1:3)
        RBAS(:,2)=VAL(4:6)
        RBAS(:,3)=VAL(7:9)
        RBAS_SET=.true.
      ELSE IF(ID.EQ.'GBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        GBAS(:,1)=VAL(1:3)
        GBAS(:,2)=VAL(4:6)
        GBAS(:,3)=VAL(7:9)
        GBAS_SET=.true.
      ELSE IF(ID.EQ.'VOFR') THEN
        IF(.not.NRL_SET)THEN
          CALL ERROR$MSG('NRL NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.not.NDIMD_SET)THEN
          CALL ERROR$MSG('NDIMD NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(LEN.NE.NRL*NDIMD)THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        IF(.not.allocated(VOFR))allocate(VOFR(NRL,NDIMD))
        VOFR=RESHAPE(VAL,(/NRL,NDIMD/))
      ELSE IF(ID.EQ.'R') THEN
        IF(.not.NAT_SET)THEN
          CALL ERROR$MSG('NAT NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(LEN.NE.3*NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(R)) ALLOCATE(R(3,NAT))
        R=RESHAPE(VAL,(/3,NAT/))
      ELSE IF(ID.EQ.'PROOFG') THEN
        IF(.not.NSP_SET)THEN
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.not.LNXX_SET)THEN
          CALL ERROR$MSG('LNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.not.NG_PROTO_SET)THEN
          CALL ERROR$MSG('NG_PROTO NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(LEN.NE.NG_PROTO*LNXX*NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(PROOFG)) ALLOCATE(PROOFG(NG_PROTO,LNXX,NSP))
        PROOFG=RESHAPE(VAL,(/NG_PROTO,LNXX,NSP/))
      ELSE IF(ID.EQ.'DO') THEN
        IF(.not.LMNXX_SET)THEN
          CALL ERROR$MSG('LMNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.not.NDIMD_SET)THEN
          CALL ERROR$MSG('NDIMD NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.not.NAT_SET)THEN
          CALL ERROR$MSG('NAT NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(LEN.NE.LMNXX*LMNXX*NDIMD*NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(DO)) ALLOCATE(DO(LMNXX,LMNXX,NDIMD,NAT))
        DO=RESHAPE(VAL,(/LMNXX,LMNXX,NDIMD,NAT/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$SETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$SETC8A(ID,LEN,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      COMPLEX(8)  ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'DH') THEN
        IF(.not.LMNXX_SET)THEN
          CALL ERROR$MSG('LMNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETC8A')
        ENDIF
        IF(.not.NDIMD_SET)THEN
          CALL ERROR$MSG('NDIMD NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETC8A')
        ENDIF
        IF(.not.NAT_SET)THEN
          CALL ERROR$MSG('NAT NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETC8A')
        ENDIF
        IF(LEN.NE.LMNXX*LMNXX*NDIMD*NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETC8A')
        END IF
        IF(.NOT.ALLOCATED(DH)) ALLOCATE(DH(LMNXX,LMNXX,NDIMD,NAT))
        DH=RESHAPE(VAL,(/LMNXX,LMNXX,NDIMD,NAT/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$SETC8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$SETCHA(ID,LEN,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.eq.'ATOMID')THEN
        IF(.not.NAT_SET)THEN
          CALL ERROR$MSG('NAT NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETCHA')
        ENDIF
        IF(LEN.NE.NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETCHA')
        END IF
        IF(.not.ALLOCATED(ATOMID))allocate(ATOMID(NAT))
        ATOMID(:)=VAL(:)
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$SETCHA')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$WRITEFILE
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: NFIL
!     ..................................................................
                             CALL TRACE$PUSH('BANDDATA$WRITEFILE')
      CALL FILEHANDLER$UNIT('BANDDATA',NFIL)
      REWIND NFIL
      CALL BANDDATA_WRITE(NFIL)
      CALL LIB$FLUSHFILE(NFIL)
      CALL FILEHANDLER$CLOSE('BANDDATA')
                             CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDDATA$WRITEFILE
!
!     ..................................................................
      SUBROUTINE BANDDATA$READFILE
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: NFIL
!     ..................................................................
                             CALL TRACE$PUSH('BANDDATA$READFILE')
      CALL FILEHANDLER$UNIT('BANDDATAIN',NFIL)
      REWIND NFIL
      CALL BANDDATA_READ(NFIL)
      CALL FILEHANDLER$CLOSE('BANDDATAIN')
                             CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDDATA$READFILE
!     ..................................................................
    
!
!     ..................................................................
      SUBROUTINE BANDDATA_WRITE(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: I,IKPT,ISPIN,IB
      INTEGER(4)             :: LNX1,NB
      LOGICAL(4)             :: SET
      LOGICAL(4)             :: TPRINT=.FALSE.
!     ******************************************************************
                             CALL TRACE$PUSH('BANDDATA_WRITE')
      SET=.true.
      SET=SET.and.NDIM_SET.and.NSPIN_SET.and.NDIMD_SET
      IF(TPRINT)PRINT*,NDIM_SET,NSPIN_SET,NDIMD_SET
      SET=SET.and.EPW_SET.and.GBAS_SET.and.RBAS_SET
      IF(TPRINT)PRINT*,EPW_SET,GBAS_SET,RBAS_SET
      SET=SET.and.NR1_SET.and.NR2_SET.and.NR3_SET.and.NRL_SET.and.allocated(VOFR)
      IF(TPRINT)PRINT*,NR1_SET,NR2_SET,NR3_SET,NRL_SET,allocated(VOFR)
      SET=SET.and.NAT_SET.and.NSP_SET.and.NPRO_SET.and.LNXX_SET.and.LMNXX_set
      IF(TPRINT)PRINT*,NAT_SET,NSP_SET,NPRO_SET,LNXX_SET,LMNXX_set
      SET=SET.and.NBAREPRO_SET.and.LMX_SET
      IF(TPRINT)PRINT*,NBAREPRO_SET,LMX_SET
      SET=SET.and.NG_PROTO_SET.and.TYPE_PROTO_SET.and.GMAX_PROTO_SET
      SET=SET.and.G1_PROTO_SET.and.DEX_PROTO_SET
      IF(TPRINT)PRINT*,NG_PROTO_SET,TYPE_PROTO_SET,GMAX_PROTO_SET,G1_PROTO_SET,DEX_PROTO_SET
      SET=SET.and.ALLOCATED(ISPECIES).and.allocated(LNX)
      IF(TPRINT)PRINT*,ALLOCATED(ISPECIES),allocated(LNX)
      SET=SET.and.ALLOCATED(LOX).and.allocated(LMNX)
      IF(TPRINT)PRINT*,ALLOCATED(LOX),allocated(LMNX)
      SET=SET.and.ALLOCATED(R).and.allocated(ATOMID).and.allocated(PROOFG)
      IF(TPRINT)PRINT*,ALLOCATED(R),allocated(ATOMID),allocated(PROOFG)
      SET=SET.and.NEL_SET.and.SPIN_SET
      IF(TPRINT)PRINT*,NEL_SET,SPIN_SET
      IF(.not.set)THEN
        CALL ERROR$MSG('ERROR WRITING BANDDATA FILE')
        CALL ERROR$MSG('NOT ALL QUANTITIES SET')
        CALL ERROR$STOP('BANDDATA$WRITE')
      ENDIF
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      WRITE(NFIL)NDIM,NSPIN,NDIMD
!
!     ==================================================================
!     == LATTICE AND G-GRID                                           ==
!     ==================================================================
      WRITE(NFIL)EPW,GBAS(:,:),RBAS(:,:)
!
!     ==================================================================
!     == POTENTIAL                                                    ==
!     ==================================================================
      WRITE(NFIL)NR1,NR2,NR3,NRL
      WRITE(NFIL)VOFR(:,:)
!
!     ==================================================================
!     == SETUPS                                                       ==
!     ==================================================================
      WRITE(NFIL)NAT,NSP,NPRO,LNXX,LMNXX,LMX,NBAREPRO
      WRITE(NFIL)NG_PROTO,TYPE_PROTO,GMAX_PROTO,G1_PROTO,DEX_PROTO
      WRITE(NFIL)ISPECIES(:),LNX(:),LOX(:,:),LMNX(:),R(:,:),ATOMID(:)
      WRITE(NFIL)PROOFG(:,:,:)
      WRITE(NFIL)DH(:,:,:,:)
      WRITE(NFIL)DO(:,:,:,:)
      WRITE(NFIL)NEL,SPIN
                             CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA_READ(NFIL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: ISPIN,IB
      INTEGER(4)             :: LNX1,NB
      INTEGER(4)             :: IOS
      CHARACTER(82)          :: IOSTATMSG
      CHARACTER(32)          :: FLAG   ! DATE SPECIFYING A VERSION
      LOGICAL(4)             :: TPRINT=.FALSE.
!     ******************************************************************
                             CALL TRACE$PUSH('BANDDATA_READ')
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      READ(NFIL)NDIM,NSPIN,NDIMD
      IF(TPRINT)PRINT*,NDIM,NSPIN,NDIMD
      
      READ(NFIL)EPW,GBAS,RBAS
      IF(TPRINT)PRINT*,"EPW,GBAS,RBAS",EPW,GBAS,RBAS
      
      READ(NFIL)NR1,NR2,NR3,NRL
      IF(TPRINT)PRINT*,"NR123, NRL",NR1,NR2,NR3,NRL
      IF(.not.allocated(VOFR))allocate(VOFR(NRL,NDIMD))
      READ(NFIL)VOFR(:,:)
      IF(TPRINT)PRINT*,'POTENTIAL',NRL,VOFR(1,1)
      
      READ(NFIL)NAT,NSP,NPRO,LNXX,LMNXX,LMX,NBAREPRO
      IF(TPRINT)PRINT*,"NAT...",NAT,NSP,NPRO,LNXX,LMNXX,LMX,NBAREPRO
      
      READ(NFIL)NG_PROTO,TYPE_PROTO,GMAX_PROTO,G1_PROTO,DEX_PROTO
      IF(TPRINT)PRINT*,"NG_PROTO...",NG_PROTO,TYPE_PROTO,GMAX_PROTO,G1_PROTO,DEX_PROTO
      
      IF(.not.allocated(ISPECIES))allocate(ISPECIES(NAT))
      IF(.not.allocated(LNX))     allocate(LNX(NSP))
      IF(.not.allocated(LMNX))    allocate(LMNX(NSP))
      IF(.not.allocated(LOX))     allocate(LOX(NSP,LNXX))
      IF(.not.allocated(R))       allocate(R(3,NAT))
      IF(.not.allocated(ATOMID))  allocate(ATOMID(NAT))
      IF(.not.allocated(PROOFG))  allocate(PROOFG(NG_PROTO,LNXX,NSP))
      IF(.not.allocated(DH))      allocate(DH(LMNXX,LMNXX,NDIMD,NAT))
      IF(.not.allocated(DO))      allocate(DO(LMNXX,LMNXX,NDIMD,NAT))

      READ(NFIL)ISPECIES(:),LNX(:),LOX(:,:),LMNX(:),R(:,:),ATOMID(:)
      READ(NFIL)PROOFG(:,:,:)
      READ(NFIL)DH(:,:,:,:)
      READ(NFIL)DO(:,:,:,:)
      READ(NFIL)NEL,SPIN
      IF(TPRINT)PRINT*,"ISPECIES ",ISPECIES(:)
      IF(TPRINT)PRINT*,"LNX ",LNX(:)
      IF(TPRINT)PRINT*,"LOX ",LOX(:,:)
      IF(TPRINT)PRINT*,"LMNX ",LMNX(:)
      IF(TPRINT)PRINT*,"R ",R(1,1)
      IF(TPRINT)PRINT*,"ATOMID ",ATOMID(:)
      IF(TPRINT)PRINT*,"PROOFG ",PROOFG(1,1,:)
      IF(TPRINT)PRINT*,"DH ",DH(1,1,1,1)
      IF(TPRINT)PRINT*,"DO ",DO(1,1,1,1)
      IF(TPRINT)PRINT*,"NEL ",NEL
      IF(TPRINT)PRINT*,"SPIN ",SPIN
      CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDDATA_READ
