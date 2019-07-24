!
!.......................................................................
      MODULE BANDDATA_MODULE
        INTEGER(4)             :: NDIM
        LOGICAL(4)             :: NDIM_SET=.FALSE.
        INTEGER(4)             :: NSPIN
        LOGICAL(4)             :: NSPIN_SET=.FALSE.
        INTEGER(4)             :: NDIMD
        LOGICAL(4)             :: NDIMD_SET=.FALSE.

        REAL(8)                :: EPW !IN HARTREE
        LOGICAL(4)             :: EPW_SET=.FALSE.
        REAL(8)                :: RBAS(3,3)
        LOGICAL(4)             :: RBAS_SET=.FALSE.
        REAL(8)                :: GBAS(3,3)
        LOGICAL(4)             :: GBAS_SET=.FALSE. !PLANEWAVE$INITIALIZE

        INTEGER(4)             :: NRG
        LOGICAL(4)             :: NRG_SET=.FALSE.
        INTEGER(4)             :: NRL
        LOGICAL(4)             :: NRL_SET=.FALSE.
        INTEGER(4)             :: NR1
        LOGICAL(4)             :: NR1_SET=.FALSE.
        INTEGER(4)             :: NR2
        LOGICAL(4)             :: NR2_SET=.FALSE.
        INTEGER(4)             :: NR3
        LOGICAL(4)             :: NR3_SET=.FALSE.
        REAL(8),ALLOCATABLE    :: VOFRL(:,:) ! ONLY LOCAL POINTS 
        REAL(8),ALLOCATABLE    :: VOFR(:,:)  ! COMPLETE GLOBAL ARRAY
     
        INTEGER(4)             :: NAT
        LOGICAL(4)             :: NAT_SET=.FALSE.
        INTEGER(4)             :: NSP
        LOGICAL(4)             :: NSP_SET=.FALSE.
        INTEGER(4)             :: NPRO
        LOGICAL(4)             :: NPRO_SET=.FALSE.
        INTEGER(4)             :: LNXX
        LOGICAL(4)             :: LNXX_SET=.FALSE.
        INTEGER(4)             :: LMNXX
        LOGICAL(4)             :: LMNXX_SET=.FALSE.
        INTEGER(4)             :: LMX
        LOGICAL(4)             :: LMX_SET=.FALSE.
        INTEGER(4)             :: NBAREPRO
        LOGICAL(4)             :: NBAREPRO_SET=.FALSE.
        
        !THE FOLLOWING IS INCLUDED IN CASE OF GRID CHANGES
        INTEGER(4)             :: NG_PROTO
        LOGICAL(4)             :: NG_PROTO_SET=.FALSE.
        CHARACTER(6)           :: TYPEID_PROTO
        LOGICAL(4)             :: TYPEID_PROTO_SET=.FALSE.
        REAL(8)                :: GMAX_PROTO
        LOGICAL(4)             :: GMAX_PROTO_SET=.FALSE.
        REAL(8)                :: G1_PROTO
        LOGICAL(4)             :: G1_PROTO_SET=.FALSE.
        REAL(8)                :: DEX_PROTO
        LOGICAL(4)             :: DEX_PROTO_SET=.FALSE.
        
        INTEGER(4),ALLOCATABLE :: ISPECIES(:) !NAT
        INTEGER(4),ALLOCATABLE :: LNX(:) !NSP
        INTEGER(4),ALLOCATABLE :: LMNX(:) !NSP
        INTEGER(4),ALLOCATABLE :: LOX(:,:) !NSP,LNXX
        REAL(8)   ,ALLOCATABLE :: R(:,:) !NAT
        CHARACTER(16),ALLOCATABLE :: ATOMID(:) !NAT
        REAL(8)   ,ALLOCATABLE :: PROOFG(:,:,:) !
        COMPLEX(8),ALLOCATABLE :: DH(:,:,:,:)!LMNXX,LMNXX,NDIMD,NAT
        REAL(8)   ,ALLOCATABLE :: DO(:,:,:,:)!LMNXX,LMNXX,NDIMD,NAT

        REAL(8)                :: NEL
        LOGICAL(4)             :: NEL_SET=.FALSE.
        REAL(8)                :: SPIN
        LOGICAL(4)             :: SPIN_SET=.FALSE.
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
!        IF(LEN.NE.NRG*NDIMD) THEN
!          CALL ERROR$MSG('INCONSISTENT SIZE')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        IF(.NOT.ALLOCATED(VOFR)) THEN
!          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
!          CALL ERROR$CHVAL('ID',ID)
!          CALL ERROR$STOP('BANDDATA$GETR8A')
!        END IF
!        VAL=RESHAPE(VOFR,(/NRG*NDIMD/))
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
        NDIM_SET=.TRUE.
      ELSE IF(ID.EQ.'NSPIN') THEN
        NSPIN=VAL
        NSPIN_SET=.TRUE.
      ELSE IF(ID.EQ.'NDIMD') THEN
        NDIMD=VAL
        NDIMD_SET=.TRUE.
      ELSE IF(ID.EQ.'NR1') THEN
        NR1=VAL
        NR1_SET=.TRUE.
      ELSE IF(ID.EQ.'NR2') THEN
        NR2=VAL
        NR2_SET=.TRUE.
      ELSE IF(ID.EQ.'NR3') THEN
        NR3=VAL
        NR3_SET=.TRUE.
      ELSE IF(ID.EQ.'NRL') THEN
        NRL=VAL
        NRL_SET=.TRUE.
      ELSE IF(ID.EQ.'NRG') THEN
        NRG=VAL
        NRG_SET=.TRUE.
      ELSE IF(ID.EQ.'NAT') THEN
        NAT=VAL
        NAT_SET=.TRUE.
      ELSE IF(ID.EQ.'NSP') THEN
        NSP=VAL
        NSP_SET=.TRUE.
      ELSE IF(ID.EQ.'NPRO') THEN
        NPRO=VAL
        NPRO_SET=.TRUE.
      ELSE IF(ID.EQ.'LNXX') THEN
        LNXX=VAL
        LNXX_SET=.TRUE.
      ELSE IF(ID.EQ.'LMNXX') THEN
        LMNXX=VAL
        LMNXX_SET=.TRUE.
      ELSE IF(ID.EQ.'LMX') THEN
        LMX=VAL
        LMX_SET=.TRUE.
      ELSE IF(ID.EQ.'NBAREPRO') THEN
        NBAREPRO=VAL
        NBAREPRO_SET=.TRUE.
      ELSE IF(ID.EQ.'NG_PROTO') THEN
        NG_PROTO=VAL
        NG_PROTO_SET=.TRUE.
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
        IF(.NOT.NAT_SET)THEN
          CALL ERROR$MSG('NAT NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        ENDIF
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
        IF(.NOT.NSP_SET)THEN
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        ENDIF
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
        IF(.NOT.NSP_SET)THEN
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        ENDIF
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
        IF(.NOT.NSP_SET)THEN
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4AA')
        ENDIF
        IF(.NOT.LNXX_SET)THEN
          CALL ERROR$MSG('LNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETI4A')
        ENDIF
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
        EPW_SET=.TRUE.
      ELSE IF(ID.EQ.'GMAX_PROTO') THEN
        GMAX_PROTO=VAL
        GMAX_PROTO_SET=.TRUE.
      ELSE IF(ID.EQ.'G1_PROTO') THEN
        G1_PROTO=VAL
        G1_PROTO_SET=.TRUE.
      ELSE IF(ID.EQ.'DEX_PROTO') THEN
        DEX_PROTO=VAL
        DEX_PROTO_SET=.TRUE.
      ELSE IF(ID.EQ.'NEL') THEN
        NEL=VAL
        NEL_SET=.TRUE.
      ELSE IF(ID.EQ.'SPIN') THEN
        SPIN=VAL
        SPIN_SET=.TRUE.
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
        RBAS_SET=.TRUE.
      ELSE IF(ID.EQ.'GBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        GBAS(:,1)=VAL(1:3)
        GBAS(:,2)=VAL(4:6)
        GBAS(:,3)=VAL(7:9)
        GBAS_SET=.TRUE.
      ELSE IF(ID.EQ.'VOFRL') THEN
        IF(.NOT.NRL_SET)THEN
          CALL ERROR$MSG('NRL NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.NOT.NDIMD_SET)THEN
          CALL ERROR$MSG('NDIMD NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(LEN.NE.NRL*NDIMD)THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(VOFRL))ALLOCATE(VOFRL(NRL,NDIMD))
        VOFRL=RESHAPE(VAL,(/NRL,NDIMD/))
      ELSE IF(ID.EQ.'R') THEN
        IF(.NOT.NAT_SET)THEN
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
        IF(.NOT.NSP_SET)THEN
          CALL ERROR$MSG('NSP NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.NOT.LNXX_SET)THEN
          CALL ERROR$MSG('LNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.NOT.NG_PROTO_SET)THEN
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
        IF(.NOT.LMNXX_SET)THEN
          CALL ERROR$MSG('LMNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.NOT.NDIMD_SET)THEN
          CALL ERROR$MSG('NDIMD NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETR8A')
        ENDIF
        IF(.NOT.NAT_SET)THEN
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
        IF(.NOT.LMNXX_SET)THEN
          CALL ERROR$MSG('LMNXX NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETC8A')
        ENDIF
        IF(.NOT.NDIMD_SET)THEN
          CALL ERROR$MSG('NDIMD NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETC8A')
        ENDIF
        IF(.NOT.NAT_SET)THEN
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
      SUBROUTINE BANDDATA$SETCH(ID,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     *****************************************************************
      IF(ID.EQ.'TYPEID_PROTO')THEN
        TYPEID_PROTO=VAL
        TYPEID_PROTO_SET=.TRUE.
      ELSE 
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('BANDDATA$SETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BANDDATA$SETCHA(ID,LEN,VAL)
      USE BANDDATA_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
!     *****************************************************************
      IF(ID.EQ.'ATOMID')THEN
        IF(.NOT.NAT_SET)THEN
          CALL ERROR$MSG('NAT NOT SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETCHA')
        ENDIF
        IF(LEN.NE.NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('BANDDATA$SETCHA')
        END IF
        IF(.NOT.ALLOCATED(ATOMID))ALLOCATE(ATOMID(NAT))
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
      SUBROUTINE BANDDATA$COLLECTBANDDATA
!     ******************************************************************
!     **  THIS FUNCTION COLLECTS MOST DATA NEEDED FOR THE BANDDATA    **
!     **  OBJECT (SEE ALSO WAVES$ETOT, BECAUSE SOME STUFF IS ONLY     **
!     **  ACCESSIBLE THERE)
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)            :: IAT,ISP,LN
      INTEGER(4)            :: NDIM_
      INTEGER(4)            :: NSPIN_
      INTEGER(4)            :: NDIMD_
      REAL(8)               :: EPW_
      REAL(8)               :: RBAS_(3,3)
      REAL(8)               :: GBAS_(3,3)
      LOGICAL(4)            :: SET_
      REAL(8),ALLOCATABLE   :: R_(:,:)
      REAL(8),ALLOCATABLE   :: PROOFG_(:,:,:)
      CHARACTER(16),ALLOCATABLE :: ATOMID_(:)
      INTEGER(4)            :: NG_PROTO_
      INTEGER(4)            :: NR1_
      INTEGER(4)            :: NR2_
      INTEGER(4)            :: NR3_
      INTEGER(4)            :: NRG_
      REAL(8)               :: NEL_
      REAL(8)               :: SPIN_
      INTEGER(4)            :: NAT_
      INTEGER(4)            :: NSP_
      INTEGER(4)            :: LNXX_
      INTEGER(4)            :: LMX_
      INTEGER(4)            :: NBAREPRO_
      INTEGER(4)            :: NPRO_
      INTEGER(4),ALLOCATABLE:: ISPECIES_(:)
      INTEGER(4),ALLOCATABLE:: LNX_(:)
      INTEGER(4),ALLOCATABLE:: LMNX_(:)
      INTEGER(4),ALLOCATABLE:: LOX_(:,:)
!     ******************************************************************
      !GENERAL QUANTITIES
      CALL WAVES$GETI4('NDIM',NDIM_)
      CALL WAVES$GETI4('NSPIN',NSPIN_)
      CALL WAVES$GETI4('NDIMD',NDIMD_)
      CALL BANDDATA$SETI4('NDIM' ,NDIM_)
      CALL BANDDATA$SETI4('NSPIN',NSPIN_)
      CALL BANDDATA$SETI4('NDIMD',NDIMD_)
      
      CALL DYNOCC$GETR8('NEL',NEL_)
      CALL BANDDATA$SETR8('NEL',NEL_)
      CALL DYNOCC$GETR8('SPIN',SPIN_)
      CALL BANDDATA$SETR8('SPIN',SPIN_)

      CALL WAVES$GETR8('EPWPSI',EPW_)
      CALL BANDDATA$SETR8('EPW',EPW_)
      CALL WAVES$GETI4('NR1',NR1_)
      CALL WAVES$GETI4('NR2',NR2_)
      CALL WAVES$GETI4('NR3',NR3_)
      CALL BANDDATA$SETI4('NR1',NR1_)
      CALL BANDDATA$SETI4('NR2',NR2_)
      CALL BANDDATA$SETI4('NR3',NR3_)
      
      !LATTICE PARAMETERS
      CALL CELL$GETR8A('T0',9,RBAS_)
      CALL BANDDATA$SETR8A('RBAS',9,RBAS_)
      
      CALL PLANEWAVE$GETR8A('GBAS',9,GBAS_)
      CALL BANDDATA$SETR8A('GBAS',9,GBAS_)

      !POTENTIAL
      NRG_=NR1_*NR2_*NR3_
      CALL BANDDATA$SETI4('NRG',NRG_)
      CALL BANDDATA_COLLECTPOTENTIAL
      
      !SETUPS (MODIFIED CODE FROM WAVES$INITIALIZE)
      CALL SETUP$GETI4('NSP',NSP_)
      CALL BANDDATA$SETI4('NSP',NSP_)
      CALL SETUP$GETI4('LNXX',LNXX_) 
      CALL BANDDATA$SETI4('LNXX',LNXX_)
      CALL ATOMLIST$NATOM(NAT_)
      ALLOCATE(ISPECIES_(NAT_))
      ALLOCATE(LNX_(NSP_))
      ALLOCATE(LMNX_(NSP_))
      ALLOCATE(LOX_(LNXX_,NSP_))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT_,ISPECIES_)

      LMX_=0
      LOX_(:,:)=0
      NBAREPRO_=0
      DO ISP=1,NSP_
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX_(ISP))
        CALL SETUP$GETI4A('LOX',LNX_(ISP),LOX_(1:LNX_(ISP),ISP))
        LMNX_(ISP)=0
        DO LN=1,LNX_(ISP)
          LMNX_(ISP)=LMNX_(ISP)+2*LOX_(LN,ISP)+1
          LMX_=MAX(LMX_,(LOX_(LN,ISP)+1)**2)
        ENDDO
        NBAREPRO_=NBAREPRO_+LNX_(ISP)
        CALL SETUP$UNSELECT()
      ENDDO
      NPRO_=0
      DO IAT=1,NAT_
        ISP=ISPECIES_(IAT)
        NPRO_=NPRO_+LMNX_(ISP)
      ENDDO
      CALL BANDDATA$SETI4('NPRO',NPRO_)
      CALL BANDDATA$SETI4('LMNXX',MAXVAL(LMNX_))
      
      ALLOCATE(R_(3,NAT_))
      ALLOCATE(ATOMID_(NAT_))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT_,R_)
      CALL ATOMLIST$GETCHA('NAME',0,NAT_,ATOMID_)
      CALL BANDDATA$SETR8A('R',3*NAT_,R_)
      CALL BANDDATA$SETCHA('ATOMID',NAT_,ATOMID_)
      DEALLOCATE(R_)
      DEALLOCATE(ATOMID_)

      CALL BANDDATA$SETI4('LMX',LMX_)
      CALL BANDDATA$SETI4('NBAREPRO',NBAREPRO_)
      CALL BANDDATA$SETI4A('LNX',NSP_,LNX_)
      CALL BANDDATA$SETI4A('LMNX',NSP_,LMNX_)
      CALL BANDDATA$SETI4A('LOX',LNXX_*NSP_,LOX_)
      CALL BANDDATA$SETI4A('ISPECIES',NAT_,ISPECIES_)
      
      CALL BANDDATA$GETI4('NG_PROTO',NG_PROTO_)
      ALLOCATE(PROOFG_(NG_PROTO_,LNXX_,NSP_))
      DO ISP=1,NSP_
        CALL SETUP$ISELECT(ISP)
        PROOFG_(:,:,ISP)=0.0D0
        CALL SETUP$GETR8A('PROOFG',NG_PROTO_*LNX_(ISP),&
    &         PROOFG_(:,1:LNX_(ISP),ISP))
        CALL SETUP$UNSELECT()
      ENDDO
      CALL BANDDATA$SETR8A('PROOFG',NG_PROTO_*LNXX_*NSP_,PROOFG_)
      DEALLOCATE(PROOFG_)
      DEALLOCATE(ISPECIES_)
      DEALLOCATE(LNX_)
      DEALLOCATE(LMNX_)
      DEALLOCATE(LOX_)
      RETURN
      END SUBROUTINE BANDDATA$COLLECTBANDDATA
!
!     ..................................................................
      SUBROUTINE BANDDATA_COLLECTPOTENTIAL
!     ******************************************************************
!     **  THIS FUNCTION COLLETS THE REAL SPACE POTENTIAL FROM ALL     **
!     **  PROCESSES SO THAT IS IT AVAIABLE ON THE FIRST NODE IN VOFR  **
!     ******************************************************************
      USE BANDDATA_MODULE
      INTEGER(4)            :: IDIMD
      INTEGER(4)            :: NTASKS,THISTASK
!     *****************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(.NOT.NDIMD_SET)THEN
        CALL ERROR$MSG('NDIMD NOT SET')
        CALL ERROR$STOP('BANDDATA_COLLECTPOTENTIAL')
      ENDIF
      IF(.NOT.NRL_SET)THEN
        CALL ERROR$MSG('NRL NOT SET')
        CALL ERROR$STOP('BANDDATA_COLLECTPOTENTIAL')
      ENDIF
      IF(.NOT.NRG_SET)THEN
        CALL ERROR$MSG('NRG NOT SET')
        CALL ERROR$STOP('BANDDATA_COLLECTPOTENTIAL')
      ENDIF
      IF(.NOT.ALLOCATED(VOFRL))THEN
        CALL ERROR$MSG('VOFRL NOT ALLOCATED')
        CALL ERROR$STOP('BANDDATA_COLLECTPOTENTIAL')
      ENDIF
      
      ALLOCATE(VOFR(NRG,NDIMD))
      DO IDIMD=1,NDIMD
        CALL PLANEWAVE$RSPACECOLLECTR8(NRL,VOFRL(:,IDIMD),NRG,VOFR(:,IDIMD))
      ENDDO
      RETURN
      END SUBROUTINE BANDDATA_COLLECTPOTENTIAL
!
!     ..................................................................
      SUBROUTINE BANDDATA$WRITEFILE
!     ******************************************************************
!     **  THIS FUNCTION CREATES AN OUTPUT-FILEHANDLER FOR THE         **
!     **  BANDDATA-FILE AND AND CALLS BANDDATA_WRITE WHICH DOES THE   **
!     **  ACTUAL WRITING                                              **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: NTASKS,THISTASK
!     *****************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
                             CALL TRACE$PUSH('BANDDATA$WRITEFILE')
      CALL BANDDATA$COLLECTBANDDATA
      
      IF(.NOT.NRG_SET)THEN
        CALL ERROR$MSG('NRG NOT SET')
        CALL ERROR$STOP('BANDDATA$WRITEFILE')
      ENDIF
      
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('BANDDATA',NFIL)
        REWIND NFIL
        CALL BANDDATA_WRITE(NFIL)
        FLUSH(NFIL)   ! FLUSH IS STANDARD FORTRAN SINCER F2003)
        CALL FILEHANDLER$CLOSE('BANDDATA')
      ENDIF
      DEALLOCATE(VOFR) !ALLOCATED IN BANDDATA_COLLECTPOTENTIAL
                             CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDDATA$WRITEFILE
!
!     ..................................................................
      SUBROUTINE BANDDATA$READFILE
!     ******************************************************************
!     **  THIS FUNCTION CREATES AN INPUT-FILEHANDLER FOR THE          **
!     **  BANDDATA-FILE AND AND CALLS BANDDATA_READ WHICH DOES THE    **
!     **  ACTUAL READING                                              **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: NFIL
!     ******************************************************************
                             CALL TRACE$PUSH('BANDDATA$READFILE')
      CALL FILEHANDLER$UNIT('BANDDATAIN',NFIL)
      REWIND NFIL
      CALL BANDDATA_READ(NFIL)
      CALL FILEHANDLER$CLOSE('BANDDATAIN')
                             CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDDATA$READFILE
!
!     ..................................................................
      SUBROUTINE BANDDATA_WRITE(NFIL)
!     ******************************************************************
!     ** THIS FUNCTION WRITES THE BANDDATA-FILE TO NFIL               **
!     ** FOR THAT THE *_SET-VARIABLES IN THE BANDDATA-MODULE ARE      **
!     ** CHECKED AND WHEN SOME VARIABLE, WHICH IS SUPPOSED TO BE      **
!     ** WRITTEN TO THE FILE HAS NOT BEEN SET, IT WILL STOP WITH AN   ** 
!     ** ERROR                                                        **
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
      SET=.TRUE.
      SET=SET.AND.NDIM_SET.AND.NSPIN_SET.AND.NDIMD_SET
      SET=SET.AND.EPW_SET.AND.GBAS_SET.AND.RBAS_SET
      SET=SET.AND.NR1_SET.AND.NR2_SET.AND.NR3_SET.AND.NRG_SET
      SET=SET.AND.ALLOCATED(VOFR)
      SET=SET.AND.NAT_SET.AND.NSP_SET.AND.NPRO_SET
      SET=SET.AND.LNXX_SET.AND.LMNXX_SET
      SET=SET.AND.NBAREPRO_SET.AND.LMX_SET
      SET=SET.AND.NG_PROTO_SET.AND.TYPEID_PROTO_SET.AND.GMAX_PROTO_SET
      SET=SET.AND.G1_PROTO_SET.AND.DEX_PROTO_SET
      SET=SET.AND.ALLOCATED(ISPECIES).AND.ALLOCATED(LNX)
      SET=SET.AND.ALLOCATED(LOX).AND.ALLOCATED(LMNX)
      SET=SET.AND.ALLOCATED(R).AND.ALLOCATED(ATOMID)
      SET=SET.AND.ALLOCATED(PROOFG)
      SET=SET.AND.NEL_SET.AND.SPIN_SET
      IF(.NOT.SET)THEN
        CALL ERROR$MSG('ERROR WRITING BANDDATA FILE')
        CALL ERROR$MSG('NOT ALL QUANTITIES SET')
        CALL ERROR$STOP('BANDDATA_WRITE')
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
      WRITE(NFIL)NR1,NR2,NR3,NRG
      WRITE(NFIL)VOFR(:,:)
!
!     ==================================================================
!     == SETUPS                                                       ==
!     ==================================================================
      WRITE(NFIL)NAT,NSP,NPRO,LNXX,LMNXX,LMX,NBAREPRO
      WRITE(NFIL)NG_PROTO,TYPEID_PROTO,GMAX_PROTO,G1_PROTO,DEX_PROTO
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
!     ** THIS FUNCTION READS THE BANDDATA-FILE FROM NFIL              **
!     ******************************************************************
      USE BANDDATA_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: ISPIN,IB
      INTEGER(4)             :: LNX1,NB
      INTEGER(4)             :: IOS
!     ******************************************************************
                             CALL TRACE$PUSH('BANDDATA_READ')
      READ(NFIL)NDIM,NSPIN,NDIMD
      READ(NFIL)EPW,GBAS,RBAS
      READ(NFIL)NR1,NR2,NR3,NRG
      IF(.NOT.ALLOCATED(VOFR))ALLOCATE(VOFR(NRG,NDIMD))
      READ(NFIL)VOFR(:,:)
      READ(NFIL)NAT,NSP,NPRO,LNXX,LMNXX,LMX,NBAREPRO
      READ(NFIL)NG_PROTO,TYPEID_PROTO,GMAX_PROTO,G1_PROTO,DEX_PROTO
      IF(.NOT.ALLOCATED(ISPECIES))ALLOCATE(ISPECIES(NAT))
      IF(.NOT.ALLOCATED(LNX))     ALLOCATE(LNX(NSP))
      IF(.NOT.ALLOCATED(LMNX))    ALLOCATE(LMNX(NSP))
      IF(.NOT.ALLOCATED(LOX))     ALLOCATE(LOX(LNXX,NSP))
      IF(.NOT.ALLOCATED(R))       ALLOCATE(R(3,NAT))
      IF(.NOT.ALLOCATED(ATOMID))  ALLOCATE(ATOMID(NAT))
      IF(.NOT.ALLOCATED(PROOFG))  ALLOCATE(PROOFG(NG_PROTO,LNXX,NSP))
      IF(.NOT.ALLOCATED(DH))      ALLOCATE(DH(LMNXX,LMNXX,NDIMD,NAT))
      IF(.NOT.ALLOCATED(DO))      ALLOCATE(DO(LMNXX,LMNXX,NDIMD,NAT))
      READ(NFIL)ISPECIES(:),LNX(:),LOX(:,:),LMNX(:),R(:,:),ATOMID(:)
      READ(NFIL)PROOFG(:,:,:)
      READ(NFIL)DH(:,:,:,:)
      READ(NFIL)DO(:,:,:,:)
      READ(NFIL)NEL,SPIN
                            CALL TRACE$POP
      RETURN
      END SUBROUTINE BANDDATA_READ
