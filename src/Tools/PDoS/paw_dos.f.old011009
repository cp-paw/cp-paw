!************************************************************************
!***********************************************************************
!**                                                                   **
!**  NAME: PDOS                                                       **
!**                                                                   **
!**  PURPOSE: ANALYSIS TOOL FOR DENSITY OF STATES                     **
!**                                                                   **
!***********************************************************************
!***********************************************************************
      PROGRAM PDOS
      USE PDOS_MODULE, ONLY: STATE,STATEARR
      IMPLICIT NONE
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NAT
      INTEGER(4)                :: NB
      INTEGER(4)                :: NKPT
      INTEGER(4)                :: NSPIN
      INTEGER(4)                :: LENG
      INTEGER(4)                :: NSET
      INTEGER(4)   ,ALLOCATABLE :: LMX(:)
      REAL(8)      ,ALLOCATABLE :: RPOS(:,:)
      REAL(8)      ,ALLOCATABLE :: EIG(:,:,:)
      REAL(8)      ,ALLOCATABLE :: SET(:,:,:,:)
      CHARACTER(32),ALLOCATABLE :: LEGEND(:)
      REAL(8)                   :: EMIN
      REAL(8)                   :: EMAX
      INTEGER(4)                :: NE
      REAL(8)                   :: EBROAD
      INTEGER(4)                :: NFILIN
      INTEGER(4)                :: NPRO
      INTEGER(4)                :: IKPT,ISPIN,IB,IPRO,IDIM
      INTEGER(4)   ,ALLOCATABLE :: NBARR(:,:)
!     ******************************************************************
      CALL TRACE$PUSH('MAIN')
!
!     ==================================================================
!     ==  RESOLVE ARGUMENTLIST AND INITIALIZE FILE HANDLER            ==
!     ==================================================================
      CALL INITIALIZEFILEHANDLER
!
!     ==================================================================
!     ==  ANALYZE CONTROL FILE                                        ==
!     ==================================================================
      CALL READCNTL
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     ==  WRITE HEADER                                                ==
!     ==================================================================
      WRITE(NFILO,FMT='(72("*"))')
      WRITE(NFILO,FMT='(72("*"),T15 &
     &             ,"           PDOS ANALYSIS TOOL                ")')
      WRITE(NFILO,FMT='(72("*"),T15 &
     &             ,"    FOR THE PROJECTOR-AUGMENTED WAVE METHOD  ")')
      WRITE(NFILO,FMT='(72("*"))')
      WRITE(NFILO,FMT='(T28 &
     &           ," P.E. BLOECHL, IBM ZURICH RESEARCH LABORATORY ")')
      WRITE(NFILO,FMT='(T10 &
     &      ,"(C) IBM, 1995-1997 * ANY USE REQUIRES WRITTEN LICENSE FROM IBM")')
      WRITE(NFILO,*)
!
!     ==================================================================
!     ==  READ PDOSFILE                                               ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('PDOS',NFILIN)
      REWIND(NFILIN)
      CALL PDOS$READ(NFILIN)
      CALL PDOS$GETI4('NAT',NAT)
      CALL PDOS$GETI4('NKPT',NKPT)
      CALL PDOS$GETI4('NSPIN',NSPIN)
      CALL PDOS$GETI4('NPRO',NPRO)
      ALLOCATE(NBARR(NKPT,NSPIN))
      CALL PDOS$GETI4A('NB',nkpt*nspin,nbarr)
print*,'nspin',nspin
print*,'nbarr',nbarr
      NB=MAXVAL(NBARR)
      DEALLOCATE(NBARR)
      LENG=NPRO
      ALLOCATE(LMX(NAT))
      ALLOCATE(RPOS(3,NAT))
      CALL PDOS$GETR8A('R',3*NAT,RPOS)
      ALLOCATE(EIG(NB,NKPT,NSPIN))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          DO IB=1,NB
            EIG(IB,IKPT,ISPIN)=STATE%EIG(IB)
          ENDDO
        ENDDO
      ENDDO
                            CALL TRACE$PASS('AFTER READPDOS')
!
!     ==================================================================
!     ==  READ PREDEFINED ORBITALS                                    ==
!     ==================================================================
      CALL READCNTL$ORBITAL(LENG,NAT,LMX,RPOS)
                            CALL TRACE$PASS('AFTER READCNTL$ORBITAL')
!
!     ==================================================================
!     ==  SELECT MATRIXELEMENTS                                       ==
!     ==================================================================
      CALL READCNTL$SETNUMBER(NSET)
      ALLOCATE(SET(NB,NKPT,NSPIN,NSET))
      ALLOCATE(LEGEND(NSET))
      CALL READCNTL$SETS(NB,NKPT,NSPIN,NSET &
     &                  ,NAT,LMX,RPOS,LENG &
     &                  ,SET,LEGEND)
      DEALLOCATE(LMX)
                            CALL TRACE$PASS('AFTER READCNTL$SETS')
!
!     ==================================================================
!     ==  MAKE PLOTS                                                  ==
!     ==================================================================
      CALL READCNTL$GRID(EMIN,EMAX,NE,EBROAD)
                            CALL TRACE$PASS('AFTER READCNTL$GRID')
      CALL READCNTL$OUTPUT(EMIN,EMAX,NE,EBROAD &
     &                    ,NB,NKPT,NSPIN,EIG,NSET,SET,LEGEND)
                            CALL TRACE$PASS('AFTER READCNTL$OUTPUT')
!
!     ==================================================================
!     ==  CLOSING                                                     ==
!     ==================================================================
      CALL FILEHANDLER$REPORT(NFILO,'USED')
      WRITE(NFILO,FMT='(72("="))')
      WRITE(NFILO,FMT='(72("="),T20,"  PAW_DOS TOOL FINISHED  ")')
      WRITE(NFILO,FMT='(72("="))')
                            CALL TRACE$PASS('AFTER CLOSING')
!
!     ==================================================================
!     ==  CLOSE FILES                                                 ==
!     ==================================================================
      CALL FILEHANDLER$CLOSEALL
                            CALL TRACE$PASS('AFTER FILEHANDLER$CLOSEALL')
      CALL TRACE$POP
      CALL ERROR$NORMALSTOP
      STOP
      END
!      
!     ...................................................................
      SUBROUTINE INITIALIZEFILEHANDLER
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: PDOSINNAME
      INTEGER(4)     :: ISVAR
!     ******************************************************************
      IF(IARGC().LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE PDOS TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GETARG(1,PDOSINNAME)
      ISVAR=INDEX(PDOSINNAME,-'.DCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=PDOSINNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('DCNTL',.FALSE.,PDOSINNAME)
      RETURN
      END
!
!     ..................................................................     
      SUBROUTINE STANDARDFILES
!     *****************************************************************
!     **                                                             **
!     *****************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: T=.TRUE.
      LOGICAL(4),PARAMETER :: F=.FALSE.
      CHARACTER(32)        :: CH32SVAR1
      CHARACTER(32)        :: ID
      INTEGER(4)           :: NFILO
!     *****************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!  
!     ==================================================================
!     == SET STANDARD FILENAMES                                       ==
!     ==================================================================
!
!     ==  ERROR FILE ===================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,T,-'.DERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOLL FILE================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.DPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =============================================
      ID=+'DCNTL'
      CALL FILEHANDLER$SETFILE(ID,T,-'.DCNTL')
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
!
!     ==  STRUCTURE FILE   =============================================
      ID=+'PDOSOUT'
      CALL FILEHANDLER$SETFILE(ID,T,-'.PDOSOUT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      MODULE ORBITALS_MODULE
      PRIVATE
      PUBLIC ORBITALS$SETORB
      PUBLIC ORBITALS$GETORB
      LOGICAL(4)                :: TINI=.FALSE.
      INTEGER(4)                :: NORB=0      ! #(ORBITALS)
      INTEGER(4)                :: NORBX=0     ! MAX#(ORBITALS)
      INTEGER(4)    ,PARAMETER  :: NORBSTEP=10 ! STEP IN #(ORBITALS)
      INTEGER(4)                :: LENG=0      ! LENGTH OF ORBITALVECTOR
      complex(8)   ,ALLOCATABLE :: ORBITAL(:,:)
      CHARACTER(21),ALLOCATABLE :: ORBITALNAME(:)
      CONTAINS
!       ................................................................
        SUBROUTINE ORBITALS$SETORB(NAME_,LENG_,ORBITAL_)
        CHARACTER(*),INTENT(IN) :: NAME_
        INTEGER(4)  ,INTENT(IN) :: LENG_
        complex(8)  ,INTENT(IN) :: ORBITAL_(LENG_)
!       ****************************************************************
        IF(.NOT.TINI) THEN
          LENG=LENG_
          TINI=.TRUE.
        END IF
        IF(LENG.NE.LENG_) THEN
          CALL ERROR$MSG('LENGTH OF ORBITAL VECTOR INCONSISTENT')
          CALL ERROR$I4VAL('LENG_',LENG_)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$STOP('ORBITALS$SETORB')
        END IF
        CALL RESIZE
        NORB=NORB+1
        ORBITAL(:,NORB)=ORBITAL_(:)
        ORBITALNAME(NORB)=NAME_
        RETURN
        END SUBROUTINE ORBITALS$SETORB
!
!       ................................................................
        SUBROUTINE ORBITALS$GETORB(NAME_,LENG_,ORBITAL_)
        CHARACTER(*) ,INTENT(IN) :: NAME_
        INTEGER(4)   ,INTENT(IN) :: LENG_
        complex(8)   ,INTENT(OUT):: ORBITAL_(LENG_)
!       ****************************************************************
        IF(.NOT.TINI) THEN
          CALL ERROR$MSG('ORBITALS MODULE NOT YET INITIALIZED')
          CALL ERROR$STOP('ORBITALS$GETORB')
        END IF
        IF(LENG.NE.LENG_) THEN
          CALL ERROR$MSG('LENGTH OF ORBITAL VECTOR INCONSISTENT')
          CALL ERROR$I4VAL('LENG_',LENG_)
          CALL ERROR$I4VAL('LENG',LENG)
          CALL ERROR$STOP('ORBITALS$GETORB')
        END IF
        DO IORB=1,NORB
          IF(TRIM(ORBITALNAME(IORB)).EQ.TRIM(NAME_)) THEN
            ORBITAL_(:)=ORBITAL(:,IORB)
            RETURN
          END IF
        ENDDO
        CALL ERROR$MSG('ORBITAL NAME NOT RECOGNIZED')
        CALL ERROR$CHVAL('NAME_',NAME_)
        CALL ERROR$STOP('ORBITALS$GETORB')
        RETURN
        END SUBROUTINE ORBITALS$GETORB
!       ...............................................................
        SUBROUTINE RESIZE
        REAL(8)          ,ALLOCATABLE :: TMPORBITAL(:,:)
        CHARACTER(32),ALLOCATABLE :: TMPNAME(:)
        IF(NORB+1.LT.NORBX) RETURN
        IF(NORB.GT.0) THEN
          ALLOCATE(TMPORBITAL(LENG,NORB))
          ALLOCATE(TMPNAME(NORB))
          TMPORBITAL(:,:)=ORBITAL(:,1:NORB)
          TMPNAME(:)     =ORBITALNAME(1:NORB)
          DEALLOCATE(ORBITAL)
          DEALLOCATE(ORBITALNAME)
        END IF
        NORBX=NORB+NORBSTEP
        ALLOCATE(ORBITAL(LENG,NORBX))
        ALLOCATE(ORBITALNAME(NORBX))
        IF(NORB.GT.0) THEN
          ORBITAL(:,1:NORB)=TMPORBITAL(:,:)
          ORBITALNAME(1:NORB)=TMPNAME(:)
          DEALLOCATE(TMPORBITAL)
          DEALLOCATE(TMPNAME)
        END IF
        ORBITAL(:,NORB+1:NORBX)=0.D0
        ORBITALNAME(NORB+1:NORBX)=' '
        END SUBROUTINE RESIZE
      END MODULE ORBITALS_MODULE
!
!     ..................................................................
      SUBROUTINE READONEORB(LL_CNTL,NAT,LMX,RPOS,npro,ORBITAL)
!     ******************************************************************
!     **                                                              **
!     ** READ  AN ORBITAL BLOCK FROM LIST LL_CNTL AND RETURN          **
!     ** A VECTOR DEFINING THAT ORBITAL                               **
!     **                                                              **
!     ******************************************************************
      USE ORBITALS_MODULE
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NAT
      INTEGER(4)   ,INTENT(IN) :: LMX(NAT)
      REAL(8)      ,INTENT(IN) :: RPOS(3,NAT)
      INTEGER(4)   ,INTENT(IN) :: npro
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL
      COMPLEX(8)   ,INTENT(OUT):: ORBITAL(NPRO)
      INTEGER(4)   ,PARAMETER  :: LMXX=16
      REAL(8)                  :: ORB(LMXX)
      INTEGER(4)               :: IAT,I,iat2
      INTEGER(4)               :: I1,I2
      CHARACTER(32)            :: ATOM1
      CHARACTER(32)            :: ORBITALNAME1
      CHARACTER(8)             :: TYPE
      REAL(8)                  :: SVAR
      REAL(8)                  :: FAC
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: DRz(3)
      REAL(8)                  :: DRx(3)
      REAL(8)                  :: ROT(3,3)
      REAL(8)      ,ALLOCATABLE:: YLMROT(:,:)
!     ******************************************************************
      ORBITAL(:)=0.D0
!
!     ==================================================================
!     ==  GET PREFACTOR                                               ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FAC',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'FAC',1,FAC)
      ELSE
        FAC=1.D0
      END IF
!
!     ==================================================================
!     ==  SEARCH PREDEFINED ORBITALS                                  ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,ORBITALNAME1)
        CALL ORBITALS$GETORB(ORBITALNAME1,npro,ORBITAL)
      END IF
!
!     ==================================================================
!     ==  BUILD NEW ORBITAL COMPONENT FROM BASIC BUILDING BLOCKS      ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK)
      IF(TCHK) THEN   
        CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,ATOM1)
        call RESOLVEATOM(ATOM1,IAT)
        CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
        CALL RESOLVETYPE(LMXX,TYPE,FAC,ORB)
!       
!       ============================================================
!       ==  FIND NEARES NEIGHBOR DIRECTIONS                       ==
!       ============================================================
        drz=0.d0
        drx=0.d0
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NNZ',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NNZ',1,ATOM1)
          CALL RESOLVEATOM(ATOM1,IAT2)
          DRZ(:)=RPOS(:,IAT2)-RPOS(:,IAT)
          CALL LINKEDLIST$SET(LL_CNTL,'Z',0,DRZ(:))
        END IF
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NNX',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NNX',1,ATOM1)
          CALL RESOLVEATOM(ATOM1,IAT2)
          DRX(:)=RPOS(:,IAT2)-RPOS(:,IAT)
          CALL LINKEDLIST$SET(LL_CNTL,'X',0,DRX(:))
        END IF
!       
!       ============================================================
!       ==  make orbital                                          ==
!       ============================================================
        CALL RESOLVEROTATION(DRZ,DRX,ROT)
        ALLOCATE(YLMROT(LMXx,LMXx))
        CALL ROTATEYLM(LMXX,ROT,YLMROT)
        ORB=MATMUL(YLMROT,ORB)
        DEALLOCATE(YLMROT)
        call MAKEORBITAL(ATOM1,LMXX,ORB,NPRO,ORBITAL)
      END IF
      RETURN
      CONTAINS
!       ................................................................
        SUBROUTINE RESOLVETYPE(LMX,TYPE,FAC,ORBITAL)
        IMPLICIT NONE
        INTEGER(4)  ,INTENT(IN)   :: LMX
        CHARACTER(*),INTENT(IN)   :: TYPE
        REAL(8)     ,INTENT(IN)   :: FAC
        REAL(8)     ,INTENT(INOUT):: ORBITAL(LMX)
        REAL(8)                   :: ORB(9)
        INTEGER(4)                :: LM
!       ****************************************************************
!
!       ================================================================
!       ==  S-ONLY                                                    ==
!       ================================================================
        ORB(:)=0.D0
        IF(TRIM(TYPE).EQ.'S') THEN
          ORB(1)=FAC
        ELSE IF(TRIM(TYPE).EQ.'PX') THEN
          ORB(2)=FAC
        ELSE IF(TRIM(TYPE).EQ.'PZ') THEN
          ORB(3)=FAC
        ELSE IF(TRIM(TYPE).EQ.'PY') THEN
          ORB(4)=FAC
        ELSE IF(TRIM(TYPE).EQ.'SP1') THEN
          ORB(1)=FAC*DSQRT(1.D0/2.D0)
          ORB(3)=FAC*DSQRT(1.D0/2.D0)
        ELSE IF(TRIM(TYPE).EQ.'SP2') THEN
          ORB(1)=FAC*DSQRT(1.D0/3.D0)
          ORB(3)=FAC*DSQRT(2.D0/3.D0)
        ELSE IF(TRIM(TYPE).EQ.'SP3') THEN
          ORB(1)=FAC*DSQRT(1.D0/4.D0)
          ORB(3)=FAC*DSQRT(3.D0/4.D0)
        ELSE IF(TRIM(TYPE).EQ.'DX2-Y2') THEN
          ORB(5)=FAC
        ELSE IF(TRIM(TYPE).EQ.'DXZ') THEN
          ORB(6)=FAC
        ELSE IF(TRIM(TYPE).EQ.'D3Z2-R2') THEN
          ORB(7)=FAC
        ELSE IF(TRIM(TYPE).EQ.'DYZ') THEN
          ORB(8)=FAC
        ELSE IF(TRIM(TYPE).EQ.'DXY') THEN
          ORB(9)=FAC
        ELSE
          CALL ERROR$MSG('TYPE NOT IDENTIFIED')
          CALL ERROR$CHVAL('TYPE',TRIM(TYPE))
          CALL ERROR$STOP('RESOLVETYPE')
        END IF
        DO lm=LMX+1,9
          IF(ORB(lm).NE.0) THEN
            CALL ERROR$MSG('LMX TOO SMALL')
            CALL ERROR$STOP('RESOLVETYPE')
          END IF
        ENDDO
        ORBITAL(1:9)=ORB
        ORBITAL(10:)=0.D0
        RETURN
        END SUBROUTINE RESOLVETYPE
      END
!
!     ..................................................................
      SUBROUTINE MAKEORBITAL(ATOM,LMXX,ORB,NPRO_,ORBITAL)
!     **                                                              **
!     **  ONLY THE FIRST PARTIAL WAVE PER ANGULAR MOMENTUM IS         **
!     **  CONSIDERED                                                  **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ATOM
      INTEGER(4)  ,INTENT(IN) :: LMXX
      REAL(8)     ,INTENT(IN) :: ORB(LMXX)
      INTEGER(4)  ,INTENT(IN) :: NPRO_
      COMPLEX(8)  ,INTENT(OUT):: ORBITAL(npro_)
      INTEGER(4)              :: IPRO,IAT,ISP,LN,L,LM,M
      LOGICAL(4)              :: TCHK
      logical(4)              :: lchk(10)
!     ******************************************************************
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(Iat)
        TCHK=(ATOM.EQ.ATOMID(IAT))
        lchk=.true.
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(TCHK.AND.LCHK(L+1)) THEN
            LCHK(L+1)=.false.
            LM=L**2
            DO M=1,2*L+1
              IPRO=IPRO+1
              LM=LM+1
              ORBITAL(IPRO)=ORB(LM)
            ENDDO
          ELSE
            IPRO=IPRO+2*L+1
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE RESOLVEROTATION(DZ_,DX_,ROT)
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: DZ_(3)
      REAL(8)   ,INTENT(IN)  :: DX_(3)
      REAL(8)   ,INTENT(OUT) :: ROT(3,3)
      REAL(8)                :: DZ(3)
      REAL(8)                :: DX(3)
      REAL(8)                :: DY(3)
      REAL(8)                :: DZLEN,DXLEN
!     ******************************************************************
      dx=dx_
      dz=dz_
!     
!     ==================================================================
!     == NORMALIZE AND COMPLETE VECTORS                               ==
!     ==================================================================
!     == SET DZ ========================================================
      DZLEN=SQRT(DZ(1)**2+DZ(2)**2+DZ(3)**2)
      IF(DZLEN.EQ.0.D0) THEN
        DZ=(/0.D0,0.D0,1.D0/)
      ELSE
        DZ=DZ/DZLEN
      END IF
!     == SET DX ========================================================
      DXLEN=SQRT(DX(1)**2+DX(2)**2+DX(3)**2)
      IF(DxLEN.EQ.0.D0) THEN
        DX=(/1.D0,0.D0,0.D0/)
      ELSE
        DX=DX/DXLEN
      END IF
      DX=DX-DZ*DOT_PRODUCT(DZ,DX)
      DXLEN=SQRT(DX(1)**2+DX(2)**2+DX(3)**2)
      IF(DXLEN.EQ.0.D0) THEN
        DX=(/0.D0,1.D0,0.D0/)
      END IF
      DX=DX/DXLEN
!     == SET DY ========================================================
      DY(1)=Dz(2)*Dx(3)-Dz(3)*Dx(2)
      DY(2)=Dz(3)*Dx(1)-Dz(1)*Dx(3)
      DY(3)=Dz(1)*Dx(2)-Dz(2)*Dx(1)
!     
!     ==================================================================
!     == NORMALIZE AND COMPLETE VECTORS                               ==
!     ==================================================================
      ROT(:,3)=DZ(:)
      ROT(:,2)=DY(:)
      ROT(:,1)=DX(:)
!     WRITE(*,FMT='("ROT",3F10.5)')ROT(1,:)
!     WRITE(*,FMT='("ROT",3F10.5)')ROT(2,:)
!     WRITE(*,FMT='("ROT",3F10.5/)')ROT(3,:)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE normalizeorbital(npro_,orbital)
!     USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(in)    :: npro_
      complex(8)  ,intent(inout) :: orbital(npro_)
      INTEGER(4)                 :: I
      real(8)                    :: sum
!     ******************************************************************
      SUM=0.D0
      DO I=1,NPRO_
        SUM=SUM+REAL(CONJG(ORBITAL(I))*ORBITAL(I))
      ENDDO
      SUM=1.D0/SQRT(SUM)
      DO I=1,NPRO_
        ORBITAL(I)=ORBITAL(I)*SUM
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE RESOLVEATOM(ATOM,IAT)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ATOM
      INTEGER(4)  ,INTENT(OUT) :: IAT
      INTEGER(4)               :: I
!     ******************************************************************
      DO I=1,NAT
        IF(ATOM.EQ.ATOMID(I)) THEN
          IAT=I
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ATOM NAME NOT FOUND')
      CALL ERROR$CHVAL('ATOM',ATOM)
      CALL ERROR$STOP('RESOLVEATOM')
      RETURN
      END
!
!     ..................................................................
MODULE READCNTL_MODULE
USE LINKEDLIST_MODULE
TYPE(LL_TYPE)   :: LL_CNTL
SAVE
END MODULE READCNTL_MODULE
!
!     ..................................................................
      SUBROUTINE READCNTL
      USE LINKEDLIST_MODULE
      USE READCNTL_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
      LOGICAL(4)           :: TCHK
      INTEGER(4)           :: NFIL
      CHARACTER(32)        :: ID 
      CHARACTER(256)       :: FILENAME
      INTEGER(4)           :: ITH
      INTEGER(4)           :: NUM
      INTEGER(4)           :: NFILO
!     ****************************************************************
                          CALL TRACE$PUSH('READCNTL')
!
!     ==================================================================
!     ==  READ CONTROL FILE                                           ==
!     ==================================================================
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('DCNTL',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL)
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO) 
        CALL LINKEDLIST$REPORT(LL_CNTL,NFILO)
      END IF
!
!     ==================================================================
!     ==  !PDOSIN!FILES!FILE                                          ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILES',1,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',ITH)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
!       ==  READ ACTUAL VALUES  ======================================
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILENAME)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        CALL FILEHANDLER$SETFILE(ID,TCHK,FILENAME)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE READCNTL$GRID(EMIN,EMAX,NE,EBROAD)
      USE READCNTL_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: EMIN
      REAL(8)   ,INTENT(OUT) :: EMAX
      REAL(8)   ,INTENT(OUT) :: EBROAD
      INTEGER(4),INTENT(OUT) :: NE
      REAL(8)                :: DE
      REAL(8)                :: EV
      INTEGER(4)             :: NUM
      LOGICAL(4)             :: TCHK
!     ******************************************************************
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GRID')
!     ==  READ ACTUAL VALUES  ==========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMIN[EV]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EMIN[EV]',0,-20.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'EMIN[EV]',1,EMIN)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMAX[EV]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EMAX[EV]',0,5.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'EMAX[EV]',1,EMAX)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DE[EV]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'DE[EV]',0,1.D-3)
      CALL LINKEDLIST$GET(LL_CNTL,'DE[EV]',1,DE)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'BROADENING[EV]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'BROADENING[EV]',0,1.D-1)
      CALL LINKEDLIST$GET(LL_CNTL,'BROADENING[EV]',1,EBROAD)
!
      CALL CONSTANTS('EV',EV)
      NE=INT((EMAX-EMIN)/DE)+1
      DE=(EMAX-EMIN)/(NUM-1)
      EMIN=EMIN*EV
      EMAX=EMAX*EV
      DE=DE*EV
      EBROAD=EBROAD*EV
      RETURN
      END      
!
!     ..................................................................
      SUBROUTINE READCNTL$ORBITAL(NPRO,NAT,LMX,RPOS)
!     ******************************************************************
!     ******************************************************************
      USE READCNTL_MODULE
      USE ORBITALS_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NPRO
      INTEGER(4)   ,INTENT(IN) :: NAT
      INTEGER(4)   ,INTENT(IN) :: LMX(NAT)
      REAL(8)      ,INTENT(IN) :: RPOS(3,NAT)
      CHARACTER(32)            :: ORBITALNAME
      complex(8)               :: ORBITAL(NPRO)
      complex(8)               :: ORBITAL1(NPRO)
      REAL(8)                  :: SVAR
      INTEGER(4)               :: IORB,ITH,I
      INTEGER(4)               :: NUM
      INTEGER(4)               :: NORB
!     ******************************************************************
                          CALL TRACE$PUSH('READCNTL$ORBITAL')
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'ORBITAL',NORB)
      DO IORB=1,NORB
        CALL LINKEDLIST$SELECT(LL_CNTL,'ORBITAL',IORB)
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,ORBITALNAME)
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB',NUM)
        ORBITAL(:)=0.D0
        DO ITH=1,NUM
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB',ITH)
          CALL READONEORB(LL_CNTL,NAT,LMX,RPOS,NPRO,ORBITAL1)
          ORBITAL(:)=ORBITAL(:)+ORBITAL1(:)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL NORMALIZEORBITAL(NPRO,ORBITAL)
        CALL ORBITALS$SETORB(ORBITALNAME,NPRO,ORBITAL)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE READCNTL$SETNUMBER(NSET)
!     ******************************************************************
!     ******************************************************************
      USE READCNTL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NSET
      INTEGER(4)             :: I
!     ******************************************************************
                          CALL TRACE$PUSH('READCNTL$SETNUMBER')
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'COOP',I)
      NSET=I
      CALL LINKEDLIST$NLISTS(LL_CNTL,'WEIGHT',I)
      NSET=NSET+I
                          CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE READCNTL$SETS(NB,NKPT,NSPIN,NSET &
     &                        ,NAT,LMX,RPOS,LENG,SET,LEGEND)
      USE READCNTL_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)  :: NB
      INTEGER(4)   ,INTENT(IN)  :: NKPT
      INTEGER(4)   ,INTENT(IN)  :: NSPIN
      INTEGER(4)   ,INTENT(IN)  :: NSET
      INTEGER(4)   ,INTENT(IN)  :: NAT
      INTEGER(4)   ,INTENT(IN)  :: LMX(NAT)
      INTEGER(4)   ,INTENT(IN)  :: LENG
      REAL(8)      ,INTENT(IN)  :: RPOS(3,NAT)
      REAL(8)      ,INTENT(OUT) :: SET(NB,NKPT,NSPIN,NSET)
      CHARACTER(32),INTENT(OUT) :: LEGEND(NSET)
      INTEGER(4)                :: ISET   ! SET COUNTER
      CHARACTER(32)             :: ORBITALNAME
      CHARACTER(8)              :: TYPE
      CHARACTER(32)             :: NAME
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: ITH, NUM
      INTEGER(4)                :: JTH, NUMJTH
      INTEGER(4)                :: IORB1,IORB2,IORB
      INTEGER(4)                :: NORB1,NORB2,NORB
      INTEGER(4)                :: I,IB,IKPT,ISPIN
      COMPLEX(8)                :: ORBITAL1(LENG)
      COMPLEX(8)                :: ORBITAL2(LENG)
      COMPLEX(8)                :: ORBITALI(LENG)
      REAL(8)                   :: SVAR1,SVAR2
      INTEGER(4)                :: IAT,IPOS
!     ******************************************************************
                          CALL TRACE$PUSH('READCNTL$SETS')
      SET(:,:,:,:)=0.D0
      ISET=0
!
!     ==================================================================
!     ==================================================================
!     ==  SCAN COOPS                                                  ==
!     ==================================================================
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'COOP',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNTL,'COOP',ITH)
        ISET=ISET+1
!
!       == LOOK UP ORBITALS ============================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB1',NORB1)
        ORBITAL1=0.D0 
        DO IORB1=1,NORB1
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB1',IORB1) 
          CALL READONEORB(LL_CNTL,NAT,LMX,RPOS,LENG,ORBITALI)
          ORBITAL1(:)=ORBITAL1(:)+ORBITALI(:)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB2',NORB2)
        ORBITAL2=0.D0 
        DO IORB2=1,NORB2
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB2',IORB2)
          CALL READONEORB(LL_CNTL,NAT,LMX,RPOS,LENG,ORBITALI)
          ORBITAL2(:)=ORBITAL2(:)+ORBITALI(:)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!       PRINT*,'SET ',SET(:,:,:,ISET)

        CALL SET$PROJECT(LENG,NB,NKPT,NSPIN,ORBITAL1,ORBITAL2,SET(:,:,:,ISET))
        WRITE(LEGEND(ISET),FMT='("COOP:SET",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ID',0,LEGEND(ISET))
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,LEGEND(ISET))
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$PASS('COOPS FINISHED')

!
!     ==================================================================
!     ==================================================================
!     ==  SCAN WEIGHTS                                                ==
!     ==================================================================
!     ==================================================================
      CALL LINKEDLIST$NLISTS(LL_CNTL,'WEIGHT',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNTL,'WEIGHT',ITH)
        ISET=ISET+1
        SET(:,:,:,ISET)=0.D0
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TYPE',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
!         ==============================================================
!         ==  'TOTAL' = TOTAL DENSITY OF STATES                       ==
!         ==============================================================
          IF(TRIM(TYPE).EQ.'TOTAL') THEN
             SET(:,:,:,ISET)=1.D0
!         ==============================================================
!         ==  'ALL' = ALL PROJECTED DENSITY OF STATES                 ==
!         ==============================================================
          ELSE IF(TRIM(TYPE).EQ.'ALL') THEN
             CALL SET$WEIGHT('ALL','',NB,NKPT,NSPIN,SET(1,1,1,ISET))
!         ==============================================================
!         ==  'EMPTY' = VACCUM DENSITY OF STATES                      ==
!         ==============================================================
          ELSE IF(TRIM(TYPE).EQ.'EMPTY') THEN
             CALL SET$WEIGHT('ALL','',NB,NKPT,NSPIN,SET(1,1,1,ISET))
             SET(:,:,:,ISET)=1.D0-SET(:,:,:,ISET)
!         ==============================================================
!         ==  'ALL' ALL PRORJECTED DENSITY OF STATES                  ==
!         ==============================================================
          ELSE 
            CALL ERROR$MSG('TYPE NOT RECOGNIZED')
            CALL ERROR$MSG('MUST BE EITHER TOTAL,ALL OR EMPTY')
            CALL ERROR$CHVAL('TYPE ',TYPE)
            CALL ERROR$STOP('')
          END IF
        END IF
!
!       ================================================================
!       == COLLECT CONTRIBUTIONS FROM INDIVIDUAL ATOMS                ==
!       ================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ATOM',NUMJTH)
        DO JTH=1,NUMJTH
          CALL LINKEDLIST$SELECT(LL_CNTL,'ATOM',JTH)
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,NAME)
!
!         == SELECT TYPE ============================================
          CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
          IF(SCAN(TYPE,'ALL').NE.0) THEN
            CALL SET$WEIGHT(NAME,'ALL',NB,NKPT,NSPIN,SET(1,1,1,ISET))
          ELSE
            IF(SCAN(TYPE,'S').NE.0) THEN
              CALL SET$WEIGHT(NAME,'S',NB,NKPT,NSPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'P').NE.0) THEN
              CALL SET$WEIGHT(NAME,'P',NB,NKPT,NSPIN,SET(1,1,1,ISET)) 
            END IF
            IF(SCAN(TYPE,'D').NE.0) THEN
              CALL SET$WEIGHT(NAME,'D',NB,NKPT,NSPIN,SET(1,1,1,ISET))
            END IF
            IF(SCAN(TYPE,'F').NE.0) THEN
              CALL SET$WEIGHT(NAME,'F',NB,NKPT,NSPIN,SET(1,1,1,ISET))
            END IF
          END IF
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
!       ================================================================
!       == ADD CONTRIBUTION FROM PREDEFINED ORBITALS ===================
!       ================================================================
!
!       == LOOK UP ORBITALS ============================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ORB',NORB)
        DO IORB=1,NORB
          CALL LINKEDLIST$SELECT(LL_CNTL,'ORB',IORB)
          CALL READONEORB(LL_CNTL,NAT,LMX,RPOS,LENG,ORBITALI)
          call SET$PROJECT(LENG,NB,NKPT,NSPIN,ORBITALI,ORBITALI,SET(1,1,1,iset))
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        WRITE(LEGEND(ISET),FMT='("WEIGHT",I5)')ISET
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ID',0,LEGEND(ISET))
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,LEGEND(ISET))
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SET$WEIGHT(ATOMID_,ORBITALID,NB_,NKPT_,NSPIN_,SET)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ATOMID_
      CHARACTER(*),INTENT(IN) :: ORBITALID
      INTEGER(4)  ,INTENT(IN) :: NB_
      INTEGER(4)  ,INTENT(IN) :: NKPT_
      INTEGER(4)  ,INTENT(IN) :: NSPIN_
      REAL(8)     ,INTENT(OUT):: SET(NB_,NKPT_,NSPIN_)
      INTEGER(4)              :: ISPIN,IKPT,IB,IDIM
      INTEGER(4)              :: IPRO,IPRO0,IPRO1,IPRO2
      INTEGER(4)              :: IAT,IAT0,ISP
      INTEGER(4)              :: L,L1,L2,M,LN,LN1,LN2
      REAL(8)                 :: SUM
!     ******************************************************************
                                 CALL TRACE$PUSH('SET$WEIGHT')
      SET(:,:,:)=0.D0
!
!     ==================================================================
!     ==  SELECT ATOM                                                 ==
!     ==================================================================
      IF(ATOMID_.EQ.'ALL'.OR.ATOMID_.EQ.'') THEN
        IAT0=-1
      ELSE
        IAT0=-1
        DO IAT=1,NAT
          IF(ATOMID_.EQ.ATOMID(IAT)) THEN
            IAT0=IAT
            EXIT
          END IF
        ENDDO
        IF(IAT0.EQ.-1) THEN
          CALL ERROR$MSG('ATOMID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ATOMID_',ATOMID_)
          CALL ERROR$STOP('SET$WEIGHT')
        END IF
      END IF
!
!     ==================================================================
!     ==  SELECT ORBITAL                                              ==
!     ==================================================================
      IF(ORBITALID.EQ.'S') THEN
        L=0
      ELSE IF(ORBITALID.EQ.'P') THEN
        L=1
      ELSE IF(ORBITALID.EQ.'D') THEN
        L=2
      ELSE IF(ORBITALID.EQ.'F') THEN
        L=3
      ELSE IF(ORBITALID.EQ.'ALL'.OR.ORBITALID.EQ.'') THEN
        L=-1
      ELSE
        CALL ERROR$MSG('ORBITALID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ORBITALID',ORBITALID)
        CALL ERROR$STOP('SET$WEIGHT')
      END IF
!
!     ==================================================================
!     ==  NOW EVALUATE WEIGHTS                                        ==
!     ==================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          IPRO0=0
          DO IAT=1,NAT
            ISP=ISPECIES(IAT)
            IF(IAT.EQ.IAT0.OR.IAT0.EQ.-1) THEN
              IPRO1=IPRO0
              DO LN1=1,LNX(ISP)
                L1=LOX(LN1,ISP)
                IPRO2=IPRO0
                DO LN2=1,LNX(ISP)
                  L2=LOX(LN2,ISP)
                  IF(L1.NE.L2.OR.(L.GE.0.AND.L1.NE.L)) THEN
                    IPRO2=IPRO2+2*L2+1
                    CYCLE
                  END IF
                  DO IB=1,STATE%NB
                    SUM=0.D0
                    DO M=1,2*L1+1
                      DO IDIM=1,NDIM
                        SUM=SUM+REAL(CONJG(STATE%VEC(IDIM,IPRO1+M,IB)) &
     &                                    *STATE%VEC(IDIM,IPRO2+M,IB))
                      ENDDO
                    ENDDO
                    SET(IB,IKPT,ISPIN)=SET(IB,IKPT,ISPIN)+SUM*OV(LN1,LN2,ISP)
                  ENDDO
                  IPRO2=IPRO2+2*L2+1
                ENDDO            
                IPRO1=IPRO1+2*L1+1
              ENDDO
            END IF
            DO LN=1,LNX(ISP)
              IPRO0=IPRO0+2*LOX(LN,ISP)+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                                 CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SET$PROJECT(NPRO_,NB_,NKPT_,NSPIN_,ORBITAL1,ORBITAL2,SET)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NPRO_
      INTEGER(4),INTENT(IN) :: NB_
      INTEGER(4),INTENT(IN) :: NKPT_
      INTEGER(4),INTENT(IN) :: NSPIN_
      COMPLEX(8),INTENT(IN) :: ORBITAL1(NPRO)
      COMPLEX(8),INTENT(IN) :: ORBITAL2(NPRO)
      REAL(8)   ,INTENT(OUT):: SET(NB_,NKPT_,NSPIN_)
      INTEGER(4)            :: ISPIN,IKPT,IB,Ipro,IDIM
      COMPLEX(8)            :: CSVAR1,csvar2
!     ******************************************************************
                                 CALL TRACE$PUSH('SET$PROJECT')
      SET(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          STATE=>STATEARR(IKPT,ISPIN)
          DO IB=1,NB_
            CSVAR1=0.D0
            CSVAR2=0.D0
            DO IPRO=1,NPRO
              DO IDIM=1,NDIM
                CSVAR1=CSVAR1+CONJG(ORBITAL1(IPRO))*STATE%VEC(IDIM,ipro,IB)
                CSVAR2=CSVAR2+CONJG(ORBITAL2(IPRO))*STATE%VEC(IDIM,ipro,IB)
              ENDDO                      
            ENDDO
            SET(IB,IKPT,ISPIN)=REAL(CONJG(CSVAR1)*CSVAR2)
          ENDDO
        ENDDO
      ENDDO
                                 CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE READCNTL$OUTPUT(EMIN,EMAX,NE,EBROAD &
     &                    ,NB,NKPT,NSPIN,EIG,NSET,SET,LEGEND)
!     ******************************************************************
!     ******************************************************************
      USE READCNTL_MODULE
      USE ORBITALS_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NE
      REAL(8)      ,INTENT(IN) :: EMIN
      REAL(8)      ,INTENT(IN) :: EMAX
      REAL(8)      ,INTENT(IN) :: EBROAD
      INTEGER(4)   ,INTENT(IN) :: NB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      INTEGER(4)   ,INTENT(IN) :: NSPIN
      REAL(8)      ,INTENT(IN) :: EIG(NB,NKPT,NSPIN)
      INTEGER(4)   ,INTENT(IN) :: NSET
      REAL(8)      ,INTENT(IN) :: SET(NB,NKPT,NSPIN,NSET)
      CHARACTER(32),INTENT(IN) :: LEGEND(NSET)
      CHARACTER(32)        :: LEGEND1
      CHARACTER(256)       :: FILE
      LOGICAL(4)           :: TIB,TE,TIK,TIS,TCHK
      INTEGER(4)           :: IB,IKPT,ISPIN,I,IOUT,ITH
      INTEGER(4)           :: IK,IS
      INTEGER(4)           :: ISET
      INTEGER(4)           :: NOUT
      INTEGER(4)           :: NFIL
      INTEGER(4)           :: NFILO
      INTEGER(4)           :: IB0,IK0,IS0
      REAL(8)              :: EV
      REAL(8)              :: ENERGY
      REAL(8)              :: SUM,SUMS,SVAR 
!     ******************************************************************
                         CALL TRACE$PUSH('READCNTL$OUTPUT')
      CALL CONSTANTS('EV',EV)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'OUTPUT',NOUT)
      DO IOUT=1,NOUT
        CALL LINKEDLIST$SELECT(LL_CNTL,'OUTPUT',IOUT)
                          CALL TRACE$PASS('NEXT IOUT')
!
!       ==  SPECIFY SET ===============================================
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,LEGEND1)
        ISET=0
        DO I=1,NSET
          IF(TRIM(LEGEND1).EQ.TRIM(LEGEND(I))) THEN
            ISET=I ; EXIT
          END IF
        ENDDO
        IF(ISET.EQ.0) THEN
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          DO I=1,NSET
            PRINT*,'LEGEND ',I,TRIM(LEGEND(I)),'::'
          ENDDO
          CALL ERROR$CHVAL('LEGEND1',LEGEND1)
          CALL ERROR$STOP('READCNTL$OUTPUT')
        END IF
!
!       ==  SPECIFY OUTPUT FILE =======================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
        IF(TCHK) THEN        
          CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FILE)
          CALL FILEHANDLER$SETFILE('TEMP',.FALSE.,FILE)
          CALL FILEHANDLER$SETSPECIFICATION('TEMP','FORM','FORMATTED')
          CALL FILEHANDLER$SETSPECIFICATION('TEMP','POSITION','REWIND')
          CALL FILEHANDLER$UNIT('TEMP',NFIL)
          WRITE(NFILO,FMT='("OUTPUT WRITTEN FOR SET",A &
       &           ," IS WRITTEN TO FILE:"/A)') &
       &           TRIM(LEGEND(ISET)),TRIM(FILE)
        ELSE
          CALL FILEHANDLER$UNIT('PDOSOUT',NFIL)
        ENDIF
!
!       ==  OUTPUT TYPE ===============================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'B',1,TIB)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'E[EV]',1,TE)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'K',1,TIK)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'S',1,TIS)
        IK0=0
        IS0=0
        IF(TIK)CALL LINKEDLIST$GET(LL_CNTL,'K',1,IK0)
        IF(TIS)CALL LINKEDLIST$GET(LL_CNTL,'S',1,IS0)
        IF(TIB.AND.TE) THEN
          CALL ERROR$MSG('NB AND E ARE INCOMPATIBLE')
          CALL ERROR$STOP('READCNTL$OUTPUT')
        END IF
!
!       ================================================================
!       ==  WRITE PROJECTION OF A GIVEN STATE ON FILE                 ==
!       ================================================================
        IF(TIB) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'B',1,IB0)
          DO ISPIN=1,NSPIN
            IF(TIS.AND.(ISPIN.NE.IS0)) CYCLE
            DO IKPT=1,NKPT
              IF(TIK.AND.(IKPT.NE.IK0)) CYCLE
              WRITE(NFIL,FMT='(A32,"; B=",I3,"; K=",I2,"; S=",I1 &
     &                  ,"; E[EV]=",F10.5,"; PRO=",F10.5)') &
     &                LEGEND(ISET),IB0,IKPT,ISPIN &
     &               ,EIG(IB0,IKPT,ISPIN)/EV &
     &               ,SET(IB0,IKPT,ISPIN,ISET)
            ENDDO
          ENDDO
        END IF
!
!       ================================================================
!       ==  WRITE INTEGRATED DENSITY OF STATES AT A GIVEN ENERGY      ==
!       ================================================================
        IF(TIB) THEN
          ENERGY=-1.D+10
          SUM=0.D0
          SUMS=0.D0
          DO ISPIN=1,NSPIN
            IF(TIS.AND.(ISPIN.NE.IS0)) CYCLE
            IF(ISPIN.EQ.1) SVAR=1.D0
            IF(ISPIN.EQ.2) SVAR=-1.D0
            DO IKPT=1,NKPT
              IF(TIK.AND.(IKPT.NE.IK0)) CYCLE
              DO IB=1,IB0
                SUM=SUM+SET(IB,IKPT,ISPIN,ISET)
                SUMS=SUMS+SET(IB,IKPT,ISPIN,ISET)*SVAR
                ENERGY=MAX(ENERGY,EIG(IB,IKPT,ISPIN))
              ENDDO
            ENDDO        
          ENDDO
          IF(.NOT.TIS) SUM=SUM*2.D0/DBLE(NSPIN)
          WRITE(NFIL,FMT='(A32,"; B=",I3,"; K=",I2,"; S=",I1 &
     &                  ,"; E[EV]=",F10.5,"; NOS=",F10.5)') &
     &                  LEGEND(ISET),IB0,IK0,IS0,ENERGY/EV,SUM
!         IF(NSPIN.EQ.2.AND.(.NOT.TIS)) THEN
!           WRITE(NFIL,FMT='(A32,"; B=",I3,"; K=",I2,"; S=",I1 &
!    &                  ,";       ",10X,  "; SOS=",F10.5)') &
!    &                  LEGEND(ISET),IB0,IK0,IS0,SUMS
!         END IF
        END IF
!
!       ================================================================
!       ==  WRITE INTEGRATED DENSITY OF STATES AT A GIVEN ENERGY      ==
!       ================================================================
        IF(TE) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'E[EV]',1,ENERGY)
          ENERGY=ENERGY*EV
          SUM=0.D0
          SUMS=0.D0
          DO ISPIN=1,NSPIN
            IF(TIS.AND.(ISPIN.NE.IS0)) CYCLE
            IF(ISPIN.EQ.1) SVAR=1.D0
            IF(ISPIN.EQ.2) SVAR=-1.D0
            DO IKPT=1,NKPT
              IF(TIK.AND.(IKPT.NE.IK0)) CYCLE
              DO IB=1,NB
                IF(EIG(IB,IKPT,ISPIN).LT.ENERGY) THEN
                  SUM=SUM+SET(IB,IKPT,ISPIN,ISET)
                  SUMS=SUMS+SET(IB,IKPT,ISPIN,ISET)*SVAR
                END IF
              ENDDO
            ENDDO        
          ENDDO
          IF(.NOT.TIS) SUM=SUM*2.D0/DBLE(NSPIN)
          WRITE(NFIL,FMT='(A32,8X," K=",I2,"; S=",I1 &
     &                  ,"; E[EV]=",F10.5,"; NOS=",F10.5)') &
     &                  LEGEND(ISET),IK0,IS0,ENERGY/EV,SUM
          IF(NSPIN.EQ.2.AND.(.NOT.TIS)) THEN
            WRITE(NFIL,FMT='(A32,8X,"; K=",I2,"; S=",I1 &
     &                  ,"; E[EV]=",F10.5,"; SOS=",F10.5)') &
     &                  LEGEND(ISET),IK0,IS0,ENERGY/EV,SUMS
          END IF
        END IF
!
!       ================================================================
!       ==  WRITE DOS AND INTEGRATED DOS ON FILE                      ==
!       ================================================================
        IF(.NOT.(TIB.OR.TE)) THEN
          CALL PUTONGRID(NFIL,EMIN,EMAX,NE,EBROAD,NB,NKPT,NSPIN,EIG &
     &                    ,SET(:,:,:,ISET),LEGEND(ISET))
        END IF
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PUTONGRID(NFIL,EMIN,EMAX,NE,EBROAD,NB,NKPT,NSPIN,EIG &
     &                    ,SET,LEGEND)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4)   ,INTENT(IN) :: NE
      REAL(8)      ,INTENT(IN) :: EMIN
      REAL(8)      ,INTENT(IN) :: EMAX
      REAL(8)      ,INTENT(IN) :: EBROAD
      INTEGER(4)   ,INTENT(IN) :: NB
      INTEGER(4)   ,INTENT(IN) :: NKPT
      INTEGER(4)   ,INTENT(IN) :: NSPIN
      REAL(8)      ,INTENT(IN) :: EIG(NB,NKPT,NSPIN)
      REAL(8)      ,INTENT(IN) :: SET(NB,NKPT,NSPIN)
      INTEGER(4)   ,INTENT(IN) :: NFIL
      CHARACTER(32),INTENT(IN) :: LEGEND
      REAL(8)              :: DE
      INTEGER(4)           :: IE1,IE2,IDE
      INTEGER(4)           :: ND
      REAL(8)              :: NOS(NE,NSPIN)
      REAL(8)              :: DOS(NE,NSPIN)
      REAL(8)              :: EV
      REAL(8)              :: W1,W2,X,FAC
      REAL(8)              :: YMIN,YMAX
      REAL(8)              :: XNOS,XDOS
      CHARACTER(256)       :: CMD
!     ******************************************************************
                                 CALL TRACE$PUSH('PUTONGRID')
      CALL CONSTANTS('EV',EV)
      DE=(EMAX-EMIN)/DBLE(NE-1)
      ND=20
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      NOS(:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NB
            X=(EIG(IB,IKPT,ISPIN)-EMIN)/DE
            IE1=INT(X)
            IE2=IE1+1
            W2=(X-DBLE(IE1))
            W1=1.D0-W2
            IE1=MAX(IE1,1)
            IE2=MAX(IE2,1)
            IF(IE1.LE.NE)NOS(IE1,ISPIN)=NOS(IE1,ISPIN)+W1*SET(IB,IKPT,ISPIN)
            IF(IE2.LE.NE)NOS(IE2,ISPIN)=NOS(IE2,ISPIN)+W2*SET(IB,IKPT,ISPIN)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  RENORMALIZE WITH K-POINT WEIGHT AND SPIN MULTIPLICITY       ==
!     ==================================================================
      NOS(:,:)=NOS(:,:)/DBLE(NKPT)*2.D0/DBLE(NSPIN)
!
!     ==================================================================
!     ==  CALCULATE DOS                                               ==
!     ==================================================================
      FAC=0.D0
      DO IDE=-ND,ND
        FAC=FAC+EXP(-(DE*DBLE(IDE)/EBROAD)**2)
      ENDDO
      FAC=1.D0/FAC
      DOS(:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IDE=-ND,ND
          IE1=MAX(1,1-IDE)
          IE2=MIN(NE,NE-IDE)
          W1=FAC*EXP(-(DE*DBLE(IDE)/EBROAD)**2)
          DO IE=IE1,IE2
            DOS(IE,ISPIN)=DOS(IE,ISPIN)+NOS(IE+IDE,ISPIN)*W1
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  INTEGRATE NOS                                               ==
!     ==================================================================
      DO ISPIN=1,NSPIN
        DO IE=2,NE
          NOS(IE,ISPIN)=NOS(IE,ISPIN)+NOS(IE-1,ISPIN)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ROUND NOS AND DOS TO OBTAIN PROPER EDITING                  ==
!     ==================================================================
      DO ISPIN=1,NSPIN
        DO IE=1,NE
          IF(DABS(NOS(IE,ISPIN)).LE.1.D-99)NOS(IE,ISPIN)=0.D0
          IF(DABS(DOS(IE,ISPIN)).LE.1.D-99)DOS(IE,ISPIN)=0.D0
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ROUND NOS AND DOS TO OBTAIN PROPER EDITING                  ==
!     ==================================================================
      XDOS=0.D0
      XNOS=0.D0
      DO ISPIN=1,NSPIN
        DO IE=1,NE
          XDOS=MAX(XDOS,DABS(DOS(IE,ISPIN)))
          XNOS=MAX(XNOS,DABS(NOS(IE,ISPIN)))
        ENDDO
      ENDDO
      IF(XDOS.GT.0.D0) THEN
        FAC=XNOS/XDOS
      ELSE
        FAC=1.D0
      END IF
      DO ISPIN=1,NSPIN
        DO IE=1,NE
          DOS(IE,ISPIN)=DOS(IE,ISPIN)*FAC
          IF(DABS(NOS(IE,ISPIN)).LE.1.D-99)NOS(IE,ISPIN)=0.D0
          IF(DABS(DOS(IE,ISPIN)).LE.1.D-99)DOS(IE,ISPIN)=0.D0
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  WRITE RESULT ON PSODOUT                                     ==
!     ==================================================================
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)*DBLE(IE-1)/DBLE(NE-1)
        IF(NSPIN.EQ.1) THEN
          WRITE(NFIL,FMT='(F10.5,2E15.5)') &
     &          E/EV,DOS(IE,1),NOS(IE,1)
        ELSE
          WRITE(NFIL,FMT='(F10.5,4E15.5)') &
     &        E/EV,DOS(IE,1),NOS(IE,1) &
     &        ,-DOS(IE,2),-NOS(IE,2)
        END IF
      ENDDO
      WRITE(NFIL,FMT='("THIS WAS: ",A)')LEGEND
                                 CALL TRACE$POP
      RETURN
      END
