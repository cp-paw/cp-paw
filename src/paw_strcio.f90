!PROGRAM TEST
!IMPLICIT NONE
!INTEGER(4)                :: NAT
!REAL(8)                   :: RBAS(3,3)
!REAL(8)      ,ALLOCATABLE :: R(:,:)
!CHARACTER(10),ALLOCATABLE :: NAME(:)
!!***********************************************************************
!NAT=2
!ALLOCATE(NAME(NAT))
!ALLOCATE(R(3,NAT))
!R(:,1)=(/1.D0,0.D0,0.D0/)
!R(:,2)=(/-1.D0,0.D0,0.D0/)
!NAME=(/'ZR1','SI3'/)
!CALL STRCIO$CLEAR
!CALL STRCIO$SETR8A('R',3*NAT,R)
!CALL STRCIO$SETCHA('NAME',NAT,NAME)
!CALL STRCIO$WRITE(6)
!DEALLOCATE(NAME)
!DEALLOCATE(R)
!PRINT*,'=========================================='
!!
!CALL STRCIO$CLEAR
!CALL STRCIO$READ(5)
!CALL STRCIO$GETI4('NAT',NAT)
!PRINT*,'NAT ',NAT
!ALLOCATE(NAME(NAT))
!ALLOCATE(R(3,NAT))
!CALL STRCIO$GETCHA('NAME',NAT,NAME)
!CALL STRCIO$GETR8A('R',3*NAT,R)
!CALL STRCIO$GETR8A('RBAS',9,RBAS)
!!
!CALL STRCIO$CLEAR
!CALL STRCIO$SETR8A('RBAS',9,RBAS)
!CALL STRCIO$SETR8A('R',3*NAT,R)
!CALL STRCIO$SETCHA('NAME',NAT,NAME)
!CALL STRCIO$WRITE(6)
!STOP
!END
!***********************************************************************
!***********************************************************************
!**  OBJECT STRCIO                                                    **
!**                                                                   **
!**  WRITES AND READS ATOMIC STRUCTURES IN VARIOUS FILE FORMATS       **
!**  COMPLETES INFORMATION WHEN NECCESARY                             **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    CLEAR        DEALLOCATES INTERNAL DATA                         **
!**    READ(NFIL)   READS A DATA FILE                                 **
!**    WRITE(NFIL)  WRITES A DATA FILE                                **
!**    SETI4(ID,VAL)  ID=(NAT)                                        **
!**    SETI4A(ID,LEN,VAL) /GETI4A(ID,LEN,VAL)                         **
!**      ID='IZ'                                                      **
!**      ID='NEIGHBORS'                                               **
!**      ID='IT(NEIGHBORS)'                                           **
!**    SETI8A(ID,LEN,VAL) /GETR8A(ID,LEN,VAL)                         **
!**      ID='RBAS'     LATTICE VECTORS                                **
!**      ID='R'        ATOMIC COORIDNATES                             **
!**      ID='Q'        ATOMIC CHARGES                                 **
!**    SETCHA(ID,LEN,VAL) /GETCHA(ID,LEN,VAL)                         **
!**      ID='NAME'     ATOM NAME                                      **
!**                                                                   **
!**  REMARKS:                                                         **
!**      NOT FULLY COMPLETED. CURRENTLY WORKS ONLY WITH CSSR FILES    **
!**                                                                   **
!***********************************************************************
!***********************************************************************
!.......................................................................
MODULE STRCIO_MODULE
INTERFACE OPERATOR (.DYAD.)
  MODULE PROCEDURE DYADISCHES_PRODUCT
END INTERFACE 
REAL(8)                   :: RBAS(3,3)
INTEGER(4)                :: NAT
CHARACTER(32),ALLOCATABLE :: NAME(:)
INTEGER(4)   ,ALLOCATABLE :: IZ(:)
REAL(8)      ,ALLOCATABLE :: R(:,:)
REAL(8)      ,ALLOCATABLE :: Q(:)
INTEGER(4)                :: NNEIGH
INTEGER(4)   ,ALLOCATABLE :: NEIGH(:,:)
INTEGER(4)   ,ALLOCATABLE :: ITNEIGH(:,:,:)
LOGICAL(4)                :: TCRYSTAL
CHARACTER(16)             :: FORMAT
CHARACTER(128)            :: OBJECTNAME
CONTAINS
!     ..................................................................
      FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
        REAL(8), INTENT(IN) :: R1(3)
        REAL(8), INTENT(IN) :: R2(3)
        REAL(8)             :: R3(3)
        R3(1)=R1(2)*R2(3)-R1(3)*R2(2)
        R3(2)=R1(3)*R2(1)-R1(1)*R2(3)
        R3(3)=R1(1)*R2(2)-R1(2)*R2(1)
      END FUNCTION DYADISCHES_PRODUCT
END MODULE STRCIO_MODULE
!
!     ..................................................................
      SUBROUTINE STRCIO$CLEAR
      USE STRCIO_MODULE
      IMPLICIT NONE
!     ******************************************************************
      TCRYSTAL=.FALSE.
      FORMAT='CSSR'
      OBJECTNAME='NONAME'
      IF(ALLOCATED(NAME)) DEALLOCATE(NAME)
      IF(ALLOCATED(IZ))   DEALLOCATE(IZ)
      IF(ALLOCATED(R))    DEALLOCATE(R)
      IF(ALLOCATED(Q))    DEALLOCATE(Q)
      IF(ALLOCATED(NEIGH))DEALLOCATE(NEIGH)
      IF(ALLOCATED(ITNEIGH))DEALLOCATE(ITNEIGH)
      NAT=0
      NNEIGH=0
      RBAS(:,:)=0.D0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$SETI4(ID,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NAT') THEN
        NAT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$SETI4A(ID,LEN,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'IZ') THEN
        IF(NAT.EQ.0) NAT=LEN
        IF(NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$SETI4A')
        END IF
        IF(ALLOCATED(IZ))DEALLOCATE(IZ)
        ALLOCATE(IZ(NAT))
        IZ(:)=VAL(:)
      ELSE IF(ID.EQ.'NEIGHBORS') THEN
        IF(NAT.EQ.0) THEN
          CALL ERROR$MSG('#(ATOMS) NOT KNOWN')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('STRCIO$SETI4A')
        END IF
        NNEIGH=LEN/NAT
        IF(NAT*NNEIGH.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$SETI4A')
        END IF
        IF(ALLOCATED(NEIGH))DEALLOCATE(NEIGH)
        ALLOCATE(NEIGH(NNEIGH,NAT))
        NEIGH(:,:)=RESHAPE(VAL,(/NNEIGH,NAT/))
      ELSE IF(ID.EQ.'IT(NEIGHBORS)') THEN
        IF(NAT.EQ.0) THEN
          CALL ERROR$MSG('#(ATOMS) NOT KNOWN')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('STRCIO$SETI4A')
        END IF
        IF(NNEIGH.EQ.0) THEN
          CALL ERROR$MSG('X#(NEIGHBORS) NOT KNOWN')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('STRCIO$SETI4A')
        END IF
        IF(3*NAT*NNEIGH.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$SETI4A')
        END IF
        IF(ALLOCATED(ITNEIGH))DEALLOCATE(ITNEIGH)
        ALLOCATE(ITNEIGH(3,NNEIGH,NAT))
        ITNEIGH(:,:,:)=RESHAPE(VAL,(/3,NNEIGH,NAT/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$SETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$GETI4A(ID,LEN,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)              :: NNEIGH1,IAT,IN,I
!     ******************************************************************
      IF(ID.EQ.'IZ') THEN
        IF(NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$GETI4A')
        END IF
        VAL(:)=IZ(:)
      ELSE IF(ID.EQ.'NEIGHBORS') THEN
        NNEIGH1=LEN/NAT
        IF(NAT*NNEIGH1.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$GETI4A')
        END IF
        VAL(:)=0
        DO IAT=1,NAT
          DO IN=1,NNEIGH
            IF(IN.GT.NNEIGH1) CYCLE
            I=IN+NNEIGH1*(IAT-1)
            VAL(I)=NEIGH(IN,IAT)
          ENDDO
        ENDDO
      ELSE IF(ID.EQ.'IT(NEIGHBORS)') THEN
        NNEIGH1=LEN/(3*NAT)
        IF(3*NAT*NNEIGH1.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$GETI4A')
        END IF
        VAL(1:LEN)=0
        DO IAT=1,NAT
          DO IN=1,NNEIGH
            IF(IN.GT.NNEIGH1) CYCLE
            I=3*(IN-1+NNEIGH1*(IAT-1))
            VAL(I+1)=ITNEIGH(1,IN,IAT)
            VAL(I+2)=ITNEIGH(2,IN,IAT)
            VAL(I+3)=ITNEIGH(3,IN,IAT)
          ENDDO
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$GETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$SETR8A(ID,LEN,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'RBAS') THEN
        IF(9.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('STRCIO$SETR8A')
        END IF
        RBAS(:,:)=RESHAPE(VAL,(/3,3/))
        TCRYSTAL=.TRUE.
      ELSE IF(ID.EQ.'R') THEN
        IF(NAT.EQ.0) NAT=LEN/3
        IF(3*NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$SETR8A')
        END IF
        IF(ALLOCATED(R))DEALLOCATE(R)
        ALLOCATE(R(3,NAT))
        R(:,:)=RESHAPE(VAL,(/3,NAT/))
      ELSE IF(ID.EQ.'Q') THEN
        IF(NAT.EQ.0) NAT=LEN
        IF(NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$SETR8A')
        END IF
        IF(ALLOCATED(Q))DEALLOCATE(Q)
        ALLOCATE(Q(NAT))
        Q(:)=VAL(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$SETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$GETR8A(ID,LEN,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'RBAS') THEN
        IF(9.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('STRCIO$GETR8A')
        END IF
        VAL=RESHAPE(RBAS,(/9/))
      ELSE IF(ID.EQ.'R') THEN
        IF(3*NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(R)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('STRCIO$GETR8A')
        END IF
        VAL=RESHAPE(R,(/3*NAT/))
      ELSE IF(ID.EQ.'Q') THEN
        IF(NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(Q)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('STRCIO$GETR8A')
        END IF
        VAL(:)=Q(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$GETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$SETCHA(ID,LEN,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'NAME') THEN
        IF(NAT.EQ.0) NAT=LEN
        IF(NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$SETI4A')
        END IF
        IF(ALLOCATED(NAME))DEALLOCATE(NAME)
        ALLOCATE(NAME(NAT))
        NAME(:)=VAL(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$SETCHA')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$GETCHA(ID,LEN,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(OUT):: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'NAME') THEN
        IF(NAT.NE.LEN) THEN 
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('STRCIO$GETI4A')
        END IF
        IF(.NOT.ALLOCATED(Q)) THEN
          CALL ERROR$MSG('NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('STRCIO$GETCHA')
        END IF
        VAL(:)=NAME(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$GETCHA')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$GETI4(ID,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'NAT') THEN
        VAL=NAT
      ELSE IF(ID.EQ.'NN') THEN
        VAL=NNEIGH
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$SETCH(ID,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'FORMAT') THEN
        FORMAT=VAL
      ELSE IF(ID.EQ.'OBJECTID') THEN
        OBJECTNAME=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$SETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$GETCH(ID,VAL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'FORMAT') THEN
        VAL=FORMAT
      ELSE IF(ID.EQ.'OBJECTID') THEN
        VAL=OBJECTNAME
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('STRCIO$SETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$WRITE(NFIL)
      USE STRCIO_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN):: NFIL
!     ******************************************************************
      CALL STRCIO_COMPLETE
      IF(FORMAT.EQ.'XYZ') THEN
        CALL ERROR$MSG('XYZ FORMAT NOT YET IMPLEMENTED')
        CALL ERROR$MSG('STRCIO$WRITE')
      ELSE
        CALL STRCIO_WRITECSSR(NFIL)
      END IF
      RETURN 
      END
!
!     ..................................................................
      SUBROUTINE STRCIO$READ(NFIL)
!     ******************************************************************
!     **                                                              **
!     **  READS A STRUCTURE FILE                                      **
!     **                                                              **
!     ******************************************************************
      USE STRCIO_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN):: NFIL
!     ******************************************************************
!
!     ==================================================================
!     == CLEAR ALL LOCAL ARRAYS                                       ==
!     ==================================================================
      CALL STRCIO$CLEAR
!
!     ==================================================================
!     == SELECT FILE AND READ                                         ==
!     ==================================================================
      IF(FORMAT.EQ.'XYZ') THEN
        CALL ERROR$MSG('XYZ FORMAT NOT YET IMPLEMENTED')
        CALL ERROR$MSG('STRCIO$WRITE')
      ELSE
        CALL STRCIO_READCSSR(NFIL)
      END IF
!
!     ==================================================================
!     == COMPLETES MISSING INFORMATION                                ==
!     ==================================================================
      CALL STRCIO_COMPLETE
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO_COMPLETE
      USE PERIODICTABLE_MODULE
      USE STRCIO_MODULE
      IMPLICIT NONE
      INTEGER(4)                :: IAT
      CHARACTER(32)             :: STRING
      REAL(8)      ,ALLOCATABLE :: RAD(:)
      REAL(8)                   :: T(3)
      INTEGER(4)                :: DI
      REAL(8)                   :: DR(3)
      REAL(8)                   :: DIS
      REAL(8)                   :: DISMAX
      INTEGER(4)   ,ALLOCATABLE :: INN(:)
      INTEGER(4)                :: I1,I2,I3,IAT1,IAT2,IN
      REAL(8)                   :: ANGSTROM
      REAL(8)                   :: SQR2      !SQRT(2)
      REAL(8)                   :: ALAT
      LOGICAL(4)                :: TCHK
!     ******************************************************************
      SQR2=SQRT(2.D0)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
!
!     ==================================================================
!     ==  CHECK IF MANDATORY DATA ARE PRESENT                         ==
!     ==================================================================
      IF(NAT.EQ.0) THEN
        CALL ERROR$MSG('#(ATOMS) UNKNOWN')
        CALL ERROR$STOP('STRCIO_COMPLETE')
      END IF
!
      IF(.NOT.ALLOCATED(R)) THEN
        CALL ERROR$MSG('ATOMIC POSITIONS UNKNOWN')
        CALL ERROR$STOP('STRCIO_COMPLETE')
      END IF
!
      IF((.NOT.ALLOCATED(IZ)).AND.(.NOT.ALLOCATED(NAME))) THEN
        CALL ERROR$MSG('INFORMATION ON ELEMENTS MUST BE GIVEN EITHER')
        CALL ERROR$MSG('BY ATOMIC NUMBER OR NAME')
        CALL ERROR$STOP('STRCIO_COMPLETE')
      END IF   
!
!     ==================================================================
!     == SET DEFAULT ATOM NAMES                                       ==
!     ==================================================================
      IF(.NOT.ALLOCATED(NAME)) THEN
        ALLOCATE(NAME(NAT))
        DO IAT=1,NAT
          WRITE(STRING,FMT=*)IAT
          CALL PERIODICTABLE$GET(IZ(IAT),'SYMBOL',NAME(IAT))
          NAME(IAT)=NAME(IAT)(1:2)//ADJUSTL(STRING)
        ENDDO
      END IF 
!
!     ==================================================================
!     == SET DEFAULT ATOMIC NUMBER                                    ==
!     ==================================================================
      IF(.NOT.ALLOCATED(IZ)) THEN
        ALLOCATE(IZ(NAT))
        DO IAT=1,NAT
          CALL PERIODICTABLE$GET(NAME(IAT)(1:2),'Z',IZ(IAT))
        ENDDO
      END IF
!
!     ==================================================================
!     == SET DEFAULT CHARGES                                          ==
!     ==================================================================
      IF(.NOT.ALLOCATED(Q)) THEN
        ALLOCATE(Q(NAT))
        Q(:)=0.D0
      END IF
!
!     ==================================================================
!     == SET DEFAULT LATTICE CONSTANTS                                ==
!     ==================================================================
      IF(.NOT.TCRYSTAL) THEN
        DISMAX=0.D0
        DO IAT1=1,NAT
          DO IAT2=IAT1+1,NAT
            DR(:)=R(:,IAT2)+T(:)-R(:,IAT1)
            DIS=SQRT(DR(1)**2+DR(2)**2+DR(3)**2)
            IF(DIS.LT.DISMAX) CYCLE
            DISMAX=DIS
            ALAT=SQR2*(DIS+6.D0*ANGSTROM)
          ENDDO
        ENDDO
        RBAS(:,:)=0.5D0*ALAT
        RBAS(1,1)=0.D0
        RBAS(2,2)=0.D0
        RBAS(3,3)=0.D0
      END IF
!
!     ==================================================================
!     == SET DEFAULT NEIGHBORLIST                                     ==
!     ==================================================================
!     == RESET IF NEIGBORLIST IS PRESENT BUT EMPTY =====================
      IF(ALLOCATED(NEIGH)) THEN
        TCHK=.FALSE.
        DO IAT=1,NAT
          TCHK=NEIGH(1,IAT).NE.0
          IF(TCHK) EXIT
        ENDDO
        IF(.NOT.TCHK) THEN
          DEALLOCATE(NEIGH)
          IF(ALLOCATED(ITNEIGH))DEALLOCATE(ITNEIGH)
        END IF
      END IF

      IF(.NOT.ALLOCATED(NEIGH)) THEN
        NNEIGH=12
        ALLOCATE(NEIGH(NNEIGH,NAT))
        NEIGH(:,:)=0
        IF(ALLOCATED(ITNEIGH))DEALLOCATE(ITNEIGH)
        ALLOCATE(ITNEIGH(3,NNEIGH,NAT))
        ITNEIGH(:,:,:)=0
        ALLOCATE(INN(NAT))
        INN(:)=0
        ALLOCATE(RAD(NAT))
        DO IAT=1,NAT
          CALL PERIODICTABLE$GET(IZ(IAT),'R(COV)',RAD(IAT))
          RAD(IAT)=RAD(IAT)*1.2D0
        ENDDO   
        DI=0
        IF(TCRYSTAL) DI=1
        DO I1=-DI,DI
          DO I2=-DI,DI
            DO I3=-DI,DI
              T(:)=RBAS(:,1)*DBLE(I1)+RBAS(:,2)*DBLE(I2)+RBAS(:,3)*DBLE(I3)
              DO IAT1=1,NAT
                IAT=IAT1
                IF(I1.EQ.0.AND.I2.EQ.0.AND.I3.EQ.0) IAT=IAT1+1
                DO IAT2=IAT,NAT
                  DISMAX=RAD(IAT1)+RAD(IAT2)
                  DR(:)=R(:,IAT2)+T(:)-R(:,IAT1)
                  DIS=SQRT(DR(1)**2+DR(2)**2+DR(3)**2)
                  IF(DIS.GT.DISMAX) CYCLE
                  INN(IAT1)=INN(IAT1)+1
                  IF(INN(IAT1).GT.NNEIGH) THEN
                    CALL ERROR$MSG('TOO MANY NEIGHBORS')
                    CALL ERROR$I4VAL('#(NEIGHBORS) ALLOWED',NNEIGH)
                    CALL ERROR$I4VAL('IAT',IAT1)
                    CALL ERROR$STOP('STRCIO_COMPLETE')
                  END IF
                  INN(IAT2)=INN(IAT2)+1
                  IF(INN(IAT2).GT.NNEIGH) THEN
                    CALL ERROR$MSG('TOO MANY NEIGHBORS')
                    CALL ERROR$I4VAL('#(NEIGHBORS) ALLOWED',NNEIGH)
                    CALL ERROR$I4VAL('IAT',IAT2)
                    CALL ERROR$STOP('STRCIO_COMPLETE')
                  END IF
                  NEIGH(INN(IAT1),IAT1)=IAT2            
                  ITNEIGH(1,INN(IAT1),IAT1)=I1
                  ITNEIGH(2,INN(IAT1),IAT1)=I2
                  ITNEIGH(3,INN(IAT1),IAT1)=I3
                  NEIGH(INN(IAT2),IAT2)=IAT1            
                  ITNEIGH(1,INN(IAT2),IAT2)=-I1
                  ITNEIGH(2,INN(IAT2),IAT2)=-I2
                  ITNEIGH(3,INN(IAT2),IAT2)=-I3
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(RAD)
        DEALLOCATE(INN)
      ELSE 
        IF(.NOT.ALLOCATED(ITNEIGH)) THEN
          ALLOCATE(ITNEIGH(3,NNEIGH,NAT))
          ITNEIGH(:,:,:)=0
          IF(TCRYSTAL) THEN
            DO IAT=1,NAT
              DO IN=1,NNEIGH
                IAT2=NEIGH(IN,IAT)
                IF(IAT2.EQ.0) CYCLE
                DR(:)=R(:,IAT2)-R(:,IAT)
                DISMAX=1.D+10
                DO I1=-1,1
                  DO I2=-1,1
                    DO I3=-1,1
                      T(:)=RBAS(:,1)*DBLE(I1)+RBAS(:,2)*DBLE(I2)+RBAS(:,3)*DBLE(I3)
                      DIS=SQRT((DR(1)+T(1))**2+(DR(2)+T(2))**2+(DR(3)+T(3))**2)
                      IF(DIS.LT.DISMAX) THEN
                        ITNEIGH(1,IN,IAT)=I1
                        ITNEIGH(2,IN,IAT)=I2
                        ITNEIGH(3,IN,IAT)=I3
                        DISMAX=DIS
                      END IF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          END IF
        END IF
      END IF
      RETURN
      END

!
!     ..................................................................
      SUBROUTINE STRCIO_WRITECSSR(NFIL)
      USE STRCIO_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: RBASINV(3,3) ! RBAS**(-1)
      REAL(8)               :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)               :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      REAL(8)               :: DET
      INTEGER(4)            :: IAT
      REAL(8)               :: VEC(3),RBASNEU(3,3)
      REAL(8)               :: ANGSTROM
      INTEGER(4)            :: NEIGH1(8)
!     ******************************************************************
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
!
!     ==================================================================
!     == CONVERT DATA TO UNITS OF LATTICE VECTORS                     ==
!     ==================================================================
      IF(TCRYSTAL) THEN
        RBASINV(1,:)=RBAS(:,2).DYAD.RBAS(:,3)
        RBASINV(2,:)=RBAS(:,3).DYAD.RBAS(:,1)
        RBASINV(3,:)=RBAS(:,1).DYAD.RBAS(:,2)
        DET=DOT_PRODUCT(RBAS(:,1),RBASINV(1,:))
        RBASINV=RBASINV/DET
        A=SQRT(DOT_PRODUCT(RBAS(:,1),RBAS(:,1)))
        B=SQRT(DOT_PRODUCT(RBAS(:,2),RBAS(:,2)))
        C=SQRT(DOT_PRODUCT(RBAS(:,3),RBAS(:,3)))
        GAMMA =ACOS(DOT_PRODUCT(RBAS(:,1),RBAS(:,2))/(A*B))
        ALPHA =ACOS(DOT_PRODUCT(RBAS(:,2),RBAS(:,3))/(B*C))
        BETA  =ACOS(DOT_PRODUCT(RBAS(:,3),RBAS(:,1))/(C*A))
        ALPHA=180.D0/PI*ALPHA
        BETA =180.D0/PI*BETA
        GAMMA=180.D0/PI*GAMMA
        RBASNEU(:,:)=0.D0
        RBASNEU(3,3)=C
        RBASNEU(3,2)=B*COS(ALPHA*PI/180.D0)
        RBASNEU(2,2)=SQRT(B**2-RBASNEU(3,2)**2)
        RBASNEU(3,1)=A*COS(BETA*PI/180.D0)
        RBASNEU(2,1)=(A*B*COS(GAMMA*PI/180.D0)-RBASNEU(3,1)*RBASNEU(3,2))/RBASNEU(2,2)
        RBASNEU(1,1)=SQRT(A**2-RBASNEU(2,1)**2-RBASNEU(3,1)**2)
      END IF
!
!     ==================================================================
!     == WRITE CSSR FILE                                              ==
!     ==================================================================
      IF(TCRYSTAL) THEN
        WRITE(NFIL,FMT='(T39,3F8.3' &
     &   //'/T22,3F8.3,T50,"SPGR = 1 P 1",T72,"OPT = 1"' &
      &  //' /I4,"   1 CREATED BY PAW    "' &
     &   //'/"     0 ",A4,": ",A4)') &
     &   A/ANGSTROM,B/ANGSTROM,C/ANGSTROM,ALPHA,BETA,GAMMA &
     &   ,NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      ELSE
        WRITE(NFIL,FMT='(//I4,"   1 CREATED BY PAW    "' &
     &   //'/"     0 ",A4,": ",A4)')NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      END IF
      DO IAT=1,NAT
        IF(NNEIGH.GT.8) THEN
          NEIGH1(:)=NEIGH(1:8,IAT)
        ELSE
          NEIGH1(:)=0
          NEIGH1(1:NNEIGH)=NEIGH(:,IAT)
        END IF
        IF(NAME(IAT)(2:2).EQ.' ')NAME(IAT)(2:)=NAME(IAT)(3:)
        IF(TCRYSTAL) THEN
          VEC=R(:,IAT)
          VEC=MATMUL(RBASINV,VEC)+100.D0
          VEC=MOD(VEC,1.D0)
          VEC=MATMUL(RBASNEU,VEC)
          WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
     &           IAT,NAME(IAT),VEC(:)/ANGSTROM,NEIGH1,Q(IAT)

        ELSE
           WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
      &           IAT,NAME(IAT),R(:,IAT)/ANGSTROM,NEIGH1,Q(IAT)
         END IF
      ENDDO

      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO_READCSSR(NFIL)
      USE STRCIO_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)               :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      INTEGER(4)            :: IAT
      REAL(8)               :: ANGSTROM
      INTEGER(4)            :: SPCGRPNR
      CHARACTER(11)         :: SPCGRPNM
      INTEGER(4)            :: IORTH,IATNR,ISVAR
      CHARACTER(60)         :: TITLE
      CHARACTER(53)         :: STRING
      CHARACTER(1)          :: CH1
!     ******************************************************************
      NEIGH(:,:)=0
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
!
!     ==================================================================
!     == WRITE CSSR FILE                                              ==
!     ==================================================================
      READ(NFIL,FMT='(38X,3F8.3)')A,B,C
      READ(NFIL,FMT='(21X,3F8.3,4X,6X,I2,1X,A11)') &
     &               ALPHA,BETA,GAMMA,SPCGRPNR,SPCGRPNM
      READ(NFIL,FMT='(2I4,1X,A60)')NAT,IORTH,TITLE
      READ(NFIL,FMT='(A53)')STRING
      OBJECTNAME=STRING
      ALLOCATE(R(3,NAT))
      ALLOCATE(Q(NAT))
      ALLOCATE(NAME(NAT))
      NNEIGH=8
      ALLOCATE(NEIGH(NNEIGH,NAT))
      DO IAT=1,NAT
        READ(NFIL,FMT='(I4,1X,A4,2X,3(F9.5,1X),8I4,1X,F7.3)') &
     &           IATNR,NAME(IAT),R(:,IAT),NEIGH(:,IAT),Q(IAT)
        CH1=+NAME(IAT)(2:2)
        ISVAR=IACHAR(CH1)
        IF(ISVAR.LT.65.OR.ISVAR.GT.90) THEN
          NAME(IAT)(3:)=NAME(IAT)(2:)
          NAME(IAT)(2:2)=' '
        END IF
      ENDDO
      RBAS(:,:)=0.D0
      RBAS(3,3)=C
      RBAS(3,2)=B*COS(ALPHA*PI/180.D0)
      RBAS(2,2)=SQRT(B**2-RBAS(3,2)**2)
      RBAS(3,1)=A*COS(BETA*PI/180.D0)
      RBAS(2,1)=(A*B*COS(GAMMA*PI/180.D0)-RBAS(3,1)*RBAS(3,2))/RBAS(2,2)
      RBAS(1,1)=SQRT(A**2-RBAS(2,1)**2-RBAS(3,1)**2)
      TCRYSTAL=(A*B*C.NE.0) 
      IF(IORTH.EQ.0) THEN
        DO IAT=1,NAT
          R(:,IAT)=MATMUL(RBAS,R(:,IAT))
        ENDDO
      END IF
      RBAS(:,:)=RBAS(:,:)*ANGSTROM
      R(:,:)=R(:,:)*ANGSTROM
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STRCIO_WRITEXYZ(NFIL,FRAME,NAT,ID,R)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)  ,INTENT(IN) :: FRAME
      INTEGER(4)  ,INTENT(IN) :: NAT
      CHARACTER(2),INTENT(IN) :: ID(NAT)
      REAL(8)     ,INTENT(IN) :: R(3,NAT)
      INTEGER(4)              :: IAT
      REAL(8)                 :: ANGSTROM
!     ******************************************************************
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      IF(FRAME.EQ.1) REWIND NFIL
      WRITE(NFIL,*)NAT
      WRITE(NFIL,FMT='(A10,I10)')'NONAME',FRAME
      DO IAT=1,NAT
        WRITE(NFIL,FMT='(A2,2X,3(F10.5,1X))')ID(IAT),R(:,IAT)/ANGSTROM
      ENDDO
      RETURN
      END
