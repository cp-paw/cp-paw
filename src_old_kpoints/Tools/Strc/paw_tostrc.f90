!     ..................................................................
      PROGRAM PAW_STRC
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT none
      TYPE(LL_TYPE)               :: LL_STRC
      INTEGER(4)                  :: NFIL
      INteGER(4)                  :: NFILSTRC,nfilpdb
      CHARACTER(128)              :: ROOTNAME     ! COMMON ROOT OF THE FILENAMES
      CHARACTER(128)              :: OBJECTNAME
      LOGICAL                     :: TCHK
      CHARACTER(32)               :: STRING
      logical(4)                  :: tcrystal
      INTEGER(4)                  :: i
      CALL LINKEDLIST$NEW(LL_STRC)
!
!     ==================================================================
!     == GET FILE NAME ROOT FROM THE ARGUMENT LIST AND CONSTRUCT      ==
!     == FILE NAMES                                                   ==
!     ==================================================================
      CALL GETARG(1,ROOTNAME)
      IF(ROOTNAME(1:1).EQ.'-') THEN
        TCRYSTAL=(+ROOTNAME(2:2).EQ.+'C')
        CALL GETARG(2,ROOTNAME)
      END IF
      IF(LEN(TRIM(ROOTNAME)).EQ.0) THEN
        STOP 'NO ROOTNAME SUPPLIED'
      END IF
      I=INDEX(ROOTNAME,'/',BACK=.TRUE.)
      OBJECTNAME=ROOTNAME(I+1:)
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.PRESTRC')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('CSSR',.TRUE.,-'.CSSR')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','FORM','FORMATTED')
!
!     ==================================================================
!     == READ CSSR FILE                                               ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('CSSR',NFIL)
      CALL FILEHANDLER$UNIT('STRC',NFILSTRC)
      CALL readCSSR(NFIL,NFILSTRC)
      WRITE(*,FMT='("======= TASK FINISHED ========")')
      STOP
      END
!     ..................................................................
      SUBROUTINE READCSSR(NFIL,NFILSTRC)
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE BOND_TYPE
        INTEGER(4) :: ATOM1
        INTEGER(4) :: ATOM2
        REAL(8)    :: BO
      END TYPE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      integer(4)  ,INTENT(IN) :: NFILSTRC
      INTEGER(4)              :: NAT
      REAL(8)                 :: RBAS(3,3)
      CHARACTER(8),ALLOCATABLE:: NAME(:)
      REAL(8)     ,ALLOCATABLE:: R(:,:)
      REAL(8)     ,ALLOCATABLE:: Q(:)
      INTEGER(4)  ,ALLOCATABLE:: NN(:,:)
      INTEGER(4)              :: NBOND
      TYPE(BOND_TYPE),ALLOCATABLE :: BOND(:)
      REAL(8)                 :: PI
      REAL(8)                 :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)                 :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      INTEGER(4)              :: IAT,I,J,IBOND,isvar,iat1,iat2
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
!
!     ==================================================================
!     == READ HEADER                                                  ==
!     ==================================================================
      READ(NFIL,FMT='(T39,3F8.3/T22,3F8.3)')A,B,C,ALPHA,BETA,GAMMA
      READ(NFIL,FMT='(I4)')NAT
      READ(NFIL,*)
!
!     ==================================================================
!     == CALCULATE LATTICE VECTORS                                    ==
!     ==================================================================
      RBAS(:,:)=0.D0
      RBAS(3,3)=C
      RBAS(3,2)=B*COS(ALPHA*PI/180.D0)
      RBAS(2,2)=SQRT(B**2-RBAS(3,2)**2)
      RBAS(3,1)=A*COS(BETA*PI/180.D0)
      RBAS(2,1)=(A*B*COS(GAMMA*PI/180.D0)-RBAS(3,1)*RBAS(3,2))/RBAS(2,2)
      RBAS(1,1)=SQRT(A**2-RBAS(2,1)**2-RBAS(3,1)**2)
!
!     ==================================================================
!     == READ POSITIONS                                              ==
!     ==================================================================
      ALLOCATE(R(3,NAT))
      ALLOCATE(NAME(NAT))
      ALLOCATE(Q(NAT))
      ALLOCATE(NN(8,NAT))
      DO IAT=1,NAT
        READ(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
     &           ISVAR,NAME(IAT),R(:,IAT),NN(:,IAT),Q(IAT)
       i=iachar(name(iat)(2:2))
       if(i.ge.48.and.i.le.57) then
         name(iat)=name(iat)(1:1)//'_'//name(iat)(2:)
       end if
      ENDDO
!
!     ==================================================================
!     == REMOVE BONDS THAT ARE COUNTED TWICE                          ==
!     ==================================================================
      NBOND=0
      DO IAT=1,NAT
        DO I=1,8
          IF(NN(I,IAT).LT.IAT) NN(I,IAT)=0
          IF(NN(i,IAT).NE.0) NBOND=NBOND+1
        ENDDO
      ENDDO
      ALLOCATE(BOND(NBOND))
      IBOND=0
      DO IAT=1,NAT
        DO I=1,8
          IF(NN(i,IAT).EQ.0) CYCLE
          IBOND=IBOND+1
          BOND(IBOND)%ATOM1=IAT
          BOND(IBOND)%ATOM2=NN(I,IAT)
          BOND(IBOND)%BO=0.D0
          BOND(IBOND)%BO=1.D0
          NN(I,IAT)=0
          DO J=I+1,8
            IF(NN(J,IAT).EQ.BOND(IBOND)%ATOM2) THEN
              BOND(IBOND)%BO=BOND(IBOND)%BO+1.D0
              NBOND=NBOND-1
              NN(J,IAT)=0
            END IF              
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == WRITE CSSR FILE                                              ==
!     ==================================================================
      DO IAT=1,NAT
        WRITE(NFILSTRC,FMT='(A,3F10.5,A,F10.5,A)')"!ATOM NAME='"//TRIM(NAME(IAT))//"' R=",R(:,IAT) &
     &                   ," Q=",Q(IAT)," FFTYPE='"//name(iat)(1:2)//"' !END"
      ENDDO
      DO IBOND=1,NBOND
        IAT1=BOND(IBOND)%ATOM1
        IAT2=BOND(IBOND)%ATOM2
        WRITE(NFILSTRC,fmt='(a,f5.2,a)')"!BOND ATOM1='"//TRIM(NAME(IAT1)) &
     &                  //"' ATOM2='"//TRIM(NAME(IAT2)) &
     &                  //"' BO=",BOND(IBOND)%BO," !END"
      ENDDO        
      RETURN
      END
