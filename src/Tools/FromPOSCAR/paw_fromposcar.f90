       PROGRAM MAIN
!      *************************************************************************
!      ** READS A STRUCTURE FILE IN THE POSCAR FORMAT OF VASP AND WRITES      **
!      ** THE STRUCTURE IN A FORMAT THAT CAN BE INSERTED INTO CP-PAW'S        **
!      ** STRUCTURE FILE                                                      **
!      **                                                                     **
!      ** USAGE:                                                              **
!      **   SUPPLY NAME OF THE POSCAR FILE AS ARGUMENT                        **
!      **   RESULT IS WRITTEN TO STANDARD OUT SO THAT IT CAN BE PIPED INTO    **
!      **   THE STRUCTURE FILE                                                **
!      **                                                                     **
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)     ,PARAMETER   :: ANGSTROM=1.8897259886D0
       REAL(8)                  :: LUNIT
       REAL(8)                  :: RBAS(3,3)
       INTEGER(4)               :: NSP
       INTEGER(4)  ,ALLOCATABLE :: NATOFISP(:)
       CHARACTER(2),ALLOCATABLE :: SYOFISP(:)
       CHARACTER(2),ALLOCATABLE :: SY(:)
       REAL(8)     ,ALLOCATABLE :: R0(:,:)
       CHARACTER(:),ALLOCATABLE :: POSCARFILE
       INTEGER(4)               :: NARG
       INTEGER(4)               :: ARGLEN
       INTEGER(4)               :: STAT
       CHARACTER(256)           :: FIRSTLINE
       CHARACTER(256)           :: LINE,LINE1
       INTEGER(4)               :: NFIL=20
       LOGICAL(4)               :: TCARTESIAN
       INTEGER(4)               :: ISVAR,IAT,ISP,J,NAT
       CHARACTER(256)           :: COMMENT
       INTEGER(4),PARAMETER     :: NATX=500       
       REAL(8)                  :: R(3,NATX)
       CHARACTER(20)            :: ID(NATX)
!      *************************************************************************
!
!      =========================================================================
!      == RESOLVE  ARGUMENT LIST                                              ==
!      =========================================================================
       NARG=COMMAND_ARGUMENT_COUNT()
       IF(NARG.EQ.0) THEN
         WRITE(*,*)'ERROR IN PAW_FROMPOSCAR TOOL'
         WRITE(*,*)'NO ARGUMENTS SUPPLIED'
         WRITE(*,*)'SUPPLY NAME OF THE POSCAR FILE'
         STOP 'MAIN OF PAW_FROMPOSCAR'
       ELSE 
         CALL GET_COMMAND_ARGUMENT(NUMBER=1,LENGTH=ARGLEN)
         ALLOCATE(CHARACTER(ARGLEN):: POSCARFILE)
         CALL GET_COMMAND_ARGUMENT(NUMBER=1,VALUE=POSCARFILE,STATUS=STAT)
         IF(STAT.NE.0) THEN
           WRITE(*,*)'ERROR IN PAW_FROMPOSCAR TOOL'
           WRITE(*,*)'ERROR ',STAT,' RETRIEVING ARGUMENT'
           STOP 'MAIN OF PAW_FROMPOSCAR'
         END IF
       END IF
!
!      =========================================================================
!      ==  READ STRUCTURE FORM POSCAR FILE                                    ==
!      =========================================================================
       OPEN(NFIL,FILE=POSCARFILE)
       CALL READPOSCAR(NFIL,COMMENT,RBAS,NATX,NAT,ID,R)
       ALLOCATE(R0(3,NAT))
       R0(:,:)=R(:,:NAT)
!
!      =========================================================================
!      =========================================================================
!      =========================================================================
       WRITE(*,FMT='(T3,"!GENERIC LUNIT[AA]=1.0 !END")') 
       WRITE(*,FMT='(T3,"!LATTICE T=",3F15.9)') RBAS(:,1)/ANGSTROM
       WRITE(*,FMT='(T3,"           ",3F15.9)') RBAS(:,2)/ANGSTROM
       WRITE(*,FMT='(T3,"           ",3F15.9," !END")') RBAS(:,3)/ANGSTROM
       DO IAT=1,NAT
         WRITE(*,FMT='(T3,"!ATOM NAME=",A," R=",3F15.9," !END")') &
      &                 "'"//TRIM(ID(IAT))//"'",R0(:,IAT)/ANGSTROM
       ENDDO
       STOP
       END


!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE READPOSCAR(NFIL,COMMENT,RBAS,NATX,NAT,ID,R)
!      *************************************************************************
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN)  :: NFIL
       INTEGER(4)  ,INTENT(IN)  :: NATX
       INTEGER(4)  ,INTENT(OUT) :: NAT
       CHARACTER(*),INTENT(OUT) :: COMMENT
       REAL(8)     ,INTENT(OUT) :: RBAS(3,3)
       CHARACTER(*),INTENT(OUT) :: ID(NATX)
       REAL(8)     ,INTENT(OUT) :: R(3,NATX)
       REAL(8)     ,PARAMETER   :: ANGSTROM=1.8897259886D0
       REAL(8)                  :: LUNIT
       INTEGER(4)               :: NSP
       INTEGER(4)  ,ALLOCATABLE :: NOFISP(:)
       CHARACTER(2),ALLOCATABLE :: SY(:)
       CHARACTER(256)           :: LINE,LINE1
       LOGICAL(4)               :: TCARTESIAN
       INTEGER(4)               :: ISVAR,IAT,ISP,J
!      *************************************************************************
       R(:,:)=0.D0
       ID(:)=' '
       READ(NFIL,FMT='(A)')COMMENT
!      == PER DEFAULTS A POSCAR FILE IS IN ANGSTROM ============================
       READ(NFIL,*)LUNIT
       LUNIT=LUNIT*ANGSTROM

       READ(NFIL,*)RBAS(:,1)
       READ(NFIL,*)RBAS(:,2)
       READ(NFIL,*)RBAS(:,3)
!
!      == READ ELEMENT SYMBOLS =================================================
       READ(NFIL,FMT='(A)')LINE
       LINE1=LINE
       NSP=0
       DO
         LINE1=ADJUSTL(LINE1)
         IF(LINE1(1:2).EQ.'  ') EXIT
         NSP=NSP+1
         LINE1=ADJUSTL(LINE1(3:))
       ENDDO
       ALLOCATE(SY(NSP))
       ALLOCATE(NOFISP(NSP))
       READ(LINE,*)SY
       DO ISP=1,NSP
         IF(SY(ISP)(2:2).EQ.' ')SY(ISP)(2:2)='_'
       ENDDO
!
!      == READ #(ATOMS) FOR EACH ELEMENT-SYMBOL ================================
       READ(NFIL,*)NOFISP
       NAT=SUM(NOFISP)
       IF(NAT.GT.NATX) THEN
         STOP '#(ATOMS) EXCEEDS MAX'
       END IF
!      == RESOLVE WHETHER COORDINATES ARE FRACTIONAL OR CARTESIAN ==============
       READ(NFIL,FMT='(A)')LINE
       LINE=ADJUSTL(LINE)
!      == THE ONLY KEY CHARACTERS RECOGNIZED BY VASP ARE 'C', 'C', 'K' OR 'K' 
!      == FOR SWITCHING TO THE CARTESIAN MODE.
       ISVAR=IACHAR(LINE(1:1))
       TCARTESIAN=(ISVAR.EQ.67.OR.ISVAR.EQ.75.OR.ISVAR.EQ.99.OR.ISVAR.EQ.107)
       IF(.NOT.TCARTESIAN.AND. &
      &   .NOT.(ISVAR.EQ.100.OR.ISVAR.EQ.102.OR.ISVAR.EQ.68 &
      &                                     .OR.ISVAR.EQ.70)) THEN
         STOP 'ERROR: NEITHER FRACTIONAL NOR CARTESIAN'
       END IF
!
!      == READ COORDINATES =====================================================
       IAT=0
       DO ISP=1,NSP
         DO J=1,NOFISP(ISP)
           IAT=IAT+1
           WRITE(LINE,*)IAT
           ID(IAT)=SY(ISP)//TRIM(ADJUSTL(LINE))
           READ(NFIL,*)R(:,IAT)
           IF(.NOT.TCARTESIAN) R(:,IAT)=MATMUL(RBAS,R(:,IAT))
         ENDDO
       ENDDO
       RBAS=RBAS*LUNIT
       R=R*LUNIT
       RETURN 
       END
