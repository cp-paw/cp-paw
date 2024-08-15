       PROGRAM MAIN
!      *************************************************************************
!      ** READS A STRUCTURE FILE IN THE POSCAR FORMAT OF VASP AND WRITES      **
!      ** THE STRUCTURE IN A FORMAT THAT CAN BE INSERTED INTO CP-PAWS         **
!      ** STRUCTURE FILE                                                      **
!      **                                                                     **
!      ** USAGE:                                                              **
!      **   SUPPLY NAME OF THE POSCAR FILE AS ARGUMENT                        **
!      **   SUPPLY OPTIONALLY THE LENGTH UNIT IN ANGSTROM OF THE OUTPUT FILE  **
!      **   RESULT IS WRITTEN TO STANDARD OUT SO THAT IT CAN BE PIPED INTO    **
!      **   THE STRUCTURE FILE                                                **
!      **                                                                     **
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)   ,PARAMETER     :: ANGSTROM=1.8897259886D0
       INTEGER(4),PARAMETER     :: NATX=500       
       REAL(8)                  :: LUNIT      ! LENGTH UNIT FOR RBAS AND R
       REAL(8)                  :: RBAS(3,3)  ! UNIT CELL
       CHARACTER(20)            :: NAME(NATX)   ! ATOM NAME
       REAL(8)                  :: R(3,NATX)  ! ATOM POSITION
       CHARACTER(:),ALLOCATABLE :: POSCARFILE
       CHARACTER(:),ALLOCATABLE :: LUNITSTRING
       INTEGER(4)               :: NARG
       INTEGER(4)               :: ARGLEN
       INTEGER(4)               :: STAT
       INTEGER(4)               :: NFIL=20   ! FILE UNIT FOR POSCAR INPUT FILE
       INTEGER(4)               :: NFILOUT=6 ! STANDARD OUT
       REAL(8)     ,ALLOCATABLE :: R0(:,:)
       LOGICAL(4)               :: TCARTESIAN
       INTEGER(4)               :: ISVAR,IAT,ISP,J,NAT
       CHARACTER(256)           :: COMMENT
       CHARACTER(2)             :: MINUSH  ! LOWERCASE -H
!      *************************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!      =========================================================================
!      == RESOLVE  ARGUMENT LIST                                              ==
!      =========================================================================
       NARG=COMMAND_ARGUMENT_COUNT()
       IF(NARG.EQ.1.OR.NARG.EQ.2) THEN
!        == COLLECT FILE NAME FROM ARGUMENT LIST ===============================
         CALL GET_COMMAND_ARGUMENT(NUMBER=1,LENGTH=ARGLEN)
         ALLOCATE(CHARACTER(ARGLEN):: POSCARFILE)
         CALL GET_COMMAND_ARGUMENT(NUMBER=1,VALUE=POSCARFILE,STATUS=STAT)
         IF(STAT.NE.0) THEN
           WRITE(*,*)'ERROR IN PAW_FROMPOSCAR TOOL'
           WRITE(*,*)'ERROR ',STAT,' RETRIEVING ARGUMENT 1'
           STOP 'MAIN OF PAW_FROMPOSCAR'
         END IF
!
!        == CHECK FOR HELP REQUEST AND PRINT USAGE =============================
         MINUSH='-H'
         MINUSH(2:2)=ACHAR(104)
         IF(POSCARFILE.EQ.MINUSH) THEN
           CALL PRUSAGE()
           STOP 'MAIN OF PAW_FROMPOSCAR AFTER PRINTING USAGE'
         END IF
!        
!        == READ OPTIONAL ARGUMENT "LENGTH UNIT IN ANGSTROM" ===================
         IF(NARG.EQ.2) THEN
           CALL GET_COMMAND_ARGUMENT(NUMBER=2,LENGTH=ARGLEN)
           ALLOCATE(CHARACTER(ARGLEN):: LUNITSTRING)
         CALL GET_COMMAND_ARGUMENT(NUMBER=2,VALUE=lunitstring,STATUS=STAT)
         IF(STAT.NE.0) THEN
           WRITE(*,*)'ERROR IN PAW_FROMPOSCAR TOOL'
           WRITE(*,*)'ERROR ',STAT,' RETRIEVING ARGUMENT 2'
           STOP 'MAIN OF PAW_FROMPOSCAR'
         END IF
           READ(LUNITSTRING,*)LUNIT
           LUNIT=LUNIT*ANGSTROM  !PROCEED WITH HARTREE ATOMIC UNITS
         ELSE
           LUNIT=1.D0*ANGSTROM   !PROCEED WITH HARTREE ATOMIC UNITS
         END IF
       ELSE
         WRITE(*,*)'ERROR IN PAW_FROMPOSCAR TOOL'
         WRITE(*,*)'NO OR MORE THAN TWO ARGUMENTS SUPPLIED'
         CALL PRUSAGE()
         STOP 'MAIN OF PAW_FROMPOSCAR'
       END IF
!
!      =========================================================================
!      ==  READ STRUCTURE FORM POSCAR FILE                                    ==
!      ==  STRUCTURE IS RETURNED CARTESIAN AND IN HARTREE ATOMIC UNITS        ==
!      =========================================================================
       OPEN(NFIL,FILE=POSCARFILE)
       CALL READPOSCAR(NFIL,COMMENT,RBAS,NATX,NAT,NAME,R)
 
!      =========================================================================
!      =========================================================================
!      =========================================================================
       ALLOCATE(R0(3,NAT))
       R0(:,:)=R(:,:NAT)
       CALL WRITECPPAW(NFILOUT,COMMENT,NAT,LUNIT,NAME,RBAS,R0)
!
       CALL ERROR$NORMALSTOP()
       STOP
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE PRUSAGE()
!      *************************************************************************
!      ** WRITES USAGE INFORMATION OF THIS TOOL TO STANDARD OUT               **
!      *************************************************************************
       WRITE(*,*)
       WRITE(*,*)'USAGE:'
       WRITE(*,*)
       WRITE(*,*)'       PAW_FROMPOSCAR.X IN.POSCAR [LUNIT] > CP-PAWFILE'
       WRITE(*,*)
       WRITE(*,*)'IN.POSCAR IS THE NAME OF THE POSCAR FILE'
       WRITE(*,*)'CP-PAWFILE IS THE NAME OF THE FILE IN CP-PAW STYLE'
       WRITE(*,*)'LUNIT(OPTIONAL): LENGTH UNIT IN ANGSTROM FOR OUTPUT'
       WRITE(*,*)
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE WRITECPPAW(NFIL,COMMENT,NAT,LUNIT,NAME,RBAS,R)
!      *************************************************************************
!      **  WRITES ATOMIC STRUCTURE IN THE FORMAT OF THE CP-PAW CODE           **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)              :: NFIL
       CHARACTER(*),INTENT(IN) :: COMMENT
       INTEGER(4)  ,INTENT(IN) :: NAT
       REAL(8)     ,INTENT(IN) :: LUNIT ! LENGTH UNIT IN ANGSTROM
       REAL(8)     ,INTENT(IN) :: RBAS(3,3)
       CHARACTER(*),INTENT(IN) :: NAME(NAT)
       REAL(8)     ,INTENT(IN) :: R(3,NAT)
       REAL(8)     ,PARAMETER  :: ANGSTROM=1.8897259886D0
       INTEGER(4)              :: I
!      *************************************************************************
       WRITE(NFIL,FMT='(80("<"))')
       WRITE(NFIL,FMT='(T3,"!GENERIC LUNIT[AA]= ",F10.5," !END")')LUNIT/ANGSTROM
      
       WRITE(NFIL,FMT='("# ",A)')COMMENT
       WRITE(NFIL,FMT='(T3,"!LATTICE",T30,"T= ",3F10.5)')RBAS(:,1)/LUNIT
       WRITE(NFIL,FMT='(T33,3F10.5)')RBAS(:,2)/LUNIT
       WRITE(NFIL,FMT='(T33,3F10.5," !END ")')RBAS(:,3)/LUNIT
       DO I=1,NAT
         WRITE(NFIL,FMT='(T3,"!ATOM NAME= ",A,T30," R=",3F10.5," !END")') &
      &                           "'"//TRIM(ADJUSTL(NAME(I)))//"'",R(:,I)/LUNIT
       ENDDO
       WRITE(NFIL,FMT='(80(">"))')
       RETURN
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
       CHARACTER(*),INTENT(OUT) :: ID(NATX)     ! ATOM NAME
       REAL(8)     ,INTENT(OUT) :: R(3,NATX)
       REAL(8)     ,PARAMETER   :: ANGSTROM=1.8897259886D0
       REAL(8)                  :: LUNIT
       INTEGER(4)               :: NSP
       INTEGER(4)  ,ALLOCATABLE :: NOFISP(:)
       CHARACTER(2),ALLOCATABLE :: SY(:)
       CHARACTER(256)           :: LINE,LINE1
       LOGICAL(4)               :: TCARTESIAN
       INTEGER(4)               :: ISVAR,IAT,ISP,J
       REAL(8)                  :: SVAR
!      *************************************************************************
       R(:,:)=0.D0
       ID(:)=' '
!
!      == LINE 1: NAME =========================================================
       READ(NFIL,FMT='(A)')COMMENT
!
!      == LINE 2: LENGTH UNIT OR VOLUME (IF NEGATIVE)===========================
!      == PER DEFAULT A POSCAR FILE IS IN ANGSTROM ============================
       READ(NFIL,*)LUNIT
!
!      == LINES 3-5: LATTICE VECTORS IN UNITS OF LUNIT =========================
       READ(NFIL,*)RBAS(:,1)
       READ(NFIL,*)RBAS(:,2)
       READ(NFIL,*)RBAS(:,3)
       IF(LUNIT.LT.0.D0) THEN 
         SVAR=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
      &      +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
      &      +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3))
         LUNIT=(-LUNIT/SVAR)**(1.D0/3.D0)
       END IF
       LUNIT=LUNIT*ANGSTROM
!
!      == LINE 6: ELEMENT SYMBOLS OF ATOM TYPES ================================
       READ(NFIL,FMT='(A)')LINE
       LINE1=LINE
       NSP=0
       DO      ! COUNT ATOM TYPES BEFORE READING
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
!      == LINE 7: NUMBER OF ATOMS PER ATOM TYPE ================================
       READ(NFIL,*)NOFISP
       NAT=SUM(NOFISP)
       IF(NAT.GT.NATX) THEN
         STOP '#(ATOMS) EXCEEDS MAX'
       END IF
!
!      == LINE 8: RESOLVE WHETHER COORDINATES ARE FRACTIONAL OR CARTESIAN ======
       READ(NFIL,FMT='(A)')LINE
       LINE=ADJUSTL(LINE)
!      __THE ONLY KEY CHARACTERS RECOGNIZED BY VASP ARE UPPER- AND LOWERCASE____
!      __'C' OR 'K' FOR SWITCHING TO THE CARTESIAN MODE AND_____________________
!      __'D' OR 'F' FOR FRACTIONAL COORDINATES__________________________________
!      __67=UPPERCASE C, 75=UPPERCASE K, 99=LOWERCASE C, 107=LOWERCASE K________
!      __68=UPPERCASE D, 70=UPPERCASE F, 100=LOWERCASE D,102=LOWERCASE F________
       ISVAR=IACHAR(LINE(1:1))
       TCARTESIAN=(ISVAR.EQ.67.OR.ISVAR.EQ.75.OR.ISVAR.EQ.99.OR.ISVAR.EQ.107)
       IF(.NOT.TCARTESIAN.AND. &
      &   .NOT.(ISVAR.EQ.100.OR.ISVAR.EQ.102.OR.ISVAR.EQ.68 &
      &                                     .OR.ISVAR.EQ.70)) THEN
         STOP 'ERROR: NEITHER FRACTIONAL NOR CARTESIAN'
       END IF
!
!      == LINES 9+ : ATOMIC POSITIONS ==========================================
       IAT=0
       DO ISP=1,NSP
         DO J=1,NOFISP(ISP)
           IAT=IAT+1
           WRITE(LINE,*)J
           ID(IAT)=SY(ISP)//TRIM(ADJUSTL(LINE))
           READ(NFIL,*)R(:,IAT)
           IF(.NOT.TCARTESIAN) R(:,IAT)=MATMUL(RBAS,R(:,IAT))
         ENDDO
       ENDDO
       RBAS=RBAS*LUNIT
       R=R*LUNIT
       RETURN 
       END
