!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)               :: LL_STRC
      INTEGER(4)                  :: NFIL
      LOGICAL                     :: TCRYSTAL=.FALSE.
      REAL(8)                     :: ANGSTROM
      CHARACTER(256)              :: ROOTNAME  ! COMMON ROOT OF THE FILENAMES
      CHARACTER(256)              :: OBJECTNAME
      LOGICAL                     :: TCHK
      REAL(8)                     :: RUNIT     ! LENGTH UNIT ON STRUCTURE FILE
      REAL(8)                     :: RBAS(3,3) ! LATTICE VECTORS
      CHARACTER(32),ALLOCATABLE   :: NAME(:)   ! ATOM NAMES
      INTEGER(4)                  :: NAT       ! #(ATOMS IN THE QM PART)
      REAL(8),      ALLOCATABLE   :: R(:,:)    ! ATOMIC POSITIONS
      REAL(8),      ALLOCATABLE   :: Q(:)      ! CHARGES
      INTEGER(4)                  :: NATMM     ! #(ATOMS OF THE MM PART)
      CHARACTER(32),ALLOCATABLE   :: MMNAME(:) ! ATOM NAMES
      REAL(8),      ALLOCATABLE   :: MMR(:,:)  ! POSITIONS FOR MM PART OF QM-MM
      REAL(8),      ALLOCATABLE   :: MMQ(:)    ! ATOMIC CHARGES MM PART OF QM-MM
      LOGICAL(4)                  :: TQMMM
      INTEGER(4)                  :: IAT,IAT1
      CHARACTER(256)              :: STRING
      INTEGER(4)                  :: I,ISVAR
      INTEGER(4)                  :: NARGS
      INTEGER(4)                  :: NFILO
      INTEGER(4)                  :: NDUP(3)
      LOGICAL                     :: TINPUT
      LOGICAL                     :: THELP
      LOGICAL                     :: TCM
!     **************************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
      NDUP(:)=1
      TINPUT=.FALSE.
      CALL LINKEDLIST$NEW(LL_STRC)
!
!     ==========================================================================
!     == GET FILE NAME ROOT FROM THE ARGUMENT LIST AND CONSTRUCT              ==
!     == FILE NAMES                                                           ==
!     ==========================================================================
      NARGS=COMMAND_ARGUMENT_COUNT()
      
!     ==  DETECT HELP REQUESTS =================================================
      THELP=.FALSE.
      TCM=.FALSE.     ! CRYSTAL OR MOLECULE SPECIFIED
      DO I=1,NARGS
        CALL GET_COMMAND_ARGUMENT(I,STRING)
        STRING=+STRING
        IF(STRING(1:2).EQ.'-H'.OR.STRING.EQ.'?') THEN
          THELP=.TRUE.
        END IF
        IF(STRING(1:2).EQ.'-C'.OR.STRING(1:2).EQ.'-M') THEN
          IF(TCM) THEN
            CALL ERROR$MSG('ONE ON OF THE OPTIONS -C AND -M IS ALLOWED')
            CALL ERROR$STOP('MAIN')
          END IF
          TCM=.TRUE.
        END IF
        TCHK=STRING(1:2).EQ.'?'.OR.STRING(1:2).EQ.'-H'.OR.STRING(1:2).EQ.'-I' &
     &                         .OR.STRING(1:2).EQ.'-C'.OR.STRING(1:2).EQ.'-M'
        IF(.NOT.TCHK.AND.I.LT.NARGS) THEN
          WRITE(*,'("ILLEGAL ARGUMENT ",A)')TRIM(STRING)
          WRITE(*,*)
          THELP=.TRUE.
        END IF
      ENDDO
      THELP=THELP.OR.(.NOT.TCM)
      IF(.NOT.TCM) THEN
        WRITE(*,'("ONE OF THE OPTIONS -M, -C OR -CIJK MUST BE SPECIFIED")')
        WRITE(*,*)
      END IF
!
!     == WRITE INFO  ===========================================================
      IF(THELP) THEN
        WRITE(*,'(A)')-"CALLING SEQUENCE: PAW_STRC.X ARGS ROOTNAME"
        WRITE(*,'(A)')-"ROOTNAME IS THE ROOT NAME OF THE STRC_OUT FILE"
        WRITE(*,'(A)')-"ARGUMENTS CAN BE:"
        WRITE(*,'(T2,A,T10,A)')-'?','PRINT HELP MESSAGE'
        WRITE(*,'(T2,A,T10,A)')-'-H','PRINT HELP MESSAGE'
        WRITE(*,'(T2,A,T10,A)') &
    &         -'-I','USE INPUT STRUCTURE FILE INSTEAD OF STRC_OUT'
        WRITE(*,'(T2,A,T10,A)')-'-M','CONSIDER AS MOLECULE'
        WRITE(*,'(T2,A,T10,A)')-'-C','CONSIDER AS CRYSTAL'
        WRITE(*,'(T2,A,T10,A)')-'-CIJK','CONSIDER AS CRYSTAL'
        WRITE(*,'(T10,A)') &
    &         -'AND MULTIPLY UNIT CELL BY FACTORS I,J,K' &
    &       //-' ALONG THE THREE LATTICE VECTORS'
        WRITE(*,'(T10,A)')-'I,J,K ARE POSITIVE SINGLE-DIGIT INTEGERS'
        WRITE(*,'("OUTPUT:")')
        WRITE(*,'(T4,"ROOTNAME",A,T20,"PROTOCOL FILE ")')-'.SPROT'
        WRITE(*,'(T4,"ROOTNAME",A,T20,A)')-'.CML' &
    &                                ,"CRYSTAL STRUCTURE FILE IN THE CML FORMAT"
        WRITE(*,'(T4,"ROOTNAME",A,T20,A)')-'.XYZ' &
    &                                  ,-"MOLECULAR OUTPUT FOR VIEWING PURPOSES"
        WRITE(*,'(T4,"ROOTNAME",A,T20,A)')-'.CSSR' &
    &                                  ,-"MOLECULAR OUTPUT FOR VIEWING PURPOSES"
        WRITE(*,'("REMARKS:")')
        WRITE(*,'(T4,A)')-"REQUIRES ATOMNAMES TO START WITH THE ELEMENT SYMBOL"
        STOP
      END IF
!
!     == RESOLVE ARGUMENTS =====================================================
      TCRYSTAL=.FALSE.
      DO I=1,NARGS-1
        CALL GET_COMMAND_ARGUMENT(I,STRING)
        WRITE(*,FMT='("ARGUMENT(",I2,"):",A)')I,TRIM(STRING)
        IF(STRING(1:1).NE.'-') THEN
          CALL ERROR$MSG('ARGUMENT DOES NOT HAVE PRECEDING "-"')
          CALL ERROR$I4VAL('ARGUMENT NR:',I)
          CALL ERROR$CHVAL('ARGUMENT',TRIM(STRING))
          CALL ERROR$STOP('MAIN')
        END IF
        STRING(1:)=+STRING(2:) ! REMOVE PRECEEDING DASH AND MAKE UPPERCASE
        IF(STRING(1:1).EQ.'C') THEN
          TCRYSTAL=.TRUE.
          IF(STRING(2:2).NE.' ') THEN
            STRING(1:)=STRING(2:)
            READ(STRING,FMT=*)ISVAR
            IF(ISVAR.GT.999.OR.ISVAR.LT.111) THEN
              CALL ERROR$MSG('UNIT CELL MULTIPLICATION ARGUMENT OUT OF RANGE')
              CALL ERROR$MSG('ARGUMENT MUST BE IJK WHERE I,J,K ARE SINGLE-DIGIT INTEGERS')
              CALL ERROR$CHVAL('ARGUMENT',TRIM(STRING))
              CALL ERROR$STOP('MAIN')
            END IF
            NDUP(1)=INT(ISVAR/100)
            ISVAR=ISVAR-100*NDUP(1)
            NDUP(2)=INT(ISVAR/10)
            ISVAR=ISVAR-10*NDUP(2)
            NDUP(3)=ISVAR
            IF(NDUP(1)*NDUP(2)*NDUP(3).EQ.0) THEN
              CALL ERROR$MSG('LATTICE DISPLACEMENT MUST BE GREATER THAN ZERO')
              CALL ERROR$I4VAL('NDUP(1)',NDUP(1))
              CALL ERROR$I4VAL('NDUP(2)',NDUP(2))
              CALL ERROR$I4VAL('NDUP(3)',NDUP(3))
              CALL ERROR$CHVAL('ARGUMENT',TRIM(STRING))
              CALL ERROR$STOP('MAIN')
            END IF
          END IF
        ELSE IF(STRING(1:1).EQ.'I') THEN
          TINPUT=.TRUE.
        ELSE IF(STRING(1:1).EQ.'M') THEN
          TCRYSTAL=.FALSE.
        ELSE
          CALL ERROR$MSG('ARGUMENT NOT RECOGNIZED')
          CALL ERROR$MSG('OBTAIN ARGUMENT LIST USING -H ARGUMENT')
          CALL ERROR$CHVAL('ARGUMENT',TRIM(STRING))
          CALL ERROR$I4VAL('ARGUMENT NUMBER',I)
          CALL ERROR$STOP('MAIN')
        END IF
      ENDDO
!
!     == ROOTNAME ==============================================================
      CALL GET_COMMAND_ARGUMENT(NARGS,ROOTNAME) !LAST ARGUMENT IS THE ROOT NAME
      WRITE(*,FMT='("ROOTNAME: ",A)')TRIM(ROOTNAME)
      IF(LEN(TRIM(ROOTNAME)).EQ.0) THEN
        STOP 'NO ROOTNAME SUPPLIED'
      END IF
!
      I=INDEX(ROOTNAME,'/',BACK=.TRUE.)
      OBJECTNAME=ROOTNAME(I+1:)
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('PROT',.TRUE.,-'.SPROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
!
      IF(TINPUT) THEN
        CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC')
      ELSE
        CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC_OUT')
      END IF
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('CSSR',.TRUE.,-'.CSSR')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('CSSR','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('MMCSSR',.TRUE.,-'_MM.CSSR')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('MMCSSR','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('SHCSSR',.TRUE.,-'_SH.CSSR')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('SHCSSR','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('CML',.TRUE.,-'.CML')
      CALL FILEHANDLER$SETSPECIFICATION('CML','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('CML','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('CML','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('CML','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('XYZ',.TRUE.,-'.XYZ')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('XYZ','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('POSCAR',.TRUE.,-'.POSCAR.VASP')
      CALL FILEHANDLER$SETSPECIFICATION('POSCAR','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('POSCAR','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('POSCAR','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('POSCAR','FORM','FORMATTED')
!
!     ==========================================================================
!     == READ STRUCTURE FILE                                                  ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
!     -- MARK ELEMENTS AS READ FROM INPUT IF READING FROM STRC BUT NOT
!     -- WHEN READING FROM STRC_OUT
      IF(TINPUT)CALL LINKEDLIST$MARK(LL_STRC,1)
!
!     ==========================================================================
!     == GET LENGTH UNIT                                                      ==
!     ==========================================================================
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')

      CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT[AA]',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_STRC,'LUNIT[AA]',1,RUNIT)
        RUNIT=RUNIT*ANGSTROM
      ELSE
        CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,RUNIT)
      END IF
!
!     ==========================================================================
!     == GET LATTICE VECTORS                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'LATTICE')
      CALL LINKEDLIST$GET(LL_STRC,'T',1,RBAS)
      RBAS=RBAS*RUNIT
!
!     ==========================================================================
!     ==  READ ATOM DATA FROM STRC FILE                                       ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT)
      ALLOCATE(NAME(NAT))
      ALLOCATE(R(3,NAT))
      ALLOCATE(Q(NAT)) ; Q(:)=0.D0
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
        CALL LINKEDLIST$GET(LL_STRC,'R',1,R(:,IAT))
        R(:,IAT)=R(:,IAT)*RUNIT
        CALL LINKEDLIST$EXISTD(LL_STRC,'Q',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'Q',1,Q(IAT))
        ELSE
          Q(IAT)=0.D0
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,NAME(IAT))
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
!
!     ==========================================================================
!     ==  READ MM ATOM DATA FROM STRC FILE                                    ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'QM-MM',1,TQMMM)
      IF(TQMMM) THEN
        CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
        CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NATMM)
        ALLOCATE(MMNAME(NATMM))
        ALLOCATE(MMR(3,NATMM)) ;MMR(:,:)=0.D0
        ALLOCATE(MMQ(NATMM)) ;MMQ(:)=0.D0
        DO IAT=1,NATMM
          CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
          CALL LINKEDLIST$GET(LL_STRC,'NAME',1,MMNAME(IAT))
          CALL LINKEDLIST$EXISTD(LL_STRC,'R',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'R',1,MMR(:,IAT))
            MMR(:,IAT)=MMR(:,IAT)*RUNIT
          ELSE
            CALL LINKEDLIST$GET(LL_STRC,'QMATOM',1,STRING)
            DO IAT1=1,NAT
              IF(STRING.EQ.NAME(IAT1)) THEN
                MMR(:,IAT)=R(:,IAT1)
                EXIT
              END IF
            ENDDO
          END IF
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO  
      END IF
!
!     ==========================================================================
!     == WRITE HEADER TO PROTOCOL FILE .SPROT                                 ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='()')
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(80("*"),T15,A)') &
     &                   "          STRUCTURE ANALYSIS                "
      WRITE(NFILO,FMT='(80("*"),T15,A)') &
     &                   "   FOR THE PROJECTOR-AUGMENTED WAVE METHOD   "
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(T10,A)') &
     &                  "P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY"
      WRITE(NFILO,FMT='(T10,A)') &
     &                  "DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3"
      WRITE(NFILO,*)
      IF(TINPUT) THEN
        WRITE(NFILO,FMT='("INPUT STRUCTURE FILE READ")')
      END IF
      IF(TINPUT) THEN
        WRITE(NFILO,FMT='("STRUCTURE FILE READ:",A)')TRIM(ROOTNAME)//-'.STRC'
      ELSE
        WRITE(NFILO,FMT='("STRUCTURE FILE READ:",A)') &
    &                    TRIM(ROOTNAME)//-'.STRC_OUT'
      END IF
      IF(TCRYSTAL) THEN
        WRITE(NFILO,FMT='("STRUCTURE INPRETED AS CRYSTAL")')
      ELSE
        WRITE(NFILO,FMT='("STRUCTURE INPRETED AS MOLECULE")')
      END IF
   
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$REPORT_UNUSED(LL_STRC,NFILO)
!
!     ==========================================================================
!     == WRITE BOND LENGTH TO PROTOCOL                                        ==
!     ==========================================================================
      CALL REPORTLOCALSTRUCTURE(RBAS,NAT,NAME,R,TCRYSTAL)
!
!     ==========================================================================
!     == WRITE SECTION FOR STRC FILE                                          ==
!     ==========================================================================
      CALL WRITESTRCFILE(RBAS,NAT,NAME,R,TCRYSTAL,RUNIT)
!
!     ==========================================================================
!     ==  WRITE CML FILE AS INPUT FOR AVOGADRO VIEWER                         ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CML',NFIL)
      CALL WRITEAVOGADRO(NFIL,OBJECTNAME,RBAS,NAT,NAME,R,Q,TCRYSTAL)
!
!     ==========================================================================
!     ==  WRITE XYZ FILE                                                      ==
!     ==========================================================================
      CALL WRITEXYZ(RBAS,NAT,NAME,R,TCRYSTAL,NDUP,ROOTNAME)
!
!     ==========================================================================
!     ==  WRITE STRUCTURE IN POSCAR FORMAT (OF THE VASP PACKAGE)              ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('POSCAR',NFIL)
      CALL WRITEPOSCAR(NFIL,OBJECTNAME,RBAS,NAT,NAME,R)
!
!     ==========================================================================
!     == CONVERT DATA TO ANGSTROM AND ELECTRON CHARGES                        ==
!     ==========================================================================
      RBAS=RBAS/ANGSTROM
      R=R/ANGSTROM
      Q=-Q
      IF(TQMMM) THEN
        MMR=MMR/ANGSTROM
      END IF
!
!     ==========================================================================
!     == WRITE CSSR FILE                                                      ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CSSR',NFIL)
      CALL WRITECSSR(NFIL,OBJECTNAME,RBAS,NAT,NAME,R,Q,TCRYSTAL)
      IF(TQMMM) THEN
        CALL FILEHANDLER$UNIT('MMCSSR',NFIL)
        CALL WRITECSSR(NFIL,TRIM(OBJECTNAME)//-'MM',RBAS,NATMM,MMNAME,MMR &
     &                                                            ,MMQ,TCRYSTAL)
      END IF

      IF(TCRYSTAL) THEN
        WRITE(*,FMT='("CRYSTAL OUTPUT PRODUCED")')
      ELSE
        WRITE(*,FMT='("MOLECULAR OUTPUT PRODUCED")')
      END IF
      CALL FILEHANDLER$CLOSEALL
      WRITE(*,FMT='("======= TASK FINISHED ========")')
      CALL ERROR$NORMALSTOP()
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITESTRCFILE(RBAS,NAT,NAME,R,TCRYSTAL,RUNIT)
!     **************************************************************************
!     ** REPORT BOND ANGLES AND BOND LENGTHS                                  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      CHARACTER(*),INTENT(IN):: NAME(NAT)
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(IN) :: RUNIT
      LOGICAL(4),INTENT(IN) :: TCRYSTAL
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: IAT
      REAL(8)                     :: ANGSTROM
!     **************************************************************************
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='(80("="))')
      WRITE(NFILO,FMT='(80("="),T20," SEGMENT FOR STRUCTURE INPUT FILE  ")')
      WRITE(NFILO,FMT='(80("="),T20,"          IN UNITS OF LUNIT        ")')
      WRITE(NFILO,FMT='(80("="))')
      WRITE(NFILO,FMT=*)
      WRITE(NFILO,FMT='(T3,"LUNIT=    ",F10.5)')RUNIT
      WRITE(NFILO,FMT='(T3,"LUNIT[AA]=",F10.5)')RUNIT/ANGSTROM
      WRITE(NFILO,FMT=*)
      WRITE(NFILO,FMT='(T3,"!LATTICE",T25," T= ",3F12.5)')RBAS(:,1)/RUNIT
      WRITE(NFILO,FMT='(T29,3F12.5)')RBAS(:,2)/RUNIT
      WRITE(NFILO,FMT='(T29,3F12.5," !END")')RBAS(:,3)/RUNIT
      WRITE(NFILO,FMT=*)
      DO IAT=1,NAT
        WRITE(NFILO,FMT='(T3,"!ATOM NAME= ",T15,A,T25," R= ",3F12.5," !END")') &
     &                   "'"//TRIM(NAME(IAT))//"'",R(:,IAT)/RUNIT
      ENDDO
      WRITE(NFILO,FMT=*)
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEXYZ(RBAS,NAT,NAME,R,TCRYSTAL,NDUP,TITLE)
!     **************************************************************************
!     **  WRITES AN XYZ FILE FOR THE STRUCTURE                                **
!     **                                                                      **
!     **  CAUTION: ASSUMES THAT THE ATOM NAME STARTS WITH THE ELEMENT SYMBOL!!**
!     **                                                                      **
!     ** see "extended XYZ" format:                                           **
!     ** https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz       **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN):: NAT
      REAL(8)     ,INTENT(IN):: RBAS(3,3)
      CHARACTER(*),INTENT(IN):: NAME(NAT)
      CHARACTER(*),INTENT(IN):: TITLE
      REAL(8)     ,INTENT(IN):: R(3,NAT)
      INTEGER(4)  ,INTENT(IN):: NDUP(3)
      LOGICAL(4)  ,INTENT(IN):: TCRYSTAL
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: IAT
      INTEGER(4)            :: IT1,IT2,IT3
      REAL(8)               :: T(3)
      REAL(8)               :: ANGSTROM
      CHARACTER(2)          :: EL
      CHARACTER(160)        :: STRING
!     **************************************************************************
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL FILEHANDLER$UNIT('XYZ',NFIL)
      REWIND(NFIL)
      IF(TCRYSTAL) THEN
        WRITE(NFIL,FMT='(I10)')NAT*NDUP(1)*NDUP(2)*NDUP(3)
!!$        WRITE(STRING,FMT='(I1,A,I1,A,I1)')NDUP(1),-'X',NDUP(2),-'X',NDUP(3)
!!$        STRING=TRIM(STRING)//' AS CLUSTER'        

        WRITE(STRING,FMT='(9F10.5)')RBAS(:,1)*NDUP(1)/ANGSTROM, &
     &        RBAS(:,2)*NDUP(2)/ANGSTROM,RBAS(:,3)*NDUP(3)/ANGSTROM
        STRING=+'L'//-'ATTICE="'//TRIM(ADJUSTL(STRING))//'" '//+'P'// &
     &         -'ROPERTIES=SPECIES:'//+'S'//-':1:POS:'//+'R'//':3'
        WRITE(NFIL,FMT='(A)')STRING

      ELSE
        WRITE(NFIL,FMT='(I10)')NAT
        STRING=' '
        WRITE(NFIL,FMT='(A)')TRIM(TITLE)//' '//TRIM(STRING)
      END IF

      DO IT1=1,NDUP(1)
        DO IT2=1,NDUP(2)
           DO IT3=1,NDUP(3)
             T(:)=RBAS(:,1)*REAL(IT1-1) &
     &           +RBAS(:,2)*REAL(IT2-1) &
     &           +RBAS(:,3)*REAL(IT3-1)
             DO IAT=1,NAT
               EL=NAME(IAT)(1:2)
               IF(EL(2:2).EQ.'_')EL(2:2)=' '
               WRITE(NFIL,FMT='(A2,2X,3(F10.5,1X),2X,A)')EL, &
     &              (R(:,IAT)+T(:))/ANGSTROM
             ENDDO
           ENDDO
         ENDDO
      ENDDO
      CALL FILEHANDLER$CLOSE('XYZ')
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE REPORTLOCALSTRUCTURE(RBAS,NAT,NAME,R,TCRYSTAL)
!     **************************************************************************
!     ** REPORT BOND ANGLES AND BOND LENGTHS                                  **
!     **************************************************************************
      IMPLICIT NONE
      TYPE BOND_TYPE
        INTEGER(4) :: IAT1
        INTEGER(4) :: IAT2
        INTEGER(4) :: IT(3)
        REAL(8)    :: DR(3)
      END TYPE BOND_TYPE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      CHARACTER(*),INTENT(IN):: NAME(NAT)
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      LOGICAL(4),INTENT(IN) :: TCRYSTAL
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: NT
      REAL(8)               :: DISX
      INTEGER(4)            :: IAT,IAT1,IAT2,IT1,IT2,IT3
      REAL(8)               :: T1(3),T2(3),T3(3)
      REAL(8)               :: DR(3),DIS
      INTEGER(4),PARAMETER  :: NBONDX=100
      TYPE(BOND_TYPE)       :: BOND(NBONDX)
      TYPE(BOND_TYPE)       :: SVARBOND
      REAL(8)               :: DISARR(NBONDX)
      REAL(8)               :: ANGLE(NBONDX,NBONDX)
      REAL(8)               :: ANGLEN
      REAL(8)               :: SVAR
      INTEGER(4)            :: NBOND,I,J,K
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: IT(3)
      REAL(8)               :: ANGSTROM
      REAL(8)               :: RBASIN(3,3)
      CHARACTER(64)         :: EXTENDEDNAME1,EXTENDEDNAME2
!     **************************************************************************
      ANGLEN=PI/180.D0*59.D0    ! MIN ANGLE=60 DEG
      DISX=6.5D0                 ! MAX DISTANCE
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
!
!     ==========================================================================
!     ==  ANALYZE LATTICE
!     ==========================================================================
      IF(TCRYSTAL) THEN
        WRITE(NFILO,FMT='(80("="))')
        WRITE(NFILO,FMT='(80("="),T20,"      LATTICE CONSTANTS AND        ")')
        WRITE(NFILO,FMT='(80("="),T20," POSITIONS IN RELATIVE COORDINATES ")')
        WRITE(NFILO,FMT='(80("="))')
        WRITE(NFILO,FMT=*)
!
!       == WRITE LATTICE VECTORS IN CARTESIAN COORDINATES ======================
        WRITE(NFILO,FMT='("LATTICE VECTORS IN CARTESIAN COORDINATES:")')
        WRITE(NFILO,FMT='(A,T30,3F10.5," ANGSTROM")')'T_A',RBAS(:,1)/ANGSTROM
        WRITE(NFILO,FMT='(A,T30,3F10.5," ANGSTROM")')'T_B',RBAS(:,2)/ANGSTROM
        WRITE(NFILO,FMT='(A,T30,3F10.5," ANGSTROM")')'T_C',RBAS(:,3)/ANGSTROM
        WRITE(NFILO,FMT=*)
!
!       == CONSTRACT LATTICE CONSTANTS AND ANGLES ==============================
        DO I=1,3
          DISARR(I)=SQRT(SUM(RBAS(:,I)**2))
        ENDDO
        DO I=1,3
          J=1+MODULO(I+1-1,3)
          K=1+MODULO(I+2-1,3)
          SVAR=DOT_PRODUCT(RBAS(:,I),RBAS(:,J))/(DISARR(I)*DISARR(J))
          SVAR=MIN(1.D0,MAX(-1.D0,SVAR))
          SVAR=ACOS(SVAR)
          DISARR(3+K)=SVAR
        ENDDO
        WRITE(NFILO,FMT='("LATTICE CONSTANTS:")')
        WRITE(NFILO,FMT='("A,B,C = 1ST, 2ND AND 3RD LATTICE VECTOR")')
        WRITE(NFILO,FMT='("ALPHA = ANGLE BETWEEN B AND C")')
        WRITE(NFILO,FMT='("BETA  = ANGLE BETWEEN C AND A")')
        WRITE(NFILO,FMT='("GAMMA = ANGLE BETWEEN A AND B")')
        WRITE(NFILO &
     &            ,FMT='("LATTICE CONSTANTS A,B,C: ",T40,3F10.5," ANGSTROM")') &
     &                                                      DISARR(1:3)/ANGSTROM
        WRITE(NFILO &
     &      ,FMT='("LATTICE ANGLES ALPHA,BETA,GAMMA: ",T40,3F10.5," DEGREE")') &
     &                                             DISARR(4:6)/PI*180.D0
        WRITE(NFILO,FMT=*)
!
!       == WRITE ATOMIC POSITIONS IN RELATIVE COORDINATES ======================
        CALL LIB$INVERTR8(3,RBAS,RBASIN)
        WRITE(NFILO,FMT='("ATOMIC POSITIONS IN RELATIVE COORDINATES:")')
        WRITE(NFILO,FMT='("(IMAGE FROM THE FIRST UNIT CELL)")')
        DO IAT=1,NAT
          DR(:)=MATMUL(RBASIN,R(:,IAT))
          T1(:)=DR(:)
          DO I=1,3
            DR(I)=MODULO(DR(I)+1.D-5,1.D0)-1.D-5
          ENDDO
          T1(:)=DR(:)-T1(:)
          IT(:)=NINT(T1)
          CALL EXTENDNAME(NAME(IAT),IT,EXTENDEDNAME1)
          WRITE(NFILO,FMT='(A,T30,3F10.5)')TRIM(EXTENDEDNAME1),DR
        ENDDO
        WRITE(NFILO,FMT=*)
      END IF
!
!     ==========================================================================
!     == 
!     ==========================================================================
      WRITE(NFILO,FMT='(80("="))')
      WRITE(NFILO,FMT='(80("="),T20,"  BOND LENGTHS AND ANGLES    ")')
      WRITE(NFILO,FMT='(80("="))')
      WRITE(NFILO,FMT='("NEIGHBORS WITH BOND LENGTHS BEYOND ",F5.1," ANGSTROM ARE DROPPED")') &
     &                DISX/ANGSTROM
      WRITE(NFILO,FMT='("NEIGHBORS WITH ANGLES BELOW ",F5.0," DEGREE ARE DROPPED")') &
     &                ANGLEN/PI*180.D0
      WRITE(NFILO,FMT='("CAUTION: ANGLE REQUIREMENT MAY RESULT IN INCOMPLETE SHELLS")')
      WRITE(NFILO,FMT=*)
      NT=0
      IF(TCRYSTAL) NT=2
      DO IAT1=1,NAT
!       ========================================================================
!       ==  COLLECT BONDS TO ATOM IAT1                                        ==
!       ========================================================================
        NBOND=0
        DO IAT2=1,NAT
          DO IT1=-NT,NT 
            T1(:)=RBAS(:,1)*REAL(IT1,KIND=8)
            DO IT2=-NT,NT 
              T2(:)=RBAS(:,2)*REAL(IT2,KIND=8)
              DO IT3=-NT,NT 
                T3(:)=RBAS(:,3)*REAL(IT3,KIND=8)
                DR(:)=R(:,IAT2)-R(:,IAT1)+T1+T2+T3
                DIS=SQRT(SUM(DR**2))
                IF(IAT1.EQ.IAT2.AND.IT1.EQ.0.AND.IT2.EQ.0.AND.IT3.EQ.0) CYCLE
                IF(DIS.GT.DISX) CYCLE
                NBOND=NBOND+1
                IF(NBOND.GT.NBONDX) THEN
                  CALL ERROR$MSG('ARRAY SIZE EXCEEDED')
                  CALL ERROR$I4VAL('NBOND',NBOND)
                  CALL ERROR$I4VAL('NBONDX',NBONDX)
                  CALL ERROR$STOP('REPORTLOCALSTRUCTURE')
                END IF
                BOND(NBOND)%IAT1=IAT1
                BOND(NBOND)%IAT2=IAT2
                BOND(NBOND)%IT(1)=IT1
                BOND(NBOND)%IT(2)=IT2
                BOND(NBOND)%IT(3)=IT3
                BOND(NBOND)%DR(:)=DR(:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  DETERMINE BOND LENGTH ARRAY DISARR                                ==
!       ========================================================================
        DO I=1,NBOND
          DISARR(I)=SQRT(SUM(BOND(I)%DR**2))
        ENDDO
!
!       ========================================================================
!       ==  ORDER BONDS ACCORDING TO  INCREASING LENGTH                       ==
!       ========================================================================
        DO I=1,NBOND
          J=I
          DO K=I+1,NBOND
            IF(DISARR(K).LT.DISARR(J)) J=K
          ENDDO
          IF(J.EQ.I) CYCLE
          SVAR=DISARR(J)
          DISARR(J)=DISARR(I)
          DISARR(I)=SVAR
          SVARBOND=BOND(J)
          BOND(J)=BOND(I)
          BOND(I)=SVARBOND
        ENDDO
        IF(DISARR(1).LT.3.D-3) THEN
          CALL ERROR$MSG('UNPHYSICAL BOND DISTANCE')
          CALL ERROR$R8VAL('MIN(DIS)',DISARR(1))
          CALL ERROR$I4VAL('IAT1',BOND(1)%IAT1)
          CALL ERROR$I4VAL('IAT2',BOND(1)%IAT2)
          CALL ERROR$I4VAL('IT1',BOND(1)%IT(1))
          CALL ERROR$I4VAL('IT2',BOND(1)%IT(2))
          CALL ERROR$I4VAL('IT3',BOND(1)%IT(3))
          CALL ERROR$STOP('REPORTLOCALSTRUCTURE')
        END IF
!
!       ========================================================================
!       == DETERMINE BOND ANGLES                                              ==
!       ========================================================================
        DO I=1,NBOND
          ANGLE(I,I)=0.D0
          DO J=I+1,NBOND
            SVAR=DOT_PRODUCT(BOND(I)%DR,BOND(J)%DR)/DISARR(I)/DISARR(J)
            SVAR=MIN(1.D0,MAX(-1.D0,SVAR))  ! AVOID NANS
            SVAR=ACOS(SVAR)
            ANGLE(I,J)=SVAR
            ANGLE(J,I)=SVAR
          ENDDO
        ENDDO
!
!       ========================================================================
!       == DROP NEIGHBORS WITH TOO SMALL BOND ANGLES                          ==
!       ========================================================================
        DO I=1,NBOND
          IF(DISARR(I).LT.0.D0) CYCLE
          DO J=I+1,NBOND
            IF(DISARR(J).LT.0.D0) CYCLE
            IF(ANGLE(I,J).LT.ANGLEN) DISARR(J)=-1.D0
          ENDDO
        ENDDO
        J=0
        DO I=1,NBOND
          IF(DISARR(I).LT.0.D0) CYCLE
          J=J+1
          IF(J.EQ.I) CYCLE
          ANGLE(:,J)=ANGLE(:,I) 
          ANGLE(J,:)=ANGLE(I,:) 
          BOND(J)   =BOND(I)
          DISARR(J)=DISARR(I)
        ENDDO
        NBOND=J
!
!       ========================================================================
!       == REPORT TO FILE BOND-LENGTHS AND ANGLES                             ==
!       ========================================================================
        WRITE(NFILO,FMT='(80("-"),T30,"  ",A,"   ")')TRIM(NAME(IAT1))
        DO I=1,NBOND
          CALL EXTENDNAME(NAME(BOND(I)%IAT2),BOND(I)%IT,EXTENDEDNAME1)
          WRITE(NFILO,FMT='(A,T19,A,T57,"LENGTH: ",F5.3," ANGSTROM")') &
     &                  TRIM(NAME(BOND(I)%IAT1)),TRIM(EXTENDEDNAME1) &
     &                 ,DISARR(I)/ANGSTROM 
        ENDDO
        DO I=1,NBOND
          CALL EXTENDNAME(NAME(BOND(I)%IAT2),BOND(I)%IT,EXTENDEDNAME1)
          DO J=I+1,NBOND
            CALL EXTENDNAME(NAME(BOND(J)%IAT2),BOND(J)%IT,EXTENDEDNAME2)
            WRITE(NFILO,FMT='(A,T19,A,T38,A,T57,"ANGLE: ",F5.1," DEGREE")') &
     &              TRIM(EXTENDEDNAME1),TRIM(NAME(BOND(I)%IAT1)) &
     &             ,TRIM(EXTENDEDNAME2) &
     &             ,ANGLE(I,J)/PI*180.
          ENDDO
        ENDDO
        DO I=1,NBOND
          CALL EXTENDNAME(NAME(BOND(I)%IAT2),BOND(I)%IT,EXTENDEDNAME1)
          WRITE(NFILO,FMT='(A,T19,A,T38,"DIRECTION: ",3F10.5)') &
     &                  TRIM(NAME(BOND(I)%IAT1)),TRIM(EXTENDEDNAME1) &
     &                 ,BOND(I)%DR/SQRT(SUM(BOND(I)%DR**2))
        ENDDO
      ENDDO

      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTENDNAME(NAME,IT,EXTENDEDNAME)
!     **************************************************************************
!     ** CONSTRUCTS THE ATOM NAME IN EXTENDED NOTATION, THAT IS INCLUDING     **
!     ** THE LATTICE TRANSLATIONS                                             **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: IT(3)
      CHARACTER(*),INTENT(OUT):: EXTENDEDNAME
      INTEGER(4)              :: I
      CHARACTER(8)            :: STRING
!     **************************************************************************
      IF(SUM(IT**2).EQ.0) THEN
        EXTENDEDNAME=NAME
        RETURN
      END IF
!
      EXTENDEDNAME=''
      DO I=3,1,-1
        WRITE(STRING,FMT='(I8)',ERR=100)IT(I)
        EXTENDEDNAME=TRIM(ADJUSTL(STRING))//ADJUSTL(EXTENDEDNAME)
      ENDDO
      EXTENDEDNAME=TRIM(ADJUSTL(NAME))//':'//TRIM(EXTENDEDNAME)
      RETURN
100   CONTINUE
      CALL ERROR$MSG('ERROR CONSTRUCTING EXTENDED ATOM NAME')
      CALL ERROR$CHVAL('NAME',NAME)
      CALL ERROR$I4VAL('LEN(NAME)',LEN(NAME))
      CALL ERROR$I4VAL('LEN(EXTENDEDNAME)',LEN(EXTENDEDNAME))
      CALL ERROR$STOP('EXTENDEDNAME')
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITECSSR(NFIL,OBJECTNAME,RBAS,NAT,NAME,R,Q,TCRYSTAL)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTERFACE OPERATOR (.DYAD.)
        FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
          REAL(8), INTENT(IN) :: R1(3)
          REAL(8), INTENT(IN) :: R2(3)
          REAL(8)             :: R3(3)
        END FUNCTION DYADISCHES_PRODUCT
      END INTERFACE 
      INTEGER(4)  ,INTENT(IN) :: NFIL
      LOGICAL(4)  ,INTENT(IN) :: TCRYSTAL
      CHARACTER(*),INTENT(IN) :: OBJECTNAME
      INTEGER(4)  ,INTENT(IN) :: NAT
!     REAL(8)     ,INTENT(IN) :: RBAS(3,3)
      REAL(8)                 :: RBAS(3,3)
      CHARACTER(*),INTENT(IN) :: NAME(NAT)
      REAL(8)                 :: R(3,NAT)
!     REAL(8)     ,INTENT(IN) :: R(3,NAT)
      REAL(8)     ,INTENT(IN) :: Q(NAT)
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)              :: NEIGH(8,NAT)
      REAL(8)                 :: RBASINV(3,3) ! RBAS**(-1)
      REAL(8)                 :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)                 :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      REAL(8)                 :: DET
      INTEGER(4)              :: IAT
      REAL(8)                 :: VEC(3),RBASNEU(3,3)
!     **************************************************************************
      NEIGH(:,:)=0
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
        WRITE(NFIL,FMT='(T39,3F8.3 &
     &   /T22,3F8.3,T50,"SPGR = 1 P 1",T72,"OPT = 1" &
      &   /I4,"   1 CREATED BY PAW    " &
     &   /"     0 ",A4,": ",A4)') &
     &   A,B,C,ALPHA,BETA,GAMMA,NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      ELSE
        WRITE(NFIL,FMT='(//I4,"   1 CREATED BY PAW    " &
     &   /"     0 ",A4,": ",A4)')NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      END IF
      DO IAT=1,NAT
        IF(TCRYSTAL) THEN
          VEC=R(:,IAT)
          VEC=MATMUL(RBASINV,VEC)+100.D0
          VEC=MOD(VEC,1.D0)
          VEC=MATMUL(RBASNEU,VEC)
          WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
     &           IAT,NAME(IAT),VEC(:),NEIGH(:,IAT),Q(IAT)

        ELSE
           WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
      &           IAT,NAME(IAT),R(:,IAT),NEIGH(:,IAT),Q(IAT)
         END IF
      ENDDO

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEPOSCAR(NFIL,OBJECTNAME,RBAS,NAT,NAME,R)
!     **************************************************************************
!     ** WRITE STRUCTURE TO THE POSCAR FORMAT FILE USED BY VASP               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: OBJECTNAME
      INTEGER(4)  ,INTENT(IN) :: NAT
      REAL(8)     ,INTENT(IN) :: RBAS(3,3)
      CHARACTER(*),INTENT(IN) :: NAME(NAT)
      REAL(8)     ,INTENT(IN) :: R(3,NAT)
      INTEGER(4)  ,PARAMETER  :: NSPX=100
      INTEGER(4)              :: NSP=0
      CHARACTER(2)            :: SPNAME(NSPX)
      INTEGER(4)              :: SPNUM(NSPX)
      REAL(8)                 :: ANGSTROM
      INTEGER(4)              :: IAT
!     **************************************************************************
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)

      WRITE(NFIL,*)TRIM(ADJUSTL(OBJECTNAME))
      WRITE(NFIL,*)1.D0     ! SCALE FACTOR TO ANGSTROM
      WRITE(NFIL,FMT='(3F20.8)')RBAS/ANGSTROM
      NSP=1
      SPNAME(NSP)=NAME(1)(1:2)
      SPNUM(NSP)=1
      DO IAT=2,NAT
        IF(NAME(IAT)(1:2).EQ.SPNAME(NSP)) THEN
          SPNUM(NSP)=SPNUM(NSP)+1
        ELSE
          NSP=NSP+1
          IF(NSP.GT.NSPX) THEN
            CALL ERROR$MSG('NUMBER OF ATOM TYPES EXCEEDS HARDCODED LIMIT')
            CALL ERROR$I4VAL('NSPX',NSPX)
            CALL ERROR$STOP('WRITEPOSCAR')
          END IF
          SPNAME(NSP)=NAME(IAT)(1:2)
          SPNUM(NSP)=1
        END IF
      ENDDO
      WRITE(NFIL,FMT='(100(A," "))')SPNAME(:NSP)
      WRITE(NFIL,*)SPNUM(:NSP)
      WRITE(NFIL,FMT='("CARTESIAN")')
      DO IAT=1,NAT
        WRITE(NFIL,FMT='(3F10.5)')R(:,IAT)/ANGSTROM
      ENDDO

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEAVOGADRO(NFIL,OBJECTNAME,RBAS,NAT,NAME,R,Q,TCRYSTAL)
      USE STRINGS_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTERFACE OPERATOR (.DYAD.)
        FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
          REAL(8), INTENT(IN) :: R1(3)
          REAL(8), INTENT(IN) :: R2(3)
          REAL(8)             :: R3(3)
        END FUNCTION DYADISCHES_PRODUCT
      END INTERFACE 
      INTEGER(4)  ,INTENT(IN) :: NFIL
      LOGICAL(4)  ,INTENT(IN) :: TCRYSTAL
      CHARACTER(*),INTENT(IN) :: OBJECTNAME
      INTEGER(4)  ,INTENT(IN) :: NAT
      REAL(8)     ,INTENT(IN) :: RBAS(3,3)
      CHARACTER(*),INTENT(IN) :: NAME(NAT)
      REAL(8)     ,INTENT(IN) :: R(3,NAT)
      REAL(8)     ,INTENT(IN) :: Q(NAT)
      REAL(8)     ,PARAMETER  :: TOLBOND=5.D-2
      REAL(8)     ,PARAMETER  :: TOLBOX=1.D-1
      REAL(8)                 :: X(3,NAT)
      CHARACTER(2)            :: EL(NAT)
      REAL(8)                 :: PI
      REAL(8)                 :: RBASINV(3,3) ! RBAS**(-1)
      REAL(8)                 :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)                 :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECTORS
      REAL(8)                 :: DET
      INTEGER(4)              :: IAT,I,J,K,L
      REAL(8)                 :: SVAR,RBASNEU(3,3)
      REAL(8)                 :: ANGSTROM,DEGREE
      INTEGER(4) ,ALLOCATABLE :: MAP(:,:)
      INTEGER(4)              :: NMAP,IMAP,J1,J2,J3,NDUP1,NBOND
      INTEGER(4)              :: NBONDX
      CHARACTER(8)           :: DISTID(-1:1,-1:1,-1:1)
      INTEGER(4) ,ALLOCATABLE:: BOND(:,:)
      REAL(8)                 :: RCOV(NAT)
      REAL(8)                 :: RCOV1,RCOV2
      REAL(8)                 :: R1(3),R2(3)
!     ******************************************************************
      CALL CONSTANTS('PI',PI)
      CALL CONSTANTS('DEGREE',DEGREE)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      DO I=-1,1
        DO J=-1,1
          DO K=-1,1
            WRITE(DISTID(I,J,K),FMT='("(",3I2,")")')I,J,K
            DO L=2,6,2
              IF(DISTID(I,J,K)(L:L).EQ.' ')DISTID(I,J,K)(L:L)='+'
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == CONVERT DATA TO UNITS OF LATTICE VECTORS                     ==
!     ==================================================================
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
      RBASNEU(:,:)=0.D0
      RBASNEU(3,3)=C
      RBASNEU(3,2)=B*COS(ALPHA)
      RBASNEU(2,2)=SQRT(B**2-RBASNEU(3,2)**2)
      RBASNEU(3,1)=A*COS(BETA)
      RBASNEU(2,1)=(A*B*COS(GAMMA)-RBASNEU(3,1)*RBASNEU(3,2))/RBASNEU(2,2)
      RBASNEU(1,1)=SQRT(A**2-RBASNEU(2,1)**2-RBASNEU(3,1)**2)
!
!     ==========================================================================
!     == WRITE AVOGADRO FILE                                                  ==
!     ==========================================================================
!
!     ==========================================================================
!     == WRITE UNIT CELL                                                      ==
!     ==========================================================================
      WRITE(NFIL,FMT='(A)')-'<MOLECULE>'
      WRITE(NFIL,FMT='(A)')-'<CRYSTAL>'
      WRITE(NFIL,FMT='(A,F10.5,A)')-'<SCALAR TITLE="A" UNITS="UNITS:ANGSTROM">',A/ANGSTROM,-'</SCALAR>'
      WRITE(NFIL,FMT='(A,F10.5,A)')-'<SCALAR TITLE="B" UNITS="UNITS:ANGSTROM">',B/ANGSTROM,-'</SCALAR>'
      WRITE(NFIL,FMT='(A,F10.5,A)')-'<SCALAR TITLE="C" UNITS="UNITS:ANGSTROM">',C/ANGSTROM,-'</SCALAR>'
      WRITE(NFIL,FMT='(A,F10.5,A)')-'<SCALAR TITLE="ALPHA" UNITS="UNITS:DEGREE">',ALPHA/DEGREE,-'</SCALAR>'
      WRITE(NFIL,FMT='(A,F10.5,A)')-'<SCALAR TITLE="BETA" UNITS="UNITS:DEGREE">',BETA/DEGREE,-'</SCALAR>'
      WRITE(NFIL,FMT='(A,F10.5,A)')-'<SCALAR TITLE="GAMMA" UNITS="UNITS:DEGREE">',GAMMA/DEGREE,-'</SCALAR>'
      WRITE(NFIL,FMT='(A)')-'</CRYSTAL>'
!
!     ==========================================================================
!     == CONSTRUCT RELATIVE COORDINATES AND MAP INTO FIRST UNIT CELL          ==
!     ==========================================================================
      NMAP=0
      DO IAT=1,NAT
         X(:,IAT)=MATMUL(RBASINV,R(:,IAT))
         IF(TCRYSTAL) THEN
           X(:,IAT)=MOD(X(:,IAT)+100.D0,1.D0)
           NDUP1=0
           DO I=1,3
             IF(X(I,IAT).LT.1.D-1.OR.X(I,IAT).GT.1.D0-1.D-1) NDUP1=NDUP1+1
           ENDDO
           NMAP=NMAP+2**NDUP1
         ELSE
           NMAP=NMAP+1
         END IF
         EL(IAT)=NAME(IAT)(1:2)
         IF(EL(IAT)(2:2).EQ.'_')EL(IAT)(2:2)=' '
         CALL PERIODICTABLE$GET(EL(IAT),'R(COV)',RCOV(IAT))
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT PERIODIC IMAGES                                            ==
!     ==========================================================================
      ALLOCATE(MAP(4,NMAP))
      MAP(:,:)=0
      IF(TCRYSTAL) THEN
        IMAP=1
        DO IAT=1,NAT
          MAP(1,IMAP)=IAT
          NDUP1=1
          DO I=1,3
            IF(X(I,IAT).LT.TOLBOX.OR.X(I,IAT).GT.1.D0-TOLBOX) THEN
              J1=IMAP
              J2=IMAP+NDUP1
              J3=IMAP+2*NDUP1
              MAP(:,J2:J3-1)=MAP(:,J1:J2-1)
              IF(X(I,IAT).LT.0.5D0) THEN
                MAP(1+I,J2:J3-1)=1
              ELSE
                MAP(1+I,J2:J3-1)=-1
              END IF
              NDUP1=2*NDUP1
            END IF
          ENDDO
          IMAP=IMAP+NDUP1
        ENDDO
      ELSE
        DO IAT=1,NAT
          MAP(1,IAT)=IAT
        ENDDO
      END IF
!
!     ==========================================================================
!     == WRITE ATOMIC POSITIONS                                               ==
!     ==========================================================================
      WRITE(NFIL,FMT='("<",A)')-'ATOM'//+'A'//-'RRAY'
      WRITE(NFIL,FMT='(A)')-'ATOM'//+'ID'//'="'
!      WRITE(NFIL,*)(TRIM(NAME(MAP(1,IMAP)))//DISTID(MAP(2,IMAP),MAP(3,IMAP),MAP(4,IMAP))//' ',IMAP=1,NMAP)
      DO IMAP=1,NMAP
        WRITE(NFIL,*)TRIM(NAME(MAP(1,IMAP)))//DISTID(MAP(2,IMAP),MAP(3,IMAP),MAP(4,IMAP))//' '
      ENDDO
      WRITE(NFIL,FMT='(A)')'"'
      WRITE(NFIL,FMT='(A)')-'ELEMENT'//+'T'//-'YPE="'
      WRITE(NFIL,FMT='(10(A2," "))')EL(MAP(1,:))
!      WRITE(NFIL,*)(EL(MAP(1,I))//' ',I=1,SIZE(MAP(1,:))) !THIS DID NOT WORK
      WRITE(NFIL,FMT='(A)')'"'
      WRITE(NFIL,FMT='(A)')-'X'//+'F'//-'RACT="'
      WRITE(NFIL,FMT='(10F10.5)')X(1,MAP(1,:))+REAL(MAP(2,:))
      WRITE(NFIL,FMT='(A)')'"'
      WRITE(NFIL,FMT='(A)')-'Y'//+'F'//-'RACT="'
      WRITE(NFIL,FMT='(10F10.5)')X(2,MAP(1,:))+REAL(MAP(3,:))
      WRITE(NFIL,FMT='(A)')'"'
      WRITE(NFIL,FMT='(A)')-'Z'//+'F'//-'RACT="'
      WRITE(NFIL,FMT='(10F10.5)')X(3,MAP(1,:))+REAL(MAP(4,:))
      WRITE(NFIL,FMT='(A)')'"/>'
!
!     ==========================================================================
!     == CALCULATE BONDS                                                      ==
!     ==========================================================================
      NBONDX=16*NMAP
      ALLOCATE(BOND(2,NBONDX))
      NBOND=0
      DO I=1,NMAP
        R1(:)=MATMUL(RBASNEU,X(:,MAP(1,I))+REAL(MAP(2:4,I)))
        RCOV1=RCOV(MAP(1,I))
        DO J=I+1,NMAP
          R2(:)=MATMUL(RBASNEU,X(:,MAP(1,J))+REAL(MAP(2:4,J)))
          RCOV2=RCOV(MAP(1,J))
          SVAR=SQRT(SUM((R2-R1)**2))/(RCOV1+RCOV2)-1.D0
          IF(SVAR.LT.TOLBOND) THEN
            NBOND=NBOND+1
            IF(NBOND.GT.NBONDX) THEN
              CALL ERROR$MSG('ARRAY BOUNDS EXCEEDED')
              CALL ERROR$STOP('WRITEAVOGADRO')
            END IF
            BOND(1,NBOND)=I
            BOND(2,NBOND)=J
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WRITE BONDS                                                          ==
!     ==========================================================================
      WRITE(NFIL,FMT='(A)')-'<BOND'//+'A'//-'RRAY'
      WRITE(NFIL,FMT='(A)')-'ATOM'//+'R'//-'EF1="'
!      WRITE(NFIL,*)(TRIM(NAME(MAP(1,BOND(1,I))))//DISTID(MAP(2,BOND(1,I)),MAP(3,BOND(1,I)),MAP(4,BOND(1,I)))//' ',I=1,NBOND)
      DO I=1,NBOND
        WRITE(NFIL,*)TRIM(NAME(MAP(1,BOND(1,I))))//DISTID(MAP(2,BOND(1,I)),MAP(3,BOND(1,I)),MAP(4,BOND(1,I)))//' '
      ENDDO
      WRITE(NFIL,FMT='(A)')'"'
      WRITE(NFIL,FMT='(A)')-'ATOM'//+'R'//-'EF2="'
!      WRITE(NFIL,*)(TRIM(NAME(MAP(1,BOND(2,I))))//DISTID(MAP(2,BOND(2,I)),MAP(3,BOND(2,I)),MAP(4,BOND(2,I)))//' ',I=1,NBOND)
      DO I=1,NBOND
        WRITE(NFIL,*)TRIM(NAME(MAP(1,BOND(2,I))))//DISTID(MAP(2,BOND(2,I)),MAP(3,BOND(2,I)),MAP(4,BOND(2,I)))//' '
      ENDDO
      WRITE(NFIL,FMT='(A)')'"'
      WRITE(NFIL,FMT='(A)')-'ORDER="'
      WRITE(NFIL,FMT='(10I7)')(1,I=1,NBOND)
      WRITE(NFIL,FMT='(A)')'"/>'
      WRITE(NFIL,FMT='(A)')-'</MOLECULE>'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
        REAL(8), INTENT(IN) :: R1(3)
        REAL(8), INTENT(IN) :: R2(3)
        REAL(8)             :: R3(3)
        R3(1)=R1(2)*R2(3)-R1(3)*R2(2)
        R3(2)=R1(3)*R2(1)-R1(1)*R2(3)
        R3(3)=R1(1)*R2(2)-R1(2)*R2(1)
      END FUNCTION DYADISCHES_PRODUCT




