!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      SAVE
      TYPE(LL_TYPE)   :: LL_CNTL
      INTEGER(4)      :: NFILO,NFIL
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
!     **************************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT

      CALL TRACE$SETL4('ON',.FALSE.)
      CALL TRACE$PUSH('MAIN')
 !
!     ==========================================================================
!     ==  RESOLVE ARGUMENTLIST AND INITIALIZE FILE HANDLER                    ==
!     ==========================================================================
      CALL INITIALIZEFILEHANDLER
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('DPCNTL',NFIL)
!      CALL FILEHANDLER$REPORT(6,'ALL')
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)

      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO) 
        CALL LINKEDLIST$REPORT(LL_CNTL,NFILO)
      END IF
!
!     ==========================================================================
!     ==  ANALYZE CONTROL FILE                                                ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==========================================================================
!     ==  WRITE HEADER                                                        ==
!     ==========================================================================
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"           DOS PLOT TOOL                ")')
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"    FOR THE PROJECTOR-AUGMENTED WAVE METHOD  ")')
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(T10 &
     &           ,"P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY")')
      WRITE(NFILO,FMT='(T10 &
     &            ,"DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3")')
      WRITE(NFILO,*)
!
!     ==========================================================================
!     ==  RESOLVE CONTROL FILE: RESULT IS ENCODED IN DOSSET_MODULE
!     ==========================================================================
      CALL READCNTL(LL_CNTL)
!
!     ==========================================================================
!     ==  RESOLVE CONTROL FILE: RESULT IS ENCODED IN DOSSET_MODULE
!     ==========================================================================
      CALL CUMMULATEDDOS()
      
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DPCNTL')
      CALL LINKEDLIST$REPORT_UNUSED(LL_CNTL,NFILO)
      CALL TRACE$POP()
      CALL ERROR$NORMALSTOP()
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
MODULE DOSSETS_MODULE
TYPE DOSSET_TYPE
  INTEGER(4)    :: IG
  INTEGER(4)    :: IS
  CHARACTER(256) :: PREFIX
  CHARACTER(64) :: ID
  CHARACTER(32) :: FILLCOLOR
  LOGICAL(4)    :: STACK
  REAL(8)       :: SCALE
  REAL(8)       :: XZERO
  LOGICAL(4)    :: LEGEND
  TYPE(DOSSET_TYPE),POINTER :: NEXT
END TYPE DOSSET_TYPE
TYPE(DOSSET_TYPE)  ,TARGET  :: FIRSTDOSSET
TYPE(DOSSET_TYPE)  ,POINTER :: DOSSET
REAL(8), ALLOCATABLE  :: XMIN(:) !(NGRAPH)
REAL(8), ALLOCATABLE  :: XMAX(:) !(NGRAPH)
REAL(8), ALLOCATABLE  :: YMIN(:) !(NGRAPH)
REAL(8), ALLOCATABLE  :: YMAX(:) !(NGRAPH)
END MODULE DOSSETS_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL(LL_CNTL)
!     **************************************************************************
!     ** RESOLVE CONTROLE FILE                                                **
!     **************************************************************************
      USE DOSSETS_MODULE , ONLY : FIRSTDOSSET,DOSSET &
     &                           ,XMIN,XMAX,YMIN,YMAX
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)   :: LL_CNTL
      LOGICAL(4)      :: TCHK
      INTEGER(4)      :: NGRAPH
      INTEGER(4)      :: NSET
      INTEGER(4)      :: IG,IS,IS1,IS2,I
      LOGICAL(4),ALLOCATABLE :: TARR(:)
      INTEGER(4),ALLOCATABLE :: ISARR(:)
      CHARACTER(256)  :: PREFIX
      REAL(8)         :: SCALE
      REAL(8)         :: XZEROGLOB
      REAL(8)         :: XZERO
      REAL(8)         :: YMINGLOB
      REAL(8)         :: YMAXGLOB
      REAL(8)         :: XMINGLOB
      REAL(8)         :: XMAXGLOB
!     **************************************************************************
                                     CALL TRACE$PUSH('READCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DPCNTL')
!
!     ==========================================================================
!     == READ FILES                                                           ==
!     ==========================================================================
                                     CALL TRACE$PASS('READ FILES')
      XZEROGLOB=0.D0
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EZERO[EV]',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'EZERO[EV]',0,XZEROGLOB)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMIN[EV]',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!DPCNTL:EMIN[EV] IS MANDATORY INPUT')
        CALL ERROR$STOP('READCNTL')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'EMIN[EV]',0,XMINGLOB)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EMAX[EV]',1,TCHK)
      IF(.NOT.TCHK) THEN 
        CALL ERROR$MSG('!DPCNTL:EMAX[EV] IS MANDATORY INPUT')
        CALL ERROR$STOP('READCNTL')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'EMAX[EV]',0,XMAXGLOB)
!
      IF(XMAXGLOB.LE.XMINGLOB) THEN
        CALL ERROR$MSG('!DPCNTL:EMAX[EV] IS BELOW EMIN[EV]')
        CALL ERROR$R8VAL('XMINGLOB',XMINGLOB)
        CALL ERROR$R8VAL('XMAXGLOB',YMAXGLOB)
        CALL ERROR$STOP('READCNTL')
      END IF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'YMIN',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!DPCNTL:YMIN IS MANDATORY INPUT')
        CALL ERROR$STOP('READCNTL')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'YMIN',0,YMINGLOB)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'YMAX',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!DPCNTL:YMAX IS MANDATORY INPUT')
        CALL ERROR$STOP('READCNTL')
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'YMAX',0,YMAXGLOB)
!
!     ==========================================================================
!     == READ GRAPHS AND SETS                                                 ==
!     ==========================================================================
                                     CALL TRACE$PASS('READ GRAPHS')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'GRAPH',NGRAPH)
      ALLOCATE(YMIN(NGRAPH))
      ALLOCATE(YMAX(NGRAPH))
      ALLOCATE(XMIN(NGRAPH))
      ALLOCATE(XMAX(NGRAPH))
      XMIN=XMINGLOB
      XMAX=XMAXGLOB
      DO IG=1,NGRAPH
        CALL LINKEDLIST$SELECT(LL_CNTL,'GRAPH',IG)
!
!       == GET DEFAULT PREFIX FOR THIS GRAPH ===================================
        PREFIX=' '
        CALL LINKEDLIST$EXISTD(LL_CNTL,'PREFIX',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'PREFIX',0,PREFIX)
!
!       == GET DEFAULT SCALE FOR THIS GRAPH ===================================
        SCALE=1.D0
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SCALE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'SCALE',0,SCALE)
!
!       == GET DEFAULT SCALE FOR THIS GRAPH ===================================
        XZERO=XZEROGLOB
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EZERO[EV]',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'EZERO[EV]',0,XZERO)
!
        YMIN(IG)=YMINGLOB
        CALL LINKEDLIST$EXISTD(LL_CNTL,'YMIN',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'YMIN',0,YMIN(IG))
!
        YMAX(IG)=YMAXGLOB
        CALL LINKEDLIST$EXISTD(LL_CNTL,'YMAX',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'YMAX',0,YMAX(IG))

        IF(YMAX(IG).LE.YMIN(IG)) THEN
          CALL ERROR$MSG('YMAX IS SMALLER THAN YMIN')
          CALL ERROR$I4VAL('IGRAPH',IG)
          CALL ERROR$R8VAL('YMIN(IG)',YMIN(IG))
          CALL ERROR$R8VAL('YMAX(IG)',YMAX(IG))
          CALL ERROR$R8VAL('YMINGLOB',YMINGLOB)
          CALL ERROR$R8VAL('YMAXGLOB',YMAXGLOB)
          CALL ERROR$STOP('READCNTL')
        END IF
!
!       ========================================================================
!       == LOOP OVER SETS                                                     ==
!       ========================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'SET',NSET)
!
!       == ORDERING TO THAT STACKED ARRAYS APPEAR IN REVERSE ORDER =============
        ALLOCATE(TARR(NSET))
        TARR=.FALSE.
        DO IS=1,NSET
          CALL LINKEDLIST$SELECT(LL_CNTL,'SET',IS)
          CALL LINKEDLIST$EXISTD(LL_CNTL,'STACK',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'STACK',1,TARR(IS))
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO      
        ALLOCATE(ISARR(NSET))
        DO IS=1,NSET
          ISARR(IS)=IS
        ENDDO
!
        IS=1
        DO WHILE (IS.LT.NSET)
          IF(TARR(IS+1)) THEN
            IS1=IS
            IS2=IS+1
            DO I=IS1+2,NSET
              IF(.NOT.TARR(I)) EXIT
              IS2=I
            ENDDO
            DO I=IS1,IS2
              ISARR(I)=IS2+IS1-I
            ENDDO
            IS=IS2+1
          ELSE
            IS=IS+1
          END IF
        ENDDO
        DEALLOCATE(TARR)
!
!       == COLLECT DATA IN THE NEW ORDER =======================================
        DO IS=1,NSET
          CALL LINKEDLIST$SELECT(LL_CNTL,'SET',ISARR(IS))
          IF(ASSOCIATED(DOSSET)) THEN
            ALLOCATE(DOSSET%NEXT)
            DOSSET=>DOSSET%NEXT
          ELSE
            DOSSET => FIRSTDOSSET
          END IF
          DOSSET%IG=IG
          DOSSET%IS=IS
          DOSSET%FILLCOLOR=' '
          DOSSET%STACK=.FALSE.
          DOSSET%LEGEND=.FALSE.
          NULLIFY(DOSSET%NEXT)

!         == GET PREFIX ========================================================
          DOSSET%PREFIX=PREFIX
          CALL LINKEDLIST$EXISTD(LL_CNTL,'PREFIX',1,TCHK)
          IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'PREFIX',1,DOSSET%PREFIX)
!
!         == GET ID ============================================================
          CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!DPCNTL!GRAPH!SET:ID IS MANDATORY INPUT')
            CALL ERROR$STOP('READCNTL')
          END IF
          CALL LINKEDLIST$GET(LL_CNTL,'ID',1,DOSSET%ID)
!
!         == GET STACK === =====================================================
          DOSSET%STACK=.FALSE.
          CALL LINKEDLIST$EXISTD(LL_CNTL,'STACK',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'STACK',1,DOSSET%STACK)
!
!         == GET LEGEND === ====================================================
          DOSSET%LEGEND=.FALSE.
          CALL LINKEDLIST$EXISTD(LL_CNTL,'LEGEND',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'LEGEND',1,DOSSET%LEGEND)
!
!         == GET SCALE =========================================================
          DOSSET%SCALE=SCALE
          CALL LINKEDLIST$EXISTD(LL_CNTL,'SCALE',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'SCALE',1,DOSSET%SCALE)
!
!         == GET XZERO =========================================================
          DOSSET%XZERO=XZERO
          CALL LINKEDLIST$EXISTD(LL_CNTL,'EZERO[EV]',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'EZERO[EV]',1,DOSSET%XZERO)
!
!         == GET FILLCOLOR =====================================================
          CALL LINKEDLIST$EXISTD(LL_CNTL,'COLOR',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FILLCOLOR',1,TCHK)
            IF(.NOT.TCHK) THEN
              DOSSET%FILLCOLOR=' '
            ELSE
              CALL LINKEDLIST$GET(LL_CNTL,'FILLCOLOR',1,DOSSET%FILLCOLOR)
            END IF
          ELSE
            CALL LINKEDLIST$GET(LL_CNTL,'COLOR',1,DOSSET%FILLCOLOR)
          END IF
!
          CALL LINKEDLIST$SELECT(LL_CNTL,'..') ! LEAVE SET
        ENDDO ! END FOR IS
        DEALLOCATE(ISARR)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')  !LEAVE GRAPH
      ENDDO ! END FOR IG
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CUMMULATEDDOS()
      USE STRINGS_MODULE
      USE DOSSETS_MODULE , ONLY : FIRSTDOSSET,DOSSET,DOSSET_TYPE &
     &                            ,XMIN,XMAX,YMIN,YMAX
      IMPLICIT NONE
      TYPE(DOSSET_TYPE),POINTER:: DOSSET1
      INTEGER(4)               :: NGRAPH
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: IG,IS,IS1,IS2,I
      REAL(8)                  :: LW !LINEWIDTH
      REAL(8)                  :: TICKSPACING
      REAL(8)                  :: SVAR
      CHARACTER(32)            :: FC
      CHARACTER(32)            :: LC
      CHARACTER(256)           :: STRING
      INTEGER(4),ALLOCATABLE   :: SETSINGRAPH(:)
!     *************************************************************************
                                     CALL TRACE$PUSH('CUMULATEDDOS')
!
!     =========================================================================
!     == START XMGRACE FILE, SET NGRAPH                                      ==
!     =========================================================================
                                     CALL TRACE$PASS('START XMGRACE FILE')
      NGRAPH=0
      DOSSET=>FIRSTDOSSET
      DO 
!        PRINT*,DOSSET%IG,DOSSET%IS,DOSSET%SCALE,DOSSET%STACK,TRIM(DOSSET%ID)
        NGRAPH=MAX(NGRAPH,DOSSET%IG)
        IF(.NOT.ASSOCIATED(DOSSET%NEXT)) EXIT
        DOSSET=>DOSSET%NEXT
      ENDDO
      CALL FILEHANDLER$UNIT('GRACEBATCH',NFIL)
      CALL GRACE_INITIALIZE(NFIL)
      CALL GRACE_ARRANGE(NFIL,NGRAPH)
!
!     =========================================================================
!     == READ DOS DATA FILES                                                 ==
!     =========================================================================
                                     CALL TRACE$PASS('READ DOS DATA FILES')
      DOSSET=>FIRSTDOSSET
      DO 
        STRING=TRIM(DOSSET%PREFIX)//TRIM(DOSSET%ID)//-'.DOS'
        CALL GRACE_READ(NFIL,DOSSET%IG-1,STRING)
        DO IS=2*(DOSSET%IS-1),2*(DOSSET%IS-1)+1  !IS= 0,1, 2,3, 4,5, ...
!         == SCALE SET =========================================================
          CALL GRACE_SCALESET(NFIL,DOSSET%IG-1,IS,DOSSET%SCALE)
!         == POSITION ZERO AT XZERO ============================================
          CALL GRACE_SHIFTSET(NFIL,DOSSET%IG-1,IS,-DOSSET%XZERO,0.D0)
        ENDDO
        IF(.NOT.ASSOCIATED(DOSSET%NEXT)) EXIT
        DOSSET=>DOSSET%NEXT
      ENDDO
!
!     =========================================================================
!     == READ SETTINGS                                                       ==
!     =========================================================================
                                     CALL TRACE$PASS('READ SETTINGS')
      ALLOCATE(SETSINGRAPH(NGRAPH))
      DOSSET=>FIRSTDOSSET
      IG=0
      DO 
        IF(DOSSET%IG.GT.IG) THEN
!         == NEW GRAPH
          IF(IG.NE.0)SETSINGRAPH(IG)=IS
          IG=DOSSET%IG
          CALL GRACE_WORLD(NFIL,IG-1,XMIN(IG),XMAX(IG),YMIN(IG),YMAX(IG))
          CALL GRACE_SETZEROAXIS(NFIL,IG-1,.FALSE.,.FALSE.)
          IS=0
        END IF
!
!       == ADD STACK TO THIS LINE ==============================================
        IF(DOSSET%STACK) THEN
          DOSSET1=>DOSSET
          DO WHILE (ASSOCIATED(DOSSET1%NEXT))
            DOSSET1=>DOSSET1%NEXT
            IF(DOSSET1%IG.NE.DOSSET%IG) EXIT
            IS1=2*(DOSSET%IS-1)
            IS2=2*(DOSSET1%IS-1)
            DO I=0,1
              CALL GRACE_ADDSET(NFIL,IG-1,IS1+I,IG-1,IS2+I,1.D0)
            ENDDO
            IF(.NOT.DOSSET1%STACK) EXIT
          ENDDO
        END IF
!
!       == DRAW LINES ==========================================================
        LW=0.5D0 !LINEWIDTH
        LC='BLACK'   ! LINECOLOR
        FC=DOSSET%FILLCOLOR
!       __FILLED AND EMPTY______________________________________________________
        CALL GRACE_LINE(NFIL,IG-1,IS+0,LW,LC,'L'//TRIM(FC))
!       __FILLED ONLY___________________________________________________________
        CALL GRACE_LINE(NFIL,IG-1,IS+1,LW,LC,FC)
        IS=IS+2
        IF(.NOT.ASSOCIATED(DOSSET%NEXT)) EXIT
        DOSSET=>DOSSET%NEXT
      ENDDO
      IF(IG.NE.0)SETSINGRAPH(IG)=IS
!     
!     ==========================================================================
!     == WRITE LEGENDS                                                        ==
!     ==========================================================================
      DOSSET=>FIRSTDOSSET
      IG=0
      DO 
        IF(DOSSET%IG.GT.IG) THEN
!         == NEW GRAPH
          IG=DOSSET%IG
          !CALL GRACE_WORLD(NFIL,IG-1,XMIN(IG),XMAX(IG),YMIN(IG),YMAX(IG))
          !CALL GRACE_SETZEROAXIS(NFIL,IG-1,.FALSE.,.FALSE.)
          IS=0
        END IF
        IF(DOSSET%LEGEND)THEN
          FC=DOSSET%FILLCOLOR
!         __FILLED ONLY_________________________________________________________
          CALL GRACE_COPYSET(NFIL,IG-1,IS+SETSINGRAPH(IG),IG-1,IS+1)
          CALL GRACE_SHIFTSET(NFIL,IG-1,IS+SETSINGRAPH(IG),1000.0D0,1000.0D0)
          CALL GRACE_LINE(NFIL,IG-1,IS+SETSINGRAPH(IG),LW*10.0D0,FC,FC)
          CALL GRACE_LEGEND(NFIL,IG-1,IS+SETSINGRAPH(IG) &
     &                               ,TRIM(DOSSET%ID)//' OCC')
!         __FILLED AND EMPTY____________________________________________________
          CALL GRACE_COPYSET(NFIL,IG-1,IS+SETSINGRAPH(IG)+1,IG-1,IS+0)
          CALL GRACE_SHIFTSET(NFIL,IG-1,IS+SETSINGRAPH(IG)+1,1000.0D0,1000.0D0)
          CALL GRACE_LINE(NFIL,IG-1,IS+SETSINGRAPH(IG)+1 &
     &                             ,LW*10.0D0,'L'//TRIM(FC),'L'//TRIM(FC))
          CALL GRACE_LEGEND(NFIL,IG-1,IS+SETSINGRAPH(IG)+1 &
     &                               ,TRIM(DOSSET%ID)//' UNOCC')
        ENDIF
        IS=IS+2
        
        IF(.NOT.ASSOCIATED(DOSSET%NEXT)) EXIT
        DOSSET=>DOSSET%NEXT
      ENDDO

!
!     ==========================================================================
!     == WRITE AXES AND TICKS                                                 ==
!     ==========================================================================
                                     CALL TRACE$PASS('WRITE AXES AND TICKS')
      DO IG=0,NGRAPH-1
        WRITE(STRING,*)IG
        WRITE(NFIL,FMT='("FOCUS G",A)')TRIM(ADJUSTL(STRING))
        WRITE(NFIL,FMT='("XAXIS TICK ON")')
!       == DISTANCE OF MAJOR TICKS (WITH LABELS) ==============================
!       == CALCULATE SUITABLE TICKSPACING FROM THE INTERVAL BOUNDS
        SVAR=LOG10((XMAX(IG+1)-XMIN(IG+1))/5.D0)   
        TICKSPACING=10.D0**FLOOR(SVAR)
        SVAR=MODULO(SVAR,1.D0)
        IF(SVAR.LT.LOG10(2.D0)) THEN
          SVAR=1.D0
        ELSE IF(SVAR.GT.LOG10(5.D0)) THEN
          SVAR=5.D0
        ELSE
          SVAR=2.D0
        END IF
        TICKSPACING=TICKSPACING*SVAR
        WRITE(NFIL,FMT='("XAXIS TICK MAJOR ",F20.5)')TICKSPACING
!       == LENGTH OF TICKS
        WRITE(NFIL,FMT='("XAXIS TICK MAJOR SIZE ",F10.3)')0.3D0
        WRITE(NFIL,FMT='("XAXIS TICKLABEL OFF")')
!       == SWITCH TICKS ON THE Y AXIS OFF
        WRITE(NFIL,FMT='("YAXIS TICK OFF")')
!       == SWITCH TICK LABELS ON THE Y AXIS OFF
        WRITE(NFIL,FMT='("YAXIS TICKLABEL OFF")')
!       == SET LEGEND POSITION         
        CALL GRACE_LEGEND_POSITION(NFIL,IG,0.0D0,0.9D0-(IG)/REAL(NGRAPH,KIND=8)*0.8D0)
      ENDDO
      WRITE(NFIL,FMT='("XAXIS TICKLABEL ON")')
      WRITE(NFIL,FMT='("XAXIS LABEL ",A)')+'"E('//-'E'//+'V)"'
!
!     =========================================================================
!     == EXECUTE XMGRACE                                                     ==
!     =========================================================================
                                     CALL TRACE$PASS('BEFORE EXECUTE')
      CALL GRACE_EXECUTE(NFIL)
      CALL FILEHANDLER$CLOSE('GRACEBATCH')
                                     CALL TRACE$POP()
      RETURN
      END
!      
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE INITIALIZEFILEHANDLER
      USE STRINGS_MODULE
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: PDOSINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
!     ******************************************************************
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE OF THE GRACE TOOL IS MANDATORY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,PDOSINNAME)
      ISVAR=INDEX(PDOSINNAME,-'.DPCNTL',BACK=.TRUE.)
      IF(ISVAR.NE.0) THEN
        ROOTNAME=PDOSINNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('DPCNTL',.FALSE.,PDOSINNAME)
      RETURN
      END SUBROUTINE INITIALIZEFILEHANDLER
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STANDARDFILES
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: ID
!     **************************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!  
!     ==================================================================
!     == SET STANDARD FILENAMES                                       ==
!     ==================================================================
!
!     ==  ERROR FILE ===================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DPERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOL FILE ================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DPPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =============================================
      ID=+'DPCNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.DPCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  COMMAND STACK FOR XMGRACE =========================================
      ID=+'GRACEBATCH'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.BAT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  EPS GRAPHICS FILE   =========================================
      ID=+'EPS'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.EPS')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PDF GRAPHICS FILE   =========================================
      ID=+'PDF'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PDF')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END SUBROUTINE STANDARDFILES
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!******                                                                  *******
!******                      GRACEBIB LIBRARY                            *******
!******                                                                  *******
!******                 FORTRAN INTERFACE TO GRACE.                      *******
!******                                                                  *******
!*******************************************************************************
!*******************************************************************************
!**  LIBRARY TO ACCESS XMGRACE FROM WITHIN FORTRAN                            **
!**
!**  GRACE_ARRANGE(NFIL,NGRAPHS)
!**  GRACE_WORLD(NFIL,IGRAPH,XMIN,XMAX,YMIN,YMAX)
!**  GRACE_VIEW(NFIL,XMIN,XMAX,YMIN,YMAX)
!**  GRACE_READ(NFIL,IGRAPH,FILE)
!**  GRACE_SCALESET(NFIL,IGRAPH,ISET,FACTOR)
!**  GRACE_HIDESET(NFIL,IGRAPH,ISET)
!**  GRACE_LINE(NFIL,IGRAPH,ISET,LINEWIDTH,LINECOLOR,FILLCOLOR)
!**  GRACE_OUTFILE(NFIL,FILE)
!**  GRACE_EXECUTE(NFIL,TSHOW)
!**  
!*******************************************************************************
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_ITOSET(IGRAPH,ISET,SET)
!      *************************************************************************
!      ** CONVERTS (IGRAPH,ISET) INTO THE IDENTIFIER SET SUCH AS 'G0.S1'      **
!      *************************************************************************
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN)  :: IGRAPH
       INTEGER(4)  ,INTENT(IN)  :: ISET
       CHARACTER(*),INTENT(OUT) :: SET
       CHARACTER(16)            :: STRING
!      *************************************************************************
       WRITE(STRING,*)IGRAPH
       SET='G'//ADJUSTL(STRING)
       WRITE(STRING,*)ISET
       SET=TRIM(SET)//'.S'//ADJUSTL(STRING)
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_MAPCOLOR(NFIL,I,R,G,B,NAME)
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: I
       INTEGER(4)  ,INTENT(IN) :: R,G,B
       CHARACTER(*),INTENT(IN) :: NAME
!      *************************************************************************
       WRITE(NFIL,FMT='("MAP COLOR ",I5," TO (",I3,", ",I3,", ",I3,"), ",A)') &
      &      I,R,G,B,'"'//-TRIM(ADJUSTL(NAME))//'"'
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_INITIALIZE(NFIL)
!      *************************************************************************
!
!      *************************************************************************
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)    ,INTENT(IN) :: NFIL
       CHARACTER(512)            :: PALETTE=' '
       INTEGER(4)                :: I
       CHARACTER(64)             :: FMT
!      *************************************************************************
       I=0
       CALL GRACE_MAPCOLOR(NFIL,I,255,255,255,'WHITE') ; I=I+1
       CALL GRACE_MAPCOLOR(NFIL,I,  0,  0,  0,'BLACK') ; I=I+1
!
       IF(PALETTE.EQ.'THANKYOU!') THEN
         CALL GRACE_MAPCOLOR(NFIL,I,235, 24, 83,'RED')    ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,140, 12,108,'BLUE')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,235,227, 24,'YELLOW') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I, 24,235,156,'GREEN')  ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,235,167, 24,'ORANGE') ; I=I+1
!
       ELSE IF(PALETTE.EQ.'COLORCOMBO178') THEN
!        == HTTP://WWW.COLORCOMBOS.COM/COLOR-SCHEMES/178/COLORCOMBO178.HTML
         CALL GRACE_MAPCOLOR(NFIL,I,239, 89,123,'RED')    ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,109, 49,'ORANGE') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,115,182,107,'GREEN')  ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,203, 24,'YELLOW') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I, 41,162,198,'BLUE')   ; I=I+1
!
       ELSE IF(PALETTE.EQ.'COLORCOMBO168') THEN
!        == HTTP://WWW.COLORCOMBOS.COM/COLOR-SCHEMES/168/COLORCOMBO168.HTML
         CALL GRACE_MAPCOLOR(NFIL,I,  0,158,206,'BLUE')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,215,  0,'ORANGE') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,247,215,  8,'YELLOW') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,206,  0,  0,'RED')    ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,256,207, 49,'GREEN')   ; I=I+1
!      
       ELSE
         CALL GRACE_MAPCOLOR(NFIL,I,255,  0,  0,'RED')    ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,150,150,'LRED')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,  0,230,  0,'GREEN')  ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,150,255,150,'LGREEN') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,  0,  0,255,'BLUE')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,150,150,255,'LBLUE')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,230,230,  0,'YELLOW') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,255,150,'LYELLOW') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,  0,230,230,'CYAN')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,150,255,255,'LCYAN')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,  0,255,'MAGENTA'); I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,150,255,'LMAGENTA'); I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,165,  0,'ORANGE') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,210,120,'LORANGE') ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,200,150,150,'BROWN')  ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,255,180,180,'LBROWN')  ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,170,130,  0,'GOLD')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,210,180,100,'LGOLD')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,  0,155,  0,'DARKGREEN')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,120,220,120,'LDARKGREEN')   ; I=I+1
         CALL GRACE_MAPCOLOR(NFIL,I,150,150,150,'LBLACK')   ; I=I+1
       CALL GRACE_MAPCOLOR(NFIL,I,240,240,240,'LGREY')   ; I=I+1
       END IF
       CALL GRACE_MAPCOLOR(NFIL,I,200,200,200,'GREY')
!
!      =========================================================================
!      == SWITCH BACKGROUND FILLING OFF                                       ==
!      =========================================================================
       FMT=-'("PAGE BACKGROUND FILL OFF")'
       WRITE(NFIL,FMT=FMT)
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_EXECUTE(NFIL)
!      *************************************************************************
!
!      *************************************************************************
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL  ! FILE UNIT OF THE BATCHFILE
       CHARACTER(250)          :: GRBFILE
       CHARACTER(250)          :: ROOT  
       INTEGER(4)              :: IPOS
!      *************************************************************************
                                     CALL TRACE$PUSH('GRACE_EXECUTE')
       CALL FILEHANDLER$UNIT('GRACEBATCH',NFIL)
       CALL FILEHANDLER$FILENAME('GRACEBATCH',GRBFILE)
       IPOS=INDEX(GRBFILE,-'.BAT',.TRUE.)
       ROOT=GRBFILE(1:IPOS-1)
       PRINT*,'EXECUTE WITH THE FOLLOWING COMMANDS:'
       PRINT*,-'XMGRACE -NOSAFE -BATCH ',TRIM(GRBFILE)
       PRINT*,-'GRACEBAT -NOSAFE -HDEVICE '//+'EPS'//-' -BATCH ',TRIM(GRBFILE) &
      &      ,-' -HARDCOPY -PRINTFILE ',TRIM(ROOT)//-'_DOS.EPS'
       PRINT*,-'GRACEBAT -NOSAFE -HDEVICE '//+'PDF'//-' -BATCH ',TRIM(GRBFILE) &
      &         ,-' -HARDCOPY -PRINTFILE ',TRIM(ROOT)//-'_DOS.PDF'
       PRINT*,-'GRACEBAT -NOSAFE -HDEVICE '//+'SVG'//-' -BATCH ',TRIM(GRBFILE) &
      &         ,-' -HARDCOPY -PRINTFILE ',TRIM(ROOT)//-'_DOS.SVG'
!
!      =========================================================================
!      == CALL SYSTEM TO PRODUCE A PLOT FILE                                  ==
!      =========================================================================
!!$       WRITE(NFIL,FMT='("PRINT TO ",A)')'"'//TRIM(ADJUSTL(PLOTFILE))//'"'
!!$       WRITE(NFIL,FMT='("HARDCOPY DEVICE ",A)')'"'//+TRIM(ADJUSTL(TYPE))//'"'
!!$       WRITE(NFIL,FMT='("PRINT")')
!!$       COMMAND=-'GRACEBAT -BATCH '//TRIM(GRBFILE)//-' -NOASK -NOSAFE'
!!$       COMMAND=TRIM(COMMAND)//-' -HARDCOPY'
!!$       CALL SYSTEM(COMMAND)
                                     CALL TRACE$POP()
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_ARRANGE(NFIL,NGRAPHS)
!      *************************************************************************
!      **  ARRANGE(NROWS, NCOLS, OFFSET, HGAP, VGAP)                          **
!      **  ARRANGE(NROWS, NCOLS, OFFSET, HGAP, VGAP, HVINV, HINV, VINV)       **
!      **  ARRANGE(NROWS, NCOLS, OFFSET, HGAP, VGAP, HVINV, HINV, VINV, SNAKE)**
!      **  OFFSET IS THE DISTANCE OF THE DRAWING BOX FROM THE PAGE BOUNDARY   **
!      **  IT MUST BE SUFFICIENTLY LARGE TO ACCOMODATE THE AXIS LABELS ETC.   **
!      *************************************************************************
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: NGRAPHS
       CHARACTER(64)           :: FMT
       CHARACTER(128)          :: STRING
       CHARACTER(128)          :: GRAPH
       INTEGER(4)              :: IGRAPH
!      *************************************************************************
       FMT=-'("ARRANGE(",I5,",1,0.1,0.0,0.0)")'
       WRITE(NFIL,FMT=FMT)NGRAPHS
!
!      =========================================================================
!      =========================================================================
!      =========================================================================
       DO IGRAPH=1,NGRAPHS
         WRITE(STRING,*)IGRAPH-1
         GRAPH='G'//TRIM(ADJUSTL(STRING))
         WRITE(NFIL,FMT=-'("WITH ",A)')TRIM(GRAPH)
         FMT=-'("FRAME BACKGROUND COLOR 0")'
         WRITE(NFIL,FMT=FMT)
         FMT=-'("FRAME BACKGROUND PATTERN 1")'
         WRITE(NFIL,FMT=FMT)
       ENDDO
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_READ(NFIL,IGRAPH,FILE)
!      *************************************************************************
!      ** READ DATA FROM A X,Y1,Y2,Y3... FILE                                 **
!      **                                                                     **
!      ** THE BATCH FILE MAY CONTAIN EITHER A READ COMMAND OR THE FILE CAN BE **
!      ** INSERTED DIRECTLY. DECIDE WITH HARD-WIRED PARAMETER 'TINSERT'       **
!      ** CAUTION: THE OPTION 'TINSERT=.TRUE.' DOES NOT WORK!!!               **
!      *************************************************************************
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       CHARACTER(*),INTENT(IN) :: FILE
       LOGICAL(4)  ,PARAMETER  :: TINSERT=.FALSE.
       INTEGER(4)              :: NFILFROM
       INTEGER(4)              :: I
       LOGICAL(4)              :: TOPEN
       CHARACTER(128)          :: STRING
       CHARACTER(128)          :: GRAPH
       CHARACTER(1024)          :: LINE
!      *************************************************************************
       WRITE(STRING,*)IGRAPH
       STRING=ADJUSTL(STRING)
       GRAPH='G'//TRIM(STRING)
       WRITE(NFIL,FMT=-'("WITH ",A)')TRIM(GRAPH)
       IF(TINSERT) THEN
         WRITE(NFIL,FMT=-'("TYPE NXY")')
!        == FIND FILE UNIT FOR THE FILE TO INSERT ==============================
         NFILFROM=-1
         DO I=1000,2000
           INQUIRE(I,OPENED=TOPEN)
           IF(TOPEN) CYCLE
           NFILFROM=I
           EXIT
         ENDDO
         IF(NFILFROM.EQ.-1) THEN
           CALL ERROR$MSG('NO OPEN FILE UNIT FOUND')
           CALL ERROR$STOP('GRACE_READ')
         END IF
!        == OPEN FILE, COPY LINE-BY-LINE AND CLOSE FILE AGAIN ==================
         OPEN(UNIT=NFILFROM,FILE=TRIM(FILE))
         DO 
           READ(NFILFROM,FMT='(A)',END=1000)LINE
           IF(LEN_TRIM(LINE).EQ.0) CYCLE ! DISREGARD EMPTY LINES
           WRITE(NFIL,FMT='(A)')TRIM(LINE)
         ENDDO
1000     CONTINUE
         CLOSE(NFILFROM)
!        == ADD TERMINATION SYMBOL =============================================
         WRITE(NFIL,FMT=-'("&")')
       ELSE
         WRITE(NFIL,FMT=-'("READ NXY ",A)')'"'//TRIM(ADJUSTL(FILE))//'"'
       END IF
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_SCALESET(NFIL,IGRAPH,ISET,FACTOR)
!      *************************************************************************
!      ** READ DATA FROM A X,Y1,Y2,Y3... FILE                                 **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       INTEGER(4)  ,INTENT(IN) :: ISET
       REAL(8)     ,INTENT(IN) :: FACTOR
       CHARACTER(128)          :: STRING
       CHARACTER(128)          :: SET
       CHARACTER(128)          :: COMMAND
!      *************************************************************************
       WRITE(STRING,*)IGRAPH
       STRING=ADJUSTL(STRING)
       SET='G'//TRIM(STRING)
       WRITE(STRING,*)ISET
       STRING=ADJUSTL(STRING)
       SET=TRIM(SET)//'.S'//TRIM(STRING)//'.Y'
       COMMAND=TRIM(SET)//' = '//TRIM(SET)//' * '
       WRITE(NFIL,FMT='(A,F10.5)')TRIM(COMMAND),FACTOR
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_SHIFTSET(NFIL,IGRAPH,ISET,DX,DY)
!      *************************************************************************
!      ** READ DATA FROM A X,Y1,Y2,Y3... FILE                                 **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       INTEGER(4)  ,INTENT(IN) :: ISET
       REAL(8)     ,INTENT(IN) :: DX,DY
       CHARACTER(128)          :: STRING
       CHARACTER(128)          :: SET
       CHARACTER(128)          :: COMMAND
!      *************************************************************************
       WRITE(STRING,*)IGRAPH
       STRING=ADJUSTL(STRING)
       SET='G'//TRIM(STRING)
       WRITE(STRING,*)ISET
       STRING=ADJUSTL(STRING)
       SET=TRIM(SET)//'.S'//TRIM(STRING)
       IF(DY.GT.0.D0) THEN
         COMMAND=TRIM(SET)//'.Y'//' = '//TRIM(SET)//'.Y'//' + '
       ELSE
         COMMAND=TRIM(SET)//'.Y'//' = '//TRIM(SET)//'.Y'//' - '
       END IF
       WRITE(NFIL,FMT='(A,F10.5)')TRIM(COMMAND),ABS(DY)
       IF(DX.GT.0.D0) THEN
         COMMAND=TRIM(SET)//'.X'//' = '//TRIM(SET)//'.X'//' + '
       ELSE
         COMMAND=TRIM(SET)//'.X'//' = '//TRIM(SET)//'.X'//' - '
       END IF
       WRITE(NFIL,FMT='(A,F10.5)')TRIM(COMMAND),ABS(DX)
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_ADDSET(NFIL,IGRAPH,ISET,IGRAPH2,ISET2,FACTOR)
!      *************************************************************************
!      ** ADD FACTOR*SET2 TO SET1                                             **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       INTEGER(4)  ,INTENT(IN) :: ISET
       INTEGER(4)  ,INTENT(IN) :: IGRAPH2
       INTEGER(4)  ,INTENT(IN) :: ISET2
       REAL(8)     ,INTENT(IN) :: FACTOR
       CHARACTER(128)          :: Y1,Y2
       CHARACTER(128)          :: COMMAND
!      *************************************************************************
       CALL GRACE_ITOSET(IGRAPH,ISET,Y1)
       Y1=TRIM(Y1)//'.Y'
       CALL GRACE_ITOSET(IGRAPH2,ISET2,Y2)
       Y2=TRIM(Y2)//'.Y'
!
       COMMAND=TRIM(Y1)//' = '//TRIM(Y1)//' + '//TRIM(Y2)//' * '
       WRITE(NFIL,FMT='(A,F10.5)')TRIM(COMMAND),FACTOR
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_COPYSET(NFIL,IGRAPH,ISET,IGRAPH2,ISET2)
!      *************************************************************************
!      ** DEFINES A NEW SET ISET WITH THE SAME LENGTH AS ISET2                **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       INTEGER(4)  ,INTENT(IN) :: ISET
       INTEGER(4)  ,INTENT(IN) :: IGRAPH2
       INTEGER(4)  ,INTENT(IN) :: ISET2
       CHARACTER(128)          :: Y1,Y2
       CHARACTER(128)          :: COMMAND
!      *************************************************************************
       CALL GRACE_ITOSET(IGRAPH,ISET,Y1)
       CALL GRACE_ITOSET(IGRAPH2,ISET2,Y2)
!
       COMMAND=TRIM(Y1)//' LENGTH '//TRIM(Y2)//'.LENGTH'
       WRITE(NFIL,FMT='(A)')TRIM(COMMAND)
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_LEGEND(NFIL,IGRAPH,ISET,TEXT)
!      *************************************************************************
!      ** DEFINES A NEW SET ISET WITH THE SAME LENGTH AS ISET2                **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       INTEGER(4)  ,INTENT(IN) :: ISET
       CHARACTER(*)            :: TEXT
       CHARACTER(128)          :: Y1
       CHARACTER(128)          :: COMMAND
!      *************************************************************************
       CALL GRACE_ITOSET(IGRAPH,ISET,Y1)
!
       COMMAND=TRIM(Y1)//' LEGEND "'//TRIM(TEXT)//'"'
       WRITE(NFIL,FMT='(A)')TRIM(COMMAND)
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_LEGEND_POSITION(NFIL,IGRAPH,X,Y)
!      *************************************************************************
!      ** DEFINES A NEW SET ISET WITH THE SAME LENGTH AS ISET2                **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       REAL(8)     ,INTENT(IN) :: X,Y
       CHARACTER(128)          :: COMMAND
       CHARACTER(16)           :: STRING
!      *************************************************************************
       WRITE(STRING,*)IGRAPH
       COMMAND='WITH G'//ADJUSTL(STRING)
       WRITE(NFIL,FMT='(A)')TRIM(COMMAND)
       COMMAND='LEGEND '
       WRITE(NFIL,FMT='(A,F12.5,A,F12.5)')TRIM(COMMAND),X,',',Y
       WRITE(NFIL,FMT='(A)')'LEGEND CHAR SIZE 1.0000'
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_HIDESET(NFIL,IGRAPH,ISET)
!      *************************************************************************
!      ** READ DATA FROM A X,Y1,Y2,Y3... FILE                                 **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       INTEGER(4)  ,INTENT(IN) :: ISET
       CHARACTER(128)          :: STRING
       CHARACTER(128)          :: SET
!      *************************************************************************
       WRITE(STRING,*)IGRAPH
       STRING=ADJUSTL(STRING)
       SET='G'//TRIM(STRING)
       WRITE(STRING,*)ISET
       STRING=ADJUSTL(STRING)
       SET=TRIM(SET)//'.S'//TRIM(STRING)
       WRITE(NFIL,FMT='(A," HIDDEN TRUE")')TRIM(SET)
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_WORLD(NFIL,IGRAPH,XMIN,XMAX,YMIN,YMAX)
!      *************************************************************************
!      ** COLORS MAY BE : BLACK, GREY, RED, ORANGE, MAGENTA, VIOLET, BLUE     **
!      **                ,BROWN, GREEN4, ORANGE, TURQUOISE                    **
!      *************************************************************************
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       REAL(8)     ,INTENT(IN) :: XMIN,XMAX
       REAL(8)     ,INTENT(IN) :: YMIN,YMAX
       CHARACTER(128)          :: STRING
       CHARACTER(128)          :: GRAPH
!      *************************************************************************
       WRITE(STRING,*)IGRAPH
       STRING=ADJUSTL(STRING)
       GRAPH='G'//TRIM(STRING)
       WRITE(NFIL,FMT=-'("WITH ",A)')GRAPH
       WRITE(NFIL,FMT=-'("WORLD XMIN ",F10.3)')XMIN
       WRITE(NFIL,FMT=-'("WORLD XMAX ",F10.3)')XMAX
       WRITE(NFIL,FMT=-'("WORLD YMIN ",F10.3)')YMIN
       WRITE(NFIL,FMT=-'("WORLD YMAX ",F10.3)')YMAX
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_VIEW(NFIL,XMIN,XMAX,YMIN,YMAX)
!      *************************************************************************
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       REAL(8)     ,INTENT(IN) :: XMIN,XMAX
       REAL(8)     ,INTENT(IN) :: YMIN,YMAX
!      *************************************************************************
       WRITE(NFIL,FMT='("VIEW XMIN ",F10.3)')XMIN
       WRITE(NFIL,FMT='("VIEW XMAX ",F10.3)')XMAX
       WRITE(NFIL,FMT='("VIEW YMIN ",F10.3)')YMIN
       WRITE(NFIL,FMT='("VIEW YMAX ",F10.3)')YMAX
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_SETZEROAXIS(NFIL,IGRAPH,TXAXIS,TYAXIS)
!      *************************************************************************
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       LOGICAL(4)  ,INTENT(IN) :: TXAXIS
       LOGICAL(4)  ,INTENT(IN) :: TYAXIS
!      *************************************************************************
       IF(TXAXIS) THEN
         WRITE(NFIL,FMT='("XAXIS TYPE ZERO TRUE")')
       ELSE
         WRITE(NFIL,FMT='("XAXIS TYPE ZERO FALSE")')
       END IF
       IF(TYAXIS) THEN
         WRITE(NFIL,FMT='("YAXIS TYPE ZERO TRUE")')
       ELSE
         WRITE(NFIL,FMT='("YAXIS TYPE ZERO FALSE")')
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE GRACE_LINE(NFIL,IGRAPH,ISET,LINEWIDTH,LINECOLOR,FILLCOLOR)
!      *************************************************************************
!      ** COLORS MAY BE : NONE,BLACK, GREY, RED, ORANGE, MAGENTA, VIOLET      **
!      **                ,BLUE ,BROWN, GREEN4, ORANGE, TURQUOISE              **
!      *************************************************************************
       USE STRINGS_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: IGRAPH
       INTEGER(4)  ,INTENT(IN) :: ISET
       REAL(8)     ,INTENT(IN) :: LINEWIDTH
       CHARACTER(*),INTENT(IN) :: LINECOLOR
       CHARACTER(*),INTENT(IN) :: FILLCOLOR
       CHARACTER(64)            :: GRAPH
       CHARACTER(64)            :: SET
!      *************************************************************************
       WRITE(GRAPH,*)IGRAPH
       GRAPH='G'//TRIM(ADJUSTL(GRAPH))
       WRITE(SET,*)ISET
       SET=TRIM(GRAPH)//'.S'//TRIM(ADJUSTL(SET))
!
       WRITE(NFIL,FMT='(A," LINE TYPE 1")')TRIM(SET)
       WRITE(NFIL,FMT='(A," LINE LINEWIDTH ",F10.5)')TRIM(SET),LINEWIDTH
       WRITE(NFIL,FMT='(A," LINE COLOR ",A)')TRIM(SET) &
     &                 ,'"'//-TRIM(ADJUSTL(LINECOLOR))//'"'
       IF(FILLCOLOR.NE.' ') THEN
         WRITE(NFIL,FMT='(A," FILL TYPE 2 ")')TRIM(SET)
         WRITE(NFIL,FMT='(A," FILL RULE 0")')TRIM(SET)
         WRITE(NFIL,FMT='(A," FILL COLOR ",A)')TRIM(SET) &
     &                 ,'"'//-TRIM(ADJUSTL(FILLCOLOR))//'"'
         WRITE(NFIL,FMT='(A," FILL PATTERN 1 ")')TRIM(SET)
       END IF
       RETURN
       END
