!***********************************************************************
!***********************************************************************
!
!    do NOT run f90pp over it before compilation. instead use 
!    the following bash script
! 
!    #!/bin/bash
!    export TMP=$(mktemp).f90
!    sed 's/[$]/__/g'  src/Buildtools/F90PP/f90pp_tmplts.f90 > $TMP
!    gfortran -o f90_tmplts.x $TMP
!    rm $TMP
!
!***********************************************************************
!***********************************************************************
!
!     ..................................................................
      SUBROUTINE UPPERCASE1(OLD,NEW)
!      USE STRINGS
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: OLD
      CHARACTER(*),INTENT(OUT):: NEW
      INTEGER(4)              :: I,ISVAR
      NEW=OLD 
      DO I=1,LEN(TRIM(NEW))
        ISVAR=IACHAR(NEW(I:I))
        IF(ISVAR.GE.97.AND.ISVAR.LE.122) NEW(I:I)=ACHAR(ISVAR-32)
      ENDDO
      RETURN
      END
!.......................................................................
MODULE ERROR_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: ERROR                                                      **
!**                                                                   **
!**  PURPOSE: COLLECTS AND REPORTS ERROR MESSAGES AND                 **
!**    STOPS THE PROGRAM                                              **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    ERROR$MSG(MSG)                                                 **
!**    ERROR$R8VAL(NAME,VAL)                                          **
!**    ERROR$I4VAL(NAME,VAL)                                          **
!**    ERROR$L4VAL(NAME,VAL)                                          **
!**    ERROR$CHVAL(NAME,VAL)                                          **
!**    ERROR$STOP(ROUTINENAME)                                        **
!**    ERROR$NORMALSTOP                                               **
!**    ERROR$OVERFLOW  (SHOUD NOT BE USED)                            **
!**                                                                   **
!**  DEPENDECIES:                                                     **
!**    MPE(MPE$STOPALL,MPE$QUERY)                                     **
!**    FILEHANDLER(FILEHANDLER$UNIT,FILEHANDLER$CLOSEALL)             **
!**                                                                   **
!**  REMARKS:                                                         **
!**    THE DEPENDENCY WITH THE FILEHANDLER CAN CAUSE PROBLEMS         **
!**                                                                   **
!***********************************************************************
INTEGER(4),PARAMETER :: IMESSAGEX=50                                
CHARACTER(82)        :: MESSAGE(IMESSAGEX)                            
INTEGER(4)           :: IMESSAGE=0                                  
INTEGER(4)           :: NFILERR=0                                   
INTEGER(4)           :: NCODE=2                                     
INTEGER(4)           :: ITASK=0                                     
INTEGER(4)           :: NSTOP=0                                     
INTEGER(4)           :: NCALL=0                                     
END MODULE ERROR_MODULE                                             
!                                                                   
!     ..................................................................
      SUBROUTINE ERROR$MSG(MESSAGE_)
!     ******************************************************************
!     ** STORE ERROR MESSAGE                                          **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
!     ******************************************************************
      IMESSAGE=IMESSAGE+1
      IF(IMESSAGE.GT.IMESSAGEX) THEN
        IMESSAGE=IMESSAGE-1
        MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
      ELSE 
        MESSAGE(IMESSAGE)=MESSAGE_
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ERROR$OVERFLOW(MESSAGE_,IACTUAL,ITARGET)
!     ******************************************************************
!     ** STORE OVERFLOW MESSAGE                                       **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      INTEGER(4)  ,INTENT(IN) :: IACTUAL
      INTEGER(4)  ,INTENT(IN) :: ITARGET
!     ******************************************************************
      IMESSAGE=IMESSAGE+1
      IF(IMESSAGE.GT.IMESSAGEX-1) THEN
        IMESSAGE=IMESSAGE-1
        MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
      ELSE 
        WRITE(MESSAGE(IMESSAGE),FMT='("VARIABLE ",A," OUT OF RANGE")') &
     &         MESSAGE_
        IMESSAGE=IMESSAGE+1
        WRITE(MESSAGE(IMESSAGE),FMT='("ACTUAL VALUE ",I10' &
     &         //'," MAXIMUM ALLOWED ",I10)')IACTUAL,ITARGET
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ERROR$R8VAL(MESSAGE_,R8VAL)
!     ******************************************************************
!     ** STORE MESSAGE WITH DOUBLE PRECISION VALUE                    **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      REAL(8)     ,INTENT(IN) :: R8VAL
!     ******************************************************************
      IMESSAGE=IMESSAGE+1
      IF(IMESSAGE.GT.IMESSAGEX-1) THEN
        IMESSAGE=IMESSAGE-1
        MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
      ELSE 
        WRITE(MESSAGE(IMESSAGE),FMT='("VARIABLE ",A' &
     &        //'," HAS THE VALUE",E20.10)')MESSAGE_,R8VAL
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ERROR$I4VAL(MESSAGE_,I4VAL)
!     ******************************************************************
!     ** STORE MESSAGE WITH INTEGER VALUE                             **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      INTEGER(4)  ,INTENT(IN) :: I4VAL
!     ******************************************************************
      IMESSAGE=IMESSAGE+1
      IF(IMESSAGE.GT.IMESSAGEX-1) THEN
        IMESSAGE=IMESSAGE-1
        MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
      ELSE 
        WRITE(MESSAGE(IMESSAGE),FMT='("VARIABLE ",A' &
     &       //'," HAS THE VALUE",I10)')MESSAGE_,I4VAL
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ERROR$L4VAL(MESSAGE_,L4VAL)
!     ******************************************************************
!     ** STORE MESSAGE WITH INTEGER VALUE                             **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      LOGICAL(4)  ,INTENT(IN) :: L4VAL
!     ******************************************************************
      IMESSAGE=IMESSAGE+1
      IF(IMESSAGE.GT.IMESSAGEX-1) THEN
        IMESSAGE=IMESSAGE-1
        MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
      ELSE 
        WRITE(MESSAGE(IMESSAGE),FMT='("VARIABLE ",A' &
     &        //'," HAS THE VALUE",L7)')MESSAGE_,L4VAL
      END IF
      RETURN
      END
! 
!     ..................................................................
      SUBROUTINE ERROR$CHVAL(MESSAGE_,CHVAL)
!     ******************************************************************
!     ** STORE MESSAGE WITH CHARACTER VALUE                           **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      CHARACTER(*),INTENT(IN) :: CHVAL
      INTEGER(4)              :: MAXLEN
      INTEGER(4)              :: I1,I2
      LOGICAL(4)              :: TFIRST
!     ******************************************************************
      TFIRST=.TRUE.
      MAXLEN=LEN(TRIM(CHVAL))+LEN(TRIM(MESSAGE_))+25
      IF(MAXLEN.LE.82) THEN
        IMESSAGE=IMESSAGE+1
        WRITE(MESSAGE(IMESSAGE),FMT='("VARIABLE ",A' &
     &        //'," HAS THE VALUE ",A)')TRIM(MESSAGE_),TRIM(CHVAL)
      ELSE
        MAXLEN=LEN(TRIM(CHVAL))
        I1=1
        I2=0
        DO WHILE(I2.LT.MAXLEN)
          IMESSAGE=IMESSAGE+1
          IF(IMESSAGE.GT.IMESSAGEX-1) THEN
            IMESSAGE=IMESSAGE-1
            MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
            RETURN
          ELSE 
            IF(TFIRST) THEN
              TFIRST=.FALSE.
              I2=82-LEN(TRIM(MESSAGE_))-25
              I2=MIN(I2,MAXLEN)
              I2=MAX(0,I2)
              MESSAGE(IMESSAGE)='VARIABLE '//TRIM(MESSAGE_)//' HAS THE VALUE '//CHVAL(I1:I2)
              I1=I2+1
            ELSE
              I2=I1-1+82
              I2=MIN(I2,MAXLEN)
              MESSAGE(IMESSAGE)=CHVAL(I1:I2)
              I1=I2+1
            END IF
          END IF
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ERROR$STOP(MESSAGE_)
!     ******************************************************************
!     **  WRITE ERROR MESSAGE, FLUSH FILES AND STOP                   **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      INTEGER(4)              :: NTASKNUM=1
      INTEGER(4)              :: NTASKID=1
      INTEGER(4)              :: NFILO=6
      INTEGER(4)              :: I
!     ******************************************************************
!     __ CATCHES ADDITIONAL ERRORS OCCURING DURING STOPPING
      NCALL=NCALL+1
!
!     ==================================================================
!     == WRITE ERROR MESSAGE ON STANDARD ERROR                        ==
!     ==================================================================
      IF(NCALL.LE.2) THEN
        NFILERR=0
        DO I=1,IMESSAGE
          WRITE(NFILERR,FMT='(A)')TRIM(MESSAGE(I))
        ENDDO
        WRITE(NFILERR,FMT='("STOP IN ",A)')TRIM(MESSAGE_)
        IF(NTASKNUM.GT.1) THEN
          WRITE(NFILERR,FMT='("PROCESSOR REQUESTING STOP: ",I3' &
     &       //'," OF ",I3)')NTASKID,NTASKNUM
        END IF
      END IF
!
!     ================================================================
!     == WRITE ERROR MESSAGE ON PROTOCOL                            ==
!     ================================================================
      IF(NCALL.LE.1) THEN
        DO I=1,IMESSAGE
          WRITE(NFILO,FMT='(A)')TRIM(MESSAGE(I))
        ENDDO
        WRITE(NFILO,FMT='("STOP IN ",A)')TRIM(MESSAGE_)
        IF(NTASKNUM.GT.1) THEN
          WRITE(NFILO,FMT='("PROCESSOR REQUESTING STOP: ",I3' &
     &                  //'," OF ",I3)')NTASKID,NTASKNUM
        END IF
      END IF
!
!     ================================================================
!     == FLUSH ALL FILES                                            ==
!     ================================================================
!     CALL FLUSH_(6)
!     DO I=1000,1100
!       CALL FLUSH_(I)
!      ENDDO
!
!     ================================================================
!     == STOP                                                       ==
!     ================================================================
      STOP
      END
!
!     ..................................................................
      SUBROUTINE ERROR$NORMALSTOP
!     ******************************************************************
!     **  WRITE ERROR MESSAGE, FLUSH FILES AND STOP                   **
!     ******************************************************************
      STOP
      END
!....................................................................
MODULE TEMPLATE_MODULE
TYPE REP_TYPE
  SEQUENCE
  CHARACTER(64) :: OUT
  INTEGER(4)    :: OUTLEN
  CHARACTER(64) :: IN
  INTEGER(4)    :: INLEN
END TYPE
END MODULE TEMPLATE_MODULE
!     ..................................................................
      PROGRAM MAIN
!      USE STRINGS
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: LINELENG=256
      INTEGER(4),PARAMETER :: NLINEX=100000
      CHARACTER(LINELENG),ALLOCATABLE  :: BUFFER(:) !(NLINEX)
      INTEGER(4)           :: NLINE
      INTEGER(4)           :: I
!     ******************************************************************
      ALLOCATE(BUFFER(NLINEX))
!
!     ==================================================================
!     == READ INPUT INTO BUFFER     
!     ==================================================================
      NLINE=0
      DO I=1,NLINEX 
        READ(*,FMT='(A)',END=1000)BUFFER(I)
        NLINE=I
      ENDDO
 1000 CONTINUE
!
!     ==================================================================
!     == RESOLVE TEMPLATES    
!     ==================================================================
      I=0
      DO 
        I=I+1
        IF(I.GT.NLINE) EXIT
!       == SELECT BEGINNING OF THE TEMPLATE, I.E. A LINE '#...TEMPLATE ...'
        IF(BUFFER(I)(1:1).NE.'#')CYCLE
        CALL UPPERCASE1(BUFFER(I),BUFFER(I))
!       ___SKIP C-PREPROCESSOR DIRECTIVES_______________________________________
        IF(INDEX(BUFFER(I),'#INCLUDE').NE.0) CYCLE
        IF(INDEX(BUFFER(I),'#IFDEF').NE.0) CYCLE
        IF(INDEX(BUFFER(I),'#ELSE').NE.0) CYCLE
        IF(INDEX(BUFFER(I),'#ENDIF').NE.0) CYCLE
!
!       __KEEP ONLY THE BEGINNIG OF A NEW TEMPLATE DEFINITION___________________
!       __ #  TEMPLATE NAME_____________________________________________________
        IF(INDEX(BUFFER(I),'MODULE TEMPLATE').NE.0) CYCLE
        IF(INDEX(BUFFER(I),'END TEMPLATE').NE.0) CYCLE
!
        IF(INDEX(BUFFER(I),'TEMPLATE').EQ.0) CYCLE
        CALL RESOLVE_TEMPLATE(NLINEX,BUFFER,I,NLINE)
      ENDDO
!
!     ==================================================================
!     ==  WRITE BUFFER OUT    
!     ==================================================================
      DO I=1,NLINE
        WRITE(*,FMT='(A)')TRIM(BUFFER(I))
      ENDDO
      STOP
      END
!
!     ...............................................................  
      SUBROUTINE RESOLVE_TEMPLATE(NLINEX,BUFFER,IBEGIN,NLINE)
!     ****************************************************************
!     **  RESOLVETEMPLATE                                           **
!     **                                                            **
!     **  PURPOSE: RESOLVE TEMPLATE STRUCTURES                      **
!     **                                                            **
!     **  SYNTAX:                                                   **
!     **    THE FOLLOWING EXAMPLE PRODUCES A GENERIC INTERFACE      **
!     **    EXAMPLE WHICH TAKES A REAL(8) OR COMPLEX(8) ARGUMENT    **
!     **    WHICH CAN BE SCALAR OR ARRAY VALUED WITH RANK 1 OR 2    **
!     **                                                            **
!     **    INTERFACE EXAMPLE                                       **
!     **    # MODULE TEMPLATE NAME                                  **
!     **    END EXAMPLE                                             **
!     **    ..                                                      **
!     **    #TEMPLATE NAME                                          **
!     **    # (<TYPEID><TYPE>)([R8][REAL(8)])([C8][COMPLEX(8)])     **
!     **    # (<RANKID><RANK>)([S][])([R1][(:)])([R2][(:,:)])       **
!     **    #BODY                                                   **
!     **    !................................................       **
!     **        SUBROUTINE EXAMPLE<TYPEID><RANKID>(IN)              **
!     **        <TYPE>,INTENT(IN) :: IN<RANK>                       **   
!     **        END SUBROUTINE EXAMPLE<TYPEID><RANKID>              **
!     **    #END TEMPLATE NAME                                      **
!     **                                                            **
!     **  REMARKS:                                                  **
!     **    IT IS IMPORTANT TO INCLUD INTO THIS ROUTINE ALL         **
!     **    LINE=REPLACEMENTS,WHICH ARE NECCESARY                   **
!     **                                                            **
!     ****************************************************************
      USE TEMPLATE_MODULE
!      USE STRINGS
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)   :: NLINEX
      INTEGER(4)   ,INTENT(INOUT):: IBEGIN
      INTEGER(4)   ,INTENT(INOUT):: NLINE
      CHARACTER(*) ,INTENT(INOUT):: BUFFER(NLINEX)
      INTEGER(4)     ,PARAMETER  :: NREPX=20     ! X#(REPLACEMENTS)/INSTANCE
      INTEGER(4)     ,PARAMETER  :: NMEMX=20      
      INTEGER(4)     ,PARAMETER  :: NSETX=20  
      INTEGER(4)     ,PARAMETER  :: NINSTX=200  
      CHARACTER(LEN(BUFFER))     :: INSTANCE0
      CHARACTER(LEN(BUFFER))     :: INSTANCE(NINSTX)
      CHARACTER(128)             :: NAME         ! TEMPLATE NAME
      INTEGER(4)                 :: NREP(NSETX)  ! #(REPLACEMENTS)
      INTEGER(4)                 :: NSET
      INTEGER(4)                 :: NINST
      INTEGER(4)                 :: NMEM(NSETX)  ! #(REPLACEMENTS)
      INTEGER(4)                 :: IBODY,IEND,IINSTANCE
      TYPE(REP_TYPE),allocatable :: REP(:,:,:) !(NREPX,NMEMX,NSETX)
      CHARACTER(LEN(BUFFER))     :: LINE
      LOGICAL(4)                 :: TSET,TREP
      INTEGER(4)                 :: I1,I2,I3,I,ILINE,ILINE2
      INTEGER(4)                 :: IREP,IMEM,ISET,ISHIFT,INST
      INTEGER(4)                 :: IMEMP(NSETX)
!     ****************************************************************
      ALLOCATE(REP(NREPX,NMEMX,NSETX))
!
!     ===============================================================
!     == EXTRACT NAME AND DETERMIN END OF TEMPLATE                 ==
!     ===============================================================
      CALL UPPERCASE1(BUFFER(IBEGIN),LINE)
!     == extract template name from "#template name" ===========================
      I1=INDEX(LINE,'TEMPLATE')
      LINE=LINE(I1:)
      I1=INDEX(LINE,' ')
      LINE=LINE(I1:)
      LINE=ADJUSTL(LINE)
      I2=INDEX(LINE,' ')
      NAME=LINE(1:I2-1)
!
!     == IDENTIFY CORRESPONDING END OF TEMPLATE "#END TEMPLATE NAME" ===========
      IEND=0
      DO ILINE=IBEGIN+1,NLINE
        IF(BUFFER(ILINE)(1:1).NE.'#') CYCLE
        CALL UPPERCASE1(BUFFER(ILINE)(2:),LINE)
        I1=INDEX(LINE,'END ')
        IF(I1.EQ.0) CYCLE
        LINE=LINE(I1:)
        I1=INDEX(LINE,' ')
        LINE=ADJUSTL(LINE(I1:))
        I2=INDEX(LINE,' ')
        IF(LINE(1:I2).NE.'TEMPLATE') CYCLE
        I1=INDEX(LINE,' ')
        LINE=ADJUSTL(LINE(I1:))
        I2=INDEX(LINE,' ')
        IF(LINE(1:I2-1).EQ.TRIM(NAME)) THEN
          IEND=ILINE
          EXIT
        END IF
      ENDDO
      IF(IEND.EQ.0) THEN
        CALL ERROR$MSG('NO END TEMPLATE FOUND')
        CALL ERROR$CHVAL('TEMPLATE NAME',TRIM(NAME))
        CALL ERROR$CHVAL('FIRSTLINE',TRIM(BUFFER(IBEGIN)))
        CALL ERROR$CHVAL('LASTLINE',TRIM(BUFFER(NLINE)))
        CALL ERROR$STOP('RESOLVE_TEMPLATES')
      END IF
!
!     ===============================================================
!     == DETERMINE COMMAND REGION AND BODY OF THE PROTOTYPE        ==
!     ===============================================================
      IBODY=0
      DO ILINE=IBEGIN,IEND
        LINE=BUFFER(ILINE)
        IF(LINE(1:5).EQ.'#BODY') THEN
          IBODY=ILINE
          EXIT
        END IF
      ENDDO
      IF(IBODY.EQ.0) THEN
        CALL ERROR$MSG('NO #BODY STATEMENT FOUND')
        CALL ERROR$STOP('RESOLVE_TEMPLATES')
      END IF
!
!     ===============================================================
!     == READ REPLACEMENTS                                         ==
!     ===============================================================
      NSET=0
      NREP(:)=0
      ISET=0
      TSET=.FALSE.
      TREP=.FALSE.
      DO ILINE=IBEGIN+1,IBODY-1
        CALL UPPERCASE1(BUFFER(ILINE),LINE)
        I1=1
        CALL NEXT(LINE,I1,I2)
        DO WHILE (I1.NE.0)
!         ============================================================
!         == CHECK FOR BEGINNING OF NEW ENVIRONMENT
!         ============================================================
          IF(LINE(I1-1:I1-1).EQ.'(') THEN
            IF(TSET.OR.TREP) THEN
              CALL ERROR$MSG('SYNTAX ERROR: SET ALREADY ON')
              CALL ERROR$CHVAL('LINE',TRIM(LINE))
              CALL ERROR$STOP('RESOLVETEMPLATE')
            END IF
            IF(LINE(I1:I1).EQ.'<') THEN
              TSET=.TRUE.
              ISET=ISET+1
              IMEM=0
              NSET=ISET
              IF(ISET.GT.NSETX) THEN
                CALL ERROR$MSG('#(REPLACEMENT SETS TOO LARGE')
                CALL ERROR$CHVAL('LINE',TRIM(LINE))
                CALL ERROR$STOP('RESOLVETEMPLATE')
              END IF
            ELSE IF(LINE(I1:I1).EQ.'[') THEN
              IF(ISET.EQ.0) THEN
                CALL ERROR$MSG('#(REPLACEMENT SETS TOO LARGE')
                CALL ERROR$CHVAL('LINE',TRIM(LINE))
                CALL ERROR$STOP('RESOLVETEMPLATE')
              END IF
              TREP=.TRUE.
              IMEM=IMEM+1
              NMEM(ISET)=IMEM
              IF(IMEM.GT.NMEMX) THEN
                CALL ERROR$MSG('#(REPLACEMENT MEMBERS) TOO LARGE')
                CALL ERROR$CHVAL('LINE',TRIM(LINE))
                CALL ERROR$STOP('RESOLVETEMPLATE')
              END IF
            ELSE
              CALL ERROR$MSG('CONFUSED')
              CALL ERROR$CHVAL('LINE',TRIM(LINE))
              CALL ERROR$STOP('RESOLVETEMPLATE')
            END IF
            IREP=0
          ENDIF
!
!         ============================================================
!         == CHECK ENVIRONMENT
!         ============================================================
          IF((LINE(I1:I1).EQ.'<').AND.(.NOT.TSET)) THEN
            CALL ERROR$MSG('OUT REPLACEMENT MUST BE IN OUT ENVIRONMENT')
            CALL ERROR$CHVAL('LINE',TRIM(LINE))
            CALL ERROR$STOP('RESOLVETEMPLATE')
          END IF
          IF((LINE(I1:I1).EQ.']').AND.(.NOT.TREP)) THEN
            CALL ERROR$MSG('IN REPLACEMENT MUST BE IN IN ENVIRONMENT')
            CALL ERROR$CHVAL('LINE',TRIM(LINE))
            CALL ERROR$STOP('RESOLVETEMPLATE')
          END IF
!
!         ============================================================
!         == COPY INTO REPLACEMENT RULES             
!         ============================================================
          IF(LINE(I1:I1).EQ.'<') THEN
            IREP=IREP+1
            IF(IREP.GT.NREPX) THEN
              CALL ERROR$MSG('SYNTAX ERROR')
              CALL ERROR$CHVAL('LINE',TRIM(LINE))
              CALL ERROR$STOP('RESOLVETEMPLATE')
            END IF
            REP(IREP,1,ISET)%OUT=LINE(I1:I2)
            REP(IREP,1,ISET)%OUTLEN=I2-I1+1
          ELSE IF(LINE(I1:I1).EQ.'[') THEN
            IREP=IREP+1
            IF(IREP.GT.NREP(ISET)) THEN
              CALL ERROR$MSG('SYNTAX ERROR')
              CALL ERROR$CHVAL('LINE',TRIM(LINE))
              CALL ERROR$STOP('RESOLVETEMPLATE')
            END IF
            REP(IREP,IMEM,ISET)%IN=LINE(I1+1:I2-1)
            REP(IREP,IMEM,ISET)%INLEN=I2-I1-1
            REP(IREP,IMEM,ISET)%OUT=REP(IREP,1,ISET)%OUT
            REP(IREP,IMEM,ISET)%OUTLEN=REP(IREP,1,ISET)%OUTLEN
          ELSE
            CALL ERROR$MSG('CONFUSED')
            CALL ERROR$CHVAL('LINE',TRIM(LINE))
            CALL ERROR$STOP('RESOLVETEMPLATE')
          END IF      
!
!         ============================================================
!         == CHECK CLOSING ENVIRONMENT  
!         ============================================================
          IF(LINE(I2+1:I2+1).EQ.')') THEN
            IF(TSET) THEN
              TSET=.FALSE.
              NREP(ISET)=IREP
            ELSE IF(TREP) THEN
              IF(IREP.NE.NREP(ISET)) THEN
                CALL ERROR$MSG('NUMBER OF REPLACEMENTS MUST BE SAME IN A SET')
                CALL ERROR$STOP('RESOLVETEMPLATE')
              END IF
              TREP=.FALSE.
            END IF
          ENDIF
!
!         ============================================================
!         == SELECT NEW WORD                         
!         ============================================================
          I1=I2+1
          CALL NEXT(LINE,I1,I2)
        ENDDO
      ENDDO
      IF(TREP.OR.TSET) THEN
        CALL ERROR$MSG('ENVIRONMENT NOT PROPERLY CLOSED')
        CALL ERROR$STOP('RESOLVETEMPLATE')
      END IF
!
!     ===============================================================
!     == CREATE SPACE FOR INSTANCES                                ==
!     ===============================================================
      NINST=1
      DO ISET=1,NSET
        NINST=NINST*NMEM(ISET)
      ENDDO
      ISHIFT=(IEND-IBODY-1)*NINST+1
      IF(NLINE+ISHIFT.GT.NLINEX) THEN
        CALL ERROR$MSG('BUFFER SIZE EXCEEDED')
        CALL ERROR$STOP('RESOLVE_TEMPLATE')
      END IF
!
      DO ILINE=NLINE,IEND,-1
        BUFFER(ILINE+ISHIFT)=BUFFER(ILINE)
      ENDDO
      BUFFER(IEND)='#INSTANCES'
      IINSTANCE=IEND
      IEND=IEND+ISHIFT
      NLINE=NLINE+ISHIFT
!
!     ===============================================================
!     ==  CREATE CALLS TO INSTANCES                                ==
!     ===============================================================
      DO I=IBODY+1,IINSTANCE-1
        I1=INDEX(BUFFER(I),'SUBROUTINE')
        IF(I1.EQ.0) CYCLE
        CALL UPPERCASE1(ADJUSTL(BUFFER(I)),LINE)
        IF(INDEX(LINE,'SUBROUTINE').NE.1) CYCLE
        I1=INDEX(LINE,' ')
        LINE=ADJUSTL(LINE(I1:))
        I2=INDEX(LINE,'(')
        I3=INDEX(LINE,' ')
        IF(I3.NE.0.AND.I3.LE.I2)I2=I3
        IF(I2.EQ.0) THEN
          CALL ERROR$MSG('COULD NOT LOCATE SUBROUTINE NAME')
          CALL ERROR$STOP('RESOLVE_TEMPLATE')
        END IF
        I2=I2-1
        INSTANCE0=LINE(1:I2)
        EXIT
      ENDDO
!
!     ===============================================================
!     == CREATE INSTANCES                                          ==
!     ===============================================================
      INST=0
      IMEMP(:)=1    ! SELECTS ONE MEMBER FROM EACH SET
      ILINE2=IINSTANCE
      DO 
        INST=INST+1
        INSTANCE(INST)=INSTANCE0
        DO ISET=1,NSET
          IMEM=IMEMP(ISET)
          DO IREP=1,NREP(ISET)
            CALL REPLACE(REP(IREP,IMEM,ISET),INSTANCE(INST))
          ENDDO
        ENDDO
        INSTANCE(INST)='  MODULE PROCEDURE '//INSTANCE(INST)
        DO ILINE=IBODY+1,IINSTANCE-1
          LINE=BUFFER(ILINE)
          DO ISET=1,NSET
            IMEM=IMEMP(ISET)
            DO IREP=1,NREP(ISET)
              CALL REPLACE(REP(IREP,IMEM,ISET),LINE)
            ENDDO
          ENDDO
          ILINE2=ILINE2+1
          BUFFER(ILINE2)=LINE
        ENDDO
!
!       == INCREMENT INSTANCE =====================================
        I=1
        DO 
          IMEMP(I)=IMEMP(I)+1
          IF(IMEMP(I).GT.NMEM(I)) THEN
            IMEMP(I)=1
            I=I+1
            IF(I.GT.NSET) GOTO 2000
          ELSE
            EXIT
          END IF
        ENDDO
      ENDDO
 2000 CONTINUE
!
!     ===============================================================
!     == COMMENT PROTOTYPE AND COMMANDS                            ==
!     ===============================================================
      DO ILINE=IBEGIN,IINSTANCE
        BUFFER(ILINE)='!'//BUFFER(ILINE)
      ENDDO
      BUFFER(IEND)='!'//BUFFER(IEND)
      
!
!     ===============================================================
!     == DETERMINE CALLS OF TEMPLATE                               ==
!     ===============================================================
      DO ILINE=1,NLINE
        IF(BUFFER(ILINE)(1:1).NE.'#') CYCLE
        CALL UPPERCASE1(BUFFER(ILINE),LINE)
        IF(INDEX(LINE,'MODULE').EQ.0) CYCLE
        IF(INDEX(LINE,TRIM(NAME)).EQ.0) CYCLE
        IF(INDEX(LINE,'TEMPLATE').EQ.0) CYCLE
!
!       == CREATE SPACE         
        ISHIFT=NINST-1
        DO I=NLINE,ILINE-1,-1
          BUFFER(I+ISHIFT)=BUFFER(I)
        ENDDO
        NLINE=NLINE+ISHIFT
        IF(IBEGIN.LT.ILINE) IBEGIN=IBEGIN+ISHIFT
!
!       == INSERT
        DO I=1,NINST
          BUFFER(ILINE+I-1)=INSTANCE(I)
        ENDDO
      ENDDO
      RETURN
      END
!     ................................................................
      SUBROUTINE NEXT(LINE,I1,I2)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)      :: LINE
      INTEGER(4)  ,INTENT(INOUT)   :: I1
      INTEGER(4)  ,INTENT(OUT)     :: I2
      INTEGER(4)                   :: LENG,IA,IB,IC
!     ****************************************************************
      LENG=LEN(LINE)
      IA=INDEX(LINE(I1:),'<')
      IB=INDEX(LINE(I1:),'[')
      IF(IA.EQ.0) IA=LENG+1
      IF(IB.EQ.0) IB=LENG+1
      IC=MIN(IA,IB)
      I1=IC+I1-1
!     == NOT FOUND ======
      IF(IC.GE.LENG+1) THEN
        I1=0
        I2=0
        RETURN
      END IF
!     == FIND END OF WORD
      IF(LINE(I1:I1).EQ.'<') THEN
        I2=INDEX(LINE(I1:),'>')
      ELSE 
        I2=INDEX(LINE(I1:),']')
      END IF
!     == NOT FOUND ======
      IF(I2.EQ.0) THEN
        I1=0
        I2=0
        RETURN
      END IF
      I2=I2+I1-1
      RETURN
      END        
!
!     ...............................................................  
      SUBROUTINE COLLECT_REPLACEMENT(NREPX,NREP,REP,LINE)
!     ****************************************************************
!     **  READS THE REPLACEMENT RULES FROM A LINE                   **
!     ****************************************************************
      USE TEMPLATE_MODULE
!      USE STRINGS
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN)    :: NREPX
      INTEGER(4)    ,INTENT(INOUT) :: NREP
      TYPE(REP_TYPE),INTENT(INOUT) :: REP(NREPX)
      CHARACTER(*)  ,INTENT(IN)    :: LINE
      CHARACTER(LEN(LINE))         :: COPY
      INTEGER(4)                   :: I1,I2
!     ****************************************************************
      COPY=LINE
      DO WHILE(INDEX(COPY,'<<').NE.0)
        NREP=NREP+1
        IF(NREP.GT.NREPX) THEN
          CALL ERROR$MSG('TOO MANY REPLACEMENTS')
          CALL ERROR$CHVAL('LINE',TRIM(COPY))
          CALL ERROR$STOP('COLLECTREPLACEMENTS')
        END IF
!     
!       ==  REPLACEMENT%OUT
        I1=INDEX(COPY,'<<')
        I2=INDEX(COPY,'>>')
        IF(I2.LE.I1) THEN
          CALL ERROR$MSG('SYNTAX ERROR')
          CALL ERROR$CHVAL('COPY',TRIM(COPY))
          CALL ERROR$STOP('COLLECTREPLACEMENTS')
        END IF
        REP(NREP)%OUT=COPY(I1:I2+1)
        REP(NREP)%OUTLEN=I2-I1+2
        COPY=COPY(I2+2:)
!     
!       == REPLACEMENT%IN
        I1=INDEX(COPY,'[[')
        I2=INDEX(COPY,']]')
        IF(I2.LE.I1) THEN
          CALL ERROR$MSG('SYNTAX ERROR')
          CALL ERROR$CHVAL('COPY',TRIM(COPY))
          CALL ERROR$STOP('RESOLVETEMPLATE')
        END IF
        REP(NREP)%IN=COPY(I1+2:I2-1)
        REP(NREP)%INLEN=I2-I1-2
        COPY=COPY(I2+2:)
      ENDDO
      RETURN
      END
!
!     ...............................................................  
      SUBROUTINE REPLACE(REP,LINE)
!     ****************************************************************
!     **  APPLY A REPLACEMENT RULE TO A LINE                        **
!     ****************************************************************
      USE TEMPLATE_MODULE
!      USE STRINGS
      IMPLICIT NONE
      TYPE(REP_TYPE),INTENT(IN)   :: REP
      CHARACTER(*)  ,INTENT(INOUT):: LINE
      INTEGER(4)                  :: I1,I2
      CHARACTER(LEN(REP%OUT))     :: TEST1
      CHARACTER(LEN(LINE))        :: LINE1
!     ****************************************************************
      IF(REP%OUTLEN.EQ.0) THEN
        CALL ERROR$STOP('REPLACE')
      END IF
      IF(INDEX(REP%IN(1:REP%INLEN),REP%OUT(1:REP%OUTLEN)).NE.0) THEN
        CALL ERROR$STOP('REPLACE')
      END IF
      CALL UPPERCASE1(LINE,LINE1)
      LINE=LINE1
      CALL UPPERCASE1(REP%OUT(1:REP%OUTLEN),TEST1)
      I1=INDEX(LINE,TEST1(1:REP%OUTLEN))
!     I1=INDEX(+LINE,+REP%OUT(1:REP%OUTLEN))
      DO WHILE (I1.NE.0)
        I2=I1+REP%OUTLEN-1
        LINE=LINE(1:I1-1)//REP%IN(1:REP%INLEN)//LINE(I2+1:)
!       I1=INDEX(LINE,REP%OUT(1:REP%OUTLEN))
        I1=INDEX(LINE,TEST1(1:REP%OUTLEN))
      ENDDO
      RETURN
      END
