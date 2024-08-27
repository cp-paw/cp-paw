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
     &        //'," MAXIMUM ALLOWED ",I10)')IACTUAL,ITARGET
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ERROR$C8VAL(MESSAGE_,C8VAL)
!     ******************************************************************
!     ** STORE MESSAGE WITH DOUBLE PRECISION VALUE                    **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      COMPLEX(8)  ,INTENT(IN) :: C8VAL
!     ******************************************************************
      IMESSAGE=IMESSAGE+1
      IF(IMESSAGE.GT.IMESSAGEX-1) THEN
        IMESSAGE=IMESSAGE-1
        MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
      ELSE 
        WRITE(MESSAGE(IMESSAGE),FMT='("VARIABLE ",A' &
     &       //'," HAS THE VALUE (",E20.10,",",E20.10,")")')MESSAGE_,C8VAL
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
     &        //'," HAS THE VALUE ",E20.10)')MESSAGE_,R8VAL
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
     &        //'," HAS THE VALUE ",I10)')MESSAGE_,I4VAL
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
     &        //'," HAS THE VALUE ",L7)')MESSAGE_,L4VAL
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
        IF(IMESSAGE.GT.IMESSAGEX-1) THEN
          IMESSAGE=IMESSAGE-1
          MESSAGE(IMESSAGE)='TOO MANY ERROR MESSAGE IN ERROR'
          RETURN
        END IF
        WRITE(MESSAGE(IMESSAGE),FMT='("VARIABLE ",A' &
     &          //'," HAS THE VALUE ",A)')TRIM(MESSAGE_),TRIM(CHVAL)
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
      RECURSIVE SUBROUTINE ERROR$STOP(MESSAGE_)
!     ******************************************************************
!     **  WRITE ERROR MESSAGE, FLUSH FILES AND STOP                   **
!     ******************************************************************
      USE ERROR_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MESSAGE_
      INTEGER(4)              :: NTASKNUM
      INTEGER(4)              :: NTASKID
      INTEGER(4)              :: NFILO
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
        CALL MPE$QUERY('~',NTASKNUM,NTASKID)
        IF(NTASKNUM.GT.1) THEN
          WRITE(NFILERR,FMT='("PROCESSOR REQUESTING STOP: ",I3' &
     &                    //'," OF ",I3)')NTASKID,NTASKNUM
        END IF
      END IF
!
!     ================================================================
!     == WRITE ERROR MESSAGE ON PROTOCOL                            ==
!     ================================================================
      IF(NCALL.LE.1) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO) 
        DO I=1,IMESSAGE
          WRITE(NFILO,FMT='(A)')TRIM(MESSAGE(I))
        ENDDO
        WRITE(NFILO,FMT='("STOP IN ",A)')TRIM(MESSAGE_)
        IF(NTASKNUM.GT.1) THEN
          WRITE(NFILO,FMT='("PROCESSOR REQUESTING STOP: ",I3' &
     &         //'," OF ",I3)')NTASKID,NTASKNUM
        END IF
      END IF
!
!     ================================================================
!     == FLUSH ALL FILES                                            ==
!     ================================================================
!      FLUSH(6)
!      DO I=1000,1100
!        FLUSH(I)
!      ENDDO
!
!     ================================================================
!     == STOP                                                       ==
!     ================================================================
      IF(NCALL.LE.3) THEN
        CALL FILEHANDLER$CLOSEALL
      END IF
      CALL MPE$STOPALL(NCODE)
!     == JOB STOPS IN MPE$STOPALL WITH ERROR CODE NCODE ========================
      ERROR STOP 'IN ERROR$STOP' ! THIS STATEMENT IS SUPERFLUOUS
      END
!
!     ..................................................................
      SUBROUTINE ERROR$NORMALSTOP
!     ******************************************************************
!     **  WRITE ERROR MESSAGE, FLUSH FILES AND STOP                   **
!     ******************************************************************
      CALL FILEHANDLER$CLOSEALL
!
!     ==========================================================================
!     == MPI$STOPALL WILL STOP THE EXECUTION ALSO IN THE SCALAR VERSION
!     ==========================================================================
      WRITE(*,FMT='(A)')'NORMAL STOP: CALLING MPE$STOPALL TO CLOSE DOWN'
      CALL MPE$STOPALL(0) ! BRING DOWN ALL MPI TASKS RELATED TO MPI_WORLD_COMM
!
!     ==========================================================================
!     == THE FOLLOWING STOP IS SUPERFLUOUS
!     ==========================================================================
      STOP 'NORMAL STOP'  ! THIS STATEMENT IS SUPERFLUOUS
      END
