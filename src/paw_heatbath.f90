MODULE HEATBATH_MODULE
!***********************************************************************
!***********************************************************************
!**  DATA TO BE SET BEFORE STARTING                                   **
!**    TON                                                            **
!**    NPULSE                                                         **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!***********************************************************************
!***********************************************************************
LOGICAL(4)            :: TON=.FALSE.
INTEGER(4)            :: ITOTAL         ! TIMESTEP COUNTER
LOGICAL(4)            :: TINI=.FALSE.
REAL(8)   ,PARAMETER  :: TOL=1.D-10
REAL(8)   ,ALLOCATABLE:: DELTAV(:,:)    !(3,NAT) PULSE DIRECTION
INTEGER(4)            :: NPULSE=30      ! #PULSES
INTEGER(4)            :: NSTEP=10       ! #STEPS PER PULSE
INTEGER(4)            :: NAT=0          ! #ATOMS
REAL(8)               :: TFINAL=0.D0   ! #TARGET TEMPERATURE IN HARTREE!!
REAL(8)               :: RKB=0.D0       ! BOLTZMANN CONSTANT
END MODULE HEATBATH_MODULE
!
!     ..................................................................
      SUBROUTINE HEATBATH_INITIALIZE
!     ******************************************************************
!     ******************************************************************
      USE HEATBATH_MODULE
      IMPLICIT NONE
      INTEGER(4)   :: NFILO
!     ******************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
      ITOTAL=0
      ALLOCATE(DELTAV(3,NAT))
      IF(NPULSE.EQ.0)NPULSE=3*NAT-6
      CALL CONSTANTS('KB',RKB)
!     
!     ================================================================
!     ==  SET DEAULT FINAL TEMPERATURE                              ==
!     ================================================================
      IF(TFINAL.LE.0.D0)TFINAL=300.D0*RKB 
!     
!     ================================================================
!     ==  WRITE TO PROTOCOLL                                        ==
!     ================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,*)
      WRITE(NFILO,FMT='("HEATBATH    ")')
      WRITE(NFILO,FMT='("==============================")')
      WRITE(NFILO,FMT='("NUMBER OF PULSES:     ",I4)')NPULSE
      WRITE(NFILO,FMT='("NUMBER OF STEPS/PULSE:",I4)')NSTEP
      WRITE(NFILO,FMT='("TOTAL PUMPED DELTA T: ",F10.5)')TFINAL/RKB
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE HEATBATH_DESTRUCT
!     ******************************************************************
!     ******************************************************************
      USE HEATBATH_MODULE
      IMPLICIT NONE
      INTEGER(4)   :: NFILO
!     ******************************************************************
      TON=.FALSE.
      DEALLOCATE(DELTAV)
!     
!     ================================================================
!     ==  WRITE TO PROTOCOLL                                        ==
!     ================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,*)'HEATBATH: FINISHING HEATBATH '
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE HEATBATH$SETL4(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE HEATBATH_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'ON') THEN
        TON=VAL
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$STOP('HEATBATH$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE HEATBATH$SETI4(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE HEATBATH_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NPULSE') THEN
        NPULSE=VAL
      ELSE IF(ID.EQ.'NSTEP') THEN
        NSTEP=VAL
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$STOP('HEATBATH$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE HEATBATH$SETR8(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE HEATBATH_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'TFINAL') THEN
        TFINAL=VAL
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$STOP('HEATBATH$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE ATOMS_HEATBATHAPPLY(NAT_,RM,R0,RP,FORCE,RMASS,DELT)
!     ******************************************************************
!     **                                                              **
!     **  ADDS DELTA T TO ATOMIC VELOCITIES (RM)                      **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    RM,R0     POSITIONS OF THE ATOMS, TFINAL DTERMINES HOW    **
!     **      MUCH ENERGY IS FINALLY GIVEN TO THE ATOMS               **
!     **      NPULSE IS THE NUMBER OF SINUSOIDAL PULSES               **
!     **      NSTEP DETERMINES HOW MANY TIMESTEP ONE PULSE LASTS      **
!     **                                                              **
!     **  OUTPUT:                                                     **
!     **    RM        CHANGED R(-DELTA T)                             **
!     **                                                              **
!     ******************************************************************
      USE HEATBATH_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_ 
      REAL(8)   ,INTENT(IN) :: R0(3,NAT_)
      REAL(8)   ,INTENT(IN) :: RM(3,NAT_)
      REAL(8)   ,INTENT(IN) :: RMASS(NAT_)
      REAL(8)   ,INTENT(OUT):: RP(3,NAT_)
      REAL(8)   ,INTENT(IN) :: FORCE(3,NAT_)
      REAL(8)   ,INTENT(IN) :: DELT
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: IPULSE       ! PULSE COUNT
      INTEGER(4)            :: ISTEP        ! STEPS COUNT FOR THIS PULSE
      LOGICAL(4)            :: TDEBUG=.FALSE.
      REAL(8)               :: PI,TWOPI
      INTEGER(4)            :: NG           !=3*NAT
      INTEGER(4)            :: I,IAT 
      REAL(8)               :: PRO(3,NAT_)   
      REAL(8)               :: FORCE1(3,NAT_)   
      REAL(8)               :: SVAR
      REAL(8)               :: RKAPPA
      REAL(8)               :: A,B,C,X
      REAL(8)               :: TEMPER,TARGETTEMP
!     ******************************************************************
      NAT=NAT_
      PI=4.D0*DATAN(1.D0)
      TWOPI=2.D0*PI
      NG=3*NAT
!
      IF(.NOT.TON) RETURN
      CALL HEATBATH_INITIALIZE
      ITOTAL=ITOTAL+1               ! TOTAL NUMBER OF TIMESTEPS IN HEATBATH
      IPULSE=(ITOTAL-1)/NSTEP+1     ! CURRENT PULSE
      ISTEP=ITOTAL-(IPULSE-1)*NSTEP ! # TIMESTEPS SPENT IN CURRENT PULSE
      IF(IPULSE.GT.NPULSE) THEN
        CALL HEATBATH_DESTRUCT
        RETURN
      ENDIF
!
!     ==================================================================
!     ==  CREATE NEW MODE FOR THIS PULSE                              ==
!     ==================================================================
      IF(ISTEP.EQ.1) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'APPLYING NEW HEATBATH VECTOR'
! 
        CALL ATOMS_HEATBATHVECTOR(IPULSE,NAT,RM,RP,RMASS,DELTAV)
      ENDIF
!
!     ==================================================================
!     == PRO= M|V>/SQRT(<V|M**2|V>)                                   ==
!     ==================================================================
! --- NORMALIZE FOR CONSTRUCTION OF PROJECTOR
! --- WHICH IS USED TO REMOVE ALL THE FORCES WHICH POINT
! --- IN THE DIRECTION OF |M|DELTAV>
      SVAR=0.D0
      DO IAT=1,NAT
        DO I=1,3
          PRO(I,IAT)=RMASS(IAT)*DELTAV(I,IAT)
          SVAR=SVAR+PRO(I,IAT)
        ENDDO
      ENDDO
      SVAR=1.D0/DSQRT(SVAR)
      DO IAT=1,NAT
        DO I=1,3
          PRO(I,IAT)=PRO(I,IAT)*SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     == PROJECT OUT FORCE COMPONENT PARALLEL TO PRO                  ==
!     ==================================================================
      SVAR=0.D0
      DO IAT=1,NAT
        DO I=1,3
          SVAR=SVAR+PRO(I,IAT)*FORCE(I,IAT)
        ENDDO
      ENDDO
      DO IAT=1,NAT
        DO I=1,3
          FORCE1(I,IAT)=FORCE(I,IAT)-PRO(I,IAT)*SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     == PROPAGATE ATOMIC POSITIONS                                   ==
!     ==================================================================
      DO IAT=1,NAT
        DO I=1,3
          RP(I,IAT)=2.D0*R0(I,IAT)-RM(I,IAT)+FORCE1(I,IAT)*DELT**2/RMASS(IAT)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == CALCULATE TEMPERATURE                                        ==
!     ==================================================================
      TEMPER=0.D0
      DO IAT=1,NAT
        DO I=1,3
          SVAR     =(RP(I,IAT)-R0(I,IAT))/DELT
          TEMPER   =TEMPER+RMASS(IAT)*SVAR**2**2
        ENDDO
      ENDDO
      TEMPER   =TEMPER/DBLE(NG)
!
!     ==================================================================
!     == CALCULATE REQUIRED TEMPERATURE                               ==
!     ==================================================================
      RKAPPA=TFINAL/DBLE(NPULSE*NSTEP) &
     &       *(1.D0+COS(TWOPI*DBLE(ITOTAL)/DBLE(NSTEP)))
      TARGETTEMP=TEMPER+RKAPPA
!
      IF(TDEBUG) THEN
        PRINT*,'TEMPERATURE AT BEGINNING OF LOOP FROM (+)',TEMPER/RKB
        PRINT*,'TARGET TEMPERATURE                       ',TARGETTEMP/RKB
      END IF
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
!     == R(+)=RBAR+X*DELTAV
!     == FINAL TEMPERATURE IS A+B*X+C*X**2
      A=0.D0
      B=0.D0
      C=0.D0
      DO IAT=1,NAT
        DO I=1,3
          SVAR=(RP(I,IAT)-R0(I,IAT))/DELT
          A=A+RMASS(IAT)*SVAR**2
          B=B+2.D0*RMASS(IAT)*DELTAV(I,IAT)*SVAR
          C=C+RMASS(IAT)*DELTAV(I,IAT)**2
        ENDDO
      ENDDO
      A=A/DBLE(NG)
      B=B/DBLE(NG)
      C=C/DBLE(NG)
      A=A-TARGETTEMP
      X=-0.5D0*B/C*(1.D0-SQRT(1.D0-4.D0*A*C/B**2))
!
!     ==================================================================
!     == UPDATE NEXT POSITIONS                                        ==
!     ==================================================================
      SVAR=0.D0
      DO IAT=1,NAT
        DO I=1,3
          RP(I,IAT)=RP(I,IAT)+DELTAV(I,IAT)*X
          SVAR=SVAR+RMASS(IAT)*((RP(I,IAT)-R0(I,IAT))/DELT)**2
        ENDDO
      ENDDO
      SVAR=SVAR/DBLE(NG)
      IF (TDEBUG)PRINT*,'TEMPERATURE AT END OF LOOP (1)',SVAR/RKB
!
      IF (DABS(SVAR-TARGETTEMP).GT.1.D-5) THEN
        CALL ERROR$MSG('HEATBATH: DESIRED DELTA T NOT EFFECTED')
        CALL ERROR$STOP('IN HEATBATH_APPLY')
      ENDIF
!
      RETURN
      END
!
!      
!     .................................................................
      SUBROUTINE ATOMS_HEATBATHVECTOR(IPULSE,NAT,RM,RP,RMASS,DELTAV)
!     *****************************************************************
!     ** RETURNS A VELOCITY VECTOR DELTAV                            **
!     ** DELTAV IS CHOSEN RANDOMLY IN THE FIRST PULSE, ALL SUBSEQUENT**
!     ** DELTAV'S ARE GENERATED IN DIRECTIONS NORMAL TO THE CURRENT  **
!     ** VELOCITIES                                                  **
!     *****************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IPULSE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RM(3,NAT)
      REAL(8)   ,INTENT(IN) :: RP(3,NAT)
      REAL(8)   ,INTENT(IN) :: RMASS(NAT)
      REAL(8)   ,INTENT(OUT):: DELTAV(3,NAT)
      INTEGER(4)            :: NFILI !FORTRAN FILE NUMBER FOR INPUT FILES
      INTEGER(4)            :: K,I,J
      REAL(8)               :: FACTOR,SVAR     
      REAL(8)               :: VEC(3,NAT)
      REAL(8)               :: RCHECK(3,NAT)
!     *****************************************************************
      IF (IPULSE .EQ. 1) THEN
        CALL RANDOM_NUMBER(DELTAV(:,:))
        DELTAV(:,:)=DELTAV(:,:)-0.5D0
      ELSE
!       ==================================================================
!       ==    CALCULATE VELOCITIES WEIGHTED WITH SQRT(M)
!       ==    BE AWARE THAT DELTAV CONTAINS NOW TEMPORARILY A VELOCITY
!       ==================================================================
        SVAR=0.D0
        DO I=1,NAT
          FACTOR=DSQRT(RMASS(I))
          DO J=1,3
            VEC(J,I)=(RP(J,I)-RM(J,I))*FACTOR
            SVAR=SVAR+VEC(J,I)**2
          ENDDO
        ENDDO
        SVAR=DSQRT(SVAR)
        DO I=1,NAT
          DO J=1,3
            VEC(J,I)=VEC(J,I)/SVAR
            RCHECK(J,I)=VEC(J,I)
          ENDDO
        ENDDO
!     
!       ==================================================================
!       == INVENT SOME RANDOM VECTOR AND NORMALIZE IT
!       ==================================================================
        CALL RANDOM_NUMBER(DELTAV(:,:))
        DELTAV(:,:)=DELTAV(:,:)-0.5D0
        SVAR=0.D0
        DO I=1,NAT
          DO J=1,3
            SVAR=SVAR+DELTAV(J,I)**2
          ENDDO
        ENDDO
        SVAR=DSQRT(SVAR)
        DO I=1,NAT
          DO J=1,3
            DELTAV(J,I)=DELTAV(J,I)/SVAR
          ENDDO
        ENDDO
!       ==================================================================
!       ==  PROJECT UPON PREVIOUS VECTORS
!       ==================================================================
        SVAR=0.D0
        DO I=1,NAT
          DO J=1,3
            SVAR=SVAR+VEC(J,I)*DELTAV(J,I)
          ENDDO
        ENDDO
        DELTAV(:,:)=DELTAV(:,:)-VEC(:,:)*SVAR
      
        SVAR=0.D0
        DO I=1,NAT
          DO J=1,3
            SVAR=SVAR+DELTAV(J,I)**2
          ENDDO
        ENDDO
        IF (DABS(SVAR).GT.1.D-10) THEN
          CALL ERROR$MSG('FAILED TO FIND ORTHOGONAL VELOCITIES!')
          CALL ERROR$STOP('IN HEATBATH')
        ENDIF
      ENDIF
!
!     ==================================================================
!     == DIVIDE BY SQRT(MASS) BEFORE HANDING IT BACK
!     ==================================================================
      DO I=1,NAT
        FACTOR=1.D0/DSQRT(RMASS(I))
        DO J=1,3
          DELTAV(J,I)=DELTAV(J,I)*FACTOR
        ENDDO
      ENDDO
!
      RETURN
      END











