MODULE CORE_MODULE
!**                                                                   **
!**  THE COMPLEX POINTER LIST THIS IS NEEDED FOR THE PARALLELIZATION  **
!**  BECAUSE EVERY NODE WILL CALCULATE THE CORE LEVELS ONLY FOR ITS   **
!**  OWN SET OF ATOMS AND ONLY THE FIRST NODE CAN WRITE               **
!**                                                                   **
!**  COMMUNIATION IN CORE$REPORT NOT DONE YET.                        **
!**                                                                   **
TYPE CORESHIFT_TYPE
INTEGER(4)           :: IAT
INTEGER(4)           :: N
CHARACTER(8),POINTER :: TYPE(:)
REAL(8)     ,POINTER :: E(:)
REAL(8)     ,POINTER :: EATOM(:)
TYPE(CORESHIFT_TYPE),POINTER :: NEXT
END TYPE CORESHIFT_TYPE
LOGICAL(4)                  :: TCORESHIFTS=.TRUE.
LOGICAL(4)                  :: DEFAULT=.FALSE.
INTEGER(4)                  :: NATOMS=0
CHARACTER(32),ALLOCATABLE   :: ATOMS(:)
TYPE(CORESHIFT_TYPE),TARGET :: FIRST
TYPE(CORESHIFT_TYPE),POINTER:: THIS
LOGICAL(4),SAVE             :: TINI=.FALSE.
END MODULE CORE_MODULE
!
!     ..................................................................
      SUBROUTINE CORE$SETL4(ID,VAL)
      USE CORE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'ON') THEN
        TCORESHIFTS=VAL
      ELSE IF(ID.EQ.'DEFAULT') THEN
        DEFAULT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CORE$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CORE$SETCHA(ID,LEN,VAL)
      USE CORE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
      INTEGER(4)              :: I
!     ******************************************************************
      IF(ID.EQ.'ATOMS') THEN
        IF(ALLOCATED(ATOMS))DEALLOCATE(ATOMS)
        NATOMS=LEN
        ALLOCATE(ATOMS(NATOMS))
        DO I=1,LEN
          ATOMS(I)=VAL(I)
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CORE$SETCHA')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE$REPORT(NFIL)
      USE CORE_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: NFIL
      REAL(8)                  :: EV
      CHARACTER(32)            :: NAME
      INTEGER(4)               :: I
      INTEGER(4)               :: N
      REAL(8)                  :: SVAR
      INTEGER(4)               :: THISTASK,NTASKS
      CHARACTER(128)           :: STRING
      INTEGER(4)               :: NAT
      INTEGER(4)               :: IAT
      INTEGER(4)  ,ALLOCATABLE :: TASKARR(:)
      INTEGER(4)  ,ALLOCATABLE :: TASKARR1(:)
      REAL(8)     ,ALLOCATABLE :: E(:)
      REAL(8)     ,ALLOCATABLE :: EATOM(:)
      CHARACTER(8),ALLOCATABLE :: XTYPE(:)
      INTEGER(4)               :: SNDTASK,RCVTASK
!     ******************************************************************
      IF(.NOT.TINI) RETURN
                             CALL TRACE$PUSH('CORE$REPORT')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL CONSTANTS('EV',EV)

!     == CREATE A MAPPING ARRAY FROM ATOMS TO TASKS FOR EACH NODE.          ==
!     == ATOMS NOT PRESENT ON THIS TASK OBTAIN 0 FOR THE ATOM NUMBER        ==
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(TASKARR(NAT))
      TASKARR(:)=0
!
      IF(FIRST%IAT.NE.0) THEN 
        THIS=>FIRST
        DO 
          IAT=THIS%IAT
          TASKARR(IAT)=THISTASK
          IF(.NOT.ASSOCIATED(THIS%NEXT)) EXIT
          THIS=>THIS%NEXT
        ENDDO
      ENDIF
!
!     == COMMUNICATE A UNIQUE MAPPING KNOWN TO ALL TASKS =====================
! THIS COMMUNICATION IS UNNECCESARILY COMPLICATED.
      ALLOCATE(TASKARR1(NAT))
      TASKARR1(:)=0
      DO I=1,NTASKS
        SNDTASK=THISTASK+1
        SNDTASK=MODULO(SNDTASK-1,NTASKS)+1
        RCVTASK=THISTASK-1
        RCVTASK=MODULO(RCVTASK-1,NTASKS)+1
        !== MPE$COMBINE WITH THE MAX FUNCTION MAY HAVE BEEN BETTER
        TASKARR1=TASKARR
        CALL MPE$SENDRECEIVE('MONOMER',THISTASK,SNDTASK,TASKARR1)
        CALL MPE$SENDRECEIVE('MONOMER',RCVTASK,THISTASK,TASKARR1)
        DO IAT=1,NAT
          TASKARR(IAT)=MAX(TASKARR(IAT),TASKARR1(IAT))
        ENDDO
      ENDDO       
      DEALLOCATE(TASKARR1)
!     == TASKARRAY IS NOW COMMON TO ALL TASKS        

      DO IAT=1,NAT
        IF(TASKARR(IAT).EQ.0) CYCLE ! NO INFORMATION FOR THIS ATOM AVAILABLE

        IF(THISTASK.EQ.TASKARR(IAT)) THEN
          THIS=>FIRST
          DO 
            IF(THIS%IAT.EQ.IAT) EXIT
            IF(.NOT.ASSOCIATED(THIS%NEXT)) EXIT
            THIS=>THIS%NEXT
          ENDDO
          N=THIS%N
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',TASKARR(IAT),1,N)
        ALLOCATE(E(N))
        ALLOCATE(EATOM(N))
        ALLOCATE(XTYPE(N))
        IF(THISTASK.EQ.TASKARR(IAT)) THEN
          E=THIS%E
          EATOM=THIS%EATOM
          XTYPE=THIS%TYPE
        END IF
        CALL MPE$SENDRECEIVE('MONOMER',TASKARR(IAT),1,E)
        CALL MPE$SENDRECEIVE('MONOMER',TASKARR(IAT),1,EATOM)
        CALL MPE$SENDRECEIVE('MONOMER',TASKARR(IAT),1,XTYPE)
!       == TASK 1 RECEIVED DATA ==============================================
!
!       == PRINT INFORMATION ON THE FIRST TASK ===============================
        IF(THISTASK.EQ.1) THEN
          CALL ATOMLIST$GETCH('NAME',IAT,NAME)
          CALL REPORT$TITLE(NFIL,'EIGENVALUES OF CORE STATES FROM ATOM '& 
     &                     //TRIM(NAME))
          WRITE(NFIL,*)'     CURRENT SYSTEM          ISOLATED ATOM' 
          STRING='(T3,"ENERGY[H]",T14,"ENERGY[EV]"'
          STRING=TRIM(ADJUSTL(STRING))//',T28,"ENERGY[H]",T39,"ENERGY[EV]"'
          STRING=TRIM(ADJUSTL(STRING))//',T53,"SHIFT[H]",T64,"SHIFT[EV]"'
          STRING=TRIM(ADJUSTL(STRING))//')'
          WRITE(NFIL,STRING)
          DO I=1,N
            WRITE(NFIL,FMT='(A5,6F12.5)')TRIM(XTYPE(I)) &
     &                           ,E(I),E(I)/EV &
     &                           ,EATOM(I),EATOM(I)/EV &
     &                           ,E(I)-EATOM(I),(E(I)-EATOM(I))/EV
          ENDDO
        END IF
        DEALLOCATE(E)
        DEALLOCATE(EATOM)
        DEALLOCATE(XTYPE)
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END
!
! SANTOS040617 BEGIN
!     ..................................................................
      SUBROUTINE CORE_CORESHIFTS(IAT,ISP,GID,NR,LMRXX,AEPOT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE EIGENVALUES OF CORE HAMILTONIAN              **
!     **                                                              **
!     ******************************************************************
!      USE ATOMS_MODULE
      USE CORE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: ISP
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRXX
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRXX)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: ATPOT(NR)
      REAL(8)   ,ALLOCATABLE:: AEPOT1(:,:) !(NR,LMRXX)
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC
      INTEGER(4),ALLOCATABLE:: LB(:)
!     REAL(8)   ,ALLOCATABLE:: FB(:)
      REAL(8)   ,ALLOCATABLE:: EB(:)
      REAL(8)   ,ALLOCATABLE:: AEPSI(:,:)
      REAL(8)               :: R(NR)
      CHARACTER(32)         :: NAME
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: LMRX

      INTEGER(4)            :: I
      REAL(8)               :: AUX2

      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LM1,LM2,LM3
      INTEGER(4)            :: L1,L2,IM1,IM2
      REAL(8)               :: AEDMU(NR)
      REAL(8)               :: DWORK1(NR)
      REAL(8)               :: CG
      REAL(8)               :: SVAR

      REAL(8)   ,ALLOCATABLE:: HAMIL(:,:)
      REAL(8)   ,ALLOCATABLE:: EIGENVAL(:)
      REAL(8)   ,ALLOCATABLE:: EIGENVEC(:,:)
      INTEGER(4)            :: NFILO
      REAL(8)               :: EV             ! ELECTRON VOLT
      LOGICAL(4)            :: TCHK
!     ******************************************************************
      IF(.NOT.TCORESHIFTS) RETURN
      CALL RADIAL$R(GID,NR,R)
      CALL ATOMLIST$GETCH('NAME',IAT,NAME)
      TCHK=DEFAULT
      DO I=1,NATOMS
        IF(NAME.EQ.ATOMS(I)) THEN
          TCHK=.NOT.DEFAULT
        END IF
      ENDDO
      IF(.NOT.TCHK) RETURN
!
!     ===================================================================
!     ==  SELECT THE PROPER ENTRY IN THE TABLE                         ==
!     ===================================================================
      THIS=>FIRST
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        THIS%IAT=IAT
        THIS%N=0
        NULLIFY(THIS%E)
        NULLIFY(THIS%EATOM)
        NULLIFY(THIS%TYPE)
        NULLIFY(THIS%NEXT)
      ELSE
        DO
          IF(THIS%IAT.EQ.IAT) EXIT
          IF(ASSOCIATED(THIS%NEXT)) THEN
            THIS=>THIS%NEXT
          ELSE
            ALLOCATE(THIS%NEXT)
            THIS=>THIS%NEXT
            NULLIFY(THIS%NEXT)
            THIS%IAT=IAT
            THIS%N=0
            NULLIFY(THIS%E)
            NULLIFY(THIS%EATOM)
            NULLIFY(THIS%TYPE)
          END IF
        ENDDO
      END IF
!
!     ==================================================================
!     ==  COLLECT CORE PARTIAL WAVES FROM SETUP OBJECT                ==
!     ==================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('NB',NB)
      CALL SETUP$GETI4('NC',NC)
      IF (NC.EQ.0) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        CALL REPORT$TITLE(NFILO,'NO CORE STATES FOR ATOM '//TRIM(NAME))
        RETURN
      END IF
      CALL SETUP$GETR8A('AEPOT',NR,ATPOT)
      ALLOCATE(LB(NB))
      CALL SETUP$GETI4A('LB',NB,LB)
      ALLOCATE(EB(NB))
      CALL SETUP$GETR8A('EB',NB,EB)
      ALLOCATE(AEPSI(NR,NB))
      CALL SETUP$GETR8A('AEPSI',NR*NB,AEPSI)
      CALL SETUP$UNSELECT()
!
!     ==================================================================
!     ==  CONSTANTS                                                   ==
!     ==================================================================
      LMNX=0
      LMRX=0
      DO I=1,NC
        LMNX=LMNX+2*LB(I)+1
        LMRX=MAX(LMRX,(2*LB(I)+1)**2)
      ENDDO
      LMRX=MIN(LMRX,LMRXX)
      IF(THIS%N.NE.LMNX) THEN
        THIS%N=LMNX
        ALLOCATE(THIS%TYPE(LMNX))
        ALLOCATE(THIS%E(LMNX))
        ALLOCATE(THIS%EATOM(LMNX))
      END IF

!     ==================================================================
!     == SUBTRACTS ATOMIC AE POTENTIAL FROM AE TOTAL POTENTIAL        ==
!     ==================================================================
      ALLOCATE(AEPOT1(NR,LMRXX))
      AEPOT1(:,:)=AEPOT(:,:)
      AEPOT1(:,1)=AEPOT(:,1)-ATPOT(:)

!     ==================================================================
!     ==   CALCULATE HAMILTONIAN                                      ==
!     ==================================================================
      ALLOCATE(HAMIL(LMNX,LMNX))      
      HAMIL(:,:)=0.D0
!
      LMN1=0
      DO LN1=1,NC
        L1=LB(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          LM1=L1**2+IM1
          DO LN2=1,NC
            L2=LB(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
            
              IF(LMN1.EQ.LMN2) THEN
                HAMIL(LMN1,LMN2)=EB(LN1)
              END IF
!
              AEDMU(:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  AEDMU(:)=AEDMU(:)+CG*AEPOT1(:,LM3)
                END IF
              ENDDO
!
              DWORK1(:)=(AEDMU(:)*AEPSI(:,LN1)*AEPSI(:,LN2))*R(:)**2
              CALL RADIAL$INTEGRAL(GID,NR,DWORK1,SVAR)
              HAMIL(LMN1,LMN2)=HAMIL(LMN1,LMN2)+SVAR
!
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(AEPOT1)
!
!     ==================================================================
!     ==  DIAGONALIZATION OF THE HAMILTONIAN                          ==
!     ==================================================================
      ALLOCATE(EIGENVAL(LMNX))
      ALLOCATE(EIGENVEC(LMNX,LMNX))
      CALL LIB$DIAGR8(LMNX,HAMIL,EIGENVAL,EIGENVEC)
      DEALLOCATE(EIGENVEC)
!
!     ==================================================================
!     ==  WRITE INTO TABLE                                            ==
!     ==================================================================
      LMN1=0
      DO LN1=1,NC
        L1=LB(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          THIS%E(LMN1)=EIGENVAL(LMN1)
          THIS%EATOM(LMN1)=EB(LN1)
          I=L1
          DO LN2=1,LN1
            IF(LB(LN2).NE.L1) CYCLE
            I=I+1
          ENDDO
          WRITE(THIS%TYPE(LMN1),FMT='(I2)')I
          IF(L1.EQ.0) THEN
            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'S'
          ELSE IF(L1.EQ.1) THEN
            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'P'
          ELSE IF(L1.EQ.2) THEN
            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'D'
          ELSE IF(L1.EQ.3) THEN
            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'F'
          ELSE
            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'?'
          END IF
        ENDDO
      ENDDO
      

!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(LB)
!     DEALLOCATE(FB)
      DEALLOCATE(EB)
      DEALLOCATE(AEPSI)
      DEALLOCATE(HAMIL)
      DEALLOCATE(EIGENVAL)

      RETURN
      END
!
! SANTOS040617 END
