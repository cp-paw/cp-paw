MODULE CORE_MODULE
!**                                                                   **
!**  THE COMPLEX POINTER LIST THIS IS NEEDED FOR THE PARALLELIZATION  **
!**  BECAUSE EVERY NODE WILL CALCULATE THE CORE LEVELS ONLY FOR ITS   **
!**  OWN SET OF ATOMS AND ONLY THE FIRST NODE CAN WRITE               **
!**                                                                   **
!**  COMMUNIATION IN CROE$REPORT NOT DONE YET.                        **
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
!     ..................................................................
      SUBROUTINE CORE$REPORT(NFIL)
      USE CORE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      REAL(8)                 :: EV
      CHARACTER(32)           :: NAME
      INTEGER(4)              :: I,IAT
      REAL(8)                 :: SVAR
      INTEGER(4)              :: THISTASK,NTASKS
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(thistask.NE.1) RETURN
      CALL CONSTANTS('EV',EV)
      THIS=>FIRST
      DO 
        IF(THIS%IAT.EQ.0) EXIT
        CALL ATOMLIST$GETCH('NAME',THIS%IAT,NAME)
        CALL REPORT$TITLE(NFIL,'EIGENVALUES OF CORE STATES FROM ATOM '& 
     &                       //TRIM(NAME))
        WRITE(NFIL,*)'     CURRENT SYSTEM          ISOLATED ATOM' 
        WRITE(NFIL,*)'  ENERGY[H]  ENERGY[EV]   ENERGY[H]  ENERGY[EV]&
     &    SHIFT[H]   SHIFT[EV]'
        DO I=1,THIS%N
          SVAR = THIS%E(I) - THIS%EATOM(I)
          WRITE(NFIL,FMT='(A5,6F12.5)')TRIM(THIS%TYPE(I)),THIS%E(I),THIS%E(I)/EV &
    &          ,THIS%EATOM(I),THIS%EATOM(I)/EV,SVAR,SVAR/EV
        ENDDO
        IF(.NOT.ASSOCIATED(THIS%NEXT)) EXIT
        THIS=>THIS%NEXT
      ENDDO  
      RETURN
      END
!
! SANTOS040617 BEGIN
!     ..................................................................
      SUBROUTINE CORE_CORESHIFTS(IAT,ISP,gid,NR,LMRXX,AEPOT)
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
      INTEGER(4),INTENT(IN) :: gid
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRXX
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRXX)

      REAL(8)               :: ATPOT(NR)
      REAL(8)               :: AEPOT1(NR,LMRXX)
      INTEGER(4)            :: NC
      INTEGER(4),ALLOCATABLE:: LB(:)
!     REAL(8)   ,ALLOCATABLE:: FB(:)
      REAL(8)   ,ALLOCATABLE:: EB(:)
      REAL(8)   ,ALLOCATABLE:: AEPSI(:,:)
      REAL(8)               :: r(nr)
      CHARACTER(32)         :: NAME
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: LMRX

      REAL(8)               :: PI
      INTEGER(4)            :: I
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: AUX2
      CHARACTER(LEN=82)     :: STRING

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
      call radial$r(gid,nr,r)
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
      CALL SETUP$GETI4('NC',NC)
      IF (NC.EQ.0) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        CALL REPORT$TITLE(NFILO,'NO CORE STATES FOR ATOM '//TRIM(NAME))
        RETURN
      END IF
      CALL SETUP$GETR8A('ATOMICAEPOT',NR,ATPOT)
      ALLOCATE(LB(NC))
      CALL SETUP$GETI4A('LB',NC,LB)
!     ALLOCATE(FB(NC))
!     CALL SETUP$FB(ISP,NC,FB)
      ALLOCATE(EB(NC))
      CALL SETUP$GETR8A('EB',NC,EB)
      ALLOCATE(AEPSI(NR,NC))
      CALL SETUP$GETR8A('AECOREPSI',NR*NC,AEPSI)
!
!     ==================================================================
!     ==  CONSTANTS                                                   ==
!     ==================================================================
      PI=4.D0*DATAN(1.D0)

      LMNX=0
      LMRX=0
      DO I=1,NC
        LMNX=LMNX+2*LB(I)+1
        LMRX=MAX(LMRX,(2*LB(I)+1)**2)
      ENDDO
      IF(THIS%N.NE.LMNX) THEN
        THIS%N=LMNX
        ALLOCATE(THIS%TYPE(LMNX))
        ALLOCATE(THIS%E(LMNX))
        ALLOCATE(THIS%EATOM(LMNX))
      END IF

!     ==================================================================
!     == SUBTRACTS ATOMIC AE POTENTIAL FROM AE TOTAL POTENTIAL        ==
!     ==================================================================
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
              CALL RADIAL$INTEGRAL(gid,NR,DWORK1,SVAR)
              HAMIL(LMN1,LMN2)=HAMIL(LMN1,LMN2)+SVAR
!
            ENDDO
          ENDDO
        ENDDO
      ENDDO
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
            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'??'
          END IF
        ENDDO
      ENDDO
      
!
!!$!     ==================================================================
!!$!     ==  PRINTOUT                                                    ==
!!$!     ==================================================================
!!$      CALL CONSTANTS('EV',EV)
!!$!     ==  PRINTOUT
!!$      CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$      CALL REPORT$TITLE(NFILO,'EIGENVALUES OF CORE STATES FROM ATOM '& 
!!$        &//TRIM(NAME))
!!$      WRITE(NFILO,*)'     CURRENT SYSTEM          ISOLATED ATOM' 
!!$      WRITE(NFILO,*)'  ENERGY[H]  ENERGY[EV]   ENERGY[H]  ENERGY[EV]&
!!$        &    SHIFT[H]   SHIFT[EV]'
!!$
!!$      LMN1=0
!!$      DO LN1=1,NC
!!$        L1=LB(LN1)
!!$        DO IM1=1,2*L1+1
!!$          LMN1=LMN1+1
!!$          AUX2 = EIGENVAL(LMN1) - EB(LN1)
!!$          WRITE(NFILO,FMT='(6F12.5)')EIGENVAL(LMN1),EIGENVAL(LMN1)/EV&
!!$            &,EB(LN1),EB(LN1)/EV,AUX2,AUX2/EV
!!$        ENDDO
!!$      ENDDO
!!$
!      CALL FILEHANDLER$UNIT('PROT',NFILO)
!      CALL CORE$REPORT(NFILO)
!     DO I=1,LMNX
!       AUX2 = EIGENVAL(I) - EB(I)
!       WRITE(NFILO,FMT='(6F12.5)')EIGENVAL(I),EIGENVAL(I)/EV,EB(I)&
!         &,EB(I)/EV,AUX2,AUX2/EV
!     ENDDO

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
