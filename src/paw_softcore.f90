MODULE CORE_MODULE
!**                                                                           **
!**  THE COMPLEX POINTER LIST THIS IS NEEDED FOR THE PARALLELIZATION          **
!**  BECAUSE EVERY NODE WILL CALCULATE THE CORE LEVELS ONLY FOR ITS           **
!**  OWN SET OF ATOMS AND ONLY THE FIRST NODE CAN WRITE                       **
!**                                                                           **
!**  COMMUNIATION IN CORE$REPORT NOT DONE YET.                                **
!**                                                                           **
TYPE CORESHIFT_TYPE
INTEGER(4)           :: IAT        !ATOM INDEX
INTEGER(4)           :: N          !#(CORE WAVE FUNCTIONS PER SPIN (L,M,N))
CHARACTER(8),POINTER :: TYPE(:)    !ONE OF 'S','P','D','F','?'
REAL(8)     ,POINTER :: E(:)       !CRYSTAL CORE ENERGY LEVELS
REAL(8)     ,POINTER :: EATOM(:)   !ATOMIC ENERGY EIGENVALUE
TYPE(CORESHIFT_TYPE),POINTER :: NEXT
END TYPE CORESHIFT_TYPE
LOGICAL(4)                  :: TCORESHIFTS=.TRUE.
LOGICAL(4)                  :: DEFAULT=.FALSE.
INTEGER(4)                  :: NATOMS=0
CHARACTER(32),ALLOCATABLE   :: ATOMS(:)
TYPE(CORESHIFT_TYPE),TARGET :: FIRST
TYPE(CORESHIFT_TYPE),POINTER:: THIS
LOGICAL(4),SAVE             :: TINI=.FALSE.
CHARACTER(32)  ,PARAMETER   :: FILEHFID='HFDAT'
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXTRAPOLATE(NP,RI,FI_,R0,F0)
!     **************************************************************************
!     **  POLYNOMIAL EXTRAPOLATION OF ORDER NP FROM NP POINTS (RI,FI_)        **
!     **  TO THE POINT (R0,F0)                                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP
      REAL(8)   ,INTENT(IN) :: RI(NP)
      REAL(8)   ,INTENT(IN) :: FI_(NP)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: FI(NP)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J,IP
!     **************************************************************************
      FI(:)=FI_(:)
      F0=0.D0
      DO I=1,NP
        SVAR=1.D0
        DO J=1,I-1
          SVAR=SVAR*(R0-RI(J))/(RI(I)-RI(J))
        ENDDO
        F0=F0+SVAR*FI(I)
        DO IP=I+1,NP
          SVAR=1.D0
          DO J=1,I-1
            SVAR=SVAR*(RI(IP)-RI(J))/(RI(I)-RI(J))
          ENDDO
          FI(IP)=FI(IP)-SVAR*FI(I)
        ENDDO
        FI(I)=0.D0
      ENDDO
      RETURN
      END SUBROUTINE EXTRAPOLATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DTOXYZ(VLM_,V2)
!     ****************************************************************
!     **                                                            **
!     **  GIVEN THE VALUE OF F(R)/R**2 FOR THE 5 D-TYPE RADIAL      **
!     **  FUNCTIONS, CALCULATED THE MATRIX OF THE SECOND            **
!     **  DERIVATIVES AT THE ORIGIN                                 **
!     **                                                            **
!     **  THE FIRST NINE REAL SPHERICAL HARMONICS ARE:              **
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2    **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2    **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2    **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2    **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2    **
!     ****************************************************************
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: VLM_(5)
      REAL(8), INTENT(OUT):: V2(3,3)
      REAL(8)             :: VLM(5)
      REAL(8)             :: PI
      REAL(8)             :: SQ15
      REAL(8)             :: SQ5
      REAL(8)             :: SQ60
      REAL(8)             :: SQ16PI
!     ****************************************************************
      PI=4.D0*ATAN(1.D0)
      SQ15=SQRT(15.D0)
      SQ5=SQRT(5.D0)
      SQ60=SQRT(60.D0)
      SQ16PI=SQRT(16.D0*PI)
      VLM(1)=VLM_(1)*SQ15/SQ16PI
      VLM(2)=VLM_(2)*SQ60/SQ16PI
      VLM(3)=VLM_(3)*SQ5/SQ16PI
      VLM(4)=VLM_(4)*SQ60/SQ16PI
      VLM(5)=VLM_(5)*SQ60/SQ16PI
      V2(1,1)=+2.D0*VLM(1)-2.D0*VLM(3)
      V2(1,2)=+VLM(5)
      V2(1,3)=+VLM(2)
      V2(2,2)=-2.D0*VLM(1)-2.D0*VLM(3)
      V2(2,3)=+VLM(4)
      V2(3,3)=+4.D0*VLM(3)
      V2(2,1)=V2(1,2)
      V2(3,1)=V2(1,3)
      V2(3,2)=V2(2,3)
      RETURN
      END SUBROUTINE DTOXYZ
END MODULE CORE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE$SETL4(ID,VAL)
      USE CORE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE$SETCHA(ID,LEN,VAL)
      USE CORE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
      INTEGER(4)              :: I
!     **************************************************************************
      IF(ID.EQ.'ATOMS') THEN
!       == LIST OF ATOMS FOR WHICH THE CORE STATES SHALL BE EVALUATED ==========
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
!     **************************************************************************
      IF(.NOT.TINI) RETURN
                             CALL TRACE$PUSH('CORE$REPORT')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL CONSTANTS('EV',EV)

!     == CREATE A MAPPING ARRAY FROM ATOMS TO TASKS FOR EACH NODE.            ==
!     == ATOMS NOT PRESENT ON THIS TASK OBTAIN 0 FOR THE ATOM NUMBER          ==
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(TASKARR(NAT))
      TASKARR(:)=0   ! TASK ID FOR EACH ATOM
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
!     == COMMUNICATE A UNIQUE MAPPING KNOWN TO ALL TASKS =======================
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

!
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
!       == TASK 1 RECEIVED DATA ================================================
!
!       == PRINT INFORMATION ON THE FIRST TASK =================================
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
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE CORE_CORESHIFTS(IAT,ISP,GID,NR,LMRXX,NDIMD,AEPOT)
!!$!     **************************************************************************
!!$!     **  CALCULATES THE EIGENVALUES OF CORE HAMILTONIAN                      **
!!$!     **                                                                      **
!!$!     **  IS CALLED FROM PAW_AUGMENTATION, WHICH EXECUTES ONLY ON ONE TASK    **
!!$!     **                                                                      **
!!$!     **                                                                      **
!!$!     **************************************************************************
!!$      USE CORE_MODULE, ONLY : CORESHIFT_TYPE &
!!$     &                       ,TCORESHIFTS &
!!$     &                       ,DEFAULT &
!!$     &                       ,TINI &
!!$     &                       ,THIS,FIRST &
!!$     &                       ,NATOMS &
!!$     &                       ,ATOMS 
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: IAT     ! ATOM INDEX
!!$      INTEGER(4),INTENT(IN) :: ISP     ! ATOM TYPE
!!$      INTEGER(4),INTENT(IN) :: GID     ! GRID ID
!!$      INTEGER(4),INTENT(IN) :: NR      ! #(GRID POINTS)
!!$      INTEGER(4),INTENT(IN) :: NDIMD   ! #(DENSITY SPIN COMPONENTS)
!!$      INTEGER(4),INTENT(IN) :: LMRXX
!!$      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRXX,NDIMD) ! 1C-AE POTENTIAL
!!$      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
!!$      REAL(8)               :: ATPOT(NR)  !RADIAL ATOM POTENTIAL 
!!$      REAL(8)   ,ALLOCATABLE:: AEPOT1(:,:,:)!(NR,LMRX,NDIMD)
!!$      INTEGER(4)            :: NB         !#(ATOMIC CORE AND VALENCE STATES)
!!$      INTEGER(4)            :: NC         !#(ATOMIC CORE STATES)
!!$      INTEGER(4),ALLOCATABLE:: LB(:)      !MAIN ANGULAR MOMENTUM
!!$!     REAL(8)   ,ALLOCATABLE:: FB(:)  
!!$      REAL(8)   ,ALLOCATABLE:: EB(:)      !ENERGY LEVEL
!!$      REAL(8)   ,ALLOCATABLE:: AEPSI(:,:) !AE WAVE FUNCTION
!!$      REAL(8)               :: R(NR)      !RADIAL GRID
!!$      CHARACTER(32)         :: NAME
!!$      INTEGER(4)            :: LMNX
!!$      INTEGER(4)            :: LMRX
!!$
!!$      INTEGER(4)            :: I
!!$      REAL(8)               :: AUX2
!!$
!!$      INTEGER(4)            :: LMN1,LMN2
!!$      INTEGER(4)            :: LN1,LN2
!!$      INTEGER(4)            :: LM1,LM2,LM3
!!$      INTEGER(4)            :: L1,L2,IM1,IM2
!!$      REAL(8)               :: AEDMU(NR)
!!$      REAL(8)               :: DWORK1(NR)
!!$      REAL(8)               :: CG          !GAUNT COEFFICIENT
!!$      REAL(8)               :: SVAR
!!$
!!$      REAL(8)   ,ALLOCATABLE:: HAMIL(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: EIGENVAL(:)
!!$      REAL(8)   ,ALLOCATABLE:: EIGENVEC(:,:)
!!$      INTEGER(4)            :: NFILO
!!$      REAL(8)               :: EV             ! ELECTRON VOLT
!!$      LOGICAL(4)            :: TCHK
!!$!     **************************************************************************
!!$      IF(.NOT.TCORESHIFTS) RETURN
!!$      CALL RADIAL$R(GID,NR,R)
!!$      CALL ATOMLIST$GETCH('NAME',IAT,NAME)
!!$!     == CHECK WHETHER CORE LEVELS SHALL BE CALCULATED FOR THIS ATOM (IAT)  ====
!!$      TCHK=DEFAULT
!!$      DO I=1,NATOMS
!!$        IF(NAME.EQ.ATOMS(I)) THEN
!!$          TCHK=.NOT.DEFAULT
!!$        END IF
!!$      ENDDO
!!$      IF(.NOT.TCHK) RETURN
!!$!
!!$!     ==========================================================================
!!$!     ==  SELECT THE PROPER ENTRY IN THE TABLE                                ==
!!$!     ==========================================================================
!!$      THIS=>FIRST
!!$      IF(.NOT.TINI) THEN
!!$        TINI=.TRUE.
!!$        THIS%IAT=IAT
!!$        THIS%N=0
!!$        NULLIFY(THIS%E)
!!$        NULLIFY(THIS%EATOM)
!!$        NULLIFY(THIS%TYPE)
!!$        NULLIFY(THIS%NEXT)
!!$      ELSE
!!$        DO WHILE (THIS%IAT.NE.IAT) 
!!$          IF(ASSOCIATED(THIS%NEXT)) THEN
!!$            THIS=>THIS%NEXT
!!$          ELSE
!!$            ALLOCATE(THIS%NEXT)
!!$            THIS=>THIS%NEXT
!!$            NULLIFY(THIS%NEXT)
!!$            THIS%IAT=IAT
!!$            THIS%N=0
!!$            NULLIFY(THIS%E)
!!$            NULLIFY(THIS%EATOM)
!!$            NULLIFY(THIS%TYPE)
!!$          END IF
!!$        ENDDO
!!$      END IF
!!$!     == THE POINTER "THIS" REFERS NOW TO THE ATOM AT HAND =====================
!!$!
!!$!     ==========================================================================
!!$!     ==  COLLECT ATOMIC CORE WAVE FUNCTIONS AEPSI FROM SETUP OBJECT          ==
!!$!     ==========================================================================
!!$      CALL SETUP$ISELECT(ISP)
!!$      CALL SETUP$GETI4('NB',NB)
!!$      CALL SETUP$GETI4('NC',NC)
!!$      IF (NC.EQ.0) THEN
!!$        CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$        CALL REPORT$TITLE(NFILO,'NO CORE STATES FOR ATOM '//TRIM(NAME))
!!$        RETURN
!!$      END IF
!!$      ALLOCATE(LB(NB))
!!$      ALLOCATE(EB(NB))
!!$      ALLOCATE(AEPSI(NR,NB))
!!$      CALL SETUP$GETR8A('AEPOT',NR,ATPOT)    !POTENTIAL OF THE ATOM
!!$      CALL SETUP$GETI4A('LB',NB,LB)          !MAIN ANGULAR MOMENTUM
!!$      CALL SETUP$GETR8A('EB',NB,EB)          !ATOMIC ENERGY LEVEL
!!$      CALL SETUP$GETR8A('AEPSI',NR*NB,AEPSI) !ATOMIC WAVE FUNCTION
!!$      CALL SETUP$UNSELECT()
!!$!
!!$!     ==========================================================================
!!$!     ==  CONSTANTS                                                           ==
!!$!     ==========================================================================
!!$      LMNX=0
!!$      LMRX=0
!!$      DO I=1,NC
!!$        LMNX=LMNX+2*LB(I)+1
!!$        LMRX=MAX(LMRX,(2*LB(I)+1)**2)
!!$      ENDDO
!!$      LMRX=MIN(LMRX,LMRXX)
!!$      IF(THIS%N.NE.LMNX) THEN
!!$        THIS%N=LMNX
!!$        ALLOCATE(THIS%TYPE(LMNX))
!!$        ALLOCATE(THIS%E(LMNX))
!!$        ALLOCATE(THIS%EATOM(LMNX))
!!$      END IF
!!$
!!$!     ==========================================================================
!!$!     == SUBTRACTS ATOMIC AE POTENTIAL FROM AE TOTAL POTENTIAL                ==
!!$!     ==========================================================================
!!$      ALLOCATE(AEPOT1(NR,LMRXX,NDIMD))
!!$      AEPOT1(:,LMRX,:)=AEPOT(:,LMRX,:)
!!$      AEPOT1(:,1)=AEPOT(:,1)-ATPOT(:)
!!$!
!!$!     ==========================================================================
!!$!     ==   CALCULATE HAMILTONIAN                                              ==
!!$!     ==========================================================================
!!$      ALLOCATE(HAMIL(LMNX,LMNX))      
!!$      HAMIL(:,:)=0.D0
!!$!
!!$      LMN1=0
!!$      DO LN1=1,NC
!!$        L1=LB(LN1)
!!$        DO IM1=1,2*L1+1
!!$          LMN1=LMN1+1
!!$          LM1=L1**2+IM1
!!$!
!!$          LMN2=0
!!$          DO LN2=1,NC
!!$            L2=LB(LN2)
!!$            DO IM2=1,2*L2+1
!!$              LMN2=LMN2+1
!!$              LM2=L2**2+IM2
!!$!
!!$              IF(LMN1.EQ.LMN2) THEN
!!$                HAMIL(LMN1,LMN2)=EB(LN1)
!!$              END IF
!!$!
!!$              AEDMU(:)=0.D0
!!$              DO LM3=1,LMRX
!!$                CALL CLEBSCH(LM1,LM2,LM3,CG)
!!$                IF(CG.NE.0.D0) THEN
!!$                  AEDMU(:)=AEDMU(:)+CG*AEPOT1(:,LM3)
!!$                END IF
!!$              ENDDO
!!$!
!!$              DWORK1(:)=(AEDMU(:)*AEPSI(:,LN1)*AEPSI(:,LN2))*R(:)**2
!!$              CALL RADIAL$INTEGRAL(GID,NR,DWORK1,SVAR)
!!$              HAMIL(LMN1,LMN2)=HAMIL(LMN1,LMN2)+SVAR
!!$!
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      DEALLOCATE(AEPOT1)
!!$!
!!$!     ==========================================================================
!!$!     ==  DIAGONALIZATION OF THE HAMILTONIAN                                  ==
!!$!     ==========================================================================
!!$      ALLOCATE(EIGENVAL(LMNX))
!!$      ALLOCATE(EIGENVEC(LMNX,LMNX))
!!$      CALL LIB$DIAGR8(LMNX,HAMIL,EIGENVAL,EIGENVEC)
!!$      DEALLOCATE(EIGENVEC)
!!$!
!!$!     ==========================================================================
!!$!     ==  WRITE INTO TABLE                                                    ==
!!$!     ==========================================================================
!!$      LMN1=0
!!$      DO LN1=1,NC
!!$        L1=LB(LN1)
!!$        DO IM1=1,2*L1+1
!!$          LMN1=LMN1+1
!!$          THIS%E(LMN1)=EIGENVAL(LMN1)
!!$          THIS%EATOM(LMN1)=EB(LN1)
!!$!
!!$!         == MAIN QUANTUM NUMBER I =============================================
!!$          I=L1
!!$          DO LN2=1,LN1
!!$            IF(LB(LN2).NE.L1) CYCLE
!!$            I=I+1
!!$          ENDDO
!!$!
!!$!         == COMPOSE STRING 1S,2S,2P,3S,3P,3D,... ==============================
!!$          WRITE(THIS%TYPE(LMN1),FMT='(I2)')I
!!$          IF(L1.EQ.0) THEN
!!$            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'S'
!!$          ELSE IF(L1.EQ.1) THEN
!!$            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'P'
!!$          ELSE IF(L1.EQ.2) THEN
!!$            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'D'
!!$          ELSE IF(L1.EQ.3) THEN
!!$            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'F'
!!$          ELSE
!!$            THIS%TYPE(LMN1)=TRIM(THIS%TYPE(LMN1))//'?'
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$      
!!$
!!$!     ==========================================================================
!!$!     ==  CLOSE DOWN                                                          ==
!!$!     ==========================================================================
!!$      DEALLOCATE(LB)
!!$!     DEALLOCATE(FB)
!!$      DEALLOCATE(EB)
!!$      DEALLOCATE(AEPSI)
!!$      DEALLOCATE(HAMIL)
!!$      DEALLOCATE(EIGENVAL)
!!$
!!$      RETURN
!!$      END
!
! SANTOS040617 END
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_BEYONDFROZENCORE(IAT,ISP,GID,NR,LMRXX,NDIMD,AEPOT,AERHO &
     &                                ,LMNX,DENMAT,SELFENERGY)
!     **************************************************************************
!     **  WORK OUT HAMILTONIAN BETWEEN FROZEN CORE STATES AND ALL-ELECTRON    **
!     **  PARTIAL WAVES                                                       **
!     **                                                                      **
!     **  DENMAT IS IN (T), (T,SZ), OR (T,SX,SY,SZ) REPRESENTATION              **
!     **                                                                      **
!     **  IS CALLED FROM PAW_AUGMENTATION, WHICH EXECUTES ONLY ON ONE TASK    **
!     **                                                                      **
! TODO: THE HVV LACKS THE ATOMIC POTENTIAL AND THE KINETIC ENERGY
! TODO: INCLUDE SPIN-ORBIT COUPLING
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT     ! ATOM INDEX
      INTEGER(4),INTENT(IN) :: ISP     ! ATOM TYPE
      INTEGER(4),INTENT(IN) :: GID     ! GRID ID
      INTEGER(4),INTENT(IN) :: NR      ! #(GRID POINTS)
      INTEGER(4),INTENT(IN) :: LMRXX
      INTEGER(4),INTENT(IN) :: NDIMD   ! #(DENSITY SPIN COMPONENTS)
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRXX,NDIMD) ! 1C-AE POTENTIAL
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRXX,NDIMD) ! 1C-AE DENSITY
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(OUT):: SELFENERGY(LMNX,LMNX,NDIMD) 
      REAL(8)   ,INTENT(OUT):: DENMAT(LMNX,LMNX,NDIMD) !
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  !COMPLEX SQRT(-1)
      COMPLEX(8),PARAMETER  :: CZERO=(0.D0,0.D0)  !COMPLEX ZERO
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: EPSILONNU=0.D0
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      LOGICAL(8),PARAMETER  :: TPR=.TRUE.  ! PRINT INFORMATION
      LOGICAL(8),PARAMETER  :: TSO=.FALSE. ! INCLUDE SPIN-ORBIT COUPLING
      LOGICAL(8),PARAMETER  :: TNOSELFENERGY=.TRUE. 
      REAL(8)               :: ATPOT(NR)  !RADIAL ATOM POTENTIAL 
      REAL(8)   ,ALLOCATABLE:: AEPOT1(:,:,:)!(NR,LMRXX,NDIMD)
      INTEGER(4)            :: LNX        !
      INTEGER(4),ALLOCATABLE:: LOX(:)     !MAIN ANGULAR MOMENTUM OF PARTIAL W.
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:) !(NR,LNX) AE PARTIAL WAVES
      INTEGER(4)            :: NB         !#(ATOMIC CORE AND VALENCE STATES)
      INTEGER(4)            :: NC         !#(ATOMIC CORE STATES)
      INTEGER(4)            :: LMNCX
      INTEGER(4)            :: LMNTX
      INTEGER(4),ALLOCATABLE:: LB(:)      !MAIN ANGULAR MOMENTUM
      REAL(8)   ,ALLOCATABLE:: EB(:)      !ENERGY LEVEL
      REAL(8)   ,ALLOCATABLE:: AEPSI(:,:) !(NR,NB) AE WAVE FUNCTION
      REAL(8)               :: R(NR)      !RADIAL GRID
      CHARACTER(32)         :: NAME
      INTEGER(4)            :: LMRX

      INTEGER(4)            :: IC

      INTEGER(4)            :: LN
      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LM1,LM2,LM3
      INTEGER(4)            :: L1,L2,IM1,IM2
      INTEGER(4)            :: I1A,I1B,I2A,I2B,I1,I2,IS1,IS2
      INTEGER(4)            :: IDIMD
      REAL(8)               :: AEDMU(NR,NDIMD)
      REAL(8)               :: DREL(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: DWORK1(NR) ! SHALL BE REPLACED BY AUX
      REAL(8)               :: CG          !GAUNT COEFFICIENT
      REAL(8)               :: SVAR
      REAL(8)               :: X(NDIMD)
      COMPLEX(8)            :: CMAT22(2,2)
      COMPLEX(8),ALLOCATABLE:: HT(:,:) ! CORE AND VALENCE HAMILTONIAN
      COMPLEX(8),ALLOCATABLE:: HCC(:,:)
      COMPLEX(8),ALLOCATABLE:: HCV(:,:)
      COMPLEX(8),ALLOCATABLE:: HVV(:,:)
      COMPLEX(8),ALLOCATABLE:: OCC(:,:)
      COMPLEX(8),ALLOCATABLE:: OCV(:,:)
      COMPLEX(8),ALLOCATABLE:: OVV(:,:)
      COMPLEX(8),ALLOCATABLE:: OT(:,:)   ! CORE AND VALENCE OVERLAP
      COMPLEX(8),ALLOCATABLE:: RHO1VV(:,:)   
      REAL(8)   ,ALLOCATABLE:: EIGENVAL(:)
      COMPLEX(8),ALLOCATABLE:: EIGENVEC(:,:)
      COMPLEX(8),ALLOCATABLE:: SELFENERGY1(:,:)
      COMPLEX(8),ALLOCATABLE:: cmat2(:,:) !(2*lmncx,2*lmncx)
      COMPLEX(8),ALLOCATABLE:: GREEN(:,:)
      INTEGER(4)             :: LMX
      COMPLEX(8),ALLOCATABLE :: CLS(:,:,:,:)
      INTEGER(4)            :: NFILO
      REAL(8)               :: E1,ERHO
      REAL(8)               :: Ev
      CHARACTER(128)        :: FMTSTRING='(40("."),":",F20.10,T1,A)'
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      CALL ATOMLIST$GETCH('NAME',IAT,NAME)
WRITE(*,FMT='(80("="),T20,A)')'CORE_BEYONDFROZENCORE FOR ATOM '//TRIM(NAME)
      CALL AUGMENTATION_FROZENCORENERGY(ISP,SVAR)
WRITE(*,FMT='("CORE ENERGY=",10F20.8)')SVAR
      CALL CORE_HYPERFINE(IAT,GID,NR,NDIMD,LMRXX,AERHO)
!
!     ==========================================================================
!     ==  COLLECT ATOMIC CORE WAVE FUNCTIONS AEPSI FROM SETUP OBJECT          ==
!     ==========================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETI4('NB',NB)
      CALL SETUP$GETI4('NC',NC)
      IF (NC.EQ.0) THEN
!       == NO CORE STATES, NOTHING TO DO, RETURN ===============================
        CALL FILEHANDLER$UNIT('PROT',NFILO)
!        CALL REPORT$TITLE(NFILO,'NO CORE STATES FOR ATOM '//TRIM(NAME))
        CALL SETUP$UNSELECT()
        RETURN
      END IF
      ALLOCATE(LB(NB))
      ALLOCATE(EB(NB))
      ALLOCATE(AEPSI(NR,NB))
      CALL SETUP$GETR8A('AEPOT',NR,ATPOT)    !POTENTIAL OF THE ATOM
      CALL SETUP$GETI4A('LB',NB,LB)          !MAIN ANGULAR MOMENTUM
      CALL SETUP$GETR8A('EB',NB,EB)          !ATOMIC ENERGY LEVEL
      CALL SETUP$GETR8A('AEPSI',NR*NB,AEPSI) !ATOMIC WAVE FUNCTION
      CALL SETUP$UNSELECT()

WRITE(*,FMT='("LNX=",10I5)')LNX
WRITE(*,FMT='("LOX=",10I5)')LOX
WRITE(*,FMT='("NB=",10I5)')NB
WRITE(*,FMT='("NC=",10I5)')NC
WRITE(*,FMT='("LB=",10I5)')LB(:NC)
WRITE(*,FMT='("EB=",5F12.3)')EB(:NC)
!
!     ==========================================================================
!     ==  CONSTANTS                                                           ==
!     ==========================================================================
      LMRX=0
      LMNCX=0
      DO IC=1,NC
        LMNCX=LMNCX+2*LB(IC)+1
        LMRX=MAX(LMRX,(2*LB(IC)+1)**2)
      ENDDO
      DO LN=1,LNX
        LMRX=MAX(LMRX,(2*LOX(LN)+1)**2)
      ENDDO
      LMRX=MIN(LMRX,LMRXX)
      IF(LMNX.NE.SUM(2*LOX(:)+1)) THEN
        CALL ERROR$MSG('CONSISTENCY CHECK FOR LMNX FAILED')
        CALL ERROR$MSG('LMNX MUST BE EQUAL TO SUM(2*LOX(:)+1)')
        CALL ERROR$I4VAL('LMNX',LMNX)
        CALL ERROR$I4VAL('SUM(2*LOX+1)',SUM(2*LOX(:)+1))
        CALL ERROR$STOP('CORE_BEYONDFROZENCORE')
      END IF
WRITE(*,FMT='("NDIMD=",10I5)')NDIMD
WRITE(*,FMT='("LMNX=",10I5)')LMNX
WRITE(*,FMT='("LMNCX=",10I5)')LMNCX
WRITE(*,FMT='("LMRX=",10I5)')LMRX
!
!     == L*S MATRIX ELEMENTS FOR SPIN ORBIT COUPLING ===========================
      LMX=(MAX(MAXVAL(LOX),MAXVAL(LB(:NC)))+1)**2
      ALLOCATE(CLS(LMX,2,LMX,2))
      CALL SCHROEDINGER_LS(LMX,CLS)
      CLS=0.5D0*CLS   ! CONVERT L*SIGMA INTO L*S
!
!     ==========================================================================
!     == SUBTRACT ATOMIC AE POTENTIAL FROM AE TOTAL POTENTIAL                 ==
!     ==========================================================================
      ALLOCATE(AEPOT1(NR,LMRX,NDIMD))
      AEPOT1(:,:,:)=AEPOT(:,:LMRX,:)
      AEPOT1(:,1,1)=AEPOT1(:,1,1)-ATPOT(:)
!
!     ==========================================================================
!     ==   CALCULATE CORE-CORE HAMILTONIAN                                    ==
!     ==========================================================================
      ALLOCATE(HCC(2*LMNCX,2*LMNCX))      
      ALLOCATE(OCC(2*LMNCX,2*LMNCX))      
      HCC(:,:)=CZERO
      OCC(:,:)=CZERO
!
      LMN1=0
      DO LN1=1,NC
        L1=LB(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LM1=L1**2+IM1
          I1A=2*(LMN1-1)+1
          I1B=I1A+1
!
          LMN2=0
          DO LN2=1,NC
            L2=LB(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
              I2A=2*(LMN2-1)+1
              I2B=I2A+1
!
              AEDMU(:,:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  AEDMU(:,:)=AEDMU(:,:)+CG*AEPOT1(:,LM3,:)
                END IF
              ENDDO
!
!             == X(IDIMD)=<PSI(LMN1)|AEPOT(IDIMD)|PSI(LMN2)>
              DO IDIMD=1,NDIMD
                AUX(:)=AEDMU(:,IDIMD)*AEPSI(:,LN1)*AEPSI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,AUX,X(IDIMD))
              ENDDO
PRINT*,'X ',LMN1,LMN2,X
              CALL CORE_RHOMAG2UPDN(NDIMD,X,CMAT22)
              HCC(I1A:I1B,I2A:I2B)=HCC(I1A:I1B,I2A:I2B)+CMAT22(:,:)

!OCC IS THE UNIT MATRIX BY CONSTRUCTION
              IF(LM1.EQ.LM2) THEN
                DWORK1(:)=AEPSI(:,LN1)*AEPSI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,SVAR)
                OCC(I1A,I2A)=OCC(I1A,I2A)+SVAR
                OCC(I1B,I2B)=OCC(I1B,I2B)+SVAR
                IF(LMN1.EQ.LMN2) THEN
                  HCC(I1A,I2A)=HCC(I1A,I2A)+SVAR*EB(LN1)
                  HCC(I1B,I2B)=HCC(I1B,I2B)+SVAR*EB(LN1)
                END IF
              END IF
!
!             == SPIN ORBIT MATRIX ELEMENTS ====================================
              IF(TSO) THEN
                SVAR=0.5D0*(EB(LN1)+EB(LN2))
                CALL SCHROEDINGER$DREL(GID,NR,AEPOT(:,1,1),SVAR,DREL)
                CALL RADIAL$DERIVE(GID,NR,DREL,AUX)
                AUX(:)=AUX(:)*R(:)*AEPSI(:,LN1)*AEPSI(:,LN2)  ! (AUX/R)*R**2
                CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
                DO IS1=1,2
                  DO IS2=1,2
                    I1=I1A-1+IS1
                    I2=I2A-1+IS2
                    HCC(I1,I2)=HCC(I1,I2)+CLS(LM1,IS1,LM2,IS2)*SVAR
                  ENDDO
                ENDDO
              END IF !TSO
!
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(TPR) THEN
        WRITE(*,FMT='(80("="),T20,A)')' CORE HAMILTONIAN '
        DO LMN1=1,2*LMNCX
          DO LMN2=LMN1,2*LMNCX
            IF(ABS(HCC(LMN1,LMN2)).GT.TOL) THEN
              WRITE(*,FMT='("HCC(",I3,",",I3,")= ",2F12.5)') &
    &                LMN1,LMN2,HCC(LMN1,LMN2)
            END IF      
          ENDDO
        ENDDO
!
!       == TEST ORTHONORMALITY OF CORE STATES ==================================
        DO I1=1,LMNCX
          IF(ABS(OCC(I1,I1)-(1.D0,0.D0)).GT.TOL.OR. &
     &        MAXVAL(ABS(OCC(I1,I1+1:))).GT.TOL.OR. &
     &        MAXVAL(ABS(OCC(I1+1:,I1))).GT.TOL) THEN
            WRITE(*,FMT='(80("="),T20,A)')' CORE OVERLAP '
            WRITE(*,FMT='(I5,"     ",100I10)')0,(LMN2,LMN2=1,LMNCX)
            DO LMN1=1,2*LMNCX
              WRITE(*,FMT='(I5,"OCC: ",100F10.5)')LMN1,REAL(OCC(LMN1,:))
            ENDDO
            CALL ERROR$MSG('OCC IS NOT UNITY')
            CALL ERROR$L4VAL('TSO',TSO)
            CALL ERROR$I4VAL('I1',I1)
            CALL ERROR$C8VAL('OCC(I1,I1)',OCC(I1,I1))
            CALL ERROR$R8VAL('MAX|OCC(I1,I1+1:)|',MAXVAL(ABS(OCC(I1,I1+1:))))
            CALL ERROR$R8VAL('MAX|OCC(I1+1:,I1)|',MAXVAL(ABS(OCC(I1+1:,I1))))
            CALL ERROR$STOP('CORE_BEYONDFROZENCORE')
          END IF     
        ENDDO
      END IF !TPR
!
!     ==========================================================================
!     ==   CALCULATE CORE-VALENCE HAMILTONIAN                                 ==
!     ==========================================================================
      ALLOCATE(HCV(2*LMNCX,2*LMNX))      
      ALLOCATE(OCV(2*LMNCX,2*LMNX))      
      HCV(:,:)=CZERO
      OCV(:,:)=CZERO
      LMN1=0
      DO LN1=1,NC
        L1=LB(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LM1=L1**2+IM1
          I1A=2*(LMN1-1)+1
          I1B=I1A+1
!
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
              I2A=2*(LMN2-1)+1
              I2B=I2A+1
!
              AEDMU(:,:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  AEDMU(:,:)=AEDMU(:,:)+CG*AEPOT1(:,LM3,:)
                END IF
              ENDDO
!
!             == X(IDIMD)=<PSI(LMN1)|AEPOT(IDIMD)|PSI(LMN2)>
              DO IDIMD=1,NDIMD
                DWORK1(:)=AEDMU(:,IDIMD)*AEPSI(:,LN1)*AEPHI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,X(IDIMD))
              ENDDO
              CALL CORE_RHOMAG2UPDN(NDIMD,X,CMAT22)
              HCV(I1A:I1B,I2A:I2B)=HCV(I1A:I1B,I2A:I2B)+CMAT22(:,:)
! 
!OCV IS ZERO BY CONSTRUCTION
              IF(LM1.EQ.LM2) THEN
                DWORK1(:)=AEPSI(:,LN1)*AEPHI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,X(1))
                CALL CORE_RHOMAG2UPDN(1,X,CMAT22)
                OCV(I1A:I1B,I2A:I2B)=OCV(I1A:I1B,I2A:I2B)+CMAT22(:,:)
                HCV(I1A:I1B,I2A:I2B)=HCV(I1A:I1B,I2A:I2B)+CMAT22(:,:)*EB(LN1)
              END IF
!
!             == SPIN ORBIT MATRIX ELEMENTS ====================================
              IF(TSO) THEN
                SVAR=0.5D0*(EB(LN1)+0.D0)
                CALL SCHROEDINGER$DREL(GID,NR,AEPOT(:,1,1),SVAR,DREL)
                CALL RADIAL$DERIVE(GID,NR,DREL,AUX)
                AUX(:)=AUX(:)*R(:)*AEPSI(:,LN1)*AEPHI(:,LN2)  ! (AUX/R)*R**2
                CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
                DO IS1=1,2
                  DO IS2=1,2
                    I1=I1A-1+IS1
                    I2=I2A-1+IS2
                    HCV(I1,I2)=HCV(I1,I2)+CLS(LM1,IS1,LM2,IS2)*SVAR
                  ENDDO
                ENDDO
              END IF !TSO 

            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      IF(TPR) THEN
        WRITE(*,FMT='(80("="),T20,A)')' CORE-VALENCE HAMILTONIAN '
        DO LMN1=1,2*LMNCX
          DO LMN2=1,2*LMNX
            IF(ABS(HCV(LMN1,LMN2)).GT.1.D-5) THEN
              WRITE(*,FMT='("HCV(",I3,",",I3,")= ",2F12.5)') &
                   LMN1,LMN2,HCV(LMN1,LMN2)
            END IF      
          ENDDO
        ENDDO
        WRITE(*,FMT='(80("="),T20,A)')' CORE-VALENCE OVERLAP '
        DO LMN1=1,2*LMNCX
          DO LMN2=1,2*LMNX
            IF(ABS(OCV(LMN1,LMN2)).GT.1.D-5) THEN
              WRITE(*,FMT='("OCV(",I3,",",I3,")= ",2F12.5)') &
                   LMN1,LMN2,OCV(LMN1,LMN2)
            END IF      
          ENDDO
        ENDDO
!
!       == TEST CORE-VALENCE ORTHOGONALITY =====================================
        IF(MAXVAL(ABS(OCV)).GT.TOL) THEN
            CALL ERROR$MSG('OCV IS NOT UNITY')
            CALL ERROR$L4VAL('TSO',TSO)
            CALL ERROR$R8VAL('MAX|OCV|',MAXVAL(ABS(OCV)))
            CALL ERROR$STOP('CORE_BEYONDFROZENCORE')
        END IF
      END IF

!!$WRITE(*,FMT='(80("="),T20,A)')' CORE-VALENCE OVERLAP '
!!$WRITE(*,FMT='(I5,"     ",100I10)')0,(LMN2,LMN2=1,2*LMNX)
!!$DO LMN1=1,2*LMNCX
!!$  WRITE(*,FMT='(I5,"OCV: ",100F10.5)')LMN1,REAL(OCV(LMN1,:))
!!$ENDDO
!!$WRITE(*,FMT='(80("="),T20,A)')' CORE-VALENCE HAMILTONIAN '
!!$WRITE(*,FMT='(I5,"     ",100I10)')0,(LMN2,LMN2=1,2*LMNX)
!!$DO LMN1=1,2*LMNCX
!!$  WRITE(*,FMT='(I5,"HCV: ",100F10.2)')LMN1,REAL(HCV(LMN1,:))
!!$ENDDO

!     ==========================================================================
!     ==  VALENCE HAMILTONIAN AND OVERLAP
!     ==========================================================================
      ALLOCATE(HVV(2*LMNX,2*LMNX))      
      ALLOCATE(OVV(2*LMNX,2*LMNX))      
      HVV(:,:)=CZERO
      OVV(:,:)=CZERO
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LM1=L1**2+IM1
          I1A=2*(LMN1-1)+1
          I1B=I1A+1
!
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
              I2A=2*(LMN2-1)+1
              I2B=I2A+1
!
!
              AEDMU(:,:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  AEDMU(:,:)=AEDMU(:,:)+CG*AEPOT1(:,LM3,:)
                END IF
              ENDDO
!
              DO IDIMD=1,NDIMD
                DWORK1(:)=AEDMU(:,IDIMD)*AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,X(IDIMD))
              ENDDO
!             == TRANSFORM FROM (T,X,Y,Z) TO (UP,DN) REPRESENTATION ============
              CALL CORE_RHOMAG2UPDN(NDIMD,X,CMAT22)
              HVV(I1A:I1B,I2A:I2B)=HVV(I1A:I1B,I2A:I2B)+CMAT22(:,:)

!OCV IS ZERO BY CONSTRUCTION
              IF(LM1.EQ.LM2) THEN
                DWORK1(:)=AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,X(1))
                CALL CORE_RHOMAG2UPDN(1,X,CMAT22)
                OVV(I1A:I1B,I2A:I2B)=OVV(I1A:I1B,I2A:I2B)+CMAT22(:,:)
              END IF
!
!             == SPIN ORBIT MATRIX ELEMENTS ====================================
              IF(TSO) THEN
                SVAR=0.5D0*(0.D0+0.D0)
                CALL SCHROEDINGER$DREL(GID,NR,AEPOT(:,1,1),SVAR,DREL)
                CALL RADIAL$DERIVE(GID,NR,DREL,AUX)
                AUX(:)=AUX(:)*R(:)*AEPHI(:,LN1)*AEPHI(:,LN2)  ! (AUX/R)*R**2
                CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
                DO IS1=1,2
                  DO IS2=1,2
                    I1=I1A-1+IS1
                    I2=I2A-1+IS2
                    HVV(I1,I2)=HVV(I1,I2)+CLS(LM1,IS1,LM2,IS2)*SVAR
                  ENDDO
                ENDDO
              END IF
! 
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(TPR) THEN
        WRITE(*,FMT='(80("="),T20,A)')' VALENCE OVERLAP '
        WRITE(*,FMT='(I5,"     ",100I10)')0,(LMN2,LMN2=1,2*LMNX)
        DO LMN1=1,2*LMNX
          WRITE(*,FMT='(I5,"OVV: ",100F10.5)')LMN1,REAL(OVV(LMN1,:))
        ENDDO
        WRITE(*,FMT='(80("="),T20,A)')' VALENCE HAMILTONIAN '
        WRITE(*,FMT='(I5,"     ",100I10)')0,(LMN2,LMN2=1,2*LMNX)
        DO LMN1=1,2*LMNX
          WRITE(*,FMT='(I5,"HVV: ",100F10.5)')LMN1,REAL(HVV(LMN1,:))
        ENDDO
      END IF

!!$WRITE(*,*)'CORE DENSITY AT NUCLEUS    ', AEPSI(1,1)**2
!!$WRITE(*,*)'VALENCE DENSITY AT NUCLEUS ', AEPHI(1,1)**2
!!$CALL SETUP_WRITEPHI('AEPOT.DAT',GID,NR,LMRX*NDIMD,AEPOT)
!!$CALL SETUP_WRITEPHI('AEPOT1.DAT',GID,NR,LMRX*NDIMD,AEPOT1)
!!$CALL SETUP_WRITEPHI('ATPOT.DAT',GID,NR,1,ATPOT)
!!$CALL ERROR$STOP('FORCED')
      DEALLOCATE(AEPOT1)
!
!     ==========================================================================
!     ==  DIAGONALIZATION OF THE CORE HAMILTONIAN                             ==
!     ==========================================================================
      ALLOCATE(EIGENVAL(2*LMNCX))
      ALLOCATE(EIGENVEC(2*LMNCX,2*LMNCX))
PRINT*,'MARKE 3A',2*LMNCX
PRINT*,'MARKE 3B',SHAPE(HCC)
PRINT*,'MARKE 3C',SHAPE(EIGENVAL)
PRINT*,'MARKE 3D',SHAPE(EIGENVEC)
      CALL LIB$DIAGC8(2*LMNCX,HCC,EIGENVAL,EIGENVEC)
PRINT*,'MARKE 4'

      IF(TPR) THEN  
        CALL CONSTANTS$GET('EV',EV)
        WRITE(*,FMT='("EIGENVALUES/ev ",5F12.3)')EIGENVAL/ev
      END IF
      DEALLOCATE(EIGENVEC)
      DEALLOCATE(EIGENVAL)

!
!     ==========================================================================
!     ==  DIAGONALIZATION OF THE VALENCE HAMILTONIAN                          ==
!     ==========================================================================
      IF(TPR) THEN  
        ALLOCATE(EIGENVAL(2*LMNX))
        ALLOCATE(EIGENVEC(2*LMNX,2*LMNX))
        CALL LIB$GENERALEIGENVALUEC8(2*LMNX,HVV,OVV,EIGENVAL,EIGENVEC)
        WRITE(*,FMT='("===== VALENCE EIGENVALUES =====")')
        WRITE(*,FMT='("EIGENVALUES ",5F12.3)')EIGENVAL
        DEALLOCATE(EIGENVEC)
        DEALLOCATE(EIGENVAL)
      END IF
!
!     ==========================================================================
!     ==  DIAGONALIZATION OF THE CORE-AND-VALENCE HAMILTONIAN                 ==
!     ==========================================================================
! CAUTION IT DOES NOT CONTAIN THE KINETIC AND ATOMIC POTENTIAL ENERGY!
      LMNTX=LMNCX+LMNX
      ALLOCATE(EIGENVAL(2*LMNTX))
      ALLOCATE(EIGENVEC(2*LMNTX,2*LMNTX))
      ALLOCATE(HT(2*LMNTX,2*LMNTX))
      ALLOCATE(OT(2*LMNTX,2*LMNTX))
      HT(:2*LMNCX,:2*LMNCX)=HCC
      HT(2*LMNCX+1:,2*LMNCX+1:)=HVV
      HT(:2*LMNCX,2*LMNCX+1:)=HCV
      HT(2*LMNCX+1:,:2*LMNCX)=TRANSPOSE(CONJG(HCV))
      OT(:2*LMNCX,:2*LMNCX)=OCC
      OT(2*LMNCX+1:,2*LMNCX+1:)=OVV
      OT(:2*LMNCX,2*LMNCX+1:)=OCV
      OT(2*LMNCX+1:,:2*LMNCX)=TRANSPOSE(CONJG(OCV))
      CALL LIB$GENERALEIGENVALUEC8(2*LMNTX,HT,OT,EIGENVAL,EIGENVEC)
!      CALL LIB$DIAGC8(2*LMNTX,HT,EIGENVAL,EIGENVEC)

      IF(TPR) THEN  
        WRITE(*,FMT='("===== CORE AND VALENCE EIGENVALUES =====")')
        WRITE(*,FMT='("EIGENVALUES ",5F12.3)')EIGENVAL
      END IF
      DEALLOCATE(EIGENVEC)
      DEALLOCATE(EIGENVAL)
      DEALLOCATE(HT)
!
!     ==========================================================================
!     ==  convert one-center density matrix into up-dn-representation         ==
!     ==========================================================================
      ALLOCATE(RHO1VV(2*LMNX,2*LMNX))
      RHO1VV=(0.D0,0.D0)
      DO LMN1=1,LMNX
        DO LMN2=1,LMNX
          I1A=2*(LMN1-1)+1
          I1B=2*(LMN1-1)+2
          I2A=2*(LMN2-1)+1
          I2B=2*(LMN2-1)+2
          RHO1VV(I1A,I2A)=0.5d0*DENMAT(LMN1,LMN2,1)
          RHO1VV(I1B,I2B)=0.5d0*DENMAT(LMN1,LMN2,1)
          IF(NDIMD.EQ.2) THEN
            RHO1VV(I1A,I2A)=RHO1VV(I1A,I2A)+0.5d0*DENMAT(LMN1,LMN2,2)
            RHO1VV(I1B,I2B)=RHO1VV(I1B,I2B)-0.5d0*DENMAT(LMN1,LMN2,2)
          ELSE IF (NDIMD.EQ.4) THEN
            RHO1VV(I1A,I2B)=RHO1VV(I1A,I2B)+0.5d0*DENMAT(LMN1,LMN2,2)
            RHO1VV(I1B,I2A)=RHO1VV(I1B,I2A)+0.5d0*DENMAT(LMN1,LMN2,2)
            RHO1VV(I1A,I2B)=RHO1VV(I1A,I2B)-0.5d0*CI*DENMAT(LMN1,LMN2,3)
            RHO1VV(I1B,I2A)=RHO1VV(I1B,I2A)+0.5d0*CI*DENMAT(LMN1,LMN2,3)
            RHO1VV(I1A,I2A)=RHO1VV(I1A,I2A)+0.5d0*DENMAT(LMN1,LMN2,4)
            RHO1VV(I1B,I2B)=RHO1VV(I1B,I2B)-0.5d0*DENMAT(LMN1,LMN2,4)
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ENERGY CORRECTION BEYOND FROZEN CORE
!     ==========================================================================
!
!     ===== CALCULATE SELF ENERGY ==============================================
      ALLOCATE(GREEN(2*LMNCX,2*LMNCX))
      CALL LIB$INVERTC8(2*LMNCX,EPSILONNU*OCC-HCC,GREEN)

      ALLOCATE(cmat2(2*LMNCX,2*LMNCX))
      CALL CORE_MAXCORECORR(GID,NR,NC,LB,AEPSI,ATPOT,NDIMD,LMRX,AEPOT &
     &                     ,lmncx,cmat2)
!!$print*,'green',green
!!$print*,'cmat2',cmat2
!!$print*,'cmat2*green',cmat2*green
      e1=-sum(real(green*cmat2))
      WRITE(*,FMT=FMTSTRING)E1,'CORE - NONCORE CONTRIBUTION (new)'



      SELFENERGY1=MATMUL(TRANSPOSE(HCV-EPSILONNU*OCV) &
     &                 ,MATMUL(GREEN,HCV-EPSILONNU*OCV))
      DEALLOCATE(GREEN)
      ERHO=SUM(REAL(SELFENERGY1*RHO1VV))
!!$print*,'selfenergy1 ',selfenergy1
!!$print*,'rho1vv      ',rho1vv
print*,'erho        ',erho



      ALLOCATE(GREEN(2*LMNX,2*LMNX))
      CALL LIB$INVERTC8(2*LMNX,OVV,GREEN)
      E1=SUM(REAL(SELFENERGY1*GREEN))

      WRITE(*,FMT='("FROZEN CORE CORRECTION :")')
      WRITE(*,FMT=FMTSTRING)E1,'CORE - NONCORE CONTRIBUTION'
      WRITE(*,FMT=FMTSTRING)ERHO,'CORE - FILLED NONCORE CONTRIBUTION'
      WRITE(*,FMT=FMTSTRING)E1-ERHO,'CORE - EMPTY NONCORE CONTRIBUTION'


      DEALLOCATE(GREEN)
!
!     ==========================================================================
!     ==  SELF ENERGY AT ENERGY EPSILONNU                                     ==
!     ==========================================================================

      IF(TPR) THEN
        WRITE(*,FMT='(80("="),T20,A)')' SELF ENERGY '
        WRITE(*,FMT='(I5,"     ",100I10)')0,(LMN2,LMN2=1,2*LMNX)
        DO LMN1=1,2*LMNX
          WRITE(*,FMT='(I5,"SIGMA: ",100F10.5)')LMN1,REAL(SELFENERGY1(LMN1,:))
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CONVERT SELF ENERGY INTO (TOT,X,Y,Z) 
!     ==========================================================================
      SELFENERGY(:,:,:)=(0.D0,0.D0)
      DO LMN1=1,LMNX
        DO LMN2=1,LMNX
          I1=2*(LMN1-1)
          I2=2*(LMN2-1)
          SELFENERGY(LMN1,LMN2,1)=REAL(SELFENERGY1(I1+1,I2+1) &
     &                                +SELFENERGY1(I1+2,I2+2),KIND=8)
          IF(NDIMD.EQ.2) THEN
            SELFENERGY(LMN1,LMN2,2)=REAL(SELFENERGY1(I1+1,I2+1) &
     &                            -SELFENERGY1(I1+2,I2+2),KIND=8)
          ELSE IF(NDIMD.EQ.4) THEN
            SELFENERGY(LMN1,LMN2,2)=REAL(SELFENERGY1(I1+1,I2+2) &
     &                            +SELFENERGY1(I1+2,I2+1),KIND=8)
            SELFENERGY(LMN1,LMN2,3)=AIMAG(SELFENERGY1(I1+1,I2+2) &
     &                             -SELFENERGY1(I1+2,I2+1))
            SELFENERGY(LMN1,LMN2,4)=REAL(SELFENERGY1(I1+1,I2+1) &
     &                            -SELFENERGY1(I1+2,I2+2),KIND=8)
          END IF
        ENDDO
      ENDDO



      IF(TPR) THEN
        WRITE(*,FMT='(80("="),T20,A)')' SELF ENERGY IN (T,X,Y,Z) REPRESENTATION'
        WRITE(*,FMT='(I5,"     ",100I10)')0,(LMN2,LMN2=1,2*LMNX)
        DO IDIMD=1,NDIMD
          IF(IDIMD.EQ.1) THEN 
            WRITE(*,FMT='(80("="),T20,A)')' TOTAL SELF ENERGY'
          ELSE IF(NDIMD.EQ.2.AND.IDIMD.EQ.2) THEN
            WRITE(*,FMT='(80("="),T20,A)')' SPIN SELF ENERGY'
          ELSE IF(NDIMD.EQ.4) THEN
            IF(IDIMD.EQ.2) THEN
              WRITE(*,FMT='(80("="),T20,A)')' SELF ENERGY X-SPIN'
            ELSE IF(IDIMD.EQ.3) THEN
              WRITE(*,FMT='(80("="),T20,A)')' SELF ENERGY Y-SPIN'
            ELSE IF(IDIMD.EQ.4) THEN
              WRITE(*,FMT='(80("="),T20,A)')' SELF ENERGY Z-SPIN'
            END IF
          END IF
          DO LMN1=1,LMNX
            WRITE(*,FMT='(I5,"SIGMA: ",100F12.5)')LMN1,SELFENERGY(LMN1,:,IDIMD)
          ENDDO
        ENDDO
      END IF




IF(TNOSELFENERGY) SELFENERGY(:,:,:)=(0.D0,0.D0)

!
!     ==========================================================================
!     ==  ENERGY GAIN FROM CORE-VALENCE INTERACTION                           ==
!     ==========================================================================
!!$      ALLOCATE(OVVIN(LMNX,LMNX))      
!!$      CALL LIB$INVERTR8(LMNX,OVV,OVVIN)
!!$         MATMUL(MATMUL(HCV,OVVIN),TRANSPOSE(HCV))

!
!     ==========================================================================
!     ==  CLOSE DOWN                                                          ==
!     ==========================================================================
      DEALLOCATE(LOX)
      DEALLOCATE(AEPHI)
      DEALLOCATE(LB)
      DEALLOCATE(EB)
      DEALLOCATE(AEPSI)
      DEALLOCATE(HCC)
      DEALLOCATE(HCV)

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_maxcorecorr(gid,nr,nc,lb,aepsi,atpot,ndimd,lmrx,aepot &
     &                           ,lmncx,cmat2)
!     **************************************************************************
!     ** calculates matrix elements between the frozen core states            **
!     **    <psi_1|(v-v_at) (1-P_c) (v-v_at)|psi_2>                           **
!     **   = <psi_1|(v-v_at)^2|psi_2>                                         **
!     **    -sum_3 <psi_1|v-v_at|psi_3><psi_3|v-v_at|psi_2>                   **
!     **  P_c is the projection onto the core states                          **
!     **                                                                      **
!     ** the potentials are provided in a (total,x,y,z) representation        **
!     ** the components are the prefactors of Pauli-matrices.                 **
!     **************************************************************************
      implicit none
      integer(4),intent(in) :: gid            ! grid id
      integer(4),intent(in) :: nr             ! #(radial grid points)
      integer(4),intent(in) :: nc             ! #(radial core states)
      integer(4),intent(in) :: lb(nc)         ! azimuthal angular momentum
      real(8)   ,intent(in) :: aepsi(nr,nc)   ! core states
      real(8)   ,intent(in) :: atpot(nr)      ! atom potential
      integer(4),intent(in) :: ndimd          ! #(density components)
      integer(4),intent(in) :: lmrx           ! x#(density angular momenta)
      real(8)   ,intent(in) :: aepot(nr,lmrx,ndimd) ! potential
      integer(4),intent(in) :: lmncx          ! #(core states w/o spin)
      complex(8),intent(out):: cmat2(2*lmncx,2*lmncx)
      real(8)   ,parameter  :: dsmall=1.d-8
      real(8)   ,parameter  :: pi=4.d0*atan(1.d0)
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      complex(8),parameter  :: ci=(0.d0,1.d0)  !sqrt(-1)
      REAL(8)   ,ALLOCATABLE:: DPOT(:,:,:) !(NR,LMR2X,NDIMD)
      REAL(8)   ,ALLOCATABLE:: DPOT2(:,:,:) !(NR,LMR2X,NDIMD)
      INTEGER(4)            :: LRX
      integer(4)            :: lmr2x
      integer(4)            :: lm1,lm2,lm
      integer(4)            :: lmn1,lmn2
      integer(4)            :: ln1,ln2
      integer(4)            :: l1,l2
      integer(4)            :: idimd
      integer(4)            :: i1a,i1b,i2a,i2b
      real(8)               :: cg ! gaunt coefficient
      real(8)               :: aux(nr)
      real(8)               :: val
      real(8)               :: x,x2
      real(8)               :: r(nr)
      complex(8),allocatable:: cmat(:,:) !(2*lmnc,2*lmnc)
      real(8)   ,allocatable:: mat(:,:,:) !(nc,nc,ndimd)
      real(8)   ,allocatable:: mat2(:,:,:) !(nc,nc,ndimd)
!     **************************************************************************
      call radial$r(gid,nr,r)
!
!     ==========================================================================
!     == difference potential and squared difference potential
!     ==========================================================================
      allocate(dpot(nr,lmrx,ndimd))
      dpot(:,:,:)=aepot(:,:,:)
      dpot(:,1,1)=dpot(:,1,1)-atpot(:)
!
!     ==========================================================================
!     == squared difference potential
!     ==========================================================================
      lrx=int(sqrt(real(lmrx-1)+dsmall))
      lmr2x=(2*lrx+1)**2
      allocate(dpot2(nr,lmr2x,ndimd))
      dpot2=0.d0
      do lm1=1,lmrx
        do lm2=1,lmrx
          do lm=1,lmr2x
            CALL CLEBSCH(LM1,LM2,LM,CG) !gaunt coefficient via table lookup
            if(abs(cg).lt.dsmall) cycle
            do idimd=1,ndimd
              dpot2(:,lm,1)=dpot2(:,lm,1)+cg*dpot(:,lm1,idimd)*dpot(:,lm2,idimd)
            enddo
            do idimd=2,ndimd
              dpot2(:,lm,idimd)=dpot2(:,lm,idimd) &
                               +cg*(dpot(:,lm1,1)*dpot(:,lm2,idimd) &
     &                             +dpot(:,lm1,idimd)*dpot(:,lm2,1))
            enddo
          enddo
        enddo
      enddo
!
!     ==========================================================================
!     == matrix elements of the difference potential and the squared diffpot
!     ==========================================================================
      IF(LMNCX.NE.SUM(2*LB+1)) THEN
        CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED')
        CALL ERROR$STOP('CORE_MAXCORECORR')
      END IF
      allocate(mat(lmncx,lmncx,ndimd))
      allocate(mat2(lmncx,lmncx,ndimd))
      mat(:,:,:)=0.d0
      mat2(:,:,:)=0.d0
      do lm=1,lmr2x
        do idimd=1,ndimd
          do ln1=1,nc
            do ln2=1,nc
              if(lm.le.lmrx) then
                aux=r**2*aepsi(:,ln1)*aepsi(:,ln2)*dpot(:,lm,idimd)
                call radial$integral(gid,nr,aux,x)
              else
                x=0.d0
              end if
              aux=r**2*aepsi(:,ln1)*aepsi(:,ln2)*dpot2(:,lm,idimd)
              call radial$integral(gid,nr,aux,x2)
              l1=lb(ln1)
              l2=lb(ln2)
              lmn1=sum(2*lb(:ln1-1)+1)
              do lm1=l1**2+1,(l1+1)**2
                lmn1=lmn1+1
                lmn2=sum(2*lb(:ln2-1)+1)
                do lm2=l2**2+1,(l2+1)**2
                  lmn2=lmn2+1
                  CALL CLEBSCH(LM1,LM2,LM,CG) !gaunt coeff. via table lookup
                  mat(lmn1,lmn2,idimd) =mat(lmn1,lmn2,idimd) +cg*x
                  mat2(lmn1,lmn2,idimd)=mat2(lmn1,lmn2,idimd)+cg*x2
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
!     ==========================================================================
!     == 
!     ==========================================================================
      allocate(cmat(2*lmncx,2*lmncx))
      cmat=(0.d0,0.d0)
      cmat2=(0.d0,0.d0)
      do lmn1=1,lmncx
        i1a=2*(lmn1-1)+1
        i1b=i1a+1
        do lmn2=1,lmncx
          i2a=2*(lmn2-1)+1
          i2b=i2a+1
          cmat(i1a,i2a) =cmat(i1a,i2a) +mat(lmn1,lmn2,1)
          cmat(i1b,i2b) =cmat(i1b,i2b) +mat(lmn1,lmn2,1)
          cmat2(i1a,i2a)=cmat2(i1a,i2a)+mat2(lmn1,lmn2,1)
          cmat2(i1b,i2b)=cmat2(i1b,i2b)+mat2(lmn1,lmn2,1)
          if(ndimd.gt.1) then
            cmat(i1a,i2a) =cmat(i1a,i2a) +mat(lmn1,lmn2,ndimd)
            cmat(i1b,i2b) =cmat(i1b,i2b) -mat(lmn1,lmn2,ndimd)
            cmat2(i1a,i2a)=cmat2(i1a,i2a)+mat2(lmn1,lmn2,ndimd)
            cmat2(i1b,i2b)=cmat2(i1b,i2b)-mat2(lmn1,lmn2,ndimd)
            if(ndimd.eq.4) then
              cmat(i1a,i2b) =cmat(i1a,i2b)+mat(lmn1,lmn2,2)-ci*mat(lmn1,lmn2,3)
              cmat(i1b,i2a) =cmat(i1b,i2a)+mat(lmn1,lmn2,2)+ci*mat(lmn1,lmn2,3)
              cmat2(i1a,i2b)=cmat2(i1a,i2b) &
     &                      +mat2(lmn1,lmn2,2)-ci*mat2(lmn1,lmn2,3)
              cmat2(i1b,i2a)=cmat2(i1b,i2a) &
     &                      +mat2(lmn1,lmn2,2)+ci*mat2(lmn1,lmn2,3)
            end if
          end if
        enddo
      enddo
!
!     ==========================================================================
!     == subtract the projection onto the core states
!     ==========================================================================
      cmat2=cmat2-matmul(cmat,cmat)
!!
!     ==========================================================================
!     == write result                                                         ==
!     ==========================================================================
      if(tpr) then
        do i1a=1,2*lmncx
          write(*,fmt='("mat---- ",20f10.6)')real(cmat(i1a,:))
        enddo
        do i1a=1,2*lmncx
          write(*,fmt='("mat2---- ",20f10.6)')real(cmat2(i1a,:))
        enddo
      end if
      return
      end

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_RHOMAG2UPDN(NDIMD,X,H)
!     **************************************************************************
!     ** TRANSFORMS MATRIX ELEMENTS INTO THE UPDN REPRESENTATION              **
!     **   NDIMD=1 H=SIGMA_0*X(1)                                             **
!     **   NDIMD=2 H=SIGMA_0*X(1)+SIGMAZ*X(2)                                 **
!     **   NDIMD=4 H=SIGMA_0*X(1)+SIGMAX*X(2)+SIGMAY*X(3)+SIGMAZ*X(4)         **
!     ** WHERE THE SIGMA ARE PAULI MATRICES                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: X(NDIMD)
      COMPLEX(8),INTENT(OUT):: H(2,2)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8),PARAMETER  :: CZERO=(0.D0,0.D0)
!     **************************************************************************
      H(1,1)=X(1)
      H(2,1)=CZERO
      H(1,2)=CZERO
      H(2,2)=X(1)
      IF(NDIMD.EQ.2) THEN
        H(1,1)=H(1,1)+X(2)
        H(2,2)=H(2,2)-X(2)
      ELSE IF(NDIMD.EQ.4) THEN
        H(1,1)=H(1,1)+X(4)
        H(2,2)=H(2,2)-X(4)
        H(1,2)=H(1,2)+X(2)-CI*X(3)
        H(2,1)=H(2,1)+X(2)+CI*X(3)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_HYPERFINE(IAT,GID,NR,NDIMD,LMRX,AERHO)
!     **************************************************************************
!     ** DRIVER FOR THE CALCULATION OF HYPERFINE PARAMETERS                   **
!     ** (ISOMERSHIFT, ELECTRIC FIELD GRADIENTS, FERMI CONTACT, ANISOTROPIC)  **
!     **************************************************************************
      use strings_module
      USE CORE_MODULE, ONLY : FILEHFID      
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT  ! ATOM INDEX
      INTEGER(4),INTENT(IN) :: GID  ! GRID ID
      INTEGER(4),INTENT(IN) :: NR   ! #(GRID POINTS)
      INTEGER(4),INTENT(IN) :: NDIMD 
      INTEGER(4),INTENT(IN) :: LMRX ! #(ANGULAR MOMENTA) FOR 1C-POTENTIAL)
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX,NDIMD)  ! 1C-AE-DENSITY
      CHARACTER(32)        :: ID
      logical(4)           :: tini=.false.
!     **************************************************************************
PRINT*,'================================================================'
PRINT*,'==  CAUTION!  THIS ROUTINE IS UNDER CONSTRUCTION!!!   =========='
PRINT*,'================================================================'
      ID=filehfid
      if(.not.tini) then
        tini=.true.
        CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.HFDAT')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
      end if
!
!     ==========================================================================
!     ==  ANALYSIS (ELECTRIC FIELD GRADIENTS, ETC.)                           ==
!     ==========================================================================
print*,'marke 1'
      CALL CORE_ISOMERSHIFT(IAT,GID,NR,NDIMD,LMRX,AERHO)
print*,'marke 2'
      CALL CORE_EFG(IAT,GID,NR,NDIMD,LMRX,AERHO)
print*,'marke 3'
      CALL CORE_FERMICONTACT(IAT,GID,NR,NDIMD,LMRX,AERHO)
print*,'marke 4'
      CALL CORE_ANISOTROPIC(IAT,GID,NR,NDIMD,LMRX,AERHO)
print*,'marke 5'
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_EFG(IAT,GID,NR,NDIMD,LMRX,AERHO)
!     **************************************************************************
!     ** CALCULATE ELECTRIC FIELD GRADIENTS FROM THE 1C-ALL ELECTRON POTENTIAL**
!     **                                                                      **
!     ** CAUTION! THE ELECTRIC FIELD GRADIENT IS DETERMINED FROM THE          **
!     ** ALL-ELECTRON 1C DENSITY FROM WITHIN THE COVALENT RADIUS ONLY.        **
!     ** THE POTENTIAL FROM THE DENSITY BEYOND THE COVALENT RADIUS IS MISSING!**
!     **
!     ** DO NOT USE AEPOT DIRECTLY BECAUSE (1) IT INCLUDES ERRONESOULY        **
!     ** THE XC POTENTIAL  AND (2) INCLUDES WEIRD NUMERICAL BEHAVIOR AT THE   **
!     **  ORIGIN WHICH I DO NOT YET UNDERSTAND. PB                            **
!     **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE CORE_MODULE, ONLY: EXTRAPOLATE &
     &                      ,DTOXYZ &
     &                      ,FILEHFID
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT  ! ATOM INDEX
      INTEGER(4),INTENT(IN) :: GID  ! GRID ID
      INTEGER(4),INTENT(IN) :: NR   ! #(GRID POINTS)
      INTEGER(4),INTENT(IN) :: NDIMD 
      INTEGER(4),INTENT(IN) :: LMRX ! #(ANGULAR MOMENTA) FOR 1C-POTENTIAL)
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX,NDIMD)  ! 1C-AE-DENSITY
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4),PARAMETER  :: NP=5
      INTEGER(4)            :: NFIL
      REAL(8)               :: V2(5)     ! 
      REAL(8)               :: R(NR)     ! RADIAL GRID
      REAL(8)               :: EFG(3,3)  ! EFT MATRIX D2 POT/(DX_I DX_J)
      REAL(8)               :: EIG(3)    ! EIGENVALUES OF EFG MATRIX
      REAL(8)               :: U(3,3)    ! EIGENVECTORS OF EG MATRIX
      CHARACTER(32)         :: ATOMNAME
      INTEGER(4)            :: I,LM,ISP
      REAL(8)               :: SVAR,SVAR1
      REAL(8)               :: AEPOT(NR,5) !TRASH
      REAL(8)               :: WORK1(NR),WORK2(NR)
      REAL(8)               :: AEZ  ! ATOMIC NUMBER
      REAL(8)               :: RCOV ! COVALENT RADIUS
      REAL(8)               :: RNUC ! RADIUS OF THE NUCLEUS
!     **************************************************************************
      CALL FILEHANDLER$UNIT(FILEHFID,NFIL)
      CALL RADIAL$R(GID,NR,R)
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
      CALL ATOMLIST$GETR8('Z',IAT,AEZ)
      CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RCOV)
      CALL PERIODICTABLE$GET(NINT(AEZ),'RNUC',RNUC)
      
      IF(LMRX.GE.9) THEN
        DO LM=5,9
          CALL RADIAL$POISSON(GID,NR,2,AERHO(:,LM,1),WORK1)
!         ==  REMOVE THE POTENTIAL OF CHARGEDENSITY BEYOND COVALENT RADIUS =====
          CALL RADIAL$VALUE(GID,NR,WORK1,RCOV,SVAR)
          WORK1(:)=WORK1(:)-SVAR*(R(:)/RCOV)**2

!         == EXTRACT EFG FROM THE ORIGIN =======================================
          CALL EXTRAPOLATE(NP,R(1:NP),WORK1(1:NP)/R(1:NP)**2,0.D0,V2(LM-4))
!AEPOT(:,LM-4)=WORK1
        ENDDO
!!$PRINT*,'========= AERHO FOR L=2 ==============='
!!$DO I=1,NR
!!$!PRINT*,R(I),AERHO(I,5:9,1)
!!$PRINT*,R(I)-RCOV,AEPOT(I,:)
!!$ENDDO
!!$PRINT*,'========= AERHO DONE ==============='
      ELSE
        V2(:)=0.D0
      END IF
!
!     __TRANSFORM LM REPRESENTATION INTO CARTESIAN COORDINATES______________
!     __VLM=V(R)/R**2___VXYZ=D2V/(DI*DJ)____________________________________
      CALL DTOXYZ(V2,EFG)
!
!     ==========================================================================
!     == WRITE RESULT 
!     ==========================================================================
      WRITE(NFIL,FMT='(/80("-"))')
      SVAR=1.D-21
      CALL CONSTANTS('VOLT',SVAR1) ; SVAR=SVAR/SVAR1 
      CALL CONSTANTS('METER',SVAR1) ; SVAR=SVAR*SVAR1**2 
      SVAR=-SVAR   ! THE ELEMENTARY CHARGE IS -1 ELECTRON CHARGE
      PRINT*,SVAR,'SHOULD BE EQUAL TO ',-9.7175D0
!
!     == TOTAL ELECTIC FIELD GRADIENT ====================================
      CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
      CALL REPORT$STRING(NFIL,'ELECTRIC FIELD GRADIENT FOR ATOM ' &
     &                      //TRIM(ATOMNAME)//' IN [10**21*V/M**2]')
      CALL LIB$DIAGR8(3,EFG,EIG,U)

      WRITE(NFIL,FMT='(T20,"VALUE",T40,"DIRECTION")')
      DO I=1,3
        WRITE(NFIL,FMT='(T5,F20.5,T30,3F10.5)')SVAR*EIG(I),U(:,I)
      ENDDO
  
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_ISOMERSHIFT(IAT,GID,NR,NDIMD,LMRX,AERHO)
!     **************************************************************************
!     ** CALCULATE ELECTRIC FIELD GRADIENTS FROM THE 1C-ALL ELECTRON POTENTIAL**
!     **************************************************************************
      USE CORE_MODULE, ONLY: EXTRAPOLATE &
     &                      ,FILEHFID
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT  ! ATOM INDEX
      INTEGER(4),INTENT(IN) :: GID  ! GRID ID
      INTEGER(4),INTENT(IN) :: NR   ! #(GRID POINTS)
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: LMRX ! #(ANGULAR MOMENTA) FOR 1C-POTENTIAL)
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX)  ! 1C-AE-DENSITY
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4),PARAMETER  :: NP=5
      INTEGER(4)            :: NFIL
      REAL(8)               :: R(NR)     ! RADIAL GRID
      REAL(8)               :: aecore(NR)     ! core density
      REAL(8)               :: RHO0      ! ELECTRON DENSITY AT THE ORIGIN
      REAL(8)               :: svar
      CHARACTER(32)         :: ATOMNAME
      integer(4)            :: isp !atom type 
!     **************************************************************************
      CALL FILEHANDLER$UNIT(FILEHFID,NFIL)
!
!     ==  collect core density =================================================
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
      CALL SETUP$unSELECT()

!     ==  core density at the nuclear site =====================================
      CALL RADIAL$R(GID,NR,R)
      CALL EXTRAPOLATE(NP,R(1:NP),AERHO(1:NP,1),0.D0,RHO0)
      CALL EXTRAPOLATE(NP,R(1:NP),AEcore(1:NP),0.D0,svar)
      rho0=rho0+svar
      RHO0=RHO0*Y0
      CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
      WRITE(NFIL,FMT='(/80("-"))')
      CALL REPORT$R8VAL(NFIL &
     &                      ,'ELECTRON DENSITY AT THE NUCLEUS OF ATOM ' &
     &                       //TRIM(ATOMNAME),RHO0,'1/ABOHR^3')
      CALL REPORT$R8VAL(NFIL,'ISOMERSHIFT FOR ATOM '//TRIM(ATOMNAME) &
     &                            ,RHO0,'E/A_0^3')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_FERMICONTACT(IAT,GID,NR,NDIMD,LMRX,AERHO)
!     **************************************************************************
!     ** THE FERMI CONTACT HYPERFINE INTERACTION TERM
!     ** SEE BUCHER00_EURJPHYS21_19 AND SOLIVEREZ80_JPCSS13_L1017             **
!     **************************************************************************
      USE CORE_MODULE, ONLY: EXTRAPOLATE &
     &                      ,FILEHFID
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT  ! ATOM INDEX
      INTEGER(4),INTENT(IN) :: GID  ! GRID ID
      INTEGER(4),INTENT(IN) :: NR   ! #(GRID POINTS)
      INTEGER(4),INTENT(IN) :: LMRX ! #(ANGULAR MOMENTA) FOR 1C-POTENTIAL)
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX,NDIMD)  ! 1C-AE-DENSITY
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4),PARAMETER  :: NP=5
      INTEGER(4)            :: NFIL
      REAL(8)               :: R(NR)     ! RADIAL GRID
      REAL(8)               :: WORK1(NR) ! WORK ARRAY
      REAL(8)               :: RHOS0     ! ELECTRON DENSITY AT THE ORIGIN
      REAL(8)               :: SVAR,SVAR1
      CHARACTER(32)         :: ATOMNAME
      INTEGER(4)            :: IDIMD
!     **************************************************************************
      IF(NDIMD.EQ.1) RETURN
      CALL FILEHANDLER$UNIT(FILEHFID,NFIL)
      CALL RADIAL$R(GID,NR,R)
!     == DETERMINE SPIN DENSITY AT THE NUCLEAR SITE ============================
      WORK1=0.D0
      DO IDIMD=2,NDIMD
        WORK1(:)=WORK1(:)+AERHO(:,1,IDIMD)**2
      ENDDO
      WORK1=SQRT(WORK1)
      CALL EXTRAPOLATE(NP,R(1:NP),WORK1(1:NP),0.D0,RHOS0)
      RHOS0=RHOS0*Y0

      SVAR=2.D0/3.D0
      CALL CONSTANTS('MU0',SVAR1)          ; SVAR=SVAR*SVAR1 
!     == MULTIPLY WITH MAGNETIC MOMENT OF THE ELECTRON
      CALL CONSTANTS('BOHRMAGNETON',SVAR1) ; SVAR=SVAR*SVAR1
      CALL CONSTANTS('GE',SVAR1)           ; SVAR=SVAR*SVAR1
      CALL CONSTANTS('HBAR',SVAR1)         ; SVAR=SVAR*0.5D0*SVAR1
!     ==   
      CALL CONSTANTS('TESLA',SVAR1)        ; SVAR=SVAR/SVAR1
      WRITE(NFIL,FMT='(/80("-"))')
      PRINT*,'THIS NUMBER SHOULD BE 104.98 : ',SVAR
      CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
      CALL REPORT$R8VAL(NFIL,'SPIN DENSITY AT THE NUCLEUS OF ATOM ' &
     &                            //TRIM(ATOMNAME),RHOS0,'1/ABOHR^3')
      CALL REPORT$R8VAL(NFIL,'FERMI CONTACT HYPERFINE FIELD FOR ATOM ' &
     &                            //TRIM(ATOMNAME),SVAR*RHOS0,'TESLA')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE_ANISOTROPIC(IAT,GID,NR,NDIMD,LMRX,AERHO)
!     **************************************************************************
!     ** THE ANISOTROPIC  HYPERFINE INTERACTIONS                              **

!     **************************************************************************
      USE CORE_MODULE, ONLY: EXTRAPOLATE &
     &                      ,FILEHFID
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT  ! ATOM INDEX
      INTEGER(4),INTENT(IN) :: GID  ! GRID ID
      INTEGER(4),INTENT(IN) :: NR   ! #(GRID POINTS)
      INTEGER(4),INTENT(IN) :: LMRX ! #(ANGULAR MOMENTA) FOR 1C-POTENTIAL)
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: AERHO(NR,LMRX,NDIMD)  ! 1C-AE-DENSITY
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4),PARAMETER  :: NP=5
      INTEGER(4)            :: NFIL
      REAL(8)               :: R(NR)     ! RADIAL GRID
      REAL(8)               :: WORK1(NR),WORK2(NR)    
      REAL(8)               :: RHOS0      ! ELECTRON DENSITY AT THE ORIGIN
      REAL(8)               :: V2(5)
      REAL(8)               :: ANIS(3,3)
      REAL(8)               :: EIG(3)
      REAL(8)               :: U(3,3)
      REAL(8)               :: SVAR,SVAR1
      CHARACTER(32)         :: ATOMNAME
      INTEGER(4)            :: LM,I,IDIMD
!     **************************************************************************
      IF(NDIMD.EQ.1) RETURN
      CALL FILEHANDLER$UNIT(FILEHFID,NFIL)
      CALL ATOMLIST$GETCH('NAME',IAT,ATOMNAME)
      CALL RADIAL$R(GID,NR,R)

      IF(LMRX.GE.9) THEN
        DO LM=5,9
          WORK1=0.D0
          DO IDIMD=2,NDIMD
            WORK1(:)=WORK1(:)+AERHO(:,LM,IDIMD)**2
          ENDDO
          WORK1=SQRT(WORK1)
          CALL RADIAL$POISSON(GID,NR,2,WORK1,WORK2)
          CALL EXTRAPOLATE(NP,R(1:NP),WORK2(1:NP)/R(1:NP)**2,0.D0,V2(LM-4))
        ENDDO
      ELSE
        V2(:)=0.D0
      END IF

!     __TRANSFORM LM REPRESENTATION INTO CARTESIAN COORDINATES__
      CALL DTOXYZ(V2,ANIS)

      SVAR=1.D0/(4.D0*PI)
      CALL CONSTANTS('MU0',SVAR1)          ; SVAR=SVAR*SVAR1 
!     ====================================================================
      CALL CONSTANTS('BOHRMAGNETON',SVAR1) ; SVAR=SVAR*SVAR1
      CALL CONSTANTS('GE',SVAR1)           ; SVAR=SVAR*SVAR1
      CALL CONSTANTS('HBAR',SVAR1)         ; SVAR=SVAR*0.5D0*SVAR1
!     ====================================================================
      CALL CONSTANTS('TESLA',SVAR1)        ; SVAR=SVAR/SVAR1
      WRITE(NFIL,FMT='(/80("-"))')
      PRINT*,'THIS NUMBER SHOULD BE 12.531 : ',SVAR
      CALL REPORT$STRING(NFIL,'ANISOTROPIC HYPERFINE FIELD FOR ATOM ' &
     &                //TRIM(ATOMNAME)//' IN UNITS OF TESLA')
      CALL REPORT$STRING(NFIL,'MULTIPLY MATRIX WITH A UNITY VECTOR ' &
     &                //'INTO THE DIRECTION OF THE MAGNETIC FIELD')
      CALL LIB$DIAGR8(3,ANIS,EIG,U)              
!      __ FACTOR 0.5 IS THE ELECTRON SPIN____________
      WRITE(NFIL,FMT='(T5,"VALUE",T20,"DIRECTION")')
      DO I=1,3
        WRITE(NFIL,FMT='(T5,F10.5,T20,3F10.5)')SVAR*EIG(I),U(:,I)
      ENDDO

      RETURN
      END
