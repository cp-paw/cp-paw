MODULE CORE_MODULE
!**                                                                   **
!**  the complex pointer list this is needed for the parallelization  **
!**  because every node will calculate the core levels only for its   **
!**  own set of atoms and only the first node can write               **
!**                                                                   **
!**  COMMUNIATION IN CROE$REPORT NOT DONE YET.                        **
!**                                                                   **
type coreshift_type
integer(4)           :: iat
integer(4)           :: n
character(8),pointer :: type(:)
real(8)     ,pointer :: e(:)
real(8)     ,pointer :: eatom(:)
type(coreshift_type),pointer :: next
end type coreshift_type
LOGICAL(4)                  :: TCORESHIFTS=.TRUE.
LOGICAL(4)                  :: DEFAULT=.FALSE.
INTEGER(4)                  :: NATOMS=0
CHARACTER(32),ALLOCATABLE   :: ATOMS(:)
type(coreshift_type),target :: first
type(coreshift_type),pointer:: this
logical(4),save             :: tini=.false.
END MODULE CORE_MODULE
!
!     ..................................................................
      SUBROUTINE CORE$SETL4(ID,VAL)
      USE CORE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(8)  ,INTENT(IN) :: VAL
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
      SUBROUTINE CORE$report(nfil)
      USE CORE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: nfil
      real(8)                 :: ev
      CHARACTER(32)           :: NAME
      INTEGER(4)              :: I,IAT
      REAL(8)                 :: SVAR
      integer(4)              :: thistask,ntasks
!     ******************************************************************
      call mpe$query(ntasks,thistask)
if(ntasks.ne.1) return
      CALL CONSTANTS('EV',EV)
      THIS=>FIRST
      DO 
        if(this%iat.eq.0) exit
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
      SUBROUTINE CORE_CORESHIFTS(IAT,ISP,R1,DEX,NR,LMRXX,AEPOT)
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
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
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
      CHARACTER(32)         :: NAME
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: LMRX

      REAL(8)               :: PI
      REAL(8)               :: XEXP,RI
      INTEGER(4)            :: I,IR
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
!     ==  select the proper entry in the table                         ==
!     ===================================================================
      this=>first
      if(.not.tini) then
        tini=.true.
        this%iat=iat
        this%n=0
        nullify(this%e)
        nullify(this%eatom)
        nullify(this%type)
        nullify(this%next)
      else
        do
          if(this%iat.eq.iat) exit
          if(associated(this%next)) then
            this=>this%next
          else
            allocate(this%next)
            this=>this%next
            nullify(this%next)
            this%iat=iat
            this%n=0
            nullify(this%e)
            nullify(this%eatom)
            nullify(this%type)
          end if
        enddo
      end if
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
      XEXP=DEXP(DEX)

      LMNX=0
      LMRX=0
      DO I=1,NC
        LMNX=LMNX+2*LB(I)+1
        LMRX=MAX(LMRX,(2*LB(I)+1)**2)
      ENDDO
      if(this%n.ne.lmnx) then
        this%n=lmnx
        allocate(this%type(lmnx))
        allocate(this%e(lmnx))
        allocate(this%eatom(lmnx))
      end if

!     ==================================================================
!     == SUBTRACTS ATOMIC AE POTENTIAL FROM AE TOTAL POTENTIAL        ==
!     ==================================================================
      AEPOT1(:,:)=AEPOT(:,:)
      DO IR=1,NR
        AEPOT1(IR,1)=AEPOT(IR,1)-ATPOT(IR)
      ENDDO

!     ==================================================================
!     ==   CALCULATE HAMILTONIAN                                      ==
!     ==================================================================
      XEXP=DEXP(DEX)
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
                    DO IR=1,NR
                      AEDMU(IR)=AEDMU(IR)+CG*AEPOT1(IR,LM3)
                    ENDDO
                END IF
              ENDDO
!
              RI=R1/XEXP
              DO IR=1,NR
                RI=RI*XEXP
                DWORK1(IR)=(AEDMU(IR)*AEPSI(IR,LN1)*AEPSI(IR,LN2))*RI**2
              ENDDO
              CALL RADIAL$INTEGRAL(R1,DEX,NR,DWORK1,SVAR)
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
!     ==  write into table                                            ==
!     ==================================================================
      LMN1=0
      DO LN1=1,NC
        L1=LB(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          this%e(lmn1)=eigenval(lmn1)
          this%eatom(lmn1)=eb(ln1)
          i=l1
          do ln2=1,ln1
            if(lb(ln2).ne.l1) cycle
            i=i+1
          enddo
          write(this%type(lmn1),fmt='(i2)')i
          if(l1.eq.0) then
            this%type(lmn1)=trim(this%type(lmn1))//'s'
          else if(l1.eq.1) then
            this%type(lmn1)=trim(this%type(lmn1))//'p'
          else if(l1.eq.2) then
            this%type(lmn1)=trim(this%type(lmn1))//'d'
          else if(l1.eq.3) then
            this%type(lmn1)=trim(this%type(lmn1))//'f'
          else
            this%type(lmn1)=trim(this%type(lmn1))//'??'
          end if
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
!      CALL CORE$report(nfilo)
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
