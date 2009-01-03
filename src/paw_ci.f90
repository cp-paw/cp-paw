!
!        1         2         3         4         5         6         7         8
!*******************************************************************************
!**  NAME CI                                                                  **
!**                                                                           **
!**    $DELETE                                                                **
!**    $SET(PHI,N,ID,C)                                                       **
!**    _ORDER(PHI)                                                            **
!**    _CLEAN(PHI)                                                            **
!**    $NORMALIZE(PHI)                                                        **
!**    $SCALARPRODUCT(PHI1,PHI2,VAL)                                          **
!**    $ADD(PHI1,PHI2)                                                        **
!**    $SCALE(PHI,C)                                                          **
!**    $CREATE(PHI,IORB)                                                      **
!**    $ANNIHILATE(PHI,IORB)                                                  **
!**    $ORBOCC(PHI,IORB,OCC)                                                  **
!**    $COPY(PHI1,PHI2)                                                       **
!**    $NUMBER(PHI,IORB,OCC)                                                  **
!*******************************************************************************
MODULE CI_MODULE
TYPE CISTATE_TYPE
  INTEGER(4)             :: NX   !X (# SLATER DETERMINANTS)
  INTEGER(4)             :: N    !ACTUAL # SLATER-DETERMINANTS
  COMPLEX(8),POINTER     :: C(:) !COEFFICIENTS
  INTEGER   ,POINTER     :: ID(:)!NUMBER REPRESENTATION IN BIT FORMAT
END TYPE CISTATE_TYPE
TYPE U_TYPE
  INTEGER(4)             :: NX=0   !
  INTEGER(4)             :: N=0    !
  COMPLEX(8),POINTER     :: C(:) !
  INTEGER   ,POINTER     :: I(:)!
  INTEGER   ,POINTER     :: J(:)!
  INTEGER   ,POINTER     :: K(:)!
  INTEGER   ,POINTER     :: L(:)!
END TYPE U_TYPE
TYPE H_TYPE
  INTEGER(4)             :: NX   !
  INTEGER(4)             :: N    !
  COMPLEX(8),POINTER     :: C(:) !
  INTEGER   ,POINTER     :: I(:)!
  INTEGER   ,POINTER     :: J(:)!
END TYPE H_TYPE
TYPE CIHAMIL_TYPE
  TYPE(U_TYPE)            :: U
  TYPE(H_TYPE)            :: H
END TYPE CIHAMIL_TYPE
! MINC IS USED BY CI_COMPACTPSI TO REMOVE SLATER-DETERMINANTS WITH TOO 
! SMALL WEIGHT. 
REAL(8),SAVE              :: CI_MINC=1.D-10  ! MINIMUM ACCEPTABLE COEFFICIENT
real(8),pointer           :: cimat(:,:)
END MODULE CI_MODULE
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$setr8(id,val)
!     **************************************************************************
!     **  CI$setr8                                                            **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      character(*),intent(in) :: id
      REAL(8)     ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'MIN(PSI)') THEN
        CI_MINC=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('CI$SETR8')
      END IF
      RETURN
      END
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$DELETEPSI(PHI)
!     **************************************************************************
!     **  CI$DELETEPSI                                                     **
!     **  DEALLOCATES STATE PHI                                               **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
!     **************************************************************************
      PHI%NX=0
      PHI%N=0
      IF(ASSOCIATED(PHI%C)) THEN
        DEALLOCATE(PHI%C)
        DEALLOCATE(PHI%ID)
      END IF
      RETURN
      END SUBROUTINE CI$DELETEPSI
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$zeroPSI(PHI)
!     **************************************************************************
!     **  CI$SETPSI                                                           **
!     **  ADD A SLATER DETERMINANT TO A STATE                                 **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PHI
      INTEGER(4)                       :: N
!     **************************************************************************
      phi%nx=0
      CALL CI_EXPANDPSI(PHI,20)
      RETURN
      END SUBROUTINE CI$ZEROPSI
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$SETPSI(PHI,ID,C)
!     **************************************************************************
!     **  CI$SETPSI                                                           **
!     **  ADD A SLATER DETERMINANT TO A STATE                                 **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PHI
      INTEGER(4)        ,INTENT(IN)    :: ID
      COMPLEX(8)        ,INTENT(IN)    :: C
      INTEGER(4)                       :: N
!     **************************************************************************
      N=PHI%N
      IF(N+1.GE.PHI%NX) THEN
        CALL CI_EXPANDPSI(PHI,20)
      END IF
      PHI%N=PHI%N+1
      N=PHI%N
      PHI%C(N)=C
      PHI%ID(N)=ID
      RETURN
      END SUBROUTINE CI$SETPSI
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_EXPANDPSI(PHI,NFURTHER)
!     **************************************************************************
!     **  CI_EXPANDPSI                                                        **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)    :: NFURTHER
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PHI
      COMPLEX(8)        ,ALLOCATABLE   :: C(:)
      INTEGER(4)        ,ALLOCATABLE   :: ID(:)
      INTEGER(4)                       :: N
!     **************************************************************************
      IF(PHI%NX.EQ.0) THEN
        PHI%NX=NFURTHER
        PHI%N=0
        ALLOCATE(PHI%ID(NFURTHER))
        ALLOCATE(PHI%C(NFURTHER))
        RETURN
      END IF
      N=PHI%N
      ALLOCATE(C(N))
      ALLOCATE(ID(N))
      C(:)=PHI%C(:N)
      ID(:)=PHI%ID(:N)
      DEALLOCATE(PHI%C)
      DEALLOCATE(PHI%ID)
      PHI%NX=PHI%NX+NFURTHER
      ALLOCATE(PHI%C(PHI%NX))
      ALLOCATE(PHI%ID(PHI%NX))
      PHI%C(:N)=C
      PHI%ID(:N)=ID
      DEALLOCATE(C)
      DEALLOCATE(ID)
      RETURN
      END SUBROUTINE CI_EXPANDPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$CLEANPSI(PHI)
!     **************************************************************************
!     **  CI$CLEANPSI                                                      **
!     **  SORT THE SLATER DETERMINANTS ACCORDING TO ID                        **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT) :: PHI
      INTEGER(4)                        :: N
      INTEGER(4)                        :: ID
      COMPLEX(8)                        :: C
      INTEGER(4)                        :: I,FROM,TO
      REAL(8)           ,ALLOCATABLE    :: CRIT(:)
!     **************************************************************************
      C=(0.D0,0.D0)   ! JUST TO MAKE THE COMPILER HAPPY
      ID=0
!
!     ==========================================================================
!     ==  REMOVE SLATER DETERMINANTS  WITH ZERO WEIGHT                        ==
!     ==========================================================================
      CALL CI_COMPACTPSI(PHI) 
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      N=PHI%N
      ALLOCATE(CRIT(N))
      DO I=1,N
        CRIT(I)=REAL(PHI%ID(I))   ! CRITERION TO ORDER
      ENDDO
      CALL SORT__SET(N,CRIT)
      CALL SORT__RESTART
      CALL SORT__FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          C=PHI%C(FROM)
          ID=PHI%ID(FROM)
        ELSE IF (FROM.EQ.0) THEN
          PHI%C(TO)=C
          PHI%ID(TO)=ID
        ELSE
          PHI%ID(TO)=PHI%ID(FROM)
          PHI%C(TO)=PHI%C(FROM)
        END IF
        CALL SORT__FLIP(FROM,TO)
      ENDDO
      CALL SORT__UNSET
!
!     ==========================================================================
!     ==  REMOVE IDENTICAL SLATER DETERMINANTS                                ==
!     ==========================================================================
      CALL CI_COMPACTPSI(PHI) 
      RETURN
      END SUBROUTINE CI$CLEANPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_COMPACTPSI(PHI)
!     **************************************************************************
!     **  CI_COMPACTPSI                                                      **
!     **  REMOVES STATES WITH ZERO COEFFICIENTS                               **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PHI
      INTEGER(4)               :: N
      INTEGER(4)               :: I,J
!     **************************************************************************
      if(phi%n.eq.0) return
!
!     ==========================================================================
!     == REMOVE SLATER DETERMINANTS WITH ZERO COEFFICIENTS                    ==
!     ==========================================================================
      N=PHI%N
      J=0
      DO I=1,N
        IF(ABS(PHI%C(I)).LE.CI_MINC) CYCLE
        J=J+1
        PHI%C(J)=PHI%C(I)
        PHI%ID(J)=PHI%ID(I)
      ENDDO
      PHI%C(J+1:N)=(0.D0,0.D0)
      PHI%ID(J+1:N)=0
      PHI%N=J
!
!     ==========================================================================
!     == CONTRACT IDENTICAL SLATER DETERMINANTS IN SUBSEQUENT POSITIONS       ==
!     ==========================================================================
      N=PHI%N
      IF(N.LT.2) RETURN   ! CONTRACTION NEEDS AT LEAST TWO DETERMINANTS
      J=1
      DO I=2,N
        IF(PHI%ID(I).EQ.PHI%ID(J)) THEN
          PHI%C(J)=PHI%C(J)+PHI%C(I)
          PHI%C(I)=(0.D0,0.D0)
        ELSE
          J=J+1
          PHI%C(J)=PHI%C(I)
          PHI%ID(J)=PHI%ID(I)
        END IF
      ENDDO
      PHI%C(J+1:N)=(0.D0,0.D0)
      PHI%N=J
      RETURN
      END SUBROUTINE CI_COMPACTPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$NORMALIZE(PHI)
!     **************************************************************************
!     **  CI$NORMALIZE                                                  **
!     **  NORMALIZES THE STATE                                                **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT) :: PHI
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: N
!     **************************************************************************
      CALL CI$CLEANPSI(PHI)
      N=PHI%N
      SVAR=SUM(ABS(PHI%C(1:N))**2)
      SVAR=1.D0/SQRT(SVAR)
      PHI%C(1:N)=PHI%C(1:N)*SVAR
      RETURN
      END SUBROUTINE CI$NORMALIZE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$SCALARPRODUCT(PHI1,PHI2,VAL)
!     **************************************************************************
!     **  CI$SCALARPRODUCT                                                    **
!     **  SCALARPRODUCT <PHI1|PHI2> OF TWO STATES                             **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI1
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI2
      COMPLEX(8)  ,INTENT(OUT)  :: VAL
      INTEGER(4)                :: N1,N2
      INTEGER(4)                :: I1,I2
!     **************************************************************************
      CALL CI$CLEANPSI(PHI1)
      CALL CI$CLEANPSI(PHI2)
      N1=PHI1%N
      N2=PHI2%N
      VAL=(0.D0,0.D0)
      I1=1
      I2=1
      IF(N1.LT.1.OR.N2.LT.1) THEN
        VAL=(0.D0,0.D0)
        RETURN
      END IF
      DO 
        IF(PHI1%ID(I1).EQ.PHI2%ID(I2)) THEN
          VAL=VAL+CONJG(PHI1%C(I1))*PHI2%C(I2)
          I1=I1+1
          I2=I2+1
        ELSE
          IF(PHI1%ID(I1).LT.PHI2%ID(I2)) THEN
            I1=I1+1
          ELSE
            I2=I2+1
          END IF
        END IF
        IF(I1.GT.N1) EXIT
        IF(I2.GT.N2) EXIT
      ENDDO
      RETURN
      END SUBROUTINE CI$SCALARPRODUCT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$ADDPSI(PHI1,PHI2)
!     **************************************************************************
!     **  CI$ADDPSI                                                        **
!     **  ADD PHI2 TO PHI1                                                    **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI1
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI2
      INTEGER(4)                :: I1,I2
      INTEGER(4)                :: N1,N2
      INTEGER(4)                :: NX,N
      INTEGER(4)  ,ALLOCATABLE  :: ID(:)              
      COMPLEX(8)  ,ALLOCATABLE  :: C(:)              
!     **************************************************************************
      CALL CI$CLEANPSI(PHI1)
      CALL CI$CLEANPSI(PHI2)
      N1=PHI1%N
      N2=PHI2%N
!
!     ==========================================================================
!     ==  TREAT SPECIAL CASES, WHEN ON OF THE DETERMINANTS IS ZERO            ==
!     ==========================================================================
      IF(N2.EQ.0) THEN
        RETURN
      END IF
      IF(N1.EQ.0) THEN
        CALL CI$COPYPSI(PHI2,PHI1)
        RETURN
      END IF
!
!     ==========================================================================
!     == COUNT NUMBER OF DETERMINANTS IN RESULTING STATE                      ==
!     ==========================================================================
      I1=1
      I2=1
      NX=0
      DO 
        NX=NX+1
        IF(PHI2%ID(I2).EQ.PHI1%ID(I1)) THEN
          I1=I1+1
          I2=I2+1
        ELSE
          IF(PHI1%ID(I1).LT.PHI2%ID(2)) THEN
            I1=I1+1
          ELSE
            I2=I2+1
          END IF
        END IF
        IF(I1.GT.N1) THEN
          NX=NX+N2-I2+1
          EXIT
        ELSE IF(I2.GT.N2) THEN
          NX=NX+N1-I1+1
          EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     == DETERMINE RESULT
!     ==========================================================================
      ALLOCATE(C(NX))
      ALLOCATE(ID(NX))
      I1=1
      I2=1
      N=0
      DO 
        N=N+1
        IF(PHI2%ID(I2).EQ.PHI1%ID(I1)) THEN
          ID(N)=PHI1%ID(I1)
          C(N)=PHI1%C(I1)+PHI2%C(I2)
          I1=I1+1
          I2=I2+1
        ELSE
          IF(PHI1%ID(I1).LT.PHI2%ID(2)) THEN
            ID(N)=PHI1%ID(I1)
            C(N)=PHI1%C(I1)
            I1=I1+1
          ELSE
            ID(N)=PHI2%ID(I2)
            C(N)=PHI2%C(I2)
            I2=I2+1
          END IF
        END IF
        IF(I1.GT.N1) THEN
          ID(N+1:NX)=PHI2%ID(I2:N2)
          C(N+1:NX)=PHI2%C(I2:N2)
          EXIT
        ELSE IF(I2.GT.N2) THEN
          ID(N+1:NX)=PHI1%ID(I1:N1)
          C(N+1:NX)=PHI1%C(I1:N1)
          EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     == MAP RESULT INTO PHI1                                                 ==
!     ==========================================================================
      IF(NX.GT.PHI1%NX) THEN
        CALL CI_EXPANDPSI(PHI1,NX-PHI1%NX)
      END IF
      PHI1%N=NX
      PHI1%C(1:NX)=C
      PHI1%ID(1:NX)=ID
      DEALLOCATE(ID)
      DEALLOCATE(C)
      RETURN
      END SUBROUTINE CI$ADDPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$SCALEPSI(PHI,C)
!     **************************************************************************
!     **  CI$DSCALEPSI                                                        **
!     **  MULTIPLY PSI WITH A COMPLEX NUMBER                                  **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
      COMPLEX(8)  ,INTENT(IN)   :: C
      INTEGER(4)                :: N
!     **************************************************************************
      N=PHI%N
      PHI%C(1:N)=PHI%C(1:N)*C
      RETURN
      END SUBROUTINE CI$SCALEPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$CREATOR(PHI,IORB)
!     **************************************************************************
!     **  CI$CREATOR                                                          **
!     **  APPLY A CREATION OPERATOR FOR ORBITAL IORB TO PHI                   **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
      INTEGER(4)  ,INTENT(IN)   :: IORB
      INTEGER(4)                :: N,I,j
      LOGICAL(4)                :: TOCC
      LOGICAL(4)                :: Tminus
!     **************************************************************************
      IF(IORB.GE.BIT_SIZE(PHI%ID(1))) THEN
        STOP 'IORB OUT OF RANGE'
      END IF
      N=PHI%N
      DO I=1,N
        TOCC=BTEST(PHI%ID(I)-1,IORB-1) 
        IF(TOCC) THEN
          PHI%C(I)=(0.D0,0.D0)
        ELSE
          PHI%ID(I)=1+IBSET(PHI%ID(I)-1,IORB-1)
!         == DETERMINE SIGN CHANGE BY THE NUMBER OF PERMUTATIONS ===============
          TMINUS=.FALSE.
          DO J=1,IORB-1
            IF(BTEST(PHI%ID(I)-1,J-1))TMINUS=.NOT.TMINUS
          ENDDO
          IF(TMINUS)PHI%C(I)=-PHI%C(I)
        END IF
      ENDDO
      CALL CI_COMPACTPSI(PHI)
      RETURN
      END SUBROUTINE CI$CREATOR
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$ANNIHILATOR(PHI,IORB)
!     **************************************************************************
!     **  CI$ANNIHILATOR                                                      **
!     **  APPLY ANNIHILATION NOPERATOR FOR ORBITAL IORB TO PHI                **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
      INTEGER(4)  ,INTENT(IN)   :: IORB
      INTEGER(4)                :: N,I,j
      LOGICAL(4)                :: TOCC
      LOGICAL(4)                :: Tminus
!     **************************************************************************
      IF(IORB.GE.BIT_SIZE(PHI%ID(1))) THEN
        STOP 'IORB OUT OF RANGE'
      END IF
      N=PHI%N
      DO I=1,N
        TOCC= BTEST(PHI%ID(I)-1,IORB-1) 
        IF(TOCC) THEN
          PHI%ID(I)=1+IBCLR(PHI%ID(I)-1,IORB-1)
!         == DETERMINE SIGN CHANGE BY THE NUMBER OF PERMUTATIONS ===============
          TMINUS=.FALSE.
          DO J=1,IORB-1
            IF(BTEST(PHI%ID(I)-1,J-1)) TMINUS=.NOT.TMINUS
          ENDDO
          IF(TMINUS)PHI%C(I)=-PHI%C(I)
        ELSE
          PHI%C(I)=(0.D0,0.D0)
        END IF
      ENDDO
      CALL CI_COMPACTPSI(PHI)
      RETURN
      END SUBROUTINE CI$ANNIHILATOR
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$COPYPSI(PHI1,PHI2)
!     **************************************************************************
!     **  CI$COPYPSI                                                          **
!     **  COPY PHI1 INTO PHI2                                                 **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI1
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI2
!     **************************************************************************
      CALL CI$DELETEPSI(PHI2)
      PHI2%NX=PHI1%NX
      PHI2%N=PHI1%N
      ALLOCATE(PHI2%ID(PHI1%NX))
      ALLOCATE(PHI2%C(PHI1%NX))
      PHI2%ID=PHI1%ID
      PHI2%C=PHI1%C
      RETURN
      END SUBROUTINE CI$COPYPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$GETORBOCC(PHI,IORB,OCC)
!     **************************************************************************
!     **  CI$GETORBOCC                                                        **
!     **  EXPECTATION VALUE OF THE OCUPATION OF A ONE-PARTICLE ORBITAL        **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
      INTEGER(4)  ,INTENT(IN)   :: IORB
      REAL(8)     ,INTENT(OUT)  :: OCC 
      INTEGER(4)                :: I,N
      LOGICAL(4)                :: TOCC
      REAL(8)                   :: RNORM
!     **************************************************************************
      IF(IORB.GE.BIT_SIZE(PHI%ID(1))) THEN
        STOP 'IORB OUT OF RANGE'
      END IF
      CALL CI$CLEANPSI(PHI)
      N=PHI%N
      OCC=0.D0
      RNORM=0.D0
      DO I=1,N
        TOCC=BTEST(PHI%ID(I)-1,IORB-1)
        IF(TOCC)OCC=OCC+ABS(PHI%C(I))**2
        RNORM=RNORM+ABS(PHI%C(I))**2
      ENDDO
      OCC=OCC/RNORM
      RETURN
      END SUBROUTINE CI$GETORBOCC
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$GETNPARTICLE(PHI,OCC)
!     **************************************************************************
!     **  CI$GETNPARTICLE                                                     **
!     **  DETERMINES THE EXPECTATION VALUE OF THE PARTICLE NUMBER             **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
      REAL(8)     ,INTENT(OUT)  :: OCC 
      INTEGER(4)                :: NORB,ICOUNT,I,J,N
      REAL(8)                   :: RNORM
!     **************************************************************************
      CALL CI$CLEANPSI(PHI)
      N=PHI%N
      NORB=BIT_SIZE(PHI%ID(1))
      OCC=0.D0
      RNORM=0.D0
      DO I=1,N
        ICOUNT=0
        DO J=1,NORB
          IF(BTEST(PHI%ID(I)-1,J-1))ICOUNT=ICOUNT+1
        ENDDO
        OCC=OCC+REAL(ICOUNT)*ABS(PHI%C(I))**2
        RNORM=RNORM+ABS(PHI%C(I))**2
      ENDDO
      OCC=OCC/RNORM
      RETURN
      END SUBROUTINE CI$GETNPARTICLE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CISTATES$NTOTPSI(PHI)
!     **************************************************************************
!     **  CISTATES$OCC                                                        **
!     **  DETERMINES THE EXPECTATION VALUE OF THE PARTILE NUMBER              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
      INTEGER(4)                       :: N,NORB,I,J,ICOUNT
!     **************************************************************************
      N=PHI%N
      NORB=BIT_SIZE(PHI%ID(1))
      DO I=1,N
        ICOUNT=0
        DO J=1,NORB
          IF(BTEST(PHI%ID(I)-1,J-1))ICOUNT=ICOUNT+1
        ENDDO
        PHI%C(I)=REAL(ICOUNT,KIND=8)*PHI%C(I)
      ENDDO
      RETURN
      END SUBROUTINE CISTATES$NTOTPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$WRITEPSI(PHI,NFIL)
!     **************************************************************************
!     **  CI$WRITEPSI                                                      **
!     **  WRITE STATE TO FILE                                                 **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI
      INTEGER(4)  ,INTENT(IN)   :: NFIL
      INTEGER(4)                :: N
      INTEGER(4)  ,ALLOCATABLE  :: NUMREP(:)
      INTEGER(4)                :: LENGTH
      INTEGER(4)                :: I,J
!     **************************************************************************
      LENGTH=BIT_SIZE(PHI%ID(1))
      IF(LENGTH.NE.32) THEN
        CALL ERROR$STOP('CI$WRITEPSI')
      END IF
      ALLOCATE(NUMREP(LENGTH))
      N=PHI%N
      WRITE(nfil,*)'========================== PSI ============================'
      DO I=1,N
        NUMREP(:)=0
        DO J=1,LENGTH
          IF(BTEST(PHI%ID(I)-1,J-1)) NUMREP(J)=1
        ENDDO
        WRITE(NFIL,FMT='(I5,"|",32I1,"> * ",2F10.7)')I,NUMREP,PHI%C(I)
      ENDDO
      DEALLOCATE(NUMREP)
      RETURN
      END SUBROUTINE CI$WRITEPSI
!
!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**   HAMILTON OPERATIONS                                                     **
!**                                                                           **
!*******************************************************************************
!*******************************************************************************
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$SETU(HAM,I,J,K,L,VAL)
!     **************************************************************************
!     **  CIHAMI$ADDTOU                                                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
      INTEGER(4)  ,INTENT(IN)   :: I
      INTEGER(4)  ,INTENT(IN)   :: J
      INTEGER(4)  ,INTENT(IN)   :: K
      INTEGER(4)  ,INTENT(IN)   :: L
      COMPLEX(8)  ,INTENT(IN)   :: VAL
!     **************************************************************************
      CALL CI_SETUELEMENT(HAM%U,I,J,K,L,VAL)
      RETURN 
      END SUBROUTINE CI$SETU
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$SETH(HAM,I,J,VAL)
!     **************************************************************************
!     **  CIHAMI$ADDTOU                                                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
      INTEGER(4)  ,INTENT(IN)   :: I
      INTEGER(4)  ,INTENT(IN)   :: J
      COMPLEX(8)  ,INTENT(IN)   :: VAL
!     **************************************************************************
      CALL CI_SETHELEMENT(HAM%H,I,J,VAL)
      RETURN 
      END SUBROUTINE CI$SETH
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$WRITEHAMILTONIAN(HAM,NFIL)
!     **************************************************************************
!     **  CIHAMI$ADDTOU                                                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
      INTEGER(4)  ,INTENT(IN)   :: NFIL
!     **************************************************************************
      CALL CI_WRITEH(HAM%H,NFIL)
      CALL CI_WRITEU(HAM%U,NFIL)
      RETURN 
      END SUBROUTINE CI$WRITEHAMILTONIAN
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$CLEANHAMILTONIAN(HAM)
!     **************************************************************************
!     **  CI$CLEANHAMILTONIAN                                                 **
!     **  BRING MATRIX ELEMENTS INTO OPTIMUM ORDER                            **
!     **  REMOVE ZERO ELEMENTS                                                **
!     **  COMBINE IDENTICAL ELEMENTS                                          **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
!     **************************************************************************
CALL TIMING$CLOCKON('CI$CLEANH')
      CALL CI_CLEANU(HAM%U)
      CALL CI_CLEANH(HAM%H)
CALL TIMING$CLOCKOFF('CI$CLEANH')
      RETURN 
      END SUBROUTINE CI$CLEANHAMILTONIAN
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$deleteHAMILTONIAN(HAM)
!     **************************************************************************
!     **  CI$deleteHAMILTONIAN                                                **
!     **  deallocate all arrays                                               **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
!     **************************************************************************
      if(associated(ham%u%c)) then
        deallocate(ham%u%c)
        deallocate(ham%u%i)
        deallocate(ham%u%j)
        deallocate(ham%u%k)
        deallocate(ham%u%l)
      end if
      ham%u%nx=0
      ham%u%n=0
      if(associated(ham%H%c)) then
        deallocate(ham%H%c)
        deallocate(ham%H%i)
        deallocate(ham%H%j)
      end if
      ham%H%nx=0
      ham%H%n=0
      RETURN 
      END SUBROUTINE CI$DELETEHAMILTONIAN
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_SETUELEMENT(U,I,J,K,L,VAL)
!     **************************************************************************
!     **  CI_SETUELEMENT                                                      **
!     **  SET ONE ELEMENT FOR THE U-TENSOR                                    **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(U_TYPE),INTENT(INOUT):: U
      INTEGER(4)  ,INTENT(IN)   :: I
      INTEGER(4)  ,INTENT(IN)   :: J
      INTEGER(4)  ,INTENT(IN)   :: K
      INTEGER(4)  ,INTENT(IN)   :: L
      COMPLEX(8)  ,INTENT(IN)   :: VAL
      INTEGER(4)                :: N
!     **************************************************************************
      IF(U%N.GE.U%NX) CALL CI_CLEANU(U)
      IF(U%N.GE.U%NX) CALL CI_EXPANDU(U,20)
      U%N=U%N+1
      N=U%N
      U%I(N)=I
      U%J(N)=J
      U%K(N)=K
      U%L(N)=L
      U%C(N)=VAL
      RETURN
      END SUBROUTINE CI_SETUELEMENT
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_WRITEU(U,NFIL)
!     **************************************************************************
!     **  CI_WRITEU                                                           **
!     **  WRITE U-TENSOR TO FILE                                              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(U_TYPE),INTENT(INOUT):: U
      INTEGER(4)  ,INTENT(IN)   :: NFIL
      INTEGER(4)                :: N,I
!     **************************************************************************
      N=U%N
      WRITE(NFIL,FMT='("============== U  =======================")')
      DO I=1,N
        WRITE(NFIL,FMT='(4I5,2F17.10)')U%I(I),U%J(I),U%K(I),U%L(I),U%C(I)
      ENDDO
      RETURN
      END SUBROUTINE CI_WRITEU
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_WRITEH(H,NFIL)
!     **************************************************************************
!     **  WRITE ONE PARTICLE HAMILTONIAN TO FILE                              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(H_TYPE),INTENT(INOUT):: H
      INTEGER(4)  ,INTENT(IN)   :: NFIL
      INTEGER(4)                :: N,I
!     **************************************************************************
      N=H%N
      WRITE(NFIL,FMT='("============== H  =======================")')
      DO I=1,N
        WRITE(NFIL,FMT='(2I5,2F17.10)')H%I(I),H%J(I),H%C(I)
      ENDDO
      RETURN
      END SUBROUTINE CI_WRITEH
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_SETHELEMENT(H,I,J,VAL)
!     **************************************************************************
!     **  SET ONE MATRIX ELEMENT OF ONE-PARTICLE HAMILTONIAN                  **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(H_TYPE),INTENT(INOUT):: H
      INTEGER(4)  ,INTENT(IN)   :: I
      INTEGER(4)  ,INTENT(IN)   :: J
      COMPLEX(8)  ,INTENT(IN)   :: VAL
      INTEGER(4)                :: N
!     **************************************************************************
      IF(H%N.GE.H%NX) CALL CI_CLEANH(H)
      IF(H%N.GE.H%NX) CALL CI_EXPANDH(H,20)
      H%N=H%N+1
      N=H%N
      H%I(N)=I
      H%J(N)=J
      H%C(N)=VAL
      RETURN
      END SUBROUTINE CI_SETHELEMENT
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_EXPANDU(U,NFURTHER)
!     **************************************************************************
!     **  ADD NFURTHER ELEMENTS FOR THE U-TENSOR                              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(U_TYPE) ,INTENT(INOUT):: U
      INTEGER(4)   ,INTENT(IN)   :: NFURTHER
      COMPLEX(8)   ,ALLOCATABLE  :: C(:)
      INTEGER(4)   ,ALLOCATABLE  :: IND(:,:)
      INTEGER(4)                 :: N,NX
!     **************************************************************************
      IF(U%NX.EQ.0) THEN
        U%NX=NFURTHER
        U%N=0
        NX=NFURTHER
        ALLOCATE(U%I(NX))
        ALLOCATE(U%J(NX))
        ALLOCATE(U%K(NX))
        ALLOCATE(U%L(NX))
        ALLOCATE(U%C(NX))
        RETURN
      END IF
      N=U%N
      ALLOCATE(IND(N,4))
      ALLOCATE(C(N))
      IND(:,1)=U%I(:N)
      IND(:,2)=U%J(:N)
      IND(:,3)=U%K(:N)
      IND(:,4)=U%L(:N)
      C(:)=U%C(:N)
      DEALLOCATE(U%I)
      DEALLOCATE(U%J)
      DEALLOCATE(U%K)
      DEALLOCATE(U%L)
      DEALLOCATE(U%C)
      U%NX=U%NX+NFURTHER
      NX=U%NX
      ALLOCATE(U%I(NX))
      ALLOCATE(U%J(NX))
      ALLOCATE(U%K(NX))
      ALLOCATE(U%L(NX))
      ALLOCATE(U%C(NX))
      U%I(:N)=IND(:,1)
      U%J(:N)=IND(:,2)
      U%K(:N)=IND(:,3)
      U%L(:N)=IND(:,4)
      U%C(:N)=C(:)
      DEALLOCATE(IND)
      DEALLOCATE(C)
      RETURN
      END SUBROUTINE CI_EXPANDU
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_EXPANDH(H,NFURTHER)
!     **************************************************************************
!     **  ADD NFURTHER ELEMENTS FOR THE ONE-PARTICLE HAMILTONIAN              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(H_TYPE) ,INTENT(INOUT):: H
      INTEGER(4)   ,INTENT(IN)   :: NFURTHER
      COMPLEX(8)   ,ALLOCATABLE  :: C(:)
      INTEGER(4)   ,ALLOCATABLE  :: IND(:,:)
      INTEGER(4)                 :: N,NX
!     **************************************************************************
      IF(H%NX.EQ.0) THEN
        H%NX=NFURTHER
        H%N=0
        NX=NFURTHER
        ALLOCATE(H%I(NX))
        ALLOCATE(H%J(NX))
        ALLOCATE(H%C(NX))
        RETURN
      END IF
      N=H%N
      ALLOCATE(IND(N,2))
      ALLOCATE(C(N))
      IND(:,1)=H%I(:N)
      IND(:,2)=H%J(:N)
      C(:)=H%C(:N)
      DEALLOCATE(H%I)
      DEALLOCATE(H%J)
      DEALLOCATE(H%C)
      H%NX=H%NX+NFURTHER
      NX=H%NX
      ALLOCATE(H%I(NX))
      ALLOCATE(H%J(NX))
      ALLOCATE(H%C(NX))
      H%I(:N)=IND(:,1)
      H%J(:N)=IND(:,2)
      H%C(:N)=C(:)
      DEALLOCATE(IND)
      DEALLOCATE(C)
      RETURN
      END SUBROUTINE CI_EXPANDH
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_CLEANU(U)
!     **************************************************************************
!     **  CI_CLEANU                                                           **
!     **  REMOVES ZERO ELEMENTS OF THE U-TENSOR, ORDER THE ELEMENTS           **
!     **  AND ADD IDENTICAL ELEMENTS                                          **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(U_TYPE) ,INTENT(INOUT):: U
      INTEGER(4)                 :: N
      INTEGER(4)                 :: I,J,K,L,I1
      INTEGER(4)                 :: ISVAR
      COMPLEX(8)                 :: C
      REAL(8)  ,ALLOCATABLE      :: CRIT(:)
      INTEGER(4)                 :: FROM,TO
      INTEGER(4)                 :: ISTART,IEND,LENGTH
      LOGICAL(4)                 :: TDO,TDIFF
!     **************************************************************************
      IF(U%N.EQ.0) RETURN
      N=U%N
      C=(0.D0,0.D0) ! JUST TO MAKE THE COMPILER HAPPY
      L=0
      K=0
!
!     ==========================================================================
!     ==  ENSURE I>J AND K>L                                                  ==
!     ==========================================================================
      DO I=1,N
        IF(U%J(I).GT.U%I(I)) THEN
          ISVAR=U%I(I)
          U%I(I)=U%J(I)          
          U%J(I)=ISVAR
          U%C(I)=-U%C(I)
        END IF
        IF(U%L(I).GT.U%K(I)) THEN
          ISVAR=U%K(I)
          U%K(I)=U%L(I)          
          U%L(I)=ISVAR
          U%C(I)=-U%C(I)
        END IF
        IF(U%I(I).EQ.U%J(I))U%C(I)=(0.D0,0.D0) !CREATOR SQUARED IS ZERO
        IF(U%K(I).EQ.U%L(I))U%C(I)=(0.D0,0.D0) !ANNIHILATOR SQUARED IS ZERO
      ENDDO
!
!     ==========================================================================
!     ==  REMOVE ZERO ELEMENTS                                                ==
!     ==========================================================================
      J=0
      DO I=1,N
        IF(ABS(U%C(I)).EQ.0.D0) CYCLE
        J=J+1
        U%C(J)=U%C(I)
        U%I(J)=U%I(I)
        U%J(J)=U%J(I)
        U%K(J)=U%K(I)
        U%L(J)=U%L(I)
      END DO
      U%N=J
      N=U%N
!
!     ==========================================================================
!     ==  ORDER ENTRIES                                                       ==
!     ==========================================================================
!     ==  ORDER WITH RESPECT TO LAST INDEX
      ALLOCATE(CRIT(N))
      DO I=1,N
        CRIT(I)=REAL(U%L(I))
      ENDDO
      CALL SORT__SET(N,CRIT)
      CALL SORT__RESTART
      CALL SORT__FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          I=U%I(FROM)
          J=U%J(FROM)
          K=U%K(FROM)
          L=U%L(FROM)
          C=U%C(FROM)
        ELSE IF (FROM.EQ.0) THEN
          U%I(TO)=I
          U%J(TO)=J
          U%K(TO)=K
          U%L(TO)=L
          U%C(TO)=C
        ELSE
          U%I(TO)=U%I(FROM)
          U%J(TO)=U%J(FROM)
          U%K(TO)=U%K(FROM)
          U%L(TO)=U%L(FROM)
          U%C(TO)=U%C(FROM)
        END IF
        CALL SORT__FLIP(FROM,TO)
      ENDDO
      CALL SORT__UNSET
!
!     == ORDER WITH RESPECT TO SECOND INDEX
      ISTART=1
      DO I1=1,N
        CRIT(I1-ISTART+1)=REAL(U%K(I1))
        TDO=(I1.EQ.N)
        IF(.NOT.TDO) TDO=(U%L(I1+1).NE.U%L(I1))
        IF(TDO) THEN
          IEND=I1
          LENGTH=IEND-ISTART+1
          CALL SORT__SET(LENGTH,CRIT)
          CALL SORT__RESTART
          CALL SORT__FLIP(FROM,TO)
          DO WHILE (FROM.NE.0.OR.TO.NE.0)
            IF(TO.EQ.0) THEN
              from=from+istart-1
              I=U%I(FROM)
              J=U%J(FROM)
              K=U%K(FROM)
              L=U%L(FROM)
              C=U%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              to=to+istart-1
              U%I(TO)=I
              U%J(TO)=J
              U%K(TO)=K
              U%L(TO)=L
              U%C(TO)=C
            ELSE
              from=from+istart-1
              to=to+istart-1
              U%I(TO)=U%I(FROM)
              U%J(TO)=U%J(FROM)
              U%K(TO)=U%K(FROM)
              U%L(TO)=U%L(FROM)
              U%C(TO)=U%C(FROM)
            END IF
            CALL SORT__FLIP(FROM,TO)
          ENDDO
          CALL SORT__UNSET
          ISTART=IEND+1
        END IF
      ENDDO
!
!     == ORDER WITH RESPECT TO THIRD INDEX
      ISTART=1
      DO I1=1,N
        CRIT(I1-ISTART+1)=REAL(U%J(I1))
        TDO=(I1.EQ.N)
        IF(.NOT.TDO)TDO=(U%L(I1+1).NE.U%L(I1)).OR.(U%K(I1+1).NE.U%K(I1))
        IF(TDO) THEN
          IEND=I1
          LENGTH=IEND-ISTART+1
          CALL SORT__SET(LENGTH,CRIT)
          CALL SORT__RESTART
          CALL SORT__FLIP(FROM,TO)
          DO WHILE (FROM.NE.0.OR.TO.NE.0)
            IF(TO.EQ.0) THEN
              from=from+istart-1
              I=U%I(FROM)
              J=U%J(FROM)
              K=U%K(FROM)
              L=U%L(FROM)
              C=U%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              to=to+istart-1
              U%I(TO)=I
              U%J(TO)=J
              U%K(TO)=K
              U%L(TO)=L
              U%C(TO)=C
            ELSE
              from=from+istart-1
              to=to+istart-1
              U%I(TO)=U%I(FROM)
              U%J(TO)=U%J(FROM)
              U%K(TO)=U%K(FROM)
              U%L(TO)=U%L(FROM)
              U%C(TO)=U%C(FROM)
            END IF
            CALL SORT__FLIP(FROM,TO)
          ENDDO
          CALL SORT__UNSET
          ISTART=IEND+1
        END IF
      ENDDO
!
!     == ORDER WITH RESPECT TO LAST INDEX
      ISTART=1
      DO I1=1,N
        CRIT(I1-ISTART+1)=REAL(U%I(I1))
        TDO=(I1.EQ.N)
        IF(.NOT.TDO)TDO=(U%L(I1+1).NE.U%L(I1)) &
     &              .OR.(U%K(I1+1).NE.U%K(I1)) &
     &              .OR.(U%J(I1+1).NE.U%J(I1))
        IF(TDO) THEN
          IEND=I1
          LENGTH=IEND-ISTART+1
          CALL SORT__SET(LENGTH,CRIT)
          CALL SORT__RESTART
          CALL SORT__FLIP(FROM,TO)
          DO WHILE (FROM.NE.0.OR.TO.NE.0)
            IF(TO.EQ.0) THEN
              from=from+istart-1
              I=U%I(FROM)
              J=U%J(FROM)
              K=U%K(FROM)
              L=U%L(FROM)
              C=U%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              to=to+istart-1
              U%I(TO)=I
              U%J(TO)=J
              U%K(TO)=K
              U%L(TO)=L
              U%C(TO)=C
            ELSE
              from=from+istart-1
              to=to+istart-1
              U%I(TO)=U%I(FROM)
              U%J(TO)=U%J(FROM)
              U%K(TO)=U%K(FROM)
              U%L(TO)=U%L(FROM)
              U%C(TO)=U%C(FROM)
            END IF
            CALL SORT__FLIP(FROM,TO)
          ENDDO
          CALL SORT__UNSET
          ISTART=IEND+1
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  COMBINE IDENTICAL ENTRIES                                           ==
!     ==========================================================================
      J=1
      DO I=2,N
        TDIFF=(J.EQ.0)
        IF(.NOT.TDIFF) THEN
          TDIFF=(U%I(I).NE.U%I(J))
          IF(.NOT.TDIFF) THEN
            TDIFF=(U%J(I).NE.U%J(J))
            IF(.NOT.TDIFF) THEN
              TDIFF=(U%K(I).NE.U%K(J))
              IF(.NOT.TDIFF) THEN
                TDIFF=(U%L(I).NE.U%L(J))
              END IF
            END IF
          END IF
        END IF
        IF(TDIFF) THEN
          J=J+1
          U%I(J)=U%I(I)
          U%J(J)=U%J(I)
          U%K(J)=U%K(I)
          U%L(J)=U%L(I)
          U%C(J)=U%C(I)
        ELSE
          U%C(J)=U%C(J)+U%C(I)
          U%C(I)=(0.D0,0.D0)
          IF(ABS(U%C(J)).EQ.0.D0) THEN
            J=J-1
          END IF
        END IF
      ENDDO
      U%N=J
      RETURN
      END SUBROUTINE CI_CLEANU
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_CLEANH(H)
!     **************************************************************************
!     **  CI_CLEANH                                                           **
!     **  REMOVES ZERO ELEMENTS OF THE ONE-PARTICLE HAMILTONIAN               **
!     **  ORDER THE ELEMENTS AND ADD IDENTICAL ELEMENTS                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(H_TYPE) ,INTENT(INOUT):: H
      INTEGER(4)                 :: N
      INTEGER(4)                 :: I,J,I1
      COMPLEX(8)                 :: C
      REAL(8)  ,ALLOCATABLE      :: CRIT(:)
      INTEGER(4)                 :: FROM,TO
      INTEGER(4)                 :: ISTART,IEND,LENGTH
      LOGICAL(4)                 :: TDO,TDIFF
!     **************************************************************************
      IF(H%N.EQ.0) RETURN
      N=H%N
      C=(0.D0,0.D0)
!
!     ==========================================================================
!     ==  REMOVE ZERO ELEMENTS                                                ==
!     ==========================================================================
      J=0
      DO I=1,N
        IF(ABS(H%C(I)).EQ.0.D0) CYCLE
        J=J+1
        H%C(J)=H%C(I)
        H%I(J)=H%I(I)
        H%J(J)=H%J(I)
      END DO
      H%N=J
      N=H%N
!
!     ==========================================================================
!     ==  ORDER ENTRIES                                                       ==
!     ==========================================================================
!     ==  ORDER WITH RESPECT TO LAST INDEX
      ALLOCATE(CRIT(N))
      DO I=1,N
        CRIT(I)=REAL(H%J(I))
      ENDDO
      CALL SORT__SET(N,CRIT)
      CALL SORT__RESTART
      CALL SORT__FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          I=H%I(FROM)
          J=H%J(FROM)
          C=H%C(FROM)
        ELSE IF (FROM.EQ.0) THEN
          H%I(TO)=I
          H%J(TO)=J
          H%C(TO)=C
        ELSE
          H%I(TO)=H%I(FROM)
          H%J(TO)=H%J(FROM)
          H%C(TO)=H%C(FROM)
        END IF
        CALL SORT__FLIP(FROM,TO)
      ENDDO
      CALL SORT__UNSET
!
!     == ORDER WITH RESPECT TO SECOND INDEX
      ISTART=1
      DO I1=1,N
        CRIT(I1-ISTART+1)=REAL(H%I(I1))
        TDO=(I1.EQ.N)
        IF(.NOT.TDO) TDO=(H%J(I1+1).NE.H%J(I1))
        IF(TDO) THEN
          IEND=I1
          LENGTH=IEND-ISTART+1
          CALL SORT__SET(LENGTH,CRIT)
          CALL SORT__RESTART
          CALL SORT__FLIP(FROM,TO)
          DO WHILE (FROM.NE.0.OR.TO.NE.0)
            IF(TO.EQ.0) THEN
              from=from+istart-1
              I=H%I(FROM)
              J=H%J(FROM)
              C=H%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              to=to+istart-1
              H%I(TO)=I
              H%J(TO)=J
              H%C(TO)=C
            ELSE
              from=from+istart-1
              to=to+istart-1
              H%I(TO)=H%I(FROM)
              H%J(TO)=H%J(FROM)
              H%C(TO)=H%C(FROM)
            END IF
            CALL SORT__FLIP(FROM,TO)
          ENDDO
          CALL SORT__UNSET
          ISTART=IEND+1
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  COMBINE IDENTICAL ENTRIES                                           ==
!     ==========================================================================
      J=1
      DO I=2,N
        TDIFF=(J.EQ.0)
        IF(.NOT.TDIFF) THEN
          TDIFF=(H%I(I).NE.H%I(J))
          IF(.NOT.TDIFF) THEN
            TDIFF=(H%J(I).NE.H%J(J))
          END IF
        END IF
        IF(TDIFF) THEN
          J=J+1
          H%I(J)=H%I(I)
          H%J(J)=H%J(I)
          H%C(J)=H%C(I)
        ELSE
          H%C(J)=H%C(J)+H%C(I)
          H%C(I)=(0.D0,0.D0)
          IF(ABS(H%C(J)).EQ.0.D0) THEN
            J=J-1
          END IF
        END IF
      ENDDO
      H%N=J
      RETURN
      END SUBROUTINE CI_CLEANH
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$HPSI(HAM,PSI)
!     **************************************************************************
!     **  CI$HPSI                                                             **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(IN)    :: HAM
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PSI
      TYPE(CISTATE_TYPE)               :: PSI1
      TYPE(CISTATE_TYPE)               :: PSI2
      TYPE(CISTATE_TYPE)               :: PSI3
      TYPE(CISTATE_TYPE)               :: PSI4
      TYPE(CISTATE_TYPE)               :: HPSI
      INTEGER(4)                       :: N
      INTEGER(4)                       :: I,J,K,L,I1
!     **************************************************************************
!
!     ==========================================================================
!     == APPLY U-TENSOR                                                       ==
!     == W=0.5d0*sum_{i,j,k,l} W_{i,j,k,l} cdagger_i cdagger_j c_l c_k        ==
!     ==========================================================================
      call ci$zeropsi(hpsi)
      N=HAM%U%N
      I=0
      J=0
      K=0
      L=0
      DO I1=1,N
        IF(HAM%U%k(I1).NE.k) THEN
          I=0
          J=0
          l=0
          k=HAM%U%k(I1)
          CALL CI$COPYPSI(PSI,PSI1)
          CALL CI$ANNIHILATOR(PSI1,k)    ! |psi1>=c_k|psi>
        END IF
        IF(HAM%U%l(I1).NE.l) THEN
          I=0
          J=0
          l=HAM%U%l(I1)
          CALL CI$COPYPSI(PSI1,PSI2)
          CALL CI$ANNIHILATOR(PSI2,l)    ! |psi2>=c_l c_k |psi>
        END IF
        IF(HAM%U%J(I1).NE.J) THEN
          I=0
          J=HAM%U%J(I1)
          CALL CI$COPYPSI(PSI2,PSI3)
          CALL CI$CREATOR(PSI3,J)        ! |psi3>=cdagger_j c_l c_k |psi>
        END IF
        I=HAM%U%I(I1)
        CALL CI$COPYPSI(PSI3,PSI4)
        CALL CI$CREATOR(PSI4,I)       ! |psi4>=cdagger_i cdagger_j c_l c_k |psi>
        CALL CI$SCALEPSI(PSI4,0.5d0*HAM%U%C(I1))  
        CALL CI$ADDPSI(HPSI,PSI4)
      ENDDO
!
!     ==========================================================================
!     == APPLY ONE-PARTICLE HAMILTONIAN                                       ==
!     ==========================================================================
      N=HAM%H%N
      I=0
      J=0
      DO I1=1,N
        IF(HAM%H%J(I1).NE.J) THEN
          I=0
          J=HAM%H%J(I1)
          CALL CI$COPYPSI(PSI,PSI1)
          CALL CI$ANNIHILATOR(PSI1,J)
        END IF
        I=HAM%H%I(I1)
        CALL CI$COPYPSI(PSI1,PSI2)
        CALL CI$CREATOR(PSI2,I)
        CALL CI$SCALEPSI(PSI2,HAM%H%C(I1))
        CALL CI$ADDPSI(HPSI,PSI2)
      ENDDO
!
!     ==========================================================================
!     == APPLY ONE-PARTICLE HAMILTONIAN                                       ==
!     ==========================================================================
      CALL CI$COPYPSI(HPSI,PSI)
!
!     ==========================================================================
!     == DESTROY INTERMEDIATE STATES                                          ==
!     ==========================================================================
      CALL CI$DELETEPSI(PSI1)
      CALL CI$DELETEPSI(PSI2)
      CALL CI$DELETEPSI(PSI3)
      CALL CI$DELETEPSI(PSI4)
      CALL CI$DELETEPSI(HPSI)
      RETURN
      END
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$1POCCPSI(PSI,IORB)
!     **************************************************************************
!     **  CIHAMI_ADDTOU                                                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PSI
      INTEGER(4)        ,INTENT(IN)    :: IORB
!     **************************************************************************
      CALL CI$ANNIHILATOR(PSI,IORB)
      CALL CI$CREATOR(PSI,IORB)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$STATEFROMDENSITYMATRIX(NCHI,RHO,PSI)
!     **************************************************************************
!     ** CONSTRUCTS A MANY-PARTICLE WAVE FUNCTION SUCH THAT ITS ONE-PARTICLE  **
!     ** DENSITY MATRIX IS EQUAL TO A SPECIFIED ONE (RHO)                     **
!     **                                                                      **
!     ** THE WAVE FUNCTIONS CONSISTS OF A STATISTICAL AVERAGE OF N-PARTICLE   **
!     ** SLATER DETERMINANTS FOR DIFFERENT N                                  **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4),        INTENT(IN)   :: NCHI
      COMPLEX(8),        INTENT(IN)   :: RHO(NCHI,NCHI) ! DENSITY MATRIX
      TYPE(CISTATE_TYPE),INTENT(OUT)  :: PSI
      COMPLEX(8)                      :: U(NCHI,NCHI)
      REAL(8)                         :: Fn(0:NCHI+1)
      REAL(8)                         :: EIG(NCHI)
      TYPE(CISTATE_TYPE)              :: PSIN
      TYPE(CISTATE_TYPE)              :: PSIADD2
      TYPE(CISTATE_TYPE)              :: PSICOPY
      complex(8)                      :: csvar
      real(8)                         :: diff
      integer(4)                      :: i,n,j,k
      logical(4)       ,parameter     :: ttest=.true.
!     **************************************************************************
!
!     ==========================================================================
!     == DIAGONALIZE DENSITY MATRIX                                           ==
!     ==========================================================================
      CALL LIB$DIAGC8(NCHI,RHO,EIG,U)
      DO I=1,NCHI
        FN(I)=EIG(NCHI+1-I)
      ENDDO
      FN(NCHI+1)=0.D0
      FN(0)=1.D0
print*,' fn ',fn(1:nchi)
!
!     ==========================================================================
!     == CONSTRUCT WAVE FUNCTION                                              ==
!     ==========================================================================
      CALL CI$DELETEPSI(PSIN)
      CALL CI$SETPSI(PSIN,1,(1.D0,0.D0))  ! VACUUM STATE
!
      CALL CI$DELETEPSI(PSI)
      CSVAR=CMPLX(SQRT(FN(0)-FN(1)),kind=8)
      CALL CI$SETPSI(PSI,1,CSVAR)
      DO N=1,NCHI
!       ========================================================================
!       == APPLY NEXT CREATOR IN THE BASIS OF EIGENSTATES
!       ========================================================================
        CALL CI$DELETEPSI(PSIADD2)
        CALL CI$DELETEPSI(PSICOPY)
        DO K=1,NCHI
          CALL CI$COPYPSI(PSIN,PSICOPY)  
          CALL CI$CREATOR(PSICOPY,K)
          CALL CI$SCALEPSI(PSICOPY,CONJG(U(K,Nchi+1-n)))
          CALL CI$ADDPSI(PSIADD2,PSICOPY)
        ENDDO
        CALL CI$COPYPSI(PSIADD2,PSIN)       
        CALL CI$DELETEPSI(PSIADD2)
        CALL CI$CLEANPSI(PSIN)
!       ========================================================================
!       ==  ADD CONMTRIBUTION OF N-PARTICLE STATE TO WAVE FUNCTION            ==
!       ========================================================================
        CALL CI$COPYPSI(PSIN,PSICOPY)       
        CSVAR=CMPLX(SQRT(FN(N)-FN(N+1)),kind=8)
        CALL CI$SCALEPSI(PSICOPY,CSVAR)
        CALL CI$ADDPSI(PSI,PSICOPY)
        CALL CI$CLEANPSI(PSI)
      ENDDO       
      CALL CI$DELETEPSI(PSIN)
      CALL CI$DELETEPSI(PSIADD2)
      CALL CI$DELETEPSI(PSICOPY)
!
!     ==========================================================================
!     == test                                                                 ==
!     ==========================================================================
!!$      if(ttest) then
!!$        DO I=1,NCHI
!!$          DO J=1,NCHI
!!$            CALL CI$COPYPSI(PSI,PSICOPY)
!!$            CALL CI$IJOCC(I,J,PSICOPY,Csvar)
!!$            diff=abs(csvar-rho(i,j))
!!$            if(diff.gt.1.d-7) then
!!$              WRITE(*,'(2i3,"dev.:",e10.2," from wv.",2F10.7,"from input",2f10.7)') &
!!$     &                                                   i,j,diff,rho(i,j),csvar
!!$            end if
!!$          ENDDO
!!$        ENDDO
!!$      end if
      RETURN
      END 
!
!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**   OPERATIONS FOR DYNAMICS                                                 **
!**                                                                           **
!*******************************************************************************
!*******************************************************************************
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$PROPAGATE(DT,ANNE,MPSI,PSI0,PSIM,HPSI,PSIP)
!     **************************************************************************
!     **                                                                      **
!     **  PSI(+)=PSI(0)*2/(1+A) + PSI(-)/2+A - HPSI*DT**2/M/(1+A)             **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      REAL(8)           ,INTENT(IN) :: DT
      REAL(8)           ,INTENT(IN) :: ANNE
      REAL(8)           ,INTENT(IN) :: MPSI
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSI0
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSIM
      TYPE(CISTATE_TYPE),INTENT(OUT):: PSIP
      TYPE(CISTATE_TYPE),INTENT(IN) :: HPSI
      TYPE(CISTATE_TYPE)            :: PSI1
      COMPLEX(8)                    :: CSVAR1,CSVAR2,CSVAR3
!     **************************************************************************
      CSVAR1=CMPLX(2.D0/(1.D0+ANNE),0.D0,kind=8)
      CSVAR2=CMPLX(1.D0,0.D0,kind=8)-CSVAR1
      CSVAR3=CMPLX(-DT**2/MPSI/(1.D0+ANNE),0.D0,kind=8)
      CALL CI$COPYPSI(PSI0,PSIP)
      CALL CI$SCALEPSI(PSIP,CSVAR1)
      CALL CI$COPYPSI(PSIM,PSI1)
      CALL CI$SCALEPSI(PSI1,CSVAR2)
      CALL CI$ADDPSI(PSIP,PSI1)
      CALL CI$COPYPSI(HPSI,PSI1)
      CALL CI$SCALEPSI(PSI1,CSVAR3)
      CALL CI$ADDPSI(PSIP,PSI1)
      CALL CI$DELETEPSI(PSI1)
      RETURN
      END SUBROUTINE CI$PROPAGATE
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$EPOT(PSI,HPSI,EPOT)
!     **************************************************************************
!     **                                                                      **
!     **  EPOT=<PSI|HPSI>/<PSI|PSI>                                           **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSI
      TYPE(CISTATE_TYPE),INTENT(IN) :: HPSI
      REAL(8)           ,INTENT(OUT):: EPOT
      COMPLEX(8)                    :: CSVAR
!     **************************************************************************
      CALL CI$SCALARPRODUCT(PSI,HPSI,CSVAR)
      EPOT=REAL(CSVAR)
      CALL CI$SCALARPRODUCT(PSI,PSI,CSVAR)
      EPOT=EPOT/REAL(CSVAR)
      RETURN
      END SUBROUTINE CI$EPOT
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$EKIN(DT,MPSI,PSIP,PSIM,EKIN)
!     **************************************************************************
!     **                                                                      **
!     **  Ekin=(<psi(+)|-<psi(-)|) (|\psi(+)>-psi(-)> *mpsi/(2*dt)**2         **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      REAL(8)           ,INTENT(IN) :: DT    ! TIME STEP
      REAL(8)           ,INTENT(IN) :: mpsi  ! wave function mass
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSIP  ! |psi(+)>
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSIM  ! |psi(-)>
      REAL(8)           ,INTENT(OUT):: EKIN
      TYPE(CISTATE_TYPE)            :: PSIcopy
      COMPLEX(8)        ,parameter  :: Cminus=(-1.D0,0.D0)
      COMPLEX(8)                    :: Csvar
!     **************************************************************************
      CALL CI$COPYPSI(PSIM,PSICOPY)
      CALL CI$SCALEPSI(PSICOPY,CMINUS)
      CALL CI$ADDPSI(PSICOPY,PSIP)
      CALL CI$SCALARPRODUCT(PSICOPY,PSICOPY,CSVAR)
      EKIN=MPSI*REAL(CSVAR)/(2.D0*DT)**2
      CALL CI$DELETEPSI(PSICOPY)
      RETURN
      END SUBROUTINE CI$EKIN
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_LAGRANGE(NCHI,PSI0,PSIP,rho,V,E)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4),        INTENT(IN)   :: NCHI
      TYPE(CISTATE_TYPE),INTENT(IN)   :: PSI0
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSIP
      COMPLEX(8)        ,INTENT(IN)   :: RHO(NCHI,NCHI)
      COMPLEX(8)        ,INTENT(OUT)  :: V(NCHI,NCHI)
      real(8)           ,INTENT(OUT)  :: E
      INTEGER(4)                      :: NC   !#(constraints)
      REAL(8)                         :: MAT(NCHI**2+1,NCHI**2+1)
      REAL(8)                         :: MATINV(NCHI**2+1,NCHI**2+1)
      REAL(8)                         :: VEC(NCHI**2+1)
      REAL(8)                         :: VECSUM(NCHI**2+1)
      INTEGER(4)        ,PARAMETER    :: NITER=100  
      REAL(8)           ,PARAMETER    :: TOL=1.D-5
      integer(4)                      :: iter
      integer(4)                      :: ic,icbar,i,j
      real(8)                         :: svar,svar1
      logical(4)        ,parameter    :: ttest=.true.
      complex(8)                      :: csvar
      logical(4)                      :: convg
!     **************************************************************************
      NC=NCHI**2+1
!
!     ==========================================================================
!     == construct matrix for iteration                                       ==
!     == MAT IS THE REAL (NC*NC) MATRIX WITH ELEMENTS                         ==
!     ==    <PSI0| cONSTR(I)*CONSTR(J) |PSI0>                                 ==
!     ==========================================================================
      if(associated(cimat)) then
print*,'1st request cimat (associated)'
        if(size(cimat).ne.size(mat)) then
          deallocate(cimat)
print*,'deallocating cimat'
        end if
      end if
      if(associated(cimat)) then
print*,'2nd request cimat: associated'
        mat=cimat
      else
print*,'2nd request cimat (not associated)'
        CALL CI_LAGRANGEMAT(NCHI,PSI0,NC,MAT)
        allocate(cimat(nc,nc))
        cimat=mat
      end if

!     ==========================================================================
!     == invert matrix for loop                                               ==
!     ==========================================================================
!     == A LITTLE OF A UNITY MATRIX TO AVOID SINGULAR VALUES
      DO I=1,NC
        MAT(I,I)=MAT(I,I)+1.D-8
      ENDDO
      CALL LIB$INVERTR8(NC,MAT,MATINV)
!
!     ==========================================================================
!     == loop to enforce constraints                                          ==
!     ==========================================================================
      VECSUM=0.D0
      DO ITER=1,NITER
!       == CALCULATE DEVIATION OF |PSI(+)> FROM CONSTRAINT =====================
!       == <PSI(+)|CONSTR(I)|PSI(+)>-VALUE(I) ==================================
        CALL CI_CONSTRAINT(NCHI,PSIP,RHO,NC,VEC)
        CONVG=MAXVAL(ABS(VEC)).LT.TOL
        IF(CONVG) EXIT
!       == DETERMINE CORRECTION FOR LAGRANGE MULTIPLIERS =======================
print*,'marke 4',iter,maxval(abs(vec))
        VEC(:)=-MATMUL(MATINV,VEC)
        VECSUM=VECSUM+VEC    ! ADD CORRECTIONS UP TO OBTAIN TOTAL
!       == ADD CORRECTION TO |PSI(+)> == =======================================
        CALL CI_ADDCONSTRAINT(NCHI,PSIP,PSI0,NC,VEC)
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$STOP('CI_LAGRANGE')
      END IF
!
!     ==========================================================================
!     ==  calculates Lagrange multiplicators                                  ==
!     ==  friction and time step not included                                 ==
!     ==========================================================================
      DO J=1,NCHI
        DO I=1,NCHI
! indices i,j has been exchanged due to christians remark
! must be consistent with similar expression in ci_addconstraint
          IC=i+(j-1)*NCHI
          ICBAR=j+(i-1)*NCHI
          IF(J.GT.I) THEN
            V(I,J)=0.5D0*CMPLX(VECSUM(IC),VECSUM(ICBAR),kind=8)
          ELSE IF(J.EQ.I) THEN
            V(I,J)=CMPLX(VECSUM(IC),0.D0,kind=8)
          ELSE IF(J.LT.I) THEN
            V(I,J)=0.5D0*CMPLX(VECSUM(ICBAR),-VECSUM(IC),kind=8)
          END IF
        ENDDO
      ENDDO
      E=VEC(NC)
      RETURN
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ci_lagrangemat(NCHI,psi0,nc,mat)
!     **************************************************************************
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)  :: NCHI
      TYPE(CISTATE_TYPE),INTENT(IN)  :: PSI0
      INTEGER(4)        ,INTENT(IN)  :: nc
      real(8)           ,INTENT(out) :: mat(nc,nc)
      TYPE(CISTATE_TYPE)             :: PSI1,psi2,psi3,psi4
      complex(8)                     :: cmat4(nchi,nchi,nchi,nchi)
      complex(8)                     :: cmat2(nchi,nchi)
      complex(8)                     :: cmat0
      real(8)                        :: rmat2(nchi,nchi)
      real(8)                        :: rmat0
      integer(4)                     :: iindex(2,nc)
      integer(4)                     :: i,j,k,l,ic1,ic2
      complex(8)                     :: csvar1,csvar2
      complex(8)        ,parameter   :: ci=(0.d0,1.d0)
!     **************************************************************************
      IF(NC.NE.NCHI**2+1) THEN
        CALL ERROR$STOP('CI_LAGRANGEMAT')
      END IF
!
!     ==========================================================================
!     == DETERMINE EXPECTATION VALUES OF OPERATORS                            ==
!     ==========================================================================
      CALL CI$SCALARPRODUCT(PSI0,PSI0,CMAT0)
      CMAT4(:,:,:,:)=(0.D0,0.D0)
!print*,'before matrix elements in lagrangemat'
      DO L=1,NCHI
        CALL CI$COPYPSI(PSI0,PSI1)
        CALL CI$ANNIHILATOR(PSI1,L)
        DO K=1,NCHI
          IF(K.GT.L) CYCLE
          CALL CI$COPYPSI(PSI1,PSI2)
          CALL CI$CREATOR(PSI2,K)
          CALL CI$SCALARPRODUCT(PSI0,PSI2,CMAT2(K,L))
          DO J=1,NCHI
            IF(J.GT.L)  CYCLE 
            IF(J.EQ.L.AND.L.NE.K) cycle   ! matrix elements are zero
            CALL CI$COPYPSI(PSI2,PSI3)
            CALL CI$ANNIHILATOR(PSI3,J)
            DO I=1,NCHI
              IF(I.GT.K) CYCLE
              IF(I.EQ.K.AND.J.NE.K) cycle ! matrix elements are zero
              CALL CI$COPYPSI(PSI3,PSI4)
              CALL CI$CREATOR(PSI4,I)
              CALL CI$SCALARPRODUCT(PSI0,PSI4,CMAT4(I,J,K,L))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!print*,'after matrix elements in lagrangemat'
!
!     ==========================================================================
!     ==  complete matrices                                                   ==
!     ==========================================================================
      do i=1,nchi
        do j=1,i-1
          cmat2(i,j)=conjg(cmat2(j,i))
        enddo
      enddo
      do i=1,nchi
        do j=1,nchi
          do k=1,nchi
            if(i.gt.k)  cycle 
            do l=1,nchi
              if(j.gt.l) cycle
              if(k.gt.l) cycle
              if(i.eq.k.and.j.ne.k) cycle !exclude zeros
              if(j.eq.l.and.j.ne.k) cycle ! exclude zeros
              cmat4(l,k,j,i)=conjg(cmat4(i,j,k,l))
!             == permute creators ===============================================
              csvar1=-cmat4(i,j,k,l)
              if(i.eq.j)csvar1=csvar1+cmat2(k,l)
              if(j.eq.k)csvar1=csvar1+cmat2(i,l)
              CMAT4(k,j,i,l)=csvar1
              CMAT4(l,i,j,k)=conjg(CMAT4(k,j,i,l))
!
!             == permute annihilators============================================
              csvar1=-cmat4(i,j,k,l)
              if(k.eq.l)csvar1=csvar1+cmat2(i,j)
              if(j.eq.k)csvar1=csvar1+cmat2(i,l)
              CMAT4(i,l,k,j)=csvar1
              CMAT4(j,k,l,i)=conjg(CMAT4(i,l,k,j))
!
!             == permute creators and annihilators ==============================
              csvar1=cmat4(i,j,k,l)
              if(j.eq.k)csvar1=csvar1-cmat2(i,l)
              if(i.eq.l)csvar1=csvar1+cmat2(k,j)
              CMAT4(k,l,i,j)=csvar1
              CMAT4(j,i,l,k)=conjg(CMAT4(k,l,i,j))
            enddo
          enddo
        enddo
      enddo
!
!     ==========================================================================
!     ==  Symmetrize                                                       ==
!     ==========================================================================
      rmat0=real(cmat0)
      do i=1,nchi
        rmat2(i,i)=real(cmat2(i,i))
        do j=i+1,nchi
          rmat2(i,j)=0.5d0*real(cmat2(i,j)+cmat2(j,i))
          rmat2(j,i)=-0.5d0*aimag(cmat2(i,j)-cmat2(j,i))
         enddo
      enddo
      do l=1,nchi
        do k=1,nchi
          do i=1,nchi
            do j=i,nchi
              csvar1=cmat4(i,j,k,l)
              csvar2=cmat4(j,i,k,l)
              cmat4(i,j,k,l)=0.5d0*(csvar1+csvar2)
              if(i.eq.j) cycle
              cmat4(j,i,k,l)=0.5d0*ci*(csvar1-csvar2)
            enddo
          enddo
        enddo
      enddo
      do k=1,nchi
        do l=k,nchi
          do j=1,nchi
            do i=1,nchi
              csvar1=cmat4(i,j,k,l)
              csvar2=cmat4(i,j,l,k)
              cmat4(i,j,k,l)=0.5d0*(csvar1+csvar2)
              if(k.eq.l) cycle
              cmat4(i,j,l,k)=0.5d0*ci*(csvar1-csvar2)
            enddo
          enddo
        enddo
      enddo
!
!     ==========================================================================
!     ==  mapping array                                                       ==
!     ==========================================================================
      ic1=0
      do j=1,nchi
        do i=1,nchi
          ic1=ic1+1
          iindex(1,ic1)=i
          iindex(2,ic1)=j
        enddo
      enddo
!
!     ==========================================================================
!     ==  construct linear dependence of conmstraints on lagrange multipliers ==
!     ==========================================================================
      mat(:,:)=0.d0
      mat(nc,nc)=rmat0
      do ic1=1,nc-1
        i=iindex(1,ic1)
        j=iindex(2,ic1)
        mat(ic1,nc)=rmat2(i,j)
        mat(nc,ic1)=rmat2(i,j)
      enddo
      do ic1=1,nc-1
        i=iindex(1,ic1)
        j=iindex(2,ic1)
        do ic2=1,nc-1
          k=iindex(1,ic2)
          l=iindex(2,ic2)
          mat(ic1,ic2)=real(cmat4(i,j,k,l))          
        enddo
      enddo
!
!     ==========================================================================
!     ==  apply factor of 2                                                   ==
!     ==========================================================================
      mat=mat+transpose(mat)
!
!     ==========================================================================
!     ==  wrap up                                                             ==
!     ==========================================================================
      call ci$deletepsi(psi1)
      call ci$deletepsi(psi2)
      call ci$deletepsi(psi3)
      call ci$deletepsi(psi4)
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ci_constraint(NCHI,psip,rho,nc,vec)
!     **************************************************************************
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4),        INTENT(IN)  :: NCHI
      INTEGER(4),        INTENT(IN)  :: nc
      TYPE(CISTATE_TYPE),INTENT(IN)  :: PSIp
      complex(8),        INTENT(IN)  :: rho(nchi,nchi) ! target density matrix
      real(8),           INTENT(out) :: vec(nc)
      TYPE(CISTATE_TYPE)             :: PSI1,psi2
      complex(8)                     :: cmat4(nchi,nchi,nchi,nchi)
      complex(8)                     :: cmat2(nchi,nchi)
      complex(8)                     :: cmat0
      real(8)                        :: rmat4(nchi,nchi,nchi,nchi)
      real(8)                        :: rmat2(nchi,nchi)
      real(8)                        :: rmat0
      integer(4)                     :: i,j,k,l,ic1,ic2
      complex(8)                     :: csvar1,csvar2
!     **************************************************************************
      if(nc.ne.nchi**2+1) then
        call error$stop('...')
      end if
!
!     ==========================================================================
!     == DETERMINE EXPECTATION VALUES OF OPERATORS==============================
!     ==========================================================================
      CALL CI$SCALARPRODUCT(PSIp,PSIp,cmAT0)
      DO I=1,NCHI
        CALL CI$COPYPSI(PSIp,PSI1)
        CALL CI$ANNIHILATOR(PSI1,I)
        DO J=1,NCHI
          CALL CI$COPYPSI(PSI1,PSI2)
          CALL CI$CREATOR(PSI2,J)
          CALL CI$SCALARPRODUCT(PSIp,PSI2,cmAT2(J,I))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  subtract target constraint value                                 ==
!     ==========================================================================
      cmat2(:,:)=cmat2(:,:)-rho(:,:)
      cmat0=cmat0-(1.d0,0.d0)
!
!     ==========================================================================
!     ==  Symmetrize                                                       ==
!     ==========================================================================
      rmat0=real(cmat0)
      do i=1,nchi
        rmat2(i,i)=real(cmat2(i,i))
        do j=i+1,nchi
          rmat2(i,j)=0.5d0*real(cmat2(i,j)+cmat2(j,i))
          rmat2(j,i)=-0.5d0*aimag(cmat2(i,j)-cmat2(j,i))
        enddo
      enddo
!
!     ==========================================================================
!     ==  mapping array                                                       ==
!     ==========================================================================
      ic1=0
      do j=1,nchi
        do i=1,nchi
          ic1=ic1+1
          vec(ic1)=rmat2(i,j)
        enddo
      enddo
      VEC(NC)=RMAT0
!
!     ==========================================================================
!     ==  WRAP UP                                                             ==
!     ==========================================================================
      CALL CI$DELETEPSI(PSI1)
      CALL CI$DELETEPSI(PSI2)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_ADDCONSTRAINT(NCHI,PSIP,PSI0,NC,VEC)
!     **************************************************************************
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4),        INTENT(IN)   :: NCHI
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSIP
      TYPE(CISTATE_TYPE),INTENT(IN)   :: PSI0
      INTEGER(4),        INTENT(IN)   :: NC
      REAL(8),           INTENT(IN)   :: VEC(NC)
      TYPE(CISTATE_TYPE)              :: PSI1
      TYPE(CISTATE_TYPE)              :: PSI2
      REAL(8)                         :: SVAR
      INTEGER(4)                      :: I,J,IC,ICBAR
      COMPLEX(8)                      :: LAMBDA(NCHI,NCHI)
      COMPLEX(8)                      :: CE
!     **************************************************************************
      IF(NC.NE.NCHI**2+1) THEN
        CALL ERROR$STOP('...')
      END IF
!
!     ==========================================================================
!     ==  DESYMMETRIZE                                                       ==
!     ==========================================================================
      DO J=1,NCHI
        DO I=1,NCHI
          IC=I+(J-1)*NCHI
          ICBAR=J+(I-1)*NCHI
          IF(J.GT.I) THEN
            LAMBDA(I,J)=0.5D0*CMPLX(VEC(IC),VEC(ICBAR),kind=8)
          ELSE IF(J.EQ.I) THEN
            LAMBDA(I,J)=CMPLX(VEC(IC),0.D0,kind=8)
          ELSE IF(J.LT.I) THEN
!print*,'ic',icbar,ic
!print*,'vec',vec(icbar),vec(ic)
!print*,'cmplx',CMPLX(VEC(ICBAR),-VEC(IC),kind=8)
            LAMBDA(I,J)=0.5D0*CMPLX(VEC(ICBAR),-VEC(IC),kind=8)
!print*,'over'
          END IF
        ENDDO
      ENDDO
      CE=CMPLX(VEC(NC),0.D0,kind=8)
!
!     ==========================================================================
!     == DETERMINE EXPECTATION VALUES OF OPERATORS==============================
!     ==========================================================================
      CALL CI$COPYPSI(PSI0,PSI1)
      CALL CI$SCALEPSI(PSI1,CE)
      CALL CI$ADDPSI(PSIP,PSI1)
!
      DO J=1,NCHI
        CALL CI$COPYPSI(PSI0,PSI1)
        CALL CI$ANNIHILATOR(PSI1,J)
        DO I=1,NCHI
          CALL CI$COPYPSI(PSI1,PSI2)
          CALL CI$CREATOR(PSI2,I)
          CALL CI$SCALEPSI(PSI2,LAMBDA(I,J))
          CALL CI$ADDPSI(PSIP,PSI2)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  WRAP UP                                                             ==
!     ==========================================================================
      CALL CI$DELETEPSI(PSI1)
      CALL CI$DELETEPSI(PSI2)
      RETURN
      END SUBROUTINE CI_ADDCONSTRAINT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$DYNWITHFIXEDDENMAT(NCHI,RHO2,HAMILTON,PSI0,V,E)
!     **************************************************************************
!     **  DETERMINE MANY PARTICLE WAVE FUNCTION CONSISTENT WITH THE           **
!     **  ONE-PARTICLE DENSITY MATRIX RHO, WHICH MINIMIZES THE TOTAL ENERGY.  **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN) :: NCHI            ! #(1P ORBITALS)
      COMPLEX(8)        ,INTENT(IN) :: RHO2(NCHI,NCHI) ! 1P DENSITY MATRIX
      TYPE(CIHAMIL_TYPE),INTENT(IN) :: HAMILTON        ! MANY-P HAMILTONIAN
      TYPE(CISTATE_TYPE),INTENT(OUT):: PSI0            ! CURRENT WAVE FUNCTION
      COMPLEX(8)        ,INTENT(OUT):: V(NCHI,NCHI)    ! DE/DRHO2
      REAL(8)           ,INTENT(OUT):: E               ! ENERGY
      INTEGER(4)        ,PARAMETER  :: NOITER=10000     ! X#(ITERATIONS)
      REAL(8)           ,PARAMETER  :: DELTA=1.D-1     ! TIME STEP IN LOOP
      REAL(8)                       :: ALPHA=0.D0      ! FRICTION (1.D-1)
      REAL(8)           ,PARAMETER  :: MASS=1.D0       ! MASS
      real(8)           ,parameter  :: maxekin=0.1d0   ! x(for hpsi contribution)
      real(8)           ,parameter  :: maxhpsi=1.d-3   ! x(for hpsi contribution)
      real(8)           ,parameter  :: maxpsi=1.d-4    ! x(for psi contribution)
      real(8)           ,parameter  :: targetekindot=1.d-2    ! x(for psi contribution)
      integer(4)        ,parameter  :: nexpandbasis=2  ! 
      TYPE(CISTATE_TYPE)            :: PSIM            ! PREVIOUS WAVE FUNCTION
      TYPE(CISTATE_TYPE)            :: PSIP            ! NEXT WAVE FUNCTION
      TYPE(CISTATE_TYPE)            :: PSIBAR          ! PROPAGATED W/O CONSTR.
      TYPE(CISTATE_TYPE)            :: HPSI            ! H|PSI>
      TYPE(CISTATE_TYPE)            :: PSIB
      REAL(8)                       :: EPOT,epotlast
      REAL(8)                       :: EKIN
      INTEGER(4)                    :: ITER
      TYPE(CISTATE_TYPE)            :: PSItest         ! CURRENT WAVE FUNCTION
      integer(4)                    :: i,j
      logical(4)                    :: convg
      real(8)          ,parameter   :: tol=1.d-5
      integer(4)       ,parameter   :: nwait=20
      integer(4)                    :: iwait
      real(8)                       :: ewait
      REAL(8)                       :: SVAR      
!     **************************************************************************
!     OPEN(123, FILE="CIINFO", STATUS="OLD")
      OPEN(123, FILE="CIINFO")
      REWIND 123
      nullify(cimat)
!     ===========================================================================
!     == CHECK IF THE DENSITY MATRIX IS HERMITEAN                              ==
!     ===========================================================================
      DO I=1,NCHI
        DO J=I,NCHI
          IF(ABS(RHO2(I,J)-CONJG(RHO2(J,I))).GT.1.D-8) THEN
            CALL ERROR$MSG('DENSITY MATRIX NOT HERMITEAN. NOT ALLOWED...')
            CALL ERROR$STOP('CI$DYNWITHFIXEDDENMAT')
          END IF
        END DO
      ENDDO
!     ===========================================================================
!     == CONSTRUCT STARTING WAVE FUNCTION                                      ==
!     ===========================================================================
      CALL CI$STATEFROMDENSITYMATRIX(NCHI,RHO2,PSI0)
      CALL CI$LIMITSIZE(MAXPSI,PSI0) ! REDUCE THE NUMBER OF SLATER DETERMINANTS
      CALL CI$COPYPSI(PSI0,PSIM) ! ZERO INITIAL VELOCITY 
!
!CALL CI$SETR8('MIN(PSI)',1.D-5)
CALL CI$CLEANPSI(PSI0)
CALL CI$WRITEPSI(PSI0,6)
PRINT*,'PSI%N ',PSI0%N
!
!     ===========================================================================
!     == ITERATE TO FIND GROUND STATE                                          ==
!     ===========================================================================
      EPOT=HUGE(EPOT)
      EKIN=0.D0
      iwait=0
      ewait=epot
      DO ITER=1,NOITER
        CALL CI$CLEANPSI(PSI0)
        CALL CI$COPYPSI(PSI0,HPSI)
        CALL CI$HPSI(HAMILTON,HPSI)  ! CONSTRUCT H|PSI(0)>
        CALL CI$LIMITSIZE(maxhpsi,HPSI) ! REDUCE THE NUMBER OF SLATER DETERMINANTS
!PRINT*,'hPSI%N ',hPSI%N
        EPOTLAST=EPOT
        CALL CI$EPOT(PSI0,HPSI,EPOT)
!
!       == SET FRICTION VALUE ==================================================
!        IF(EPOT.GT.EPOTLAST.AND.EKIN.GT.1.D-6) THEN
        IF(EPOT.GT.EPOTLAST+1.d-6) THEN
          CALL CI$LIMITSIZE(MAXPSI,PSIM) ! REDUCE THE NUMBER OF SLATER DETERMINANTS
          CALL CI$COPYPSI(PSIM,PSI0)
          PRINT*,'PSI0%N AFTER LIMITSIZE',PSI0%N,epot-epotlast
          ekin=1.d-10
          EPOT=EPOTLAST
          cycle
        END IF
        IF(EKIN.GT.MAXEKIN) THEN
!         == keep kinetic energy about constant if ekin>maxekin ===============
          ALPHA=-(EPOT-EPOTLAST)/(4.D0*MAXEKIN)
        ELSE
!         == CHOOSE FRICTION SO THAT EKIN GROWS WITH RATE TARGETEKINDOT
          ALPHA=-(TARGETEKINDOT+EPOT-EPOTLAST)/(4.D0*EKIN)
          ALPHA=MAX(ALPHA,-0.1D0)   ! AVOID TOO LARGE ACCELERATIONS (INSTABLE)
          ALPHA=MIN(ALPHA,0.D0)     ! DO NOT SLOW DOWN 
        END IF
!
!       == PROPAGATE WITHOUT CONSTRAINTS =======================================
        CALL CI$PROPAGATE(DELTA,ALPHA,MASS,PSI0,PSIM,HPSI,PSIBAR)
!
!       == APPLY CONSTRAINTS
        if(iter.lt.20.or.mod(iter,1).eq.0) then
          if(associated(cimat)) deallocate(cimat)
        end if
        CALL CI_LAGRANGE(NCHI,PSI0,PSIBAR,RHO2,V,E)
        CALL CI$COPYPSI(PSIBAR,PSIP)
        CALL CI$DELETEPSI(PSIBAR)
!
!       ========================================================================
!       == LIMIT THE GROWTH OF THE NUMBER OF SLATER DETERMINANTS              ==
!       == THE MAXIMUM ALLOWED VALUE IS NEXPANDBASIS*NCHI                     ==
!       ========================================================================
!PRINT*,'PSI0%N BEFORE GROWTH CHECK',PSIP%N
        CALL CI$GROWTHCHECK(nexpandbasis*NCHI,PSIP,PSI0)
!PRINT*,'PSI0%N AFTER GROWTH CHECK',PSIP%N
!
!       ========================================================================
!       == CHECK CONVERGENCE                                                  ==
!       ========================================================================
        IF(ABS(EPOT-EWAIT).LT.TOL) THEN
          IWAIT=IWAIT+1
        ELSE
          EWAIT=EPOT
          IWAIT=0
        END IF
        CONVG=(IWAIT.GT.NWAIT) 
        IF(CONVG) EXIT
!
!       ========================================================================
!       == KINETIC ENERGY AND ENERGY REPORT                                   ==
!       ========================================================================
        CALL CI$EKIN(DELTA,MASS,PSIP,PSIM,EKIN)
!
        WRITE(*,'(I5,3F30.15,i10)') ITER,EKIN,EPOT,EKIN+EPOT,psi0%n
        WRITE(123,'(I5,3F40.20)') ITER,EPOT,EKIN,EKIN+EPOT
!
!       ========================================================================
!       == SWITCH TO NEXT TIME STEP                                           ==
!       ========================================================================
        CALL CI$COPYPSI(PSI0,PSIM)
        CALL CI$COPYPSI(PSIP,PSI0) 
      END DO
      if(.not.convg) then
        call error$msg('self-consistency loop not converged')
        call error$stop('ci$dynwithfixeddenmat')
      end if
      if(associated(cimat)) deallocate(cimat)
      CALL CI$CLEANPSI(PSI0)
!
!     ===========================================================================
!     == CONSTRUCT CONSTRAINING FORCES                                         ==
!     ===========================================================================
      SVAR=delta**2/mass/(1.d0+alpha)
!      e=e/svar  energy is not the lagrange parameter
      e=epot    ! E=<psi|H|psi>
      v(:,:)=v(:,:)/svar
      close(123)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ci$growthcheck(nmax,psip,psi0)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      integer(4)        ,intent(in)    :: nmax
      TYPE(CISTATE_TYPE),INTENT(in)    :: PSI0         
      TYPE(CISTATE_TYPE),INTENT(inout) :: PSIp         
      integer(4)                       :: i0,ip
      integer(4)                       :: id0,idp
      integer(4)                       :: np
      integer(4)                       :: count
      real(8)           ,allocatable   :: test(:)      
      real(8)                          :: svar
      real(8)                          :: factor
      real(8)                          :: xtest
!     **************************************************************************
      NP=PSIP%N
      ALLOCATE(TEST(NP))
      TEST(:)=-1.d0
      I0=1
      ID0=PSI0%ID(I0)
      DO IP=1,NP
        IDP=PSIP%ID(IP)
        IF(IDP.LT.ID0) THEN  ! NO CORRESPONDING
!          PRINT*,'LONE SLATER DETERMINANT ',IP,PSIP%C(IP) 
          TEST(IP)=ABS(PSIP%C(IP))
        ELSE IF(IDP.GE.ID0) THEN
          DO WHILE (ID0.LE.IDP)
            I0=I0+1
            IF(I0.LE.PSI0%N) THEN
              ID0=PSI0%ID(I0)
            ELSE
              ID0=PSIP%ID(PSIP%N)+1
            END IF
          ENDDO      
        END IF
      ENDDO
!
!     ==========================================================================
      if(psi0%n.ge.nmax) then
        do ip=1,np
          if(test(ip).gt.0.d0)test(ip)=0.d0
        enddo
      else
        xtest=MAXVAL(TEST)
        factor=0.5d0
1000   continue
        SVAR=xtest*(1.d0-factor)
        count=0
        DO IP=1,NP
          IF(TEST(IP).LT.0.d0) cycle
          IF(TEST(IP).LT.SVAR) test(ip)=0.d0
          if(test(ip).gt.0.d0) count=count+1
        enddo
        if(psi0%n+count.gt.nmax) then
          factor=0.5d0*factor
          goto 1000
        end if
      end if
      do ip=1,np
        IF(TEST(IP).EQ.0.D0) PSIP%C(IP)=(0.D0,0.D0)
      ENDDO
      CALL CI_COMPACTPSI(PSIP)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ci$limitsize(xlimit,pHi)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      real(8)           ,intent(in)    :: xlimit
      TYPE(CISTATE_TYPE),INTENT(inOUT) :: PHI         
      integer(4)                       :: i,J
      integer(4)                       :: N
!     **************************************************************************
      N=PHI%N
      J=0
      do i=1,n
        IF(ABS(PHI%C(I)).LT.XLIMIT) CYCLE
        J=J+1
        PHI%C(J)=PHI%C(I)
        PHI%ID(J)=PHI%ID(I)
      ENDDO
      PHI%C(J+1:)=(0.D0,0.D0)
      PHI%ID(J+1:)=0
      PHI%N=J
      RETURN
      END
