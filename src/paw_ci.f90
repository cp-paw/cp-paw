!
!........1.........2.........3.........4.........5.........6.........7.........8
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
!**                                                                           **
!**   remarks:                                                                **
!**     the following points need to be implemented!!!                        **
!**                                                                           **
!**     ci$creator and ci$annihilator do not destroy the order of the array   **
!**     because it is either an addition with a fixed number or a deletion of **
!**     the corresponding element. Is there any operation that destroys the   **
!**     order                                                                 **
!**                                                                           **
!**     avoid deallocate and allocate for ci$copypsi, if not necessary        **
!**                                                                           **
!**     do not overwrite elements wit n>nx with zeros                         **
!**                                                                           **
!*******************************************************************************
MODULE CI_MODULE
TYPE CISTATE_TYPE
  INTEGER(4)             :: NX     ! X (# SLATER DETERMINANTS)
  INTEGER(4)             :: N      ! ACTUAL # SLATER-DETERMINANTS
  logical                :: tclean ! is array ordered?
  COMPLEX(8),POINTER     :: C(:)   ! COEFFICIENTS
  INTEGER   ,POINTER     :: ID(:)  ! NUMBER REPRESENTATION IN BIT FORMAT
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
REAL(8),POINTER           :: CIMAT(:,:)
END MODULE CI_MODULE
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$SETR8(ID,VAL)
!     **************************************************************************
!     **  CI$SETR8                                                            **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
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
      PHI%TCLEAN=.TRUE.
      IF(ASSOCIATED(PHI%C)) THEN
        DEALLOCATE(PHI%C)
        DEALLOCATE(PHI%ID)
      END IF
      RETURN
      END SUBROUTINE CI$DELETEPSI
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$ZEROPSI(PHI)
!     **************************************************************************
!     **  CI$SETPSI                                                           **
!     **  EXPAND ARRAY OF SLATER DETERMINANTS WITH ZERO ENTRIES               **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PHI
!     **************************************************************************
      PHI%NX=0
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
      phi%tclean=.false.
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
        phi%tclean=.true.
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
      PHI%C(N+1:)=(0.d0,0.d0)
      PHI%ID(n+1:)=0
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
      if(phi%tclean) return
!
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
      IF(PHI%N.EQ.0) RETURN
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
      INTEGER(4)                :: N,I,J
      LOGICAL(4)                :: TOCC
      LOGICAL(4)                :: TMINUS
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
      INTEGER(4)                :: N,I,J
      LOGICAL(4)                :: TOCC
      LOGICAL(4)                :: TMINUS
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
      WRITE(NFIL,*)'========================== PSI ============================'
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
      SUBROUTINE CI$DELETEHAMILTONIAN(HAM)
!     **************************************************************************
!     **  CI$DELETEHAMILTONIAN                                                **
!     **  DEALLOCATE ALL ARRAYS                                               **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
!     **************************************************************************
      IF(ASSOCIATED(HAM%U%C)) THEN
        DEALLOCATE(HAM%U%C)
        DEALLOCATE(HAM%U%I)
        DEALLOCATE(HAM%U%J)
        DEALLOCATE(HAM%U%K)
        DEALLOCATE(HAM%U%L)
      END IF
      HAM%U%NX=0
      HAM%U%N=0
      IF(ASSOCIATED(HAM%H%C)) THEN
        DEALLOCATE(HAM%H%C)
        DEALLOCATE(HAM%H%I)
        DEALLOCATE(HAM%H%J)
      END IF
      HAM%H%NX=0
      HAM%H%N=0
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
              FROM=FROM+ISTART-1
              I=U%I(FROM)
              J=U%J(FROM)
              K=U%K(FROM)
              L=U%L(FROM)
              C=U%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              TO=TO+ISTART-1
              U%I(TO)=I
              U%J(TO)=J
              U%K(TO)=K
              U%L(TO)=L
              U%C(TO)=C
            ELSE
              FROM=FROM+ISTART-1
              TO=TO+ISTART-1
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
              FROM=FROM+ISTART-1
              I=U%I(FROM)
              J=U%J(FROM)
              K=U%K(FROM)
              L=U%L(FROM)
              C=U%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              TO=TO+ISTART-1
              U%I(TO)=I
              U%J(TO)=J
              U%K(TO)=K
              U%L(TO)=L
              U%C(TO)=C
            ELSE
              FROM=FROM+ISTART-1
              TO=TO+ISTART-1
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
              FROM=FROM+ISTART-1
              I=U%I(FROM)
              J=U%J(FROM)
              K=U%K(FROM)
              L=U%L(FROM)
              C=U%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              TO=TO+ISTART-1
              U%I(TO)=I
              U%J(TO)=J
              U%K(TO)=K
              U%L(TO)=L
              U%C(TO)=C
            ELSE
              FROM=FROM+ISTART-1
              TO=TO+ISTART-1
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
              FROM=FROM+ISTART-1
              I=H%I(FROM)
              J=H%J(FROM)
              C=H%C(FROM)
            ELSE IF (FROM.EQ.0) THEN
              TO=TO+ISTART-1
              H%I(TO)=I
              H%J(TO)=J
              H%C(TO)=C
            ELSE
              FROM=FROM+ISTART-1
              TO=TO+ISTART-1
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
!     == W=0.5D0*SUM_{I,J,K,L} W_{I,J,K,L} CDAGGER_I CDAGGER_J C_L C_K        ==
!     ==========================================================================
      CALL CI$ZEROPSI(HPSI)
      N=HAM%U%N
      I=0
      J=0
      K=0
      L=0
      DO I1=1,N
        IF(HAM%U%K(I1).NE.K) THEN
          I=0
          J=0
          L=0
          K=HAM%U%K(I1)
          CALL CI$COPYPSI(PSI,PSI1)
          CALL CI$ANNIHILATOR(PSI1,K)    ! |PSI1>=C_K|PSI>
        END IF
        IF(HAM%U%L(I1).NE.L) THEN
          I=0
          J=0
          L=HAM%U%L(I1)
          CALL CI$COPYPSI(PSI1,PSI2)
          CALL CI$ANNIHILATOR(PSI2,L)    ! |PSI2>=C_L C_K |PSI>
        END IF
        IF(HAM%U%J(I1).NE.J) THEN
          I=0
          J=HAM%U%J(I1)
          CALL CI$COPYPSI(PSI2,PSI3)
          CALL CI$CREATOR(PSI3,J)        ! |PSI3>=CDAGGER_J C_L C_K |PSI>
        END IF
        I=HAM%U%I(I1)
        CALL CI$COPYPSI(PSI3,PSI4)
        CALL CI$CREATOR(PSI4,I)       ! |PSI4>=CDAGGER_I CDAGGER_J C_L C_K |PSI>
        CALL CI$SCALEPSI(PSI4,0.5D0*HAM%U%C(I1))  
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
      SUBROUTINE CI$1Pdenmat(iorb1,iorb2,PSI,cval)
!     **************************************************************************
!     **  CIHAMI_ADDTOU                                                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSI
      INTEGER(4)        ,INTENT(IN) :: IORB1
      INTEGER(4)        ,INTENT(IN) :: IORB2
      complex(8)        ,intent(out):: cval
      TYPE(CISTATE_TYPE)            :: PSI1
!     **************************************************************************
      CALL CI$COPYPSI(PSI,PSI1)
      CALL CI$ANNIHILATOR(PSI1,IORB2)
      CALL CI$CREATOR(PSI1,IORB1)
      call ci$scalarproduct(psi,psi1,cval)
      CALL CI$DELETEPSI(PSI1)
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
      REAL(8)                         :: FN(0:NCHI+1)
      REAL(8)                         :: EIG(NCHI)
      TYPE(CISTATE_TYPE)              :: PSIN
      TYPE(CISTATE_TYPE)              :: PSIADD2
      TYPE(CISTATE_TYPE)              :: PSICOPY
      COMPLEX(8)                      :: CSVAR
      REAL(8)                         :: DIFF
      INTEGER(4)                      :: I,j,N,K
      LOGICAL(4)       ,PARAMETER     :: TTEST=.FALSE.
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
!
!     ==========================================================================
!     == CONSTRUCT WAVE FUNCTION                                              ==
!     ==========================================================================
      CALL CI$DELETEPSI(PSIN)
      CALL CI$SETPSI(PSIN,1,(1.D0,0.D0))  ! VACUUM STATE
!
      CALL CI$DELETEPSI(PSI)
      CSVAR=CMPLX(SQRT(FN(0)-FN(1)),KIND=8)
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
          CALL CI$SCALEPSI(PSICOPY,CONJG(U(K,NCHI+1-N)))
          CALL CI$ADDPSI(PSIADD2,PSICOPY)
        ENDDO
        CALL CI$COPYPSI(PSIADD2,PSIN)       
        CALL CI$DELETEPSI(PSIADD2)
        CALL CI$CLEANPSI(PSIN)
!       ========================================================================
!       ==  ADD CONMTRIBUTION OF N-PARTICLE STATE TO WAVE FUNCTION            ==
!       ========================================================================
        CALL CI$COPYPSI(PSIN,PSICOPY)       
        CSVAR=CMPLX(SQRT(FN(N)-FN(N+1)),KIND=8)
        CALL CI$SCALEPSI(PSICOPY,CSVAR)
        CALL CI$ADDPSI(PSI,PSICOPY)
        CALL CI$CLEANPSI(PSI)
      ENDDO       
      CALL CI$DELETEPSI(PSIN)
      CALL CI$DELETEPSI(PSIADD2)
      CALL CI$DELETEPSI(PSICOPY)
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        DO I=1,NCHI
          DO J=1,NCHI
            CALL CI$1PDENMAT(I,J,PSI,CSVAR)
            DIFF=ABS(CSVAR-RHO(I,J))
            IF(DIFF.GT.1.D-7) THEN
              WRITE(*,'(2I3,"DEV.:",E10.2," FROM WV.",2F10.7,"FROM INPUT",2F10.7)') &
     &                                                 I,J,DIFF,RHO(I,J),CSVAR
            END IF
          ENDDO
        ENDDO
      END IF
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
      CSVAR1=CMPLX(2.D0/(1.D0+ANNE),0.D0,KIND=8)
      CSVAR2=CMPLX(1.D0,0.D0,KIND=8)-CSVAR1
      CSVAR3=CMPLX(-DT**2/MPSI/(1.D0+ANNE),0.D0,KIND=8)
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
!     **  EKIN=(<PSI(+)|-<PSI(-)|) (|\PSI(+)>-PSI(-)> *MPSI/(2*DT)**2         **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      REAL(8)           ,INTENT(IN) :: DT    ! TIME STEP
      REAL(8)           ,INTENT(IN) :: MPSI  ! WAVE FUNCTION MASS
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSIP  ! |PSI(+)>
      TYPE(CISTATE_TYPE),INTENT(IN) :: PSIM  ! |PSI(-)>
      REAL(8)           ,INTENT(OUT):: EKIN
      TYPE(CISTATE_TYPE)            :: PSICOPY
      COMPLEX(8)        ,PARAMETER  :: CMINUS=(-1.D0,0.D0)
      COMPLEX(8)                    :: CSVAR
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
      SUBROUTINE CI_LAGRANGE(NCHI,PSI0,PSIP,RHO,V,E)
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
      REAL(8)           ,INTENT(OUT)  :: E
      INTEGER(4)                      :: NC   !#(CONSTRAINTS)
      REAL(8)                         :: MAT(NCHI**2+1,NCHI**2+1)
      REAL(8)                         :: MATINV(NCHI**2+1,NCHI**2+1)
      REAL(8)                         :: VEC(NCHI**2+1)
      REAL(8)                         :: VECSUM(NCHI**2+1)
      INTEGER(4)        ,PARAMETER    :: NITER=100  
      REAL(8)           ,PARAMETER    :: TOL=1.D-5
      INTEGER(4)                      :: ITER
      INTEGER(4)                      :: IC,ICBAR,I,J
      LOGICAL(4)                      :: CONVG
!     **************************************************************************
      NC=NCHI**2+1
!
!     ==========================================================================
!     == CONSTRUCT MATRIX FOR ITERATION                                       ==
!     == MAT IS THE REAL (NC*NC) MATRIX WITH ELEMENTS                         ==
!     ==    <PSI0| CONSTR(I)*CONSTR(J) |PSI0>                                 ==
!     ==========================================================================
      IF(ASSOCIATED(CIMAT)) THEN
        IF(SIZE(CIMAT).NE.SIZE(MAT)) THEN
          DEALLOCATE(CIMAT)
        END IF
      END IF
      IF(ASSOCIATED(CIMAT)) THEN
        MAT=CIMAT
      ELSE
        CALL CI_LAGRANGEMAT(NCHI,PSI0,NC,MAT)
        ALLOCATE(CIMAT(NC,NC))
        CIMAT=MAT
      END IF

!     ==========================================================================
!     == INVERT MATRIX FOR LOOP                                               ==
!     ==========================================================================
!     == A LITTLE OF A UNITY MATRIX TO AVOID SINGULAR VALUES
      DO I=1,NC
        MAT(I,I)=MAT(I,I)+1.D-8
      ENDDO
      CALL LIB$INVERTR8(NC,MAT,MATINV)
!
!     ==========================================================================
!     == LOOP TO ENFORCE CONSTRAINTS                                          ==
!     ==========================================================================
      VECSUM=0.D0
      DO ITER=1,NITER
!       == CALCULATE DEVIATION OF |PSI(+)> FROM CONSTRAINT =====================
!       == <PSI(+)|CONSTR(I)|PSI(+)>-VALUE(I) ==================================
        CALL CI_CONSTRAINT(NCHI,PSIP,RHO,NC,VEC)
        CONVG=MAXVAL(ABS(VEC)).LT.TOL
        IF(CONVG) EXIT
!       == DETERMINE CORRECTION FOR LAGRANGE MULTIPLIERS =======================
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
!     ==  CALCULATES LAGRANGE MULTIPLICATORS                                  ==
!     ==  FRICTION AND TIME STEP NOT INCLUDED                                 ==
!     ==========================================================================
      DO J=1,NCHI
        DO I=1,NCHI
! INDICES I,J HAS BEEN EXCHANGED DUE TO CHRISTIANS REMARK
! MUST BE CONSISTENT WITH SIMILAR EXPRESSION IN CI_ADDCONSTRAINT
          IC=I+(J-1)*NCHI
          ICBAR=J+(I-1)*NCHI
          IF(J.GT.I) THEN
            V(I,J)=0.5D0*CMPLX(VECSUM(IC),VECSUM(ICBAR),KIND=8)
          ELSE IF(J.EQ.I) THEN
            V(I,J)=CMPLX(VECSUM(IC),0.D0,KIND=8)
          ELSE IF(J.LT.I) THEN
            V(I,J)=0.5D0*CMPLX(VECSUM(ICBAR),-VECSUM(IC),KIND=8)
          END IF
        ENDDO
      ENDDO
      E=VEC(NC)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_LAGRANGEMAT(NCHI,PSI0,NC,MAT)
!     **************************************************************************
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)  :: NCHI
      TYPE(CISTATE_TYPE),INTENT(IN)  :: PSI0
      INTEGER(4)        ,INTENT(IN)  :: NC
      REAL(8)           ,INTENT(OUT) :: MAT(NC,NC)
      TYPE(CISTATE_TYPE)             :: PSI1,PSI2,PSI3,PSI4
      COMPLEX(8)                     :: CMAT4(NCHI,NCHI,NCHI,NCHI)
      COMPLEX(8)                     :: CMAT2(NCHI,NCHI)
      COMPLEX(8)                     :: CMAT0
      REAL(8)                        :: RMAT2(NCHI,NCHI)
      REAL(8)                        :: RMAT0
      INTEGER(4)                     :: IINDEX(2,NC)
      INTEGER(4)                     :: I,J,K,L,IC1,IC2
      COMPLEX(8)                     :: CSVAR1,CSVAR2
      COMPLEX(8)        ,PARAMETER   :: CI=(0.D0,1.D0)
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
!PRINT*,'BEFORE MATRIX ELEMENTS IN LAGRANGEMAT'
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
            IF(J.EQ.L.AND.L.NE.K) CYCLE   ! MATRIX ELEMENTS ARE ZERO
            CALL CI$COPYPSI(PSI2,PSI3)
            CALL CI$ANNIHILATOR(PSI3,J)
            DO I=1,NCHI
              IF(I.GT.K) CYCLE
              IF(I.EQ.K.AND.J.NE.K) CYCLE ! MATRIX ELEMENTS ARE ZERO
              CALL CI$COPYPSI(PSI3,PSI4)
              CALL CI$CREATOR(PSI4,I)
              CALL CI$SCALARPRODUCT(PSI0,PSI4,CMAT4(I,J,K,L))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!PRINT*,'AFTER MATRIX ELEMENTS IN LAGRANGEMAT'
!
!     ==========================================================================
!     ==  COMPLETE MATRICES                                                   ==
!     ==========================================================================
      DO I=1,NCHI
        DO J=1,I-1
          CMAT2(I,J)=CONJG(CMAT2(J,I))
        ENDDO
      ENDDO
      DO I=1,NCHI
        DO J=1,NCHI
          DO K=1,NCHI
            IF(I.GT.K)  CYCLE 
            DO L=1,NCHI
              IF(J.GT.L) CYCLE
              IF(K.GT.L) CYCLE
              IF(I.EQ.K.AND.J.NE.K) CYCLE !EXCLUDE ZEROS
              IF(J.EQ.L.AND.J.NE.K) CYCLE ! EXCLUDE ZEROS
              CMAT4(L,K,J,I)=CONJG(CMAT4(I,J,K,L))
!             == PERMUTE CREATORS ==============================================
              CSVAR1=-CMAT4(I,J,K,L)
              IF(I.EQ.J)CSVAR1=CSVAR1+CMAT2(K,L)
              IF(J.EQ.K)CSVAR1=CSVAR1+CMAT2(I,L)
              CMAT4(K,J,I,L)=CSVAR1
              CMAT4(L,I,J,K)=CONJG(CMAT4(K,J,I,L))
!
!             == PERMUTE ANNIHILATORS===========================================
              CSVAR1=-CMAT4(I,J,K,L)
              IF(K.EQ.L)CSVAR1=CSVAR1+CMAT2(I,J)
              IF(J.EQ.K)CSVAR1=CSVAR1+CMAT2(I,L)
              CMAT4(I,L,K,J)=CSVAR1
              CMAT4(J,K,L,I)=CONJG(CMAT4(I,L,K,J))
!
!             == PERMUTE CREATORS AND ANNIHILATORS =============================
              CSVAR1=CMAT4(I,J,K,L)
              IF(J.EQ.K)CSVAR1=CSVAR1-CMAT2(I,L)
              IF(I.EQ.L)CSVAR1=CSVAR1+CMAT2(K,J)
              CMAT4(K,L,I,J)=CSVAR1
              CMAT4(J,I,L,K)=CONJG(CMAT4(K,L,I,J))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  SYMMETRIZE                                                          ==
!     ==========================================================================
      RMAT0=REAL(CMAT0)
      DO I=1,NCHI
        RMAT2(I,I)=REAL(CMAT2(I,I))
        DO J=I+1,NCHI
          RMAT2(I,J)=0.5D0*REAL(CMAT2(I,J)+CMAT2(J,I))
          RMAT2(J,I)=-0.5D0*AIMAG(CMAT2(I,J)-CMAT2(J,I))
         ENDDO
      ENDDO
      DO L=1,NCHI
        DO K=1,NCHI
          DO I=1,NCHI
            DO J=I,NCHI
              CSVAR1=CMAT4(I,J,K,L)
              CSVAR2=CMAT4(J,I,K,L)
              CMAT4(I,J,K,L)=0.5D0*(CSVAR1+CSVAR2)
              IF(I.EQ.J) CYCLE
              CMAT4(J,I,K,L)=0.5D0*CI*(CSVAR1-CSVAR2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO K=1,NCHI
        DO L=K,NCHI
          DO J=1,NCHI
            DO I=1,NCHI
              CSVAR1=CMAT4(I,J,K,L)
              CSVAR2=CMAT4(I,J,L,K)
              CMAT4(I,J,K,L)=0.5D0*(CSVAR1+CSVAR2)
              IF(K.EQ.L) CYCLE
              CMAT4(I,J,L,K)=0.5D0*CI*(CSVAR1-CSVAR2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  MAPPING ARRAY                                                       ==
!     ==========================================================================
      IC1=0
      DO J=1,NCHI
        DO I=1,NCHI
          IC1=IC1+1
          IINDEX(1,IC1)=I
          IINDEX(2,IC1)=J
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  CONSTRUCT LINEAR DEPENDENCE OF CONMSTRAINTS ON LAGRANGE MULTIPLIERS ==
!     ==========================================================================
      MAT(:,:)=0.D0
      MAT(NC,NC)=RMAT0
      DO IC1=1,NC-1
        I=IINDEX(1,IC1)
        J=IINDEX(2,IC1)
        MAT(IC1,NC)=RMAT2(I,J)
        MAT(NC,IC1)=RMAT2(I,J)
      ENDDO
      DO IC1=1,NC-1
        I=IINDEX(1,IC1)
        J=IINDEX(2,IC1)
        DO IC2=1,NC-1
          K=IINDEX(1,IC2)
          L=IINDEX(2,IC2)
          MAT(IC1,IC2)=REAL(CMAT4(I,J,K,L))          
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  APPLY FACTOR OF 2                                                   ==
!     ==========================================================================
      MAT=MAT+TRANSPOSE(MAT)
!
!     ==========================================================================
!     ==  WRAP UP                                                             ==
!     ==========================================================================
      CALL CI$DELETEPSI(PSI1)
      CALL CI$DELETEPSI(PSI2)
      CALL CI$DELETEPSI(PSI3)
      CALL CI$DELETEPSI(PSI4)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_CONSTRAINT(NCHI,PSIP,RHO,NC,VEC)
!     **************************************************************************
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4),        INTENT(IN)  :: NCHI
      INTEGER(4),        INTENT(IN)  :: NC
      TYPE(CISTATE_TYPE),INTENT(IN)  :: PSIP
      COMPLEX(8),        INTENT(IN)  :: RHO(NCHI,NCHI) ! TARGET DENSITY MATRIX
      REAL(8),           INTENT(OUT) :: VEC(NC)
      TYPE(CISTATE_TYPE)             :: PSI1,PSI2
      COMPLEX(8)                     :: CMAT2(NCHI,NCHI)
      COMPLEX(8)                     :: CMAT0
      REAL(8)                        :: RMAT4(NCHI,NCHI,NCHI,NCHI)
      REAL(8)                        :: RMAT2(NCHI,NCHI)
      REAL(8)                        :: RMAT0
      INTEGER(4)                     :: I,J,K,L,IC1,IC2
!     **************************************************************************
      IF(NC.NE.NCHI**2+1) THEN
        CALL ERROR$STOP('...')
      END IF
!
!     ==========================================================================
!     == DETERMINE EXPECTATION VALUES OF OPERATORS==============================
!     ==========================================================================
      CALL CI$SCALARPRODUCT(PSIP,PSIP,CMAT0)
      DO I=1,NCHI
        CALL CI$COPYPSI(PSIP,PSI1)
        CALL CI$ANNIHILATOR(PSI1,I)
        DO J=1,NCHI
          CALL CI$COPYPSI(PSI1,PSI2)
          CALL CI$CREATOR(PSI2,J)
          CALL CI$SCALARPRODUCT(PSIP,PSI2,CMAT2(J,I))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  SUBTRACT TARGET CONSTRAINT VALUE                                 ==
!     ==========================================================================
      CMAT2(:,:)=CMAT2(:,:)-RHO(:,:)
      CMAT0=CMAT0-(1.D0,0.D0)
!
!     ==========================================================================
!     ==  SYMMETRIZE                                                       ==
!     ==========================================================================
      RMAT0=REAL(CMAT0)
      DO I=1,NCHI
        RMAT2(I,I)=REAL(CMAT2(I,I))
        DO J=I+1,NCHI
          RMAT2(I,J)=0.5D0*REAL(CMAT2(I,J)+CMAT2(J,I))
          RMAT2(J,I)=-0.5D0*AIMAG(CMAT2(I,J)-CMAT2(J,I))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  MAPPING ARRAY                                                       ==
!     ==========================================================================
      IC1=0
      DO J=1,NCHI
        DO I=1,NCHI
          IC1=IC1+1
          VEC(IC1)=RMAT2(I,J)
        ENDDO
      ENDDO
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
            LAMBDA(I,J)=0.5D0*CMPLX(VEC(IC),VEC(ICBAR),KIND=8)
          ELSE IF(J.EQ.I) THEN
            LAMBDA(I,J)=CMPLX(VEC(IC),0.D0,KIND=8)
          ELSE IF(J.LT.I) THEN
!PRINT*,'IC',ICBAR,IC
!PRINT*,'VEC',VEC(ICBAR),VEC(IC)
!PRINT*,'CMPLX',CMPLX(VEC(ICBAR),-VEC(IC),KIND=8)
            LAMBDA(I,J)=0.5D0*CMPLX(VEC(ICBAR),-VEC(IC),KIND=8)
!PRINT*,'OVER'
          END IF
        ENDDO
      ENDDO
      CE=CMPLX(VEC(NC),0.D0,KIND=8)
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
      REAL(8)           ,PARAMETER  :: MAXEKIN=0.1D0  ! X(FOR HPSI CONTRIBUTION)
      REAL(8)           ,PARAMETER  :: MAXHPSI=1.D-3  ! X(FOR HPSI CONTRIBUTION)
      REAL(8)           ,PARAMETER  :: MAXPSI=1.D-4   ! X(FOR PSI CONTRIBUTION)
      REAL(8)           ,PARAMETER  :: TARGETEKINDOT=1.D-2 !X(FOR PSI CONTRIB.)
      INTEGER(4)        ,PARAMETER  :: NEXPANDBASIS=2  ! 
      TYPE(CISTATE_TYPE)            :: PSIM            ! PREVIOUS WAVE FUNCTION
      TYPE(CISTATE_TYPE)            :: PSIP            ! NEXT WAVE FUNCTION
      TYPE(CISTATE_TYPE)            :: PSIBAR          ! PROPAGATED W/O CONSTR.
      TYPE(CISTATE_TYPE)            :: HPSI            ! H|PSI>
      REAL(8)                       :: EPOT,EPOTLAST
      REAL(8)                       :: EKIN
      INTEGER(4)                    :: ITER
      INTEGER(4)                    :: I,J
      LOGICAL(4)                    :: CONVG
      REAL(8)          ,PARAMETER   :: TOL=1.D-5
      INTEGER(4)       ,PARAMETER   :: NWAIT=20
      INTEGER(4)                    :: IWAIT
      REAL(8)                       :: EWAIT
      REAL(8)                       :: SVAR      
!     **************************************************************************
!     OPEN(123, FILE="CIINFO", STATUS="OLD")
      OPEN(123, FILE="CIINFO")
      REWIND 123
      NULLIFY(CIMAT)
!     ==========================================================================
!     == CHECK IF THE DENSITY MATRIX IS HERMITEAN                             ==
!     ==========================================================================
      DO I=1,NCHI
        DO J=I,NCHI
          IF(ABS(RHO2(I,J)-CONJG(RHO2(J,I))).GT.1.D-8) THEN
            CALL ERROR$MSG('DENSITY MATRIX NOT HERMITEAN. NOT ALLOWED...')
            CALL ERROR$STOP('CI$DYNWITHFIXEDDENMAT')
          END IF
        END DO
      ENDDO
!     ==========================================================================
!     == CONSTRUCT STARTING WAVE FUNCTION                                     ==
!     ==========================================================================
      CALL CI$STATEFROMDENSITYMATRIX(NCHI,RHO2,PSI0)
      CALL CI$LIMITSIZE(MAXPSI,PSI0) ! REDUCE THE NUMBER OF SLATER DETERMINANTS
      CALL CI$COPYPSI(PSI0,PSIM) ! ZERO INITIAL VELOCITY 
!
!CALL CI$SETR8('MIN(PSI)',1.D-5)
CALL CI$CLEANPSI(PSI0)
CALL CI$WRITEPSI(PSI0,6)
PRINT*,'PSI%N ',PSI0%N
!
!     ==========================================================================
!     == ITERATE TO FIND GROUND STATE                                         ==
!     ==========================================================================
      EPOT=HUGE(EPOT)
      EKIN=0.D0
      IWAIT=0
      EWAIT=EPOT
      DO ITER=1,NOITER
        CALL CI$CLEANPSI(PSI0)
        CALL CI$COPYPSI(PSI0,HPSI)
        CALL CI$HPSI(HAMILTON,HPSI)  ! CONSTRUCT H|PSI(0)>
        CALL CI$LIMITSIZE(MAXHPSI,HPSI) !REDUCE NUMBER OF SLATER DETERMINANTS
!PRINT*,'HPSI%N ',HPSI%N
        EPOTLAST=EPOT
        CALL CI$EPOT(PSI0,HPSI,EPOT)
!
!       == SET FRICTION VALUE ==================================================
!        IF(EPOT.GT.EPOTLAST.AND.EKIN.GT.1.D-6) THEN
        IF(EPOT.GT.EPOTLAST+1.D-6) THEN
          CALL CI$LIMITSIZE(MAXPSI,PSIM) ! REDUCE NUMBER OF SLATER DETERMINANTS
          CALL CI$COPYPSI(PSIM,PSI0)
          PRINT*,'PSI0%N AFTER LIMITSIZE',PSI0%N,EPOT-EPOTLAST
          EKIN=1.D-10
          EPOT=EPOTLAST
          CYCLE
        END IF
        IF(EKIN.GT.MAXEKIN) THEN
!         == KEEP KINETIC ENERGY ABOUT CONSTANT IF EKIN>MAXEKIN ================
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
        IF(ITER.LT.20.OR.MOD(ITER,1).EQ.0) THEN
          IF(ASSOCIATED(CIMAT)) DEALLOCATE(CIMAT)
        END IF
        CALL CI_LAGRANGE(NCHI,PSI0,PSIBAR,RHO2,V,E)
        CALL CI$COPYPSI(PSIBAR,PSIP)
        CALL CI$DELETEPSI(PSIBAR)
!
!       ========================================================================
!       == LIMIT THE GROWTH OF THE NUMBER OF SLATER DETERMINANTS              ==
!       == THE MAXIMUM ALLOWED VALUE IS NEXPANDBASIS*NCHI                     ==
!       ========================================================================
!PRINT*,'PSI0%N BEFORE GROWTH CHECK',PSIP%N
        CALL CI$GROWTHCHECK(NEXPANDBASIS*NCHI,PSIP,PSI0)
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
        WRITE(*,'(I5,3F30.15,I10)') ITER,EKIN,EPOT,EKIN+EPOT,PSI0%N
        WRITE(123,'(I5,3F40.20)') ITER,EPOT,EKIN,EKIN+EPOT
!
!       ========================================================================
!       == SWITCH TO NEXT TIME STEP                                           ==
!       ========================================================================
        CALL CI$COPYPSI(PSI0,PSIM)
        CALL CI$COPYPSI(PSIP,PSI0) 
      END DO
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('SELF-CONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$STOP('CI$DYNWITHFIXEDDENMAT')
      END IF
      IF(ASSOCIATED(CIMAT)) DEALLOCATE(CIMAT)
      CALL CI$CLEANPSI(PSI0)
!
!     ==========================================================================
!     == CONSTRUCT CONSTRAINING FORCES                                        ==
!     ==========================================================================
      SVAR=DELTA**2/MASS/(1.D0+ALPHA)
!      E=E/SVAR  ENERGY IS NOT THE LAGRANGE PARAMETER
      E=EPOT    ! E=<PSI|H|PSI>
      V(:,:)=V(:,:)/SVAR
      CLOSE(123)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$GROWTHCHECK(NMAX,PSIP,PSI0)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)    :: NMAX
      TYPE(CISTATE_TYPE),INTENT(IN)    :: PSI0         
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PSIP         
      INTEGER(4)                       :: I0,IP
      INTEGER(4)                       :: ID0,IDP
      INTEGER(4)                       :: NP
      INTEGER(4)                       :: COUNT
      REAL(8)           ,ALLOCATABLE   :: TEST(:)      
      REAL(8)                          :: SVAR
      REAL(8)                          :: FACTOR
      REAL(8)                          :: XTEST
!     **************************************************************************
      NP=PSIP%N
      ALLOCATE(TEST(NP))
      TEST(:)=-1.D0
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
      IF(PSI0%N.GE.NMAX) THEN
        DO IP=1,NP
          IF(TEST(IP).GT.0.D0)TEST(IP)=0.D0
        ENDDO
      ELSE
        XTEST=MAXVAL(TEST)
        FACTOR=0.5D0
1000   CONTINUE
        SVAR=XTEST*(1.D0-FACTOR)
        COUNT=0
        DO IP=1,NP
          IF(TEST(IP).LT.0.D0) CYCLE
          IF(TEST(IP).LT.SVAR) TEST(IP)=0.D0
          IF(TEST(IP).GT.0.D0) COUNT=COUNT+1
        ENDDO
        IF(PSI0%N+COUNT.GT.NMAX) THEN
          FACTOR=0.5D0*FACTOR
          GOTO 1000
        END IF
      END IF
      DO IP=1,NP
        IF(TEST(IP).EQ.0.D0) PSIP%C(IP)=(0.D0,0.D0)
      ENDDO
      CALL CI_COMPACTPSI(PSIP)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$LIMITSIZE(XLIMIT,PHI)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      REAL(8)           ,INTENT(IN)    :: XLIMIT
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PHI         
      INTEGER(4)                       :: I,J
      INTEGER(4)                       :: N
!     **************************************************************************
      N=PHI%N
      J=0
      DO I=1,N
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
