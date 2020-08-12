! TODO: 
!      - FORTRAN 2008 HAS A NUMBER OF NEW INTRINSIC BIT FUNCTIONS SUCH AS 
!        POPCNT AND POPPAR, WHICH IMPLEMENT SOME OF THE FUNCTIONALITY DIRECTLY
!
!      - USE INTEGER(8) FOR ID TO ALLOW LARGER 1-PARTICLE BASISSETS
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
!**  HAMILTONIAN                                                              **
!**    H = SUM A_I CDAGGER_I + SUM CONJG(A_I) C_I                             **
!**      + SUM H_IJ CDAGGER_I C_J                                             **
!**      + SUM U_IJKL CDAGGER_I CDAGGER_J C_K C_L                             **
!**                                                                           **
!**    H = SUM_N CDAGGER(I_N) * A%C_N + C(I_N) CONJG(A%C_N)                   **
!**      + SUM_N CDAGGER(I_N) C(J_N) * H%C_N                                  **
!**      + SUM_N CDAGGER(I_N) CDAGGER(J_N) C(K_N) C(L_N) * U%C_N              **
!**                                                                           **
!**  ONE-PARTICLE REDUCED DENSITY MATRIX                                      **
!**    RHO(I,J)=<PSI| CDAGGER_J C_I |PSI>                                     **
!**                                                                           **
!**  ONE-PARTICLE ENERGY: E_1 = SUM H(I,J)*RHO(J,I)                           **
!**                                                                           **
!**   REMARKS:                                                                **
!**     1) BEFORE USING A CI-STATE IT MUST BE INITIALIZED BY CI$NEWPSI.       **
!**        SAME FOR THE HAMILTONIAN USING CI$NEWHAMILTONIAN                   **
!**     2) PSI%NX=0 IMPLIES DEALLOCATED ARRAYS. ONLY CI$DELETEPSI AND         **
!**        CI$NEWPSI CAN PRODUCE A STATE WITH DEALLOCATED ARRAYS. ONLY        **
!**        CI$NEWPSI AND CI_EXPANDPSI CAN LEAD FROM A STATE WITH              **
!**        DEALLOCATED ARRAYS TO ONE WITH ALLOCATED ARRAYS.                   **
!**     3) A STATE WITH DEALLOCATED ARRAYS DOES NOT EVEN CONTAIN THE          **
!**        VACUUM STATE                                                       **
!**     4) STATES INTERNAL TO A SUBROUTINE MUST BE DELETED WITH CI$DELETEPSI  **
!**        TO AVOID MEMORY LEAKS                                              **
!**                                                                           **
!**     THE FOLLOWING POINTS NEED TO BE IMPLEMENTED!!!                        **
!**     DO NOT OVERWRITE ELEMENTS WIT N>NX WITH ZEROS                         **
!**                                                                           **
!**     MAKE INITIAL TRANSFORM ONTO ORTHONORMAL NATURAL ORBITALS              **
!**     AND FINAL BACK TRANSFORM                                              **
!**                                                                           **
!**     CONSTRAINT ROUTINE CAN ONLY WORK IF THE WAVE FUNCTION CAN HAVE MORE   **
!**     COEFFICIENTS THAN CONSTRAINTS                                         **
!**                                                                           **
!**     ARE ALL OPERATIONS COMPATIBLE WITH A STATE WITH DEALLOCATED ARRAYS?   **
!**                                                                           **
!**     W=0.5D0*SUM_{I,J,K,L} W_{I,J,K,L} CDAGGER_I CDAGGER_J C_L C_K         **
!**     THIS CONVENTION IS CONSISTENT                                         **
!**     --WITH SZABO OSTLUND  W_{I,J,K,L}=<IJ|KL>                             **
!**     --FRANZ WEGENER SKRIPT "THEORETISCHE FESTKOPERPHYSIK"                 **
!**                                                                           **
!*******************************************************************************
MODULE CI_MODULE
TYPE CISTATE_TYPE
  INTEGER(4)             :: NX=0   ! X (# SLATER DETERMINANTS)
  INTEGER(4)             :: N=0    ! ACTUAL # SLATER-DETERMINANTS
  LOGICAL                :: TCLEAN ! IS ARRAY ORDERED?
  COMPLEX(8),POINTER     :: C(:)   ! COEFFICIENTS
  INTEGER   ,POINTER     :: ID(:)  ! NUMBER REPRESENTATION IN BIT FORMAT
END TYPE CISTATE_TYPE
!== INTERACTION PART OF THE HAMILTONIAN A^+_I A^+_J A_K A_L 
TYPE U_TYPE
  INTEGER(4)             :: NX=0  !
  INTEGER(4)             :: N=0   !
  COMPLEX(8),POINTER     :: C(:) !
  INTEGER   ,POINTER     :: I(:)!
  INTEGER   ,POINTER     :: J(:)!
  INTEGER   ,POINTER     :: K(:)!
  INTEGER   ,POINTER     :: L(:)!
END TYPE U_TYPE
!== ONE-PARTICLE PART OF THE HAMILTONIAN
TYPE H_TYPE
  INTEGER(4)             :: NX=0   !
  INTEGER(4)             :: N=0    !
  COMPLEX(8),POINTER     :: C(:) !
  INTEGER   ,POINTER     :: I(:) !
  INTEGER   ,POINTER     :: J(:) !
END TYPE H_TYPE
!== SOURCE TERM OF THE HAMILTONIAN ====================================
TYPE A_TYPE
  INTEGER(4)             :: NX=0   !
  INTEGER(4)             :: N=0    !
  COMPLEX(8),POINTER     :: C(:) !
  INTEGER   ,POINTER     :: I(:) !
END TYPE A_TYPE
!== HAMILTONIAN INCLUDING ONE-PARTICLE TERM, INTERACTION AND SOURCES ==
TYPE CIHAMIL_TYPE
  TYPE(U_TYPE)            :: U  ! INTERACTION 
  TYPE(H_TYPE)            :: H  ! ONE-PARTICLE TERM
  TYPE(A_TYPE)            :: A  ! SOURCE AND SINK TERMS
END TYPE CIHAMIL_TYPE
! MINC IS USED BY CI_COMPACTPSI TO REMOVE SLATER-DETERMINANTS WITH TOO 
! SMALL WEIGHT. 
REAL(8),SAVE              :: CI_MINC=1.D-10  ! MINIMUM ACCEPTABLE COEFFICIENT
END MODULE CI_MODULE
!
!     ..........................................................................
      SUBROUTINE CI_ODDPARITY(IVAL,POS,TMINUS)
!     **************************************************************************
!     ** COUNTS THE NUMBER OF BITS FROM POSITION ONE TO POSITION POS          **
!     ** AND RETURNS TRUE OF THAT NUMBER IS ODD.                              **
!     ** THIS IS THE NUMBER OF PERMUTATIONS OF AN ANNIHILATOR, UNTIL IT IS    **
!     ** PLACED IN FRONT OF THE CREATOR AT POSITION POS+1.
!     **************************************************************************
      IMPLICIT NONE
      INTEGER   ,INTENT(IN) :: IVAL
      INTEGER(4),INTENT(IN) :: POS     !BITS ARE COUNTED FROM 1 TO POS
      LOGICAL(4),INTENT(OUT):: TMINUS
      INTEGER   ,PARAMETER  :: ONE=1
      INTEGER               :: JVAL
      INTEGER               :: MASK
      INTEGER(4),PARAMETER  :: NSHIFTS=8
      INTEGER(4),PARAMETER  :: ISHIFTS(NSHIFTS)=(/128,64,32,16,8,4,2,1/)
      INTEGER(4)            :: I
!     **************************************************************************
      JVAL=IVAL
      MASK=IBSET(0,POS)-1   ! =2**POS-1
      JVAL=IAND(JVAL,MASK)  ! SET ALL BITS BEYOND POS TO ZERO
      DO I=1,NSHIFTS
        IF(POS.LT.ISHIFTS(I)) CYCLE   !THIS STATEMENT IS NOT NECESSARY AND VERY! COSTLY
       JVAL=IEOR(JVAL,ISHFT(JVAL,-ISHIFTS(I)))
      ENDDO
!!$!     THIS IS AN ALTERNATIVE THAT IS FASTER BUT REQUIRES 8-BIT INTEGER
!!$      IF(REALLYFAST) THEN
!!$        JVAL=IAND(JVAL,MASK)  ! SET ALL BITS BEYOND POS TO ZERO
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-128))
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-64))
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-32))
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-16))
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-8))
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-4))
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-2))
!!$        JVAL=IEOR(JVAL,ISHFT(JVAL,-1))
!!$      END IF
      TMINUS=IAND(JVAL,ONE).EQ.1
      RETURN
      END
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
!       == ONLY CI_EXPANDPSI CAN LEAD OUT OF THE NULLIFIED STATE ===============
        DEALLOCATE(PHI%C)
        DEALLOCATE(PHI%ID)
      END IF
      RETURN
      END SUBROUTINE CI$DELETEPSI
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$ZEROPSI(PHI)
!     **************************************************************************
!     **  CI$ZEROPSI                                                          **
!     **  SET WAVE FUNCTION TO ZERO                                           **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PHI
!     **************************************************************************
      PHI%N=0
      PHI%TCLEAN=.TRUE.
      RETURN
      END SUBROUTINE CI$ZEROPSI
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$NEWPSI(PHI,NX)
!     **************************************************************************
!     **  CI$NEWPSI                                                           **
!     **  INITIALIZES A WAVE FUNCTION THAT HAS NOT BEEN USED BEFORE           **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PHI
      INTEGER(4)        ,INTENT(IN)    :: NX
!     **************************************************************************
      PHI%N=0
      PHI%NX=NX
      PHI%TCLEAN=.TRUE.
      IF(NX.GT.0) THEN
        ALLOCATE(PHI%ID(NX))
        ALLOCATE(PHI%C(NX))
      ELSE
!       == ONLY CI_EXPANDPSI CAN LEAD OUT OF THE NULLIFIED STATE ===============
        NULLIFY(PHI%ID)
        NULLIFY(PHI%C)
      END IF
      RETURN
      END SUBROUTINE CI$NEWPSI
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$SETPSI(PHI,ID,C)
!     **************************************************************************
!     **  ADD A SLATER DETERMINANT TO A STATE.                                **
!     **  IF THIS SLATER DETERMINANT IS ALREADY PRESENT, THE VALUE IS ADDED   **
!     **  TO THE EXISTING COEFFICIENT                                         **
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
      PHI%TCLEAN=.FALSE.
      RETURN
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
      COMPLEX(8)         ,POINTER      :: C(:)
      INTEGER            ,POINTER      :: ID(:)
      INTEGER(4)                       :: N
      INTEGER(4)                       :: NPLUS
!     **************************************************************************
      IF(PHI%NX.EQ.0) THEN
        PHI%NX=NFURTHER
        PHI%N=0
        PHI%TCLEAN=.TRUE.
        ALLOCATE(PHI%ID(NFURTHER))
        ALLOCATE(PHI%C(NFURTHER))
        RETURN
      END IF
      NPLUS=PHI%NX+NFURTHER
      N=PHI%N
      ALLOCATE(C(NPLUS))
      ALLOCATE(ID(NPLUS))
      C(:N)=PHI%C(:N)
      ID(:N)=PHI%ID(:N)
      DEALLOCATE(PHI%C)
      DEALLOCATE(PHI%ID)
      PHI%NX=NPLUS
      PHI%ID=>ID 
      PHI%C=>C 
      NULLIFY(ID)
      NULLIFY(C)
      RETURN
      END SUBROUTINE CI_EXPANDPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$CLEANPSI(PHI)
!     **************************************************************************
!     **  CI$CLEANPSI                                                         **
!     **  SORT THE SLATER DETERMINANTS ACCORDING TO ID                        **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE) ,INTENT(INOUT) :: PHI
      INTEGER(4)                        :: N
      INTEGER                           :: ID
      COMPLEX(8)                        :: C
      INTEGER(4)                        :: I,FROM,TO
      REAL(8)           ,ALLOCATABLE    :: CRIT(:)
!     **************************************************************************
      IF(PHI%TCLEAN) RETURN
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
      PHI%TCLEAN=.TRUE.
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
      IF(N.LE.0) THEN
        CALL ERROR$MSG('STATE IS THE ZERO STATE')
        CALL ERROR$MSG('THE ZERO STATE CANNOT BE NORMALIZED')
        CALL ERROR$STOP('CI$NORMALIZE')
      END IF
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
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PHI1
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PHI2
      COMPLEX(8)        ,INTENT(OUT)  :: VAL
      INTEGER(4)                      :: N1,N2
      INTEGER(4)                      :: I1,I2
!     **************************************************************************
      CALL CI$CLEANPSI(PHI1)
      CALL CI$CLEANPSI(PHI2)
      N1=PHI1%N
      N2=PHI2%N
      VAL=(0.D0,0.D0)
      IF(N1.LE.0.OR.N2.LE.0) THEN
        VAL=(0.D0,0.D0)
        RETURN
      END IF
      I1=1
      I2=1
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
          IF(PHI1%ID(I1).LT.PHI2%ID(I2)) THEN
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
!     == EXPAND ARRAY IF NECESSARY =============================================
      IF(NX.GT.PHI1%NX) THEN
        CALL CI_EXPANDPSI(PHI1,NX-PHI1%NX)
      END IF
!     == COPY ARRAY TO THE END SO THAT THE ADDITION CAN BE DONE IN PLACE =======
      PHI1%ID(NX-N1+1:NX)=PHI1%ID(1:N1)
      PHI1%C(NX-N1+1:NX)=PHI1%C(1:N1)
      PHI1%N=NX
!
!     ==========================================================================
!     == DETERMINE RESULT                                                     ==
!     ==========================================================================
      I1=NX-N1+1
      I2=1
      N=0
      DO 
        N=N+1
        IF(PHI2%ID(I2).EQ.PHI1%ID(I1)) THEN
          PHI1%ID(N)=PHI1%ID(I1)
          PHI1%C(N)=PHI1%C(I1)+PHI2%C(I2)
          I1=I1+1
          I2=I2+1
        ELSE
          IF(PHI1%ID(I1).LT.PHI2%ID(I2)) THEN
            PHI1%ID(N)=PHI1%ID(I1)
            PHI1%C(N)=PHI1%C(I1)
            I1=I1+1
          ELSE
            PHI1%ID(N)=PHI2%ID(I2)
            PHI1%C(N)=PHI2%C(I2)
            I2=I2+1
          END IF
        END IF
        IF(I1.GT.NX) THEN
          PHI1%ID(N+1:NX)=PHI2%ID(I2:N2)
          PHI1%C(N+1:NX)=PHI2%C(I2:N2)
          EXIT
        ELSE IF(I2.GT.N2) THEN
          EXIT
        END IF
      ENDDO
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
      IF(N.EQ.0) RETURN   !TAKE CARE OF STATE WITH DEALLOCATED ARRAYS
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
      INTEGER(4)                :: N,I,J,K
      LOGICAL(4)                :: TOCC
      LOGICAL(4)                :: TMINUS
      LOGICAL(4)  ,PARAMETER    :: TNEWPARITY=.TRUE.
!     **************************************************************************
      IF(IORB.GE.BIT_SIZE(PHI%ID(1))) THEN
        STOP 'IORB OUT OF RANGE'
      END IF
      CALL CI$CLEANPSI(PHI)
      N=PHI%N
      J=0
      DO I=1,N
        TOCC=BTEST(PHI%ID(I)-1,IORB-1) 
        IF(.NOT.TOCC) THEN
          J=J+1
!         == CREATE A PARTICLE AT POSITION IORB ================================
          PHI%ID(J)=1+IBSET(PHI%ID(I)-1,IORB-1)
!
!         == DETERMINE SIGN CHANGE BY THE NUMBER OF PERMUTATIONS ===============
          IF(TNEWPARITY) THEN
            CALL CI_ODDPARITY(PHI%ID(J)-1,IORB-1,TMINUS)
          ELSE
            TMINUS=.FALSE.
            DO K=1,IORB-1
              IF(BTEST(PHI%ID(J)-1,K-1))TMINUS=.NOT.TMINUS
            ENDDO
          END IF
          IF(TMINUS) THEN
            PHI%C(J)=-PHI%C(I)
          ELSE
            PHI%C(J)=+PHI%C(I)
          END IF
        END IF
      ENDDO
      PHI%N=J
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
      INTEGER(4)                :: N,I,J,K
      LOGICAL(4)                :: TOCC
      LOGICAL(4)                :: TMINUS
      LOGICAL(4)  ,PARAMETER    :: TNEWPARITY=.TRUE.
!     **************************************************************************
      IF(IORB.GE.BIT_SIZE(PHI%ID(1))) THEN
        STOP 'IORB OUT OF RANGE'
      END IF
      CALL CI$CLEANPSI(PHI)
      N=PHI%N
      J=0
      DO I=1,N
        TOCC= BTEST(PHI%ID(I)-1,IORB-1) 
        IF(TOCC) THEN
          J=J+1
          PHI%ID(J)=1+IBCLR(PHI%ID(I)-1,IORB-1)
!
!         == DETERMINE SIGN CHANGE BY THE NUMBER OF PERMUTATIONS ===============
          IF(TNEWPARITY) THEN
             CALL CI_ODDPARITY(PHI%ID(J)-1,IORB-1,TMINUS)
          ELSE
            TMINUS=.FALSE.
            DO K=1,IORB-1
              IF(BTEST(PHI%ID(J)-1,K-1)) TMINUS=.NOT.TMINUS
            ENDDO
          END IF
!
          IF(TMINUS)THEN
            PHI%C(J)=-PHI%C(I)
          ELSE
            PHI%C(J)=+PHI%C(I)
          END IF
        END IF
      ENDDO
      PHI%N=J
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
      TYPE(CISTATE_TYPE) ,INTENT(IN)   :: PHI1
      TYPE(CISTATE_TYPE) ,INTENT(INOUT):: PHI2
      INTEGER(4)                       :: N
!     **************************************************************************
      N=PHI1%N
      IF(N.EQ.0) THEN ! TAKE CARE OF STATE WITH DEALLOCATED ARRAYS ==========
        CALL CI$ZEROPSI(PHI2)
        RETURN  
      END IF
      IF(N.GT.PHI2%NX) THEN
        CALL CI_EXPANDPSI(PHI2,N-PHI2%NX)
      END IF
      PHI2%N=N
      PHI2%ID(1:N)=PHI1%ID(1:N)
      PHI2%C(1:N) =PHI1%C(1:N)
      PHI2%TCLEAN =PHI1%TCLEAN
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
      IF(N.EQ.0) RETURN ! TAKE CARE OF ZERO STATES
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
      SUBROUTINE CI$NEWHAMILTONIAN(HAM)
!     **************************************************************************
!     **  NULLIFIES A HAMILTONIAN BEFORE FIRST USE                            **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
!     **************************************************************************
      HAM%A%N=0
      HAM%A%NX=0
      NULLIFY(HAM%A%I)
      NULLIFY(HAM%A%C)
      HAM%H%N=0
      HAM%H%NX=0
      NULLIFY(HAM%H%I)
      NULLIFY(HAM%H%J)
      NULLIFY(HAM%H%C)
      HAM%U%N=0
      HAM%U%NX=0
      NULLIFY(HAM%U%I)
      NULLIFY(HAM%U%J)
      NULLIFY(HAM%U%K)
      NULLIFY(HAM%U%L)
      NULLIFY(HAM%U%C)
      RETURN 
      END SUBROUTINE CI$NEWHAMILTONIAN
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
      SUBROUTINE CI$SETA(HAM,I,VAL)
!     **************************************************************************
!     **  CIHAMI$ADDTOU                                                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM
      INTEGER(4)  ,INTENT(IN)   :: I
      COMPLEX(8)  ,INTENT(IN)   :: VAL
!     **************************************************************************
      CALL CI_SETAELEMENT(HAM%A,I,VAL)
      RETURN 
      END SUBROUTINE CI$SETA
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
      CALL CI_WRITEA(HAM%A,NFIL)
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
!CALL TIMING$CLOCKON('CI$CLEANH')
      CALL CI_CLEANU(HAM%U)
      CALL CI_CLEANH(HAM%H)
      CALL CI_CLEANA(HAM%A)
!CALL TIMING$CLOCKOFF('CI$CLEANH')
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
!
      IF(ASSOCIATED(HAM%H%C)) THEN
        DEALLOCATE(HAM%H%C)
        DEALLOCATE(HAM%H%I)
        DEALLOCATE(HAM%H%J)
      END IF
      HAM%H%NX=0
      HAM%H%N=0
!
      IF(ASSOCIATED(HAM%A%C)) THEN
        DEALLOCATE(HAM%A%C)
        DEALLOCATE(HAM%A%I)
      END IF
      HAM%A%NX=0
      HAM%A%N=0
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
      SUBROUTINE CI_WRITEA(A,NFIL)
!     **************************************************************************
!     **  WRITE  TO FILE                              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(A_TYPE),INTENT(INOUT):: A
      INTEGER(4)  ,INTENT(IN)   :: NFIL
      INTEGER(4)                :: N,I
!     **************************************************************************
      N=A%N
      WRITE(NFIL,FMT='("============== A  =======================")')
      DO I=1,N
        WRITE(NFIL,FMT='(I5,2F17.10)')A%I(I),A%C(I)
      ENDDO
      RETURN
      END SUBROUTINE CI_WRITEA
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
      SUBROUTINE CI_SETAELEMENT(A,I,VAL)
!     **************************************************************************
!     **  SET ONE MATRIX ELEMENT OF ONE-PARTICLE HAMILTONIAN                  **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(A_TYPE),INTENT(INOUT):: A
      INTEGER(4)  ,INTENT(IN)   :: I
      COMPLEX(8)  ,INTENT(IN)   :: VAL
      INTEGER(4)                :: N
!     **************************************************************************
      IF(A%N.GE.A%NX) CALL CI_CLEANA(A)
      IF(A%N.GE.A%NX) CALL CI_EXPANDA(A,20)
      A%N=A%N+1
      N=A%N
      A%I(N)=I
      A%C(N)=VAL
      RETURN
      END SUBROUTINE CI_SETAELEMENT
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
      SUBROUTINE CI_EXPANDA(A,NFURTHER)
!     **************************************************************************
!     **  ADD NFURTHER ELEMENTS FOR THE ONE-PARTICLE HAMILTONIAN              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(A_TYPE) ,INTENT(INOUT):: A
      INTEGER(4)   ,INTENT(IN)   :: NFURTHER
      COMPLEX(8)   ,ALLOCATABLE  :: C(:)
      INTEGER(4)   ,ALLOCATABLE  :: IND(:)
      INTEGER(4)                 :: N,NX
!     **************************************************************************
      IF(A%NX.EQ.0) THEN
        A%NX=NFURTHER
        A%N=0
        NX=NFURTHER
        ALLOCATE(A%I(NX))
        ALLOCATE(A%C(NX))
        RETURN
      END IF
      N=A%N
      ALLOCATE(IND(N))
      ALLOCATE(C(N))
      IND(:)=A%I(:N)
      C(:)=A%C(:N)
      DEALLOCATE(A%I)
      DEALLOCATE(A%C)
      A%NX=A%NX+NFURTHER
      NX=A%NX
      ALLOCATE(A%I(NX))
      ALLOCATE(A%C(NX))
      A%I(:N)=IND(:)
      A%C(:N)=C(:)
      DEALLOCATE(IND)
      DEALLOCATE(C)
      RETURN
      END SUBROUTINE CI_EXPANDA
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
      SUBROUTINE CI_CLEANA(A)
!     **************************************************************************
!     **  CI_CLEANA                                                           **
!     **  REMOVES ZERO ELEMENTS OF THE CREATION ANNIHILATION PART             **
!     **  ORDER THE ELEMENTS AND ADD IDENTICAL ELEMENTS                       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(A_TYPE) ,INTENT(INOUT):: A
      INTEGER(4)                 :: N
      INTEGER(4)                 :: I,J
      COMPLEX(8)                 :: C
      REAL(8)  ,ALLOCATABLE      :: CRIT(:)
      INTEGER(4)                 :: FROM,TO
      LOGICAL(4)                 :: TDIFF
!     **************************************************************************
      IF(A%N.EQ.0) RETURN
      N=A%N
      C=(0.D0,0.D0)
!
!     ==========================================================================
!     ==  REMOVE ZERO ELEMENTS                                                ==
!     ==========================================================================
      J=0
      DO I=1,N
        IF(ABS(A%C(I)).EQ.0.D0) CYCLE
        J=J+1
        A%C(J)=A%C(I)
        A%I(J)=A%I(I)
      END DO
      A%N=J
      N=A%N
!
!     ==========================================================================
!     ==  ORDER ENTRIES                                                       ==
!     ==========================================================================
!     ==  ORDER WITH RESPECT TO LAST INDEX
      ALLOCATE(CRIT(N))
      DO I=1,N
        CRIT(I)=REAL(A%I(I))
      ENDDO
      CALL SORT__SET(N,CRIT)
      CALL SORT__RESTART
      CALL SORT__FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          I=A%I(FROM)
          C=A%C(FROM)
        ELSE IF (FROM.EQ.0) THEN
          A%I(TO)=I
          A%C(TO)=C
        ELSE
          A%I(TO)=A%I(FROM)
          A%C(TO)=A%C(FROM)
        END IF
        CALL SORT__FLIP(FROM,TO)
      ENDDO
      CALL SORT__UNSET
!
!     ==========================================================================
!     ==  COMBINE IDENTICAL ENTRIES                                           ==
!     ==========================================================================
      J=1
      DO I=2,N
        TDIFF=(J.EQ.0)
        IF(.NOT.TDIFF) THEN
          TDIFF=(A%I(I).NE.A%I(J))
        END IF
        IF(TDIFF) THEN
          J=J+1
          A%I(J)=A%I(I)
          A%C(J)=A%C(I)
        ELSE
          A%C(J)=A%C(J)+A%C(I)
          A%C(I)=(0.D0,0.D0)
          IF(ABS(A%C(J)).EQ.0.D0) THEN
            J=J-1
          END IF
        END IF
      ENDDO
      A%N=J
      RETURN
      END SUBROUTINE CI_CLEANA
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$HAMILTONCOPY(HAM1,HAM2)
!     **************************************************************************
!     **  EXPRESSES THE HAMILTONIAN IN A NEW BASIS                            **
!     **            |CHI-PRIME_I>=SUM_J |CHI_J>*U(J,I)                        **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CIHAMIL_TYPE),INTENT(IN)    :: HAM1           
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM2           
      INTEGER(4)                       :: N1
!     **************************************************************************
      N1=HAM1%A%N
      IF(N1.GT.0) THEN
        IF(HAM2%A%NX.LT.N1)CALL CI_EXPANDA(HAM2%A,N1-HAM2%A%NX)
        HAM2%A%N=N1
        HAM2%A%I(:N1)=HAM1%A%I(:N1)
        HAM2%A%C(:N1)=HAM1%A%C(:N1)
      END IF
!
      N1=HAM1%H%N
      IF(N1.GT.0) THEN
        IF(HAM2%H%NX.LT.N1)CALL CI_EXPANDH(HAM2%H,N1-HAM2%H%NX)
        HAM2%H%N=HAM1%H%N
        HAM2%H%I(:N1)=HAM1%H%I(:N1)
        HAM2%H%J(:N1)=HAM1%H%J(:N1)
        HAM2%H%C(:N1)=HAM1%H%C(:N1)
      END IF
!
      N1=HAM1%U%N
      IF(N1.GT.0) THEN
        IF(HAM2%U%NX.LT.N1)CALL CI_EXPANDU(HAM2%U,N1-HAM2%U%NX)
        HAM2%U%N=HAM1%U%N
        HAM2%U%I(:N1)=HAM1%U%I(:N1)
        HAM2%U%J(:N1)=HAM1%U%J(:N1)
        HAM2%U%K(:N1)=HAM1%U%K(:N1)
        HAM2%U%L(:N1)=HAM1%U%L(:N1)
        HAM2%U%C(:N1)=HAM1%U%C(:N1)
      END IF
      RETURN
      END
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$HAMILTONTRANSFORM(NCHI,U,HAM)
!     **************************************************************************
!     **  EXPRESSES THE HAMILTONIAN IN A NEW BASIS                            **
!     **               |CHI-PRIME_I>=SUM_J |CHI_J>*U(J,I)
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)   :: NCHI
      COMPLEX(8)  ,INTENT(IN)   :: U(NCHI,NCHI)
      TYPE(CIHAMIL_TYPE),INTENT(INOUT) :: HAM           
      REAL(8)     ,PARAMETER    :: CMIN=1.D-8
      TYPE(CIHAMIL_TYPE)        :: HAMOLD
      INTEGER(4)                :: I,J,K,L,I1,J1,K1,L1,M
      COMPLEX(8)                :: CSVAR
!     **************************************************************************
      CALL CI$HAMILTONCOPY(HAM,HAMOLD)
      HAM%A%N=0
      HAM%H%N=0
      HAM%U%N=0
!
!     ==========================================================================
!     == TRANSFORM SOURCE TERMS OF THE HAMILTONIAN                            ==
!     ==========================================================================
      DO I=1,NCHI
        CSVAR=0.D0
        DO M=1,HAMOLD%A%N
          I1=HAMOLD%A%I(M)
          CSVAR=CSVAR+U(I,I1)*HAMOLD%A%C(M)
        ENDDO
        IF(ABS(CSVAR).LT.CMIN) CYCLE
        CALL CI_SETAELEMENT(HAM%A,I,CSVAR)
      ENDDO
      CALL CI_CLEANA(HAM%A)
!
!     ==========================================================================
!     == TRANSFORM ONE-CENTER PART OF THE  HAMILTONIAN                        ==
!     ==========================================================================
      DO I=1,NCHI
        DO J=1,NCHI
          CSVAR=0.D0
          DO M=1,HAMOLD%H%N
            I1=HAMOLD%H%I(M)
            J1=HAMOLD%H%J(M)
            CSVAR=CSVAR+CONJG(U(I,I1))*U(J,J1)*HAMOLD%H%C(M)
          ENDDO
          IF(ABS(CSVAR).LT.CMIN) CYCLE
          CALL CI_SETHELEMENT(HAM%H,I,J,CSVAR)
        ENDDO
      ENDDO
      CALL CI_CLEANH(HAM%H)
!
!     ==========================================================================
!     == TRANSFORM INTERACTION PART OF THE  HAMILTONIAN                       ==
!     ==========================================================================
      DO I=1,NCHI
        DO J=1,NCHI
          DO K=1,NCHI
            DO L=1,NCHI
              CSVAR=(0.D0,0.D0)
              DO M=1,HAMOLD%U%N
                I1=HAMOLD%U%I(M)
                J1=HAMOLD%U%J(M)
                K1=HAMOLD%U%K(M)
                L1=HAMOLD%U%L(M)
                CSVAR=CSVAR+CONJG(U(I,I1)*U(J,J1))*U(K,K1)*U(L,L1)*HAMOLD%U%C(M)
              ENDDO
              IF(ABS(CSVAR).LT.CMIN) CYCLE
              CALL CI_SETUELEMENT(HAM%U,I,J,K,L,CSVAR)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL CI_CLEANU(HAM%U)
      RETURN 
      END SUBROUTINE CI$HAMILTONTRANSFORM
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
      CALL CI$NEWPSI(PSI1,PSI%N)
      CALL CI$NEWPSI(PSI2,PSI%N)
      CALL CI$NEWPSI(PSI3,PSI%N/2)
      CALL CI$NEWPSI(PSI4,PSI%N/4)
!
!     ==========================================================================
!     == APPLY U-TENSOR                                                       ==
!     == W=0.5D0*SUM_{I,J,K,L} W_{I,J,K,L} CDAGGER_I CDAGGER_J C_L C_K        ==
!     ==========================================================================
      CALL CI$NEWPSI(HPSI,PSI%N)
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
!     == APPLY CREATORS AND ANNIHILATORS                                      ==
!     ==========================================================================
      N=HAM%A%N
      DO I1=1,N
        I=HAM%A%I(I1)
        CALL CI$COPYPSI(PSI,PSI1)
        CALL CI$CREATOR(PSI1,I)
        CALL CI$SCALEPSI(PSI1,HAM%A%C(I1))
        CALL CI$ADDPSI(HPSI,PSI1)
!
        CALL CI$COPYPSI(PSI,PSI1)
        CALL CI$ANNIHILATOR(PSI1,I)
        CALL CI$SCALEPSI(PSI1,CONJG(HAM%A%C(I1)))
        CALL CI$ADDPSI(HPSI,PSI1)
      ENDDO
!
!     ==========================================================================
!     == COPY RESULT INTO FINAL ARRAY                                         ==
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
      SUBROUTINE CI$1PDENMAT(IORB1,IORB2,PSI,CVAL)
!     **************************************************************************
!     **  MATRIX ELEMENT OF THE REDUCED ONE-PARTICLE DENSITY MATRIX           **
!     **       CVAL=RHO(I,J)=<PSI|CDAGGER_J C_I |PSI>                         **
!     **                                                                      **
!     **  PSI IS INTENT(INOUT) BECAUSE IT MAY BE CLEANED UP                   **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSI  ! MODIFIED BY CLEANING
      INTEGER(4)        ,INTENT(IN)   :: IORB1
      INTEGER(4)        ,INTENT(IN)   :: IORB2
      COMPLEX(8)        ,INTENT(OUT)  :: CVAL
      TYPE(CISTATE_TYPE)              :: PSI1
!     **************************************************************************
      CALL CI$NEWPSI(PSI1,PSI%N)
      CALL CI$COPYPSI(PSI,PSI1)
      CALL CI$ANNIHILATOR(PSI1,IORB1)
      CALL CI$CREATOR(PSI1,IORB2)
      CALL CI$SCALARPRODUCT(PSI,PSI1,CVAL)
      CALL CI$DELETEPSI(PSI1)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$STATEFROMDENSITYMATRIX(NCHI,NPSI,RHO,P,PSI)
!     **************************************************************************
!     ** CONSTRUCTS A MANY-PARTICLE WAVE FUNCTION SUCH THAT ITS ONE-PARTICLE  **
!     ** DENSITY MATRIX IS EQUAL TO A SPECIFIED ONE (RHO)                     **
!     **                                                                      **
!     ** THE WAVE FUNCTIONS CONSISTS OF A STATISTICAL AVERAGE OF N-PARTICLE   **
!     ** SLATER DETERMINANTS FOR DIFFERENT N                                  **
!     **                                                                      **
!     ** PSI MUST BE CREATED BY NEWPSI BEFORE CALLING THIS ROUTINE            **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)   :: NCHI           ! #(1-P ORBITALS)
      INTEGER(4)        ,INTENT(IN)   :: NPSI           ! #(WAVE FUNCTIONS)
      COMPLEX(8)        ,INTENT(IN)   :: RHO(NCHI,NCHI) ! DENSITY MATRIX
      REAL(8)           ,INTENT(OUT)  :: P(NPSI)        ! PROBABILITIES
      TYPE(CISTATE_TYPE),INTENT(OUT)  :: PSI(NPSI)      ! WAVE FUNCTIONS
      COMPLEX(8)                      :: U(NCHI,NCHI)
      REAL(8)                         :: FN(0:NCHI+1)
      REAL(8)                         :: EIG(NCHI)
      TYPE(CISTATE_TYPE)              :: PSIN
      TYPE(CISTATE_TYPE)              :: PSIADD2
      TYPE(CISTATE_TYPE)              :: PSICOPY
      COMPLEX(8)                      :: CSVAR,CSVAR1
      REAL(8)                         :: DIFF
      REAL(8)                         :: C(0:NCHI,NPSI)
      REAL(8)                         :: FI(0:NCHI+1,NPSI)
      INTEGER(4)                      :: MAP(0:NCHI+1,NPSI)
      INTEGER(4)                      :: I,J,N,K,N1
      LOGICAL(4)       ,PARAMETER     :: TTEST=.TRUE.
!     **************************************************************************
      CALL CI$NEWPSI(PSIN,NCHI*NCHI)
      CALL CI$NEWPSI(PSIADD2,NCHI*NCHI)
      CALL CI$NEWPSI(PSICOPY,NCHI*NCHI)
!
!     ==========================================================================
!     == DIAGONALIZE DENSITY MATRIX                                           ==
!     ==========================================================================
      CALL LIB$DIAGC8(NCHI,RHO,EIG,U)
      DO I=1,NCHI
        FN(I)=EIG(NCHI+1-I)
        FN(I)=MAX(0.D0,FN(I))
        FN(I)=MIN(1.D0,FN(I))
      ENDDO
      FN(NCHI+1)=0.D0
      FN(0)=1.D0
!
!     ==========================================================================
!     == CONSTRUCT WAVE FUNCTION                                              ==
!     == PURE SLATER DETERMINANTS WITH IDEAL PARTICLE NUMBER FOR I=0,NPSI-1   ==
!     == UNCORRELATED MULTIPARTICLE NUMBER STATE FOR I=NPSI                   ==
!     ==========================================================================
      FI(:,:)=0.D0
      DO I=0,NCHI+1
        MAP(I,:)=I
      ENDDO
      DO I=1,NPSI-1
        N=NINT(SUM(FN(1:NCHI)))
        P(I)=MIN(FN(MAP(N,I)),1.D0-FN(MAP(N+1,I)))
        FI(0,I)=1.D0
        DO J=1,N
          FI(MAP(J,I),I)=1.D0
        ENDDO
        FN(:)=FN(:)-FI(:,I)*P(I)
        IF(1.D0-P(I).LE.0.D0) THEN
          CALL ERROR$STOP('CI$STATEFROMDENSITYMATRIX')
        END IF
        FN(:)=FN(:)/(1.D0-P(I))
        FN(:)=MIN(FN(:),1.D0)
        FN(:)=MAX(FN(:),0.D0)
!       == ORDER FN(MAP(I)) IN DECREASING ORDER 
        DO J=1,NCHI
          DO K=J+1,NCHI
            IF(FN(MAP(K,I+1)).GT.FN(MAP(J,I+1))) THEN
              N=MAP(J,I+1)
              MAP(J,I+1)=MAP(K,I+1)
              MAP(K,I+1)=N
            END IF
          ENDDO
        ENDDO           
      ENDDO
      P(NPSI)=1.D0
      FI(:,NPSI)=FN(:)
      DO I=2,NPSI
        P(I)=P(I)*(1.D0-SUM(P(:I-1)))
      ENDDO
!
PRINT*,' SUM OF PROBABILITIES ',SUM(P)
DO I=1,NPSI
  WRITE(*,FMT='("P=",F7.5," N=",F7.3," FI=",20F7.4)')P(I),SUM(FI(1:NCHI,I)),FI(:,I)
ENDDO
!     == CONVERT OOCUPATIONS INTO CONSTANTS ====================================
      DO I=1,NPSI
        DO J=0,NCHI
          C(MAP(J,I),I)=SQRT(FI(MAP(J,I),I)-FI(MAP(J+1,I),I))
        ENDDO
      ENDDO
!
DO I=1,NPSI
  WRITE(*,FMT='("P=",F7.5," N=",F7.3," C=",20F7.4)')P(I),SUM(C(:,I)**2),C(:,I)
ENDDO
!
!     ==========================================================================
!     == CONSTRUCT WAVE FUNCTION                                              ==
!     ==========================================================================
      DO I=1,NPSI
        CALL CI$ZEROPSI(PSI(I))
        CALL CI$ZEROPSI(PSIN)
        CALL CI$SETPSI(PSIN,1,(1.D0,0.D0))  ! VACUUM STATE
!       == ADD CONTRIBUTION OF THE VACCUM STATE
        IF(ABS(C(0,I)).GT.1.D-6) THEN
          CALL CI$COPYPSI(PSIN,PSICOPY)       
          CSVAR=CMPLX(C(0,I),KIND=8)
          CALL CI$SCALEPSI(PSICOPY,CSVAR)
          CALL CI$ADDPSI(PSI(I),PSICOPY)
        END IF
        DO N1=1,NCHI
          N=MAP(N1,I)
!         ======================================================================
!         == APPLY NEXT CREATOR IN THE BASIS OF EIGENSTATES
!         ======================================================================
          CALL CI$ZEROPSI(PSIADD2)
          DO K=1,NCHI
            IF(ABS(U(K,NCHI+1-N)).LT.1.D-6) CYCLE
            CALL CI$COPYPSI(PSIN,PSICOPY)  
            CALL CI$CREATOR(PSICOPY,K)
            CALL CI$SCALEPSI(PSICOPY,CONJG(U(K,NCHI+1-N)))
            CALL CI$ADDPSI(PSIADD2,PSICOPY)
          ENDDO
          CALL CI$COPYPSI(PSIADD2,PSIN)       
!
!         ======================================================================
!         ==  ADD CONTRIBUTION OF N-PARTICLE STATE TO WAVE FUNCTION           ==
!         ======================================================================
          IF(ABS(C(N,I)).LT.1.D-6) CYCLE
          CALL CI$COPYPSI(PSIN,PSICOPY)       
          CSVAR=CMPLX(C(N,I),KIND=8)
          CALL CI$SCALEPSI(PSICOPY,CSVAR)
          CALL CI$ADDPSI(PSI(I),PSICOPY)
        ENDDO
        CALL CI$CLEANPSI(PSI(I))
      ENDDO       
      CALL CI$DELETEPSI(PSIN)
      CALL CI$DELETEPSI(PSIADD2)
      CALL CI$DELETEPSI(PSICOPY)
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        DO I=1,NPSI
          CALL CI$WRITEPSI(PSI(I),6)
        ENDDO
        DO I=1,NCHI
          DO J=1,NCHI
            CSVAR=(0.D0,0.D0)
            DO K=1,NPSI
              CALL CI$1PDENMAT(I,J,PSI(K),CSVAR1)
              CSVAR=CSVAR+P(K)*CSVAR1
            ENDDO
            DIFF=ABS(CSVAR-RHO(I,J))
            IF(DIFF.GT.1.D-7) THEN
              WRITE(*,'(2I3,"DEV.:",E10.2," FROM WV.",2F10.7,"FROM INPUT",2F10.7)') &
     &                                                 I,J,DIFF,RHO(I,J),CSVAR
            END IF
          ENDDO
        ENDDO
        DO I=1,NPSI
          DO J=1,NPSI
            CALL CI$SCALARPRODUCT(PSI(I),PSI(J),CSVAR)
            IF(I.EQ.J) CSVAR=CSVAR-(1.D0,0.D0)
            DIFF=ABS(CSVAR)
            IF(DIFF.GT.1.D-7) THEN
              WRITE(*,'(2I3,"DEV.:",E10.2," <PSI_I|PSI_J>-DELTA_IJ=",2E10.3)')I,J,DIFF,CSVAR
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
      CALL CI$NEWPSI(PSI1,PSI0%N)
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
!     **  EPOT=<PSI|HPSI>/<PSI|PSI>                                           **
!     **                                                                      **
!     **  REMARK: PSI AND HPSI ARE INTENT(INOUT) ONLY BECAUSE OF THE CLEAN    **
!     **          FUNCTION IN CI$SCALARPRODUCT. VALUES REMAIN UNCHANGED       **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PSI
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: HPSI
      REAL(8)           ,INTENT(OUT)   :: EPOT
      COMPLEX(8)                       :: CSVAR
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
!     **  PSIP AND PSIM ARE INTENT(INOUT) BECAUSE THEY MAY BE CLEANED UP      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      REAL(8)           ,INTENT(IN)   :: DT    ! TIME STEP
      REAL(8)           ,INTENT(IN)   :: MPSI  ! WAVE FUNCTION MASS
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSIP  ! |PSI(+)>
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSIM  ! |PSI(-)>
      REAL(8)           ,INTENT(OUT)  :: EKIN
      TYPE(CISTATE_TYPE)              :: PSICOPY
      COMPLEX(8)        ,PARAMETER    :: CMINUS=(-1.D0,0.D0)
      COMPLEX(8)                      :: CSVAR
!     **************************************************************************
      CALL CI$NEWPSI(PSICOPY,PSIM%N)
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
      SUBROUTINE CI_LAGRANGE(NCHI,NPSI,X0,PSI0,XP,PSIP,QX,QPSI,RHO,H,LAMBDA,IERR)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)   :: NCHI
      INTEGER(4)        ,INTENT(IN)   :: NPSI
      REAL(8)           ,INTENT(IN)   :: X0(NPSI)  ! CURRENT PROBABILITIES
      REAL(8)           ,INTENT(INOUT):: XP(NPSI)  ! CURRENT PROBABILITIES
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSI0(NPSI) !(ACTUALLY INTENT(IN)
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSIP(NPSI)
      REAL(8)           ,INTENT(IN)   :: QX
      REAL(8)           ,INTENT(IN)   :: QPSI(NPSI)
      COMPLEX(8)        ,INTENT(IN)   :: RHO(NCHI,NCHI)
      COMPLEX(8)        ,INTENT(OUT)  :: H(NCHI,NCHI)
      COMPLEX(8)        ,INTENT(OUT)  :: LAMBDA(NPSI,NPSI)
      INTEGER(4)        ,INTENT(OUT)  :: IERR
      INTEGER(4)                      :: NC   !#(CONSTRAINTS)
      REAL(8)                         :: MAT(NCHI**2+NPSI**2+1,NCHI**2+NPSI**2+1)
      REAL(8)                         :: MATINV(NCHI**2+NPSI**2+1,NCHI**2+NPSI**2+1)
      REAL(8)                         :: U(NCHI**2+NPSI**2+1,NCHI**2+NPSI**2+1)
      REAL(8)                         :: S(NCHI**2+NPSI**2+1)
      REAL(8)                         :: V(NCHI**2+NPSI**2+1,NCHI**2+NPSI**2+1)
      REAL(8)                         :: CNSTRNTVEC(NCHI**2+NPSI**2+1)
      REAL(8)                         :: LAGRNGEVEC(NCHI**2+NPSI**2+1)
      REAL(8)                         :: VEC(NCHI**2+NPSI**2+1)
      REAL(8)                         :: VECSUM(NCHI**2+NPSI**2+1)
      INTEGER(4)        ,PARAMETER    :: NITER=1000
      REAL(8)           ,PARAMETER    :: TOL=1.D-6
      REAL(8)           ,PARAMETER    :: DEVMAX=1.D+1
      LOGICAL(4)        ,PARAMETER    :: TTEST=.FALSE.
      REAL(8)                         :: DEV,DEVLAST
      INTEGER(4)                      :: ITER
      INTEGER(4)                      :: I,J,IND
      LOGICAL(4)                      :: CONVG
      TYPE(CISTATE_TYPE)              :: PSITEST(NPSI)
REAL(8)  :: XTEST(NPSI)
REAL(8)  :: STEP
REAL(8)  :: VECX(NCHI**2+NPSI**2+1)
REAL(8)  :: MATTEST(NCHI**2+NPSI**2+1,NCHI**2+NPSI**2+1)
REAL(8)  :: MAT1(NCHI**2+1,NCHI**2+1)
REAL(8)  :: QPSI1(NPSI),QX1
!     **************************************************************************
      IERR=0
      NC=NCHI**2+NPSI**2+1
!
!     ==========================================================================
!     == CONSTRUCT MATRIX FOR ITERATION                                       ==
!     == MAT IS THE REAL (NC*NC) MATRIX WITH ELEMENTS                         ==
!     ==    <PSI0| CONSTR(I)*CONSTR(J) |PSI0>                                 ==
!     ==========================================================================
PRINT*,'BEFORE LAGRANGEMAT'
      CALL CI_LAGRANGEMAT(NCHI,NPSI,QX,X0,QPSI,PSI0,NC,MAT)
PRINT*,'AFTER LAGRANGEMAT'
!
!     ==========================================================================
!     == TEST LAGRANGEMAT ======================================================
!     ==========================================================================
      IF(TTEST) THEN
        PRINT*,'NCHI ',NCHI,' NPSI=',NPSI,' NC=',NC
        PRINT*,'X0=',X0,' X0^2=',X0**2
        QPSI1(:)=1.D0
        QX1=1.D0
        CALL CI_LAGRANGEMAT(NCHI,NPSI,QX1,X0,QPSI1,PSI0,NC,MAT)
!       == PRINT MATRIX ========================================================
        DO I=1,NC
          WRITE(*,FMT='("MAT ",100F10.6)')MAT(I,:)
        ENDDO
        WRITE(*,*)
!       == WRITE CONSTRAINTS OF REFERENCE STATE ================================
        CALL CI_CONSTRAINT(NCHI,NPSI,X0,PSI0,RHO,NC,CNSTRNTVEC)
        WRITE(*,FMT='("REF.",I3,100F10.6)')0,CNSTRNTVEC
        PRINT*,' '
!       == COMPARE WITH NUMERICAL DIFFERENTIATION ==============================
        WRITE(*,FMT='("<>",I3,100I10)')0,(I,I=1,NC)
        DO IND=1,NC
          STEP=1.D-6
          LAGRNGEVEC(:)=0.D0
          LAGRNGEVEC(IND)=STEP
          DO I=1,NPSI
            CALL CI$ZEROPSI(PSITEST(I))      
            CALL CI$COPYPSI(PSI0(I),PSITEST(I))
            XTEST(I)=X0(I)
          ENDDO
          CALL CI_ADDCONSTRAINT(NCHI,NPSI,QX1,QPSI1,XTEST,PSITEST,X0,PSI0,NC,LAGRNGEVEC)
          CALL CI_CONSTRAINT(NCHI,NPSI,XTEST,PSITEST,RHO,NC,VEC)
          VECX = CNSTRNTVEC + MATMUL(MAT,LAGRNGEVEC)
          WRITE(*,FMT='("NUM",I3,100F10.6)')IND,VEC/STEP
          WRITE(*,FMT='("AN.",I3,100F10.6)')IND,VECX/STEP,MAXVAL(ABS(VEC-VECX))/STEP
          WRITE(*,*)
          IF(IND.EQ.16)WRITE(*,*)
        ENDDO
!
        MAT1=0.D0
        CALL CI_LAGRANGEMATOLD(NCHI,PSI0(1),NCHI**2+1,MAT1)
        DO I=1,NCHI**2+1
          WRITE(*,FMT='("MAT1 ",100F10.6)')MAT1(I,:)
        ENDDO
        WRITE(*,*)
!       == CLEAN UP ============================================================
        DO I=1,NPSI
          CALL CI$ZEROPSI(PSITEST(I))      
        ENDDO
        STOP 'FORCED AFTER TEST IN CI$LAGRANGE'
      END IF
!
!     ==========================================================================
!     == INVERT MATRIX FOR LOOP                                               ==
!     ==========================================================================
      CALL LIB$SVDR8(NC,NC,MAT,U,S,V)
      DO I=1,NC
        IF(S(I).GT.1.D-4.AND.S(I).LT.1.D+6) THEN
          S(I)=1.D0/S(I)
        ELSE
          S(I)=0.D0
        END IF
      ENDDO
      DO I=1,NC
        U(:,I)=U(:,I)*S(I)
      ENDDO
      MATINV=MATMUL(U,V)
!
!     ==========================================================================
!     == CORRECT THE PSEUDO INVERSE SO THAT IT DOES NOT MIX IMAGINARY STATES ==
!     ==========================================================================
      CALL CI_FILTERIMAGINARY(NCHI,NPSI,MATINV)
!!$DO I=1,NC
!!$ WRITE(*,FMT='("MATINV ",30F7.2)')MATINV(I,:)
!!$ENDDO
!!$WRITE(*,*)
!!$MATTEST=MATMUL(MATINV,MAT)
!!$DO I=1,NC
!!$ WRITE(*,FMT='("MATINV*MAT ",30F7.2)')MATTEST(I,:)
!!$ENDDO
!!$WRITE(*,*)
!!$STOP
!
!     ==========================================================================
!     == LOOP TO ENFORCE CONSTRAINTS                                          ==
!     ==========================================================================
PRINT*,'BEFORE ITERATION'
      DEVLAST=HUGE(DEVLAST)
      VECSUM=0.D0
      DO ITER=1,NITER
!       == CALCULATE DEVIATION OF |PSI(+)> FROM CONSTRAINT =====================
!       == <PSI(+)|CONSTR(I)|PSI(+)>-VALUE(I) ==================================
        CALL CI_CONSTRAINT(NCHI,NPSI,XP,PSIP,RHO,NC,CNSTRNTVEC)
!
!       == DETERMINE CORRECTION FOR LAGRANGE MULTIPLIERS =======================
        LAGRNGEVEC(:)=-MATMUL(MATINV,CNSTRNTVEC)
!
!       == IF CONVERGENCE FAILS, JUMP OUT OF THE LOOP ==========================
!       == IF IT FAILS, IT DOES SO EXPONENTIALLY, SO THAT THE CHOICE OF THE   ==
!       == THRESHHOLD IS NOT SIGNIFICANT =======================================
        DEV=SUM(ABS(CNSTRNTVEC))
!        IF(DEV.GT.DEVLAST) EXIT
        IF(DEV.GT.DEVMAX) THEN
          IERR=1
          CALL ERROR$MSG('DEVIATION UNACCEPTABLE')
          CALL ERROR$STOP('CI$LAGRANGE')
!          RETURN
        END IF
        CONVG=MAXVAL(ABS(CNSTRNTVEC)).LT.TOL
        IF(CONVG) EXIT
        IF(MAXVAL(ABS(LAGRNGEVEC)).LT.1.D-10) EXIT
        DEVLAST=DEV
!
!       == ADD CORRECTION TO |PSI(+)> == =======================================
        CALL CI_ADDCONSTRAINT(NCHI,NPSI,QX,QPSI,XP,PSIP,X0,PSI0,NC,LAGRNGEVEC)
        VECSUM=VECSUM+LAGRNGEVEC(:)       ! ADD CORRECTIONS UP TO OBTAIN TOTAL
!
!!$VECX=CNSTRNTVEC+MATMUL(MAT,LAGRNGEVEC)
!!$CALL CI_CONSTRAINT(NCHI,NPSI,XP,PSIP,RHO,NC,VEC)
!!$WRITE(*,FMT='("ITER ",I5," MAXDEV=",4F20.10," DLAGR=",2F20.10)')ITER &
!!$,MAXVAL(ABS(VEC)),MAXVAL(ABS(VECX)),MAXVAL(ABS(VEC))-MAXVAL(ABS(VECX)) &
!!$,MAXVAL(ABS(CNSTRNTVEC)),MAXVAL(ABS(LAGRNGEVEC)),MAXVAL(ABS(VECSUM))
      ENDDO
      IF(.NOT.CONVG) THEN
!PRINT*,'LOOP NOT CONVERGED ',ITER,CONVG
!!$        CALL ERROR$MSG('LOOP NOT CONVERGED')
!!$        CALL ERROR$STOP('CI_LAGRANGE')
      END IF
PRINT*,'AFTER ITERATION'
!
!     ==========================================================================
!     ==  CALCULATES LAGRANGE MULTIPLICATORS                                  ==
!     ==  FRICTION AND TIME STEP NOT INCLUDED                                 ==
!     ==========================================================================
      H(:,:)=0.D0
      LAMBDA(:,:)=0.D0
      IND=0
      DO J=1,NCHI    ! BUGFIX MAR.22, 20012: INDICES I,J REVERSED. PB
        DO I=1,NCHI
          IND=IND+1
          IF(J.GE.I) THEN
            H(I,J)=H(I,J)+CMPLX(VECSUM(IND),0.D0,KIND=8)
            IF(I.NE.J)H(J,I)=H(J,I)+CMPLX(VECSUM(IND),0.D0,KIND=8)
          ELSE 
            H(I,J)=H(I,J)+CMPLX(0.D0,VECSUM(IND),KIND=8)
            H(J,I)=H(J,I)-CMPLX(0.D0,VECSUM(IND),KIND=8)
          END IF
        ENDDO
      ENDDO
      DO J=1,NPSI  ! BUGFIX MAR.22, 20012: INDICES I,J REVERSED. PB 
        DO I=1,NPSI
          IND=IND+1
          IF(J.GE.I) THEN
            LAMBDA(J,I)=LAMBDA(J,I)+CMPLX(VECSUM(IND),0.D0,KIND=8)
            IF(I.NE.J)LAMBDA(I,J)=LAMBDA(I,J)+CMPLX(VECSUM(IND),0.D0,KIND=8)
          ELSE 
            LAMBDA(J,I)=LAMBDA(J,I)+CMPLX(0.D0,VECSUM(IND),KIND=8)
            LAMBDA(I,J)=LAMBDA(I,J)-CMPLX(0.D0,VECSUM(IND),KIND=8)
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_LAGRANGEMAT(NCHI,NPSI,QX,X0,QPSI,PSI0,NC,MAT)
!     **************************************************************************
!     **  DETERMINE LINEAR TERM OF EQUATION FOR LAGRANGE MULTIPLIERS          **
!     **                                                                      **
!     **  PSI0 IS INTENT(INOUT) BECAUSE IT MAY BECOME CLEANED UP              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)  :: NCHI
      INTEGER(4)        ,INTENT(IN)  :: NPSI
      REAL(8)           ,INTENT(IN)  :: QX
      REAL(8)           ,INTENT(IN)  :: X0(NPSI) 
      REAL(8)           ,INTENT(IN)  :: QPSI(NPSI)
      TYPE(CISTATE_TYPE),INTENT(INOUT) :: PSI0(NPSI)
      INTEGER(4)        ,INTENT(IN)  :: NC
      REAL(8)           ,INTENT(OUT) :: MAT(NC,NC)
      REAL(8)                        :: P0(NPSI)
      TYPE(CISTATE_TYPE)             :: PSI1,PSI2,PSI3,PSI4
      COMPLEX(8)                     :: CMAT4(NCHI,NCHI,NCHI,NCHI)
      COMPLEX(8)                     :: CMAT2(NCHI,NCHI,NPSI,NPSI)
      COMPLEX(8)                     :: CMAT2A(NCHI,NCHI)
      COMPLEX(8)                     :: CMAT2B(NCHI,NCHI)
      COMPLEX(8)                     :: CMAT3(NPSI,NPSI,NPSI,NPSI)
      COMPLEX(8)                     :: CMAT0(NPSI,NPSI)
      INTEGER(4)                     :: I,J,K,L,M,N,IND1,IND2
      REAL(8)                        :: S1,S2
      COMPLEX(8)                     :: CSVAR,CSVAR1
      COMPLEX(8)        ,PARAMETER   :: CI=(0.D0,1.D0)
      INTEGER(4)                     :: IINDEX(2,NC)
!     **************************************************************************
      IF(NC.NE.NCHI**2+NPSI**2+1) THEN
        CALL ERROR$STOP('CI_LAGRANGEMAT')
      END IF
      P0(:)=X0(:)**2
!
!     ==========================================================================
!     == DETERMINE EXPECTATION VALUES OF OPERATORS                            ==
!     ==========================================================================
!     == CMAT0(N,M)=<PSI(N)|PSI(M)> ============================================
      CMAT0(:,:)=(0.D0,0.D0)
      DO N=1,NPSI
        DO M=N,NPSI
          CALL CI$SCALARPRODUCT(PSI0(N),PSI0(M),CSVAR)
          CMAT0(N,M)=CSVAR
          CMAT0(M,N)=CONJG(CSVAR)
        ENDDO
      ENDDO
!
!     ==  CMAT2(I,J,M,N)=<PSI(M)|CDAGGER(I)C(J)|\PSI(N)> =======================
!     ==  CMAT4(I,J,K,L)=SUM_N P(N)<PSI(N)|CDAGGER(I)C(J)CDAGGER(K)C(L)|\PSI(N)>
      N=MAXVAL(PSI0(:)%N)
      CALL CI$NEWPSI(PSI1,N)
      CALL CI$NEWPSI(PSI2,N)
      CALL CI$NEWPSI(PSI3,N/2)
      CALL CI$NEWPSI(PSI4,N/4)
      CMAT4(:,:,:,:)=(0.D0,0.D0)
      CMAT2(:,:,:,:)=(0.D0,0.D0)
PRINT*,'MARKE 1'
      DO N=1,NPSI
        DO L=1,NCHI
          CALL CI$COPYPSI(PSI0(N),PSI1)
          CALL CI$ANNIHILATOR(PSI1,L)
          DO K=1,NCHI
            IF(K.GT.L) CYCLE
            CALL CI$COPYPSI(PSI1,PSI2)
            CALL CI$CREATOR(PSI2,K)
            DO M=1,NPSI
              CALL CI$SCALARPRODUCT(PSI0(M),PSI2,CMAT2(K,L,M,N))
            ENDDO
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
                CALL CI$SCALARPRODUCT(PSI0(N),PSI4,CSVAR)
                CMAT4(I,J,K,L)=CMAT4(I,J,K,L)+P0(N)*QPSI(N)*CSVAR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
PRINT*,'MARKE 2'
      CALL CI$DELETEPSI(PSI1)
      CALL CI$DELETEPSI(PSI2)
      CALL CI$DELETEPSI(PSI3)
      CALL CI$DELETEPSI(PSI4)
!
!     ==========================================================================
!     ==  COMPLETE MATRICES                                                   ==
!     ==========================================================================
!     == <PSI_M|CDAGGER_KC_L|PSI_N>=<PSI_N|CDAGGER_LC_K|PSI_N>* ================
      DO N=1,NPSI
        DO M=1,NPSI
          DO I=1,NCHI
            DO J=1,I-1
              CMAT2(I,J,N,M)=CONJG(CMAT2(J,I,M,N))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      CMAT2A(:,:)=(0.D0,0.D0)
      CMAT2B(:,:)=(0.D0,0.D0)
      DO N=1,NPSI
        CMAT2A(:,:)=CMAT2A(:,:)+QPSI(N)*P0(N)*CMAT2(:,:,N,N)
        CMAT2B(:,:)=CMAT2B(:,:)+4.D0*QX*P0(N)*CMAT2(:,:,N,N)
      ENDDO
!
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
              IF(I.EQ.J)CSVAR1=CSVAR1+CMAT2A(K,L)
              IF(J.EQ.K)CSVAR1=CSVAR1+CMAT2A(I,L)
              CMAT4(K,J,I,L)=CSVAR1
              CMAT4(L,I,J,K)=CONJG(CMAT4(K,J,I,L))
!
!             == PERMUTE ANNIHILATORS===========================================
              CSVAR1=-CMAT4(I,J,K,L)
              IF(K.EQ.L)CSVAR1=CSVAR1+CMAT2A(I,J)
              IF(J.EQ.K)CSVAR1=CSVAR1+CMAT2A(I,L)
              CMAT4(I,L,K,J)=CSVAR1
              CMAT4(J,K,L,I)=CONJG(CMAT4(I,L,K,J))
!
!             == PERMUTE CREATORS AND ANNIHILATORS =============================
              CSVAR1=CMAT4(I,J,K,L)
              IF(J.EQ.K)CSVAR1=CSVAR1-CMAT2A(I,L)
              IF(I.EQ.L)CSVAR1=CSVAR1+CMAT2A(K,J)
              CMAT4(K,L,I,J)=CSVAR1
              CMAT4(J,I,L,K)=CONJG(CMAT4(K,L,I,J))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  CONSTRUCT K4                                                        ==
!     ==========================================================================
      DO I=1,NCHI
        DO J=1,NCHI
          DO K=1,NCHI
            DO L=1,NCHI
              CSVAR=0.5D0*(CMAT4(I,J,K,L)+CMAT4(K,L,I,J))
              CMAT4(I,J,K,L)=CSVAR
              CMAT4(K,L,I,J)=CSVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO I=1,NCHI
        DO J=1,NCHI
          DO K=1,NCHI
            DO L=1,NCHI
              CSVAR=(0.D0,0.D0)
              DO N=1,NPSI
                CSVAR=CSVAR+P0(N)*CMAT2(I,J,N,N)*CMAT2(K,L,N,N) 
              ENDDO
              CMAT4(I,J,K,L)=2.D0*CMAT4(I,J,K,L) +4.D0*QX*CSVAR 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      DO N=1,NPSI
        DO M=1,NPSI
          CMAT2(:,:,M,N)=CMAT2(:,:,M,N)*(QPSI(M)+QPSI(N))
        ENDDO
      ENDDO
!
      CMAT3(:,:,:,:)=(0.D0,0.D0)
      DO I=1,NPSI
        DO J=1,NPSI
          DO K=1,NPSI
            L=I
            CMAT3(I,J,K,L)=CMAT3(I,J,K,L)+CMAT0(K,J)*QPSI(I)/P0(I) 
          ENDDO
          DO L=1,NPSI
            K=J
            CMAT3(I,J,K,L)=CMAT3(I,J,K,L)+CMAT0(I,L)*QPSI(J)/P0(J)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  MAP ONTO REAL MATRIX                                                ==
!     ==========================================================================
      MAT(:,:)=0.D0
      IND1=0
      DO J=1,NCHI
        DO I=1,NCHI
          IND1=IND1+1
          S1=1.D0
          IF(J.LT.I) S1=-1.D0 !IMAGINARY PART OF X(I,J)
          IND2=0
          DO L=1,NCHI
            DO K=1,NCHI
              IND2=IND2+1
              S2=1.D0
              IF(L.LT.K) S2=-1.D0
              CSVAR=0.25D0*(    (CMAT4(I,J,K,L)+S2*CMAT4(I,J,L,K) &
     &                      +S1*(CMAT4(J,I,K,L)+S2*CMAT4(J,I,L,K))))
              IF(S1.LT.0.D0) CSVAR=-CI*CSVAR   
              IF(S2.LT.0.D0) CSVAR=+CI*CSVAR
              MAT(IND1,IND2)=REAL(CSVAR)
            ENDDO
          ENDDO
!
!         ======================================================================
          DO N=1,NPSI
            DO M=1,NPSI
              IND2=IND2+1
              S2=1.D0
              IF(N.LT.M) S2=-1.D0
              CSVAR=0.25D0*(    (CMAT2(I,J,M,N)+S2*CMAT2(I,J,N,M)) &
    &                       +S1*(CMAT2(J,I,M,N)+S2*CMAT2(J,I,N,M)))
              IF(S1.LT.0.D0) CSVAR=-CI*CSVAR    
              IF(S2.LT.0.D0) CSVAR=+CI*CSVAR    
              MAT(IND1,IND2)=REAL(CSVAR)
              MAT(IND2,IND1)=REAL(CSVAR)
            ENDDO
          ENDDO
!
!         ======================================================================
          IND2=IND2+1
          CSVAR=0.5D0*(CMAT2B(I,J)+S1*CMAT2B(J,I))
          IF(S1.LT.0.D0) CSVAR=-CI*CSVAR        
          MAT(IND1,IND2)=REAL(CSVAR)
          MAT(IND2,IND1)=REAL(CSVAR)
        ENDDO
      ENDDO
!
      IND1=NCHI**2
      DO J=1,NPSI
        DO I=1,NPSI
          IND1=IND1+1
          S1=1.D0
          IF(J.LT.I) S1=-1.D0
          IND2=NCHI**2 
          DO L=1,NPSI
            DO K=1,NPSI
              IND2=IND2+1
              S2=1.D0
              IF(L.LT.K) S2=-1.D0
              CSVAR=0.25D0*(    (CMAT3(I,J,K,L)+S2*CMAT3(I,J,L,K)) &
      &                     +S1*(CMAT3(J,I,K,L)+S2*CMAT3(J,I,L,K)))
              IF(S1.LT.0.D0) CSVAR=-CI*CSVAR   
              IF(S2.LT.0.D0) CSVAR=+CI*CSVAR
              MAT(IND1,IND2)=REAL(CSVAR)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IND1=NCHI**2+NPSI**2+1
      MAT(IND1,IND1)=4.D0*QX*SUM(P0)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_FILTERIMAGINARY(NCHI,NPSI,MAT)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: NPSI
      REAL(8)   ,INTENT(INOUT) :: MAT(NCHI**2+NPSI**2+1,NCHI**2+NPSI**2+1)
      INTEGER(4)               :: I,J,K,L,IND1,IND2
      REAL(8)                  :: S1,S2
      LOGICAL(4)               :: TPR=.FALSE.
!     **************************************************************************
      IND1=0
      DO J=1,NCHI
        DO I=1,NCHI
          IND1=IND1+1
          S1=1.D0
          IF(J.LT.I) S1=-1.D0
          IND2=0
          DO L=1,NCHI
            DO K=1,NCHI
              IND2=IND2+1
              S2=1.D0
              IF(L.LT.K) S2=-1.D0
!              IF((S1.LT.0.D0.OR.S2.LT.0.D0).AND.MAT(IND1,IND2).NE.0.D0) THEN
              IF(S1*S2.LT.0.D0.AND.MAT(IND1,IND2).NE.0.D0) THEN
                 IF(TPR)PRINT*,'IMAGINARY ',IND1,IND2,MAT(IND1,IND2)
                 MAT(IND1,IND2)=0.D0
              END IF
            ENDDO
          ENDDO
          DO L=1,NPSI
            DO K=1,NPSI
              IND2=IND2+1
              S2=1.D0
              IF(L.LT.K) S2=-1.D0  
!              IF((S1.LT.0.D0.OR.S2.LT.0.D0).AND.MAT(IND1,IND2).NE.0.D0) THEN
              IF(S1*S2.LT.0.D0.AND.MAT(IND1,IND2).NE.0.D0) THEN
                 IF(TPR)PRINT*,'IMAGINARY ',IND1,IND2,MAT(IND1,IND2)
                 MAT(IND1,IND2)=0.D0
              END IF
            ENDDO
          ENDDO
          IND2=IND2+1 
!          IF(S1.LT.0.D0.AND.MAT(IND1,IND2).NE.0.D0) THEN
          IF(S1.LT.0.D0.AND.MAT(IND1,IND2).NE.0.D0) THEN
            IF(TPR)PRINT*,'IMAGINARY ',IND1,IND2,MAT(IND1,IND2)
            MAT(IND1,IND2)=0.D0
          END IF
         ENDDO
      ENDDO
      DO J=1,NPSI
        DO I=1,NPSI
          IND1=IND1+1
          S1=1.D0
          IF(J.LT.I) S1=-1.D0
          IND2=0
          DO L=1,NCHI
            DO K=1,NCHI
              IND2=IND2+1
              S2=1.D0
              IF(L.LT.K) S2=-1.D0
!              IF((S1.LT.0.D0.OR.S2.LT.0.D0).AND.MAT(IND1,IND2).NE.0.D0) THEN
              IF(S1*S2.LT.0.D0.AND.MAT(IND1,IND2).NE.0.D0) THEN
                 IF(TPR)PRINT*,'IMAGINARY ',IND1,IND2,MAT(IND1,IND2)
                 MAT(IND1,IND2)=0.D0
              END IF
            ENDDO
          ENDDO
          DO L=1,NPSI
            DO K=1,NPSI
              IND2=IND2+1
              S2=1.D0
              IF(L.LT.K) S2=-1.D0
!              IF((S1.LT.0.D0.OR.S2.LT.0.D0).AND.MAT(IND1,IND2).NE.0.D0) THEN
              IF(S1*S2.LT.0.D0.AND.MAT(IND1,IND2).NE.0.D0) THEN
                 IF(TPR)PRINT*,'IMAGINARY ',IND1,IND2,MAT(IND1,IND2)
                 MAT(IND1,IND2)=0.D0
              END IF
            ENDDO
          ENDDO
          IND2=IND2+1
          IF(S1.LT.0.D0.AND.MAT(IND1,IND2).NE.0.D0) THEN
             IF(TPR)PRINT*,'IMAGINARY ',IND1,IND2,MAT(IND1,IND2)
             MAT(IND1,IND2)=0.D0
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_CONSTRAINT(NCHI,NPSI,X,PSI,RHO,NC,VEC)
!     **************************************************************************
!     ** EVALUATES THE DEVIATION FROM THE CONSTRAINT CONDITIONS               **
!     **                                                                      **
!     ** PSI IS INTENT(INOUT) BECAUSE IT MAY BE CLEANED UP                    **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)   :: NCHI
      INTEGER(4)        ,INTENT(IN)   :: NPSI
      REAL(8)           ,INTENT(IN)   :: X(NPSI)
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSI(NPSI)
      COMPLEX(8)        ,INTENT(IN)   :: RHO(NCHI,NCHI) ! TARGET DENSITY MATRIX
      INTEGER(4)        ,INTENT(IN)   :: NC
      REAL(8)           ,INTENT(OUT)  :: VEC(NC)
      TYPE(CISTATE_TYPE)              :: PSI1,PSI2
      COMPLEX(8)                      :: CMAT2(NCHI,NCHI)
      COMPLEX(8)                      :: CMAT0(NPSI,NPSI)
      INTEGER(4)                      :: I,J,N,IND
      REAL(8)                         :: S
      REAL(8)                         :: P(NPSI)
      COMPLEX(8)                      :: CSVAR
      COMPLEX(8)       ,PARAMETER     :: CI=(0.D0,1.D0)
!     **************************************************************************
      IF(NC.NE.NCHI**2+NPSI**2+1) THEN
        CALL ERROR$STOP('CI_CONSTRAINT')
      END IF
      P(:)=X(:)**2
!
!     ==========================================================================
!     == DETERMINE EXPECTATION VALUES OF OPERATORS==============================
!     ==========================================================================
      N=MAXVAL(PSI(:)%N)
!     == CMAT0(I,J)=<PSI_I|PSI_J> ==============================================
      CALL CI$NEWPSI(PSI1,N)
      CALL CI$NEWPSI(PSI2,N)
      DO I=1,NPSI
        DO J=I,NPSI
          CALL CI$SCALARPRODUCT(PSI(I),PSI(J),CSVAR)
          CMAT0(I,J)=CSVAR
          CMAT0(J,I)=CONJG(CSVAR)
        ENDDO
      ENDDO
!
!     == RHO(J,I)=SUM_N P_N <PSI_N|CDAGGER_I C_J |PSI_N> =======================
      CMAT2(:,:)=(0.D0,0.D0)
      DO N=1,NPSI
        DO J=1,NCHI
          CALL CI$COPYPSI(PSI(N),PSI1)
          CALL CI$ANNIHILATOR(PSI1,J)
          DO I=1,NCHI
            CALL CI$COPYPSI(PSI1,PSI2)
            CALL CI$CREATOR(PSI2,I)
            CALL CI$SCALARPRODUCT(PSI(N),PSI2,CSVAR)
            CMAT2(I,J)=CMAT2(I,J)+P(N)*CSVAR
          ENDDO
        ENDDO
      ENDDO
      CALL CI$DELETEPSI(PSI1)
      CALL CI$DELETEPSI(PSI2)
!
!     ==========================================================================
!     ==  SUBTRACT TARGET CONSTRAINT VALUE                                    ==
!     ==========================================================================
      CMAT2(:,:)=CMAT2(:,:)-TRANSPOSE(RHO(:,:))
      DO I=1,NPSI
        CMAT0(I,I)=CMAT0(I,I)-(1.D0,0.D0)
      ENDDO

!     ==========================================================================
!     ==  MAP ONTO REAL VECTOR                                                ==
!     ==========================================================================
      IND=0
      DO J=1,NCHI
        DO I=1,NCHI
          IND=IND+1
          S=1.D0
          IF(J.LT.I) S=-1.D0
          CSVAR=0.5D0*(CMAT2(I,J)+S*CMAT2(J,I))
          IF(S.LT.0.D0) CSVAR=-CI*CSVAR   
          VEC(IND)=REAL(CSVAR)
        ENDDO
      ENDDO
      DO J=1,NPSI
        DO I=1,NPSI
          IND=IND+1
          S=1.D0
          IF(J.LT.I) S=-1.D0
          CSVAR=0.5D0*(CMAT0(I,J)+S*CMAT0(J,I))
          IF(S.LT.0.D0) CSVAR=-CI*CSVAR   
          VEC(IND)=REAL(CSVAR)
        ENDDO
      ENDDO
      IND=IND+1
      VEC(IND)=SUM(P)-1.D0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_ADDCONSTRAINT(NCHI,NPSI,QX,QPSI,XP,PSIP,X0,PSI0,NC,VEC)
!     **************************************************************************
!     ** ADD CONSTRAINTS FORCES TO PSIP                                       **
!     ** CONSTRAINT FORCES ARE EVALUATED FOR THE CURRENT WAVE FUNCTIONS AND   **
!     ** THE CURRENT PROBABILITIES P0=P(0)                                    **
!     **                                                                      **
!     ** PSI0 IS INTENT(OINOUT) BECAUSE IT MAY BECOME CLEANED UP              **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)   :: NCHI
      INTEGER(4)        ,INTENT(IN)   :: NPSI
      REAL(8)           ,INTENT(IN)   :: QX
      REAL(8)           ,INTENT(IN)   :: QPSI(NPSI)
      REAL(8)           ,INTENT(IN)   :: X0(NPSI)
      REAL(8)           ,INTENT(INOUT):: XP(NPSI)
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSIP(NPSI)
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSI0(NPSI)
      INTEGER(4)        ,INTENT(IN)   :: NC
      REAL(8)           ,INTENT(IN)   :: VEC(NC)
      TYPE(CISTATE_TYPE)              :: PSI1
      TYPE(CISTATE_TYPE)              :: PSI2
      REAL(8)                         :: P0(NPSI)
      REAL(8)                         :: Y
      INTEGER(4)                      :: I,J,N,M,IND
      COMPLEX(8)                      :: H(NCHI,NCHI)
      COMPLEX(8)                      :: LAMBDA(NPSI,NPSI)
      COMPLEX(8)                      :: CSVAR
      COMPLEX(8)        ,PARAMETER    :: CI=(0.D0,1.D0)
!     **************************************************************************
      IF(NC.NE.NCHI**2+NPSI**2+1) THEN
        CALL ERROR$STOP('CI_ADDCONSTRAINT')
      END IF
      P0=X0**2
!
!     ==========================================================================
!     ==  EXTRACT LAGRANGE MULTIPLIERS IN THEIR COMPLEX FORM                  ==
!     ==========================================================================
      H(:,:)=(0.D0,0.D0)
      LAMBDA(:,:)=(0.D0,0.D0)
      IND=0
      DO J=1,NCHI
        DO I=1,NCHI
          IND=IND+1
          IF(J.LT.I) THEN
            H(I,J)=H(I,J)+CI*VEC(IND)
            H(J,I)=H(J,I)-CI*VEC(IND)
          ELSE
            H(I,J)=H(I,J)+VEC(IND)
            IF(I.NE.J)H(J,I)=H(J,I)+VEC(IND)
          END IF
        ENDDO
      ENDDO
      DO J=1,NPSI
        DO I=1,NPSI
          IND=IND+1
          IF(J.LT.I) THEN
            LAMBDA(I,J)=LAMBDA(I,J)+CI*VEC(IND)
            LAMBDA(J,I)=LAMBDA(J,I)-CI*VEC(IND)
          ELSE
            LAMBDA(I,J)=LAMBDA(I,J)+VEC(IND)
            IF(I.NE.J)LAMBDA(J,I)=LAMBDA(J,I)+VEC(IND)
          END IF
        ENDDO
      ENDDO
      IND=IND+1
      Y=VEC(IND)
!
!     ==========================================================================
!     == FUDGE                                                                ==
!     ==========================================================================
      DO I=1,NCHI
        DO J=I+1,NCHI
          H(I,J)=0.5D0*H(I,J)    
          H(J,I)=0.5D0*H(J,I)
        ENDDO
      ENDDO
      DO I=1,NPSI
        DO J=I+1,NPSI
          LAMBDA(I,J)=0.5D0*LAMBDA(I,J)
          LAMBDA(J,I)=0.5D0*LAMBDA(J,I)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == ADD CONSTRAINT FORCES TO PSIP                                        ==
!     ==========================================================================
      N=MAXVAL(PSI0(:)%N)
      CALL CI$NEWPSI(PSI1,N)
      CALL CI$NEWPSI(PSI2,N)
      DO N=1,NPSI
        DO J=1,NCHI
          CALL CI$COPYPSI(PSI0(N),PSI1)
          CALL CI$ANNIHILATOR(PSI1,J)
          DO I=1,NCHI
            CALL CI$COPYPSI(PSI1,PSI2)
            CALL CI$CREATOR(PSI2,I)
            CALL CI$SCALARPRODUCT(PSI0(N),PSI2,CSVAR)
            XP(N)=XP(N)+2.D0*QX*X0(N)*REAL(CSVAR)*H(I,J)
            CALL CI$SCALEPSI(PSI2,H(I,J)*QPSI(N))
            CALL CI$ADDPSI(PSIP(N),PSI2)
          ENDDO
        ENDDO
      ENDDO

      DO N=1,NPSI
        DO M=1,NPSI
          CALL CI$COPYPSI(PSI0(M),PSI1)
          CALL CI$SCALEPSI(PSI1,LAMBDA(N,M)*QPSI(N)/(P0(N)+1.D-12))
          CALL CI$ADDPSI(PSIP(N),PSI1)
        ENDDO
      ENDDO
!
      XP(:)=XP(:)+2.D0*QX*X0(:)*Y
      CALL CI$DELETEPSI(PSI1)
      CALL CI$DELETEPSI(PSI2)
      RETURN
      END SUBROUTINE CI_ADDCONSTRAINT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$DYNWITHFIXEDDENMAT(NCHI,RHO,HAMILTON,NPSI,P0,PSI0,H,E)
!     **************************************************************************
!     **  DETERMINE ENSEMBLE OF MANY PARTICLE WAVE FUNCTION CONSISTENT WITH   **
!     **  THE ONE-PARTICLE DENSITY MATRIX RHO, WHICH MINIMIZES THE TOTAL      **
!     **  ENERGY SUM P_I<PSI_I|HAMILTON|PSI_I>                                **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN) :: NCHI           ! #(1P ORBITALS)
      INTEGER(4)        ,INTENT(IN) :: NPSI           ! #(MP WAVE FUNCTIONS)
      COMPLEX(8)        ,INTENT(IN) :: RHO(NCHI,NCHI) ! 1P DENSITY MATRIX
      TYPE(CIHAMIL_TYPE),INTENT(INOUT):: HAMILTON       ! MANY-P HAMILTONIAN
      REAL(8)           ,INTENT(OUT):: P0(NPSI)       ! PROBABILITIES
      TYPE(CISTATE_TYPE),INTENT(OUT):: PSI0(NPSI)     ! CURRENT WAVE FUNCTION
      COMPLEX(8)        ,INTENT(OUT):: H(NCHI,NCHI)   ! DE/DRHO2
      REAL(8)           ,INTENT(OUT):: E              ! ENERGY
      INTEGER(4)        ,PARAMETER  :: NOITER=20000     ! X#(ITERATIONS)
      REAL(8)                       :: DELTA=1.D-2  ! TIME STEP IN LOOP
      REAL(8)           ,PARAMETER  :: M0=1.D0        ! MASS 
      REAL(8)           ,PARAMETER  :: ANNE0=1.D-2
      REAL(8)                       :: MPSI(NPSI)
      REAL(8)                       :: ANNEPSI(NPSI)
      REAL(8)                       :: QPSI(NPSI)
      REAL(8)           ,PARAMETER  :: MX=1.D+2      ! MASS PROBABILITIES
      REAL(8)           ,PARAMETER  :: ANNEX=1.D-2     ! FRICTION (1.D-1)
      REAL(8)           ,PARAMETER  :: QX=M0*(1+ANNE0)/(MX*(1+ANNEX))
      REAL(8)           ,PARAMETER  :: MAXEKIN=0.02D0 ! X(KINETIC ENERGY)
      REAL(8)           ,PARAMETER  :: MAXHPSI=1.D-5  ! X(FOR HPSI CONTRIBUTION)
      REAL(8)           ,PARAMETER  :: MAXPSI=1.D-5   ! X(FOR PSI CONTRIBUTION)
      REAL(8)           ,PARAMETER  :: TARGETEKINDOT=1.D-2 !X(FOR PSI CONTRIB.)
      INTEGER(4)        ,PARAMETER  :: NEXPANDBASIS=10 ! 
      TYPE(CISTATE_TYPE)            :: PSIM(NPSI)     ! PREVIOUS WAVE FUNCTION
      TYPE(CISTATE_TYPE)            :: PSIP(NPSI)     ! NEXT WAVE FUNCTION
      TYPE(CISTATE_TYPE)            :: HPSI(NPSI)     ! H|PSI>
      COMPLEX(8)                    :: RHO2(NCHI,NCHI) ! DIAGONAL DENSITY MATRIX
      REAL(8)                       :: F(NCHI)        ! OCCUPATIONS
      COMPLEX(8)                    :: U(NCHI,NCHI)   ! NATURAL ORBITALS
      COMPLEX(8)                    :: LAMBDA(NPSI,NPSI) ! ENERGY
      REAL(8)                       :: X0(NPSI),XM(NPSI),XP(NPSI)
      REAL(8)                       :: EPOT,EPOTI(NPSI),EPOTLAST
      REAL(8)                       :: EKIN,EKINI(NPSI)
      REAL(8)                       :: EKINP
      REAL(8)                       :: FX
      INTEGER(4)                    :: ITER
      INTEGER(4)                    :: I,J,N
      LOGICAL(4)                    :: CONVG
      REAL(8)          ,PARAMETER   :: TOL=1.D-8
      INTEGER(4)       ,PARAMETER   :: NWAIT=20
      LOGICAL(4)       ,PARAMETER   :: TNATURAL=.FALSE.
      INTEGER(4)                    :: IWAIT
      REAL(8)                       :: EWAIT
      REAL(8)                       :: SVAR,SVAR1,SVAR2,SVAR3
      INTEGER(4)                    :: IERR
      INTEGER(4)                    :: NC
      LOGICAL(4)                    :: TCYCLED=.FALSE.
      LOGICAL(4)                    :: TSTOP
      INTEGER(4)                    :: ITERSTOP
      REAL(8) :: VEC(NCHI**2+NPSI**2+1)
      COMPLEX(8)                   :: CSVAR
!     **************************************************************************
      NC=NCHI**2+NPSI**2+1
      DO I=1,NPSI
        CALL CI$NEWPSI(PSIM(I),200)
        CALL CI$NEWPSI(PSIP(I),200)
        CALL CI$NEWPSI(HPSI(I),200)
      ENDDO

      OPEN(123, FILE="CIINFO")
      REWIND 123
!
!     ==========================================================================
!     == SET INDIVIDUAL MASSES AND FRICTIONS                                  ==
!     ==========================================================================
      MPSI(:)=M0
      ANNEPSI(:)=ANNE0
      QPSI(:)=M0*(1+ANNE0)/(MPSI(:)*(1.D0+ANNEPSI(:)))
!
!     ==========================================================================
!     == CHECK IF THE DENSITY MATRIX IS HERMITEAN                             ==
!     ==========================================================================
      DO I=1,NCHI
        DO J=I,NCHI
          IF(ABS(RHO(I,J)-CONJG(RHO(J,I))).GT.1.D-8) THEN
            CALL ERROR$MSG('DENSITY MATRIX NOT HERMITEAN. NOT ALLOWED...')
            CALL ERROR$STOP('CI$DYNWITHFIXEDDENMAT')
          END IF
        END DO
      ENDDO
!
!     ==========================================================================
!     == DIAGONALIZE DENSITY MATRIX                                           ==
!     ==========================================================================
      IF(TNATURAL) THEN
        CALL LIB$DIAGC8(NCHI,RHO,F,U)   ! SUM_K RHO(I,K)*U(K,N)= U(I,N)*EIG(N)
!       == REORDER OCCUPATIONS IN DECENDING ORDER ================================
        DO N=1,NCHI/2
          SVAR=F(N)
          F(N)=F(NCHI+1-N)
          F(NCHI+1-N)=SVAR
        ENDDO
        RHO2=U  ! RHO2 IS HERE ONLY A PLACE HOLDER
        DO N=1,NCHI/2+1
          U(:,N)       =RHO2(:,NCHI+1-N)
          U(:,NCHI+1-N)=RHO2(:,N)
        ENDDO
!       ==  MAKE DENSITY MATRIX N-REPRESENTABLE ================================
        DO N=1,NCHI
          F(N)=MIN(1.D0,MAX(0.D0,F(N)))
        ENDDO
!       ==  CONSTRUCT DIAGONAL DENSITY MATRIX ==================================
        RHO2(:,:)=(0.D0,0.D0)
        DO I=1,NCHI
          RHO2(I,I)=F(NCHI+1-I)
        ENDDO
!       == TRANSFORM HAMILTONIAN INTO THE NEW BASIS ============================
!       == THIS IS VERY TIME CONSUMING
PRINT*,'TRANSFORM HAMILTONIAN TO NATURAL ORBITALS....'
        CALL CI$HAMILTONTRANSFORM(NCHI,U,HAMILTON)      
PRINT*,'..........................TRANSFORMATION DONE'
      ELSE
        RHO2=RHO
        U(:,:)=(0.D0,0.D0)
        DO N=1,NCHI
          U(N,N)=(1.D0,0.D0)
        ENDDO
      END IF
!
!     ==========================================================================
!     == CONSTRUCT STARTING WAVE FUNCTION                                     ==
!     ==========================================================================
      CALL CI$STATEFROMDENSITYMATRIX(NCHI,NPSI,RHO2,P0,PSI0)
      DO I=1,NPSI
        CALL CI$LIMITSIZE(MAXPSI,PSI0(I)) ! REDUCE NUMBER OF SLATER DETERMINANTS
        CALL CI$COPYPSI(PSI0(I),PSIM(I))  ! ZERO INITIAL VELOCITY 
      ENDDO
      X0(:)=SQRT(P0(:))
      XM(:)=X0(:)
! HIER UNTERSCHEIDET SICH CHRISTIANS VERSION STARK VON MEINER. PB
!
!     ==========================================================================
!     == ITERATE TO FIND GROUND STATE                                         ==
!     ==========================================================================
      EPOTLAST=HUGE(EPOTLAST)
      EPOT=HUGE(EPOT)
      EKIN=0.D0
      EKINI(:)=0.D0
      IWAIT=0
      EWAIT=EPOT
      TCYCLED=.FALSE.
      ITERSTOP=0
      DO ITER=1,NOITER
        TSTOP=.FALSE.
        P0(:)=X0(:)**2
        EPOTLAST=EPOT
        EPOT=0.D0
        DO I=1,NPSI
          CALL CI$CLEANPSI(PSI0(I))
          CALL CI$COPYPSI(PSI0(I),HPSI(I))
          CALL CI$HPSI(HAMILTON,HPSI(I))     ! CONSTRUCT H|PSI(0)>
!          CALL CI$LIMITSIZE(MAXHPSI,HPSI(I)) ! REDUCE #(SLATER DETERMINANTS)
          CALL CI$EPOT(PSI0(I),HPSI(I),EPOTI(I))
          EPOT=EPOT+P0(I)*EPOTI(I)
        ENDDO
!
!       == SET FRICTION VALUE ==================================================
!!$PRINT*,'EPOT ',EPOT,EPOTLAST
!!$        IF(EPOT.GT.EPOTLAST+1.D-6) THEN
!!$          DO I=1,NPSI
!!$            CALL CI$LIMITSIZE(MAXPSI,PSIM(I))!REDUCE NUM. OF SLATER DETERMINANTS
!!$            CALL CI$COPYPSI(PSIM(I),PSI0(I))
!!$          ENDDO
!!$          PRINT*,'PSI0%N AFTER LIMITSIZE',PSI0(:)%N
!!$          EKIN=1.D-10
!!$          EPOT=EPOTLAST
!!$          CYCLE
!!$        END IF

!
!       == PROPAGATE WITHOUT CONSTRAINTS =======================================
        DO I=1,NPSI
          CALL CI$PROPAGATE(DELTA,ANNEPSI(I),MPSI(I),PSI0(I),PSIM(I),HPSI(I),PSIP(I))
        ENDDO
!
!       ========================================================================
!       == PROPAGATE PROBABILITIES                                            ==
!       ========================================================================
        SVAR1=2.D0/(1.D0+ANNEX)
        SVAR2=1.D0-SVAR1
        SVAR3=DELTA**2/MX/(1.D0+ANNEX)
        DO N=1,NPSI
          FX=-2.D0*X0(N)*(EPOTI(N)-EKINI(N)) !KINETIC ENERGY IGNORED
          XP(N)=X0(N)*SVAR1+XM(N)*SVAR2+FX*SVAR3
        ENDDO         
!
!       ========================================================================
!       == APPLY CONSTRAINTS                                                  ==
!       ========================================================================
!       == APPLY CONSTRAINTS ===================================================
        CALL CI_LAGRANGE(NCHI,NPSI,X0,PSI0,XP,PSIP,QX,QPSI,RHO2,H,LAMBDA,IERR)
!!$PRINT*,'++++++++++++++++++++++++++'
!!$DO I=1,NPSI
!!$  CALL CI$WRITEPSI(PSIP(I),6)
!!$  CALL CI$SCALARPRODUCT(PSIP(I),PSIP(I),CSVAR)
!!$  PRINT*,'NORM ',CSVAR
!!$ENDDO
!!$        IF(IERR.EQ.1) THEN
!!$          DO N=1,NPSI
!!$            CALL CI$COPYPSI(PSI0(N),PSIM(N))
!!$            X0(N)=XM(N)
!!$            EKINI(N)=0.D0
!!$          ENDDO
!!$          DELTA=0.5D0*DELTA  
!!$          PRINT*,'CYCLED TWICE'
!!$          IF(TCYCLED) EXIT
!!$          TCYCLED=.TRUE.
!!$          CYCLE
!!$        END IF
!!$        TCYCLED=.FALSE.
!
!       ========================================================================
!       == KINETIC ENERGY AND ENERGY REPORT                                   ==
!       ========================================================================
        EKIN=0.D0
        DO I=1,NPSI
          CALL CI$EKIN(DELTA,MPSI(I),PSIP(I),PSIM(I),EKINI(I))
          EKIN=EKIN+P0(I)*EKINI(I)
        ENDDO
        EKINP=0.5D0*MX*SUM((XP(:)-XM(:))**2)/(2.D0*DELTA)**2

!       == SET VELOCITY TO ZERO IF KINETIC ENERGY EXCEEDS LIMIT ================
        TSTOP=TSTOP.OR.(EKIN+EKINP.GT.MAXEKIN)
!        TSTOP=TSTOP.OR.(EPOT.GT.EPOTLAST)
        IF(TSTOP.AND.ITER.GT.ITERSTOP+50) THEN
          ITERSTOP=ITER
          DO I=1,NPSI
            CALL CI$COPYPSI(PSIM(I),PSI0(I))
            XM(I)=X0(I)
          ENDDO
          CYCLE
        END IF
!!$!
!!$        IF(EKIN.GT.MAXEKIN) THEN
!!$!         == KEEP KINETIC ENERGY ABOUT CONSTANT IF EKIN>MAXEKIN ================
!!$          ANNEPSI=-(EPOT-EPOTLAST)/(4.D0*MAXEKIN)
!!$        ELSE
!!$!         == CHOOSE FRICTION SO THAT EKIN GROWS WITH RATE TARGETEKINDOT
!!$          ANNEPSI=-(TARGETEKINDOT+EPOT-EPOTLAST)/(4.D0*EKIN)
!!$          ANNEPSI=MAX(ANNEPSI,-0.1D0)   ! AVOID TOO LARGE ACCELERATIONS (INSTABLE)
!!$          ANNEPSI=MIN(ANNEPSI,0.D0)     ! DO NOT SLOW DOWN 
!!$        END IF
!
!       ========================================================================
!       == LIMIT THE GROWTH OF THE NUMBER OF SLATER DETERMINANTS              ==
!       == THE MAXIMUM ALLOWED VALUE IS NEXPANDBASIS*NCHI                     ==
!       ========================================================================
PRINT*,'PSI0%N BEFORE GROWTH CHECK',PSIP%N
        DO I=1,NPSI
          CALL CI$GROWTHCHECK(NEXPANDBASIS*NCHI,PSIP(I),PSI0(I))
        ENDDO
PRINT*,'PSI0%N AFTER GROWTH CHECK',PSIP%N
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
!PRINT*,'-----',CONVG,IWAIT,EWAIT,EPOT
        IF(CONVG) EXIT
!
CALL CI_CONSTRAINT(NCHI,NPSI,X0,PSI0,RHO2,NC,VEC)
        WRITE(*,'("!>",I5,4F30.15,20F10.6)')ITER,EKINP,EKIN,EPOT,(EKIN+EPOT),SUM(ABS(VEC)),X0**2
        WRITE(123,'(I5,5F40.20)') ITER,EPOT,EKIN+EPOT,EKINP+EKIN+EPOT
!
!       ========================================================================
!       == SWITCH TO NEXT TIME STEP                                           ==
!       ========================================================================
        DO I=1,NPSI
          CALL CI$COPYPSI(PSI0(I),PSIM(I)) !COPIES TO THE RIGHT
          CALL CI$COPYPSI(PSIP(I),PSI0(I)) 
        ENDDO
        XM(:)=X0(:)
        X0(:)=XP(:)
      END DO
PRINT*,'LOOP OVER'
PRINT*,'-----',CONVG,IWAIT,EWAIT,EPOT
      CLOSE(123)
      IF(.NOT.CONVG) THEN
!!$        CALL ERROR$MSG('SELF-CONSISTENCY LOOP NOT CONVERGED')
!!$        CALL ERROR$STOP('CI$DYNWITHFIXEDDENMAT')
      END IF
      DO I=1,NPSI
        CALL CI$CLEANPSI(PSI0(I))
      ENDDO
!
!     ==========================================================================
!     == TRANSFORM LAGRANGE MULTIPLIERS BACK TO ORIGINAL BASIS                ==
!     ==========================================================================
      H=MATMUL(CONJG(TRANSPOSE(U)),MATMUL(H,U))    
!
!     ==========================================================================
!     == CONSTRUCT CONSTRAINING FORCES                                        ==
!     ==========================================================================
      SVAR=DELTA**2/(M0*(1.D0+ANNE0))
      H=-H/SVAR
      LAMBDA=TRANSPOSE(LAMBDA)/SVAR
      E=EPOT    ! E=<PSI|HAMILTON|PSI>
      DO I=1,NPSI
        CALL CI$DELETEPSI(PSIM(I))
        CALL CI$DELETEPSI(PSIP(I))
        CALL CI$DELETEPSI(HPSI(I))
      ENDDO
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
      INTEGER                          :: ID0,IDP
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI$TESTCI()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NSITE=8
      INTEGER(4),PARAMETER :: NCHI=2*NSITE
      INTEGER(4),PARAMETER :: NPSI=2
      REAL(8)   ,PARAMETER :: NPARTICLE=8.2D0 !REAL(NSITE)+0.2D0
      COMPLEX(8),PARAMETER :: U=(1.D0,0.D0)
      COMPLEX(8),PARAMETER :: T=-1.D0
      REAL(8)              :: H(NCHI,NCHI)
      REAL(8)              :: EI(NCHI)
      COMPLEX(8)           :: RHO2(NCHI,NCHI)
      COMPLEX(8)           :: V(NCHI,NCHI)
      TYPE(CIHAMIL_TYPE)   :: HAM          ! MANY-P HAMILTONIAN
      REAL(8)              :: P(NPSI)
      TYPE(CISTATE_TYPE)   :: PSI(NPSI)
      REAL(8)              :: E
      REAL(8)              :: ORB(NCHI,NCHI)
      REAL(8)              :: F(NCHI)
      INTEGER(4)           :: I,J,K,ISITE,ISPIN
      REAL(8)              :: SVAR       ! SUPPORT VARIABLE
      REAL(8)              :: DELTA
!     **************************************************************************
                              CALL TRACE$PUSH('CI$TESTCI')
!     ==========================================================================
!     == CONSTRUCT NATURAL ORBITALS FOR A FINITE HUBBARD CHAIN =================
!     ==========================================================================
      H(:,:)=0.D0
      DO ISPIN=1,2
        DO ISITE=1,NSITE-1
          I=2*(ISITE-1)+ISPIN
          J=2*ISITE+ISPIN
          H(I,J)=T
          H(J,I)=T
        ENDDO
        DO ISITE=1,NSITE
          I=2*(ISITE-1)+ISPIN
          H(I,I)=1.D-3*REAL(2*ISPIN-3,KIND=8)
        ENDDO
      ENDDO
      CALL LIB$DIAGR8(NCHI,H,EI,ORB)
      DO I=1,NCHI
        WRITE(*,FMT='("E",F10.3," ORB",20F10.3)')EI(I),ORB(:,I)
      ENDDO
!     ==========================================================================
!     == OCCUPATIONS ===========================================================
!     ==========================================================================
      SVAR=NPARTICLE
      F(:)=0.D0
      DO I=1,NCHI
        IF(SVAR.GT.0.D0) THEN
          F(I)=MIN(1.D0,SVAR)
          SVAR=SVAR-F(I)
        END IF
      ENDDO    
!     == RANDOMIZE OCCUPATIONS A LITTLE
      DO K=1,30
        DELTA=1.D-2
        CALL RANDOM_NUMBER(SVAR)
        SVAR=0.005D0+0.99D0*SVAR
        I=NINT(0.5D0+REAL(NCHI)*SVAR)
        IF(F(I)-DELTA.LT.0.D0) CYCLE
        CALL RANDOM_NUMBER(SVAR)
        SVAR=0.005D0+0.99D0*SVAR
        J=NINT(0.5D0+REAL(NCHI)*SVAR)
        IF(F(J)+DELTA.GT.1.D0) CYCLE
        F(I)=F(I)-DELTA
        F(J)=F(J)+DELTA
      ENDDO
      WRITE(*,FMT='("N=",F10.5," F=",20F10.5)')SUM(F),F
!
!     ==========================================================================
!     == CONSTRUCT DENSITY MATRIX FROM NATURAL ORBITALS ========================
!     ==========================================================================
      DO I=1,NCHI
        DO J=1,NCHI
          SVAR=0.D0
          DO K=1,NCHI
            SVAR=SVAR+ORB(I,K)*F(K)*ORB(J,K)
          ENDDO
          RHO2(I,J)=CMPLX(SVAR,KIND=8)
        ENDDO
      ENDDO
      WRITE(*,*)
      DO I=1,NCHI
        WRITE(*,FMT='(" RHO",20F10.3)')REAL(RHO2(:,I))
      ENDDO
!
!     ==========================================================================
!     == SET UP HAMILTONIAN ====================================================
!     ==========================================================================
      DO ISITE=1,NSITE
        I=2*ISITE-1
        J=2*ISITE
        CALL CI$SETU(HAM,I,J,I,J,U)
        CALL CI$SETU(HAM,J,I,J,I,U)
      ENDDO
      CALL CI$WRITEHAMILTONIAN(HAM,6)
!
!     ==========================================================================
!     == RUN DENSITY MATRIX FUNCTIONAL =========================================
!     ==========================================================================
      DO I=1,NPSI
        CALL CI$NEWPSI(PSI(I),200)
      ENDDO
PRINT*,'STARTING DENSITY-MATRIX FUNCTIONAL'
      CALL CI$DYNWITHFIXEDDENMAT(NCHI,RHO2,HAM,NPSI,P,PSI,V,E)
!
!     ==========================================================================
!     == REPORT RESULT                                                        ==
!     ==========================================================================
      DO I=1,NPSI
        CALL CI$LIMITSIZE(1.D-2,PSI(I)) ! REDUCE NUMBER OF SLATER DETERMINANTS
        CALL CI$WRITEPSI(PSI(I),6)
      ENDDO
STOP 'FORCED'
                                        CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CI_LAGRANGEMATOLD(NCHI,PSI0,NC,MAT)
!     **************************************************************************
!     ** PSI0 IS INTENT(INOUT) BECAUSE IT MAY BE CLEANED UP
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)   :: NCHI
      TYPE(CISTATE_TYPE),INTENT(INOUT):: PSI0
      INTEGER(4)        ,INTENT(IN)   :: NC
      REAL(8)           ,INTENT(OUT)  :: MAT(NC,NC)
      TYPE(CISTATE_TYPE)              :: PSI1,PSI2,PSI3,PSI4
      COMPLEX(8)                      :: CMAT4(NCHI,NCHI,NCHI,NCHI)
      COMPLEX(8)                      :: CMAT2(NCHI,NCHI)
      COMPLEX(8)                      :: CMAT0
      REAL(8)                         :: RMAT2(NCHI,NCHI)
      REAL(8)                         :: RMAT0
      INTEGER(4)                      :: IINDEX(2,NC)
      INTEGER(4)                      :: I,J,K,L,IC1,IC2
      COMPLEX(8)                      :: CSVAR1,CSVAR2
      COMPLEX(8)        ,PARAMETER    :: CI=(0.D0,1.D0)
!     **************************************************************************
      IF(NC.NE.NCHI**2+1) THEN
        CALL ERROR$STOP('CI_LAGRANGEMAT')
      END IF
      CALL CI$NEWPSI(PSI1,PSI0%N)
      CALL CI$NEWPSI(PSI2,PSI0%N)
      CALL CI$NEWPSI(PSI3,PSI0%N/2)
      CALL CI$NEWPSI(PSI4,PSI0%N/4)
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
