!*******************************************************************************
!**                                                                           **
!**  NAME: SPHERICAL                                                          **
!**                                                                           **
!**  PURPOSE: REAL SPHERICAL HARMONICS AND CLEBSCH GORDON COEFFICIENTS        **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    GETYLM                                                                 **
!**    ROTATEYLM                                                              **
!**    CLEBSCH                                                                **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    1) FUNCTIONS NEED TO BE PUT INTO THE FORM SPHERICAL$FUNCTION           **
!**    2) INTERNAL ROUTINES SHALL BE INCLUDED INTO A MODULE                   **
!**                                                                           **
!*******************************************************************************
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SPHERICAL_MODULE
INTEGER(4)          :: LMXX=0
REAL(8),ALLOCATABLE :: GAUNTMAT(:,:,:)
END MODULE SPHERICAL_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ROTATEYLM(LMX,ROT,YLMROT)
!     **************************************************************************
!     ** OLD INTERFACE FOR SPHERICAL$ROTATEYLM                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(IN) :: ROT(3,3)
      REAL(8)   ,INTENT(OUT):: YLMROT(LMX,LMX)
!     **************************************************************************
      CALL SPHERICAL$ROTATEYLM(LMX,ROT,YLMROT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GETYLM(LMX,R,YLM)
!     **************************************************************************
!     ** OLD INTERFACE FOR SPHERICAL$YLM                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: LMX
      REAL(8),   INTENT(IN)   :: R(3)
      REAL(8),   INTENT(OUT)  :: YLM(LMX)
!     **************************************************************************
      CALL SPHERICAL$YLM(LMX,R,YLM)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLEBSCH(LM1,LM2,LM3,CG)
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM1
      INTEGER(4),INTENT(IN) :: LM2
      INTEGER(4),INTENT(IN) :: LM3
      REAL(8)   ,INTENT(OUT):: CG
!     **************************************************************************
      CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$YLM(LMX,R,YLM)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATE REAL SPHERICAL HARMONICS YLM AT POINT R                   **
!     **  THE SPHERICAL HARMONICS ARE ORDERD WITH INDICES                     **
!     **        LM(L,M)=1+L**2+L-M                                            **
!     **                                                                      **
!     **  THE FIRST NINE REAL SPHERICAL HARMONICS ARE:                        **
!     **      YLM(1)=SQRT( 1/( 4*PI))    * 1                                  **
!     **      YLM(2)=SQRT( 3/( 4*PI))    * X / R                              **
!     **      YLM(3)=SQRT( 3/( 4*PI))    * Z / R                              **
!     **      YLM(4)=SQRT( 3/( 4*PI))    * Y / R                              **
!     **      == D-TYPE SPHERICAL HARMONICS 
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2              **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2              **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2              **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2              **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2              **
!     **      == F-TYPE REAL HARMONICS (PLEASE CHECK!)                        **
!     **      YLM(10)=SQRT( 35/(32*PI))  * ( (X^3-3*X*Y^2 ) /R**3             **
!     **      YLM(11)=SQRT(210/(32*PI))  * ( Z*(X^2-Y^2)    /R**3             **
!     **      YLM(12)=SQRT( 21/(32*PI))  * ( X*(5Z^2-R^2)   /R**3             **
!     **      YLM(13)=SQRT( 14/(32*PI))  * ( Z*(5Z^2-3*R^2) /R**3             **
!     **      YLM(14)=SQRT( 21/(32*PI))  * ( Y*(5Z^2-R^2)   /R**3             **
!     **      YLM(15)=SQRT(840/(32*PI))  * ( X*Y*Z          /R**3             **
!     **      YLM(16)=SQRT( 35/(32*PI))  * ( Y*(3X^2-Y^2)   /R**3             **
!     ** 
!     **    REMARK: THE F-TYPE ORBITALS IN THIS COMMENT HAVE BEEN TAKEN FROM  **
!     **    HTTPS://EN.WIKIPEDIA.ORG/WIKI/TABLE_OF_SPHERICAL_HARMONICS.       **
!     **    THE CODE WORKS INDEPENDENTLY.                                     **
!     **                                                                      **
!     **                                         P.E. BLOECHL, (1991,2018)    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: LMX
      REAL(8),   INTENT(IN)   :: R(3)
      REAL(8),   INTENT(OUT)  :: YLM(LMX)
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      COMPLEX(8)              :: EIPHI,EIMPHI
      REAL(8)                 :: FPI,SQ2
      INTEGER(4)              :: LX
      REAL(8)                 :: DIS,DISXY
      REAL(8)                 :: COSTHE,SINPHI,COSPHI
      INTEGER(4)              :: L,LM0,M,LM
      REAL(8)                 :: FAC,SVAR
      INTEGER(4)              :: LMM,LMP
!     **************************************************************************
      FPI=4.D0*PI
      SQ2=SQRT(2.D0)
      IF(LMX.LE.1) THEN
        IF(LMX.LE.0) RETURN
        YLM(1)=1.D0/SQRT(FPI)
        RETURN
      END IF
      LX=INT(SQRT(REAL(LMX))-1.D0)
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('GETYLM')
      END IF
      DISXY=SQRT(R(1)**2+R(2)**2)
      DIS=SQRT(DISXY**2+R(3)**2)
      IF(ABS(DIS).LT.1.D-12) THEN
        YLM(1)=1.D0/SQRT(FPI)
        DO LM=2,LMX
          YLM(LM)=0.D0
        ENDDO
        RETURN
      ENDIF
      COSTHE=R(3)/DIS
!     == SPHERICAL_PLGNDR RETURNS THE ASSOCIATED LEGENDRE POLYNOMIALS ==========
!     == MULTIPLIED WITH  SQRT( (L-M)! / (L+M)! ) FOR LM(L,M)=L(L+1)+1-M =======
      CALL SPHERICAL_PLGNDR(LMX,LX,COSTHE,YLM)
      IF(DISXY.NE.0.D0) THEN
        SINPHI=R(2)/DISXY
        COSPHI=R(1)/DISXY
        EIPHI=CMPLX(COSPHI,SINPHI,KIND=8)
      ELSE
        EIPHI=(1.D0,0.D0)
      END IF
      DO L=0,LX
        LM0=L*(L+1)+1
        FAC=SQRT(DBLE(2*L+1)/FPI)
        YLM(LM0)=YLM(LM0)*FAC
        EIMPHI=(1.D0,0.D0)
        DO M=1,L
          EIMPHI=EIMPHI*EIPHI
          LMM=LM0+M
          LMP=LM0-M
          SVAR=FAC*YLM(LMM)*SQ2
          YLM(LMP)=SVAR*REAL(EIMPHI,KIND=8)
          YLM(LMM)=SVAR*AIMAG(EIMPHI)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$YLMPOLYNOMIALS(LX,INDX,LMX,YLMPOL)
!     **************************************************************************
!     ** REAL SPHERICAL HARMONICS (TIMES R^L) AS POWER SERIES                 **
!     **                                                                      **
!     ** REFERS TO OBJECT PAW_GAUSSIAN REGARDING THE ORDERING OF POWERS       **
!     **                                                                      **
!     ********************************PETER BLOECHL, GOSLAR 2011****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LX     ! MAX(ANGULAR MOMENTUM QUANTUM NUMBER)
      INTEGER(4),INTENT(IN) :: INDX   ! 
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(OUT):: YLMPOL(INDX,LMX) 
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: FACTORIAL(0:2*LX)
      INTEGER(4)            :: L,MPOS,IX,JY,KZ,P,I,J,K,J1,K1,N,IND,LM,M
      REAL(8)               :: SVAR1,SVAR2,SVAR3,SVAR4,SVAR5,RES
      INTEGER(4)            :: ISVAR
      LOGICAL(4),PARAMETER  :: TWRITE=.FALSE.
!     **************************************************************************
!
!     ==========================================================================
!     == CHECK CONSISTENCY OF INDICES                                         ==
!     ==========================================================================
      IF(LX.LT.0) THEN
        CALL ERROR$MSG('INDEX LX MUST NOT BE NEGATIVE')
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('SPHERICAL$YLMPOLYNOMIALS')
      END IF
      IF(LMX.LT.(LX+1)**2) THEN
        CALL ERROR$MSG('INDEX LMX TOO SMALL FOR SPECIFIED MAX ANGULAR MOMENTUM')
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('REQUIRED MINIMUM LMX',(LX+1)**2)
        CALL ERROR$STOP('SPHERICAL$YLMPOLYNOMIALS')
      END IF
      IF(INDX.LT.(LX+1)*(LX+2)*(LX+3)/6) THEN
        CALL ERROR$MSG('INDEX INDX TOO SMALL FOR SPECIFIED MAX ANG. MOMENTUM')
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$I4VAL('INDX',INDX)
        CALL ERROR$I4VAL('REQUIRED MINIMUM INDX',(LX+1)*(LX+2)*(LX+3)/6)
        CALL ERROR$STOP('SPHERICAL$YLMPOLYNOMIALS')
      END IF
!
!     ==========================================================================
!     == PREPARATION: PI AND FACTORIAL                                        ==
!     ==========================================================================
      ISVAR=0
      FACTORIAL(0)=1.D0
      DO I=1,2*LX
        FACTORIAL(I)=FACTORIAL(I-1)*REAL(I,KIND=8)
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT SPHERICAL HARMONICS                                        ==
!     ==========================================================================
      YLMPOL(:,:)=0.D0
      DO L=0,LX
        DO MPOS=0,L
          SVAR1=SQRT(REAL(2*L+1)/(4.D0*PI)*FACTORIAL(L-MPOS)/FACTORIAL(L+MPOS))
          SVAR1=SVAR1*(-1.D0)**MPOS/(2.D0**L*FACTORIAL(L))  !FROM ASS. LEG. POL.
          DO P=0,MPOS
            CALL BINOMIALCOEFFICIENT(MPOS,P,RES)
            SVAR2=SVAR1*RES
!           == CONVERT TO REAL SPHERICAL HARMONICS =============================
            SVAR2=SVAR2*(-1.D0)**MPOS
            IF(MPOS.NE.0)SVAR2=SVAR2*SQRT(2.D0)  
            IF(MODULO(P,2).EQ.0) THEN    ! P IS EVEN => M=+MPOS
              LM=L**2+L+1-MPOS
              ISVAR=P/2
              SVAR2=SVAR2*(-1.D0)**ISVAR ! SIGN RESULTING FROM I^P
            ELSE                         ! P IS ODD => M=-MPOS 
              LM=L**2+L+1+MPOS 
              ISVAR=(P-1)/2              
              SVAR2=SVAR2*(-1.D0)**ISVAR ! SIGN RESULTING FROM I^P
            END IF
!
            J1=INT(0.5D0*REAL(L+MPOS+1))
            DO J=J1,L
              CALL BINOMIALCOEFFICIENT(L,J,RES)
              SVAR3=SVAR2*FACTORIAL(2*J)/FACTORIAL(2*J-L-MPOS)*RES
              IF(MODULO(L-J,2).EQ.1) SVAR3=-SVAR3   ! (-1)^(L-J)
              DO K1=0,L-J
                CALL BINOMIALCOEFFICIENT(L-J,K1,RES)
                SVAR4=SVAR3*RES
                DO N=0,K1
                  CALL BINOMIALCOEFFICIENT(K1,N,RES)
                  SVAR5=SVAR4*RES
                  IX=2*N+MPOS-P
                  JY=2*K1-2*N+P
                  KZ=L-MPOS-2*K1
                  CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND,IX,JY,KZ)
                  YLMPOL(IND,LM)=YLMPOL(IND,LM)+SVAR5
                ENDDO
              ENDDO
            ENDDO
          ENDDO              
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WRITE POWER SERIES EXPANSION                                         ==
!     ==========================================================================
      IF(TWRITE) THEN
        LM=0
        DO L=0,LX
          DO M=-L,L
            LM=LM+1
            WRITE(*,*)'================ L=',L,' M=',M,' LM=',LM,'=============='
            DO IND=1,INDX
              IF(YLMPOL(IND,LM).EQ.0.D0) CYCLE
              CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
              WRITE(*,*)I,J,K,YLMPOL(IND,LM)
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$LMBYNAME(NAME,LM)
!     **************************************************************************
!     ** RETURNS THE COMBINED ANGULAR MOMENTUM INDEX FOR A NAME               **
!     ** OF REAL SPHERICAL HARMONICS                                          **
!     ** 
!     **************************************************************************
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(OUT):: LM
!     **************************************************************************
      IF(NAME.EQ.'S') THEN 
        LM=1
      ELSE IF(NAME.EQ.'PX') THEN
        LM=2
      ELSE IF(NAME.EQ.'PZ') THEN
        LM=3
      ELSE IF(NAME.EQ.'PY') THEN
        LM=4
      ELSE IF(NAME.EQ.'DX2-Y2') THEN
        LM=5
      ELSE IF(NAME.EQ.'DXZ') THEN
        LM=6
      ELSE IF(NAME.EQ.'D3Z2-R2'.OR.NAME.EQ.'DZ2') THEN
        LM=7
      ELSE IF(NAME.EQ.'DYZ') THEN
        LM=8
      ELSE IF(NAME.EQ.'DXY') THEN
        LM=9
      ELSE IF(NAME.EQ.'FX(X2-3Y2)'.OR.NAME.EQ.'F1') THEN
        LM=10
      ELSE IF(NAME.EQ.'FZ(X2-Y2)'.OR.NAME.EQ.'F2') THEN
        LM=11
      ELSE IF(NAME.EQ.'FX(5Z2-R2)'.OR.NAME.EQ.'FXZ2'.OR.NAME.EQ.'F3') THEN
        LM=12
      ELSE IF(NAME.EQ.'FZ(5Z2-3R2)'.OR.NAME.EQ.'FZ3'.OR.NAME.EQ.'F4') THEN
        LM=13
      ELSE IF(NAME.EQ.'FY(5Z2-R2)'.OR.NAME.EQ.'FYZ2'.OR.NAME.EQ.'F5') THEN
        LM=14
      ELSE IF(NAME.EQ.'FXYZ'.OR.NAME.EQ.'F6') THEN
        LM=15
      ELSE IF(NAME.EQ.'FY(3X2-Y2)'.OR.NAME.EQ.'F7') THEN
        LM=16
      ELSE
        CALL ERROR$MSG('ORBITAL NAME NOT RECOGNIZED')
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('SPHERICAL$LMBYNAME')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$YLMNAME(LM,NAME)
!     **************************************************************************
!     ** RETURNS THE NAME FOR A REAL HARMONICS FOR A GIVEN INDEX              **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: LM
      CHARACTER(*),INTENT(OUT) :: NAME
      CHARACTER(8)             :: STRING
      INTEGER(4)   ,PARAMETER  :: LMX=16
      CHARACTER(11),PARAMETER  :: NAMES(LMX)= &
!     ** '12345678901','12345678901','12345678901','12345678901','12345678901'
     & (/'S          ' &
     &  ,'PX         ','PZ         ','PY         ' &
     &  ,'DX2-Y2     ','DXZ        ','D3Z2-R2    ','DYZ        ','DXY        ' &
     &  ,'FX(X2-3Y2) ','FZ(X2-Y2)  ','FX(5Z2-R2) ','FZ(5Z2-3R2)' &
     &                ,'FY(5Z2-R2) ','FXYZ       ','FY(3X2-Y2) '/)
!     **************************************************************************
      IF(LM.LE.0) THEN
        CALL ERROR$MSG('INDEX FOR CUBIC HARMONIC MUST BE POSITIVE')
        CALL ERROR$STOP('SPHERICAL$YLMNAME')
      ELSE IF(LM.GT.LMX) THEN
        WRITE(STRING,*)LM
        NAME='LM='//TRIM(STRING)
      ELSE
        NAME=TRIM(NAMES(LM))
      END IF       
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL_PLGNDR(LMX,LX,X,PLM)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATE THE ASSOCIATED LEGENDRE POLYNOMIALS MULTIPLIED WITH       **
!     **  SQRT((L-M)!/(L+M)!).                                                **
!     **                                                                      **
!     **  REMARK                                                              **
!     **    1) LM(L,M)=L(L+1)-M+1                                             **
!     **    2) THE ROUTINE PLAYS WITH                                         **
!     **                                                                      **
!     **  THE ASSOCIATED LEGENDRE POLYNOMIALS ARE DEFINED  FOR NON-NEGATIVE M **
!     **    P_{L,M}=(-1)^M/(2^L*L!) (1-X^2)^{M/2) (D/DX)^{M+L} (1-X^2)^L      **
!     **  AND FOR NEGATIVE M AS                                               **
!     **    P_{L,-M}=(-1)^M  (L-M)!/(L+M)! P_{L,M}                            **
!     **                                                                      **
!     **                                                                      **
!     **  P_{L,L}=(-1)^L (2L-1)!! (1-X^2)^{L/2}                               **
!     **                                                                      **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: LMX
      INTEGER(4)  ,INTENT(IN)  :: LX
      REAL(8)     ,INTENT(IN)  :: X
      REAL(8)     ,INTENT(OUT) :: PLM(LMX)
      REAL(8)                  :: FACT,FAC2
      REAL(8)                  :: SVAR
      INTEGER(4)               :: LM,L,M,I1,LMM,LMP
      REAL(8)                  :: PMM,CM,C0,CP
!     **************************************************************************
      FACT=1.D0-X**2
      FACT=MAX(FACT,0.D0)  ! AVOID NUMERICAL PROBLEM THAT 1-1**2 < 0
      FACT=SQRT(FACT)
      IF(FACT.EQ.0.D0) THEN
        SVAR=-X
        LM=0
        DO L=0,LX
          SVAR=-SVAR
          DO M=1,2*L+1
            LM=LM+1
            PLM(LM)=SVAR
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == DETERMINE ASSOCIATED LEGENDRE POLYNOMIALS FOR M.GE.0                 ==
!     == USE P_{M,M}(X)=(-1)^M  (2*M-1)!!  [SQRT(1-X^2)]^M                    ==
!     == AND THE RECURRENCE RELATION                                          ==
!     ==    (L-M)*P_{L,M}(X)=X*(2L-1)*P_{L-1,M}(X)-(L+M-1)*P_{L-2,M}(X)       ==
!     == THE INDICES ARE MAPPED ACCORDING TO                                  ==
!     ==    LM(L,M)=L(L+1)-M+1                                                ==
!     ==========================================================================
      PMM=1.D0
      DO M=0,LX
        CM=0.D0
        C0=PMM
        LM=M**2+1   ! LM=1,2,5,10   
        PLM(LM)=C0
        DO L=M+1,LX
          LM=LM+2*L    
!         == RECURRENCE RELATION OF ASSOCIATED LEGENDRE POLYNOMIALS ============
          CP=(X*DBLE(2*L-1)*C0-DBLE(L+M-1)*CM)/DBLE(L-M)
          PLM(LM)=CP
          CM=C0
          C0=CP
        ENDDO
!       ==  DETERMINE P_{M,M}(X) ITERATIVELY ===================================
        PMM=-DBLE(2*M+1)*FACT*PMM
      ENDDO
!
!     ==========================================================================
!     ==  MULTIPLY FACTOR SQRT[(L+M)!/(L-M)!]
!     ==  AND COMPLETE FOR M<0                                                ==
!     ==========================================================================
      DO L=1,LX
        I1=L*(L+1)+1
        FACT=1.D0
        FAC2=1.D0
        DO M=1,L
          LMM=I1+M
          LMP=I1-M
          FACT=FACT*SQRT(DBLE((L-M+1)*(L+M)))  ! FACT=SQRT[(L+M)!/(L-M)!]
          FAC2=-FAC2                           ! FAC2=(-1)^M
          PLM(LMP)=PLM(LMP)/FACT
          PLM(LMM)=FAC2*PLM(LMP)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$AXISUROT(DR,U,DU)
!     **************************************************************************
!     ** DETERMINES A 3X3 ROTATION MATRIX U, WHICH TRANSFORMS THE VECTOR DR   **
!     ** INTO THE Z-DIRECTION, AND ITS DERIVATIVE DU                          **
!     **     E_Z|DR|=U*DR, WHERE E_Z IS THE UNIT VECTOR IN Z-DIRECTION        **
!     **     DU(I,J,K)=DU_IJ/DR_K                                             **
!     **                                                                      **
!     ** REMARK: THE MATRIX U CHANGES DISCONTINUOUSLY, WHEN THE TWO ABSOLUTE  **
!     **   SMALLEST COMPONENTS OF DR INTERCHANGE. THIS CAN CONFUSE TESTS      **
!     **   BASED ON NUMERICAL DIFFERENCES                                     **
!     **************************************PETER BLOECHL, GOSLAR 2020**********
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: DR(3)      ! VECTOR WILL BE ROTATED IN Z-DIRECTION
      REAL(8),INTENT(OUT) :: U(3,3)     ! UNITARY ROTATION MATRIX U(I,J)
      REAL(8),INTENT(OUT) :: DU(3,3,3)  ! DERIVATIVE DU(I,J)/DR(K)
      REAL(8)             :: AVEC(3)
      REAL(8)             :: Y(3,3)
      REAL(8)             :: Y1L,Y2L,Y3L
      REAL(8)             :: DYDR(3,3,3)
      REAL(8)             :: ADR
      REAL(8)             :: DR2
      REAL(8)             :: SVAR
      INTEGER(4)          :: I,J,N
!     **************************************************************************
!     ==========================================================================
!     == CHOOSE VECTOR AVEC, WHICH IS NOT LINEAR DEPENDENT WITH DR            ==
!     ==========================================================================
      AVEC(:)=0.D0
      AVEC(MINLOC(ABS(DR),1))=1.D0
      ADR=SUM(AVEC*DR)
      DR2=SUM(DR*DR)
!
!     ==========================================================================
!     == DETERMINE ORTHOGONAL SET OF VECTORS WITH Y3 LINEAR DEPENDENT WITH DR ==
!     ==========================================================================
      Y(:,3)=DR(:)
      Y(1,2)=AVEC(2)*DR(3)-AVEC(3)*DR(2)
      Y(2,2)=AVEC(3)*DR(1)-AVEC(1)*DR(3)
      Y(3,2)=AVEC(1)*DR(2)-AVEC(2)*DR(1)
      Y(:,1)=DR(:)*ADR-AVEC(:)*DR2   ! Y2 TIMES Y3 (SEE BAC-CAB RULE)
!
!     ==========================================================================
!     == DETERMINE GRADIENT OF THE VECTORS Y                                  ==
!     ==========================================================================
      DYDR(:,:,:)=0.D0
      DO I=1,3
        DYDR(I,3,I)=1.D0
      ENDDO
      DYDR(1,2,3)= AVEC(2)
      DYDR(1,2,2)=-AVEC(3)
      DYDR(2,2,1)= AVEC(3)
      DYDR(2,2,3)=-AVEC(1)
      DYDR(3,2,2)= AVEC(1)
      DYDR(3,2,1)=-AVEC(2)
      DO I=1,3
        DYDR(:,1,I)=DR(:)*AVEC(I)-AVEC(:)*2.D0*DR(I)
        DYDR(I,1,I)=DYDR(I,1,I)+ADR
      ENDDO
!
!     ==========================================================================
!     == NORMALIZE Y-VECTORS                                                  ==
!     ==========================================================================
      Y1L=SQRT(SUM(Y(:,1)**2))
      Y2L=SQRT(SUM(Y(:,2)**2))
      Y3L=SQRT(SUM(Y(:,3)**2))
      Y(:,1)=Y(:,1)/Y1L
      Y(:,2)=Y(:,2)/Y2L
      Y(:,3)=Y(:,3)/Y3L
      DYDR(:,1,:)=DYDR(:,1,:)/Y1L
      DYDR(:,2,:)=DYDR(:,2,:)/Y2L
      DYDR(:,3,:)=DYDR(:,3,:)/Y3L
      DO N=1,3
        DO J=1,3
          SVAR=SUM(Y(:,J)*DYDR(:,J,N))
          DYDR(:,J,N)=DYDR(:,J,N)-Y(:,J)*SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == ROTATION MATRIX AND GRADIENT                                         ==
!     == ROTATION MATRIX U IS THE TRANSPOSE OF THE VECTOR-MATRIX Y            ==
!     ==========================================================================
      U=TRANSPOSE(Y)
      DO N=1,3
        DU(:,:,N)=TRANSPOSE(DYDR(:,:,N))
      ENDDO

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$ROTATEYLM(LMX,ROT,YLMROT)
!     **************************************************************************
!     **  PRODUCES A LMX*LMX MATRIX YLMROT                                    **
!     **  THAT TRANSFORMS A COEFFICIENT VECTOR C_LM                           **
!     **  OF A FUNCTION EXPRESSED AS F(R)=SUM_LM Y_LM(R)*C_LM                 **
!     **  INTO A COORDINATE SYSTEM RPRIME=ROT*R                               **
!     **  WITH YPRIME_LM(RPRIME)=Y_LM(R)                                      **
!     **  SUCH THAT    F(R)=SUM_LM YPRIME_LM(R)*CPRIME_LM                     **
!     **  WITH CPRIME_LM1=SUM_LM2 YLMROT_LM1,LM2 C_LM2                        **
!     **  REMAINS INVARIANT UNDER ROTATION                                    **
!     **                                                                      **
!     **  WORKS ONLY FOR REAL HARMONICS  WITH Y_2=C*X;Y_3=C*Z;Y_4=C*Y         **
!     **  USES THE CLEBSCH GORDAN ROUTINE CLEBSCH                             **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX             ! #(L,M-ANGULAR MOMENTA)
      REAL(8)   ,INTENT(IN) :: ROT(3,3)        ! ROTATION MATRIX FOR POSITIONS
      REAL(8)   ,INTENT(OUT):: YLMROT(LMX,LMX) !ROTATION MATRIX FOR COEFFICIENTS
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.     ! PRINTS RESULT
      LOGICAL(4),PARAMETER  :: TTEST=.TRUE.    ! CHECKS WHETHER RO IS UNITARY
      INTEGER(4)            :: LX
      INTEGER(4)            :: I1,I2,J1,J2
      INTEGER(4)            :: INDEX(3)
      INTEGER(4)            :: LM1A,LM1B,LM2A,LM2B
      INTEGER(4)            :: LM1,LM2,LM3,LM4,LM5,LM6,L
      REAL(8)               :: SVAR,SVAR1,SVAR2,CG
!     **************************************************************************
!
!      =========================================================================
!      == TEST WHETHER ROTATION MATRIX IS UNITARY                             ==
!      =========================================================================
       IF(TTEST) THEN
         SVAR=ROT(1,1)*(ROT(2,2)*ROT(3,3)-ROT(2,3)*ROT(3,2)) &
      &      +ROT(1,2)*(ROT(2,3)*ROT(3,1)-ROT(2,1)*ROT(3,3)) &
      &      +ROT(1,3)*(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))
         IF(ABS(SVAR-1.D0).GT.1.D-12) THEN
           CALL ERROR$MSG('ROTATION MATRIX NOT UNITARY')
           CALL ERROR$STOP('SPHERICAL$ROTATYLM')
         END IF
       END IF
!
!      =========================================================================
!      == INITIALIZE YLMROT FOR L=0 AND L=1                                   ==
!      =========================================================================
       LX=INT(SQRT(REAL(LMX))-1.D0)
       IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
         CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
         CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
         CALL ERROR$I4VAL('LMX',LMX)
         CALL ERROR$I4VAL('LX',LX)
         CALL ERROR$STOP('SPHERICAL$ROTATYLM')
       END IF
       YLMROT(:,:)=(0.D0,0.D0)
       IF(LX.GE.0)YLMROT(1,1)=1.D0
       IF(LX.GE.1) THEN
         INDEX(1)=2
         INDEX(2)=4
         INDEX(3)=3
         DO I1=1,3
           I2=INDEX(I1)
           DO J1=1,3
             J2=INDEX(J1)
             YLMROT(I2,J2)=ROT(I1,J1)
           ENDDO
         ENDDO
       END IF
!
!      =========================================================================
!      == APPLY RECURSION FOR YLMROT FOR L>1                                  ==
!      =========================================================================
       DO L=2,LX
         LM1A=L**2+1
         LM1B=(L+1)**2
         LM2A=(L-1)**2+1
         LM2B=L**2
         SVAR2=0.D0
         DO LM3=2,4
           DO LM5=LM2A,LM2B
             CALL CLEBSCH(LM1A,LM3,LM5,CG)
             SVAR2=SVAR2+CG**2
           ENDDO
         ENDDO
         
         DO LM1=LM1A,LM1B
           DO LM2=LM1A,LM1B
             SVAR=0.D0
             DO LM3=2,4
               DO LM5=LM2A,LM2B
                 SVAR1=0.D0
                 DO LM4=2,4
                   DO LM6=LM2A,LM2B
                     CALL CLEBSCH(LM2,LM4,LM6,CG)
                     SVAR1=SVAR1+YLMROT(LM3,LM4)*YLMROT(LM5,LM6)*CG
                   ENDDO
                 ENDDO
                 CALL CLEBSCH(LM1,LM3,LM5,CG)
                 SVAR=SVAR+SVAR1*CG
               ENDDO
             ENDDO
             YLMROT(LM1,LM2)=SVAR/SVAR2
           ENDDO
         ENDDO 
       ENDDO
!
!      =================================================================
!      == PRINT RESULT FOR TEST PURPOSES                              ==
!      =================================================================
       IF(TPR) THEN
         DO L=0,LX
           LM1A=L**2+1
           LM1B=(L+1)**2
           DO LM1=LM1A,LM1B
             WRITE(*,FMT='(9F10.3)')YLMROT(LM1A:LM1B,LM1)
           ENDDO
         ENDDO
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$ROTATEYLMWDER(LMX,N,ROT,DROT,YLMROT,DYLMROT)
!     **************************************************************************
!     **  SIMILAR TO SPHERICAL$ROTATEYLM, BUT WITH GRADIENTS                  **
!     **                                                                      **
!     **  PRODUCES A LMX*LMX MATRIX YLMROT AND GRADIENTS                      **
!     **  THAT TRANSFORMS A COEFFICIENT VECTOR C_LM                           **
!     **  OF A FUNCTION EXPRESSED AS F(R)=SUM_LM Y_LM(R)*C_LM                 **
!     **  INTO A COORDINATE SYSTEM RPRIME=ROT*R                               **
!     **  WITH YPRIME_LM(RPRIME)=Y_LM(R)                                      **
!     **  SUCH THAT    F(R)=SUM_LM YPRIME_LM(R)*CPRIME_LM                     **
!     **  WITH CPRIME_LM1=SUM_LM2 YLMROT_LM1,LM2 C_LM2                        **
!     **  REMAINS INVARIANT UNDER ROTATION                                    **
!     **                                                                      **
!     **  WORKS ONLY FOR REAL HARMONICS  WITH Y_2=C*X;Y_3=C*Z;Y_4=C*Y         **
!     **  USES THE CLEBSCH GORDAN ROUTINE CLEBSCH                             **
!     **                                                                      **
!     **  DROT(LMX,LMX,N) DESCRIBES N DERIVATIVES OF the rotation matrix      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX             ! #(L,M-ANGULAR MOMENTA)
      INTEGER(4),INTENT(IN) :: N               ! #(DERIVATIVES)
      REAL(8)   ,INTENT(IN) :: ROT(3,3)        ! ROTATION MATRIX FOR POSITIONS
      REAL(8)   ,INTENT(IN) :: DROT(3,3,N)     ! PERTURBATION
      REAL(8)   ,INTENT(OUT):: YLMROT(LMX,LMX) !ROTATION MATRIX FOR COEFFICIENTS
      REAL(8)   ,INTENT(OUT):: DYLMROT(LMX,LMX,N) !PERTURBATION
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.     ! PRINTS RESULT
      LOGICAL(4),PARAMETER  :: TTEST=.TRUE.    ! CHECKS WHETHER ROt IS UNITARY
      INTEGER(4)            :: LX
      INTEGER(4)            :: I1,I2,J1,J2
      INTEGER(4)            :: INDEX(3)
      INTEGER(4)            :: LM1A,LM1B,LM2A,LM2B
      INTEGER(4)            :: LM1,LM2,LM3,LM4,LM5,LM6,L
      REAL(8)               :: SVAR,SVAR1,SVAR2,CG
      REAL(8)               :: DSVAR(N),DSVAR1(N)
!     **************************************************************************
!
!      =========================================================================
!      == TEST WHETHER ROTATION MATRIX IS UNITARY                             ==
!      =========================================================================
       IF(TTEST) THEN
         SVAR=ROT(1,1)*(ROT(2,2)*ROT(3,3)-ROT(2,3)*ROT(3,2)) &
      &      +ROT(1,2)*(ROT(2,3)*ROT(3,1)-ROT(2,1)*ROT(3,3)) &
      &      +ROT(1,3)*(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))
         IF(ABS(SVAR-1.D0).GT.1.D-12) THEN
           CALL ERROR$MSG('ROTATION MATRIX NOT UNITARY')
           CALL ERROR$STOP('SPHERICAL$ROTATYLMWDER')
         END IF
       END IF
!
!      =========================================================================
!      == INITIALIZE YLMROT FOR L=0 AND L=1                                   ==
!      =========================================================================
       LX=INT(SQRT(REAL(LMX))-1.D0)
       IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
         CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
         CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
         CALL ERROR$I4VAL('LMX',LMX)
         CALL ERROR$I4VAL('LX',LX)
         CALL ERROR$STOP('SPHERICAL$ROTATYLMWDER')
       END IF
       YLMROT(:,:)=(0.D0,0.D0)
       dYLMROT(:,:,:)=(0.D0,0.D0)
       IF(LX.GE.0) THEN
         YLMROT(1,1)=1.D0
         DYLMROT(1,1,:)=0.D0
       END IF
       IF(LX.GE.1) THEN
         INDEX(1)=2
         INDEX(2)=4
         INDEX(3)=3
         DO I1=1,3
           I2=INDEX(I1)
           DO J1=1,3
             J2=INDEX(J1)
             YLMROT(I2,J2)   =ROT(I1,J1)
             DYLMROT(I2,J2,:)=DROT(I1,J1,:)
           ENDDO
         ENDDO
       END IF
!
!      =========================================================================
!      == APPLY RECURSION FOR YLMROT FOR L>1                                  ==
!      =========================================================================
       DO L=2,LX
         LM1A=L**2+1
         LM1B=(L+1)**2
         LM2A=(L-1)**2+1
         LM2B=L**2
         SVAR2=0.D0
         DO LM3=2,4
           DO LM5=LM2A,LM2B
             CALL CLEBSCH(LM1A,LM3,LM5,CG)
             SVAR2=SVAR2+CG**2
           ENDDO
         ENDDO
         
         DO LM1=LM1A,LM1B
           DO LM2=LM1A,LM1B
             SVAR=0.D0
             dSVAR=0.D0
             DO LM3=2,4
               DO LM5=LM2A,LM2B
                 SVAR1=0.D0
                 DSVAR1=0.D0
                 DO LM4=2,4
                   DO LM6=LM2A,LM2B
                     CALL CLEBSCH(LM2,LM4,LM6,CG)
                     SVAR1=SVAR1+YLMROT(LM3,LM4)*YLMROT(LM5,LM6)*CG
                     DSVAR1=DSVAR1+(DYLMROT(LM3,LM4,:)*YLMROT(LM5,LM6) &
       &                           +YLMROT(LM3,LM4)*DYLMROT(LM5,LM6,:))*CG
                   ENDDO
                 ENDDO
                 CALL CLEBSCH(LM1,LM3,LM5,CG)
                 SVAR=SVAR+SVAR1*CG
                 DSVAR(:)=DSVAR(:)+DSVAR1(:)*CG
               ENDDO
             ENDDO
             YLMROT(LM1,LM2)=SVAR/SVAR2
             DYLMROT(LM1,LM2,:)=DSVAR(:)/SVAR2
           ENDDO
         ENDDO 
       ENDDO
!
!      =================================================================
!      == PRINT RESULT FOR TEST PURPOSES                              ==
!      =================================================================
       IF(TPR) THEN
         DO L=0,LX
           LM1A=L**2+1
           LM1B=(L+1)**2
           DO LM1=LM1A,LM1B
             WRITE(*,FMT='(9F10.3)')YLMROT(LM1A:LM1B,LM1)
           ENDDO
         ENDDO
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$YLMTRANS(LMX,C)
!     **                                                                      **
!     **  TRANSFORMATION MATRIX FROM REAL SPHERICAL HARMONICS YBAR            **
!     **  TO ANGULAR MOMENTUM EIGENSTATES Y                                   **
!     **                                                                      **
!     **   1=\SUM_{LM,LMPRIME} |Y_{LM}><Y_{LM}|YBAR_{LMPRIME}><YBAR_{LMPRIME}|**
!     **    =\SUM_{LM,LMPRIME} |Y_{LM}>C_{LM,LMPRIME}<YBAR_{LMPRIME}|         **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: C(LMX,LMX)
      INTEGER(4)            :: L,M,LMP,LMM,I
      INTEGER(4)            :: LMAX
      REAL(8)               :: FAC
      REAL(8)               :: SQRTINV
      LOGICAL(4)            :: TTEST=.FALSE.
      COMPLEX(8)            :: CTEST(LMX,LMX)
      REAL(8)               :: SVAR
!     **************************************************************************
      LMAX=INT(SQRT(REAL(LMX)+1.D-8))-1
      IF((LMAX+1)**2.NE.LMX) THEN
        CALL ERROR$MSG('LMX MUST SPAN FULL L-SHELLS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$STOP('SPHERICAL_YLMTRANS')
      END IF
      SQRTINV=1.D0/SQRT(2.D0)
      C=(0.D0,0.D0)
      DO L=0,LMAX
        LMP=L**2+L+1
        LMM=LMP
        C(LMP,LMP)=(1.D0,0.D0)
        FAC=1.D0
        DO M=1,L
          FAC=-FAC
          LMP=LMP+1
          LMM=LMM-1
          C(LMP,LMP)=    (1.D0, 0.D0)*SQRTINV
          C(LMP,LMM)=    (0.D0,-1.D0)*SQRTINV
          C(LMM,LMP)=FAC*(1.D0, 0.D0)*SQRTINV
          C(LMM,LMM)=FAC*(0.D0, 1.D0)*SQRTINV
        ENDDO
      ENDDO
!
      IF(TTEST) THEN
        CTEST=MATMUL(TRANSPOSE(CONJG(C)),C)
        DO I=1,LMX
          CTEST(I,I)=CTEST(I,I)-(1.D0,0.D0)
        ENDDO
        SVAR=MAXVAL(ABS(CTEST))
        IF(SVAR.GT.1.D-10) THEN
          CALL ERROR$MSG('TRANSFORMATION IS NOT UNITARY')   
          CALL ERROR$R8VAL('DEVIATION',SVAR)
          CALL ERROR$STOP('SPHERICAL_YLMTRANS')   
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$ER(LMX,X,Y,Z)
!     **************************************************************************
!     **  MATRIX ELEMENTS OF THE UNIT VECTORS X/R, Y/R, Z/Y                   **
!     **                                                                      **
!     **   XI/R=\SUM_{LM,LMPRIME} |Y_LM><Y_LM|XI/R|Y_LMPRIME><YBAR_LMPRIME|   **
!     **       =\SUM_{LM,LMPRIME} |Y_LM> XI_{LM,LMPRIME} <YBAR_LMPRIME|       **
!     **                                                                      **
!     **       X_{LM,LMPRIME}=<LM|X/R|LMPRIME>                                **
!     **       Y_{LM,LMPRIME}=<LM|Y/R|LMPRIME>                                **
!     **       Z_{LM,LMPRIME}=<LM|Z/R|LMPRIME>                                **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(OUT):: X(LMX,LMX)
      REAL(8)   ,INTENT(OUT):: Y(LMX,LMX)
      REAL(8)   ,INTENT(OUT):: Z(LMX,LMX)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: SQ4PIBY3
      INTEGER(4)            :: LM1,LM2
!     **************************************************************************
      SQ4PIBY3=SQRT(4.D0*PI/3.D0)
      DO LM1=1,LMX
        DO LM2=1,LMX
          CALL CLEBSCH(LM1,LM2,2,X(LM1,LM2))
          CALL CLEBSCH(LM1,LM2,3,Y(LM1,LM2))
          CALL CLEBSCH(LM1,LM2,4,Z(LM1,LM2))
        ENDDO
      ENDDO
      X=X*SQ4PIBY3
      Y=Y*SQ4PIBY3
      Z=Z*SQ4PIBY3
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$L(LMX,LX,LY,LZ)
!     **************************************************************************
!     **  ANGULAR MOMENTUM MATRIX ELEMENTS IN A REPRESENTATION                **
!     **  OF REAL SPHERICAL HARMONICS                                         **
!     **                                                                      **
!     **    LX_{LM,LMPRIME}=<LM|LX|LMPRIME>                                   **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: LX(LMX,LMX)
      COMPLEX(8),INTENT(OUT):: LY(LMX,LMX)
      COMPLEX(8),INTENT(OUT):: LZ(LMX,LMX)
      COMPLEX(8)            :: C(LMX,LMX)
      INTEGER(4)            :: LMAX
      INTEGER(4)            :: LM,L,M
      INTEGER(4)            :: LOX(LMX),MOX(LMX)
      REAL(8)               :: LPFAC(LMX),LMFAC(LMX),LZFAC(LMX)
      COMPLEX(8)            :: LPMAT(LMX,LMX),LMMAT(LMX,LMX)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
!     **************************************************************************
      LMAX=INT(SQRT(REAL(LMX)+1.D-8))-1
      IF((LMAX+1)**2.NE.LMX) THEN
        CALL ERROR$STOP('SPHERICAL_L')
      END IF
      LM=0
      DO L=0,LMAX
        DO M=-L,L
          LM=LM+1
          LOX(LM)=L
          MOX(LM)=M
          LPFAC(LM)=SQRT(REAL((L-M)*(L+M+1),KIND=8))
          LMFAC(LM)=SQRT(REAL((L+M)*(L-M+1),KIND=8))
          LZFAC(LM)=REAL(M,KIND=8)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == ANGULAR MOMENTA IN ANGULAR MOMENTUM EIGENSTATES                      ==
!     ==========================================================================
      LZ=(0.D0,0.D0)
      LPMAT=(0.D0,0.D0)
      LMMAT=(0.D0,0.D0)
      DO LM=1,LMX
        LZ(LM,LM)=CMPLX(LZFAC(LM),0.D0,KIND=8)
        IF(MOX(LM).NE.+LOX(LM))LPMAT(LM+1,LM)=CMPLX(LPFAC(LM),0.D0,KIND=8)
        IF(MOX(LM).NE.-LOX(LM))LMMAT(LM-1,LM)=CMPLX(LMFAC(LM),0.D0,KIND=8)
      ENDDO
!
!     ==========================================================================
!     == TRANSFORM FROM A BASIS OF ANGULAR MOMENTUM EIGENSTATES TO A BASIS    ==
!     == OF REAL SPHERICAL HARMONICS                                          ==
!     ==========================================================================
      CALL SPHERICAL$YLMTRANS(LMX,C)
      LZ=MATMUL(LZ,C)
      LPMAT=MATMUL(LPMAT,C)
      LMMAT=MATMUL(LMMAT,C)
      C=TRANSPOSE(CONJG(C))
      LZ=MATMUL(C,LZ)
      LPMAT=MATMUL(C,LPMAT)
      LMMAT=MATMUL(C,LMMAT)
!
!     ==========================================================================
!     == CARTESIAN ANGULAR MOMENTA                                            ==
!     ==========================================================================
      LX= 0.5D0   *(LPMAT+LMMAT)
      LY=-0.5D0*CI*(LPMAT-LMMAT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL_TEST()
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: LMX=9
      COMPLEX(8)           :: LX(LMX,LMX),LY(LMX,LMX),LZ(LMX,LMX)
      COMPLEX(8)           :: C(LMX,LMX)
      COMPLEX(8)           :: CDAGGER(LMX,LMX)
      COMPLEX(8)           :: MAT(LMX,LMX)
      INTEGER(4)           :: I,LM,L,M
      INTEGER(4)           :: MOX(LMX)
      INTEGER(4)           :: LMAX
      REAL(8)              :: SVAR
      COMPLEX(8)           :: CI=(0.D0,1.D0)
!     **************************************************************************
      LMAX=INT(SQRT(REAL(LMX)+1.D-8))-1
      IF((LMAX+1)**2.NE.LMX) THEN
        CALL ERROR$STOP('SPHERICAL_L')
      END IF
      LM=0
      DO L=0,LMAX
        DO M=-L,L
          LM=LM+1
          MOX(LM)=M
         ENDDO
      ENDDO
!
      CALL SPHERICAL$L(LMX,LX,LY,LZ)
      CALL SPHERICAL$YLMTRANS(LMX,C)
      CDAGGER=TRANSPOSE(CONJG(C))
!
!     == CHECK EIGENVALUE EQUATION FOR LZ
      MAT=MATMUL(LZ,CDAGGER)
      DO I=1,LMX
        MAT(:,I)=MAT(:,I)-CDAGGER(:,I)*REAL(MOX(I))
      ENDDO
      SVAR=MAXVAL(ABS(MAT))
      PRINT*,'DEVIATION EIGENVALUE LZ ',SVAR

!     == CHECK EIGENVALUE EQUATION FOR L**2
      MAT=MATMUL(LX,LX)+MATMUL(LY,LY)+MATMUL(LZ,LZ)
      PRINT*,'L**2'
      DO I=1,LMX
        WRITE(*,FMT='(20("(",2F6.2,")"))')MAT(I,:)
      ENDDO
!
!     == PRINT TRANSFORMATION FROM REAL HARMONICS
      PRINT*,'C'
      MAT=C
      DO I=1,LMX
        WRITE(*,FMT='(20("(",2F6.2,")"))')MAT(I,:)
      ENDDO
!
!     == CHECK IF C IS UNITARY
      PRINT*,'CDAGGER*C'
      MAT=MATMUL(CDAGGER,C)
      DO I=1,LMX
        WRITE(*,FMT='(20("(",2F6.2,")"))')MAT(:,I)
      ENDDO
      PRINT*,'C*CDAGGER'
      MAT=MATMUL(C,CDAGGER)
      DO I=1,LMX
        WRITE(*,FMT='(20("(",2F6.2,")"))')MAT(:,I)
      ENDDO
!
      PRINT*,'CHECK COMMUTATOR RELATIONS'
      MAT=MATMUL(LX,LY)-MATMUL(LY,LX)-CI*LZ
      PRINT*,'[LX,LY]-I*LZ',MAXVAL(ABS(MAT))
      MAT=MATMUL(LY,LZ)-MATMUL(LZ,LY)-CI*LX
      PRINT*,'[LY,LZ]-I*LX',MAXVAL(ABS(MAT))
      MAT=MATMUL(LZ,LX)-MATMUL(LX,LZ)-CI*LY
      PRINT*,'[LZ,LX]-I*LY',MAXVAL(ABS(MAT))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
!     **************************************************************************
!     **  CALCULATE GAUNT COEFFICIENT FOR REAL SPHERICAL HARMONICS.           **
!     **    Y(LM1)*Y(LM2)= SUM(LM3|CG(LM1,LM2,LM3)*Y(LM3))                    **
!     ****************************************** P.E. BLOECHL, 1991 ************
      USE SPHERICAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM1
      INTEGER(4),INTENT(IN) :: LM2
      INTEGER(4),INTENT(IN) :: LM3
      REAL(8)   ,INTENT(OUT):: CG
      INTEGER(4)            :: LM123X
      INTEGER(4)            :: LMAX
!     **************************************************************************
      LM123X=MAX(LM1,LM2,LM3)
      IF(LM123X.GT.LMXX) THEN
        LMAX=INT(SQRT(REAL(LM123X-1)+0.01))
        LMAX=MAX(4,LMAX)
        LMXX=(LMAX+1)**2
        IF(ALLOCATED(GAUNTMAT))DEALLOCATE(GAUNTMAT)
        ALLOCATE(GAUNTMAT(LMXX,LMXX,LMXX))
        CALL SPHERICAL_GAUNTMATRIX(LMXX,LMAX,GAUNTMAT)
      END IF
      CG=GAUNTMAT(LM1,LM2,LM3)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL_GAUNTMATRIX(LMXX,LMAX,CG)
!     **************************************************************************
!     **  CALCULATE MATRIX OF GAUNT COEFFICIENTS FOR REAL SPHERICAL HARMONICS **
!     **    Y(LM1)*Y(LM2)= SUM(LM3|CG(LM1,LM2,LM3)*Y(LM3))                    **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      INTEGER(4),INTENT(IN) :: LMXX
      INTEGER(4),INTENT(IN) :: LMAX
      REAL(8)   ,INTENT(OUT):: CG(LMXX,LMXX,LMXX)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: FACT(0:170) !UPPER LIMIT IS MAX IN REAL(8)
      REAL(8)               :: SQFACT(0:170)
      REAL(8)               :: SQPI,SQ2
      INTEGER(4)            :: I,IMAX
      INTEGER(4)            :: LM1,L1,M1,N1,IS1
      INTEGER(4)            :: LM2,L2,M2,N2,IS2
      INTEGER(4)            :: LM3,L3,M3,N3,IS3
      INTEGER(4)            :: L3MIN,L3MAX
      REAL(8)               :: FAC0,FAC1,FAC2,SVAR
      INTEGER(4)            :: LMUP
      INTERFACE
        DOUBLE PRECISION FUNCTION THREEJ(L1,M1,L2,M2,L3,M3,FACT,SQFACT)
        INTEGER(4),INTENT(IN) :: L1,M1,L2,M2,L3,M3
        REAL(8)   ,INTENT(IN) :: FACT(0:170),SQFACT(0:170)
        END
      END INTERFACE
      INTERFACE
        DOUBLE PRECISION FUNCTION PRECG(L1,L2,L3,FACT,SQFACT)
        INTEGER(4),INTENT(IN) :: L1,L2,L3
        REAL(8)   ,INTENT(IN) :: FACT(0:170),SQFACT(0:170)
        END
      END INTERFACE
!     **************************************************************************
      SQ2=SQRT(2.D0)
      SQPI=SQRT(PI)
      IMAX=4*LMAX+1
      IF(IMAX.GT.150) THEN
!       == PREVISIOUSLY IT WAS TESTED FOR IMAX.GE.50 ===========================
!       == POSSIBLY THE RESULTS ARE INACCURATE BEYOND THIS VALUE ===============
        CALL ERROR$MSG('DIMENSIONS IN CLBSCH TO SMALL')
        CALL ERROR$MSG('FACTORIAL PRODUCES OVERFLOW FOR I>170')
        CALL ERROR$I4VAL('FACTORIAL REQUESTED UP TO',IMAX)
        CALL ERROR$I4VAL('MAX L REQUESTED',LMAX)
        CALL ERROR$STOP('CLBSCH')
      END IF
      FACT(0)=1.D0
      SQFACT(0)=1.D0
      DO I=1,IMAX
        FACT(I)=FACT(I-1)*DBLE(I)
        SQFACT(I)=SQRT(FACT(I))
      ENDDO
!
!     ==========================================================================
!     ==  INITIALIZE ARRAY                                                    ==
!     ==========================================================================
      DO LM3=1,LMXX
        DO LM2=1,LMXX
          DO LM1=1,LMXX
            CG(LM1,LM2,LM3)=0.D0
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE CLEBSCH GORDAN COEFFICIENTS                               ==
!     ==========================================================================
      DO L1=0,LMAX
        DO M1=-L1,L1
          N1=IABS(M1)
          IS1=0
          IF(M1.LT.0) IS1=1
!       ==
          DO L2=0,L1
            DO M2=-L2,L2
              N2=IABS(M2)
              IS2=0
              IF(M2.LT.0) IS2=1
!             ==
              L3MIN=IABS(L2-L1)
              L3MAX=MIN0(LMAX,L1+L2)
              DO L3=L3MIN,L3MAX,2
                DO M3=-L3,L3
                  N3=IABS(M3)
                  IS3=0
                  IF(M3.LT.0) IS3=1
!                 ==
                  IF(M1*M2.LT.0.AND.(M3.NE.-IABS(N1+N2) &
     &                         .AND. M3.NE.-IABS(N1-N2))) GOTO 210
                  IF(M1*M2.EQ.0.AND. M3.NE.M1+M2        ) GOTO 210
                  IF(M1*M2.GT.0.AND.(M3.NE.+IABS(N1+N2) &
     &                         .AND. M3.NE.+IABS(N1-N2))) GOTO 210
!                 ==
                  FAC0=PRECG(L1,L2,L3,FACT,SQFACT)
                  FAC1=0.D0
               IF(N1.EQ.0.AND.N2.EQ.0)FAC1=FAC1+FAC0/SQRT(DBLE(2*L3+1))
                  FAC1=FAC1+THREEJ(L1, N1,L2, N2,L3, N3,FACT,SQFACT) &
     &                     *(-1.D0)**(N3+IS3) &
     &                     +THREEJ(L1, N1,L2,-N2,L3,-N3,FACT,SQFACT) &
     &                     *(-1.D0)**(N2+IS2) &
     &                     +THREEJ(L1,-N1,L2, N2,L3,-N3,FACT,SQFACT) &
     &                     *(-1.D0)**(N1+IS1)
                  FAC2=SQRT( 0.25D0*DBLE((2*L1+1)*(2*L2+1)) ) / SQ2 &
     &                          *(-1.D0)**(N3+(IS1+IS2+IS3)/2)    !WARNING CAN IS1+IS2+IS3 BE EVEN?
                  IF(M1.EQ.0)FAC2=FAC2/SQ2
                  IF(M2.EQ.0)FAC2=FAC2/SQ2
                  IF(M3.EQ.0)FAC2=FAC2/SQ2
                  LM1=(L1+1)**2-L1-M1
                  LM2=(L2+1)**2-L2-M2
                  LM3=(L3+1)**2-L3-M3
                  CG(LM1,LM2,LM3)=FAC0*FAC1*FAC2/SQPI
                  CG(LM2,LM1,LM3)=CG(LM1,LM2,LM3)
210               CONTINUE
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      IF(.NOT.TPR) RETURN
!     ==========================================================================
!     == PRINTOUT FOR CHECK                                                   ==
!     ==========================================================================
      LMUP=(LMAX+1)**2
      DO LM1=1,LMUP
        DO LM2=1,LMUP
          DO LM3=1,LMUP
            IF(ABS(CG(LM1,LM2,LM3)).GT.1.D-6) THEN
              SVAR=CG(LM1,LM2,LM3)*SQRT(4.D0*PI)
              WRITE(*,6000)LM1,LM2,LM3,CG(LM1,LM2,LM3),SVAR
6000          FORMAT(' CG ',3I3,2F10.5)
            END IF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      DOUBLE PRECISION FUNCTION THREEJ(L1,M1,L2,M2,L3,M3,FACT,SQFACT)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATE WIEGNER 3J SYMBOLS                                        **
!     **                                                                      **
!     **  THREEJ = (-1)**L1-L2-M3 * ( L1  L2  L3 )                            **
!     **                            ( M1  M2  M3 )                            **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L1
      INTEGER(4),INTENT(IN) :: M1
      INTEGER(4),INTENT(IN) :: L2
      INTEGER(4),INTENT(IN) :: M2
      INTEGER(4),INTENT(IN) :: L3
      INTEGER(4),INTENT(IN) :: M3
      REAL(8)   ,INTENT(IN) :: FACT(0:170)
      REAL(8)   ,INTENT(IN) :: SQFACT(0:170)
      REAL(8)               :: R1,R2
      INTEGER(4)            :: KMAX,KMIN,K
      REAL(8)               :: SUM,SVAR
!     **************************************************************************
      IF(M3.NE.M1+M2) THEN
        THREEJ=0.D0
        RETURN
      END IF
!
!     ==========================================================================
!     ==  CALCULATE THE 3J SYMBOLS  TIMES ....                                ==
!     ==========================================================================
      R1=SQFACT(L1+L2-L3)*SQFACT(L3+L1-L2)*SQFACT(L3+L2-L1) &
     &  /SQFACT(L1+L2+L3+1)
      R2=SQFACT(L1+M1)*SQFACT(L1-M1)*SQFACT(L2+M2)*SQFACT(L2-M2) &
     &                              *SQFACT(L3+M3)*SQFACT(L3-M3)
      KMAX=MIN0(L1+L2-L3,L1-M1,L2+M2)
      KMIN=MAX0(0,L2-L3-M1,L1-L3+M2)
      SUM=0.D0
      DO K=KMIN,KMAX
        SVAR=FACT(K)*FACT(L1+L2-L3-K)*FACT(L1-M1-K)*FACT(L2+M2-K) &
     &      *FACT(L3-L2+M1+K)*FACT(L3-L1-M2+K)
        SUM=SUM+(-1.D0)**K/SVAR
      ENDDO
      THREEJ=R1*R2*SUM
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      DOUBLE PRECISION FUNCTION PRECG(L1,L2,L3,FACT,SQFACT)
!     **************************************************************************
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L1
      INTEGER(4),INTENT(IN) :: L2
      INTEGER(4),INTENT(IN) :: L3
      REAL(8)   ,INTENT(IN) :: FACT(0:170)
      REAL(8)   ,INTENT(IN) :: SQFACT(0:170)
      INTEGER(4)            :: LT,LTH
!     **************************************************************************
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      LT=L1+L2+L3
      LTH=LT/2
      IF(2*LTH.NE.LT) THEN
        PRECG=0.D0
      ELSE
        PRECG=SQRT(DBLE(2*L3+1)/DBLE(LT+1))*FACT(LTH)/SQFACT(LT) &
     &       *SQFACT(LT-2*L1)/FACT(LTH-L1) &
     &       *SQFACT(LT-2*L2)/FACT(LTH-L2) &
     &       *SQFACT(LT-2*L3)/FACT(LTH-L3) &
     &       *(-1.D0)**(LTH-L3)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$TESTGAUNTREL()
!     **************************************************************************
!     ** THIS ROUTINE HAS NO SPECIAL USE OTHER THAN TO TEST A MATHEMATICAL    **
!     ** IDENTITY                                                             **
!     ** ROUTINE TO TEST THE IDENTITY                                         **
!     **   4*PI*(2*L3+1) SUM_{M1,M2} C_{LM1,LM2,LM3A)C_(LM1,LM2,LM3B)         **
!     **                                   =Q(L1,L2,L3)*DELTA(LM3A,LM3B)      **
!     **  WITH                                                                **
!     **    Q(L1,L2,L3)=Q(L2,L1,L3)                                           **
!     **    Q(L1,L2,L3)=Q(L2,L3,L1)=Q(L3,L1,L2)                               **
!     **    Q(L1,L2,L3)=0 IF L1+L2+L3 IS ODD                                  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: LX1=10
      INTEGER(4),PARAMETER :: LX2=10
      INTEGER(4),PARAMETER :: LX3=10
      REAL(8)   ,PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8) ,allocatable :: COEFF(:,:,:) !(LX1+1,LX2+1,(LX3+1)**2)
      INTEGER(4)           :: L1,L2,L3,L3A,L3B
      INTEGER(4)           :: IM1,IM2,IM3,IM3A,IM3B
      INTEGER(4)           :: LM1,LM2,LM3,LM3A,LM3B
      REAL(8)              :: CG1,CG2,SVAR
!     **************************************************************************
      allocate(coeff(LX1+1,LX2+1,(LX3+1)**2))
!
!     ==========================================================================
!     == TAKE CARE OF ANGULAR PART                                            ==
!     ==========================================================================
PRINT*,'========================================================='
      COEFF(:,:,:)=0.D0
      DO L1=0,LX1
        DO L2=0,LX2
          DO L3A=0,LX3
            DO IM3A=1,2*L3A+1
              LM3A=L3A**2+IM3A
              DO L3B=0,LX3
                DO IM3B=1,2*L3B+1
                  LM3B=L3B**2+IM3B
!
                  SVAR=0.D0
                  DO IM1=1,2*L1+1
                    LM1=L1**2+IM1
                    DO IM2=1,2*L2+1
                      LM2=L2**2+IM2
                      CALL SPHERICAL$GAUNT(LM1,LM2,LM3A,CG1)
                      CALL SPHERICAL$GAUNT(LM1,LM2,LM3B,CG2)
                      SVAR=SVAR+4.D0*PI*CG1*CG2/REAL((2*L1+1)*(2*L2+1))
                    ENDDO
                  ENDDO
!
                  IF(LM3A.EQ.LM3B) THEN
                    COEFF(L1+1,L2+1,LM3A)=SVAR
                  ELSE
                    IF(ABS(SVAR).GT.1.D-10) THEN
                      WRITE(*,FMT='(4I5,F20.10)')L1,L2,LM3A,LM3B,SVAR
                      CALL ERROR$MSG('TEST 1 FAILED ')
                      CALL ERROR$MSG('NONZERO RESULT FOR LM3A.NE.LM3B')
                      CALL ERROR$STOP('SPHERICAL$TESTGAUNTREL1')
                    END IF
                    CYCLE
                  END IF
!
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  TEST DEPENDENCE ON M3                                               ==
!     ==========================================================================
      DO L3=0,LX3
        DO IM3=2,2*L3+1
          SVAR=MAXVAL(ABS(COEFF(:,:,L3**2+IM3)-COEFF(:,:,L3**2+1)))
          IF(SVAR.GT.1.D-8) THEN
            DO L2=0,LX2
              DO L1=0,LX1
                SVAR=ABS(COEFF(L1+1,L2+1,L3**2+IM3)-COEFF(L1+1,L2,L3**2+1))
                IF(SVAR.LT.1.D-8) CYCLE
                WRITE(*,FMT='(3I5,20F15.9)')L1,L2,L3 &
     &                      ,COEFF(L1+1,L2+1,L3**2+1:(L3+1)**2)
              ENDDO
            ENDDO
            CALL ERROR$MSG('TEST 2 FAILED ')
            CALL ERROR$MSG('RESULT DEPENDS ON M3')
            CALL ERROR$I4VAL('L3',L3)
            CALL ERROR$I4VAL('IM3',IM3)
            CALL ERROR$R8VAL('DEVIATION',SVAR)
            CALL ERROR$STOP('SPHERICAL$TESTGAUNTREL1')
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  TEST IF RESULT VANISHES FOR ODD L1+L2+L3                            ==
!     ==========================================================================
      DO L1=0,LX1
        DO L2=0,LX2
          DO L3=0,LX3
            IF(MOD(L1+L2+L3,2).EQ.0) CYCLE
            LM3=L3**2  
            IM3=2*L3+1
            SVAR=MAXVAL(ABS(COEFF(L1+1,L2+1,LM3+1:LM3+IM3)))
            IF(SVAR.GT.1.D-8) THEN
              CALL ERROR$MSG('TEST 3 FAILED ')
              CALL ERROR$MSG('RESULT NONZERO FOR ODD L1+L2+L3')
              CALL ERROR$STOP('SPHERICAL$TESTGAUNTREL')
            END IF
          ENDDO
        ENDDO
      ENDDO
      STOP 'TESTGAUNTREL1 COMPLETED SUCESSFULLY'
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SPHERICALCCMAT_MODULE
!*******************************************************************************
!**                                                                           **
!**  USED TO EVALUATE THE MATRIXES A0 AND A1 TO EVALUATE STRESSES             **
!**  IN A SPHERICAL HARMONICS REPRESENTATION                                  **
!**                                                                           **
!**  GI(D/DGJ)SUM_L V_L Y_L                                                   **
!**        = G_IG_J SUM_L [1/G^2 DV_L/DG -L V_L] Y_L                          **
!**        + SUM_L SUM_LPRIME V_L [P0_IJ(L,LPRIME)+PM_IJ(L,LPRIME)] Y_LPRIME  **
!**                                                                           **
!**  WHERE P_0 IS NON-ZERO ONLY IF THE MAIN ANGULAR MOMENTA ARE               **
!**  IDENTICAL AND FOR P_M THE SECOND ANGULAR MOMENTUM IS LOWER               **
!**  BY 2 THAN THE FIRST                                                      **
!**                                                                           **
!**  USE THE INTERFACES                                                       **
!**    SPHERICAL$CCMAT0(LM1,I,LM2,CC)                                         **
!**    SPHERICAL$CCMATM(LM1,I,LM2,CC)                                         **
!**                                                                           **
!*******************************************************************************
TYPE A_TYPE 
LOGICAL(4):: EXIST0(10)
LOGICAL(4):: EXISTM(10)
REAL(8)   :: A0(3,3,10)
REAL(8)   :: AM(3,3,10)
INTEGER(4):: LM20(10)
INTEGER(4):: LM2M(10)
END TYPE A_TYPE 
INTEGER(4)               :: LRXX=0
INTEGER(4)               :: LMRXX=0
TYPE(A_TYPE),ALLOCATABLE :: THIS(:)
END MODULE SPHERICALCCMAT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$CLEARCCMAT
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
!     **************************************************************************
      IF(ALLOCATED(THIS))DEALLOCATE(THIS)
      LRXX=0
      LMRXX=0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICALCCMAT_INITIALIZE(LM_)
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN):: LM_
      REAL(8)              :: A1(3,3)
      REAL(8)              :: FOURPIBY3
      REAL(8)              :: DSMALL=1.D-12
      REAL(8)              :: CG1,CG2
      INTEGER(4),PARAMETER :: MAPLM(3)=(/2,4,3/)
      INTEGER(4)           :: I,J,LMI,LMJ,LM1,LM2,LM3,LM3X,L1,L2,LM
      LOGICAL(4)           :: TCHK,T0,TM
!     **************************************************************************
      IF(ALLOCATED(THIS))DEALLOCATE(THIS)
      LRXX=MAX(INT(SQRT(REAL(LM_-1,KIND=8)+DSMALL)),2)
      LMRXX=(LRXX+1)**2
      ALLOCATE(THIS(LMRXX))
      DO LM=1,LMRXX
        THIS(LM)%A0(:,:,:)=0.D0
        THIS(LM)%AM(:,:,:)=0.D0
        THIS(LM)%EXIST0(:)=.FALSE.
        THIS(LM)%EXISTM(:)=.FALSE.
        THIS(LM)%LM20(:)=0
        THIS(LM)%LM2M(:)=0
      ENDDO
      FOURPIBY3=16.D0*ATAN(1.D0)/3.D0  !4PI/3
!
!     ==========================================================================
      DO LM1=1,LMRXX
        L1=INT(SQRT(REAL(LM1-1,KIND=8)+DSMALL))
        LM3X=L1**2
        DO LM2=1,LMRXX
          L2=INT(SQRT(REAL(LM2-1,KIND=8)+DSMALL))
          T0=L2.EQ.L1
          TM=L2.EQ.L1-2
          IF(.NOT.(T0.OR.TM)) CYCLE
!
!         ======================================================================
!         ==  EVALUATE MATRIX                                                 ==
!         ======================================================================
          A1(:,:)=0.D0
          DO LM3=1,LM3X
            DO I=1,3
              LMI=MAPLM(I)
              CALL CLEBSCH(LMI,LM1,LM3,CG1)
              DO J=1,3
                LMJ=MAPLM(J)
                CALL CLEBSCH(LMJ,LM2,LM3,CG2)
                A1(I,J)=A1(I,J)+CG1*CG2
              ENDDO
            ENDDO
          ENDDO 
          A1(:,:)=FOURPIBY3*REAL(2*L1+1,KIND=8)*A1(:,:)
!
!         ======================================================================
!         ==  CYCLE IF THERE ARE ONLY ZERO ELEMENTS                           ==
!         ======================================================================
          TCHK=.FALSE.
          DO I=1,3
            DO J=1,3
              TCHK=TCHK.OR.(A1(I,J).NE.0.D0)
            ENDDO
          ENDDO
          IF(.NOT.TCHK) CYCLE                    
!
!         ======================================================================
!         ==  STORE L2=L1 ELEMENTS                                            ==
!         ======================================================================
          IF(T0) THEN
            I=1
            DO WHILE(THIS(LM1)%EXIST0(I)) 
              I=I+1
              IF(I.GT.10) THEN 
                CALL ERROR$MSG('SUMRULE VIOLATED - OVERFLOW A0')
                CALL ERROR$STOP('SPHERICAL_CCMAT')
              END IF
            END DO
            THIS(LM1)%A0(:,:,I)=A1(:,:)
            THIS(LM1)%LM20(I)=LM2 
            THIS(LM1)%EXIST0(I)=.TRUE.
          END IF
!
!         ======================================================================
!         ==  STORE L2=L1-2 ELEMENTS                                          ==
!         ======================================================================
          IF(TM) THEN
            I=1
            DO WHILE(THIS(LM1)%EXISTM(I)) 
              I=I+1
              IF(I.GT.10) THEN 
                CALL ERROR$MSG('SUMRULE VIOLATED - OVERFLOW AM')
                CALL ERROR$STOP('SPHERICAL_CCMAT')
              END IF
            END DO
            THIS(LM1)%AM(:,:,I)=A1(:,:)
            THIS(LM1)%LM2M(I)=LM2
            THIS(LM1)%EXISTM(I)=.TRUE.
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$CCMAT0(LM1,I,LM2,CC)
!     **************************************************************************
!     **  MATRIX P0 USED TO EVALUATE THE STRESS IN A REPRESENTATION           **
!     **  OF SPHERICAL HARMONICS                                              **
!     **  SEE SPHERICALCCMAT_MODULE FOR DETAILS                               **
!     **  LM2=0 AND CC=0.D0 ON RETURN IF I IS TOO LARGE                       **
!     **************************************************************************
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM1     ! FIRST ANGULAR MOMENTUM INDEX  
      INTEGER(4),INTENT(IN) :: I       ! COUNTER FOR LM2               
      INTEGER(4),INTENT(OUT):: LM2     ! SECOND ANGULAR MOMENTUM INDEX 
      REAL(8)   ,INTENT(OUT):: CC(3,3) ! P0_IJ(LM1,LM2)                
      LOGICAL(4)            :: TCHK
!     **************************************************************************
      IF(LM1.GT.LMRXX) CALL SPHERICALCCMAT_INITIALIZE(LM1)
      TCHK=I.LE.10
      IF(TCHK) TCHK=THIS(LM1)%EXIST0(I) 
      IF(TCHK) THEN
        LM2=THIS(LM1)%LM20(I)
        CC=THIS(LM1)%A0(:,:,I)
      ELSE
        LM2=0
        CC(:,:)=0.D0
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$CCMATM(LM1,I,LM2,CC)
!     **************************************************************************
!     **  MATRIX P0 USED TO EVALUATE THE STRESS IN A REPRESENTATION           **
!     **  OF SPHERICAL HARMONICS                                              **
!     **  SEE SPHERICALCCMAT_MODULE FOR DETAILS                               **
!     **  LM2=0 AND CC=0.D0 ON RETURN IF I IS TOO LARGE                       **
!     **************************************************************************
      USE SPHERICALCCMAT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LM1      ! FIRST ANGULAR MOMENTUM INDEX  
      INTEGER(4),INTENT(IN) :: I        ! COUNTER FOR LM2               
      INTEGER(4),INTENT(OUT):: LM2      ! SECOND ANGULAR MOMENTUM INDEX 
      REAL(8)   ,INTENT(OUT):: CC(3,3)  ! PM_IJ(LM1,LM2)                
      LOGICAL(4)            :: TCHK
!     **************************************************************************
      IF(LM1.GT.LMRXX) CALL SPHERICALCCMAT_INITIALIZE(LM1)
      TCHK=I.LE.10
      IF(TCHK) TCHK=THIS(LM1)%EXISTM(I)
      IF(TCHK) THEN
        LM2=THIS(LM1)%LM2M(I)
        CC=THIS(LM1)%AM(:,:,I)
      ELSE
        LM2=0
        CC(:,:)=0.D0
      END IF
      RETURN
      END

