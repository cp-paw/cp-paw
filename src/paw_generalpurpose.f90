!
!     .....................................................................
      SUBROUTINE ABCALPHABETAGAMMA(SWITCH,A,B,C,ALPHA,BETA,GAMMA,T)
!     **                                                                 **
!     **  CONVERTS THE LATTICE REPRESENTATION OF LENTHS OF AND ANGLES    **
!     **  BETWEEN LATTICE VECTORS                                        **
!     **                                                                 **
!     **  ATTENTION!! ANGLES ARE IN DEGREE I.E. IN UNITS OF  2*PI/360    **
!     **                                                                 **
!     *********************************************************************
      IMPLICIT NONE
      LOGICAL, INTENT(IN)    :: SWITCH  ! ABC...-> T / T->ABC...
      REAL(8),  INTENT(INOUT) :: A,B,C   ! LENGTH OF LATTICE VECTORS
      REAL(8),  INTENT(INOUT) :: ALPHA,BETA,GAMMA  ! ANGLES BETWEEN LATTICE VECTORS
      REAL(8),  INTENT(INOUT) :: T(3,3)  !LATTICE VECTORS  
      REAL(8)                 :: PI
      REAL(8)                 :: DEGREE
      REAL(8)                 :: COSA,COSB,COSG,SING
      PI=4.D0*ATAN(1.D0)
      DEGREE=2.D0*PI/360.D0
      IF(SWITCH) THEN
        COSA=COS(ALPHA*DEGREE)
        COSB=COS(BETA *DEGREE)
        COSG=COS(GAMMA*DEGREE)
        SING=SQRT(1.D0-COSG**2)
        T(1,1)=A
        T(2,1)=0.D0   
        T(3,1)=0.D0   
        T(1,2)=B*COSG
        T(2,2)=B*SING
        T(3,2)=0.D0   
        T(1,3)=C*COSB
        T(2,3)=C*(COSA-COSB*COSG)/SING
        T(3,3)=C*SQRT(SING**2+2.D0*COSA*COSB*COSG-COSA**2-COSB**2)/SING
      ELSE
        A=SQRT(T(1,1)**2+T(2,1)**2+T(3,1)**2)      
        B=SQRT(T(1,2)**2+T(2,2)**2+T(3,2)**2)      
        C=SQRT(T(1,3)**2+T(2,3)**2+T(3,3)**2)      
        COSA=(T(1,2)*T(1,3)+T(2,2)*T(2,3)+T(3,2)*T(3,3))/(B*C)
        COSB=(T(1,1)*T(1,3)+T(2,1)*T(2,3)+T(3,1)*T(3,3))/(A*C)
        COSG=(T(1,1)*T(1,2)+T(2,1)*T(2,2)+T(3,1)*T(3,2))/(A*B)
        ALPHA=ACOS(COSA)/DEGREE
        BETA =ACOS(COSB)/DEGREE
        GAMMA=ACOS(COSG)/DEGREE
      END IF
      RETURN
      END
!
!     .............................................FUNCTION TFAPOT......
      SUBROUTINE TFAPOT(R,Z,V)
!     **                                                              **
!     **  GENERALIZED THOMAS FERMI ATOMIC POTENTIAL                   **
!     **                                                              **
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: R   !RADIUS
      REAL(8),INTENT(IN)  :: Z   !ATOMIC NUMBER
      REAL(8),INTENT(OUT) :: V   ! POTENTIAL
      REAL(8),PARAMETER   :: BY3=1.D0/3.D0
      REAL(8)             :: B,X,XS,T
!     ******************************************************************
      B=(0.69395656D0/Z)**BY3
      X=R/B
      XS=SQRT(X)
      T=Z/(1.0D0+XS*(0.02747D0 - X*(0.1486D0 - 0.007298D0*X)) &
     &      + X*(1.243D0 + X*(0.2302D0 + 0.006944D0*X)))
      IF(T .LT. 1.0D0) T=1.0D0
      V=-T/R
      RETURN
      END
!
!     .....................................................GBASS........
      SUBROUTINE GBASS(RBAS,GBAS,DET)
!     ******************************************************************
!     **                                                              **
!     **  GENERATES RECIPROCAL LATTICE VECTORS G1,G2,G3               **
!     **  FROM THE REAL SPACE LATTICE VECTORS                         **
!     **          RI*GI=2*PI ; GI*GJ=0 FOR I.NE.J                     **
!     **  AND THE REAL SPACE UNIT CELL VOLUME                         **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: RBAS(3,3) ! REAL SPACE LATTIC VECTORS
      REAL(8), INTENT(OUT) :: GBAS(3,3) ! RECIPROCAL SPACE LATTICE VECTORS
      REAL(8), INTENT(OUT) :: DET       ! REAL SPACE UNIT CELL VOLUME
      REAL(8)              :: PI
      REAL(8)              :: FAC
      INTEGER(4)           :: I,J
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      GBAS(1,1) = RBAS(2,2)*RBAS(3,3) - RBAS(3,2)*RBAS(2,3)
      GBAS(2,1) = RBAS(3,2)*RBAS(1,3) - RBAS(1,2)*RBAS(3,3)
      GBAS(3,1) = RBAS(1,2)*RBAS(2,3) - RBAS(2,2)*RBAS(1,3)
      GBAS(1,2) = RBAS(2,3)*RBAS(3,1) - RBAS(3,3)*RBAS(2,1)
      GBAS(2,2) = RBAS(3,3)*RBAS(1,1) - RBAS(1,3)*RBAS(3,1)
      GBAS(3,2) = RBAS(1,3)*RBAS(2,1) - RBAS(2,3)*RBAS(1,1)
      GBAS(1,3) = RBAS(2,1)*RBAS(3,2) - RBAS(3,1)*RBAS(2,2)
      GBAS(2,3) = RBAS(3,1)*RBAS(1,2) - RBAS(1,1)*RBAS(3,2)
      GBAS(3,3) = RBAS(1,1)*RBAS(2,2) - RBAS(2,1)*RBAS(1,2)
      DET =RBAS(1,1)*GBAS(1,1) +RBAS(2,1)*GBAS(2,1) +RBAS(3,1)*GBAS(3,1)
      IF(DET.EQ.0.D0) THEN
        CALL ERROR$MSG('LATTICE VECTORS LINEAR DEPENDEND')
        PRINT*,'RBAS=',RBAS
        CALL ERROR$STOP('GBASS')
      END IF
      FAC = 2.D0*PI/DET
      DO I=1,3
        DO J=1,3
          GBAS(I,J)=GBAS(I,J)*FAC
        ENDDO
      ENDDO
      DET=ABS(DET)
      RETURN
      END
!
!     .....................................................BOXSPH.......
      SUBROUTINE BOXSPH(RBAS,X0,Y0,Z0,RMAX &
     &                 ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!     **                                                              **
!     **  BOXSPH DESCRIBES A BOX AROUND AN SPHERE                     **
!     **  CENTERED AT (Z0,Y0,Z0) AND WITH RADIUS RMAX                 **
!     **  AND RETURNS THE MINIMUM AND MAXIMUM NUMBER OF DISPLACEMENTS **
!     **  IN STEPS OF THAT CREATE POINTS WITHIN THE SPHERE            **
!     **                                                              **
!     **  |R-R0| =< RMAX        ONLY IF                               **
!     **  MIN1 =< I1 =< MAX1; MIN2 =< I2 =< MAX2; MIN3 =< I3 =< MAX3  **
!     **  WHERE:  R(I) = RBAS(I,1)*I1 + RBAS(I,2)*I2 + RBAS(I,3)*I3   **
!     **  AND     R0=(X0,Y0,Z0)                                       **
!     **                                                              **
!     **  INPUT :                                                     **
!     **    RBAS        DISPLACEMENT VECTORS                          **
!     **    RMAX        RADIUS OF THE SPHERE                          **
!     **    X0,Y0,Z0    CENTER OF THE SPHERE                          **
!     **  OUTPUT :                                                    **
!     **    MIN1,MAX1,MIN2,MAX2,MIN3,MAX3     (SEE ABOVE)             **
!     **                                                              **
!     **  WARNING: THE DISPLACEMENTS MUST NOT BE SMALLER THAN -1.E+6  **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: RBAS(3,3)
      REAL(8)   ,INTENT(IN)  :: X0,Y0,Z0
      REAL(8)   ,INTENT(IN)  :: RMAX
      INTEGER(4),INTENT(OUT) :: MIN1,MAX1
      INTEGER(4),INTENT(OUT) :: MIN2,MAX2
      INTEGER(4),INTENT(OUT) :: MIN3,MAX3
      REAL(8)                :: T1(3),T2(3),T3(3)
      REAL(8)                :: G1(3),G2(3),G3(3)
      REAL(8)                :: DET
      REAL(8)                :: XN1,XN2,XN3
      REAL(8)                :: XP0,YP0,ZP0
      INTEGER(4)             :: I
      INTEGER(4)             :: IBIG
!     ******************************************************************
      DO I=1,3
        T1(I)=RBAS(I,1)
        T2(I)=RBAS(I,2)
        T3(I)=RBAS(I,3)
      ENDDO
!     ==================================================================
!     ==  CALCULATE RECIPROCAL LATTICE VECTORS G1,G2,G3               ==
!     ==  note that the factor 2pi is dropped                         ==
!     ==================================================================
      G1(1) = T2(2)*T3(3) - T2(3)*T3(2)
      G1(2) = T2(3)*T3(1) - T2(1)*T3(3)
      G1(3) = T2(1)*T3(2) - T2(2)*T3(1)
      G2(1) = T3(2)*T1(3) - T3(3)*T1(2)
      G2(2) = T3(3)*T1(1) - T3(1)*T1(3)
      G2(3) = T3(1)*T1(2) - T3(2)*T1(1)
      G3(1) = T1(2)*T2(3) - T1(3)*T2(2)
      G3(2) = T1(3)*T2(1) - T1(1)*T2(3)
      G3(3) = T1(1)*T2(2) - T1(2)*T2(1)
      DET = T1(1)*G1(1) + T1(2)*G1(2) + T1(3)*G1(3)
      IF(DET.EQ.0.D0) THEN
        CALL ERROR$MSG('DISPLACEMENT VECTORS LINAR DEPENDEND')
        CALL ERROR$STOP('BOXSPH')
      END IF
      DO I=1,3
        G1(I)=G1(I)/DET
        G2(I)=G2(I)/DET
        G3(I)=G3(I)/DET
      ENDDO
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      XN1=SQRT(G1(1)**2+G1(2)**2+G1(3)**2)*RMAX
      XN2=SQRT(G2(1)**2+G2(2)**2+G2(3)**2)*RMAX
      XN3=SQRT(G3(1)**2+G3(2)**2+G3(3)**2)*RMAX
      XP0=G1(1)*X0+G1(2)*Y0+G1(3)*Z0
      YP0=G2(1)*X0+G2(2)*Y0+G2(3)*Z0
      ZP0=G3(1)*X0+G3(2)*Y0+G3(3)*Z0
      IBIG=1000000
      MIN1=INT(XP0-XN1+1.D0+IBIG)-IBIG
      MIN2=INT(YP0-XN2+1.D0+IBIG)-IBIG
      MIN3=INT(ZP0-XN3+1.D0+IBIG)-IBIG
      MAX1=INT(XP0+XN1+IBIG)-IBIG
      MAX2=INT(YP0+XN2+IBIG)-IBIG
      MAX3=INT(ZP0+XN3+IBIG)-IBIG
!     WRITE(*,FMT='("BX ",6F10.5)')XN1,XN2,XN3,XP0,YP0,ZP0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BOXBOX(RBAS,R0BOX,TBOX &
     &                   ,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
!     ******************************************************************
!     **  CIRCUMSCRIBES A BOX                                         **
!     **  DEFINED BY THE LOWER LEFT FORDER CORNER R0BOX               **
!     **  AND THREE EDGE VECTORS TBOX                                 **
!     **  BY ANOTHER BOX DEFINED THROUGH THE GRID TRANSLATION VECTORS **
!     **  RBAS SUCH THAT ALL POINTS IN THE CIRCUMSCRIBED BOX          **
!     **  FULFILL R=RBAS*X WITH XMIN<X<XMAX FOR ALL THREE COORDINATES **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: R0BOX(3)
      REAL(8)   ,INTENT(IN) :: TBOX(3,3)
      REAL(8)   ,INTENT(OUT):: XMIN1,XMAX1
      REAL(8)   ,INTENT(OUT):: XMIN2,XMAX2
      REAL(8)   ,INTENT(OUT):: XMIN3,XMAX3
!     INTEGER(4),PARAMETER  :: NAUXX=300 ! USED FOR DEICD
!     REAL(8)               :: AUX(NAUXX)! USED FOR DEICD
!     REAL(8)               :: RCOND,DET ! USED FOR DEICD         
      REAL(8)               :: RBASIN(3,3)  !INVERSE OF RBAS
      REAL(8)               :: X(3)
      REAL(8)               :: DT(3)
      REAL(8)               :: R(3)
      INTEGER(4)            :: I,J,I1,I2,I3
      LOGICAL(4)            :: TFIRST
!     ******************************************************************
!
!     ==================================================================
!     ==  CALCULATE RBASIN, THE INVERSE OF RBAS                       ==
!     ==================================================================
      CALL LIB$INVERTR8(3,RBAS,RBASIN)
!      DO I=1,3
!        DO J=1,3
!          RBASIN(I,J)=RBAS(I,J)
!        ENDDO
!      ENDDO
!      CALL DGEICD(RBASIN,3,3,0,RCOND,DET,AUX,NAUXX)
!
!     ==================================================================
!     ==  CALCULATE THE 8 CORNERS R OF THE BOX                        ==
!     ==  AND EXPRESS THEM IN TERMS R=RBAS*X                          ==
!     ==  THE LOWEST AND LARGEST COORDINATES OD X DETERMINE           == 
!     ==  XMIN,XMAX IN EACHB DIRECTION                                ==
!     ==================================================================
      TFIRST=.TRUE.
      DO I1=0,1
        DT(1)=DBLE(I1)
        DO I2=0,1
          DT(2)=DBLE(I2)
          DO I3=0,1
            DT(3)=DBLE(I3)
!           == EVALUATE CORNER OF THE BOX ====================
            DO I=1,3
              R(I)=R0BOX(I)
              DO J=1,3
                R(I)=R(I)+TBOX(I,J)*DT(J)
              ENDDO
            ENDDO
!           == TRANSFORM CORNER INTO RELATIVE COORDINATES
            DO I=1,3
              X(I)=0.D0
              DO J=1,3
                X(I)=X(I)+RBASIN(I,J)*R(J)
              ENDDO
            ENDDO
            IF(TFIRST) THEN
              XMIN1=X(1)
              XMIN2=X(2)
              XMIN3=X(3)
              XMAX1=X(1)
              XMAX2=X(2)
              XMAX3=X(3)
              TFIRST=.FALSE.
            ELSE
              XMIN1=MIN(XMIN1,X(1))
              XMAX1=MAX(XMAX1,X(1))
              XMIN2=MIN(XMIN2,X(2))
              XMAX2=MAX(XMAX2,X(2))
              XMIN3=MIN(XMIN3,X(3))
              XMAX3=MAX(XMAX3,X(3))
            END IF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!      .................................................................
       SUBROUTINE NOMOM(STRING,NAT,R1,R2,RMASS)
!      **                                                             **
!      ** ROTATES AND TRANSLATES THE VECTORS R2 SO THAT THE ANGULAR   **
!      ** AND/OR TRANSLATIONAL MOMENTUM FROM R1 TO R2 VANISHES        **
!      **                                                             **
!      **  INPUT:                                                     **
!      **  STRING    'TR' SUBTRACT TRANSLATION AND ROTATION           **
!      **            'T'  SUBTRACT TRANSLATION                        **
!      **            'R'  SUBTRACT ROTATION                           **
!      **  NAT       NUMBER OF ATOMS                                  **
!      **  R1        REFERENCE STRUCTURE                              **
!      **  R2        STRUCTURE TO BE ROTATED                          **
!      **  RMASS     ATOMIC MASSES                                    **
!      **                                                             **
!      **  OUTPUT:                                                    **
!      **  R2        ROTATED STRUCTURE                                **
!      **                                                             **
!      **  REMARKS:                                                   **
!      **  PROGRAM WILL FAIL FOR ROTATIONS OF PI                      **
!      **  PROGRAM IS NOT PROTECTED FOR 1-D SYSTEM                    **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: STRING
       INTEGER(4),INTENT(IN)   :: NAT
       REAL(8)   ,INTENT(IN)   :: R1(3,NAT)
       REAL(8)   ,INTENT(INOUT):: R2(3,NAT)
       REAL(8)   ,INTENT(IN)   :: RMASS(NAT)
       REAL(8)                 :: T(3)
       REAL(8)                 :: RL(3)
       REAL(8)                 :: RI(3,3)
       REAL(8)                 :: ROT(3,3)
       REAL(8)                 :: W(3)
       REAL(8)                 :: C(3)
!      REAL(8)                 :: WORK(300)
       REAL(8)                 :: TOTM,WTOT,RM,SVAR
       REAL(8)                 :: X1,Y1,Z1
       REAL(8)                 :: X2,Y2,Z2
       REAL(8)                 :: DX,DY,DZ
       REAL(8)                 :: XM,YM,ZM
       INTEGER(4)              :: I,J,K,IAT,ITER
       REAL(8)   ,PARAMETER    :: TOL=1.D-6
       INTEGER(4),PARAMETER    :: ITERX=10000
!      *****************************************************************
       IF(STRING.NE.'T'.AND.STRING.NE.'R'.AND.STRING.NE.'TR') THEN
         CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED IN NOMOM')
         CALL ERROR$STOP('NOMOM')
       END IF 
!
!      =================================================================
!      == SUBTRACT TRANSLATION T                                      ==
!      =================================================================
       TOTM=0.D0
       DO I=1,3
         T(I)=0.D0
         C(I)=0.D0
       ENDDO
       DO IAT=1,NAT
         TOTM=TOTM+RMASS(IAT)
         DO I=1,3
           T(I)=T(I)+RMASS(IAT)*(R1(I,IAT)-R2(I,IAT))
           C(I)=C(I)+RMASS(IAT)*R2(I,IAT)
         ENDDO
       ENDDO
       DO I=1,3
         T(I)=T(I)/TOTM
         C(I)=C(I)/TOTM
       ENDDO 
       DO IAT=1,NAT
         DO I=1,3
           R2(I,IAT)=R2(I,IAT)+T(I)
         ENDDO
       ENDDO
       IF(STRING.EQ.'T') RETURN
       IF(NAT.EQ.1) GOTO 2000
!
!      =================================================================
!      == SUBTRACT ROTATION                                           ==
!      =================================================================
!      == R(S)=U(S)*(R1+(R2-R1))*S)
!      == WHERE U IS A ROTATION MATRIX (1+W)**N/N AND
!      ==   (   0 ; PHIZ; -PHIY)      
!      == W=(-PHIZ;   0 ;  PHIX)
!      ==   ( PHIY;-PHIX;    0 )
!
!      == CALCULATE CENTER OF GRAVITY OF R2
!      == ASSUMING THAT CENTERS OF GRAVITY OF R1 AND R2 COINCIDE
       TOTM=0.D0
       DO I=1,3
         C(I)=0.D0
       ENDDO
       DO IAT=1,NAT
         TOTM=TOTM+RMASS(IAT)
         DO I=1,3
           C(I)=C(I)+RMASS(IAT)*R2(I,IAT)
         ENDDO
       ENDDO
       DO I=1,3
         C(I)=C(I)/TOTM
       ENDDO 
       ITER=0
       DO ITER=1,ITERX 
!
!      == CALCULATE ANGULAR MOMENTUM AND
         DO I=1,3
           RL(I)=0.D0 
           DO J=1,3
             RI(I,J)=0.D0 
           ENDDO
         ENDDO
         DO IAT=1,NAT
           X1=R1(1,IAT)-C(1)
           Y1=R1(2,IAT)-C(2)
           Z1=R1(3,IAT)-C(3)
           X2=R2(1,IAT)-C(1)
           Y2=R2(2,IAT)-C(2)
           Z2=R2(3,IAT)-C(3)
           XM=0.5D0*(X1+X2)
           YM=0.5D0*(Y1+Y2)
           ZM=0.5D0*(Z1+Z2)
           DX=X2-X1
           DY=Y2-Y1
           DZ=Z2-Z1
           RM=RMASS(IAT)
           RL(1)=RL(1)+RM*(YM*DZ-ZM*DY)
           RL(2)=RL(2)+RM*(ZM*DX-XM*DZ)
           RL(3)=RL(3)+RM*(XM*DY-YM*DX)
           SVAR=-RM*(XM**2+YM**2+ZM**2)
           RI(1,1)=RI(1,1)+RM*XM**2-SVAR
           RI(1,2)=RI(1,2)+RM*XM*YM
           RI(1,3)=RI(1,3)+RM*XM*ZM
           RI(2,3)=RI(2,3)+RM*YM*ZM
           RI(2,2)=RI(2,2)+RM*YM**2-SVAR
           RI(3,3)=RI(3,3)+RM*ZM**2-SVAR
         ENDDO
         RI(2,1)=RI(1,2)
         RI(3,1)=RI(1,3)
         RI(3,2)=RI(2,3)
!
!        -- INVERT RI 
         CALL LIB$INVERTR8(3,RI,RI)
!        NAUX=300
!        CALL DGEICD(RI,3,3,0,RCOND,DET,WORK,NAUX)
!
!        ==  CALCULATE ANGULAR VELOCITY W
         DO I=1,3
           W(I)=0.D0
           DO J=1,3
             W(I)=W(I)-RI(I,J)*RL(J)
           ENDDO
         ENDDO
         WTOT=SQRT(W(1)**2+W(2)**2+W(3)**2)
!
         IF(WTOT.LT.TOL.OR.ITER.GT.ITERX) GOTO 2000
         CALL ROTATIONMATRIX(W,ROT)       
         DO IAT=1,NAT
           DO I=1,3
             SVAR=0.D0
             DO K=1,3
               SVAR=SVAR+ROT(I,K)*(R2(K,IAT)-C(K))
             ENDDO
             W(I)=SVAR+C(I)
           ENDDO
           DO I=1,3
             R2(I,IAT)=W(I)
           ENDDO
         ENDDO
       ENDDO
       CALL ERROR$MSG('ITERATION IN NOMOM NOT CONVERGED')
       CALL ERROR$R8VAL('WTOT',WTOT)
       CALL ERROR$STOP('NOMOM')
2000   CONTINUE
       IF(STRING.EQ.'R') THEN
         DO IAT=1,NAT
           DO I=1,3
             R2(I,IAT)=R2(I,IAT)-T(I)
           ENDDO
         ENDDO
       END IF 
       RETURN
       END
!
!      .................................................................
       SUBROUTINE EULERANGLE(PHI,THETA,PSI,R)
!      **                                                             **
!      ** ROTATION MATRIX FROM EULER ANGLES (SEE GOLDSTEIN)           **
!      **                                                             **
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: PHI
       REAL(8)   ,INTENT(IN) :: THETA
       REAL(8)   ,INTENT(IN) :: PSI
       REAL(8)   ,INTENT(OUT):: R(3,3)
       REAL(8)               :: C1,S1
       REAL(8)               :: C2,S2
       REAL(8)               :: C3,S3
       INTEGER(4)            :: I,J,K
       REAL(8)               :: SVAR
       LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
!      ******************************************************************
       C1=COS(PHI)
       S1=SIN(PHI)
       C2=COS(THETA)
       S2=SIN(THETA)
       C3=COS(PSI)
       S3=SIN(PSI)
       R(1,1)=C3*C1-C2*S1*S3
       R(1,2)=-S3*C1-C2*S1*C3
       R(1,3)=S2*S1
       R(2,1)=-S3*C1-C2*S1*C3
       R(2,2)=-S3*S1+C2*C1*C3
       R(2,3)=C3*S2
       R(3,1)=S2*S1
       R(3,2)=-S2*C1
       R(3,3)=C2
!
!      =================================================================
!      == TEST
!      =================================================================
       IF(TTEST) THEN
         DO I=1,3
           DO J=1,3
             SVAR=0.D0
             IF(I.EQ.J) SVAR=-1.D0
             DO K=1,3
               SVAR=SVAR+R(I,K)*R(J,K)
             ENDDO
             IF(ABS(SVAR).GT.1.D-8) THEN
               CALL ERROR$STOP('EULERANGLE')
             END IF
           ENDDO
         ENDDO
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE ROTATIONMATRIX(PHI,R)
!      **                                                             **
!      **  CONSTRUCTS ROTATION MATRIX FOR A GIVEN ANGLE VECTOR        **
!      **  THE ANGLE IS |PHI| THE AXIS IS PHI/|PHI|                   **
!      **  THE ROTATION IS COUNTER CLOCKWISE                          **
!      **  THE TRANSFORMED VECTOR IS XPRIME=R*X                       **
!      **                                                             **
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN)  :: PHI(3)
       REAL(8)   ,INTENT(OUT) :: R(3,3)
       REAL(8)   ,PARAMETER   :: TOL=1.D-6
       REAL(8)                :: R1(3,3)
       REAL(8)                :: R2(3,3)
       REAL(8)                :: ABSPHI
       REAL(8)                :: COSPHI,SINPHI
       REAL(8)                :: SVAR
       INTEGER(4)             :: I,J,K
!      *****************************************************************
       ABSPHI=SQRT(PHI(1)**2+PHI(2)**2+PHI(3)**2)
!
!      =================================================================
!      == RETURN IDENTITY IF ANGLE=0                                  ==
!      =================================================================
       IF(ABSPHI.LT.TOL) THEN
         DO I=1,3
           DO J=1,3
             R(I,J)=0.D0
           ENDDO
           R(I,I)=1.D0
         ENDDO
         RETURN
       END IF
       COSPHI=COS(ABSPHI)
       SINPHI=SIN(ABSPHI)
       DO I=1,3
         R1(I,1)=PHI(I)/ABSPHI
       ENDDO
!
!      =================================================================
!      == TREAT ROTATION ABOUT Z-AXIS EXTRA                           ==
!      =================================================================
       IF(ABS(R1(3,1)-1.D0).LT.TOL) THEN
         R(1,1)=COSPHI         
         R(1,2)=-SINPHI
         R(1,3)=0.D0
         R(2,1)=SINPHI
         R(2,2)=COSPHI
         R(2,3)=0.D0
         R(3,1)=0.D0
         R(3,2)=0.D0
         R(3,3)=1.D0
         RETURN
       ELSE IF(ABS(R1(3,1)+1.D0).LT.TOL) THEN
         R(1,1)=COSPHI
         R(1,2)=SINPHI
         R(1,3)=0.D0
         R(2,1)=-SINPHI
         R(2,2)=COSPHI
         R(2,3)=0.D0
         R(3,1)=0.D0
         R(3,2)=0.D0
         R(3,3)=1.D0
         RETURN
       END IF
!
!      =================================================================
!      == SET UP ORTHONORMAL TRIAL VECTORS                            ==
!      =================================================================
       SVAR=1.D0/SQRT(R1(1,1)**2+R1(2,1)**2)
       R1(1,2)=R1(2,1)*SVAR
       R1(2,2)=-R1(1,1)*SVAR
       R1(3,2)=0.D0
       R1(1,3)=R1(3,1)*R1(1,1)*SVAR
       R1(2,3)=R1(3,1)*R1(2,1)*SVAR
       R1(3,3)=(R1(3,1)*R1(3,1)-1.D0)*SVAR
!
!      =================================================================
!      == ROTATE TRIAL VECTORS                                        ==
!      =================================================================
       DO I=1,3
         R2(I,1)=R1(I,1)
         R2(I,2)=R1(I,2)*COSPHI+R1(I,3)*SINPHI
         R2(I,3)=-R1(I,2)*SINPHI+R1(I,3)*COSPHI      
       ENDDO
!
!      =================================================================
!      == EXTRACT ROTATION MATRIX                                     ==
!      =================================================================
       DO I=1,3
         DO J=1,3
           R(I,J)=0.D0
           DO K=1,3
             R(I,J)=R(I,J)+R2(I,K)*R1(J,K)
           ENDDO
         ENDDO
       ENDDO
       RETURN
       END
!
!     ..................................................................
      SUBROUTINE BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
!     ******************************************************************
!     **                                                              **
!     **  FINDS THE ZERO OF A MONOTONIC FUNCTION Y(X)                 **
!     **  BY DOUBLING STEP SIZE AND SUBSEQUENT BYSECTION.             **
!     **  THE ROUTINE MUST BE CALLED IN A LOOP                        **
!     **  THAT STOPS WHEN CONVERGENCE IS OBTAINED                     **
!     **  AND SUPPLIES NEW FUNCTION VALUES Y0 FOR THE VALUE X0        **
!     **  THAT IS SUPPLIED BY THIS ROUTINE.                           **
!     **  REMARK: CALL THIS ROUTINE ONCE BEFORE THE LOOP              **
!     **                                                              **
!     **   INITIALIZE BEFORE CALL WITH                                **
!     **   ISTART=1 FOR MONOTONICALLY INCREASING OR WITH              **
!     **   ISTART=-1 FOR MONOTONICALLY DECREASING FUNCTIONS:          **
!     **   AND WITH X0 THE STARTING ARGUMENT AND WITH DX THE          **
!     **   STARTING STEP SIZE                                         **
!     **                                                              **
!     **   DO NOT CHANGE THE VALUES FOR IBI,DX,XM,Y DURING THE        **
!     **   ITERATION                                                  **
!     **                                                              **
!     **   X0=??                                                      **
!     **   DX=??                                                      **
!     **   CALL BISEC(1,IBI,X0,Y0,DX,XM,YM)                           **
!     **   DO I=1,MAX                                                 **
!     **     CALCULATE Y0 FOR VALUE X0                                **
!     **     CALL BISEC(0,IBI,X0,Y0,DX,XM,YM)                         **
!     **     IF(ABS(Y0).LT.TOL) EXIT                                  **
!     **   ENDDO                                                      **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: ISTART  !=1 BEFORE ITERATION/ =0 OTHERWISE
      INTEGER(4),INTENT(INOUT) :: IBI     ! SWITCH BETWEEN EXPANSION AND CONTRACTION 
      REAL(8)   ,INTENT(INOUT) :: X0      ! CURRENT ARGUMENT
      REAL(8)   ,INTENT(IN)    :: Y0      ! FUNCTION VALUE AT X0
      REAL(8)   ,INTENT(INOUT) :: DX      ! STEP WIDTH
      REAL(8)   ,INTENT(INOUT) :: YM      ! FUNCTION VALUE AT XM
      REAL(8)   ,INTENT(INOUT) :: XM      ! PREVIOUS ARGUMENT
      REAL(8)   ,save          :: SLOPE
      REAL(8)                  :: XP
!     ******************************************************************
!
!     ==   STARTUP
      IF(ISTART.ne.0) THEN
        slope=-1.d0
        if(istart.gt.0) slope=1.d0
        ISTART=0
        IBI=0
        YM=0.D0
        XM=X0-DX
        RETURN
      END IF
!
!     ==  SWITCH TO BISECTING
      IF(IBI.EQ.0.AND.YM*Y0.LT.0.D0) THEN
!       PRINT*,'Y0,YM ',Y0,YM
        IBI=1
        DX=0.25D0*DX
      END IF
!
!     == PROPAGATE
      XP=X0-DSIGN(DX,Y0*SLOPE)
      IF(IBI.EQ.0)DX=DX*2.D0
      IF(IBI.EQ.1)DX=DX/2.D0
      XM=X0
      YM=Y0
      X0=XP
      RETURN
      END
!
!...............................................................................
MODULE BROYDEN_MODULE
LOGICAL(4)         :: TON=.FALSE.
INTEGER(4)         :: NSTEPX=0
INTEGER(4)         :: NSTEP=0
INTEGER(4)         :: NX=0
REAL(8)            :: ALPHA
REAL(8),ALLOCATABLE :: XPREV(:,:)
REAL(8),ALLOCATABLE :: YPREV(:,:)
REAL(8),ALLOCATABLE :: WEIGHT(:)
END MODULE BROYDEN_MODULE
!      .........................................................................
       SUBROUTINE BROYDEN$NEW(NX_,NSTEPX_,ALPHA_)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NX_
       INTEGER(4),INTENT(IN)    :: NSTEPX_
       REAL(8)   ,INTENT(IN)    :: ALPHA_
!      *************************************************************************
       IF(TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT ALREADY IN USE')
         CALL ERROR$STOP('BROYDEN$NEW')
       END IF
       TON=.TRUE.
       NSTEP=0
       NX=NX_
       NSTEPX=NSTEPX_
       ALPHA=ALPHA_
       ALLOCATE(XPREV(NX,NSTEPX))
       ALLOCATE(YPREV(NX,NSTEPX))
       ALLOCATE(WEIGHT(NX))
       WEIGHT(:)=1.D0
       XPREV(:,:)=0.D0
       YPREV(:,:)=0.D0
       RETURN
       END
!      .........................................................................
       SUBROUTINE BROYDEN$CLEAR
       USE BROYDEN_MODULE
       IMPLICIT NONE
!      *************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$CLEAR')
       END IF
       TON=.FALSE.
       NSTEPX=0
       NSTEP=0
       NX=0
       ALPHA=0.D0
       DEALLOCATE(XPREV)
       DEALLOCATE(YPREV)
       DEALLOCATE(WEIGHT)
       RETURN
       END
!
!      .........................................................................
       SUBROUTINE BROYDEN$SETWEIGHT(NX_,WEIGHT_)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NX_
       REAL(8)   ,INTENT(IN)  :: WEIGHT_(NX_)
!      *************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$SETWEIGHT')
       END IF
       IF(NX_.NE.NX) THEN
         CALL ERROR$MSG('SIZE INCONSISTENT')
         CALL ERROR$STOP('BROYDEN$SETWEIGHT')
       END IF
       WEIGHT(:)=WEIGHT_(:)
       RETURN
       END
!
!      .........................................................................
       SUBROUTINE BROYDEN$STEP(NX_,X,Y)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NX_
       REAL(8)   ,INTENT(INOUT) :: X(NX_)
       REAL(8)   ,INTENT(IN)    :: Y(NX_)
       REAL(8)   ,ALLOCATABLE   :: DX(:,:)
       REAL(8)   ,ALLOCATABLE   :: DY(:,:)
       REAL(8)   ,ALLOCATABLE   :: B(:,:)
       REAL(8)   ,ALLOCATABLE   :: BINV(:,:)
       REAL(8)                  :: WY(NX_)
       INTEGER(4)               :: I
!      *************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
       IF(NX_.NE.NX) THEN
         CALL ERROR$MSG('SIZE INCONSISTENT')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
!PRINT*,'NSTEP',NSTEP
!
!      =================================================================
!      == APPLY WEIGHTING                                             ==
!      =================================================================
       WY(:)=WEIGHT(:)*Y(:)
!
!      =================================================================
!      == SIMPLE MIXING IN THE FIRST STEP                             ==
!      =================================================================
       IF(NSTEP.EQ.0) THEN
         IF(NSTEPX.GT.0)THEN
           NSTEP=1
           XPREV(:,1)=X(:)     
           YPREV(:,1)=WY(:)     
         END IF
         X=X+ALPHA*WY
         RETURN
       END IF
!
!      =================================================================
!      == DETERMINE INVERSE HESSIAN ALPHA+DX OTIMES DY                ==
!      =================================================================
       ALLOCATE(DX(NX,NSTEP))
       ALLOCATE(DY(NX,NSTEP))
       DO I=1,NSTEP
         DY(:,I)=YPREV(:,I)-WY(:)  
         DX(:,I)=XPREV(:,I)-X(:)+ALPHA*DY(:,I)
       ENDDO
       ALLOCATE(B(NSTEP,NSTEP))
       ALLOCATE(BINV(NSTEP,NSTEP))
       B=MATMUL(TRANSPOSE(DY),DY)   !OVERLAP MATRIX OF DY
!PRINT*,'B',B
       CALL LIB$INVERTR8(NSTEP,B,BINV)           
!ALLOCATE(W(NX,NSTEP))
!W=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
!PRINT*,'W ',MATMUL(TRANSPOSE(W),DY)
!DEALLOCATE(W)
       DY=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
       DEALLOCATE(B)
       DEALLOCATE(BINV)
!
!      =================================================================
!      == STORE HISTORY                                               ==
!      =================================================================
       IF(NSTEP.LT.NSTEPX)NSTEP=NSTEP+1
       DO I=NSTEP,2,-1
         YPREV(:,I)=YPREV(:,I-1)
         XPREV(:,I)=XPREV(:,I-1)
       ENDDO
       XPREV(:,1)=X(:)     
       YPREV(:,1)=WY(:)     
!
!      =================================================================
!      == PREDICT NEW VECTOR                                          ==
!      =================================================================
       X=X+ALPHA*WY-MATMUL(DX,MATMUL(TRANSPOSE(DY),WY))
       DEALLOCATE(DX)
       DEALLOCATE(DY)
       RETURN
       END
!
!     .....................................................GAUSSN.......
      SUBROUTINE GAUSSN(L,ALPHA,C)
!     ******************************************************************
!     **                                                              **
!     ** CALCULATES THE NORMALIZATION CONSTANT C OF GENERALIZED,      **
!     ** ANGULAR MOMENTUM DEPENDENT GAUSSIANS.                        **
!     **                                                              **
!     **   G_L(R) = C * EXP(-ALPHA*R**2) * R**L * Y_LM(R)             **
!     **                                                              **
!     ** A NORMALIZED GAUSSIAN HAS A MULTIPOLE MOMENT Q EQUAL TO ONE: **
!     **                                                              **
!     **   Q = INT(DR**3): G_L(R) * R**L * Y_LM(R) !=! 1              **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: ALPHA
      REAL(8)   ,INTENT(OUT):: C
      REAL(8)               :: PI
      REAL(8)               :: RINT
      INTEGER(4)            :: K,IFAC,I
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
!     ==================================================================
!     ==  CALCULATE INT(DR): R**(2*L+2) *EXP(-ALPHA * R**2)           ==
!     ==  SEE BRONSTEIN P66.                                          ==
!     ==================================================================
      K=L+1
      IFAC=1
      DO I=2,K
        IFAC=IFAC*(2*I-1)
      ENDDO
      RINT=SQRT(PI/(4.D0*ALPHA))*DBLE(IFAC)/(2.D0*ALPHA)**K
      C=1.D0/RINT
      RETURN
      END
!
      subroutine cg$test()
!     **                                                              **
!     ** test routine for conjugate gradient                          **
!     **                                                              **
!     **                                                              **
      implicit none
      INTEGER(4), PARAMETER :: N=2
      REAL(8)               :: R(N),F1(N),F2(N),D1(N),D2(N)
      REAL(8)               :: e
      INTEGER(4),PARAMETER  :: NITER=10
      INTEGER(4)            :: ITER,inner
      real(8)               :: lambda
      logical(4)            :: tconv
!     ******************************************************************
      R=(/10.D0,2.D0/)
      CALL ETOT(N,R,E,F1)
write(*,fmt='(i5," e ",f20.15," r=",2f10.5," f= ",2e10.3)')0,E,r,f1
      D1=F1
      DO ITER=1,NITER
        lambda=1.d-2/SQRT(DOT_PRODUCT(F1,F1))
        do inner=1,100
          CALL ETOT(N,R+lambda*d1,E,F2)
          call cg$linesearch(n,f1,d1,f2,lambda,tconv)
          if(tconv) exit
        enddo
        if(.not.tconv) stop 'not converged'
        R=R+D1*lambda
write(*,fmt='(i5," e ",f20.15," r=",2f10.5," f= ",2e10.3)')ITER,E,r,f2
        CALL cg$NEWDIR(n,f1,d1,f2,d2)
        F1=F2
        D1=D2
      ENDDO         
      stop
      contains
!     ....................................................................
        SUBROUTINE ETOT(N,R,E,F)
        integer(4),intent(in) :: n
        real(8)   ,intent(in) :: R(n)  
        real(8)   ,intent(OUT):: E
        real(8)   ,intent(OUT):: F(n)  
        INTEGER(4),PARAMETER  :: NLOC=2
        REAL(8)   ,PARAMETER  :: C=2.D0
        REAL(8)   ,PARAMETER  :: B(NLOC)=(/0.D0,0.D0/)
        REAL(8)               :: A(NLOC,NLOC)
!       *********************************************************************
        a(:,1)=(/2.D+1,0.D0/)
        a(:,2)=(/0.D0,2.D-3/)
        IF(N.NE.NLOC) STOP 'ERROR IN ETOT'
        E=C-DOT_PRODUCT(B,R)+0.5D0*DOT_PRODUCT(R,MATMUL(A,R))
        F=B-MATMUL(A,R)
        RETURN
        END SUBROUTINE ETOT
      end
!
!     ..................................................................
      subroutine cg$linesearch(n,f1,d1,f2,dlambda,tconv)
!     **                                                              **
!     ** conjugate gradient line search                               **
!     **                                                              **
!     **  adjusts lambda such that the force at x+d1*lambda           **
!     **  parallel to the search direction d1 converges to zer        **
!     **                                                              **
!     **  x(lambda)=x1+d1*lambda                                      **
!     **  f(lambda)=f1+lambda (f2-f1)/lambda_in                       **
!     **  f2=f(lambda_in)                                             **
!     **  f1=f(lambda=0)                                              **
!     **                                                              **
      implicit none
      integer(4),intent(in)   :: n
      real(8)   ,intent(in)   :: f1(n)  ! force at x1
      real(8)   ,intent(in)   :: d1(n)  ! x(LAMBDA)=x1+lambda*d1
      real(8)   ,intent(in)   :: f2(n)  ! f(LAMBDA) WITH F(LAMBDA)*D1=0
      real(8)   ,intent(inout):: dlAMBDA  ! NEXT DIRECTION FOR LINE SEARCH
      logical(4),intent(out)  :: tconv
      reAL(8)                 :: SVAR1,svar2
!     ******************************************************************
      svar1=dot_product(d1,f2)
      svar2=dot_product(d1,f2-f1)
      if(dlambda*svar2.lt.0.d0) then
        dlambda=-svar1/svar2*dlambda
      else
print*,'svar2 ',svar2,dlambda
print*,'warning! hessian not positive definite; switch to stepping'
!       == correction for wrong curvature
        if(svar2.gt.0.d0) then
          dlambda=1.d-2/sqrt(dot_product(d1,d1))
        else
          dlambda=-1.d-2/sqrt(dot_product(d1,d1))
        end if
      end if
      tconv=(abs(svar2/svar1).lt.1.d-4)
      RETURN
      END
!
!     ..................................................................
      subroutine cg$linesearch_old(n,f1,d1,f2,lambda,tconv)
!     **                                                              **
!     ** conjugate gradient line search                               **
!     **                                                              **
!     **  adjusts lambda such that the force at x+d1*lambda           **
!     **  parallel to the search direction d1 converges to zer        **
!     **                                                              **
!     **  x(lambda)=x1+d1*lambda                                      **
!     **  f(lambda)=f1+lambda (f2-f1)/lambda_in                       **
!     **  f2=f(lambda_in)                                             **
!     **  f1=f(lambda=0)                                              **
!     **                                                              **
      implicit none
      integer(4),intent(in)   :: n
      real(8)   ,intent(in)   :: f1(n)  ! force at x1
      real(8)   ,intent(in)   :: d1(n)  ! x(LAMBDA)=x1+lambda*d1
      real(8)   ,intent(in)   :: f2(n)  ! f(LAMBDA) WITH F(LAMBDA)*D1=0
      real(8)   ,intent(inout):: lAMBDA  ! NEXT DIRECTION FOR LINE SEARCH
      logical(4),intent(out)  :: tconv
      reAL(8)                 :: SVAR1,svar2
!     ******************************************************************
      svar1=dot_product(d1,f1)
      svar2=dot_product(d1,f2-f1)
      if(lambda*svar2.lt.0.d0) then
        lambda=-svar1/svar2*lambda
      else
print*,'warning! hessian not positive definite; switch to stepping'
!       == correction for wrong curvature
        if(svar1+svar2.gt.0.d0) then
          lambda=lambda+1.d-2/sqrt(dot_product(d1,d1))
        else
          lambda=lambda-1.d-2/sqrt(dot_product(d1,d1))
        end if
      end if
      tconv=(abs((svar2+svar1)/svar1).lt.1.d-4)
      RETURN
      END
!
!     ..................................................................
      subroutine cg$newdir(n,f1,d1,f2,d2)
!     **                                                              **
!     ** conjugate gradient new search direction                      **
!     **                                                              **
      implicit none
      integer(4),intent(in) :: n
      real(8)   ,intent(in) :: f1(n)  ! force at x1
      real(8)   ,intent(in) :: d1(n)  ! x(LAMBDA)=x1+lambda*d1
      real(8)   ,intent(in) :: f2(n)  ! f(LAMBDA) WITH F(LAMBDA)*D1=0
      real(8)   ,intent(out):: d2(n)  ! NEXT DIRECTION FOR LINE SEARCH
!     ******************************************************************
      d2=f2+d1*dot_product(f2-f1,f2)/dot_product(f1,f1)
      RETURN
      END
!!$!
!!$!     .....................................................INITDC.......
!!$      SUBROUTINE INITDC(N,CARRAY)
!!$      COMPLEX(8) CARRAY(N)
!!$      CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$      WRITE(NFILO,FMT='("THE ROUTINE INITDC SHOULD NOT BE USED")')
!!$      DO I=1,N
!!$        CARRAY(I)=(0.D0,0.D0)
!!$      ENDDO
!!$      RETURN
!!$      END
!
!      .................................................................
       SUBROUTINE GAUSS_RANDOM_NUMBER(HARVEST)
!      *****************************************************************
!      **                                                             **
!      **  GENERATES GAUSS DISTRIBUTED RANDOM NUMBERS                 **
!      **  WHICH APPROXIMATE THE DISTRIBUTION P(X)=EXP(-X^2)/SQRT(PI) **
!      **  THE DISTRIBUTION IS NORMALIZED AND HAS THE VARIANCE        **
!      **            INT:P(X)X^2=0.5                                  **
!      **                                                             **
!      **  REMARKS:                                                   **
!      **  -- THE DISTRIBUTION CAN BE IMPROVED BY INCREASING N        **
!      **  -- N=6 AND N=24 AVOID CALCULATING THE SQUARE ROOT          **
!      **  -- TO OBTAIN A DISTRIBUTION WITH MEAN SQUARE 1, MULTIPLY   **
!      **     HARVEST SQRT(2)                                         **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),PARAMETER   :: N=6
       LOGICAL(4),PARAMETER   :: T1=(N.EQ.6)
       LOGICAL(4),PARAMETER   :: T2=(N.EQ.24)
       REAL(8)   ,INTENT(OUT) :: HARVEST
       REAL(8)                :: DEL
       INTEGER(4)             :: I
!      *****************************************************************
       HARVEST=0.D0
       DO I=1,N
         CALL RANDOM_NUMBER(DEL)
         HARVEST=HARVEST+DEL-0.5D0
       ENDDO
       IF(T1) THEN
         RETURN
       ELSE IF(T2) THEN
         HARVEST=0.5D0*HARVEST
         RETURN
       ELSE
         HARVEST=HARVEST*SQRT(6.D0/DBLE(N))
         RETURN 
       END IF 
       END
!
!.......................................................................
MODULE SORT_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: SORT                                                       **
!**                                                                   **
!**  SORT SUCH THAT A GIVEN  R8 ARRAY INCREASES MONOTONICALLY         **
!**                                                                   **
!**  INTERFACE:                                                       **
!**    SORT$SET                                                       **
!**    SORT$FLIP                                                      **
!**    SORT$ORDER                                                     **
!**    SORT$RESTART                                                   **
!**    SORT$UNSET                                                     **
!**                                                                   **
!**  USE THE FOLLOWING LOOP:                                          **
!**     __CREAT INDEX ARRAY FOR ORDERING_________________________     **
!**     CALL SORT$SET(N,F)                                            **
!**                                                                   **
!**     __REORDER ARRAY_THIS OPERATION CAN BE REPEATED SEVERAL TIMES_ **
!**     CALL SORT$RESTART                                             **
!**     CALL SORT$FLIP(FROM,TO)                                       **
!**     DO WHILE (FROM.NE.0.OR.TO.NE.0)                               **
!**       IF(TO.EQ.0) THEN                                            **
!**         SVAR=F(FROM,ISTEP)                                        **
!**       ELSE IF(FROM.EQ.0) THEN                                     **
!**         F(TO)=SVAR                                                **
!**      ELSE                                                         **
!**         F(TO)=F(FROM)                                             **
!**       END IF                                                      **
!**       CALL SORT$FLIP(FROM,TO)                                     **
!**     ENDDO                                                         **
!**                                                                   **
!**     __CLEAR MEMORY OF SORT OBJECT TO ALLOW SETTING OF NEW ARRAY   **
!**     CALL SORT$UNSET                                               **
!**                                                                   **
!***********************************************************************
LOGICAL(4)             :: TSET=.FALSE.
INTEGER(4)             :: LEN=0
INTEGER(4),ALLOCATABLE :: IND(:)
INTEGER(4),ALLOCATABLE :: RANK(:)
INTEGER(4)             :: I0=0
INTEGER(4)             :: IHOLE=0
INTEGER(4)             :: RANK0
CONTAINS
!      ..................................................................
       SUBROUTINE HEAPSORT(LEN,X_,IND)
!      ******************************************************************
!      **                                                              **
!      **  HEAPSORT SORTING ALGORITHM                                  **
!      **  COMPUTES IND SUCH X_(IND(I)) INCREASES WITH INCREASING I    **
!      **                                                              **
!      **  SEE  "NUMERICAL RECIPES; THE ART OF SCIENTIFIC COMPUTING"   **
!      **        W.H. PRESS, B.P. FLANNERY, S.A TEUKOLSKY AND          **
!      **                                   W.T. VETTERLING            **
!      **        CAMBRIDGE UNIVERSITY PRESS 1986                       **
!      **                                                              **
!      **  1) DEFINE HIRARCHICAL STRUCTURE:                            **
!      **  EACH BOSS (I) IS CONNECTED TO TWO WORKERS (2I; 2I+1)        **
!      **  THE FIRST LEN/2 ELEMENTS ARE BOSSES                         **
!      **                                                              **
!      **  2) HIRING PHASE:                                            **
!      **  BEGINNING WITH THE "LOWEST" BOSS (LEN/2),                   **
!      **  EACH BOSS IS COMPARED TO THE LARGER OF HIS TWO WORKERS      **
!      **  IF THE BOSS IS SMALLER, WORKER AND BOSS TRADE PLACES        **
!      **                                                              **
!      **  THE FORMER BOSS HAS TO COMPETE NOW IMMEDIATELY              **
!      **  WITH THE WORKERS OF HIS NEW POSITION, UNTIL HE FOUND HIS    **
!      **  PLACE BEING LARGER THAN HIS WORKERS                         **
!      **                                                              **
!      **  3) RETIRING PHASE:                                          **
!      **                                                              **
!      **                                                              **
!      **                                                              **
!      ******************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)  :: LEN
       REAL(8)   ,INTENT(IN)  :: X_(LEN)
       INTEGER(4),INTENT(OUT) :: IND(LEN)
       REAL(8)                :: X(LEN)
       INTEGER(4)             :: L,IR,I,J
       INTEGER(4)             :: INDA
       REAL(8)                :: XA
!      ******************************************************************
       X(:)=X_(:)
       DO I=1,LEN
         IND(I)=I
       ENDDO
       L=LEN/2+1  ! POSITION OF LOWEST BOSS +1 (1 WILL BE SUBTRACTED LATER)
       IR=LEN
 10    CONTINUE
       IF(L.GT.1) THEN      ! HIRING PHASE
         L=L-1              ! SELECT BOSS FOR COMPETITION
         XA=X(L)            ! PLACE BOSS INTO XA
         INDA=IND(L)  
       ELSE 
         XA=X(IR)           ! CLEAR SPACE AT THE END OF THE ARRAY
         INDA=IND(IR)       
         X(IR)=X(1)         ! RETIRE BOSS OF ALL BOSSES INTO IT
         IND(IR)=IND(1)       
         IR=IR-1            ! DECREASE SIZE OF CORPORATION
         IF(IR.EQ.1) THEN   ! SHUT DOWN
           X(1)=XA
           IND(1)=INDA  
           RETURN
         END IF
       END IF
       I=L                           ! BOSS IS CALLED I AND IS LOCATED IN XA
       J=2*L                         ! SELECT FIRST OF HIS TWO WORKERS
       DO WHILE (J.LE.IR)              
         IF(J.LT.IR) THEN     
           IF(X(J).LT.X(J+1))J=J+1   !  SELECT THE BETTER WORKER
         END IF
         IF(XA.LT.X(J)) THEN         ! DECIDE ON DEMOTION
           X(I)=X(J)                 ! PROMOTE WORKER
           IND(I)=IND(J)
           I=J                       ! DEMOTE BOSS
           J=2*J                     ! SELECT ONE OF HIS NEW WORKERS
         ELSE
           EXIT                      ! TERMINATE DEMOTION
         END IF
       ENDDO
       X(I)=XA                       ! PLACE DEMOTED BOSS FROM XA IN THE NEW POSITION
       IND(I)=INDA  
       GOTO 10                       ! INITIATE NEW CYCLE
       END SUBROUTINE HEAPSORT
!
!      ..................................................................
       SUBROUTINE SORTRANK(LEN,IND,RANK)
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: LEN
       INTEGER(4),INTENT(IN) :: IND(LEN)
       INTEGER(4),INTENT(OUT):: RANK(LEN)
       INTEGER(4)            :: I
       DO I=1,LEN
         RANK(IND(I))=I
       ENDDO
       RETURN
       END SUBROUTINE SORTRANK
END MODULE SORT_MODULE
!
!     ..................................................................
      SUBROUTINE SORT$indexarray(lEN_,x,ind_)
!     **                                                              **
!     ** direct interface for heapsort                                **
!     **    x(ind(i)) increases with inreasing i                      **
!     **                                                              **
      USE SORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN_
      real(8)   ,intent(in) :: x(len_)
      integer(4),intent(out):: ind_(len_)
!     *******************************************************************
      call HEAPSORT(LEN_,x,IND_)
      return 
      end
!
!     ..................................................................
      SUBROUTINE SORT$SET(LEN_,CRIT)
!     ******************************************************************
!     ******************************************************************
      USE SORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN_
      REAL(8)   ,INTENT(IN) :: CRIT(LEN_)
      INTEGER(4)            :: ISVAR,I
!     ******************************************************************
      IF(TSET) THEN
        CALL ERROR$MSG('SORT OBJECT ALREADY IN USE')
        CALL ERROR$STOP('SORT$SET')
      END IF
      TSET=.TRUE.
      LEN=LEN_
if(len.le.1) then
  if(len.eq.1) len=0
  return
end if
      ALLOCATE(IND(LEN))
      ALLOCATE(RANK(LEN))
      CALL HEAPSORT(LEN,CRIT,IND)
!     WRITE(*,FMT='("CRIT",9F9.5)')CRIT
!     WRITE(*,FMT='("IND ",9I9)')IND
      DO I=1,LEN-1
        IF(IND(I+1)+1.EQ.IND(I)) THEN
          IF(ABS(CRIT(IND(I+1))-CRIT(IND(I))).LE.0) THEN
            ISVAR=IND(I)
            IND(I)=IND(I+1)
            IND(I+1)=ISVAR
          END IF
        END IF
      ENDDO
!     WRITE(*,FMT='("IND ",9I9)')IND
      CALL SORTRANK(LEN,IND,RANK)
      I0=1
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SORT$RESTART
!     ******************************************************************
!     ******************************************************************
      USE SORT_MODULE
      IMPLICIT NONE
!     ******************************************************************
      IF(.NOT.TSET) THEN
        CALL ERROR$MSG('SORT OBJECT NOT INITIALIZED')        
        CALL ERROR$STOP('SORT$RESTART')
      END IF
if(len.eq.0) return
      CALL SORTRANK(LEN,IND,RANK)
      I0=1
      IHOLE=0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SORT$UNSET
!     ******************************************************************
!     ******************************************************************
      USE SORT_MODULE
      IMPLICIT NONE
!     ******************************************************************
      TSET=.FALSE.
i0=0
if(len.eq.0)return
len=0
      DEALLOCATE(IND)
      DEALLOCATE(RANK)
      LEN=1
      I0=0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SORT$FLIP(FROM,TO)
!     ******************************************************************
!     **                                                              **
!     **  ORDER ARRAY ACCORDING TO AN ORDERING SCHEME DEFINED IN      **
!     **  SORT$SET.                                                   **
!     **  USE THE FOLLOWING LOOP:                                     **
!     **                                                              **
!     **     CALL SORT$RESTART                                        **
!     **     CALL SORT$FLIP(FROM,TO)                                  **
!     **     DO WHILE (FROM.NE.0.OR.TO.NE.0)                          **
!     **       IF(TO.EQ.0) THEN                                       **
!     **         SVAR=FIOFT(FROM,ISTEP)                               **
!     **       ELSE IF(FROM.EQ.0) THEN                                **
!     **         FIOFT(TO)=SVAR                                       **
!     **      ELSE                                                    **
!     **         FIOFT(TO)=FIOFT(FROM)                                **
!     **       END IF                                                 **
!     **       CALL SORT$FLIP(FROM,TO)                                **
!     **     ENDDO                                                    **
!     **                                                              **
!     ******************************************************************
      USE SORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: FROM
      INTEGER(4),INTENT(OUT) :: TO
!     ******************************************************************
      IF(.NOT.TSET) THEN
        CALL ERROR$MSG('SORT OBJECT NOT INITIALIZED')        
        CALL ERROR$STOP('SORT$FLIP')
      END IF
if(len.eq.0) then
  from=0
  to=0
  return
end if
      IF(IHOLE.EQ.0) THEN
!       == FIND ELEMENT WHICH IS ON THE WRONG POSITION ================
        DO WHILE(RANK(I0).EQ.I0)
          IF(I0.EQ.LEN) THEN
!           == FINISH SORTING =========
            FROM=0
            TO=0
            RETURN
          END IF
          I0=I0+1
        ENDDO
!       ==  COPY ELEMENT INTO BUFFER ===================================
        FROM=I0
        TO=0
        RANK0=RANK(FROM)
        RANK(FROM)=0
        IHOLE=FROM
      ELSE
        TO=IHOLE
        FROM=IND(TO)
        IF(FROM.EQ.I0) THEN
!         == COPY ELEMENT FROM BUFFER INTO THE HOLE ====================
          FROM=0
          RANK(TO)=RANK0
          RANK0=0
        ELSE
          RANK(TO)=RANK(FROM)
          RANK(FROM)=0
        END IF
        IHOLE=FROM       
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SORT$ORDERC16(LEN,NB,ARRAY)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LEN
      COMPLEX(8),INTENT(INOUT) :: ARRAY(LEN,NB)
      INTEGER(4)               :: FROM,TO
      COMPLEX(8)               :: TMPARRAY(LEN)
!     ******************************************************************
      CALL SORT$RESTART
      CALL SORT$FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          TMPARRAY(:)=ARRAY(:,FROM)
        ELSE IF(FROM.EQ.0) THEN
          ARRAY(:,TO)=TMPARRAY(:)
        ELSE
          ARRAY(:,TO)=ARRAY(:,FROM)
        END IF
        CALL SORT$FLIP(FROM,TO)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SORT$ORDERR8(LEN,NB,ARRAY)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LEN
      REAL(8)   ,INTENT(INOUT) :: ARRAY(LEN,NB)
      INTEGER(4)               :: FROM,TO
      REAL(8)                  :: TMPARRAY(LEN)
!     ******************************************************************
      CALL SORT$RESTART
      CALL SORT$FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          TMPARRAY(:)=ARRAY(:,FROM)
        ELSE IF(FROM.EQ.0) THEN
          ARRAY(:,TO)=TMPARRAY(:)
        ELSE
          ARRAY(:,TO)=ARRAY(:,FROM)
        END IF
        CALL SORT$FLIP(FROM,TO)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SORT$ORDERI4(LEN,NB,ARRAY)
!     ******************************************************************
!     ******************************************************************
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LEN
      INTEGER(4),INTENT(INOUT) :: ARRAY(LEN,NB)
      INTEGER(4)               :: FROM,TO
      INTEGER(4)               :: TMPARRAY(LEN)
!     ******************************************************************
      CALL SORT$RESTART
      CALL SORT$FLIP(FROM,TO)
      DO WHILE (FROM.NE.0.OR.TO.NE.0)
        IF(TO.EQ.0) THEN
          TMPARRAY(:)=ARRAY(:,FROM)
        ELSE IF(FROM.EQ.0) THEN
          ARRAY(:,TO)=TMPARRAY(:)
        ELSE
          ARRAY(:,TO)=ARRAY(:,FROM)
        END IF
        CALL SORT$FLIP(FROM,TO)
      ENDDO
      RETURN
      END
!      
!     ......................................................MADELUNG....
      SUBROUTINE MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
!     ******************************************************************
!     **                                                              **
!     ** EVALUATES MADELUNG ENERGY, POTENTIAL AND FORCES              **
!     **                                                              **
!     ** USES: MPE$QUERY                                              **
!     **       MPE$COMBINE                                            **
!     **       GBASS                                                  **
!     **       BOXSPH                                                 **
!     **                                                              **
!     ******************************************************************
!     USE MPE_MODULE
      IMPLICIT NONE
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8),   PARAMETER  :: TOL=1.D-8
      REAL(8),   INTENT(IN) :: RBAS(3,3)     ! LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: NBAS          ! NUMBER OF CHARGES
      REAL(8),   INTENT(IN) :: BAS(3,NBAS)   ! POSITIONS OF CHARGES
      REAL(8),   INTENT(IN) :: Q(NBAS)       ! CHARGES
      REAL(8),   INTENT(OUT):: EMAD          ! MADELUNG ENERGY   
      REAL(8),   INTENT(OUT):: VMAD(NBAS)    ! POTENTIAL AT CHARGE POSITIONS
      REAL(8),   INTENT(OUT):: FMAD(3,NBAS)  ! FORCES ON CHARGES
      REAL(8)               :: GBAS(3,3)     ! RECIPROCAL LATTICE VECTORS
      COMPLEX(8)            :: EIGR(NBAS)    ! FORM FACTOR
      COMPLEX(8)            :: EIGR12
      COMPLEX(8)            :: EIGR1
      REAL(8)               :: PI,FOURPI,ROOT2
      REAL(8)               :: VOL
      REAL(8)               :: RC,X,Y,C1,C2,SVAR,FAC,GFAC
      INTEGER(4)            :: K,I,IR,IR1,IR2
      REAL(8)               :: GMAX,RMAX
      INTEGER(4)            :: IG1MIN,IG1MAX
      INTEGER(4)            :: IG2MIN,IG2MAX
      INTEGER(4)            :: IG3MIN,IG3MAX
      INTEGER(4)            :: IG1,IG2,IG3
      REAL(8)               :: T1,T2,T3
      REAL(8)               :: G1,G2,G3,GSQUARE
      REAL(8)               :: SINFAC,COSFAC
      REAL(8)               :: Q12
      REAL(8)               :: DR1,DR2,DR3
      INTEGER(4)            :: IT1MIN,IT1MAX
      INTEGER(4)            :: IT2MIN,IT2MAX
      INTEGER(4)            :: IT3MIN,IT3MAX
      INTEGER(4)            :: IT1,IT2,IT3
      REAL(8)               :: DX,DY,DZ,DLEN
      REAL(8)               :: RFAC1,RFAC2,DV,QTOT
      REAL(8)               :: GR,G2MAX
      REAL(8)               :: ERFCX
!     INTEGER(4)            :: NTASKNUM,ITASK,ICOUNT
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      ROOT2=SQRT(2.D0)
!
!     ==================================================================
!     == CHECKS FOR PARALLEIZATION                                    ==
!     ==================================================================
!     CALL MPE$QUERY('NONE',NTASKNUM,ITASK)
!
!     ==================================================================
!     == CALCULATE RECIPROCAL TRANSLATION VECTORS                     ==
!     ==================================================================
      CALL GBASS(RBAS,GBAS,VOL)
!
!     ==================================================================
!     == CALCULATE RANGE OF R- AND G-SPACE SUMMATIONS                 ==
!     ==================================================================
      RC=1.D0
      DO K=1,10
        SVAR=PI*RC*TOL/2.D0**(1.5D0)
        DO I=1,10000
          X=0.1D0*DBLE(I)
          CALL LIB$ERFCR8(X,ERFCX)   !ERFC=1-ERF
          Y=PI/2.D0*ERFCX
!         PRINT*,'R: X=',X,' Y=',Y
          IF(Y.LT.SVAR) THEN
            C1=2.D0*X
            GOTO 100
          END IF
        ENDDO
        CALL ERROR$MSG('G-cutoff NOT FOUND')
        CALL ERROR$STOP('MADELUNG')
 100    CONTINUE
        SVAR=VOL*TOL/(8.D0*PI*RC**2)
        DX=0.1D0
        Y=0.D0
        DO I=1000,1,-1
          X=DX*DBLE(I)
          CALL LIB$ERFCR8(X,ERFCX)   !ERFC=1-ERF
          Y=Y+DX*X*ERFCX
!         PRINT*,'G: X=',X,' Y=',Y
          IF(Y.GT.SVAR) THEN
            X=X+DX          
            C2=X*SQRT(2.D0)
            GOTO 200
          END IF
        ENDDO
        CALL ERROR$MSG('R-cutoff NOT FOUND')
        CALL ERROR$STOP('MADELUNG')
 200    CONTINUE
        RC=SQRT(C1/(2.D0*PI*C2))*VOL**(1.D0/3.D0)
        GMAX=C1/RC
        RMAX=C2*RC
!       WRITE(*,FMT='(''RC='',F10.5,'' RMAX='',F10.5,'' GMAX='',F10.5)')
!    &       RC,RMAX,GMAX
      ENDDO
      GMAX=C1/RC*2.D0
      RMAX=C2*RC*2.D0
!     WRITE(*,FMT='(''RC='',F10.5,'' RMAX='',F10.5,'' GMAX='',F10.5)')
!    &       RC,RMAX,GMAX
!
!     ==================================================================
!     == DETERMINE CUTOFF RADII FOR REAL AND RECIPROCAL SPACE         ==
!     ==================================================================
      DO IR=1,NBAS
        VMAD(IR)=0.D0
        FMAD(1,IR)=0.D0
        FMAD(2,IR)=0.D0
        FMAD(3,IR)=0.D0
      ENDDO
!
!     ==================================================================
!     == G-SPACE SUM                                                  ==
!     ==================================================================
!     == G(R)=(SQRT(PI)*RC)**3 EXP(-(R/RC)**2)
!     == G(G)=1/V INT[V;DR:SUM G(R-T) EXP(-G*R)]=
!     ==     =INT[INFTY;DR:G(R)EXP(-G*R)]
!     ==     =1/V EXP(-(G*RC/2)**2) 
!     == E=Q1*Q2*SUM{G|G(G)**2*4*PI/G**2*COS(G*(R1-R2))               ==
!
      CALL BOXSPH(GBAS,0.D0,0.D0,0.D0,GMAX &
     &           ,IG1MIN,IG1MAX,IG2MIN,IG2MAX,IG3MIN,IG3MAX)
      G2MAX=GMAX**2
!     ICOUNT=0
      DO IG1=IG1MIN,IG1MAX
        T1=DBLE(IG1)
!       IF(IG1.EQ.0) THEN
!         IG2MIN=0
!       ELSE
!         IG2MIN=-IG2MAX
!       END IF  
        DO IG2=IG2MIN,IG2MAX
          T2=DBLE(IG2)
!         IF(IG2.EQ.0) THEN
!           IG3MIN=0
!         ELSE
!           IG3MIN=-IG3MAX
!         END IF  
          DO IG3=IG3MIN,IG3MAX
!           ICOUNT=ICOUNT+1
!           __ SELECTION FOR PARALLEL PROCESSING________________________
!           IF(MOD(ICOUNT-1,NTASKNUM).ne.ITASK-1) cycle
!
            T3=DBLE(IG3)
            G1=GBAS(1,1)*T1+GBAS(1,2)*T2+GBAS(1,3)*T3  
            G2=GBAS(2,1)*T1+GBAS(2,2)*T2+GBAS(2,3)*T3  
            G3=GBAS(3,1)*T1+GBAS(3,2)*T2+GBAS(3,3)*T3  
            GSQUARE=G1*G1+G2*G2+G3*G3
            IF(GSQUARE.LE.G2MAX.AND.GSQUARE.GT.1.D-7) THEN 
              FAC=2.D0*FOURPI/VOL*0.5D0
              SVAR=-0.5D0*GSQUARE*RC**2
              GFAC=FAC*EXP(SVAR)/GSQUARE
!             ========================================================
!             == this is the first time-critical part 
!             == can be streamlined:
!             ========================================================
              DO IR=1,NBAS
                GR=G1*BAS(1,IR)+G2*BAS(2,IR)+G3*BAS(3,IR) 
                EIGR(IR)=EXP(-CI*GR)
              ENDDO  
              DO IR1=1,NBAS
                EIGR1=CONJG(EIGR(IR1)) 
                SINFAC=0.D0
                COSFAC=0.D0
!               == better blas2: sinfac=aimag(eigr1*sum(eigr*q))
!               ==               cosfac= real(eigr1*sum(eigr*q))
                DO IR2=1,NBAS
                  EIGR12=EIGR1*EIGR(IR2)
                  SINFAC=SINFAC-AIMAG(EIGR12)*Q(IR2)
                  COSFAC=COSFAC+REAL(EIGR12,KIND=8)*Q(IR2)
                ENDDO
                SINFAC=SINFAC*GFAC
                COSFAC=COSFAC*GFAC 
                VMAD(IR1)=VMAD(IR1)+COSFAC
                SVAR=Q(IR1)*SINFAC           
                FMAD(1,IR1)=FMAD(1,IR1)+SVAR*G1
                FMAD(2,IR1)=FMAD(2,IR1)+SVAR*G2
                FMAD(3,IR1)=FMAD(3,IR1)+SVAR*G3
              ENDDO
            END IF
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == R-SPACE SUM                                                  ==
!     ==================================================================
!     ICOUNT=0
      FAC=1.D0/(ROOT2*RC)
      DO IR1=1,NBAS
!       ICOUNT=ICOUNT+1
!       __ SELECTION FOR PARALLEL PROCESSING__________________________
!       IF(MOD(ICOUNT-1,NTASKNUM).ne.ITASK-1) cycle
        DO IR2=1,NBAS
          Q12=Q(IR1)*Q(IR2)
          DR1=BAS(1,IR2)-BAS(1,IR1)
          DR2=BAS(2,IR2)-BAS(2,IR1)
          DR3=BAS(3,IR2)-BAS(3,IR1)
!         == calling boxsph is too complicated
          CALL BOXSPH(RBAS,-DR1,-DR2,-DR3,RMAX &
     &           ,IT1MIN,IT1MAX,IT2MIN,IT2MAX,IT3MIN,IT3MAX)
          DO IT1=IT1MIN,IT1MAX
            T1=DBLE(IT1)
            DO IT2=IT2MIN,IT2MAX
              T2=DBLE(IT2)
              DO IT3=IT3MIN,IT3MAX
                T3=DBLE(IT3)
                DX=DR1+RBAS(1,1)*T1+RBAS(1,2)*T2+RBAS(1,3)*T3  
                DY=DR2+RBAS(2,1)*T1+RBAS(2,2)*T2+RBAS(2,3)*T3  
                DZ=DR3+RBAS(3,1)*T1+RBAS(3,2)*T2+RBAS(3,3)*T3  
                DLEN=SQRT(DX*DX+DY*DY+DZ*DZ)
                IF(DLEN.LT.RMAX) THEN
!                 == this is time critical
                  IF(IR1.EQ.IR2 &
     &                   .AND.IT1.EQ.0.AND.IT2.EQ.0.AND.IT3.EQ.0) THEN
                    RFAC1=-SQRT(2.D0/PI)/RC
                    RFAC2=0.D0
                  ELSE
!                   == table lookup may be faster
!                   == store spline of p(x):=erfc(x)/x
!                   ==    rfac1=fac*p(dlen*fac)
!                   ==    rfac2=[ac**2*dp(dlen*fac)/d(dlen*fac)]/dlen
                    CALL LIB$ERFCR8(DLEN*FAC,ERFCX)
                    RFAC1=ERFCX/DLEN
                    RFAC2=-(RFAC1+FAC*2.D0/SQRT(PI) &
     &                               *EXP(-(FAC*DLEN)**2))/DLEN**2
                  END IF
                  RFAC1=0.5D0*RFAC1
                  RFAC2=0.5D0*RFAC2*Q12
                  VMAD(IR1)=VMAD(IR1)+RFAC1*Q(IR2)
                  VMAD(IR2)=VMAD(IR2)+RFAC1*Q(IR1)
                  FMAD(1,IR1)=FMAD(1,IR1)-RFAC2*DX     
                  FMAD(2,IR1)=FMAD(2,IR1)-RFAC2*DY     
                  FMAD(3,IR1)=FMAD(3,IR1)-RFAC2*DZ     
                  FMAD(1,IR2)=FMAD(1,IR2)+RFAC2*DX     
                  FMAD(2,IR2)=FMAD(2,IR2)+RFAC2*DY     
                  FMAD(3,IR2)=FMAD(3,IR2)+RFAC2*DZ     
                END IF              
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     CALL MPE$COMBINE('NONE','+',FMAD)
!     CALL MPE$COMBINE('NONE','+',VMAD)

!     == FROM YIN/COHEN???NIO
      QTOT=0.D0
      DO IR=1,NBAS
        QTOT=QTOT+Q(IR)
      ENDDO
      DV=-2.D0*RC**2*PI/VOL*QTOT
      DO IR=1,NBAS
        VMAD(IR)=VMAD(IR)+DV
      ENDDO
!
!
!     ==================================================================
!     == MADELUNG ENERGY                                              ==
!     ==================================================================
      EMAD=0.D0
      DO IR=1,NBAS
        EMAD=EMAD+0.5D0*Q(IR)*VMAD(IR)
      ENDDO
!
!     ==================================================================
!     == TURN DERIVATIVES INTO FORCES                                 ==
!     ==================================================================
      DO IR=1,NBAS
        DO I=1,3
          FMAD(I,IR)=-FMAD(I,IR)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      subroutine slaterkoster(r,ov,h)
!     ******************************************************************
!     ** slater-Koster energy integrals                               **
!     **  j.c. slater and G.F. Koster, Phys. Rev. 94, 1498 (1954)     **
!     **                                                              **
!     **  THE FIRST NINE REAL SPHERICAL HARMONICS ARE:                **
!     **      YLM(1)=SQRT( 1/( 4*PI))    * 1                          **
!     **      YLM(2)=SQRT( 3/( 4*PI))    * X / R                      **
!     **      YLM(3)=SQRT( 3/( 4*PI))    * Z / R                      **
!     **      YLM(4)=SQRT( 3/( 4*PI))    * Y / R                      **
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2      **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2      **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2      **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2      **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2      **
!     ******************************************************************
      implicit none
      real(8),intent(in)  :: r(3)
      real(8),intent(in)  :: ov(10)
      real(8),intent(out) :: h(9,9)
      real(8)             :: sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd
      real(8)             :: l,m,n,l2,m2,n2
      real(8)             :: svar,sq3,p1,p
      integer(4)          :: i,j
!     ******************************************************************
      sq3=sqrt(3.d0)
      svar=sqrt(dot_product(r,r))
      if(svar.lt.1.d-20) then
        stop 'in slaterkoster: distance=0'
      end if
      l=r(1)/svar
      m=r(2)/svar
      n=r(3)/svar
      l2=l**2
      m2=m**2
      n2=n**2
      sss=ov(1)
      sps=ov(2)
      pps=ov(3)
      ppp=ov(4)
      sds=ov(5)
      pds=ov(6)
      pdp=ov(7)
      dds=ov(8)
      ddp=ov(9)
      ddd=ov(10)
      h(:,:)=0.d0
!     == s-s block=================================================
      h(1,1)=sss
!     == s-p block=================================================
      h(1,2)=l*sps        
         h(1,3)=n*sps        
         h(1,4)=m*sps        
!     == p-p block=================================================
      h(2,2)=l2*pps+(1.d0-l2)*ppp
         h(3,3)=n2*pps+(1.d0-n2)*ppp
         h(4,4)=m2*pps+(1.d0-m2)*ppp
      h(2,4)=l*m*(pps-ppp)
      h(2,3)=l*n*(pps-ppp)
         h(3,4)=m*n*(pps-ppp)
!     == s-d block=================================================
      h(1,9)=sq3*l*m*sds
         h(1,8)=sq3*m*n*sds
         h(1,6)=sq3*n*l*sds
      h(1,5)=0.5d0*sq3*(l2-m2)*sds
      h(1,7)=(n2-0.5d0*(l2+m2))*sds
!     == p-d block=================================================
      h(2,9)=sq3*l2*m*pds+m*(1.d0-2.d0*l2)*pdp        
        h(4,8)=sq3*m2*n*pds+n*(1.d0-2.d0*m2)*pdp
        h(3,6)=sq3*n2*l*pds+l*(1.d0-2.d0*n2)*pdp
      h(2,8)=sq3*l*m*n*pds-2.d0*l*m*n*pdp
        h(4,6)=sq3*l*m*n*pds-2.d0*l*m*n*pdp
        h(3,9)=sq3*l*m*n*pds-2.d0*l*m*n*pdp
      h(2,6)=sq3*l2*n*pds+n*(1.d0-2.d0*l2)*pdp   !x,xz
        h(4,9)=sq3*m2*l*pds+l*(1.d0-2.d0*m2)*pdp     !y,yx
        h(3,8)=sq3*n2*m*pds+m*(1.d0-2.d0*n2)*pdp     !z,zy
      h(2,5)=0.5d0*sq3*l*(l2-m2)*pds+l*(1.d0-l2+m2)*pdp
      h(4,5)=0.5d0*sq3*m*(l2-m2)*pds-m*(1.d0+l2-m2)*pdp
      h(3,5)=0.5d0*sq3*n*(l2-m2)*pds-n*(l2-m2)*pdp
      h(2,7)=l*(n2-0.5d0*(l2+m2))*pds-sq3*l*n2*pdp
      h(4,7)=m*(n2-0.5d0*(l2+m2))*pds-sq3*m*n2*pdp
      h(3,7)=n*(n2-0.5d0*(l2+m2))*pds+sq3*n*(l2+m2)*pdp
!     == d-d block=================================================
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2      **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2      **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2      **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2      **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2      **
      h(9,9)=3.d0*l2*m2*dds+(l2+m2-4.d0*l2*m2)*ddp+(n2+l2*m2)*ddd   !xy,xy
      h(8,8)=3.d0*m2*n2*dds+(m2+n2-4.d0*m2*n2)*ddp+(l2+m2*n2)*ddd   !yz,yz cyc.perm
      h(6,6)=3.d0*n2*l2*dds+(n2+l2-4.d0*n2*l2)*ddp+(m2+n2*l2)*ddd   !zx,zx cyc.perm

      h(9,8)=3.d0*l*m2*n*dds+l*n*(1-4.d0*m2)*ddp+l*n*(m2-1.d0)*ddd  !xy,yz
      h(8,6)=3.d0*m*n2*l*dds+m*l*(1-4.d0*n2)*ddp+m*l*(n2-1.d0)*ddd  !yz,zx cyc.perm
!     h(6,9)=3.d0*n*l2*m*dds+n*m*(1-4.d0*l2)*ddp+n*m*(l2-1.d0)*ddd  

      h(9,6)=3.d0*l2*m*n*dds+m*n*(1-4.d0*l2)*ddp+m*n*(l2-1.d0)*ddd  !xy,xz
!     h(8,9)=3.d0*m2*n*l*dds+n*l*(1-4.d0*m2)*ddp+n*l*(m2-1.d0)*ddd  !
!     h(6,8)=3.d0*n2*l*m*dds+l*m*(1-4.d0*n2)*ddp+l*m*(n2-1.d0)*ddd  !
      h(9,5)=1.5d0*l*m*(l2-m2)*dds+2.d0*l*m*(m2-l2)*ddp &           !xy,x2-y2 
     &                            +0.5d0*l*m*(l2-m2)*ddd            
      h(8,5)=1.5d0*m*n*(l2-m2)*dds-m*n*(1.d0+2.d0*(l2-m2))*ddp &    !yz,x2-y2
     &      +m*n*(1.d0+0.5d0*(l2-m2))*ddd
      h(6,5)=1.5d0*n*l*(l2-m2)*dds+n*l*(1.d0-2.d0*(l2-m2))*ddp &    !xz,x2-y2
     &      -n*l*(1.d0-0.5d0*(l2-m2))*ddd
      h(9,7)=sq3*l*m*(n2-0.5d0*(l2+m2))*dds-2.d0*sq3*l*m*n2*ddp &   !xy,3z2-r2
     &      +0.5d0*sq3*l*m*(1.d0+n2)*ddd
      h(8,7)=sq3*m*n*(n2-0.5d0*(l2+m2))*dds+sq3*m*n*(l2+m2-n2)*ddp & !yz,3z2-r2
     &      -0.5d0*sq3*m*n*(l2+m2)*ddd
      h(6,7)=sq3*l*n*(n2-0.5d0*(l2+m2))*dds+sq3*l*n*(l2+m2-n2)*ddp & !xz,3z2-r2
     &      -0.5d0*sq3*l*n*(l2+m2)*ddd
      h(5,5)=0.75d0*(l2-m2)**2*dds+(l2+m2-(l2-m2)**2)*ddp &          !x2-y2,x2-y2
     &      +(n2+0.25d0*(l2-m2)**2)*ddd
      h(5,7)=0.5d0*sq3*(l2-m2)*(n2-0.5d0*(l2+m2))*dds &              !x2-y2,3z2-r2
     &      +sq3*n2*(m2-l2)*ddp+0.25d0*sq3*(1.d0+n2)*(l2-m2)*ddd
      h(7,7)=(n2-0.5d0*(l2+m2))**2*dds+3.d0*n2*(l2+m2)*ddp &         !3z2-r2,3z2-r2
     &      +0.75d0*(l2+m2)**2*ddd
!
!     =================================================================
!     == make h hermitean =============================================
!     =================================================================
      do i=1,9
        p1=(-1.d0)**int(sqrt(real(i)+1.d-6-1.d0))
        do j=i+1,9
          p=p1*(-1.d0)**int(sqrt(real(j)+1.d-6-1.d0))
          h(j,i)=p*h(i,j)+h(j,i)
          h(i,j)=p*h(j,i)
        enddo
      enddo
      return
      end
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!****                                                                         ****
!****   test routines                                                         ****
!****                                                                         ****
!****                                                                         ****
!*********************************************************************************
!*********************************************************************************
!*********************************************************************************
!
!     ...............................................................................
      subroutine test_madelung()
!     **                                                                        **
!     **  CALCULATE VARIOUS MADELUNG CONSTANTS AND COMPARE WITH                 **
!     **  "CONDENSED MATTER PHYSICS" BY M.P.MARDER                              **
!     **                                                                        **
      integer(4),parameter :: nbas=2
      real(8)              :: rbas(3,3)
      real(8)              :: bas(3,nbas)
      reAL(8)              :: Q(NBAS)
      reAL(8)              :: EMAD
      reAL(8)              :: VMAD(NBAS)
      reAL(8)              :: FMAD(3,NBAS)
      reAL(8)              :: D,RS,DET
      reAL(8)              :: PI
!     ******************************************************************************
      PI=4.D0*ATAN(1.D0)
      q(1)=1.d0
      q(2)=-1.d0
!
!     == sodium chloride: 1.74757 =============================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      bas(:,1) =(/0.0d0,0.0d0,0.0d0/)
      bas(:,2) =(/0.5d0,0.0d0,0.0d0/)
      d=sqrt(sum((bas(:,2)-bas(:,1))**2))
      call MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      write(*,fmt='("Nacl structure",t30,"e=",f10.5,"  dev=",f10.5)')-emad*d,-EMAD*D-1.74757
!
!     == CsCl structure: 1.76268 ==============================================
      RBAS(1,:)=(/1.0D0,0.0D0,0.0D0/)
      RBAS(2,:)=(/0.0D0,1.0D0,0.0D0/)
      RBAS(3,:)=(/0.0D0,0.0D0,1.0D0/)
      bas(:,1) =(/0.0d0,0.0d0,0.0d0/)
      bas(:,2) =(/0.5d0,0.5d0,0.5d0/)
      d=sqrt(sum((bas(:,2)-bas(:,1))**2))
      call MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      write(*,fmt='("cscl structure",t30,"e=",f10.5,"  dev=",f10.5)')-emad*d,-EMAD*D-1.76268
!
!     == Zns structure: 1.63806 ===============================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      bas(:,1) =(/0.0d0,0.0d0,0.0d0/)
      bas(:,2) =(/0.25d0,0.25d0,0.25d0/)
      d=sqrt(sum((bas(:,2)-bas(:,1))**2))
      call MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      write(*,fmt='("zns structure",t30,"e=",f10.5,"  dev=",f10.5)')-emad*d,-EMAD*D-1.63806
!
!     ========================================================================
!     == MADELUNG CONSTANTS FOR METALS                                      ==
!     == ATTENTION THE DEFINITION IS DIFFERENT!!!!                          ==
!     ========================================================================
!
!     == FCC: 1.79186 ========================================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      bas(:,1) =(/0.0d0,0.0d0,0.0d0/)
      DET=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)) 
      RS=(3.D0*DET/(4.D0*PI))**(1.D0/3.D0)
      call MADELUNG(1,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      write(*,fmt='("fcc structure",t30,"e=",f10.5,"  dev=",f10.5)')-2.D0*emad*RS,-2.D0*EMAD*RS-1.79186
!
!     == SIC: 1.76012 ========================================================
      RBAS(1,:)=(/1.0D0,0.0D0,0.0D0/)
      RBAS(2,:)=(/0.0D0,1.0D0,0.0D0/)
      RBAS(3,:)=(/0.0D0,0.0D0,1.0D0/)
      bas(:,1) =(/0.0d0,0.0d0,0.0d0/)
      DET=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)) 
      RS=(3.D0*DET/(4.D0*PI))**(1.D0/3.D0)
      call MADELUNG(1,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      write(*,fmt='("sic structure",t30,"e=",f10.5,"  dev=",f10.5)')-2.D0*emad*RS,-2.D0*EMAD*RS-1.76012
!
!     == DIAMOND: 1.67085 ===================================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      bas(:,1) =(/0.0d0,0.0d0,0.0d0/)
      bas(:,2) =(/0.25d0,0.25d0,0.25d0/)
      DET=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)) 
      DET=0.5D0*DET
      Q(:)=1.D0
      RS=(3.D0*DET/(4.D0*PI))**(1.D0/3.D0)
      call MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      write(*,fmt='("diamond structure",t30,"e=",f10.5,"  dev=",f10.5)')-emad*RS,-EMAD*RS-1.67085
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RANDOM_MSLNG(RAN)
!     **************************************************************************
!     ** MINIMAL STANDARD LINEAR CONGRUENTIAL RANDOM NUMBER GENERATOR         **
!     ** S.K.PARK AND K.W.MILLER, COMMUNICATIONS OF THE ACM, 31, 1192 (1988)  **
!     **                                                                      **
!     ** THIS VERSION ONLY WORKS IF HUGE(SEED).GE.2147483647                  **
!     **                                                                      **
!     ** IF IT WORKS CORRECTLY THE SEED IN STEP 10000 IS 1043618065           **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT):: RAN
      INTEGER(4),PARAMETER  :: M=2147483647  !=2_8**31-1_8
      INTEGER(4),PARAMETER  :: A=16807
      INTEGER(4),PARAMETER  :: Q=127773    !=INT(M/A)
      INTEGER(4),PARAMETER  :: R=2836      !=MOD(M,A)
      INTEGER(4),SAVE       :: SEED=1
      INTEGER(4)            :: HI,LO,TEST
!     **********************************************************
      HI=INT(SEED/Q)
      LO=MOD(SEED,Q)
      TEST=A*LO-R*HI
      IF(TEST.GT.0) THEN
        SEED=TEST
      ELSE
        SEED=TEST+M
      END IF
      RAN=REAL(SEED,8)/REAL(M,8)
      RETURN
      END
