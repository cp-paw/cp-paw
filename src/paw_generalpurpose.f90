!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ABCALPHABETAGAMMA(SWITCH,A,B,C,ALPHA,BETA,GAMMA,T)
!     **************************************************************************
!     **  CONVERTS THE LATTICE REPRESENTATION OF LENTHS OF AND ANGLES         **
!     **  BETWEEN LATTICE VECTORS                                             **
!     **                                                                      **
!     **  ATTENTION!! ANGLES ARE IN DEGREE I.E. IN UNITS OF  2*PI/360         **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)    :: SWITCH  ! ABC...-> T / T->ABC...
      REAL(8)   ,INTENT(INOUT) :: A,B,C   ! LENGTH OF LATTICE VECTORS
      REAL(8)   ,INTENT(INOUT) :: ALPHA,BETA,GAMMA !ANGLES BETWEEN LAT. VECTORS
      REAL(8)   ,INTENT(INOUT) :: T(3,3)  !LATTICE VECTORS  
      REAL(8)   ,PARAMETER     :: PI=4.D0*ATAN(1.D0)
      REAL(8)                  :: DEGREE
      REAL(8)                  :: COSA,COSB,COSG,SING
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TFAPOT(R,Z,V)
!     **************************************************************************
!     **  GENERALIZED THOMAS FERMI ATOMIC POTENTIAL                           **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: R   !RADIUS
      REAL(8),INTENT(IN)  :: Z   !ATOMIC NUMBER
      REAL(8),INTENT(OUT) :: V   ! POTENTIAL
      REAL(8),PARAMETER   :: BY3=1.D0/3.D0
      REAL(8)             :: B,X,XS,T
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GBASS(RBAS,GBAS,DET)
!     **************************************************************************
!     **                                                                      **
!     **  GENERATES RECIPROCAL LATTICE VECTORS G1,G2,G3                       **
!     **  FROM THE REAL SPACE LATTICE VECTORS                                 **
!     **          RI*GI=2*PI ; GI*GJ=0 FOR I.NE.J                             **
!     **  AND THE REAL SPACE UNIT CELL VOLUME                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: RBAS(3,3) ! REAL SPACE LATTIC VECTORS
      REAL(8),INTENT(OUT) :: GBAS(3,3) ! RECIPROCAL SPACE LATTICE VECTORS
      REAL(8),INTENT(OUT) :: DET       ! REAL SPACE UNIT CELL VOLUME
      REAL(8),PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)             :: FAC
      INTEGER(4)          :: I,J
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BOXSPH(RBAS,X0,Y0,Z0,RMAX &
     &                 ,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!     **************************************************************************
!     **                                                                      **
!     **  BOXSPH DESCRIBES A BOX AROUND AN SPHERE                             **
!     **  CENTERED AT (Z0,Y0,Z0) AND WITH RADIUS RMAX                         **
!     **  AND RETURNS THE MINIMUM AND MAXIMUM NUMBER OF DISPLACEMENTS         **
!     **  IN STEPS OF THAT CREATE POINTS WITHIN THE SPHERE                    **
!     **                                                                      **
!     **  |R-R0| =< RMAX        ONLY IF                                       **
!     **  MIN1 =< I1 =< MAX1; MIN2 =< I2 =< MAX2; MIN3 =< I3 =< MAX3          **
!     **  WHERE:  R(I) = RBAS(I,1)*I1 + RBAS(I,2)*I2 + RBAS(I,3)*I3           **
!     **  AND     R0=(X0,Y0,Z0)                                               **
!     **                                                                      **
!     **  INPUT :                                                             **
!     **    RBAS        DISPLACEMENT VECTORS                                  **
!     **    RMAX        RADIUS OF THE SPHERE                                  **
!     **    X0,Y0,Z0    CENTER OF THE SPHERE                                  **
!     **  OUTPUT :                                                            **
!     **    MIN1,MAX1,MIN2,MAX2,MIN3,MAX3     (SEE ABOVE)                     **
!     **                                                                      **
!     **  WARNING: THE DISPLACEMENTS MUST NOT BE SMALLER THAN -1.E+6          **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
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
!     **************************************************************************
      DO I=1,3
        T1(I)=RBAS(I,1)
        T2(I)=RBAS(I,2)
        T3(I)=RBAS(I,3)
      ENDDO
!     ==================================================================
!     ==  CALCULATE RECIPROCAL LATTICE VECTORS G1,G2,G3               ==
!     ==  NOTE THAT THE FACTOR 2PI IS DROPPED                         ==
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BOXBOX(RBAS,R0BOX,TBOX &
     &                   ,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
!     **************************************************************************
!     **  CIRCUMSCRIBES A BOX                                                 **
!     **  DEFINED BY THE LOWER LEFT FORDER CORNER R0BOX                       **
!     **  AND THREE EDGE VECTORS TBOX                                         **
!     **  BY ANOTHER BOX DEFINED THROUGH THE GRID TRANSLATION VECTORS         **
!     **  RBAS SUCH THAT ALL POINTS IN THE CIRCUMSCRIBED BOX                  **
!     **  FULFILL R=RBAS*X WITH XMIN<X<XMAX FOR ALL THREE COORDINATES         **
!     **************************************************************************
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
!     **************************************************************************
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
        DT(1)=REAL(I1,KIND=8)
        DO I2=0,1
          DT(2)=REAL(I2,KIND=8)
          DO I3=0,1
            DT(3)=REAL(I3,KIND=8)
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NOMOM(STRING,NAT,R1,R2,RMASS)
!     **************************************************************************
!     **                                                                      **
!     ** ROTATES AND TRANSLATES THE VECTORS R2 SO THAT THE ANGULAR            **
!     ** AND/OR TRANSLATIONAL MOMENTUM FROM R1 TO R2 VANISHES                 **
!     **                                                                      **
!     **  INPUT:                                                              **
!     **  STRING    'TR' SUBTRACT TRANSLATION AND ROTATION                    **
!     **            'T'  SUBTRACT TRANSLATION                                 **
!     **            'R'  SUBTRACT ROTATION                                    **
!     **  NAT       NUMBER OF ATOMS                                           **
!     **  R1        REFERENCE STRUCTURE                                       **
!     **  R2        STRUCTURE TO BE ROTATED                                   **
!     **  RMASS     ATOMIC MASSES                                             **
!     **                                                                      **
!     **  OUTPUT:                                                             **
!     **  R2        ROTATED STRUCTURE                                         **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **  PROGRAM WILL FAIL FOR ROTATIONS OF PI                               **
!     **  PROGRAM IS NOT PROTECTED FOR 1-D SYSTEM                             **
!     **                                                                      **
!     **************************************************************************
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
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE EULERANGLE(PHI,THETA,PSI,R)
!     **************************************************************************
!     ** ROTATION MATRIX FROM EULER ANGLES (SEE GOLDSTEIN)                    **
!     **************************************************************************
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
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ROTATIONMATRIX(PHI,R)
!     **************************************************************************
!     **  CONSTRUCTS ROTATION MATRIX FOR A GIVEN ANGLE VECTOR                 **
!     **  THE ANGLE IS |PHI| THE AXIS IS PHI/|PHI|                            **
!     **  THE ROTATION IS COUNTER CLOCKWISE                                   **
!     **  THE TRANSFORMED VECTOR IS XPRIME=R*X                                **
!     **************************************************************************
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
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ROTATION$ANGLETOMATRIX(PHI,TINV,U)
!     **************************************************************************
!     **  COMPOSES ROTATION MATRIX FROM AXIS AND ANGLE                        **
!     **  LENGTH OF PHI IS THE ANGLE, DIRECTION OF PHI IS THE AXIS            **
!     **  THE ANGLE IS COUNTERCLOCKWISE IF VIEWED AGAINST THE AXIS            **
!     **  SOURCE: HTTP://EN.WIKIPEDIA.ORG/WIKI/ROTATION_MATRIX                **
!     **                                                                      **
!     ********************************PETER E. BLOECHL, GOSLAR 2011 ************
      REAL(8)   ,INTENT(IN)  :: PHI(3)
      REAL(8)   ,INTENT(OUT) :: U(3,3)
      LOGICAL(4),INTENT(IN)  :: TINV ! INCLUDES INVERSION?
      REAL(8)                :: C,S
      REAL(8)                :: ANGLE,AXIS(3)
      INTEGER(4)             :: J
!     **************************************************************************
      ANGLE=SQRT(SUM(PHI**2))
      AXIS(:)=PHI(:)/ANGLE
      C=COS(ANGLE)
      S=SIN(ANGLE)
      DO J=1,3
        U(:,J)=AXIS(:)*AXIS(J)*(1.D0-C)
        U(J,J)=U(J,J)+C
      ENDDO
      U(1,2)=U(1,2)-AXIS(3)*S
      U(2,3)=U(2,3)-AXIS(1)*S
      U(3,1)=U(3,1)-AXIS(2)*S
      U(3,2)=U(3,2)+AXIS(1)*S
      U(2,1)=U(2,1)+AXIS(3)*S
      U(1,3)=U(1,3)+AXIS(2)*S
      IF(TINV)U=-U
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ROTATION$MATRIXTOANGLE(U,PHI,TINV)
!     **************************************************************************
!     **  EXTRACTS AXIS AND ANGLE FROM A ROTATION MATRIX                      **
!     **                                                                      **
!     ********************************PETER E. BLOECHL, GOSLAR 2011 ************
      REAL(8),INTENT(IN)     :: U(3,3)
      REAL(8),INTENT(OUT)    :: PHI(3)
      LOGICAL(4),INTENT(OUT) :: TINV
      REAL(8)                :: TRACE
      REAL(8)                :: DET
      REAL(8)                :: AXIS(3)
      REAL(8)                :: ROT(3,3)
      REAL(8)                :: ANGLE
!     **************************************************************************
      TRACE=U(1,1)+U(2,2)+U(3,3)
      DET=U(1,1)*(U(2,2)*U(3,3)-U(3,2)*U(2,3)) &
     &   +U(1,2)*(U(2,3)*U(3,1)-U(3,3)*U(2,1)) &
     &   +U(1,3)*(U(2,1)*U(3,2)-U(3,1)*U(2,2)) 
      TINV=DET.LT.0.D0
      IF(TINV) THEN
         ROT=-U
         TRACE=-TRACE
      ELSE
         ROT=U
      END IF
      IF(TRACE.EQ.3.D0) THEN
        ANGLE=0.D0
        AXIS(:)=0.D0
        AXIS(3)=1.D0
      ELSE IF(TRACE.EQ.-1.D0) THEN
        ANGLE=4.D0*ATAN(1.D0)  !=PI
        AXIS(1)=SQRT(0.5D0*(ROT(1,1)+1.D0))
        AXIS(2)=SQRT(0.5D0*(ROT(2,2)+1.D0))
        AXIS(3)=SQRT(0.5D0*(ROT(3,3)+1.D0))
      ELSE  
        ANGLE=ACOS(0.5D0*(TRACE-1.D0))
        AXIS(1)=ROT(3,2)-ROT(2,3)
        AXIS(2)=ROT(1,3)-ROT(3,1)
        AXIS(3)=ROT(2,1)-ROT(1,2)
        AXIS(:)=AXIS/SQRT(SUM(AXIS**2))
      END IF
      PHI=AXIS*ANGLE
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RGB$RAINBOW(XINV,R,G,B)
!     **************************************************************************
!     **  MAPS [0,1] ONTO THE RGB VALUES OF A RAINBOW                         **
!     **  XINV=0 MAPS ONTO VIOLET AND XINV=1 MAPS TO RED                      **
!     **  RAINBOW COLORS FROM HTTP://SIMPLE.WIKIPEDIA.ORG/WIKI/RAINBOW        **
!     **                                                                      **
!     ********************************PETER E. BLOECHL, GOSLAR 2011 ************
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: XINV  ! IN[0,1]
      INTEGER(4),INTENT(OUT) :: R
      INTEGER(4),INTENT(OUT) :: G
      INTEGER(4),INTENT(OUT) :: B
      INTEGER(4),PARAMETER   :: NC=7
      INTEGER(4)             :: COLOR(3,NC)
      REAL(8)                :: X1(NC)=(/0.D0,0.15D0,0.3D0,0.45D0,0.6D0 &
     &                                  ,0.75D0,0.9D0/)
      INTEGER(4)             :: I
      REAL(8)                :: SVAR
      REAL(8)                :: X
!     **************************************************************************
      X=1.D0-XINV
      COLOR(:,1)=(/255,0,0/)   !RED
      COLOR(:,2)=(/255,127,0/) !ORANGE
      COLOR(:,3)=(/255,255,0/) !YELLOW
      COLOR(:,4)=(/0,255,0/)   !GREEN
      COLOR(:,5)=(/255,0,0/)   !BLUE
      COLOR(:,6)=(/111,0,255/) !INDIGO
      COLOR(:,7)=(/143,0,255/) !VIOLET
!     ==  RED 255,0,0 ORANGE 255,127,0 YELLOW 255,255,0 GREEN 0,255,0 ==========
!     ==  BLUE 0,0,255 INDIGO 111,0,255 VIOLET 143,0,255
      R=COLOR(1,1)
      G=COLOR(2,1)
      B=COLOR(3,1)
      DO I=1,NC-1
        IF(X1(I).LE.X.AND.X.LE.X1(I+1)) THEN
          SVAR=(X-X1(I))/(X1(I+1)-X1(I))
          R=NINT(COLOR(1,I)*(1.D0-SVAR)+COLOR(1,I+1)*SVAR)
          G=NINT(COLOR(2,I)*(1.D0-SVAR)+COLOR(2,I+1)*SVAR)
          B=NINT(COLOR(3,I)*(1.D0-SVAR)+COLOR(3,I+1)*SVAR)
        END IF
      ENDDO
      IF(X.GT.X1(NC)) THEN
        R=COLOR(1,NC)
        G=COLOR(2,NC)
        B=COLOR(3,NC)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RGB$BLACKBODY(T,R,G,B)
!     **************************************************************************
!     **  RGB VALUES OF A BLACK BODY WITH TEMPERATURE T KELVIN                **
!     ** SOURCE: HTTP://WWW.PHYSICS.SFASU.EDU/ASTRO/COLOR/BLACKBODY.HTML      **
!     **         MITCHELL CHARITY                                             **
!     ** HEIGHT PROFILE FROM 1000-15000 K                                     **
!     **                                                                      **
!     ********************************PETER E. BLOECHL, GOSLAR 2011 ************
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: T
      INTEGER(4),INTENT(OUT) :: R
      INTEGER(4),INTENT(OUT) :: G
      INTEGER(4),INTENT(OUT) :: B
!     **************************************************************************
      R=NINT(148.D0+56100000.D0*T**(-1.5D0))
      G=NINT(100.04D0*LOG(T)-623.6D0)
      IF(T.GT.6500.D0) G=NINT(184.D0+35200000.D0*T**(-1.5D0))
      B=NINT(194.18D0*LOG(T)-1448.6D0)
      R=MIN(MAX(R,0),255)
      G=MIN(MAX(G,0),255)
      B=MIN(MAX(B,0),255)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
!     **************************************************************************
!     **                                                                      **
!     **  FINDS THE ZERO OF A MONOTONIC FUNCTION Y(X)                         **
!     **  BY DOUBLING STEP SIZE AND SUBSEQUENT BYSECTION.                     **
!     **  THE ROUTINE MUST BE CALLED IN A LOOP                                **
!     **  THAT STOPS WHEN CONVERGENCE IS OBTAINED                             **
!     **  AND SUPPLIES NEW FUNCTION VALUES Y0 FOR THE VALUE X0                **
!     **  THAT IS SUPPLIED BY THIS ROUTINE.                                   **
!     **  REMARK: CALL THIS ROUTINE ONCE BEFORE THE LOOP                      **
!     **                                                                      **
!     **   INITIALIZE BEFORE CALL WITH                                        **
!     **   ISTART=1 FOR MONOTONICALLY INCREASING OR WITH                      **
!     **   ISTART=-1 FOR MONOTONICALLY DECREASING FUNCTIONS:                  **
!     **   AND WITH X0 THE STARTING ARGUMENT AND WITH DX THE                  **
!     **   STARTING STEP SIZE                                                 **
!     **                                                                      **
!     **   DO NOT CHANGE THE VALUES FOR IBI,DX,XM,Y DURING THE                **
!     **   ITERATION                                                          **
!     **                                                                      **
!     **   X0=??                                                              **
!     **   DX=??                                                              **
!     **   CALL BISEC(1,IBI,X0,Y0,DX,XM,YM)                                   **
!     **   DO I=1,MAX                                                         **
!     **     CALCULATE Y0 FOR VALUE X0                                        **
!     **     CALL BISEC(0,IBI,X0,Y0,DX,XM,YM)                                 **
!     **     IF(ABS(Y0).LT.TOL) EXIT                                          **
!     **   ENDDO                                                              **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: ISTART  !=1 BEFORE ITERATION/ =0 OTHERWISE
      INTEGER(4),INTENT(INOUT) :: IBI     ! SWITCH BETWEEN EXPANSION AND CONTRACTION 
      REAL(8)   ,INTENT(INOUT) :: X0      ! CURRENT ARGUMENT
      REAL(8)   ,INTENT(IN)    :: Y0      ! FUNCTION VALUE AT X0
      REAL(8)   ,INTENT(INOUT) :: DX      ! STEP WIDTH
      REAL(8)   ,INTENT(INOUT) :: YM      ! FUNCTION VALUE AT XM
      REAL(8)   ,INTENT(INOUT) :: XM      ! PREVIOUS ARGUMENT
      REAL(8)   ,SAVE          :: SLOPE
      REAL(8)                  :: XP
!     **************************************************************************
!
!     ==   STARTUP
      IF(ISTART.NE.0) THEN
        SLOPE=-1.D0
        IF(ISTART.GT.0) SLOPE=1.D0
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
!........1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LBFGS(N,NMEMX,NMEM,AMIX,YK,SK,X0,XM,F0,FM,XP)
!     **************************************************************************
!     ** LIMITED-MEMORY BROYDEN-FLETSCHER-GOLDFARB-SHANNO  (L-BFGS) ALGORITHM **
!     ** FOR OPTIMIZING NEARLY QUADRATIC FUNCTIONS
!     ** 
!     ** SEE P 779 OF UPDATING QUASI-NEWTON MATRICES IWTH LIMITED STORAGE,    **
!     ** JORGE NOCEDAL, MATHEMATICS OF COMPUTATION, 35. P773 (1980)           **
!     **                                                                      **
!     ** INITIALIZE NMEM=0 BEFORE THE FIRST CALL OR TO RESTART HISTORY.       **
!     ** DURING OPTIMIZATION LEAVE NMEM,YK,SK UNTOUCHED. THEY ARE UPDATED     **
!     ** INSIDE THIS ROUTINE.                                                 **
!     ** FOLLOWING THE ROUTINE,                                               **
!     **  (1) SHIFT POSITIONS AND FORCES, I.E. XM=X0,FM=X0,X0=XP              **
!     **  (2) PERFORM A LINE SEARCH ALONG X(S)=X0+(X0-XM)*S                   **
!     **      UNTIL F(S)*(X0-XM)=0                                            **
!     **                                                                      **
!     **************************************PETER BLOECHL GOSLAR 2015***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N     ! DIMENSION OD THE SEARCH SPACE
      INTEGER(4),INTENT(IN)    :: NMEMX ! MAX #(STEPS STORED IN HISTORY)
      INTEGER(4),INTENT(INOUT) :: NMEM  ! ACTUAL #(STEPS STORED IN HISTORY)
      REAL(8)   ,INTENT(INOUT) :: YK(N,NMEMX)
      REAL(8)   ,INTENT(INOUT) :: SK(N,NMEMX)
      REAL(8)   ,INTENT(IN)    :: AMIX  ! MIXING FACTOR
      REAL(8)   ,INTENT(IN)    :: X0(N) ! ACTUAL COORDINATES
      REAL(8)   ,INTENT(IN)    :: XM(N) ! PREVIOUS COORDINATES
      REAL(8)   ,INTENT(IN)    :: F0(N) ! ACTUAL FORCE
      REAL(8)   ,INTENT(IN)    :: FM(N) ! PREVIOUS FORCE
      REAL(8)   ,INTENT(OUT)   :: XP(N) ! NEXT COORDINATES
      REAL(8)                  :: ALPHA(NMEMX)
      REAL(8)                  :: BETA(NMEMX)
      REAL(8)                  :: RHO(NMEMX)
      REAL(8)                  :: Q(N)
      INTEGER(4)               :: I
!     **************************************************************************
!
!     ==========================================================================
!     ==  UPDATE MEMORY                                                       ==
!     ==========================================================================
      Q(:)=X0(:)-XM(:)
      IF(DOT_PRODUCT(Q,Q).GT.1.D-10) THEN
        NMEM=MIN(NMEMX,NMEM+1)
        NMEM=MIN(N,NMEM)
        DO I=NMEM-1,1,-1
          SK(:,I+1)=SK(:,I)
          YK(:,I+1)=YK(:,I)
        ENDDO
        SK(:,1)=X0(:)-XM(:)
        YK(:,1)=F0(:)-FM(:)
      END IF
!
!     ==========================================================================
!     ==  PREDICT NEXT POSITION                                               ==
!     ==========================================================================
      Q(:)=F0(:)
      DO I=1,NMEM
        RHO(I)=1.D0/DOT_PRODUCT(YK(:,I),SK(:,I))
        ALPHA(I)=RHO(I)*DOT_PRODUCT(SK(:,I),Q)
        Q(:)=Q(:)-YK(:,I)*ALPHA(I)
      ENDDO
      Q=-AMIX*Q
      DO I=NMEM,1,-1
        BETA(I)=RHO(I)*DOT_PRODUCT(YK(:,I),Q)
        Q(:)=Q(:)+SK(:,I)*(ALPHA(I)-BETA(I))
      ENDDO
      XP(:)=X0(:)-Q(:)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: RINT
      INTEGER(4)            :: K,IFAC,I
!     ******************************************************************
!     ==================================================================
!     ==  CALCULATE INT(DR): R**(2*L+2) *EXP(-ALPHA * R**2)           ==
!     ==  SEE BRONSTEIN P66.                                          ==
!     ==================================================================
      K=L+1
      IFAC=1
      DO I=2,K
        IFAC=IFAC*(2*I-1)
      ENDDO
      RINT=SQRT(PI/(4.D0*ALPHA))*REAL(IFAC,KIND=8)/(2.D0*ALPHA)**K
      C=1.D0/RINT
      RETURN
      END
!
      SUBROUTINE CG$TEST()
!     **                                                              **
!     ** TEST ROUTINE FOR CONJUGATE GRADIENT                          **
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4), PARAMETER :: N=2
      REAL(8)               :: R(N),F1(N),F2(N),D1(N),D2(N)
      REAL(8)               :: E
      INTEGER(4),PARAMETER  :: NITER=10
      INTEGER(4)            :: ITER,INNER
      REAL(8)               :: LAMBDA
      LOGICAL(4)            :: TCONV
!     ******************************************************************
      R=(/10.D0,2.D0/)
      CALL ETOT(N,R,E,F1)
WRITE(*,FMT='(I5," E ",F20.15," R=",2F10.5," F= ",2E10.3)')0,E,R,F1
      D1=F1
      DO ITER=1,NITER
        LAMBDA=1.D-2/SQRT(DOT_PRODUCT(F1,F1))
        DO INNER=1,100
          CALL ETOT(N,R+LAMBDA*D1,E,F2)
          CALL CG$LINESEARCH(N,F1,D1,F2,LAMBDA,TCONV)
          IF(TCONV) EXIT
        ENDDO
        IF(.NOT.TCONV) STOP 'NOT CONVERGED'
        R=R+D1*LAMBDA
WRITE(*,FMT='(I5," E ",F20.15," R=",2F10.5," F= ",2E10.3)')ITER,E,R,F2
        CALL CG$NEWDIR(N,F1,D1,F2,D2)
        F1=F2
        D1=D2
      ENDDO         
      STOP
      CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE ETOT(N,R,E,F)
        INTEGER(4),INTENT(IN) :: N
        REAL(8)   ,INTENT(IN) :: R(N)  
        REAL(8)   ,INTENT(OUT):: E
        REAL(8)   ,INTENT(OUT):: F(N)  
        INTEGER(4),PARAMETER  :: NLOC=2
        REAL(8)   ,PARAMETER  :: C=2.D0
        REAL(8)   ,PARAMETER  :: B(NLOC)=(/0.D0,0.D0/)
        REAL(8)               :: A(NLOC,NLOC)
!       ************************************************************************
        A(:,1)=(/2.D+1,0.D0/)
        A(:,2)=(/0.D0,2.D-3/)
        IF(N.NE.NLOC) STOP 'ERROR IN ETOT'
        E=C-DOT_PRODUCT(B,R)+0.5D0*DOT_PRODUCT(R,MATMUL(A,R))
        F=B-MATMUL(A,R)
        RETURN
        END SUBROUTINE ETOT
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CG$LINESEARCH(N,F1,D1,F2,DLAMBDA,TCONV)
!     **                                                              **
!     ** CONJUGATE GRADIENT LINE SEARCH                               **
!     **                                                              **
!     **  ADJUSTS LAMBDA SUCH THAT THE FORCE AT X+D1*LAMBDA           **
!     **  PARALLEL TO THE SEARCH DIRECTION D1 CONVERGES TO ZER        **
!     **                                                              **
!     **  X(LAMBDA)=X1+D1*LAMBDA                                      **
!     **  F(LAMBDA)=F1+LAMBDA (F2-F1)/LAMBDA_IN                       **
!     **  F2=F(LAMBDA_IN)                                             **
!     **  F1=F(LAMBDA=0)                                              **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      REAL(8)   ,INTENT(IN)   :: F1(N)  ! FORCE AT X1
      REAL(8)   ,INTENT(IN)   :: D1(N)  ! X(LAMBDA)=X1+LAMBDA*D1
      REAL(8)   ,INTENT(IN)   :: F2(N)  ! F(LAMBDA) WITH F(LAMBDA)*D1=0
      REAL(8)   ,INTENT(INOUT):: DLAMBDA  ! NEXT DIRECTION FOR LINE SEARCH
      LOGICAL(4),INTENT(OUT)  :: TCONV
      REAL(8)                 :: SVAR1,SVAR2
!     ******************************************************************
      SVAR1=DOT_PRODUCT(D1,F2)
      SVAR2=DOT_PRODUCT(D1,F2-F1)
      IF(DLAMBDA*SVAR2.LT.0.D0) THEN
        DLAMBDA=-SVAR1/SVAR2*DLAMBDA
      ELSE
PRINT*,'SVAR2 ',SVAR2,DLAMBDA
PRINT*,'WARNING! HESSIAN NOT POSITIVE DEFINITE; SWITCH TO STEPPING'
!       == CORRECTION FOR WRONG CURVATURE
        IF(SVAR2.GT.0.D0) THEN
          DLAMBDA=1.D-2/SQRT(DOT_PRODUCT(D1,D1))
        ELSE
          DLAMBDA=-1.D-2/SQRT(DOT_PRODUCT(D1,D1))
        END IF
      END IF
      TCONV=(ABS(SVAR2/SVAR1).LT.1.D-4)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CG$LINESEARCH_OLD(N,F1,D1,F2,LAMBDA,TCONV)
!     **************************************************************************
!     ** CONJUGATE GRADIENT LINE SEARCH                                       **
!     **                                                                      **
!     **  ADJUSTS LAMBDA SUCH THAT THE FORCE AT X+D1*LAMBDA                   **
!     **  PARALLEL TO THE SEARCH DIRECTION D1 CONVERGES TO ZER                **
!     **                                                                      **
!     **  X(LAMBDA)=X1+D1*LAMBDA                                              **
!     **  F(LAMBDA)=F1+LAMBDA (F2-F1)/LAMBDA_IN                               **
!     **  F2=F(LAMBDA_IN)                                                     **
!     **  F1=F(LAMBDA=0)                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      REAL(8)   ,INTENT(IN)   :: F1(N)  ! FORCE AT X1
      REAL(8)   ,INTENT(IN)   :: D1(N)  ! X(LAMBDA)=X1+LAMBDA*D1
      REAL(8)   ,INTENT(IN)   :: F2(N)  ! F(LAMBDA) WITH F(LAMBDA)*D1=0
      REAL(8)   ,INTENT(INOUT):: LAMBDA  ! NEXT DIRECTION FOR LINE SEARCH
      LOGICAL(4),INTENT(OUT)  :: TCONV
      REAL(8)                 :: SVAR1,SVAR2
!     ******************************************************************
      SVAR1=DOT_PRODUCT(D1,F1)
      SVAR2=DOT_PRODUCT(D1,F2-F1)
      IF(LAMBDA*SVAR2.LT.0.D0) THEN
        LAMBDA=-SVAR1/SVAR2*LAMBDA
      ELSE
PRINT*,'WARNING! HESSIAN NOT POSITIVE DEFINITE; SWITCH TO STEPPING'
!       == CORRECTION FOR WRONG CURVATURE
        IF(SVAR1+SVAR2.GT.0.D0) THEN
          LAMBDA=LAMBDA+1.D-2/SQRT(DOT_PRODUCT(D1,D1))
        ELSE
          LAMBDA=LAMBDA-1.D-2/SQRT(DOT_PRODUCT(D1,D1))
        END IF
      END IF
      TCONV=(ABS((SVAR2+SVAR1)/SVAR1).LT.1.D-4)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CG$NEWDIR(N,F1,D1,F2,D2)
!     **                                                              **
!     ** CONJUGATE GRADIENT NEW SEARCH DIRECTION                      **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: F1(N)  ! FORCE AT X1
      REAL(8)   ,INTENT(IN) :: D1(N)  ! X(LAMBDA)=X1+LAMBDA*D1
      REAL(8)   ,INTENT(IN) :: F2(N)  ! F(LAMBDA) WITH F(LAMBDA)*D1=0
      REAL(8)   ,INTENT(OUT):: D2(N)  ! NEXT DIRECTION FOR LINE SEARCH
!     ******************************************************************
      D2=F2+D1*DOT_PRODUCT(F2-F1,F2)/DOT_PRODUCT(F1,F1)
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
         HARVEST=HARVEST*SQRT(6.D0/REAL(N,KIND=8))
         RETURN 
       END IF 
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SORT$INDEXARRAY(LEN_,X,IND_)
!     **************************************************************************
!     ** DIRECT INTERFACE FOR HEAPSORT                                        **
!     **    X(IND(I)) INCREASES WITH INREASING I                              **
!     **************************************************************************
      USE SORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN_
      REAL(8)   ,INTENT(IN) :: X(LEN_)
      INTEGER(4),INTENT(OUT):: IND_(LEN_)
!     **************************************************************************
      CALL HEAPSORT(LEN_,X,IND_)
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
IF(LEN.LE.1) THEN
  IF(LEN.EQ.1) LEN=0
  RETURN
END IF
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
IF(LEN.EQ.0) RETURN
      CALL SORTRANK(LEN,IND,RANK)
      I0=1
      IHOLE=0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SORT$UNSET
!     ******************************************************************
!     ******************************************************************
      USE SORT_MODULE
      IMPLICIT NONE
!     ******************************************************************
      TSET=.FALSE.
I0=0
IF(LEN.EQ.0)RETURN
LEN=0
      DEALLOCATE(IND)
      DEALLOCATE(RANK)
      LEN=1
      I0=0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SORT$FLIP(FROM,TO)
!     **************************************************************************
!     **  ORDER ARRAY ACCORDING TO AN ORDERING SCHEME DEFINED IN SORT$SET.    **
!     **  USE THE FOLLOWING LOOP:                                             **
!     **                                                                      **
!     **     CALL SORT$RESTART                                                **
!     **     CALL SORT$FLIP(FROM,TO)                                          **
!     **     DO WHILE (FROM.NE.0.OR.TO.NE.0)                                  **
!     **       IF(TO.EQ.0) THEN                                               **
!     **         SVAR=FIOFT(FROM,ISTEP)                                       **
!     **       ELSE IF(FROM.EQ.0) THEN                                        **
!     **         FIOFT(TO)=SVAR                                               **
!     **      ELSE                                                            **
!     **         FIOFT(TO)=FIOFT(FROM)                                        **
!     **       END IF                                                         **
!     **       CALL SORT$FLIP(FROM,TO)                                        **
!     **     ENDDO                                                            **
!     **************************************************************************
      USE SORT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: FROM
      INTEGER(4),INTENT(OUT) :: TO
!     **************************************************************************
      IF(.NOT.TSET) THEN
        CALL ERROR$MSG('SORT OBJECT NOT INITIALIZED')        
        CALL ERROR$STOP('SORT$FLIP')
      END IF
IF(LEN.EQ.0) THEN
  FROM=0
  TO=0
  RETURN
END IF
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SORT$ORDERC16(LEN,NB,ARRAY)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LEN
      COMPLEX(8),INTENT(INOUT) :: ARRAY(LEN,NB)
      INTEGER(4)               :: FROM,TO
      COMPLEX(8)               :: TMPARRAY(LEN)
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SORT$ORDERR8(LEN,NB,ARRAY)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LEN
      REAL(8)   ,INTENT(INOUT) :: ARRAY(LEN,NB)
      INTEGER(4)               :: FROM,TO
      REAL(8)                  :: TMPARRAY(LEN)
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SORT$ORDERI4(LEN,NB,ARRAY)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LEN
      INTEGER(4),INTENT(INOUT) :: ARRAY(LEN,NB)
      INTEGER(4)               :: FROM,TO
      INTEGER(4)               :: TMPARRAY(LEN)
!     **************************************************************************
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
!     **************************************************************************
!     ** EVALUATES MADELUNG ENERGY, POTENTIAL AND FORCES                      **
!     **                                                                      **
!     ** USES: MPE$QUERY                                                      **
!     **       MPE$COMBINE                                                    **
!     **       GBASS                                                          **
!     **       BOXSPH                                                         **
!     **                                                                      **
!     **************************************************************************
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
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: GBAS(3,3)     ! RECIPROCAL LATTICE VECTORS
      COMPLEX(8)            :: EIGR(NBAS)    ! FORM FACTOR
      COMPLEX(8)            :: EIGR12
      COMPLEX(8)            :: EIGR1
      REAL(8)               :: FOURPI,ROOT2
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
!     **************************************************************************
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
          X=0.1D0*REAL(I,KIND=8)
          CALL LIB$ERFCR8(X,ERFCX)   !ERFC=1-ERF
          Y=PI/2.D0*ERFCX
!         PRINT*,'R: X=',X,' Y=',Y
          IF(Y.LT.SVAR) THEN
            C1=2.D0*X
            GOTO 100
          END IF
        ENDDO
        CALL ERROR$MSG('G-CUTOFF NOT FOUND')
        CALL ERROR$STOP('MADELUNG')
 100    CONTINUE
        SVAR=VOL*TOL/(8.D0*PI*RC**2)
        DX=0.1D0
        Y=0.D0
        DO I=1000,1,-1
          X=DX*REAL(I,KIND=8)
          CALL LIB$ERFCR8(X,ERFCX)   !ERFC=1-ERF
          Y=Y+DX*X*ERFCX
!         PRINT*,'G: X=',X,' Y=',Y
          IF(Y.GT.SVAR) THEN
            X=X+DX          
            C2=X*SQRT(2.D0)
            GOTO 200
          END IF
        ENDDO
        CALL ERROR$MSG('R-CUTOFF NOT FOUND')
        CALL ERROR$STOP('MADELUNG')
 200    CONTINUE
        RC=SQRT(C1/(2.D0*PI*C2))*VOL**(1.D0/3.D0)
        GMAX=C1/RC
        RMAX=C2*RC
!       WRITE(*,FMT='("RC=",F10.5," RMAX=",F10.5," GMAX=",F10.5)')
!    &       RC,RMAX,GMAX
      ENDDO
      GMAX=C1/RC*2.D0
      RMAX=C2*RC*2.D0
!     WRITE(*,FMT='("RC=",F10.5," RMAX=",F10.5," GMAX=",F10.5)')
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
        T1=REAL(IG1,KIND=8)
!       IF(IG1.EQ.0) THEN
!         IG2MIN=0
!       ELSE
!         IG2MIN=-IG2MAX
!       END IF  
        DO IG2=IG2MIN,IG2MAX
          T2=REAL(IG2,KIND=8)
!         IF(IG2.EQ.0) THEN
!           IG3MIN=0
!         ELSE
!           IG3MIN=-IG3MAX
!         END IF  
          DO IG3=IG3MIN,IG3MAX
!           ICOUNT=ICOUNT+1
!           __ SELECTION FOR PARALLEL PROCESSING________________________
!           IF(MOD(ICOUNT-1,NTASKNUM).NE.ITASK-1) CYCLE
!
            T3=REAL(IG3,KIND=8)
            G1=GBAS(1,1)*T1+GBAS(1,2)*T2+GBAS(1,3)*T3  
            G2=GBAS(2,1)*T1+GBAS(2,2)*T2+GBAS(2,3)*T3  
            G3=GBAS(3,1)*T1+GBAS(3,2)*T2+GBAS(3,3)*T3  
            GSQUARE=G1*G1+G2*G2+G3*G3
            IF(GSQUARE.LE.G2MAX.AND.GSQUARE.GT.1.D-7) THEN 
              FAC=2.D0*FOURPI/VOL*0.5D0
              SVAR=-0.5D0*GSQUARE*RC**2
              GFAC=FAC*EXP(SVAR)/GSQUARE
!             ========================================================
!             == THIS IS THE FIRST TIME-CRITICAL PART 
!             == CAN BE STREAMLINED:
!             ========================================================
              DO IR=1,NBAS
                GR=G1*BAS(1,IR)+G2*BAS(2,IR)+G3*BAS(3,IR) 
                EIGR(IR)=EXP(-CI*GR)
              ENDDO  
              DO IR1=1,NBAS
                EIGR1=CONJG(EIGR(IR1)) 
                SINFAC=0.D0
                COSFAC=0.D0
!               == BETTER BLAS2: SINFAC=AIMAG(EIGR1*SUM(EIGR*Q))
!               ==               COSFAC= REAL(EIGR1*SUM(EIGR*Q))
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
!       IF(MOD(ICOUNT-1,NTASKNUM).NE.ITASK-1) CYCLE
        DO IR2=1,NBAS
          Q12=Q(IR1)*Q(IR2)
          DR1=BAS(1,IR2)-BAS(1,IR1)
          DR2=BAS(2,IR2)-BAS(2,IR1)
          DR3=BAS(3,IR2)-BAS(3,IR1)
!         == CALLING BOXSPH IS TOO COMPLICATED
          CALL BOXSPH(RBAS,-DR1,-DR2,-DR3,RMAX &
     &           ,IT1MIN,IT1MAX,IT2MIN,IT2MAX,IT3MIN,IT3MAX)
          DO IT1=IT1MIN,IT1MAX
            T1=REAL(IT1,KIND=8)
            DO IT2=IT2MIN,IT2MAX
              T2=REAL(IT2,KIND=8)
              DO IT3=IT3MIN,IT3MAX
                T3=REAL(IT3,KIND=8)
                DX=DR1+RBAS(1,1)*T1+RBAS(1,2)*T2+RBAS(1,3)*T3  
                DY=DR2+RBAS(2,1)*T1+RBAS(2,2)*T2+RBAS(2,3)*T3  
                DZ=DR3+RBAS(3,1)*T1+RBAS(3,2)*T2+RBAS(3,3)*T3  
                DLEN=SQRT(DX*DX+DY*DY+DZ*DZ)
                IF(DLEN.LT.RMAX) THEN
!                 == THIS IS TIME CRITICAL
                  IF(IR1.EQ.IR2 &
     &                   .AND.IT1.EQ.0.AND.IT2.EQ.0.AND.IT3.EQ.0) THEN
                    RFAC1=-SQRT(2.D0/PI)/RC
                    RFAC2=0.D0
                  ELSE
!                   == TABLE LOOKUP MAY BE FASTER
!                   == STORE SPLINE OF P(X):=ERFC(X)/X
!                   ==    RFAC1=FAC*P(DLEN*FAC)
!                   ==    RFAC2=[AC**2*DP(DLEN*FAC)/D(DLEN*FAC)]/DLEN
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SLATERKOSTER(R,OV,H)
!     **************************************************************************
!     ** SLATER-KOSTER ENERGY INTEGRALS                                       **
!     **  J.C. SLATER AND G.F. KOSTER, PHYS. REV. 94, 1498 (1954)             **
!     **                                                                      **
!     **  THE FIRST NINE REAL SPHERICAL HARMONICS ARE:                        **
!     **      YLM(1)=SQRT( 1/( 4*PI))    * 1                                  **
!     **      YLM(2)=SQRT( 3/( 4*PI))    * X / R                              **
!     **      YLM(3)=SQRT( 3/( 4*PI))    * Z / R                              **
!     **      YLM(4)=SQRT( 3/( 4*PI))    * Y / R                              **
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2              **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2              **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2              **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2              **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: R(3)
      REAL(8),INTENT(IN)  :: OV(10)
      REAL(8),INTENT(OUT) :: H(9,9)
      REAL(8)             :: SSS,SPS,PPS,PPP,SDS,PDS,PDP,DDS,DDP,DDD
      REAL(8)             :: L,M,N,L2,M2,N2
      REAL(8)             :: SVAR,SQ3,P1,P
      INTEGER(4)          :: I,J
!     **************************************************************************
      SQ3=SQRT(3.D0)
      SVAR=SQRT(DOT_PRODUCT(R,R))
      IF(SVAR.LT.1.D-20) THEN
        STOP 'IN SLATERKOSTER: DISTANCE=0'
      END IF
      L=R(1)/SVAR
      M=R(2)/SVAR
      N=R(3)/SVAR
      L2=L**2
      M2=M**2
      N2=N**2
      SSS=OV(1)
      SPS=OV(2)
      PPS=OV(3)
      PPP=OV(4)
      SDS=OV(5)
      PDS=OV(6)
      PDP=OV(7)
      DDS=OV(8)
      DDP=OV(9)
      DDD=OV(10)
      H(:,:)=0.D0
!     == S-S BLOCK=================================================
      H(1,1)=SSS
!     == S-P BLOCK=================================================
      H(1,2)=L*SPS        
         H(1,3)=N*SPS        
         H(1,4)=M*SPS        
!     == P-P BLOCK=================================================
      H(2,2)=L2*PPS+(1.D0-L2)*PPP
         H(3,3)=N2*PPS+(1.D0-N2)*PPP
         H(4,4)=M2*PPS+(1.D0-M2)*PPP
      H(2,4)=L*M*(PPS-PPP)
      H(2,3)=L*N*(PPS-PPP)
         H(3,4)=M*N*(PPS-PPP)
!     == S-D BLOCK=================================================
      H(1,9)=SQ3*L*M*SDS
         H(1,8)=SQ3*M*N*SDS
         H(1,6)=SQ3*N*L*SDS
      H(1,5)=0.5D0*SQ3*(L2-M2)*SDS
      H(1,7)=(N2-0.5D0*(L2+M2))*SDS
!     == P-D BLOCK=================================================
      H(2,9)=SQ3*L2*M*PDS+M*(1.D0-2.D0*L2)*PDP        
        H(4,8)=SQ3*M2*N*PDS+N*(1.D0-2.D0*M2)*PDP
        H(3,6)=SQ3*N2*L*PDS+L*(1.D0-2.D0*N2)*PDP
      H(2,8)=SQ3*L*M*N*PDS-2.D0*L*M*N*PDP
        H(4,6)=SQ3*L*M*N*PDS-2.D0*L*M*N*PDP
        H(3,9)=SQ3*L*M*N*PDS-2.D0*L*M*N*PDP
      H(2,6)=SQ3*L2*N*PDS+N*(1.D0-2.D0*L2)*PDP   !X,XZ
        H(4,9)=SQ3*M2*L*PDS+L*(1.D0-2.D0*M2)*PDP     !Y,YX
        H(3,8)=SQ3*N2*M*PDS+M*(1.D0-2.D0*N2)*PDP     !Z,ZY
      H(2,5)=0.5D0*SQ3*L*(L2-M2)*PDS+L*(1.D0-L2+M2)*PDP
      H(4,5)=0.5D0*SQ3*M*(L2-M2)*PDS-M*(1.D0+L2-M2)*PDP
      H(3,5)=0.5D0*SQ3*N*(L2-M2)*PDS-N*(L2-M2)*PDP
      H(2,7)=L*(N2-0.5D0*(L2+M2))*PDS-SQ3*L*N2*PDP
      H(4,7)=M*(N2-0.5D0*(L2+M2))*PDS-SQ3*M*N2*PDP
      H(3,7)=N*(N2-0.5D0*(L2+M2))*PDS+SQ3*N*(L2+M2)*PDP
!     == D-D BLOCK=================================================
!     **      YLM(5)=SQRT(15/(16*PI))    * (  X**2-Y**2  ) /R**2      **
!     **      YLM(6)=SQRT(60/(16*PI))    * (     X*Z     ) /R**2      **
!     **      YLM(7)=SQRT( 5/(16*PI))    * ( 3*Z**2-R**2 ) /R**2      **
!     **      YLM(8)=SQRT(60/(16*PI))    * (      Y*Z    ) /R**2      **
!     **      YLM(9)=SQRT(60/(16*PI))    * (      X*Y    ) /R**2      **
      H(9,9)=3.D0*L2*M2*DDS+(L2+M2-4.D0*L2*M2)*DDP+(N2+L2*M2)*DDD   !XY,XY
      H(8,8)=3.D0*M2*N2*DDS+(M2+N2-4.D0*M2*N2)*DDP+(L2+M2*N2)*DDD   !YZ,YZ CYC.PERM
      H(6,6)=3.D0*N2*L2*DDS+(N2+L2-4.D0*N2*L2)*DDP+(M2+N2*L2)*DDD   !ZX,ZX CYC.PERM

      H(9,8)=3.D0*L*M2*N*DDS+L*N*(1-4.D0*M2)*DDP+L*N*(M2-1.D0)*DDD  !XY,YZ
      H(8,6)=3.D0*M*N2*L*DDS+M*L*(1-4.D0*N2)*DDP+M*L*(N2-1.D0)*DDD  !YZ,ZX CYC.PERM
!     H(6,9)=3.D0*N*L2*M*DDS+N*M*(1-4.D0*L2)*DDP+N*M*(L2-1.D0)*DDD  

      H(9,6)=3.D0*L2*M*N*DDS+M*N*(1-4.D0*L2)*DDP+M*N*(L2-1.D0)*DDD  !XY,XZ
!     H(8,9)=3.D0*M2*N*L*DDS+N*L*(1-4.D0*M2)*DDP+N*L*(M2-1.D0)*DDD  !
!     H(6,8)=3.D0*N2*L*M*DDS+L*M*(1-4.D0*N2)*DDP+L*M*(N2-1.D0)*DDD  !
      H(9,5)=1.5D0*L*M*(L2-M2)*DDS+2.D0*L*M*(M2-L2)*DDP &           !XY,X2-Y2 
     &                            +0.5D0*L*M*(L2-M2)*DDD            
      H(8,5)=1.5D0*M*N*(L2-M2)*DDS-M*N*(1.D0+2.D0*(L2-M2))*DDP &    !YZ,X2-Y2
     &      +M*N*(1.D0+0.5D0*(L2-M2))*DDD
      H(6,5)=1.5D0*N*L*(L2-M2)*DDS+N*L*(1.D0-2.D0*(L2-M2))*DDP &    !XZ,X2-Y2
     &      -N*L*(1.D0-0.5D0*(L2-M2))*DDD
      H(9,7)=SQ3*L*M*(N2-0.5D0*(L2+M2))*DDS-2.D0*SQ3*L*M*N2*DDP &   !XY,3Z2-R2
     &      +0.5D0*SQ3*L*M*(1.D0+N2)*DDD
      H(8,7)=SQ3*M*N*(N2-0.5D0*(L2+M2))*DDS+SQ3*M*N*(L2+M2-N2)*DDP & !YZ,3Z2-R2
     &      -0.5D0*SQ3*M*N*(L2+M2)*DDD
      H(6,7)=SQ3*L*N*(N2-0.5D0*(L2+M2))*DDS+SQ3*L*N*(L2+M2-N2)*DDP & !XZ,3Z2-R2
     &      -0.5D0*SQ3*L*N*(L2+M2)*DDD
      H(5,5)=0.75D0*(L2-M2)**2*DDS+(L2+M2-(L2-M2)**2)*DDP &          !X2-Y2,X2-Y2
     &      +(N2+0.25D0*(L2-M2)**2)*DDD
      H(5,7)=0.5D0*SQ3*(L2-M2)*(N2-0.5D0*(L2+M2))*DDS &              !X2-Y2,3Z2-R2
     &      +SQ3*N2*(M2-L2)*DDP+0.25D0*SQ3*(1.D0+N2)*(L2-M2)*DDD
      H(7,7)=(N2-0.5D0*(L2+M2))**2*DDS+3.D0*N2*(L2+M2)*DDP &         !3Z2-R2,3Z2-R2
     &      +0.75D0*(L2+M2)**2*DDD
!
!     =================================================================
!     == MAKE H HERMITEAN =============================================
!     =================================================================
      DO I=1,9
        P1=(-1.D0)**INT(SQRT(REAL(I)+1.D-6-1.D0))
        DO J=I+1,9
          P=P1*(-1.D0)**INT(SQRT(REAL(J)+1.D-6-1.D0))
          H(J,I)=P*H(I,J)+H(J,I)
          H(I,J)=P*H(J,I)
        ENDDO
      ENDDO
      RETURN
      END
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!****                                                                       ****
!****   TEST ROUTINES                                                       ****
!****                                                                       ****
!****                                                                       ****
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!
!     ..........................................................................
      SUBROUTINE TEST_MADELUNG()
!     **                                                                      **
!     **  CALCULATE VARIOUS MADELUNG CONSTANTS AND COMPARE WITH               **
!     **  "CONDENSED MATTER PHYSICS" BY M.P.MARDER                            **
!     **                                                                      **
      INTEGER(4),PARAMETER :: NBAS=2
      REAL(8)   ,PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8)              :: RBAS(3,3)
      REAL(8)              :: BAS(3,NBAS)
      REAL(8)              :: Q(NBAS)
      REAL(8)              :: EMAD
      REAL(8)              :: VMAD(NBAS)
      REAL(8)              :: FMAD(3,NBAS)
      REAL(8)              :: D,RS,DET
!     **************************************************************************
      Q(1)=1.D0
      Q(2)=-1.D0
!
!     == SODIUM CHLORIDE: 1.74757 ==============================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      BAS(:,1) =(/0.0D0,0.0D0,0.0D0/)
      BAS(:,2) =(/0.5D0,0.0D0,0.0D0/)
      D=SQRT(SUM((BAS(:,2)-BAS(:,1))**2))
      CALL MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      WRITE(*,FMT='("NACL STRUCTURE",T30,"E=",F10.5,"  DEV=",F10.5)') &
     &             -EMAD*D,-EMAD*D-1.74757
!
!     == CSCL STRUCTURE: 1.76268 ===============================================
      RBAS(1,:)=(/1.0D0,0.0D0,0.0D0/)
      RBAS(2,:)=(/0.0D0,1.0D0,0.0D0/)
      RBAS(3,:)=(/0.0D0,0.0D0,1.0D0/)
      BAS(:,1) =(/0.0D0,0.0D0,0.0D0/)
      BAS(:,2) =(/0.5D0,0.5D0,0.5D0/)
      D=SQRT(SUM((BAS(:,2)-BAS(:,1))**2))
      CALL MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      WRITE(*,FMT='("CSCL STRUCTURE",T30,"E=",F10.5,"  DEV=",F10.5)') &
     &            -EMAD*D,-EMAD*D-1.76268
!
!     == ZNS STRUCTURE: 1.63806 ================================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      BAS(:,1) =(/0.0D0,0.0D0,0.0D0/)
      BAS(:,2) =(/0.25D0,0.25D0,0.25D0/)
      D=SQRT(SUM((BAS(:,2)-BAS(:,1))**2))
      CALL MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      WRITE(*,FMT='("ZNS STRUCTURE",T30,"E=",F10.5,"  DEV=",F10.5)') &
     &            -EMAD*D,-EMAD*D-1.63806
!
!     ==========================================================================
!     == MADELUNG CONSTANTS FOR METALS                                        ==
!     == ATTENTION THE DEFINITION IS DIFFERENT!!!!                            ==
!     ==========================================================================
!
!     == FCC: 1.79186 ==========================================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      BAS(:,1) =(/0.0D0,0.0D0,0.0D0/)
      DET=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)) 
      RS=(3.D0*DET/(4.D0*PI))**(1.D0/3.D0)
      CALL MADELUNG(1,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      WRITE(*,FMT='("FCC STRUCTURE",T30,"E=",F10.5,"  DEV=",F10.5)') &
     &             -2.D0*EMAD*RS,-2.D0*EMAD*RS-1.79186
!
!     == SIC: 1.76012 ==========================================================
      RBAS(1,:)=(/1.0D0,0.0D0,0.0D0/)
      RBAS(2,:)=(/0.0D0,1.0D0,0.0D0/)
      RBAS(3,:)=(/0.0D0,0.0D0,1.0D0/)
      BAS(:,1) =(/0.0D0,0.0D0,0.0D0/)
      DET=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)) 
      RS=(3.D0*DET/(4.D0*PI))**(1.D0/3.D0)
      CALL MADELUNG(1,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      WRITE(*,FMT='("SIC STRUCTURE",T30,"E=",F10.5,"  DEV=",F10.5)') &
     &             -2.D0*EMAD*RS,-2.D0*EMAD*RS-1.76012
!
!     == DIAMOND: 1.67085 ======================================================
      RBAS(1,:)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(2,:)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(3,:)=(/0.5D0,0.5D0,0.0D0/)
      BAS(:,1) =(/0.0D0,0.0D0,0.0D0/)
      BAS(:,2) =(/0.25D0,0.25D0,0.25D0/)
      DET=RBAS(1,1)*(RBAS(2,2)*RBAS(3,3)-RBAS(3,2)*RBAS(2,3)) &
     &   +RBAS(2,1)*(RBAS(3,2)*RBAS(1,3)-RBAS(1,2)*RBAS(3,3)) &
     &   +RBAS(3,1)*(RBAS(1,2)*RBAS(2,3)-RBAS(2,2)*RBAS(1,3)) 
      DET=0.5D0*DET
      Q(:)=1.D0
      RS=(3.D0*DET/(4.D0*PI))**(1.D0/3.D0)
      CALL MADELUNG(NBAS,RBAS,BAS,Q,EMAD,VMAD,FMAD)
      WRITE(*,FMT='("DIAMOND STRUCTURE",T30,"E=",F10.5,"  DEV=",F10.5)') &
     &             -EMAD*RS,-EMAD*RS-1.67085
      RETURN
      END
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BINOMIALCOEFFICIENT(I,J,RES)
!     **************************************************************************
!     **  CALCULATES THE BINOMIAL COEFFICIENT (I OVER J)=I!/J!/(I-J)!         **
!     **                                                                      **
!     **  USES A LOOKUP TABLE                                                 **
!     **  USES RECURSIVE FORMULA (N OVER K)=(N-1 OVER K-1) (N-1 OVER K)       **
!     **  WITH INITIAL VALUES   (N OVER 0)=1                                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)       :: I
      INTEGER(4),INTENT(IN)       :: J
      REAL(8)   ,INTENT(OUT)      :: RES
      INTEGER(4),SAVE             :: NX=-1
      INTEGER(4),ALLOCATABLE,SAVE :: B(:,:)
      INTEGER(4)                  :: N,K     
!     **************************************************************************
!
!     ==========================================================================
!     == RETURN RESULT, IF ALREADY AVAILABLE                                  ==
!     ==========================================================================
      IF(I.LT.NX) THEN
        RES=REAL(B(I,J),KIND=8)
        RETURN
      END IF
!
!     ==========================================================================
!     == REALLOCATE PERMANENT ARRAY                                           ==
!     ==========================================================================
      IF(ALLOCATED(B)) DEALLOCATE(B)
      NX=MAX(I,NX+20)
      ALLOCATE(B(0:NX,0:NX))
!
!     ==========================================================================
!     == CALCULATE TABLE OF BINOMIAL COEFFICIENTS                             ==
!     ==========================================================================
      B(:,:)=0
      DO N=0,NX
        B(N,0)=1     !INITIAL VALUE (N OVER 0)=1
        DO K=1,N-1
          B(N,K)=B(N-1,K-1)+B(N-1,K)   ! RECURSIVE FORMULA
        ENDDO
        B(N,N)=1     ! INITIAL VALUE (N OVER N)=1
      ENDDO
!
!     ==========================================================================
!     == RETURN RESULT                                                        ==
!     ==========================================================================
      RES=REAL(B(I,J),KIND=8)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEBITSTRINGI4(NFIL,NAME,VAL)
!     **************************************************************************
!     ** WRITE THE BIT STRING OF A INTEGER(4) VARIABLE VAL                    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: VAL
      CHARACTER(64)           :: BITSTRING
      INTEGER                 :: I
!     **************************************************************************
      BITSTRING=" "
      DO I=1,BIT_SIZE(VAL)
        IF(BTEST(VAL,I-1)) THEN
          BITSTRING(I:I)='1'
        ELSE 
         BITSTRING(I:I)='0'
        END IF
      ENDDO
      WRITE(NFIL,FMT='(" BITSTRING=",A," VALUE=",I10," VARIABLE NAME=",A)') &
                     BITSTRING(1:BIT_SIZE(VAL)),VAL,TRIM(NAME)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEBITSTRINGL4(NFIL,NAME,VAL)
!     **************************************************************************
!     ** WRITE THE BIT STRING OF A LOGICAL(4) VARIABLE VAL                    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      LOGICAL(4)  ,INTENT(IN) :: VAL
      INTEGER                 :: IVAL
      INTEGER                 :: LEN
      CHARACTER(64)           :: BITSTRING
      INTEGER                 :: I
!     **************************************************************************
      LEN=STORAGE_SIZE(VAL)  ! NUMBER OF BITS TO HOLD VALUE OF VAL IN MEMORY
      LEN=MIN(LEN,BIT_SIZE(IVAL))
      IVAL=TRANSFER(VAL,IVAL)
      BITSTRING=" "
      DO I=1,LEN
        IF(BTEST(IVAL,I-1)) THEN
          BITSTRING(I:I)='1'
        ELSE 
         BITSTRING(I:I)='0'
        END IF
      ENDDO
      WRITE(NFIL,FMT='(" BITSTRING=",A," VALUE=",L10," VARIABLE NAME=",A)') &
                     BITSTRING(1:LEN),VAL,TRIM(NAME)
      RETURN
      END
