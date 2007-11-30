!
!     ..................................................................
      SUBROUTINE SHADOW
      IMPLICIT INTEGER ($)                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TCHK,TCHK1,TINI
      LOGICAL TCHK_
      CHARACTER*(*) IDENT_
      DIMENSION RMASS_(NAT_),FORCE_(3,NAT_),R0_(3,NAT_)
!     == POINTER ARRAYS
      POINTER ($VAL,VAL)
      POINTER ($D2E,D2E),($W2,W2),($U,U)
      DIMENSION D2E(3*NAT,3*NAT),W2(3*NAT),U(3*NAT,3*NAT)
      POINTER ($WORK,WORK),($WORK2,WORK2)
      DIMENSION WORK(*),WORK2(3*NAT)
!     == FIXED DIMENSION ARRAYS
      CHARACTER*32 IDENT
      DIMENSION SV(6)
      DATA TINI/.FALSE./
      CALL ERROR$STOP('SHADOW')
!
!     ******************************************************************
!     **   SET                                                        **
!     ******************************************************************
      ENTRY SHADOW$SET(IDENT_,NBYTE_,VAL_)
        CALL LINKEDLIST$SET($LIST,IDENT_,NBYTE_,VAL_)
        RETURN
!
!     ******************************************************************
!     **   GET                                                        **
!     ******************************************************************
      ENTRY SHADOW$GET(IDENT_,NBYTE_,VAL_)
        CALL LINKEDLIST$EXIST($LIST,IDENT_,TCHK)
        IF(.NOT.TCHK) THEN
          IDENT=IDENT_
          IF(IDENT(1:8).EQ.'OPTIMIZE') THEN
            CALL LINKEDLIST$SET($LIST,'OPTIMIZE',4,.FALSE.)
          ELSE IF(IDENT(1:12).EQ.'PRECONDITION') THEN
            CALL LINKEDLIST$SET($LIST,'PRECONDITION',4,.FALSE.)
          END IF
        END IF
        CALL LINKEDLIST$GET($LIST,IDENT_,NBYTE_,VAL_)
        RETURN
!
!     ******************************************************************
!     **   INITIALIZE                                                 **
!     ******************************************************************
      ENTRY SHADOW$INITIALIZE
        IF(TINI) THEN
          CALL ERROR$MSG('ATTEMPT TO INITIALIZE SHADOW TWICE')
          CALL ERROR$STOP('SHADOW$INITIALIZE')
        END IF
!
!       == CHECK WHETHER SHADOW WILL BE USED ===========================
        TCHK1=.FALSE.
        CALL LINKEDLIST$EXIST($LIST,'OPTIMIZE',TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET($LIST,'OPTIMIZE',4,TCHK)
          TCHK1=TCHK
        ELSE
          CALL LINKEDLIST$SET($LIST,'OPTIMIZE',4,.FALSE.)
        END IF
        CALL LINKEDLIST$EXIST($LIST,'PRECONDITION',TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET($LIST,'PRECONDITION',4,TCHK)
          TCHK1=TCHK
        ELSE
          CALL LINKEDLIST$SET($LIST,'PRECONDITION',4,.FALSE.)
        END IF
        IF(.NOT.TCHK1) RETURN
        TINI=.TRUE.
!
!       == SET DATA FOR CLASSICAL OBJECT ===============================
        CALL CLASSICAL$SELECT('SHADOW')
        CALL ATOMLIST$NATOM(NAT)
        CALL CLASSICAL$SET('NAT',4,NAT)
        CALL STACK$ALLOCATE(8,$VAL,3*NAT)
        CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,VAL)
        CALL CLASSICAL$SET('R(0)',8*3*NAT,VAL)
        CALL ATOMLIST$GETR8A('R(-)',0,3*NAT,VAL)
        CALL CLASSICAL$SET('R(-)',8*3*NAT,VAL)
        CALL ATOMLIST$GETR8A('MASS',0,NAT,VAL)
        CALL CLASSICAL$SET('MASS',8*NAT,RMASS)
        CALL ATOMLIST$GETR8A('Q',0,NAT,VAL)
        CALL CLASSICAL$SET('QEL',8*NAT,QEL)
        CALL ATOMLIST$GETCHA('UFFTYPE',0,NAT,VAL)
        CALL CLASSICAL$SET('TYPE',5*NAT,VAL)
        CALL STACK$FREE($VAL)
!
        CALL LINKEDLIST$GET($LIST,'NBOND',4,NBOND)
        CALL CLASSICAL$SET('NBOND',4,NBOND)
        CALL STACK$ALLOCATE(MAX(4*3,8)*NBOND,$VAL,1)
        CALL LINKEDLIST$GET($LIST,'BONDORDER',8*NBOND,VAL)
        CALL CLASSICAL$SET('BONDORDER',8*NBOND,VAL)
        CALL LINKEDLIST$GET($LIST,'BOND',4*3*NBOND,VAL)
        CALL CLASSICAL$SET('INDEX2',4*3*NBOND,VAL)
        CALL LINKEDLIST$GET($LIST,'TLONGRANGE',4,VAL)
        CALL CLASSICAL$SET('TLONGRANGE',4,VAL)
        CALL STACK$FREE($VAL)
!
        CALL LINKEDLIST$DELETE($LIST,'BOND')
        CALL LINKEDLIST$DELETE($LIST,'BONDORDER')
!
        CALL CLASSICAL$SELECT('SHADOW')
        CALL CLASSICAL$INITIALIZE
        CALL CLASSICAL$SELECT('SHADOW')
!
!       CALL FILEHANDLER$UNIT('PROT',NFILO)
!       CALL CLASSICAL$REPORT(NFILO)
        RETURN
!
!     ******************************************************************
!     **   PRECONVERGE                                                **
!     ******************************************************************
      ENTRY SHADOW$OPTIMIZE(TCHK_)
        IF(.NOT.TINI) RETURN
        CALL LINKEDLIST$GET($LIST,'OPTIMIZE',4,TCHK)
        IF(.NOT.TCHK) RETURN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,FMT='("PREOPTIMIZE WITH SHADOW")')
        CALL CLASSICAL$SELECT('SHADOW')
        CALL LINKEDLIST$GET($LIST,'TOL',8,TOL)
        CALL LINKEDLIST$GET($LIST,'NSTEP',4,NITER)
        WRITE(NFILO,FMT='("PREOPTIMIZE WITH SHADOW")')
        CALL CLASSICAL$MINIMIZE(TOL,NSTEP,TCHK)
        IF(TCHK) THEN
          WRITE(NFILO,FMT='("PREOPTIMIZE CONVERGED")')
        ELSE
          WRITE(NFILO,FMT='("PREOPTIMIZE NOT CONVERGED")')
        END IF
!       CALL CONSTRAINTS$APPLY(NAT,R0,RP,REDRMASS,FORCE,DELT)
!       CALL SYMMETRIZE$POSITION(NAT,RP)
!
!       == COPY NEW ATOMIC POSITIONS INTO ATOMLIST =====================
        CALL STACK$ALLOCATE(8,$WORK,3*NAT)
        CALL CLASSICAL$GET('R(0)',8*3*NAT,WORK)
        CALL ATOMLIST$SETR8A('R(0)',0,3*NAT,WORK)
        CALL STACK$FREE($WORK)
        TCHK_=TCHK
        RETURN
!
!     ******************************************************************
!     **   PRECONDITION                                               **
!     ******************************************************************
      ENTRY SHADOW$PRECONDITION(NAT_,RMASS_,R0_,FORCE_)
        IF(.NOT.TINI) RETURN
        CALL LINKEDLIST$GET($LIST,'PRECONDITION',4,TCHK)
        IF(.NOT.TCHK) RETURN
        CALL CLASSICAL$SELECT('SHADOW')
        CALL CLASSICAL$GET('NAT',4,NAT)
        IF(NAT.NE.NAT_) THEN
          CALL ERROR$MSG('NAT NOT CONSISTENT')
          CALL ERROR$STOP('SHADOW$PRECONDITION')
        END IF
        CALL LINKEDLIST$GET($LIST,'OMEGA1',8,OMEGA1)
        CALL LINKEDLIST$GET($LIST,'OMEGA2',8,OMEGA2)
!       WRITE(*,FMT='("OLD FORCE ",12E10.3)') &
!    &       ((FORCE_(I,IAT),I=1,3),IAT=1,NAT)
!
!       ================================================================
!       == EVALUATE MATRIX OF SECOND DERIVATIVES                      ==
!       ================================================================
        CALL CLASSICAL$SET('R(0)',8*3*NAT_,R0_)
        CALL STACK$ALLOCATE(8,VAL,NAT)
        CALL ATOMLIST$GETR8A('QEL',0,NAT,VAL)
        CALL CLASSICAL$SET('QEL',8*NAT_,VAL)
        CALL STACK$FREE($VAL)
        CALL STACK$ALLOCATE(8,$D2E,3*NAT*3*NAT)
        CALL CLASSICAL$D2EDX2(NAT,D2E)
!       PRINT*,'D2E'
!       DO I=1,3*NAT
!         WRITE(*,FMT='(12E10.3)')(D2E(J,I),J=1,3*NAT)
!       ENDDO
!
!       ================================================================
!       == EXTRACT ROTATIONS AND TRANSLATIONS                         ==
!       ================================================================
        CALL CLASSICAL_D2ECORRECT(NAT,R0_,D2E)
!       PRINT*,'D2E'
!       DO I=1,3*NAT
!         WRITE(*,FMT='(12E10.3)')(D2E(J,I),J=1,3*NAT)
!       ENDDO
!
!       ================================================================
!       == EVALUATE DYNAMICAL MATRIX                                  ==
!       ================================================================
        CALL STACK$ALLOCATE(8,$WORK,3*NAT)
        CALL STACK$ALLOCATE(8,$WORK2,3*NAT)
        II=0
        DO IAT=1,NAT
          SVAR=1.D0/SQRT(RMASS_(IAT))
          DO I=1,3
            II=II+1
            WORK(II)=SVAR
            WORK2(II)=FORCE_(I,IAT)*SVAR
          ENDDO
        ENDDO
        DO I1=1,3*NAT
          SVAR1=WORK(I1)
          DO I2=1,3*NAT
            SVAR=SVAR1*WORK(I2)
            D2E(I1,I2)= D2E(I1,I2)*SVAR 
          ENDDO
        ENDDO
        CALL STACK$FREE($WORK)
!
!       PRINT*,'DYNAMICAL MATRIX'
!       DO I=1,3*NAT
!         WRITE(*,FMT='(12E10.3)')(D2E(J,I),J=1,3*NAT)
!       ENDDO
!       PRINT*,'FORCE/SQRT(M)'
!       WRITE(*,FMT='(12E10.3)')(WORK2(J),J=1,3*NAT)
!    
!       ================================================================
!       == DIAGONALIZE DYNAMICAL MATRIX                               ==
!       ================================================================
        CALL STACK$ALLOCATE(8,$W2,3*NAT)
        CALL STACK$ALLOCATE(8,$U,9*NAT*NAT)
        CALL LIB$DIAGR8(3*NAT,D2E,W2,U)
!CALL DIAG(3*NAT,3*NAT,D2E,W2,U)
!       WRITE(*,FMT='("EIG",12E10.3)')(W2(J),J=1,3*NAT)
        W2IN=2.D0*OMEGA1**2
        W2OUT=OMEGA2**2
        DO I=1,3*NAT
          SVAR2=ABS(W2(I)/W2IN)
          IF(SVAR2.LT.1.D-10) THEN
            W2(I)=0.D0
          ELSE IF(SVAR2.LT.1.D-2) THEN
            W2(I)=W2OUT/W2IN
          ELSE IF(SVAR2.GT.1.D+2) THEN
            W2(I)=2.D0*W2OUT/ABS(W2(I))
          ELSE 
            SVAR=EXP(W2(I)/W2IN)
            W2(I)=2.D0*W2OUT/W2(I)*(SVAR-1.D0)/(SVAR+1.D0)
          END IF
        ENDDO
!       WRITE(*,FMT='("FAC",12E10.3)')(W2(J),J=1,3*NAT)
        CALL STACK$ALLOCATE(8,$WORK,3*NAT)
        CALL DGEMV('T',3*NAT,3*NAT,1.D0,U,3*NAT,WORK2,1,0.D0,WORK,1)
!       WRITE(*,FMT='("<U|F>",12E10.3)')(WORK(J),J=1,3*NAT)
        DO I=1,3*NAT
          WORK(I)=WORK(I)*W2(I)
        ENDDO
!       WRITE(*,FMT='("C<U|F>",12E10.3)')(WORK(J),J=1,3*NAT)
!       WRITE(*,FMT='("U",12E10.3)')(U(J,3*NAT),J=1,3*NAT)

        CALL DGEMV('N',3*NAT,3*NAT,1.D0,U,3*NAT,WORK,1,0.D0,FORCE_,1)

!       WRITE(*,FMT='("NEW FORCE ",12E10.3)') &
!    &       ((FORCE_(I,IAT),I=1,3),IAT=1,NAT)
!       DO I=1,3*NAT
!         DO J=1,3*NAT
!           SUM=0.D0
!           DO K=1,3*NAT
!             SUM=SUM+U(I,K)*U(J,K)*W2(K)
!           ENDDO
!         ENDDO
!       ENDDO
        CALL STACK$FREE($WORK)
        CALL STACK$FREE($WORK2)
        CALL STACK$FREE($U)
        CALL STACK$FREE($W2)
        CALL STACK$FREE($D2E)
!    
!       ================================================================
!       == SINGULAR..I DYNAMICAL MATRIX                               ==
!       ================================================================
!       CALL STACK$ALLOCATE(8,$W2,3*NAT)
!       CALL STACK$ALLOCATE(8,$WORK,9*NAT)
!       CALL DGESVF(12,D2E,3*NAT,WORK2,3*NAT,1,W2,3*NAT &
!    &             ,3*NAT,WORK,9*NAT)
!       WRITE(*,FMT='("SV ",12E10.3)')(W2(J),J=1,3*NAT)
!       W2IN=2.D0*WIN_**2
!       W2OUT=WOUT_**2
!       DO I=1,3*NAT
!         SVAR1=W2(I)
!         SVAR2=ABS(W2(I)/W2IN)
!         IF(SVAR2.LT.1.D-2) THEN
!           W2(I)=W2OUT/W2IN
!         ELSE IF(SVAR2.GT.1.D+2) THEN
!           W2(I)=2.D0*W2OUT/W2(I)
!         ELSE 
!           SVAR=EXP(W2(I)/W2IN)
!           W2(I)=2.D0*W2OUT/W2(I)*(SVAR-1.D0)/(SVAR+1.D0)
!         END IF
!         WRITE(*,FMT='("W2BEF",E15.5,"W2AFT",2E15.5)')SVAR1,W2(I),SVAR2
!         W2(I)=1.D0/W2(I)
!       ENDDO
!       CALL DGESVS(D2E,3*NAT,WORK2,3*NAT,1,W2,FORCE_ &
!    &             ,3*NAT,3*NAT,3*NAT,1.D-7)
!       CALL STACK$FREE($W2)
!       CALL STACK$FREE($WORK)
!       CALL STACK$FREE($WORK2)
!       CALL STACK$FREE($D2E)
!
!       ================================================================
!       == |F>= M**{-1/2)|U_J> MAX(1,W2MAX/W2_J) <U_J|M**(+1/2)|F>    ==
!       ================================================================
        DO IAT=1,NAT
          SVAR=SQRT(RMASS_(IAT))
          DO I=1,3
            FORCE_(I,IAT)=FORCE_(I,IAT)*SVAR
          ENDDO
        ENDDO
!       WRITE(*,FMT='("NEW FORCE ",12E10.3)') &
!    &       ((FORCE_(I,IAT),I=1,3),IAT=1,NAT)
!       STOP
        RETURN
      END
!
!     ..................................................................
      SUBROUTINE CLASSICAL_D2ECORRECT(NAT,R0,D2E)
      IMPLICIT INTEGER ($)                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D2E(3*NAT,3*NAT),R0(3,NAT)        
!     == POINTER ARRAYS ================================================
      POINTER ($VEC,VEC),($D2EVEC,D2EVEC)
      DIMENSION VEC(3*NAT,6),D2EVEC(3*NAT,6)
      POINTER ($VECVEC,VECVEC),($VECD2EVEC,VECD2EVEC)
      DIMENSION VECVEC(6,6),VECD2EVEC(6,6)
      POINTER ($WORK1,WORK1),($WORK2,WORK2),($WORK3,WORK3)
!     == FIXED DIMENSION ARRAYS ========================================
      DIMENSION SV(6)
!
      CALL STACK$ALLOCATE(8,$VEC,3*NAT*6) 
      CALL STACK$ALLOCATE(8,$D2EVEC,3*NAT*6)
      DO IVEC=1,6
        DO I=1,3*NAT
          VEC(I,IVEC)=0.D0
          D2EVEC(I,IVEC)=0.D0
        ENDDO
      ENDDO
      II=0
      DO IAT=1,NAT
        VEC(II+1,1)=1.D0
        VEC(II+2,2)=1.D0
        VEC(II+3,3)=1.D0
        VEC(II+2,4)=+R0(3,IAT)
        VEC(II+3,4)=-R0(2,IAT)
        VEC(II+3,5)=+R0(1,IAT)
        VEC(II+1,5)=-R0(3,IAT)
        VEC(II+1,6)=+R0(2,IAT)
        VEC(II+2,6)=-R0(1,IAT)
        II=II+3
      ENDDO
!
      DO J=1,6
        SUM=0.D0
        DO I=1,3*NAT
          SUM=SUM+VEC(I,J)**2
        ENDDO
        IF(SUM.NE.0.D0) SUM=1.D0/SQRT(SUM)
        DO I=1,3*NAT
          VEC(I,J)=VEC(I,J)*SUM
        ENDDO
      ENDDO
!
      DO I=1,3*NAT
        DO J=1,3*NAT
          SVAR=D2E(I,J)
          D2EVEC(I,1)=D2EVEC(I,1)+SVAR*VEC(J,1)
          D2EVEC(I,2)=D2EVEC(I,2)+SVAR*VEC(J,2)
          D2EVEC(I,3)=D2EVEC(I,3)+SVAR*VEC(J,3)
          D2EVEC(I,4)=D2EVEC(I,4)+SVAR*VEC(J,4)
          D2EVEC(I,5)=D2EVEC(I,5)+SVAR*VEC(J,5)
          D2EVEC(I,6)=D2EVEC(I,6)+SVAR*VEC(J,6)
        ENDDO
      ENDDO
!     DO I=1,6
!       PRINT*,'VEC ',(VEC(J,I),J=1,3*NAT)
!       PRINT*,'D2EVEC',(D2EVEC(J,I),J=1,3*NAT)
!     ENDDO
!
!     == VECVEC=<VEC|VEC>; VECD2EVEC=<VEC|D2E|VEC>  ==================
      CALL STACK$ALLOCATE(8,$VECVEC,6*6)
      CALL STACK$ALLOCATE(8,$VECD2EVEC,6*6)
      DO I=1,6
        DO J=I,6
          SVAR1=0.D0
          SVAR2=0.D0
          DO K=1,3*NAT
            SVAR1=SVAR1+VEC(K,I)*VEC(K,J)
            SVAR2=SVAR2+VEC(K,I)*D2EVEC(K,J)
          ENDDO
          VECVEC(I,J)=SVAR1
          VECVEC(J,I)=SVAR1
          VECD2EVEC(I,J)=SVAR2
          VECD2EVEC(J,I)=SVAR2
        ENDDO
      ENDDO
!     DO I=1,6
!       WRITE(*,FMT='("VECD2EVEC",6F10.5)')(VECD2EVEC(I,J),J=1,6)
!     ENDDO
!     DO I=1,6
!       WRITE(*,FMT='("VECVEC",6F10.5)')(VECVEC(I,J),J=1,6)
!     ENDDO
!     == |VEC>=|VEC><VEC|VEC>**-1  ===================================
      CALL STACK$ALLOCATE(8,$WORK1,3*NAT*6)
!     __WORK=TRANSPOSE(VEC)___________________________________________
      CALL DGETMO(VEC,3*NAT,3*NAT,6,WORK1,6)
!     __SINGULAR VALUE DECOMPOSITION__________________________________
      NAUX=2*6+MAX(6,3*NAT)
      CALL STACK$ALLOCATE(8,$WORK2,NAUX)
      CALL DGESVF(12,VECVEC,6,WORK1,6,3*NAT,SV,6,6,WORK2,NAUX)
      CALL STACK$FREE($WORK2)
!     PRINT*,'SINGVAL',SV
!     __ INVERSION_____________________________________________________
      CALL STACK$ALLOCATE(8,$WORK2,3*NAT*6)
      CALL DGESVS(VECVEC,6,WORK1,6,3*NAT,SV,WORK2,6,6,6,1.D-1)
!     __ VEC=TRANSPOSE(WORK2)
      CALL DGETMO(WORK2,6,6,3*NAT,VEC,3*NAT)
      CALL STACK$FREE($WORK2)
      CALL STACK$FREE($WORK1)
!     == D2E=D2E-|VEC><D2EVEC|-|D2EVEC><VEC|+|VEC><VECD2EVEC><VEC|
      DO J=1,3*NAT
        DO K=1,6
          SVAR1=-D2EVEC(J,K)
          SVAR2=-VEC(J,K)
          DO L=1,6
            SVAR1=SVAR1+VECD2EVEC(K,L)*VEC(J,L)
          ENDDO
          DO I=1,3*NAT
              D2E(I,J)=D2E(I,J)+VEC(I,K)*SVAR1+D2EVEC(I,K)*SVAR2
          ENDDO          
        ENDDO
      ENDDO
      CALL STACK$FREE($VECVEC)
      CALL STACK$FREE($VECD2EVEC)
      CALL STACK$FREE($VEC)
      CALL STACK$FREE($D2EVEC)
      RETURN
      END
      

