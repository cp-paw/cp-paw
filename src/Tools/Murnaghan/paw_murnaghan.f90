     PROGRAM MAIN
!    ***************************************************************************
!    ** DETERMINE EQUILIBRIUM VOLUME, EQUILIBRIUM ENERGY, BULK MODULUS        **
!    ** FROM A MURNAGHAN FIT TO A SET OF (VOLUME,ENERGY) DATA RECEIVED FROM   **
!    ** STANDARD INPUT                                                        **
!    ***************************MATTHE UITTEVAAL AND PETER E. BLOECHL***********
     IMPLICIT NONE
     REAL(8)    :: V(100), E(100)     ! INPUT VOLUMES, ENERGIES
     INTEGER(4) :: HIGH,LOW,I, NP=0   ! INDICES, NO. DATA POINTS (INIT 0)
     INTEGER(4) :: IARRAY(1)          ! 1D ARRAY FOR MIN INDEX
     REAL(8)    :: PARMS(4),GRAD(4),X ! PARMS, GRAD AND MEAN SQR DEV OF MURN FIT
     REAL(8)    :: VI, EFIT           ! VOLUME IN, ENERGY OUT
!    ***************************************************************************
!    ==========================================================================
!    == START WITH PRINTING HEADER
!    ==========================================================================
     WRITE(*,FMT='(80("="))')
     WRITE(*,FMT='(80("="),T15," FIT MURNAGHAN EQUATION OF STATE ")')
     WRITE(*,FMT='(80("="))')
     WRITE(*,FMT='(A80)') &
    &     'MURNAGHAN EQUATION OF STATE IS BASED ON THE ASSUMPTION THAT THE' &
    &    ,'BULK MODULUS DEPENDS LINEARLY ON PRESSURE.' &
    &    ,'SOURCE: F.D. MURNAGHAN, PNAS 30, 244 (1934)'
!
!    ==========================================================================
!    == READ DATA UNTIL ERROR THEN RESET NP
!    ==========================================================================
     DO 
       NP=NP+1
       IF (NP.EQ.101) THEN
         WRITE(*,*)"TOO MUCH INPUT DATA, ONLY FIRST 100 USED"
         EXIT
       END IF
       READ(*,*,END=100,ERR=100)V(NP),E(NP)
     ENDDO
100  CONTINUE
     NP=NP-1
     IF (NP.LT.4) THEN
       WRITE(*,*)"TOO FEW INPUT DATA, AT LEAST 4 REQUIRED!"
       STOP
     END IF
!    == REPORT INPUT DATA =====================================================
     DO I=1,NP
       WRITE(*,FMT='(I5," V=",F10.5," E=",F10.5)')I,V(I),E(I)
     ENDDO
!
!    ==========================================================================
!    == PRINT DATA 
!    ==========================================================================
     IARRAY=MINLOC(V(1:NP))
     LOW=IARRAY(1) 
     IARRAY=MAXLOC(V(1:NP))
     HIGH=IARRAY(1) 
!     WRITE(*,FMT='(" VMIN=",F10.5," VMAX=",F10.5)')V(LOW),V(HIGH)
     IARRAY=MINLOC(E(1:NP)) 
!    I=MINLOC(V(1:NP))(1) ??
     PARMS(1)=E(IARRAY(1)) ! MINE
     PARMS(2)=V(IARRAY(1)) ! VOL MINE
!     WRITE(*,FMT='(" V(EMIN)=",F10.5," EMIN=",F10.5)')PARMS(2),PARMS(1)
     PARMS(3)=2*(E(HIGH)+E(LOW)-2*PARMS(1))/(V(HIGH)-V(LOW))  !INITIAL B0
     PARMS(4)=3.5D0  ! APPROX CONST
!     WRITE(*,FMT='(" BINI=",F10.5," BPRIMEINI=",F10.5)')PARMS(3),PARMS(4)
!
!    ==========================================================================
!    == DO FIT
!    ==========================================================================
     CALL CG(NP,V,E,PARMS,X)
!
!    ==========================================================================
!    == PRINT PARAMETERS
!    ==========================================================================
     WRITE(*,*)
     WRITE(*,FMT='(80("="),T10,"FIT PARAMETERS OF MURNGAHAN EQUATION OF STATE")')
     WRITE(*,FMT='("EQUILIBRIUM ENERGY E0",T40,F10.5)')PARMS(1)
     WRITE(*,FMT='("EQUILIBRIUM VOLUME V0",T40,F10.5)')PARMS(2)
     WRITE(*,FMT='("EQUILIBRIUM BULK MODULUS B0",T40,F10.5)')PARMS(3)
     WRITE(*,FMT='("PRESSURE DERIVATIVE OF BULK MODULUS BPRIME",T40,F10.5)') &
   &       PARMS(4)
!
!    ===========================================================================
!    == PRINT ENERGY VOLUME CURVE WITH INPUT FOR COMPARISON                   ==
!    ===========================================================================
     WRITE(*,*)
     WRITE(*,FMT='(80("="),T10,"COMPARE ORIGINAL AND INTERPOLATED DATA")')
     DO I=1,NP
       CALL MURNAGHAN(PARMS,V(I),EFIT,GRAD)
       WRITE(*,FMT='(I5," V=",F10.5," E(IN)=",F10.5," E(FIT)=",F10.5," E(FIT)-E(IN))=",F10.5)') &
    &        I,V(I),E(I),EFIT,EFIT-E(I)
     ENDDO
!
!    ===========================================================================
!    == WRITE EQUIDISTANT ENERGY VOLUME CURVE (INPUT RANGE +10%)              ==
!    ===========================================================================
     WRITE(*,FMT='(80("="),T10,"INTERPOLATED EQUATION OF STATE IS WRITTEN TO FILE MURN.DAT")')
     OPEN(UNIT=8,FILE='MURN.DAT')
     DO I=-10,110
       VI=V(LOW)+(V(HIGH)-V(LOW))/REAL(100)*REAL(I-1)
       CALL MURNAGHAN(PARMS,VI,EFIT,GRAD)
       WRITE(8,FMT='(2F10.5)')VI,EFIT
     ENDDO
     STOP
     END
!
!    ...........................................................................
     SUBROUTINE CG(NP,V,E,PARMS,X)
!    ***************************************************************************
!    ** CONJUGATE GRADIENT FIT OF THE MURNAGHAN CURVE (ITS PARAMETERS)        **
!    ** RETURNS THE MURNAGHAN PARAMETERS AND THE QUALITY OF THE FIT           **
!    **************************MATTHE UITTEVAAL*********************************
     IMPLICIT NONE
     INTEGER(4),INTENT(IN)   :: NP
     REAL(8)   ,INTENT(IN)   :: V(NP)        ! VOLUMES OF INPUT DATA SET
     REAL(8)   ,INTENT(IN)   :: E(NP)        ! ENERGIES OF THE INPUT DATA SET
     REAL(8)   ,INTENT(INOUT):: PARMS(4)     ! MURNAGHAN PARAMETERS
     REAL(8)   ,INTENT(OUT)  :: X            ! QUALITY OF THE FIT
     REAL(8)   ,PARAMETER    :: STEP=1.D-4
     REAL(8)   ,PARAMETER    :: FACT=.13
     REAL(8)   ,PARAMETER    :: TOL=7.D-5
     REAL(8)                 :: PARMS1(4),PARMS2(4)
     REAL(8)                 :: GRAD1(4),GRAD2(4)
     REAL(8)                 :: LAPL(4)
     LOGICAL(4)              :: TCONV
     INTEGER                 :: LOOP=0
     INTEGER   ,PARAMETER    :: LOOPMAX=1000
     INTEGER(4),PARAMETER    :: NFIL=6
!    ***************************************************************************
!     OPEN(NFIL,FILE='ITER.DAT')
     WRITE(NFIL,*)
     WRITE(NFIL,FMT='(80("="),T10," OPTIMIZING PARAMETERS ...")')
     CALL ONESTEP(NP,V,E,PARMS,X,GRAD1)
     WRITE(NFIL,FMT='(I4," E=",F10.5," V=",F10.5," B=",F10.5," DB/DV=",F10.5," PENALTY=",E15.7)') &
    &     LOOP,PARMS,X
     PARMS1(:)=PARMS(:)*(1.D0+STEP) !PERCENTAGE (DIMENSIONS,ENERGY 0??)
     CALL ONESTEP(NP,V,E,PARMS1,X,GRAD2)
     LAPL(:)=(GRAD2(:)-GRAD1(:))/PARMS(:)/STEP
!     WRITE(NFIL,FMT='("          GRADIENT AND LAPLACIAN OF PENALTY FUNCTION ")')
!     WRITE(NFIL,FMT='(4E10.2,"::",4E10.2)')GRAD1, LAPL
!
     PARMS1=PARMS
     DO LOOP=1,LOOPMAX
       PARMS2=PARMS1-GRAD1/LAPL*FACT  !FACT TOWARDS MINIMUM
       TCONV=(SQRT(SUM((GRAD1/LAPL/PARMS1)**2)).LT.TOL) 
       IF(TCONV) EXIT
       CALL ONESTEP(NP,V,E,PARMS2,X,GRAD2)
       WRITE(NFIL,FMT='(I4," E=",F10.5," V=",F10.5," B=",F10.5," DB/DV=",F10.5," PENALTY=",E15.7)') &
    &                 LOOP,PARMS2,X
       LAPL(:)=(GRAD2(:)-GRAD1(:))/(PARMS2-PARMS1)
!       WRITE(NFIL,FMT='("        GRADIENT AND LAPLACIAN OF PENALTY FUNCTION ___")')
!       WRITE(NFIL,FMT='(4E10.2,"::"4E10.2)')GRAD2, LAPL
       GRAD1=GRAD2
       PARMS1=PARMS2
     ENDDO
     IF (.NOT.TCONV) THEN
       WRITE(NFIL,FMT='("-----CONJUGATE GRADIENT NOT CONVERGED---------")')
     END IF
     CLOSE(NFIL)
     PARMS=PARMS2
     RETURN
     END
!
!    ...........................................................................
     SUBROUTINE ONESTEP(NP,V,E,PARMS,X,GRAD)
!    ***************************************************************************
!    ** PENALTY FUNCTION X (MEAN SQUARE DISPLACEMENT) OF THE FIT              **
!    ** AND GRADIENT WITH RESPECT TO FIT PARAMETERS                           **
!    ***************************************************************************
     IMPLICIT NONE
     INTEGER(4),INTENT(IN) :: NP           !#(DATA)
     REAL(8)   ,INTENT(IN) :: V(NP),E(NP)  !VOLUMES, ENERGIES
     REAL(8)   ,INTENT(IN) :: PARMS(4)     !CURRENT FIT PARAMETERS  
     REAL(8)   ,INTENT(OUT):: X,GRAD(4)    !VALUE AND GRADIENT OF PENALTY FUNC.
     REAL(8)               :: DEDP(4)      !GRADIENT OF MURNAGHAN CURVE
     REAL(8)               :: EFIT         !FITTING ENERGY
     INTEGER(4)            :: IP           !COUNTER
!    ***************************************************************************
     X=0.D0
     GRAD(:)=0.D0
     DO IP=1,NP
       CALL MURNAGHAN(PARMS,V(IP),EFIT,DEDP)
       X=X+(EFIT-E(IP))**2
       GRAD(:)=GRAD(:)+2.D0*(EFIT-E(IP))*DEDP(:)
     ENDDO
     X=X/REAL(NP)
     GRAD=GRAD/REAL(NP)
     RETURN
     END
!
!    ...........................................................................
     SUBROUTINE MURNAGHAN(PARMS,V,E,GRAD)
!    ***************************************************************************
!    **   MURNAGHAN EQUATION OF STATE:  MURNAGHAN44_PNAS30_244                **
!    **                                                                       **
!    **   MURNAGHAN'S EQUATION OF STATE IS BASED ON THE ASSUMPTION            **
!    **   THAT THE BULK MODULUS DEPENDS LINEARLY ON PRESSURE                  **
!    **   $P=-\FRAC{\PARTIAL E}{\PARTIAL V}$                                  **
!    **   $B=\FRAC{1}{\KAPPA_T}=-V\FRAC{\PARTIAL P}{\PARTIAL V}$              **
!    **   $B'=\FRAC{\PARTIAL B}{\PARTIAL P}\EQUIV B_0'\APPROX 3.5$            **
!    **   SUBROUTINE RETURNS ENERGY AND GRADIENT AT REQUESTED VOLUME          **
!    **                                                                       **
!    **   SEE BLOECHL LECTURE NOTES FOR THE PAW PRACTIKUM FOR A DERIVATION    **
!    *************************PETER BLOECHL, GOSLAR 2010************************
     IMPLICIT NONE
     REAL(8),INTENT(IN) :: PARMS(4) ! PARAMETERS OF MURN CURVE: E0, V0, B0, B'
     REAL(8),INTENT(IN) :: V        ! INPUT VOLUME
     REAL(8),INTENT(OUT):: E        ! OUTPUT ENERGY
     REAL(8),INTENT(OUT):: GRAD(4)  ! GRADIENT OF MURN CURVE AT INPUT VOLUME
     REAL(8)            :: VREL     ! V/V0
     REAL(8)            :: E0,V0,B0,BPRIME
     REAL(8)            :: DEDE0,DEDV0,DEDB0,DEDBPRIME ! DERIVATIVES
     REAL(8)            :: SVAR1,SVAR2
!    ***************************************************************************
     E0=PARMS(1)
     V0=PARMS(2)
     B0=PARMS(3)
     BPRIME=PARMS(4)
     SVAR1=B0*V0/BPRIME
     SVAR2=SVAR1/(BPRIME-1.D0)
     VREL=V/V0
     E        =E0+SVAR2*(VREL**(-BPRIME+1.D0)-1.D0)+SVAR1*(VREL-1.D0)
     DEDE0    =1.D0
     DEDV0    =(E-E0)/V0+(SVAR2*(-BPRIME+1.D0)*VREL**(-BPRIME)+SVAR1)*(-VREL/V0)
     DEDB0    =(E-E0)/B0
     DEDBPRIME=-(E-E0)/BPRIME &
    &          -SVAR2/(BPRIME-1.D0)*(VREL**(-BPRIME+1.D0)-1.D0) &
    &          -SVAR2*VREL**(-BPRIME+1.D0)*LOG(VREL)
     GRAD(1)=DEDE0
     GRAD(2)=DEDV0
     GRAD(3)=DEDB0
     GRAD(4)=DEDBPRIME
     RETURN
     END
