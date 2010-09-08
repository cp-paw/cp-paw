     PROGRAM MAIN
!    ***************************************************************************
!    ** DETERMINE EQUILIBRIUM VOLUME, EQUILIBRIUM ENERGY, BULK MODULUS        **
!    ** FROM A MURNAGHAN FIT TO A SET OF (VOLUME,ENERGY) DATA RECEIVED FROM   **
!    ** STANDARD INPUT                                                        **
!    ***************************MATTHE UITTEVAAL AND PETER E. BLOECHL***********
     IMPLICIT NONE
     integer(4),parameter :: npx=100
     REAL(8)    :: V(npx), E(npx)     ! INPUT VOLUMES, ENERGIES
     INTEGER(4) :: np                 ! #(datasets)
     INTEGER(4) :: HIGH,LOW,I         ! INDICES, 
     INTEGER(4) :: IARRAY(1)          ! 1D ARRAY FOR MIN INDEX
     REAL(8)    :: x                  ! penalty
     REAL(8)    :: PARMS(4)           ! fit paramers (E0,V0,B0,Bprime)
     REAL(8)    :: GRAD(4)            ! gradients of penalty function
     REAL(8)    :: VI,EFIT            ! VOLUME IN, ENERGY OUT
!    ***************************************************************************
!    ==========================================================================
!    == START WITH PRINTING HEADER                                           ==
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
!    == READ DATA UNTIL ERROR THEN RESET NP                                  ==
!    ==========================================================================
     np=0
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
!
!    == REPORT INPUT DATA =====================================================
     DO I=1,NP
       WRITE(*,FMT='(I5," V=",F10.5," E=",F10.5)')I,V(I),E(I)
     ENDDO
!
!    ==========================================================================
!    == ESTIMATE FIT PARAMETERS AS STARTING VALUES FOR OPTIMIZATION          ==
!    ==========================================================================
     IARRAY=MINLOC(V(1:NP))
     LOW=IARRAY(1) 
     IARRAY=MAXLOC(V(1:NP))
     HIGH=IARRAY(1) 
     IARRAY=MINLOC(E(1:NP)) 
     PARMS(1)=E(IARRAY(1)) ! MINE
     PARMS(2)=V(IARRAY(1)) ! VOL MINE
     PARMS(3)=2*(E(HIGH)+E(LOW)-2*PARMS(1))/(V(HIGH)-V(LOW))  !INITIAL B0
     PARMS(4)=3.5D0  ! APPROX CONST
!
!    ==========================================================================
!    == DO FIT USING QUENCHED DYNAMICS                                       ==
!    ==========================================================================
     CALL MD(NP,V,E,PARMS,X)
!
!    ==========================================================================
!    == PRINT PARAMETERS                                                     ==
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
     SUBROUTINE MD(NP,V,E,PARMS,X)
!    ***************************************************************************
!    ** CONJUGATE GRADIENT FIT OF THE MURNAGHAN CURVE (ITS PARAMETERS)        **
!    ** RETURNS THE MURNAGHAN PARAMETERS AND THE QUALITY OF THE FIT           **
!    **************************PETER BLOECHL, GOSLAR 2010***********************
     IMPLICIT NONE
     INTEGER(4),INTENT(IN)   :: NP
     REAL(8)   ,INTENT(IN)   :: V(NP)        ! VOLUMES OF INPUT DATA SET
     REAL(8)   ,INTENT(IN)   :: E(NP)        ! ENERGIES OF THE INPUT DATA SET
     REAL(8)   ,INTENT(INOUT):: PARMS(4)     ! MURNAGHAN PARAMETERS
     REAL(8)   ,INTENT(OUT)  :: X            ! QUALITY OF THE FIT
     INTEGER(4),PARAMETER    :: NITER=1000000
     REAL(8)   ,PARAMETER    :: ANNE=1.D-2
     REAL(8)   ,PARAMETER    :: DT=1.D-2
     REAL(8)   ,PARAMETER    :: TOL=1.D-8
     REAL(8)   ,PARAMETER    :: STEP=1.D-2
     REAL(8)                 :: M(4)
     REAL(8)                 :: Y0(4),YM(4),YP(4),GRAD(4),D2(4,4)
     REAL(8)                 :: XLAST
     INTEGER(4)              :: ITER,I
     REAL(8)                 :: SVAR1,SVAR2,SVAR3
     REAL(8)                 :: EKIN
     INTEGER(4)              :: NWRITE
!    ***************************************************************************
     Y0(:)=PARMS(:)
     CALL ONESTEP(NP,V,E,Y0,X,GRAD)
     DO I=1,4
       Y0(:)=PARMS(:)
       Y0(I)=Y0(I)+STEP
       CALL ONESTEP(NP,V,E,Y0,X,D2(:,I))
       D2(:,I)=(D2(:,I)-GRAD(:))/STEP
       M(I)=ABS(D2(I,I))
     ENDDO
!
     WRITE(*,FMT='(80("="),T10,"OPTIMIZE FIT PARAMETERS WITH MD")')
     Y0(:)=PARMS(:)
     YM=Y0
     XLAST=HUGE(X)
     NWRITE=10
     DO ITER=1,NITER
       CALL ONESTEP(NP,V,E,Y0,X,GRAD)
       IF(SUM(GRAD**2).LT.TOL) EXIT
       IF(X.GT.XLAST) YM=Y0
       XLAST=X
       SVAR1=2.D0/(1.D0+ANNE)
       SVAR2=1.D0-SVAR1
       SVAR3=DT**2/(1+ANNE)
       YP(:)=Y0(:)*SVAR1+YM(:)*SVAR2-GRAD(:)/M*SVAR3
       EKIN=0.5D0*SUM(M*(YP-YM)**2)/(2.D0*DT)**2
       IF(ITER.GT.10*NWRITE)NWRITE=10*NWRITE
       IF(MODULO(ITER,NWRITE).EQ.0) THEN
         WRITE(*,FMT='(I8," EKIN=",F12.5," PENALTY=",F12.7," CONSERVED=",F12.7)') &
      &      ITER,EKIN,X,EKIN+X
       END IF
       YM=Y0
       Y0=YP
     ENDDO
     WRITE(*,FMT='(I8," EKIN=",F12.5," PENALTY=",F12.7," CONSERVED=",F12.7)') &
    &      ITER,EKIN,XLAST,EKIN+XLAST
     PARMS=Y0
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
!    **   SEE BLOECHL LECTURE NOTES FOR THE PAW hands-on FOR DERIVATION       **
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
