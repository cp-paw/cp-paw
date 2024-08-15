     PROGRAM MAIN
!    ***************************************************************************
!    ** DETERMINE EQUILIBRIUM VOLUME, EQUILIBRIUM ENERGY, BULK MODULUS        **
!    ** FROM A MURNAGHAN FIT TO A SET OF (VOLUME,ENERGY) DATA RECEIVED FROM   **
!    ** STANDARD INPUT                                                        **
!    ***************************MATTHE UITTEVAAL AND PETER E. BLOECHL***********
     USE STRINGS_MODULE
     IMPLICIT NONE
     INTEGER(4),PARAMETER :: NPX=100
     REAL(8)    :: V(NPX), E(NPX)     ! INPUT VOLUMES, ENERGIES
     INTEGER(4) :: NP                 ! #(DATASETS)
     INTEGER(4) :: HIGH,LOW,I         ! INDICES, 
     INTEGER(4) :: IARRAY(1)          ! 1D ARRAY FOR MIN INDEX
     REAL(8)    :: X                  ! PENALTY
     REAL(8)    :: PARMS(4)           ! FIT PARAMERS (E0,V0,B0,BPRIME)
     REAL(8)    :: GRAD(4)            ! GRADIENTS OF PENALTY FUNCTION
     REAL(8)    :: VI,EFIT            ! VOLUME IN, ENERGY OUT
     REAL(8)    :: EUNIT              ! ENERGY UNIT IN HARTREE
     REAL(8)    :: LUNIT              ! LENGTH UNIT IN ATOMIC UNITS
     REAL(8)    :: VBYL3              ! VOLUME IS VBYL3*L**3
     REAL(8)    :: SVAR               ! AUXILIARY VARIABLE
     REAL(8)    :: SCALE              ! SCALES THE RESULT
     CHARACTER(64) :: ARG     
     CHARACTER(250) :: STRING     
     CHARACTER(250) :: LINE
     LOGICAL(4) :: TL                 ! INPUT IS LENGTH INSTEAD OF VOLUME
     REAL(8)   ,PARAMETER :: ANGSTROM=1.D0/0.52917721092D0
     REAL(8)   ,PARAMETER :: JOULE=1.D0/4.35974434D-18
     REAL(8)   ,PARAMETER :: METER=ANGSTROM*1.D+10
     REAL(8)   ,PARAMETER :: GPASCAL=1.D+9*JOULE/METER**3
!    ***************************************************************************
!    ===========================================================================
!    == START WITH PRINTING HEADER                                            ==
!    ===========================================================================
     WRITE(*,FMT='(80("="))')
     WRITE(*,FMT='(80("="),T15," FIT MURNAGHAN EQUATION OF STATE ")')
     WRITE(*,FMT='(80("="))')
     WRITE(*,FMT='(A80)') &
    &     'MURNAGHAN EQUATION OF STATE IS BASED ON THE ASSUMPTION THAT THE' &
    &    ,'BULK MODULUS DEPENDS LINEARLY ON PRESSURE.' &
    &    ,'SOURCE: F.D. MURNAGHAN, PNAS 30, 244 (1934)'
!
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!    ==========================================================================
!    == COLLECT UNITS FROM COMMAND LINE                                       ==
!    ==========================================================================
     EUNIT=1.D0
     LUNIT=1.D0
     SCALE=1.D0
     TL=.FALSE.
     VBYL3=1.D0
     I=0
     DO 
       I=I+1
       CALL GET_COMMAND_ARGUMENT(I,ARG)
       IF(LEN_TRIM(ARG).EQ.0) EXIT
       ARG=ADJUSTL(ARG)
       IF(ARG.EQ.-'-H') THEN
         WRITE(*,FMT=-'("PAW_MURNAGHAN.X OPTIONS < INPUT")')
         WRITE(*,FMT='("OPTIONS:")')
         WRITE(*,FMT=-'(T5,"-L",T30 &
    &                  ,"FIRST COLUMN OF INPUT IS A LATTICE CONSTANT")')
         WRITE(*,FMT=-'(T5,"-V",T30 &
    &                 ,"FIRST COLUMN OF INPUT IS A VOLUME (DEFAULT)")')
         WRITE(*,FMT=-'(T5,"-EU VALUE",T30 &
    &                    ,"ENERGY UNIT OF SECOND COLUMN IN HARTREE")')
         WRITE(*,FMT=-'(T5,"-LU VALUE",T30 &
    &                    ,"LENGTH UNIT OF FIRST COLUMN IN ABOHR")')
         WRITE(*,FMT=-'(T5,"-VBL VALUE",T30 &
    &                    ,"VOLUME / LATTICE CONSTANT^3")')
         WRITE(*,FMT='("INPUT CONTAINS TWO COLUMNS WITH DATA:")')
         WRITE(*,FMT=-'(T5,"FIRST COLUMN:",T30,"LATTICE CONSTANT OR VOLUME")')
         WRITE(*,FMT=-'(T5,"SECOND COLUMN:",T30,"ENERGY")')
         STOP 'STOPPING AFTER REPORTING HELP INFORMATION'
       ELSE IF(ARG.EQ.-'-L') THEN
         TL=.TRUE.
       ELSE IF(ARG.EQ.-'-V') THEN
         TL=.FALSE.
       ELSE IF(ARG.EQ.-'-EU') THEN
         I=I+1
         CALL GET_COMMAND_ARGUMENT(I,ARG)
         READ(ARG,*) EUNIT
       ELSE IF(ARG.EQ.-'-LU') THEN
         I=I+1
         CALL GET_COMMAND_ARGUMENT(I,ARG)
         READ(ARG,*) LUNIT
       ELSE IF(ARG.EQ.-'-VBL') THEN
         I=I+1
         CALL GET_COMMAND_ARGUMENT(I,ARG)
         READ(ARG,*)VBYL3
       ELSE IF(ARG.EQ.-'-SCALE') THEN
         I=I+1
         CALL GET_COMMAND_ARGUMENT(I,ARG)
         READ(ARG,*)SCALE
       ELSE
         WRITE(*,*) I,'TH ARGUMENT NOT RECOGNIZED'
         WRITE(*,*) 'ARGUMENT VALUE: ',TRIM(ARG)
         STOP 'IN MAIN'
       END IF
     ENDDO
     WRITE(*,FMT='(50("."),T1,"ENERGY UNIT IN HARTREE",T50,F10.5)')EUNIT
     WRITE(*,FMT='(50("."),T1,"LENGTH UNIT IN BOHR RADII",T50,F10.5)')LUNIT
     IF(TL) THEN
       WRITE(*,FMT='("FIRST COLUMN INTERPRETED AS LENGTH")')
     ELSE
       WRITE(*,FMT='("FIRST COLUMN INTERPRETED AS VOLUME")')
     END IF
     WRITE(*,FMT='(50("."),T1,"CELL-VOLUME / LENGTH**3",T50,F10.5)')VBYL3
!
!    ==========================================================================
!    == READ DATA UNTIL ERROR, THEN RESET NP                                 ==
!    ==========================================================================
     NP=0
     DO 
       NP=NP+1
       IF (NP.EQ.101) THEN
         WRITE(*,*)"TOO MANY INPUT DATA, ONLY FIRST 100 USED"
         EXIT
       END IF
       READ(*,FMT='(A)',END=100,ERR=100)LINE
       LINE=ADJUSTL(LINE)
       IF(LINE(1:1).EQ.'#') THEN ! SKIP COMMENT LINES
         NP=NP-1
         CYCLE
       END IF
       READ(LINE,*,END=100,ERR=100)V(NP),E(NP)
       IF(TL) V(NP)=VBYL3*V(NP)**3   ! FIRST, V=ALAT, THEN V=VOLUME OPER CELL
     ENDDO
100  CONTINUE
     NP=NP-1
     IF (NP.LT.4) THEN
       WRITE(*,*)"TOO FEW INPUT DATA, AT LEAST 4 ARE REQUIRED!"
       STOP
     END IF
!
!    ==========================================================================
!    == CONVERT UNITS ACCORDING TO COMMAND LINE ARGUMENTS                    ==
!    ==========================================================================
     E(:)=E(:)*EUNIT
     V(:)=V(:)*LUNIT**3
!
!    ==========================================================================
!    == SCALE RESULTS
!    ==========================================================================
     E(:)=E(:)*SCALE
     V(:)=V(:)*SCALE
     VBYL3=VBYL3*SCALE
!
!    ==========================================================================
!    == REPORT INPUT DATA =====================================================
!    ==========================================================================
     DO I=1,NP
       SVAR=(V(I)/VBYL3)**(1.D0/3.D0)
       WRITE(*,FMT='(I5," L=",F15.5," A0 V=",F15.5," A0^3 E=",F10.5," H")') &
    &          I,SVAR,V(I),E(I)
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
     WRITE(*,FMT='("EQUILIBRIUM ENERGY E0 IN HARTREE",T50,F10.5)')PARMS(1)
     WRITE(*,FMT='("EQUILIBRIUM VOLUME V0 IN ABOHR^3",T50,F10.5)')PARMS(2)
     WRITE(*,FMT='("EQUILIBRIUM VOLUME V0 IN ANGSTROM^3",T50,F10.5)') &
    &                                    PARMS(2)/ANGSTROM**3
     SVAR=(PARMS(2)/VBYL3)**(1.D0/3.D0)/ANGSTROM
     WRITE(*,FMT='("EQUILIBRIUM LATTICE CONSTANT IN ANGSTROM",T50,F10.5)')SVAR
     WRITE(*,FMT='("EQUILIBRIUM BULK MODULUS B0 IN A.U.",T50,F10.5)')PARMS(3)
     SVAR=PARMS(3)/GPASCAL
     WRITE(*,FMT='("EQUILIBRIUM BULK MODULUS B0 IN GPA",T50,F10.5)')SVAR
     WRITE(*,FMT='("PRESSURE DERIVATIVE OF BULK MODULUS BPRIME",T50,F10.5)') &
   &       PARMS(4)
!
!    ===========================================================================
!    == PRINT ENERGY VOLUME CURVE WITH INPUT FOR COMPARISON                   ==
!    ===========================================================================
     WRITE(*,*)
     WRITE(*,FMT='(80("="),T10,"COMPARE ORIGINAL AND INTERPOLATED DATA")')
     DO I=1,NP
       CALL MURNAGHAN(PARMS,V(I),EFIT,GRAD)
       WRITE(*,FMT='(I5," V=",F10.5," E(IN)=",F10.5," E(FIT)=",F10.5' &
    &            //'," E(FIT)-E(IN))=",F10.5)') &
    &        I,V(I),E(I),EFIT,EFIT-E(I)
     ENDDO
!
!    ===========================================================================
!    == WRITE EQUIDISTANT ENERGY VOLUME CURVE (INPUT RANGE +10%)              ==
!    ===========================================================================
     STRING=-'INTERPOLATED EQUATION OF STATE IS WRITTEN TO FILE MURN.DAT'
     STRING='(80("="),T10,"'//TRIM(STRING)//'")'
     WRITE(*,FMT=TRIM(STRING))
!     WRITE(*,FMT=-'(80("="),T10,"INTERPOLATED EQUATION OF STATE IS WRITTEN' &
!    &                         //-' TO FILE MURN.DAT")')
     OPEN(UNIT=8,FILE=-'MURN.DAT')
!     OPEN(UNIT=9,FILE=-'EOFV.DAT')
     DO I=-10,110
       VI=V(LOW)+(V(HIGH)-V(LOW))/REAL(100)*REAL(I-1)
       CALL MURNAGHAN(PARMS,VI,EFIT,GRAD)
!       WRITE(9,FMT='(2F20.5)')VI/LUNIT**3,EFIT/EUNIT
!      __ CONVERT CONSISTENT WITH THE INPUT DATA________________________________
       EFIT=EFIT/EUNIT
       VI=VI/LUNIT**3
!      __ WRITE_________________________________________________________________
       IF(TL) THEN
!        __FIRST COLUMN IS THE LATTICE CONSTANT_________________________________
         WRITE(8,FMT='(3F15.10)')(VI/VBYL3)**(1.D0/3.D0),EFIT
       ELSE
!        __FIRST COLUMN IS THE VOLUME PER UNIT CELL_____________________________
!        __SECOND COLUMN IS THE ENERGY PER UNIT CELL____________________________
!        __THIRD COLUMN IS THE LATTICE CONSTANT_________________________________
         WRITE(8,FMT='(3F15.10)')VI,EFIT,(VI/VBYL3)**(1.D0/3.D0)
       END IF
     ENDDO
     CALL ERROR$NORMALSTOP()
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
!    **   MURNAGHAN EQUATION OF STATE IS BASED ON THE ASSUMPTION              **
!    **   THAT THE BULK MODULUS DEPENDS LINEARLY ON PRESSURE                  **
!    **   $P=-\FRAC{\PARTIAL E}{\PARTIAL V}$                                  **
!    **   $B=\FRAC{1}{\KAPPA_T}=-V\FRAC{\PARTIAL P}{\PARTIAL V}$              **
!    **   $B'=\FRAC{\PARTIAL B}{\PARTIAL P}\EQUIV B_0'\APPROX 3.5$            **
!    **   SUBROUTINE RETURNS ENERGY AND GRADIENT AT REQUESTED VOLUME          **
!    **                                                                       **
!    **   SEE BLOECHL LECTURE NOTES FOR THE PAW HANDS-ON FOR DERIVATION       **
!    *************************PETER BLOECHL, GOSLAR 2010************************
     IMPLICIT NONE
     REAL(8),INTENT(IN) :: PARMS(4) ! PARAMETERS OF MURN CURVE: E0,V0,B0,BPRIME
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
