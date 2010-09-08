     PROGRAM MAIN
!    ***************************************************************************
!    ** determine equilibrium volume, equilibrium energy, bulk modulus        **
!    ** from a murnaghan fit to a set of (Volume,energy) data received from   **
!    ** standard input                                                        **
!    **                                                                       **
!    ***************************************************************************
!    ***************************Matthe Uittevaal and Peter E. Bloechl***********
     IMPLICIT NONE
     REAL(8)    :: V(100), E(100)     ! input volumes, energies
     INTEGER(4) :: high,low,I, NP=0   ! indices, no. data points (init 0)
     INTEGER(4) :: IARRAY(1)          ! 1D array for min index
     REAL(8)    :: PARMS(4),GRAD(4),X ! parms, grad and mean sqr dev of murn fit
     REAL(8)    :: VI, EFIT           ! volume in, energy out
!    ***************************************************************************
!    ==========================================================================
!    == start with printing header
!    ==========================================================================
     WRITE(*,FMT='(60("="))')
     WRITE(*,FMT='(60("="),T15," FIT MURNAGHAN EQUATION OF STATE ")')
     WRITE(*,FMT='(60("="))')
!
!    ==========================================================================
!    == read data until error then reset NP
!    ==========================================================================
     DO 
       NP=NP+1
       IF (NP.EQ.101) THEN
         WRITE(*,*)"Too much input data, only first 100 used"
         EXIT
       END IF
       READ(*,*,END=100,ERR=100)V(NP),E(NP)
     ENDDO
100  CONTINUE
     NP=NP-1
     IF (NP.LT.4) THEN
       WRITE(*,*)"Too few input data, at least 4 required!"
       STOP
     END IF
!
!    ==========================================================================
!    == print data 
!    ==========================================================================
     DO I=1,NP
       WRITE(*,FMT='(I5," V=",F10.5," E=",F10.5)')I,V(I),E(I)
     ENDDO
     IARRAY=MINLOC(V(1:NP))
     low=IARRAY(1) 
     IARRAY=MAXLOC(V(1:NP))
     high=IARRAY(1) 
     WRITE(*,FMT='(" Vmin=",F10.5," Vmax=",F10.5)')V(low),V(high)
     IARRAY=MINLOC(E(1:NP)) 
!    I=MINLOC(V(1:NP))(1) ??
     PARMS(1)=E(IARRAY(1)) ! minE
     PARMS(2)=V(IARRAY(1)) ! vol minE
     WRITE(*,FMT='(" V(Emin)=",F10.5," Emin=",F10.5)')PARMS(2),PARMS(1)
     PARMS(3)=2*(E(high)+E(low)-2*PARMS(1))/(V(high)-V(low))  !initial B0
     PARMS(4)=3.5D0  ! approx const
     WRITE(*,FMT='(" Bini=",F10.5," BPRIMEini=",F10.5)')PARMS(3),PARMS(4)
!
!    ==========================================================================
!    == do fit
!    ==========================================================================
     CALL CG(NP,V,E,PARMS,X)
!
!    ==========================================================================
!    == print parameters
!    ==========================================================================
     WRITE(*,FMT='("FIT PARAMETERS OF MURNGAHAN EQUATION OF STATE")')
     WRITE(*,FMT='("EQUILIBRIUM ENERGY E0",T40,F10.5)')PARMS(1)
     WRITE(*,FMT='("EQUILIBRIUM VOLUME V0",T40,F10.5)')PARMS(2)
     WRITE(*,FMT='("EQUILIBRIUM BULK MODULUS B0",T40,F10.5)')PARMS(3)
     WRITE(*,FMT='("PRESSURE DERIVATIVE OF BULK MODULUS BPRIME",T40,F10.5)') &
   &       PARMS(4)
!
!    ===========================================================================
!    == print energy volume curve with input for comparison                   ==
!    ===========================================================================
     DO I=1,NP
       CALL MURNAGHAN(PARMS,V(I),EFIT,GRAD)
       WRITE(*,FMT='(I5," V=",F10.5," E=",F10.5," EFIT ",F10.5)') &
    &        I,V(I),E(I),EFIT
     ENDDO
!
!    ===========================================================================
!    == write equidistant energy volume curve (input range +10%)              ==
!    ===========================================================================
     OPEN(UNIT=8,FILE='murn.dat')
     DO I=-10,110
       VI=V(low)+(V(high)-V(low))/REAL(100)*REAL(I-1)
       CALL MURNAGHAN(PARMS,VI,EFIT,GRAD)
       WRITE(8,FMT='(2F10.5)')VI,EFIT
     ENDDO
     STOP
     END
!
!    ...........................................................................
     SUBROUTINE CG(NP,V,E,PARMS,X)
!    **                                                                       **
!    ** SUBROUTINE does CG fitting of the murnaghan curve (its parameters)    **
!    ** returns the murnaghan parameters and the quality of the fit           **
!    **                                                                       **
     IMPLICIT NONE
     INTEGER(4),INTENT(IN)   :: NP
     REAL(8)   ,INTENT(IN)   :: V(NP)
     REAL(8)   ,INTENT(IN)   :: E(NP)
     REAL(8)   ,INTENT(INOUT):: PARMS(4)     !murnaghan parameters
     REAL(8)   ,INTENT(OUT)  :: X            !quality of the fit
     REAL(8)   ,PARAMETER    :: STEP=1.D-4, FACT=.13, TOL=7.D-5
     REAL(8)   ,DIMENSION(4) :: PARMS1,PARMS2,GRAD1,GRAD2,LAPL
     LOGICAL(4)              :: TCONV
     INTEGER                 :: LOOP=0,LOOPMAX=1000
!    ***************************************************************************
     OPEN(7,FILE='ITER.DAT')
     CALL ONESTEP(NP,V,E,PARMS,X,GRAD1)
     WRITE(7,FMT='("========= PARAMETERS AND PENALTY FUNCTION OF FIT =======")')
     WRITE(7,FMT='("ITER",I4,":",4F10.5,"::",E10.2)')LOOP, PARMS, X
     PARMS1(:)=PARMS(:)*(1.+STEP) !percentage (DIMENSIONS,energy 0??)
     CALL ONESTEP(NP,V,E,PARMS1,X,GRAD2)
     LAPL(:)=(GRAD2(:)-GRAD1(:))/PARMS(:)/STEP
     WRITE(7,FMT='("          gradient and laplacian of penalty function ")')
     WRITE(7,FMT='(4E10.2,"::",4E10.2)')GRAD1, LAPL
!
     PARMS1=PARMS
     DO LOOP=1,LOOPMAX
       PARMS2=PARMS1-GRAD1/LAPL*FACT  !FACT towards minimum
       TCONV=(SQRT(SUM((GRAD1/LAPL/PARMS1)**2)).LT.TOL) 
       IF(TCONV) EXIT
       CALL ONESTEP(NP,V,E,PARMS2,X,GRAD2)
       WRITE(7,FMT='("======= PARAMETERS AND PENALTY FUNCTION OF FIT =======")')
       WRITE(7,FMT='("ITER",I4,":",4F10.5,"::",E10.5)')LOOP, PARMS2, X
       LAPL(:)=(GRAD2(:)-GRAD1(:))/(PARMS2-PARMS1)
       WRITE(7,FMT='("        gradient and laplacian of penalty function ___")')
       WRITE(7,FMT='(4E10.2,"::"4E10.2)')GRAD2, LAPL
       GRAD1=GRAD2
       PARMS1=PARMS2
     ENDDO
     IF (.NOT.TCONV) THEN
       WRITE(7,FMT='("-----CONJUGATE GRADIENT NOT CONVERGED---------")')
     END IF
     CLOSE(7)
     PARMS=PARMS2
     RETURN
     END
!
!    ...........................................................................
     SUBROUTINE ONESTEP(NP,V,E,PARMS,X,GRAD)
!    **                                                                       **
!    ** SUBROUTINE calculates how the fit changes by a change of parameters   **
!    ** returns penalty function and derivative                               **
!    **                                                                       **
     IMPLICIT NONE
     INTEGER(4),INTENT(IN) :: NP           !#(data)
     REAL(8)   ,INTENT(IN) :: V(NP),E(NP)  !volumes, energies
     REAL(8)   ,INTENT(IN) :: PARMS(4)     !current fit parameters  
     REAL(8)   ,INTENT(OUT):: X,GRAD(4)    !value and gradient of penalty func.
     REAL(8)               :: DEDP(4)      !gradient of murnaghan curve
     REAL(8)               :: EFIT         !fitting energy
     INTEGER(4)            :: IP           !counter
!    ***************************************************************************
     X=0.D0
     GRAD(:)=0.D0
     DO IP=1,NP
       CALL MURNAGHAN(PARMS,V(IP),EFIT,DEDP)
       X=X+(EFIT-E(IP))**2
       GRAD(:)=GRAD(:)+(EFIT-E(IP))*DEDP(:)
     ENDDO
     X=X/REAL(NP)
     GRAD=GRAD/REAL(NP)
     RETURN
     END
!
!    ...........................................................................
     SUBROUTINE MURNAGHAN(PARMS,V,E,GRAD)
!    **                                                                       **
!    **   MURNAGHAN EQUATION OF STATE:  MURNAGHAN44_PNAS30_244                **
!    **                                                                       **
!    **   MURNAGHAN'S EQUATION OF STATE IS BASED ON THE ASSUMPTION            **
!    **   THAT THE BULK MODULUS DEPENDS LINEARLY ON PRESSURE                  **
!    **   $p=-\frac{\partial E}{\partial V}$                                  **
!    **   $B=\frac{1}{\kappa_T}=-V\frac{\partial p}{\partial V}$              **
!    **   $B'=\frac{\partial B}{\partial p}\equiv B_0'\approx 3.5$            **
!    **   SUBROUTINE returns energy and gradient at requested volume          **
!    **                                                                       **
     IMPLICIT NONE
     REAL(8),INTENT(IN) :: PARMS(4) ! parameters of murn curve: E0, V0, B0, B'
     REAL(8),INTENT(IN) :: V        ! input volume
     REAL(8),INTENT(OUT):: E        ! output energy
     REAL(8),INTENT(OUT):: GRAD(4)  ! gradient of murn curve at input volume
     REAL(8)            :: IB,IB1,F1,G1  !tmp var
!    ***************************************************************************
     F1=(PARMS(2)/V)**PARMS(4)
     IB=1.D0/PARMS(4)         !inverse B'
     IB1=1.D0/(PARMS(4)-1.D0) !inverse B'-1
     G1=PARMS(3)*V*IB*IB1
     E=PARMS(1) + G1 * (F1+PARMS(4)-1.D0) - PARMS(3)*PARMS(2)*IB1
!
     GRAD(1)=1.D0                                    !dE/dE0
     GRAD(2)=PARMS(3)*IB1 * ((PARMS(2)/V)**(PARMS(4)-1.D0) - 1.D0) !dE/dV0
     GRAD(3)=(E-PARMS(1)) / PARMS(3)                 !dE/dB0
     GRAD(4)=G1 * (IB+IB1+F1 * (LOG(PARMS(2)/V) - IB1))    !dE/dB'
     RETURN
     END
