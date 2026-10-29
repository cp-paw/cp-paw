! Native SCAN kernel contributed by Konstantin Tamoev.
! Based on Sun, Ruzsinszky and Perdew, Phys. Rev. Lett. 115, 036402 (2015).
!*******************************************************************************
!*******************************************************************************
!****  SCAN META-GGA FUNCTIONAL (NATIVE IMPLEMENTATION)                     ****
!****                                                                       ****
!****  SUN, RUZSINSZKY, PERDEW, PHYS.REV.LETT. 115, 036402 (2015)           ****
!****  FORMULAS FROM THE PAPER AND ITS SUPPLEMENTAL MATERIAL.               ****
!****  THE LSDA CORRELATION EPSILON_C^LSDA1 IS THE PERDEW-WANG 1992         ****
!****  PARAMETRIZATION (PHYS.REV.B 45, 13244 (1992)).                       ****
!****                                                                       ****
!****  STATE VARIABLES (META-GGA SET OF THE PAW_DFT OBJECT):                ****
!****    VAL(1)=RHOT  =TOTAL DENSITY                                        ****
!****    VAL(2)=RHOS  =SPIN DENSITY (NUP-NDN)                               ****
!****    VAL(3)=GRHOT2=(GRAD RHOT)**2                                       ****
!****    VAL(4)=GRHOS2=(GRAD RHOS)**2                                       ****
!****    VAL(5)=GRHOST=(GRAD RHOT)*(GRAD RHOS)                              ****
!****    VAL(6)=LAPLT =LAPLACIAN OF RHOT (NOT USED BY SCAN)                 ****
!****    VAL(7)=LAPLS =LAPLACIAN OF RHOS (NOT USED BY SCAN)                 ****
!****    VAL(8)=TAUT  =TOTAL KINETIC ENERGY DENSITY SUM_N |GRAD PSI_N|**2/2 ****
!****    VAL(9)=TAUS  =SPIN KINETIC ENERGY DENSITY (TAUUP-TAUDN)            ****
!****                                                                       ****
!****  SCAN$EVAL1 RETURNS EXC=RHOT*EPSILON_XC AND ITS ANALYTIC FIRST        ****
!****  DERIVATIVES. SCAN$EVAL2/SCAN$EVAL3 OBTAIN SECOND AND THIRD           ****
!****  DERIVATIVES BY FINITE DIFFERENCES OF THE ANALYTIC FIRST DERIVATIVES. ****
!****                                                                       ****
!****  THE EXCHANGE PART USES THE EXACT SPIN-SCALING RELATION               ****
!****     EX[NUP,NDN]=(EX[2*NUP]+EX[2*NDN])/2                               ****
!****  (OLIVER,PERDEW, PHYS.REV.A 20,397 (1979))                            ****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE SCAN_MODULE
!     **************************************************************************
!     ** PARAMETERS OF THE SCAN FUNCTIONAL                                    **
!     ** EXCHANGE CONSTANTS FROM PRL 115, 036402 (2015), EQS.(5)-(9)          **
!     ** CORRELATION CONSTANTS FROM THE SUPPLEMENTAL MATERIAL                 **
!     **************************************************************************
      IMPLICIT NONE
!     == CONSTANT INITIALIZATION ALSO MAKES FIRST CALLS THREAD SAFE =============
      LOGICAL(4),PARAMETER:: TINI=.TRUE.
      REAL(8)   ,PARAMETER:: PI=4.D0*ATAN(1.D0)
!     == MINIMUM DENSITY: BELOW THIS THE FUNCTIONAL RETURNS ZERO ===============
      REAL(8)   ,PARAMETER:: RHOMIN=1.D-12
!     == EXCHANGE ==============================================================
      REAL(8)   ,PARAMETER:: K1=0.065D0
      REAL(8)   ,PARAMETER:: H0X=1.174D0
      REAL(8)   ,PARAMETER:: A1=4.9479D0
      REAL(8)   ,PARAMETER:: C1X=0.667D0
      REAL(8)   ,PARAMETER:: C2X=0.8D0
      REAL(8)   ,PARAMETER:: DX=1.24D0
      REAL(8)   ,PARAMETER:: MUAK=10.D0/81.D0
      REAL(8)   ,PARAMETER:: B3=0.5D0
      REAL(8)   ,PARAMETER:: B2=SQRT(5913.D0/405000.D0)
      REAL(8)   ,PARAMETER:: B1=(511.D0/13500.D0)/(2.D0*B2)
      REAL(8)   ,PARAMETER:: B4=MUAK**2/K1-1606.D0/18225.D0-B1**2
!     == CORRELATION ===========================================================
      REAL(8)   ,PARAMETER:: B1C=0.0285764D0
      REAL(8)   ,PARAMETER:: B2C=0.0889D0
      REAL(8)   ,PARAMETER:: B3C=0.125541D0
      REAL(8)   ,PARAMETER:: C1C=0.64D0
      REAL(8)   ,PARAMETER:: C2C=1.5D0
      REAL(8)   ,PARAMETER:: DC=0.7D0
      REAL(8)   ,PARAMETER:: GAMMAC=(1.D0-LOG(2.D0))/PI**2
!     == USE LIBXC'S SCAN/PW92 PARAMETER CONVENTION, INCLUDING GC=2.363 ===========
      REAL(8)   ,PARAMETER:: BETA0=0.066724550603149220D0
      REAL(8)   ,PARAMETER:: GCPREF=2.363D0
!     == CHI_INFINITY OF G_INFINITY(ZETA=0,S), SEE SUPPL.MAT. OF SCAN PAPER ====
      REAL(8)   ,PARAMETER:: CHIINF=0.12802585262625815D0
!     == DERIVED CONSTANTS =====================================================
      REAL(8)   ,PARAMETER:: C23=(3.D0*PI**2)**(2.D0/3.D0)
      REAL(8)   ,PARAMETER:: CEX=(3.D0/(4.D0*PI))*(3.D0*PI**2)**(1.D0/3.D0)
      REAL(8)   ,PARAMETER:: CT2=(3.D0*PI**2/16.D0)**(2.D0/3.D0)
      REAL(8)   ,PARAMETER:: CRS=(3.D0/(4.D0*PI))**(1.D0/3.D0)
!     == PERDEW-WANG 92 PARAMETERS (A,ALPHA1,BETA1,BETA2,BETA3,BETA4) ==========
!     == FOR EC(RS,0), EC(RS,1) AND -ALPHAC(RS)                               ==
      REAL(8)   ,PARAMETER:: PWU(6)=(/0.0310907D0,0.21370D0,7.5957D0 &
     &                               ,3.5876D0,1.6382D0,0.49294D0/)
      REAL(8)   ,PARAMETER:: PWP(6)=(/0.01554535D0,0.20548D0,14.1189D0 &
     &                               ,6.1977D0,3.3662D0,0.62517D0/)
      REAL(8)   ,PARAMETER:: PWA(6)=(/0.0168869D0,0.11125D0,10.357D0 &
     &                               ,3.6231D0,0.88026D0,0.49671D0/)
      REAL(8)   ,PARAMETER:: FDD0=4.D0/(9.D0*(2.D0**(1.D0/3.D0)-1.D0))
      END MODULE SCAN_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN_INITIALIZE()
!     **************************************************************************
!     ** COMPATIBILITY ENTRY POINT. ALL CONSTANTS ARE INITIALIZED ABOVE.     **
!     **************************************************************************
      USE SCAN_MODULE
      IMPLICIT NONE
!     **************************************************************************
      RETURN
      END SUBROUTINE SCAN_INITIALIZE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN_PW92(RS,ZETA,EC,DEC_DRS,DEC_DZ)
!     **************************************************************************
!     ** PERDEW-WANG 1992 PARAMETRIZATION OF THE LSDA CORRELATION ENERGY      **
!     ** PER PARTICLE, PHYS.REV.B 45, 13244 (1992)                            **
!     ** RETURNS EC(RS,ZETA) [HARTREE] AND ITS DERIVATIVES                    **
!     **************************************************************************
      USE SCAN_MODULE, ONLY : PWU,PWP,PWA,FDD0
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RS
      REAL(8)   ,INTENT(IN) :: ZETA
      REAL(8)   ,INTENT(OUT):: EC
      REAL(8)   ,INTENT(OUT):: DEC_DRS
      REAL(8)   ,INTENT(OUT):: DEC_DZ
      REAL(8)               :: EU,DEU,EP,DEP,ALFM,DALFM
      REAL(8)               :: F,DF,Z3,Z4,OPZ,OMZ
!     **************************************************************************
      CALL SCAN_PW92G(RS,PWU,EU,DEU)
      CALL SCAN_PW92G(RS,PWP,EP,DEP)
      CALL SCAN_PW92G(RS,PWA,ALFM,DALFM)  ! ALFM=-ALPHAC (SPIN STIFFNESS)
      OPZ=1.D0+ZETA
      OMZ=1.D0-ZETA
      F=(OPZ**(4.D0/3.D0)+OMZ**(4.D0/3.D0)-2.D0) &
     &  /(2.D0**(4.D0/3.D0)-2.D0)
      DF=(4.D0/3.D0)*(OPZ**(1.D0/3.D0)-OMZ**(1.D0/3.D0)) &
     &  /(2.D0**(4.D0/3.D0)-2.D0)
      Z3=ZETA**3
      Z4=ZETA**4
      EC     =EU-ALFM*F/FDD0*(1.D0-Z4)+(EP-EU)*F*Z4
      DEC_DRS=DEU-DALFM*F/FDD0*(1.D0-Z4)+(DEP-DEU)*F*Z4
      DEC_DZ =-ALFM/FDD0*(DF*(1.D0-Z4)-4.D0*Z3*F) &
     &        +(EP-EU)*(DF*Z4+4.D0*Z3*F)
      RETURN
      END SUBROUTINE SCAN_PW92
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN_PW92G(RS,P,G,DG_DRS)
!     **************************************************************************
!     ** AUXILIARY FUNCTION G(RS) OF THE PERDEW-WANG 92 PARAMETRIZATION       **
!     ** P=(A,ALPHA1,BETA1,BETA2,BETA3,BETA4)                                 **
!     **   G=-2A(1+ALPHA1*RS)*LN(1+1/(2A*(B1*SQRT(RS)+B2*RS                  **
!     **                                 +B3*RS**1.5+B4*RS**2)))              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RS
      REAL(8)   ,INTENT(IN) :: P(6)
      REAL(8)   ,INTENT(OUT):: G
      REAL(8)   ,INTENT(OUT):: DG_DRS
      REAL(8)               :: SQRS,Q0,Q1,DQ1,ARG
!     **************************************************************************
      SQRS=SQRT(RS)
      Q0=-2.D0*P(1)*(1.D0+P(2)*RS)
      Q1=2.D0*P(1)*(P(3)*SQRS+P(4)*RS+P(5)*RS*SQRS+P(6)*RS**2)
      DQ1=P(1)*(P(3)/SQRS+2.D0*P(4)+3.D0*P(5)*SQRS+4.D0*P(6)*RS)
      ARG=1.D0+1.D0/Q1
      G=Q0*LOG(ARG)
      DG_DRS=-2.D0*P(1)*P(2)*LOG(ARG)-Q0*DQ1/(Q1**2+Q1)
      RETURN
      END SUBROUTINE SCAN_PW92G
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN_INTERPOLATE(ALPHA,C1,C2,D,F,DF_DA)
!     **************************************************************************
!     ** SCAN INTERPOLATION/EXTRAPOLATION FUNCTION OF ALPHA                   **
!     **    F(ALPHA)=EXP(-C1*ALPHA/(1-ALPHA))          FOR ALPHA<1            **
!     **    F(ALPHA)=-D*EXP(C2/(1-ALPHA))              FOR ALPHA>1            **
!     ** USED FOR FX (C1X,C2X,DX) AND FC (C1C,C2C,DC).                        **
!     ** THE FUNCTION AND ALL ITS DERIVATIVES VANISH AT ALPHA=1.              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: ALPHA
      REAL(8)   ,INTENT(IN) :: C1
      REAL(8)   ,INTENT(IN) :: C2
      REAL(8)   ,INTENT(IN) :: D
      REAL(8)   ,INTENT(OUT):: F
      REAL(8)   ,INTENT(OUT):: DF_DA
      REAL(8)   ,PARAMETER  :: EPS=1.D-8
      REAL(8)               :: OMA
!     **************************************************************************
      OMA=1.D0-ALPHA
      IF(ALPHA.LT.1.D0-EPS) THEN
        F=EXP(-C1*ALPHA/OMA)
        DF_DA=-F*C1/OMA**2
      ELSE IF(ALPHA.GT.1.D0+EPS) THEN
        F=-D*EXP(C2/OMA)
        DF_DA=F*C2/OMA**2
      ELSE
!       == ESSENTIAL SINGULARITY AT ALPHA=1: FUNCTION AND DERIVATIVES VANISH ===
        F=0.D0
        DF_DA=0.D0
      END IF
      RETURN
      END SUBROUTINE SCAN_INTERPOLATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN_FX(P,ALPHA,FX,DFX_DP,DFX_DA)
!     **************************************************************************
!     ** SCAN EXCHANGE ENHANCEMENT FACTOR, EQS.(5)-(9) OF PRL 115,036402      **
!     **    FX={H1X(P,ALPHA)+FXA(ALPHA)*[H0X-H1X(P,ALPHA)]}*GX(P)             **
!     ** WITH P=S**2 AND ANALYTIC DERIVATIVES WITH RESPECT TO P AND ALPHA     **
!     **************************************************************************
      USE SCAN_MODULE, ONLY : K1,H0X,A1,C1X,C2X,DX,MUAK,B1,B2,B3,B4
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: P
      REAL(8)   ,INTENT(IN) :: ALPHA
      REAL(8)   ,INTENT(OUT):: FX
      REAL(8)   ,INTENT(OUT):: DFX_DP
      REAL(8)   ,INTENT(OUT):: DFX_DA
      REAL(8)               :: OMA          ! 1-ALPHA
      REAL(8)               :: EB3          ! EXP(-B3*(1-ALPHA)**2)
      REAL(8)               :: U,DU_DP,DU_DA! B1*P+B2*(1-ALPHA)*EB3
      REAL(8)               :: EB4          ! EXP(-B4*P/MUAK)
      REAL(8)               :: X,DX_DP,DX_DA
      REAL(8)               :: H1X,DH1X_DX
      REAL(8)               :: GX,DGX_DP
      REAL(8)               :: FXA,DFXA_DA
      REAL(8)               :: EXPO
!     **************************************************************************
!     ==========================================================================
!     == X OF EQ.(6): X=MUAK*P*[1+(B4*P/MUAK)*EXP(-|B4|*P/MUAK)]              ==
!     ==              +{B1*P+B2*(1-ALPHA)*EXP[-B3*(1-ALPHA)**2]}**2           ==
!     == NOTE: B4>0 SO THAT |B4|=B4                                           ==
!     ==========================================================================
      OMA=1.D0-ALPHA
      EB3=EXP(-B3*OMA**2)
      U=B1*P+B2*OMA*EB3
      DU_DP=B1
      DU_DA=B2*EB3*(2.D0*B3*OMA**2-1.D0)
      EB4=EXP(-B4*P/MUAK)
      X=MUAK*P*(1.D0+(B4*P/MUAK)*EB4)+U**2
      DX_DP=MUAK+B4*P*EB4*(2.D0-P*B4/MUAK)+2.D0*U*DU_DP
      DX_DA=2.D0*U*DU_DA
!
!     ==========================================================================
!     == H1X OF EQ.(5): H1X=1+K1-K1/(1+X/K1)                                  ==
!     ==========================================================================
      H1X=1.D0+K1-K1/(1.D0+X/K1)
      DH1X_DX=1.D0/(1.D0+X/K1)**2
!
!     ==========================================================================
!     == GX OF EQ.(7): GX=1-EXP(-A1/P**(1/4))                                 ==
!     == FOR P->0 GX->1 AND DGX/DP->0 (ESSENTIAL SINGULARITY)                 ==
!     ==========================================================================
      IF(P.GT.1.D-12) THEN
        EXPO=A1/P**0.25D0
        IF(EXPO.LT.500.D0) THEN
          GX=1.D0-EXP(-EXPO)
          DGX_DP=-0.25D0*A1*EXP(-EXPO)/P**1.25D0
        ELSE
          GX=1.D0
          DGX_DP=0.D0
        END IF
      ELSE
        GX=1.D0
        DGX_DP=0.D0
      END IF
!
!     ==========================================================================
!     == FXA OF EQ.(9) AND ASSEMBLY OF FX (EQ.(8))                            ==
!     ==========================================================================
      CALL SCAN_INTERPOLATE(ALPHA,C1X,C2X,DX,FXA,DFXA_DA)
      FX=(H1X+FXA*(H0X-H1X))*GX
      DFX_DP=DH1X_DX*DX_DP*(1.D0-FXA)*GX+(H1X+FXA*(H0X-H1X))*DGX_DP
      DFX_DA=(DH1X_DX*DX_DA*(1.D0-FXA)+DFXA_DA*(H0X-H1X))*GX
      RETURN
      END SUBROUTINE SCAN_FX
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN_XUNP(N,GN2,TAU,AX,DAX_DN,DAX_DGN2,DAX_DTAU)
!     **************************************************************************
!     ** SCAN EXCHANGE ENERGY PER VOLUME OF THE SPIN-UNPOLARIZED SYSTEM       **
!     **    AX=N*EPSX_UNIF(N)*FX(P,ALPHA)                                     **
!     ** WITH EPSX_UNIF=-(3/(4PI))*(3PI**2*N)**(1/3)                          **
!     ** INPUT:  N=DENSITY, GN2=(GRAD N)**2, TAU=KINETIC ENERGY DENSITY       **
!     ** OUTPUT: AX AND ITS DERIVATIVES WITH RESPECT TO N, GN2 AND TAU        **
!     **************************************************************************
      USE SCAN_MODULE, ONLY : C23,CEX
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: GN2
      REAL(8)   ,INTENT(IN) :: TAU
      REAL(8)   ,INTENT(OUT):: AX
      REAL(8)   ,INTENT(OUT):: DAX_DN
      REAL(8)   ,INTENT(OUT):: DAX_DGN2
      REAL(8)   ,INTENT(OUT):: DAX_DTAU
      REAL(8)               :: P,DP_DN,DP_DGN2
      REAL(8)               :: TAUW,TAUUNIF
      REAL(8)               :: ALPHA,DA_DN,DA_DGN2,DA_DTAU
      REAL(8)               :: FX,DFX_DP,DFX_DA
      REAL(8)               :: EXUNIF,DEXUNIF_DN
!     **************************************************************************
!     == DIMENSIONLESS GRADIENT P=S**2=(GRAD N)**2/(4*(3PI**2)**(2/3)*N**(8/3))
      DP_DGN2=1.D0/(4.D0*C23*N**(8.D0/3.D0))
      P=GN2*DP_DGN2
      DP_DN=-8.D0*P/(3.D0*N)
!
!     == ALPHA=(TAU-TAUW)/TAUUNIF (DIRECT FORM, AVOIDS 0/0 OF Z=TAUW/TAU) ======
      TAUW=GN2/(8.D0*N)
      TAUUNIF=0.3D0*C23*N**(5.D0/3.D0)
      ALPHA=(TAU-TAUW)/TAUUNIF
      DA_DN=(TAUW/N)/TAUUNIF-ALPHA*(5.D0/3.D0)/N
      DA_DGN2=-1.D0/(8.D0*N*TAUUNIF)
      DA_DTAU=1.D0/TAUUNIF
!     == CLAMP ALPHA TO THE PHYSICAL RANGE ALPHA>=0. (NUMERICALLY TAU<TAUW    ==
!     == CAN OCCUR.) THE DERIVATIVES OF ALPHA ARE ZEROED CONSISTENTLY.        ==
      IF(ALPHA.LT.0.D0) THEN
        ALPHA=0.D0
        DA_DN=0.D0
        DA_DGN2=0.D0
        DA_DTAU=0.D0
      END IF
!
!     == ENHANCEMENT FACTOR AND ASSEMBLY =======================================
      CALL SCAN_FX(P,ALPHA,FX,DFX_DP,DFX_DA)
      EXUNIF=-CEX*N**(1.D0/3.D0)
      DEXUNIF_DN=EXUNIF/(3.D0*N)
      AX=N*EXUNIF*FX
      DAX_DN=(EXUNIF+N*DEXUNIF_DN)*FX &
     &      +N*EXUNIF*(DFX_DP*DP_DN+DFX_DA*DA_DN)
      DAX_DGN2=N*EXUNIF*(DFX_DP*DP_DGN2+DFX_DA*DA_DGN2)
      DAX_DTAU=N*EXUNIF*DFX_DA*DA_DTAU
      RETURN
      END SUBROUTINE SCAN_XUNP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN_CORR(RUP,RDN,GNT2,TAU,AC,DAC_DRUP,DAC_DRDN &
     &                    ,DAC_DGNT2,DAC_DTAU)
!     **************************************************************************
!     ** SCAN CORRELATION ENERGY PER VOLUME                                   **
!     **    AC=N*EPSC WITH EPSC=EPSC1+FC(ALPHA)*(EPSC0-EPSC1)                 **
!     ** FROM THE SUPPLEMENTAL MATERIAL OF PRL 115, 036402 (2015).            **
!     **   EPSC1=EPSC_LSDA1(PW92)+H1  (ALPHA APPROX 1 LIMIT)                  **
!     **   EPSC0=(EPSC_LDA0+H0)*GC(ZETA) (ALPHA=0 LIMIT)                      **
!     ** INPUT:  RUP,RDN = SPIN DENSITIES                                     **
!     **         GNT2    = (GRAD (RUP+RDN))**2                                **
!     **         TAU     = TOTAL KINETIC ENERGY DENSITY                       **
!     ** OUTPUT: AC AND ITS DERIVATIVES                                       **
!     **************************************************************************
      USE SCAN_MODULE, ONLY : C23,CRS,CT2,GAMMAC,B1C,B2C,B3C,C1C,C2C,DC,CHIINF,BETA0,GCPREF
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RUP
      REAL(8)   ,INTENT(IN) :: RDN
      REAL(8)   ,INTENT(IN) :: GNT2
      REAL(8)   ,INTENT(IN) :: TAU
      REAL(8)   ,INTENT(OUT):: AC
      REAL(8)   ,INTENT(OUT):: DAC_DRUP
      REAL(8)   ,INTENT(OUT):: DAC_DRDN
      REAL(8)   ,INTENT(OUT):: DAC_DGNT2
      REAL(8)   ,INTENT(OUT):: DAC_DTAU
      REAL(8)   ,PARAMETER  :: ZETAMAX=1.D0-1.D-10
      REAL(8)               :: N,ZETA,DZ_DRUP,DZ_DRDN
      REAL(8)               :: RS,DRS_DN
      REAL(8)               :: P,DP_DN,DP_DG
      REAL(8)               :: OPZ,OMZ
      REAL(8)               :: DS,DDS_DZ                 ! DS(ZETA)
      REAL(8)               :: PHI,DPHI_DZ
      REAL(8)               :: TAUW,TAUUNIF
      REAL(8)               :: ALPHA,DA_DN,DA_DZ,DA_DG,DA_DTAU
      REAL(8)               :: ECLSDA,DECLSDA_DRS,DECLSDA_DZ
      REAL(8)               :: PHI3,W1,DW1_DRS,DW1_DZ
      REAL(8)               :: BETA,DBETA_DRS
      REAL(8)               :: T2,DT2_DRS,DT2_DZ,DT2_DP
      REAL(8)               :: AA,DAA_DRS,DAA_DZ         ! A=BETA/(GAMMAC*W1)
      REAL(8)               :: Y,DY_DRS,DY_DZ,DY_DP      ! Y=A*T**2
      REAL(8)               :: G,DG_DY                   ! G=(1+4Y)**(-1/4)
      REAL(8)               :: DD                        ! 1+W1*(1-G)
      REAL(8)               :: H1,DH1_DRS,DH1_DZ,DH1_DP
      REAL(8)               :: EC1,DEC1_DRS,DEC1_DZ,DEC1_DP
      REAL(8)               :: SQRS,QLDA0
      REAL(8)               :: ECLDA0,DECLDA0_DRS
      REAL(8)               :: W0,DW0_DRS
      REAL(8)               :: GINF,DGINF_DP
      REAL(8)               :: DD0                       ! 1+W0*(1-GINF)
      REAL(8)               :: H0,DH0_DRS,DH0_DP
      REAL(8)               :: DXZ,DDXZ_DZ,GC,DGC_DZ
      REAL(8)               :: EC0,DEC0_DRS,DEC0_DZ,DEC0_DP
      REAL(8)               :: FC,DFC_DA
      REAL(8)               :: EC,DEC_DRS,DEC_DZ,DEC_DP,DEC_DA
      REAL(8)               :: DEC_DN,DEC_DZTOT
!     **************************************************************************
      N=RUP+RDN
!
!     ==========================================================================
!     == BASIC VARIABLES: ZETA, RS, P, DS(ZETA), PHI(ZETA), ALPHA             ==
!     ==========================================================================
      ZETA=(RUP-RDN)/N
      IF(ABS(ZETA).GE.ZETAMAX) THEN
        ZETA=MAX(-ZETAMAX,MIN(ZETAMAX,ZETA))
        DZ_DRUP=0.D0
        DZ_DRDN=0.D0
      ELSE
        DZ_DRUP=(1.D0-ZETA)/N
        DZ_DRDN=-(1.D0+ZETA)/N
      END IF
      OPZ=1.D0+ZETA
      OMZ=1.D0-ZETA
!
      RS=CRS/N**(1.D0/3.D0)
      DRS_DN=-RS/(3.D0*N)
!
      DP_DG=1.D0/(4.D0*C23*N**(8.D0/3.D0))
      P=GNT2*DP_DG
      DP_DN=-8.D0*P/(3.D0*N)
!
      DS=(OPZ**(5.D0/3.D0)+OMZ**(5.D0/3.D0))/2.D0
      DDS_DZ=(5.D0/6.D0)*(OPZ**(2.D0/3.D0)-OMZ**(2.D0/3.D0))
      PHI=(OPZ**(2.D0/3.D0)+OMZ**(2.D0/3.D0))/2.D0
      DPHI_DZ=(OPZ**(-1.D0/3.D0)-OMZ**(-1.D0/3.D0))/3.D0
!
!     == ALPHA=(TAU-TAUW)/(TAUUNIF*DS(ZETA)) ===================================
      TAUW=GNT2/(8.D0*N)
      TAUUNIF=0.3D0*C23*N**(5.D0/3.D0)*DS
      ALPHA=(TAU-TAUW)/TAUUNIF
      DA_DN=(TAUW/N)/TAUUNIF-ALPHA*(5.D0/3.D0)/N   ! AT FIXED ZETA,GNT2,TAU
      DA_DZ=-ALPHA*DDS_DZ/DS
      DA_DG=-1.D0/(8.D0*N*TAUUNIF)
      DA_DTAU=1.D0/TAUUNIF
      IF(ALPHA.LT.0.D0) THEN
        ALPHA=0.D0
        DA_DN=0.D0
        DA_DZ=0.D0
        DA_DG=0.D0
        DA_DTAU=0.D0
      END IF
!
!     ==========================================================================
!     == EPSC1=EPSC_LSDA1(PW92)+H1  (SEE SUPPL.MAT., PBE-LIKE WITH            ==
!     == G(A*T**2)=(1+4*A*T**2)**(-1/4) INSTEAD OF THE PBE RATIONAL FORM)     ==
!     ==========================================================================
      CALL SCAN_PW92(RS,ZETA,ECLSDA,DECLSDA_DRS,DECLSDA_DZ)
      PHI3=PHI**3
      W1=EXP(-ECLSDA/(GAMMAC*PHI3))-1.D0
      DW1_DRS=-(1.D0+W1)*DECLSDA_DRS/(GAMMAC*PHI3)
      DW1_DZ=(1.D0+W1)*(3.D0*ECLSDA*DPHI_DZ/(GAMMAC*PHI3*PHI) &
     &                 -DECLSDA_DZ/(GAMMAC*PHI3))
!
      BETA=BETA0*(1.D0+0.1D0*RS)/(1.D0+0.1778D0*RS)
      DBETA_DRS=BETA0*(0.1D0-0.1778D0)/(1.D0+0.1778D0*RS)**2
!
!     == T**2=(3PI**2/16)**(2/3)*P/(PHI**2*RS) =================================
      DT2_DP=CT2/(PHI**2*RS)
      T2=P*DT2_DP
      DT2_DRS=-T2/RS
      DT2_DZ=-2.D0*T2*DPHI_DZ/PHI
!
      AA=BETA/(GAMMAC*W1)
      DAA_DRS=DBETA_DRS/(GAMMAC*W1)-AA*DW1_DRS/W1
      DAA_DZ=-AA*DW1_DZ/W1
!
      Y=AA*T2
      DY_DRS=DAA_DRS*T2+AA*DT2_DRS
      DY_DZ=DAA_DZ*T2+AA*DT2_DZ
      DY_DP=AA*DT2_DP
      G=(1.D0+4.D0*Y)**(-0.25D0)
      DG_DY=-(1.D0+4.D0*Y)**(-1.25D0)
!
      DD=1.D0+W1*(1.D0-G)
      H1=GAMMAC*PHI3*LOG(DD)
      DH1_DRS=GAMMAC*PHI3*(DW1_DRS*(1.D0-G)-W1*DG_DY*DY_DRS)/DD
      DH1_DZ=GAMMAC*PHI3*(DW1_DZ*(1.D0-G)-W1*DG_DY*DY_DZ)/DD &
     &      +3.D0*GAMMAC*PHI**2*DPHI_DZ*LOG(DD)
      DH1_DP=-GAMMAC*PHI3*W1*DG_DY*DY_DP/DD
!
      EC1=ECLSDA+H1
      DEC1_DRS=DECLSDA_DRS+DH1_DRS
      DEC1_DZ=DECLSDA_DZ+DH1_DZ
      DEC1_DP=DH1_DP
!
!     ==========================================================================
!     == EPSC0=(EPSC_LDA0+H0)*GC(ZETA)  (ALPHA=0 LIMIT)                       ==
!     == G_INFINITY IS EVALUATED AT ZETA=0 (SEE SUPPL.MAT.): CHIINF=0.128026  ==
!     ==========================================================================
      SQRS=SQRT(RS)
      QLDA0=1.D0+B2C*SQRS+B3C*RS
      ECLDA0=-B1C/QLDA0
      DECLDA0_DRS=B1C*(0.5D0*B2C/SQRS+B3C)/QLDA0**2
!
      W0=EXP(-ECLDA0/B1C)-1.D0
      DW0_DRS=-(1.D0+W0)*DECLDA0_DRS/B1C
!
      GINF=(1.D0+4.D0*CHIINF*P)**(-0.25D0)
      DGINF_DP=-CHIINF*(1.D0+4.D0*CHIINF*P)**(-1.25D0)
!
      DD0=1.D0+W0*(1.D0-GINF)
      H0=B1C*LOG(DD0)
      DH0_DRS=B1C*DW0_DRS*(1.D0-GINF)/DD0
      DH0_DP=-B1C*W0*DGINF_DP/DD0
!
      DXZ=(OPZ**(4.D0/3.D0)+OMZ**(4.D0/3.D0))/2.D0
      DDXZ_DZ=(2.D0/3.D0)*(OPZ**(1.D0/3.D0)-OMZ**(1.D0/3.D0))
      GC=(1.D0-GCPREF*(DXZ-1.D0))*(1.D0-ZETA**12)
      DGC_DZ=-GCPREF*DDXZ_DZ*(1.D0-ZETA**12) &
     &      -12.D0*ZETA**11*(1.D0-GCPREF*(DXZ-1.D0))
!
      EC0=(ECLDA0+H0)*GC
      DEC0_DRS=(DECLDA0_DRS+DH0_DRS)*GC
      DEC0_DZ=(ECLDA0+H0)*DGC_DZ
      DEC0_DP=DH0_DP*GC
!
!     ==========================================================================
!     == INTERPOLATION: EPSC=EPSC1+FC(ALPHA)*(EPSC0-EPSC1)                    ==
!     ==========================================================================
      CALL SCAN_INTERPOLATE(ALPHA,C1C,C2C,DC,FC,DFC_DA)
      EC=EC1+FC*(EC0-EC1)
      DEC_DRS=DEC1_DRS+FC*(DEC0_DRS-DEC1_DRS)
      DEC_DZ=DEC1_DZ+FC*(DEC0_DZ-DEC1_DZ)
      DEC_DP=DEC1_DP+FC*(DEC0_DP-DEC1_DP)
      DEC_DA=DFC_DA*(EC0-EC1)
!
!     ==========================================================================
!     == COLLECT TOTAL DERIVATIVES OF AC=N*EPSC                               ==
!     ==========================================================================
      DEC_DN=DEC_DRS*DRS_DN+DEC_DP*DP_DN+DEC_DA*DA_DN  ! AT FIXED ZETA,GNT2,TAU
      DEC_DZTOT=DEC_DZ+DEC_DA*DA_DZ
      AC=N*EC
      DAC_DRUP=EC+N*(DEC_DN+DEC_DZTOT*DZ_DRUP)
      DAC_DRDN=EC+N*(DEC_DN+DEC_DZTOT*DZ_DRDN)
      DAC_DGNT2=N*(DEC_DP*DP_DG+DEC_DA*DA_DG)
      DAC_DTAU=N*DEC_DA*DA_DTAU
      RETURN
      END SUBROUTINE SCAN_CORR
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN$EVAL1(VAL,EXC,DEXC)
!     **************************************************************************
!     ** SCAN EXCHANGE-CORRELATION ENERGY PER VOLUME EXC=RHOT*EPSILON_XC      **
!     ** AND ITS ANALYTIC FIRST DERIVATIVES WITH RESPECT TO THE META-GGA      **
!     ** STATE VARIABLES                                                      **
!     **    VAL=(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST,LAPLT,LAPLS,TAUT,TAUS)        **
!     ** SCAN DOES NOT DEPEND ON THE LAPLACIANS: DEXC(6)=DEXC(7)=0.           **
!     **************************************************************************
      USE SCAN_MODULE, ONLY : TINI,RHOMIN
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL(9)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: DEXC(9)
      REAL(8)               :: RUP,RDN,SUU,SDD,TUP,TDN
      REAL(8)               :: RHOS
      REAL(8)               :: AX,XU_N,XU_S,XU_T,XD_N,XD_S,XD_T
      REAL(8)               :: AC,C_UP,C_DN,C_G,C_T
      REAL(8)               :: SVAR
!     **************************************************************************
      IF(.NOT.TINI) CALL SCAN_INITIALIZE()
      EXC=0.D0
      DEXC(:)=0.D0
      IF(VAL(1).LT.RHOMIN) RETURN
!
!     ==========================================================================
!     == TRANSFORM TO SPIN-CHANNEL QUANTITIES                                 ==
!     ==   RUP/RDN: SPIN DENSITIES;  SUU/SDD: (GRAD RUP)**2,(GRAD RDN)**2     ==
!     ==   TUP/TDN: SPIN KINETIC ENERGY DENSITIES                             ==
!     ==========================================================================
      RHOS=VAL(2)
      RHOS=MAX(-VAL(1)+RHOMIN,MIN(VAL(1)-RHOMIN,RHOS))
      RUP=0.5D0*(VAL(1)+RHOS)
      RDN=0.5D0*(VAL(1)-RHOS)
      SUU=MAX(0.D0,0.25D0*(VAL(3)+VAL(4)+2.D0*VAL(5)))
      SDD=MAX(0.D0,0.25D0*(VAL(3)+VAL(4)-2.D0*VAL(5)))
      TUP=MAX(0.D0,0.5D0*(VAL(8)+VAL(9)))
      TDN=MAX(0.D0,0.5D0*(VAL(8)-VAL(9)))
!
!     ==========================================================================
!     == EXCHANGE VIA SPIN SCALING: EX[NUP,NDN]=(EX[2*NUP]+EX[2*NDN])/2       ==
!     == CHAIN RULE: D(0.5*AX[2*RUP,4*SUU,2*TUP])/D(RUP)=DAX_DN ETC.          ==
!     ==========================================================================
      CALL SCAN_XUNP(2.D0*RUP,4.D0*SUU,2.D0*TUP,AX,XU_N,XU_S,XU_T)
      EXC=EXC+0.5D0*AX
      XU_S=2.D0*XU_S
      CALL SCAN_XUNP(2.D0*RDN,4.D0*SDD,2.D0*TDN,AX,XD_N,XD_S,XD_T)
      EXC=EXC+0.5D0*AX
      XD_S=2.D0*XD_S
!
!     ==========================================================================
!     == CORRELATION (DEPENDS ON THE TOTAL GRADIENT AND TOTAL TAU ONLY)       ==
!     ==========================================================================
      CALL SCAN_CORR(RUP,RDN,MAX(0.D0,VAL(3)),TUP+TDN,AC,C_UP,C_DN,C_G,C_T)
      EXC=EXC+AC
!     Differentiate the actual clipped energy, including numerical PAW tails.
      XU_N=XU_N+C_UP
      XD_N=XD_N+C_DN
      XU_T=XU_T+C_T
      XD_T=XD_T+C_T
      IF(VAL(8)+VAL(9).LT.0.D0) XU_T=0.D0
      IF(VAL(8)-VAL(9).LT.0.D0) XD_T=0.D0
      IF(VAL(3)+VAL(4)+2.D0*VAL(5).LT.0.D0) XU_S=0.D0
      IF(VAL(3)+VAL(4)-2.D0*VAL(5).LT.0.D0) XD_S=0.D0
      IF(VAL(3).LT.0.D0) C_G=0.D0
!
!     ==========================================================================
!     == CHAIN RULE BACK TO THE CP-PAW STATE VARIABLES                        ==
!     ==   RUP=(RHOT+RHOS)/2, SUU=(GRHOT2+GRHOS2+2*GRHOST)/4, ETC.            ==
!     ==========================================================================
      DEXC(1)=0.5D0*(XU_N+XD_N)
      DEXC(2)=0.5D0*(XU_N-XD_N)
      IF(VAL(2).GT.VAL(1)-RHOMIN) THEN
        DEXC(1)=XU_N
        DEXC(2)=0.D0
      ELSE IF(VAL(2).LT.-VAL(1)+RHOMIN) THEN
        DEXC(1)=XD_N
        DEXC(2)=0.D0
      END IF
      SVAR=0.25D0*(XU_S+XD_S)
      DEXC(3)=SVAR+C_G
      DEXC(4)=SVAR
      DEXC(5)=0.5D0*(XU_S-XD_S)
      DEXC(6)=0.D0
      DEXC(7)=0.D0
      DEXC(8)=0.5D0*(XU_T+XD_T)
      DEXC(9)=0.5D0*(XU_T-XD_T)
      RETURN
      END SUBROUTINE SCAN$EVAL1
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN$EVAL2(VAL,EXC,DEXC,D2EXC)
!     **************************************************************************
!     ** SCAN ENERGY WITH FIRST AND SECOND DERIVATIVES.                       **
!     ** THE SECOND DERIVATIVES ARE OBTAINED BY CENTRAL FINITE DIFFERENCES    **
!     ** OF THE ANALYTIC FIRST DERIVATIVES.                                   **
!     ** CAUTION: NEAR THE NON-ANALYTIC POINTS OF SCAN (ALPHA=1, ALPHA=0      **
!     ** CLAMP, ZETA=+-1) THE NUMERICAL DERIVATIVES ARE APPROXIMATE.          **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL(9)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: DEXC(9)
      REAL(8)   ,INTENT(OUT):: D2EXC(9,9)
      REAL(8)   ,PARAMETER  :: STEP=1.D-5
      REAL(8)               :: VAL1(9),H(9)
      REAL(8)               :: EP,EM,DP(9),DM(9)
      INTEGER(4)            :: I,J
!     **************************************************************************
      CALL SCAN$EVAL1(VAL,EXC,DEXC)
      DO I=1,9
        H(I)=STEP*MAX(ABS(VAL(I)),1.D-2)
      ENDDO
      D2EXC(:,:)=0.D0
      DO I=1,9
        VAL1=VAL
        VAL1(I)=VAL(I)+H(I)
        CALL SCAN$EVAL1(VAL1,EP,DP)
        VAL1(I)=VAL(I)-H(I)
        CALL SCAN$EVAL1(VAL1,EM,DM)
        D2EXC(I,:)=(DP(:)-DM(:))/(2.D0*H(I))
      ENDDO
!     == SYMMETRIZE ============================================================
      DO I=1,9
        DO J=I+1,9
          D2EXC(I,J)=0.5D0*(D2EXC(I,J)+D2EXC(J,I))
          D2EXC(J,I)=D2EXC(I,J)
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE SCAN$EVAL2
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCAN$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
!     **************************************************************************
!     ** SCAN ENERGY WITH FIRST, SECOND AND THIRD DERIVATIVES.                **
!     ** SECOND AND THIRD DERIVATIVES ARE OBTAINED BY CENTRAL FINITE          **
!     ** DIFFERENCES OF THE ANALYTIC FIRST DERIVATIVES:                       **
!     **    D3EXC(I,J,K)=D2/(DVI*DVJ) [DEXC(K)]                               **
!     ** NOTE: THIS COSTS O(100) EVALUATIONS OF SCAN$EVAL1 PER CALL AND IS    **
!     ** INTENDED FOR THE ONE-CENTER (RADIAL) EXPANSIONS ONLY.                **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL(9)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: DEXC(9)
      REAL(8)   ,INTENT(OUT):: D2EXC(9,9)
      REAL(8)   ,INTENT(OUT):: D3EXC(9,9,9)
      REAL(8)   ,PARAMETER  :: STEP=1.D-4
      REAL(8)               :: VAL1(9),H(9)
      REAL(8)               :: E0,EP,EM,D0(9),DP(9),DM(9)
      REAL(8)               :: DPP(9),DPM(9),DMP(9),DMM(9)
      INTEGER(4)            :: I,J,K
!     **************************************************************************
      CALL SCAN$EVAL2(VAL,EXC,DEXC,D2EXC)
      D3EXC(:,:,:)=0.D0
      DO I=1,9
        H(I)=STEP*MAX(ABS(VAL(I)),1.D-2)
      ENDDO
      CALL SCAN$EVAL1(VAL,E0,D0)
!     == DIAGONAL BLOCKS I=J ===================================================
      DO I=1,9
        VAL1=VAL
        VAL1(I)=VAL(I)+H(I)
        CALL SCAN$EVAL1(VAL1,EP,DP)
        VAL1(I)=VAL(I)-H(I)
        CALL SCAN$EVAL1(VAL1,EM,DM)
        D3EXC(I,I,:)=(DP(:)-2.D0*D0(:)+DM(:))/H(I)**2
      ENDDO
!     == OFF-DIAGONAL BLOCKS I<J ===============================================
      DO I=1,9
        DO J=I+1,9
          VAL1=VAL
          VAL1(I)=VAL(I)+H(I)
          VAL1(J)=VAL(J)+H(J)
          CALL SCAN$EVAL1(VAL1,EP,DPP)
          VAL1(J)=VAL(J)-H(J)
          CALL SCAN$EVAL1(VAL1,EP,DPM)
          VAL1(I)=VAL(I)-H(I)
          VAL1(J)=VAL(J)+H(J)
          CALL SCAN$EVAL1(VAL1,EM,DMP)
          VAL1(J)=VAL(J)-H(J)
          CALL SCAN$EVAL1(VAL1,EM,DMM)
          D3EXC(I,J,:)=(DPP(:)-DPM(:)-DMP(:)+DMM(:))/(4.D0*H(I)*H(J))
          D3EXC(J,I,:)=D3EXC(I,J,:)
        ENDDO
      ENDDO
!     == SYMMETRIZE OVER ALL INDEX PERMUTATIONS ================================
      DO I=1,9
        DO J=I,9
          DO K=J,9
            D3EXC(I,J,K)=(D3EXC(I,J,K)+D3EXC(I,K,J)+D3EXC(J,I,K) &
     &                   +D3EXC(J,K,I)+D3EXC(K,I,J)+D3EXC(K,J,I))/6.D0
            D3EXC(I,K,J)=D3EXC(I,J,K)
            D3EXC(J,I,K)=D3EXC(I,J,K)
            D3EXC(J,K,I)=D3EXC(I,J,K)
            D3EXC(K,I,J)=D3EXC(I,J,K)
            D3EXC(K,J,I)=D3EXC(I,J,K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE SCAN$EVAL3

      SUBROUTINE SCAN$BLOCK(NPOINT,RHO,GRAD,TAU,WEIGHT,ENERGY,DRHO,DGRAD,DTAU,DWEIGHT)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NPOINT
      REAL(8),INTENT(IN) :: RHO(NPOINT,2),GRAD(NPOINT,3,2),TAU(NPOINT,2),WEIGHT(NPOINT)
      REAL(8),INTENT(OUT) :: ENERGY,DRHO(NPOINT,2),DGRAD(NPOINT,3,2),DTAU(NPOINT,2),DWEIGHT(NPOINT)
      REAL(8) :: V(9),E,D(9),GT(3),GS(3),SIGN
      INTEGER :: I,SPIN
      IF(ANY(.NOT.IEEE_IS_FINITE(RHO)).OR.ANY(.NOT.IEEE_IS_FINITE(GRAD)) &
     &   .OR.ANY(.NOT.IEEE_IS_FINITE(TAU)).OR.ANY(.NOT.IEEE_IS_FINITE(WEIGHT))) &
     &   ERROR STOP 'Nonfinite SCAN quadrature input'
      ENERGY=0.D0
      DO I=1,NPOINT
        GT=GRAD(I,:,1)+GRAD(I,:,2)
        GS=GRAD(I,:,1)-GRAD(I,:,2)
        V=[SUM(RHO(I,:)),RHO(I,1)-RHO(I,2),SUM(GT**2),SUM(GS**2),DOT_PRODUCT(GT,GS) &
     &     ,0.D0,0.D0,SUM(TAU(I,:)),TAU(I,1)-TAU(I,2)]
        CALL SCAN$EVAL1(V,E,D)
        IF(.NOT.IEEE_IS_FINITE(E).OR.ANY(.NOT.IEEE_IS_FINITE(D))) &
     &    ERROR STOP 'Nonfinite SCAN quadrature result'
        ENERGY=ENERGY+WEIGHT(I)*E
        DWEIGHT(I)=E
        DO SPIN=1,2
          SIGN=REAL(3-2*SPIN,8)
          DRHO(I,SPIN)=WEIGHT(I)*(D(1)+SIGN*D(2))
          DGRAD(I,:,SPIN)=WEIGHT(I)*(2.D0*D(3)*GT+2.D0*SIGN*D(4)*GS &
     &                                          +D(5)*(GS+SIGN*GT))
          DTAU(I,SPIN)=WEIGHT(I)*(D(8)+SIGN*D(9))
        END DO
      END DO
      END SUBROUTINE SCAN$BLOCK
