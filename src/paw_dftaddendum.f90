!*******************************************************************************
!*******************************************************************************
!****  TEMPLATE FOR NEW DENSITY FUNCTIONALS IN THE PAW_DFT OBJECT           ****
!****                                                                       ****
!****  THE PAW METHOD REQUIRES TWO CALLS EVAL1 AND EVAL3                    ****
!****  - EVAL1 PROVIDES THE ENERGY AND ITS XC POTENTIAL FOR THE PLANE WAVE  ****
!****          PART                                                         ****
!****  - EVAL3 PROVIDES THE ENERGY AND ITS FIRST THREE DERIVATIVES WITH     ****
!****          RESPECT TO THE ENTERING VARIABLES FOR THE ONE-CENTER EXPANS. ****
!****          THE NON-SPHERICAL PARTS ARE INCLUDED BY A SECOND ORDER       ****
!****          TAYLOR EXPANSION IN THE NON-SPHERICAL CONTRIBUTIONS ABOUT    ****
!****          THE SPHERICAL PARTS                                          ****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
MODULE NEWDFT_MODULE
LOGICAL(4) :: TINI=.FALSE.
END MODULE NEWDFT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWDFT_INITIALIZE()
      USE NEWDFT_MODULE, ONLY : TINI
      IMPLICIT NONE
      IF(TINI) RETURN
      TINI=.TRUE.
      END SUBROUTINE NEWDFT_INITIALIZE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWDFT$EVAL1(VAL,EXC,DEXC)
!     **************************************************************************
!     **  BRIEF DESCRIPTION                                           **
!     **************************************************************************
      USE NEWDFT_MODULE, ONLY : TINI
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: VAL(5)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)
!     ******************************************************************
!EXC=RHOT*EPSILON_XC
!VAL(1)=TOTAL DENSITY=RHOT
!VAL(2)=SPINDENSITY=RHOS= (NUP-NDOWN)
!VAL(3)= (GRAD*RHOT)**2
!VAL(4)= (GRAD*RHOS)**2
!VAL(3)= (GRAD*RHOT)*(GRAD*RHOS)
      CALL NEWDFT_INITIALIZE()
      EXC=0.D0
      DEXC(:)=0.D0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWDFT$EVAL2(VAL,EXC,DEXC,D2EXC)
!     **************************************************************************
!     **  BRIEF DESCRIPTION                                           **
!     **************************************************************************
      USE NEWDFT_MODULE, ONLY : TINI
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: VAL(5)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)
      REAL(8),INTENT(OUT) :: D2EXC(5,5)
!     ******************************************************************
      CALL NEWDFT_INITIALIZE()
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEWDFT$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
!     **************************************************************************
!     **  BRIEF DESCRIPTION                                           **
!     **************************************************************************
      USE NEWDFT_MODULE
      REAL(8),INTENT(IN)  :: VAL(5)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)
      REAL(8),INTENT(OUT) :: D2EXC(5,5)
      REAL(8),INTENT(OUT) :: D3EXC(5,5,5)
!     ******************************************************************
      CALL NEWDFT_INITIALIZE()
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      D3EXC(:,:,:)=0.D0
      RETURN
      END
!*******************************************************************************
!*******************************************************************************
!*****  INTERFACE TO LIBXC                                                ******
!*******************************************************************************
!*******************************************************************************
#IF DEFINED(CPPVAR_NOLIBXC)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GETL4(ID,VAL)
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
      CALL ERROR$MSG('LIBXC NOT AVAILABLE. INSTALL WITH LIBXC')
      CALL ERROR$STOP('PAWLIBXC$GETL4')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$REPORT(NFIL)
      INTEGER(4)       ,INTENT(IN) :: NFIL
      CALL ERROR$MSG('LIBXC NOT AVAILABLE. INSTALL WITH LIBXC')
      CALL ERROR$STOP('PAWLIBXC$REPORT')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$SETCH(ID,VAL)
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
      CALL ERROR$MSG('LIBXC NOT AVAILABLE. INSTALL WITH LIBXC')
      CALL ERROR$STOP('PAWLIBXC$SETCH')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$SETCHA(ID,LEN,VAL)
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
      CALL ERROR$MSG('LIBXC NOT AVAILABLE. INSTALL WITH LIBXC')
      CALL ERROR$STOP('PAWLIBXC$SETCHA')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GGA1(VAL,EXC,DER)
      REAL(8)   ,INTENT(IN) :: VAL(5)     ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC         ! SPIN DENSITY
      REAL(8)   ,INTENT(OUT):: DER(5)      ! 
      CALL ERROR$MSG('LIBXC NOT AVAILABLE. INSTALL WITH LIBXC')
      CALL ERROR$STOP('PAWLIBXC$GGA1')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GGA2(VAL,EXC,DER,DER2)
      REAL(8)   ,INTENT(IN) :: VAL(5)     ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC         ! SPIN DENSITY
      REAL(8)   ,INTENT(OUT):: DER(5)      ! 
      REAL(8)   ,INTENT(OUT):: DER2(5,5)   ! 
      CALL ERROR$MSG('LIBXC NOT AVAILABLE. INSTALL WITH LIBXC')
      CALL ERROR$STOP('PAWLIBXC$GGA2')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GGA3(VAL,EXC,DER,DER2,DER3)
      REAL(8)   ,INTENT(IN) :: VAL(5)     ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC         ! SPIN DENSITY
      REAL(8)   ,INTENT(OUT):: DER(5)      ! 
      REAL(8)   ,INTENT(OUT):: DER2(5,5)   ! 
      REAL(8)   ,INTENT(OUT):: DER3(5,5,5) ! 
      CALL ERROR$MSG('LIBXC NOT AVAILABLE. INSTALL WITH LIBXC')
      CALL ERROR$STOP('PAWLIBXC$GGA3')
      RETURN
      END
#ELSE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE PAWLIBXC_MODULE
!     **************************************************************************
!     ** COMMON INFORMATION FOR INTERFACING CP-PAW WITH LIBXC                 **
!     **                                                                      **
!     **   RHO(1)  =RHOUP    =0.5*(RHOT+RHOS)
!     **   RHO(2)  =RHODN    =0.5*(RHOT-RHOS)
!     **   SIGMA(1)=GRHO2UPUP=0.25*(GRHOT2+GRHOS2+2*GRHOST)
!     **   SIGMA(2)=GRHO2UPDN=0.25*(GRHOT2-GRHOS2)
!     **   SIGMA(3)=GRHO2DNDN=0.25*(GRHOT2+GRHOS2-2*GRHOST)
!     **   LAPL(1) =LAPLUP    =0.5*(LAPLT+LAPLS)
!     **   LAPL(2) =LAPLDN    =0.5*(LAPLT-LAPLS)
!     **   TAU(1)  =TAUUP     =0.5*(TAUT+TAUS)
!     **   TAU(2)  =TAUDN     =0.5*(TAUT-TAUS)
!     **                                                                      **
!     **  HTTPS://WWW.TDDFT.ORG/PROGRAMS/LIBXC/MANUAL/LIBXC-5.1.X/            **
!     **************************************************************************
!     **  FIND THE INTERFACES TO LIBXC SUBROUTINES IN FILE                    **
!     **     LIBXC/LIBXC-MASTER/SRC/LIBXC_MASTER.F90                          **
!     **************************************************************************
!     **  FOR THE FUNCTIONAL IDENTIFIERS OF LIBXC SEE                         **
!     **  MARQUES, OLIVIERA, BURNUS, ARXIV:1203.1739                          **
!     **************************************************************************
!     **  (RHO,SIGMA,LAPL,TAU)=MATINV* (RHOT,RHOS,GRHOT2,GRHOS2,GRHOST,LAPLT2..)
!     **   RHO:    MATINV=0.25 * ( 2  2 )      MAT=0.25 * ( 4  4 )
!     **                         ( 2 -2 )                 ( 4 -4 )
!     **   SIGMA:  MATINV=0.25 * ( 1  1  2 )   MAT=0.25 * (  1  2  1 ) 
!     **                         ( 1 -1  0 )              (  1 -2  1 )
!     **                         ( 1  1 -2 )              (  1  0 -1 )
!     **************************************************************************
      USE XC_F03_LIB_M, ONLY : XC_F03_FUNC_T
      IMPLICIT NONE
      LOGICAL(4) :: TINI=.FALSE.
!
!     ==========================================================================
!     == DEFINITION OF XC FUNCTIONAL                                          ==
!     == SOME FUNCTIONALS ARE COMPOSED OF SEVERAL PARTS TO BE ADDED UP        ==
!     ==========================================================================
      TYPE(XC_F03_FUNC_T),ALLOCATABLE :: XC_FUNC(:)
      INTEGER(4)                      :: NXC=0
!
!     ==========================================================================
!     == AVAILABLE IDS FOR DENSITY FUNCTIONALS IN LIBXC                       ==
!     ==========================================================================
      INTEGER(4),PARAMETER :: NIDX=734
      CHARACTER(32)        :: XCID(NIDX)
!
!     ==========================================================================
!     == TRANSFORMATION MATRIX FROM THE PARAMETERS USED IN CPPAW AND LIBXC
!     ==========================================================================
!     == THE TRANSPOSITION IS REQUIRED BECAUSE THE RUNNING (FIRST) INDEX IS
!     == TYPED HORIZONTALLY RATHER THAN VERTICALLY
!     ___LDA: (RHO)_____________________________________________________________
      REAL(8)    ,PARAMETER :: MAT_LDA(2,2)=0.25D0*REAL(TRANSPOSE(RESHAPE( &
     &                                   [+2,+2 &
     &                                   ,+2,-2] &
     &                                   ,[2,2])),KIND=8)
      REAL(8)    ,PARAMETER :: MATINV_LDA(2,2)=0.25D0*REAL(TRANSPOSE(RESHAPE( &
     &                                   [+4,+4 &
     &                                   ,+4,-4] &
     &                                   ,[2,2])),KIND=8)
!     ___GGA: (RHO,SIGMA)_______________________________________________________
      REAL(8)    ,PARAMETER :: MAT_GGA(5,5)=0.25D0*REAL(TRANSPOSE(RESHAPE( &
     &                                   [+2,+2,+0,+0,+0 &
     &                                   ,+2,-2,+0,+0,+0 &
     &                                   ,+0,+0,+1,+1,+2 &
     &                                   ,+0,+0,+1,-1,+0 &
     &                                   ,+0,+0,+1,+1,-2] &
     &                                   ,[5,5])),KIND=8)
      REAL(8)    ,PARAMETER :: MATINV_GGA(5,5)=0.25D0*REAL(TRANSPOSE(RESHAPE( &
     &                                   [+4,+4,+0,+0,+0 &
     &                                   ,+4,-4,+0,+0,+0 &
     &                                   ,+0,+0,+1,+2,+1 &
     &                                   ,+0,+0,+1,-2,+1 &
     &                                   ,+0,+0,+1,+0,-1] &
     &                                   ,[5,5])),KIND=8)
!     ___META GGA: (RHO,SIGMA,LAPL,TAU)_________________________________________
      REAL(8)    ,PARAMETER :: MAT_MGGA(9,9)=0.25D0*REAL(TRANSPOSE(RESHAPE( &
     &                                   [+2,+2,+0,+0,+0,+0,+0,+0,+0 &
     &                                   ,+2,-2,+0,+0,+0,+0,+0,+0,+0 &
     &                                   ,+0,+0,+1,+1,+2,+0,+0,+0,+0 &
     &                                   ,+0,+0,+1,-1,+0,+0,+0,+0,+0 &
     &                                   ,+0,+0,+1,+1,-2,+0,+0,+0,+0 &
     &                                   ,+0,+0,+0,+0,+0,+2,+2,+0,+0 &
     &                                   ,+0,+0,+0,+0,+0,+2,-2,+0,+0 &
     &                                   ,+0,+0,+0,+0,+0,+0,+0,+2,+2 &
     &                                   ,+0,+0,+0,+0,+0,+0,+0,+2,-2] &
     &                                   ,[9,9])),KIND=8)
      REAL(8)    ,PARAMETER :: MATINV_MGGA(9,9)=0.25D0*REAL(TRANSPOSE(RESHAPE( &
     &                                   [+2,+2,+0,+0,+0,+0,+0,+0,+0 &
     &                                   ,+2,-2,+0,+0,+0,+0,+0,+0,+0 &
     &                                   ,+0,+0,+1,+1,+1,+0,+0,+0,+0 &
     &                                   ,+0,+0,+1,-2,+1,+0,+0,+0,+0 &
     &                                   ,+0,+0,+1,+0,-1,+0,+0,+0,+0 &
     &                                   ,+0,+0,+0,+0,+0,+4,+4,+0,+0 &
     &                                   ,+0,+0,+0,+0,+0,+4,-4,+0,+0 &
     &                                   ,+0,+0,+0,+0,+0,+0,+0,+4,+4 &
     &                                   ,+0,+0,+0,+0,+0,+0,+0,+4,-4] &
     &                                   ,[9,9])),KIND=8)
      END MODULE PAWLIBXC_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC_INITIALIZE()
!     **************************************************************************
!     **************************************************************************
      USE PAWLIBXC_MODULE, ONLY: TINI &
     &                          ,XC_FUNC
      IMPLICIT NONE
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!
!     ==========================================================================
!     == SANITY CHECK OF THE PAWLIBXC MODULE                                  ==
!     ==========================================================================
!      CALL PAWLIBXC_TEST1()
!
!     ==========================================================================
!     == CHECK WHETHER FUNCTIONAL IS SET                                      ==
!     ==========================================================================
      IF(.NOT.ALLOCATED(XC_FUNC)) THEN
        PRINT*,'ERROR: XC_FUNC IS NOT SET'
        STOP 'IN PAWLIBXC_INITIALIZE'
      END IF
      ! THIS STATEMENT IS USED IN THE TEST PHASE
      CALL PAWLIBXC$REPORT(6)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC_TEST1()
!     **************************************************************************
!     ** SANITY CHECK: WHETHER THE MATRICES ARE INVERSES OF EACH OTHER        **
!     **************************************************************************
      USE PAWLIBXC_MODULE
      IMPLICIT NONE
      REAL(8) :: MAT2(2,2),MAT5(5,5),MAT9(9,9)
      INTEGER(4) :: I
!     **************************************************************************
      MAT2=MATMUL(MAT_LDA,MATINV_LDA)
      DO I=1,2
        MAT2(I,I)=MAT2(I,I)-1.D0
      ENDDO
      PRINT*,'MAXVAL(MAT2) ',MAXVAL(MAT2)
      MAT5=MATMUL(MAT_GGA,MATINV_GGA)
      DO I=1,5
        MAT5(I,I)=MAT5(I,I)-1.D0
      ENDDO
      PRINT*,'MAXVAL(MAT5) ',MAXVAL(MAT5)
      MAT9=MATMUL(MAT_MGGA,MATINV_MGGA)
      DO I=1,9
        MAT9(I,I)=MAT9(I,I)-1.D0
      ENDDO
      PRINT*,'MAXVAL(MAT9) ',MAXVAL(MAT9)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC_IDLIST()
      USE PAWLIBXC_MODULE, ONLY : XCID
      IMPLICIT NONE
      LOGICAL(4),SAVE :: TLISTINI=.FALSE.
!     **************************************************************************
      IF(TLISTINI) RETURN
      TLISTINI=.TRUE.
      XCID(:)=' '
      XCID(  1)='XC_LDA_X'                              
      XCID(  2)='XC_LDA_C_WIGNER'
      XCID(  3)='XC_LDA_C_RPA'
      XCID(  4)='XC_LDA_C_HL'
      XCID(  5)='XC_LDA_C_GL'
      XCID(  6)='XC_LDA_C_XALPHA'
      XCID(  7)='XC_LDA_C_VWN'
      XCID(  8)='XC_LDA_C_VWN_RPA'
      XCID(  9)='XC_LDA_C_PZ'
      XCID( 10)='XC_LDA_C_PZ_MOD'
      XCID( 11)='XC_LDA_C_OB_PZ'
      XCID( 12)='XC_LDA_C_PW'
      XCID( 13)='XC_LDA_C_PW_MOD'
      XCID( 14)='XC_LDA_C_OB_PW'
      XCID( 15)='XC_LDA_C_2D_AMGB'
      XCID( 16)='XC_LDA_C_2D_PRM'
      XCID( 17)='XC_LDA_C_VBH'
      XCID( 18)='XC_LDA_C_1D_CSC'
      XCID( 19)='XC_LDA_X_2D'
      XCID( 20)='XC_LDA_XC_TETER93'
      XCID( 21)='XC_LDA_X_1D_SOFT'
      XCID( 22)='XC_LDA_C_ML1'
      XCID( 23)='XC_LDA_C_ML2'
      XCID( 24)='XC_LDA_C_GOMBAS'
      XCID( 25)='XC_LDA_C_PW_RPA'
      XCID( 26)='XC_LDA_C_1D_LOOS'
      XCID( 27)='XC_LDA_C_RC04'
      XCID( 28)='XC_LDA_C_VWN_1'
      XCID( 29)='XC_LDA_C_VWN_2'
      XCID( 30)='XC_LDA_C_VWN_3'
      XCID( 31)='XC_LDA_C_VWN_4'
      XCID( 32)='XC_GGA_X_GAM'
      XCID( 33)='XC_GGA_C_GAM'
      XCID( 34)='XC_GGA_X_HCTH_A'
      XCID( 35)='XC_GGA_X_EV93'
      XCID( 36)='XC_HYB_MGGA_X_DLDF'
      XCID( 37)='XC_MGGA_C_DLDF'
      XCID( 38)='XC_GGA_X_BCGP'
      XCID( 39)='XC_GGA_C_ACGGA'
      XCID( 40)='XC_GGA_X_LAMBDA_OC2_N'
      XCID( 41)='XC_GGA_X_B86_R'
      XCID( 42)='XC_MGGA_XC_ZLP'
      XCID( 43)='XC_LDA_XC_ZLP'
      XCID( 44)='XC_GGA_X_LAMBDA_CH_N'
      XCID( 45)='XC_GGA_X_LAMBDA_LO_N'
      XCID( 46)='XC_GGA_X_HJS_B88_V2'
      XCID( 47)='XC_GGA_C_Q2D'
      XCID( 48)='XC_GGA_X_Q2D'
      XCID( 49)='XC_GGA_X_PBE_MOL'
      XCID( 50)='XC_LDA_K_TF'
      XCID( 51)='XC_LDA_K_LP'
      XCID( 52)='XC_GGA_K_TFVW'
      XCID( 53)='XC_GGA_K_REVAPBEINT'
      XCID( 54)='XC_GGA_K_APBEINT'
      XCID( 55)='XC_GGA_K_REVAPBE'
      XCID( 56)='XC_GGA_X_AK13'
      XCID( 57)='XC_GGA_K_MEYER'
      XCID( 58)='XC_GGA_X_LV_RPW86'
      XCID( 59)='XC_GGA_X_PBE_TCA'
      XCID( 60)='XC_GGA_X_PBEINT'
      XCID( 61)='XC_GGA_C_ZPBEINT'
      XCID( 62)='XC_GGA_C_PBEINT'
      XCID( 63)='XC_GGA_C_ZPBESOL'
      XCID( 64)='XC_MGGA_XC_OTPSS_D'
      XCID( 65)='XC_GGA_XC_OPBE_D'
      XCID( 66)='XC_GGA_XC_OPWLYP_D'
      XCID( 67)='XC_GGA_XC_OBLYP_D'
      XCID( 68)='XC_GGA_X_VMT84_GE'
      XCID( 69)='XC_GGA_X_VMT84_PBE'
      XCID( 70)='XC_GGA_X_VMT_GE'
      XCID( 71)='XC_GGA_X_VMT_PBE'
      XCID( 72)='XC_MGGA_C_CS'
      XCID( 73)='XC_MGGA_C_MN12_SX'
      XCID( 74)='XC_MGGA_C_MN12_L'
      XCID( 75)='XC_MGGA_C_M11_L'
      XCID( 76)='XC_MGGA_C_M11'
      XCID( 77)='XC_MGGA_C_M08_SO'
      XCID( 78)='XC_MGGA_C_M08_HX'
      XCID( 79)='XC_GGA_C_N12_SX'
      XCID( 80)='XC_GGA_C_N12'
      XCID( 81)='XC_HYB_GGA_X_N12_SX'
      XCID( 82)='XC_GGA_X_N12'
      XCID( 83)='XC_GGA_C_REGTPSS'
      XCID( 84)='XC_GGA_C_OP_XALPHA'
      XCID( 85)='XC_GGA_C_OP_G96'
      XCID( 86)='XC_GGA_C_OP_PBE'
      XCID( 87)='XC_GGA_C_OP_B88'
      XCID( 88)='XC_GGA_C_FT97'
      XCID( 89)='XC_GGA_C_SPBE'
      XCID( 90)='XC_GGA_X_SSB_SW'
      XCID( 91)='XC_GGA_X_SSB'
      XCID( 92)='XC_GGA_X_SSB_D'
      XCID( 93)='XC_GGA_XC_HCTH_407P'
      XCID( 94)='XC_GGA_XC_HCTH_P76'
      XCID( 95)='XC_GGA_XC_HCTH_P14'
      XCID( 96)='XC_GGA_XC_B97_GGA1'
      XCID( 97)='XC_GGA_C_HCTH_A'
      XCID( 98)='XC_GGA_X_BPCCAC'
      XCID( 99)='XC_GGA_C_REVTCA'
      XCID(100)='XC_GGA_C_TCA'
      XCID(101)='XC_GGA_X_PBE'
      XCID(102)='XC_GGA_X_PBE_R'
      XCID(103)='XC_GGA_X_B86'
      XCID(105)='XC_GGA_X_B86_MGC'
      XCID(106)='XC_GGA_X_B88'
      XCID(107)='XC_GGA_X_G96'
      XCID(108)='XC_GGA_X_PW86'
      XCID(109)='XC_GGA_X_PW91'
      XCID(110)='XC_GGA_X_OPTX'
      XCID(111)='XC_GGA_X_DK87_R1'
      XCID(112)='XC_GGA_X_DK87_R2'
      XCID(113)='XC_GGA_X_LG93'
      XCID(114)='XC_GGA_X_FT97_A'
      XCID(115)='XC_GGA_X_FT97_B'
      XCID(116)='XC_GGA_X_PBE_SOL'
      XCID(117)='XC_GGA_X_RPBE'
      XCID(118)='XC_GGA_X_WC'
      XCID(119)='XC_GGA_X_MPW91'
      XCID(120)='XC_GGA_X_AM05'
      XCID(121)='XC_GGA_X_PBEA'
      XCID(122)='XC_GGA_X_MPBE'
      XCID(123)='XC_GGA_X_XPBE'
      XCID(124)='XC_GGA_X_2D_B86_MGC'
      XCID(125)='XC_GGA_X_BAYESIAN'
      XCID(126)='XC_GGA_X_PBE_JSJR'
      XCID(127)='XC_GGA_X_2D_B88'
      XCID(128)='XC_GGA_X_2D_B86'
      XCID(129)='XC_GGA_X_2D_PBE'
      XCID(130)='XC_GGA_C_PBE'
      XCID(131)='XC_GGA_C_LYP'
      XCID(132)='XC_GGA_C_P86'
      XCID(133)='XC_GGA_C_PBE_SOL'
      XCID(134)='XC_GGA_C_PW91'
      XCID(135)='XC_GGA_C_AM05'
      XCID(136)='XC_GGA_C_XPBE'
      XCID(137)='XC_GGA_C_LM'
      XCID(138)='XC_GGA_C_PBE_JRGX'
      XCID(139)='XC_GGA_X_OPTB88_VDW'
      XCID(140)='XC_GGA_X_PBEK1_VDW'
      XCID(141)='XC_GGA_X_OPTPBE_VDW'
      XCID(142)='XC_GGA_X_RGE2'
      XCID(143)='XC_GGA_C_RGE2'
      XCID(144)='XC_GGA_X_RPW86'
      XCID(145)='XC_GGA_X_KT1'
      XCID(146)='XC_GGA_XC_KT2'
      XCID(147)='XC_GGA_C_WL'
      XCID(148)='XC_GGA_C_WI'
      XCID(149)='XC_GGA_X_MB88'
      XCID(150)='XC_GGA_X_SOGGA'
      XCID(151)='XC_GGA_X_SOGGA11'
      XCID(152)='XC_GGA_C_SOGGA11'
      XCID(153)='XC_GGA_C_WI0'
      XCID(154)='XC_GGA_XC_TH1'
      XCID(155)='XC_GGA_XC_TH2'
      XCID(156)='XC_GGA_XC_TH3'
      XCID(157)='XC_GGA_XC_TH4'
      XCID(158)='XC_GGA_X_C09X'
      XCID(159)='XC_GGA_C_SOGGA11_X'
      XCID(160)='XC_GGA_X_LB'
      XCID(161)='XC_GGA_XC_HCTH_93'
      XCID(162)='XC_GGA_XC_HCTH_120'
      XCID(163)='XC_GGA_XC_HCTH_147'
      XCID(164)='XC_GGA_XC_HCTH_407'
      XCID(165)='XC_GGA_XC_EDF1'
      XCID(166)='XC_GGA_XC_XLYP'
      XCID(167)='XC_GGA_XC_KT1'
      XCID(168)='XC_GGA_X_LSPBE'
      XCID(169)='XC_GGA_X_LSRPBE'
      XCID(170)='XC_GGA_XC_B97_D'
      XCID(171)='XC_GGA_X_OPTB86B_VDW'
      XCID(172)='XC_MGGA_C_REVM11'
      XCID(173)='XC_GGA_XC_PBE1W'
      XCID(174)='XC_GGA_XC_MPWLYP1W'
      XCID(175)='XC_GGA_XC_PBELYP1W'
      XCID(176)='XC_GGA_C_ACGGAP'
      XCID(177)='XC_HYB_LDA_XC_LDA0'
      XCID(178)='XC_HYB_LDA_XC_CAM_LDA0'
      XCID(179)='XC_GGA_X_B88_6311G'
      XCID(180)='XC_GGA_X_NCAP'
      XCID(181)='XC_GGA_XC_NCAP'
      XCID(182)='XC_GGA_X_LBM'
      XCID(183)='XC_GGA_X_OL2'
      XCID(184)='XC_GGA_X_APBE'
      XCID(185)='XC_GGA_K_APBE'
      XCID(186)='XC_GGA_C_APBE'
      XCID(187)='XC_GGA_K_TW1'
      XCID(188)='XC_GGA_K_TW2'
      XCID(189)='XC_GGA_K_TW3'
      XCID(190)='XC_GGA_K_TW4'
      XCID(191)='XC_GGA_X_HTBS'
      XCID(192)='XC_GGA_X_AIRY'
      XCID(193)='XC_GGA_X_LAG'
      XCID(194)='XC_GGA_XC_MOHLYP'
      XCID(195)='XC_GGA_XC_MOHLYP2'
      XCID(196)='XC_GGA_XC_TH_FL'
      XCID(197)='XC_GGA_XC_TH_FC'
      XCID(198)='XC_GGA_XC_TH_FCFO'
      XCID(199)='XC_GGA_XC_TH_FCO'
      XCID(200)='XC_GGA_C_OPTC'
      XCID(201)='XC_MGGA_X_LTA'
      XCID(202)='XC_MGGA_X_TPSS'
      XCID(203)='XC_MGGA_X_M06_L'
      XCID(204)='XC_MGGA_X_GVT4'
      XCID(205)='XC_MGGA_X_TAU_HCTH'
      XCID(206)='XC_MGGA_X_BR89'
      XCID(207)='XC_MGGA_X_BJ06'
      XCID(208)='XC_MGGA_X_TB09'
      XCID(209)='XC_MGGA_X_RPP09'
      XCID(210)='XC_MGGA_X_2D_PRHG07'
      XCID(211)='XC_MGGA_X_2D_PRHG07_PRP10'
      XCID(212)='XC_MGGA_X_REVTPSS'
      XCID(213)='XC_MGGA_X_PKZB'
      XCID(214)='XC_MGGA_X_BR89_1'
      XCID(215)='XC_GGA_X_ECMV92'
      XCID(216)='XC_GGA_C_PBE_VWN'
      XCID(217)='XC_GGA_C_P86_FT'
      XCID(218)='XC_GGA_K_RATIONAL_P'
      XCID(219)='XC_GGA_K_PG1'
      XCID(220)='XC_MGGA_K_PGSL025'
      XCID(221)='XC_MGGA_X_MS0'
      XCID(222)='XC_MGGA_X_MS1'
      XCID(223)='XC_MGGA_X_MS2'
      XCID(224)='XC_HYB_MGGA_X_MS2H'
      XCID(225)='XC_MGGA_X_TH'
      XCID(226)='XC_MGGA_X_M11_L'
      XCID(227)='XC_MGGA_X_MN12_L'
      XCID(228)='XC_MGGA_X_MS2_REV'
      XCID(229)='XC_MGGA_XC_CC06'
      XCID(230)='XC_MGGA_X_MK00'
      XCID(231)='XC_MGGA_C_TPSS'
      XCID(232)='XC_MGGA_C_VSXC'
      XCID(233)='XC_MGGA_C_M06_L'
      XCID(234)='XC_MGGA_C_M06_HF'
      XCID(235)='XC_MGGA_C_M06'
      XCID(236)='XC_MGGA_C_M06_2X'
      XCID(237)='XC_MGGA_C_M05'
      XCID(238)='XC_MGGA_C_M05_2X'
      XCID(239)='XC_MGGA_C_PKZB'
      XCID(240)='XC_MGGA_C_BC95'
      XCID(241)='XC_MGGA_C_REVTPSS'
      XCID(242)='XC_MGGA_XC_TPSSLYP1W'
      XCID(243)='XC_MGGA_X_MK00B'
      XCID(244)='XC_MGGA_X_BLOC'
      XCID(245)='XC_MGGA_X_MODTPSS'
      XCID(246)='XC_GGA_C_PBELOC'
      XCID(247)='XC_MGGA_C_TPSSLOC'
      XCID(248)='XC_HYB_MGGA_X_MN12_SX'
      XCID(249)='XC_MGGA_X_MBEEF'
      XCID(250)='XC_MGGA_X_MBEEFVDW'
      XCID(251)='XC_MGGA_C_TM'
      XCID(252)='XC_GGA_C_P86VWN'
      XCID(253)='XC_GGA_C_P86VWN_FT'
      XCID(254)='XC_MGGA_XC_B97M_V'
      XCID(255)='XC_GGA_XC_VV10'
      XCID(256)='XC_MGGA_X_JK'
      XCID(257)='XC_MGGA_X_MVS'
      XCID(258)='XC_GGA_C_PBEFE'
      XCID(259)='XC_LDA_XC_KSDT'
      XCID(260)='XC_MGGA_X_MN15_L'
      XCID(261)='XC_MGGA_C_MN15_L'
      XCID(262)='XC_GGA_C_OP_PW91'
      XCID(263)='XC_MGGA_X_SCAN'
      XCID(264)='XC_HYB_MGGA_X_SCAN0'
      XCID(265)='XC_GGA_X_PBEFE'
      XCID(266)='XC_HYB_GGA_XC_B97_1P'
      XCID(267)='XC_MGGA_C_SCAN'
      XCID(268)='XC_HYB_MGGA_X_MN15'
      XCID(269)='XC_MGGA_C_MN15'
      XCID(270)='XC_GGA_X_CAP'
      XCID(271)='XC_GGA_X_EB88'
      XCID(272)='XC_GGA_C_PBE_MOL'
      XCID(273)='XC_HYB_GGA_XC_PBE_MOL0'
      XCID(274)='XC_HYB_GGA_XC_PBE_SOL0'
      XCID(275)='XC_HYB_GGA_XC_PBEB0'
      XCID(276)='XC_HYB_GGA_XC_PBE_MOLB0'
      XCID(277)='XC_GGA_K_ABSP3'
      XCID(278)='XC_GGA_K_ABSP4'
      XCID(279)='XC_HYB_MGGA_X_BMK'
      XCID(280)='XC_GGA_C_BMK'
      XCID(281)='XC_GGA_C_TAU_HCTH'
      XCID(282)='XC_HYB_MGGA_X_TAU_HCTH'
      XCID(283)='XC_GGA_C_HYB_TAU_HCTH'
      XCID(284)='XC_MGGA_X_B00'
      XCID(285)='XC_GGA_X_BEEFVDW'
      XCID(286)='XC_GGA_XC_BEEFVDW'
      XCID(287)='XC_LDA_C_CHACHIYO'
      XCID(288)='XC_MGGA_XC_HLE17'
      XCID(289)='XC_LDA_C_LP96'
      XCID(290)='XC_HYB_GGA_XC_PBE50'
      XCID(291)='XC_GGA_X_PBETRANS'
      XCID(292)='XC_MGGA_C_SCAN_RVV10'
      XCID(293)='XC_MGGA_X_REVM06_L'
      XCID(294)='XC_MGGA_C_REVM06_L'
      XCID(295)='XC_HYB_MGGA_X_M08_HX'
      XCID(296)='XC_HYB_MGGA_X_M08_SO'
      XCID(297)='XC_HYB_MGGA_X_M11'
      XCID(298)='XC_GGA_X_CHACHIYO'
      XCID(299)='XC_MGGA_X_RTPSS'
      XCID(300)='XC_MGGA_X_MS2B'
      XCID(301)='XC_MGGA_X_MS2BS'
      XCID(302)='XC_MGGA_X_MVSB'
      XCID(303)='XC_MGGA_X_MVSBS'
      XCID(304)='XC_HYB_MGGA_X_REVM11'
      XCID(305)='XC_HYB_MGGA_X_REVM06'
      XCID(306)='XC_MGGA_C_REVM06'
      XCID(307)='XC_LDA_C_CHACHIYO_MOD'
      XCID(308)='XC_LDA_C_KARASIEV_MOD'
      XCID(309)='XC_GGA_C_CHACHIYO'
      XCID(310)='XC_HYB_MGGA_X_M06_SX'
      XCID(311)='XC_MGGA_C_M06_SX'
      XCID(312)='XC_GGA_X_REVSSB_D'
      XCID(313)='XC_GGA_C_CCDF'
      XCID(314)='XC_HYB_GGA_XC_HFLYP'
      XCID(315)='XC_HYB_GGA_XC_B3P86_NWCHEM'
      XCID(316)='XC_GGA_X_PW91_MOD'
      XCID(317)='XC_LDA_C_W20'
      XCID(318)='XC_LDA_XC_CORRKSDT'
      XCID(319)='XC_MGGA_X_FT98'
      XCID(320)='XC_GGA_X_PBE_MOD'
      XCID(321)='XC_GGA_X_PBE_GAUSSIAN'
      XCID(322)='XC_GGA_C_PBE_GAUSSIAN'
      XCID(323)='XC_MGGA_C_TPSS_GAUSSIAN'
      XCID(324)='XC_GGA_X_NCAPR'
      XCID(325)='XC_HYB_GGA_XC_RELPBE0'
      XCID(327)='XC_GGA_XC_B97_3C'
      XCID(387)='XC_MGGA_C_CC'
      XCID(388)='XC_MGGA_C_CCALDA'
      XCID(389)='XC_HYB_MGGA_XC_BR3P86'
      XCID(390)='XC_HYB_GGA_XC_CASE21'
      XCID(391)='XC_MGGA_C_RREGTM'
      XCID(392)='XC_HYB_GGA_XC_PBE_2X'
      XCID(393)='XC_HYB_GGA_XC_PBE38'
      XCID(394)='XC_HYB_GGA_XC_B3LYP3'
      XCID(395)='XC_HYB_GGA_XC_CAM_O3LYP'
      XCID(396)='XC_HYB_MGGA_XC_TPSS0'
      XCID(397)='XC_MGGA_C_B94'
      XCID(398)='XC_HYB_MGGA_XC_B94_HYB'
      XCID(399)='XC_HYB_GGA_XC_WB97X_D3'
      XCID(400)='XC_HYB_GGA_XC_LC_BLYP'
      XCID(401)='XC_HYB_GGA_XC_B3PW91'
      XCID(402)='XC_HYB_GGA_XC_B3LYP'
      XCID(403)='XC_HYB_GGA_XC_B3P86'
      XCID(404)='XC_HYB_GGA_XC_O3LYP'
      XCID(405)='XC_HYB_GGA_XC_MPW1K'
      XCID(406)='XC_HYB_GGA_XC_PBEH'
      XCID(407)='XC_HYB_GGA_XC_B97'
      XCID(408)='XC_HYB_GGA_XC_B97_1'
      XCID(409)='XC_HYB_GGA_XC_APF'
      XCID(410)='XC_HYB_GGA_XC_B97_2'
      XCID(411)='XC_HYB_GGA_XC_X3LYP'
      XCID(412)='XC_HYB_GGA_XC_B1WC'
      XCID(413)='XC_HYB_GGA_XC_B97_K'
      XCID(414)='XC_HYB_GGA_XC_B97_3'
      XCID(415)='XC_HYB_GGA_XC_MPW3PW'
      XCID(416)='XC_HYB_GGA_XC_B1LYP'
      XCID(417)='XC_HYB_GGA_XC_B1PW91'
      XCID(418)='XC_HYB_GGA_XC_MPW1PW'
      XCID(419)='XC_HYB_GGA_XC_MPW3LYP'
      XCID(420)='XC_HYB_GGA_XC_SB98_1A'
      XCID(421)='XC_HYB_GGA_XC_SB98_1B'
      XCID(422)='XC_HYB_GGA_XC_SB98_1C'
      XCID(423)='XC_HYB_GGA_XC_SB98_2A'
      XCID(424)='XC_HYB_GGA_XC_SB98_2B'
      XCID(425)='XC_HYB_GGA_XC_SB98_2C'
      XCID(426)='XC_HYB_GGA_X_SOGGA11_X'
      XCID(427)='XC_HYB_GGA_XC_HSE03'
      XCID(428)='XC_HYB_GGA_XC_HSE06'
      XCID(429)='XC_HYB_GGA_XC_HJS_PBE'
      XCID(430)='XC_HYB_GGA_XC_HJS_PBE_SOL'
      XCID(431)='XC_HYB_GGA_XC_HJS_B88'
      XCID(432)='XC_HYB_GGA_XC_HJS_B97X'
      XCID(433)='XC_HYB_GGA_XC_CAM_B3LYP'
      XCID(434)='XC_HYB_GGA_XC_TUNED_CAM_B3LYP'
      XCID(435)='XC_HYB_GGA_XC_BHANDH'
      XCID(436)='XC_HYB_GGA_XC_BHANDHLYP'
      XCID(437)='XC_HYB_GGA_XC_MB3LYP_RC04'
      XCID(438)='XC_HYB_MGGA_X_M05'
      XCID(439)='XC_HYB_MGGA_X_M05_2X'
      XCID(440)='XC_HYB_MGGA_XC_B88B95'
      XCID(441)='XC_HYB_MGGA_XC_B86B95'
      XCID(442)='XC_HYB_MGGA_XC_PW86B95'
      XCID(443)='XC_HYB_MGGA_XC_BB1K'
      XCID(444)='XC_HYB_MGGA_X_M06_HF'
      XCID(445)='XC_HYB_MGGA_XC_MPW1B95'
      XCID(446)='XC_HYB_MGGA_XC_MPWB1K'
      XCID(447)='XC_HYB_MGGA_XC_X1B95'
      XCID(448)='XC_HYB_MGGA_XC_XB1K'
      XCID(449)='XC_HYB_MGGA_X_M06'
      XCID(450)='XC_HYB_MGGA_X_M06_2X'
      XCID(451)='XC_HYB_MGGA_XC_PW6B95'
      XCID(452)='XC_HYB_MGGA_XC_PWB6K'
      XCID(453)='XC_HYB_GGA_XC_MPWLYP1M'
      XCID(454)='XC_HYB_GGA_XC_REVB3LYP'
      XCID(455)='XC_HYB_GGA_XC_CAMY_BLYP'
      XCID(456)='XC_HYB_GGA_XC_PBE0_13'
      XCID(457)='XC_HYB_MGGA_XC_TPSSH'
      XCID(458)='XC_HYB_MGGA_XC_REVTPSSH'
      XCID(459)='XC_HYB_GGA_XC_B3LYPS'
      XCID(460)='XC_HYB_GGA_XC_QTP17'
      XCID(461)='XC_HYB_GGA_XC_B3LYP_MCM1'
      XCID(462)='XC_HYB_GGA_XC_B3LYP_MCM2'
      XCID(463)='XC_HYB_GGA_XC_WB97'
      XCID(464)='XC_HYB_GGA_XC_WB97X'
      XCID(465)='XC_HYB_GGA_XC_LRC_WPBEH'
      XCID(466)='XC_HYB_GGA_XC_WB97X_V'
      XCID(467)='XC_HYB_GGA_XC_LCY_PBE'
      XCID(468)='XC_HYB_GGA_XC_LCY_BLYP'
      XCID(469)='XC_HYB_GGA_XC_LC_VV10'
      XCID(470)='XC_HYB_GGA_XC_CAMY_B3LYP'
      XCID(471)='XC_HYB_GGA_XC_WB97X_D'
      XCID(472)='XC_HYB_GGA_XC_HPBEINT'
      XCID(473)='XC_HYB_GGA_XC_LRC_WPBE'
      XCID(474)='XC_HYB_MGGA_X_MVSH'
      XCID(475)='XC_HYB_GGA_XC_B3LYP5'
      XCID(476)='XC_HYB_GGA_XC_EDF2'
      XCID(477)='XC_HYB_GGA_XC_CAP0'
      XCID(478)='XC_HYB_GGA_XC_LC_WPBE'
      XCID(479)='XC_HYB_GGA_XC_HSE12'
      XCID(480)='XC_HYB_GGA_XC_HSE12S'
      XCID(481)='XC_HYB_GGA_XC_HSE_SOL'
      XCID(482)='XC_HYB_GGA_XC_CAM_QTP_01'
      XCID(483)='XC_HYB_GGA_XC_MPW1LYP'
      XCID(484)='XC_HYB_GGA_XC_MPW1PBE'
      XCID(485)='XC_HYB_GGA_XC_KMLYP'
      XCID(486)='XC_HYB_GGA_XC_LC_WPBE_WHS'
      XCID(487)='XC_HYB_GGA_XC_LC_WPBEH_WHS'
      XCID(488)='XC_HYB_GGA_XC_LC_WPBE08_WHS'
      XCID(489)='XC_HYB_GGA_XC_LC_WPBESOL_WHS'
      XCID(490)='XC_HYB_GGA_XC_CAM_QTP_00'
      XCID(491)='XC_HYB_GGA_XC_CAM_QTP_02'
      XCID(492)='XC_HYB_GGA_XC_LC_QTP'
      XCID(493)='XC_MGGA_X_RSCAN'
      XCID(494)='XC_MGGA_C_RSCAN'
      XCID(495)='XC_GGA_X_S12G'
      XCID(496)='XC_HYB_GGA_X_S12H'
      XCID(497)='XC_MGGA_X_R2SCAN'
      XCID(498)='XC_MGGA_C_R2SCAN'
      XCID(499)='XC_HYB_GGA_XC_BLYP35'
      XCID(500)='XC_GGA_K_VW'
      XCID(501)='XC_GGA_K_GE2'
      XCID(502)='XC_GGA_K_GOLDEN'
      XCID(503)='XC_GGA_K_YT65'
      XCID(504)='XC_GGA_K_BALTIN'
      XCID(505)='XC_GGA_K_LIEB'
      XCID(506)='XC_GGA_K_ABSP1'
      XCID(507)='XC_GGA_K_ABSP2'
      XCID(508)='XC_GGA_K_GR'
      XCID(509)='XC_GGA_K_LUDENA'
      XCID(510)='XC_GGA_K_GP85'
      XCID(511)='XC_GGA_K_PEARSON'
      XCID(512)='XC_GGA_K_OL1'
      XCID(513)='XC_GGA_K_OL2'
      XCID(514)='XC_GGA_K_FR_B88'
      XCID(515)='XC_GGA_K_FR_PW86'
      XCID(516)='XC_GGA_K_DK'
      XCID(517)='XC_GGA_K_PERDEW'
      XCID(518)='XC_GGA_K_VSK'
      XCID(519)='XC_GGA_K_VJKS'
      XCID(520)='XC_GGA_K_ERNZERHOF'
      XCID(521)='XC_GGA_K_LC94'
      XCID(522)='XC_GGA_K_LLP'
      XCID(523)='XC_GGA_K_THAKKAR'
      XCID(524)='XC_GGA_X_WPBEH'
      XCID(525)='XC_GGA_X_HJS_PBE'
      XCID(526)='XC_GGA_X_HJS_PBE_SOL'
      XCID(527)='XC_GGA_X_HJS_B88'
      XCID(528)='XC_GGA_X_HJS_B97X'
      XCID(529)='XC_GGA_X_ITYH'
      XCID(530)='XC_GGA_X_SFAT'
      XCID(531)='XC_HYB_MGGA_XC_WB97M_V'
      XCID(532)='XC_LDA_X_REL'
      XCID(533)='XC_GGA_X_SG4'
      XCID(534)='XC_GGA_C_SG4'
      XCID(535)='XC_GGA_X_GG99'
      XCID(536)='XC_LDA_XC_1D_EHWLRG_1'
      XCID(537)='XC_LDA_XC_1D_EHWLRG_2'
      XCID(538)='XC_LDA_XC_1D_EHWLRG_3'
      XCID(539)='XC_GGA_X_PBEPOW'
      XCID(540)='XC_MGGA_X_TM'
      XCID(541)='XC_MGGA_X_VT84'
      XCID(542)='XC_MGGA_X_SA_TPSS'
      XCID(543)='XC_MGGA_K_PC07'
      XCID(544)='XC_GGA_X_KGG99'
      XCID(545)='XC_GGA_XC_HLE16'
      XCID(546)='XC_LDA_X_ERF'
      XCID(547)='XC_LDA_XC_LP_A'
      XCID(548)='XC_LDA_XC_LP_B'
      XCID(549)='XC_LDA_X_RAE'
      XCID(550)='XC_LDA_K_ZLP'
      XCID(551)='XC_LDA_C_MCWEENY'
      XCID(552)='XC_LDA_C_BR78'
      XCID(553)='XC_GGA_C_SCAN_E0'
      XCID(554)='XC_LDA_C_PK09'
      XCID(555)='XC_GGA_C_GAPC'
      XCID(556)='XC_GGA_C_GAPLOC'
      XCID(557)='XC_GGA_C_ZVPBEINT'
      XCID(558)='XC_GGA_C_ZVPBESOL'
      XCID(559)='XC_GGA_C_TM_LYP'
      XCID(560)='XC_GGA_C_TM_PBE'
      XCID(561)='XC_GGA_C_W94'
      XCID(562)='XC_MGGA_C_KCIS'
      XCID(563)='XC_HYB_MGGA_XC_B0KCIS'
      XCID(564)='XC_MGGA_XC_LP90'
      XCID(565)='XC_GGA_C_CS1'
      XCID(566)='XC_HYB_MGGA_XC_MPW1KCIS'
      XCID(567)='XC_HYB_MGGA_XC_MPWKCIS1K'
      XCID(568)='XC_HYB_MGGA_XC_PBE1KCIS'
      XCID(569)='XC_HYB_MGGA_XC_TPSS1KCIS'
      XCID(570)='XC_GGA_X_B88M'
      XCID(571)='XC_MGGA_C_B88'
      XCID(572)='XC_HYB_GGA_XC_B5050LYP'
      XCID(573)='XC_LDA_C_OW_LYP'
      XCID(574)='XC_LDA_C_OW'
      XCID(575)='XC_MGGA_X_GX'
      XCID(576)='XC_MGGA_X_PBE_GX'
      XCID(577)='XC_LDA_XC_GDSMFB'
      XCID(578)='XC_LDA_C_GK72'
      XCID(579)='XC_LDA_C_KARASIEV'
      XCID(580)='XC_LDA_K_LP96'
      XCID(581)='XC_MGGA_X_REVSCAN'
      XCID(582)='XC_MGGA_C_REVSCAN'
      XCID(583)='XC_HYB_MGGA_X_REVSCAN0'
      XCID(584)='XC_MGGA_C_SCAN_VV10'
      XCID(585)='XC_MGGA_C_REVSCAN_VV10'
      XCID(586)='XC_MGGA_X_BR89_EXPLICIT'
      XCID(587)='XC_GGA_XC_KT3'
      XCID(588)='XC_HYB_LDA_XC_BN05'
      XCID(589)='XC_HYB_GGA_XC_LB07'
      XCID(590)='XC_LDA_C_PMGB06'
      XCID(591)='XC_GGA_K_GDS08'
      XCID(592)='XC_GGA_K_GHDS10'
      XCID(593)='XC_GGA_K_GHDS10R'
      XCID(594)='XC_GGA_K_TKVLN'
      XCID(595)='XC_GGA_K_PBE3'
      XCID(596)='XC_GGA_K_PBE4'
      XCID(597)='XC_GGA_K_EXP4'
      XCID(598)='XC_HYB_MGGA_XC_B98'
      XCID(599)='XC_LDA_XC_TIH'
      XCID(600)='XC_LDA_X_1D_EXPONENTIAL'
      XCID(601)='XC_GGA_X_SFAT_PBE'
      XCID(602)='XC_MGGA_X_BR89_EXPLICIT_1'
      XCID(603)='XC_MGGA_X_REGTPSS'
      XCID(604)='XC_GGA_X_FD_LB94'
      XCID(605)='XC_GGA_X_FD_REVLB94'
      XCID(606)='XC_GGA_C_ZVPBELOC'
      XCID(607)='XC_HYB_GGA_XC_APBE0'
      XCID(608)='XC_HYB_GGA_XC_HAPBE'
      XCID(609)='XC_MGGA_X_2D_JS17'
      XCID(610)='XC_HYB_GGA_XC_RCAM_B3LYP'
      XCID(611)='XC_HYB_GGA_XC_WC04'
      XCID(612)='XC_HYB_GGA_XC_WP04'
      XCID(613)='XC_GGA_K_LKT'
      XCID(614)='XC_HYB_GGA_XC_CAMH_B3LYP'
      XCID(615)='XC_HYB_GGA_XC_WHPBE0'
      XCID(616)='XC_GGA_K_PBE2'
      XCID(617)='XC_MGGA_K_L04'
      XCID(618)='XC_MGGA_K_L06'
      XCID(619)='XC_GGA_K_VT84F'
      XCID(620)='XC_GGA_K_LGAP'
      XCID(621)='XC_MGGA_K_RDA'
      XCID(622)='XC_GGA_X_ITYH_OPTX'
      XCID(623)='XC_GGA_X_ITYH_PBE'
      XCID(624)='XC_GGA_C_LYPR'
      XCID(625)='XC_HYB_GGA_XC_LC_BLYP_EA'
      XCID(626)='XC_MGGA_X_REGTM'
      XCID(627)='XC_MGGA_K_GEA2'
      XCID(628)='XC_MGGA_K_GEA4'
      XCID(629)='XC_MGGA_K_CSK1'
      XCID(630)='XC_MGGA_K_CSK4'
      XCID(631)='XC_MGGA_K_CSK_LOC1'
      XCID(632)='XC_MGGA_K_CSK_LOC4'
      XCID(633)='XC_GGA_K_LGAP_GE'
      XCID(634)='XC_MGGA_K_PC07_OPT'
      XCID(635)='XC_GGA_K_TFVW_OPT'
      XCID(636)='XC_HYB_GGA_XC_LC_BOP'
      XCID(637)='XC_HYB_GGA_XC_LC_PBEOP'
      XCID(638)='XC_MGGA_C_KCISK'
      XCID(639)='XC_HYB_GGA_XC_LC_BLYPR'
      XCID(640)='XC_HYB_GGA_XC_MCAM_B3LYP'
      XCID(641)='XC_LDA_X_YUKAWA'
      XCID(642)='XC_MGGA_C_R2SCAN01'
      XCID(643)='XC_MGGA_C_RMGGAC'
      XCID(644)='XC_MGGA_X_MCML'
      XCID(645)='XC_MGGA_X_R2SCAN01'
      XCID(646)='XC_HYB_GGA_X_CAM_S12G'
      XCID(647)='XC_HYB_GGA_X_CAM_S12H'
      XCID(648)='XC_MGGA_X_RPPSCAN'
      XCID(649)='XC_MGGA_C_RPPSCAN'
      XCID(650)='XC_MGGA_X_R4SCAN'
      XCID(651)='XC_MGGA_X_VCML'
      XCID(652)='XC_MGGA_XC_VCML_RVV10'
      XCID(653)='XC_HYB_LDA_X_ERF'
      XCID(654)='XC_LDA_C_PW_ERF'
      XCID(655)='XC_GGA_X_PBE_ERF_GWS'
      XCID(656)='XC_HYB_GGA_X_PBE_ERF_GWS'
      XCID(657)='XC_GGA_C_PBE_ERF_GWS'
      XCID(658)='XC_HYB_MGGA_XC_GAS22'
      XCID(659)='XC_HYB_MGGA_XC_R2SCANH'
      XCID(660)='XC_HYB_MGGA_XC_R2SCAN0'
      XCID(661)='XC_HYB_MGGA_XC_R2SCAN50'
      XCID(681)='XC_HYB_GGA_XC_CAM_PBEH'
      XCID(682)='XC_HYB_GGA_XC_CAMY_PBEH'
      XCID(683)='XC_LDA_C_UPW92'
      XCID(684)='XC_LDA_C_RPW92'
      XCID(685)='XC_MGGA_X_TLDA'
      XCID(686)='XC_MGGA_X_EDMGGA'
      XCID(687)='XC_MGGA_X_GDME_NV'
      XCID(688)='XC_MGGA_X_RLDA'
      XCID(689)='XC_MGGA_X_GDME_0'
      XCID(690)='XC_MGGA_X_GDME_KOS'
      XCID(691)='XC_MGGA_X_GDME_VT'
      XCID(692)='XC_LDA_X_SLOC'
      XCID(693)='XC_MGGA_X_REVTM'
      XCID(694)='XC_MGGA_C_REVTM'
      XCID(695)='XC_HYB_MGGA_XC_EDMGGAH'
      XCID(696)='XC_MGGA_X_MBRXC_BG'
      XCID(697)='XC_MGGA_X_MBRXH_BG'
      XCID(698)='XC_MGGA_X_HLTA'
      XCID(699)='XC_MGGA_C_HLTAPW'
      XCID(700)='XC_MGGA_X_SCANL'
      XCID(701)='XC_MGGA_X_REVSCANL'
      XCID(702)='XC_MGGA_C_SCANL'
      XCID(703)='XC_MGGA_C_SCANL_RVV10'
      XCID(704)='XC_MGGA_C_SCANL_VV10'
      XCID(705)='XC_HYB_MGGA_X_JS18'
      XCID(706)='XC_HYB_MGGA_X_PJS18'
      XCID(707)='XC_MGGA_X_TASK'
      XCID(708)='XC_HYB_GGA_X_LCGAU'
      XCID(709)='XC_HYB_GGA_X_LCGAU_CORE'
      XCID(710)='XC_HYB_GGA_X_LC2GAU'
      XCID(711)='XC_MGGA_X_MGGAC'
      XCID(712)='XC_GGA_C_MGGAC'
      XCID(713)='XC_HYB_GGA_XC_B2PLYP'
      XCID(714)='XC_HYB_GGA_XC_SRC1_BLYP'
      XCID(715)='XC_HYB_GGA_XC_SRC2_BLYP'
      XCID(716)='XC_MGGA_X_MBR'
      XCID(717)='XC_HYB_GGA_XC_HISS'     
      XCID(718)='XC_MGGA_X_R2SCANL'
      XCID(719)='XC_MGGA_C_R2SCANL'
      XCID(720)='XC_HYB_MGGA_XC_LC_TMLYP'
      XCID(721)='XC_HYB_GGA_XC_B2GPPLYP'
      XCID(722)='XC_HYB_GGA_XC_WB2PLYP'
      XCID(723)='XC_HYB_GGA_XC_WB2GPPLYP'
      XCID(724)='XC_MGGA_X_MTASK'
      XCID(725)='XC_HYB_GGA_XC_PBE0_DH'
      XCID(726)='XC_HYB_GGA_XC_PBE0_2'
      XCID(727)='XC_HYB_GGA_XC_PBE_QIDH'
      XCID(728)='XC_HYB_GGA_XC_LS1DH_PBE' 
      XCID(734)='XC_GGA_X_Q1D'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$REPORT(NFIL)
!     **************************************************************************
!     ** THIS IS A MODIFIED VERSION XCINFO OF LIBXC                           **
!     **************************************************************************
      USE XC_F03_LIB_M   , ONLY: XC_F03_FUNC_INFO_T      & !TYPE
     &                          ,XC_F03_FUNC_REFERENCE_T & !TYPE
     &                          ,XC_F03_FUNC_GET_INFO    &  !FUNCTION
     &                          ,XC_F03_FUNC_INFO_GET_KIND & !FUNCTION
     &                          ,XC_F03_FUNC_INFO_GET_FAMILY & !FUNCTION
     &                          ,XC_F03_FUNC_INFO_GET_REFERENCES & !FUNCTION
     &                          ,XC_F03_FUNC_REFERENCE_GET_REF & !FUNCTION
     &                          ,XC_F03_FUNC_INFO_GET_NAME & !FUNCTION
     &                          ,XC_EXCHANGE &
     &                          ,XC_CORRELATION &
     &                          ,XC_EXCHANGE_CORRELATION &
     &                          ,XC_KINETIC &
     &                          ,XC_FAMILY_LDA &
     &                          ,XC_FAMILY_HYB_LDA &
     &                          ,XC_FAMILY_GGA &
     &                          ,XC_FAMILY_HYB_GGA &
     &                          ,XC_FAMILY_MGGA &
     &                          ,XC_FAMILY_HYB_MGGA &
     &                          ,XC_FAMILY_LCA &
     &                          ,XC_MAX_REFERENCES 
      USE PAWLIBXC_MODULE, ONLY: NXC &
     &                          ,XC_FUNC
      IMPLICIT NONE
      INTEGER(4)       ,INTENT(IN) :: NFIL
      TYPE(XC_F03_FUNC_INFO_T)     :: XC_INFO
      TYPE(XC_F03_FUNC_REFERENCE_T):: XC_REF
      INTEGER(4)                   :: IXC
      INTEGER(4)                   :: NUMBER,I
      CHARACTER(LEN=128)           :: S1, S2
!     **************************************************************************
      DO IXC=1,NXC
        XC_INFO = XC_F03_FUNC_GET_INFO(XC_FUNC(IXC))

        SELECT CASE(XC_F03_FUNC_INFO_GET_KIND(XC_INFO))
          CASE(XC_EXCHANGE)
            WRITE(NFIL, '(A)') 'EXCHANGE'
          CASE(XC_CORRELATION)
            WRITE(NFIL, '(A)') 'CORRELATION'
          CASE(XC_EXCHANGE_CORRELATION)
            WRITE(NFIL, '(A)') 'EXCHANGE-CORRELATION'
          CASE(XC_KINETIC)
            WRITE(NFIL, '(A)') 'KINETIC'
        END SELECT

        S1 = XC_F03_FUNC_INFO_GET_NAME(XC_INFO)
        SELECT CASE(XC_F03_FUNC_INFO_GET_FAMILY(XC_INFO))
          CASE (XC_FAMILY_LDA);       WRITE(S2,'(A)') "LDA"
          CASE (XC_FAMILY_HYB_LDA);   WRITE(S2,'(A)') "HYBRID LDA"
          CASE (XC_FAMILY_GGA);       WRITE(S2,'(A)') "GGA"
          CASE (XC_FAMILY_HYB_GGA);   WRITE(S2,'(A)') "HYBRID GGA"
          CASE (XC_FAMILY_MGGA);      WRITE(S2,'(A)') "MGGA"
          CASE (XC_FAMILY_HYB_MGGA);  WRITE(S2,'(A)') "HYBRID MGGA"
          CASE (XC_FAMILY_LCA);       WRITE(S2,'(A)') "LCA"
        END SELECT
        WRITE(NFIL, '(4A)') TRIM(S1), ' (', TRIM(S2), ')'

        DO I=0, XC_MAX_REFERENCES-1
!         __GET INFO ON REFERENCE NUMBER I; CANT REUSE THE LOOP VARIABLE
!         __DUE TO FORTRAN RESTRICTIONS FOR INTENT(INOUT)
          NUMBER = I
          XC_REF = XC_F03_FUNC_INFO_GET_REFERENCES(XC_INFO, NUMBER)
          IF( NUMBER < 0 ) THEN
            EXIT
          END IF
          WRITE(NFIL, '(A,I1,2A)') '[', I+1, '] ' &
    &         ,TRIM(XC_F03_FUNC_REFERENCE_GET_REF(XC_REF))
        END DO
      END DO
      RETURN
      END


!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE PAWLIBXC_MODULE, ONLY: NXC &
     &                          ,XC_FUNC
      USE XC_F03_LIB_M   , ONLY: XC_F03_FUNC_INFO_T      & !TYPE
     &                          ,XC_FAMILY_LDA &
     &                          ,XC_FAMILY_HYB_LDA &
     &                          ,XC_FAMILY_GGA &
     &                          ,XC_FAMILY_HYB_GGA &
     &                          ,XC_FAMILY_MGGA &
     &                          ,XC_FAMILY_HYB_MGGA &
     &                          ,XC_FAMILY_UNKNOWN &
     &                          ,XC_FAMILY_LCA &
     &                          ,XC_F03_FUNC_GET_INFO &
     &                          ,XC_F03_FUNC_INFO_GET_FAMILY
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
      TYPE(XC_F03_FUNC_INFO_T):: XC_INFO
      INTEGER(4)              :: I
!     **************************************************************************
      IF(ID.EQ.'GRADIENT') THEN 
        IF(NXC.EQ.0) THEN
          CALL ERROR$MSG('FUNCTIONALS ARE NOT YET SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PAWLIBXC$GETL4')
        END IF
        VAL=.FALSE.
        DO I=1,NXC
          XC_INFO=XC_F03_FUNC_GET_INFO(XC_FUNC(I))
          SELECT CASE (XC_F03_FUNC_INFO_GET_FAMILY(XC_INFO))
          CASE(XC_FAMILY_LDA,XC_FAMILY_HYB_LDA) 
            VAL=VAL
          CASE(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA &
       &       ,XC_FAMILY_MGGA,XC_FAMILY_HYB_MGGA)
            VAL=.TRUE.
          CASE(XC_FAMILY_UNKNOWN) 
            CALL ERROR$MSG('XC_FAMILY UNKOWN')
            CALL ERROR$STOP('PAWLIBXC$GETL4')
          CASE DEFAULT
            CALL ERROR$MSG('XC_FAMILY NOT IMPLEMENTED')
            CALL ERROR$STOP('PAWLIBXC$GETL4')
          END SELECT
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PAWLIBXC$SETL4')
      END IF
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE PAWLIBXC$GETCH(ID,VAL)
!!$!     **************************************************************************
!!$!     **************************************************************************
!!$      IMPLICIT NONE
!!$      CHARACTER(*),INTENT(IN) :: ID
!!$      CHARACTER(*),INTENT(IN) :: VAL
!!$!     **************************************************************************
!!$      IF(ID.EQ.'TGRA') THEN 
!!$
!!$      ELSE
!!$        CALL ERROR$MSG('ID NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID',ID)
!!$        CALL ERROR$STOP('PAWLIBXC$SETL4')
!!$      END IF
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$SETCH(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE XC_F03_LIB_M   , ONLY: XC_POLARIZED &
     &                          ,XC_LDA_X &
     &                          ,XC_LDA_C_VWN &
     &                          ,XC_GGA_X_PBE &
     &                          ,XC_GGA_C_PBE &
     &                          ,XC_MGGA_X_R2SCAN &
     &                          ,XC_MGGA_C_R2SCAN &
     &                          ,XC_F03_FUNC_INIT ! FUNCTION
      USE PAWLIBXC_MODULE, ONLY: NXC &
     &                          ,XC_FUNC
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'FUNCTIONAL') THEN 
!       ========================================================================
!       == VOSKO WILK NUSAIR
!       ========================================================================
        IF(VAL.EQ.'VWN') THEN
          IF(ALLOCATED(XC_FUNC)) THEN
            PRINT*,'ERROR: XC_FUNC IS ALREADY SET'
            STOP 'DFT_LIBXC_SETCH'
          END IF
          NXC=2
          ALLOCATE(XC_FUNC(NXC))          
          CALL XC_F03_FUNC_INIT(XC_FUNC(1),XC_LDA_X,XC_POLARIZED)
          CALL XC_F03_FUNC_INIT(XC_FUNC(2),XC_LDA_C_VWN,XC_POLARIZED)

!       ========================================================================
!       == PERDEW-BURKE ERNZERHOF (PBE) FUNCTIONAL                            ==
!       ========================================================================
        ELSE IF(VAL.EQ.'PBE') THEN
          IF(ALLOCATED(XC_FUNC)) THEN
            PRINT*,'ERROR: XC_FUNC IS ALREADY SET'
            STOP 'DFT_LIBXC_SETCH'
          END IF
          NXC=2
          ALLOCATE(XC_FUNC(NXC))          
          CALL XC_F03_FUNC_INIT(XC_FUNC(1),XC_GGA_X_PBE,XC_POLARIZED)
          CALL XC_F03_FUNC_INIT(XC_FUNC(2),XC_GGA_C_PBE,XC_POLARIZED)
!
!       ========================================================================
!       == SCAN-R2 FUNCTIONAL                                                 ==
!       == FURNESS,KAPLAN,NING,PERDEW,SUN, J.PHYS.CHEM.LETT.11, 8208 (2020)   ==
!       ==                                 (DOI: 10.1021/ACS.JPCLETT.0C02405) ==
!       == FURNESS,KAPLAN,NING,PERDEW,SUN, J.PHYS.CHEM.LETT.11, 9248 (2020)   ==
!       ==                                (DOI: 10.1021/ACS.JPCLETT.0C03077)  ==
!       ========================================================================
        ELSE IF(VAL.EQ.'SCANR2') THEN
          IF(ALLOCATED(XC_FUNC)) THEN
            PRINT*,'ERROR: XC_FUNC IS ALREADY SET'
            STOP 'DFT_LIBXC_SETCH'
          END IF
          NXC=2
          ALLOCATE(XC_FUNC(NXC))          
          CALL XC_F03_FUNC_INIT(XC_FUNC(1),XC_MGGA_X_R2SCAN,XC_POLARIZED)
          CALL XC_F03_FUNC_INIT(XC_FUNC(2),XC_MGGA_C_R2SCAN,XC_POLARIZED)
!
       ELSE
          PRINT*,'FUNCTIONAL ID NOT RECOGNIZED'
          PRINT*,'VAL=',VAL
          STOP 'IN DFT_LIBXC1_ADDFUNC'
        END IF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PAWLIBXC$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$SETCHA(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE, INTRINSIC :: ISO_C_BINDING
      USE XC_F03_LIB_M   , ONLY: XC_POLARIZED &
     &                          ,XC_LDA_X &
     &                          ,XC_LDA_C_VWN &
     &                          ,XC_GGA_X_PBE &
     &                          ,XC_GGA_C_PBE &
     &                          ,XC_MGGA_X_R2SCAN &
     &                          ,XC_MGGA_C_R2SCAN &
     &                          ,XC_F03_FUNC_INIT &! FUNCTION
     &                          ,XC_F03_FUNCTIONAL_GET_NUMBER &
     &                          ,XC_F03_FUNCTIONAL_GET_NAME
      USE PAWLIBXC_MODULE, ONLY: NXC &
     &                          ,XC_FUNC &
     &                          ,NIDX &
     &                          ,XCID
      use strings_module
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
      INTEGER(4)              :: I
      INTEGER(C_INT)          :: J
      LOGICAL(4)              :: TCHK
      CHARACTER(128)          :: NAME
!     **************************************************************************
!PRINT*,'STARTING PAWLIBXC$SETCHA'
      IF(ID.EQ.'FUNCTIONAL') THEN 
        IF(ALLOCATED(XC_FUNC)) THEN
          PRINT*,'ERROR: XC_FUNC IS ALREADY SET'
          STOP 'DFT_LIBXC_SETCHA'
        END IF
        NXC=LEN
        ALLOCATE(XC_FUNC(NXC))          
        DO I=1,NXC
!
!         ======================================================================
!         == RESOLVE FUNCTIONAL IDS AND INITIALIZE LIBXC                      ==
!         == FORTRAN ROUTINES ARE IN LIBXC/LIBXC-MASTER/SRC/LIBXC_MASTER.F90  ==
!         == FUNCTIONAL IDS ARE IN LIBXC/LIBXC-MASTER/SRC/LIBXC_INC.F90       ==
!         ======================================================================
          IF(+VAL(I)(:7).NE.+'XC_LDA_'.AND.+VAL(I)(:7).NE.+'XC_GGA_') THEN
            CALL ERROR$MSG('LIBXC SELECTION IS NOT YET SUPPORTED')
            CALL ERROR$CHVAL('SELECTION',TRIM(ADJUSTL(VAL(I))))
            CALL ERROR$MSG('USE SELECTION BEGINNING WITH XC_LDA_ OR XC_GGA_')
            CALL ERROR$STOP('PAWLIBXC$SETCHA')
          END IF
          J=XC_F03_FUNCTIONAL_GET_NUMBER(VAL(I))
          IF(J.LT.0) THEN
            CALL ERROR$MSG('FUNCTIONAL ID NOT RECOGNIZED OR INCOMPATIBLE')
            CALL ERROR$CHVAL('FUNCTIONAL ID',VAL(I))
            CALL ERROR$STOP('PAWLIBXC$SETCHA')
          END IF
          NAME=XC_F03_FUNCTIONAL_GET_NAME(J)
!          PRINT*,'INDEX ',J,' VAL=',TRIM(VAL(I)),' NAME=',TRIM(NAME)
          CALL XC_F03_FUNC_INIT(XC_FUNC(I),J,XC_POLARIZED)
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PAWLIBXC$SETCHA')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GGA1(VAL,EXC,DER)
      USE PAWLIBXC_MODULE, ONLY : NXC &
     &                          ,XC_FUNC
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL(5)     ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC         ! SPIN DENSITY
      REAL(8)   ,INTENT(OUT):: DER(5)      ! 
      INTEGER(4)            :: I
      REAL(8)               :: X_EXC
      REAL(8)               :: X_DER(5)
!     **************************************************************************
      CALL PAWLIBXC_INITIALIZE()
      EXC=0.D0
      DER=0.D0
      DO I=1,NXC
        CALL PAWLIBXC_GGA1(XC_FUNC(I),VAL,X_EXC,X_DER)
        EXC=EXC+X_EXC*VAL(1)
        DER=DER+X_DER
      ENDDO
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE PAWLIBXC$MGGA1(VAL,EXC,DER)
!!$      USE PAWLIBXC_MODULE, ONLY : NXC &
!!$     &                          ,XC_FUNC
!!$      IMPLICIT NONE
!!$      REAL(8)   ,INTENT(IN) :: VAL(9)     ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
!!$      REAL(8)   ,INTENT(OUT):: EXC         ! SPIN DENSITY
!!$      REAL(8)   ,INTENT(OUT):: DER(9)      ! 
!!$      INTEGER(4)            :: I
!!$      REAL(8)               :: X_EXC
!!$      REAL(8)               :: X_DER(9)
!!$!     **************************************************************************
!!$      CALL PAWLIBXC_INITIALIZE()
!!$      EXC=0.D0
!!$      DER=0.D0
!!$      DO I=1,NXC
!!$        CALL PAWLIBXC_MGGA1(XC_FUNC(I),VAL,X_EXC,X_DER)
!!$        EXC=EXC+X_EXC*VAL(1)
!!$        DER=DER+X_DER
!!$      ENDDO
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GGA2(VAL,EXC,DER,DER2)
      USE PAWLIBXC_MODULE, ONLY : NXC &
     &                           ,XC_FUNC
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL(5)     ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC         ! SPIN DENSITY
      REAL(8)   ,INTENT(OUT):: DER(5)      ! 
      REAL(8)   ,INTENT(OUT):: DER2(5,5)   ! 
      REAL(8)               :: X_EXC
      REAL(8)               :: X_DER(5),X_DER2(5,5),X_DER3(5,5,5)=0.D0
      INTEGER(4)            :: I
!     **************************************************************************
      CALL PAWLIBXC_INITIALIZE()
      EXC=0.D0
      DER=0.D0
      DER2=0.D0
      DO I=1,NXC
        CALL PAWLIBXC_GGA3(XC_FUNC(I),VAL,X_EXC,X_DER,X_DER2,X_DER3)
        EXC=EXC+X_EXC*VAL(1)
        DER=DER+X_DER
        DER2=DER2+X_DER2
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC$GGA3(VAL,EXC,DER,DER2,DER3)
      USE PAWLIBXC_MODULE, ONLY : NXC &
     &                           ,XC_FUNC
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL(5)     ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC         ! SPIN DENSITY
      REAL(8)   ,INTENT(OUT):: DER(5)      ! 
      REAL(8)   ,INTENT(OUT):: DER2(5,5)   ! 
      REAL(8)   ,INTENT(OUT):: DER3(5,5,5) ! 
      REAL(8)               :: X_EXC
      REAL(8)               :: X_DER(5),X_DER2(5,5),X_DER3(5,5,5)
      INTEGER(4)            :: I
!     **************************************************************************
      CALL PAWLIBXC_INITIALIZE()
      EXC=0.D0
      DER=0.D0
      DER2=0.D0
      DER3=0.D0
      DO I=1,NXC
        CALL PAWLIBXC_GGA3(XC_FUNC(I),VAL,X_EXC,X_DER,X_DER2,X_DER3)
        EXC=EXC+X_EXC*VAL(1)
        DER=DER+X_DER
        DER2=DER2+X_DER2
        DER3=DER3+X_DER3
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC_GGA1(XC_FUNC,VAL,EXC,DER)
!     **************************************************************************
!     ** DFT INTERFACE TO LIBXC CALL FOR THE GGA FAMILY                      **
!     **************************************************************************
      USE XC_F03_LIB_M   , ONLY : XC_F03_FUNC_T  &        !TYPE
     &                           ,XC_F03_FUNC_INFO_T &    !TYPE
     &                           ,XC_F03_LDA_EXC_VXC &    !FUNCTION
     &                           ,XC_F03_GGA_EXC_VXC &    !FUNCTION
     &                           ,XC_F03_FUNC_GET_INFO &  !FUNCTION
     &                           ,XC_F03_FUNC_INFO_GET_FAMILY & !FUNCTION
     &                           ,XC_FAMILY_LDA &         !VALUE
     &                           ,XC_FAMILY_GGA &         !VALUE
     &                           ,XC_FAMILY_HYB_GGA       !VALUE
      USE PAWLIBXC_MODULE, ONLY : MAT => MAT_GGA
      IMPLICIT NONE
      TYPE(XC_F03_FUNC_T),INTENT(IN) :: XC_FUNC
      REAL(8)   ,INTENT(IN) :: VAL(5)      ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC         ! XC ENERGY DENSITY 
      REAL(8)   ,INTENT(OUT):: DER(5)      ! FIRST DERIVARIVES OF EXC
      INTEGER(8),PARAMETER  :: NP=1
      TYPE(XC_F03_FUNC_INFO_T):: XC_INFO
      REAL(8)               :: EXCARR(1)
      REAL(8)               :: RHO(2)         
      REAL(8)               :: SIGMA(3)       
      REAL(8)               :: VRHO(2)         ! 0,1
      REAL(8)               :: VSIGMA(3)       ! 0,1,2
      INTEGER(4)            :: I1,I2,I3
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
!     **************************************************************************
      RHO(1)=0.5D0*(VAL(1)+VAL(2))                 !RHOUP
      RHO(2)=0.5D0*(VAL(1)-VAL(2))                 !RHODN
      SIGMA(1)=0.25D0*(VAL(3)+VAL(4)+2.D0*VAL(5))  !GRHOUP*GRHOUP
      SIGMA(2)=0.25D0*(VAL(3)-VAL(4))              !GRHOUP*GRHODN
      SIGMA(3)=0.25D0*(VAL(3)+VAL(4)-2.D0*VAL(5))  !GRHODN*GRHODN
!
!     ==========================================================================
!     == TEST TRANSFORMATION MATRIX                                           ==
!     ==========================================================================
      IF(TPR) THEN
        PRINT*,'SHAPE(MAT)',SHAPE(MAT)
        WRITE(*,*)'=== TRANSFORMATION MATRIX ==='
        DO I1=1,5
          WRITE(*,FMT='(5F10.5)')MAT(I1,:)
        ENDDO
        WRITE(*,FMT='("VAL0=",10F10.5)')VAL
        WRITE(*,FMT='("VAL1=",10F10.5)')RHO,SIGMA
        WRITE(*,FMT='("VAL2=",10F10.5)')MATMUL(MAT,VAL)
!        STOP
      END IF
!
!     ==========================================================================
!     == EVALUATE DENSITY FUNCTIONAL AND DERIVATIVES                          ==
!     ==========================================================================
      EXCARR=0.D0
      DER=0.D0
      XC_INFO=XC_F03_FUNC_GET_INFO(XC_FUNC)
      SELECT CASE (XC_F03_FUNC_INFO_GET_FAMILY(XC_INFO))
        CASE(XC_FAMILY_LDA)
          CALL XC_F03_LDA_EXC_VXC(XC_FUNC,NP,RHO,EXCARR,VRHO)
          VSIGMA=0.D0
        CASE(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
          CALL XC_F03_GGA_EXC_VXC(XC_FUNC,NP,RHO,SIGMA,EXCARR,VRHO,VSIGMA)
      END SELECT
      EXC=EXCARR(1)
      DER(1:2)=VRHO(:)
      DER(3:5)=VSIGMA(:)      
!
!     ==========================================================================
!     == TRANSFORM BACK                                                       ==
!     == (RHO,SIGMA) = MAT * (RHOT,RHOS,GRHOT2,GRHOS2,GRHOST)
!     ==========================================================================
      DER=MATMUL(DER,MAT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWLIBXC_GGA3(XC_FUNC,VAL,EXC,DER,DER2,DER3)
!     **************************************************************************
!     ** DFT3 INTERFACE TO LIBXC CALL FOR THE GGA FAMILY                      **
!     **************************************************************************
      USE XC_F03_LIB_M    , ONLY : XC_F03_FUNC_T &               !TYPE
     &                            ,XC_F03_FUNC_INFO_T &          !TYPE
     &                            ,XC_F03_LDA_EXC_VXC_FXC_KXC &  !FUNCTION
     &                            ,XC_F03_GGA_EXC_VXC_FXC_KXC &  !FUNCTION
     &                            ,XC_F03_FUNC_GET_INFO &        !FUNCTION
     &                            ,XC_F03_FUNC_INFO_GET_FAMILY & !FUNCTION
     &                            ,XC_FAMILY_LDA &         !VALUE
     &                            ,XC_FAMILY_GGA &         !VALUE
     &                            ,XC_FAMILY_HYB_GGA       !VALUE
      USE PAWLIBXC_MODULE, ONLY : MAT => MAT_GGA
      IMPLICIT NONE
      TYPE(XC_F03_FUNC_T),INTENT(IN) :: XC_FUNC
      REAL(8)   ,INTENT(IN) :: VAL(5)      ! (RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS)
      REAL(8)   ,INTENT(OUT):: EXC          ! XC ENERGY DENSITY 
      REAL(8)   ,INTENT(OUT):: DER(5)      ! FIRST DERIVARIVES OF EXC
      REAL(8)   ,INTENT(OUT):: DER2(5,5)   ! SECOND DERIVARIVES OF EXC
      REAL(8)   ,INTENT(OUT):: DER3(5,5,5) ! THIRD DERIVARIVES OF EXC
      INTEGER(8),PARAMETER  :: NP=1
      TYPE(XC_F03_FUNC_INFO_T):: XC_INFO
      REAL(8)               :: EXCARR(1)
      REAL(8)               :: RHO(2)         
      REAL(8)               :: SIGMA(3)       
      REAL(8)               :: VRHO(2)         ! 0,1
      REAL(8)               :: VSIGMA(3)       ! 0,1,2
      REAL(8)               :: V2RHO2(3)       ! 00,01,11
      REAL(8)               :: V2RHOSIGMA(6)   ! 00,01,02,10,11,12
      REAL(8)               :: V2SIGMA2(6)     ! 00,01,02,11,12,22
      REAL(8)               :: V3RHO3(4)       ! 000,001,011,111
      REAL(8)               :: V3RHO2SIGMA(9) 
    !                          ! 000,001,0002,010,011,012,110,111,112
      REAL(8)               :: V3RHOSIGMA2(12)
     !                         ! 000,001,002,011,012,022,100,101,102,111,112,122
      REAL(8)               :: V3SIGMA3(10)
                               ! 000,001,002,011,012,022,111,112,122,222
      REAL(8)               :: V4RHO4(5)
      REAL(8)               :: V4RHO3SIGMA(12)
      REAL(8)               :: V4RHO2SIGMA2(15)
      REAL(8)               :: V4RHOSIGMA3(20)
      REAL(8)               :: V4SIGMA4(15)
      INTEGER(4)            :: I1,I2,I3
      REAL(8)               :: VAL1(5)
!     **************************************************************************
      RHO(1)=0.5D0*(VAL(1)+VAL(2))                 !RHOUP
      RHO(2)=0.5D0*(VAL(1)-VAL(2))                 !RHODN
      SIGMA(1)=0.25D0*(VAL(3)+VAL(4)+2.D0*VAL(5))  !GRHOUP*GRHOUP
      SIGMA(2)=0.25D0*(VAL(3)-VAL(4))              !GRHOUP*GRHODN
      SIGMA(3)=0.25D0*(VAL(3)+VAL(4)-2.D0*VAL(5))  !GRHODN*GRHODN
!
!     ==========================================================================
!     == EVALUATE DENSITY FUNCTIONAL AND DERIVATIVES                          ==
!     ==========================================================================
!PRINT*,'BEFORE XC_F03_GGA'
!PRINT*,'BEFORE XC_F03_GGA_EXC_VXC_FXC_KXC'
! F03_GGA_EXC_VXC_FXC_KXC?

      XC_INFO=XC_F03_FUNC_GET_INFO(XC_FUNC)
      SELECT CASE (XC_F03_FUNC_INFO_GET_FAMILY(XC_INFO))
        CASE(XC_FAMILY_LDA)
          CALL XC_F03_LDA_EXC_VXC_FXC_KXC(XC_FUNC,NP,RHO &
     &           ,EXCARR,VRHO,V2RHO2,V3RHO3)
          VSIGMA=0.D0
          V2RHOSIGMA=0.D0
          V2SIGMA2=0.D0
          V3RHO2SIGMA=0.D0
          V3RHOSIGMA2=0.D0
          V3SIGMA3=0.D0
        CASE(XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
          CALL XC_F03_GGA_EXC_VXC_FXC_KXC(XC_FUNC,NP,RHO,SIGMA &
     &           ,EXCARR,VRHO,VSIGMA &
     &           ,V2RHO2,V2RHOSIGMA ,V2SIGMA2 &
     &           ,V3RHO3,V3RHO2SIGMA,V3RHOSIGMA2, V3SIGMA3)
      END SELECT

      EXC=EXCARR(1)
      DER=0.D0
      DER2=0.D0
      DER3=0.D0
      DER(1:2)=VRHO(:)
      DER(3:5)=VSIGMA(:)      
      DER2(1,1)=V2RHO2(1)
      DER2(1,2)=V2RHO2(2)
      DER2(2,2)=V2RHO2(3)
      DER2(1,3)=V2RHOSIGMA(1)
      DER2(1,4)=V2RHOSIGMA(2)
      DER2(1,5)=V2RHOSIGMA(3)
      DER2(2,3)=V2RHOSIGMA(4)
      DER2(2,4)=V2RHOSIGMA(5)
      DER2(2,5)=V2RHOSIGMA(6)
      DER2(3,3)=V2SIGMA2(1)
      DER2(3,4)=V2SIGMA2(2)
      DER2(3,5)=V2SIGMA2(3)
      DER2(4,4)=V2SIGMA2(4)
      DER2(4,5)=V2SIGMA2(5)
      DER2(5,5)=V2SIGMA2(6)

      DER3(1,1,1)=V3RHO3(1)
      DER3(1,1,2)=V3RHO3(2)
      DER3(1,2,2)=V3RHO3(3)
      DER3(2,2,2)=V3RHO3(4)
      DER3(1,1,3)=V3RHO2SIGMA(1)
      DER3(1,1,4)=V3RHO2SIGMA(2)
      DER3(1,1,5)=V3RHO2SIGMA(3)
      DER3(1,2,3)=V3RHO2SIGMA(4)
      DER3(1,2,4)=V3RHO2SIGMA(5)
      DER3(1,2,5)=V3RHO2SIGMA(6)
      DER3(2,2,3)=V3RHO2SIGMA(7)
      DER3(2,2,4)=V3RHO2SIGMA(8)
      DER3(2,2,5)=V3RHO2SIGMA(9)
      DER3(1,3,3)=V3RHOSIGMA2(1)
      DER3(1,3,4)=V3RHOSIGMA2(2)
      DER3(1,3,5)=V3RHOSIGMA2(3)
      DER3(1,4,4)=V3RHOSIGMA2(4)
      DER3(1,4,5)=V3RHOSIGMA2(5)
      DER3(1,5,5)=V3RHOSIGMA2(6)
      DER3(2,3,3)=V3RHOSIGMA2(7)
      DER3(2,3,4)=V3RHOSIGMA2(8)
      DER3(2,3,5)=V3RHOSIGMA2(9)
      DER3(2,4,4)=V3RHOSIGMA2(10)
      DER3(2,4,5)=V3RHOSIGMA2(11)
      DER3(2,5,5)=V3RHOSIGMA2(12)
      DER3(3,3,3)=V3SIGMA3(1)
      DER3(3,3,4)=V3SIGMA3(2)
      DER3(3,3,5)=V3SIGMA3(3)
      DER3(3,4,4)=V3SIGMA3(4)
      DER3(3,4,5)=V3SIGMA3(5)
      DER3(3,5,5)=V3SIGMA3(6)
      DER3(4,4,4)=V3SIGMA3(7)
      DER3(4,4,5)=V3SIGMA3(8)
      DER3(4,5,5)=V3SIGMA3(9)
      DER3(5,5,5)=V3SIGMA3(10)
!
!     ==========================================================================
!     == COMPLETE MATRIX ELEMENTS                                             ==
!     ==========================================================================
      DO I2=1,5
        DO I1=1,I2-1
          DER2(I2,I1)=DER2(I1,I2)
        ENDDO
      ENDDO
!
      DO I3=1,5
        DO I2=1,I3-1
          DO I1=1,I2-1
            DER3(I2,I3,I1)=DER3(I1,I2,I3)
            DER3(I3,I2,I1)=DER3(I1,I2,I3)
            DER3(I1,I3,I2)=DER3(I1,I2,I3)
            DER3(I3,I1,I2)=DER3(I1,I2,I3)
            DER3(I2,I1,I3)=DER3(I1,I2,I3)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == TRANSFORM BACK                                                       ==
!     ==========================================================================
      DER=MATMUL(DER,MAT)
      DER2=MATMUL(DER2,MAT)
      DO I1=1,5
        DER2(:,I1)=MATMUL(DER2(:,I1),MAT)
      ENDDO
      DO I1=1,5
        DO I2=1,5
          DER3(I1,I2,:)=MATMUL(DER3(I1,I2,:),MAT)
        ENDDO
      ENDDO
      DO I1=1,5
        DO I2=1,5
          DER3(:,I1,I2)=MATMUL(DER3(:,I1,I2),MAT)
        ENDDO
      ENDDO
      DO I1=1,5
        DO I2=1,5
          DER3(I1,:,I2)=MATMUL(DER3(I1,:,I2),MAT)
        ENDDO
      ENDDO
      RETURN
      END
#ENDIF
