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
!     **                                                                      **
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
      CHARACTER(LEN=120)           :: S1, S2
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
!         __GET INFO ON REFERENCE NUMBER I; CAN'T REUSE THE LOOP VARIABLE
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
        PRINT*,'ID NOT RECOGNIZED'
        PRINT*,'ID=',ID
        STOP 'IN DFT_LIBXC1_ADDFUNC'
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
      REAL(8)               :: X_DER(5),X_DER2(5,5),X_DER3(5,5,5)=0.d0
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
        DO I1=1,I1-1
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
