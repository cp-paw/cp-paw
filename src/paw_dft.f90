!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE DFT_MODULE
!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**  NAME:  DFT                                                               **
!**                                                                           **
!**  PURPOSE :  CALCULATES THE EXCHANGE CORRELATION ENERGY                    **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**                                                                           **
!**    DFT$SETI4(ID,VAL)                                                      **
!**    DFT$GETI4(ID,VAL)                                                      **
!**    DFT$SETL4(ID,VAL)                                                      **
!**    DFT$GETL4(ID,VAL)                                                      **
!**    DFT$REPORT(NFIL)                                                       **
!**    DFT(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST,EXC,VT,VS,VGT2,VGS2,VGTS)           **
!**    DFT2(X,E,DE,D2E)                                                       **
!**    DFT3(X,E,DE,D2E,D3E)                                                   **
!**                                                                           **
!**  USAGE:                                                                   **
!**    1) SELECT A FUNCTIONAL: DFT$SETI4('TYPE',ITYPE)                        **
!**    2) CHOSE SPIN (UN) RESTRICTED: DFT$SETL4('SPIN',TSPIN)                 **
!**         (DEFAULT IS TSPIN=.TRUE.)                                         **
!**    3A) EVALUATE WITH GRADIENTS :                                          **
!**        DFT(RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS,EXC &                           **
!**    3B) EVALUATE WITH FIRST AND SECOND DERIVATIVES                         **
!**            DFT2(VAL,EXC,DEXC,D2EXC)                                       **
!**    3C) EVALUATE WITH FIRST, SECOND AND THIRD DERIVATIVES                  **
!**            DFT3(VAL,EXC,DEXC,D2EXC,D3EXC)                                 **
!**                                                                           **
!**  NOTES:                                                                   **
!**    THE PARAMETER RHOTMIN MODIFIES THE FUNCTIONAL IF THE                   **
!**    DENSITIES BECOME TOO SMALL, IN ORDER TO AVOID INSTABILITIES            **
!**    AS SOME QUANTITIES GO TO INFINITY.                                     **
!**                                                                           **
!**  DEPENDENCIES:                                                            **
!**    EXCHANGE                                                               **
!**    PERDEWZUNGER                                                           **
!**    PERDEWWANG91L                                                          **
!**    PERDEWWANG91G  (NO THIRD DERIVATIVES)                                  **
!**    PBE96                                                                  **
!**    RPBE                                                                   **
!**    PERDEW    (NO THIRD DERIVATIVES)                                       **
!**                                                                           **
!**  SPIN INVERSION SYMMETRY:                                                 **
!**          RHOS <-> -RHOS AND GRHOS <-> -GRHOS                              **
!**  SPACE INVERSION SYMMETRY                                                 **
!**          GRHOT <-> -GRHOT AND GRHOS <-> -GRHOS                            **
!**                                                                           **
!*******************************************************************************
!*******************************************************************************
!     ==  SWITCHES =============================================================
      INTEGER(4)           :: IT=0   ! FUNCTIONAL SELECTOR SEE DFT$REPORT
      LOGICAL(4),PARAMETER :: TSAFE  =.TRUE.  ! NO DENSITY-GRID INTERPOLATION
      LOGICAL(4)           :: TSPIN  =.TRUE.  ! CAN BE SET W.O. INITIALIZATION
      LOGICAL(4)           :: TGRA   =.FALSE. ! GRADIENTS ARE USED (OUTPUT)
      LOGICAL(4)           :: TSIC   =.FALSE. ! NOT USED
      LOGICAL(4)           :: TPLUSU =.FALSE. ! NOT USED
!     == SWITCHES FOR INTERNAL BOOKKEEPING =====================================
      REAL(8)   ,PARAMETER :: RHOTMIN=1.D-6   ! MINIMUM DENSITY
      LOGICAL(4)           :: TINI   =.FALSE. ! INITIALIZATION DONE
      LOGICAL(4)           :: TX     =.FALSE. ! USE EXCHANGE
      LOGICAL(4)           :: TCORRELATION =.TRUE. ! SWITCH FOR CORRELATION
      LOGICAL(4)           :: TPZ    =.FALSE. ! USE PERDEW ZUNGER
      LOGICAL(4)           :: TBH    =.FALSE. ! USE BARTH HEDIN
      LOGICAL(4)           :: TPERDEW=.FALSE. ! USE PERDEW GC
      LOGICAL(4)           :: TPW91L =.FALSE. ! USE PERDEW WANG LOCAL
      LOGICAL(4)           :: TPW91G =.FALSE. ! USE PERDEW WANG91 GC
      LOGICAL(4)           :: TPBE96 =.FALSE. ! USE PERDEW BURKE ERNZERHOF GC
      LOGICAL(4)           :: TLYP88 =.FALSE. ! USE LEE-YANG-PARR 88 CORRELATION
      LOGICAL(4)           :: TLIBXC =.FALSE. ! USE LIBXC INTERFACE
      REAL(8)              :: SCALEX=1.D0     ! SCALES EXCHANGE CONTRIBUTION
!                                             ! USED FOR HYBRID FUNCTIONALS
      INTEGER(4),PARAMETER :: NDESCRIPTION=5
      CHARACTER(128)       :: DESCRIPTION(NDESCRIPTION)
      END MODULE DFT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT_INITIALIZE
!     **************************************************************************
!     **                                                                      **
!     **  DEFINES THE FUNCTIONAL SELECTORS                                    **
!     **                                                                      **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      LOGICAL(4)     :: TGRATARGET=.FALSE.
!     **************************************************************************
      IF(TINI) RETURN
!
!     ==========================================================================
!     ==  RESET FUNCTIONAL SELECTORS                                          ==
!     ==========================================================================
      TPZ    =.FALSE.  ! USE PERDEW ZUNGER
      TBH    =.FALSE.  ! USE BARTH HEDIN
      TX     =.FALSE.  ! USE EXCHANGE MODULE
      TPERDEW=.FALSE.  ! USE PERDEW GC
      TPW91L =.FALSE.  ! USE PERDEW WANG LOCAL
      TPW91G =.FALSE.  ! USE PERDEW WANG91 GC
      TPBE96 =.FALSE.  ! USE PERDEW BURKE ERNZERHOF GC
      TLYP88 =.FALSE.  ! USE LEE YANG PARR GC
      DESCRIPTION(:)=' '
      DESCRIPTION(1)='NO FUNCTIONAL DESCRIPTION GIVEN'
!
!     ==========================================================================
!     ==  HANDLE LIBXC                                                        ==
!     ==========================================================================
      IF(TLIBXC) THEN
        RETURN
      END IF
!
!     ==========================================================================
!     ==  NOW SET THE FUNCTIONAL SELECTORS                                    ==
!     ==========================================================================
!
!     ==========================================================================
!     == TYPE 1:  LSD (PERDEW ZUNGER PARAMETERIZATION)                        ==
!     ==========================================================================
      IF(IT.EQ.0) THEN
        TGRATARGET=.FALSE.
        DESCRIPTION(1)='NO EXCHANGE OR CORRELATION'
      ELSE IF(IT.EQ.1) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',1)  ! SELECT LOCAL EXCHANGE
        TPZ=.TRUE.        
        TGRATARGET=.FALSE.
        DESCRIPTION(1)='PERDEW-ZUNGER PARAMETRIZATION: ' &
     &                          //'(PHYS.REV.B23,P5048(1981)'
!
!     ==========================================================================
!     == TYPE 2:  LSD (PERDEW WANG PARAMETERIZATION)                          ==
!     ==========================================================================
      ELSE IF(IT.EQ.2) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',1)  ! SELECT LOCAL EXCHANGE
        TPW91L=.TRUE.        
        TGRATARGET=.FALSE.
        DESCRIPTION(1)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
!
!     ==========================================================================
!     == TYPE 3:  LOCAL HARTREE FOCK                                          ==
!     ==========================================================================
      ELSE IF(IT.EQ.3) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',1)   ! SELECT LOCAL EXCHANGE
        TGRATARGET=.FALSE.
        DESCRIPTION(1)='LOCAL EXCHANGE (X-ALPHA=2/3)'
!
!     ==========================================================================
!     == TYPE 4:  X_ALPHA WITH ALPHA=0.7 (=1.05 TIMES LOCAL HF)               ==
!     ==========================================================================
      ELSE IF(IT.EQ.4) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',12)  ! SELECT LOCAL EXCHANGE TIMES 1.5*0.7
        TGRATARGET=.FALSE.
        DESCRIPTION(1)='X-ALPHA = 0.7'
!
!     ==========================================================================
!     == TYPE 6:  LOCAL HARTREE FOCK + BECKE GC                               ==
!     ==========================================================================
      ELSE IF (IT.EQ.6) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',2)
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='LOCAL HARTREE FOCK'
        DESCRIPTION(2)='BECKE GRADIENT CORRECTION FOR EXCHANGE ' &
     &                          //'(J.CHEM.PHYS.96,P2155(1992))'
!
!     ==========================================================================
!     == TYPE 7:  LSD + BECKE GC                                              ==
!     ==========================================================================
      ELSE IF (IT.EQ.7) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',2)
        TPW91L=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='BECKE GRADIENT CORRECTION FOR EXCHANGE ' &
     &                          //'(J.CHEM.PHYS.96,P2155(1992))'
!
!     ==========================================================================
!     == TYPE 71:  LSD(PZ) + BECKE GC                                         ==
!     ==========================================================================
      ELSE IF (IT.EQ.71) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',2)
        TPZ=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW-ZUNGER PARAMETRIZATION: ' &
     &                          //'(PHYS.REV.B23,P5048(1981)'
        DESCRIPTION(2)='BECKE GRADIENT CORRECTION FOR EXCHANGE ' &
     &                          //'(J.CHEM.PHYS.96,P2155(1992))'
!
!     ==========================================================================
!     == TYPE 8:  LSD + BECKE GC + PERDEW86                                   ==
!     ==========================================================================
      ELSE IF (IT.EQ.8) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',2)
        TPW91L=.TRUE.
        TPERDEW=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='BECKE GRADIENT CORRECTION FOR EXCHANGE ' &
     &                          //'(J.CHEM.PHYS.96,P2155(1992))'
        DESCRIPTION(3)='PERDEW GRADIENT CORRECTION FOR CORRELATION ' &
     &                          //'(PHYS.REV.B33,P8822(1986))'
!
!     ==========================================================================
!     == TYPE 81:  LSD + BECKE GC + PERDEW86                                  ==
!     ==========================================================================
      ELSE IF (IT.EQ.81) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',2)
        TPZ=.TRUE.
        TPERDEW=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW-ZUNGER PARAMETRIZATION: ' &
     &                         //'(PHYS.REV.B23,P5048(1981)'
        DESCRIPTION(2)='BECKE GRADIENT CORRECTION FOR EXCHANGE ' &
     &                          //'(J.CHEM.PHYS.96,P2155(1992))'
        DESCRIPTION(3)='PERDEW GRADIENT CORRECTION FOR CORRELATION ' &
     &                          //'(PHYS.REV.B33,P8822(1986))'
!
!     ==========================================================================
!     == TYPE 9:  LSD + BECKE GC + PBEC                                       ==
!     ==========================================================================
      ELSE IF (IT.EQ.9) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',2)
        TPBE96=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='BECKE GRADIENT CORRECTION FOR EXCHANGE ' &
     &                          //'(J.CHEM.PHYS.96,P2155(1992))'
        DESCRIPTION(3)='PERDEW-BURKE-ERNZERHOF GGA FOR CORRELATION ' &
     &                          //'(PHYS.REV.B 46, 6671 (1992-I))'
        DESCRIPTION(4)='(PBE96 SHOULD BE SIMILAR TO PERDEW WANG GGA 91) ' &
     &                          //'(PHYS.REV.B 46, 6671 (1992-I)))'
!
!     ==========================================================================
!     == TYPE 10:  LSD + PBE96-GC                                             ==
!     ==========================================================================
      ELSE IF (IT.EQ.10.OR.IT.EQ.-10) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',3)
        TPBE96=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='PERDEW-BURKE-ERNZERHOF GGA FOR EXCHANGE AND' &
     &                   //' CORRELATION (PHYS.REV.LETT 77, 3865 (1996))'
        IF(IT.EQ.-10) THEN
          DESCRIPTION(3)= &
     &     'CAUTION!!! PARAMETER TFAC DIFFERS FROM REFERENCE ABOVE'
        END IF
!
!     ==========================================================================
!     == TYPE 11:  LSD + RPBE-GC                                              ==
!     ==========================================================================
      ELSE IF (IT.EQ.11) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',4)
        TPBE96=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW-WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='PERDEW-BURKE-ERNZERHOF GGA FOR EXCHANGE AND' &
     &                //' CORRELATION (PHYS.REV.LETT.77, 3865 (1992))'
        DESCRIPTION(3)='HAMMER-HANSEN-NORSKOV RPBE-GGA FOR EXCHANGE AND' &
     &                //' CORRELATION (PHYS.REV.B 59, 7413 (1999))'
!
!     ==========================================================================
!     == TYPE 12:  LSD + PW91-GC                                              ==
!     ==========================================================================
      ELSE IF (IT.EQ.12) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',3)
        TPW91G=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW-WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='PERDEW-WANG 91 GGA FOR EXCHANGE AND CORRELATION ' &
     &                          //'(PHYS.REV.B 46, 6671 (1992))'
!
!     ==========================================================================
!     == TYPE 13:  LSD + REVPBE-GC                                            ==
!     ==========================================================================
      ELSE IF (IT.EQ.13) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',5)
        TPBE96=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='REVPBE OF Y. ZHANG AND W. YANG GGA FOR EXCHANGE AND' &
     &                //' CORRELATION (PHYS.REV.LETT 80, 890 (1998))'
!
!     ==========================================================================
!     == TYPE 5001:  LOCAL HF+ BECKE XC-GC+ LYP88 CORRELATION                 ==
!     ==========================================================================
      ELSE IF (IT.EQ.5001) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',2)
        TLYP88=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='BECKE GRADIENT CORRECTION FOR EXCHANGE ' &
     &                          //'(J.CHEM.PHYS.96,P2155(1992))'
        DESCRIPTION(2)='LEE-YANG-PARR-88 CORRELATION' &
     &                          //'(PHYS.REV.B 37, 785 (1988-I))'
!
!     ==========================================================================
!     == TYPE 5002:  LSD + PBE96-GC WITH SCALED EXCHANGE                      ==
!     ==========================================================================
      ELSE IF (IT.EQ.5002) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',3)
        SCALEX=0.75D0
        TPBE96=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
        DESCRIPTION(2)='PERDEW-BURKE-ERNZERHOF GGA FOR EXCHANGE AND CORRELATION ' &
     &                          //'(PHYS.REV.LETT 77, 3865 (1996))'
        DESCRIPTION(3)='EXCHANGE SCALED BY 0.75 FOR USE IN THE PBE0 HYBRID FUNCTIONAL'
!
!     ==========================================================================
!     == TYPE 10001:  CORRELATION (PERDEW ZUNGER PARAMETERIZATION)            ==
!     ==========================================================================
      ELSE IF(IT.EQ.10001) THEN
        TPZ=.TRUE.        
        TGRATARGET=.FALSE.
        DESCRIPTION(1)='NO EXCHANGE! (TEST OPTION) '
        DESCRIPTION(2)='PERDEW-ZUNGER PARAMETRIZATION: ' &
     &                          //'(PHYS.REV.B23,P5048(1981)'
!
!     ==========================================================================
!     == TYPE 10002:  CORRELATION (PERDEW WANG PARAMETERIZATION)              ==
!     ==========================================================================
      ELSE IF(IT.EQ.10002) THEN
        TPW91L=.TRUE.        
        TGRATARGET=.FALSE.
        DESCRIPTION(1)='NO EXCHANGE! (TEST OPTION) '
        DESCRIPTION(2)='PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
!
!     ==========================================================================
!     == TYPE 10003: PERDEW-86 GRADIENT CORRECTION FOR CORRELATION            ==
!     ==========================================================================
      ELSE IF(IT.EQ.10003) THEN
        TPERDEW=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='NO EXCHANGE AND NO LOCAL CORRELATION! (TEST OPTION) '
        DESCRIPTION(2)='PERDEW GRADIENT CORRECTION FOR CORRELATION ' &
     &                          //'(PHYS.REV.B33,P8822(1986))'
!
!     ==========================================================================
!     == TYPE 10004:CORRELATION (RPBE-TYPE11 OR PBE-TYPE10 OR REVPBE-TYPE-13) ==
!     ==========================================================================
      ELSE IF (IT.EQ.10004) THEN
        TPBE96=.TRUE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='NO EXCHANGE !  (TEST OPTION) ' 
        DESCRIPTION(2)='PERDEW-BURKE-ERNZERHOF GGA FOR EXCHANGE AND' &
     &               //' CORRELATION (PHYS.REV.B 46, 6671 (1992-I))'
        DESCRIPTION(3)='HAMMER-HANSEN-NORSKOV RPBE-GGA FOR EXCHANGE AND' &
     &               //' CORRELATION (PHYS.REV.B 59, 7413 (1999))'
!
!     ==========================================================================
!     == TYPE 10005:  PBE EXCHANGE ONLY                                       ==
!     ==========================================================================
      ELSE IF (IT.EQ.10005) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',3)
        TPBE96=.FALSE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='PERDEW-BURKE-ERNZERHOF GGA FOR EXCHANGE ONLY ' &
     &                          //'(PHYS.REV.LETT 77, 3865 (1996))'
!
!     ==========================================================================
!     == TYPE 10006:  CORRELATION - PW91-GC(PERDEW WANG PARAMETERIZATION)     ==
!     ==========================================================================
      ELSE IF(IT.EQ.10006) THEN
        TPW91G=.TRUE.        
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='NO EXCHANGE! (TEST OPTION) '
        DESCRIPTION(2)='PERDEW-WANG PARAMETERIZATION OF LOCAL CORRELATION ' &
     &                          //'(PHYS.REV.B 45, 13244 (1992-I))'
!   
!     ==========================================================================
!     == TYPE 10007:  RPBE EXCHANGE ONLY                                      ==
!     ==========================================================================
      ELSE IF (IT.EQ.10007) THEN
        TX=.TRUE.
        CALL EXCHANGE$SETI4('TYPE',4)
        TPBE96=.FALSE.
        TGRATARGET=.TRUE.
        DESCRIPTION(1)='HAMMER-HANSEN-NORSKOV RPBE-GGA FOR EXCHANGE AND ' &
     &                          //'CORRELATION (PHYS.REV.B 59, 7413 (1999))'
!
!     ==========================================================================
!     == TYPE 20:  LIBXC-GGA INTERFACE                                    ==
!     ==========================================================================
      ELSE IF (IT.EQ.20) THEN
        TLIBXC=.TRUE.
        TGRATARGET=.FALSE.
        CALL PAWLIBXC$SETCH('FUNCTIONAL','VWN')
!       == THE LIBXC INTERFACE MUST USE NOT USE OUR OWN IMPLEMENTATIONS FOR   ==
!       == EXCHANGE (HENCE TX=.FALSE), BUT THE EXCHANGE TERM                  ==
!       == IS INTEGRATED WITH THE CORRELATION TERM                            ==
        TX=.FALSE.

        DESCRIPTION(1)='LIBXC INTERFACE CURRENTLY VWN' &
     &               //'(---------------------)'
!
!     ==========================================================================
!     == TYPE 111111:  LIBXC-GGA INTERFACE                                    ==
!     ==========================================================================
      ELSE IF (IT.EQ.111111) THEN
        TLIBXC=.TRUE.
        TGRATARGET=.TRUE.
        CALL PAWLIBXC$SETCH('FUNCTIONAL','PBE')
!       == THE LIBXC INTERFACE MUST USE NOT USE OUR OWN IMPLEMENTATIONS FOR   ==
!       == (HENCE TX=.FALSE), BUT THE EXCHANGE TERM OF THE SCAN FUNCTIONAL    ==
!       == IS INTEGRATED WITH THE CORRELATION TERM                            ==
        TX=.FALSE.

        DESCRIPTION(1)='LIBXC INTERFACE CURRENTLY PBE' &
     &               //'(---------------------)'
!
!     ==========================================================================
!     == ILLEGAL SELECTION                                                    ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('CHOICE OF FUNCTIONAL NOT DEFINED')
        CALL ERROR$I4VAL('IT',IT)
        CALL ERROR$STOP('DFT_INITIALIZE')
      END IF
!     
!     ==========================================================================
!     == SET SWITCH FOR INITALIALIZATION                                      ==
!     ==========================================================================
      IF(TGRA.NEQV.TGRATARGET) THEN
        CALL ERROR$MSG('FLAG FOR GRADIENTS IS NOT SET PROPERLY')
        CALL ERROR$I4VAL('IT',IT)
        CALL ERROR$L4VAL('TGRA',TGRA)
        CALL ERROR$L4VAL('TGRATARGET',TGRATARGET)
        CALL ERROR$STOP('DFT_INITIALIZE')
      END IF
!     
!     ==========================================================================
!     == SET SWITCH FOR INITALIALIZATION                                      ==
!     ==========================================================================
      TINI=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$REPORT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **  WRITES A DESCRIPTION OF THE CURRENT (FUNCTIONAL) SETTING            **
!     **  TO A FILE                                                           **
!     **                                                                      **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: I
!     **************************************************************************
      IF(.NOT.TINI) CALL DFT_INITIALIZE
      CALL REPORT$TITLE(NFIL,'DENSITY FUNCTIONAL')
!
!     ==========================================================================
!     == DELEGATE REPORTING TO PAWLIBXC OBJECT                                ==
!     ==========================================================================
      IF(TLIBXC) THEN
        CALL PAWLIBXC$REPORT(NFIL)
        RETURN
      END IF
!
!     ==========================================================================
!     == DO REPORTING FOR CPPAW INTRINSIC FUNCTIONAL IMPLEMENTATIONS          ==
!     ==========================================================================
      IF(TSPIN) THEN
        CALL REPORT$STRING(NFIL,'SPIN POLARIZED (LSD) FUNCTIONAL')
      ELSE
        CALL REPORT$STRING(NFIL,'NON-SPIN POLARIZED (LDA) FUNCTIONAL')
      END IF        
      IF(.NOT.TCORRELATION) THEN
        CALL REPORT$STRING(NFIL,'CORRELATIONS ARE SWITCHED OFF')
      END IF
      IF(SCALEX.NE.1.D0) THEN
         CALL REPORT$R8VAL(NFIL,'EXCHANGE CONTRIBUTION SCALED BY FACTOR ' &
     &                         ,SCALEX,' ')
      END IF
      DO I=1,NDESCRIPTION
        IF(LEN(TRIM(DESCRIPTION(I))).EQ.0) CYCLE
        CALL REPORT$STRING(NFIL,TRIM(DESCRIPTION(I)))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$SETI4(ID,VALUE)
!     **************************************************************************
!     ** PASS INTEGER VALUE TO DFT OBJECT                                     **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VALUE
      LOGICAL(4)              :: TCHK
!     **************************************************************************
!
!     ==========================================================================
!     == SELECT DENSITY FUNCTIONAL TYPE PER INTEGER ID                        ==
!     ==========================================================================
      IF(ID.EQ.'TYPE') THEN
        IT=VALUE
        TGRA=.FALSE.
        TCHK=.FALSE.
        IF(IT.EQ.0)  THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
        IF(IT.EQ.1)  THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
        IF(IT.EQ.2)  THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
        IF(IT.EQ.3)  THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
        IF(IT.EQ.4)  THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
        IF(IT.EQ.6)  THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.7)  THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.71) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.8)  THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.81) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.9)  THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.10) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.-10) THEN ;TCHK=.TRUE. ;TGRA=.TRUE. 
          CALL PBE$SETL4('BUGFIX1',.FALSE.)
        END IF
!       == 111111 LIBXC-PBE   ==============================================
        IF(IT.EQ.20)THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
!       == 111111 LIBXC-PBE   ==============================================
        IF(IT.EQ.111111)THEN ;TCHK=.TRUE. ;TGRA=.TRUE. ;END IF
        IF(IT.EQ.11) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.12) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.13) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.5001) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.5002) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.10001) THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
        IF(IT.EQ.10002) THEN ;TCHK=.TRUE. ;TGRA=.FALSE. ;END IF
        IF(IT.EQ.10003) THEN ;TCHK=.TRUE. ;TGRA=.TRUE.  ;END IF
        IF(IT.EQ.10004) THEN ;TCHK=.TRUE. ;TGRA=.TRUE. ;END IF
        IF(IT.EQ.10005) THEN ;TCHK=.TRUE. ;TGRA=.TRUE. ;END IF
        IF(IT.EQ.10006) THEN ;TCHK=.TRUE. ;TGRA=.TRUE. ;END IF
        IF(IT.EQ.10007) THEN ;TCHK=.TRUE. ;TGRA=.TRUE. ;END IF
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('DFT FUNCTIONAL SELECTION INVALID')
          CALL ERROR$I4VAL('TYPE',VALUE)
          CALL ERROR$STOP('DFT$SETI4')
        END IF
        TINI=.FALSE.
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DFT$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$GETI4(ID,VALUE)
!     **************************************************************************
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VALUE
!     **************************************************************************
      IF(ID.EQ.'TYPE') THEN
        VALUE=IT
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DFT$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$SETCH(ID,VAL)
!     **************************************************************************
!     ** PASS STRING VALUE TO OBJECT                                          **
!     **                                                                      **
!     ** E.G. SELECT DENSITY FUNCTIONAL BY NAME                               **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
      LOGICAL(4)              :: TCHK
!     **************************************************************************
!
!     ==========================================================================
!     == SET FUNCTIONAL TYPE (INTERFACE FOR SETI4)                            ==
!     ==========================================================================
      IF(ID.EQ.'TYPE') THEN   
        TCHK=.TRUE.
        IF(VAL.EQ.'XC_') THEN
!         == CAUTION: PASSING SCALAR (VAL) INTO ARRAY ARGUMENT. IS THIS LEGAL? =
          CALL DFT$SETCHA(ID,1,VAL)
        ELSE IF(VAL.EQ.'PZ') THEN
          CALL  DFT$SETI4('TYPE',1)
        ELSE IF(VAL.EQ.'PERDEW-WANG:LOCAL') THEN
          CALL  DFT$SETI4('TYPE',2)
        ELSE IF(VAL.EQ.'HF:LOCAL') THEN
          CALL  DFT$SETI4('TYPE',3)
        ELSE IF(VAL.EQ.'PBE') THEN
          CALL  DFT$SETI4('TYPE',10)
        ELSE IF(VAL.EQ.'RPBE') THEN
          CALL  DFT$SETI4('TYPE',11)
        ELSE IF(VAL.EQ.'PW91-GGA') THEN
          CALL  DFT$SETI4('TYPE',12)
        ELSE IF(VAL.EQ.'BLYP') THEN
          CALL  DFT$SETI4('TYPE',5001)
        ELSE IF(VAL.EQ.'LIBXC-PBE') THEN
          CALL  DFT$SETI4('TYPE',111111)
        ELSE
          TCHK=.FALSE.
        END IF
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('INVALID DENSITY FUNCTIONAL TYPE SELECTION')
          CALL ERROR$CHVAL('DENSITY FUNCTIONAL TYPE',VAL)
          CALL ERROR$STOP('DFT$SETCH')
        END IF
        TINI=.FALSE.
!
!     ==========================================================================
!     == UNRECOGNIZED ID                                                      ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DFT$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$SETCHA(ID,LEN,VAL)
!     **************************************************************************
!     ** PASS STRING VALUES TO OBJECT                                         **
!     **                                                                      **
!     ** E.G. SELECT DENSITY FUNCTIONAL BY NAME                               **
!     **************************************************************************
      USE DFT_MODULE, ONLY : TLIBXC &
     &                      ,TGRA
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(*),INTENT(IN) :: VAL(LEN)
      LOGICAL(4)              :: TCHK
      INTEGER(4)              :: I
!     **************************************************************************
      IF(ID.EQ.'TYPE') THEN   
        TCHK=.TRUE.
        TLIBXC=(VAL(1)(1:3).EQ.'XC_')
        IF(TLIBXC) THEN
!         == THIS IS FOR THE FUNCTIONALS IN LIBXC ==============================
          CALL PAWLIBXC$SETCHA('FUNCTIONAL',LEN,VAL)
          CALL PAWLIBXC$GETL4('GRADIENT',TGRA)
        ELSE
!         == THIS IS FOR THE CPPAW INTRINSIC FUNCTIONALS =======================
          IF(LEN.EQ.1) THEN
            CALL DFT$SETCH(ID,VAL(1))
          ELSE
            CALL ERROR$MSG('UNKNOWN COMBINATION OF FUNCTIONAL ID AND LEN')
            DO I=1,LEN
              CALL ERROR$CHVAL('FUNCTIONAL ID',VAL(I))
            ENDDO
            CALL ERROR$STOP('DFT$SETCHA')
          END IF
        END IF  
!
!     ==========================================================================
!     == UNRECOGNIZED ID                                                      ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DFT$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$SETL4(ID,VAL)
!     **************************************************************************
!     ** SET LOGICAL SCALAR VARIABLE INSIDE THIS OBJECT                       **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'SPIN') THEN
        TSPIN=VAL
      ELSE IF(ID.EQ.'XCONLY') THEN
        TCORRELATION=.NOT.VAL
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DFT$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$GETL4(ID,VAL)
!     **************************************************************************
!     ** REQUEST LOGICAL SCALAR VARIABLE FROM THIS OBJECT                     **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'GC') THEN
        VAL=TGRA
      ELSE IF(ID.EQ.'SPIN') THEN
        VAL=TSPIN    ! SPIN-POLARIZATION CONSIDERED
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DFT$GETL4')
      END IF
      RETURN
      END
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$GRADIENTSWITCH(TGRA_)
!     **************************************************************************
!     == SET GRADIENT CORRECTION SWITCH                                       ==
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(OUT) :: TGRA_
!     **************************************************************************
      CALL ERROR$MSG('THIS SUBROUTINE IS MARKED FOR DELETION')
      CALL ERROR$MSG('AND MUST NOT BE USED.')
      CALL ERROR$MSG('PLEASE INFORM DEVELOPERS!')
      CALL ERROR$STOP('DFT$GRADIENTSWITCH')
      TGRA_=TGRA
      RETURN          
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE DFT$SETR8(ID,VAL)
!!$!     **************************************************************************
!!$!     **************************************************************************
!!$      USE DFT_MODULE
!!$      IMPLICIT NONE
!!$      CHARACTER(*),INTENT(IN) :: ID
!!$      LOGICAL(4)  ,INTENT(IN) :: VAL
!!$!     **************************************************************************
!!$      IF(ID.EQ.'SCALEX') THEN
!!$        SCALEX=VAL
!!$      ELSE
!!$        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
!!$        CALL ERROR$CHVAL('ID',ID)
!!$        CALL ERROR$STOP('DFT$SETR8')
!!$      END IF
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST &
     &              ,EXC,VXCT,VXCS,GVXCT2,GVXCS2,GVXCST)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATE EXCHANGE CORRELATION ENERGY AND POTENTIAL                  **
!     **                                                                      **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RHOT
      REAL(8)   ,INTENT(IN) :: RHOS
      REAL(8)   ,INTENT(IN) :: GRHOT2
      REAL(8)   ,INTENT(IN) :: GRHOS2
      REAL(8)   ,INTENT(IN) :: GRHOST
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: VXCT
      REAL(8)   ,INTENT(OUT):: VXCS
      REAL(8)   ,INTENT(OUT):: GVXCT2
      REAL(8)   ,INTENT(OUT):: GVXCS2
      REAL(8)   ,INTENT(OUT):: GVXCST
      REAL(8)               :: VAL(5)
      REAL(8)               :: DEXC(5)
      REAL(8)               :: EXC1
      REAL(8)               :: DEXC1(5)
!     **************************************************************************
      IF(.NOT.TINI) CALL DFT_INITIALIZE
      EXC=0.D0
      VXCT=0.D0
      VXCS=0.D0
      GVXCT2=0.D0
      GVXCS2=0.D0
      GVXCST=0.D0
!
!     ==================================================================
!     == INITIALIZE VALUES                                            ==
!     ==================================================================
      VAL(1)=RHOT
      VAL(3)=GRHOT2
      IF(TSPIN) THEN
        VAL(2)=RHOS
        VAL(4)=GRHOS2
        VAL(5)=GRHOST
      ELSE
        VAL(2)=0.D0
        VAL(4)=0.D0
        VAL(5)=0.D0
      END IF
      EXC    =0.D0
      DEXC(:)=0.D0
!
!     ==================================================================
!     == AVOID VERY SMALL DENSITIES                                   ==
!     ==================================================================
      IF(VAL(1).LT.RHOTMIN) THEN
!RETURN !#
        VAL(1)=RHOTMIN
      END IF
!
!     == SPIN DENSITY MUST NOT BE LARGER THAN THE TOTAL DENSITY ========
      IF(ABS(VAL(2)).GE.VAL(1)-RHOTMIN) THEN
        VAL(2)=MAX(VAL(2),-VAL(1)+RHOTMIN)
        VAL(2)=MIN(VAL(2),VAL(1)-RHOTMIN)
      END IF
!
!     ==================================================================
!     ==  EXCHANGE FUNCTIONAL                                         ==
!     ==================================================================
      IF(TX) THEN 
        CALL EXCHANGE$EVAL1(VAL,EXC1,DEXC1)
        EXC=EXC+EXC1*SCALEX
        DEXC=DEXC+DEXC1*SCALEX
      END IF
!
!     ==================================================================
!     ==  PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION           ==
!     ==================================================================
      IF(TPW91L.AND.TCORRELATION) THEN
        CALL PERDEWWANG91L$EVAL1(VAL,EXC1,DEXC1)
        EXC=EXC+EXC1
        DEXC(:)=DEXC(:)+DEXC1(:)
      END IF
!
!     ==================================================================
!     ==  PERDEW-ZUNGER PARAMETERIZATION OF LOCAL CORRELATION         ==
!     ==================================================================
      IF(TPZ.AND.TCORRELATION) THEN
        CALL PERDEWZUNGER$EVAL1(VAL,EXC1,DEXC1)
        EXC=EXC+EXC1
        DEXC(:)=DEXC(:)+DEXC1(:)
      END IF
!
!     ==================================================================
!     ==  PERDEW GRADIENT CORRECTION TO CORRELATION                   ==
!     ==================================================================
      IF(TPERDEW.AND.TCORRELATION) THEN
        CALL PERDEW$EVAL1(VAL,EXC1,DEXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
      END IF
!
!     ==================================================================
!     ==  PERDEW WANG 91 GRADIENT CORRECTION TO CORRELATION           ==
!     ==  ( INCLUDES ALSO LOCAL CORRELATION)                          ==
!     ==================================================================
      IF(TPW91G.AND.TCORRELATION) THEN
        DEXC1=0.D0
        CALL PERDEWWANG91G$EVAL(VAL(1),VAL(2),VAL(3),EXC1 &
     &                         ,DEXC1(1),DEXC1(2),DEXC1(3))
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
      END IF
!
!     ==================================================================
!     ==  PBE96 EXCHANGE CORRELATION FUNCTIONAL                       ==
!     ==================================================================
      IF(TPBE96.AND.TCORRELATION) THEN
        CALL PBE$EVAL1(VAL,EXC1,DEXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
      END IF
!
!     ==================================================================
!     ==  LYP88 CORRELATION FUNCTIONAL                                ==
!     ==================================================================
      IF(TLYP88.AND.TCORRELATION) THEN
        CALL LYP88$EVAL1(VAL,EXC1,DEXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
      END IF
!
!     ==========================================================================
!     ==  LIBXC INTERFACE 
!     ==========================================================================
!     == THE LIBXC FUNCTIONAL DOES NOT (YET) ALLOW TO SWITCH OFF CORRELATION  ==
!     == THEREFORE THIS TERM IS NOT SWITCHED OFF BY TCORRELATION=.FALSE.      ==
!     == TCORRELATION IS USED FOR HYBRID FUNCTIONALS TO SELECT ONLY EXCHANGE. ==
      IF(TLIBXC) THEN
!       == EXCHANGE (TYPE 3) IS CALLED BEFOREHAND AS REQUIRED BY PBE
!       == NO!!! EXCHANGE IS NOT SEPARATED OUT IN LIBXC
!CAUTION!!!
!       
        EXC=0.D0
        DEXC=0.D0
        CALL PAWLIBXC$GGA1(VAL,EXC1,DEXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
      END IF
!
!     ==================================================================
!     ==  WRAP-UP  (CHECKS ETC.)                                      ==
!     ==================================================================
     IF(RHOT.LT.RHOTMIN) THEN
        DEXC(1)=EXC/RHOTMIN
        DEXC(2:5)=DEXC(2:5)*RHOT/RHOTMIN
        EXC=EXC*RHOT/RHOTMIN
      END IF
!
      VXCT  =DEXC(1)
      VXCS  =DEXC(2)
      GVXCT2=DEXC(3)
      GVXCS2=DEXC(4)
      GVXCST=DEXC(5)
!
!     ==========================================================================
!     == CATCH NAN-S (NAN=NOT-A-NUMBER)                                       ==
!     ==========================================================================
      IF(EXC.NE.EXC.OR.VXCT.NE.VXCT.OR.VXCS.NE.VXCS &
     &             .OR.GVXCT2.NE.GVXCT2.OR.GVXCS2.NE.GVXCS2 &
                   .OR.GVXCST.NE.GVXCST) THEN
        CALL ERROR$MSG('ERROR AFTER DFT')
        CALL ERROR$L4VAL('TSPIN ',TSPIN)
        CALL ERROR$R8VAL('EXC ',EXC)
        CALL ERROR$R8VAL('RHOT ',RHOT)
        CALL ERROR$R8VAL('RHOS ',RHOS)
        CALL ERROR$R8VAL('GRHOT2 ',GRHOT2)
        CALL ERROR$R8VAL('GRHOS2 ',GRHOS2)
        CALL ERROR$R8VAL('GRHOST ',GRHOST)
        CALL ERROR$R8VAL('VXCT ',VXCT)
        CALL ERROR$R8VAL('VXCS ',VXCS)
        CALL ERROR$R8VAL('GVXCT2 ',GVXCT2)
        CALL ERROR$R8VAL('GVXCS2 ',GVXCS2)
        CALL ERROR$R8VAL('GVXCST ',GVXCST)
        CALL ERROR$STOP('DFT')
      END IF
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT2(VAL_,EXC,DEXC,D2EXC)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE EXCHANGE AND CORRELATION ENERGY                    **
!     **  AND ITS FIRST AND SECOND DERIVATIVES                        **
!     **                                                              **
!     ******************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL_(5)  !(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: DEXC(5)     ! DEXC/DVAL(I)
      REAL(8)   ,INTENT(OUT):: D2EXC(5,5)  ! D2EXC/(DVAL(I)DVAL(J))
      REAL(8)               :: VAL(5)
      REAL(8)               :: EXC1,DEXC1(5),D2EXC1(5,5)
!     ******************************************************************
      IF(.NOT.TINI) CALL DFT_INITIALIZE
!
!     ==================================================================
!     == INITIALIZE VALUES                                            ==
!     ==================================================================
      VAL(:)=VAL_(:)
      IF(.NOT.TSPIN) THEN
        VAL(2)=0.D0
        VAL(4)=0.D0
        VAL(5)=0.D0
      END IF
      EXC       =0.D0
      DEXC(:)   =0.D0
      D2EXC(:,:)=0.D0
!
!     ==================================================================
!     == RETURN IF TOTAL DENSITY BELOW ZERO                           ==
!     ==================================================================
      IF(VAL(1).LE.0.D0) RETURN
!
!     ==================================================================
!     == AVOID VERY SMALL DENSITIES                                   ==
!     ==================================================================
      IF(VAL(1).LT.RHOTMIN) THEN
        VAL(1)=RHOTMIN
      END IF

      IF(ABS(VAL(2)).GE.VAL(1)-RHOTMIN) THEN
        VAL(2)=MAX(VAL(2),-VAL(1)+RHOTMIN)
        VAL(2)=MIN(VAL(2),VAL(1)-RHOTMIN)
      END IF
!
!     ==================================================================
!     ==  EXCHANGE ENERGY                                             ==
!     ==================================================================
      IF(TX) THEN
        CALL EXCHANGE$EVAL2(VAL,EXC1,DEXC1,D2EXC1)
        EXC=EXC+EXC1*SCALEX
        DEXC(:)=DEXC(:)+DEXC1(:)*SCALEX
        D2EXC(:,:)=D2EXC(:,:)+D2EXC1(:,:)*SCALEX
      END IF
!
!     ==================================================================
!     ==  PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION           ==
!     ==================================================================
      IF(TPW91L.AND.TCORRELATION) THEN
        CALL PERDEWWANG91L$EVAL2(VAL,EXC1,DEXC1,D2EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
      END IF
!
!     ==================================================================
!     ==  PERDEW ZUNGER PARAMETERIZATION OF LOCAL CORRELATION         ==
!     ==================================================================
      IF(TPZ.AND.TCORRELATION) THEN
        CALL PERDEWZUNGER$EVAL2(VAL,EXC1,DEXC1,D2EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
      END IF
!
!     ==================================================================
!     ==  PERDEW GRADIENT CORRECTION TO CORRELATION                   ==
!     ==================================================================
      IF(TPERDEW.AND.TCORRELATION) THEN
        CALL PERDEW$EVAL2(VAL,EXC1,DEXC1,D2EXC1)
        EXC=EXC+EXC1
        DEXC(:)=DEXC(:)+DEXC1(:)
        D2EXC(:,:)=D2EXC(:,:)+D2EXC1(:,:)
      END IF
!
!     ==================================================================
!     ==  PERDEW WANG 91 GRADIENT CORRECTION TO CORRELATION           ==
!     ==  ( INCLUDES ALSO LOCAL CORRELATION)                          ==
!     ==================================================================
      IF(TPW91G.AND.TCORRELATION) THEN
        CALL PERDEWWANG91G$EVAL2(VAL(1),VAL(2),VAL(3),EXC1 &
     &       ,DEXC1(1),DEXC1(2),DEXC1(3) &
     &       ,D2EXC1(1,1),D2EXC1(1,2),D2EXC1(1,3) &
     &       ,D2EXC1(2,2),D2EXC1(2,3),D2EXC1(3,3))
        D2EXC1(2,1)=D2EXC1(1,2)
        D2EXC1(3,1)=D2EXC1(1,3)
        D2EXC1(3,2)=D2EXC1(2,3)
        EXC=EXC+EXC1
        DEXC(1:3)=DEXC(1:3)+DEXC1(1:3)        
        D2EXC(1:3,1:3)=D2EXC(1:3,1:3)+D2EXC(1:3,1:3)
      END IF
!
!     ==================================================================
!     ==  PBE96 EXCHANGE CORRELATION FUNCTIONAL                       ==
!     ==================================================================
      IF(TPBE96.AND.TCORRELATION) THEN
        CALL PBE$EVAL2(VAL,EXC1,DEXC1,D2EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
      END IF
!
!     ==================================================================
!     ==  LEE YANG-PARR 88 CORRELATION                                ==
!     ==================================================================
      IF(TLYP88.AND.TCORRELATION) THEN
        CALL LYP88$EVAL2(VAL,EXC1,DEXC1,D2EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
      END IF
!
!     ==========================================================================
!     ==  LIBXC INTERFACE
!     ==========================================================================
!     == THE LIBXC INTERFACE DOES NOT (YET) ALLOW TO SWITCH OFF CORRELATION   ==
!     == THEREFORE THIS TERM IS NOT SWITCHED OFF BY TCORRELATION=.FALZSE.     ==
!     == TCORRELATION IS USED FOR HYBRID FUNCTIONALS TO SELECT ONLY EXCHANGE. ==
      IF(TLIBXC) THEN
        CALL PAWLIBXC$GGA2(VAL,EXC1,DEXC1,D2EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
      END IF
!
!     ==================================================================
!     ==  WRAPUP AND CHECKS                                           ==
!     ==================================================================
      IF(VAL_(1).LT.RHOTMIN) THEN
        D2EXC(1,1)=0.D0
        D2EXC(1,2:5)=DEXC(2:5)/RHOTMIN
        D2EXC(2:5,1)=DEXC(2:5)/RHOTMIN
        D2EXC(2:5,2:5)=D2EXC(2:5,2:5)*VAL_(1)/RHOTMIN
        DEXC(1)=EXC/RHOTMIN
        DEXC(2:5)=DEXC(2:5)*VAL_(1)/RHOTMIN
        EXC=EXC*VAL_(1)/RHOTMIN
      END IF
      IF(EXC.NE.EXC.OR.DEXC(1).NE.DEXC(1).OR.DEXC(2).NE.DEXC(2) &
            .OR.DEXC(3).NE.DEXC(3).OR.D2EXC(1,1).NE.D2EXC(1,1) &
            .OR.D2EXC(1,2).NE.D2EXC(1,2).OR.D2EXC(1,3).NE.D2EXC(1,3) &
            .OR.D2EXC(2,2).NE.D2EXC(2,2).OR.D2EXC(2,3).NE.D2EXC(2,3) &
            .OR.D2EXC(3,3).NE.D2EXC(3,3)) THEN
        CALL ERROR$R8VAL('VAL(1)',VAL(1))
        CALL ERROR$R8VAL('VAL(2)',VAL(2))
        CALL ERROR$R8VAL('VAL(3)',VAL(3))
        CALL ERROR$R8VAL('EXC ',EXC)
        CALL ERROR$R8VAL('DEXC(1)',DEXC(1))
        CALL ERROR$R8VAL('DEXC(2)',DEXC(2))
        CALL ERROR$R8VAL('DEXC(3)',DEXC(3))
        CALL ERROR$R8VAL('D2EXC(1,1)',D2EXC(1,1))
        CALL ERROR$R8VAL('D2EXC(1,2)',D2EXC(1,2))
        CALL ERROR$R8VAL('D2EXC(1,3)',D2EXC(1,3))
        CALL ERROR$R8VAL('D2EXC(2,2)',D2EXC(2,2))
        CALL ERROR$R8VAL('D2EXC(2,3)',D2EXC(2,3))
        CALL ERROR$R8VAL('D2EXC(3,3)',D2EXC(3,3))
        CALL ERROR$STOP('DFT2')
      END IF
      RETURN
    END SUBROUTINE DFT2
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT3(VAL_,EXC,DEXC,D2EXC,D3EXC)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATE EXCHANGE AND CORRELATION ENERGY                            **
!     **  AND ITS FIRST, SECOND AND THIRD DERIVATIVES                         **
!     **                                                                      **
!     **************************************************************************
      USE DFT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: VAL_(5)  !(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST)
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)   ,INTENT(OUT):: DEXC(5)     ! DEXC/DVAL(I)
      REAL(8)   ,INTENT(OUT):: D2EXC(5,5)  ! D2EXC/(DVAL(I)DVAL(J))
      REAL(8)   ,INTENT(OUT):: D3EXC(5,5,5)! D3EXC/(DVAL(I)DVAL(J),DVAL(K))
      REAL(8)               :: VAL(5)
      REAL(8)               :: EXC1,DEXC1(5),D2EXC1(5,5),D3EXC1(5,5,5)
!     **************************************************************************
      IF(.NOT.TINI) CALL DFT_INITIALIZE
!
!     ==========================================================================
!     == INITIALIZE VALUES                                                    ==
!     ==========================================================================
      VAL(:)=VAL_(:)
      IF(.NOT.TSPIN) THEN
       VAL(2)=0.D0
       VAL(4)=0.D0
       VAL(5)=0.D0
      END IF
      EXC       =0.D0
      DEXC(:)   =0.D0
      D2EXC(:,:)=0.D0
      D3EXC(:,:,:)=0.D0
!
!     ==========================================================================
!     == RETURN IF TOTAL DENSITY BELOW ZERO                                   ==
!     ==========================================================================
      IF(VAL(1).LE.0.D0) RETURN
!
!     ==========================================================================
!     == AVOID VERY SMALL DENSITIES                                           ==
!     ==========================================================================
      IF(VAL(1).LT.RHOTMIN) THEN
        VAL(1)=RHOTMIN
      END IF

      IF(ABS(VAL(2)).GE.VAL(1)-RHOTMIN) THEN
        VAL(2)=MAX(VAL(2),-VAL(1)+RHOTMIN)
        VAL(2)=MIN(VAL(2),VAL(1)-RHOTMIN)
      END IF
!
!     ==========================================================================
!     ==  EXCHANGE ENERGY                                                     ==
!     ==========================================================================
      IF(TX) THEN
        CALL EXCHANGE$EVAL3(VAL,EXC1,DEXC1,D2EXC1,D3EXC1)
        EXC=EXC+EXC1*SCALEX
        DEXC(:)=DEXC(:)+DEXC1(:)*SCALEX
        D2EXC(:,:)=D2EXC(:,:)+D2EXC1(:,:)*SCALEX
        D3EXC(:,:,:)=D3EXC(:,:,:)+D3EXC1(:,:,:)*SCALEX
      END IF
!
!     ==========================================================================
!     ==  PERDEW WANG PARAMETERIZATION OF LOCAL CORRELATION                   ==
!     ==========================================================================
      IF(TPW91L.AND.TCORRELATION) THEN
        CALL PERDEWWANG91L$EVAL3(VAL,EXC1,DEXC1,D2EXC1,D3EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
        D3EXC=D3EXC+D3EXC1
      END IF
!
!     ==========================================================================
!     ==  PERDEW ZUNGER PARAMETERIZATION OF LOCAL CORRELATION                 ==
!     ==========================================================================
      IF(TPZ.AND.TCORRELATION) THEN
        CALL PERDEWZUNGER$EVAL3(VAL,EXC1,DEXC1,D2EXC1,D3EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
        D3EXC=D3EXC+D3EXC1
      END IF
!
!     ==========================================================================
!     ==  PERDEW GRADIENT CORRECTION TO CORRELATION                           ==
!     ==========================================================================
      IF(TPERDEW.AND.TCORRELATION) THEN
        CALL PERDEW$EVAL3(VAL,EXC1,DEXC1,D2EXC1,D3EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
        D3EXC=D3EXC+D3EXC1
      END IF
!
!     ==================================================================
!     ==  PERDEW WANG 91 GRADIENT CORRECTION TO CORRELATION           ==
!     ==  ( INCLUDES ALSO LOCAL CORRELATION)                          ==
!     ==================================================================
      IF(TPW91G.AND.TCORRELATION) THEN
        CALL PERDEWWANG91G$EVAL2(VAL(1),VAL(2),VAL(3),EXC1 &
     &       ,DEXC1(1),DEXC1(2),DEXC1(3) &
     &       ,D2EXC1(1,1),D2EXC1(1,2),D2EXC1(1,3) &
     &       ,D2EXC1(2,2),D2EXC1(2,3),D2EXC1(3,3))
        D2EXC1(2,1)=D2EXC1(1,2)
        D2EXC1(3,1)=D2EXC1(1,3)
        D2EXC1(3,2)=D2EXC1(2,3)
        EXC=EXC+EXC1
        DEXC(1:3)=DEXC(1:3)+DEXC1(1:3)        
        D2EXC(1:3,1:3)=D2EXC(1:3,1:3)+D2EXC(1:3,1:3)
      END IF
!
!     ==========================================================================
!     ==  PBE96 EXCHANGE CORRELATION FUNCTIONAL                               ==
!     ==========================================================================
      IF(TPBE96.AND.TCORRELATION) THEN
        CALL PBE$EVAL3(VAL,EXC1,DEXC1,D2EXC1,D3EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
        D3EXC=D3EXC+D3EXC1
      END IF
!
!     ==========================================================================
!     ==  LEE YANG-PARR 88 CORRELATION FUNCTIONAL                             ==
!     ==========================================================================
      IF(TLYP88.AND.TCORRELATION) THEN
        CALL LYP88$EVAL3(VAL,EXC1,DEXC1,D2EXC1,D3EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
        D3EXC=D3EXC+D3EXC1
      END IF
!
!     ==========================================================================
!     ==  LIBXC INTERFACE
!     ==========================================================================
!     == THE LIBXC INTERFACE DOES NOT (YET) ALLOW TO SWITCH OFF CORRELATION   ==
!     == THEREFORE THIS TERM IS NOT SWITCHED OFF BY TCORRELATION=.FALZSE.     ==
!     == TCORRELATION IS USED FOR HYBRID FUNCTIONALS TO SELECT ONLY EXCHANGE. ==
      IF(TLIBXC) THEN
        CALL PAWLIBXC$GGA3(VAL,EXC1,DEXC1,D2EXC1,D3EXC1)
        EXC=EXC+EXC1
        DEXC=DEXC+DEXC1
        D2EXC=D2EXC+D2EXC1
        D3EXC=D3EXC+D3EXC1
      END IF
!
!     ==========================================================================
!     ==  WRAPUP AND CHECKS                                                   ==
!     ==========================================================================
      IF(VAL_(1).LT.RHOTMIN) THEN
        D3EXC(1,1,1)=0.D0
        D3EXC(1,1,2:5)=0.D0
        D3EXC(1,2:5,1)=0.D0
        D3EXC(2:5,1,1)=0.D0
        D3EXC(1,2:5,2:5)=D2EXC(2:5,2:5)/RHOTMIN
        D3EXC(2:5,1,2:5)=D2EXC(2:5,2:5)/RHOTMIN
        D3EXC(2:5,2:5,1)=D2EXC(2:5,2:5)/RHOTMIN
        D3EXC(2:5,2:5,2:5)=D3EXC(2:5,2:5,2:5)*VAL_(1)/RHOTMIN
        D2EXC(1,1)=0.D0
        D2EXC(1,2:5)=DEXC(2:5)/RHOTMIN
        D2EXC(2:5,1)=DEXC(2:5)/RHOTMIN
        D2EXC(2:5,2:5)=D2EXC(2:5,2:5)*VAL_(1)/RHOTMIN
        DEXC(1)=EXC/RHOTMIN
        DEXC(2:5)=DEXC(2:5)*VAL_(1)/RHOTMIN
        EXC=EXC*VAL_(1)/RHOTMIN
      END IF
      IF(EXC.NE.EXC.OR.DEXC(1).NE.DEXC(1).OR.DEXC(2).NE.DEXC(2) &
            .OR.DEXC(3).NE.DEXC(3).OR.D2EXC(1,1).NE.D2EXC(1,1) &
            .OR.D2EXC(1,2).NE.D2EXC(1,2).OR.D2EXC(1,3).NE.D2EXC(1,3) &
            .OR.D2EXC(2,2).NE.D2EXC(2,2).OR.D2EXC(2,3).NE.D2EXC(2,3) &
            .OR.D2EXC(3,3).NE.D2EXC(3,3)) THEN
        CALL ERROR$R8VAL('VAL(1)',VAL(1))
        CALL ERROR$R8VAL('VAL(2)',VAL(2))
        CALL ERROR$R8VAL('VAL(3)',VAL(3))
        CALL ERROR$R8VAL('EXC ',EXC)
        CALL ERROR$R8VAL('DEXC(1)',DEXC(1))
        CALL ERROR$R8VAL('DEXC(2)',DEXC(2))
        CALL ERROR$R8VAL('DEXC(3)',DEXC(3))
        CALL ERROR$R8VAL('D2EXC(1,1)',D2EXC(1,1))
        CALL ERROR$R8VAL('D2EXC(1,2)',D2EXC(1,2))
        CALL ERROR$R8VAL('D2EXC(1,3)',D2EXC(1,3))
        CALL ERROR$R8VAL('D2EXC(2,2)',D2EXC(2,2))
        CALL ERROR$R8VAL('D2EXC(2,3)',D2EXC(2,3))
        CALL ERROR$R8VAL('D2EXC(3,3)',D2EXC(3,3))
        CALL ERROR$STOP('DFT3')
      END IF
      RETURN
    END SUBROUTINE DFT3
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DFT$MACDONALD(RHO,E,POT)
!     **************************************************************************
!     ** RELATIVISTIC CORRECTION TO ENERGY DENSITY                            **
!     ** OF MACDONALD AND VOSKO, MACDONALD79_JPCSS12_2977                     **
!     ** MACDONALD AND VOSKO, J.PHYS. SOL.ST. PHYS 12, 2977 (1979)            **
!     **                                                                      **
!     ** PHIREL= MACDONALD EQ 3.4 (CAUTION TYPO: N <-> ETA)                   **
!     **                                                                      **
!     ** ATTENTION: I USE THE LOCAL DENSITY EXCHANGE ENERGY OF THE            **
!     ** NON-MAGNETIC ELECTRON GAS.                                           **
!     **                                                                      **
!     ** THE CORRECTION IS FORMULATED AS ADDITIVE CORRECTION RATHER THAN AS   **
!     ** MULTIPLICATIVE CORRECTION AS IN MACDONALD AND VOSKO.                 **
!     *************************PETER BLOECHL, GOSLAR 24.12.2023*****************
      IMPLICIT NONE 
      REAL(8)    ,INTENT(IN) :: RHO
      REAL(8)    ,INTENT(OUT):: E
      REAL(8)    ,INTENT(OUT):: POT
      REAL(8)    ,PARAMETER  :: ONETHIRD=1.D0/3.D0
      REAL(8)    ,PARAMETER  :: FOURTHIRD=4.D0/3.D0
      REAL(8)    ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)    ,PARAMETER  :: SPEEDOFLIGHT=137.035999084D0
      REAL(8)    ,PARAMETER  :: EXFAC=-(3.d0/4.d0)*( 3.D0/PI )**(1.D0/3.D0)
      REAL(8)                :: PHIREL,DPHIREL
      REAL(8)                :: BETA,DBETA
      REAL(8)                :: ETA,DETA
      REAL(8)                :: SVAR,DSVAR
      REAL(8)                :: EX,DEX
!     **************************************************************************
!
!     ==========================================================================
!     == avoid divide-by-zero for small densities/ non-relativistic limit ======
!     ==========================================================================
      if(rho.lt.1.d-12) then
        ex=0.d0
        dex=0.d0
        return
      end if
!
!     ==========================================================================
!     == EXCHANGE ENERGY DENSITY OF HOMOGENEOUS ELECTRON GAS ===================
!     ==========================================================================
!     == RS=(3.D0/(4.D0*PI*RHO))**(1/3)
!     == EX=-3.D0*( 9.D0/(32.D0*PI**2)**ONETHIRD / RS * RHO
!     ==   =-1.5D0*( 3.D0/PI )**(1/3) * RHO**(4/3)
!     ==========================================================================
      EX=EXFAC*RHO**FOURTHIRD
      DEX=FOURTHIRD*EXFAC*RHO**ONETHIRD
!
!     ==========================================================================
!     == RELATIVISTIC CORRECTION FACTOR FROM MACDONALD AND VOSKO              ==
!     == PHIREL=PHIC+PHIT                                                     ==
!     ==========================================================================
      BETA=(3.D0*PI**2*RHO)**ONETHIRD / SPEEDOFLIGHT
      DBETA=ONETHIRD*BETA/RHO
!
      ETA=SQRT(1.D0+BETA**2)
      DETA=BETA/ETA*DBETA
!
      SVAR=(BETA*ETA-LOG(BETA+ETA)) / BETA**2 
      DSVAR=(DBETA*ETA+BETA*DETA-(DBETA+DETA)/(BETA+ETA) ) / BETA**2 &
     &     -2.D0*SVAR/BETA*DBETA
!
      PHIREL=1.D0-1.5D0*SVAR**2
      DPHIREL=-3.D0*SVAR*DSVAR
!
!     ==========================================================================
!     == CONVERT INTO A CORRECTION TO THE ENERGY DENSITY                      ==
!     ==========================================================================
      E=EX*(PHIREL-1.D0)
      POT=EX*DPHIREL+DEX*(PHIREL-1.D0)
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE TABLE1D_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: TABLE1D                                                            **
!**                                                                           **
!**  PURPOSE: INTERPOLATES ON A EQUISPACED GRID                               **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    TABLE1D$MAKE(XTABLE,XMIN,XMAX,NX                                       **
!**    TABLE1D$XOFI(XTABLE,I,X)                                               **
!**    TABLE1D$SETX(XTABLE,X,TCHK)                                            **
!**    TABLE1D$GETF(XTABLE,F,FOFX)                                            **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    ALL ROUTINES ARE INTERNAL TO THE MODULE, ALWAYS USE THE MODULE         **
!**                                                                           **
!**  DEPENDENCIES:                                                            **
!**    ERROR_MODULE                                                           **
!**                                                                           **
!*******************************************************************************
TYPE XTABLE_TYPE
  INTEGER(4) :: IX
  REAL(8)    :: W(4)
  INTEGER(4) :: NX
  REAL(8)    :: XMIN
  REAL(8)    :: XMAX
  REAL(8)    :: DX
  INTEGER(4) :: NHITS
  INTEGER(4) :: NNONHITS
END TYPE XTABLE_TYPE
!***********************************************************************
CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TABLE1D$MAKE(XTABLE,XMIN,XMAX,NX)
!     ******************************************************************
!     **  F(X) = W1*F(I1)+W2*F(I1+1)+W3*F(I1+2)+W4*F(I1+3)            **
!     ******************************************************************
      IMPLICIT NONE
      TYPE (XTABLE_TYPE),INTENT(OUT):: XTABLE
      REAL(8)           ,INTENT(IN) :: XMIN
      REAL(8)           ,INTENT(IN) :: XMAX
      INTEGER(4)        ,INTENT(IN) :: NX
!     ******************************************************************
      XTABLE%IX=0
      XTABLE%W(:)=0.D0
      XTABLE%NX=NX
      XTABLE%DX=(XMAX-XMIN)/DBLE(NX-1)
      XTABLE%XMIN=XMIN
      XTABLE%XMAX=XMAX
!
      XTABLE%NHITS=0
      XTABLE%NNONHITS=0
      RETURN
      END SUBROUTINE TABLE1D$MAKE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TABLE1D$XOFI(XTABLE,I,X)
!     ******************************************************************
!     **  F(X) = W1*F(I1)+W2*F(I1+1)+W3*F(I1+2)+W4*F(I1+3)            **
!     ******************************************************************
      IMPLICIT NONE
      TYPE (XTABLE_TYPE),INTENT(IN)  :: XTABLE
      INTEGER(4)        ,INTENT(IN)  :: I
      REAL(8)           ,INTENT(OUT) :: X
!     ******************************************************************
      X=XTABLE%XMIN+XTABLE%DX*REAL(I-1,KIND=8)
      RETURN
      END SUBROUTINE TABLE1D$XOFI
!     ..................................................................
      SUBROUTINE TABLE1D$SETX(XTABLE,X,TCHK)
!     ******************************************************************
!     **  F(X) = W1*F(I1)+W2*F(I1+1)+W3*F(I1+2)+W4*F(I1+3)            **
!     ******************************************************************
      IMPLICIT NONE
      TYPE (XTABLE_TYPE),INTENT(INOUT) :: XTABLE
      REAL(8)           ,INTENT(IN)    :: X
      LOGICAL(4)        ,INTENT(OUT)   :: TCHK  ! X LIES INSIDE TABLE
      REAL(8)           ,PARAMETER     :: ONEBY6=1.D0/6.D0
      REAL(8)                          :: XREL,XREL2,XREL3
      INTEGER(4)                       :: IREL
!     ******************************************************************
!
!     ==========================================================================
!     ==  RETURN WITH I1=0 IF X OUTSIDE OF GRID                       ==
!     ==========================================================================
      TCHK=(X.GE.XTABLE%XMIN.AND.X.LE.XTABLE%XMAX) 
      IF(.NOT.TCHK) THEN
        XTABLE%IX=0
        XTABLE%NNONHITS=XTABLE%NNONHITS+1
        RETURN
      ELSE
        XTABLE%NHITS=XTABLE%NHITS+1
      END IF
!
!     ==========================================================================
!     ==  OBTAINE POSITION ON THE GRID                                ==
!     ==========================================================================
      XREL=(X-XTABLE%XMIN)/XTABLE%DX+1.D0
      IREL=INT(XREL)
      IREL=MIN(XTABLE%NX-2,IREL)
      IREL=MAX(2,IREL)
      XTABLE%IX=IREL-1
!
!     ==========================================================================
!     ==  OBTAIN INTERPOLATION WEIGHTS                                ==
!     ==========================================================================
      XREL=XREL-DBLE(IREL)
      XREL2=XREL*XREL
      XREL3=XREL2*XREL
      XTABLE%W(1)=ONEBY6*     (-2.D0*XREL+3.D0*XREL2-XREL3)
      XTABLE%W(3)=ONEBY6*3.D0*(+2.D0*XREL     +XREL2-XREL3)
      XTABLE%W(4)=ONEBY6*          (-XREL           +XREL3)
      XTABLE%W(2)=1.D0-XTABLE%W(1)-XTABLE%W(3)-XTABLE%W(4)
      RETURN
      END SUBROUTINE TABLE1D$SETX
!     ..................................................................
      SUBROUTINE TABLE1D$GETF(XTABLE,NX,F,FOFX)
!     ******************************************************************
!     **  F(X) = W1*F(I1)+W2*F(I1+1)+W3*F(I1+2)+W4*F(I1+3)            **
!     ******************************************************************
      IMPLICIT NONE
      TYPE (XTABLE_TYPE),INTENT(IN) :: XTABLE
      INTEGER(4)        ,INTENT(IN) :: NX
      REAL(8)           ,INTENT(IN) :: F(NX)
      REAL(8)           ,INTENT(OUT):: FOFX
      INTEGER(4)                    :: IX
!     ******************************************************************
      IX=XTABLE%IX
      IF(IX.EQ.0) THEN
        CALL ERROR$MSG('XTABLE IS NOT POSITIONED')
        CALL ERROR$STOP('TABLE1D$GETF')
      END IF
      IF(NX.NE.XTABLE%NX) THEN
        CALL ERROR$MSG('LENGTH OF ARRAY NOT CONSISTENT WITH XTABLE')
        CALL ERROR$STOP('TABLE1D$GETF')
      END IF
!
      FOFX=DOT_PRODUCT(XTABLE%W(:),F(IX:IX+3))
      RETURN
      END SUBROUTINE TABLE1D$GETF
END MODULE TABLE1D_MODULE
!
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE EXCHANGE_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: EXCHANGE                                                   **
!**                                                                   **
!**  PURPOSE: CALCULATE THE EXCHANGE ENERGY                           **
!**    (WITHOUT OR WITH GRADIENT CORRECTION OF BECKE OR               **
!**     PERDEW-BURKE-ERNZERHOF)                                       **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    EXCHANGE$SETI4(ID,VAL)                                         **
!**    EXCHANGE$EVAL1(VAL,EXC,DEXC)                                   **
!**    EXCHANGE$EVAL2(VAL,EXC,DEXC,D2EXC)                             **
!**    EXCHANGE$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)                       **
!**                                                                   **
!**  METHODS:                                                         **
!**                                                                   **
!***********************************************************************
USE DFT_MODULE, ONLY: TSAFE,TSPIN
USE TABLE1D_MODULE
LOGICAL(4)         :: TINI=.FALSE.
!=======================================================================
!== THIS PARAMETER MUST BE DEFINED BEFORE THE FIRST CALL              ==
!=======================================================================
INTEGER(4)         :: ITYPE=-1
!=======================================================================
!== PARAMETERS AND SWITCHES CALCULATED DURING INITIALIZATION============
!=======================================================================
LOGICAL(4)         :: TGRA=.TRUE.
CHARACTER(8)       :: FXTYPE='NONE' ! ID FOR FX: EXC_GGA=EXC_LDA*FX
LOGICAL(4)         :: TXALPHA=.FALSE.
REAL(8)            :: S2FAC   ! S**2= S2FAC*(GRHO/RHO**(4/3))**2
REAL(8)            :: EXFAC   ! LOCAL EX=EXFAC*RHO**(4/3)
REAL(8)            :: MU
REAL(8)            :: KAPPA     
REAL(8)            :: MUBYKAPPA
REAL(8)            :: ALPHAB2  ! USED TO OBTAIN E_X=RHO/2R SCALING
REAL(8)            :: XALPHA=2.D0/3.D0
!===============================================================================
!== ARRAYS FOR TABLE LOOKUP OF THE GRADIENT ENHANCEMENT FACTOR        ==
!===============================================================================
INTEGER(4),PARAMETER   :: NS2=10000
REAL(8)   ,ALLOCATABLE :: FXARRAY(:)   !(NS2)
REAL(8)   ,ALLOCATABLE :: DFXARRAY(:)  !(NS2)
REAL(8)   ,ALLOCATABLE :: D2FXARRAY(:) !(NS2)
REAL(8)   ,ALLOCATABLE :: D3FXARRAY(:) !(NS2)
TYPE(XTABLE_TYPE)      :: S2TABLE
CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE_INITIALIZE
      IMPLICIT NONE
      REAL(8)           :: PI
      REAL(8)           :: KFFAC
      REAL(8)           :: NU
      REAL(8),PARAMETER :: CC0=4.235D-3
      REAL(8)           :: BETA    ! BETA=NU*CC0 (PWG:EQ.13)
      INTEGER(4)        :: IPAR
      REAL(8)           :: S2
      INTEGER(4)        :: I
!     *****************************************************************
!
!     ==========================================================================
!     ==  ALLOCATE ARRAYS FOR EXCHANGE                                        ==
!     ==  ALLOCATED HERE RATHER THAN IN THE MODULE TO NOT OVERLOAD THE STACK  ==
!     ==========================================================================
      ALLOCATE(FXARRAY(NS2))
      ALLOCATE(DFXARRAY(NS2))
      ALLOCATE(D2FXARRAY(NS2))
      ALLOCATE(D3FXARRAY(NS2))
!
!     =========================================================================
!     ==  INITIALIZE PARAMETERS TO LOCAL EXCHANGE ETC.               ==
!     =========================================================================
      IPAR=2
      TXALPHA=.FALSE.
      TGRA=.FALSE.
!
!     =========================================================================
!     ==  RESOLVE SELECTION FOR EXCHANGE FUNCTIONAL                  ==
!     =========================================================================
      IF(ITYPE.EQ.1) THEN
!       == LOCAL EXCHANGE =============================================        
      ELSE IF(ITYPE.EQ.12) THEN
!       == X_ALPHA WITH THE COMMONLY USED VALUE X_ALPHA=0.7 ===========
        TXALPHA=.TRUE.
        XALPHA=0.7D0
      ELSE IF(ITYPE.EQ.2) THEN
!       == LOCAL EXCHANGE +BECKE 88 GRADIENT CORRECTION ================        
        TGRA=.TRUE.
        FXTYPE='BECKE'
        IPAR=2       ! SELECT PARAMETERS CONSISTENT WITH BECKE88
      ELSE IF(ITYPE.EQ.23) THEN
!       == LOCAL EXCHANGE =============================================
!       == + BECKE 88 GRADIENT CORRECTION WITHOUT PROPER LONG-RANGE SCALING
!       == ACCORDING TO PBE
        TGRA=.TRUE.
        FXTYPE='PBE'
        IPAR=2       ! SELECT PARAMETERS CONSISTENT WITH BECKE88
      ELSE IF(ITYPE.EQ.3) THEN
!       == LOCAL EXCHANGE + PBE96 GRADIENT CORRECTION ================        
        TGRA=.TRUE.
        FXTYPE='PBE'
        IPAR=3       ! SELECT PARAMETERS CONSISTENT WITH PBE EXCHANGE
      ELSE IF(ITYPE.EQ.4) THEN
!       == LOCAL EXCHANGE + PBE96 GRADIENT CORRECTION ================        
        TGRA=.TRUE.
        FXTYPE='RPBE'
        IPAR=3       ! SELECT PARAMETERS CONSISTENT WITH PBE EXCHANGE
      ELSE IF(ITYPE.EQ.5) THEN
!       == LOCAL EXCHANGE + REVPBE GRADIENT CORRECTION ================        
        TGRA=.TRUE.
        FXTYPE='PBE'
        IPAR=4       ! SELECT PARAMETERS CONSISTENT WITH REVPBE EXCHANGE
      ELSE IF(ITYPE.EQ.32) THEN
!       == LOCAL EXCHANGE ==============================================
!       == + PBE96 GRADIENT CORRECTION CORRECT LONG-RANGE SCALING ======
        TGRA=.TRUE.
        FXTYPE='BECKE'
        IPAR=3       ! SELECT PARAMETERS CONSISTENT WITH PBE EXCHANGE
      ELSE IF(ITYPE.EQ.111111) THEN
!       == LIBXC INTERFACE =====================================================
        TGRA=.TRUE.
! ATOMIC (SPHERICAL,NON-SPIN-POLARIZED) CALCULATIONS GIVE CORRECT NUMBERS,
! BUIT FAIL TO CONVERGE, FOR EXAMPLE, FOR MG
        CALL ERROR$MSG('THIS CODE (LIBXC) NEEDS CHECKING')
        CALL ERROR$STOP('EXCHANGE_INITIALIZE')
      ELSE
        CALL ERROR$MSG('FUNCTIONAL SELECTION FOR EXCHANGE NOT RECOGNIZED')
        CALL ERROR$I4VAL('ITYPE',ITYPE)
        CALL ERROR$STOP('EXCHANGE_INITIALIZE (SEE DFT OBJECT)')
      END IF
!
!     =========================================================================
!     ==  CALCULATE PARAMETERS                                       ==
!     =================================================================
      PI=4.D0*ATAN(1.D0)
      KFFAC=(3.D0*PI**2)**(1.D0/3.D0)     ! KF=KFFAC/RS
      EXFAC=-3.D0/(4.D0*PI)*KFFAC         ! EX_LOC(NT)=EXFAC*RHO**(4/3)
      S2FAC=0.25D0/KFFAC**2               ! PREFACTOR FOR THE REDUCED GRADIENT
      IF(IPAR.EQ.2) THEN
!       ================================================================
!       ==  BECKE88 EXCHANGE (PRA38,3098(1988) IS OBTAINED            ==
!       ==  WITH THE FOLLOWING PARAMETERS                             ==
!       ==  MU   =16.D0*PI/3.D0*(6.D0*PI**2)**(1.D0/3.D0)*0.0042D0    ==
!       ==       =0.2742931D0                                         ==
!       ==  KAPPA=2.D0*PI/9.D0/(6.D0*PI**2)**(1.D0/3.D0)              ==
!       ==       =0.1791102D0                                         ==
!       ================================================================
        MU=16.D0*PI/3.D0*(6.D0*PI**2)**(1.D0/3.D0)*0.0042D0
        KAPPA=2.D0*PI/9.D0/(6.D0*PI**2)**(1.D0/3.D0)
      ELSE IF(IPAR.EQ.3) THEN
!       ================================================================
!       == PARAMETERS FROM PBE PAPER                                  ==
!       == THE PARAMETERS MAKE SENS ONLY IN CONNECTION WITH THE       ==
!       == PBE CORRELATION FUNCTIONAL                                 ==
!       == KAPPA MAKES SENSE ONLY WITH THE SIMPLE FORM                ==
!       ================================================================
!       == MU IS CHOSEN TO CANCEL THE GRADIENT DEPENDENCE OF ===========
!       == THE PBE CORRELATION IN THE SMALL GRADIENT LIMIT AND RHOS=0 ==
        NU=16.D0/PI*(3*PI**2)**(1.D0/3.D0)
        BETA=NU*CC0
        MU=BETA*PI**2/3.D0
!       == KAPPA IS CHOSEN TO SATISFY THE LIEB OXFORD BOUND  ============
!       == EX(GRHO->INFTY)=EX(GRHO=0)*KAPPA =============================
        KAPPA=0.804D0
      ELSE IF(IPAR.EQ.4) THEN
!       ================================================================
!       == PARAMETERS FROM REVPBE PAPER                               ==
!       == THE PARAMETERS MAKE SENS ONLY IN CONNECTION WITH THE       ==
!       == PBE CORRELATION FUNCTIONAL                                 ==
!       == KAPPA MAKES SENSE ONLY WITH THE SIMPLE FORM                ==
!       ================================================================
!       == MU IS CHOSEN TO CANCEL THE GRADIENT DEPENDENCE OF ===========
!       == THE PBE CORRELATION IN THE SMALL GRADIENT LIMIT AND RHOS=0 ==
        NU=16.D0/PI*(3*PI**2)**(1.D0/3.D0)
        BETA=NU*CC0
        MU=BETA*PI**2/3.D0
!       == KAPPA IS CHOSEN TO SATISFY THE LIEB OXFORD BOUND  ============
!       == EX(GRHO->INFTY)=EX(GRHO=0)*KAPPA =============================
        KAPPA=1.245D0
      ELSE
        CALL ERROR$MSG('AT LEAST ONE PARAMETER SET MUST BE SELECTED')
        CALL ERROR$STOP('EXCHANGE_INITIALIZE')
      END IF
      MUBYKAPPA=MU/KAPPA
!     == ALPHAB2 IS NEEDED FOR THE FORM WITH PROPER LONG-RANGE BEHAVIOR ========
      ALPHAB2=(4.D0*PI/(9.D0*KAPPA))**2
!
!     ==========================================================================
!     ==  DEFINE ARRAYS FOR FAST TABLE LOOKUP                        ==
!     ==========================================================================
      IF(TGRA) THEN
        CALL TABLE1D$MAKE(S2TABLE,1.D-3,10.D0,NS2)
        DO I=1,NS2
          CALL TABLE1D$XOFI(S2TABLE,I,S2)
          CALL EXCHANGE_FX3(S2,FXARRAY(I),DFXARRAY(I),D2FXARRAY(I),D3FXARRAY(I))
        ENDDO
      END IF
!
      TINI=.TRUE.
      END SUBROUTINE EXCHANGE_INITIALIZE
END MODULE EXCHANGE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE$SETI4(ID,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE EXCHANGE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID.EQ.'TYPE') THEN
        ITYPE=VAL_
        TINI=.FALSE.
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE$EVAL1(VAL,EXC,DEXC)
      USE EXCHANGE_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: VAL(5)     !RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS
      REAL(8),INTENT(OUT):: EXC        ! EXCHANGE ENERGY (RHOT*EPS_X)
      REAL(8),INTENT(OUT):: DEXC(5)    ! DEXC(I)=D[EXC]/D[VAL(I)]
      REAL(8)            :: RHO(2),GRHO2(2)
      REAL(8)            :: EX(2),EX_R(2),EX_G(2)
      INTEGER(4)         :: I
      INTEGER(4)         :: NSPIN
!     ******************************************************************
      IF(.NOT.TINI) CALL EXCHANGE_INITIALIZE
      EXC=0.D0
      DEXC(:)=0.D0
!
!     ==========================================================================
!     == TRANSFORM INTO SPIN UP AND SPIN DOWN COMPONENTS              ==
!     ==========================================================================
      IF(TSPIN) THEN
        RHO(1)=VAL(1)+VAL(2)               !=2*RHO(UP)
        RHO(2)=VAL(1)-VAL(2)               !=2*RHO(DOWN)
        GRHO2(1)=VAL(3)+2.D0*VAL(5)+VAL(4) !=2*GRAD[RHO(UP)]
        GRHO2(2)=VAL(3)-2.D0*VAL(5)+VAL(4) !=2*GRAD[RHO[DOWN]
        NSPIN=2
      ELSE
        RHO(1)=VAL(1)
        GRHO2(1)=VAL(3)
        NSPIN=1
      END IF
!
!     ==================================================================
!     == EXCHANGE_X2 ASSUMES A NON-SPIN POLARIZED TOTAL DENSITY       ==
!     == E_X(RHO_UP,RHO_DOWN)=0.5*(E_X(2*RHO_UP,0)+E_X(2*RHO_DOWN,0))/2=
!     ==================================================================
      DO I=1,NSPIN
        CALL EXCHANGE_X1(RHO(I),GRHO2(I),EX(I),EX_R(I),EX_G(I))
      ENDDO
!
!     ==========================================================================
!     == MODIFY FOR XALPHA                                            ==
!     ==========================================================================
      IF(TXALPHA) THEN
        EX(:)=1.5D0*XALPHA*EX(:)
        EX_R(:)=1.5D0*XALPHA*EX_R(:)
        EX_G(:)=1.5D0*XALPHA*EX_G(:)
      END IF
!
!     ==========================================================================
!     == TRANSFORM BACK TO TOTAL AND SPIN DENSITIES                   ==
!     ==========================================================================
      IF(TSPIN) THEN
        EX(:)   =0.5D0*EX(:)
        EX_R(:) =0.5D0*EX_R(:)
        EXC       =      EX(1)+EX(2)
        DEXC(1)   =      EX_R(1)+EX_R(2)
        DEXC(2)   =      EX_R(1)-EX_R(2)
        IF(.NOT.TGRA) THEN
          RETURN
        END IF
        EX_G(:) =0.5D0*EX_G(:)
        DEXC(3)   =      EX_G(1)+EX_G(2)
        DEXC(4)   =      EX_G(1)+EX_G(2)
        DEXC(5)   =2.D0*(EX_G(1)-EX_G(2))
      ELSE
        EXC       = EX(1)
        DEXC(1)   = EX_R(1)
        IF(.NOT.TGRA) RETURN
        DEXC(3)   = EX_G(1)
      END IF
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE$EVAL2(VAL,EXC,DEXC,D2EXC)
      USE EXCHANGE_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: VAL(5)     !RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS
      REAL(8),INTENT(OUT):: EXC        ! EXCHANGE ENERGY (RHOT*EPS_X)
      REAL(8),INTENT(OUT):: DEXC(5)    ! DEXC(I)=D[EXC]/D[VAL(I)]
      REAL(8),INTENT(OUT):: D2EXC(5,5) ! D2EXC(I,J)=D2[EXC]/(D[VAL(I)]D[VAL(J)])
      REAL(8)            :: RHO(2),GRHO2(2)
      REAL(8)            :: EX(2),EX_R(2),EX_G(2)
      REAL(8)            :: EX_RR(2),EX_RG(2),EX_GG(2)
      INTEGER(4)         :: I,J
      INTEGER(4)         :: NSPIN
!     ******************************************************************
      IF(.NOT.TINI) CALL EXCHANGE_INITIALIZE
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
!
!     ==========================================================================
!     == TRANSFORM INTO SPIN UP AND SPIN DOWN COMPONENTS                      ==
!     ==========================================================================
      IF(TSPIN) THEN
        RHO(1)=VAL(1)+VAL(2)               !=2*RHO(UP)
        RHO(2)=VAL(1)-VAL(2)               !=2*RHO(DOWN)
        GRHO2(1)=VAL(3)+2.D0*VAL(5)+VAL(4) !=2*GRAD[RHO(UP)]
        GRHO2(2)=VAL(3)-2.D0*VAL(5)+VAL(4) !=2*GRAD[RHO[DOWN]
        NSPIN=2
      ELSE
        RHO(1)=VAL(1)
        GRHO2(1)=VAL(3)
        NSPIN=1
      END IF
!
!     ==================================================================
!     == EXCHANGE_X2 ASSUMES A NON-SPIN POLARIZED TOTAL DENSITY       ==
!     == E_X(RHO_UP,RHO_DOWN)=0.5*(E_X(2*RHO_UP,0)+E_X(2*RHO_DOWN,0))/2=
!     ==================================================================
      DO I=1,NSPIN
        CALL EXCHANGE_X2(RHO(I),GRHO2(I),EX(I),EX_R(I),EX_G(I) &
     &            ,EX_RR(I),EX_RG(I),EX_GG(I))
      ENDDO
!
!     ==========================================================================
!     == MODIFY FOR XALPHA                                                    ==
!     ==========================================================================
      IF(TXALPHA) THEN
        EX(:)=1.5D0*XALPHA*EX(:)
        EX_R(:)=1.5D0*XALPHA*EX_R(:)
        EX_G(:)=1.5D0*XALPHA*EX_G(:)
        EX_RR(:)=1.5D0*XALPHA*EX_RR(:)
        EX_RG(:)=1.5D0*XALPHA*EX_RG(:)
        EX_GG(:)=1.5D0*XALPHA*EX_GG(:)
      END IF
!
!     ==========================================================================
!     == TRANSFORM BACK TO TOTAL AND SPIN DENSITIES                           ==
!     ==========================================================================
      IF(TSPIN) THEN
        EX(:)   =0.5D0*EX(:)
        EX_R(:) =0.5D0*EX_R(:)
        EX_RR(:)=0.5D0*EX_RR(:)
        EXC       =      EX(1)+EX(2)
        DEXC(1)   =      EX_R(1)+EX_R(2)
        DEXC(2)   =      EX_R(1)-EX_R(2)
        D2EXC(1,1)=      EX_RR(1)+EX_RR(2) 
        D2EXC(1,2)=      EX_RR(1)-EX_RR(2) 
        D2EXC(2,2)=      EX_RR(1)+EX_RR(2) 
        IF(.NOT.TGRA) THEN
          D2EXC(2,1)=D2EXC(1,2)
          RETURN
        END IF
        EX_G(:) =0.5D0*EX_G(:)
        EX_RG(:)=0.5D0*EX_RG(:)
        EX_GG(:)=0.5D0*EX_GG(:)
        DEXC(3)   =      EX_G(1)+EX_G(2)
        DEXC(4)   =      EX_G(1)+EX_G(2)
        DEXC(5)   =2.D0*(EX_G(1)-EX_G(2))
        D2EXC(1,3)=      EX_RG(1)+EX_RG(2) 
        D2EXC(1,4)=      EX_RG(1)+EX_RG(2) 
        D2EXC(1,5)=2.D0*(EX_RG(1)-EX_RG(2))
        D2EXC(2,3)=      EX_RG(1)-EX_RG(2) 
        D2EXC(2,4)=      EX_RG(1)-EX_RG(2) 
        D2EXC(2,5)=2.D0*(EX_RG(1)+EX_RG(2))
        D2EXC(3,3)=      EX_GG(1)+EX_GG(2) 
        D2EXC(3,4)=      EX_GG(1)+EX_GG(2) 
        D2EXC(3,5)=2.D0*(EX_GG(1)-EX_GG(2))
        D2EXC(4,4)=      EX_GG(1)+EX_GG(2) 
        D2EXC(4,5)=2.D0*(EX_GG(1)-EX_GG(2))
        D2EXC(5,5)=4.D0*(EX_GG(1)+EX_GG(2))
        DO I=1,4
          DO J=I+1,5
            D2EXC(J,I)=D2EXC(I,J)
          ENDDO
        ENDDO
      ELSE
        EXC       = EX(1)
        DEXC(1)   = EX_R(1)
        D2EXC(1,1)= EX_RR(1)
        IF(.NOT.TGRA) RETURN
        DEXC(3)   = EX_G(1)
        D2EXC(1,3)= EX_RG(1)
        D2EXC(3,3)= EX_GG(1)
        D2EXC(3,1)= D2EXC(1,3)
      END IF
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
      USE EXCHANGE_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: VAL(5)     !RHOT,RHOS,GRHOT2,GRHOS2,GRHOTS
      REAL(8),INTENT(OUT):: EXC        ! EXCHANGE ENERGY (RHOT*EPS_X)
      REAL(8),INTENT(OUT):: DEXC(5)    ! DEXC(I)=D[EXC]/D[VAL(I)]
      REAL(8),INTENT(OUT):: D2EXC(5,5) ! D2EXC(I,J)=D2[EXC]/(D[VAL(I)]D[VAL(J)])
      REAL(8),INTENT(OUT):: D3EXC(5,5,5)! D3EXC(I,J,K)=D2[EXC]/(D[VAL(I)]D[VAL(J)D[VAL[K]])
      REAL(8)            :: RHO(2),GRHO2(2)
      REAL(8)            :: EX(2),EX_R(2),EX_G(2)
      REAL(8)            :: EX_RR(2),EX_RG(2),EX_GG(2)
      REAL(8)            :: EX_RRR(2),EX_RRG(2),EX_RGG(2),EX_GGG(2)
      INTEGER(4)         :: I,J,K
      INTEGER(4)         :: NSPIN
!     ******************************************************************
      IF(.NOT.TINI) CALL EXCHANGE_INITIALIZE
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      D3EXC(:,:,:)=0.D0
!
!     ==================================================================
!     == TRANSFORM INTO SPIN UP AND SPIN DOWN COMPONENTS              ==
!     ==================================================================
      IF(TSPIN) THEN
        RHO(1)=VAL(1)+VAL(2)               !=2*RHO(UP)
        RHO(2)=VAL(1)-VAL(2)               !=2*RHO(DOWN)
        GRHO2(1)=VAL(3)+2.D0*VAL(5)+VAL(4) !=2*GRAD[RHO(UP)]
        GRHO2(2)=VAL(3)-2.D0*VAL(5)+VAL(4) !=2*GRAD[RHO[DOWN]
        NSPIN=2
      ELSE
        RHO(1)=VAL(1)
        GRHO2(1)=VAL(3)
        NSPIN=1
      END IF
!
!     ==================================================================
!     == EXCHANGE_X2 ASSUMES A NON-SPIN POLARIZED TOTAL DENSITY       ==
!     == E_X(RHO_UP,RHO_DOWN)=0.5*(E_X(2*RHO_UP,0)+E_X(2*RHO_DOWN,0))/2=
!     ==================================================================
      DO I=1,NSPIN
        CALL EXCHANGE_X3(RHO(I),GRHO2(I),EX(I),EX_R(I),EX_G(I) &
     &  ,EX_RR(I),EX_RG(I),EX_GG(I),EX_RRR(I),EX_RRG(I),EX_RGG(I),EX_GGG(I))
      ENDDO
!
!     ==================================================================
!     == MODIFY FOR XALPHA                                            ==
!     ==================================================================
      IF(TXALPHA) THEN
        EX(:)=1.5D0*XALPHA*EX(:)
        EX_R(:)=1.5D0*XALPHA*EX_R(:)
        EX_G(:)=1.5D0*XALPHA*EX_G(:)
        EX_RR(:)=1.5D0*XALPHA*EX_RR(:)
        EX_RG(:)=1.5D0*XALPHA*EX_RG(:)
        EX_GG(:)=1.5D0*XALPHA*EX_GG(:)
        EX_RRR(:)=1.5D0*XALPHA*EX_RRR(:)
        EX_RRG(:)=1.5D0*XALPHA*EX_RRG(:)
        EX_RGG(:)=1.5D0*XALPHA*EX_RGG(:)
        EX_GGG(:)=1.5D0*XALPHA*EX_GGG(:)
      END IF
!
!     ==================================================================
!     == TRANSFORM BACK TO TOTAL AND SPIN DENSITIES                   ==
!     ==================================================================
      IF(TSPIN) THEN
        EX(:)   =0.5D0*EX(:)
        EX_R(:) =0.5D0*EX_R(:)
        EX_RR(:)=0.5D0*EX_RR(:)
        EX_RRR(:)=0.5D0*EX_RRR(:)
        EXC       =      EX(1)+EX(2)
        DEXC(1)   =      EX_R(1)+EX_R(2)
        DEXC(2)   =      EX_R(1)-EX_R(2)
        D2EXC(1,1)=      EX_RR(1)+EX_RR(2) 
        D2EXC(1,2)=      EX_RR(1)-EX_RR(2) 
        D2EXC(2,2)=      EX_RR(1)+EX_RR(2) 
        D3EXC(1,1,1)=    EX_RRR(1)+EX_RRR(2) 
        D3EXC(1,1,2)=    EX_RRR(1)-EX_RRR(2) 
        D3EXC(1,2,2)=    EX_RRR(1)+EX_RRR(2) 
        D3EXC(2,2,2)=    EX_RRR(1)-EX_RRR(2) 
        IF(.NOT.TGRA) THEN
          D2EXC(2,1)=D2EXC(1,2)
          D3EXC(1,2,1)=D3EXC(1,1,2)
          D3EXC(2,1,1)=D3EXC(1,1,2)
          D3EXC(2,1,2)=D3EXC(1,2,2)
          D3EXC(2,2,1)=D3EXC(1,2,2)
          RETURN
        END IF
        EX_G(:)  =0.5D0*EX_G(:)
        EX_RG(:) =0.5D0*EX_RG(:)
        EX_GG(:) =0.5D0*EX_GG(:)
        EX_RRG(:)=0.5D0*EX_RRG(:)
        EX_RGG(:)=0.5D0*EX_RGG(:)
        EX_GGG(:)=0.5D0*EX_GGG(:)
        DEXC(3)   =      EX_G(1)+EX_G(2)
        DEXC(4)   =      EX_G(1)+EX_G(2)
        DEXC(5)   =2.D0*(EX_G(1)-EX_G(2))
        D2EXC(1,3)=      EX_RG(1)+EX_RG(2) 
        D2EXC(1,4)=      EX_RG(1)+EX_RG(2) 
        D2EXC(1,5)=2.D0*(EX_RG(1)-EX_RG(2))
        D2EXC(2,3)=      EX_RG(1)-EX_RG(2) 
        D2EXC(2,4)=      EX_RG(1)-EX_RG(2) 
        D2EXC(2,5)=2.D0*(EX_RG(1)+EX_RG(2))
        D2EXC(3,3)=      EX_GG(1)+EX_GG(2) 
        D2EXC(3,4)=      EX_GG(1)+EX_GG(2) 
        D2EXC(3,5)=2.D0*(EX_GG(1)-EX_GG(2))
        D2EXC(4,4)=      EX_GG(1)+EX_GG(2) 
        D2EXC(4,5)=2.D0*(EX_GG(1)-EX_GG(2))
        D2EXC(5,5)=4.D0*(EX_GG(1)+EX_GG(2))
        D3EXC(1,1,3)=     EX_RRG(1)+EX_RRG(2)
        D3EXC(1,1,4)=D3EXC(1,1,3)
        D3EXC(1,1,5)=2.D0*(EX_RRG(1)-EX_RRG(2))
        D3EXC(1,2,3)=     EX_RRG(1)-EX_RRG(2)
        D3EXC(1,2,4)=D3EXC(1,2,3)
        D3EXC(1,2,5)=2.D0*(EX_RRG(1)+EX_RRG(2))
        D3EXC(1,3,3)=      EX_RGG(1)+EX_RGG(2)
        D3EXC(1,3,4)=D3EXC(1,3,3)
        D3EXC(1,3,5)=2.D0*(EX_RGG(1)-EX_RGG(2))
        D3EXC(1,4,4)=D3EXC(1,3,3)
        D3EXC(1,4,5)=D3EXC(1,3,5)
        D3EXC(1,5,5)=4.D0*(EX_RGG(1)+EX_RGG(2))
        D3EXC(2,2,3)=      EX_RRG(1)+EX_RRG(2)
        D3EXC(2,2,4)=D3EXC(2,2,3)
        D3EXC(2,2,5)=2.D0*(EX_RRG(1)-EX_RRG(2))
        D3EXC(2,3,3)=      EX_RGG(1)-EX_RGG(2)
        D3EXC(2,3,4)=D3EXC(2,3,3)
        D3EXC(2,3,5)=2.D0*(EX_RGG(1)+EX_RGG(2))
        D3EXC(2,4,4)=D3EXC(2,3,4)
        D3EXC(2,4,5)=D3EXC(2,3,5)
        D3EXC(2,5,5)=4.D0*(EX_RGG(1)-EX_RGG(2))
        D3EXC(3,3,3)=      EX_GGG(1)+EX_GGG(2)
        D3EXC(3,3,4)=D3EXC(3,3,3)
        D3EXC(3,3,5)=2.D0*(EX_GGG(1)-EX_GGG(2))
        D3EXC(3,4,4)=D3EXC(3,3,3)
        D3EXC(3,4,5)=D3EXC(3,3,5)
        D3EXC(3,5,5)=4.D0*(EX_GGG(1)+EX_GGG(2))
        D3EXC(4,4,4)=D3EXC(3,3,3)
        D3EXC(4,4,5)=D3EXC(3,3,5)
        D3EXC(4,5,5)=D3EXC(3,5,5)
        D3EXC(5,5,5)=8.D0*(EX_GGG(1)-EX_GGG(2))
        DO I=1,4
          DO J=I+1,5
            D2EXC(J,I)=D2EXC(I,J)
            D3EXC(I,J,I)=D3EXC(I,I,J)
            D3EXC(J,I,I)=D3EXC(I,I,J)
            D3EXC(J,I,J)=D3EXC(I,J,J)
            D3EXC(J,J,I)=D3EXC(I,J,J)
            DO K=J+1,5
              D3EXC(I,K,J)=D3EXC(I,J,K)
              D3EXC(J,I,K)=D3EXC(I,J,K)
              D3EXC(J,K,I)=D3EXC(I,J,K)
              D3EXC(K,I,J)=D3EXC(I,J,K)
              D3EXC(K,J,I)=D3EXC(I,J,K)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        EXC       = EX(1)
        DEXC(1)   = EX_R(1)
        D2EXC(1,1)= EX_RR(1)
        D3EXC(1,1,1)= EX_RRR(1)
        IF(.NOT.TGRA) RETURN
        DEXC(3)   = EX_G(1)
        D2EXC(1,3)= EX_RG(1)
        D2EXC(3,3)= EX_GG(1)
        D2EXC(3,1)= D2EXC(1,3)
        D3EXC(1,1,3)=EX_RRG(1)
        D3EXC(1,3,3)=EX_RGG(1)
        D3EXC(3,3,3)=EX_GGG(1)
        D3EXC(1,3,1)=D3EXC(1,1,3)
        D3EXC(3,1,1)=D3EXC(1,1,3)
        D3EXC(3,1,3)=D3EXC(1,3,3)
        D3EXC(3,3,1)=D3EXC(1,3,3)
      END IF
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE_X1(RHO,GRHO2,EX,EX_R,EX_G)
!     *****************************************************************
!     **                                                             **
!     **  EXCHANGE ENERGY OF A NON-SPIN POLARIZED ELECTRON GAS       **
!     **                                                             **
!     **  THE PBE EXCHANGE HAS BEEN MODIFIED IN ORDER TO SATISFY     **
!     **  THE CORRECT 1/2R BEHAVIOR OF THE EXCHANGE ENERGY PER       **
!     **  PARTICLE IN THE TAIL REGION OF THE WAVE FUNCTIONS.          **
!     **  IT DEVIATES FROM THE PBE EXCHANGE ONLY IN ORDERS X**6      **
!     **  AND HIGHER, WHERE X IS THE DIMENSIONLESS GRADIENT          **
!     **  GRHO/RHO**(4/3).                                            **
!     **                                                             **
!     *****************************************************************
!     *****************************************************************
      USE TABLE1D_MODULE
      USE EXCHANGE_MODULE, ONLY : EXFAC,S2FAC,TGRA &
     &             ,TSAFE,S2TABLE,NS2,FXARRAY,DFXARRAY
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: RHO
      REAL(8),INTENT(IN) :: GRHO2      ! GRAD(RHO)**2
      REAL(8),INTENT(OUT):: EX
      REAL(8),INTENT(OUT):: EX_R,EX_G 
!     
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER  :: ONEBY6=1.D0/6.D0
      REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
      REAL(8),PARAMETER  :: EIGHTBY3=8.D0/3.D0
      REAL(8),PARAMETER  :: ELEVENBY3=11.D0/3.D0
      REAL(8)            :: RHO13,RHO23,RHO43
      REAL(8)            :: SVAR
      REAL(8)            :: S2,S2_R,S2_G 
      REAL(8)            :: FX,FX_S,FX_R,FX_G 
      REAL(8)            :: FX_SS
      LOGICAL(4)         :: TTABLE
!     *****************************************************************
!     
!     =================================================================
!     == EXIT FOR SMALL DENSITIES                                    ==
!     =================================================================
      EX=0.D0
      EX_R=0.D0
      EX_G=0.D0
      IF(ABS(RHO).LT.1.D-50) RETURN  !#
!     
!     =================================================================
!     == POWERS OF RHO**(1/3)                                        ==
!     =================================================================
      RHO13=RHO**ONEBY3  
      RHO23=RHO13*RHO13
      RHO43=RHO23*RHO23
!     
!     =================================================================
!     == LOCAL EXCHANGE                                              ==
!     =================================================================
      EX=EXFAC*RHO43
      EX_R=EX*FOURBY3/RHO        !#
!     
!     =================================================================
!     ==  RETURN POINT FOR LOCAL EXCHANGE                            ==
!     =================================================================
      IF(.NOT.TGRA) THEN
        EX_G =0.D0
        RETURN
      END IF
!     
!     =================================================================
!     ==  CALCULATE REDUCED GRADIENT S2=(0.5*GRHO/KFRHO)**2          ==
!     =================================================================
      SVAR=S2FAC/RHO43**2     !PROPORTIONAL RHO**(-8/3)
      S2=SVAR*GRHO2           
      S2_R=-EIGHTBY3*S2/RHO   !PROPORTIONAL RHO**(-11/3)
      S2_G=SVAR
!     
!     
!     =================================================================
!     ==  CALCULATE GRADIENT ENHANCEMENT FACTOR FX                   ==
!     ==  (TABLE LOOKUP IF PERMITTED BY TSAFE)                       ==
!     =================================================================
      TTABLE=.NOT.TSAFE
      IF(TTABLE) CALL TABLE1D$SETX(S2TABLE,S2,TTABLE)
      IF(TTABLE) THEN
        CALL TABLE1D$GETF(S2TABLE,NS2,FXARRAY,FX)
        CALL TABLE1D$GETF(S2TABLE,NS2,DFXARRAY,FX_S)
      ELSE
        CALL EXCHANGE_FX2(S2,FX,FX_S,FX_SS)
      END IF
      FX_R=FX_S*S2_R
      FX_G=FX_S*S2_G
!     
!     =================================================================
!     == MULTIPLY ENHANCEMENT FACTOR WITH LOCAL EXCHANGE             ==
!     == WATCH ORDER OF STATEMENTS!!                                 ==
!     =================================================================
      EX_R =EX_R*FX+EX*FX_R
      EX_G =EX*FX_G
      EX   =EX*FX
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE_X2(RHO,GRHO2,EX,EX_R,EX_G,EX_RR,EX_RG,EX_GG)
!     *****************************************************************
!     **                                                             **
!     **  EXCHANGE ENERGY OF A NON-SPIN POLARIZED ELECTRON GAS       **
!     **                                                             **
!     **  THE PBE EXCHANGE HAS BEEN MODIFIED IN ORDER TO SATISFY     **
!     **  THE CORRECT 1/2R BEHAVIOR OF THE EXCHANGE ENERGY PER       **
!     **  PARTICLE IN THE TAIL REGION OF THE WAVE FUNCTIONS.          **
!     **  IT DEVIATES FROM THE PBE EXCHANGE ONLY IN ORDERS X**6      **
!     **  AND HIGHER, WHERE X IS THE DIMENSIONLESS GRADIENT          **
!     **  GRHO/RHO**(4/3).                                            **
!     **                                                             **
!     *****************************************************************
!     *****************************************************************
      USE TABLE1D_MODULE
      USE EXCHANGE_MODULE, ONLY : EXFAC,S2FAC,TGRA &
     &             ,TSAFE,S2TABLE,NS2,FXARRAY,DFXARRAY,D2FXARRAY 
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: RHO
      REAL(8),INTENT(IN) :: GRHO2      ! GRAD(RHO)**2
      REAL(8),INTENT(OUT):: EX
      REAL(8),INTENT(OUT):: EX_R,EX_G 
      REAL(8),INTENT(OUT):: EX_RR,EX_RG,EX_GG
!     
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER  :: ONEBY6=1.D0/6.D0
      REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
      REAL(8),PARAMETER  :: FOURBY9=4.D0/9.D0
      REAL(8),PARAMETER  :: EIGHTBY3=8.D0/3.D0
      REAL(8),PARAMETER  :: ELEVENBY3=11.D0/3.D0
      REAL(8)            :: RHO13,RHO23,RHO43
      REAL(8)            :: SVAR
      REAL(8)            :: S2,S2_R,S2_G 
      REAL(8)            :: S2_RR,S2_RG
      REAL(8)            :: FX,FX_S,FX_R,FX_G 
      REAL(8)            :: FX_SS,FX_RR,FX_RG,FX_GG
      LOGICAL(4)         :: TTABLE
!     *****************************************************************
!     
!     =================================================================
!     == EXIT FOR SMALL DENSITIES                                   ==
!     =================================================================
      IF(RHO.LT.1.D-10) THEN
        EX=0.D0
        EX_R=0.D0
        EX_G=0.D0
        EX_RR=0.D0
        EX_RG=0.D0
        EX_GG=0.D0
      END IF
!     
!     =================================================================
!     == POWERS OF RHO**(1/3)                                        ==
!     =================================================================
      RHO13=RHO**ONEBY3
      RHO23=RHO13*RHO13
      RHO43=RHO23*RHO23
!     
!     =================================================================
!     == LOCAL EXCHANGE                                              ==
!     =================================================================
      EX=EXFAC*RHO43
      EX_R=EXFAC*FOURBY3*RHO13
      EX_RR=EXFAC*FOURBY9/RHO23
!     
!     =================================================================
!     ==  RETURN POINT FOR LOCAL EXCHANGE                            ==
!     =================================================================
      IF(.NOT.TGRA) THEN
        EX_G =0.D0
        EX_RG=0.D0
        EX_GG=0.D0
        RETURN
      END IF
!     
!     =================================================================
!     ==  CALCULATE REDUCED GRADIENT S2=(0.5*GRHO/KFRHO)**2         ==
!     =================================================================
      SVAR=S2FAC/RHO43**2     !PROPORTIONAL RHO**(-8/3)
      S2=SVAR*GRHO2           
      S2_R=-EIGHTBY3*S2/RHO   !PROPORTIONAL RHO**(-11/3)
      S2_G=SVAR
      S2_RR=-ELEVENBY3*S2_R/RHO  !PROPORTIONAL RHO**(-14/3)
      S2_RG=-EIGHTBY3*SVAR/RHO
!     
!     
!     =================================================================
!     ==  CALCULATE GRADIENT ENHANCEMENT FACTOR FX                   ==
!     ==  (TABLE LOOKUP IF PERMITTED BY TSAFE)                       ==
!     =================================================================
      TTABLE=.NOT.TSAFE
      IF(TTABLE) CALL TABLE1D$SETX(S2TABLE,S2,TTABLE)
      IF(TTABLE) THEN
        CALL TABLE1D$GETF(S2TABLE,NS2,FXARRAY,FX)
        CALL TABLE1D$GETF(S2TABLE,NS2,DFXARRAY,FX_S)
        CALL TABLE1D$GETF(S2TABLE,NS2,D2FXARRAY,FX_SS)
      ELSE
        CALL EXCHANGE_FX2(S2,FX,FX_S,FX_SS)
      END IF
      FX_R=FX_S*S2_R
      FX_G=FX_S*S2_G
      FX_RR=FX_SS*S2_R*S2_R+FX_S*S2_RR
      FX_RG=FX_SS*S2_R*S2_G+FX_S*S2_RG
      FX_GG=FX_SS*S2_G*S2_G
!     
!     =================================================================
!     == MULTIPLY ENHANCEMENT FACTOR WITH LOCAL EXCHANGE             ==
!     == WATCH ORDER OF STATEMENTS!!                                 ==
!     =================================================================
      EX_RR=EX_RR*FX+2.D0*EX_R*FX_R+EX*FX_RR
      EX_RG=EX_R*FX_G+EX*FX_RG
      EX_GG=EX*FX_GG 
      EX_R =EX_R*FX+EX*FX_R
      EX_G =EX*FX_G
      EX   =EX*FX
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE_X3(RHO,GRHO2,EX,EX_R,EX_G,EX_RR,EX_RG,EX_GG &
     &           ,EX_RRR,EX_RRG,EX_RGG,EX_GGG)
!     *****************************************************************
!     **                                                             **
!     **  EXCHANGE ENERGY OF A NON-SPIN POLARIZED ELECTRON GAS       **
!     **                                                             **
!     **  THE PBE EXCHANGE HAS BEEN MODIFIED IN ORDER TO SATISFY     **
!     **  THE CORRECT 1/2R BEHAVIOR OF THE EXCHANGE ENERGY PER       **
!     **  PARTICLE IN THE TAIL REGION OF THE WAVE FUNCTIONS.          **
!     **  IT DEVIATES FROM THE PBE EXCHANGE ONLY IN ORDERS X**6      **
!     **  AND HIGHER, WHERE X IS THE DIMENSIONLESS GRADIENT          **
!     **  GRHO/RHO**(4/3).                                            **
!     **                                                             **
!     *****************************************************************
!     *****************************************************************
      USE TABLE1D_MODULE
      USE EXCHANGE_MODULE, ONLY : EXFAC,S2FAC,TGRA &
     &    ,TSAFE,S2TABLE,NS2,FXARRAY,DFXARRAY,D2FXARRAY,D3FXARRAY
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: RHO
      REAL(8),INTENT(IN) :: GRHO2      ! GRAD(RHO)**2
      REAL(8),INTENT(OUT):: EX
      REAL(8),INTENT(OUT):: EX_R,EX_G 
      REAL(8),INTENT(OUT):: EX_RR,EX_RG,EX_GG
      REAL(8),INTENT(OUT):: EX_RRR,EX_RRG,EX_RGG,EX_GGG
!     
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER  :: TWOBY3=2.D0/3.D0
      REAL(8),PARAMETER  :: ONEBY6=1.D0/6.D0
      REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
      REAL(8),PARAMETER  :: FOURBY9=4.D0/9.D0
      REAL(8),PARAMETER  :: EIGHTBY3=8.D0/3.D0
      REAL(8),PARAMETER  :: ELEVENBY3=11.D0/3.D0
      REAL(8)            :: RHO13,RHO23,RHO43
      REAL(8)            :: SVAR
      REAL(8)            :: S2,S2_R,S2_G 
      REAL(8)            :: S2_RR,S2_RG
      REAL(8)            :: S2_RRR,S2_RRG
      REAL(8)            :: FX,FX_S,FX_R,FX_G 
      REAL(8)            :: FX_SS,FX_RR,FX_RG,FX_GG
      REAL(8)            :: FX_SSS,FX_RRR,FX_RRG,FX_RGG,FX_GGG
      LOGICAL(4)         :: TTABLE
!     *****************************************************************
!     
!     =================================================================
!     == EXIT FOR SMALL DENSITIES                                   ==
!     =================================================================
      IF(RHO.LT.1.D-10) THEN
        EX=0.D0
        EX_R=0.D0
        EX_G=0.D0
        EX_RR=0.D0
        EX_RG=0.D0
        EX_GG=0.D0
        EX_RRR=0.D0
        EX_RRG=0.D0
        EX_RGG=0.D0
        EX_GGG=0.D0
      END IF
!     
!     =================================================================
!     == POWERS OF RHO**(1/3)                                        ==
!     =================================================================
      RHO13=RHO**ONEBY3
      RHO23=RHO13*RHO13
      RHO43=RHO23*RHO23
!     
!     =================================================================
!     == LOCAL EXCHANGE                                              ==
!     =================================================================
      EX=EXFAC*RHO43
      EX_R=EXFAC*FOURBY3*RHO13
      EX_RR=EXFAC*FOURBY9/RHO23
      EX_RRR=-TWOBY3*EX_RR/RHO
!     
!     =================================================================
!     ==  RETURN POINT FOR LOCAL EXCHANGE                            ==
!     =================================================================
      IF(.NOT.TGRA) THEN
        EX_G =0.D0
        EX_RG=0.D0
        EX_GG=0.D0
        EX_RRG=0.D0
        EX_RGG=0.D0
        EX_GGG=0.D0
        RETURN
      END IF
!     
!     =================================================================
!     ==  CALCULATE REDUCED GRADIENT S2=(0.5*GRHO/KFRHO)**2         ==
!     =================================================================
      SVAR=S2FAC/RHO43**2     !PROPORTIONAL RHO**(-8/3)
      S2=SVAR*GRHO2           
      S2_R=-EIGHTBY3*S2/RHO   !PROPORTIONAL RHO**(-11/3)
      S2_G=SVAR
      S2_RR=-ELEVENBY3*S2_R/RHO  !PROPORTIONAL RHO**(-14/3)
      S2_RG=-EIGHTBY3*SVAR/RHO
      S2_RRR=-(14.D0/3.D0)*S2_RR/RHO
      S2_RRG=ELEVENBY3*EIGHTBY3*SVAR/(RHO**2)
!     
!     =================================================================
!     ==  CALCULATE GRADIENT ENHANCEMENT FACTOR FX                   ==
!     ==  (TABLE LOOKUP IF PERMITTED BY TSAFE)                       ==
!     =================================================================
      TTABLE=.NOT.TSAFE
      IF(TTABLE) CALL TABLE1D$SETX(S2TABLE,S2,TTABLE)
      IF(TTABLE) THEN
        CALL TABLE1D$GETF(S2TABLE,NS2,FXARRAY,FX)
        CALL TABLE1D$GETF(S2TABLE,NS2,DFXARRAY,FX_S)
        CALL TABLE1D$GETF(S2TABLE,NS2,D2FXARRAY,FX_SS)
        CALL TABLE1D$GETF(S2TABLE,NS2,D3FXARRAY,FX_SSS)
      ELSE
        CALL EXCHANGE_FX3(S2,FX,FX_S,FX_SS,FX_SSS)
      END IF
      FX_R   = FX_S*S2_R
      FX_G   = FX_S*S2_G
      FX_RR  = FX_SS*S2_R*S2_R + FX_S*S2_RR
      FX_RG  = FX_SS*S2_R*S2_G + FX_S*S2_RG
      FX_GG  = FX_SS*S2_G*S2_G
      FX_RRR = FX_SSS*S2_R**3 + 3.D0*FX_SS*S2_R*S2_RR + FX_S*S2_RRR
      FX_RRG = FX_SSS*S2_R**2*S2_G + 2.D0*FX_SS*S2_R*S2_RG + FX_SS*S2_G*S2_RR +FX_S*S2_RRG
      FX_RGG = FX_SSS*S2_G**2*S2_R + 2.D0*FX_SS*S2_G*S2_RG
      FX_GGG = FX_SSS*S2_G**3 
!     
!     =================================================================
!     == MULTIPLY ENHANCEMENT FACTOR WITH LOCAL EXCHANGE             ==
!     == WATCH ORDER OF STATEMENTS!!                                 ==
!     =================================================================
      EX_RRR=EX_RRR*FX + 3.D0*EX_RR*FX_R + 3.D0*EX_R*FX_RR + EX*FX_RRR      
      EX_RRG=EX_RR*FX_G + 2.D0*EX_R*FX_RG + EX*FX_RRG      
      EX_RGG=EX_R*FX_GG + EX*FX_RGG
      EX_GGG=EX*FX_GGG
      EX_RR=EX_RR*FX + 2.D0*EX_R*FX_R + EX*FX_RR
      EX_RG=EX_R*FX_G + EX*FX_RG
      EX_GG=EX*FX_GG 
      EX_R =EX_R*FX + EX*FX_R
      EX_G =EX*FX_G
      EX   =EX*FX
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE_FX1(S2,FX,FX_S)
      USE EXCHANGE_MODULE,ONLY : FXTYPE,MU,KAPPA,MUBYKAPPA,ALPHAB2
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: S2
      REAL(8),INTENT(OUT):: FX
      REAL(8),INTENT(OUT):: FX_S
      REAL(8),PARAMETER  :: ONEBY6=1.D0/6.D0
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8)            :: BFAC,BFAC2,BFAC_S
      REAL(8)            :: CS2,CS,CS22,CS23,CS24
      REAL(8)            :: CSH    !
      REAL(8)            :: ASNH   !ARCUS SINUS HYPERBOLICUS (CS)
      REAL(8)            :: SASNH,SASNH_S
!     *****************************************************************
      IF(FXTYPE.EQ.'PBE') THEN
!       ===============================================================
!       == THIS IS THE EXCHANGE AS USED IN PBES PAPER                ==
!       == THE SIMPLER FORM                                          ==
!       == THAT GOES TO A CONSTANT FOR S2->INFTY                     ==
!       == IT REQUIRES THE FACTORS MU,KAPPA,MUBYKAPPA                ==
!       ===============================================================
        BFAC=1.D0/(1.D0+MUBYKAPPA*S2)
        BFAC2=BFAC*BFAC
        FX=1.D0+KAPPA*(1.D0-BFAC)
        FX_S=MU*BFAC2
      ELSE IF(FXTYPE.EQ.'RPBE') THEN
!       ===============================================================
!       == THIS IS THE RPBE FORM OF EXCHANGE AS USED IN EQ. 15 OF    ==
!       == B. HAMMER, L.B. HANSEN AND J.K. NORSKOV, PRB59,7413(1999) ==
!       ===============================================================
        BFAC=EXP(-MUBYKAPPA*S2)
        FX=1.D0+KAPPA*(1.D0-BFAC)
        FX_S  =MU*BFAC
      ELSE IF(FXTYPE.EQ.'BECKE') THEN
!       ===============================================================
!       == THIS IS THE EXCHANGE AS USED IN BECKES PAPER              ==
!       == IT CONTAINS THE ARCUS SINUS HYPERBOLICUS                  ==
!       == AND FULFILLS E_X(R)->0.5*RHO/R FOR THE EXPONENTIAL        ==
!       == TAIL OF THE WAVE FUNCTION RHO(R)=A*EXP(-B*R)              ==
!       == SASNH = ALPHAB*S*ASNH(ALPHAB*S)/ALPHAB**2  REPLACES S2    ==
!       == IN THE DENOMINATOR OF FX=1+KAPPA*MU*S2/(1+MU/KAPPA*S2)    ==
!       == SASNH(S->0)    =S2;     SASNH(S->INFTY)=S*LN(S)/ALPHAB    ==
!       ===============================================================
        CS2=ALPHAB2*S2
        IF(CS2.GT.1.D-20) THEN
          CS=SQRT(CS2)
          CSH=SQRT(1.D0+CS2)
          ASNH=LOG(CS+CSH)    ! =S+O(S**3)
          SASNH   =CS*ASNH     ! 
          SASNH_S =0.5D0*(ASNH/CS+1.D0/CSH)
          SASNH=SASNH/ALPHAB2
        ELSE
          CS22=CS2*CS2
          CS23=CS22*CS2
          CS24=CS23*CS2
          SASNH=CS2-ONEBY6*CS22+0.075D0*CS23-5.D0/112.D0*CS24
          SASNH_S=1.D0-ONEBY3*CS2+0.225D0*CS22-5.D0/28.D0*CS23
          SASNH=SASNH/ALPHAB2
        END IF
        SASNH   =MUBYKAPPA*SASNH
        SASNH_S =MUBYKAPPA*SASNH_S
        BFAC  =1.D0/(1.D0+SASNH)
        BFAC2 =BFAC*BFAC
        BFAC_S=-BFAC2*SASNH_S
        FX=1.D0+MU*S2*BFAC
        FX_S=MU*(BFAC+S2*BFAC_S)
      ELSE
        CALL ERROR$MSG('FXTYPE NOT RECOGNIZED (ALLOWED VALUES ARE "PBE", "RPBE", "BECKE")')
        CALL ERROR$CHVAL('FXTYPE',FXTYPE)
        CALL ERROR$STOP('EXCHANGE_FX1')
      END IF
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE_FX2(S2,FX,FX_S,FX_SS)
      USE EXCHANGE_MODULE,ONLY : FXTYPE,MU,KAPPA,MUBYKAPPA,ALPHAB2
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: S2
      REAL(8),INTENT(OUT):: FX
      REAL(8),INTENT(OUT):: FX_S
      REAL(8),INTENT(OUT):: FX_SS
      REAL(8),PARAMETER  :: ONEBY6=1.D0/6.D0
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8)            :: BFAC,BFAC2,BFAC_S,BFAC_SS
      REAL(8)            :: CS2,CS,CS22,CS23,CS24
      REAL(8)            :: CSH    !
      REAL(8)            :: ASNH   !ARCUS SINUS HYPERBOLICUS (CS)
      REAL(8)            :: SASNH,SASNH_S,SASNH_SS
!     *****************************************************************
      IF(FXTYPE.EQ.'PBE') THEN
!       ===============================================================
!       == THIS IS THE EXCHANGE AS USED IN PBES PAPER                ==
!       == THE SIMPLER FORM                                          ==
!       == THAT GOES TO A CONSTANT FOR S2->INFTY                     ==
!       == IT REQUIRES THE FACTORS MU,KAPPA,MUBYKAPPA                ==
!       ===============================================================
        BFAC=1.D0/(1.D0+MUBYKAPPA*S2)
        BFAC2=BFAC*BFAC
        FX=1.D0+KAPPA*(1.D0-BFAC)
        FX_S=MU*BFAC2
        FX_SS=-2.D0*MU*BFAC2*BFAC*MUBYKAPPA
      ELSE IF(FXTYPE.EQ.'RPBE') THEN
!       ===============================================================
!       == THIS IS THE RPBE FORM OF EXCHANGE AS USED IN EQ. 15 OF    ==
!       == B. HAMMER, L.B. HANSEN AND J.K. NORSKOV, PRB59,7413(1999) ==
!       ===============================================================
        BFAC=EXP(-MUBYKAPPA*S2)
        FX=1.D0+KAPPA*(1.D0-BFAC)
        FX_S  =MU*BFAC
        FX_SS =-MUBYKAPPA*FX_S
      ELSE IF(FXTYPE.EQ.'BECKE') THEN
!       ===============================================================
!       == THIS IS THE EXCHANGE AS USED IN BECKES PAPER              ==
!       == IT CONTAINS THE ARCUS SINUS HYPERBOLICUS                  ==
!       == AND FULFILLS E_X(R)->0.5*RHO/R FOR THE EXPONENTIAL        ==
!       == TAIL OF THE WAVE FUNCTION RHO(R)=A*EXP(-B*R)              ==
!       == SASNH = ALPHAB*S*ASNH(ALPHAB*S)/ALPHAB**2  REPLACES S2    ==
!       == IN THE DENOMINATOR OF FX=1+KAPPA*MU*S2/(1+MU/KAPPA*S2)    ==
!       == SASNH(S->0)    =S2;     SASNH(S->INFTY)=S*LN(S)/ALPHAB    ==
!       ===============================================================
        CS2=ALPHAB2*S2
        IF(CS2.GT.1.D-3) THEN
          CS=SQRT(CS2)
          CSH=SQRT(1.D0+CS2)
          ASNH=LOG(CS+CSH)    ! =S+O(S**3)
          SASNH   =CS*ASNH     ! 
          SASNH_S =0.5D0*(ASNH/CS+1.D0/CSH)
          SASNH_SS=0.25D0/CS2*(1.D0/CSH**3-ASNH/CS) 
          SASNH=SASNH/ALPHAB2
          SASNH_SS=SASNH_SS*ALPHAB2
        ELSE
          CS22=CS2*CS2
          CS23=CS22*CS2
          CS24=CS23*CS2
          SASNH=CS2-ONEBY6*CS22+0.075D0*CS23-5.D0/112.D0*CS24
          SASNH_S=1.D0-ONEBY3*CS2+0.225D0*CS22-5.D0/28.D0*CS23
          SASNH_SS=-ONEBY3+0.45D0*CS2-15.D0/28.D0*CS22
          SASNH=SASNH/ALPHAB2
          SASNH_SS=SASNH_SS*ALPHAB2
        END IF
        SASNH   =MUBYKAPPA*SASNH
        SASNH_S =MUBYKAPPA*SASNH_S
        SASNH_SS=MUBYKAPPA*SASNH_SS
        BFAC  =1.D0/(1.D0+SASNH)
        BFAC2 =BFAC*BFAC
        BFAC_S=-BFAC2*SASNH_S
        BFAC_SS=BFAC2*(2.D0*BFAC*SASNH_S**2-SASNH_SS)
        FX=1.D0+MU*S2*BFAC
        FX_S=MU*(BFAC+S2*BFAC_S)
        FX_SS=MU*(2.D0*BFAC_S+S2*BFAC_SS)
      ELSE
        CALL ERROR$MSG('FXTYPE NOT RECOGNIZED (ALLOWED VALUES ARE "PBE", "RPBE", "BECKE")')
        CALL ERROR$CHVAL('FXTYPE',FXTYPE)
        CALL ERROR$STOP('EXCHANGE_FX2')
      END IF
      RETURN
      END
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXCHANGE_FX3(S2,FX,FX_S,FX_SS,FX_SSS)
      USE EXCHANGE_MODULE,ONLY : FXTYPE,MU,KAPPA,MUBYKAPPA,ALPHAB2
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: S2
      REAL(8),INTENT(OUT):: FX
      REAL(8),INTENT(OUT):: FX_S
      REAL(8),INTENT(OUT):: FX_SS
      REAL(8),INTENT(OUT):: FX_SSS
      REAL(8),PARAMETER  :: ONEBY6=1.D0/6.D0
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8)            :: BFAC,BFAC2,BFAC3,BFAC4,BFAC_S,BFAC_SS,BFAC_SSS
      REAL(8)            :: CS2,CS,CS22,CS23,CS24
      REAL(8)            :: CSH    !
      REAL(8)            :: ASNH   !ARCUS SINUS HYPERBOLICUS (CS)
      REAL(8)            :: SASNH,SASNH_S,SASNH_SS,SASNH_SSS
!     *****************************************************************
      IF(FXTYPE.EQ.'PBE') THEN
!       ===============================================================
!       == THIS IS THE EXCHANGE AS USED IN PBES PAPER                ==
!       == THE SIMPLER FORM                                          ==
!       == THAT GOES TO A CONSTANT FOR S2->INFTY                     ==
!       == IT REQUIRES THE FACTORS MU,KAPPA,MUBYKAPPA                ==
!       ===============================================================
        BFAC=1.D0/(1.D0+MUBYKAPPA*S2)
        BFAC2=BFAC*BFAC
        FX=1.D0+KAPPA*(1.D0-BFAC)
        FX_S  =MU*BFAC2
        FX_SS =-2.D0*MUBYKAPPA*BFAC*FX_S
        FX_SSS=-3.D0*MUBYKAPPA*BFAC*FX_SS
      ELSE IF(FXTYPE.EQ.'RPBE') THEN
!       ===============================================================
!       == THIS IS THE RPBE FORM OF EXCHANGE AS USED IN EQ. 15 OF    ==
!       == B. HAMMER, L.B. HANSEN AND J.K. NORSKOV, PRB59,7413(1999) ==
!       ===============================================================
        BFAC=EXP(-MUBYKAPPA*S2)
        FX=1.D0+KAPPA*(1.D0-BFAC)
        FX_S  =MU*BFAC
        FX_SS =-MUBYKAPPA*FX_S
        FX_SSS=-MUBYKAPPA*FX_SS
      ELSE IF(FXTYPE.EQ.'BECKE') THEN
!       ===============================================================
!       == THIS IS THE EXCHANGE AS USED IN BECKES PAPER              ==
!       == IT CONTAINS THE ARCUS SINUS HYPERBOLICUS                  ==
!       == AND FULFILLS E_X(R)->0.5*RHO/R FOR THE EXPONENTIAL        ==
!       == TAIL OF THE WAVE FUNCTION RHO(R)=A*EXP(-B*R)              ==
!       == SASNH = ALPHAB*S*ASNH(ALPHAB*S)/ALPHAB**2  REPLACES S2    ==
!       == IN THE DENOMINATOR OF FX=1+KAPPA*MU*S2/(1+MU/KAPPA*S2)    ==
!       == SASNH(S->0)    =S2;     SASNH(S->INFTY)=S*LN(S)/ALPHAB    ==
!       ===============================================================
        CS2=ALPHAB2*S2
       IF(CS2.GT.1.D-3) THEN
          CS=SQRT(CS2)
          CSH=SQRT(1.D0+CS2)
          ASNH=LOG(CS+CSH)    ! =S+O(S**3)
          SASNH   =CS*ASNH     ! 
          SASNH_S =0.5D0*(ASNH/CS+1.D0/CSH)
          SASNH_SS=0.25D0/CS2*(1.D0/CSH**3-ASNH/CS) 
          SASNH_SSS=-.125D0/CS2**2*((3D0+7.D0*CS2+CS2**2)/CSH**5-3.D0*ASNH/CS)
          SASNH=SASNH/ALPHAB2
          SASNH_SS=SASNH_SS*ALPHAB2
          SASNH_SSS=SASNH_SSS*ALPHAB2**2
        ELSE
          CS22=CS2*CS2
          CS23=CS22*CS2
          CS24=CS23*CS2
          SASNH    =CS2 -ONEBY6*CS22+0.075D0*CS23-5.D0/112.D0*CS24
          SASNH_S  =1.D0-ONEBY3*CS2 +0.225D0*CS22- 5.D0/28.D0*CS23
          SASNH_SS =    -ONEBY3     +0.45D0 *CS2 -15.D0/28.D0*CS22
          SASNH_SSS=                +0.45D0      -15.D0/14.D0*CS2
          SASNH=SASNH/ALPHAB2
          SASNH_SS=SASNH_SS*ALPHAB2
          SASNH_SSS=SASNH_SSS*ALPHAB2**2
        END IF
        SASNH    =MUBYKAPPA*SASNH
        SASNH_S  =MUBYKAPPA*SASNH_S
        SASNH_SS =MUBYKAPPA*SASNH_SS
        SASNH_SSS=MUBYKAPPA*SASNH_SSS
        BFAC  =1.D0/(1.D0+SASNH)
        BFAC2 =BFAC*BFAC
        BFAC3 =BFAC2*BFAC
        BFAC4 =BFAC2*BFAC2
        BFAC_S=-BFAC2*SASNH_S
        BFAC_SS=BFAC2*(2.D0*BFAC*SASNH_S**2-SASNH_SS)
        BFAC_SSS=-6.D0*BFAC4*SASNH_S**3+6.D0*BFAC3*SASNH_S*SASNH_SS &
       &         -BFAC2*SASNH_SSS
        FX    =1.D0+MU*S2*BFAC
        FX_S  =MU*(BFAC+S2*BFAC_S)
        FX_SS =MU*(2.D0*BFAC_S+S2*BFAC_SS)
        FX_SSS=MU*(3.D0*BFAC_SS+S2*(BFAC_SSS))
      ELSE
        CALL ERROR$MSG('FXTYPE NOT RECOGNIZED')
        CALL ERROR$MSG('(ALLOWED VALUES ARE "PBE", "RPBE", "BECKE")')
        CALL ERROR$CHVAL('FXTYPE',FXTYPE)
        CALL ERROR$STOP('EXCHANGE_FX3')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
MODULE PERDEWZUNGER_MODULE
!***********************************************************************
!**                                                                   **
!**  PURPOSE: EVALUATE THE PERDEW ZUNGER PARAMETERIZATION             **
!**    OF THE CORRELATION ENERGY OF A HOMOGENEOUS ELECTRON GAS        **
!**    AS OBTAINED FROM CEPERLEY ALDERS MONTE CARLO SIMULATIONS       **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    PERDEWZUNGER$EVAL1(VAL,EXC,DEXC)                               **
!**    PERDEWZUNGER$EVAL2(VAL,EXC,DEXC,D2EXC)                         **
!**    PERDEWZUNGER$EVAL2(VAL,EXC,DEXC,D2EXC,D3EXC)                   **
!**                                                                   **
!**  DEPENDENCIES:                                                    **
!**    DFT_MODULE(ONLY:TSPIN)                                         **
!**                                                                   **
!**  REMARKS : NO TABLE LOOKUP YET                                    **
!**                                                                   **
!**                                                                   **
!***********************************************************************
USE DFT_MODULE,ONLY : TSPIN
LOGICAL(4)        :: TINI=.FALSE.
REAL(8),PARAMETER :: GAMMAU=-142.3D-3
REAL(8),PARAMETER :: BETA1U=1052.9D-3
REAL(8),PARAMETER :: BETA2U= 333.4D-3
REAL(8),PARAMETER :: GAMMAP= -84.3D-3
REAL(8),PARAMETER :: BETA1P=1398.1D-3
REAL(8),PARAMETER :: BETA2P= 261.1D-3
REAL(8),PARAMETER :: AU=31.10D-3
REAL(8),PARAMETER :: BU=-48.0D-3
REAL(8),PARAMETER :: CU=  2.0D-3
REAL(8),PARAMETER :: DU=-11.6D-3
REAL(8),PARAMETER :: AP= 15.55D-3
REAL(8),PARAMETER :: BP=-26.9D-3
REAL(8),PARAMETER :: CP=  0.7D-3
REAL(8),PARAMETER :: DP= -4.8D-3
!     ==  DENSITY GRID =================================================
!     == THIS INTERPOATION WITH DRHO 5.E-5 IS ACCURATE FOR RHO>5.D-5
!     == SPIN GRID =====================================================
!== CONSTANTS =====================================================
REAL(8)              :: RSFAC
REAL(8)              :: FACF,FACDF,FACD2F,FACD3F
CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERDEWZUNGER_INITIALIZE
!     ******************************************************************
!     ** INITIALIZE ARRAYS
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)    :: PI
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0)
      FACF=1.D0/(2.D0**(4.D0/3.D0)-2.D0)
      FACDF=FACF*4.D0/3.D0
      FACD2F=FACDF/3.D0
      FACD3F=-FACD2F*2.D0/3.D0
      RETURN
      END SUBROUTINE PERDEWZUNGER_INITIALIZE
END MODULE PERDEWZUNGER_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERDEWZUNGER_RHODEP(RHOT,ECU,VCU,ECP,VCP)
!     **************************************************************
!     **                                                          **
!     ** CALCULATE EC =RHOT*EPSILON_C                             **
!     **       AND VC =(RHOT*EPSILON_C)/DRHOT                     **
!     **                                                          **
!     ** FOR THE FULLY NON-SPINPOLARIZED ELECTRON GAS (ECU,VCU)   **
!     ** FOR THE FULLY     SPINPOLARIZED ELECTRON GAS (ECP,VCP    **
!     ** (ONLY CORRELATION!)                                      **
!     **                                                          **
!     ** J.P. PERDEW AND A. ZUNGER, PRB23, 5048 (1981)            **
!     ** USE EQ.C3 FOR EPSILON_C AND RS >1                        **
!     ** USE EQ.C4 FOR V_C       AND RS >1                        **
!     ** USE EQ.C5 FOR EPXILON_C AND RS<1                         **
!     ** USE EQ.C6 FOR V_C       AND RS<1                         **
!     **                                                          **
!     ** REMARK: REQUIRES THAT RSFAC IS SET IN DFT_MODULE!!       **
!     **************************************************************
      USE PERDEWZUNGER_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)    :: RHOT   ! TOTAL DENSITY
      REAL(8),INTENT(OUT)   :: ECU    ! CORRELATION ENERGY (RHOS=0)
      REAL(8),INTENT(OUT)   :: VCU    ! POTENTIAL (RHOS=0)
      REAL(8),INTENT(OUT)   :: ECP    ! CORRELATION ENERGY (RHOS=RHOT)
      REAL(8),INTENT(OUT)   :: VCP    ! POTENTIAL (RHOS=RHOT)
      REAL(8),PARAMETER     :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER     :: ONEBY9=1.D0/9.D0
      REAL(8),PARAMETER     :: ONEBY27=1.D0/27.D0
      REAL(8),PARAMETER     :: TWOBY3=2.D0/3.D0
      REAL(8),PARAMETER     :: FORBY3=4.D0/3.D0
      REAL(8),PARAMETER     :: SEVBY6=7.D0/6.D0
      REAL(8)               :: RS
      REAL(8)               :: RSLN
      REAL(8)               :: RSLNRS
      REAL(8)               :: SQRS
      REAL(8)               :: SVAR,SVAR1
!     **************************************************************
      IF(RHOT.EQ.0.D0) THEN
        ECU=0.D0
        VCU=0.D0
        ECP=0.D0
        VCP=0.D0
        RETURN
      END IF
      RS=RSFAC/RHOT**ONEBY3 
!     ==================================================================
!     == THE HIGH DENSITY LIMIT  (RS<1)                               ==
!     ==================================================================
      IF(RS.LT.1.D0) THEN
        RSLN  =LOG(RS)
        RSLNRS=RS*RSLN
        ECU = RHOT*(BU + DU*RS + (AU+CU*RS)*RSLN)
        VCU = ONEBY3 *(-AU+3.D0*BU - (CU-2.D0*DU)*RS &
     &                 + 3.D0*AU*RSLN + 2.D0*CU*RSLNRS)
!
        IF(TSPIN) THEN
          ECP = RHOT*(BP + DP*RS + (AP+CP*RS)*RSLN)
          VCP = ONEBY3 *(-AP+3.D0*BP - (CP-2.D0*DP)*RS &
     &                 + 3.D0*AP*RSLN + 2.D0*CP*RSLNRS)
        ELSE
          ECP=0.D0
          VCP=0.D0
        END IF
!
!     ==================================================================
!     == THE LOW DENSITY LIMIT (RS>1)                                 ==
!     ==================================================================
      ELSE
        SQRS=SQRT(RS)
        SVAR=1.D0/(1.D0+SQRS*(BETA1U+BETA2U*SQRS))
        SVAR1=SVAR/(6.D0*RHOT)
        ECU=GAMMAU*RHOT*SVAR
        VCU=ECU*SVAR1*(6.D0+SQRS*(7.D0*BETA1U+8.D0*BETA2U*SQRS)) 
!
        IF(TSPIN) THEN
          SVAR=1.D0/(1.D0+SQRS*(BETA1P+BETA2P*SQRS))
          SVAR1=SVAR/(6.D0*RHOT)
          ECP =GAMMAP*RHOT*SVAR
          VCP =ECP*SVAR1*(6.D0+SQRS*(7.D0*BETA1P+8.D0*BETA2P*SQRS)) 
        ELSE
          ECP=0.D0
          VCP=0.D0
        END IF
      END IF                     
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERDEWZUNGER_RHODEP2(RHOT,ECU,VCU,V2CU &
     &                                    ,ECP,VCP,V2CP)
!     **************************************************************
!     **                                                          **
!     ** CALCULATE EC =RHOT*EPSILON_C                             **
!     **       AND VC =(RHOT*EPSILON_C)/DRHOT                     **
!     **                                                          **
!     ** FOR THE FULLY NON-SPINPOLARIZED ELECTRON GAS (ECU,VCU)   **
!     ** FOR THE FULLY     SPINPOLARIZED ELECTRON GAS (ECP,VCP    **
!     ** (ONLY CORRELATION!)                                      **
!     **                                                          **
!     ** J.P. PERDEW AND A. ZUNGER, PRB23, 5048 (1981)            **
!     ** USE EQ.C3 FOR EPSILON_C AND RS >1                        **
!     ** USE EQ.C4 FOR V_C       AND RS >1                        **
!     ** USE EQ.C5 FOR EPXILON_C AND RS<1                         **
!     ** USE EQ.C6 FOR V_C       AND RS<1                         **
!     **                                                          **
!     ** REMARK: REQUIRES THAT RSFAC IS SET IN DFT_MODULE!!       **
!     **************************************************************
      USE PERDEWZUNGER_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)    :: RHOT   ! TOTAL DENSITY
      REAL(8),INTENT(OUT)   :: ECU    ! CORRELATION ENERGY (RHOS=0)
      REAL(8),INTENT(OUT)   :: VCU    ! POTENTIAL (RHOS=0)
      REAL(8),INTENT(OUT)   :: V2CU   ! D2[ECU]/D[RHOT]**2
      REAL(8),INTENT(OUT)   :: ECP    ! CORRELATION ENERGY (RHOS=RHOT)
      REAL(8),INTENT(OUT)   :: VCP    ! POTENTIAL (RHOS=RHOT)
      REAL(8),INTENT(OUT)   :: V2CP   ! D2[ECP]/D[RHOT]**2
      REAL(8),PARAMETER     :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER     :: ONEBY9=1.D0/9.D0
      REAL(8),PARAMETER     :: ONEBY27=1.D0/27.D0
      REAL(8),PARAMETER     :: TWOBY3=2.D0/3.D0
      REAL(8),PARAMETER     :: FORBY3=4.D0/3.D0
      REAL(8),PARAMETER     :: SEVBY6=7.D0/6.D0
      REAL(8)               :: RS
      REAL(8)               :: RSLN
      REAL(8)               :: RSLNRS
      REAL(8)               :: SQRS
      REAL(8)               :: SVAR,SVAR1
!     **************************************************************
      IF(RHOT.EQ.0.D0) THEN
        ECU=0.D0
        VCU=0.D0
        V2CU=0.D0
        ECP=0.D0
        VCP=0.D0
        V2CP=0.D0
        RETURN
      END IF
      RS=RSFAC/RHOT**ONEBY3 
!     ==================================================================
!     == THE HIGH DENSITY LIMIT  (RS<1)                               ==
!     ==================================================================
      IF(RS.LT.1.D0) THEN
        RSLN  =LOG(RS)
        RSLNRS=RS*RSLN
        ECU = RHOT*(BU + DU*RS + (AU+CU*RS)*RSLN)
        VCU = ONEBY3 *(-AU+3.D0*BU - (CU-2.D0*DU)*RS &
     &                 + 3.D0*AU*RSLN + 2.D0*CU*RSLNRS)
        V2CU=-ONEBY9/RHOT*(3.D0*AU+(CU+2.D0*DU)*RS+2.D0*CU*RSLNRS)
!
        IF(TSPIN) THEN
          ECP = RHOT*(BP + DP*RS + (AP+CP*RS)*RSLN)
          VCP = ONEBY3 *(-AP+3.D0*BP - (CP-2.D0*DP)*RS &
     &                 + 3.D0*AP*RSLN + 2.D0*CP*RSLNRS)
          V2CP=-ONEBY9/RHOT*( 3.D0*AP+(CP+2.D0*DP)*RS+2.D0*CP*RSLNRS)
        ELSE
          ECP=0.D0
          VCP=0.D0
          V2CP=0.D0
        END IF
!
!     ==================================================================
!     == THE LOW DENSITY LIMIT (RS>1)                                 ==
!     ==================================================================
      ELSE
        SQRS=SQRT(RS)
        SVAR=1.D0/(1.D0+SQRS*(BETA1U+BETA2U*SQRS))
        SVAR1=SVAR/(6.D0*RHOT)
        ECU=GAMMAU*RHOT*SVAR
        VCU=ECU*SVAR1*(6.D0+SQRS*(7.D0*BETA1U+8.D0*BETA2U*SQRS)) 
        V2CU=ECU*SVAR1**2*SQRS &
       &     *( 5.D0*BETA1U &
       &        + SQRS*((7.D0*BETA1U**2+8.D0*BETA2U) &
       &          + SQRS*BETA2U*(21.D0*BETA1U &
       &            + SQRS*BETA2U*16.D0)))
!
        IF(TSPIN) THEN
          SVAR=1.D0/(1.D0+SQRS*(BETA1P+BETA2P*SQRS))
          SVAR1=SVAR/(6.D0*RHOT)
          ECP =GAMMAP*RHOT*SVAR
          VCP =ECP*SVAR1*(6.D0+SQRS*(7.D0*BETA1P+8.D0*BETA2P*SQRS)) 
          V2CP=ECP*SVAR1**2*SQRS &
       &       *( 5.D0*BETA1P &
       &          + SQRS*((7.D0*BETA1P**2+8.D0*BETA2P) &
       &            + SQRS*BETA2P*(21.D0*BETA1P &
       &              + SQRS*BETA2P*16.D0)))
        ELSE
          ECP=0.D0
          VCP=0.D0
          V2CP=0.D0
        END IF
      END IF                     
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERDEWZUNGER_RHODEP3(RHOT,ECU,VCU,V2CU,V3CU &
     &                                    ,ECP,VCP,V2CP,V3CP)
!     **************************************************************
!     **                                                          **
!     ** CALCULATE EC =RHOT*EPSILON_C                             **
!     **       AND VC =(RHOT*EPSILON_C)/DRHOT                     **
!     **                                                          **
!     ** FOR THE FULLY NON-SPINPOLARIZED ELECTRON GAS (ECU,VCU)   **
!     ** FOR THE FULLY     SPINPOLARIZED ELECTRON GAS (ECP,VCP    **
!     ** (ONLY CORRELATION!)                                      **
!     **                                                          **
!     ** J.P. PERDEW AND A. ZUNGER, PRB23, 5048 (1981)            **
!     ** USE EQ.C3 FOR EPSILON_C AND RS >1                        **
!     ** USE EQ.C4 FOR V_C       AND RS >1                        **
!     ** USE EQ.C5 FOR EPXILON_C AND RS<1                         **
!     ** USE EQ.C6 FOR V_C       AND RS<1                         **
!     **                                                          **
!     ** REMARK: REQUIRES THAT RSFAC IS SET IN DFT_MODULE!!       **
!     **************************************************************
      USE PERDEWZUNGER_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)    :: RHOT   ! TOTAL DENSITY
      REAL(8),INTENT(OUT)   :: ECU    ! CORRELATION ENERGY (RHOS=0)
      REAL(8),INTENT(OUT)   :: VCU    ! POTENTIAL (RHOS=0)
      REAL(8),INTENT(OUT)   :: V2CU   ! D2[ECU]/D[RHOT]**2
      REAL(8),INTENT(OUT)   :: V3CU   ! D3[ECU]/D[RHOT]**3
      REAL(8),INTENT(OUT)   :: ECP    ! CORRELATION ENERGY (RHOS=RHOT)
      REAL(8),INTENT(OUT)   :: VCP    ! POTENTIAL (RHOS=RHOT)
      REAL(8),INTENT(OUT)   :: V2CP   ! D2[ECP]/D[RHOT]**2
      REAL(8),INTENT(OUT)   :: V3CP   ! D3[ECP]/D[RHOT]**3
      REAL(8),PARAMETER     :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER     :: ONEBY9=1.D0/9.D0
      REAL(8),PARAMETER     :: ONEBY27=1.D0/27.D0
      REAL(8),PARAMETER     :: TWOBY3=2.D0/3.D0
      REAL(8),PARAMETER     :: FORBY3=4.D0/3.D0
      REAL(8),PARAMETER     :: SEVBY6=7.D0/6.D0
      REAL(8)               :: RS
      REAL(8)               :: RSLN
      REAL(8)               :: RSLNRS
      REAL(8)               :: SQRS
      REAL(8)               :: SVAR,SVAR1
!     **************************************************************
      IF(RHOT.EQ.0.D0) THEN
        ECU=0.D0
        VCU=0.D0
        V2CU=0.D0
        V3CU=0.D0
        ECP=0.D0
        VCP=0.D0
        V2CP=0.D0
        V3CP=0.D0
        RETURN
      END IF
      RS=RSFAC/RHOT**ONEBY3 
!     ==================================================================
!     == THE HIGH DENSITY LIMIT  (RS<1)                               ==
!     ==================================================================
      IF(RS.LT.1.D0) THEN
        RSLN  =LOG(RS)
        RSLNRS=RS*RSLN
        ECU = RHOT*(BU + DU*RS + (AU+CU*RS)*RSLN)
        VCU = ONEBY3 *(-AU+3.D0*BU - (CU-2.D0*DU)*RS &
     &                 + 3.D0*AU*RSLN + 2.D0*CU*RSLNRS)
        V2CU=-ONEBY9/RHOT*(3.D0*AU+(CU+2.D0*DU)*RS+2.D0*CU*RSLNRS)
        V3CU= ONEBY27/RHOT**2*(9.D0*AU+(6.D0*CU+8.D0*DU)*RS+8.D0*CU*RSLNRS)
!
        IF(TSPIN) THEN
          ECP = RHOT*(BP + DP*RS + (AP+CP*RS)*RSLN)
          VCP = ONEBY3 *(-AP+3.D0*BP - (CP-2.D0*DP)*RS &
     &                 + 3.D0*AP*RSLN + 2.D0*CP*RSLNRS)
          V2CP=-ONEBY9/RHOT*( 3.D0*AP+(CP+2.D0*DP)*RS+2.D0*CP*RSLNRS)
          V3CP=ONEBY27/RHOT**2*(9.D0*AP+(6.D0*CP+8.D0*DP)*RS+8.D0*CP*RSLNRS)
        ELSE
          ECP=0.D0
          VCP=0.D0
          V2CP=0.D0
          V3CP=0.D0
        END IF
!
!     ==================================================================
!     == THE LOW DENSITY LIMIT (RS>1)                                 ==
!     ==================================================================
      ELSE
        SQRS=SQRT(RS)
        SVAR=1.D0/(1.D0+SQRS*(BETA1U+BETA2U*SQRS))
        SVAR1=SVAR/(6.D0*RHOT)
        ECU=GAMMAU*RHOT*SVAR
        VCU=ECU*SVAR1*(6.D0+SQRS*(7.D0*BETA1U+8.D0*BETA2U*SQRS)) 
        V2CU=ECU*SVAR1**2*SQRS &
       &     *( 5.D0*BETA1U &
       &        + SQRS*((7.D0*BETA1U**2+8.D0*BETA2U) &
       &          + SQRS*BETA2U*(21.D0*BETA1U &
       &            + SQRS*BETA2U*16.D0)))
        V3CU=-ECU*SVAR1**3*SQRS &
       &     *( 35.D0*BETA1U &
       &       + SQRS*((76.D0*BETA1U**2+64.D0*BETA2U) &
       &         + SQRS*(BETA1U*(35.D0*BETA1U**2+234.D0*BETA2U) &
       &           + SQRS*BETA2U*((140.D0*BETA1U**2+176.D0*BETA2U) &
       &             + SQRS*BETA2U*(175.D0*BETA1U &
       &               + SQRS*BETA2U*64.D0)))))
!
        IF(TSPIN) THEN
          SVAR=1.D0/(1.D0+SQRS*(BETA1P+BETA2P*SQRS))
          SVAR1=SVAR/(6.D0*RHOT)
          ECP =GAMMAP*RHOT*SVAR
          VCP =ECP*SVAR1*(6.D0+SQRS*(7.D0*BETA1P+8.D0*BETA2P*SQRS)) 
          V2CP=ECP*SVAR1**2*SQRS &
       &       *( 5.D0*BETA1P &
       &          + SQRS*((7.D0*BETA1P**2+8.D0*BETA2P) &
       &            + SQRS*BETA2P*(21.D0*BETA1P &
       &              + SQRS*BETA2P*16.D0)))
          V3CP=-ECP*SVAR1**3*SQRS &
       &       *( 35.D0*BETA1P &
       &         + SQRS*((76.D0*BETA1P**2+64.D0*BETA2P) &
       &           + SQRS*(BETA1P*(35.D0*BETA1P**2+234.D0*BETA2P) &
       &             + SQRS*BETA2P*((140.D0*BETA1P**2+176.D0*BETA2P) &
       &               + SQRS*BETA2P*(175.D0*BETA1P &
       &                 + SQRS*BETA2P*64.D0)))))
        ELSE
          ECP=0.D0
          VCP=0.D0
          V2CP=0.D0
          V3CP=0.D0
        END IF
      END IF                     
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PERDEWZUNGER$EVAL1(VAL,EXC,DEXC)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE PERDEWZUNGER_MODULE
      IMPLICIT NONE
      REAL(8) ,INTENT(IN)  :: VAL(5)
      REAL(8) ,INTENT(OUT) :: EXC
      REAL(8) ,INTENT(OUT) :: DEXC(5)
      REAL(8) ,PARAMETER   :: ONEBY3=1.D0/3.D0
      REAL(8)              :: RHOT,RHOS
      REAL(8)              :: ECU,ECP
      REAL(8)              :: DECU,DECP
      REAL(8)              :: SIGMA
      REAL(8)              :: SIGMAP13,SIGMAP23,SIGMAP43
      REAL(8)              :: SIGMAM13,SIGMAM23,SIGMAM43
      REAL(8)              :: F,DF
      REAL(8)              :: F_R,F_S
      REAL(8)              :: ECS,DECS
!     ******************************************************************
      IF(.NOT.TINI) CALL PERDEWZUNGER_INITIALIZE
!
!     ==================================================================
!     ==  RESET  OUTPUT VALUES                                        ==
!     ==================================================================
      RHOT=VAL(1)
      RHOS=VAL(2)
      DEXC(:)=0.D0
!
!     ==================================================================
!     ==  CALCULATE RHO DEPENDENCES                                   ==
!     ==================================================================
      CALL PERDEWZUNGER_RHODEP(RHOT,ECU,DECU,ECP,DECP)
!
!     ==================================================================
!     ==  SHORTCUT FOR NON-SPIN POLARIZED CALCULATIONS                ==
!     ==================================================================
      IF(.NOT.TSPIN) THEN
        EXC         =ECU
        DEXC(1)     =DECU
        RETURN
      END IF
!
!     ==================================================================
!     ==  CALCULATE SPIN DEPENDENCE=====================================
!     ==================================================================
      SIGMA=RHOS/RHOT
      SIGMAP13=(1.D0+SIGMA)**ONEBY3
      SIGMAP23=SIGMAP13*SIGMAP13
      SIGMAP43=SIGMAP23*SIGMAP23
      SIGMAM13=(1.D0-SIGMA)**ONEBY3
      SIGMAM23=SIGMAM13*SIGMAM13
      SIGMAM43=SIGMAM23*SIGMAM23
      F  =FACF  *(SIGMAP43+SIGMAM43-2.D0)
      DF =FACDF *(SIGMAP13-SIGMAM13)
      F_R=-DF*SIGMA/RHOT
      F_S=DF/RHOT
!     
      ECS=ECP-ECU
      DECS=DECP-DECU
      EXC=ECU+F*ECS
      DEXC(1)=DECU+F*DECS+F_R*ECS
      DEXC(2)=F_S*ECS
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PERDEWZUNGER$EVAL2(VAL,EXC,DEXC,D2EXC)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE PERDEWZUNGER_MODULE
      IMPLICIT NONE
      REAL(8) ,INTENT(IN)  :: VAL(5)
      REAL(8) ,INTENT(OUT) :: EXC
      REAL(8) ,INTENT(OUT) :: DEXC(5)
      REAL(8) ,INTENT(OUT) :: D2EXC(5,5)
      REAL(8) ,PARAMETER   :: ONEBY3=1.D0/3.D0
      REAL(8)              :: RHOT,RHOS
      REAL(8)              :: ECU,ECP
      REAL(8)              :: DECU,DECP
      REAL(8)              :: D2ECU,D2ECP
      REAL(8)              :: SIGMA
      REAL(8)              :: SIGMAP13,SIGMAP23,SIGMAP43
      REAL(8)              :: SIGMAM13,SIGMAM23,SIGMAM43
      REAL(8)              :: F,DF,D2F
      REAL(8)              :: SIG_R,SIG_S,SIG_RR,SIG_RS
      REAL(8)              :: F_R,F_S,F_RR,F_RS,F_SS
      REAL(8)              :: ECS,DECS,D2ECS
!     ******************************************************************
      IF(.NOT.TINI) CALL PERDEWZUNGER_INITIALIZE
!
!     ==================================================================
!     ==  RESET  OUTPUT VALUES                                        ==
!     ==================================================================
      RHOT=VAL(1)
      RHOS=VAL(2)
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
!
!     ==================================================================
!     ==  CALCULATE RHO DEPENDENCES                                   ==
!     ==================================================================
      CALL PERDEWZUNGER_RHODEP2(RHOT,ECU,DECU,D2ECU,ECP,DECP,D2ECP)
!
!     ==================================================================
!     ==  SHORTCUT FOR NON-SPIN POLARIZED CALCULATIONS                ==
!     ==================================================================
      IF(.NOT.TSPIN) THEN
        EXC         =ECU
        DEXC(1)     =DECU
        D2EXC(1,1)  =D2ECU
        RETURN
      END IF
!
!     ==================================================================
!     ==  CALCULATE SPIN DEPENDENCE=====================================
!     ==================================================================
      SIGMA=RHOS/RHOT
      SIGMAP13=(1.D0+SIGMA)**ONEBY3
      SIGMAP23=SIGMAP13*SIGMAP13
      SIGMAP43=SIGMAP23*SIGMAP23
      SIGMAM13=(1.D0-SIGMA)**ONEBY3
      SIGMAM23=SIGMAM13*SIGMAM13
      SIGMAM43=SIGMAM23*SIGMAM23
      F  =FACF  *(SIGMAP43+SIGMAM43-2.D0)
      DF =FACDF *(SIGMAP13-SIGMAM13)
      D2F=FACD2F*(1.D0/SIGMAP23+1.D0/SIGMAM23)
      SIG_R=-SIGMA/RHOT
      SIG_S=1.D0/RHOT
      SIG_RR=-2.D0*SIG_R/RHOT
      SIG_RS=-SIG_S/RHOT
      F_R=DF*SIG_R
      F_S=DF*SIG_S
      F_RR=D2F*SIG_R**2+DF*SIG_RR
      F_RS=D2F*SIG_R*SIG_S+DF*SIG_RS
      F_SS=D2F*SIG_S**2
!     
      ECS=ECP-ECU
      DECS=DECP-DECU
      D2ECS=D2ECP-D2ECU
      EXC=ECU+F*ECS
      DEXC(1)=DECU+F*DECS+F_R*ECS
      DEXC(2)=F_S*ECS
      D2EXC(1,1)=D2ECU + F*D2ECS +2.D0*F_R*DECS + F_RR*ECS 
      D2EXC(2,1)=F_RS*ECS+F_S*DECS
      D2EXC(1,2)=D2EXC(2,1)
      D2EXC(2,2)=F_SS*ECS
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PERDEWZUNGER$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE PERDEWZUNGER_MODULE
      IMPLICIT NONE
      REAL(8) ,INTENT(IN)  :: VAL(5)
      REAL(8) ,INTENT(OUT) :: EXC
      REAL(8) ,INTENT(OUT) :: DEXC(5)
      REAL(8) ,INTENT(OUT) :: D2EXC(5,5)
      REAL(8) ,INTENT(OUT) :: D3EXC(5,5,5)
      REAL(8) ,PARAMETER   :: ONEBY3=1.D0/3.D0
      REAL(8)              :: RHOT,RHOS
      REAL(8)              :: ECU,ECP
      REAL(8)              :: DECU,DECP
      REAL(8)              :: D2ECU,D2ECP
      REAL(8)              :: D3ECU,D3ECP
      REAL(8)              :: SIGMA
      REAL(8)              :: SIGMAP13,SIGMAP23,SIGMAP43,SIGMAP53
      REAL(8)              :: SIGMAM13,SIGMAM23,SIGMAM43,SIGMAM53
      REAL(8)              :: F,DF,D2F,D3F
      REAL(8)              :: SIG_R,SIG_S,SIG_RR,SIG_RS
      REAL(8)              :: SIG_RRR,SIG_RRS
      REAL(8)              :: F_R,F_S,F_RR,F_RS,F_SS
      REAL(8)              :: F_RRR,F_RRS,F_RSS,F_SSS
      REAL(8)              :: ECS,DECS,D2ECS,D3ECS
!     ******************************************************************
      IF(.NOT.TINI) CALL PERDEWZUNGER_INITIALIZE
!
!     ==================================================================
!     ==  RESET  OUTPUT VALUES                                        ==
!     ==================================================================
      RHOT=VAL(1)
      RHOS=VAL(2)
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      D3EXC(:,:,:)=0.D0
!
!     ==================================================================
!     ==  CALCULATE RHO DEPENDENCES                                   ==
!     ==================================================================
      CALL PERDEWZUNGER_RHODEP3(RHOT,ECU,DECU,D2ECU,D3ECU,ECP,DECP,D2ECP,D3ECP)
!
!     ==================================================================
!     ==  SHORTCUT FOR NON-SPIN POLARIZED CALCULATIONS                ==
!     ==================================================================
      IF(.NOT.TSPIN) THEN
        EXC         =ECU
        DEXC(1)     =DECU
        D2EXC(1,1)  =D2ECU
        D3EXC(1,1,1)=D3ECU
        RETURN
      END IF
!
!     ==================================================================
!     ==  CALCULATE SPIN DEPENDENCE=====================================
!     ==================================================================
      SIGMA=RHOS/RHOT
      SIGMAP13=(1.D0+SIGMA)**ONEBY3
      SIGMAP23=SIGMAP13*SIGMAP13
      SIGMAP43=SIGMAP23*SIGMAP23
      SIGMAP53=SIGMAP13*SIGMAP43
      SIGMAM13=(1.D0-SIGMA)**ONEBY3
      SIGMAM23=SIGMAM13*SIGMAM13
      SIGMAM43=SIGMAM23*SIGMAM23
      SIGMAM53=SIGMAM13*SIGMAM43
      F  =FACF  *(SIGMAP43+SIGMAM43-2.D0)
      DF =FACDF *(SIGMAP13-SIGMAM13)
      D2F=FACD2F*(1.D0/SIGMAP23+1.D0/SIGMAM23)
      D3F=FACD3F*(1.D0/SIGMAP53-1.D0/SIGMAM53)
      SIG_R=-SIGMA/RHOT
      SIG_S=1.D0/RHOT
      SIG_RR=-2.D0*SIG_R/RHOT
      SIG_RS=-SIG_S/RHOT
      SIG_RRR=-3.D0*SIG_RR/RHOT
      SIG_RRS=-2.D0*SIG_RS/RHOT
      F_R=DF*SIG_R
      F_S=DF*SIG_S
      F_RR=D2F*SIG_R**2+DF*SIG_RR
      F_RS=D2F*SIG_R*SIG_S+DF*SIG_RS
      F_SS=D2F*SIG_S**2
      F_RRR=D3F*SIG_R**3+3.D0*D2F*SIG_R*SIG_RR+DF*SIG_RRR
      F_RRS=D3F*SIG_R**2*SIG_S+2.D0*D2F*SIG_R*SIG_RS+D2F*SIG_RR*SIG_S+DF*SIG_RRS
      F_RSS=D3F*SIG_R*SIG_S**2+2.D0*D2F*SIG_S*SIG_RS
      F_SSS=D3F*SIG_S**3
!     
      ECS=ECP-ECU
      DECS=DECP-DECU
      D2ECS=D2ECP-D2ECU
      D3ECS=D3ECP-D3ECU
      EXC         =ECU+F*ECS
      DEXC(1)     =DECU+F*DECS+F_R*ECS
      DEXC(2)     =F_S*ECS
      D2EXC(1,1)  =D2ECU + F*D2ECS +2.D0*F_R*DECS + F_RR*ECS 
      D2EXC(1,2)  =F_RS*ECS+F_S*DECS
      D2EXC(2,2)  =F_SS*ECS
      D3EXC(1,1,1)=D3ECU+F*D3ECS+3.D0*F_R*D2ECS+3.D0*F_RR*DECS+F_RRR*ECS
      D3EXC(1,1,2)=F_RRS*ECS+2.D0*F_RS*DECS+F_S*D2ECS
      D3EXC(1,2,2)=F_RSS*ECS+F_SS*DECS
      D3EXC(2,2,2)=F_SSS*ECS
      D2EXC(2,1)  =D2EXC(1,2)
      D3EXC(1,2,1)=D3EXC(1,1,2)
      D3EXC(2,1,1)=D3EXC(1,1,2)
      D3EXC(2,1,2)=D3EXC(1,2,2)
      D3EXC(2,2,1)=D3EXC(1,2,2)
      RETURN
      END
!
!     ..................................................................
      MODULE BARTHHEDIN_MODULE
!     ==  INITIALIZE BARTH HEDIN FORM OF CORRELATION                ==
!     ==  U.VON BARTH; L. HEDIN, J.PHYS.C:SOL.ST.PHYS. 51629(1972)  == 
!     ==  CAN BE USED IN DIFFERENT PARAMETRIZATIONS ACCORDING TO THE==
!     ==  PARAMETERS GIVEN ABOVE
      LOGICAL(4) :: TBH
      LOGICAL(4) :: TRSK
      REAL(8)    :: CPBH
      REAL(8)    :: RPBH
      REAL(8)    :: CFBH
      REAL(8)    :: RFBH
      REAL(8)   :: RSFAC=0.D0 
      END MODULE BARTHHEDIN_MODULE 
!
!     ...................................................................
      SUBROUTINE BARTHHEDIN$INITIALIZE(ID_)
!     ********************************************************************
!     ********************************************************************
      USE DFT_MODULE
      USE BARTHHEDIN_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_    ! EITHER BH OR RSK
!     ******************************************************************
      IF(ID_.EQ.'BH') THEN
!       == PARAMETERS FOR BARTH HEDIN (J.PHYS.C5,P1629(1972)) =================    
        CPBH=0.0504D0
        RPBH=30.D0
        CFBH=0.0254D0
        RFBH=75.D0
      ELSE IF(ID_.EQ.'RSK') THEN
!       == PARAMETERS FOR RAJAGOPAL, SINGHAL,KIMBALL
        CPBH=0.04612D0
        RPBH=39.7D0
        CFBH=0.02628D0
        RFBH=70.6D0
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$STOP('BARTHHEDIN_INITIALIZE')
      END IF
!     ********************************************************************
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE BARTHHEDIN$EVAL(RHOT,ARRU,DARRU,ARRP,DARRP)
!     ********************************************************************
!     ********************************************************************
      USE BARTHHEDIN_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)    :: RHOT
      REAL(8),INTENT(INOUT) :: ARRU
      REAL(8),INTENT(INOUT) :: DARRU
      REAL(8),INTENT(INOUT) :: ARRP
      REAL(8),INTENT(INOUT) :: DARRP 
      REAL(8),PARAMETER     :: ONEBY3=1.D0/3.D0
      REAL(8)   :: RHO13
      REAL(8)   :: RS
      REAL(8)   :: ZU,ZU2,ZU3,ULN,FOFZU
      REAL(8)   :: ZP,ZP2,ZP3,PLN,FOFZP
!     ********************************************************************
      RHO13=RHOT**ONEBY3
      RS=RSFAC/RHO13
      ZU=RS/RPBH
      ZP=RS/RFBH
      ZU2=ZU*ZU
      ZP2=ZP*ZP
      ZU3=ZU2*ZU
      ZP3=ZP2*ZP
      ULN=LOG(1.D0+1.D0/ZU)
      PLN=LOG(1.D0+1.D0/ZP)
      FOFZU=(1.D0+ZU3)*ULN+0.5D0*ZU-ZU2-ONEBY3
      FOFZP=(1.D0+ZP3)*PLN+0.5D0*ZP-ZP2-ONEBY3
      ARRU=ARRU-0.5D0*CPBH*RHOT*FOFZU
      ARRP=ARRP-0.5D0*CFBH*RHOT*FOFZP
      DARRU=DARRU-0.5D0*CPBH*ULN
      DARRP=DARRP-0.5D0*CFBH*PLN
      RETURN
      END
!
!     ..................................................................
MODULE PERDEW_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: PERDEW                                                     **
!**                                                                   **
!**  PURPOSE: CALCULATE THE PERDEW86 GRADIENT CORRECTION FOR          **
!**    CORRELATION                                                    **
!**    J.P. PERDEW, PHYS.REV.B,33,8822(1986)                          **
!**                                                                   **
!**  FUNCTIONS::                                                      **
!**    PERDEW$EVAL                                                    **
!**    PERDEW$EVAL2(VAL,EXC,DEXC,D2EXC)                               **
!**    PERDEW$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)                         **
!**                                                                   **
!***********************************************************************
USE DFT_MODULE ,ONLY : TSPIN,TSAFE
USE TABLE1D_MODULE
LOGICAL(4)        :: TINI=.FALSE.
!=======================================================================
!==  CONSTANT PARAMETERS                                              ==
!=======================================================================
REAL(8),PARAMETER :: C1    = 1.667D-3
REAL(8),PARAMETER :: C2    = 2.568D-3
REAL(8),PARAMETER :: ALPHA =23.266D-3
REAL(8),PARAMETER :: BETA  = 7.389D-6
REAL(8),PARAMETER :: GAMMA = 8.723D0
REAL(8),PARAMETER :: DELTA = 0.472D0
REAL(8),PARAMETER :: CINFTY=C1+C2
REAL(8),PARAMETER :: F1745 =1.745D0*0.11D0*CINFTY
!=======================================================================
!==  CONSTANTS CALCULATED DURING INITIALIZATION                       ==
!=======================================================================
REAL(8)           :: RSFAC
!=======================================================================
!==  ARRAYS FOR FAST TABLE LOOKUP                                     ==
!=======================================================================
INTEGER(4),PARAMETER   :: NRHO=1000
TYPE (XTABLE_TYPE)     :: RHOTABLE
REAL(8)   ,ALLOCATABLE :: FPHIARRAY(:)   !(NRHO)
REAL(8)   ,ALLOCATABLE :: DFPHIARRAY(:)  !(NRHO)
REAL(8)   ,ALLOCATABLE :: D2FPHIARRAY(:) !(NRHO)
REAL(8)   ,ALLOCATABLE :: FEXPARRAY(:)   !(NRHO)
REAL(8)   ,ALLOCATABLE :: DFEXPARRAY(:)  !(NRHO)
REAL(8)   ,ALLOCATABLE :: D2FEXPARRAY(:) !(NRHO)
INTEGER(4),PARAMETER   :: NSIG=1000
TYPE (XTABLE_TYPE)     :: SIGTABLE
REAL(8)   ,ALLOCATABLE :: ONEBYDARRAY(:)   !(NSIG)
REAL(8)   ,ALLOCATABLE :: DONEBYDARRAY(:)  !(NSIG)
REAL(8)   ,ALLOCATABLE :: D2ONEBYDARRAY(:) !(NSIG)
!***********************************************************************
CONTAINS
!     ..................................................................
      SUBROUTINE PERDEW_INITIALIZE
!     ================================================================
!     ==  INITIALIZE PERDEW CORRECTION                              ==
!     ==  J.P. PERDEW: PHYS.REV.B33,P8822(1986)                     ==
!     ================================================================
      IMPLICIT NONE
      REAL(8)        :: RHOT
      REAL(8)        :: SIGMA
      INTEGER(4)     :: I
      REAL(8)        :: PI
!     ****************************************************************
      PI=4.D0*ATAN(1.D0)
      RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0)
!
!     ================================================================
!     ==  CREATE TABLE FOR RHODEPENDENCE                            ==
!     ================================================================
      ALLOCATE(FPHIARRAY(NRHO))
      ALLOCATE(DFPHIARRAY(NRHO))
      ALLOCATE(D2FPHIARRAY(NRHO))
      ALLOCATE(FEXPARRAY(NRHO))
      ALLOCATE(DFEXPARRAY(NRHO))
      ALLOCATE(D2FEXPARRAY(NRHO))
      CALL TABLE1D$MAKE(RHOTABLE,1.D-3,1.D0,NRHO)
      DO I=1,NRHO
        CALL TABLE1D$XOFI(RHOTABLE,I,RHOT)
        CALL PERDEW_RHODEP2(RHOT,FPHIARRAY(I),DFPHIARRAY(I),D2FPHIARRAY(I) &
     &                          ,FEXPARRAY(I),DFEXPARRAY(I),D2FEXPARRAY(I))
      ENDDO
!
!     ================================================================
!     ==  CREATE TABLE FOR SIGMA DEPENDENCE                         ==
!     ================================================================
      ALLOCATE(ONEBYDARRAY(NSIG))
      ALLOCATE(DONEBYDARRAY(NSIG))
      ALLOCATE(D2ONEBYDARRAY(NSIG))
      CALL TABLE1D$MAKE(SIGTABLE,0.D0,1.D0-1.D-12,NSIG)
      DO I=1,NRHO
        CALL TABLE1D$XOFI(SIGTABLE,I,SIGMA)
        CALL PERDEW_SPINDEP2(SIGMA,ONEBYDARRAY(I),DONEBYDARRAY(I),D2ONEBYDARRAY(I))
      ENDDO
      TINI=.TRUE.
      RETURN 
      END SUBROUTINE PERDEW_INITIALIZE
END MODULE PERDEW_MODULE
!
!     ..................................................................
      SUBROUTINE PERDEW_RHODEP1(RHOT,FPHI,DFPHI &
     &                              ,FEXP,DFEXP)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES TWO DENSITY DEPENDENT FUNCTIONS                  **
!     **  OF THE PERDEW 86 GRADIENT-CORRECTION FOR CORRELATION        **
!     **                                                              **
!     **  FPHI=-PHI/|GRAD(RHOT)| WITH PHI DEFINED                     **
!     **  AS IN EQ.9 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **                                                              **
!     **  FEXP= C(RHOT)/RHOT**(4/3) WITH C(N) DEFINED AS              **
!     **  AS IN EQ.6 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **                                                              **
!     **  THE GRADIENT CORRECTION TO THE CORRELATION ENERGY IS        **
!     **     E = (1/D) * FEXP * EXP(FPHI*|GRHOT|) * GRHOT**2          **
!     **  WHERE GRHOT IS THE GRADEINT OF THE TOTAL DENSITY            **
!     **                                                              **
!     **  REMARK: THIS ROUTINE IS A KLON OF PERDEW_RHODEP3            **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)   :: RHOT
      REAL(8),INTENT(OUT)  :: FPHI      !OLD ARRPE: -PHI=FPHI*GRHO
      REAL(8),INTENT(OUT)  :: DFPHI
      REAL(8),INTENT(OUT)  :: FEXP      !OLD ARRPN1: EC=FEXP*EXP(-PHI)*GRHO2/D(SIGMA)
      REAL(8),INTENT(OUT)  :: DFEXP
      REAL(8)   ,PARAMETER :: ONEBY3=1.D0/3.D0
      REAL(8)   ,PARAMETER :: FOURBY3=4.D0/3.D0
      REAL(8)   ,PARAMETER :: SEVENBY6=7.D0/6.D0
      REAL(8)   ,PARAMETER :: ONEBY6=1.D0/6.D0
      REAL(8)              :: RHO13,RHO43,RHO16,RHO76
      REAL(8)              :: RS,RS_T,RS2,RS3
      REAL(8)              :: BFAC1,BFAC1_R
      REAL(8)              :: BFAC2,BFAC2_R
      REAL(8)              :: CN,DCN
      REAL(8)              :: CN_R
      REAL(8)              :: SVAR,DSVAR
!     ******************************************************************
      RHO13=RHOT**ONEBY3
      RHO43=RHO13**4
      RS=RSFAC/RHO13
      RS_T=-RS/(3.D0*RHOT)       
      RS2=RS*RS
      RS3=RS2*RS
!
      BFAC1   = C2 + ALPHA*RS + BETA*RS2
      BFAC1_R =      ALPHA    + 2.D0*BETA*RS
!
      BFAC2   =1.D0 + GAMMA*RS + DELTA*RS2     + 1.D+4*BETA*RS3
      BFAC2_R =       GAMMA     +2.D0*DELTA*RS + 3.D+4*BETA*RS2
!
      CN      =BFAC1/BFAC2
      CN_R    =(BFAC1_R-CN*BFAC2_R)/BFAC2
      CN=CN+C1
!
      DCN =CN_R*RS_T
!
      RHO16=RHOT**ONEBY6
      RHO76=RHO16*RHOT
      SVAR=CN*RHO76
      DSVAR=DCN*RHO76+SEVENBY6*CN*RHO16
!     ==================================================================
!     == -PHI=ARRPE(RHO)*GRHO                                         ==
!     ==================================================================
      FPHI=-F1745/SVAR
      DFPHI=- FPHI*DSVAR/SVAR
!     ==================================================================
!     == EC=1/D *EXP(-PHI)*ARRPN*GRHO**2 ===============================
!     ==================================================================
      FEXP=CN/RHO43
      DFEXP=DCN/RHO43-FOURBY3*FEXP/RHOT
      RETURN
      END SUBROUTINE PERDEW_RHODEP1
!
!     ..................................................................
      SUBROUTINE PERDEW_RHODEP2(RHOT,FPHI,DFPHI,D2FPHI &
     &                              ,FEXP,DFEXP,D2FEXP)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES TWO DENSITY DEPENDENT FUNCTIONS                  **
!     **  OF THE PERDEW 86 GRADIENT-CORRECTION FOR CORRELATION        **
!     **                                                              **
!     **  FPHI=-PHI/|GRAD(RHOT)| WITH PHI DEFINED                     **
!     **  AS IN EQ.9 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **                                                              **
!     **  FEXP= C(RHOT)/RHOT**(4/3) WITH C(N) DEFINED AS              **
!     **  AS IN EQ.6 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **                                                              **
!     **  THE GRADIENT CORRECTION TO THE CORRELATION ENERGY IS        **
!     **     E = (1/D) * FEXP * EXP(FPHI*|GRHOT|) * GRHOT**2          **
!     **  WHERE GRHOT IS THE GRADEINT OF THE TOTAL DENSITY            **
!     **                                                              **
!     **  REMARK: THIS ROUTINE IS A KLON OF PERDEW_RHODEP3            **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)   :: RHOT
      REAL(8),INTENT(OUT)  :: FPHI      !OLD ARRPE: -PHI=FPHI*GRHO
      REAL(8),INTENT(OUT)  :: DFPHI
      REAL(8),INTENT(OUT)  :: D2FPHI
      REAL(8),INTENT(OUT)  :: FEXP !OLD ARRPN1: EC=FEXP*EXP(-PHI)*GRHO2/D(SIGMA)
      REAL(8),INTENT(OUT)  :: DFEXP
      REAL(8),INTENT(OUT)  :: D2FEXP
      REAL(8)   ,PARAMETER :: ONEBY3=1.D0/3.D0
      REAL(8)   ,PARAMETER :: FOURBY3=4.D0/3.D0
      REAL(8)   ,PARAMETER :: SEVENBY6=7.D0/6.D0
      REAL(8)   ,PARAMETER :: ONEBY6=1.D0/6.D0
      REAL(8)              :: RHO13,RHO43,RHO16,RHO76
      REAL(8)              :: RS,RS_T,RS_TT,RS2,RS3
      REAL(8)              :: BFAC1,BFAC1_R,BFAC1_RR
      REAL(8)              :: BFAC2,BFAC2_R,BFAC2_RR
      REAL(8)              :: CN,DCN,D2CN
      REAL(8)              :: CN_R,CN_RR
      REAL(8)              :: SVAR,DSVAR,D2SVAR
!     ******************************************************************
      RHO13=RHOT**ONEBY3
      RHO43=RHO13**4
      RS=RSFAC/RHO13
      RS_T=-RS/(3.D0*RHOT)       
      RS_TT=-4.D0*RS_T/(3.D0*RHOT)
      RS2=RS*RS
      RS3=RS2*RS
!
      BFAC1   = C2 + ALPHA*RS + BETA*RS2
      BFAC1_R =      ALPHA    + 2.D0*BETA*RS
      BFAC1_RR=                 2.D0*BETA
!
      BFAC2   =1.D0 + GAMMA*RS + DELTA*RS2     + 1.D+4*BETA*RS3
      BFAC2_R =       GAMMA     +2.D0*DELTA*RS + 3.D+4*BETA*RS2
      BFAC2_RR=                  2.D0*DELTA    + 6.D+4*BETA*RS
!
      CN      =BFAC1/BFAC2
      CN_R    =(BFAC1_R-CN*BFAC2_R)/BFAC2
      CN_RR   =(BFAC1_RR-2.D0*CN_R*BFAC2_R-CN*BFAC2_RR)/BFAC2
      CN=CN+C1
!
      DCN =CN_R*RS_T
      D2CN=CN_RR*RS_T**2 + CN_R*RS_TT      
!
      RHO16=RHOT**ONEBY6
      RHO76=RHO16*RHOT
      SVAR=CN*RHO76
      DSVAR=DCN*RHO76+SEVENBY6*CN*RHO16
      D2SVAR=D2CN*RHO76+SEVENBY6*RHO16*(2.D0*DCN+CN/(6.D0*RHOT))
!     ==================================================================
!     == -PHI=ARRPE(RHO)*GRHO                                         ==
!     ==================================================================
      FPHI=-F1745/SVAR
      DFPHI=- FPHI*DSVAR/SVAR
      D2FPHI=-(2.D0*DFPHI*DSVAR+FPHI*D2SVAR)/SVAR
!     ==================================================================
!     == EC=1/D *EXP(-PHI)*ARRPN*GRHO**2 ===============================
!     ==================================================================
      FEXP=CN/RHO43
      DFEXP=DCN/RHO43-FOURBY3*FEXP/RHOT
      D2FEXP=D2CN/RHO43-FOURBY3*(2.D0*DFEXP+ONEBY3*FEXP/RHOT)/RHOT
      RETURN
      END SUBROUTINE PERDEW_RHODEP2
!
!     ..................................................................
      SUBROUTINE PERDEW_RHODEP3(RHOT,FPHI,DFPHI,D2FPHI,D3FPHI &
     &                              ,FEXP,DFEXP,D2FEXP,D3FEXP)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES TWO DENSITY DEPENDENT FUNCTIONS                  **
!     **  OF THE PERDEW 86 GRADIENT-CORRECTION FOR CORRELATION        **
!     **                                                              **
!     **  FPHI=-PHI/|GRAD(RHOT)| WITH PHI DEFINED                     **
!     **  AS IN EQ.9 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **                                                              **
!     **  FEXP= C(RHOT)/RHOT**(4/3) WITH C(N) DEFINED AS              **
!     **  AS IN EQ.6 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **                                                              **
!     **  THE GRADIENT CORRECTION TO THE CORRELATION ENERGY IS        **
!     **     E = (1/D) * FEXP * EXP(FPHI*|GRHOT|) * GRHOT**2          **
!     **  WHERE GRHOT IS THE GRADEINT OF THE TOTAL DENSITY            **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)   :: RHOT
      REAL(8),INTENT(OUT)  :: FPHI      !OLD ARRPE: -PHI=FPHI*GRHO
      REAL(8),INTENT(OUT)  :: DFPHI
      REAL(8),INTENT(OUT)  :: D2FPHI
      REAL(8),INTENT(OUT)  :: D3FPHI
      REAL(8),INTENT(OUT)  :: FEXP      !OLD ARRPN1: EC=FEXP*EXP(-PHI)*GRHO2/D(SIGMA)
      REAL(8),INTENT(OUT)  :: DFEXP
      REAL(8),INTENT(OUT)  :: D2FEXP
      REAL(8),INTENT(OUT)  :: D3FEXP
      REAL(8)   ,PARAMETER :: ONEBY3=1.D0/3.D0
      REAL(8)   ,PARAMETER :: FOURBY3=4.D0/3.D0
      REAL(8)   ,PARAMETER :: SEVENBY6=7.D0/6.D0
      REAL(8)   ,PARAMETER :: ONEBY6=1.D0/6.D0
      REAL(8)              :: RHO13,RHO43,RHO16,RHO76
      REAL(8)              :: RS,RS_T,RS_TT,RS_TTT,RS2,RS3
      REAL(8)              :: BFAC1,BFAC1_R,BFAC1_RR
      REAL(8)              :: BFAC2,BFAC2_R,BFAC2_RR,BFAC2_RRR
      REAL(8)              :: CN,DCN,D2CN,D3CN 
      REAL(8)              :: CN_R,CN_RR,CN_RRR
      REAL(8)              :: SVAR,DSVAR,D2SVAR,D3SVAR
!     ******************************************************************
      RHO13=RHOT**ONEBY3
      RHO43=RHO13**4
      RS=RSFAC/RHO13
      RS_T=-RS/(3.D0*RHOT)       
      RS_TT=-4.D0*RS_T/(3.D0*RHOT)
      RS_TTT=-7.D0*RS_TT/(3.D0*RHOT)
      RS2=RS*RS
      RS3=RS2*RS
!
!     ==================================================================
!     == CALCULATE THE FUNCTION C(RS) OF EQ.6                         ==
!     ==================================================================
      BFAC1   = C2 + ALPHA*RS + BETA*RS2
      BFAC1_R =      ALPHA    + 2.D0*BETA*RS
      BFAC1_RR=                 2.D0*BETA
!
      BFAC2   =1.D0 + GAMMA*RS + DELTA*RS2     + 1.D+4*BETA*RS3
      BFAC2_R =       GAMMA     +2.D0*DELTA*RS + 3.D+4*BETA*RS2
      BFAC2_RR=                  2.D0*DELTA    + 6.D+4*BETA*RS
      BFAC2_RRR=                                 6.D+4*BETA
!
      CN      =BFAC1/BFAC2
      CN_R    =(BFAC1_R-CN*BFAC2_R)/BFAC2
      CN_RR   =(BFAC1_RR-2.D0*CN_R*BFAC2_R-CN*BFAC2_RR)/BFAC2
      CN_RRR  =(-3.D0*CN_RR*BFAC2_R-3.D0*CN_R*BFAC2_RR-CN*BFAC2_RRR)/BFAC2
      CN=CN+C1
!
      DCN =CN_R*RS_T
      D2CN=CN_RR*RS_T**2 + CN_R*RS_TT      
      D3CN=CN_RRR*RS_T**3 + 3.D0*CN_RR*RS_T*RS_TT + CN_R*RS_TTT
!
!     ==================================================================
!     == CALCULATE FPHI=-PHI/|GRAD(RHOT)| OF EQ.9                     ==
!     ==================================================================
      RHO16=RHOT**ONEBY6
      RHO76=RHO16*RHOT
      SVAR=CN*RHO76
      DSVAR=DCN*RHO76+SEVENBY6*CN*RHO16
      D2SVAR=D2CN*RHO76+SEVENBY6*RHO16*(2.D0*DCN+CN/(6.D0*RHOT))
      D3SVAR=D3CN*RHO76+SEVENBY6*RHO16*(3.D0*D2CN+3.D0*DCN/(6.D0*RHOT) &
     &                                 -5.D0*CN/(6.D0*RHOT)**2)
!
      FPHI=-F1745/SVAR
      DFPHI=- FPHI*DSVAR/SVAR
      D2FPHI=-(2.D0*DFPHI*DSVAR+FPHI*D2SVAR)/SVAR
      D3FPHI=-(3.D0*D2FPHI*DSVAR+3.D0*DFPHI*D2SVAR+FPHI*D3SVAR)/SVAR
!
!     ==================================================================
!     ==  CALCULATE FEXP=CN/RHO**(4/3)                                ==
!     ==================================================================
      FEXP=CN/RHO43
      DFEXP=DCN/RHO43-FOURBY3*FEXP/RHOT
      D2FEXP=D2CN/RHO43-FOURBY3*(2.D0*DFEXP+ONEBY3*FEXP/RHOT)/RHOT
      D3FEXP=D3CN/RHO43-FOURBY3 &
     &     *(3.D0*D2CN/RHO43+6.D0*D2FEXP-(5.D0*DFEXP+2.D0*FEXP/RHOT)/RHOT) &
     &     /(3.D0*RHOT)
      RETURN
      END SUBROUTINE PERDEW_RHODEP3
!
!     ..................................................................
      SUBROUTINE PERDEW_SPINDEP1(SIGMA,ONEBYD,DONEBYD)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE SPIN DEPENDENCE OF THE PERDEW 86 GRADIENT-   **
!     **  CORRECTION FOR CORRELATION                                  **
!     **                                                              **
!     **  ONEBYD IS (1/D) WITH D DEFINED AS                           **
!     **   1/D=1/SQRT[{(1+SIGMA)**(5/3) + (1-SIGMA)**(5/3)}/2]        **
!     **  AS IN EQ.4 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **  WHERE SIGMA IS RHOT/RHOS                                    **
!     **                                                              **
!     **  REMARK: THIS ROUTINE IS A KLON OF PERDEW_SPINDEP3           **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: SIGMA
      REAL(8),INTENT(OUT):: ONEBYD
      REAL(8),INTENT(OUT):: DONEBYD
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER  :: FIVEBY6=5.D0/6.D0
      REAL(8)            :: SIGP,SIGM
      REAL(8)            :: SIGP13,SIGM13
      REAL(8)            :: SIGP23,SIGM23
      REAL(8)            :: SIGP53,SIGM53
      REAL(8)            :: D,DD
!     ******************************************************************
      SIGP=(1.D0+SIGMA)
      SIGM=(1.D0-SIGMA)
      SIGP13=SIGP**ONEBY3
      SIGM13=SIGM**ONEBY3
      SIGM23=SIGM13*SIGM13
      SIGP23=SIGP13*SIGP13
      SIGM53=SIGM*SIGM23
      SIGP53=SIGP*SIGP23
!
      D   = 0.5D0*(SIGP53+SIGM53)
      DD  = FIVEBY6*(SIGP23-SIGM23)
!
!     == 1/D AND DERIVATIVES
      ONEBYD=1.D0/SQRT(D)
      DONEBYD=-0.5D0*ONEBYD*DD/D
!
      RETURN
      END SUBROUTINE PERDEW_SPINDEP1
!
!     ..................................................................
      SUBROUTINE PERDEW_SPINDEP2(SIGMA,ONEBYD,DONEBYD,D2ONEBYD)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE SPIN DEPENDENCE OF THE PERDEW 86 GRADIENT-   **
!     **  CORRECTION FOR CORRELATION                                  **
!     **                                                              **
!     **  ONEBYD IS (1/D) WITH D DEFINED AS                           **
!     **   1/D=1/SQRT[{(1+SIGMA)**(5/3) + (1-SIGMA)**(5/3)}/2]        **
!     **  AS IN EQ.4 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **  WHERE SIGMA IS RHOT/RHOS                                    **
!     **                                                              **
!     **  REMARK: THIS ROUTINE IS A KLON OF PERDEW_SPINDEP3           **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: SIGMA
      REAL(8),INTENT(OUT):: ONEBYD
      REAL(8),INTENT(OUT):: DONEBYD
      REAL(8),INTENT(OUT):: D2ONEBYD
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER  :: FIVEBY6=5.D0/6.D0
      REAL(8),PARAMETER  :: FIVEBY9=5.D0/9.D0
      REAL(8)            :: SIGP,SIGM
      REAL(8)            :: SIGP13,SIGM13
      REAL(8)            :: SIGP23,SIGM23
      REAL(8)            :: SIGP53,SIGM53
      REAL(8)            :: D,DD,D2D
!     ******************************************************************
      SIGP=(1.D0+SIGMA)
      SIGM=(1.D0-SIGMA)
      SIGP13=SIGP**ONEBY3
      SIGM13=SIGM**ONEBY3
      SIGM23=SIGM13*SIGM13
      SIGP23=SIGP13*SIGP13
      SIGM53=SIGM*SIGM23
      SIGP53=SIGP*SIGP23
!
      D   = 0.5D0*(SIGP53+SIGM53)
      DD  = FIVEBY6*(SIGP23-SIGM23)
      D2D = FIVEBY9*(1.D0/SIGP13+1.D0/SIGM13)
!
!     == 1/D AND DERIVATIVES
      ONEBYD=1.D0/SQRT(D)
      DONEBYD=-0.5D0*ONEBYD*DD/D
      D2ONEBYD=-0.5D0*(3.D0*DONEBYD*DD+ONEBYD*D2D)/D
!
      RETURN
      END SUBROUTINE PERDEW_SPINDEP2
!
!     ..................................................................
      SUBROUTINE PERDEW_SPINDEP3(SIGMA,ONEBYD,DONEBYD,D2ONEBYD,D3ONEBYD)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE SPIN DEPENDENCE OF THE PERDEW 86 GRADIENT-   **
!     **  CORRECTION FOR CORRELATION                                  **
!     **                                                              **
!     **  ONEBYD IS (1/D) WITH D DEFINED AS                           **
!     **   1/D=1/SQRT[{(1+SIGMA)**(5/3) + (1-SIGMA)**(5/3)}/2]        **
!     **  AS IN EQ.4 OF  J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)    **
!     **  WHERE SIGMA IS RHOT/RHOS                                    **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: SIGMA
      REAL(8),INTENT(OUT):: ONEBYD
      REAL(8),INTENT(OUT):: DONEBYD
      REAL(8),INTENT(OUT):: D2ONEBYD
      REAL(8),INTENT(OUT):: D3ONEBYD
      REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
      REAL(8),PARAMETER  :: FIVEBY6=5.D0/6.D0
      REAL(8),PARAMETER  :: FIVEBY9=5.D0/9.D0
      REAL(8),PARAMETER  :: FIVEBY27=5.D0/27.D0
      REAL(8)            :: SIGP,SIGM
      REAL(8)            :: SIGP13,SIGM13
      REAL(8)            :: SIGP23,SIGM23
      REAL(8)            :: SIGP43,SIGM43
      REAL(8)            :: SIGP53,SIGM53
      REAL(8)            :: D,DD,D2D,D3D
!     ******************************************************************
      SIGP=(1.D0+SIGMA)
      SIGM=(1.D0-SIGMA)
      SIGP13=SIGP**ONEBY3
      SIGM13=SIGM**ONEBY3
      SIGM23=SIGM13*SIGM13
      SIGP23=SIGP13*SIGP13
      SIGM53=SIGM*SIGM23
      SIGP53=SIGP*SIGP23
      SIGM43=SIGM*SIGM13
      SIGP43=SIGP*SIGP13
!
      D   = 0.5D0*(SIGP53+SIGM53)
      DD  = FIVEBY6*(SIGP23-SIGM23)
      D2D = FIVEBY9*(1.D0/SIGP13+1.D0/SIGM13)
      D3D = -FIVEBY27*(1.D0/SIGP43-1.D0/SIGM43)
!
!     == 1/D AND DERIVATIVES
      ONEBYD=1.D0/SQRT(D)
      DONEBYD=-0.5D0*ONEBYD*DD/D
      D2ONEBYD=-0.5D0*(3.D0*DONEBYD*DD+ONEBYD*D2D)/D
      D3ONEBYD=-0.5D0*(5.D0*D2ONEBYD*DD+4.D0*DONEBYD*D2D+ONEBYD*D3D)/D
!
      RETURN
      END SUBROUTINE PERDEW_SPINDEP3
!
!     ..................................................................
      SUBROUTINE PERDEW$EVAL1(VAL,EXC,VXC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES PERDEWS 86 GRADIENT CORRECTION FOR CORRELATION   **
!     **    J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)                 **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    VAL=(RHOT,RHOS,GRHOT**2,DUMMY,DUMMY)                      **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: VAL(5)
      REAL(8)   ,INTENT(OUT)  :: EXC
      REAL(8)   ,INTENT(OUT)  :: VXC(5)
      REAL(8)   ,PARAMETER    :: RHOMIN=1.D-6
      REAL(8)   ,PARAMETER    :: GRHOT2MIN=1.D-30
      REAL(8)                 :: RHOT
      REAL(8)                 :: GRHOT2
      REAL(8)                 :: RHOS,SIGMA
      REAL(8)                 :: SIG_R,SIG_S
      REAL(8)                 :: ARRPE1,DARRPE1
      REAL(8)                 :: ARRPN1,DARRPN1
      REAL(8)                 :: DI,DDI
      REAL(8)                 :: DI_R,DI_S
      REAL(8)                 :: G,G_W
      REAL(8)                 :: XPHI,XPHI_R,XPHI_W
      REAL(8)                 :: A,A_R,A_S
      REAL(8)                 :: B,B_R,B_W
!     *******************************************************************
      IF(.NOT.TINI) CALL PERDEW_INITIALIZE
      RHOT=VAL(1)
      RHOS=VAL(2)
      GRHOT2=VAL(3)
      EXC=0.D0
      VXC(:)=0.D0
      IF(VAL(1).LT.RHOMIN) RETURN
      IF(GRHOT2.LT.GRHOT2MIN) RETURN
!
!     ==================================================================
!     ==  DENSITY DEPENDENCE                                          ==
!     ==================================================================
      CALL PERDEW_RHODEP1(RHOT,ARRPE1,DARRPE1,ARRPN1,DARRPN1)
!
!     ==================================================================
!     ==  GRADIENT DEPENDENCE: XPHI (EQ.9)                            ==
!     ==================================================================
      G    =SQRT(GRHOT2)
      G_W  =G/(2.D0*GRHOT2)    !D(AGRHOT)/D(G2)
      XPHI=EXP(ARRPE1*G)
      XPHI_R=G*DARRPE1*XPHI
      XPHI_W=G_W*ARRPE1*XPHI
!
!     ==================================================================
!     ==  B=W*XPHI(R,W)                                               ==
!     ==================================================================
      B     = GRHOT2*XPHI
      B_R   = GRHOT2*XPHI_R
      B_W   = XPHI+GRHOT2*XPHI_W
!
!     ==================================================================
!     ==  SPIN DEPENDENCE                                             ==
!     ==================================================================
!       == SIGMA AND DERIVATIVES WITH RESPECT RO RHOT,RHOS
        SIGMA=RHOS/RHOT
        SIG_R=-SIGMA/RHOT
        SIG_S=1.D0/RHOT
!
!       == 1/D AND DERIVATIVES WITH RESPECT TO SIGMA
        CALL PERDEW_SPINDEP1(SIGMA,DI,DDI)
!
!       == 1/D AND DERIVATIVES WITH RESPECT TO RHOT,RHOS
        DI_R  =DDI*SIG_R
        DI_S  =DDI*SIG_S
!
!       ==================================================================
!       ==  A=DI(R,S)*ARRPN1(R)                                         ==
!       ==================================================================
        A    =DI*ARRPN1
        A_R  =DI_R*ARRPN1 + DI*DARRPN1
        A_S  =DI_S*ARRPN1
!
!     ==================================================================
!     ==  ENERGY EXC(R,S,W)=A(R,S)*B(R,W)                             ==
!     ==================================================================
      EXC       =A*B
      VXC(1)    =A_R*B + A*B_R
      VXC(2)    =A_S*B
      VXC(3)    =A*B_W
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PERDEW$EVAL2(VAL,EXC,VXC,V2XC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES PERDEWS 86 GRADIENT CORRECTION FOR CORRELATION   **
!     **    J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)                 **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    VAL=(RHOT,RHOS,GRHOT**2,DUMMY,DUMMY)                      **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: VAL(5)
      REAL(8)   ,INTENT(OUT)  :: EXC
      REAL(8)   ,INTENT(OUT)  :: VXC(5)
      REAL(8)   ,INTENT(OUT)  :: V2XC(5,5)
      REAL(8)   ,PARAMETER    :: RHOMIN=1.D-5
      REAL(8)   ,PARAMETER    :: GRHOT2MIN=1.D-30
      REAL(8)                 :: RHOT
      REAL(8)                 :: GRHOT2
      REAL(8)                 :: RHOS,SIGMA
      REAL(8)                 :: SIG_R,SIG_S,SIG_RR,SIG_RS
      REAL(8)                 :: ARRPE1,DARRPE1,D2ARRPE1
      REAL(8)                 :: ARRPN1,DARRPN1,D2ARRPN1
      REAL(8)                 :: DI,DDI,D2DI
      REAL(8)                 :: DI_R,DI_S,DI_RR,DI_RS,DI_SS
      REAL(8)                 :: G,G_W,G_WW
      REAL(8)                 :: XPHI,XPHI_R,XPHI_W,XPHI_RR,XPHI_RW,XPHI_WW
      REAL(8)                 :: A,A_R,A_S,A_RR,A_RS,A_SS
      REAL(8)                 :: B,B_R,B_W,B_RR,B_RW,B_WW
!     *******************************************************************
      IF(.NOT.TINI) CALL PERDEW_INITIALIZE
      RHOT=VAL(1)
      RHOS=VAL(2)
      GRHOT2=VAL(3)
      EXC=0.D0
      VXC(:)=0.D0
      V2XC(:,:)=0.D0
      IF(VAL(1).LT.RHOMIN) RETURN
      IF(GRHOT2.LT.GRHOT2MIN) RETURN
!
!     ==================================================================
!     ==  DENSITY DEPENDENCE                                          ==
!     ==================================================================
      CALL PERDEW_RHODEP2(RHOT,ARRPE1,DARRPE1,D2ARRPE1 &
     &                        ,ARRPN1,DARRPN1,D2ARRPN1)
!
!     ==================================================================
!     ==  GRADIENT DEPENDENCE: XPHI (EQ.9)                            ==
!     ==================================================================
      G    =SQRT(GRHOT2)
      G_W  =         G/(2.D0*GRHOT2)    !D(AGRHOT)/D(G2)
      G_WW =      -G_W/(2.D0*GRHOT2)
      XPHI=EXP(ARRPE1*G)
      XPHI_R=G*DARRPE1*XPHI
      XPHI_W=G_W*ARRPE1*XPHI
      XPHI_RR=G*(D2ARRPE1*XPHI+DARRPE1*XPHI_R)
      XPHI_RW=G_W*(DARRPE1*XPHI+ARRPE1*XPHI_R)
      XPHI_WW=(G_WW*XPHI+G_W*XPHI_W)*ARRPE1
!
!     ==================================================================
!     ==  B=W*XPHI(R,W)                                               ==
!     ==================================================================
      B     = GRHOT2*XPHI
      B_R   = GRHOT2*XPHI_R
      B_W   = XPHI+GRHOT2*XPHI_W
      B_RR  = GRHOT2*XPHI_RR
      B_RW  = XPHI_R+GRHOT2*XPHI_RW
      B_WW  = 2.D0*XPHI_W+GRHOT2*XPHI_WW
!
!     ==================================================================
!     ==  SPIN DEPENDENCE                                             ==
!     ==================================================================
!       == SIGMA AND DERIVATIVES WITH RESPECT RO RHOT,RHOS
        SIGMA=RHOS/RHOT
        SIG_R=-SIGMA/RHOT
        SIG_S=1.D0/RHOT
        SIG_RR=-2.D0*SIG_R/RHOT
        SIG_RS=-SIG_S/RHOT
!
!       == 1/D AND DERIVATIVES WITH RESPECT TO SIGMA
        CALL PERDEW_SPINDEP2(SIGMA,DI,DDI,D2DI)
!
!       == 1/D AND DERIVATIVES WITH RESPECT TO RHOT,RHOS
        DI_R  =DDI*SIG_R
        DI_S  =DDI*SIG_S
        DI_RR =D2DI*SIG_R**2 + DDI*SIG_RR
        DI_SS =D2DI*SIG_S**2
        DI_RS =D2DI*SIG_R*SIG_S + DDI*SIG_RS
!
!       ==================================================================
!       ==  A=DI(R,S)*ARRPN1(R)                                         ==
!       ==================================================================
        A    =DI*ARRPN1
        A_R  =DI_R*ARRPN1 + DI*DARRPN1
        A_S  =DI_S*ARRPN1
        A_RR =DI_RR*ARRPN1 + 2.D0*DI_R*DARRPN1 + DI*D2ARRPN1
        A_RS =DI_RS*ARRPN1 +      DI_S*DARRPN1
        A_SS =DI_SS*ARRPN1
!
!     ==================================================================
!     ==  ENERGY EXC(R,S,W)=A(R,S)*B(R,W)                             ==
!     ==================================================================
      EXC       =A*B
      VXC(1)    =A_R*B + A*B_R
      VXC(2)    =A_S*B
      VXC(3)    =A*B_W
      V2XC(1,1) =A_RR*B  + 2.D0*A_R*B_R + A*B_RR
      V2XC(1,2) =A_RS*B  + A_S*B_R
      V2XC(1,3) =A_R*B_W + A*B_RW
      V2XC(2,2) =A_SS*B
      V2XC(2,3) =A_S*B_W
      V2XC(3,3) =A*B_WW
!
!     ==================================================================
!     ==  SYMMETRIZE SECOND AND THIRD DERIVATIVES                     ==
!     ==================================================================
      V2XC(2,1)=V2XC(1,2)
      V2XC(3,1)=V2XC(1,3)
      V2XC(3,2)=V2XC(2,3)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PERDEW$EVAL3(VAL,EXC,VXC,V2XC,V3XC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES PERDEWS 86 GRADIENT CORRECTION FOR CORRELATION   **
!     **    J.P. PERDEW PHYS. REV. B, 33, 8822 (1986)                 **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    VAL=(RHOT,RHOS,GRHOT**2,DUMMY,DUMMY)                      **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1996 ****
      USE PERDEW_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: VAL(5)
      REAL(8)   ,INTENT(OUT)  :: EXC
      REAL(8)   ,INTENT(OUT)  :: VXC(5)
      REAL(8)   ,INTENT(OUT)  :: V2XC(5,5)
      REAL(8)   ,INTENT(OUT)  :: V3XC(5,5,5)
      REAL(8)   ,PARAMETER    :: RHOMIN=1.D-5
      REAL(8)   ,PARAMETER    :: GRHOT2MIN=1.D-30
      REAL(8)                 :: RHOT
      REAL(8)                 :: GRHOT2
      REAL(8)                 :: RHOS,SIGMA
      REAL(8)                 :: SIG_R,SIG_S,SIG_RR,SIG_RS,SIG_RRR,SIG_RRS
      REAL(8)                 :: ARRPE1,DARRPE1,D2ARRPE1,D3ARRPE1
      REAL(8)                 :: ARRPN1,DARRPN1,D2ARRPN1,D3ARRPN1
      REAL(8)                 :: DI,DDI,D2DI,D3DI
      REAL(8)                 :: DI_R,DI_S,DI_RR,DI_RS,DI_SS
      REAL(8)                 :: DI_RRR,DI_RRS,DI_RSS,DI_SSS
      REAL(8)                 :: G,G_W,G_WW,G_WWW
      REAL(8)                 :: XPHI,XPHI_R,XPHI_W,XPHI_RR,XPHI_RW,XPHI_WW
      REAL(8)                 :: XPHI_RRR,XPHI_RRW,XPHI_RWW,XPHI_WWW
      REAL(8)                 :: A,A_R,A_S,A_RR,A_RS,A_SS,A_RRR,A_RRS,A_RSS,A_SSS
      REAL(8)                 :: B,B_R,B_W,B_RR,B_RW,B_WW,B_RRR,B_RRW,B_RWW,B_WWW
!     *******************************************************************
      IF(.NOT.TINI) CALL PERDEW_INITIALIZE
      RHOT=VAL(1)
      RHOS=VAL(2)
      GRHOT2=VAL(3)
      EXC=0.D0
      VXC(:)=0.D0
      V2XC(:,:)=0.D0
      V3XC(:,:,:)=0.D0
      IF(VAL(1).LT.RHOMIN) RETURN
      IF(GRHOT2.LT.GRHOT2MIN) RETURN
!
!     ==================================================================
!     ==  DENSITY DEPENDENCE                                          ==
!     ==================================================================
      CALL PERDEW_RHODEP3(RHOT,ARRPE1,DARRPE1,D2ARRPE1,D3ARRPE1 &
     &                        ,ARRPN1,DARRPN1,D2ARRPN1,D3ARRPN1)
!
!     ==================================================================
!     ==  GRADIENT DEPENDENCE: XPHI (EQ.9)                            ==
!     ==================================================================
      G    =SQRT(GRHOT2)
      G_W  =         G/(2.D0*GRHOT2)    !D(AGRHOT)/D(G2)
      G_WW =      -G_W/(2.D0*GRHOT2)
      G_WWW=-3.D0*G_WW/(2.D0*GRHOT2)
      XPHI=EXP(ARRPE1*G)
      XPHI_R=G*DARRPE1*XPHI
      XPHI_W=G_W*ARRPE1*XPHI
      XPHI_RR=G*(D2ARRPE1*XPHI+DARRPE1*XPHI_R)
      XPHI_RW=G_W*(DARRPE1*XPHI+ARRPE1*XPHI_R)
      XPHI_WW=(G_WW*XPHI+G_W*XPHI_W)*ARRPE1
      XPHI_RRR=G*(D3ARRPE1*XPHI+2.D0*D2ARRPE1*XPHI_R+DARRPE1*XPHI_RR)
      XPHI_RRW=G_W*(D2ARRPE1*XPHI+2.D0*DARRPE1*XPHI_R+ARRPE1*XPHI_RR)
      XPHI_RWW=(G_WW*XPHI_R+G_W*XPHI_RW)*ARRPE1+DARRPE1*(G_WW*XPHI+G_W*XPHI_W)
      XPHI_WWW=(G_WWW*XPHI+2.D0*G_WW*XPHI_W+G_W*XPHI_WW)*ARRPE1
!
!     ==================================================================
!     ==  B=W*XPHI(R,W)                                               ==
!     ==================================================================
      B     = GRHOT2*XPHI
      B_R   = GRHOT2*XPHI_R
      B_W   = XPHI+GRHOT2*XPHI_W
      B_RR  = GRHOT2*XPHI_RR
      B_RW  = XPHI_R+GRHOT2*XPHI_RW
      B_WW  = 2.D0*XPHI_W+GRHOT2*XPHI_WW
      B_RRR = GRHOT2*XPHI_RRR
      B_RRW = XPHI_RR+GRHOT2*XPHI_RRW
      B_RWW = 2.D0*XPHI_RW+GRHOT2*XPHI_RWW
      B_WWW = 3.D0*XPHI_WW+GRHOT2*XPHI_WWW
!
!     ==================================================================
!     ==  SPIN DEPENDENCE                                             ==
!     ==================================================================
!       == SIGMA AND DERIVATIVES WITH RESPECT RO RHOT,RHOS
        SIGMA=RHOS/RHOT
        SIG_R=-SIGMA/RHOT
        SIG_S=1.D0/RHOT
        SIG_RR=-2.D0*SIG_R/RHOT
        SIG_RS=-SIG_S/RHOT
        SIG_RRR=-3.D0*SIG_RR/RHOT
        SIG_RRS=-2.D0*SIG_RS/RHOT
!
!       == 1/D AND DERIVATIVES WITH RESPECT TO SIGMA
        CALL PERDEW_SPINDEP3(SIGMA,DI,DDI,D2DI,D3DI)
!
!       == 1/D AND DERIVATIVES WITH RESPECT TO RHOT,RHOS
        DI_R  =DDI*SIG_R
        DI_S  =DDI*SIG_S
        DI_RR =D2DI*SIG_R**2 + DDI*SIG_RR
        DI_SS =D2DI*SIG_S**2
        DI_RS =D2DI*SIG_R*SIG_S + DDI*SIG_RS
        DI_RRR=D3DI*SIG_R**3 + 3.D0*D2DI*SIG_R*SIG_RR + DDI*SIG_RRR
        DI_RRS=D3DI*SIG_R**2*SIG_S &
     &        +D2DI*(SIG_RR*SIG_S+2.D0*SIG_R*SIG_RS) &
     &        +DDI*SIG_RRS
        DI_RSS=D3DI*SIG_R*SIG_S**2 + 2.D0*D2DI*SIG_S*SIG_RS
        DI_SSS=D3DI*SIG_S**3
!
!       ==================================================================
!       ==  A=DI(R,S)*ARRPN1(R)                                         ==
!       ==================================================================
        A    =DI*ARRPN1
        A_R  =DI_R*ARRPN1 + DI*DARRPN1
        A_S  =DI_S*ARRPN1
        A_RR =DI_RR*ARRPN1 + 2.D0*DI_R*DARRPN1 + DI*D2ARRPN1
        A_RS =DI_RS*ARRPN1 +      DI_S*DARRPN1
        A_SS =DI_SS*ARRPN1
        A_RRR=DI_RRR*ARRPN1 + 3.D0*DI_RR*DARRPN1 + 3.D0*DI_R*D2ARRPN1+DI*D3ARRPN1
        A_RRS=DI_RRS*ARRPN1 + 2.D0*DI_RS*DARRPN1 +      DI_S*D2ARRPN1
        A_RSS=DI_RSS*ARRPN1 +      DI_SS*DARRPN1
        A_SSS=DI_SSS*ARRPN1
!
!     ==================================================================
!     ==  ENERGY EXC(R,S,W)=A(R,S)*B(R,W)                             ==
!     ==================================================================
      EXC       =A*B
      VXC(1)    =A_R*B + A*B_R
      VXC(2)    =A_S*B
      VXC(3)    =A*B_W
      V2XC(1,1) =A_RR*B  + 2.D0*A_R*B_R + A*B_RR
      V2XC(1,2) =A_RS*B  + A_S*B_R
      V2XC(1,3) =A_R*B_W + A*B_RW
      V2XC(2,2) =A_SS*B
      V2XC(2,3) =A_S*B_W
      V2XC(3,3) =A*B_WW
      V3XC(1,1,1)=A_RRR*B  + 3.D0*A_RR*B_R + 3.D0*A_R*B_RR + A*B_RRR
      V3XC(1,1,2)=A_RRS*B  + 2.D0*A_RS*B_R + A_S*B_RR
      V3XC(1,1,3)=A_RR*B_W + 2.D0*A_R*B_RW + A*B_RRW
      V3XC(1,2,2)=A_RSS*B  + A_SS*B_R
      V3XC(1,2,3)=A_RS*B_W + A_S*B_RW
      V3XC(1,3,3)=A_R*B_WW + A*B_RWW
      V3XC(2,2,2)=A_SSS*B
      V3XC(2,2,3)=A_SS*B_W
      V3XC(2,3,3)=A_S*B_WW
      V3XC(3,3,3)=A*B_WWW
!
!     ==================================================================
!     ==  SYMMETRIZE SECOND AND THIRD DERIVATIVES                     ==
!     ==================================================================
      V2XC(2,1)=V2XC(1,2)
      V2XC(3,1)=V2XC(1,3)
      V2XC(3,2)=V2XC(2,3)
      V3XC(1,2,1)=V3XC(1,1,2)
      V3XC(2,1,1)=V3XC(1,1,2)
      V3XC(1,3,1)=V3XC(1,1,3)
      V3XC(3,1,1)=V3XC(1,1,3)
      V3XC(2,1,2)=V3XC(1,2,2)
      V3XC(2,2,1)=V3XC(1,2,2)
      V3XC(1,3,2)=V3XC(1,2,3)
      V3XC(2,1,3)=V3XC(1,2,3)
      V3XC(2,3,1)=V3XC(1,2,3)
      V3XC(3,1,2)=V3XC(1,2,3)
      V3XC(3,2,1)=V3XC(1,2,3)
      V3XC(3,1,3)=V3XC(1,3,3)
      V3XC(3,3,1)=V3XC(1,3,3)
      V3XC(2,3,2)=V3XC(2,2,3)
      V3XC(3,2,2)=V3XC(2,2,3)
      V3XC(3,2,3)=V3XC(2,3,3)
      V3XC(3,3,2)=V3XC(2,3,3)
      RETURN
      END
!
!.......................................................................
MODULE PERDEWWANG91L_MODULE
!***********************************************************************      
!**                                                                   **
!** NAME: PERDEWWANG91L                                               **
!**                                                                   **
!** PURPOSE: IMPLEMENTS THE PERDEW WANG PARAMETERIZATION              **
!**   OF CEPERLEY ALDERS MONTE CARLO SIMULATIONS OF THE               **
!**   HOMOGENEOUS ELECTRON GAS                                        **
!**   J.P.PERDEW AND Y. WANG, PHYS. REV. B45, 13244 (1992-I)          **
!**                                                                   **
!** FUNCTIONS:                                                        **
!**   PERDEWWANG91L$EPSVAL(RHOT,RHOS,EPS,EPS_T,EPS_S)                 **
!**   PERDEWWANG91L$EPSVAL2(RHOT,RHOS,EPS,EPS_T,EPS_S &               **
!**                        EPS_TT,EPS_TS,EPS_SS)                      **
!**   PERDEWWANG91L$EPSVAL3(RHOT,RHOS,EPS,EPS_T,EPS_S &               **
!**                        ,EPS_TT,EPS_TS,EPS_SS                      **
!**                        ,EC_TTT,EC_TTS,EC_TSS,EC_SSS)              **
!**   PERDEWWANG91L$EVAL(VAL,EXC,DEXC,D2EXC)                          **
!**   PERDEWWANG91L$EVAL2(VAL,EXC,DEXC,D2EXC)                         **
!**   PERDEWWANG91L$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)                   **
!**                                                                   **
!** DEPENDENCIES:                                                     **
!**   1) INHERITS TSPIN FROM DFT_MODULE                               **
!**   2) TABLE1D <- ERROR                                             **
!***********************************************************************
USE DFT_MODULE, ONLY : TSPIN,TSAFE
USE TABLE1D_MODULE
LOGICAL(4)         :: TINI=.FALSE.
!== PARAMETERS FOR CORRELATION WITH SIGMA=0
!== NOTE: AU HAS ONE MORE DIGIT COMPARED TO THE PUBLICATION
REAL(8),PARAMETER  :: AU     =3.10907D-2  
REAL(8),PARAMETER  :: ALPHA1U=2.1370D-1
REAL(8),PARAMETER  :: BETA1U =7.5957D0
REAL(8),PARAMETER  :: BETA2U =3.5876D0
REAL(8),PARAMETER  :: BETA3U =1.6382D0
REAL(8),PARAMETER  :: BETA4U =4.9294D-1
!== PARAMETERS FOR CORRELATION WITH SIGMA=1
!== NOTE: AP HAS ONE MORE DIGIT COMPARED TO THE PUBLICATION
REAL(8),PARAMETER  :: AP     =1.554535D-2
REAL(8),PARAMETER  :: ALPHA1P=2.0548D-1
REAL(8),PARAMETER  :: BETA1P =1.41189D+1
REAL(8),PARAMETER  :: BETA2P =6.1977D0
REAL(8),PARAMETER  :: BETA3P =3.3662D0
REAL(8),PARAMETER  :: BETA4P =6.2517D-1
REAL(8),PARAMETER  :: AX     =1.68869D-2
REAL(8),PARAMETER  :: ALPHA1X=1.1125D-1
REAL(8),PARAMETER  :: BETA1X =10.357D0
REAL(8),PARAMETER  :: BETA2X =3.6231D0
REAL(8),PARAMETER  :: BETA3X =8.8026D-1
REAL(8),PARAMETER  :: BETA4X =4.9671D-1
REAL(8),PARAMETER  :: FSIGMAPP=1.709920934D0
REAL(8)            :: FFAC
REAL(8)            :: DFFAC
REAL(8)            :: D2FFAC
REAL(8)            :: D3FFAC
REAL(8)            :: RSFAC    !RS=RSFAC/RHO**(1/3)
!======================================================================
!== TABLE TO ACCELERATE THE CALCULATION                              ==
!======================================================================
INTEGER(4),PARAMETER :: NRHO=5000
TYPE (XTABLE_TYPE)   :: RHOTABLE
REAL(8)              :: ECUARRAY(NRHO)
REAL(8)              :: DECUARRAY(NRHO)
REAL(8)              :: D2ECUARRAY(NRHO)
REAL(8)              :: D3ECUARRAY(NRHO)
REAL(8)              :: ECPARRAY(NRHO)
REAL(8)              :: DECPARRAY(NRHO)
REAL(8)              :: D2ECPARRAY(NRHO)
REAL(8)              :: D3ECPARRAY(NRHO)
REAL(8)              :: ECXARRAY(NRHO)
REAL(8)              :: DECXARRAY(NRHO)
REAL(8)              :: D2ECXARRAY(NRHO)
REAL(8)              :: D3ECXARRAY(NRHO)
!**********************************************************************
CONTAINS
!      .................................................................
       SUBROUTINE PERDEWWANG91L_INITIALIZE
!      *****************************************************************
!      ** CALCULATE SOME FACTORS FOR LATER USE                        **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8)       :: PI
       REAL(8)       :: RHOT
       REAL(8)       :: ECU,DECU,D2ECU,D3ECU
       REAL(8)       :: ECP,DECP,D2ECP,D3ECP
       REAL(8)       :: ECX,DECX,D2ECX,D3ECX
       INTEGER(4)    :: I
!      *****************************************************************
       PI=4.D0*ATAN(1.D0)
       RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0)
       FFAC=1.D0/(2.D0**(4.D0/3.D0)-2.D0)
       DFFAC=4.D0/3.D0*FFAC
       D2FFAC=DFFAC/3.D0
       D3FFAC=-2.D0*D2FFAC/3.D0
!      ================================================================
!      ==  INITIALIZE TABLES                                         ==
!      ================================================================
       CALL TABLE1D$MAKE(RHOTABLE,1.D-3,1.D-1,NRHO)
       DO I=1,NRHO
         CALL TABLE1D$XOFI(RHOTABLE,I,RHOT)
         CALL PERDEWWANG91L_RHODEP3(RHOT,ECU,DECU,D2ECU,D3ECU &
      &                         ,AU,ALPHA1U,BETA1U,BETA2U,BETA3U,BETA4U)
         ECUARRAY(I)  =ECU
         DECUARRAY(I) =DECU
         D2ECUARRAY(I)=D2ECU
         D3ECUARRAY(I)=D3ECU
         CALL PERDEWWANG91L_RHODEP3(RHOT,ECP,DECP,D2ECP,D3ECP &
      &                         ,AP,ALPHA1P,BETA1P,BETA2P,BETA3P,BETA4P)
         ECPARRAY(I)  =ECP
         DECPARRAY(I) =DECP
         D2ECPARRAY(I)=D2ECP
         D3ECPARRAY(I)=D3ECP
         CALL PERDEWWANG91L_RHODEP3(RHOT,ECX,DECX,D2ECX,D3ECX &
      &                         ,AX,ALPHA1X,BETA1X,BETA2X,BETA3X,BETA4X)
         ECXARRAY(I)  =-ECX/FSIGMAPP    !ECX=ALPHA/FSIGMAPP
         DECXARRAY(I) =-DECX/FSIGMAPP    !ECX=ALPHA/FSIGMAPP
         D2ECXARRAY(I)=-D2ECX/FSIGMAPP    !ECX=ALPHA/FSIGMAPP
         D3ECXARRAY(I)=-D3ECX/FSIGMAPP    !ECX=ALPHA/FSIGMAPP
       ENDDO
       TINI=.TRUE.
       RETURN
       END SUBROUTINE PERDEWWANG91L_INITIALIZE
!      .................................................................
       SUBROUTINE PERDEWWANG91L_RHODEP1(RHOT,EC,DEC &
      &                ,A,ALPHA1,BETA1,BETA2,BETA3,BETA4)
!      *****************************************************************
!      **  IMPLEMENTS FUNCTION DEFINED IN EQ.10 OF PERDEW-WANG PAPER **
!      **  PW,PRB45,13224(1992-I)                                     **
!      **  HERE AS FUNCTION OF RHOT AND NOT AS FUNCTION OF RS         **
!      **  NOTE: RHOT IS A MODULE ROUTINE, FROM WHICH IT USES THE     **
!      **  DEFINITION OF RSFAC                                        **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: A
       REAL(8),INTENT(IN) :: ALPHA1
       REAL(8),INTENT(IN) :: BETA1
       REAL(8),INTENT(IN) :: BETA2
       REAL(8),INTENT(IN) :: BETA3
       REAL(8),INTENT(IN) :: BETA4
       REAL(8),INTENT(OUT):: EC
       REAL(8),INTENT(OUT):: DEC
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8)            :: RS,DRS
       REAL(8)            :: SQRS,RSSQRS
       REAL(8)            :: SVAR1,DSVAR1
       REAL(8)            :: SVAR2,DSVAR2
       REAL(8)            :: SVAR3,DSVAR3
       REAL(8)            :: SVAR4,DSVAR4
!      *****************************************************************
       RS=RSFAC/RHOT**ONEBY3
!
!      =================================================================
!      == CALCULATE DERIVATIVE WITH RESPECT TO RS                    ==
!      =================================================================
       SQRS=SQRT(RS)
       RSSQRS=RS*SQRS
!      =================================================================
       SVAR1  = -2.D0*A*(1.D0+ALPHA1*RS)
       DSVAR1 = -2.D0*A*ALPHA1
!      =================================================================
       SVAR2  =2.D0*A*SQRS*(BETA1+BETA2*SQRS+BETA3*RS+BETA4*RSSQRS)
       DSVAR2 =A/SQRS*(BETA1+2.D0*BETA2*SQRS+3.D0*BETA3*RS+4.D0*BETA4*RSSQRS)
!      =================================================================
       SVAR3  =1.D0+1.D0/SVAR2
       DSVAR3 =-DSVAR2/SVAR2**2
!      =================================================================
       SVAR4  =LOG(SVAR3)
       DSVAR4 =DSVAR3/SVAR3
!      =================================================================
       EC=SVAR1*SVAR4
       DEC=DSVAR1*SVAR4+SVAR1*DSVAR4
!
!      == TRANSFORM FROM EC(RS) TO EC(RHOT)==============================
       DRS=-ONEBY3*RS/RHOT
       DEC =DEC*DRS
       RETURN
       END SUBROUTINE PERDEWWANG91L_RHODEP1
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L_RHODEP2(RHOT,EC,DEC,D2EC &
      &                ,A,ALPHA1,BETA1,BETA2,BETA3,BETA4)
!      *****************************************************************
!      **  IMPLEMENTS FUNCTION DEFINED IN EQ.10 OF PERDEW-WANG PAPER  **
!      **  PW,PRB45,13224(1992-I)                                     **
!      **  HERE AS FUNCTION OF RHOT AND NOT AS FUNCTION OF RS         **
!      **  NOTE: RHOT IS A MODULE ROUTINE, FROM WHICH IT USES THE     **
!      **  DEFINITION OF RSFAC                                        **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: A
       REAL(8),INTENT(IN) :: ALPHA1
       REAL(8),INTENT(IN) :: BETA1
       REAL(8),INTENT(IN) :: BETA2
       REAL(8),INTENT(IN) :: BETA3
       REAL(8),INTENT(IN) :: BETA4
       REAL(8),INTENT(OUT):: EC
       REAL(8),INTENT(OUT):: DEC
       REAL(8),INTENT(OUT):: D2EC
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8),PARAMETER  :: FOURBY9=4.D0/9.D0
       REAL(8)            :: RS,DRS,D2RS
       REAL(8)            :: SQRS,RSSQRS
       REAL(8)            :: SVAR1,DSVAR1
       REAL(8)            :: SVAR2,DSVAR2,D2SVAR2
       REAL(8)            :: SVAR3,DSVAR3,D2SVAR3
       REAL(8)            :: SVAR4,DSVAR4,D2SVAR4
!      *****************************************************************
       RS=RSFAC/RHOT**ONEBY3
!
!      =================================================================
!      == CALCULATE DERIVATIVE WITH RESPECT TO RS                    ==
!      =================================================================
       SQRS=SQRT(RS)
       RSSQRS=RS*SQRS
!      =================================================================
       SVAR1  = -2.D0*A*(1.D0+ALPHA1*RS)
       DSVAR1 = -2.D0*A*ALPHA1
!      =================================================================
       SVAR2  =2.D0*A*SQRS*(BETA1+BETA2*SQRS+BETA3*RS+BETA4*RSSQRS)
       DSVAR2 =A/SQRS*(BETA1+2.D0*BETA2*SQRS+3.D0*BETA3*RS+4.D0*BETA4*RSSQRS)
       D2SVAR2=0.5D0*A/RSSQRS*(-BETA1+3.D0*BETA3*RS+8.D0*BETA4*RSSQRS)
!      =================================================================
       SVAR3  =1.D0+1.D0/SVAR2
       DSVAR3 =-DSVAR2/SVAR2**2
       D2SVAR3=-D2SVAR2/SVAR2**2+2.D0*SVAR2*DSVAR3**2
!      =================================================================
       SVAR4  =LOG(SVAR3)
       DSVAR4 =DSVAR3/SVAR3
       D2SVAR4=D2SVAR3/SVAR3-DSVAR4**2
!      =================================================================
       EC=SVAR1*SVAR4
       DEC=DSVAR1*SVAR4+SVAR1*DSVAR4
       D2EC=2.D0*DSVAR1*DSVAR4+SVAR1*D2SVAR4
!
!      == TRANSFORM FROM EC(RS) TO EC(RHOT)==============================
       DRS=-ONEBY3*RS/RHOT
       D2RS=FOURBY9*RS/RHOT**2
       D2EC=D2EC*DRS*DRS+DEC*D2RS
       DEC =DEC*DRS
       RETURN
       END SUBROUTINE PERDEWWANG91L_RHODEP2
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L_RHODEP3(RHOT,EC,DEC,D2EC,D3EC &
      &                ,A,ALPHA1,BETA1,BETA2,BETA3,BETA4)
!      *****************************************************************
!      **  IMPLEMENTS FUNCTION DEFINED IN EQ.10 OF PERDEW-WANG PAPER  **
!      **  PW,PRB45,13224(1992-I)                                     **
!      **  HERE AS FUNCTION OF RHOT AND NOT AS FUNCTION OF RS         **
!      **  NOTE: RHOT IS A MODULE ROUTINE, FROM WHICH IT USES THE     **
!      **  DEFINITION OF RSFAC                                        **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: A
       REAL(8),INTENT(IN) :: ALPHA1
       REAL(8),INTENT(IN) :: BETA1
       REAL(8),INTENT(IN) :: BETA2
       REAL(8),INTENT(IN) :: BETA3
       REAL(8),INTENT(IN) :: BETA4
       REAL(8),INTENT(OUT):: EC
       REAL(8),INTENT(OUT):: DEC
       REAL(8),INTENT(OUT):: D2EC
       REAL(8),INTENT(OUT):: D3EC
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8),PARAMETER  :: FOURBY9=4.D0/9.D0
       REAL(8)            :: RS,DRS,D2RS,D3RS
       REAL(8)            :: SQRS,RSSQRS
       REAL(8)            :: SVAR1,DSVAR1
       REAL(8)            :: SVAR2,DSVAR2,D2SVAR2,D3SVAR2
       REAL(8)            :: SVAR3,DSVAR3,D2SVAR3,D3SVAR3
       REAL(8)            :: SVAR4,DSVAR4,D2SVAR4,D3SVAR4
!      *****************************************************************
       RS=RSFAC/RHOT**ONEBY3
!
!      =================================================================
!      == CALCULATE DERIVATIVE WITH RESPECT TO RS                    ==
!      =================================================================
       SQRS=SQRT(RS)
       RSSQRS=RS*SQRS
!      =================================================================
       SVAR1  = -2.D0*A*(1.D0+ALPHA1*RS)
       DSVAR1 = -2.D0*A*ALPHA1
!      =================================================================
       SVAR2  =2.D0*A*SQRS*(BETA1+BETA2*SQRS+BETA3*RS+BETA4*RSSQRS)
       DSVAR2 =A/SQRS*(BETA1+2.D0*BETA2*SQRS+3.D0*BETA3*RS+4.D0*BETA4*RSSQRS)
       D2SVAR2=0.5D0*A/RSSQRS*(-BETA1+3.D0*BETA3*RS+8.D0*BETA4*RSSQRS)
       D3SVAR2=0.75D0*A/RSSQRS*(BETA1/RS-BETA3)
!      =================================================================
       SVAR3  =1.D0+1.D0/SVAR2
       DSVAR3 =-DSVAR2/SVAR2**2
       D2SVAR3=-D2SVAR2/SVAR2**2+2.D0*SVAR2*DSVAR3**2
       D3SVAR3=-D3SVAR2/SVAR2**2+2.D0*D2SVAR2*DSVAR2/SVAR2**3 &
      &        +2.D0*DSVAR2*DSVAR3**2+4.D0*SVAR2*DSVAR3*D2SVAR3
!      =================================================================
       SVAR4  =LOG(SVAR3)
       DSVAR4 =DSVAR3/SVAR3
       D2SVAR4=D2SVAR3/SVAR3-DSVAR4**2
       D3SVAR4=D3SVAR3/SVAR3-D2SVAR3*DSVAR3/SVAR3**2-2.D0*DSVAR4*D2SVAR4
!      =================================================================
       EC=SVAR1*SVAR4
       DEC=DSVAR1*SVAR4+SVAR1*DSVAR4
       D2EC=2.D0*DSVAR1*DSVAR4+SVAR1*D2SVAR4
       D3EC=2.D0*DSVAR1*D2SVAR4+DSVAR1*D2SVAR4+SVAR1*D3SVAR4
!
!      == TRANSFORM FROM EC(RS) TO EC(RHOT)==============================
       DRS=-ONEBY3*RS/RHOT
       D2RS=FOURBY9*RS/RHOT**2
       D3RS=-28.D0/27.D0*RS/RHOT**3
!
       D3EC=D3EC*DRS**3+3.D0*D2EC*DRS*D2RS+DEC*D3RS
       D2EC=D2EC*DRS*DRS+DEC*D2RS
       DEC =DEC*DRS
       RETURN
       END SUBROUTINE PERDEWWANG91L_RHODEP3
       END MODULE PERDEWWANG91L_MODULE
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L$EPSVAL(RHOT,RHOS,EC,EC_T,EC_S)
!      *****************************************************************
!      ** J.P.PERDEW AND Y. WANG, PHYS. REV. B45, 13244 (1992-I)      **
!      *****************************************************************
       USE PERDEWWANG91L_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: RHOS
       REAL(8),INTENT(OUT):: EC
       REAL(8),INTENT(OUT):: EC_T
       REAL(8),INTENT(OUT):: EC_S
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8)            :: SIGMA
       REAL(8)            :: SIGMAP13,SIGMAP23,SIGMAP43
       REAL(8)            :: SIGMAM13,SIGMAM23,SIGMAM43
       REAL(8)            :: FSIGMA,DFSIGMA
       REAL(8)            :: ECU,DECU
       REAL(8)            :: ECP,DECP
       REAL(8)            :: ECX,DECX
       REAL(8)            :: SIGMA2,SIGMA3,SIGMA4
       REAL(8)            :: SVAR1,SVAR1_T,SVAR1_S
       REAL(8)            :: EC1,EC1_T,EC1_S
       REAL(8)            :: SIGMA_T,SIGMA_S
       LOGICAL(4)         :: TTABLE
!      *****************************************************************
       IF(.NOT.TINI) CALL PERDEWWANG91L_INITIALIZE
!
!      =================================================================
!      == DECIDE ABOUT TABLE LOOKUP OR NOT                            ==
!      =================================================================
       TTABLE=.NOT.TSAFE
       IF(TTABLE) THEN
         CALL TABLE1D$SETX(RHOTABLE,RHOT,TTABLE)
       END IF
!
!      =================================================================
!      == GET RHO-DEPENCE FOR THE NON-SPIN POLARIZED ELECTRON GAS     ==
!      =================================================================
       IF(TTABLE) THEN
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECUARRAY,ECU)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECUARRAY,DECU)
       ELSE
         CALL PERDEWWANG91L_RHODEP1(RHOT,ECU,DECU &
      &                         ,AU,ALPHA1U,BETA1U,BETA2U,BETA3U,BETA4U)
       ENDIF
!
!      =================================================================
!      == RETURN EARLY FOR NON-SPIN POLARIZED CALCULATIONS            ==
!      =================================================================
       IF(.NOT.TSPIN) THEN
         EC=ECU
         EC_T=DECU     
         EC_S=0.D0
         RETURN
       END IF
!
!      =================================================================
!      == CALCULATE EPS_XC FOR FULLY POLARIZED ELECTRON GAS           ==
!      == AND THE SPIN STIFFNESS                                      ==
!      =================================================================
       IF(TTABLE) THEN
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECPARRAY,ECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECPARRAY,DECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECXARRAY,ECX)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECXARRAY,DECX)
       ELSE
         CALL PERDEWWANG91L_RHODEP1(RHOT,ECP,DECP &
      &                         ,AP,ALPHA1P,BETA1P,BETA2P,BETA3P,BETA4P)
         CALL PERDEWWANG91L_RHODEP1(RHOT,ECX,DECX &
      &                         ,AX,ALPHA1X,BETA1X,BETA2X,BETA3X,BETA4X)
!        == NOTE THAT ECX IS -ALPHA_C FROM THE PAPER !!
         ECX =-ECX/FSIGMAPP    !ECX=ALPHA/FSIGMAPP
         DECX=-DECX/FSIGMAPP
       END IF
!      == NOW, ECX IS +ALPHA_C/FSIGMAPP ...
!
!      =================================================================
!      == CALCULATE SPIN DEPENDENCE                                   ==
!      =================================================================
       SIGMA=RHOS/RHOT
       SIGMAP13=(1.D0+SIGMA)**ONEBY3
       SIGMAP23=SIGMAP13*SIGMAP13
       SIGMAP43=SIGMAP23*SIGMAP23
       SIGMAM13=(1.D0-SIGMA)**ONEBY3
       SIGMAM23=SIGMAM13*SIGMAM13
       SIGMAM43=SIGMAM23*SIGMAM23
       FSIGMA  = FFAC *(SIGMAP43  + SIGMAM43 - 2.D0)
       DFSIGMA = DFFAC*(SIGMAP13  - SIGMAM13)
!
!      =================================================================
!      == FORM DERIVATIVES FIRST WITH RESPECT TO RHOT AND SIGMA       ==
!      =================================================================
       SIGMA2=SIGMA*SIGMA
       SIGMA3=SIGMA2*SIGMA
       SIGMA4=SIGMA2*SIGMA2
!
!      == NOTE THAT ECX IS -ALPHA_C FROM THE PAPER !!
       SVAR1  =(ECP-ECU)  *SIGMA4+ECX*(1.D0-SIGMA4)
       SVAR1_T=(DECP-DECU)*SIGMA4+DECX*(1.D0-SIGMA4)
       SVAR1_S=(ECP-ECU-ECX)*4.D0*SIGMA3
!
       EC1=ECU+FSIGMA*SVAR1
       EC1_T=DECU+FSIGMA*SVAR1_T
       EC1_S=DFSIGMA*SVAR1 + FSIGMA*SVAR1_S
!
!      =================================================================
!      == TRANSFORM DERIVATIVES TO RHOT AND RHOS                      ==
!      =================================================================
       SIGMA_S =1.D0/RHOT
       SIGMA_T =-SIGMA*SIGMA_S
!
       EC   =EC1
       EC_T =EC1_T+EC1_S*SIGMA_T
       EC_S =EC1_S*SIGMA_S
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L$EPSVAL2(RHOT,RHOS,EC,EC_T,EC_S &
      &                        ,EC_TT,EC_TS,EC_SS)
!      *****************************************************************
!      ** J.P.PERDEW AND Y. WANG, PHYS. REV. B45, 13244 (1992-I)      **
!      *****************************************************************
       USE PERDEWWANG91L_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: RHOS
       REAL(8),INTENT(OUT):: EC
       REAL(8),INTENT(OUT):: EC_T
       REAL(8),INTENT(OUT):: EC_S
       REAL(8),INTENT(OUT):: EC_TT
       REAL(8),INTENT(OUT):: EC_TS
       REAL(8),INTENT(OUT):: EC_SS
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8)            :: SIGMA
       REAL(8)            :: SIGMAP13,SIGMAP23,SIGMAP43
       REAL(8)            :: SIGMAM13,SIGMAM23,SIGMAM43
       REAL(8)            :: FSIGMA,DFSIGMA,D2FSIGMA
       REAL(8)            :: ECU,DECU,D2ECU
       REAL(8)            :: ECP,DECP,D2ECP
       REAL(8)            :: ECX,DECX,D2ECX
       REAL(8)            :: SIGMA2,SIGMA3,SIGMA4
       REAL(8)            :: SVAR1,SVAR1_T,SVAR1_S,SVAR1_TT,SVAR1_TS,SVAR1_SS
       REAL(8)            :: EC1,EC1_T,EC1_S,EC1_TT,EC1_TS,EC1_SS
       REAL(8)            :: SIGMA_T,SIGMA_S,SIGMA_TT,SIGMA_TS
       LOGICAL(4)         :: TTABLE
!      *****************************************************************
       IF(.NOT.TINI) CALL PERDEWWANG91L_INITIALIZE
!
!      =================================================================
!      == DECIDE ABOUT TABLE LOOKUP OR NOT                            ==
!      =================================================================
       TTABLE=.NOT.TSAFE
       IF(TTABLE) THEN
         CALL TABLE1D$SETX(RHOTABLE,RHOT,TTABLE)
       END IF
!
!      =================================================================
!      == GET RHO-DEPENCE FOR THE NON-SPIN POLARIZED ELECTRON GAS     ==
!      =================================================================
       IF(TTABLE) THEN
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECUARRAY,ECU)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECUARRAY,DECU)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D2ECUARRAY,D2ECU)
       ELSE
         CALL PERDEWWANG91L_RHODEP2(RHOT,ECU,DECU,D2ECU &
      &                         ,AU,ALPHA1U,BETA1U,BETA2U,BETA3U,BETA4U)
       END IF
!
!      =================================================================
!      == RETURN EARLY FOR NON-SPIN POLARIZED CALCULATIONS            ==
!      =================================================================
       IF(.NOT.TSPIN) THEN
         EC=ECU
         EC_T=DECU     
         EC_S=0.D0
         EC_TT=D2ECU
         EC_TS=0.D0
         EC_SS=0.D0
         RETURN
       END IF
!
!      =================================================================
!      == CALCULATE EPS_XC FOR FULLY POLARIZED ELECTRON GAS           ==
!      == AND THE SPIN STIFFNESS                                      ==
!      =================================================================
       IF(TTABLE) THEN
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECPARRAY,ECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECPARRAY,DECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D2ECPARRAY,D2ECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECXARRAY,ECX)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECXARRAY,DECX)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D2ECXARRAY,D2ECX)
       ELSE
         CALL PERDEWWANG91L_RHODEP2(RHOT,ECP,DECP,D2ECP &
      &                         ,AP,ALPHA1P,BETA1P,BETA2P,BETA3P,BETA4P)
!      == NOTE THAT ECX IS -ALPHA_C FOMR THE PAPER !!
         CALL PERDEWWANG91L_RHODEP2(RHOT,ECX,DECX,D2ECX &
      &                         ,AX,ALPHA1X,BETA1X,BETA2X,BETA3X,BETA4X)
         ECX=-ECX/FSIGMAPP
         DECX=-DECX/FSIGMAPP
         D2ECX=-D2ECX/FSIGMAPP
       END IF
!
!      =================================================================
!      == CALCULATE SPIN DEPENDENCE                                   ==
!      =================================================================
       SIGMA=RHOS/RHOT
       SIGMAP13=(1.D0+SIGMA)**ONEBY3
       SIGMAP23=SIGMAP13*SIGMAP13
       SIGMAP43=SIGMAP23*SIGMAP23
       SIGMAM13=(1.D0-SIGMA)**ONEBY3
       SIGMAM23=SIGMAM13*SIGMAM13
       SIGMAM43=SIGMAM23*SIGMAM23
       FSIGMA  = FFAC*(SIGMAP43  + SIGMAM43 - 2.D0)
       DFSIGMA = DFFAC*(SIGMAP13  - SIGMAM13)
       D2FSIGMA= D2FFAC*(1.D0/SIGMAP23+ 1.D0/SIGMAM23) 
!
!      =================================================================
!      == FORM DERIVATIVES FIRST WITH RESPECT TO RHOT AND SIGMA       ==
!      =================================================================
       SIGMA2=SIGMA*SIGMA
       SIGMA3=SIGMA2*SIGMA
       SIGMA4=SIGMA3*SIGMA
!
       SVAR1   =(ECP-ECU)  *SIGMA4+ECX*(1.D0-SIGMA4)
       SVAR1_T =(DECP-DECU)*SIGMA4+DECX*(1.D0-SIGMA4)
       SVAR1_S =(ECP-ECU-ECX)*4.D0*SIGMA3
       SVAR1_SS=(ECP-ECU-ECX)*12.D0*SIGMA2
       SVAR1_TS=(DECP-DECU-DECX)*4.D0*SIGMA3
       SVAR1_TT=(D2ECP-D2ECU)*SIGMA4+D2ECX*(1.D0-SIGMA4)
!
       EC1=ECU+FSIGMA*SVAR1
       EC1_T=DECU+FSIGMA*SVAR1_T
       EC1_S=DFSIGMA*SVAR1 + FSIGMA*SVAR1_S
       EC1_TT=D2ECU+FSIGMA*SVAR1_TT
       EC1_TS=DFSIGMA*SVAR1_T+FSIGMA*SVAR1_TS
       EC1_SS=D2FSIGMA*SVAR1 + 2.D0*DFSIGMA*SVAR1_S+FSIGMA*SVAR1_SS
!
!      =================================================================
!      == TRANSFORM DERIVATIVES TO RHOT AND RHOS                      ==
!      =================================================================
       SIGMA_S =1.D0/RHOT
       SIGMA_T =-SIGMA*SIGMA_S
       SIGMA_TS=-SIGMA_S**2
       SIGMA_TT=-2.D0*SIGMA_T*SIGMA_S
!
       EC   =EC1
       EC_T =EC1_T+EC1_S*SIGMA_T
       EC_S =EC1_S*SIGMA_S
       EC_TT=EC1_TT+2.D0*EC1_TS*SIGMA_T+EC1_SS*SIGMA_T**2+EC1_S*SIGMA_TT
       EC_TS=EC1_TS*SIGMA_S+EC1_SS*SIGMA_S*SIGMA_T+EC1_S*SIGMA_TS
       EC_SS=EC1_SS*SIGMA_S**2
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L$EPSVAL3(RHOT,RHOS,EC,EC_T,EC_S &
      &                  ,EC_TT,EC_TS,EC_SS,EC_TTT,EC_TTS,EC_TSS,EC_SSS)
!      *****************************************************************
!      ** J.P.PERDEW AND Y. WANG, PHYS. REV. B45, 13244 (1992-I)      **
!      *****************************************************************
       USE PERDEWWANG91L_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: RHOS
       REAL(8),INTENT(OUT):: EC
       REAL(8),INTENT(OUT):: EC_T
       REAL(8),INTENT(OUT):: EC_S
       REAL(8),INTENT(OUT):: EC_TT
       REAL(8),INTENT(OUT):: EC_TS
       REAL(8),INTENT(OUT):: EC_SS
       REAL(8),INTENT(OUT):: EC_TTT
       REAL(8),INTENT(OUT):: EC_TTS
       REAL(8),INTENT(OUT):: EC_TSS
       REAL(8),INTENT(OUT):: EC_SSS
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8)            :: SIGMA
       REAL(8)            :: SIGMAP,SIGMAP13,SIGMAP23,SIGMAP43,SIGMAP53
       REAL(8)            :: SIGMAM,SIGMAM13,SIGMAM23,SIGMAM43,SIGMAM53
       REAL(8)            :: FSIGMA,DFSIGMA,D2FSIGMA,D3FSIGMA
       REAL(8)            :: ECU,DECU,D2ECU,D3ECU
       REAL(8)            :: ECP,DECP,D2ECP,D3ECP
       REAL(8)            :: ECX,DECX,D2ECX,D3ECX
       REAL(8)            :: SIGMA2,SIGMA3,SIGMA4
       REAL(8)            :: SVAR1,SVAR1_T,SVAR1_S,SVAR1_TT,SVAR1_TS,SVAR1_SS
       REAL(8)            :: SVAR1_TTT,SVAR1_TTS,SVAR1_TSS,SVAR1_SSS
       REAL(8)            :: EC1,EC1_T,EC1_S,EC1_TT,EC1_TS,EC1_SS
       REAL(8)            :: EC1_TTT,EC1_TTS,EC1_TSS,EC1_SSS
       REAL(8)            :: SIGMA_T,SIGMA_S,SIGMA_TT,SIGMA_TS
       REAL(8)            :: SIGMA_TTT,SIGMA_TTS
       LOGICAL(4)         :: TTABLE
!      *****************************************************************
       IF(.NOT.TINI) CALL PERDEWWANG91L_INITIALIZE
!
!      =================================================================
!      == DECIDE ABOUT TABLE LOOKUP OR NOT                            ==
!      =================================================================
       TTABLE=.NOT.TSAFE
       IF(TTABLE) THEN
         CALL TABLE1D$SETX(RHOTABLE,RHOT,TTABLE)
       END IF
!
!      =================================================================
!      == GET RHO-DEPENCE FOR THE NON-SPIN POLARIZED ELECTRON GAS     ==
!      =================================================================
       IF(TTABLE) THEN
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECUARRAY,ECU)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECUARRAY,DECU)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D2ECUARRAY,D2ECU)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D3ECUARRAY,D3ECU)
       ELSE
         CALL PERDEWWANG91L_RHODEP3(RHOT,ECU,DECU,D2ECU,D3ECU &
      &                         ,AU,ALPHA1U,BETA1U,BETA2U,BETA3U,BETA4U)
       END IF
!
!      =================================================================
!      == RETURN EARLY FOR NON-SPIN POLARIZED CALCULATIONS            ==
!      =================================================================
       IF(.NOT.TSPIN) THEN
         EC=ECU
         EC_T=DECU     
         EC_S=0.D0
         EC_TT=D2ECU
         EC_TS=0.D0
         EC_SS=0.D0
         EC_TTT=D3ECU
         EC_TTS=0.D0
         EC_TSS=0.D0
         EC_SSS=0.D0
         RETURN
       END IF
!
!      =================================================================
!      == CALCULATE EPS_XC FOR FULLY POLARIZED ELECTRON GAS           ==
!      == AND THE SPIN STIFFNESS                                      ==
!      =================================================================
       IF(TTABLE) THEN
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECPARRAY,ECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECPARRAY,DECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D2ECPARRAY,D2ECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D3ECPARRAY,D3ECP)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,ECXARRAY,ECX)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,DECXARRAY,DECX)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D2ECXARRAY,D2ECX)
         CALL TABLE1D$GETF(RHOTABLE,NRHO,D3ECXARRAY,D3ECX)
       ELSE
         CALL PERDEWWANG91L_RHODEP3(RHOT,ECP,DECP,D2ECP,D3ECP &
      &                         ,AP,ALPHA1P,BETA1P,BETA2P,BETA3P,BETA4P)
!      == NOTE THAT ECX IS -ALPHA_C FOMR THE PAPER !!
         CALL PERDEWWANG91L_RHODEP3(RHOT,ECX,DECX,D2ECX,D3ECX &
      &                         ,AX,ALPHA1X,BETA1X,BETA2X,BETA3X,BETA4X)
         ECX=-ECX/FSIGMAPP
         DECX=-DECX/FSIGMAPP
         D2ECX=-D2ECX/FSIGMAPP
         D3ECX=-D3ECX/FSIGMAPP
       END IF
!
!      =================================================================
!      == CALCULATE SPIN DEPENDENCE                                   ==
!      =================================================================
       SIGMA=RHOS/RHOT
       SIGMAP=1.D0+SIGMA
       SIGMAP13=SIGMAP**ONEBY3
       SIGMAP23=SIGMAP13*SIGMAP13
       SIGMAP43=SIGMAP23*SIGMAP23
       SIGMAP53=SIGMAP23*SIGMAP
       SIGMAM=1.D0-SIGMA
       SIGMAM13=SIGMAM**ONEBY3
       SIGMAM23=SIGMAM13*SIGMAM13
       SIGMAM43=SIGMAM23*SIGMAM23
       SIGMAM53=SIGMAM23*SIGMAM
       FSIGMA  = FFAC*(SIGMAP43  + SIGMAM43 - 2.D0)
       DFSIGMA = DFFAC*(SIGMAP13  - SIGMAM13)
       D2FSIGMA= D2FFAC*(1.D0/SIGMAP23+ 1.D0/SIGMAM23) 
       D3FSIGMA= D3FFAC*(1.D0/SIGMAP53- 1.D0/SIGMAM53) 
!
!      =================================================================
!      == FORM DERIVATIVES FIRST WITH RESPECT TO RHOT AND SIGMA       ==
!      =================================================================
       SIGMA2=SIGMA*SIGMA
       SIGMA3=SIGMA2*SIGMA
       SIGMA4=SIGMA3*SIGMA
!
       SVAR1    =(ECP-ECU)  *SIGMA4+ECX*(1.D0-SIGMA4)
       SVAR1_T  =(DECP-DECU)*SIGMA4+DECX*(1.D0-SIGMA4)
       SVAR1_S  =(ECP-ECU-ECX)*4.D0*SIGMA3
       SVAR1_SS =(ECP-ECU-ECX)*12.D0*SIGMA2
       SVAR1_TS =(DECP-DECU-DECX)*4.D0*SIGMA3
       SVAR1_TT =(D2ECP-D2ECU)*SIGMA4+D2ECX*(1.D0-SIGMA4)
       SVAR1_TTT=(D3ECP-D3ECU)*SIGMA4+D3ECX*(1.D0-SIGMA4)
       SVAR1_TTS=(D2ECP-D2ECU-D2ECX)*4.D0*SIGMA3
       SVAR1_TSS=(DECP-DECU-DECX)*12.D0*SIGMA2
       SVAR1_SSS =(ECP-ECU-ECX)*24.D0*SIGMA
!
       EC1    =ECU+FSIGMA*SVAR1
       EC1_T  =DECU+FSIGMA*SVAR1_T
       EC1_S  =DFSIGMA*SVAR1 + FSIGMA*SVAR1_S
       EC1_TT =D2ECU+FSIGMA*SVAR1_TT
       EC1_TS =DFSIGMA*SVAR1_T+FSIGMA*SVAR1_TS
       EC1_SS =D2FSIGMA*SVAR1 + 2.D0*DFSIGMA*SVAR1_S+FSIGMA*SVAR1_SS
       EC1_TTT=D3ECU+FSIGMA*SVAR1_TTT
       EC1_TTS=DFSIGMA*SVAR1_TT+FSIGMA*SVAR1_TTS
       EC1_TSS=D2FSIGMA*SVAR1_T + 2.D0*DFSIGMA*SVAR1_TS+FSIGMA*SVAR1_TSS
       EC1_SSS=D3FSIGMA*SVAR1 +3.D0*D2FSIGMA*SVAR1_S+3.D0*DFSIGMA*SVAR1_SS+FSIGMA*SVAR1_SSS
!
!      =================================================================
!      == TRANSFORM DERIVATIVES TO RHOT AND RHOS                      ==
!      =================================================================
       SIGMA_S =1.D0/RHOT
       SIGMA_T =-SIGMA*SIGMA_S
       SIGMA_TS=-SIGMA_S**2
       SIGMA_TT=-2.D0*SIGMA_T*SIGMA_S
       SIGMA_TTT=-2.D0*(SIGMA_TT*SIGMA_S+SIGMA_T*SIGMA_TS)
       SIGMA_TTS=-2.D0*SIGMA_S*SIGMA_TS
!
       EC   =EC1
       EC_T =EC1_T+EC1_S*SIGMA_T
       EC_S =EC1_S*SIGMA_S
       EC_TT=EC1_TT+2.D0*EC1_TS*SIGMA_T+EC1_SS*SIGMA_T**2+EC1_S*SIGMA_TT
       EC_TS=EC1_TS*SIGMA_S+EC1_SS*SIGMA_S*SIGMA_T+EC1_S*SIGMA_TS
       EC_SS=EC1_SS*SIGMA_S**2 
       EC_TTT=EC1_TTT+3.D0*EC1_TTS*SIGMA_T+3.D0*EC1_TSS*SIGMA_T**2+EC1_SSS*SIGMA_T**3 &
      &      +3.D0*EC1_TS*SIGMA_TT+3.D0*EC1_SS*SIGMA_T*SIGMA_TT &
      &      +EC1_S*SIGMA_TTT
       EC_TTS= EC1_TTS*SIGMA_S+2.D0*EC1_TSS*SIGMA_S*SIGMA_T+EC1_SSS*SIGMA_S*SIGMA_T**2 &
      &      +2.D0*EC1_TS*SIGMA_TS+2.D0*EC1_SS*SIGMA_TS*SIGMA_T+EC1_SS*SIGMA_S*SIGMA_TT &
      &      +EC1_S*SIGMA_TTS
       EC_TSS=EC1_TSS*SIGMA_S**2+EC1_SSS*SIGMA_S**2*SIGMA_T+2.D0*EC1_SS*SIGMA_S*SIGMA_TS
       EC_SSS=EC1_SSS*SIGMA_S**3
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L$EVAL1(VAL,EXC,DEXC)
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: VAL(5)
       REAL(8),INTENT(OUT):: EXC
       REAL(8),INTENT(OUT):: DEXC(5)
       REAL(8)            :: RHOT
       REAL(8)            :: RHOS
       REAL(8)            :: EC,EC_T,EC_S
!      *****************************************************************
       DEXC(:)=0.D0
       RHOT=VAL(1)
       RHOS=VAL(2)
       CALL PERDEWWANG91L$EPSVAL(RHOT,RHOS,EC,EC_T,EC_S)
       EXC=RHOT*EC
       DEXC(1)=RHOT*EC_T+EC
       DEXC(2)=RHOT*EC_S
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L$EVAL2(VAL,EXC,DEXC,D2EXC)
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: VAL(5)
       REAL(8),INTENT(OUT):: EXC
       REAL(8),INTENT(OUT):: DEXC(5)
       REAL(8),INTENT(OUT):: D2EXC(5,5)
       REAL(8)            :: RHOT
       REAL(8)            :: RHOS
       REAL(8)            :: EC,EC_T,EC_S
       REAL(8)            :: EC_TT,EC_TS,EC_SS
!      *****************************************************************
       DEXC(:)=0.D0
       D2EXC(:,:)=0.D0
       RHOT=VAL(1)
       RHOS=VAL(2)
       CALL PERDEWWANG91L$EPSVAL2(RHOT,RHOS,EC,EC_T,EC_S,EC_TT,EC_TS,EC_SS)
       EXC=RHOT*EC
       DEXC(1)=RHOT*EC_T+EC
       DEXC(2)=RHOT*EC_S
       D2EXC(1,1)=2.D0*EC_T+RHOT*EC_TT
       D2EXC(2,1)=EC_S+RHOT*EC_TS
       D2EXC(1,2)=D2EXC(2,1)
       D2EXC(2,2)=RHOT*EC_SS
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PERDEWWANG91L$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: VAL(5)
       REAL(8),INTENT(OUT):: EXC
       REAL(8),INTENT(OUT):: DEXC(5)
       REAL(8),INTENT(OUT):: D2EXC(5,5)
       REAL(8),INTENT(OUT):: D3EXC(5,5,5)
       REAL(8)            :: RHOT
       REAL(8)            :: RHOS
       REAL(8)            :: EC,EC_T,EC_S
       REAL(8)            :: EC_TT,EC_TS,EC_SS
       REAL(8)            :: EC_TTT,EC_TTS,EC_TSS,EC_SSS
!      *****************************************************************
       DEXC(:)=0.D0
       D2EXC(:,:)=0.D0
       D3EXC(:,:,:)=0.D0
       RHOT=VAL(1)
       RHOS=VAL(2)
       CALL PERDEWWANG91L$EPSVAL3(RHOT,RHOS,EC,EC_T,EC_S &
      &                  ,EC_TT,EC_TS,EC_SS &
      &                  ,EC_TTT,EC_TTS,EC_TSS,EC_SSS)
       EXC=RHOT*EC
       DEXC(1)=RHOT*EC_T+EC
       DEXC(2)=RHOT*EC_S
       D2EXC(1,1)=2.D0*EC_T+RHOT*EC_TT
       D2EXC(2,1)=EC_S+RHOT*EC_TS
       D2EXC(2,2)=RHOT*EC_SS
       D3EXC(1,1,1)=3.D0*EC_TT+RHOT*EC_TTT
       D3EXC(1,1,2)=2.D0*EC_TS+RHOT*EC_TTS
       D3EXC(1,2,2)=EC_SS+RHOT*EC_TSS
       D3EXC(2,2,2)=RHOT*EC_SSS
       D2EXC(1,2)  =D2EXC(2,1)
       D3EXC(1,2,1)=D3EXC(1,1,2)
       D3EXC(2,1,1)=D3EXC(1,1,2)
       D3EXC(2,1,2)=D3EXC(1,2,2)
       D3EXC(2,2,1)=D3EXC(1,2,2)
       RETURN
       END
!
!      .................................................................
MODULE PERDEWWANG91G_MODULE
LOGICAL(4)         :: TINI=.FALSE.
REAL(8)            :: PI
REAL(8)            :: RSFAC   ! RS=RSFAC*RHO**(-1/3)
REAL(8)            :: KFFAC   ! KF= KFFAC/RS
REAL(8)            :: KSFAC   ! KS=KSFAC/SQRT(RS)
REAL(8)            :: TFAC    ! T2=TFAC*GRHO2/G**2*RHOT**(-7/3)
REAL(8),PARAMETER  :: CC0=4.235D-3
REAL(8)            :: NU      
REAL(8),PARAMETER  :: ALPHA=9.D-2
REAL(8)            :: BETA    ! BETA=NU*CC0 (PWG:EQ.13)
REAL(8),PARAMETER  :: CX=-1.667212D-3
CONTAINS
!
!      .................................................................
       SUBROUTINE PERDEWWANG91G_INITIALIZE
       IMPLICIT NONE 
       PI=4.D0*ATAN(1.D0)
       RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0) ! RS=RSFAC*RHOT**(-1/3)
       KFFAC=(2.25D0*PI)**(1.D0/3.D0)      ! KF=KFFAC/RS
       KSFAC=SQRT(4.D0*KFFAC/PI)          ! KS=KSFAC/SQRT(RS)
       TFAC=0.25D0*RSFAC/KSFAC**2          ! T2=TFAC*GRHO2/G**2*RHOT**(-7/3)
       NU=(16.D0/PI)*(3.D0*PI**2)**(1.D0/3.D0)
       BETA=NU*CC0
       TINI=.TRUE.
       RETURN
       END SUBROUTINE PERDEWWANG91G_INITIALIZE
!
!      .................................................................
       SUBROUTINE H1OFTRSG2(T2,RS,G,H,H_T,H_R,H_G &
      &                    ,H_TT,H_TR,H_TG,H_RR,H_RG,H_GG)
!      *****************************************************************    
!      **  EQ.15 OF PRB46,6671(1992-II)                               **    
!      **                                                             **    
!      *****************************************************************    
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: T2     ! T**2
       REAL(8),INTENT(IN) :: RS     
       REAL(8),INTENT(IN) :: G
       REAL(8),INTENT(OUT):: H
       REAL(8),INTENT(OUT):: H_T,H_R,H_G
       REAL(8),INTENT(OUT):: H_TT,H_TR,H_TG,H_RR,H_RG,H_GG
       REAL(8),PARAMETER  :: THREEBY7=3.D0/7.D0
       REAL(8)            :: CC,DCC,D2CC
       REAL(8)            :: CEXPO
       REAL(8)            :: G2,G3,G4
       REAL(8)            :: EXPO,EXPO_T,EXPO_R
       REAL(8)            :: SVAR1,SVAR3,SVAR4
!
!      =================================================================
!      ==  CONSTANTS                                                  ==
!      =================================================================
       CEXPO=-100.D0*(KSFAC/KFFAC)**2
!
!      =================================================================
!      ==  CHECK WHETHER CUTOFF VANISHES         S                    ==
!      =================================================================
       G2=G*G
       G3=G2*G
       G4=G3*G
       EXPO=CEXPO*T2*G4*RS
       IF(EXPO.LT.-100.D0) THEN
         H=0.D0
         H_T=0.D0
         H_R=0.D0
         H_G=0.D0
         H_TT=0.D0
         H_TR=0.D0
         H_TG=0.D0
         H_RR=0.D0
         H_RG=0.D0
         H_GG=0.D0
         RETURN
       END IF
!
!      =================================================================
!      ==  CC DEPENDENCE                                              ==
!      =================================================================
       CALL GETCC(RS,CC,DCC,D2CC)
!
!      =================================================================
!      ==  H1                                                         ==
!      =================================================================
       SVAR1=CC-CC0-THREEBY7*CX
       EXPO_T=CEXPO*G4*RS
       EXPO_R=CEXPO*G4*T2
       SVAR3=NU*SVAR1*EXP(EXPO)
       SVAR4=DCC/SVAR1+EXPO_R
       H   =SVAR3*G3*T2
       H_T =SVAR3*G3*(1.D0+EXPO)
       H_R =H*SVAR4
       H_G =SVAR3*G2*T2*(3.D0+4.D0*EXPO)
       H_TT=SVAR3*G3*(2.D0+EXPO)*EXPO_T
       H_TR=SVAR3*G3*((1.D0+EXPO)*SVAR4+EXPO_R)
       H_RR=H*((D2CC+2.D0*DCC*EXPO_R)/SVAR1+EXPO_R**2)
       H_TG=SVAR3*G2*((3.D0+4.D0*EXPO)*(1.D0+EXPO)+4.D0*EXPO)
       H_RG=SVAR3*G2*T2*(SVAR4*(3.D0+4.D0*EXPO)+4.D0*EXPO_R)
       H_GG=SVAR3*G*T2*(6.D0+(36.D0+16.D0*EXPO)*EXPO)
       RETURN
       END SUBROUTINE H1OFTRSG2
!
!      .................................................................
       SUBROUTINE H0OFTEG2(T2,E,G,H,H_T,H_E,H_G,H_TT,H_TE,H_TG,H_EE,H_EG,H_GG)
!      *****************************************************************    
!      **                                                             **    
!      **                                                             **    
!      **  THE PARAMATERS MUST LIE IN THE FOLLOWING RANGE:            **    
!      **    T2>0; E<0 G>0                                            **    
!      **                                                             **    
!      *****************************************************************    
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: T2
       REAL(8),INTENT(IN) :: E
       REAL(8),INTENT(IN) :: G
       REAL(8),INTENT(OUT):: H
       REAL(8),INTENT(OUT):: H_T,H_E,H_G
       REAL(8),INTENT(OUT):: H_TT,H_TE,H_TG,H_EE,H_EG,H_GG
       REAL(8)            :: FAC1,FAC2
       REAL(8)            :: V,V_E,V_G,V_EG,V_GG
       REAL(8)            :: AI,AI_V,AI_VV
       REAL(8)            :: SVAR,SVAR2,SVAR3,DSVAR,D2SVAR
       REAL(8)            :: T4,T6
       REAL(8)            :: G2,G3
       REAL(8)            :: AI2,AI3
       REAL(8)            :: SVAR1,SVAR1_A,SVAR1_T,SVAR1_AA,SVAR1_AT,SVAR1_TT
       REAL(8)            :: SVAR1_V,SVAR1_TV,SVAR1_VV
       REAL(8)            :: H_V,H_TV,H_VV,H_VG
!      *****************************************************************    
!
!      =================================================================
!      ==  CONSTANTS                                                  ==
!      =================================================================
       FAC1=2.D0*ALPHA/BETA    
       FAC2=FAC1/BETA   !(2*ALPHA/BETA**2)
!
!      =================================================================
!      ==  CALCULATE A(VAR1=E/G3)                                     ==
!      ==  G>0; EPSC=<0 ; 0=< AINV < INFTY                            ==
!      ==  AINV BECOMES LARGE FOR LARGE DENSITIES                     ==
!      =================================================================
       G2=G*G
       G3=G2*G
       V=E/G3
!      == CHECK WHETHER EXPONENTIAL OVERFLOWS
       IF(-FAC2*V.GT.100.D0) THEN
         SVAR=LOG(1.D0+FAC1*T2)/FAC2
         DSVAR=FAC1/(1.D0+FAC1*T2)/FAC2
         D2SVAR=-FAC2*DSVAR**2
         H=G3*SVAR
         H_T=G3*DSVAR
         H_E=0.D0
         H_G=3.D0*G2*SVAR
         H_TT=G3*D2SVAR
         H_TE=0.D0
         H_TG=3.D0*G2*DSVAR
         H_EE=0.D0
         H_EG=0.D0
         H_GG=6.D0*G*SVAR
PRINT*,'SHORTCUT FOR H0'
         RETURN
       END IF
       SVAR=EXP(-FAC2*V)
       AI  =(SVAR-1.D0)/FAC1
       AI_V =-FAC2/FAC1*SVAR
       AI_VV=-FAC2*AI_V
!
!      =================================================================
!      == 0 .LE. SVAR .LT.< 1.D0                                      ==       
!      =================================================================
       AI2=AI*AI
       AI3=AI2*AI
       T4=T2*T2
       T6=T4*T2
       SVAR=1.D0/(AI2+AI*T2+T4)
       SVAR2=SVAR*SVAR
       SVAR3=6.D0*AI*T2*(AI+T2)*SVAR2*SVAR
       SVAR =SVAR*FAC1
       SVAR2=SVAR2*FAC1
       SVAR3=SVAR3*FAC1
!
       SVAR1   =AI*T2*(AI+T2)*SVAR
       SVAR1_T =AI3*(2.D0*T2+AI)*SVAR2
       SVAR1_A =T6*(2.D0*AI+T2)*SVAR2
       SVAR1_TT=-SVAR3*AI2
       SVAR1_AT=+SVAR3*AI*T2
       SVAR1_AA=-SVAR3*T4

!      ==  TRANSFORM SVAR1(T2,AI) TO SVAR1(T2,V=VAR1)  ====================
       SVAR1_V   =SVAR1_A*AI_V
       SVAR1_VV  =SVAR1_AA*AI_V**2+SVAR1_A*AI_VV
       SVAR1_TV  =SVAR1_AT*AI_V
!
!      =================================================================
!      ==  CALCULATE H(T,V,G)                                        ==
!      =================================================================
       SVAR=LOG(1.D0+SVAR1)/FAC2
       DSVAR=1.D0/(1.D0+SVAR1)/FAC2
       D2SVAR=-FAC2*DSVAR*DSVAR
       H   =      G3*SVAR
       H_G = 3.D0*G2*SVAR
       H_T =      G3*DSVAR*SVAR1_T
       H_V =      G3*DSVAR*SVAR1_V
       H_TT=      G3*(D2SVAR*SVAR1_T**2+DSVAR*SVAR1_TT)
       H_TV=      G3*(D2SVAR*SVAR1_T*SVAR1_V+DSVAR*SVAR1_TV)
       H_TG= 3.D0*G2*DSVAR*SVAR1_T
       H_VV=      G3*(D2SVAR*SVAR1_V**2+DSVAR*SVAR1_VV)
       H_VG= 3.D0*G2*DSVAR*SVAR1_V
       H_GG= 6.D0*G *SVAR
!
!      =================================================================
!      ==  TRANSFORM FROM H(T2,V(E,G),G) TO H(T2,E,G)                 ==
!      =================================================================
       V_E=1.D0/G3
       V_G=-3.D0*V/G
       V_EG=-3.D0*V_E/G
       V_GG=-4.D0*V_G/G
!
       H_TE=H_TV*V_E
       H_TG=H_TG+H_TV*V_G
       H_EE=H_VV*V_E**2
       H_EG=H_VG*V_E+H_VV*V_G*V_E+H_V*V_EG
       H_GG=H_GG+2.D0*H_VG*V_G + H_VV*V_G**2+H_V*V_GG
       H_G =H_G+H_V*V_G
       H_E =H_V*V_E
       RETURN
       END SUBROUTINE H0OFTEG2
!       ................................................................
        SUBROUTINE GETCC(RS,CXC,DCXC,D2CXC)
!       ****************************************************************
!       **  C_C(RS)=C_XC(RS)-C_X                                      **
!       **  C_XC(RS) IS GIVEN IN EQ. 4D OF                            **
!       **  M. RASOLT AND D.J.W. GELDART, PHYS.REV.B 34, 1325 (1986)  **
!       **  AND C_X IS THE SHAM COEFFICIENT C_X=-1.667212E-3          **
!       ****************************************************************
        IMPLICIT NONE
        REAL(8),INTENT(IN)  :: RS
        REAL(8),INTENT(OUT) :: CXC
        REAL(8),INTENT(OUT) :: DCXC
        REAL(8),INTENT(OUT) :: D2CXC
        REAL(8),PARAMETER   :: A=23.266D0
        REAL(8),PARAMETER   :: B=7.389D-3
        REAL(8),PARAMETER   :: C=8.723D0
        REAL(8),PARAMETER   :: D=0.472D0
        REAL(8),PARAMETER   :: E=2.568D0
        REAL(8)             :: RS2,RS3
        REAL(8)             :: U,DU,D2U
        REAL(8)             :: V,DV,D2V
        REAL(8)             :: ONEBYV,ONEBYV2,ONEBYV3
!MODULE REAL(8),PARAMETER   :: CX=-1.667D-3
!       ****************************************************************
        RS2=RS*RS
        RS3=RS2*RS
        U  =E+A*RS+     B*RS2
        DU =  A   +2.D0*B*RS
        D2U=       2.D0*B
        V  =1.D0+C*RS+     D*RS2+10.D0*B*RS3
        DV =     C   +2.D0*D*RS +30.D0*B*RS2
        D2V=          2.D0*D    +60.D0*B*RS
        ONEBYV=1.D0/V
        ONEBYV2=ONEBYV*ONEBYV
        ONEBYV3=ONEBYV*ONEBYV2
        CXC=1.D-3*U*ONEBYV-CX
        DCXC=1.D-3*(DU*ONEBYV-U*DV*ONEBYV2)
        D2CXC=1.D-3*(D2U*ONEBYV-(2.D0*DU*DV+U*D2V)*ONEBYV2 &
      &              +2.D0*U*DV**2*ONEBYV3)
        END SUBROUTINE GETCC
       END MODULE PERDEWWANG91G_MODULE

!
!      .................................................................
       SUBROUTINE PERDEWWANG91G$EVAL(RHOT,RHOS,GRHO2,EPSC,EPSC_D,EPSC_S,EPSC_N)
       USE PERDEWWANG91G_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: RHOS
       REAL(8),INTENT(IN) :: GRHO2
       REAL(8),INTENT(OUT):: EPSC
       REAL(8),INTENT(OUT):: EPSC_D,EPSC_S,EPSC_N
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
       REAL(8),PARAMETER  :: ONEBY9=1.D0/9.D0
       REAL(8)            :: RHOT13
       REAL(8)            :: P13,M13,P23,M23
       REAL(8)            :: SIG,SIG_D,SIG_S
       REAL(8)            :: G,G_SIG
       REAL(8)            :: G_S,G_D
       REAL(8)            :: E,E_D,E_S,E_DD,E_DS,E_SS
       REAL(8)            :: RS,R_D
       REAL(8)            :: T,T_N,T_G,T_D,T_S
       REAL(8)            :: H0,H0_T,H0_E,H0_G
       REAL(8)            :: H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG
       REAL(8)            :: H1,H1_T,H1_R,H1_G 
       REAL(8)            :: H1_TT,H1_TR,H1_TG,H1_RR,H1_RG,H1_GG
       REAL(8)            :: H0_D,H0_S,H0_N
       REAL(8)            :: H1_D,H1_S,H1_N
!      *****************************************************************
       IF(RHOT.LT.1.D-10) THEN
         EPSC=0.D0
         EPSC_D=0.D0
         EPSC_S=0.D0
         EPSC_N=0.D0
         RETURN
       END IF
       IF(.NOT.TINI) CALL PERDEWWANG91G_INITIALIZE
!
!      ================================================================
!      == EQ.11: SPIN DEPENDENCE G                                   ==
!      ================================================================
       SIG=RHOS/RHOT
       SIG=MIN(1.D0-1.D-6,SIG)
       SIG=MAX(-1.D0+1.D-6,SIG)
       P13=(1.D0+SIG)**ONEBY3
       M13=(1.D0-SIG)**ONEBY3
       P23=P13*P13
       M23=M13*M13
       G=0.5D0*(P23+M23)
       G_SIG=ONEBY3*(1.D0/P13-1.D0/M13)
!
!      ================================================================
!      == CALCULATE LOCAL CORRELATION                                ==
!      ================================================================
       CALL PERDEWWANG91L$EPSVAL2(RHOT,RHOS,E,E_D,E_S,E_DD,E_DS,E_SS)
!
!      ================================================================
!      == T IS T**2 OF EQ. 10                                       ==
!      ================================================================
       RHOT13=RHOT**ONEBY3
       RS=RSFAC/RHOT13
       T_N=TFAC/G**2/RHOT13**7
       T  =T_N*GRHO2
       T_G=-2.D0*T/G
       T_D=-7.D0/3.D0*T/RHOT
       IF(T.LT.1.D-8) THEN
         EPSC   =RHOT*E
         EPSC_D =RHOT*E_D+E
         EPSC_S =RHOT*E_S
         EPSC_N =0.D0
         RETURN
       END IF
!
!      ================================================================
!      == CALCULATE H0                                               ==
!      ================================================================
!      == CONSIDER H0 AS FUNCTION OF T,EPSC,G 
       CALL H0OFTEG2(T,E,G,H0,H0_T,H0_E,H0_G &
      &            ,H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG)
!
!      ================================================================
!      == CALCULATE H0                                               ==
!      ================================================================
       CALL H1OFTRSG2(T,RS,G,H1,H1_T,H1_R,H1_G &
      &             ,H1_TT,H1_TR,H1_TG,H1_RR,H1_RG,H1_GG)
!
!      ================================================================
!      == CONVERT FROM RS TO RHOT                                    ==
!      ================================================================
!      == SIG(D=RHOT,S=RHOS)
       SIG_D=-SIG/RHOT
       SIG_S=1.D0/RHOT
!      == G(D=RHOT,S=RHOS)
       G_D=G_SIG*SIG_D
       G_S=G_SIG*SIG_S
!      == RS(RHOT) ====================================================
       R_D=-ONEBY3*RS/RHOT
!      == T(D=RHOT,S=RHOS,N=GRHO2) ====================================
       T_S=T_G*G_S
       T_D=T_D+T_G*G_D

!      == TRANSFORM H0 ================================================        
       H0_D  = H0_T*T_D+H0_E*E_D+H0_G*G_D
       H0_S  = H0_T*T_S+H0_E*E_S+H0_G*G_S
       H0_N  = H0_T*T_N
!      == TRANSFORM H1
       H1_D  = H1_T*T_D+H1_R*R_D+H1_G*G_D
       H1_S  = H1_T*T_S+H1_G*G_S
       H1_N  = H1_T*T_N
!
!      == EXCHANGE ENERGY ==============================================
       EPSC   =RHOT*(E+H0+H1)
       EPSC_D =RHOT*(E_D+H0_D+H1_D)+E+H0+H1
       EPSC_S =RHOT*(E_S+H0_S+H1_S)
       EPSC_N =RHOT*(H0_N+H1_N)
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PERDEWWANG91G$EVAL2(RHOT,RHOS,GRHO2,EPSC &
      &          ,EPSC_D,EPSC_S,EPSC_N &
      &          ,EPSC_DD,EPSC_DS,EPSC_DN,EPSC_SS,EPSC_SN,EPSC_NN)
       USE PERDEWWANG91G_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: RHOT
       REAL(8),INTENT(IN) :: RHOS
       REAL(8),INTENT(IN) :: GRHO2
       REAL(8),INTENT(OUT):: EPSC
       REAL(8),INTENT(OUT):: EPSC_D,EPSC_S,EPSC_N
       REAL(8),INTENT(OUT):: EPSC_DD,EPSC_DS,EPSC_DN,EPSC_SS,EPSC_SN,EPSC_NN
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8),PARAMETER  :: ONEBY9=1.D0/9.D0
       REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
       REAL(8),PARAMETER  :: SEVENBY6=7.D0/6.D0
       REAL(8)            :: RHOT13
       REAL(8)            :: P13,M13,P23,M23,P43,M43
       REAL(8)            :: SIG,SIG_D,SIG_S,SIG_DD,SIG_DS
       REAL(8)            :: G,G_SIG,G_SIGSIG
       REAL(8)            :: G_S,G_D,G_DD,G_DS,G_SS
       REAL(8)            :: E,E_D,E_S,E_DD,E_DS,E_SS
       REAL(8)            :: RS,R_D,R_DD
       REAL(8)            :: T,T_N,T_G,T_D,T_DD,T_DN,T_DG,T_GG,T_GN,T_NN
       REAL(8)            :: T_S,T_DS,T_SS,T_SN
       REAL(8)            :: H0,H0_T,H0_E,H0_G
       REAL(8)            :: H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG
       REAL(8)            :: H1,H1_T,H1_R,H1_G 
       REAL(8)            :: H1_TT,H1_TR,H1_TG,H1_RR,H1_RG,H1_GG
       REAL(8)            :: H0_D,H0_S,H0_N
       REAL(8)            :: H0_DD,H0_DS,H0_DN,H0_SS,H0_SN,H0_NN
       REAL(8)            :: H1_D,H1_S,H1_N
       REAL(8)            :: H1_DD,H1_DS,H1_DN,H1_SS,H1_SN,H1_NN
!      *****************************************************************
       IF(RHOT.LT.1.D-10) THEN
         EPSC=0.D0
         EPSC_D=0.D0
         EPSC_S=0.D0
         EPSC_N=0.D0
         EPSC_DD=0.D0
         EPSC_DS=0.D0
         EPSC_DN=0.D0
         EPSC_SS=0.D0
         EPSC_SN=0.D0
         EPSC_NN=0.D0
         RETURN
       END IF
       IF(.NOT.TINI) CALL PERDEWWANG91G_INITIALIZE
!
!      ================================================================
!      == EQ.11: SPIN DEPENDENCE G                                   ==
!      ================================================================
       SIG=RHOS/RHOT
       SIG=MIN(1.D0-1.D-6,SIG)
       SIG=MAX(-1.D0+1.D-6,SIG)
       P13=(1.D0+SIG)**ONEBY3
       M13=(1.D0-SIG)**ONEBY3
       P23=P13*P13
       M23=M13*M13
       P43=P23*P23
       M43=M23*M23
       G=0.5D0*(P23+M23)
       G_SIG=ONEBY3*(1.D0/P13-1.D0/M13)
       G_SIGSIG=-ONEBY9*(1.D0/P43+1.D0/M43)
!
!      ================================================================
!      == CALCULATE LOCAL CORRELATION                                ==
!      ================================================================
       CALL PERDEWWANG91L$EPSVAL2(RHOT,RHOS,E,E_D,E_S,E_DD,E_DS,E_SS)
!
!      ================================================================
!      == T IS T**2 OF EQ. 10                                       ==
!      ================================================================
       RHOT13=RHOT**ONEBY3
       RS=RSFAC/RHOT13
       T_N=TFAC/G**2/RHOT13**7
       T=T_N*GRHO2
       T_NN=0.D0
       T_G=-2.D0*T/G
       T_D=-7.D0/3.D0*T/RHOT
       T_DD=-10.D0/3.D0*T_D/RHOT
       T_DN=-7.D0/3.D0*T_N/RHOT
       T_DG=-2.D0*T_D/G
       T_GG=-3.D0*T_G/G
       T_GN=-2.D0*T_N/G
!
!      ================================================================
!      == SHORTCUT FOR SMALL T                                       ==
!      ================================================================
       IF(T.LT.1.D-8) THEN
         EPSC   =RHOT*E
         EPSC_D =RHOT*E_D+E
         EPSC_S =RHOT*E_S
         EPSC_N =0.D0
         EPSC_DD=RHOT*E_DD+2.D0*E_D
         EPSC_DS=RHOT*E_DS+E_S
         EPSC_DN=0.D0
         EPSC_SS=RHOT*E_SS
         EPSC_SN=0.D0
         EPSC_NN=0.D0
         PRINT*,'SHORTCUT WITH SMALL T'
         RETURN
       END IF
!
!      ================================================================
!      == CALCULATE H0                                               ==
!      ================================================================
!      == CONSIDER H0 AS FUNCTION OF T,EPSC,G 
       CALL H0OFTEG2(T,E,G,H0,H0_T,H0_E,H0_G &
      &            ,H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG)
!
!      ================================================================
!      == CALCULATE H0                                               ==
!      ================================================================
       CALL H1OFTRSG2(T,RS,G,H1,H1_T,H1_R,H1_G &
      &             ,H1_TT,H1_TR,H1_TG,H1_RR,H1_RG,H1_GG)
!
!      ================================================================
!      == CONVERT FROM RS TO RHOT                                    ==
!      ================================================================
!      == SIG(D=RHOT,S=RHOS)
       SIG_D=-SIG/RHOT
       SIG_S=1.D0/RHOT
       SIG_DD=-2.D0*SIG_D/RHOT
       SIG_DS=-SIG_S/RHOT
!      == G(D=RHOT,S=RHOS)
       G_D=G_SIG*SIG_D
       G_S=G_SIG*SIG_S
       G_DD=G_SIGSIG*SIG_D**2+G_SIG*SIG_DD
       G_DS=G_SIGSIG*SIG_D*SIG_S+G_SIG*SIG_DS
       G_SS=G_SIGSIG*SIG_S**2
!      == RS(RHOT) ====================================================
       R_D=-ONEBY3*RS/RHOT
       R_DD=-FOURBY3*R_D/RHOT
!      == T(D=RHOT,S=RHOS,N=GRHO2) ====================================
!      == WATCH ORDER OF STATEMENTS!! =================================
!      T=SQRT(GRHO2)/(2.D0*G*KSFAC/SQRT(RSFAC))*RHOT**(-SEVENBY6)
       T_DD=T_DD+2.D0*T_DG*G_D+T_GG*G_D**2+T_G*G_DD
       T_DS=T_DG*G_S+T_GG*G_S*G_D+T_G*G_DS
       T_DN=T_DN+T_GN*G_D
       T_SS=T_GG*G_S**2+T_G*G_SS
       T_SN=T_GN*G_S
       T_S=T_G*G_S
       T_D=T_D+T_G*G_D
!      == TRANSFORM H0 ================================================        
       H0_DD = H0_TT*T_D**2 + 2.D0*H0_TE*T_D*E_D + 2.D0*H0_TG*T_D*G_D +H0_T*T_DD &
      &      + H0_EE*E_D**2 + 2.D0*H0_EG*E_D*G_D +H0_E*E_DD &
      &      + H0_GG*G_D**2 + H0_G*G_DD
       H0_SS = H0_TT*T_S**2 + 2.D0*H0_TE*T_S*E_S + 2.D0*H0_TG*T_S*G_S +H0_T*T_SS &
      &      + H0_EE*E_S**2 + 2.D0*H0_EG*E_S*G_S +H0_E*E_SS &
      &      + H0_GG*G_S**2 + H0_G*G_SS
       H0_NN = H0_TT*T_N**2  + H0_T*T_NN
       H0_DS = H0_TT*T_D*T_S + H0_TE*(T_D*E_S+T_S*E_D) + H0_TG*(T_D*G_S+T_S*G_D) +H0_T*T_DS &
      &      + H0_EE*E_D*E_S + H0_EG*(E_D*G_S+E_S*G_D) + H0_E*E_DS &
      &      + H0_GG*G_D*G_S + H0_G*G_DS
       H0_DN = H0_TT*T_D*T_N + H0_TE*T_N*E_D + H0_TG*T_N*G_D +H0_T*T_DN
       H0_SN = H0_TT*T_S*T_N + H0_TE*T_N*E_S + H0_TG*T_N*G_S +H0_T*T_SN
       H0_D  = H0_T*T_D+H0_E*E_D+H0_G*G_D
       H0_S  = H0_T*T_S+H0_E*E_S+H0_G*G_S
       H0_N  = H0_T*T_N
!      == TRANSFORM H1
       H1_DD = H1_TT*T_D**2 + 2.D0*H1_TR*T_D*R_D + 2.D0*H1_TG*T_D*G_D + H1_T*T_DD &
      &      + H1_RR*R_D**2 + 2.D0*H1_RG*R_D*G_D + H1_R*R_DD &
      &      + H1_GG*G_D**2 + H1_G*G_DD
       H1_SS = H1_TT*T_S**2 + 2.D0*H1_TG*T_S*G_S + H1_T*T_SS + H1_GG*G_S**2 + H1_G*G_SS
       H1_NN = H1_TT*T_N**2 + H1_T*T_NN
       H1_DS = H1_TT*T_D*T_S + H1_TR*T_S*R_D + H1_TG*(T_D*G_S+T_S*G_D) + H1_T*T_DS &
      &      + H1_RG*R_D*G_S + H1_GG*G_D*G_S + H1_G*G_DS 
       H1_DN = H1_TT*T_D*T_N + H1_TR*T_N*R_D + H1_TG*T_N*G_D + H1_T*T_DN
       H1_SN = H1_TT*T_S*T_N + H1_TG*T_N*G_S + H1_T*T_SN

       H1_D  = H1_T*T_D+H1_R*R_D+H1_G*G_D
       H1_S  = H1_T*T_S+H1_G*G_S
       H1_N  = H1_T*T_N
!
!      == EXCHANGE ENERGY ==============================================
       EPSC   =RHOT*(E+H0+H1)
       EPSC_D =RHOT*(E_D+H0_D+H1_D)+E+H0+H1
       EPSC_S =RHOT*(E_S+H0_S+H1_S)
       EPSC_N =RHOT*(H0_N+H1_N)
       EPSC_DD=RHOT*(E_DD+H0_DD+H1_DD)+2.D0*(E_D+H0_D+H1_D)
       EPSC_DS=RHOT*(E_DS+H0_DS+H1_DS)+(E_S+H0_S+H1_S)
       EPSC_DN=RHOT*(H0_DN+H1_DN)+(H0_N+H1_N)
       EPSC_SS=RHOT*(E_SS+H0_SS+H1_SS)
       EPSC_SN=RHOT*(H0_SN+H1_SN)
       EPSC_NN=RHOT*(H0_NN+H1_NN)
       RETURN
       END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE PBE_MODULE
!*******************************************************************************
!**                                                                   **
!**  NAME: PBE                                                        **
!**                                                                   **
!**  PURPOSE: IMPLEMENTS THE CORRELATION FUNCTIONAL                   **
!**    OF PERDEW,BURKE,ERNZERHOF (PHYS. REV. LETT.SUBMITTED 1996)     **
!**    (INCLUDES LOCAL AND NONLOCAL CORRELATION, BUT NO EXCHANGE)     **
!**    THIS FUNCTIONAL IS A SIMPLIFIED VERSION OF THE PW91-GGA        **
!**                                                                   **
!**  FUNCTIONS                                                        **
!**    PBE$EVAL1                                                      **
!**    PBE$EVAL2                                                      **
!**    PBE$EVAL3                                                      **
!**                                                                   **
!**  DEPENDENCIES:                                                    **
!**    PERDEWWANG91L(PERDEWWANG91L$EPSVAL1,PERDEWWANG91L$EPSVAL2      **
!**                 ,PERDEWWANG91L$EPSVAL3                            **
!**                                                                   **
!*******************************************************************************
!*******************************************************************************
!**  BECKE88 EXCHANGE (PRA38,3098(1988) IS OBTAINED                           **
!**  WITH THE FOLLOWING PARAMETERS                                            **
!**  MU   =16.D0*PI/3.D0*(6.D0*PI**2)**(1.D0/3.D0)*0.0042D0                   **
!**       =0.2742931D0                                                        **
!**  KAPPA=2.D0*PI/9.D0/(6.D0*PI**2)**(1.D0/3.D0)                             **
!**       =0.1791102D0                                                        **
!**  MU=16.D0*PI/3.D0*(6.D0*PI**2)**(1.D0/3.D0)*0.0042D0                      **
!**  KAPPA=2.D0*PI/9.D0/(6.D0*PI**2)**(1.D0/3.D0)                             **
!*******************************************************************************
LOGICAL(4)         :: TINI=.FALSE.
REAL(8)            :: PI
REAL(8)            :: RSFAC   ! RS=RSFAC*RHO**(-1/3)
REAL(8)            :: KFFAC   ! KF= KFFAC*RHO**(4/3)
REAL(8)            :: KSFAC   ! KS=KSFAC/SQRT(RS)
REAL(8)            :: TFAC    ! T2=TFAC*GRHO2/G**2*RHOT**(-7/3)
REAL(8)            :: GAMMA   ! (1-LOG(2))/PI**2
REAL(8)            :: KAPPA=0.804D0
REAL(8)            :: ALPHAB2  ! USED TO OBTAIN E_X=RHO/2R SCALING
REAL(8),PARAMETER  :: CC0=4.235D-3
REAL(8)            :: NU      
REAL(8)            :: MU
REAL(8)            :: MUBYKAPPA
REAL(8)            :: EXFAC
REAL(8)            :: BETA    ! BETA=NU*CC0 (PWG:EQ.13)
LOGICAL(4)         :: T_BUGFIX1=.TRUE.
END MODULE PBE_MODULE
!
!      .................................................................
       SUBROUTINE PBE_INITIALIZE
       USE PBE_MODULE
       IMPLICIT NONE 
       PI=4.D0*ATAN(1.D0)
       RSFAC=(3.D0/(4.D0*PI))**(1.D0/3.D0) ! RS=RSFAC*RHOT**(-1/3)
       KFFAC=(3.D0*PI**2)**(1.D0/3.D0)     ! KF=KFFAC/RS
       KSFAC=SQRT(4.D0*KFFAC/PI)           ! KS=KSFAC/SQRT(RS)
       IF(T_BUGFIX1) THEN
         TFAC=PI/(16.D0*(3.D0*PI**2)**(1.D0/3.D0))
       ELSE 
!        == BUG FIX ON DEC.10.2023 THE NEW VERSION IS ABOVE ====================
!        == TFAC IS USED FOR THE DIMENSION-LESS GRADIENT AS DEFINED BELOW EQ.3==
!        == IN THE PBE PAPER PERDEW96_PRL77_3865                              ==
!        == THE OLD VERSION IS BELOW:                                         ==
         TFAC=0.25D0*RSFAC/KSFAC**2        ! T2=TFAC*GRHO2/G**2*RHOT**(-7/3)
       END IF
!      =========================================================================
       GAMMA=(1.D0-LOG(2.D0))/PI**2
       NU=16.D0/PI*(3*PI**2)**(1.D0/3.D0)
       BETA=NU*CC0
       MU=BETA*PI**2/3.D0
!      ==================================================================
!      ==  BECKE88 EXCHANGE (PRA38,3098(1988) IS OBTAINED              ==
!      ==  WITH THE FOLLOWING PARAMETERS                               ==
!      ==  MU   =16.D0*PI/3.D0*(6.D0*PI**2)**(1.D0/3.D0)*0.0042D0      ==
!      ==       =0.2742931D0                                           ==
!      ==  KAPPA=2.D0*PI/9.D0/(6.D0*PI**2)**(1.D0/3.D0)                ==
!      ==       =0.1791102D0                                           ==
!      MU=16.D0*PI/3.D0*(6.D0*PI**2)**(1.D0/3.D0)*0.0042D0
!      KAPPA=2.D0*PI/9.D0/(6.D0*PI**2)**(1.D0/3.D0)
!      ==================================================================
       MUBYKAPPA=MU/KAPPA
       EXFAC=-3.D0/(4.D0*PI)       !EX_LOC(NT)=EXFAC*RHO**(4/3)
       ALPHAB2=(4.D0*PI/(9.D0*KAPPA))**2
       TINI=.TRUE.
       RETURN
       END SUBROUTINE PBE_INITIALIZE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE PBE$SETL4(ID,VAL)
       USE PBE_MODULE, ONLY : T_BUGFIX1
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       LOGICAL(4)  ,INTENT(IN) :: VAL
!      *************************************************************************
       IF(ID.EQ.'BUGFIX1') THEN
         T_BUGFIX1=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('PBE$SETL4')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PBE_H0OFTEG1(T2,E,G,H,H_T,H_E,H_G)
!      *****************************************************************    
!      **                                                             **    
!      **                                                             **    
!      **  THE PARAMATERS MUST LIE IN THE FOLLOWING RANGE:            **    
!      **    T2>0; E<0 G>0                                            **    
!      **                                                             **    
!      *****************************************************************    
       USE PBE_MODULE, ONLY: GAMMA
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: T2
       REAL(8),INTENT(IN) :: E
       REAL(8),INTENT(IN) :: G
       REAL(8),INTENT(OUT):: H
       REAL(8),INTENT(OUT):: H_T,H_E,H_G
       REAL(8)            :: FAC2
       REAL(8)            :: V,V_E,V_G
       REAL(8)            :: SVAR,DSVAR
       REAL(8)            :: G2,G3
       REAL(8)            :: SVAR1,SVAR1_T
       REAL(8)            :: SVAR1_V
       REAL(8)            :: H_V
!      *****************************************************************    
!
!      =================================================================
!      ==  CONSTANTS                                                  ==
!      =================================================================
       FAC2=1.D0/GAMMA   !(2*ALPHA/BETA**2)
!
!      =================================================================
!      ==  CALCULATE A(VAR1=E/G3)                                     ==
!      ==  G>0; EPSC=<0 ; 0=< AINV < INFTY                            ==
!      ==  AINV BECOMES LARGE FOR LARGE DENSITIES                     ==
!      =================================================================
       G2=G*G
       G3=G2*G
       V=E/G3
       CALL PBE_FOFTV1(T2,V,SVAR1,SVAR1_T,SVAR1_V)
!
!      =================================================================
!      ==  CALCULATE H(T,V,G)                                        ==
!      =================================================================
       SVAR=LOG(1.D0+SVAR1)/FAC2
       DSVAR=1.D0/(1.D0+SVAR1)/FAC2
       H   =      G3*SVAR
       H_G = 3.D0*G2*SVAR
       H_T =      G3*DSVAR*SVAR1_T
       H_V =      G3*DSVAR*SVAR1_V
!
!      =================================================================
!      ==  TRANSFORM FROM H(T2,V(E,G),G) TO H(T2,E,G)                 ==
!      =================================================================
       V_E=1.D0/G3
       V_G=-3.D0*V/G
!
       H_G =H_G+H_V*V_G
       H_E =H_V*V_E
       RETURN
       END SUBROUTINE PBE_H0OFTEG1
!
!      .................................................................
       SUBROUTINE PBE_H0OFTEG2(T2,E,G,H,H_T,H_E,H_G &
      &                       ,H_TT,H_TE,H_TG,H_EE,H_EG,H_GG)
!      *****************************************************************    
!      **                                                             **    
!      **                                                             **    
!      **  THE PARAMATERS MUST LIE IN THE FOLLOWING RANGE:            **    
!      **    T2>0; E<0 G>0                                            **    
!      **                                                             **    
!      *****************************************************************    
       USE PBE_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: T2
       REAL(8),INTENT(IN) :: E
       REAL(8),INTENT(IN) :: G
       REAL(8),INTENT(OUT):: H
       REAL(8),INTENT(OUT):: H_T,H_E,H_G
       REAL(8),INTENT(OUT):: H_TT,H_TE,H_TG,H_EE,H_EG,H_GG
       REAL(8)            :: FAC2
       REAL(8)            :: V,V_E,V_G,V_EG,V_GG
       REAL(8)            :: SVAR,DSVAR,D2SVAR
       REAL(8)            :: G2,G3
       REAL(8)            :: SVAR1,SVAR1_T,SVAR1_V,SVAR1_TT,SVAR1_TV,SVAR1_VV
       REAL(8)            :: H_V,H_TV,H_VV,H_VG
!      *****************************************************************    
!
!      =================================================================
!      ==  CONSTANTS                                                  ==
!      =================================================================
       FAC2=1.D0/GAMMA   !(2*ALPHA/BETA**2)
!
!      =================================================================
!      ==  CALCULATE A(VAR1=E/G3)                                     ==
!      ==  G>0; EPSC=<0 ; 0=< AINV < INFTY                            ==
!      ==  AINV BECOMES LARGE FOR LARGE DENSITIES                     ==
!      =================================================================
       G2=G*G
       G3=G2*G
       V=E/G3
       CALL PBE_FOFTV2(T2,V,SVAR1,SVAR1_T,SVAR1_V &
      &                    ,SVAR1_TT,SVAR1_TV,SVAR1_VV)
!
!      =================================================================
!      ==  CALCULATE H(T,V,G)                                        ==
!      =================================================================
       SVAR=LOG(1.D0+SVAR1)/FAC2
       DSVAR=1.D0/(1.D0+SVAR1)/FAC2
       D2SVAR=-FAC2*DSVAR*DSVAR
       H   =      G3*SVAR
       H_G = 3.D0*G2*SVAR
       H_T =      G3*DSVAR*SVAR1_T
       H_V =      G3*DSVAR*SVAR1_V
       H_TT=      G3*(D2SVAR*SVAR1_T**2+DSVAR*SVAR1_TT)
       H_TV=      G3*(D2SVAR*SVAR1_T*SVAR1_V+DSVAR*SVAR1_TV)
       H_TG= 3.D0*G2*DSVAR*SVAR1_T
       H_VV=      G3*(D2SVAR*SVAR1_V**2+DSVAR*SVAR1_VV)
       H_VG= 3.D0*G2*DSVAR*SVAR1_V
       H_GG= 6.D0*G *SVAR
!
!      =================================================================
!      ==  TRANSFORM FROM H(T2,V(E,G),G) TO H(T2,E,G)                 ==
!      =================================================================
       V_E=1.D0/G3
       V_G=-3.D0*V/G
       V_EG=-3.D0*V_E/G
       V_GG=-4.D0*V_G/G
!
       H_TE=H_TV*V_E
       H_TG=H_TG+H_TV*V_G
       H_EE=H_VV*V_E**2
       H_EG=H_VG*V_E+H_VV*V_G*V_E+H_V*V_EG
       H_GG=H_GG+2.D0*H_VG*V_G + H_VV*V_G**2+H_V*V_GG
       H_G =H_G+H_V*V_G
       H_E =H_V*V_E
       RETURN
       END SUBROUTINE PBE_H0OFTEG2
!
!      .................................................................
       SUBROUTINE PBE_H0OFTEG3(T2,E,G,H,H_T,H_E,H_G &
      &                       ,H_TT,H_TE,H_TG,H_EE,H_EG,H_GG &
      &                       ,H_TTT,H_TTE,H_TTG,H_TEE,H_TEG,H_TGG &
      &                       ,H_EEE,H_EEG,H_EGG,H_GGG)
!      *****************************************************************    
!      **                                                             **    
!      **                                                             **    
!      **  THE PARAMATERS MUST LIE IN THE FOLLOWING RANGE:            **    
!      **    T2>0; E<0 G>0                                            **    
!      **                                                             **    
!      *****************************************************************    
       USE PBE_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: T2
       REAL(8),INTENT(IN) :: E
       REAL(8),INTENT(IN) :: G
       REAL(8),INTENT(OUT):: H
       REAL(8),INTENT(OUT):: H_T,H_E,H_G
       REAL(8),INTENT(OUT):: H_TT,H_TE,H_TG,H_EE,H_EG,H_GG
       REAL(8),INTENT(OUT):: H_TTT,H_TTE,H_TTG,H_TEE,H_TEG,H_TGG
       REAL(8),INTENT(OUT):: H_EEE,H_EEG,H_EGG,H_GGG
       REAL(8)            :: FAC2
       REAL(8)            :: V,V_E,V_G,V_EG,V_GG,V_EGG,V_GGG
       REAL(8)            :: SVAR,DSVAR,D2SVAR,D3SVAR
       REAL(8)            :: G2,G3
       REAL(8)            :: SVAR1,SVAR1_T,SVAR1_V,SVAR1_TV,SVAR1_VV,SVAR1_TT
       REAL(8)            :: SVAR1_TTT,SVAR1_TTV,SVAR1_TVV,SVAR1_VVV
       REAL(8)            :: H_V,H_TV,H_VV,H_VG
       REAL(8)            :: H_TTV,H_TVV,H_TVG,H_VVV,H_VVG,H_VGG
!      *****************************************************************    
!
!      =================================================================
!      ==  CONSTANTS                                                  ==
!      =================================================================
       FAC2=1.D0/GAMMA   !(2*ALPHA/BETA**2)
!
!      =================================================================
!      ==  CALCULATE A(VAR1=E/G3)                                     ==
!      ==  G>0; EPSC=<0 ; 0=< AINV < INFTY                            ==
!      ==  AINV BECOMES LARGE FOR LARGE DENSITIES                     ==
!      =================================================================
       G2=G*G
       G3=G2*G
       V=E/G3
       CALL PBE_FOFTV3(T2,V,SVAR1,SVAR1_T,SVAR1_V &
      &                    ,SVAR1_TT,SVAR1_TV,SVAR1_VV &
      &                    ,SVAR1_TTT,SVAR1_TTV,SVAR1_TVV,SVAR1_VVV)
!
!      =================================================================
!      ==  CALCULATE H(T,V,G)                                        ==
!      =================================================================
       SVAR=LOG(1.D0+SVAR1)/FAC2
       DSVAR=1.D0/(1.D0+SVAR1)/FAC2
       D2SVAR=-FAC2*DSVAR*DSVAR
       D3SVAR=-2.D0*FAC2*DSVAR*D2SVAR
       H    =      G3*SVAR
       H_G  = 3.D0*G2*SVAR
       H_T  =      G3*DSVAR*SVAR1_T
       H_V  =      G3*DSVAR*SVAR1_V
       H_TT =      G3*(D2SVAR*SVAR1_T**2     +DSVAR*SVAR1_TT)
       H_TV =      G3*(D2SVAR*SVAR1_T*SVAR1_V+DSVAR*SVAR1_TV)
       H_VV =      G3*(D2SVAR*SVAR1_V**2     +DSVAR*SVAR1_VV)
       H_TG = 3.D0*G2*DSVAR*SVAR1_T
       H_VG = 3.D0*G2*DSVAR*SVAR1_V
       H_GG = 6.D0*G *SVAR
       H_TTT=      G3*(D3SVAR*SVAR1_T**3 &
      &               +3.D0*D2SVAR*SVAR1_T*SVAR1_TT &
      &               +DSVAR*SVAR1_TTT)
       H_TTV=      G3*(D3SVAR*SVAR1_T**2*SVAR1_V &
      &               +D2SVAR*(SVAR1_TT*SVAR1_V+2.D0*SVAR1_T*SVAR1_TV) &
      &               +DSVAR*SVAR1_TTV)
       H_TVV=      G3*(D3SVAR*SVAR1_T*SVAR1_V**2 &
      &               +D2SVAR*(2.D0*SVAR1_TV*SVAR1_V + SVAR1_T*SVAR1_VV) &
      &               +DSVAR*SVAR1_TVV)
       H_VVV=      G3*(D3SVAR*SVAR1_V**3 &
      &               +D2SVAR*3.D0*SVAR1_V*SVAR1_VV &
      &               +DSVAR*SVAR1_VVV)
       H_TTG= 3.D0*G2*(D2SVAR*SVAR1_T**2+DSVAR*SVAR1_TT)
       H_TVG= 3.D0*G2*(D2SVAR*SVAR1_T*SVAR1_V+DSVAR*SVAR1_TV)
       H_VVG= 3.D0*G2*(D2SVAR*SVAR1_V**2+DSVAR*SVAR1_VV)
       H_TGG= 6.D0*G *DSVAR*SVAR1_T
       H_VGG= 6.D0*G *DSVAR*SVAR1_V
       H_GGG= 6.D0   *SVAR
!
!      =================================================================
!      ==  TRANSFORM FROM H(T2,V(E,G),G) TO H(T2,E,G)                 ==
!      =================================================================
       V_E=1.D0/G3
       V_G=-3.D0*V/G
       V_EG=-3.D0*V_E/G
       V_GG=-4.D0*V_G/G
       V_EGG=-4.D0*V_EG/G
       V_GGG=-5.D0*V_GG/G
!
       H_TTE= H_TTV*V_E
       H_TTG= H_TTG+H_TTV*V_G
       H_TEE= H_TVV*V_E**2    ! V_EE=0
       H_TEG= (H_TVG + H_TVV*V_G)*V_E+H_TV*V_EG
       H_TGG= H_TGG + 2.D0*H_TVG*V_G + H_TVV*V_G**2 + H_TV*V_GG
       H_EEE= H_VVV*V_E**3     !V_EE=0
       H_EEG= (H_VVG+H_VVV*V_G)*V_E**2 + H_VV*2.D0*V_E*V_EG
       H_EGG= (H_VGG + 2.D0*H_VVG*V_G + H_VVV*V_G**2 + H_VV*V_GG)*V_E &
      &       +2.D0*(H_VG + H_VV*V_G)*V_EG + H_V*V_EGG
       H_GGG = H_GGG + 2.D0*H_VGG*V_G + H_VVG*V_G**2 + H_VG*V_GG &
      &      +(H_VGG + 2.D0*H_VVG*V_G + H_VVV*V_G**2 + H_VV*V_GG)*V_G &
      &              + 2.D0*(H_VG+H_VV*V_G)*V_GG     + H_V*V_GGG
       H_TE = H_TV*V_E
       H_TG = H_TG + H_TV*V_G
       H_EE = H_VV*V_E**2
       H_EG = H_VG*V_E + H_VV*V_G*V_E + H_V*V_EG
       H_GG = H_GG + 2.D0*H_VG*V_G + H_VV*V_G**2 + H_V*V_GG
       H_G  = H_G + H_V*V_G
       H_E  = H_V*V_E
       RETURN
       END SUBROUTINE PBE_H0OFTEG3
!
!      ..............................................................
       SUBROUTINE PBE_FOFTV1(T2,V,F,F_T,F_V)
       USE PBE_MODULE, ONLY : BETA,GAMMA
       IMPLICIT NONE
       REAL(8) ,INTENT(IN) :: T2
       REAL(8) ,INTENT(IN) :: V
       REAL(8) ,INTENT(OUT):: F
       REAL(8) ,INTENT(OUT):: F_T,F_V
       REAL(8)             :: FAC1,FAC2
       REAL(8)             :: SVAR,SVAR2
       REAL(8)             :: AI,AI_V
       REAL(8)             :: AI2,AI3,T4,T6
       REAL(8)             :: F_A
!      ***************************************************************
       FAC2=1.D0/GAMMA   !(2*ALPHA/BETA**2)
       FAC1=BETA/GAMMA
!
!      =================================================================
!      == CAPTURE OVERFLOW OF THE EXPONENTIAL (1/A->INFTY)            ==       
!      =================================================================
       IF(-FAC2*V.GT.100.D0) THEN
         F=FAC1*T2
         F_T=FAC1
         F_V=0.D0
         RETURN
       END IF
!
!      =================================================================
!      == 1/A=(EXP(-FAC2*V)-1)/FAC1 AND DERIVATIVES WITH RESPECT TO V ==       
!      =================================================================
       AI    =EXP(-FAC2*V)/FAC1
       AI_V  =-FAC2*AI
       AI=AI-1.D0/FAC1
!
!      =================================================================
!      == 0<SVAR<1.D0                                                 ==       
!      =================================================================
       AI2=AI*AI
       AI3=AI2*AI
       T4=T2*T2
       T6=T4*T2
       SVAR=1.D0/(AI2+AI*T2+T4)
       SVAR2=SVAR*SVAR
!
!      ==================================================================
!      == MULTIPLY RESULT WITH FAC1. IT WILL PROPAGATE FROM HERE  =======
!      ==================================================================
       SVAR=FAC1*SVAR
       SVAR2=FAC1*SVAR2
!
!      ==================================================================
!      == CALCULATE F AND DERIVATIVES WITH RESPECT TO AI AND T2        ==
!      ==================================================================
       F   =AI*T2*(AI+T2)*SVAR
       F_T =AI3*(2.D0*T2+AI)*SVAR2
       F_A =T6*(2.D0*AI+T2)*SVAR2
!
!      ==================================================================
!      == NOW TRANSFORM F(AI,T2) INTO F(T2,V(AI))                      ==
!      ==================================================================
       F_V   =F_A*AI_V
       RETURN
       END
!
!      ..............................................................
       SUBROUTINE PBE_FOFTV2(T2,V,F,F_T,F_V,F_TT,F_TV,F_VV)
       USE PBE_MODULE, ONLY : BETA,GAMMA
       IMPLICIT NONE
       REAL(8) ,INTENT(IN) :: T2
       REAL(8) ,INTENT(IN) :: V
       REAL(8) ,INTENT(OUT):: F
       REAL(8) ,INTENT(OUT):: F_T,F_V
       REAL(8) ,INTENT(OUT):: F_TT,F_TV,F_VV
       REAL(8)             :: FAC1,FAC2
       REAL(8)             :: SVAR,SVAR2,SVAR3
       REAL(8)             :: AI,AI_V,AI_VV
       REAL(8)             :: AI2,AI3,T4,T6
       REAL(8)             :: F_A,F_AA,F_AT
       REAL(8)             :: C2
!      ***************************************************************
       FAC2=1.D0/GAMMA   !(2*ALPHA/BETA**2)
       FAC1=BETA/GAMMA
!
!      =================================================================
!      == CAPTURE OVERFLOW OF THE EXPONENTIAL (1/A->INFTY)            ==       
!      =================================================================
       IF(-FAC2*V.GT.100.D0) THEN
         F=FAC1*T2
         F_T=FAC1
         F_V=0.D0
         F_TT=0.D0
         F_TV=0.D0
         F_VV=0.D0
         RETURN
       END IF
!
!      =================================================================
!      == 1/A=(EXP(-FAC2*V)-1)/FAC1 AND DERIVATIVES WITH RESPECT TO V ==       
!      =================================================================
       AI    =EXP(-FAC2*V)/FAC1
       AI_V  =-FAC2*AI
       AI_VV =-FAC2*AI_V
       AI=AI-1.D0/FAC1
!
!      =================================================================
!      == 0<SVAR<1.D0                                                 ==       
!      =================================================================
       AI2=AI*AI
       AI3=AI2*AI
       T4=T2*T2
       T6=T4*T2
       SVAR=1.D0/(AI2+AI*T2+T4)
       SVAR2=SVAR*SVAR
       SVAR3=SVAR2*SVAR
!
!      ==================================================================
!      == MULTIPLY RESULT WITH FAC1. IT WILL PROPAGATE FROM HERE  =======
!      ==================================================================
       SVAR=FAC1*SVAR
       SVAR2=FAC1*SVAR2
       SVAR3=FAC1*SVAR3
!      == PREFACTORS FOR SECOND AND THIRD DERIVATIVES
       C2=6.D0*AI*T2*(AI+T2)*SVAR3
!
!      ==================================================================
!      == CALCULATE F AND DERIVATIVES WITH RESPECT TO AI AND T2        ==
!      ==================================================================
       F   =AI*T2*(AI+T2)*SVAR
       F_T =AI3*(2.D0*T2+AI)*SVAR2
       F_A =T6*(2.D0*AI+T2)*SVAR2
       F_TT=-C2*AI2
       F_AT=+C2*AI*T2
       F_AA=-C2*T4
!
!      ==================================================================
!      == NOW TRANSFORM F(AI,T2) INTO F(T2,V(AI))                      ==
!      ==================================================================
       F_V   =F_A*AI_V
       F_VV  =F_AA*AI_V**2+F_A*AI_VV
       F_TV  =F_AT*AI_V
       RETURN
       END
!
!      ..............................................................
       SUBROUTINE PBE_FOFTV3(T2,V,F,F_T,F_V &
      &                     ,F_TT,F_TV,F_VV &
      &                     ,F_TTT,F_TTV,F_TVV,F_VVV)
       USE PBE_MODULE, ONLY : BETA,GAMMA
       IMPLICIT NONE
       REAL(8) ,INTENT(IN) :: T2
       REAL(8) ,INTENT(IN) :: V
       REAL(8) ,INTENT(OUT):: F
       REAL(8) ,INTENT(OUT):: F_T,F_V
       REAL(8) ,INTENT(OUT):: F_TT,F_TV,F_VV
       REAL(8) ,INTENT(OUT):: F_TTT,F_TTV,F_TVV,F_VVV
       REAL(8)             :: FAC1,FAC2
       REAL(8)             :: SVAR,SVAR2,SVAR3,SVAR4
       REAL(8)             :: AI,AI_V,AI_VV,AI_VVV
       REAL(8)             :: AI2,AI3,T4,T6
       REAL(8)             :: F_A,F_AA,F_AT,F_AAA,F_AAT,F_ATT
       REAL(8)             :: C2,C2_A,C2_T
!      ***************************************************************
       FAC2=1.D0/GAMMA   !(2*ALPHA/BETA**2)
       FAC1=BETA/GAMMA
!
!      =================================================================
!      == CAPTURE OVERFLOW OF THE EXPONENTIAL                         ==       
!      == V->-INFTY; 1/A->+INFTY; A->0                                ==       
!      =================================================================
       IF(-FAC2*V.GT.100.D0) THEN
         F=FAC1*T2
         F_T=FAC1
         F_V=0.D0
         F_TT=0.D0
         F_TV=0.D0
         F_VV=0.D0
         F_TTT=0.D0
         F_TTV=0.D0
         F_TVV=0.D0
         F_VVV=0.D0
         RETURN
       END IF
!
!      =================================================================
!      == 1/A=(EXP(-FAC2*V)-1)/FAC1 AND DERIVATIVES WITH RESPECT TO V ==       
!      =================================================================
       AI    =EXP(-FAC2*V)/FAC1
       AI_V  =-FAC2*AI
       AI_VV =-FAC2*AI_V
       AI_VVV=-FAC2*AI_VV
       AI=AI-1.D0/FAC1
!
!      =================================================================
!      == 0<SVAR<1.D0                                                 ==       
!      =================================================================
       AI2=AI*AI
       AI3=AI2*AI
       T4=T2*T2
       T6=T4*T2
       SVAR=1.D0/(AI2+AI*T2+T4)
       SVAR2=SVAR*SVAR
       SVAR3=SVAR2*SVAR
       SVAR4=SVAR2*SVAR2
!
!      ==================================================================
!      == MULTIPLY RESULT WITH FAC1. IT WILL PROPAGATE FROM HERE  =======
!      ==================================================================
       SVAR=FAC1*SVAR
       SVAR2=FAC1*SVAR2
       SVAR3=FAC1*SVAR3
       SVAR4=FAC1*SVAR4
!      == PREFACTORS FOR SECOND AND THIRD DERIVATIVES
       C2=6.D0*AI*T2*(AI+T2)*SVAR3
       C2_A=6.D0*T2*(2.D0*AI+T2)*(T4 -2.D0*AI*(AI+T2))*SVAR4
       C2_T=6.D0*AI*(AI+2.D0*T2)*(AI2-2.D0*T2*(AI+T2))*SVAR4
!
!      ==================================================================
!      == CALCULATE F AND DERIVATIVES WITH RESPECT TO AI AND T2        ==
!      ==================================================================
       F   =AI*T2*(AI+T2)*SVAR
       F_T =AI3*(2.D0*T2+AI)*SVAR2
       F_A =T6*(2.D0*AI+T2)*SVAR2
       F_TT=-C2*AI2
       F_AT=+C2*AI*T2
       F_AA=-C2*T4
       F_AAA=-C2_A*T4
       F_AAT=+C2_A*AI*T2+C2*T2
       F_ATT=+C2_T*AI*T2+C2*AI
       F_TTT=-C2_T*AI2
!
!      ==================================================================
!      == NOW TRANSFORM F(AI,T2) INTO F(T2,V(AI))                      ==
!      ==================================================================
       F_V   =F_A*AI_V
       F_VV  =F_AA*AI_V**2+F_A*AI_VV
       F_TV  =F_AT*AI_V
       F_VVV =F_AAA*AI_V**3+3.D0*F_AA*AI_V*AI_VV+F_A*AI_VVV
       F_TVV =F_AAT*AI_V**2+F_AT*AI_VV
       F_TTV =F_ATT*AI_V
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PBE$EVAL1(VAL,EXC,DEXC)
       USE PBE_MODULE,ONLY : TINI &
      &                     ,TFAC
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: VAL(5)
       REAL(8),INTENT(OUT):: EXC
       REAL(8),INTENT(OUT):: DEXC(5)
       REAL(8)            :: RHOT    !TOTAL ELECTRON DENSITY
       REAL(8)            :: RHOS    !SPIN DENSITY
       REAL(8)            :: GRHOT2  !SQUARED GRADIENT OF THE TOTAL DENSITY
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8),PARAMETER  :: ONEBY9=1.D0/9.D0
       REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
       REAL(8),PARAMETER  :: SEVENBY6=7.D0/6.D0
       REAL(8)            :: RHOT13
       REAL(8)            :: P13,M13,P23,M23
       REAL(8)            :: SIG,SIG_D,SIG_S
       REAL(8)            :: G,G_SIG
       REAL(8)            :: G_S,G_D
       REAL(8)            :: E,E_D,E_S
       REAL(8)            :: T,T_N,T_G,T_D
       REAL(8)            :: T_S
       REAL(8)            :: H0,H0_T,H0_E,H0_G
       REAL(8)            :: H0_D,H0_S,H0_N
!      *****************************************************************
       IF(.NOT.TINI) CALL PBE_INITIALIZE
       EXC=0.D0
       DEXC(:)=0.D0
       IF(VAL(1).LT.1.D-10) RETURN
       RHOT=VAL(1)
       RHOS=VAL(2)
       GRHOT2=VAL(3)
!
!      =========================================================================
!      == CALCULATE LOCAL CORRELATION. E=ECUNIV IN PW92:EQ3                   ==
!      =========================================================================
       CALL PERDEWWANG91L$EPSVAL(RHOT,RHOS,E,E_D,E_S)
!
!      =========================================================================
!      == EQ.11: SPIN DEPENDENCE G=PHI(SIGMA) IN PW92 DEFINED BELOW EQ.3      ==
!      =========================================================================
       SIG=RHOS/RHOT
       P13=(1.D0+SIG)**ONEBY3
       M13=(1.D0-SIG)**ONEBY3
       P23=P13*P13
       M23=M13*M13
       G=0.5D0*(P23+M23)
       G_SIG=ONEBY3*(1.D0/P13-1.D0/M13)
!
!      =========================================================================
!      == T_N IS T**2 OF EQ. 10 (WHICH REFERENCE?)                            ==
!      == T_N=(DIMENSIONLESS DENSITY GRADIENT T)**2                           ==
!      == SEE DEFINITION OF (DIMENSIONLESS DENSITY GRADIENT T) BELOW EQ.3     ==
!      == OF PERDEW96_PRL77_3865 (PBE PAPER)                                  ==
!      =========================================================================
       RHOT13=RHOT**ONEBY3
       T_N=TFAC/G**2/RHOT13**7
       T=T_N*GRHOT2
       T_G=-2.D0*T/G
       T_D=-7.D0/3.D0*T/RHOT
!
!      ================================================================
!      == LOCAL CORRELATION FOR SMALL T                              ==
!      ================================================================
       IF(T.LT.-1.D-8) THEN
         EXC     =RHOT*E
         DEXC(1) =RHOT*E_D+E
         DEXC(2) =RHOT*E_S
         RETURN
       END IF
!
!      ================================================================
!      == CALCULATE H0                                               ==
!      ================================================================
!      == CONSIDER H0 AS FUNCTION OF T2,E,G 
       CALL PBE_H0OFTEG1(T,E,G,H0,H0_T,H0_E,H0_G)
!
!      ================================================================
!      == CONVERT FROM RS TO RHOT                                    ==
!      ================================================================
!      == SIG(D=RHOT,S=RHOS)
       SIG_D=-SIG/RHOT
       SIG_S=1.D0/RHOT
!      == G(D=RHOT,S=RHOS)
       G_D=G_SIG*SIG_D
       G_S=G_SIG*SIG_S
!      == RS(RHOT) ====================================================
!      == T(D=RHOT,S=RHOS,N=GRHOT2) ====================================
!      == WATCH ORDER OF STATEMENTS!! =================================
!      T=SQRT(GRHOT2)/(2.D0*G*KSFAC/SQRT(RSFAC))*RHOT**(-SEVENBY6)
       T_S=T_G*G_S
       T_D=T_D+T_G*G_D
!      == TRANSFORM H0 ================================================        
       H0_D  = H0_T*T_D+H0_E*E_D+H0_G*G_D
       H0_S  = H0_T*T_S+H0_E*E_S+H0_G*G_S
       H0_N  = H0_T*T_N
!
!      == EXCHANGE ENERGY ==============================================
       EXC    =RHOT*(E+H0)      !PW92:EQ.3
       DEXC(1) =RHOT*(E_D+H0_D)+E+H0
       DEXC(2) =RHOT*(E_S+H0_S)
       DEXC(3)=RHOT*H0_N
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PBE$EVAL2(VAL,EXC,DEXC,D2EXC)
       USE PBE_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: VAL(5)
       REAL(8),INTENT(OUT):: EXC
       REAL(8),INTENT(OUT):: DEXC(5)
       REAL(8),INTENT(OUT):: D2EXC(5,5)
       REAL(8)            :: RHOT
       REAL(8)            :: RHOS
       REAL(8)            :: GRHOT2
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8),PARAMETER  :: ONEBY9=1.D0/9.D0
       REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
       REAL(8),PARAMETER  :: SEVENBY6=7.D0/6.D0
       REAL(8)            :: RHOT13
       REAL(8)            :: P13,M13,P23,M23,P43,M43
       REAL(8)            :: SIG,SIG_D,SIG_S,SIG_DD,SIG_DS
       REAL(8)            :: G,G_SIG,G_SIGSIG
       REAL(8)            :: G_S,G_D,G_DD,G_DS,G_SS
       REAL(8)            :: E,E_D,E_S,E_DD,E_DS,E_SS
       REAL(8)            :: T,T_N,T_G,T_D,T_DD,T_DN,T_DG,T_GG,T_GN,T_NN
       REAL(8)            :: T_S,T_DS,T_SS,T_SN
       REAL(8)            :: H0,H0_T,H0_E,H0_G
       REAL(8)            :: H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG
       REAL(8)            :: H0_D,H0_S,H0_N
       REAL(8)            :: H0_DD,H0_DS,H0_DN,H0_SS,H0_SN,H0_NN
!      *****************************************************************
       IF(.NOT.TINI) CALL PBE_INITIALIZE
       EXC=0.D0
       DEXC(:)=0.D0
       D2EXC(:,:)=0.D0
       IF(VAL(1).LT.1.D-10) RETURN
       RHOT=VAL(1)
       RHOS=VAL(2)
       GRHOT2=VAL(3)
!
!      ================================================================
!      == CALCULATE LOCAL CORRELATION                                ==
!      ================================================================
       CALL PERDEWWANG91L$EPSVAL2(RHOT,RHOS,E,E_D,E_S,E_DD,E_DS,E_SS)
!
!      ================================================================
!      == EQ.11: SPIN DEPENDENCE G                                   ==
!      ================================================================
       SIG=RHOS/RHOT
       P13=(1.D0+SIG)**ONEBY3
       M13=(1.D0-SIG)**ONEBY3
       P23=P13*P13
       M23=M13*M13
       P43=P23*P23
       M43=M23*M23
       G=0.5D0*(P23+M23)
       G_SIG=ONEBY3*(1.D0/P13-1.D0/M13)
       G_SIGSIG=-ONEBY9*(1.D0/P43+1.D0/M43)
!
!      ================================================================
!      == T IS T**2 OF EQ. 10                                       ==
!      ================================================================
       RHOT13=RHOT**ONEBY3
       T_N=TFAC/G**2/RHOT13**7
       T=T_N*GRHOT2
       T_NN=0.D0
       T_G=-2.D0*T/G
       T_D=-7.D0/3.D0*T/RHOT
       T_DD=-10.D0/3.D0*T_D/RHOT
       T_DN=-7.D0/3.D0*T_N/RHOT
       T_DG=-2.D0*T_D/G
       T_GG=-3.D0*T_G/G
       T_GN=-2.D0*T_N/G
!
!      ================================================================
!      == LOCAL CORRELATION FOR SMALL T                              ==
!      ================================================================
       IF(T.LT.-1.D-8) THEN
         EXC     =RHOT*E
         DEXC(1) =RHOT*E_D+E
         DEXC(2) =RHOT*E_S
         D2EXC(1,1)=RHOT*E_DD+2.D0*E_D
         D2EXC(1,2)=RHOT*E_DS+E_S
         D2EXC(2,2)=RHOT*E_SS
         D2EXC(2,1)=D2EXC(1,2)
         RETURN
       END IF
!
!      ================================================================
!      == CALCULATE H0                                               ==
!      ================================================================
!      == CONSIDER H0 AS FUNCTION OF T2,E,G 
       CALL PBE_H0OFTEG2(T,E,G,H0,H0_T,H0_E,H0_G &
      &            ,H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG)
!
!      ================================================================
!      == CONVERT FROM RS TO RHOT                                    ==
!      ================================================================
!      == SIG(D=RHOT,S=RHOS)
       SIG_D=-SIG/RHOT
       SIG_S=1.D0/RHOT
       SIG_DD=-2.D0*SIG_D/RHOT
       SIG_DS=-SIG_S/RHOT
!      == G(D=RHOT,S=RHOS)
       G_D=G_SIG*SIG_D
       G_S=G_SIG*SIG_S
       G_DD=G_SIGSIG*SIG_D**2+G_SIG*SIG_DD
       G_DS=G_SIGSIG*SIG_D*SIG_S+G_SIG*SIG_DS
       G_SS=G_SIGSIG*SIG_S**2
!      == RS(RHOT) ====================================================
!      == T(D=RHOT,S=RHOS,N=GRHOT2) ====================================
!      == WATCH ORDER OF STATEMENTS!! =================================
!      T=SQRT(GRHOT2)/(2.D0*G*KSFAC/SQRT(RSFAC))*RHOT**(-SEVENBY6)
       T_DD=T_DD+2.D0*T_DG*G_D+T_GG*G_D**2+T_G*G_DD
       T_DS=T_DG*G_S+T_GG*G_S*G_D+T_G*G_DS
       T_DN=T_DN+T_GN*G_D
       T_SS=T_GG*G_S**2+T_G*G_SS
       T_SN=T_GN*G_S
       T_S=T_G*G_S
       T_D=T_D+T_G*G_D
!      == TRANSFORM H0 ================================================        
       H0_DD = H0_TT*T_D**2 + 2.D0*H0_TE*T_D*E_D + 2.D0*H0_TG*T_D*G_D +H0_T*T_DD &
      &      + H0_EE*E_D**2 + 2.D0*H0_EG*E_D*G_D +H0_E*E_DD &
      &      + H0_GG*G_D**2 + H0_G*G_DD
       H0_SS = H0_TT*T_S**2 + 2.D0*H0_TE*T_S*E_S + 2.D0*H0_TG*T_S*G_S +H0_T*T_SS &
      &      + H0_EE*E_S**2 + 2.D0*H0_EG*E_S*G_S +H0_E*E_SS &
      &      + H0_GG*G_S**2 + H0_G*G_SS
       H0_NN = H0_TT*T_N**2  + H0_T*T_NN
       H0_DS = H0_TT*T_D*T_S + H0_TE*(T_D*E_S+T_S*E_D) + H0_TG*(T_D*G_S+T_S*G_D) +H0_T*T_DS &
      &      + H0_EE*E_D*E_S + H0_EG*(E_D*G_S+E_S*G_D) + H0_E*E_DS &
      &      + H0_GG*G_D*G_S + H0_G*G_DS
       H0_DN = H0_TT*T_D*T_N + H0_TE*T_N*E_D + H0_TG*T_N*G_D +H0_T*T_DN
       H0_SN = H0_TT*T_S*T_N + H0_TE*T_N*E_S + H0_TG*T_N*G_S +H0_T*T_SN
       H0_D  = H0_T*T_D+H0_E*E_D+H0_G*G_D
       H0_S  = H0_T*T_S+H0_E*E_S+H0_G*G_S
       H0_N  = H0_T*T_N
!
!      == EXCHANGE ENERGY ==============================================
       EXC   =RHOT*(E+H0)
       DEXC(1) =RHOT*(E_D+H0_D)+E+H0
       DEXC(2) =RHOT*(E_S+H0_S)
       DEXC(3)=RHOT*H0_N
       D2EXC(1,1)=RHOT*(E_DD+H0_DD)+2.D0*(E_D+H0_D)
       D2EXC(1,2)=RHOT*(E_DS+H0_DS)+(E_S+H0_S)
       D2EXC(1,3)=RHOT*H0_DN+H0_N
       D2EXC(2,2)=RHOT*(E_SS+H0_SS)
       D2EXC(2,3)=RHOT*H0_SN
       D2EXC(3,3)=RHOT*H0_NN
       D2EXC(2,1)=D2EXC(1,2)
       D2EXC(3,1)=D2EXC(1,3)
       D2EXC(3,2)=D2EXC(2,3)
       RETURN
       END
!
!      .................................................................
       SUBROUTINE PBE$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
       USE PBE_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: VAL(5)
       REAL(8),INTENT(OUT):: EXC
       REAL(8),INTENT(OUT):: DEXC(5)
       REAL(8),INTENT(OUT):: D2EXC(5,5)
       REAL(8),INTENT(OUT):: D3EXC(5,5,5)
       REAL(8)            :: RHOT
       REAL(8)            :: RHOS
       REAL(8)            :: GRHOT2
       REAL(8),PARAMETER  :: ONEBY3=1.D0/3.D0
       REAL(8),PARAMETER  :: ONEBY9=1.D0/9.D0
       REAL(8),PARAMETER  :: FOURBY3=4.D0/3.D0
       REAL(8),PARAMETER  :: FOURBY27=4.D0/27.D0
       REAL(8)            :: RHOT13
       REAL(8)            :: P13,M13,P23,M23,P43,M43,P73,M73
       REAL(8)            :: SIG,SIG_D,SIG_S,SIG_DD,SIG_DS,SIG_DDD,SIG_DDS
       REAL(8)            :: G,G_SIG,G_SIGSIG,G_SIGSIGSIG
       REAL(8)            :: G_S,G_D,G_DD,G_DS,G_SS
       REAL(8)            :: GM2,GM2_G,GM2_GG,GM2_GGG
       REAL(8)            :: RM7BY3,RM7BY3_D,RM7BY3_DD,RM7BY3_DDD
       REAL(8)            :: G_DDD,G_DDS,G_DSS,G_SSS
       REAL(8)            :: E,E_D,E_S,E_DD,E_DS,E_SS
       REAL(8)            :: E_DDD,E_DDS,E_DSS,E_SSS
       REAL(8)            :: T,T_N,T_G,T_D,T_DD,T_DN,T_DG,T_GG,T_GN,T_NN
       REAL(8)            :: T_DDD,T_DDG,T_DDN,T_DGG,T_DGN,T_GGG,T_GGN,T_NNN
       REAL(8)            :: T_S,T_DS,T_SS,T_SN
       REAL(8)            :: T_DDS,T_DSS,T_DSN,T_SSS,T_SSN
       REAL(8)            :: H0,H0_T,H0_E,H0_G
       REAL(8)            :: H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG
       REAL(8)            :: H0_TTT,H0_TTE,H0_TTG,H0_TEE,H0_TEG,H0_TGG
       REAL(8)            ::    H0_EEE,H0_EEG,H0_EGG,H0_GGG
       REAL(8)            :: H0_D,H0_S,H0_N
       REAL(8)            :: H0_DD,H0_DS,H0_DN,H0_SS,H0_SN,H0_NN
       REAL(8)            :: H0_DDD,H0_DDS,H0_DDN,H0_DSS,H0_DSN,H0_DNN
       REAL(8)            ::     H0_SSS,H0_SSN,H0_SNN,H0_NNN
!      *************************************************************************
       IF(.NOT.TINI) CALL PBE_INITIALIZE
       EXC=0.D0
       DEXC(:)=0.D0
       D2EXC(:,:)=0.D0
       D3EXC(:,:,:)=0.D0
       IF(VAL(1).LT.1.D-10) RETURN
       RHOT=VAL(1)
       RHOS=VAL(2)
       GRHOT2=VAL(3)
!
!      ================================================================
!      == EQ.11: SPIN DEPENDENCE G                                   ==
!      ================================================================
       SIG=RHOS/RHOT
       P13=(1.D0+SIG)**ONEBY3
       M13=(1.D0-SIG)**ONEBY3
       P23=P13*P13
       M23=M13*M13
       P43=P23*P23
       M43=M23*M23
       P73=P43*(1.D0+SIG)
       M73=M43*(1.D0-SIG)
       G=0.5D0*(P23+M23)
       G_SIG=ONEBY3*(1.D0/P13-1.D0/M13)
       G_SIGSIG=-ONEBY9*(1.D0/P43+1.D0/M43)
       G_SIGSIGSIG=FOURBY27*(1.D0/P73-1.D0/M73)
!
!      ================================================================
!      == T IS T**2 OF EQ. 10                                       ==
!      ================================================================
       RHOT13=RHOT**ONEBY3
       RM7BY3=1.D0/RHOT13**7
       GM2=1.D0/G**2
       GM2_G=-2.D0*GM2/G
       GM2_GG=-3.D0*GM2_G/G
       GM2_GGG=-4.D0*GM2_GG/G
       RM7BY3_D=-7.D0*ONEBY3*RM7BY3/RHOT
       RM7BY3_DD=-10.D0*ONEBY3*RM7BY3_D/RHOT
       RM7BY3_DDD=-13.D0*ONEBY3*RM7BY3_DD/RHOT
       T    =TFAC*GM2   *RM7BY3    *GRHOT2
       T_D  =TFAC*GM2   *RM7BY3_D  *GRHOT2
       T_G  =TFAC*GM2_G *RM7BY3    *GRHOT2
       T_N  =TFAC*GM2   *RM7BY3
       T_DD =TFAC*GM2   *RM7BY3_DD *GRHOT2
       T_DG =TFAC*GM2_G *RM7BY3_D  *GRHOT2
       T_DN =TFAC*GM2   *RM7BY3_D
       T_GG =TFAC*GM2_GG*RM7BY3    *GRHOT2
       T_GN =TFAC*GM2_G *RM7BY3
       T_NN=0.D0
       T_DDD=TFAC*GM2   *RM7BY3_DDD*GRHOT2
       T_DDG=TFAC*GM2_G *RM7BY3_DD *GRHOT2
       T_DDN=TFAC*GM2   *RM7BY3_DD
       T_DGG=TFAC*GM2_GG*RM7BY3_D  *GRHOT2
       T_DGN=TFAC*GM2_G *RM7BY3_D
       T_GGG=TFAC*GM2_GGG*RM7BY3    *GRHOT2
       T_GGN=TFAC*GM2_GG *RM7BY3
       T_NNN=0.D0
!
!      ================================================================
!      == CALCULATE LOCAL CORRELATION                                ==
!      ================================================================
       IF(T.LT.-1.D-8) THEN
         CALL PERDEWWANG91L$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
         RETURN
       ELSE
         CALL PERDEWWANG91L$EPSVAL3(RHOT,RHOS,E,E_D,E_S,E_DD,E_DS,E_SS &
      &                          ,E_DDD,E_DDS,E_DSS,E_SSS)
       END IF
!
!      ================================================================
!      == CALCULATE H0                                               ==
!      ================================================================
!      == CONSIDER H0 AS FUNCTION OF T2,E,G 
       CALL PBE_H0OFTEG3(T,E,G,H0,H0_T,H0_E,H0_G &
      &            ,H0_TT,H0_TE,H0_TG,H0_EE,H0_EG,H0_GG &
      &            ,H0_TTT,H0_TTE,H0_TTG,H0_TEE,H0_TEG,H0_TGG  &
      &           ,H0_EEE,H0_EEG,H0_EGG,H0_GGG)
!
!      ================================================================
!      == TRANSFORM G(SIG(D,S)) -> G(D,S)                            ==
!      ================================================================
!      == SIG(D=RHOT,S=RHOS)
       SIG_D=-SIG/RHOT
       SIG_S=1.D0/RHOT
       SIG_DD=-2.D0*SIG_D/RHOT
       SIG_DS=-SIG_S/RHOT
       SIG_DDD=-3.D0*SIG_DD/RHOT
       SIG_DDS=-2.D0*SIG_DS/RHOT
       G_D  =G_SIG*SIG_D
       G_S  =G_SIG*SIG_S
       G_DD =G_SIGSIG*SIG_D**2+G_SIG*SIG_DD
       G_DS =G_SIGSIG*SIG_D*SIG_S+G_SIG*SIG_DS
       G_SS =G_SIGSIG*SIG_S**2
       G_DDD=G_SIGSIGSIG*SIG_D**3      +G_SIGSIG*3.D0*SIG_D*SIG_DD+G_SIG*SIG_DDD
       G_DDS=G_SIGSIGSIG*SIG_D**2*SIG_S+G_SIGSIG*(2.D0*SIG_D*SIG_DS+SIG_S*SIG_DD)+G_SIG*SIG_DDS
       G_DSS=G_SIGSIGSIG*SIG_S**2*SIG_D+G_SIGSIG*2.D0*SIG_S*SIG_DS
       G_SSS=G_SIGSIGSIG*SIG_S**3
!
!      ==================================================================
!      == TRANSFORM T(D,G(D,S),N) -> T(D,S,N)                          ==
!      ==      D=RHOT; S=RHOS; N=GRHOT2                                ==
!      ==================================================================
       T_DDD= T_DDD + 3.D0*T_DDG*G_D + 2.D0*T_DGG*G_D**2 + 3.D0*T_DG*G_DD &
      &     + T_DGG*G_D**2 + T_GGG*G_D**3 + 3.D0*T_GG*G_D*G_DD + T_G*G_DDD 
       T_DDS= T_DDG*G_S + 2.D0*T_DGG*G_D*G_S + T_GGG*G_D**2*G_S &
      &     + 2.D0*T_DG*G_DS + T_GG*(2.D0*G_D*G_DS+G_DD*G_S) + T_G*G_DDS 
       T_DSS= T_DGG*G_S**2+T_GGG*G_D*G_S**2+T_GG*(2.D0*G_S*G_DS+G_SS*G_D) +T_DG*G_SS+T_G*G_DSS
       T_SSS=T_GGG*G_S**3+T_GG*2.D0*G_S*G_SS+T_GG*G_S*G_SS+T_G*G_SSS
       T_DDN=T_DDN+T_DGN*G_D+T_DGN*G_D+T_GGN*G_D**2+T_GN*G_DD
       T_DSN=T_DGN*G_S+T_GGN*G_D*G_S+T_GN*G_DS
       T_SSN=T_GGN*G_S**2+T_GN*G_SS
       T_NNN=0.D0
       T_DD=T_DD+2.D0*T_DG*G_D+T_GG*G_D**2+T_G*G_DD
       T_DS=T_DG*G_S+T_GG*G_S*G_D+T_G*G_DS
       T_DN=T_DN+T_GN*G_D
       T_SS=T_GG*G_S**2+T_G*G_SS
       T_SN=T_GN*G_S
       T_S=T_G*G_S
       T_D=T_D+T_G*G_D
!
!      ==================================================================
!      == TRANSFORM H0(T(D,S,N),E(D,S),G(D,S)) -> H0(D,S,N)            ==
!      ==      D=RHOT; S=RHOS; N=GRHOT2                                ==
!      ==================================================================
       H0_DDD=      H0_TTT*T_D**3     + 3.D0*H0_TTE*T_D**2*E_D  + 3.D0*H0_TTG*T_D**2*G_D &
      &      + 3.D0*H0_TEE*T_D*E_D**2 + 6.D0*H0_TEG*T_D*E_D*G_D + 3.D0*H0_TGG*T_D*G_D**2 &
      &      +      H0_EEE*E_D**3     + 3.D0*H0_EEG*E_D**2*G_D  + 3.D0*H0_EGG*E_D*G_D**2 &
      &      +      H0_GGG*G_D**3 &
      &      + 3.D0*H0_TT*T_D*T_DD + 3.D0*H0_TE*(T_DD*E_D+T_D*E_DD) &
      &      + 3.D0*H0_EE*E_D*E_DD + 3.D0*H0_TG*(T_DD*G_D+T_D*G_DD) &
      &      + 3.D0*H0_GG*G_D*G_DD + 3.D0*H0_EG*(E_DD*G_D+E_D*G_DD) &
      &      + H0_T*T_DDD + H0_E*E_DDD + H0_G*G_DDD
       H0_SSS=      H0_TTT*T_S**3     + 3.D0*H0_TTE*T_S**2*E_S  + 3.D0*H0_TTG*T_S**2*G_S &
      &      + 3.D0*H0_TEE*T_S*E_S**2 + 6.D0*H0_TEG*T_S*E_S*G_S + 3.D0*H0_TGG*T_S*G_S**2 &
      &      +      H0_EEE*E_S**3     + 3.D0*H0_EEG*E_S**2*G_S  + 3.D0*H0_EGG*E_S*G_S**2 &
      &      +      H0_GGG*G_S**3 &
      &      + 3.D0*H0_TT *T_S*T_SS + 3.D0*H0_TE*(T_SS*E_S+T_S*E_SS) &
      &      + 3.D0*H0_EE *E_S*E_SS + 3.D0*H0_TG*(T_SS*G_S+T_S*G_SS) &
      &      + 3.D0*H0_GG *G_S*G_SS + 3.D0*H0_EG*(E_SS*G_S+E_S*G_SS) &
      &      +      H0_T  *T_SSS    +      H0_E*  E_SSS + H0_G*G_SSS 
       H0_NNN=H0_TTT*T_N**3 + 3.D0*H0_TT*T_N*T_NN + H0_T*T_NNN
       H0_DDS = H0_TTT*T_D**2*T_S                                        &  
              + H0_TTE*(T_D**2*E_S + 2.D0*T_D*E_D*T_S)                   &  
              + H0_TTG*(T_D**2*G_S + 2.D0*T_D*G_D*T_S)                   &
              + H0_TEE*(E_D**2*T_S + 2.D0*T_D*E_D*E_S)                   &
              + H0_TEG*2.D0*(T_D*E_D*G_S + T_D*G_D*E_S +E_D*G_D*T_S )    &
              + H0_TGG*(G_D**2*T_S + 2.D0*T_D*G_D*G_S)                   &
              + H0_EEE*E_D**2*E_S                                        &
              + H0_EEG*(E_D**2*G_S + 2.D0*E_D*G_D*E_S)                   &
              + H0_EGG*(G_D**2*E_S + 2.D0*E_D*G_D*G_S)                   &
              + H0_GGG*G_D**2*G_S                                        &
              + H0_TT*(2.D0*T_D*T_DS + T_DD*T_S)                         &
              + H0_TE*(2.D0*(T_DS*E_D + T_D*E_DS) + T_DD*E_S + E_DD*T_S) &
              + H0_TG*(2.D0*(T_DS*G_D + T_D*G_DS) + T_DD*G_S + G_DD*T_S) &
              + H0_EE*(2.D0*E_D*E_DS + E_DD*E_S)                         &
              + H0_EG*(2.D0*(E_DS*G_D+ E_D*G_DS) + E_DD*G_S + G_DD*E_S)  &
              + H0_GG*(G_DD*G_S +2.D0*G_D*G_DS)                          &
              + H0_T*T_DDS + H0_E*E_DDS + H0_G*G_DDS
       H0_DSS = H0_TTT*T_S**2*T_D                                        &  
              + H0_TTE*(T_S**2*E_D + 2.D0*T_S*E_S*T_D)                   &  
              + H0_TTG*(T_S**2*G_D + 2.D0*T_S*G_S*T_D)                   &
              + H0_TEE*(E_S**2*T_D + 2.D0*T_S*E_S*E_D)                   &
              + H0_TEG*2.D0*(T_S*E_S*G_D + T_S*G_S*E_D +E_S*G_S*T_D )    &
              + H0_TGG*(G_S**2*T_D + 2.D0*T_S*G_S*G_D)                   &
              + H0_EEE*E_S**2*E_D                                        &
              + H0_EEG*(E_S**2*G_D + 2.D0*E_S*G_S*E_D)                   &
              + H0_EGG*(G_S**2*E_D + 2.D0*E_S*G_S*G_D)                   &
              + H0_GGG*G_S**2*G_D                                        &
              + H0_TT*(2.D0*T_S*T_DS + T_SS*T_D)                         &
              + H0_TE*(2.D0*(T_DS*E_S + T_S*E_DS) + T_SS*E_D + E_SS*T_D) &
              + H0_TG*(2.D0*(T_DS*G_S + T_S*G_DS) + T_SS*G_D + G_SS*T_D) &
              + H0_EE*(2.D0*E_S*E_DS + E_SS*E_D)                         &
              + H0_EG*(2.D0*(E_DS*G_S+ E_S*G_DS) + E_SS*G_D + G_SS*E_D)  &
              + H0_GG*(G_SS*G_D +2.D0*G_S*G_DS)                          &
              + H0_T*T_DSS + H0_E*E_DSS + H0_G*G_DSS
       H0_DDN=( H0_TTT*T_D**2 + 2.D0*H0_TTE*T_D*E_D + 2.D0*H0_TTG*T_D*G_D + H0_TT*T_DD &
      &       + H0_TEE*E_D**2 + 2.D0*H0_TEG*E_D*G_D +H0_TE*E_DD &
      &       + H0_TGG*G_D**2 + H0_TG*G_DD)*T_N &
      &       + 2.D0*(H0_TT*T_D+H0_TE*E_D+H0_TG*G_D)*T_DN + H0_T*T_DDN 
       H0_SSN=( H0_TTT*T_S**2 + 2.D0*H0_TTE*T_S*E_S + 2.D0*H0_TTG*T_S*G_S + H0_TT*T_SS &
      &       + H0_TEE*E_S**2 + 2.D0*H0_TEG*E_S*G_S +H0_TE*E_SS &
      &       + H0_TGG*G_S**2 + H0_TG*G_SS)*T_N &
      &       + 2.D0*(H0_TT*T_S+H0_TE*E_S+H0_TG*G_S)*T_SN + H0_T*T_SSN 
       H0_DSN = H0_TTT*T_D*T_S*T_N + H0_TEE*T_N*E_D*E_S + H0_TGG*T_N*G_D*G_S &
              + H0_TTE*(T_S*T_N*E_D+T_D*T_N*E_S) &
              + H0_TTG*(T_S*T_N*G_D + T_D*T_N*G_S) &
              + H0_TEG*(T_N*E_S*G_D +T_N*E_D*G_S) &
              + H0_TG*(T_DN*G_S + T_N*G_DS) &
              + H0_TE*(T_DN*E_S + T_N*E_DS) &
              + H0_TT*(T_DS*T_N + T_S*T_DN) &           
              + H0_TT*T_D*T_SN+ H0_TE*T_SN*E_D+ H0_TG*T_SN*G_D+ H0_T*T_DSN
       H0_DNN =  H0_TTT*T_N**2*T_D + H0_TTE*T_N**2*E_D + H0_TTG*T_N**2*G_D &
               + H0_TT*(2.D0*T_N*T_DN +T_NN*T_D) + H0_TE*T_NN*E_D+ H0_TG*T_NN*G_D
       H0_SNN =  H0_TTT*T_N**2*T_S + H0_TTE*T_N**2*E_S + H0_TTG*T_N**2*G_S &
               + H0_TT*(2.D0*T_N*T_SN +T_NN*T_S) + H0_TE*T_NN*E_S+ H0_TG*T_NN*G_S
       H0_DD = H0_TT*T_D**2 + 2.D0*H0_TE*T_D*E_D + 2.D0*H0_TG*T_D*G_D + H0_T*T_DD &
      &      + H0_EE*E_D**2 + 2.D0*H0_EG*E_D*G_D +H0_E*E_DD &
      &      + H0_GG*G_D**2 + H0_G*G_DD
       H0_SS = H0_TT*T_S**2 + 2.D0*H0_TE*T_S*E_S + 2.D0*H0_TG*T_S*G_S + H0_T*T_SS &
      &      + H0_EE*E_S**2 + 2.D0*H0_EG*E_S*G_S +H0_E*E_SS &
      &      + H0_GG*G_S**2 + H0_G*G_SS
       H0_NN = H0_TT*T_N**2  + H0_T*T_NN
       H0_DS = H0_TT*T_D*T_S + H0_TE*(T_D*E_S+T_S*E_D) + H0_TG*(T_D*G_S+T_S*G_D) +H0_T*T_DS &
      &      + H0_EE*E_D*E_S + H0_EG*(E_D*G_S+E_S*G_D) + H0_E*E_DS &
      &      + H0_GG*G_D*G_S + H0_G*G_DS
       H0_DN = H0_TT*T_D*T_N + H0_TE*T_N*E_D + H0_TG*T_N*G_D +H0_T*T_DN
       H0_SN = H0_TT*T_S*T_N + H0_TE*T_N*E_S + H0_TG*T_N*G_S +H0_T*T_SN
       H0_D  = H0_T*T_D+H0_E*E_D+H0_G*G_D
       H0_S  = H0_T*T_S+H0_E*E_S+H0_G*G_S
       H0_N  = H0_T*T_N
!
!      ==================================================================
!      == EXCHANGE ENERGY EXC=RHOT*(E+H0)                              ==
!      ==================================================================
       EXC      =RHOT*(E+H0)
       DEXC(1)  =RHOT*(E_D+H0_D)+E+H0
       DEXC(2)  =RHOT*(E_S+H0_S)
       DEXC(3)  =RHOT*H0_N
       D2EXC(1,1)=RHOT*(E_DD+H0_DD)+2.D0*(E_D+H0_D)
       D2EXC(1,2)=RHOT*(E_DS+H0_DS)+(E_S+H0_S)
       D2EXC(1,3)=RHOT*H0_DN+H0_N
       D2EXC(2,2)=RHOT*(E_SS+H0_SS)
       D2EXC(2,3)=RHOT*H0_SN
       D2EXC(3,3)=RHOT*H0_NN
       D2EXC(3,3)=RHOT*H0_NN
       D3EXC(1,1,1)=RHOT*(E_DDD+H0_DDD)+3.D0*(E_DD+H0_DD)
       D3EXC(1,1,2)=RHOT*(E_DDS+H0_DDS)+2.D0*(E_DS+H0_DS)
       D3EXC(1,1,3)=RHOT*H0_DDN+2.D0*H0_DN
       D3EXC(1,2,2)=RHOT*(E_DSS+H0_DSS)+(E_SS+H0_SS)
       D3EXC(1,2,3)=RHOT*H0_DSN+H0_SN
       D3EXC(1,3,3)=RHOT*H0_DNN+H0_NN
       D3EXC(2,2,2)=RHOT*(E_SSS+H0_SSS)
       D3EXC(2,2,3)=RHOT*H0_SSN
       D3EXC(2,3,3)=RHOT*H0_SNN
       D3EXC(3,3,3)=RHOT*H0_NNN
!
!      ==================================================================
!      == SYMMETRIZE                                                   ==
!      ==================================================================
       D2EXC(2,1)=D2EXC(1,2)
       D2EXC(3,1)=D2EXC(1,3)
       D2EXC(3,2)=D2EXC(2,3)
       D3EXC(1,2,1)=D3EXC(1,1,2)
       D3EXC(2,1,1)=D3EXC(1,1,2)
       D3EXC(1,3,1)=D3EXC(1,1,3)
       D3EXC(3,1,1)=D3EXC(1,1,3)
       D3EXC(2,1,2)=D3EXC(1,2,2)
       D3EXC(2,2,1)=D3EXC(1,2,2)
       D3EXC(3,1,3)=D3EXC(1,3,3)
       D3EXC(3,3,1)=D3EXC(1,3,3)
       D3EXC(2,3,2)=D3EXC(2,2,3)
       D3EXC(3,2,2)=D3EXC(2,2,3)
       D3EXC(3,2,3)=D3EXC(2,3,3)
       D3EXC(3,3,2)=D3EXC(2,3,3)
       D3EXC(1,3,2)=D3EXC(1,2,3)
       D3EXC(2,1,3)=D3EXC(1,2,3)
       D3EXC(2,3,1)=D3EXC(1,2,3)
       D3EXC(3,1,2)=D3EXC(1,2,3)
       D3EXC(3,2,1)=D3EXC(1,2,3)
       RETURN
       END
!
!.......................................................................
MODULE LYP88_MODULE
!***********************************************************************
!**                                                                   **
!**    IMPLEMENTATION OF THE CORRELATION FUNCTIONAL OF                **
!**    LEE,YANG,PARR (PHYS. REV. B 37,785 (1988-I))                   **
!**    (INCLUDES LOCAL AND NONLOCAL CORRELATION)                      **
!**                                                                   **
!**  FUNCTIONS                                                        **
!**    LYP88$EVAL1                                                    **
!**    LYP88$EVAL2                                                    **
!**    LYP88$EVAL3                                                    **
!**                                                                   **
!**                       MANUEL LOUWERSE, FREE UNIVERSITY (NL) 2003  **
!***********************************************************************
LOGICAL(4)          :: TINI=.FALSE.
REAL(8),PARAMETER   :: A=4.918D-02
REAL(8),PARAMETER   :: B=0.132D0
REAL(8),PARAMETER   :: C=0.2533D0
REAL(8),PARAMETER   :: D=0.349D0
REAL(8)             :: CF,CF2,AB
CONTAINS
!.......................................................................
      SUBROUTINE LYP88_INITIALIZE
      IMPLICIT NONE
      REAL(8)       :: PI
!     ******************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
      PI=4.D0*ATAN(1.D0)
      CF=(3.D0*PI*PI)**(2.D0/3.D0)*3.D0/10.D0
      CF2=CF*2.D0**(11.D0/3.D0)*A*B
      AB=A*B/72.D0
      RETURN
      END SUBROUTINE LYP88_INITIALIZE
END MODULE LYP88_MODULE
!     ..................................................................
      SUBROUTINE LYP88$EVAL1(VAL,EXC,DEXC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE LEE-YANG-PARR CORRELATION                         **
!     **  (PHYS. REV. B 37, 785 (1988-I))                             **
!     **                                                              **
!     **  FORMULA (2) FROM CHEM. PHYS. LETT. 157, 200 (1989)          **
!     **  IS USED FOR THE IMPLEMENTATION.                             **
!     **  FOR THE GGA-PART THIS FORMULA IS WORKED AROUND TO           **
!     **  SOMETHING ONLY CONTAINING RHOT AND RHOS.                    **
!     **                                                              **
!     ******************************************************************
      USE LYP88_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: VAL(5)   !(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)  ! DEXC/DVAL(I)
      REAL(8)             :: RHOT,RHOS,GRHOT2,GRHOS2,GRHOST,RA,RB
      REAL(8)             :: RM13,DR,DDR,DR_T,RHO2,RHOS2,RHOS3
      REAL(8)             :: OMEGA,OMEGA_T,DELTA,DELTA_T
      REAL(8)             :: RA83,RB83
      REAL(8)             :: LDA1,LDA1_T,LDA1_S,LDA2,LDA2_T,LDA2_S
      REAL(8)             :: GGA1,GGA1_T,GGA1_S,GGA1_GT
      REAL(8)             :: GGA2,GGA2_T,GGA2_S,GGA2_GS
      REAL(8)             :: GGA3,GGA3_T,GGA3_S,GGA3_GST
!     ******************************************************************
      IF(.NOT.TINI) CALL LYP88_INITIALIZE
      EXC=0.D0
      DEXC(:)=0.D0
      RHOT=VAL(1)
      RHOS=VAL(2)
      GRHOT2=VAL(3)
      GRHOS2=VAL(4)
      GRHOST=VAL(5)
      RA=(RHOT+RHOS)/2.D0
      RB=(RHOT-RHOS)/2.D0
! N.B.: DRA/DRT=0.5  DRA/DRS=0.5  DRB/DRT=0.5  DRB/DRS=-0.5

      RM13=RHOT**(-1.D0/3.D0)
      DR=1.D0/(1.D0+D*RM13)
      DDR=D*RM13*DR
      DR_T=DDR/3.D0/RHOT    ! *DR
      DELTA=C*RM13+DDR
      DELTA_T=(DDR*DDR-DELTA)/3.D0/RHOT
      OMEGA=EXP(-C*RM13)*DR*RM13**11
      OMEGA_T=(DELTA-11.D0)/3.D0/RHOT    ! *OMEGA

      RA83=RA**(8.D0/3.D0)
      RB83=RB**(8.D0/3.D0)
      LDA1=-A*4.D0*DR*RA*RB/RHOT
      LDA1_T=LDA1*DR_T-LDA1/RHOT-2.D0*A*DR
      LDA1_S=-A*2.D0*DR*(RB-RA)/RHOT
      LDA2=-CF2*OMEGA*RA*RB*(RA83+RB83)
      LDA2_T=LDA2*OMEGA_T-CF2*OMEGA*(11.D0/3.D0*(RA83*RB+RB83*RA)+RB83*RB+RA83*RA)/2.D0
      LDA2_S=-CF2*OMEGA*(11.D0/3.D0*(RA83*RB-RB83*RA)+RB83*RB-RA83*RA)/2.D0

      RHO2=RHOT*RHOT
      RHOS2=RHOS*RHOS
      RHOS3=RHOS2*RHOS
      GGA1_GT=AB*OMEGA*(7.D0*(RHO2-RHOS2)*DELTA+3.D0*RHO2+39.D0*RHOS2)
      GGA1=GGA1_GT*GRHOT2
      GGA1_T=GGA1*OMEGA_T+AB*OMEGA*GRHOT2*(7.D0*(RHO2-RHOS2)*DELTA_T+14.D0*RHOT*DELTA+6.D0*RHOT)
      GGA1_S=AB*OMEGA*GRHOT2*(-14.D0*RHOS*DELTA+78.D0*RHOS)
      GGA2_GS=AB*OMEGA*(2.D0*RHO2-8.D0*RHOS2)
      GGA2=GGA2_GS*GRHOS2
      GGA2_T=GGA2*OMEGA_T+AB*OMEGA*GRHOS2*4.D0*RHOT
      GGA2_S=-AB*OMEGA*GRHOS2*16.D0*RHOS
      GGA3_GST=AB*OMEGA*(RHOT*RHOS*(DELTA-47.D0)-RHOS3/RHOT*(DELTA-11.D0))
      GGA3=GGA3_GST*GRHOST
      GGA3_T=GGA3*OMEGA_T+AB*OMEGA*GRHOST*(RHOS*(DELTA-47.D0)+RHOT*RHOS*DELTA_T &
            +RHOS3/RHO2*(DELTA-11.D0)-RHOS3/RHOT*DELTA_T)
      GGA3_S=AB*OMEGA*GRHOST*(RHOT*(DELTA-47.D0)-3.D0*RHOS2/RHOT*(DELTA-11.D0))

      EXC=LDA1+LDA2+GGA1+GGA2+GGA3
      DEXC(1)=LDA1_T+LDA2_T+GGA1_T+GGA2_T+GGA3_T
      DEXC(2)=LDA1_S+LDA2_S+GGA1_S+GGA2_S+GGA3_S
      DEXC(3)=GGA1_GT
      DEXC(4)=GGA2_GS
      DEXC(5)=GGA3_GST

      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LYP88$EVAL2(VAL,EXC,DEXC,D2EXC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE LEE-YANG-PARR CORRELATION                         **
!     **  (PHYS. REV. B 37, 785 (1988-I))                             **
!     **                                                              **
!     **  FORMULA (2) FROM CHEM. PHYS. LETT. 157, 200 (1989)          **
!     **  IS USED FOR THE IMPLEMENTATION                              **
!     **  FOR THE GGA-PART THIS FORMULA IS WORKED AROUND TO           **
!     **  SOMETHING ONLY CONTAINING RHOT AND RHOS.                    **
!     **                                                              **
!     ******************************************************************
      USE LYP88_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: VAL(5)       !(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)      ! DEXC/DVAL(I)
      REAL(8),INTENT(OUT) :: D2EXC(5,5)   ! D2EXC/DVAL(I)DVAL(J)
      REAL(8)             :: RHOT,RHOS,GRHOT2,GRHOS2,GRHOST,RA,RB
      REAL(8)             :: RM13,DR,DDR,DR_T,DR_TT,DDR2,RHO2,RHOS2,RS
      REAL(8)             :: OMEGA,OMEGA_T,DELTA,DELTA_T
      REAL(8)             :: OMEGA_TT,DELTA_TT
      REAL(8)             :: RA53,RA83,RB53,RB83,TEMP,TEMP2
      REAL(8)             :: LDA1,LDA1_T,LDA1_S,LDA2,LDA2_T,LDA2_S
      REAL(8)             :: LDA1_TT,LDA1_TS,LDA1_SS,LDA2_TT,LDA2_TS,LDA2_SS
      REAL(8)             :: GGA1,GGA1_T,GGA1_S,GGA1_G
      REAL(8)             :: GGA1_TT,GGA1_TS,GGA1_SS,GGA1_TG,GGA1_SG
      REAL(8)             :: GGA2,GGA2_T,GGA2_S,GGA2_G
      REAL(8)             :: GGA2_TT,GGA2_TS,GGA2_SS,GGA2_TG,GGA2_SG
      REAL(8)             :: GGA3,GGA3_T,GGA3_S,GGA3_G
      REAL(8)             :: GGA3_TT,GGA3_TS,GGA3_SS,GGA3_TG,GGA3_SG
      REAL(8)             :: GGA1_SSG,GGA2_SSG
      INTEGER(4)          :: I
!     ******************************************************************
      IF(.NOT.TINI) CALL LYP88_INITIALIZE
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      RHOT=VAL(1)
      RHOS=VAL(2)
      GRHOT2=VAL(3)
      GRHOS2=VAL(4)
      GRHOST=VAL(5)
      RA=(RHOT+RHOS)/2.D0
      RB=(RHOT-RHOS)/2.D0
! N.B.: DRA/DRT=0.5  DRA/DRS=0.5  DRB/DRT=0.5  DRB/DRS=-0.5

      RM13=RHOT**(-1.D0/3.D0)
      DR=1.D0/(1.D0+D*RM13)
      DDR=D*RM13*DR
      DDR2=DDR*DDR
      RHO2=RHOT*RHOT
      DR_T=DDR/3.D0/RHOT      ! *DR
      DR_TT=DR_T*(2.D0*DDR-4.D0)/3.D0/RHOT                   ! *DR
      DELTA=C*RM13+DDR
      DELTA_T=(DDR2-DELTA)/3.D0/RHOT
      DELTA_TT=(2.D0*DDR2*DDR-6.D0*DDR2+4.D0*DELTA)/9.D0/RHO2
      OMEGA=EXP(-C*RM13)*DR*RM13**11
      OMEGA_T=(DELTA-11.D0)/3.D0/RHOT    ! *OMEGA
      OMEGA_TT=(OMEGA_T*(DELTA-14.D0)+DELTA_T)/3.D0/RHOT   ! *OMEGA

!     ===================================================================
!     ==  LDA-PART                                                     ==
!     ===================================================================
      LDA1=-A*4.D0*DR*RA*RB/RHOT
      LDA1_T=LDA1*(DR_T-1.D0/RHOT)-2.D0*A*DR
      LDA1_S=-A*2.D0*DR*(RB-RA)/RHOT
      LDA1_TT=LDA1_T*(DR_T-1.D0/RHOT)+LDA1*(DR_TT-DR_T*DR_T+1.D0/RHO2)-2.D0*A*DR_T*DR
      LDA1_TS=LDA1_S*(DR_T-1.D0/RHOT)
      LDA1_SS=A*2.D0*DR/RHOT
      RA53=RA**(5.D0/3.D0)
      RA83=RA53*RA
      RB53=RB**(5.D0/3.D0)
      RB83=RB53*RB
      TEMP=OMEGA*CF2*(11.D0/3.D0*(RA83*RB+RB83*RA)+RB83*RB+RA83*RA)/2.D0
      TEMP2=OMEGA*CF2*(44.D0*(RA53*RB+RB53*RA)+33.D0*(RA83+RB83))/18.D0
      LDA2=-CF2*OMEGA*RA*RB*(RA83+RB83)
      LDA2_T=LDA2*OMEGA_T-TEMP
      LDA2_S=-CF2*OMEGA*(11.D0/3.D0*(RA83*RB-RB83*RA)+RB83*RB-RA83*RA)/2.D0
      LDA2_TT=LDA2*OMEGA_TT-2.D0*OMEGA_T*TEMP-TEMP2
      LDA2_TS=LDA2_S*OMEGA_T-CF2*OMEGA*22.D0/9.D0*(RA53*RB-RB53*RA)
      LDA2_SS=-CF2*OMEGA*(44.D0*(RA53*RB+RB53*RA)-33.D0*(RA83+RB83))/18.D0

!     ===================================================================
!     ==  GGA-PART                                                     ==
!     ===================================================================
      RHOS2=RHOS*RHOS
      TEMP=AB*OMEGA*(7.D0*(RHO2-RHOS2)*DELTA_T+14.D0*RHOT*DELTA+6.D0*RHOT)
      TEMP2=AB*OMEGA*(7.D0*(RHO2-RHOS2)*DELTA_TT+28.D0*RHOT*DELTA_T+14.D0*DELTA+6.D0)
      GGA1_G=AB*OMEGA*(7.D0*(RHO2-RHOS2)*DELTA+3.D0*RHO2+39.D0*RHOS2)
      GGA1=GGA1_G*GRHOT2
      GGA1_TG=GGA1_G*OMEGA_T+TEMP
      GGA1_T=GGA1_TG*GRHOT2
      GGA1_SSG=AB*OMEGA*(-14.D0*DELTA+78.D0)
      GGA1_S=GGA1_SSG*GRHOT2*RHOS
      GGA1_SS=GGA1_SSG*GRHOT2
      GGA1_SG=GGA1_SSG*RHOS
      GGA1_TT=GGA1*OMEGA_TT+2.D0*OMEGA_T*TEMP*GRHOT2+TEMP2*GRHOT2
      GGA1_TS=GGA1_S*OMEGA_T-AB*OMEGA*GRHOT2*14.D0*RHOS*DELTA_T

      GGA2_G=AB*OMEGA*(2.D0*RHO2-8.D0*RHOS2)
      GGA2=GGA2_G*GRHOS2
      GGA2_TG=GGA2_G*OMEGA_T+AB*OMEGA*4.D0*RHOT
      GGA2_T=GGA2_TG*GRHOS2
      GGA2_SSG=-AB*OMEGA*16.D0
      GGA2_S=GGA2_SSG*GRHOS2*RHOS
      GGA2_SS=GGA2_SSG*GRHOS2
      GGA2_SG=GGA2_SSG*RHOS
      GGA2_TT=GGA2*OMEGA_TT+AB*OMEGA*GRHOS2*(OMEGA_T*8.D0*RHOT+4.D0)
      GGA2_TS=GGA2_S*OMEGA_T

      RS=RHOS2/RHO2
      TEMP=AB*OMEGA*(RHOS*(DELTA-47.D0+RS*(DELTA-11.D0))+RHOS*RHOT*DELTA_T*(1.D0-RS))
      TEMP2=AB*OMEGA*(2.D0*RHOS*DELTA_T*(1.D0+RS)+RHOS*RHOT*DELTA_TT*(1.D0-RS) &
           -2.D0*RHOS/RHOT*RS*(DELTA-11.D0))
      GGA3_G=AB*OMEGA*RHOS*RHOT*(DELTA-47.D0-RS*(DELTA-11.D0))
      GGA3=GGA3_G*GRHOST
      GGA3_TG=GGA3_G*OMEGA_T+TEMP
      GGA3_T=GGA3_TG*GRHOST
      GGA3_SG=AB*OMEGA*(RHOT*(DELTA-47.D0)-3.D0*RS*RHOT*(DELTA-11.D0))
      GGA3_S=GGA3_SG*GRHOST
      GGA3_TT=GGA3*OMEGA_TT+2.D0*OMEGA_T*TEMP*GRHOST+TEMP2*GRHOST
      GGA3_TS=GGA3_S*OMEGA_T+AB*OMEGA*GRHOST*(DELTA-47.D0+RHOT*DELTA_T+3.D0*RS*(DELTA-11.D0 &
             -RHOT*DELTA_T))
      GGA3_SS=-AB*OMEGA*GRHOST*6.D0*RHOS/RHOT*(DELTA-11.D0)

!     ===================================================================
!     ==  EXC=LDA+GGA                                                  ==
!     ===================================================================
      EXC=LDA1+LDA2+GGA1+GGA2+GGA3
      DEXC(1)=LDA1_T+LDA2_T+GGA1_T+GGA2_T+GGA3_T
      DEXC(2)=LDA1_S+LDA2_S+GGA1_S+GGA2_S+GGA3_S
      DEXC(3)=GGA1_G
      DEXC(4)=GGA2_G
      DEXC(5)=GGA3_G
      D2EXC(1,1)=LDA1_TT+LDA2_TT+GGA1_TT+GGA2_TT+GGA3_TT
      D2EXC(1,2)=LDA1_TS+LDA2_TS+GGA1_TS+GGA2_TS+GGA3_TS
      D2EXC(1,3)=GGA1_TG
      D2EXC(1,4)=GGA2_TG
      D2EXC(1,5)=GGA3_TG
      D2EXC(2,2)=LDA1_SS+LDA2_SS+GGA1_SS+GGA2_SS+GGA3_SS
      D2EXC(2,3)=GGA1_SG
      D2EXC(2,4)=GGA2_SG
      D2EXC(2,5)=GGA3_SG

!     ============================================================
!     == SYMMETRIZE                                             ==
!     ============================================================
      D2EXC(2,1)=D2EXC(1,2)
      DO I=3,5
        D2EXC(I,1)=D2EXC(1,I)
        D2EXC(I,2)=D2EXC(2,I)
      ENDDO

      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LYP88$EVAL3(VAL,EXC,DEXC,D2EXC,D3EXC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE LEE-YANG-PARR CORRELATION                         **
!     **  (PHYS. REV. B 37, 785 (1988-I))                             **
!     **                                                              **
!     **  FORMULA (2) FROM CHEM. PHYS. LETT. 157, 200 (1989)          **
!     **  IS USED FOR THE IMPLEMENTATION                              **
!     **  FOR THE GGA-PART THIS FORMULA IS WORKED AROUND TO           **
!     **  SOMETHING ONLY CONTAINING RHOT AND RHOS.                    **
!     **                                                              **
!     ******************************************************************
      USE LYP88_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: VAL(5)       !(RHOT,RHOS,GRHOT2,GRHOS2,GRHOST)
      REAL(8),INTENT(OUT) :: EXC
      REAL(8),INTENT(OUT) :: DEXC(5)      ! DEXC/DVAL(I)
      REAL(8),INTENT(OUT) :: D2EXC(5,5)   ! D2EXC/DVAL(I)DVAL(J)
      REAL(8),INTENT(OUT) :: D3EXC(5,5,5) ! D3EXC/DVAL(I)DVAL(J)DVAL(K)
      REAL(8)             :: RHOT,RHOS,GRHOT2,GRHOS2,GRHOST,RA,RB
      REAL(8)             :: RM13,DR,DDR,DR_T,DR_TT,DR_TTT,DDR2,RHO2,DR_T2,RHOS2,RS
      REAL(8)             :: OMEGA,OMEGA_T,DELTA,DELTA_T
      REAL(8)             :: OMEGA_TT,OMEGA_TTT,DELTA_TT,DELTA_TTT
      REAL(8)             :: RA23,RA53,RA83,RB23,RB53,RB83,TEMP,TEMP2
      REAL(8)             :: LDA1,LDA1_T,LDA1_S,LDA2,LDA2_T,LDA2_S
      REAL(8)             :: LDA1_TT,LDA1_TS,LDA1_SS,LDA2_TT,LDA2_TS,LDA2_SS
      REAL(8)             :: LDA1_TTT,LDA1_TTS,LDA1_TSS,LDA1_SSS
      REAL(8)             :: LDA2_TTT,LDA2_TTS,LDA2_TSS,LDA2_SSS
      REAL(8)             :: GGA1,GGA1_T,GGA1_S,GGA1_G
      REAL(8)             :: GGA1_TT,GGA1_TS,GGA1_SS,GGA1_TG,GGA1_SG
      REAL(8)             :: GGA1_TTT,GGA1_TTS,GGA1_TSS,GGA1_SSS,GGA1_TTG,GGA1_TSG,GGA1_SSG
      REAL(8)             :: GGA2,GGA2_T,GGA2_S,GGA2_G
      REAL(8)             :: GGA2_TT,GGA2_TS,GGA2_SS,GGA2_TG,GGA2_SG
      REAL(8)             :: GGA2_TTT,GGA2_TTS,GGA2_TSS,GGA2_SSS,GGA2_TTG,GGA2_TSG,GGA2_SSG
      REAL(8)             :: GGA3,GGA3_T,GGA3_S,GGA3_G
      REAL(8)             :: GGA3_TT,GGA3_TS,GGA3_SS,GGA3_TG,GGA3_SG
      REAL(8)             :: GGA3_TTT,GGA3_TTS,GGA3_TSS,GGA3_SSS,GGA3_TTG,GGA3_TSG,GGA3_SSG
      REAL(8)             :: GGA3_SSSG
      INTEGER(4)          :: I
!     ******************************************************************
      IF(.NOT.TINI) CALL LYP88_INITIALIZE
      EXC=0.D0
      DEXC(:)=0.D0
      D2EXC(:,:)=0.D0
      D3EXC(:,:,:)=0.D0
      RHOT=VAL(1)
      RHOS=VAL(2)
      GRHOT2=VAL(3)
      GRHOS2=VAL(4)
      GRHOST=VAL(5)
      RA=(RHOT+RHOS)/2.D0
      RB=(RHOT-RHOS)/2.D0
! N.B.: DRA/DRT=0.5  DRA/DRS=0.5  DRB/DRT=0.5  DRB/DRS=-0.5

      RM13=RHOT**(-1.D0/3.D0)
      DR=1.D0/(1.D0+D*RM13)
      DDR=D*RM13*DR
      DDR2=DDR*DDR
      RHO2=RHOT*RHOT
      DR_T=DDR/3.D0/RHOT      ! *DR
      DR_TT=DR_T*(2.D0*DDR-4.D0)/3.D0/RHOT                   ! *DR
      DR_TTT=DR_T*(6.D0*DDR2-24.D0*DDR+28.D0)/9.D0/RHO2      ! *DR
      DR_T2=DR_T*DR_T
      DELTA=C*RM13+DDR
      DELTA_T=(DDR2-DELTA)/3.D0/RHOT
      DELTA_TT=(2.D0*DDR2*DDR-6.D0*DDR2+4.D0*DELTA)/9.D0/RHO2
      DELTA_TTT=(6.D0*DDR2*DDR2-30.D0*DDR2*DDR+52.D0*DDR2-28.D0*DELTA)/27.D0/RHO2/RHOT
      OMEGA=EXP(-C*RM13)*DR*RM13**11
      OMEGA_T=(DELTA-11.D0)/3.D0/RHOT    ! *OMEGA
      OMEGA_TT=(OMEGA_T*(DELTA-14.D0)+DELTA_T)/3.D0/RHOT   ! *OMEGA
      OMEGA_TTT=(OMEGA_TT*(DELTA-17.D0)+2.D0*OMEGA_T*DELTA_T+DELTA_TT)/3.D0/RHOT   ! *OMEGA

!     ===================================================================
!     ==  LDA-PART                                                     ==
!     ===================================================================
      LDA1=-A*4.D0*DR*RA*RB/RHOT
      LDA1_T=LDA1*(DR_T-1.D0/RHOT)-2.D0*A*DR
      LDA1_S=-A*2.D0*DR*(RB-RA)/RHOT
      LDA1_TT=LDA1_T*(DR_T-1.D0/RHOT)+LDA1*(DR_TT-DR_T2+1.D0/RHO2)-2.D0*A*DR_T*DR
      LDA1_TS=LDA1_S*(DR_T-1.D0/RHOT)
      LDA1_SS=A*2.D0*DR/RHOT
      LDA1_TTT=LDA1_TT*(DR_T-1.D0/RHOT)+LDA1_T*2.D0*(DR_TT-DR_T2+1.D0/RHO2) &
              +LDA1*(DR_TTT-3.D0*DR_TT*DR_T+2.D0*DR_T2*DR_T-2.D0/RHO2/RHOT)-2.D0*A*DR_TT*DR
      LDA1_TTS=LDA1_TS*(DR_T-1.D0/RHOT)+LDA1_S*(DR_TT-DR_T2+1.D0/RHO2)
      LDA1_TSS=LDA1_SS*(DR_T-1.D0/RHOT)
      LDA1_SSS=0.D0
      RA23=RA**(2.D0/3.D0)
      RA53=RA23*RA
      RA83=RA53*RA
      RB23=RB**(2.D0/3.D0)
      RB53=RB23*RB
      RB83=RB53*RB
      TEMP=OMEGA*CF2*(11.D0/3.D0*(RA83*RB+RB83*RA)+RB83*RB+RA83*RA)/2.D0
      TEMP2=OMEGA*CF2*(44.D0*(RA53*RB+RB53*RA)+33.D0*(RA83+RB83))/18.D0
      LDA2=-CF2*OMEGA*RA*RB*(RA83+RB83)
      LDA2_T=LDA2*OMEGA_T-TEMP
      LDA2_S=-CF2*OMEGA*(11.D0/3.D0*(RA83*RB-RB83*RA)+RB83*RB-RA83*RA)/2.D0
      LDA2_TT=LDA2*OMEGA_TT-2.D0*OMEGA_T*TEMP-TEMP2
      LDA2_TS=LDA2_S*OMEGA_T-CF2*OMEGA*22.D0/9.D0*(RA53*RB-RB53*RA)
      LDA2_SS=-CF2*OMEGA*(44.D0*(RA53*RB+RB53*RA)-33.D0*(RA83+RB83))/18.D0
      LDA2_TTT=LDA2*OMEGA_TTT-3.D0*OMEGA_TT*TEMP-3.D0*OMEGA_T*TEMP2 &
              -CF2*OMEGA*(55.D0*(RA23*RB+RB23*RA)+99.D0*(RA53+RB53))/27.D0
      LDA2_TTS=LDA2_S*OMEGA_TT-CF2*OMEGA*(132.D0*OMEGA_T*(RA53*RB-RB53*RA) &
              +55.D0*(RA23*RB-RB23*RA)+33.D0*(RA53-RB53))/27.D0
      LDA2_TSS=LDA2_SS*OMEGA_T-CF2*OMEGA*(55.D0*(RA23*RB+RB23*RA)-33.D0*(RA53+RB53))/27.D0
      LDA2_SSS=-CF2*OMEGA*(55.D0*(RA23*RB-RB23*RA)+99.D0*(RB53-RA53))/27.D0

!     ===================================================================
!     ==  GGA-PART                                                     ==
!     ===================================================================
      RHOS2=RHOS*RHOS
      TEMP=AB*OMEGA*(7.D0*(RHO2-RHOS2)*DELTA_T+14.D0*RHOT*DELTA+6.D0*RHOT)
      TEMP2=AB*OMEGA*(7.D0*(RHO2-RHOS2)*DELTA_TT+28.D0*RHOT*DELTA_T+14.D0*DELTA+6.D0)
      GGA1_G=AB*OMEGA*(7.D0*(RHO2-RHOS2)*DELTA+3.D0*RHO2+39.D0*RHOS2)
      GGA1=GGA1_G*GRHOT2
      GGA1_TG=GGA1_G*OMEGA_T+TEMP
      GGA1_T=GGA1_TG*GRHOT2
      GGA1_SSG=AB*OMEGA*(-14.D0*DELTA+78.D0)
      GGA1_S=GGA1_SSG*GRHOT2*RHOS
      GGA1_SS=GGA1_SSG*GRHOT2
      GGA1_SG=GGA1_SSG*RHOS
      GGA1_TTG=GGA1_G*OMEGA_TT+2.D0*OMEGA_T*TEMP+TEMP2
      GGA1_TT=GGA1_TTG*GRHOT2
      GGA1_TSG=GGA1_SG*OMEGA_T-AB*OMEGA*14.D0*RHOS*DELTA_T
      GGA1_TS=GGA1_TSG*GRHOT2
      GGA1_TTT=GGA1*OMEGA_TTT+3.D0*OMEGA_TT*TEMP*GRHOT2+3.D0*OMEGA_T*TEMP2*GRHOT2 &
              +AB*OMEGA*GRHOT2*(7.D0*(RHO2-RHOS2)*DELTA_TTT+42.D0*RHOT*DELTA_TT+42.D0*DELTA_T)
      GGA1_TTS=GGA1_S*OMEGA_TT-AB*OMEGA*GRHOT2*(2.D0*OMEGA_T*DELTA_T+DELTA_TT)*14.D0*RHOS
      GGA1_TSS=GGA1_SS*OMEGA_T-AB*OMEGA*GRHOT2*14.D0*DELTA_T
      GGA1_SSS=0.D0

      GGA2_G=AB*OMEGA*(2.D0*RHO2-8.D0*RHOS2)
      GGA2=GGA2_G*GRHOS2
      GGA2_TG=GGA2_G*OMEGA_T+AB*OMEGA*4.D0*RHOT
      GGA2_T=GGA2_TG*GRHOS2
      GGA2_SSG=-AB*OMEGA*16.D0
      GGA2_S=GGA2_SSG*GRHOS2*RHOS
      GGA2_SS=GGA2_SSG*GRHOS2
      GGA2_SG=GGA2_SSG*RHOS
      GGA2_TTG=GGA2_G*OMEGA_TT+AB*OMEGA*(OMEGA_T*8.D0*RHOT+4.D0)
      GGA2_TT=GGA2_TTG*GRHOS2
      GGA2_TS=GGA2_S*OMEGA_T
      GGA2_TSG=GGA2_SG*OMEGA_T
      GGA2_TTT=GGA2*OMEGA_TTT+AB*OMEGA*GRHOS2*12.D0*(OMEGA_TT*RHOT+OMEGA_T)
      GGA2_TTS=GGA2_S*OMEGA_TT
      GGA2_TSS=GGA2_SS*OMEGA_T
      GGA2_SSS=0.D0

      RS=RHOS2/RHO2
      TEMP=AB*OMEGA*(RHOS*(DELTA-47.D0+RS*(DELTA-11.D0))+RHOS*RHOT*DELTA_T*(1.D0-RS))
      TEMP2=AB*OMEGA*(2.D0*RHOS*DELTA_T*(1.D0+RS)+RHOS*RHOT*DELTA_TT*(1.D0-RS) &
           -2.D0*RHOS/RHOT*RS*(DELTA-11.D0))
      GGA3_G=AB*OMEGA*RHOS*RHOT*(DELTA-47.D0-RS*(DELTA-11.D0))
      GGA3=GGA3_G*GRHOST
      GGA3_TG=GGA3_G*OMEGA_T+TEMP
      GGA3_T=GGA3_TG*GRHOST
      GGA3_SG=AB*OMEGA*(RHOT*(DELTA-47.D0)-3.D0*RS*RHOT*(DELTA-11.D0))
      GGA3_S=GGA3_SG*GRHOST
      GGA3_TTG=GGA3_G*OMEGA_TT+2.D0*OMEGA_T*TEMP+TEMP2
      GGA3_TT=GGA3_TTG*GRHOST
      GGA3_TSG=GGA3_SG*OMEGA_T+AB*OMEGA*(DELTA-47.D0+RHOT*DELTA_T+3.D0*RS*(DELTA-11.D0 &
             -RHOT*DELTA_T))
      GGA3_TS=GGA3_TSG*GRHOST
      GGA3_SSSG=-AB*OMEGA*6.D0/RHOT*(DELTA-11.D0)
      GGA3_SS=GGA3_SSSG*GRHOST*RHOS
      GGA3_SSS=GGA3_SSSG*GRHOST
      GGA3_SSG=GGA3_SSSG*RHOS
      GGA3_TTT=GGA3*OMEGA_TTT+3.D0*OMEGA_TT*TEMP*GRHOST+3.D0*OMEGA_T*TEMP2*GRHOST &
              +AB*OMEGA*GRHOST*(3.D0*RHOS*DELTA_TT*(1.D0+RS)+RHOS*RHOT*DELTA_TTT*(1.D0-RS) &
              +6.D0*RHOS/RHOT*RS*((DELTA-11.D0)/RHOT-DELTA_T))
      GGA3_TTS=GGA3_S*OMEGA_TT+AB*OMEGA*GRHOST*(6.D0*RS*(DELTA-11.D0-RHOT*DELTA_T) &
              *(OMEGA_T-1.D0/RHOT)+2.D0*OMEGA_T*(DELTA+RHOT*DELTA_T-47.D0)+2.D0*DELTA_T &
              +RHOT*DELTA_TT*(1.D0-3.D0*RS))
      GGA3_TSS=GGA3_SS*OMEGA_T+AB*OMEGA*GRHOST*6.D0*RHOS/RHO2*(DELTA-11.D0-RHOT*DELTA_T)

!     ===================================================================
!     ==  EXC=LDA+GGA                                                  ==
!     ===================================================================
      EXC=LDA1+LDA2+GGA1+GGA2+GGA3
      DEXC(1)=LDA1_T+LDA2_T+GGA1_T+GGA2_T+GGA3_T
      DEXC(2)=LDA1_S+LDA2_S+GGA1_S+GGA2_S+GGA3_S
      DEXC(3)=GGA1_G
      DEXC(4)=GGA2_G
      DEXC(5)=GGA3_G
      D2EXC(1,1)=LDA1_TT+LDA2_TT+GGA1_TT+GGA2_TT+GGA3_TT
      D2EXC(1,2)=LDA1_TS+LDA2_TS+GGA1_TS+GGA2_TS+GGA3_TS
      D2EXC(1,3)=GGA1_TG
      D2EXC(1,4)=GGA2_TG
      D2EXC(1,5)=GGA3_TG
      D2EXC(2,2)=LDA1_SS+LDA2_SS+GGA1_SS+GGA2_SS+GGA3_SS
      D2EXC(2,3)=GGA1_SG
      D2EXC(2,4)=GGA2_SG
      D2EXC(2,5)=GGA3_SG
      D3EXC(1,1,1)=LDA1_TTT+LDA2_TTT+GGA1_TTT+GGA2_TTT+GGA3_TTT
      D3EXC(1,1,2)=LDA1_TTS+LDA2_TTS+GGA1_TTS+GGA2_TTS+GGA3_TTS
      D3EXC(1,1,3)=GGA1_TTG
      D3EXC(1,1,4)=GGA2_TTG
      D3EXC(1,1,5)=GGA3_TTG
      D3EXC(1,2,2)=LDA1_TSS+LDA2_TSS+GGA1_TSS+GGA2_TSS+GGA3_TSS
      D3EXC(1,2,3)=GGA1_TSG
      D3EXC(1,2,4)=GGA2_TSG
      D3EXC(1,2,5)=GGA3_TSG
      D3EXC(2,2,2)=LDA1_SSS+LDA2_SSS+GGA1_SSS+GGA2_SSS+GGA3_SSS
      D3EXC(2,2,3)=GGA1_SSG
      D3EXC(2,2,4)=GGA2_SSG
      D3EXC(2,2,5)=GGA3_SSG

!     ============================================================
!     == SYMMETRIZE                                             ==
!     ============================================================
      D2EXC(2,1)=D2EXC(1,2)
      D3EXC(1,2,1)=D3EXC(1,1,2)
      D3EXC(2,1,1)=D3EXC(1,1,2)
      D3EXC(2,1,2)=D3EXC(1,2,2)
      D3EXC(2,2,1)=D3EXC(1,2,2)
      DO I=3,5
        D2EXC(I,1)=D2EXC(1,I)
        D2EXC(I,2)=D2EXC(2,I)
        D3EXC(1,I,1)=D3EXC(1,1,I)
        D3EXC(I,1,1)=D3EXC(1,1,I)
        D3EXC(2,I,2)=D3EXC(2,2,I)
        D3EXC(I,2,2)=D3EXC(2,2,I)
        D3EXC(1,I,2)=D3EXC(1,2,I)
        D3EXC(2,1,I)=D3EXC(1,2,I)
        D3EXC(2,I,1)=D3EXC(1,2,I)
        D3EXC(I,1,2)=D3EXC(1,2,I)
        D3EXC(I,2,1)=D3EXC(1,2,I)
      ENDDO

      RETURN
      END


