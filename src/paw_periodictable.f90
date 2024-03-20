!PROGRAM MAIN;USE PERIODICTABLE_MODULE;CALL PERIODICTABLE$REPORT(6);STOP;END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE PERIODICTABLE_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: PERIODICTABLE                                                      **
!**                                                                           **
!**  PURPOSE: MAINTAINS A PERIODIC TABLE OF CHEMICAL ELEMENTS                 **
!**                                                                           **
!**  METHOD:                                                                  **
!**    PERIODICTABLE$GET(SYMBOL,ID,VALUE)                                     **
!**    PERIODICTABLE$GET(IZ,ID,VALUE)                                         **
!**    PERIODICTABLE$GET(Z,ID,VALUE)                                          **
!**                                                                           **
!**      SYMBOL IS THE TWO-LETTER ELEMENT SYMBOL                              **
!**      IZ IS THE ATOMIC NUMBER (INTEGER(4))                                 **
!**      Z IS THE ATOMIC NUMBER (REAL(8))                                     **
!**      ID IS THE IDENTIFIER OF THE DATE REQUESTED. IT MAY HAVE THE VALUES   **
!**         'Z'         ATOMIC NUMBER                                         **
!**         'SYMBOL'    ELEMENT SYMBOL                                        **
!**         'MASS'      MASS IN ATOMIC UNITS (ELECTRON MASSES)                **
!**         'R(COV)'    COVALENT RADIUS IN A0                                 **
!**         'R(ASA)'    1.10534 TIMES THE COVALENT RADIUS IN A0               **
!**                     (FACTOR CORRESPONDS TO THE RATIO BETWEEN TOUCHING AND **
!**                     VOLUME FILLING SPHERES IN AN FCC LATTICE)             **
!**         'R(VDW)'    VAN DER WAALS RADIUS                                  **
!**         'RNUC'      NUCLEAR RADIUS                                        **
!**         'EN'        PAULI ELECTRONEGATIVITY                               **
!**         'OCC(S)'    ATOMIC OCCUPATION OF THE VALENCE S-SHELL              **
!**         'OCC(P)'    ATOMIC OCCUPATION OF THE VALENCE P-SHELL              **
!**         'OCC(D)'    ATOMIC OCCUPATION OF THE VALENCE D-SHELL              **
!**         'OCC(F)'    ATOMIC OCCUPATION OF THE VALENCE F-SHELL              **
!**         '#NODES(S)' NUMBER OF NODES OF THE VALENCE S-SHELL                **
!**         '#NODES(P)' NUMBER OF NODES OF THE VALENCE P-SHELL                **
!**         '#NODES(D)' NUMBER OF NODES OF THE VALENCE D-SHELL                **
!**         '#NODES(F)' NUMBER OF NODES OF THE VALENCE F-SHELL                **
!**         'ZCORE'     ATOMIC NUMBER OF THE CORE SHELL                       **
!**         'MAGNETIC MOMENT' NUCLEAR MAGNETIC MOMENT OF THE MOST ABUNDANT    **
!**                     ISOTOPE WITH NON-ZERO NUCLEAR SPIN                    **
!**         'GYROMAGNETICRATIO' GYROMAGNETICRATIO OF THE MOST ABUNDANT        **
!**                     ISOTOPE WITH NON-ZERO NUCLEAR SPIN                    **
!**         'SPIN'      NUCLEAR SPIN OF THE MOST ABUNDANT ISOTOPE WITH        **
!**                     NON-ZERO NUCLEAR SPIN                                 **
!**                                                                           **
!**      VALUE IS THE REQUESTED DATA                                          **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    THE ADDITIONAL ATOMS FROM 106-108 ARE DUMMY ATOMS REFERRED             **
!**      TO BY FORCE FIELDS                                                   **
!**                                                                           **
!**    THE CALLS BY ATOMIC NUMBER INTERPOLATE BETWEEN THE TWO NEIGHBORING     **
!**      ATOMIC NUMBERS                                                       **
!**                                                                           **
!**                                                                           **
!*******************************************************************************
PRIVATE
PUBLIC PERIODICTABLE$GET
INTERFACE PERIODICTABLE$GET
  MODULE PROCEDURE PERIODICTABLE$GETCH
  MODULE PROCEDURE PERIODICTABLE$GETR8
  MODULE PROCEDURE PERIODICTABLE$GETI4
  MODULE PROCEDURE PERIODICTABLE$SYMBOLGETR8
  MODULE PROCEDURE PERIODICTABLE$SYMBOLGETI4
  MODULE PROCEDURE PERIODICTABLE$GETBYZR8
END INTERFACE 
TYPE ISOTOPE_TYPE
  INTEGER(4)        :: A               ! #(PROTONS)+#(NEUTRONS)
  REAL(8)           :: ABUNDANCE    
  REAL(8)           :: S               ! SPIN
  REAL(8)           :: MAGNETICMOMENT
  REAL(8)           :: QUADRUPOLEMOMENT
END TYPE ISOTOPE_TYPE
TYPE ELEMENT_TYPE
  CHARACTER(2)      :: SYMBOL        ! ELEMENT SYMBOL
  REAL(8)           :: MASS          ! ATOMIC MASS
  REAL(8)           :: RCOV          ! COVALENT RADIUS
  REAL(8)           :: RASA          ! RCOV SCALED UP TO VOUME FILLING FCC SPHERES
  REAL(8)           :: RVDW          ! VAN DER WAALS RADIUS
  REAL(8)           :: RNUC          ! NUCLEAR RADIUS
  REAL(8)           :: EN            ! PAULING ELECTRONEGATIVITY
  INTEGER(4)        :: CORE          ! ATOMIC NUMBER OF CORE
  INTEGER(4)        :: CONFIGURATION(4) ! S,P,D,F VALENCE OCCUPATION
  INTEGER(4)        :: NODES(4)      ! #(NODES FOR S,P,D,F VALENCE WAVE FUNCTIONS)
  INTEGER(4)        :: NISOTOPES
  TYPE(ISOTOPE_TYPE),POINTER :: ISOTOPE(:)
END TYPE
INTEGER,PARAMETER     :: NEL=109
INTEGER,PARAMETER     :: IFIRSTDUMMY=106
TYPE(ELEMENT_TYPE)    :: ELEMENT(0:NEL)
LOGICAL(4)            :: TINI=.FALSE.
REAL(8)     ,PARAMETER:: U=1822.8885046287D0
REAL(8)     ,PARAMETER:: ANGSTROM=1.8897259926D0
REAL(8)     ,PARAMETER:: METER=1.889726D+10 ! DOES NOT HAVE FULL ACCURACY
! THE NUCLRAR RADIUS IS RNUC=(MASS/U)**(1/3)*1.2E-15 METER, 
! WHERE U IS THE M[C12]/12. (E.G. HALLIDAY-RESNICK-WALKER, PHYSIK, WILEY)
REAL(8)     ,PARAMETER:: RNUCFAC=1.85635065215D-6
!*******************************************************************************
CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE_INITIALIZE
!     **************************************************************************
!     **                                                                      **
!     ** ATOMIC WEIGHTS: RELATIVE ATOMIC MASS                                 **
!     **    (FOR UNSTABLE NUCLEI: ISOTOPIC MASS)                              **
!     **    FROM QUANTITIES,UNITS AND SYMBOLS IN PHYSICAL CHEMISTRY           **
!     **    BLACKWELL SCIENTIFIC PUBLICATIONS, OXFORD 1993                    **
!     **                                                                      **
!     ** THERE IS ONE ELEMENT WITH ATOMIC NUMBER ZERO. IT HAS A SMALL FINITE  **
!     ** MASS TO AVOID A DIVIDE-BY-ZERO WHILE PROPAGATING ITS POSITION        **
!     **                                                                      **
!     **************************************************************************
      TYPE SET_TYPE
        CHARACTER(2)      :: SYMBOL    ! ELEMENT SYMBOL
        REAL(8)           :: MASS      ! ATOMIC MASS
        REAL(8)           :: RCOV      ! COVALENT RADIUS
        REAL(8)           :: RVDW      ! VAN DER WAALS RADIUS
        REAL(8)           :: EN        ! PAULING SUSCEPTIBILITY
        INTEGER(4)        :: CONFIGURATION(4)  ! S,P,D,F VALENCE OCCUPATION
        CHARACTER(2)      :: CORE
      END TYPE
      TYPE(SET_TYPE)      :: SET(0:NEL)
      REAL(8)  ,PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8)             :: FACASA
      INTEGER(4)          :: I,J,IC=0
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!    
!     ==  1S   =================================================================
      SET(  0)=SET_TYPE('0 ',0.1000  ,0.53,0.530,0.00,(/0,0,0,0/),'0 ')
      SET(  1)=SET_TYPE('H ',1.00794 ,0.32,2.886,2.20,(/1,0,0,0/),'0 ')
      SET(  2)=SET_TYPE('HE',4.002602,0.93,2.362,0.00,(/2,0,0,0/),'0 ')
!     ==  2SP  =================================================================
      SET(  3)=SET_TYPE('LI',6.941    ,1.23,2.451,0.98,(/1,0,0,0/),'HE')
      SET(  4)=SET_TYPE('BE',9.012182 ,0.90,2.745,1.57,(/2,0,0,0/),'HE')
      SET(  5)=SET_TYPE('B ',10.811   ,0.82,4.083,2.04,(/2,1,0,0/),'HE')
      SET(  6)=SET_TYPE('C ',12.011   ,0.77,3.851,2.55,(/2,2,0,0/),'HE')
      SET(  7)=SET_TYPE('N ',14.00674 ,0.75,3.660,3.04,(/2,3,0,0/),'HE')
      SET(  8)=SET_TYPE('O ',15.9994  ,0.73,3.500,3.44,(/2,4,0,0/),'HE')
      SET(  9)=SET_TYPE('F ',18.998403,0.72,3.364,3.98,(/2,5,0,0/),'HE')
      SET( 10)=SET_TYPE('NE',20.1797  ,0.71,3.243,0.00,(/2,6,0,0/),'HE')
!     == 3SP  ==================================================================
      SET( 11)=SET_TYPE('NA',22.989768,1.54,2.983,0.93,(/1,0,0,0/),'NE')
      SET( 12)=SET_TYPE('MG',24.3050  ,1.36,3.021,1.31,(/2,0,0,0/),'NE')
      SET( 13)=SET_TYPE('AL',26.981539,1.18,4.499,1.61,(/2,1,0,0/),'NE')
      SET( 14)=SET_TYPE('SI',28.0855  ,1.11,4.295,1.90,(/2,2,0,0/),'NE')
      SET( 15)=SET_TYPE('P ',30.973762,1.06,4.147,2.19,(/2,3,0,0/),'NE')
      SET( 16)=SET_TYPE('S ',32.066   ,1.02,4.035,2.58,(/2,4,0,0/),'NE')
      SET( 17)=SET_TYPE('CL',35.4527  ,0.99,3.947,3.16,(/2,5,0,0/),'NE')
      SET( 18)=SET_TYPE('AR',39.948   ,0.98,3.868,0.00,(/2,6,0,0/),'NE')
!     ==  4S   =================================================================
      SET( 19)=SET_TYPE('K ',39.0983  ,2.03,3.812,0.82,(/1,0,0,0/),'AR')
      SET( 20)=SET_TYPE('CA',40.078   ,1.74,3.399,1.00,(/2,0,0,0/),'AR')
!     ==  3D   =================================================================
      SET( 21)=SET_TYPE('SC',44.955910,1.44,3.295,1.36,(/2,0,1,0/),'AR')
      SET( 22)=SET_TYPE('TI',47.88    ,1.32,3.175,1.54,(/2,0,2,0/),'AR')
      SET( 23)=SET_TYPE('V ',50.9415  ,1.22,3.144,1.63,(/2,0,3,0/),'AR')
      SET( 24)=SET_TYPE('CR',51.9961  ,1.18,3.023,1.66,(/1,0,5,0/),'AR')
      SET( 25)=SET_TYPE('MN',54.93805 ,1.17,2.961,1.55,(/2,0,5,0/),'AR')
      SET( 26)=SET_TYPE('FE',55.847   ,1.17,2.912,1.83,(/2,0,6,0/),'AR')
      SET( 27)=SET_TYPE('CO',58.93320 ,1.16,2.872,1.88,(/2,0,7,0/),'AR')
      SET( 28)=SET_TYPE('NI',58.34    ,1.15,2.834,1.91,(/2,0,8,0/),'AR')
      SET( 29)=SET_TYPE('CU',63.546   ,1.17,3.495,1.90,(/1,0,10,0/),'AR')
      SET( 30)=SET_TYPE('ZN',65.39    ,1.25,2.763,1.65,(/2,0,10,0/),'AR')
!     ==  4P   =================================================================
      SET( 31)=SET_TYPE('GA',60.723   ,1.26,4.383,1.81,(/2,1,10,0/),'AR')
      SET( 32)=SET_TYPE('GE',72.61    ,1.22,4.280,2.01,(/2,2,10,0/),'AR')
      SET( 33)=SET_TYPE('AS',74.92159 ,1.20,4.230,2.18,(/2,3,10,0/),'AR')
      SET( 34)=SET_TYPE('SE',78.96    ,1.16,4.205,2.55,(/2,4,10,0/),'AR')
      SET( 35)=SET_TYPE('BR',79.904   ,1.14,4.189,2.96,(/2,5,10,0/),'AR')
      SET( 36)=SET_TYPE('KR',83.80    ,1.12,4.141,0.00,(/2,6,10,0/),'AR')
!     ==  5S   =================================================================
      SET( 37)=SET_TYPE('RB',85.4678  ,2.16,4.114,0.82,(/1,0,0,0/),'KR')
      SET( 38)=SET_TYPE('SR',87.62    ,1.91,3.641,0.95,(/2,0,0,0/),'KR')
!     ==  4D   =================================================================
      SET( 39)=SET_TYPE('Y ',88.90585 ,1.62,3.345,1.22,(/2,0,1,0/),'KR')
      SET( 40)=SET_TYPE('ZR',91.224   ,1.45,3.124,1.33,(/2,0,2,0/),'KR')
      SET( 41)=SET_TYPE('NB',92.90638 ,1.34,3.165,1.60,(/1,0,4,0/),'KR')
      SET( 42)=SET_TYPE('MO',95.94    ,1.30,3.052,2.16,(/1,0,5,0/),'KR')
      SET( 43)=SET_TYPE('TC',97.907215,1.27,2.998,1.90,(/2,0,5,0/),'KR')
      SET( 44)=SET_TYPE('RU',101.07   ,1.25,2.963,2.20,(/1,0,7,0/),'KR')
      SET( 45)=SET_TYPE('RH',102.90550,1.25,2.929,2.28,(/1,0,8,0/),'KR')
      SET( 46)=SET_TYPE('PD',106.42   ,1.28,2.899,2.20,(/0,0,10,0/),'KR')
      SET( 47)=SET_TYPE('AG',107.868  ,1.34,3.148,1.93,(/1,0,10,0/),'KR')
      SET( 48)=SET_TYPE('CD',112.411  ,1.48,2.848,1.69,(/2,0,10,0/),'KR')
!     ==  5P ===================================================================
      SET( 49)=SET_TYPE('IN',114.818  ,1.44,4.463,1.78,(/2,1,10,0/),'KR')
      SET( 50)=SET_TYPE('SN',118.710  ,1.41,4.392,1.96,(/2,2,10,0/),'KR')
      SET( 51)=SET_TYPE('SB',121.757  ,1.40,4.420,2.05,(/2,3,10,0/),'KR')
      SET( 52)=SET_TYPE('TE',127.60   ,1.36,4.470,2.10,(/2,4,10,0/),'KR')
      SET( 53)=SET_TYPE('I ',126.90447,1.33,4.500,2.66,(/2,5,10,0/),'KR')
      SET( 54)=SET_TYPE('XE',131.29   ,1.31,4.404,0.00,(/2,6,10,0/),'KR')
!     == 6S  ===================================================================
      SET( 55)=SET_TYPE('CS',132.90543,2.35,4.517,0.79,(/1,0,0,0/),'XE')
      SET( 56)=SET_TYPE('BA',137.327  ,1.98,3.703,0.89,(/2,0,0,0/),'XE')
!     == 5D TRANSITION METALS ==================================================
      SET( 57)=SET_TYPE('LA',138.9055 ,1.69,3.522,1.10,(/2,0,1,0/),'XE')
!     == LANTHANIDES - 4F RARE EARTHS ==========================================
      SET( 58)=SET_TYPE('CE',140.115  ,1.65,3.556,1.12,(/2,0,1,1/),'XE')
      SET( 59)=SET_TYPE('PR',140.90765,1.65,3.606,1.13,(/2,0,0,3/),'XE')
      SET( 60)=SET_TYPE('ND',144.24   ,1.64,3.575,1.14,(/2,0,0,4/),'XE')
      SET( 61)=SET_TYPE('PM',145.     ,1.63,3.547,1.13,(/2,0,0,5/),'XE')
      SET( 62)=SET_TYPE('SM',150.36   ,1.62,3.520,1.17,(/2,0,0,6/),'XE')
      SET( 63)=SET_TYPE('EU',151.965  ,1.85,3.493,1.20,(/2,0,0,7/),'XE')
      SET( 64)=SET_TYPE('GD',157.25   ,1.61,3.368,1.20,(/2,0,1,7/),'XE')
      SET( 65)=SET_TYPE('TB',158.92534,1.59,3.451,1.20,(/2,0,0,9/),'XE')
      SET( 66)=SET_TYPE('DY',162.50   ,1.59,3.428,1.22,(/2,0,0,10/),'XE')
      SET( 67)=SET_TYPE('HO',164.93032,1.58,3.409,1.23,(/2,0,0,11/),'XE')
      SET( 68)=SET_TYPE('ER',167.26   ,1.57,3.391,1.24,(/2,0,0,12/),'XE')
      SET( 69)=SET_TYPE('TM',168.93421,1.56,3.374,1.25,(/2,0,0,13/),'XE')
      SET( 70)=SET_TYPE('YB',173.04   ,1.74,3.355,1.10,(/2,0,0,14/),'XE')
      SET( 71)=SET_TYPE('LU',174.967  ,1.56,3.640,1.27,(/2,0,1,14/),'XE')
!     == 5D TRANSITION METALS (CONTINUATION)====================================
      SET( 72)=SET_TYPE('HF',178.49   ,1.44,3.141,1.30,(/2,0,2,14/),'XE')
      SET( 73)=SET_TYPE('TA',180.9479 ,1.34,3.170,1.50,(/2,0,3,14/),'XE')
      SET( 74)=SET_TYPE('W ',183.84   ,1.30,3.069,2.36,(/2,0,4,14/),'XE')
      SET( 75)=SET_TYPE('RE',186.207  ,1.28,2.954,1.90,(/2,0,5,14/),'XE')
      SET( 76)=SET_TYPE('OS',190.23   ,1.26,3.120,2.20,(/2,0,6,14/),'XE')
      SET( 77)=SET_TYPE('IR',192.22   ,1.27,2.840,2.20,(/2,0,7,14/),'XE')
      SET( 78)=SET_TYPE('PT',195.08   ,1.30,2.754,2.28,(/1,0,9,14/),'XE')
      SET( 79)=SET_TYPE('AU',196.96654,1.34,3.293,2.54,(/1,0,10,14/),'XE')
      SET( 80)=SET_TYPE('HG',200.59   ,1.49,2.705,2.00,(/2,0,10,14/),'XE')
!     == 6P ======================= ============================================
      SET( 81)=SET_TYPE('TL',204.3833 ,1.48,4.347,2.04,(/2,1,10,14/),'XE')
      SET( 82)=SET_TYPE('PB',207.2    ,1.47,4.297,2.33,(/2,2,10,14/),'XE')
      SET( 83)=SET_TYPE('BI',208.98037,1.46,4.370,2.02,(/2,3,10,14/),'XE')
      SET( 84)=SET_TYPE('PO',208.98240,1.46,4.709,2.00,(/2,4,10,14/),'XE')
      SET( 85)=SET_TYPE('AT',209.98713,1.45,4.750,2.20,(/2,5,10,14/),'XE')
      SET( 86)=SET_TYPE('RN',222.01757,1.44,4.765,0.00,(/2,6,10,14/),'XE')
!     ==  7S ===================================================================
      SET( 87)=SET_TYPE('FR',223.01973,2.50,4.900,0.70,(/1,0,0,0/),'RN')
      SET( 88)=SET_TYPE('RA',226.02540,2.00,3.677,0.90,(/2,0,0,0/),'RN')
!     ==  6D TRANSITION METALS =================================================
      SET( 89)=SET_TYPE('AC',227.02775,1.65,3.478,1.10,(/2,0,1,0/),'RN')
!     ==  ACTINIDES - 5F RARE EARTHS ===========================================
      SET( 90)=SET_TYPE('TH',232.0381 ,1.65,3.396,1.30,(/2,0,2,0/),'RN')
      SET( 91)=SET_TYPE('PA',231.03588,1.65,3.424,1.50,(/2,0,1,2/),'RN')
      SET( 92)=SET_TYPE('U ',238.0289 ,1.42,3.395,1.38,(/2,0,1,3/),'RN')
      SET( 93)=SET_TYPE('NP',237.048  ,1.63,3.424,1.36,(/2,0,1,4/),'RN')
      SET( 94)=SET_TYPE('PU',244.06420,1.62,3.424,1.28,(/2,0,0,6/),'RN')
      SET( 95)=SET_TYPE('AM',243.06138,1.85,3.381,1.30,(/2,0,0,7/),'RN')
      SET( 96)=SET_TYPE('CM',247.07035,1.61,3.326,1.30,(/2,0,1,7/),'RN')
      SET( 97)=SET_TYPE('BK',247.07030,1.59,3.339,1.30,(/2,0,0,9/),'RN')
      SET( 98)=SET_TYPE('CF',251.07958,1.59,3.313,1.30,(/2,0,0,10/),'RN')
      SET( 99)=SET_TYPE('ES',252.08294,1.58,3.299,1.30,(/2,0,0,11/),'RN')
      SET(100)=SET_TYPE('FM',257.09510,1.57,3.286,1.30,(/2,0,0,12/),'RN')
      SET(101)=SET_TYPE('MD',258.09857,1.56,3.274,1.30,(/2,0,0,13/),'RN')
      SET(102)=SET_TYPE('NO',259.10093,1.74,3.248,1.30,(/2,0,0,14/),'RN')
      SET(103)=SET_TYPE('LR',260.10532,1.56,3.236,1.30,(/2,0,1,14/),'RN')
!     ==  6D TRANSITION METALS (CONTINUATION, ==================================
      SET(104)=SET_TYPE('RF',261.10869,1.44,3.500,0.00,(/2,0,2,14/),'RN')
      SET(105)=SET_TYPE('HA',262.11376,1.34,3.500,0.00,(/2,0,3,14/),'RN')
!     ==  DUMMY ATOMS FOR FORCE FIELDS =========================================
      SET(106)=SET_TYPE('CP',1.0000 ,0.00,0.000,2.55,(/0,0,0,0/),'0 ')
      SET(107)=SET_TYPE('PI',1.0000 ,0.00,0.000,2.55,(/0,0,0,0/),'0 ')
      SET(108)=SET_TYPE('CI',1.0000 ,0.00,0.000,2.55,(/0,0,0,0/),'0 ')
      SET(109)=SET_TYPE('M ',1.0000 ,0.00,0.000,2.55,(/0,0,0,0/),'0 ')  !TIP4P
!
!     ==========================================================================
!     ==  MAP INTO ELEMENT AND CONVERT TO ATOMIC UNITS                        ==
!     ==========================================================================
      FACASA=4.D0/SQRT(2.D0)*(3.D0/(16.D0*PI))**(1.D0/3.D0)
      DO I=0,NEL
        ELEMENT(I)%SYMBOL=SET(I)%SYMBOL
        ELEMENT(I)%MASS=SET(I)%MASS*U
        ELEMENT(I)%RCOV=SET(I)%RCOV*ANGSTROM
!       __LISTED ARE DISTANCES, NOT RADII: CONVERT TO RADII_____________________
        ELEMENT(I)%RVDW=SET(I)%RVDW*ANGSTROM*0.5D0 
        ELEMENT(I)%CONFIGURATION=SET(I)%CONFIGURATION
        ELEMENT(I)%RASA=FACASA*ELEMENT(I)%RCOV
        ELEMENT(I)%EN=SET(I)%EN
        ELEMENT(I)%RNUC=ELEMENT(I)%MASS**(1.D0/3.D0)*RNUCFAC
      ENDDO
      ELEMENT(0)%NODES(:)=(/0,0,0,0/)
      DO I=0,NEL
        CALL PERIODICTABLE_ATOMICNUMBER(SET(I)%CORE,IC)
        ELEMENT(I)%CORE=IC
        DO J=1,4
          IF(ELEMENT(IC)%CONFIGURATION(J).GT.0) THEN
            ELEMENT(I)%NODES(J)=ELEMENT(IC)%NODES(J)+1
          END IF
        ENDDO
        ELEMENT(I)%NISOTOPES=0
        NULLIFY(ELEMENT(I)%ISOTOPE)
      ENDDO
      CALL PERIODICTABLE_ISOTOPES
      RETURN
      END SUBROUTINE PERIODICTABLE_INITIALIZE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE_ISOTOPES
!     **************************************************************************
!     **  PROVIDES INFORMATION FOR INDIVIDUAL ISOTOPES FOR EACH ELEMENT       **
!     **  THIS INFORMATION IS NOT COMPLETE!!                                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)         :: NUCLEARMAGNETON   
      REAL(8)         :: METER
      REAL(8)         :: FEMTO
      REAL(8)         :: FM2   ! (FEMTOMETER)**2
      INTEGER(4)      :: IZ,I
      INTEGER(4)      :: NISOTOPES
!     **************************************************************************
!===============================================================================
!==  ABUNDANCE IS GIVEN IN PERCENT
!==  MAGNETIC MOMENT IN NUCLEAR MAGNETONS
!==  SPIN IN HBAR
!===============================================================================
!== HYDROGEN ==========================================================
IZ=1
NISOTOPES=3
ELEMENT(IZ)%NISOTOPES=NISOTOPES
ALLOCATE(ELEMENT(IZ)%ISOTOPE(NISOTOPES))
ELEMENT(IZ)%ISOTOPE(1)=ISOTOPE_TYPE( 1,99.985D0   ,0.5D0,2.79284738663,0.0000)
ELEMENT(IZ)%ISOTOPE(2)=ISOTOPE_TYPE( 2, 0.0150D0  ,1.0D0,0.85743823024,0.2860)
ELEMENT(IZ)%ISOTOPE(3)=ISOTOPE_TYPE( 3, 0.D0      ,0.5D0,2.97896247968D0,0.D0)
!== HELIUM ===========================================================
IZ=2
NISOTOPES=2
ELEMENT(IZ)%NISOTOPES=NISOTOPES
ALLOCATE(ELEMENT(IZ)%ISOTOPE(NISOTOPES))
ELEMENT(IZ)%ISOTOPE(1)=ISOTOPE_TYPE( 4,99.999863D0,0.D0,0.D0,0.D0)
ELEMENT(IZ)%ISOTOPE(2)=ISOTOPE_TYPE( 3,0.0001373D0,0.5D0,-2.12762484866,0.D0)
!== OXYGEN ===========================================================
IZ=8
NISOTOPES=3
ELEMENT(IZ)%NISOTOPES=NISOTOPES
ALLOCATE(ELEMENT(IZ)%ISOTOPE(NISOTOPES))
ELEMENT(IZ)%ISOTOPE(1)=ISOTOPE_TYPE(16,99.6349D0  ,0.0D0,0.D0,0.D0)           
ELEMENT(IZ)%ISOTOPE(2)=ISOTOPE_TYPE(18, 0.20012D0 ,0.0D0,0.D0,0.D0)           
ELEMENT(IZ)%ISOTOPE(3)=ISOTOPE_TYPE(17, 0.0383D0  ,2.5D0,-1.8938D0,-2.55822D0)
!== SILICON ===========================================================
IZ=14
NISOTOPES=3
ELEMENT(IZ)%NISOTOPES=NISOTOPES
ALLOCATE(ELEMENT(IZ)%ISOTOPE(NISOTOPES))
ELEMENT(IZ)%ISOTOPE(1)=ISOTOPE_TYPE(28,92.23D0    ,0.0D0,0.0,0.0)          
ELEMENT(IZ)%ISOTOPE(2)=ISOTOPE_TYPE(29, 4.671D0   ,0.5D0,-0.555293D0,0.D0) 
ELEMENT(IZ)%ISOTOPE(3)=ISOTOPE_TYPE(30, 3.101D0   ,0.D0,0.D0,0.D0)         
!== IRON ==============================================================
IZ=26
NISOTOPES=3
ELEMENT(IZ)%NISOTOPES=NISOTOPES
ALLOCATE(ELEMENT(IZ)%ISOTOPE(NISOTOPES))
ELEMENT(IZ)%ISOTOPE(1)=ISOTOPE_TYPE(56,91.7230D0  ,0.D0,0.D0,0.D0)        
ELEMENT(IZ)%ISOTOPE(2)=ISOTOPE_TYPE(54, 5.81D0    ,0.D0,0.D0,0.D0)        
ELEMENT(IZ)%ISOTOPE(3)=ISOTOPE_TYPE(57, 2.21D0    ,0.5D0,0.090623009,0.D0)
!
!=======================================================================
!== CONVERT INTO ATOMIC UNITS                                         ==
!=======================================================================
CALL CONSTANTS$GET('NUCLEARMAGNETON',NUCLEARMAGNETON)
CALL CONSTANTS$GET('METER',METER)
CALL CONSTANTS$GET('FEMTO',FEMTO)
FM2=(FEMTO*METER)**2
DO IZ=1,NEL
  DO I=1,ELEMENT(IZ)%NISOTOPES
    ELEMENT(IZ)%ISOTOPE(I)%ABUNDANCE &
&     =ELEMENT(IZ)%ISOTOPE(I)%ABUNDANCE*0.01D0
    ELEMENT(IZ)%ISOTOPE(I)%MAGNETICMOMENT &
&     =ELEMENT(IZ)%ISOTOPE(I)%MAGNETICMOMENT*NUCLEARMAGNETON
    ELEMENT(IZ)%ISOTOPE(I)%QUADRUPOLEMOMENT &
&     =ELEMENT(IZ)%ISOTOPE(I)%QUADRUPOLEMOMENT*FM2
  ENDDO
ENDDO
RETURN
END SUBROUTINE PERIODICTABLE_ISOTOPES
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE_ATOMICNUMBER(SYMBOL_,IZ)
!     **************************************************************************
!     ** RETURN INTEGER ATOMIC NUMBER FOR A GIVEN SYMBOL, RESPECETIVELY       **
!     ** THE ENTRY INDEX                                                      **
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      CHARACTER(*),INTENT(IN) :: SYMBOL_
      INTEGER(4)  ,INTENT(OUT):: IZ
      CHARACTER(2)            :: SYMBOL
      INTEGER(4)              :: I
!     **************************************************************************
      SYMBOL=SYMBOL_(1:2)
      IF(SYMBOL(2:2).EQ.'_') SYMBOL(2:2)=' '
      DO I=0,NEL
        IF(+SYMBOL.EQ.ELEMENT(I)%SYMBOL) THEN
          IZ=I
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('ELEMENT-SYMBOL NOT RECOGNIZED')
      CALL ERROR$CHVAL('SYMBOL',SYMBOL_)
      CALL ERROR$STOP('PERIODICTABLE_ATOMICNUMBER')
      RETURN
      END SUBROUTINE PERIODICTABLE_ATOMICNUMBER
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE$GETCH(IZ,ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: IZ
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(OUT):: VAL
!     **************************************************************************
      CALL PERIODICTABLE_INITIALIZE
      IF(IZ.LT.0.OR.IZ.GE.NEL) THEN
        CALL ERROR$MSG('ATOMIC NUMBER OUT OF RANGE')
        CALL ERROR$I4VAL('IZ',IZ)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETCH')
      END IF
      IF(ID.EQ.'SYMBOL') THEN
        VAL=ELEMENT(IZ)%SYMBOL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETCH')
      END IF
      RETURN
      END SUBROUTINE PERIODICTABLE$GETCH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE$GETR8(IZ,ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: IZ
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
      INTEGER(4)              :: I
!     **************************************************************************
      CALL PERIODICTABLE_INITIALIZE
      IF(IZ.LT.0.OR.IZ.GT.NEL) THEN
        CALL ERROR$MSG('ATOMIC NUMBER OUT OF RANGE')
        CALL ERROR$I4VAL('IZ',IZ)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETR8')
      END IF
      IF(ID.EQ.'R(COV)') THEN
        VAL=ELEMENT(IZ)%RCOV
      ELSE IF(ID.EQ.'R(VDW)') THEN
        VAL=ELEMENT(IZ)%RVDW
      ELSE IF(ID.EQ.'R(ASA)') THEN
        VAL=ELEMENT(IZ)%RASA
      ELSE IF(ID.EQ.'MASS') THEN
        VAL=ELEMENT(IZ)%MASS
      ELSE IF(ID.EQ.'Z') THEN
        VAL=REAL(IZ)
      ELSE IF(ID.EQ.'ZCORE') THEN
        VAL=REAL(ELEMENT(IZ)%CORE)
      ELSE IF(ID.EQ.'EN') THEN
        VAL=ELEMENT(IZ)%EN
      ELSE IF(ID.EQ.'RNUC') THEN
        VAL=ELEMENT(IZ)%RNUC
      ELSE IF(ID.EQ.'MAGNETICMOMENT') THEN
        VAL=0.D0
        DO I=1,ELEMENT(IZ)%NISOTOPES
          IF(ELEMENT(IZ)%ISOTOPE(I)%S.NE.0.D0) THEN
            VAL=ELEMENT(IZ)%ISOTOPE(I)%MAGNETICMOMENT
            EXIT
          END IF
        ENDDO
      ELSE IF(ID.EQ.'SPIN') THEN
        VAL=0.D0
        DO I=1,ELEMENT(IZ)%NISOTOPES
          IF(ELEMENT(IZ)%ISOTOPE(I)%S.NE.0.D0) THEN
            VAL=ELEMENT(IZ)%ISOTOPE(I)%S
            EXIT
          END IF
        ENDDO
      ELSE IF(ID.EQ.'GYROMAGNETICRATIO') THEN
        VAL=0.D0
        DO I=1,ELEMENT(IZ)%NISOTOPES
          IF(ELEMENT(IZ)%ISOTOPE(I)%S.NE.0.D0) THEN
            VAL=ELEMENT(IZ)%ISOTOPE(I)%MAGNETICMOMENT &
     &         /ELEMENT(IZ)%ISOTOPE(I)%S
            EXIT
          END IF
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETR8')
      END IF
      RETURN
      END SUBROUTINE PERIODICTABLE$GETR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE$GETBYZR8(Z,ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)     ,INTENT(IN) :: Z
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
      INTEGER(4)              :: IZ1,IZ2
      REAL(8)                 :: C1,C2
      LOGICAL(4)              :: TSUBALKALI
!     **************************************************************************
      CALL PERIODICTABLE_INITIALIZE
      IF(Z.LT.0.OR.Z.GT.REAL(NEL)) THEN
        CALL ERROR$MSG('ATOMIC NUMBER OUT OF RANGE')
        CALL ERROR$R8VAL('Z',Z)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETBYZR8')
      END IF
      IZ1=INT(Z)
      IZ2=IZ1+1
      C2=Z-REAL(IZ1)
      C1=1.D0-C2
      IZ2=MIN(IZ2,NEL)
      TSUBALKALI=.FALSE.
!      TSUBALKALI=TSUBALKALI.OR.ELEMENT(IZ2)%SYMBOL.EQ.'H '
      TSUBALKALI=TSUBALKALI.OR.ELEMENT(IZ2)%SYMBOL.EQ.'LI'
      TSUBALKALI=TSUBALKALI.OR.ELEMENT(IZ2)%SYMBOL.EQ.'NA'
      TSUBALKALI=TSUBALKALI.OR.ELEMENT(IZ2)%SYMBOL.EQ.'K '
      TSUBALKALI=TSUBALKALI.OR.ELEMENT(IZ2)%SYMBOL.EQ.'RB'
      TSUBALKALI=TSUBALKALI.OR.ELEMENT(IZ2)%SYMBOL.EQ.'CS'
      TSUBALKALI=TSUBALKALI.OR.ELEMENT(IZ2)%SYMBOL.EQ.'FR'
!      IF(C2.LT.1.D-5) TSUBALKALI=.FALSE.
      IF(ID.EQ.'R(COV)') THEN
        IF(TSUBALKALI) THEN
          VAL=ELEMENT(IZ2)%RCOV
        ELSE
          VAL=C1*ELEMENT(IZ1)%RCOV+C2*ELEMENT(IZ2)%RCOV
        END IF
      ELSE IF(ID.EQ.'R(VDW)') THEN
        IF(TSUBALKALI) THEN
          VAL=ELEMENT(IZ2)%RVDW
        ELSE
          VAL=C1*ELEMENT(IZ1)%RVDW+C2*ELEMENT(IZ2)%RVDW
        END IF
      ELSE IF(ID.EQ.'R(ASA)') THEN
        IF(TSUBALKALI) THEN
          VAL=ELEMENT(IZ2)%RASA
        ELSE
          VAL=C1*ELEMENT(IZ1)%RASA+C2*ELEMENT(IZ2)%RASA
        END IF
      ELSE IF(ID.EQ.'MASS') THEN
        VAL=C1*ELEMENT(IZ1)%MASS+C2*ELEMENT(IZ2)%MASS
      ELSE IF(ID.EQ.'RNUC') THEN
        VAL=C1*ELEMENT(IZ1)%RNUC+C2*ELEMENT(IZ2)%RNUC
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETBYZR8')
      END IF
      RETURN
      END SUBROUTINE PERIODICTABLE$GETBYZR8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE$GETI4(IZ,ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: IZ
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      CALL PERIODICTABLE_INITIALIZE
      IF(IZ.LT.0.OR.IZ.GT.NEL) THEN
        CALL ERROR$MSG('ATOMIC NUMBER OUT OF RANGE')
        CALL ERROR$I4VAL('IZ',IZ)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETI4')
      END IF
      IF(ID.EQ.'Z') THEN
        VAL=IZ       
      ELSE IF(ID.EQ.'ZCORE') THEN
        VAL=ELEMENT(IZ)%CORE
      ELSE IF(ID.EQ.'OCC(S)') THEN
        VAL=ELEMENT(IZ)%CONFIGURATION(1)
      ELSE IF(ID.EQ.'OCC(P)') THEN
        VAL=ELEMENT(IZ)%CONFIGURATION(2)
      ELSE IF(ID.EQ.'OCC(D)') THEN
        VAL=ELEMENT(IZ)%CONFIGURATION(3)
      ELSE IF(ID.EQ.'OCC(F)') THEN
        VAL=ELEMENT(IZ)%CONFIGURATION(4)
      ELSE IF(ID.EQ.'#NODES(S)') THEN
        VAL=ELEMENT(IZ)%NODES(1)
      ELSE IF(ID.EQ.'#NODES(P)') THEN
        VAL=ELEMENT(IZ)%NODES(2)
      ELSE IF(ID.EQ.'#NODES(D)') THEN
        VAL=ELEMENT(IZ)%NODES(3)
      ELSE IF(ID.EQ.'#NODES(F)') THEN
        VAL=ELEMENT(IZ)%NODES(4)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PERIODICTABLE$GETI4')
      END IF
      RETURN
      END SUBROUTINE PERIODICTABLE$GETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE$SYMBOLGETI4(SYMBOL,ID,VAL)
!     **************************************************************************
!     **  INTERFACE FOR CALLS BY ELEMENT SYMBOL                               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: SYMBOL
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
      INTEGER(4)              :: IZ
!     **************************************************************************
      CALL PERIODICTABLE_INITIALIZE
      CALL PERIODICTABLE_ATOMICNUMBER(SYMBOL,IZ)
      CALL PERIODICTABLE$GETI4(IZ,ID,VAL)
      RETURN
      END SUBROUTINE PERIODICTABLE$SYMBOLGETI4
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE$SYMBOLGETR8(SYMBOL,ID,VAL)
!     **************************************************************************
!     **  INTERFACE FOR CALLS BY ELEMENT SYMBOL                               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: SYMBOL
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
      INTEGER(4)              :: IZ
!     **************************************************************************
      CALL PERIODICTABLE_INITIALIZE
      CALL PERIODICTABLE_ATOMICNUMBER(SYMBOL,IZ)
      CALL PERIODICTABLE$GETR8(IZ,ID,VAL)
      RETURN
      END SUBROUTINE PERIODICTABLE$SYMBOLGETR8
!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PERIODICTABLE$REPORT(NFIL)
!     **************************************************************************
!     **  REPORT SETTINGS OF PERIODIC TABLE                                   **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)              :: IZ
      CHARACTER(5)            :: CHAR
      CHARACTER(15)           :: CONFIG
      CHARACTER(15)           :: NNODES
      INTEGER(4)              :: N
      CHARACTER(1)            :: SPDF(4)=(/'S','P','D','F'/)
!     ******************************************************************
      CALL PERIODICTABLE_INITIALIZE
      WRITE(NFIL,FMT='(72("=")/72("="),T10,"  PERIODICTABLE  "/72("="))')
      WRITE(NFIL,FMT='(A2,T5,A3,T10,A10,T20,A10,T30,A10,T40,A10,T51,A13,T66,A8)') &
     &   "SY","Z","MASS[U]","R(COV)","R(ASA)","R(VDW)","CONFIGURATION","#(NODES)"
      DO IZ=1,NEL
        CONFIG='['//ELEMENT(ELEMENT(IZ)%CORE)%SYMBOL//']'
        NNODES=' '
        DO N=1,4
          IF(ELEMENT(IZ)%CONFIGURATION(N).NE.0) THEN
            WRITE(CHAR,*)ELEMENT(IZ)%CONFIGURATION(N)
            CONFIG=TRIM(CONFIG)//TRIM(SPDF(N))//ADJUSTL(CHAR)
          END IF
          IF(ELEMENT(IZ)%NODES(N).NE.0) THEN
            WRITE(CHAR,*)ELEMENT(IZ)%NODES(N)
            NNODES=TRIM(NNODES)//TRIM(SPDF(N))//ADJUSTL(CHAR)
          END IF
        ENDDO
        WRITE(NFIL,FMT='(A2,T5,I3,T10,F10.3,T20,F10.3,T30,F10.3,T40,F10.3,T51,A15,T66,A10)') &
     &  ELEMENT(IZ)%SYMBOL,IZ,ELEMENT(IZ)%MASS/U,ELEMENT(IZ)%RCOV &
     &              ,ELEMENT(IZ)%RASA,ELEMENT(IZ)%RVDW,CONFIG,NNODES
      ENDDO
      RETURN
      END SUBROUTINE PERIODICTABLE$REPORT
END MODULE PERIODICTABLE_MODULE
          
