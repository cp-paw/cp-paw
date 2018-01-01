      SUBROUTINE CONSTANTS(STRING,VAL)
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: STRING
      REAL(8)          ,INTENT(OUT):: VAL
      CALL CONSTANTS$GET(STRING,VAL)
      RETURN
      END
      SUBROUTINE CONSTANTS$LIST(NFIL)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL 
      CALL CONSTANTS$REPORT(NFIL)
      RETURN
      END
!     ==================================================================
!     ==================================================================
!     ==== GENERAL PURPOSE ROUTINES ====================================
!     ==== (NOT DEPENDEND ON PARAMETERS AND COMMON BLOCKS)==============
!     ==================================================================
!     ==================================================================
!
!     ..................................................................
MODULE CONSTANTS_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: CONSTANTS                                                  **
!**                                                                   **
!**  PURPOSE: PROVIDE FUNDAMENTAL CONSTANTS AND CONVERSION UNITS      **
!**                                                                   **
!**  SOURCE IS:                                                       **
!**  "THE FUNDAMENTAL PHYSICAL CONSTANTS"                             **
!**   E.R. COHEN AND B.N. TAYLOR,                                     **
!**   PHYSICS TODAY (BUYERS GUIDE), P.9-14 AUGUST 1993                **
!**                                                                   **
!***********************************************************************
!**                                                                   **
!**  ATOMIC UNITS:                                                    **
!**      HBAR=E=M_E=4*PI*EPSILON_0=1                                  **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!***********************************************************************
!**********************************************************************=
!** SI UNITS:                                                         **
!** METER,KILOGRAMMG,SECOND,AMPERE,KELVIN,MOLE,CANDELA                **
!**                                                                   **
!** 1 METER IS 1650763.73 TIMES THE WAVE LENGTH OF LIGHT IN           **
!**          VACUUM EMITTED BY THE 86 ISOTOPE OF KR                   **
!**          IN THE TRANSITION 5D5->2P10                              **
!**                                                                   **
!** 1 KILOGRAMM PROTOTYPE                                             **
!**                                                                   **
!** 1 SECOND IS THE 9192631770 THE PERIOD OF LIGHT EMITTED            **
!**          IN THE TRANSITION BETWEEN THE HYPERFINESTRUCTUIRE        **
!**          NIVEAUS OF THE 133 ISOTOPE OF CR                         **
!**                                                                   **
!** 1 AMPERE IS THE CURRENT WHICH RESULTS IN A FORCE OF               **
!**          0.2*10**-6 NEWTON WHEN TWO INFINITESIMAL THIN,           **
!**          1M LONG WIRES ARE PLACED IN A DISTACE OF 1 M             **
!**                                                                   **
!** 1 KELVIN IS 1/273.16 TIMES THE THERMODYNAMIC TEMPERATURE          **
!**          OF THE TRIPLE POINT OF WATER                             **
!**                                                                   **
!** 1 MOLE   IS THE NUMBER OF C^12 ATOMS IN 12G                       **
!**                                                                   **
!** CANDELA ....                                                      **
!**********************************************************************=
TYPE CONSTANT
  CHARACTER(32) :: NAME
  REAL(8)       :: VALUE
  CHARACTER(64) :: DESCRIPTION
END TYPE CONSTANT
INTEGER(4)    ,PARAMETER :: NCX=100
INTEGER(4)               :: NC   ! ACTUAL #(CONSTANTS)
TYPE(CONSTANT)           :: XXX(NCX)
LOGICAL(4)               :: TFIRST=.TRUE.
!***********************************************************************
!** FUNDAMENTAL ATOMIC UNITS                                          **
!** THIS UNIT SYSTEM IS DEFINED AS (E=ME=HBAR=4*PI*EPSILON_0=1)       **
!***********************************************************************
REAL(8),   PARAMETER :: E    =1.D0               ! ELEMENTARY CHARGE
REAL(8),   PARAMETER :: ME   =1.D0               ! ELECTRON MASS
REAL(8),   PARAMETER :: HBAR =1.D0               ! PLANCK CONSTANT
REAL(8),   PARAMETER :: ALPHA=7.2973530833D-3    ! FINESTRUCTURE CONSTANT =E**2/(4\PI\EPSILON_0*C)
REAL(8),   PARAMETER :: C=E*E/(HBAR*ALPHA)       ! SPEED OF LIGHT
REAL(8),   PARAMETER :: KELVIN=1.D0              ! TEMPERATURE
!***********************************************************************
!** CONVERSION FROM ATOMIC UNITS TO SI UNITS USING FOUR CONSTANTS     **
!** C=299792458 M/SEC                                                 **
!** HBAR=1.0545726663D-34 KG*M**2/SEC                                 **
!** ME=9.109389754D+31 KG                                             **
!** E=1.6021773349D-19 A*S                                            **
!***********************************************************************
REAL(8)   ,PARAMETER :: CBYMBYS =299792458D0     ! VELOCITY OF LIGHT IN [M/S]
REAL(8)   ,PARAMETER :: HBARBYJS=1.0545726663D-34! PLANCK CONSTANT / [JOULE*SECOND] 
REAL(8)   ,PARAMETER :: MEBYKG  =9.109389754D-31 ! ELECTRON MASS / KG
REAL(8)   ,PARAMETER :: EBYAS   =1.6021773349D-19! ELECTRON CHARGE / [AMPERE*SECOND]
!***********************************************************************
!** DEFINE NOW THE FUNDAMENTAL SI UNITS IN TERMS OF ATOMIC UNITS      **
!***********************************************************************
REAL(8)   ,PARAMETER :: KG    =ME/MEBYKG         ! KILOGRAMM
REAL(8)   ,PARAMETER :: METER =CBYMBYS/HBARBYJS*HBAR/C/KG
REAL(8)   ,PARAMETER :: SECOND=CBYMBYS/C*METER
REAL(8)   ,PARAMETER :: AMPERE=E/(EBYAS*SECOND)
REAL(8)   ,PARAMETER :: KB    =3.166678911D-6    ! BOLTZMANN CONSTANT (1HARTREE/1KELVIN)
REAL(8)   ,PARAMETER :: MOL   =6.0221367D+23     ! MOLE
!***********************************************************************
!** OTHER INDEPENDENT PARAMETERS                                      **
!***********************************************************************
REAL(8)   ,PARAMETER :: CALBYJOULE=4.184D0       ! CALORIE IN UNITS OF JOULE
REAL(8)   ,PARAMETER :: MPBYME  =1836.152701D0   ! RATIO OF PROTON AND ELECTRON MASS
REAL(8)   ,PARAMETER :: MAGNETICMOMENTANOMALY=1.159652193D-3
                     ! THE MAGNETIC MOMENT OF AN ELECTRON IS 2*(1+...)
!***********************************************************************
!** EXPONENTS                                                         **
!***********************************************************************
REAL(8)   ,PARAMETER :: FEMTO=1.D-15
REAL(8)   ,PARAMETER :: PICO =1.D-12
REAL(8)   ,PARAMETER :: NANO =1.D-9
REAL(8)   ,PARAMETER :: MICRO=1.D-6
REAL(8)   ,PARAMETER :: MILLI=1.D-3
REAL(8)   ,PARAMETER :: KILO =1.D+3
REAL(8)   ,PARAMETER :: MEGA =1.D+6
REAL(8)   ,PARAMETER :: GIGA =1.D+9
REAL(8)   ,PARAMETER :: TERA =1.D+12
REAL(8)   ,PARAMETER :: PETA =1.D+15
!***********************************************************************
CONTAINS
!     ..................................................................
      SUBROUTINE CONSTANTS_INITIALIZE
!     ******************************************************************
!     **                                                              **
!     **  ENTER CONSTANTS AND DESCRIPTIONS INTO THE LIST              **
!     **  WHICH CAN BE SEARCHED                                       **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)              :: PI
      REAL(8)              :: MU0             ! PERMEABILITY OF VACUUM
      REAL(8)              :: EPSILON0        ! PERMITIVITY OF VACUUM
      REAL(8)              :: HARTREE         ! ATOMIC ENERGY UNIT
      REAL(8)              :: TAU0            ! ATOMIC TIME UNIT
      REAL(8)              :: ABOHR           ! ATOMIC LENGTH UNIT
      REAL(8)              :: EV              ! ELECTRON VOLT
      REAL(8)              :: DEBYE           ! DEBYE    
      REAL(8)              :: NEWTON          ! SI UNIT FOR FORCE
      REAL(8)              :: JOULE           ! SI UNIT FOR ENERGY
      REAL(8)              :: COULOMB         ! SI UNIT FOR CHARGE
      REAL(8)              :: VOLT            ! SI UNIT FOR ELCTROSTATIC POTENTIAL
      REAL(8)              :: TESLA           ! SI UNIT FOR THE MAGNETIC FIELD
      REAL(8)              :: PASCAL          ! SI UNIT FOR THE PRESSURE
      REAL(8)              :: KJBYMOL         ! KILOJOULE PER MOLE
      REAL(8)              :: KCALBYMOL       ! KILOCALORIE PER MOLE
      REAL(8)              :: CGSCHARGEUNIT   ! 
      REAL(8)              :: U               ! MASS UNIT FOR NUCLEI
!     ******************************************************************
      NC=0
!     
!     ==================================================================
!     == DIMENSION LESS CONSTANTS                                     ==
!     ==================================================================
      PI=4.D0*ATAN(1.D0) ;
      NC=NC+1;XXX(NC)=CONSTANT('PI',PI,'PI')
      NC=NC+1;XXX(NC)=CONSTANT('ALPHA',ALPHA,'FINE STRUCTURE CONSTANT')
      NC=NC+1;XXX(NC)=CONSTANT('MBYTE',2.D0**20,'MEGA BYTE')
      NC=NC+1;XXX(NC)=CONSTANT('DEGREE',2.D0*PI/360.D0,'DEGREE (UNIT FOR ANGLE)')
!     
!     ==================================================================
!     == FUNDAMENTAL ATOMIC UNITS                                     ==
!     ==================================================================
      NC=NC+1;XXX(NC)=CONSTANT('HBAR',HBAR,'A.U. FOR ANGULAR MOMENTUM')
      NC=NC+1;XXX(NC)=CONSTANT('E',E,'A.U. FOR CHARGE; ELEMENTARY CHARGE')
      NC=NC+1;XXX(NC)=CONSTANT('ME',ME,'A.U.FOR MASS = ELECTRON MASS')
      EPSILON0=1.D0/(4.D0*PI)
      NC=NC+1;XXX(NC)=CONSTANT('EPSILON0',EPSILON0,'PERMITIVITY OF VACUUM')
      NC=NC+1;XXX(NC)=CONSTANT('KELVIN',KELVIN,'TEMPERATURE UNIT (A.U. AND SI)')
!     
!     ==================================================================
!     == DERIVED ATOMIC UNITS                                         ==
!     ==================================================================
      TAU0          =HBAR**3/(ME*E**4)    ! ATOMIC TIME UNIT
      ABOHR         =(HBAR/E)**2/ME       ! BOHR RADIUS
      HARTREE       =ME*(E**2/HBAR)**2    ! HARTREE
      MU0           =1.D0/(EPSILON0*C**2) ! PERMEABILITY OF VACUUM
!     
      NC=NC+1;XXX(NC)=CONSTANT('ABOHR',ABOHR,'A.U. FOR LENGTH = BOHR RADIUS')
      NC=NC+1;XXX(NC)=CONSTANT('TAU0',TAU0,'A.U. FOR TIME')
      NC=NC+1;XXX(NC)=CONSTANT('HARTREE',HARTREE,'A.U. FOR ENERGY')
      NC=NC+1;XXX(NC)=CONSTANT('C',C,'SPEED OF LIGHT')
      NC=NC+1;XXX(NC)=CONSTANT('MU0',MU0,'PERMEABILITY OF VACCUM')
!
      NC=NC+1;XXX(NC)=CONSTANT('BOHRMAGNETON' &
     &                         ,E*HBAR/(2.D0*ME),'BOHR MAGNETON')
      NC=NC+1;XXX(NC)=CONSTANT('NUCLEARMAGNETON' &
     &                         ,E*HBAR/(2.D0*ME*MPBYME),'NUCLEAR MAGNETON')
      NC=NC+1;XXX(NC)=CONSTANT('GE' &
     &                         ,2.D0*(MAGNETICMOMENTANOMALY+1.D0),'ELECRON G-FACTOR')
!     
!     ==================================================================
!     == FUNDAMENTAL SI UNITS                                         ==
!     ==================================================================
      NC=NC+1;XXX(NC)=CONSTANT('KG',KG,'SI UNIT FOR MASS=KILOGRAMM')
      NC=NC+1;XXX(NC)=CONSTANT('SECOND',SECOND,'SI UNIT FOR TIME=SECOND')
      NC=NC+1;XXX(NC)=CONSTANT('METER',METER,'SI UNIT FOR LENGTH=METER')
      NC=NC+1;XXX(NC)=CONSTANT('AMPERE',AMPERE,'SI UNIT FOR ELECTRON CURRENT=AMPERE')
      NC=NC+1;XXX(NC)=CONSTANT('MOL',MOL,'MOLE')
      NC=NC+1;XXX(NC)=CONSTANT('KB',KB,'BOLTZMANN CONSTANT IN A.U.= HARTREE/KELVIN')
!
!     ==================================================================
!     == DERIVED SI UNITS                                             ==
!     ==================================================================
      NEWTON=KG*METER/SECOND**2 
      JOULE=NEWTON*METER       
      COULOMB=AMPERE*SECOND 
      TESLA=KG/(COULOMB*SECOND)
      VOLT=JOULE/COULOMB
      PASCAL=JOULE/METER**3
      NC=NC+1;XXX(NC)=CONSTANT('NEWTON',NEWTON,'SI UNIT FOR FORCE=NEWTON')
      NC=NC+1;XXX(NC)=CONSTANT('JOULE',JOULE,'SI UNIT FOR ENERGY=JOULE')
      NC=NC+1;XXX(NC)=CONSTANT('COULOMB',COULOMB,'SI UNIT FOR CHARGE=COULOMB')
      NC=NC+1;XXX(NC)=CONSTANT('VOLT',VOLT,'SI UNIT FOR ... =VOLT')
      NC=NC+1;XXX(NC)=CONSTANT('TESLA',TESLA,'SI UNIT FOR MAGNETIC FIELD =TESLA')
      NC=NC+1;XXX(NC)=CONSTANT('GAUSS',1.D-4*TESLA,'UNIT FOR MAGNETIC FIELD =GAUSS')
      NC=NC+1;XXX(NC)=CONSTANT('PASCAL',PASCAL,'SI UNIT FOR PRESSURE =PASCAL')
      NC=NC+1;XXX(NC)=CONSTANT('BAR',1.D+5*PASCAL,'UNIT FOR PRESSURE =BAR')
      NC=NC+1;XXX(NC)=CONSTANT('MINUTE',60.D0*SECOND,'UNIT FOR TIME =MINUTE')
      NC=NC+1;XXX(NC)=CONSTANT('HOUR',60.D0*60.D0*SECOND,'UNIT FOR TIME =HOUR')
      NC=NC+1;XXX(NC)=CONSTANT('DAY',24.D0*60.D0*60.D0*SECOND,'UNIT FOR TIME =DAY')
      NC=NC+1;XXX(NC)=CONSTANT('HERTZ',1.D0/SECOND,'UNIT FOR FREQUENCY =HERTZ')
      NC=NC+1;XXX(NC)=CONSTANT('OHM',VOLT/AMPERE,'UNIT FOR ELECTRICAL RESISTANCE =OHM')
!     
!     ==================================================================
!     == DERIVED UNITS                                                ==
!     ==================================================================
      EV=JOULE*E/COULOMB            ! ELECTRON VOLT
      KJBYMOL=1.D+3*JOULE/MOL ! KILOJOULE PER MOLE
      KCALBYMOL=CALBYJOULE*KJBYMOL  ! KILOCALORIE PER MOLE
      U=1.D-3*KG/MOL                ! ATOMIC MASS UNIT 
      NC=NC+1;XXX(NC)=CONSTANT('RY',0.5D0*HARTREE,'RYDBERG ENERGY UNIT')
      NC=NC+1;XXX(NC)=CONSTANT('EV',EV,'ELECTRON VOLT')
      NC=NC+1;XXX(NC)=CONSTANT('ANGSTROM',1.D-10*METER,'ANGSTROM=1.E-10 METER')
      NC=NC+1;XXX(NC)=CONSTANT('KJ/MOL',KJBYMOL,'KILOJOULE PER MOLE')
      NC=NC+1;XXX(NC)=CONSTANT('KCAL/MOL',KCALBYMOL,'KILOCALORIE PER MOLE')
      NC=NC+1;XXX(NC)=CONSTANT('U',U,'MASS UNIT')
!     
!     ==================================================================
!     == DERIVED UNITS IN THE CGS=ESU SYSTEM                          ==
!     ==================================================================
      CGSCHARGEUNIT=0.1D0/(C*SECOND/METER)*COULOMB
      DEBYE=CGSCHARGEUNIT*(1.D-2*METER)*1.D-18 ! CGSUNIT OF THE DIPOLE MOMENT 
      NC=NC+1;XXX(NC)=CONSTANT('DEBYE',DEBYE,'DEBYE= CGS UNIT OF DIPOLE MOMENT')
      NC=NC+1;XXX(NC)=CONSTANT('FEMTO',1.D-15,'FEMTO=1.D-15')
      NC=NC+1;XXX(NC)=CONSTANT('PICO',1.D-12,'PICO=1.D-12')
      NC=NC+1;XXX(NC)=CONSTANT('NANO',1.D-9,'NANO=1.D-9')
      NC=NC+1;XXX(NC)=CONSTANT('MICRO',1.D-6,'MICRO=1.D-6')
      NC=NC+1;XXX(NC)=CONSTANT('MILLI',1.D-3,'MILLI=1.D-3')
      NC=NC+1;XXX(NC)=CONSTANT('KILO',1.D+3,'KILO=1.D+3')
      NC=NC+1;XXX(NC)=CONSTANT('MEGA',1.D+6,'MEGA=1.D+6')
      NC=NC+1;XXX(NC)=CONSTANT('GIGA',1.D+9,'GIGA=1.D+9')
      NC=NC+1;XXX(NC)=CONSTANT('TERA',1.D+12,'TERA=1.D+12')
      NC=NC+1;XXX(NC)=CONSTANT('PETA',1.D+15,'PETA=1.D+15')
      TFIRST=.FALSE.
      RETURN
      END SUBROUTINE CONSTANTS_INITIALIZE
END MODULE CONSTANTS_MODULE
!
!     ..................................................................
      SUBROUTINE CONSTANTS$GET(ID_,VAL)
!     ******************************************************************
!     **                                                              **
!     **  RETURNS THE  VALUE OF A CONSTANT IDENTIFIED BY A STRING     **
!     **  STRING_='HELP' PROVIDES A LISTING OF THE CONSTANTS          **
!     **                                                              **
!     ******************************************************************
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: ID_ 
      REAL(8)      ,INTENT(OUT):: VAL
      INTEGER(4)               :: I
!     ******************************************************************
      IF(TFIRST) CALL CONSTANTS_INITIALIZE
      DO I=1,NC
        IF(ID_.EQ.XXX(I)%NAME) THEN
          VAL=XXX(I)%VALUE
          RETURN
        END IF
      ENDDO
      IF(TRIM(ID_).EQ.'HELP') THEN
        CALL CONSTANTS$LIST(6)
        RETURN
      END IF
      CALL ERROR$MSG('STRING NOT FOUND IN ROUTINE CONSTANTS')
      CALL ERROR$CHVAL('STRING=',ID_)
      CALL ERROR$STOP('CONSTANTS')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTANTS$REPORT(NFIL)
!     ******************************************************************
!     **                                                              **
!     **  PRINTS A  LIST OF THE CONSTANTS KNOWN TO THE OBJECT         **
!     **                                                              **
!     ******************************************************************
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL ! FORTRAN FILE UNIT
      REAL(8)                :: VALUE
      INTEGER(4)             :: I
!     ******************************************************************
      IF(TFIRST) CALL CONSTANTS_INITIALIZE
      WRITE(NFIL,FMT='("UNITS AND CONSTANTS IN ATOMIC HARTREE UNITS")')
      WRITE(NFIL,FMT='(A8,"=",A20," :",A)')'NAME','VALUE','DESCRIPTION'
      DO I=1,NC
        VALUE=XXX(I)%VALUE
        IF(ABS(VALUE).GE.0.1D0.AND.ABS(VALUE).LT.1000.D0) THEN
          WRITE(NFIL,FMT='(16("."),T1,A,T16,"=",F12.6,4X," :",A)') &
     &          TRIM(XXX(I)%NAME),XXX(I)%VALUE,TRIM(XXX(I)%DESCRIPTION)
        ELSE
          WRITE(NFIL,FMT='(16("."),T1,A,T16,"=",ES16.6," :",A)') &
     &          TRIM(XXX(I)%NAME),XXX(I)%VALUE,TRIM(XXX(I)%DESCRIPTION)
        END IF
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CONSTANTS$LISTCONVERSION(NFIL)
!     ******************************************************************
!     **                                                              **
!     **  PRINTS A  LIST OF CONVERSION FACTORS                        **
!     **                                                              **
!     ******************************************************************
      USE CONSTANTS_MODULE
      INTEGER(4)   ,INTENT(IN):: NFIL    ! FORTRAN FILE UNIT
      REAL(8)                 :: VALUE1
      REAL(8)                 :: VALUE2
      REAL(8)                 :: SVAR
!     ******************************************************************
      IF(TFIRST)CALL CONSTANTS_INITIALIZE
      WRITE(NFIL,FMT='(/"CONVERSION FACTORS AN CONSTANTS:"' &
     &              //'/"================================")') 
!
!     ==================================================================
!     ==   ATOMIC UNITS                                               ==
!     ==================================================================
      CALL TITLE('ATOMIC UNITS')
      CALL DUMP('HBAR',0.D0,'A.U.',1.D0)
      CALL DUMP('ME',0.D0,'A.U.',1.D0)
      CALL DUMP('HARTREE',0.D0,'A.U.',1.D0)
      CALL DUMP('ABOHR',0.D0,'A.U.',1.D0)
        CALL CONSTANTS('ALPHA',SVAR); VALUE1=1.D0/SVAR
      CALL DUMP('1/ALPHA',VALUE1,' ',1.D0)
      CALL DUMP('EPSILON0',0.D0,' ',1.D0)
      CALL DUMP('MU0',0.D0,' ',1.D0)
        CALL CONSTANTS('ABOHR',VALUE2)
        CALL CONSTANTS('TAU0',SVAR); VALUE2=VALUE2/SVAR
      CALL DUMP('C',0.D0,'ABOHR/TAU0',VALUE2)
      CALL DUMP('KB',0.D0,'HARTREE',0.D0)
!
!     ==================================================================
!     ==   SI UNITS                                                   ==
!     ==================================================================
      CALL TITLE('CONVERSION FROM SI TO ATOMIC UNITS')
      CALL DUMP('MOL',0.D0,' ',1.D0)
      CALL DUMP('METER',0.D0,'ABOHR',0.D0)
      CALL DUMP('SECOND',0.D0,'TAU0',0.D0)
      CALL DUMP('JOULE',0.D0,'HARTREE',0.D0)
      CALL DUMP('KG',0.D0,'ME',0.D0)
        CALL CONSTANTS('METER',VALUE2)
        CALL CONSTANTS('SECOND',SVAR); VALUE2=VALUE2/SVAR
      CALL DUMP('C',0.D0,'METER/SECOND',VALUE2)
!
!     ==================================================================
!     ==   CONVERSION OF SI UNITS TO ATOMIC UNITS                     ==
!     ==================================================================
      CALL TITLE('CONVERSION FROM ATOMIC TO SI UNITS')
      CALL DUMP('ABOHR',0.D0,'METER',0.D0)
      CALL DUMP('TAU0',0.D0,'SECOND',0.D0)
      CALL DUMP('ME',0.D0,'KG',0.D0)
      CALL DUMP('E',0.D0,'COULOMB',0.D0)
        CALL CONSTANTS('JOULE',VALUE2)
        CALL CONSTANTS('SECOND',SVAR); VALUE2=VALUE2*SVAR
      CALL DUMP('HBAR',0.D0,'JOULE*SECOND',VALUE2)
!
!     ==================================================================
!     ==   LENGTH CONVERSIONS                                         ==
!     ==================================================================
      CALL TITLE('LENGTH CONVERSIONS')
      CALL DUMP('ANGSTROM',0.D0,'ABOHR',0.D0)
      CALL DUMP('ABOHR',0.D0,'ANGSTROM',0.D0)
      CALL DUMP('ABOHR',0.D0,'METER',0.D0)
      CALL DUMP('METER',0.D0,'ABOHR',0.D0)
!
!     ==================================================================
!     ==   TIME CONVERSIONS                                           ==
!     ==================================================================
      CALL TITLE('TIME CONVERSIONS')
!
      CALL DUMP('TAU0',0.D0,'SECOND',0.D0)
        CALL CONSTANTS('SECOND',VALUE2)
        CALL CONSTANTS('FEMTO',SVAR); VALUE2=VALUE2*SVAR
      CALL DUMP('TAU0',0.D0,'FEMTO SECOND',VALUE2)
        CALL CONSTANTS('SECOND',VALUE1)
        CALL CONSTANTS('PICO',SVAR); VALUE1=VALUE1*SVAR
      CALL DUMP('PICO SECOND',VALUE1,'TAU0',0.D0)
!
!     ==================================================================
!     ==   ENERGY CONVERSIONS                                         ==
!     ==================================================================
      CALL TITLE('ENERGY CONVERSIONS')
        CALL CONSTANTS('KB',VALUE1)
        VALUE1=VALUE1*273.15
      CALL DUMP('273.15 KB',VALUE1,'EV',0.D0)
      CALL DUMP('EV',0.D0,'KJ/MOL',0.D0)
      CALL DUMP('HARTREE',0.D0,'KJ/MOL',0.D0)
      CALL DUMP('EV',0.D0,'KCAL/MOL',0.D0)
      CALL DUMP('HARTREE',0.D0,'KCAL/MOL',0.D0)
      CALL DUMP('KJ/MOL',0.D0,'EV',0.D0)
      CALL DUMP('KCAL/MOL',0.D0,'EV',0.D0)
      CALL DUMP('HARTREE',0.D0,'EV',0.D0)
      CALL DUMP('RY',0.D0,'HARTREE',0.D0)
!
!     ==================================================================
!     ==   MASS CONVERSIONS                                           ==
!     ==================================================================
      CALL TITLE('MASS CONVERSIONS')
      CALL DUMP('U',0.D0,'KG',0.D0)
      CALL DUMP('U',0.D0,'ME',0.D0)
!
!     ==================================================================
!     ==   CHARGE CONVERSIONS                                         ==
!     ==================================================================
      CALL TITLE('CHARGE CONVERSIONS')
!
      CALL DUMP('E',0.D0,'COULOMB',0.D0)
      CALL DUMP('COULOMB',0.D0,'E',0.D0)
!
!     ==================================================================
!     ==   OTHER CONVERSIONS                                         ==
!     ==================================================================
      CALL TITLE('OTHER CONVERSIONS')
        CALL CONSTANTS('ALPHA',VALUE1)
        CALL CONSTANTS('E',VALUE2)
        CALL CONSTANTS('ABOHR',SVAR); VALUE2=VALUE2*SVAR
      CALL DUMP('DEBYE',0.D0,'E*ABOHR',VALUE2)
      CALL DUMP('E*ABOHR',VALUE2,'DEBYE',0.D0)
        CALL CONSTANTS('HBAR',SVAR); VALUE1=0.5D0*SVAR
        CALL CONSTANTS('PI',SVAR); VALUE1=VALUE1*2.D0*SVAR
        CALL CONSTANTS('PICO',SVAR); VALUE1=VALUE1/SVAR
        CALL CONSTANTS('SECOND',SVAR); VALUE1=VALUE1/SVAR
      CALL DUMP('E_ZPV/T[PSEC]',VALUE1,'HARTREE',0.D0)
      CALL DUMP('E_ZPV/T[PSEC]',VALUE1,'EV',0.D0)
        CALL CONSTANTS('HBAR',SVAR); VALUE1=0.5D0*SVAR
        CALL CONSTANTS('PI',SVAR); VALUE1=VALUE1*2.D0*SVAR
        CALL CONSTANTS('C',SVAR); VALUE1=VALUE1*SVAR
        CALL CONSTANTS('METER',SVAR); VALUE1=VALUE1/(SVAR*1.D-2)
      CALL DUMP('E_ZPE/(1/LAMBDA)[CM**-1]',VALUE1,'EV',0.D0)
      CALL DUMP('E_ZPE/(1/LAMBBDA)[CM**-1]',VALUE1,'HARTREE',0.D0)
        CALL CONSTANTS('C',SVAR); VALUE1=1.D0/SVAR
        CALL CONSTANTS('PICO',SVAR); VALUE1=VALUE1/SVAR
        CALL CONSTANTS('SECOND',SVAR); VALUE1=VALUE1/SVAR
        CALL CONSTANTS('METER',SVAR); VALUE2=1.D0/(1.D-2*SVAR)
      CALL DUMP('(1/LAMBDA)',VALUE1,'CM**-1/T[PSEC]',VALUE2)
      RETURN
      CONTAINS
!       . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        SUBROUTINE TITLE(STRING)
        CHARACTER(LEN=*),INTENT(IN):: STRING
        WRITE(NFIL,FMT='(72("="),T10," ",A," ")')STRING
        RETURN
        END SUBROUTINE TITLE
!       . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        SUBROUTINE DUMP(NAME1,VALUE1_,NAME2,VALUE2_)
        CHARACTER(LEN=*),INTENT(IN) :: NAME1
        CHARACTER(LEN=*),INTENT(IN) :: NAME2
        REAL(8)          ,INTENT(IN) :: VALUE1_
        REAL(8)          ,INTENT(IN) :: VALUE2_
        REAL(8)                      :: VALUE1
        REAL(8)                      :: VALUE2
        REAL(8)                      :: VALUE
        IF(VALUE1_.EQ.0.D0) THEN
          CALL CONSTANTS(NAME1,VALUE1)
        ELSE
          VALUE1=VALUE1_
        END IF
        IF(VALUE2_.EQ.0.D0) THEN
          CALL CONSTANTS(NAME2,VALUE2)
        ELSE
          VALUE2=VALUE2_
        END IF
        VALUE=VALUE1/VALUE2
        IF(ABS(VALUE).GE.0.1D0.AND.ABS(VALUE).LT.1000.D0) THEN
          WRITE(NFIL,FMT='(24("."),T1,A,T24,"=",F12.6,4X," :",A)') &
     &          TRIM(NAME1),VALUE,TRIM(NAME2)
        ELSE
          WRITE(NFIL,FMT='(24("."),T1,A,T24,"=",ES16.6," :",A)') &
     &          TRIM(NAME1),VALUE,TRIM(NAME2)
        END IF
        RETURN
        END SUBROUTINE DUMP
      END     
