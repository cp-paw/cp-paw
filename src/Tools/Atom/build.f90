      PROGRAM MAIN
!     *******************************************************************
!!****h* 
!! NAME
!!   MAIN PROGRAM
!! PURPOSE
!!   write a
!! METHODS
!! USAGE
!! EXAMPLE
!! USES
!! NOTES
!! TODO
!! AUTHOR
!! CREATION DATE
!! COPYRIGHT c<year>-<year> by company/person
!!***
      USE STRINGS_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      logical(4),parameter :: tmakedirectories=.false.
      INTEGER(4)           :: IZ
      CHARACTER(2)         :: SY
      CHARACTER(32)        :: DIR,FILE
      CHARACTER(32)        :: STRING
      REAL(8)   ,PARAMETER :: FACTOR=0.75D0  !RC=FACTOR*RCOV
      REAL(8)   ,PARAMETER :: LAMBDA=6.D0  !RC=FACTOR*RCOV
      INTEGER(4)           :: NS,NP,ND,NF,L
      INTEGER(4)           :: NFIL
      INTEGER(4)           :: ISVAR
      REAL(8)              :: FS,FP,FD,FF
      REAL(8)              :: RCOV 
!     ********************************************
      DO IZ=1,106  !1,106
        CALL PERIODICTABLE$GET(IZ,'SYMBOL',SY)
!       ================================================================
!       ==  compose directory name  PBE/014SI                         ==
!       ==  and create it                                             ==
!       ================================================================
        WRITE(DIR,*)IZ
        DIR=ADJUSTL(DIR)
        IF(IZ.LT.10)THEN
          DIR='00'//DIR
        ELSE IF (IZ.LT.100) THEN
          DIR='0'//DIR
        END IF
        DIR=TRIM(DIR)//TRIM(SY)
        DIR='PBE/'//DIR
        STRING=-'MKDIR '//DIR
        PRINT*,'STRING ',STRING
        if(tmakedirectories) then
          CALL SYSTEM(STRING) 
          CYCLE
        end if
!
!       ===============================================================
!       == DEFINE FILE NAME  si_.75_6.0                              ==
!       ===============================================================
        WRITE(STRING,FMT='(F3.2)')FACTOR
        FILE=-TRIM(SY)//'_'//ADJUSTL(STRING)
        WRITE(STRING,FMT='(F3.1)')LAMBDA
        FILE=TRIM(FILE)//'_'//ADJUSTL(STRING)
        call system('rm '//TRIM(DIR)//'/'//TRIM(FILE)//-'.prot')
        call system('rm '//TRIM(DIR)//'/'//TRIM(FILE)//-'.out') 
        call system('rm '//TRIM(DIR)//'/'//TRIM(FILE)//-'.stp')
        FILE=TRIM(DIR)//'/'//TRIM(FILE)//-'.ACNTL'
print*,'file ',file
!       ===============================================================
!       == SPECIFY FILE                                              ==
!       ===============================================================
        CALL FILEHANDLER$SETROOT(DIR)
        CALL FILEHANDLER$SETFILE('ACNTL',.false.,FILE)
        CALL FILEHANDLER$SETSPECIFICATION('ACNTL','STATUS','UNKNOWN')
        CALL FILEHANDLER$SETSPECIFICATION('ACNTL','POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION('ACNTL','ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION('ACNTL','FORM','FORMATTED')
        CALL FILEHANDLER$UNIT('ACNTL',NFIL)
!
!       ===============================================================
!       == write input file                                          ==
!       ===============================================================
        CALL PERIODICTABLE$GET(IZ,'R(COV)',RCOV)
        CALL PERIODICTABLE$GET(IZ,'OCC(S)',ISVAR) ;FS=DBLE(ISVAR)
        CALL PERIODICTABLE$GET(IZ,'OCC(P)',ISVAR) ;FP=DBLE(ISVAR)
        CALL PERIODICTABLE$GET(IZ,'OCC(D)',ISVAR) ;FD=DBLE(ISVAR)
        CALL PERIODICTABLE$GET(IZ,'OCC(F)',ISVAR) ;FF=DBLE(ISVAR)
        WRITE(NFIL,FMT='("!ACNTL")')
        WRITE(NFIL,FMT='("  !GENERIC ELEMENT=''",A2,"'' RAUG=4.D0 RBOX=20.D0 !END")')SY
        WRITE(NFIL,FMT='("  !DFT     TYPE=10 !END")')
        WRITE(NFIL,FMT='("  !GRID    R1=1.056E-4 DEX=0.05 NR=250 !END")')
!
!       =============================================================================
!       == specify core electrons                                                  ==
!       =============================================================================
        WRITE(NFIL,FMT='("  !AECORE")')
        NS=1
        NP=0
        ND=0
        ND=0
        IF(IZ.GT.2) THEN   !HE-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')0,1,2.
          NS=2
          NP=2
        END IF
        IF(IZ.GT.10) THEN  !NE-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')0,2,2.
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')1,2,6.
          NS=3
          NP=3
        END IF
        IF(IZ.GT.18) THEN  !AR-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')0,3,2.
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')1,3,6.
          NS=4
          NP=4
          ND=3
        END IF
        IF(IZ.GT.31) THEN  !AR+d-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')2,3,10.
          if(iz.le.36) fd=0.d0
          NS=4
          NP=4
          ND=4
        END IF
        IF(IZ.GT.36) THEN  !KR-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')0,4,2.
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')1,4,6.
          NS=5
          NP=5
          ND=4
        END IF
        IF(IZ.GT.49) THEN  !Kr+d-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')2,4,10.
          if(iz.le.54) fd=0.d0
          NS=6
          NP=6
          ND=6
          NF=4
        END IF
        IF(IZ.GT.54) THEN  !XE-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')0,5,2.
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')1,5,6.
          NS=6
          NP=6
          ND=5
          NF=4
        END IF
        IF(IZ.GT.71) THEN  !Xe+f-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')3,4,14.
          if(iz.le.86) ff=0.d0
          NS=6
          NP=6
          ND=5
          NF=5
        END IF
        IF(IZ.GT.81) THEN  !Xe+f+d-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')2,5,10.
          if(iz.le.86) fd=0.d0
          NS=6
          NP=6
          ND=6
          NF=5
        END IF
        IF(IZ.GT.86) THEN  !RN-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')0,6,2.
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')1,6,6.
          NS=7
          NP=7
          ND=6
          NF=5
        END IF
        IF(IZ.GT.103) THEN  !RN+f-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')3,5,14.
          ff=0.d0
          NS=7
          NP=7
          ND=6
          NF=6
        END IF
        IF(IZ.GT.104) THEN  !RN+f+d-CORE
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')2,5,10.
          fd=0.d0
          NS=7
          NP=7
          ND=7
          NF=6
        END IF
        WRITE(NFIL,FMT='("  !END")')
!
!       =============================================================================
!       == specify valence electrons                                               ==
!       =============================================================================
        WRITE(NFIL,FMT='("  !VALENCE")')
        IF(FS.GT.0.D0) THEN
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')0,NS,FS
        END IF
        IF(FP.GT.0.D0) THEN
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')1,NP,FP
        END IF
        IF(FD.GT.0D0) THEN
          WRITE(NFIL,FMT='("    !STATE  L=",I2," N=",I2," F=",F5.1," !END")')2,ND,FD
        END IF
        WRITE(NFIL,FMT='("  !END")')
!
!       =============================================================================
!       == specify pseudo potential, pseudo core and compensation density          ==
!       =============================================================================
        WRITE(NFIL,FMT='("  !VTILDE  TYPE=",A," RC=",F5.3 &
     &                          ," POWER=",i2," !END ")') &
              '''POLYNOMIAL''',factor*RCOV-0.1,3
        WRITE(NFIL,FMT='("  !PSCORE  TYPE=''POLYNOMIAL'' RC=",F5.3 &
     &                       ," POWER=2   !END ")')factor*RCOV-0.1
        WRITE(NFIL,FMT='("  !COMPENSATE RC=",F5.3," !END ")')RCOV/4.d0
!
!       =============================================================================
!       == specify pseudo partal waves                                             ==
!       =============================================================================
        IF(NS.NE.0) THEN
           WRITE(NFIL,FMT='("  !WAVE L=0 N=",I2," PSTYPE=",A," RC=",F5.3 &
     &                     ," LAMBDA=6 !END")')NS,'"HBS"',FACTOR*RCOV
        END IF
        IF(NP.NE.0) THEN
          WRITE(NFIL,FMT='("  !WAVE L=1 N=",I2," PSTYPE=''HBS'' RC=",F5.3 &
    &                     ," LAMBDA=6 !END")')NP,FACTOR*RCOV
        END IF
        IF(ND.NE.0) THEN
          WRITE(NFIL,FMT='("  !WAVE L=2 N=",I2," PSTYPE=''HBS'' RC=",F5.3 &
    &                     ," LAMBDA=6 !END")')ND,FACTOR*RCOV
        END IF
        WRITE(NFIL,FMT='("  !WAVE L=0 E=0. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")') &
    &                  FACTOR*RCOV
        WRITE(NFIL,FMT='("  !WAVE L=1 E=0. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")') &
    &                  FACTOR*RCOV
        WRITE(NFIL,FMT='("  !WAVE L=2 E=0. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")') &
    &                  FACTOR*RCOV

        WRITE(NFIL,FMT='("  !WAVE L=0 E=1. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")') &
    &                  FACTOR*RCOV
        WRITE(NFIL,FMT='("  !WAVE L=1 E=1. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")') &
    &                  FACTOR*RCOV
        WRITE(NFIL,FMT='("  !WAVE L=2 E=1. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")') &
    &                  FACTOR*RCOV
!WRITE(NFIL,FMT='("  !WAVE L=3 E=0. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")')FACTOR*RCOV
!WRITE(NFIL,FMT='("  !WAVE L=3 E=1. PSTYPE=''HBS'' RC=",F5.3," LAMBDA=6 !END")')FACTOR*RCOV
        WRITE(NFIL,FMT='("!END")')
        WRITE(NFIL,FMT='("!EOB")')
        CALL FILEHANDLER$CLOSE('ACNTL')
!
        call system('paw_atom.x '//trim(file))
      ENDDO
      STOP
      END
