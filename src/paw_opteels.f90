MODULE OPTEELS_MODULE
TYPE OPTICAL_TYPE           !CALCULATE OPTICAL MATRIX ELEMENTS         
 CHARACTER(32)  :: ATOM
 CHARACTER(512) :: FILE
 CHARACTER(256) :: GROUP  ! GROUP NAME (MAKE A DEFAULT?!)
 LOGICAL(4)     :: ORIG   ! STATE ORIGIN CENTERED?
 LOGICAL(4)     :: QUAD   ! QUADRUPOLE CONTRIBUTION
 INTEGER(4)     :: IC     ! CORE SHELL NUMBER FOR EELS
 REAL(8)        :: ZI     ! CHANGE INTO NI&NF(GRID POINTS)
 REAL(8)        :: ZF     ! CHANGE INTO NI&NF(GRID POINTS)
 REAL(8)        :: EMAX 
END TYPE OPTICAL_TYPE
TYPE EELS_TYPE           !CALCULATE OPTICAL MATRIX ELEMENTS         
 CHARACTER(32)  :: ATOM
 CHARACTER(512) :: FILE
 LOGICAL(4)     :: QUAD   ! QUADRUPOLE CONTRIBUTION
 INTEGER(4)     :: IC     ! CORE SHELL NUMBER FOR EELS
 REAL(8)        :: INSTRUMENTALFWHM !INSTR. BROADENING (FULL-WIDTH HALF-MAX)
END TYPE EELS_TYPE
LOGICAL(4)                     :: TON=.FALSE.
INTEGER(4)                     :: NOPT=0
INTEGER(4)                     :: IOPT=0
TYPE(OPTICAL_TYPE),ALLOCATABLE :: OPTICAL(:)
INTEGER(4)                     :: NEELS=0
INTEGER(4)                     :: IEELS=0
TYPE(EELS_TYPE)   ,ALLOCATABLE :: EELS(:)
END MODULE OPTEELS_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS$SETR8(ID,VAL)
!      *************************************************************************
!      *************************************************************************
       USE OPTEELS_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       REAL(8)     ,INTENT(IN) :: VAL
!      *************************************************************************
!
!      =========================================================================
!      == SPECIFT CORE STATE FOR EELS CALCULATION                             ==
!      =========================================================================
       IF(ID.EQ.'BROADENING') THEN
         IF(IEELS.EQ.0) THEN
           CALL ERROR$MSG('IEELS HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETR8')
         END IF           
         EELS(IEELS)%INSTRUMENTALFWHM=VAL
!
!      =========================================================================
       ELSE IF(ID.EQ.'EMAX') THEN
         IF(IOPT.EQ.0) THEN
           CALL ERROR$MSG('IOPT HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETR8')
         END IF           
         OPTICAL(IOPT)%EMAX=VAL
!
!      =========================================================================
       ELSE IF(ID.EQ.'ZI') THEN
         IF(IOPT.EQ.0) THEN
           CALL ERROR$MSG('IOPT HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETR8')
         END IF           
         OPTICAL(IOPT)%ZI=VAL
!
!      =========================================================================
       ELSE IF(ID.EQ.'ZF') THEN
         IF(IOPT.EQ.0) THEN
           CALL ERROR$MSG('IOPT HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETR8')
         END IF           
         OPTICAL(IOPT)%ZF=VAL
!
!      ==  WRONG ID ============================================================
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OPTEELS$SETR8')
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS$SETI4(ID,VAL)
!      *************************************************************************
!      *************************************************************************
       USE OPTEELS_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(IN) :: VAL
!      *************************************************************************
!
!      =========================================================================
       IF(ID.EQ.'NOPT') THEN
         IF(NOPT.NE.0) THEN
           CALL ERROR$MSG('NOPT HAS ALREADY BEEN SET AND MUST BE SET ONLY ONCE')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETI4')
         END IF           
         NOPT=VAL
         ALLOCATE(OPTICAL(NOPT))
!
!      =========================================================================
       ELSE IF(ID.EQ.'IOPT') THEN
         IF(NOPT.EQ.0) THEN
           CALL ERROR$MSG('NOPT HAS NOT BEEN SET BEFORE SELECTING IOPT')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETI4')
         END IF           
         IF(IEELS.NE.0) THEN
           CALL ERROR$MSG('SET IEELS TO ZERO BEFORE SELECTING IOPT')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETI4')
         END IF
         IOPT=VAL
!
!      =========================================================================
!      == NEELS IS THE NUMBER OF EELS SPECTRA                                 ==
!      =========================================================================
       ELSE IF(ID.EQ.'NEELS') THEN
         IF(NEELS.NE.0) THEN
           CALL ERROR$MSG('NEELS HAS ALREADY BEEN SET')
           CALL ERROR$MSG('IT MUST BE SET ONLY ONCE')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETI4')
         END IF           
         NEELS=VAL
         ALLOCATE(EELS(NEELS))
         EELS(:)%INSTRUMENTALFWHM=0.D0
!
!      =========================================================================
!      == IEELS SELECTS A PARTICULAR INSTANCE OF THE EELS CALCULATION         ==
!      =========================================================================
       ELSE IF(ID.EQ.'IEELS') THEN
         IF(NEELS.EQ.0) THEN
           CALL ERROR$MSG('NEELS HAS NOT BEEN SET BEFORE SELECTING IOPT')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETI4')
         END IF           
         IF(IOPT.NE.0) THEN
           CALL ERROR$MSG('SET IOPT TO ZERO BEFORE SELECTING IEELS')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETI4')
         END IF
         IEELS=VAL
!
!      =========================================================================
!      == SPECIFT CORE STATE FOR EELS CALCULATION                             ==
!      =========================================================================
       ELSE IF(ID.EQ.'IC') THEN
         IF(IEELS.EQ.0) THEN
           CALL ERROR$MSG('IEELS HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETI4')
         END IF           
         EELS(IEELS)%IC=VAL
!
!      =========================================================================
!      ==  WRONG ID                                                           ==
!      =========================================================================
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OPTEELS$SETI4')
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS$SETL4(ID,VAL)
!      *************************************************************************
!      *************************************************************************
       USE OPTEELS_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       LOGICAL(4)  ,INTENT(IN) :: VAL
!      *************************************************************************
!
!      ==  ACTIVATE OR DEACTIVATE OPTEELS OBJECT ===============================
       IF(ID.EQ.'ON') THEN
         TON=VAL
!
!      =========================================================================
       ELSE IF(ID.EQ.'ORIG') THEN
         IF(IOPT.EQ.0) THEN
           CALL ERROR$MSG('IOPT HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETL4')
         END IF           
         OPTICAL(IOPT)%ORIG=VAL
!
!      ==  SELECT CALCULATION OF QUADRUPOLE ====================================
       ELSE IF(ID.EQ.'QUAD') THEN
         IF(IOPT.EQ.0) THEN
           CALL ERROR$MSG('IOPT HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETL4')
         END IF           
         OPTICAL(IOPT)%QUAD=VAL
!
!      ==  WRONG ID ============================================================
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OPTEELS$SETI4')
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS$SETCH(ID,VAL)
!      *************************************************************************
!      *************************************************************************
       USE OPTEELS_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       CHARACTER(*),INTENT(IN) :: VAL
!      *************************************************************************
!
!      ==  OUTPUT FILE =========================================================
       IF(ID.EQ.'FILE') THEN
         IF(IOPT.NE.0) THEN
           OPTICAL(IOPT)%FILE=VAL
         ELSE IF(IEELS.NE.0) THEN
           EELS(IEELS)%FILE=VAL
         ELSE
           CALL ERROR$MSG('NEITHER IOPT NOR IEELS HAVE BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETCH')
         END IF           
!
!      ==  SELECTS GROUP FOR OPTICAL TRANSITIONS ===============================
       ELSE IF(ID.EQ.'GROUP') THEN
         IF(IOPT.EQ.0) THEN
           CALL ERROR$MSG('IOPT HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETCH')
         END IF           
         OPTICAL(IOPT)%GROUP=VAL
!
!      ==  ATOM FOR EELS CALCULATION ===========================================
       ELSE IF(ID.EQ.'ATOM') THEN
         IF(IEELS.EQ.0) THEN
           CALL ERROR$MSG('IEELS HAS NOT BEEN SET')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$STOP('OPTEELS$SETCH')
         END IF           
         EELS(IEELS)%ATOM=VAL
!
!      ==  WRONG ID ============================================================
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OPTEELS$SETCH')
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS$PLOT()
!      *************************************************************************
!      **                                                                     **
!      **               MATTHE A. UIJTTEWAAL, 2011 (ADAPTED P.BLOECHL 2012)   **
!      *************************************************************************
       USE OPTEELS_MODULE, ONLY : TON,NOPT,OPTICAL,NEELS,EELS
       IMPLICIT NONE
       INTEGER(4)     :: I
       REAL(8)        :: EMAX
       CHARACTER(32)  :: ATOM
       LOGICAL(4)     :: ORIG
       LOGICAL(4)     :: QUAD
       CHARACTER(512) :: FILE
       CHARACTER(256) :: GROUP
       REAL(8)        :: ZI
       REAL(8)        :: INSTRUMENTALFWHM
       REAL(8)        :: ZF
       INTEGER(4)     :: IC
!      *************************************************************************
       IF(.NOT.TON) RETURN
                               CALL TRACE$PUSH('OPTEELS$PLOT')
!
!      =========================================================================
!      == CALCULATE EELS SPECTRA                                              ==
!      =========================================================================
       DO I=1,NEELS
         ATOM=EELS(I)%ATOM
         IC  =EELS(I)%IC
         FILE=EELS(I)%FILE
         INSTRUMENTALFWHM=EELS(I)%INSTRUMENTALFWHM
         CALL OPTEELS_EELS(FILE,ATOM,IC,INSTRUMENTALFWHM)
       ENDDO
!
!      =========================================================================
!      == CALCULATE OPTICAL MATRIX ELEMENTS                            --
!      =========================================================================
       IF(NOPT.NE.0.D0) THEN
         CALL ERROR$MSG('THIS OPTION HAS BEEN DISABLED')
         CALL ERROR$MSG('IT NEEDS TO BE CHECKED BEFORE GIVEN FREE')
         CALL ERROR$I4VAL('NOPT',NOPT)
         CALL ERROR$STOP('OPTEELS$PLOT')
       END IF
       DO I=1,NOPT 
         EMAX=OPTICAL(I)%EMAX
         FILE=TRIM(OPTICAL(I)%FILE)//'.OPT'
         ZI=OPTICAL(I)%ZI 
         ZF=OPTICAL(I)%ZF 
         ORIG=OPTICAL(I)%ORIG
         GROUP=OPTICAL(I)%GROUP
         CALL OPTEELS_OPTICAL(FILE,EMAX,ZI,ZF,QUAD,ORIG,GROUP)
       ENDDO
                               CALL TRACE$POP()
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS_GETGRIDIDX(ORIG,ZI,ZF,NI,NF)
!      *************************************************************************
!      **  CHANGE Z-COORDINATES INTO GRID INDICES                             **
!      **                                                                     **
!      **               MATTHE A. UIJTTEWAAL, 2011 (ADAPTED P.BLOECHL 2012)   **
!      *************************************************************************
       IMPLICIT NONE
       LOGICAL(4),INTENT(IN) :: ORIG
       REAL(8)   ,INTENT(IN) :: ZI,ZF
       INTEGER(4),INTENT(OUT):: NI,NF
       INTEGER(4)            :: NR3
       REAL(8)               :: RBAS(3,3),Z1,Z2
!      *************************************************************************
       CALL CELL$GETR8A('T(0)',9,RBAS) 
       CALL WAVES$GETI4('NR3',NR3) 
       IF(ABS(RBAS(3,3)).LT.1.E-8) THEN !NOT USEFULL
         NI=1
         NF=NR3
         RETURN
       ENDIF
       Z1=ZI
       IF(.NOT.ORIG.AND.ZI.LT.0) THEN
         Z1=.0
       ELSEIF(ZI.LT.-.5*RBAS(3,3)) THEN
         Z1=-.5*RBAS(3,3)
       ENDIF
       Z2=ZF
       IF(ORIG.AND.(ZF.GT..5*RBAS(3,3).OR.ZF.LT.(Z1+.0001))) THEN
         Z2=.5*RBAS(3,3)
       ELSEIF(ZF.GT.RBAS(3,3).OR.ZF.LT.(Z1+.0001)) THEN
         Z2=RBAS(3,3)
       ENDIF
       IF(ORIG.OR.Z1.GT.1.)Z1=Z1/RBAS(3,3)  !FRACTIN ONLY WHEN NOTORIG!
       NI=NINT(NR3*Z1)+1 
       IF(ORIG.OR.Z2.GT.1.)Z2=Z2/RBAS(3,3)
       NF=INT(NR3*Z2) 
!WRITE(*,FMT='("ZI,NI..",F5.1,I4,F5.1,2I4)')ZI,NI,ZF,NF,NR3
       RETURN
       END 

!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS_TRANSFORM(EXCEN,WGHT,OPT,OSC)
!      *************************************************************************
!      **  CHANGE OPT FROM X,Y,Z TO R,THETA,PHI(UNITS OF PI), % 
!      **  UNITS TO ANGST,EV & MAKE DIM.LESS OSCIL. STRENGTH  % 
!      **                                                                     **
!      **               MATTHE A. UIJTTEWAAL, 2011 (ADAPTED P.BLOECHL 2012)   **
!      *************************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN)    :: WGHT
       REAL(8),INTENT(INOUT) :: EXCEN
       REAL(8),INTENT(INOUT) :: OPT(3)
       REAL(8),INTENT(INOUT) :: OSC
       REAL(8),PARAMETER     :: PI=4.D0*ATAN(1.D0)
       REAL(8)               :: ANGST,EV,TMP(3)
!      *************************************************************************
       TMP(:)=(/ SQRT(SUM(OPT(:)**2)), .0D0, .0D0 /)
! THE WEIGHT MUST NOT BE SQUARED!
       OSC=2./3.*EXCEN*(WGHT*TMP(1))**2 !DIMENSIONLESS OSCILLATOR STRENGTH
       IF(TMP(1).GT.1.E-8) TMP(2)=ACOS(OPT(3) / TMP(1))/PI
       IF(OPT(2).GT.1.E-8) TMP(3)=ASIN(OPT(2) / SQRT(OPT(1)**2+OPT(2)**2))/PI
       CALL CONSTANTS$GET('ANGSTROM',ANGST)
       CALL CONSTANTS$GET('EV',EV)
       OPT(:)=(/ TMP(1)/ANGST, TMP(2), TMP(3) /) !IN E-ANGSTROEM,PI
       EXCEN=EXCEN/EV !IN EV
       RETURN
       END SUBROUTINE OPTEELS_TRANSFORM
     
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS_OPTICAL(FILE,EMAX,ZI,ZF,QUAD,ORIG,GROUP)
!      *************************************************************************
!      **  CALCULATE ALL OPTICAL MATRIX ELEMENTS AND WRITE TO FILE            **
!      **                                                                     **
!      **               MATTHE A. UIJTTEWAAL, 2011 (ADAPTED P.BLOECHL 2012)   **
!      *************************************************************************
       USE MPE_MODULE
       IMPLICIT NONE
       LOGICAL(4)    ,INTENT(IN) :: QUAD !QUADRUPOLE CONTR?
       LOGICAL(4)    ,INTENT(IN) :: ORIG !STATE ORIGIN CENTERED?
       CHARACTER(512),INTENT(IN) :: FILE*512
       CHARACTER(256),INTENT(IN) :: GROUP
       REAL(8)       ,INTENT(IN) :: EMAX
       REAL(8)       ,INTENT(IN) :: ZI,ZF !Z-COORDINATES OF ELEMENT RESTRICTION
       LOGICAL(4)                :: SAVETRAWSTATES
       LOGICAL(4)                :: TKGROUP 
       LOGICAL(4)    ,ALLOCATABLE:: TGROUP(:)
       CHARACTER(132)            :: TITLE
       INTEGER(4)                :: IB,JB,NB,NI,NF,IKPT,NKPT,ISPIN,NSPIN
       INTEGER(4)                :: NFIL,IAT,NAT,LMNXX
       INTEGER(4)                :: NTASKS_K,THISTASK_K
       REAL(8)       ,PARAMETER  :: PI=4.*ATAN(1.)
       REAL(8)                   :: RBAS(3,3)
       REAL(8)                   :: WGHT
       REAL(8)                   :: OPT(3)
       REAL(8)                   :: OPAT(3)
       REAL(8)                   :: OPGR(3)
       REAL(8)                   :: OPGRAT(3)
       REAL(8)                   :: OPQ(5)
       REAL(8)                   :: OPQAT(5)
       REAL(8)                   :: NRM,NRMAT
       REAL(8)                   :: OSC,OSCG,OSCQ
       REAL(8)                   :: EXCEN
       REAL(8)     ,ALLOCATABLE :: WKPT(:)
       REAL(8)     ,ALLOCATABLE :: OCC(:,:,:)
       REAL(8)     ,ALLOCATABLE :: POS(:,:)
       REAL(8)     ,ALLOCATABLE :: EIGVAL(:)
!      *************************************************************************
                        CALL TRACE$PUSH('OPTEELS_OPTICAL')
!      -- GET GRID INDICES --
       CALL OPTEELS_GETGRIDIDX(ORIG,ZI,ZF,NI,NF)

!      =========================================================================
!      == GET OCCUPATIONS + KP-WEIGHTS OF STATES                              ==
!      =========================================================================
!      ==LATER RESET ORIGINAL WAVE OBJECT
       CALL WAVES$GETL4('RAWSTATES',SAVETRAWSTATES) 
       CALL WAVES$SETL4('RAWSTATES',.FALSE.) 
       CALL WAVES$GETI4('NKPT',NKPT) 
       CALL WAVES$GETI4('NSPIN',NSPIN) 
       CALL DYNOCC$GETI4('NB',NB) 
       ALLOCATE(WKPT(NKPT),OCC(NB,NKPT,NSPIN),EIGVAL(NB))
       CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT) 
       CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC) !0<OCC<2*WKPT(IKPT)/NSPIN
!
!      =========================================================================
!      == GET ATOM POSITIONS, GROUP MEMBERS                                   ==
!      =========================================================================
       CALL ATOMLIST$NATOM(NAT)
       ALLOCATE(POS(3,NAT),TGROUP(NAT))
       DO IAT=1,NAT 
         CALL ATOMLIST$GETR8A('R(0)',IAT,3,POS(:,IAT)) !CARTHESIAN
       ENDDO 
       CALL GROUPLIST$MEMBERS(GROUP,NAT,TGROUP)
!PRINT*,'TGROUP',TGROUP
!
!      =========================================================================
!      == OPEN FILE & WRITE HEADER                                            ==
!      =========================================================================
       CALL FILEHANDLER$SETFILE('OPTICAL',.FALSE.,FILE)
       CALL FILEHANDLER$SETSPECIFICATION('OPTICAL','FORM','FORMATTED')
       CALL FILEHANDLER$UNIT('OPTICAL',NFIL) !OPENS FILE # NFIL AT ITS END
       TITLE=" BND1/2 KP SP WGHT  EXC(EV) F(-)WGHT |OPT|(EA) THETA,PHI(PI) NRM (FGR,FQU) ORIG?"
       WRITE(NFIL,FMT='(A84,L2)')TITLE, ORIG
!
!      =========================================================================
!      == CALC & WRITE ALL OPTICAL MATRIX ELEMENTS                            ==
!      =========================================================================
!      ==  POSSIBLE SPEED-UP: PASS-ON WAV1R&WAV1I TO OPTPW + ALL DIMS 
!      ==                                                        OF ALLOARRAYS?
!      == ALSO PASS ATOMIC AE/PSPHI(ISP) + DIMS ALLOARRAYS TO OPTAT ?
       DO IKPT=1,NKPT !PARALLELLIZE!
         CALL WAVES$SETI4('IKPT',IKPT)
         DO ISPIN=1,NSPIN
           CALL WAVES$SETI4('ISPIN',ISPIN)
           CALL WAVES$GETR8A('EIGVAL',NB,EIGVAL) 
           DO IB=1,NB 
!            == SELECT STATE AND CHECK IF PRESENT ON CURRENT NODE ==============
             CALL WAVES$SETI4('IB',IB)
!            CALL PLANEWAVE$SELECT(GSET%ID) !???,FROM WAVES
             CALL WAVES$STATESELECTED(TKGROUP) !IB SUFFICES
             CALL MPE$QUERY('K',NTASKS_K,THISTASK_K) !TWICE OUTPUT
! CALL PE$QUERY('MONOMER',NTASKS,THISTASK) !NOT RELEVANT
! SENDTSK=0 !1??
             IF(THISTASK_K.EQ.1) THEN !ONLY ON PROC 1?!
               DO JB=IB+1,NB 
                 OPGR=.0D0
                 WGHT=OCC(IB,IKPT,ISPIN) &
      &              * (1.D0-OCC(JB,IKPT,ISPIN)*NSPIN/WKPT(IKPT)/2.) 
                 EXCEN=EIGVAL(JB)-EIGVAL(IB)
!                == CALC ELEMENT IF STATE ON NODE, WEIGHT > 1/20 ===============
!                ==                             OF MAX & EXCEN NOT TOO LARGE  ==
                 IF(TKGROUP.AND.WGHT.GT.(REAL(WKPT(IKPT))/NSPIN) &
      &                    .AND.EXCEN.LT.EMAX) THEN 
                   CALL OPTEELS_OPTPW(IB,JB,IKPT,ISPIN,ORIG,NI,NF,OPT,NRM)
                   DO IAT=1,NAT
                     IF(POS(3,IAT).GT.ZI.AND.POS(3,IAT).LT.ZF) THEN 
                                                    !NO SHIFT/FRACTION ALLOWED!
                       CALL OPTEELS_OPTAT(IAT,IB,JB,POS,OPAT,OPGRAT,NRMAT)
                       OPT(:)=OPT(:)+OPAT(:)
                       IF(TGROUP(IAT))OPGR(:)=OPGR+OPGRAT
                       NRM=NRM+NRMAT
                     ENDIF
                   ENDDO
!         
!                  =============================================================
!                  == TRANSFORM VARIABLES & WRITE TO FILE                     ==
!                  =============================================================
                   CALL OPTEELS_TRANSFORM(EXCEN,WGHT,OPT,OSC)
                   NRMAT=EIGVAL(JB)-EIGVAL(IB) !DUMMY, E SCALED ONLY ONCE!
                   CALL OPTEELS_TRANSFORM(NRMAT,WGHT,OPGR,OSCG) 
                   WRITE(NFIL,FMT='(4I4,F6.3,F6.2,E11.4,F9.5,3F7.3,2E11.4)') &
     &                         IB,JB,IKPT,ISPIN,WGHT,EXCEN,OSC,OPT,NRM,OSCG,OSCQ
                 ENDIF
               ENDDO
! CALL MPE$COMBINE('MONOMER','+',SENDTASK) !SENDTASK IS CHANGED, TOO LARGE!!!
! NO COMBINE: ONLY 1 NODE!? COMBINE OVER EVERYTHING, NOT JUST WITHIN K?!
! SENDTASK=THISTASK
! CALL MPE$SENDRECEIVE('MONOMER',SENDTASK,1,NR1) 
! CALL MPE$SENDRECEIVE('MONOMER',SENDTASK,1,NR2) 
! CALL MPE$SENDRECEIVE('MONOMER',SENDTASK,1,NR3) 
! CALL MPE$SENDRECEIVE('MONOMER',SENDTASK,1,OPT)
! CALL MPE$SENDRECEIVE('MONOMER',SENDTASK,1,NRM)
             ENDIF !PARALLEL
           ENDDO 
         ENDDO
       ENDDO
!
!      =========================================================================
!      == CLOSE DOWN                                                          ==
!      =========================================================================
       CALL FILEHANDLER$CLOSE('OPTICAL')
       CALL WAVES$SETL4('RAWSTATES',SAVETRAWSTATES)
       DEALLOCATE(OCC,WKPT,POS,TGROUP,EIGVAL)
                              CALL TRACE$POP
       RETURN
       END SUBROUTINE OPTEELS_OPTICAL
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS_OPTPW(IB1,IB2,IKPT,ISPIN,ORIG,NI,NF,OPT,NRM)
!      *************************************************************************
!      **  CALCULATE PLANE WAVE PART OF OPTICAL ELEMENT FROM IB1 TO IB2       **
!      **  FOR IKPT,ISPIN RESTRICTED TO NI,NF IN ORIG?                        **
!      **                                                                     **
!      **               MATTHE A. UIJTTEWAAL, 2011 (ADAPTED P.BLOECHL 2012)   **
!      *************************************************************************
       USE MPE_MODULE
       IMPLICIT NONE
       LOGICAL(4),INTENT(IN) :: ORIG !STATE ORIGIN CENTERED?
       INTEGER(4),INTENT(IN) :: IB1,IB2
       INTEGER(4),INTENT(IN) :: IKPT, ISPIN
       INTEGER(4),INTENT(IN) :: NI,NF !BNDS, KP, SPN, GRIDPTS
       REAL(8)   ,INTENT(OUT):: OPT(3),NRM
       INTEGER(4)            :: NR1,NR1L,NR2,NR3,NNRL,NR1START, IR,IR1,IR2,IR3
       REAL(8)               :: GBAS(3:3)
       REAL(8)               :: RBAS(3,3)
       REAL(8)               :: VOL, F1,F2,F3
       REAL(8)               :: RVEC(3), OP(3) 
       REAL(8)   ,ALLOCATABLE:: WAVE(:,:,:)
       REAL(8)   ,ALLOCATABLE:: WAV1R(:,:,:),WAV2R(:,:,:)
       REAL(8)   ,ALLOCATABLE:: WAV1I(:,:,:),WAV2I(:,:,:)
!      *************************************************************************
                             CALL TRACE$PUSH('OPTEELS_OPTPW')
!      =========================================================================
!      == GET WAVE FUNCTION PARAMETERS                                        ==
!      =========================================================================
       CALL WAVES$GETI4('NR1',NR1)
       CALL WAVES$GETI4('NR1L',NR1L)
       CALL WAVES$GETI4('NR2',NR2)
       CALL WAVES$GETI4('NR3',NR3)
       NNRL=NR1L*NR2*NR3
      
!      =========================================================================
!      -- GET WAV1R,WAV1I                                                     ==
!      =========================================================================
       ALLOCATE(WAV1R(NR1,NR2,NR3),WAV1I(NR1,NR2,NR3),WAVE(NR1L,NR2,NR3))
       CALL WAVES$SETL4('TIM',.TRUE.)
       CALL WAVES$GETR8A('PSPSI',NNRL,WAVE)
       CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WAV1I)!PAR?!
       CALL WAVES$SETL4('TIM',.FALSE.)
       CALL WAVES$GETR8A('PSPSI',NNRL,WAVE)
       CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WAV1R) 

!      =========================================================================
!      -- GET WAV2R,WAV2I                                                     ==
!      =========================================================================
       ALLOCATE(WAV2R(NR1,NR2,NR3),WAV2I(NR1,NR2,NR3))
       CALL WAVES$SETI4('IB',IB2)
       CALL WAVES$SETL4('TIM',.TRUE.)
       CALL WAVES$GETR8A('PSPSI',NNRL,WAVE) 
       CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WAV2I)
       CALL WAVES$SETL4('TIM',.FALSE.)
       CALL WAVES$GETR8A('PSPSI',NNRL,WAVE) 
       CALL PLANEWAVE$RSPACECOLLECTR8(NR1L*NR2*NR3,WAVE,NR1*NR2*NR3,WAV2R)
       DEALLOCATE(WAVE)

!      =========================================================================
!      = GET CELL VECTORS, VOL                                                ==
!      =========================================================================
       CALL CELL$GETR8A('T(0)',9,RBAS)
       CALL GBASS(RBAS,GBAS,VOL) 
       VOL=ABS(VOL) !IN CASE RBAS ORDERED DIFFERENTLY
!       
!      =========================================================================
!      == CALC PLANE WAVE ELEMENT & NRM 2ND STATE                             ==
!      =========================================================================
       NRM=0.D0
       OPT=0.D0
       DO IR1=1,NR1 
         F1=0.D0
         DO IR2=1,NR2
           F2=0.D0
           OP=0.D0
           DO IR3=NI,NF !RESTRICTED IN Z-DIRECTION
             IF(IR3.LE.0) THEN 
               IR=IR3+NR3
             ELSE 
               IR=IR3 
             ENDIF
             F3 = WAV1R(IR1,IR2,IR)*WAV2R(IR1,IR2,IR) &
    &           + WAV1I(IR1,IR2,IR)*WAV2I(IR1,IR2,IR)
             RVEC(:)=(/ (IR1*1.D0)/NR1, (IR2*1.D0)/NR2, (IR3*1.D0)/NR3 /)
             IF(ORIG) THEN 
               RVEC(:)=RVEC(:)+.5D0
               RVEC(:)=(/MOD(RVEC(1),1.D0),MOD(RVEC(2),1.D0),MOD(RVEC(3),1.D0)/)
               RVEC(:)=RVEC(:)-0.5D0
             ENDIF
             OP(:)=OP(:) + F3*MATMUL(RBAS,RVEC)
             F2=F2 + WAV2R(IR1,IR2,IR)**2 + WAV2I(IR1,IR2,IR)**2
           ENDDO 
           F1=F1 + F2
           OPT(:)=OPT(:) + OP(:)
         ENDDO
         NRM=NRM + F1
       ENDDO
       OPT(:)=OPT(:) * VOL/NR1/NR2/NR3  !PROPER NORMALISATION
       NRM=NRM * VOL/NR1/NR2/NR3  
       PRINT*,'OPTPW,NRM',SQRT(SUM(OPT(:)**2)),NRM,IB1,IB2,IKPT,ISPIN
       DEALLOCATE(WAV1R,WAV1I,WAV2R,WAV2I)
                              CALL TRACE$POP
       RETURN
! -- CHECK CONSISTENCY NRM FROM GSPACE --
! COMPLEX(8) :: CSVAR1,CSVAR2
! COMPLEX(8),ALLOCATABLE :: OVERLAP(:,:),AUXMAT(:,:) !GNORM
! WRITE(*,FMT='(3I3," PWNRM",F10.5)')IB2,IKPT,ISPIN,NRM 
! ALLOCATE(OVERLAP(NB,NB),AUXMAT(NB,NB)) 
! WRITE(*,*)"NGL,NDIM,NBH,NB",GSET%NGL,NDIM,THIS%NBH,NB
! CALL WAVES_1COVERLAP(MAP,NDIM,(THIS%NBH),NB,(MAP%NPRO),(THIS%PROJ),(THIS%PROJ),AUXMAT) !PAR?
! CALL WAVES_OVERLAP(.TRUE.,(GSET%NGL),NDIM,(THIS%NBH),NB,(THIS%PSI0),(THIS%PSI0),OVERLAP)!PAR?
! CSVAR1=(.0D0,.0D0)
! CSVAR2=(.0D0,.0D0)
! DO IB=1,NB
!   DO IBA=1,NB
!     CSVAR1=CSVAR1+OVERLAP(IB,IBA)*CONJG(THIS%EIGVEC(IB,IB2))*THIS%EIGVEC(IBA,IB2)
!     CSVAR2=CSVAR2+AUXMAT(IB,IBA)*CONJG(THIS%EIGVEC(IB,IB2))*THIS%EIGVEC(IBA,IB2)
!   ENDDO
! ENDDO
! WRITE(*,FMT='("NRMPWGOPT,+1C ",3I3,2F10.5)')IB2,IKPT,ISPIN,REAL(CSVAR1),REAL(CSVAR1+CSVAR2)
! DEALLOCATE(OVERLAP,AUXMAT)
       END SUBROUTINE OPTEELS_OPTPW
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS_OPTAT(IAT,IB1,IB2,POS,OPT,OPTG,NRM)
!      *************************************************************************
!      ** CALC 1-C CONTRIBUTION OF THIS ATOM TO OPTICAL ELEMENT               **
!      **                                                                     **
!      **               MATTHE A. UIJTTEWAAL, 2011 (ADAPTED P.BLOECHL 2012)   **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: IAT,IB1,IB2
       REAL(8)   ,INTENT(IN) :: POS(3)
       REAL(8)   ,INTENT(OUT):: OPT(3),OPTG(3),NRM
       INTEGER(4)            :: ISP,GID,NR
       INTEGER(4)            :: LNX,LMNXX,LMN1,LMN2,LN1,LN2,L1,L2
       INTEGER(4)            :: M1,LM1,M2,LM2 
       INTEGER(4),ALLOCATABLE:: LOX(:) !(LNX)
       REAL(8)   ,PARAMETER  :: SQPI=SQRT(16.*ATAN(1.0)/3.)  !SQRT(4*PI/3)
       REAL(8)               :: RBOX, F2,O2(3),OG2(3)
       REAL(8)               :: AEINT,AERINT,PSINT,PSRINT, CG(3), PRO
       REAL(8)   ,ALLOCATABLE:: AEPHI(:,:),PSPHI(:,:)
       REAL(8)   ,ALLOCATABLE:: PROJ1(:),PROJ2(:)
       REAL(8)   ,ALLOCATABLE:: PROI1(:),PROI2(:)
!      *************************************************************************
                          CALL TRACE$PUSH('OPTEELS_OPTAT')
!      =========================================================================
!      == GET PARTIAL WAVES (REAL,RADIAL)                                     ==
!      =========================================================================
       CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
       CALL SETUP$ISELECT(ISP)
       CALL SETUP$GETR8('RBOX',RBOX)
       CALL SETUP$GETI4('LNX',LNX)
       ALLOCATE(LOX(LNX))
       CALL SETUP$GETI4A('LOX',LNX,LOX)
       CALL SETUP$GETI4('GID',GID)
       CALL RADIAL$GETI4(GID,'NR',NR) !=SETUP$GETI4('NR',NR), BUT GID ALSO NEEDED LATER
       CALL SETUP$GETI4('LMNXX',LMNXX)
       ALLOCATE(AEPHI(NR,LNX))
       ALLOCATE(PSPHI(NR,LNX))
       ALLOCATE(PROJ1(LMNXX),PROJ2(LMNXX),PROI1(LMNXX),PROI2(LMNXX))
       CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI) !FOR ISP
       CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI) 
       CALL SETUP$UNSELECT()
       
!      =========================================================================
!      == GET OCCUPATIONS                                                     ==
!      =========================================================================
       CALL WAVES$SETI4('IAT',IAT) 
       CALL WAVES$SETI4('IB',IB1)
       CALL WAVES$SETL4('TIM',.FALSE.) 
       CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROJ1)
       CALL WAVES$SETL4('TIM',.TRUE.) 
       CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROI1) !IMAGINARY PART
       CALL WAVES$SETI4('IB',IB2)
       CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROI2) 
       CALL WAVES$SETL4('TIM',.FALSE.) !
       CALL WAVES$GETR8A('<PSPSI|PRO>',LMNXX,PROJ2)
                
!      =========================================================================
!      == CALC LOCAL PART OF OPTICAL MATRIX ELEMENT                           ==
!      =========================================================================
!      %  THE FIRST REAL SPHERICAL HARMONICS ARE:   % 
!      %     YLM(1)=SQRT( 1/( 4*PI))    * 1         % 
!      %     YLM(2)=SQRT( 3/( 4*PI))    * X / R     % 
!      %     YLM(3)=SQRT( 3/( 4*PI))    * Z / R     % 
!      %     YLM(4)=SQRT( 3/( 4*PI))    * Y / R     % 
!      -- GET NUMERICAL CORRECT 1C OVERLAP --
!      REAL(8), ALLOCATABLE :: DOVER(:,:) !(NR,LNX) (REAL FUNCT)
!      ALLOCATE(DOVER(LNX,LNX))
!      CALL SETUP$1COVERLAP(ISP,LNX,DOVER) 
       OPT=.0D0
       OPTG=.0D0
       NRM=.0D0
       LMN1=0
       DO LN1=1,LNX
         L1=LOX(LN1) 
         LMN2=0
         DO LN2=1,LNX
           O2=0.D0
           OG2=0.D0
           F2=0.D0
           L2=LOX(LN2) 
           IF(ABS(L1-L2).LE.1) THEN
             IF(L1.EQ.L2) THEN !SPEED IT UP
               CALL RADIAL$MOMENT(GID,NR,0 &
     &                           ,(AEPHI(:,LN1)*AEPHI(:,LN2)),RBOX,AEINT)
               CALL RADIAL$MOMENT(GID,NR,0 &
     &                           ,(PSPHI(:,LN1)*PSPHI(:,LN2)),RBOX,PSINT)
             ENDIF
             CALL RADIAL$MOMENT(GID,NR,1 &
     &                         ,(AEPHI(:,LN1)*AEPHI(:,LN2)),RBOX,AERINT) 
             CALL RADIAL$MOMENT(GID,NR,1 &
     &                         ,(PSPHI(:,LN1)*PSPHI(:,LN2)),RBOX,PSRINT) 
           ENDIF
           DO M1=1,2*L1+1
             LM1=L1**2+M1
             DO M2=1,2*L2+1
               LM2=L2**2+M2
               CALL SPHERICAL$GAUNT(LM1,2,LM2,CG(1)) !X/R
               CALL SPHERICAL$GAUNT(LM1,4,LM2,CG(2)) !Y/R
               CALL SPHERICAL$GAUNT(LM1,3,LM2,CG(3)) !Z/R
               PRO= PROJ1(LMN1+M1)*PROJ2(LMN2+M2) &
     &            + PROI1(LMN1+M1)*PROI2(LMN2+M2)
               O2(:)=O2(:) + SQPI*(AERINT-PSRINT)*PRO* CG(:)
               OG2(:)=OG2(:) + SQPI*AERINT*PRO* CG(:)
               IF((M1.EQ.M2).AND.(L1.EQ.L2)) THEN
                 O2(:)=O2(:) + (AEINT-PSINT)*PRO* POS(:)
                 OG2(:)=OG2(:) +  AEINT*PRO* POS(:)
!                O2=O2 + POS(:) * DOVER(LN1,LN2)*PRO 
                 PRO= PROJ2(LMN1+M1)*PROJ2(LMN2+M2) &
     &              + PROI2(LMN1+M1)*PROI2(LMN2+M2) !FORNRM
                 F2=F2 + (AEINT-PSINT)*PRO 
!                F2=F2 + DOVER(LN1,LN2)*PRO 
               ENDIF
             ENDDO
           ENDDO
           LMN2=LMN2+2*L2+1
           OPT(:)=OPT(:)+O2(:)
           OPTG(:)=OPTG(:)+OG2(:)
           NRM=NRM+F2
         ENDDO
         LMN1=LMN1+2*L1+1
       ENDDO
       PRINT*,'OPTAT,GR',SQRT(SUM(OPT(:)**2)),SQRT(SUM(OPTG(:)**2)),IAT,IB1,IB2
       DEALLOCATE(AEPHI,PSPHI,PROJ1,PROJ2,PROI1,PROI2,LOX) 
                          CALL TRACE$POP
       RETURN
       END SUBROUTINE OPTEELS_OPTAT

!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OPTEELS_EELS(FILE,ATOM,IC,INSTRUMENTALFWHM)
!      *************************************************************************
!      **   CALCULATE ALL EELS MATRIX ELEMENTS AND WRITE TO FILE              **
!      **                                                                     **
!      **               MATTHE A. UIJTTEWAAL, 2011 (ADAPTED P.BLOECHL 2012)   **
!      *************************************************************************
       USE MPE_MODULE
!       USE CONSTANTS_MODULE
       IMPLICIT NONE
       CHARACTER(512),INTENT(IN) :: FILE
       CHARACTER(32) ,INTENT(IN) :: ATOM ! ATOM NAME
       INTEGER(4)    ,INTENT(IN) :: IC 
       REAL(8)       ,INTENT(IN) :: INSTRUMENTALFWHM  ! INSTRUMENTAL BROADENING 
       INTEGER(4)    ,PARAMETER  :: NE=1000
       COMPLEX(8)    ,PARAMETER  :: CI=(0.D0,1.D0)
       REAL(8)       ,PARAMETER  :: TOL=1.D-3
       REAL(8)                   :: EMIN,EMAX ! BOUNDS OF ENERGY GRID
       REAL(8)                   :: A,B ! PARAMETERS FOR LIFETIME CALC.
       REAL(8)                   :: EV 
       REAL(8)                   :: METER   
       REAL(8)                   :: NANOMETER   
       INTEGER(4),ALLOCATABLE    :: ISPECIES(:)
       LOGICAL(4)                :: SAVETRAWSTATES
       LOGICAL(4)                :: TKGROUP
       CHARACTER(132)            :: TITLE
       INTEGER(4)                :: GID   ! GRID ID FOR RADIAL GRID
       INTEGER(4)                :: NR    ! #(RADIAL GRID POINTS)
       REAL(8)                   :: RBOX  ! OUTERMOST ZERO FOR WAVE FUNCTIONS
       REAL(8)       ,ALLOCATABLE:: R(:)  !(NR) RADIAL GRID
       REAL(8)       ,ALLOCATABLE:: WORK(:)  !(NR) RADIAL WORK ARRAY
       COMPLEX(8)                :: DIEL(NE)
       REAL(8)                   :: EELSDOS(NE)
       REAL(8)                   :: LOSS(NE)
       REAL(8)                   :: BLOSS(NE)
       REAL(8)                   :: RBAS(3,3) !LATTICE VECTORS
       REAL(8)                   :: GBAS(3,3) !REC. SPACE LATTICE VECTORS
       REAL(8)                   :: CELLVOL   !UNIT-CELL VOLUME
       REAL(8)                   :: DENSITY
       REAL(8)                   :: LSCATT ! INELASTIC SCATTERING LENGTH
       INTEGER(4)                :: MC,LMC,LC
       INTEGER(4)                :: NKPT   !#(K-POINTS)
       INTEGER(4)                :: NSPIN  !#(COLLINEAR SPIN COMPONENTS)
       INTEGER(4)                :: NDIM   !#(SPINOR COMPONENTS)
       INTEGER(4)                :: NB     !#(BANDS)
       INTEGER(4)                :: NC     !#(CORE STATES)
       INTEGER(4)                :: LNX    !#(PARTIAL WAVES, ONLY RADIAL PARTS)
       INTEGER(4)                :: LMNX   !#(PARTIAL WAVES)
       INTEGER(4)                :: NBA    !#(NUMBER OF ATOMIC ENERGY LEVELS)
       INTEGER(4)                :: NFIL,IAT,ISP,IB,IKPT,ISPIN
       INTEGER(4)                :: LN,LMN,LM,L,M,IREALIMAG,I,IE,IE2
       INTEGER(4)                :: NTASKS,THISTASK
       INTEGER(4)                :: NTASKS_K,THISTASK_K
       REAL(8)                   :: AEZ  ! ATOMIC NUMBER
       REAL(8)                   :: CORELEVEL ! ENERGY OF THE CORE LEVEL
       REAL(8)                   :: WGHT
       REAL(8)                   :: CG(3)
       REAL(8)                   :: DELTAE! EXCITATION ENERGY FROM CORE LEVEL
       REAL(8)                   :: EFERMI
       REAL(8)                   :: E
       REAL(8)                   :: GAMMACORE
       REAL(8)                   :: GAMMA
       INTEGER(4)    ,ALLOCATABLE:: LOX(:)
       INTEGER(4)    ,ALLOCATABLE:: LOFI(:)
       REAL(8)       ,ALLOCATABLE:: EOFI(:)
       REAL(8)       ,ALLOCATABLE:: AEPSI(:,:)
       REAL(8)       ,ALLOCATABLE:: AEPHI(:,:)
       REAL(8)       ,ALLOCATABLE:: WKPT(:)
       REAL(8)       ,ALLOCATABLE:: OCC(:,:,:)
       REAL(8)       ,ALLOCATABLE:: EIGVAL(:)
       REAL(8)       ,ALLOCATABLE:: XVAL(:,:,:) !(LMNX,3,2*LX+1) MATRIX ELEMENTS
       REAL(8)       ,ALLOCATABLE:: PROJ(:) !(LMNX) PROJECTION
       REAL(8)                   :: RADINT ! RADIAL INTEGRAL
       REAL(8)                   :: STRENGTH ! OSCILLATOR STRENGTH
       REAL(8)                   :: SVAR,SVAR1
       REAL(8)                   :: SADD
       REAL(8)                   :: HBAROMEGA
       REAL(8)                   :: BANDLEVEL
       REAL(8)                   :: ENBMIN
       INTEGER(4)                :: NDEL,IDEL
       REAL(8)       ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
       REAL(8)       ,PARAMETER  :: SQPI43=SQRT(4.D0*PI/3.D0)
!      *************************************************************************
                          CALL TRACE$PUSH('OPTEELS_EELS')
       CALL MPE$QUERY('K',NTASKS_K,THISTASK_K) 
       CALL MPE$QUERY('MONOMER',NTASKS,THISTASK) 
       CALL CONSTANTS('EV',EV)
       CALL CONSTANTS('METER',METER)
       NANOMETER=1.D-9*METER
!
       CALL ATOMLIST$INDEX(ATOM,IAT)
       CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
       CALL SETUP$ISELECT(ISP)
       CALL SETUP$GETR8('AEZ',AEZ)
       CALL SETUP$GETI4('GID',GID)
       CALL SETUP$GETI4('NR',NR)
!
!      =========================================================================
!      == GET CORE STATES FOR THIS ATOM                                       ==
!      =========================================================================
       CALL SETUP$GETI4('NC',NC)
       IF (IC.LT.1.OR.IC.GT.NC) THEN             !CHECK CORE LEVEL INDEX  
         CALL ERROR$MSG('WRONG CORESHELL INDEX')
         CALL ERROR$I4VAL('IC',IC)
         CALL ERROR$I4VAL('NC',NC)
         CALL ERROR$STOP('OPTEELS$EELS')
       ENDIF
       CALL SETUP$GETI4('NR',NR)
       CALL SETUP$GETI4('NB',NBA) !NB FROM ATOM
       ALLOCATE(LOFI(NBA))
       ALLOCATE(EOFI(NBA))
       ALLOCATE(AEPSI(NR,NBA))
       CALL SETUP$GETI4A('LB',NBA,LOFI) !ANGULAR MOMENTA
       LC=LOFI(IC)
       CALL SETUP$GETR8A('EOFI',NBA,EOFI) !ENERGIES ON THIS NODE???
       CORELEVEL=EOFI(IC)
       DEALLOCATE(LOFI,EOFI)
       CALL SETUP$GETR8A('AEPSI',NR*NBA,AEPSI) 
!
!      =========================================================================
!      == GET CORE-HOLE LIFE TIME                                             ==
!      =========================================================================
       IF(LC.EQ.0) THEN
!        == EQ 12 OF EGERTON07_ULTRAMICROSCOPY107_565 ==========================
         E=0.D0-CORELEVEL
         GAMMACORE=(-0.285D0+0.0216D0*(E/EV)**0.472D0)*EV
       ELSE 
!        == FORMULA 14 OF EGERTON
         E=0.D0-CORELEVEL
         GAMMACORE=(7.4D-4*(E/EV) &
      &            -0.16D0*EXP(-(E/EV-690.D0)/400.D0)**2 &
      &            -0.05D0*EXP(-(E/EV- 30.D0)/ 90.D0)**2)*EV
       END IF
PRINT*,'GAMMACORE[EV] ',GAMMACORE/EV
!
!      =========================================================================
!      == GET PARTIAL WAVES                                                   ==
!      =========================================================================
       CALL SETUP$GETI4('LNX',LNX)
       ALLOCATE(LOX(LNX))
       CALL SETUP$GETI4A('LOX',LNX,LOX)
       ALLOCATE(AEPHI(NR,LNX)) 
       CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
       LMNX=SUM(2*LOX+1)
!
!      =========================================================================
!      == CALCULATE MATRIX ELEMENTS <PSI_C|X_J|\PHI_ALPHA>                    ==
!      =========================================================================
       CALL SETUP$GETR8('RBOX',RBOX)
       CALL SETUP$UNSELECT()
       ALLOCATE(R(NR))
       CALL RADIAL$R(GID,NR,R)
       ALLOCATE(XVAL(LMNX,3,2*LC+1))
       ALLOCATE(WORK(NR))
       XVAL(:,:,:)=0.D0
       LMN=0
       DO LN=1,LNX
         L=LOX(LN) 
         RADINT=0.D0
         IF(ABS(L-LC).EQ.1) THEN
           CALL RADIAL$INTEGRATE(GID,NR,R(:)**3*AEPSI(:,IC)*AEPHI(:,LN),WORK)
           CALL RADIAL$VALUE(GID,NR,WORK,RBOX,RADINT)
           DO M=1,2*L+1
             LMN=LMN+1
             LM=L**2+M
             DO MC=1,2*LC+1
               LMC=LC**2+MC
               CALL SPHERICAL$GAUNT(LM,2,LMC,CG(1)) !X/R
               CALL SPHERICAL$GAUNT(LM,4,LMC,CG(2)) !Y/R
               CALL SPHERICAL$GAUNT(LM,3,LMC,CG(3)) !Z/R
               XVAL(LMN,:,MC)=SQPI43*CG(:)*RADINT
             ENDDO
           ENDDO
         ELSE
           LMN=LMN+2*L+1
         END IF
       ENDDO
       DEALLOCATE(WORK)
!       
!      =========================================================================
!      == CALCULATE SITE DENSITY OF THIS ATOM                                 ==
!      =========================================================================
       CALL CELL$GETR8A('T(0)',9,RBAS)
       CALL GBASS(RBAS,GBAS,CELLVOL)
       DENSITY=1.D0/CELLVOL
!       
!      =========================================================================
!      == GET OCCUPATIONS, KP-WEIGHTS + EIGVAL OF ALL STATES                  ==
!      =========================================================================
       CALL WAVES$SETI4('IAT',IAT)
       CALL WAVES$GETL4('RAWSTATES',SAVETRAWSTATES)  !RESET ORIGINAL WAVE OBJECT
       CALL WAVES$SETL4('RAWSTATES',.FALSE.)
       CALL DYNOCC$GETI4('NB',NB) 
       CALL WAVES$GETI4('NKPT',NKPT)
       CALL WAVES$GETI4('NSPIN',NSPIN)
       CALL WAVES$GETI4('SPINORDIM',NDIM)
       ALLOCATE(WKPT(NKPT))
       ALLOCATE(OCC(NB,NKPT,NSPIN))
       ALLOCATE(EIGVAL(NB))
       CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT) 
       CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC) !0<OCC<2*WKPT(IKPT)/NSPIN
!       
!      =========================================================================
!      == ESTIMATE FERMI LEVEL FOR LIFETIME CALCULATION (ONLY A ROUGH ESTIMATE)=
!      == EMIN AND EMAX ARE THE ENERGY BOUNDS FOR THE SPECTRUM                ==
!      == EMIN IS CHOSEN 5EV BELOW THE ESTIMATED FERMI LEVEL.                 ==
!      == EMAX IS THE MINIMUM OF THE TOP MOST ENERGY BAND                     ==
!      ===    (BOTH SPINS CONSIDERED)                                         ==
!      =========================================================================
       EFERMI=0.D0
       EMIN=+1.D+12
       EMAX=-1.D+12
       ENBMIN=+1.D+12
       DO IKPT=1,NKPT
         CALL WAVES$SETI4('IKPT',IKPT) !REFERS TO GLOBAL IKPT 
         DO ISPIN=1,NSPIN
           CALL WAVES$SETI4('ISPIN',ISPIN)
           CALL WAVES$SETI4('IB',1)
           CALL WAVES$STATESELECTED(TKGROUP) 
           IF(.NOT.TKGROUP) CYCLE
PRINT*,'BEFORE EIGVAL ',THISTASK,THISTASK_K,ISPIN,IKPT,TKGROUP
           CALL WAVES$GETR8A('EIGVAL',NB,EIGVAL) !1KP,SPIN!!
PRINT*,'AFTER EIGVAL '
           ENBMIN=MIN(ENBMIN,EIGVAL(NB))
           DO IB=1,NB
             EMIN=MIN(EMIN,EIGVAL(IB))
             EMAX=MAX(EMAX,EIGVAL(IB))
             IF(OCC(IB,IKPT,ISPIN)/WKPT(IKPT).GT.0.5D0) THEN
               EFERMI=MAX(EFERMI,EIGVAL(IB))
             ELSE
               EFERMI=MIN(EFERMI,EIGVAL(IB))
             END IF
           ENDDO
         ENDDO
       ENDDO
       EMAX=EMAX+5.D0*EV
       EMAX=ENBMIN
       EMIN=EFERMI-5.D0*EV
       CALL MPE$COMBINE('MONOMER','MAX',EMAX)
       CALL MPE$COMBINE('MONOMER','MIN',EMIN)
PRINT*,'EFERMI[EV]     ',EFERMI/EV
PRINT*,'CORE LEVEL[EV] ',CORELEVEL/EV
PRINT*,'THRESHHOLD[EV] ',(EFERMI-CORELEVEL)/EV
PRINT*,'UPPER LIMIT OF ENERGY RANGE[EV] ',ENBMIN/EV
!       
!      =========================================================================
!      == DETERMINE DIELECTRIC CONSTANT                                       ==
!      =========================================================================
       ALLOCATE(PROJ(LMNX)) 
!FOR PARALLELIZATION: MAKE SURE THAT NKPT IUS THE LOCAL NKPTL
       EELSDOS(:)=0.D0
       DO IKPT=1,NKPT
         CALL WAVES$SETI4('IKPT',IKPT)
         DO ISPIN=1,NSPIN 
           CALL WAVES$SETI4('ISPIN',ISPIN)
           CALL WAVES$SETI4('IB',1)
           CALL WAVES$STATESELECTED(TKGROUP)
           IF(.NOT.TKGROUP) CYCLE
           CALL WAVES$GETR8A('EIGVAL',NB,EIGVAL) !1KP,SPIN!!
           DO IB=1,NB   
!            == PARALLELIZATION: K-POINTS WILL BE SELECTED WITH ================
!            == WAVES$STATESELECTED(TKGROUP). STATES WILL BE DISTRIBUTED ON ====
!            == THE NODES FOR ONE K-POINT USING MODULO. THE RESULT IS SUMMED ===
!            == OVER ALL STATES ================================================
             IF(MODULO(IB,NTASKS_K).NE.0) CYCLE
!
!            ===================================================================
!            == WGHT=(1-F)*WKPT ====== REMARK: NSPIN=1:NON-SPINPOL:OCC=2!=====
!            ===================================================================
!            == WEIGHT CONTAINES THE OCCUPATION, I.E. 1-F, THE K-POINT WEIGHT ==
!            == AND THE SPIN MULTIPLICITY ======================================
             IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
                WGHT=2.D0*WKPT(IKPT)-OCC(IB,IKPT,ISPIN)
             ELSE
                WGHT=WKPT(IKPT)-OCC(IB,IKPT,ISPIN) 
             END IF
             IF(ABS(WGHT/WKPT(IKPT)).LT.1.D-5) CYCLE
!
!            ===================================================================
!            == CALCULATE OSCILLATOR STRENGTHS =================================
!            ===================================================================
             CALL WAVES$SETI4('IB',IB)
     
             STRENGTH=0.D0
             DO IREALIMAG=1,2
               IF(IREALIMAG.EQ.1) THEN
                 CALL WAVES$SETL4('TIM',.TRUE.) !
               ELSE
                 CALL WAVES$SETL4('TIM',.FALSE.) !
               END IF
               CALL WAVES$GETR8A('<PSPSI|PRO>',LMNX,PROJ) 
               DO I=1,3
                 DO MC=1,LC+1
                   STRENGTH=STRENGTH+SUM(XVAL(:,I,MC)*PROJ(:))**2
                 ENDDO
               ENDDO
             ENDDO
!            == FACTOR 1/3 FROM AVERAGE OVER POLARIZATIONS OF THE LIGHT
!            == FACTOR 2 FROM DEFINITION OF OSCILLATOR STRENGTH (WOOTEN EQ.3.67)
             STRENGTH=2.D0*STRENGTH/3.D0*(EIGVAL(IB)-CORELEVEL)
!
!            ===================================================================
!            == PLACE MATRIX ELEMENTS ON ENERGY GRID                          ==
!            ===================================================================
             DELTAE=EIGVAL(IB)
             SVAR=1.D0+(DELTAE-EMIN)/(EMAX-EMIN)*REAL(NE-1)
             IE=INT(SVAR)
             SVAR=SVAR-REAL(IE)
             EELSDOS(IE  )=EELSDOS(IE)  +WGHT*STRENGTH*(1.D0-SVAR)
             EELSDOS(IE+1)=EELSDOS(IE+1)+WGHT*STRENGTH*SVAR
           ENDDO
         ENDDO
       ENDDO
       CALL MPE$COMBINE('MONOMER','+',EELSDOS)
!
!      =========================================================================
!      == APPLY LIFE TIME BROADENING                                          ==
!      =========================================================================
       DIEL(:)=(0.D0,0.D0)
       DO IE=1,NE
         IF(EELSDOS(IE).EQ.0.D0) CYCLE
!        __DELTAE=EXCITATION ENERGY
         BANDLEVEL=EMIN+(EMAX-EMIN)/REAL(NE-1)*REAL(IE-1)
         DELTAE=BANDLEVEL-CORELEVEL

!        =======================================================================
!        == DETERMINE LIFE TIME                                               ==
!        =======================================================================
!        == INELASTIC MEAN FREE PATH FROM SEAH79_SURFINTERFACEANAL1_2
!        == EQ. 5 AND TABLE 1 FOR LAMBDA_N, ENTRY FOR INORGANIC COMPOUNDS.
         A=641.D0*EV**2*NANOMETER
         B=0.096D0/SQRT(EV)*NANOMETER
         E=MAX(0.D0,BANDLEVEL-EFERMI)
         LSCATT=A/E**2+B*SQRT(E)
!        == LIFE TIME DUE TO EGERTON07_ULTRAMICROSCOPY107_565 EQ.15.
         GAMMA=SQRT(2.D0*E)/LSCATT
         GAMMA=GAMMA+GAMMACORE
!
!        ===================================================================
!        == ACCUMULATE DIELECTRIC CONSTANT                                ==
!        ===================================================================
         DO IE2=1,NE
           HBAROMEGA=EMIN+(EMAX-EMIN)/REAL(NE-1)*REAL(IE2-1)-CORELEVEL
           DIEL(IE2)=DIEL(IE2)+EELSDOS(IE) &
      &                       /(DELTAE**2-HBAROMEGA**2-CI*GAMMA*HBAROMEGA)
         ENDDO
       ENDDO
       DIEL(:)=(1.D0,0.D0)+4.D0*PI*DENSITY*DIEL(:)
!
!      =========================================================================
!      == CONVERT TO LOSS FUNCTION                                            ==
!      =========================================================================
       DO IE=1,NE
         LOSS(IE)=-AIMAG(1.D0/DIEL(IE))
       ENDDO
!
!      =========================================================================
!      == INSTRUMENTAL SMEARING OF THE LOSS FUNCTION                          ==
!      =========================================================================
!      == SMEAR WITH A GAUSSIAN 
!      == WITH FULL WIDTH HALF MAXIMUM = INSTRUMENTALFWHM
PRINT*,'INSTRUMENTALFWHM ',INSTRUMENTALFWHM,EMAX,EMIN
       SVAR1=2.D0*(EMAX-EMIN)/REAL(NE-1,KIND=8)/INSTRUMENTALFWHM
       SVAR1=LOG(2.D0)*SVAR1**2
       NDEL=NINT(SQRT(-LOG(TOL)/SVAR1))
       BLOSS(:)=LOSS(:)   !IDEL=0
       SADD=1.D0
       ALLOCATE(WORK(NE))
       WORK(:)=0.D0
       DO IDEL=1,NDEL
         SVAR=EXP(-SVAR1*REAL(IDEL,KIND=8)**2)
         BLOSS(1+IDEL:NE)=BLOSS(1+IDEL:NE)+SVAR*LOSS(1:NE-IDEL)
         BLOSS(1:NE-IDEL)=BLOSS(1:NE-IDEL)+SVAR*LOSS(1+IDEL:NE)
         SADD=SADD+2.D0*SVAR
         WORK(1+IDEL:NE)=WORK(1+IDEL:NE)+SVAR
         WORK(1:NE-IDEL)=WORK(1:NE-IDEL)+SVAR
       ENDDO
!       BLOSS(:)=BLOSS(:)/SADD
       BLOSS(:)=BLOSS(:)/WORK(:)
       DEALLOCATE(WORK)
!
!      =========================================================================
!      == CONVERT TO LOSS FUNCTION AND PRINT
!      =========================================================================
       IF(THISTASK.EQ.1.D0) THEN
PRINT*,'WRITING EELS TO FILE: ',TRIM(FILE)
         CALL FILEHANDLER$SETFILE('EELS',.FALSE.,FILE)
         CALL FILEHANDLER$SETSPECIFICATION('EELS','FORM','FORMATTED')
         CALL FILEHANDLER$UNIT('EELS',NFIL) !OPENS FILE NUMBER NFIL AT ITS END
         DO IE=1,NE
           E=EMIN+(EMAX-EMIN)/REAL(NE-1)*REAL(IE-1)
           WRITE(NFIL,*)E/EV,EELSDOS(IE),LOSS(IE),BLOSS(IE)
         ENDDO
         CALL FILEHANDLER$CLOSE('EELS')
PRINT*,'EELS WRITTEN'
       END IF
!
!      =========================================================================
!      == CLOSE DOWN                                                          ==
!      =========================================================================
       CALL WAVES$SETL4('RAWSTATES',SAVETRAWSTATES)
       DEALLOCATE(LOX,AEPSI,AEPHI, WKPT,OCC,EIGVAL)
                        CALL TRACE$POP
       RETURN
       END SUBROUTINE OPTEELS_EELS
