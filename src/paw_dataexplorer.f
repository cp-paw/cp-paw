!
!     ..................................................................
      SUBROUTINE NROUNDUP (X,NX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (X.GE.0.D0) NX=INT(X+0.99999999D0)
      IF (X.LT.0.D0) NX=INT(X)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE NROUNDDOWN  (X,NX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (X.LE.0.D0) NX=INT(X-0.99999999D0)
      IF (X.GT.0.D0) NX=INT(X)
      RETURN
      END
!
!.......................................................................
MODULE BALLSTICK_MODULE
REAL(8)   ,PARAMETER :: SCALEBOND=1.25D0
REAL(8)   ,PARAMETER :: SCALESPHERE=0.5D0
LOGICAL(4)           :: TINIT=.FALSE.
LOGICAL(4)           :: TLATTICE=.FALSE.
INTEGER(4)           :: NREP
INTEGER(4),ALLOCATABLE :: IFLAG(:)       !(NATOMX)
INTEGER(4)            :: NLAT(3,2)
INTEGER(4)           :: NATOMX
END MODULE BALLSTICK_MODULE
!     
!     ..................................................................
      SUBROUTINE BALLSTICK$INITIALIZE(NFIL)
!     ==================================================================
!     ==  INITIALIZATION                                              ==
!     ==================================================================
      USE BALLSTICK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)               :: RBAS(3,3)
      REAL(8)               :: R(3)
      REAL(8)               :: TMAT(3,3)
      REAL(8)               :: ORIGIN(3)
      INTEGER(4)            :: I,J,K,IAT,IX,IY,IZ,IATOM,IREP
      LOGICAL(4)            :: TCHECK
      REAL(8)               :: XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3
      REAL(8)               :: DUMMY1,DUMMY2
      INTEGER(4),PARAMETER  :: NAUX=300
      REAL(8)               :: AUX(NAUX)
      REAL(8)   ,ALLOCATABLE:: POS(:,:)  !(3,NATOMX)
      INTEGER(4)            :: NATOM
      INTEGER(4)            :: NAT
      REAL(8)               :: TMP(3)
      REAL(8)               :: V1(3),V2(3)
      REAL(8)   ,ALLOCATABLE:: RAD(:)         !(NATOMX)
      INTEGER(4),ALLOCATABLE:: ICOLOR(:,:)    !(3,NATOMX)
      REAL(8)               :: AEZ
      INTEGER(4)            :: IZATOM
      REAL(8)               :: RAD1
!     ******************************************************************
                            CALL TRACE$PUSH('BALLSTICK$INITIALIZE')
      IF (TINIT) THEN 
        CALL ERROR$MSG('BALLSTICK HAS ALREADY BEEN INITIALIZED')
        CALL ERROR$STOP('BALLSTICK$INITIALIZE')
      END IF
!     
!     ================================================================
!     ==  SET DEFAULT VALUES                                        ==
!     ================================================================
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      DO I =1,3
        DO J = 1,3
          TMAT(I,J)=RBAS(I,J)
        ENDDO
        ORIGIN(I)=0.0D0
      ENDDO
!      
!     ================================================================
!     ==  ESTIMATE NUMBER OF ATOMS AND ALLOCATE ARRAYS              ==
!     ================================================================
      CALL WRITEDXBOX$BOX(ORIGIN,TMAT,TCHECK)
      CALL BOXBOX(RBAS,ORIGIN,TMAT,XMIN1,XMAX1,XMIN2,XMAX2 &
     &           ,XMIN3,XMAX3)
      CALL NROUNDDOWN(XMIN1,NLAT(1,1))
      CALL NROUNDDOWN(XMIN2,NLAT(2,1))
      CALL NROUNDDOWN(XMIN3,NLAT(3,1))
      CALL NROUNDUP(XMAX1,NLAT(1,2))
      CALL NROUNDUP(XMAX2,NLAT(2,2))
      CALL NROUNDUP(XMAX3,NLAT(3,2))
!     
      call lib$invertr8(3,tmat,tmat)
!CALL DGEICD (TMAT,3,3,0,DUMMY1,DUMMY2,AUX,NAUX)
!     
      NREP=1
      DO I =1,3
        NREP=NREP*(NLAT(I,2)-NLAT(I,1)+1)
      ENDDO
!     
!     ================================================================
!     ==  CALCULATE SPHERE POSITIONS                                ==
!     ================================================================
      NATOMX=NAT*NREP
      ALLOCATE(POS(3,NATOMX))
      ALLOCATE(IFLAG(NATOMX))
!     
      NATOM = 0
      DO IAT=1,NAT
        CALL ATOMLIST$GET('R(0)',8*3,IAT,R)
        DO IX = NLAT (1,1),NLAT(1,2)
          DO IY = NLAT (2,1),NLAT(2,2)
            DO IZ = NLAT (3,1),NLAT(3,2)
              TMP(1)=DBLE(IX)
              TMP(2)=DBLE(IY)
              TMP(3)=DBLE(IZ)
              NATOM = NATOM+1
              DO J = 1,3
                POS(J,NATOM)=R(J)
                DO K = 1,3
                  POS(J,NATOM)=POS(J,NATOM)+RBAS(J,K)*TMP(K)
                ENDDO
              ENDDO
              IFLAG(NATOM)=1
              DO I = 1,3
                V2(I)=0.D0
                V1(I)=POS(I,NATOM)-ORIGIN(I)
              ENDDO
              CALL DGEMV('N',3,3,1.D0,TMAT,3,V1,1,1.D0,V2,1)
              DO J = 1,3
                IF (V2(J).LT.0.D0.OR.V2(J).GT.1.D0)IFLAG(NATOM)=0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     
!     ================================================================
!     == DETERMINE ATOM-COLORS AND RADII =============================
!     ================================================================
      ALLOCATE(RAD(NATOMX))
      ALLOCATE(ICOLOR(3,NATOMX))
      IATOM=0
      NATOM=0
      DO IAT=1,NAT
        CALL ATOMLIST$GET('Z',8,IAT,AEZ)
        IZATOM=NINT(AEZ)
        CALL PERIODICTABLER8(IZATOM,'R(COV)',RAD1)
        RAD1=RAD1*SCALESPHERE
        DO IREP=1,NREP
          IATOM=IATOM+1
          IF(IFLAG(IATOM).GT.0)THEN
            NATOM=NATOM+1
            CALL ATOMCOLOR(IZATOM,ICOLOR(1,NATOM))
            RAD(NATOM)=RAD1
          ENDIF
        ENDDO
      ENDDO  
!     
!     ================================================================
!     ==  WRITE HEADER FOR BALL STICK FILE                          ==
!     ================================================================
      REWIND NFIL
      CALL WRITEDX_BALLSTICK$START(NFIL,NATOM,RAD,ICOLOR,TLATTICE)
      DEALLOCATE(ICOLOR)
      DEALLOCATE(RAD)
      DEALLOCATE(POS)
      TINIT = .TRUE.
                            CALL TRACE$POP
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE BALLSTICK$ADDFRAME(NFIL)
!     ******************************************************************
!     **  CALCULATE SPHERE POSITIONS                                  **
!     ******************************************************************
      USE BALLSTICK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)               :: R(3)
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IATOM
      INTEGER(4)            :: NATOM
      INTEGER(4)            :: IAT,IREP,INAT,JNAT,IX,IY,IZ,J,K,I1,I2
      REAL(8)               :: AEZ
      REAL(8)               :: RAD1
      INTEGER(4)            :: IZATOM
      REAL(8)   ,ALLOCATABLE:: RAD(:)   !(NATOMX)
      REAL(8)   ,ALLOCATABLE:: POS(:,:)   !(3,NATOMX)
      REAL(8)               :: TMP(3)
      INTEGER(4)            :: NBONDX
      INTEGER(4)            :: NBOND
      INTEGER(4),ALLOCATABLE:: IBOND(:,:)  !(2,IBOND)
      REAL(8)               :: DIS
!     ******************************************************************
                            CALL TRACE$PUSH('BALLSTICK$ADDFRAME')
!     
!     ================================================================
!     == DETERMINE ATOM-COLORS AND RADII =============================
!     ================================================================
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(RAD(NATOMX))
      IATOM=0
      NATOM=0
      DO IAT=1,NAT
        CALL ATOMLIST$GET('Z',8,IAT,AEZ)
        IZATOM=NINT(AEZ)
        CALL PERIODICTABLER8(IZATOM,'R(COV)',RAD1)
        RAD1=RAD1*SCALEBOND
        DO IREP=1,NREP
          IATOM=IATOM+1
          IF(IFLAG(IATOM).GT.0)THEN
            NATOM=NATOM+1
            RAD(NATOM)=RAD1
          ENDIF
        ENDDO
      ENDDO  
!     ================================================================
!     ==  REWRITE THE POS ARRAY                                     ==
!     ================================================================
                            CALL TRACE$PASS('A')
      ALLOCATE(POS(3,NATOM))
!     
      INAT = 0
      JNAT = 0
      DO IAT=1,NAT
        CALL ATOMLIST$GET('R(0)',8*3,IAT,R)
        DO IX = NLAT(1,1),NLAT(1,2)
          DO IY = NLAT(2,1),NLAT(2,2)
            DO IZ = NLAT(3,1),NLAT(3,2)
              TMP(1)=DBLE(IX)
              TMP(2)=DBLE(IY)
              TMP(3)=DBLE(IZ)
              INAT = INAT+1
              IF(IFLAG(INAT).GT.0)THEN
                JNAT=JNAT+1
                DO J = 1,3
                  POS(J,JNAT)=R(J)
                  DO K = 1,3
                    POS(J,JNAT)=POS(J,JNAT)+RBAS(J,K)*TMP(K)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO  
                            CALL TRACE$PASS('B')
!     ================================================================
!     ==  CALCULATE BONDS                                           ==
!     ================================================================
      NBONDX = 4*NATOM      ! FIRST GUESS (OVERFLOW PROTECTION LATER)
 10   CONTINUE
      ALLOCATE(IBOND(2,NBONDX))
      NBOND=0
      DO I1=1,NATOM-1 
        DO I2=I1+1,NATOM
          DIS = DSQRT( ( POS(1,I1)-POS(1,I2))**2 &
     &                +( POS(2,I1)-POS(2,I2))**2 &
     &                +( POS(3,I1)-POS(3,I2))**2)
!     
          IF(DIS.LT.(RAD(I1)+RAD(I2))) THEN
            NBOND=NBOND+1
            IF(NBOND.LE.NBONDX) THEN
              IBOND(1,NBOND)=I1
              IBOND(2,NBOND)=I2          
           ENDIF
         ENDIF
        ENDDO
      ENDDO
!     
      IF (NBOND.GT.NBONDX )THEN
        NBONDX = NBOND
        DEALLOCATE(IBOND)
        GOTO 10
      ENDIF
!     
!     ================================================================
!     ==  WRITE BALLSTICK OBJECT                                    ==
!     ================================================================
                            CALL TRACE$PASS('C')
      CALL WRITEDX_BALLSTICK$ADDFRAME(NFIL,NATOM,POS,RBAS,NBOND,IBOND)
                            CALL TRACE$PASS('D')
!     
      DEALLOCATE(RAD)
      DEALLOCATE(POS)
      DEALLOCATE(IBOND)
                          CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BALLSTICK$END(NFIL)
!     ==================================================================
!     ==  DEFINE BOX FOR WHICH ATOMS ARE SHOWN                        ==
!     ==================================================================
      USE BALLSTICK_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
!     *******************************************************************
                            CALL TRACE$PUSH('BALLSTICK$END')
      CALL WRITEDX_BALLSTICK$END(NFIL)
      DEALLOCATE(IFLAG)
      TINIT=.FALSE.
                            CALL TRACE$POP
      RETURN
      END
!
!     .....................................................BONDS .......
      SUBROUTINE WRITEDX_BALLSTICK
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                          
      REAL(8) RBAS(3,3)
      INTEGER ICOLOR(3,NATOM)
      DIMENSION RAD(NATOM)
      DIMENSION POS(3,NATOM)
      DIMENSION IBOND(2,NBOND)
!
      LOGICAL TSET,TLATTICE,TLATTICE_
      SAVE NFRAME,TLATTICE,NATOM1
      SAVE TSET
      DATA TSET/.FALSE./
!     
!     ==================================================================
!     ==  WRITE HEADER                                                ==
!     ==================================================================
      ENTRY WRITEDX_BALLSTICK$START(NFIL,NATOM,RAD,ICOLOR,TLATTICE_)
!
        TLATTICE=TLATTICE_
        NATOM1=NATOM
        IF(TSET) THEN
          CALL ERROR$MSG('BALLSTICK HEADER ALREADY WRITTEN')
          CALL ERROR$STOP('WRITEDX_BALLSTICK$HEADER')
        END IF
        TSET=.TRUE.
        NFRAME=0
        ISPHERERADIUSOBJECT=1
        ISPHERECOLOROBJECT=2
!
!       ================================================================
!       ==   DATA ARRAY: SIZE OF THE SPHERES                          ==
!       ================================================================
        WRITE(NFIL,FMT='("#"/"#",T10,"SPHERE SIZE"/"#")')
        WRITE(NFIL,FMT='("OBJECT ",I10, &
     &                  /5X,"CLASS ARRAY" &
     &                  ,5X,"TYPE FLOAT" &
     &                  ,5X,"RANK 0" &
     &                  ,5X,"ITEMS ",I10 &
     &                  /"DATA FOLLOWS")')ISPHERERADIUSOBJECT,NATOM
        WRITE(NFIL,FMT='(10F10.5)')(RAD(I),I=1,NATOM)
        WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
        WRITE(NFIL,FMT='("#")')
!
!       ================================================================
!       ==   COLORS                                                   ==
!       ================================================================
        WRITE(NFIL,FMT='("#"/"#",T10,"COLORS"/"#")')
        WRITE(NFIL,FMT='("OBJECT ",I10 &
     &              /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &              /"DATA FOLLOWS")')ISPHERECOLOROBJECT,NATOM
!
        DO IATOM=1,NATOM
          X=DBLE(ICOLOR(1,IATOM))
          Y=DBLE(ICOLOR(2,IATOM))
          Z=DBLE(ICOLOR(3,IATOM))
!         SVAR=DSQRT(X**2+Y**2+Z**2)
          SVAR=200.D0
          X=X/SVAR
          Y=Y/SVAR
          Z=Z/SVAR
          WRITE(NFIL,FMT='(3F10.4)')X,Y,Z
        ENDDO
!     
        WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
        WRITE(NFIL,FMT='("#")')
        RETURN
!     
!     ==================================================================
!     ==  WRITE FRAME                                                 ==
!     ==================================================================
      ENTRY WRITEDX_BALLSTICK$ADDFRAME(NFIL,NATOM,POS,RBAS,NBOND,IBOND)
!
        IF(NATOM.NE.NATOM1) THEN
          CALL ERROR$MSG('NUMBER OF BALLS MUST NOT CHANGE DURING MOVIE')
          CALL ERROR$STOP('WRITEDX_BALLSTICK$ADDFRAME')
        END IF
        NFRAME=NFRAME+1
        IPOSITIONSOBJECT=(NFRAME-1)*4+2+1
        IBONDSOBJECT    =(NFRAME-1)*4+2+2
        ICELLOBJECT     =(NFRAME-1)*4+2+3
        IMOLECULEOBJECT =(NFRAME-1)*4+2+4
!
!       ================================================================
!       ==   POSITIONS ARRAY: ATOMIC POSITIONS                        ==
!       ================================================================
        WRITE(NFIL,FMT='("#"/"#",T10,"ATOMIC POSITIONS"/"#")')
        WRITE(NFIL,FMT='("OBJECT ",I10 &
     &              /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &              /"DATA FOLLOWS")')IPOSITIONSOBJECT,NATOM
        DO IA=1,NATOM
            WRITE(NFIL,FMT='(3F10.5)')(POS(I,IA),I=1,3)
        ENDDO  
        WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
        WRITE(NFIL,FMT='("#")')
!
!       ================================================================
!       ==   CONNECTIONS ARRAY: BONDS                                 ==
!       ================================================================
        WRITE(NFIL,FMT='("#"/"#",T10,"BONDS"/"#")')
        WRITE(NFIL,FMT='("OBJECT ",I10 &
     &              /"CLASS ARRAY TYPE INT RANK 1 SHAPE 2 ITEMS ",I10 &
     &              /"DATA FOLLOWS")')IBONDSOBJECT,NBOND
        WRITE(NFIL,FMT='(2I5)')((IBOND(I,IB)-1,I=1,2),IB=1,NBOND)
        WRITE(NFIL,FMT='("ATTRIBUTE ""REF"" STRING ""POSITIONS""")')
        WRITE(NFIL &
     &       ,FMT='("ATTRIBUTE ""ELEMENT TYPE"" STRING ""LINES""")')
        WRITE(NFIL,FMT='("#")')
!
!       ================================================================
!       ==   BOX: LATTICE VECTORS                                     ==
!       ================================================================
        WRITE(NFIL,FMT='("#"/"#",T10,"LATTICE VECTORS"/"#")')
        WRITE(NFIL,FMT='("OBJECT ",I10 &
     &              /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &              /"DATA FOLLOWS")')ICELLOBJECT,8
        DO I3=0,1 
           T3=DBLE(I3)-0.5D0
          DO I2=0,1
            T2=DBLE(I2)-0.5D0
            DO I1=0,1
              T1=DBLE(I1)-0.5D0
              IF (TLATTICE)THEN
                 T1=T1+0.5D0
                 T2=T2+0.5D0
                 T3=T3+0.5D0
              ENDIF
              X=RBAS(1,1)*T1+RBAS(1,2)*T2+RBAS(1,3)*T3
              Y=RBAS(2,1)*T1+RBAS(2,2)*T2+RBAS(2,3)*T3
              Z=RBAS(3,1)*T1+RBAS(3,2)*T2+RBAS(3,3)*T3
              WRITE(NFIL,FMT='(3F10.5)')X,Y,Z
            ENDDO
          ENDDO
        ENDDO
        WRITE(NFIL,FMT='("#")')
!
!       ================================================================
!       ==   OBJECT MOLECULE:                                         ==
!       ================================================================
        WRITE(NFIL,FMT='("OBJECT ",I10)')IMOLECULEOBJECT
        WRITE(NFIL,FMT='("CLASS FIELD")')
        WRITE(NFIL,FMT='("COMPONENT ""DATA"" VALUE ",I10)') &
     &        ISPHERERADIUSOBJECT
        WRITE(NFIL,FMT='("COMPONENT ""POSITIONS"" VALUE ",I10)') &
     &        IPOSITIONSOBJECT
        WRITE(NFIL,FMT='("COMPONENT ""CONNECTIONS"" VALUE ",I10)') &
     &        IBONDSOBJECT
        WRITE(NFIL,FMT='("COMPONENT ""BOX"" VALUE ",I10)') &
     &        ICELLOBJECT
        WRITE(NFIL,FMT='("COMPONENT ""COLORS"" VALUE ",I10)') &
     &       ISPHERECOLOROBJECT
        WRITE(NFIL,FMT='("ATTRIBUTE ""NAME"" STRING ""ATOMS""")')
        WRITE(NFIL,FMT='("#")')
        RETURN
!     
!     ==================================================================
!     ==  WRITE END                                                   ==
!     ==================================================================
      ENTRY WRITEDX_BALLSTICK$END(NFIL)
        IF (NFRAME.NE.1) THEN 
          WRITE(NFIL,FMT='("OBJECT ""SERIES""")')
          WRITE(NFIL,FMT='("CLASS SERIES")')
          DO IFRAME=1,NFRAME
            IMOLECULEOBJECT =(IFRAME-1)*4+2+4
            WRITE(NFIL &
     &           ,FMT='("MEMBER ",I10," VALUE ",I10," POSITION ",I10)') &
     &           IFRAME-1,IMOLECULEOBJECT,IFRAME
          ENDDO
          WRITE(NFIL,FMT='("#")')
        END IF
        WRITE(NFIL,FMT='("END")')
        TSET=.FALSE.
        RETURN
!
      END
!        
!     ..................................................................
      SUBROUTINE ATOMCOLOR(IZ,ICOLOR)
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      INTEGER ICOLOR(3),ICOLORSTANDARD(3,106)
      DATA (ICOLORSTANDARD(I,  1),I=1,3)/135,125,131/
      DATA (ICOLORSTANDARD(I,  2),I=1,3)/255,228,196/!BISQUE
      DATA (ICOLORSTANDARD(I,  3),I=1,3)/240,248,255/!ALICE BLUE
      DATA (ICOLORSTANDARD(I,  4),I=1,3)/200,200,  0/!YELLOW
      DATA (ICOLORSTANDARD(I,  5),I=1,3)/192,102,624/!CORAL
      DATA (ICOLORSTANDARD(I,  6),I=1,3)/ 50, 50, 50/
      DATA (ICOLORSTANDARD(I,  7),I=1,3)/ 10,200, 10/
      DATA (ICOLORSTANDARD(I,  8),I=1,3)/255,  0,  0/!RED
      DATA (ICOLORSTANDARD(I,  9),I=1,3)/ 50,205, 50/!LIME GREEN
      DATA (ICOLORSTANDARD(I, 10),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 11),I=1,3)/ 60,  1,  1/!REDISH BLACK
      DATA (ICOLORSTANDARD(I, 12),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 13),I=1,3)/180,  1,100/!PURPLE
      DATA (ICOLORSTANDARD(I, 14),I=1,3)/  1, 20,198/!DARK BLUE
      DATA (ICOLORSTANDARD(I, 15),I=1,3)/230,171, 17/
      DATA (ICOLORSTANDARD(I, 16),I=1,3)/240,240,  0/
      DATA (ICOLORSTANDARD(I, 17),I=1,3)/ 60,180,  0/!YELLOW GREEN
      DATA (ICOLORSTANDARD(I, 18),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 19),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 20),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 21),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 22),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 23),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 24),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 25),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 26),I=1,3)/176, 48, 96/
      DATA (ICOLORSTANDARD(I, 27),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 28),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 29),I=1,3)/278,134, 34/!FIREBRICK
      DATA (ICOLORSTANDARD(I, 30),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 31),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 32),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 33),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 34),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 35),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 36),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 37),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 38),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 39),I=1,3)/ 40, 40,140/
      DATA (ICOLORSTANDARD(I, 40),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 41),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 42),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 43),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 44),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 45),I=1,3)/230, 51, 41/
      DATA (ICOLORSTANDARD(I, 46),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 47),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 48),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 49),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 50),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 51),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 52),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 53),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 54),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 55),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 56),I=1,3)/100, 40,  0/
      DATA (ICOLORSTANDARD(I, 57),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 58),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 59),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 60),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 61),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 62),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 63),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 64),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 65),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 66),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 67),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 68),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 69),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 70),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 71),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 72),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 73),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 74),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 75),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 76),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 77),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 78),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 79),I=1,3)/255,215,  0/!GOLD
      DATA (ICOLORSTANDARD(I, 80),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 81),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 82),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 83),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 84),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 85),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 86),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 87),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 88),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 89),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 90),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 91),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 92),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 93),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 94),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 95),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 96),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 97),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 98),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 99),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,100),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,101),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,102),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,103),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,104),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,105),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,106),I=1,3)/  1,  1,  1/   
      DO I=1,3
        ICOLOR(I)=ICOLORSTANDARD(I,IZ)
      ENDDO
      RETURN
      END
!                                                                       
!     .....................................................BONDS .......
      SUBROUTINE WRITEDX_DENSITY(NFIL,NR1X,NR2X,NR3X &
     &                          ,NR1,NR2,NR3,RBAS,DENSITY)
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TCHK
      REAL(8) RBAS(3,3)
      REAL(8) DENSITY(NR1X,NR2X,NR3X)
      DIMENSION BOXORIGIN(3),BOXVECTORS(3,3)
      CALL DX$SET('RBAS',8*3*3,RBAS)
      IF(NR1X.NE.NR1.OR.NR2X.NE.NR2.OR.NR3X.NE.NR3) THEN
        CALL ERROR$MSG('NRIX.NE.NRI IS NOT ALLOWED')
        CALL ERROR$STOP('WRITEDX_DENSITY')
      END IF
      CALL DX$DENSITY(NFIL,NR1,NR2,NR3,DENSITY)
      RETURN
!
      DO I=1,3
        BOXORIGIN(I)=0.D0+1.D-8
        DO J=1,3
          BOXVECTORS(I,J)=RBAS(I,J)*(1.D0-1.D-8)
        ENDDO
      ENDDO
!
      CALL WRITEDXBOX$BOX(BOXORIGIN,BOXVECTORS,TCHK)
      CALL BOXBOX(RBAS,BOXORIGIN,BOXVECTORS &
     &           ,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
!
      XMIN1=XMIN1*NR1-1.D0
      XMAX1=XMAX1*NR1+1.D0
      XMIN2=XMIN2*NR2-1.D0
      XMAX2=XMAX2*NR2+1.D0
      XMIN3=XMIN3*NR3-1.D0
      XMAX3=XMAX3*NR3+1.D0
      MIN1=INT(ABS(XMIN1))*INT(SIGN(1.D0,XMIN1))
      MAX1=INT(ABS(XMAX1))*INT(SIGN(1.D0,XMAX1))
      MIN2=INT(ABS(XMIN2))*INT(SIGN(1.D0,XMIN2))
      MAX2=INT(ABS(XMAX2))*INT(SIGN(1.D0,XMAX2))
      MIN3=INT(ABS(XMIN3))*INT(SIGN(1.D0,XMIN3))
      MAX3=INT(ABS(XMAX3))*INT(SIGN(1.D0,XMAX3))
      N1=MAX1-MIN1+1
      N2=MAX2-MIN2+1
      N3=MAX3-MIN3+1
!
!     ==================================================================
!     ==================================================================
!     ==   WRITE DX FILE                                              ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     ==   DATA ARRAY: DENSITY                                        ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 1" &
     &                /5X,"CLASS ARRAY" &
     &                ,5X,"TYPE FLOAT" &
     &                ,5X,"RANK 0" &
     &                ,5X,"ITEMS ",I10 &
     &                /"DATA FOLLOWS")')N1*N2*N3
      WRITE(NFIL,FMT='(10F12.5)') &
     &     (((DENSITY(MOD(I+1000*NR1,NR1)+1 &
     &               ,MOD(J+1000*NR2,NR2)+1 &
     &               ,MOD(K+1000*NR3,NR3)+1) &
     &               ,K=MIN3,MAX3) &
     &               ,J=MIN2,MAX2) &
     &               ,I=MIN1,MAX1)
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   POSITIONS ARRAY: ATOMIC POSITIONS                          ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 2" &
     &              /"CLASS GRIDPOSITIONS COUNTS",3I10)') &
     &              N1,N2,N3
      SH1=RBAS(1,1)/DBLE(NR1)*MIN1 &
     &   +RBAS(1,2)/DBLE(NR2)*MIN2 &
     &   +RBAS(1,3)/DBLE(NR3)*MIN3
      SH2=RBAS(2,1)/DBLE(NR1)*MIN1 &
     &   +RBAS(2,2)/DBLE(NR2)*MIN2 &
     &   +RBAS(2,3)/DBLE(NR3)*MIN3
      SH3=RBAS(3,1)/DBLE(NR1)*MIN1 &
     &   +RBAS(3,2)/DBLE(NR2)*MIN2 &
     &   +RBAS(3,3)/DBLE(NR3)*MIN3
      WRITE(NFIL,FMT='("ORIGIN",3F10.5)')SH1,SH2,SH3
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(RBAS(I,1)/DBLE(NR1),I=1,3)
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(RBAS(I,2)/DBLE(NR2),I=1,3)
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(RBAS(I,3)/DBLE(NR3),I=1,3)
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   CONNECTIONS ARRAY: CONNECTIONS                             ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 3" &
     &              /"CLASS GRIDCONNECTIONS COUNTS",3I10)') &
     &              N1,N2,N3
      WRITE(NFIL,FMT='("ATTRIBUTE ""ELEMENT TYPE"" STRING ""CUBES""")')
      WRITE(NFIL,FMT='("ATTRIBUTE ""REF"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ================================================================
!     ==   CLIPBOX                                                  ==
!     ================================================================
      WRITE(NFIL,FMT='("OBJECT ",I10 &
     &            /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &            /"DATA FOLLOWS")')4,8
      DO I3=0,1 
         T3=DBLE(I3)
        DO I2=0,1
          T2=DBLE(I2)
          DO I1=0,1
            T1=DBLE(I1)
            X=BOXORIGIN(1)+BOXVECTORS(1,1)*T1 &
     &                    +BOXVECTORS(1,2)*T2 &
     &                    +BOXVECTORS(1,3)*T3
            Y=BOXORIGIN(2)+BOXVECTORS(2,1)*T1 &
     &                    +BOXVECTORS(2,2)*T2 &
     &                    +BOXVECTORS(2,3)*T3
            Z=BOXORIGIN(3)+BOXVECTORS(3,1)*T1 &
     &                    +BOXVECTORS(3,2)*T2 &
     &                    +BOXVECTORS(3,3)*T3
            WRITE(NFIL,FMT='(3F10.5)')X,Y,Z
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   OBJECT MOLECULE:                                           ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT ""DENSITY""")')
      WRITE(NFIL,FMT='("CLASS FIELD")')
      WRITE(NFIL,FMT='("COMPONENT ""DATA"" VALUE 1")')
      WRITE(NFIL,FMT='("COMPONENT ""POSITIONS"" VALUE 2")')
      WRITE(NFIL,FMT='("COMPONENT ""CONNECTIONS"" VALUE 3")')
      WRITE(NFIL,FMT='("COMPONENT ""BOX"" VALUE 4")')
      WRITE(NFIL,FMT='("ATTRIBUTE ""NAME"" STRING ""DENSITY""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   END                                                        ==
!     ==================================================================
      WRITE(NFIL,FMT='("END")')
      RETURN
      END
!                                                                       
!     .....................................................WRITE .......
      SUBROUTINE WRITEDX_PATH(NFIL,NAT,NSTEP,POSI,FORC)
!     **                                                              **
!     **  WRITE INPUT FOR THE DATAEXPLORER FOR PLOTTING TRAJECTORIES  **
!     **  AND FORCES                                                  **
!     **                                                              **
!     **  WARNING! POSI AND FORC ARE SINGLE PRECISION ARRAYS!!!       **
!     **                                             P. MARGL         **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      REAL(4)  POSI(3,NAT,NSTEP),FORC(3,NAT,NSTEP)
      REWIND NFIL
!
!     ================================================================
!     ==  WRITE LINES                                               ==
!     ================================================================
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'#  CONNECTIONS (COMMON FOR ALL ATOMS)'
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'OBJECT',0
      WRITE(NFIL,*)'CLASS GRIDCONNECTIONS COUNTS 1'
      WRITE(NFIL,*)'ATTRIBUTE "REF" STRING "POSITIONS"'
      WRITE(NFIL,*)'ATTRIBUTE "ELEMENT TYPE" STRING "LINES"'
      WRITE(NFIL,*)'#'
!
!     ================================================================
!     ==   LOPOP OVER ATOMS                                         ==
!     ================================================================
      DO IAT=1,NAT     
        IOBJECT1=(IAT-1)*3+1
        IOBJECT2=(IAT-1)*3+2
        IOBJECT3=(IAT-1)*3+3
!
!       ================================================================
!       ==   WRITE FORCES (DATA)                                      ==
!       ================================================================
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'#  FORCES OF ATOM ',IAT
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'OBJECT ',IOBJECT1
        WRITE(NFIL,*)'CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3', &
     & ' ITEMS',NSTEP
        WRITE(NFIL,*)'DATA FOLLOWS'
        DO ISTEP=1,NSTEP
          WRITE(NFIL,FMT='(1X,3F15.9)')(FORC(I,IAT,ISTEP),I=1,3)
        ENDDO
        WRITE(NFIL,*)'ATTRIBUTE "DEP" STRING "POSITIONS"'
        WRITE(NFIL,*)'#'
!
!       ================================================================
!       ==  WRITE THE PATHPOSITIONS                                   ==
!       ================================================================
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'#  POSITIONS OF ATOM ',IAT
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'OBJECT',IOBJECT2
        WRITE(NFIL,*)'CLASS ARRAY   TYPE FLOAT  RANK 1 SHAPE 3', &
     & ' ITEMS',NSTEP
        WRITE(NFIL,*)'DATA FOLLOWS'
        DO ISTEP=1,NSTEP
          WRITE(NFIL,FMT='(1X,3F15.9)')(POSI(I,IAT,ISTEP),I=1,3)
        ENDDO
        WRITE(NFIL,*)'ATTRIBUTE "DEP" STRING "POSITIONS"'
        WRITE(NFIL,*)'#'
!
!       ================================================================
!       ==                                                            ==
!       ================================================================
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'#  FIELD OBJECT OF ATOM ',IAT
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'OBJECT',IOBJECT3
        WRITE(NFIL,*)'CLASS FIELD'
        WRITE(NFIL,*)'COMPONENT "DATA" VALUE ',IOBJECT1
        WRITE(NFIL,*)'COMPONENT "POSITIONS" VALUE ',IOBJECT2
        WRITE(NFIL,*)'COMPONENT "CONNECTIONS" VALUE ',0
        WRITE(NFIL,*)'ATTRIBUTE "NAME" STRING "PATH"  '
        WRITE(NFIL,*)'#'
      ENDDO
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'#  COMPOSITE FIELD FOR PATH'
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'OBJECT ',3*NAT+1,' COMPOSITEFIELD'
      DO IAT=1,NAT
        WRITE(NFIL,*)'MEMBER ',IAT-1,' VALUE ',3*IAT
      ENDDO
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'END'
!
      RETURN
      END
!  
!     ..................................................................
      SUBROUTINE WRITEPATHFORDX
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **  USES:                                                       **
!     **   FILEHANDLER,STACK,WRITEDX                                  **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      REAL(8),ALLOCATABLE:: POS(:)
      REAL(4),ALLOCATABLE:: POSTRA(:)
      REAL(4),ALLOCATABLE:: FORTRA(:)
      DATA NSKIP/2/
!
!     == READ POSTIONS
      CALL FILEHANDLER$UNIT('POSITION-TRAJECTORY',NFIL)    
      REWIND NFIL
!     BACKSPACE (NFIL)
      NSTEP=0
      DO WHILE (NSTEP.LT.1000000)
        NSTEP=NSTEP+1
        READ(NFIL,END=1000)NFI,NAT
      ENDDO
 1000 CONTINUE
      NSTEP=NSTEP/(NSKIP+1)
      IF(NSTEP.LT.2) THEN
        CALL FILEHANDLER$CLOSE('POSITION-TRAJECTORY')    
        RETURN
      END IF
      REWIND(NFIL)
!
      LENG=3*NAT
      ALLOCATE(POSTRA(LENG*NSTEP))
      ALLOCATE(FORTRA(LENG*NSTEP))
      ALLOCATE(POS(LENG))
      IND=0
      DO ISTEP=1,NSTEP
        READ(NFIL)NFI,NDUMMY,(POS(I),I=1,LENG)
!       PRINT*,'POS ',NFI,(POS(I),I=1,LENG)
        DO I=1,LENG
          POSTRA(IND+I)=REAL(POS(I))
        ENDDO
        IND=IND+LENG
        IF(ISTEP.LT.NSTEP) THEN
          DO I=1,NSKIP
            READ(NFIL)
          ENDDO
        END IF
      ENDDO
      CALL FILEHANDLER$CLOSE('POSITION-TRAJECTORY')    
!
!     == CHECK THE NUMBER OF POINTS ON THE FORCE TRAJECTORY ============
      CALL FILEHANDLER$UNIT('FORCE-TRAJECTORY',NFIL)    
      REWIND NFIL
      NSTEPF=0
      DO WHILE (NSTEP1.LT.1000000)
        NSTEPF=NSTEPF+1
        READ(NFIL,END=2000)NFI,NAT
      ENDDO
 2000 CONTINUE
      NSTEPF=NSTEPF/(NSKIP+1)
      IF(NSTEPF.LE.2) THEN
        CALL FILEHANDLER$CLOSE('FORCE-TRAJECTORY')
        RETURN
      END IF
!
!     == READ FORCES ===================================================
      IF(NSTEPF.EQ.NSTEP) THEN
        REWIND(NFIL)
        READ(NFIL)NFI,NAT1
        IF(NAT.NE.NAT1) THEN
          CALL ERROR$MSG('NUMBER OF ATOMS INCONSISTENT')
          CALL ERROR$STOP('WRITEPATHFORDX')
        ENDIF
        REWIND(NFIL)
        IND=0
        DO ISTEP=1,NSTEP
          READ(NFIL)NFI,NDUMMY,(POS(I),I=1,LENG)
!         PRINT*,'FOR ',NFI,(POS(I),I=1,LENG)
          DO I=1,LENG
            FORTRA(IND+I)=REAL(POS(I))
          ENDDO
          IND=IND+LENG
          IF(ISTEP.LT.NSTEP) THEN
            DO I=1,NSKIP
              READ(NFIL)
            ENDDO
          END IF
        ENDDO
        CALL FILEHANDLER$CLOSE('FORCE-TRAJECTORY')    
      ELSE
        DO ISTEP=1,NSTEP
          DO I=1,LENG
            FORTRA(IND+I)=0.0
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(POS)
!
!     == READING FINISHED
      CALL FILEHANDLER$UNIT('PATH.DX',NFIL)
      CALL DX$SPAGETTI(NFIL,NAT,NSTEP,POSTRA,FORTRA)
!     CALL WRITEDX_PATH(NFIL,NAT,NSTEP,POSTRA,FORTRA)
      CALL FILEHANDLER$CLOSE('PATH.DX')
      DEALLOCATE(FORTRA)
      DEALLOCATE(POSTRA)
 9999 CONTINUE
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WRITEDXBOX
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TCHK_,TINI
      DIMENSION BOXORIGIN(3),BOXVECTORS(3,3)
      DIMENSION BOXORIGIN_(3),BOXVECTORS_(3,3)
      SAVE BOXORIGIN,BOXVECTORS
      DATA TINI/.FALSE./
!
      ENTRY WRITEDXBOX$SET(BOXORIGIN_,BOXVECTORS_)
        DO I=1,3
          BOXORIGIN(I)=BOXORIGIN_(I)
          DO J=1,3
            BOXVECTORS(I,J)=BOXVECTORS_(I,J)
          ENDDO
        ENDDO
        TINI=.TRUE.
        RETURN

      ENTRY WRITEDXBOX$BOX(BOXORIGIN_,BOXVECTORS_,TCHK_)
        IF(.NOT.TINI) THEN
          TCHK_=.FALSE.
          RETURN
        ENDIF
        DO I=1,3
          BOXORIGIN_(I)=BOXORIGIN(I)
          DO J=1,3
            BOXVECTORS_(I,J)=BOXVECTORS(I,J)
          ENDDO
        ENDDO
        TCHK_=.TRUE.
        RETURN
      END
!.......................................................................
MODULE DX_MODULE
PRIVATE
PUBLIC DX__SET
PUBLIC DX__GET
USE LINKEDLIST_MODULE
TYPE(LL_TYPE)        :: $LIST
LOGICAL(4)           :: TSET=.FALSE.
INTEGER(4)           :: NAT1
INTEGER(4)           :: NFRAME
LOGICAL(4)           :: TLATTICE=.FALSE.
INTERFACE DX$SET
#  MODULE TEMPLATE DX$SET
END INTERFACE DX$SET
INTERFACE DX$SET
#  MODULE TEMPLATE DX$GET
END INTERFACE DX$SET
!......................................................................
#MODULE DX$SET
(<TYPEID><TYPE>)=([R8][REAL(8)])([I4][INTEGER(4)])
(<RANKID><RANK>)=([R0][])([R1][(:)])([R2][(:,:)]
#BODY
!     
!     ..................................................................
      SUBROUTINE DX$SET<TYPEID>(IDENT_,NBYTE_,VAL_)
!     ******************************************************************
!     **  SET DATA                                                    **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      <TYPE>      ,INTENT(IN) :: VAL_<RANK>
!     ******************************************************************
      CALL LINKEDLIST$SET($LIST,IDENT_,0,VAL_)
      RETURN
      END
#END MODULE DX$SET
!......................................................................
#MODULE DX$GET
(<TYPEID><TYPE>)=([R8][REAL(8)])([I4][INTEGER(4)])
(<RANKID><RANK>)=([R0][])([R1][(:)])([R2][(:,:)]
#BODY
!     
!     ..................................................................
      SUBROUTINE DX$GET<TYPEID><RANKID>(IDENT_,NBYTE_,VAL_)
!     ******************************************************************
!     **  GET DATA                                                    **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: IDENT_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      <TYPE>      ,INTENT(OUT) :: VAL_<RANK>
!     ******************************************************************
      CALL LINKEDLIST$GET($LIST,IDENT_,1,VAL_)
      RETURN
      END
#END MODULE DX$GET
END MODULE DX_MODULE
!     
!     ..................................................................
      SUBROUTINE DX$BALLSTICKSTART(NFIL)
!     ******************************************************************
!     **  WRITE HEADER FOR BALLSTICK MODEL (MOVIE)                    **
!     ******************************************************************
      USE DX_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)              :: NAT
      REAL(8)     ,ALLOCATABLE:: RAD(:)        !(NAT)
      INTEGER(4)  ,ALLOCATABLE:: ICOLOR(:,:)   !(3,NAT)
      INTEGER(4)              :: ISPHERERADIUSOBJECT
      INTEGER(4)              :: ISPHERECOLOROBJECT
      INTEGER(4)              :: I,IATOM
      REAL(8)                 :: X,Y,Z,SVAR
!     ******************************************************************
      IF(TSET) THEN
        CALL ERROR$MSG('BALLSTICK HEADER ALREADY WRITTEN')
        CALL ERROR$STOP('WRITEDX_BALLSTICK$HEADER')
      END IF
      TSET=.TRUE.
      CALL LINKEDLIST$REPORT($LIST,6)
      CALL LINKEDLIST$GET($LIST,'NAT',4,NAT)
      ALLOCATE(RAD(NAT))
      ALLOCATE(ICOLOR(3,NAT))
      CALL LINKEDLIST$GET($LIST,'RAD',8*NAT,RAD)
      CALL LINKEDLIST$GET($LIST,'ICOLOR',4*3*NAT,ICOLOR)
      
      NAT1=NAT
      NFRAME=0
      ISPHERERADIUSOBJECT=1
      ISPHERECOLOROBJECT=2
!     
!     ================================================================
!     ==   DATA ARRAY: SIZE OF THE SPHERES                          ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"SPHERE SIZE"/"#")')
      WRITE(NFIL,FMT='("OBJECT ",I10, &
     &                /5X,"CLASS ARRAY" &
     &                ,5X,"TYPE FLOAT" &
     &                ,5X,"RANK 0" &
     &                ,5X,"ITEMS ",I10 &
     &                /"DATA FOLLOWS")')ISPHERERADIUSOBJECT,NAT
      WRITE(NFIL,FMT='(10F10.5)')(RAD(I),I=1,NAT)
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   COLORS                                                   ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"COLORS"/"#")')
      WRITE(NFIL,FMT='("OBJECT ",I10 &
     &            /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &            /"DATA FOLLOWS")')ISPHERECOLOROBJECT,NAT
!     
      DO IATOM=1,NAT
        X=DBLE(ICOLOR(1,IATOM))
        Y=DBLE(ICOLOR(2,IATOM))
        Z=DBLE(ICOLOR(3,IATOM))
        SVAR=200.D0
        X=X/SVAR
        Y=Y/SVAR
        Z=Z/SVAR
        WRITE(NFIL,FMT='(3F10.4)')X,Y,Z
      ENDDO
!     
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
      DEALLOCATE(RAD)
      DEALLOCATE(ICOLOR)
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE DX$BALLSTICKADD(NFIL)
!     ******************************************************************
!     **  WRITE FRAME                                                 **
!     ******************************************************************
      USE DX_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NAT
      REAL(8)   ,ALLOCATABLE:: POS(:,:)
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: NBOND
      INTEGER(4),ALLOCATABLE:: IBOND(:,:)  !(2,NBOND)
      REAL(8)               :: BOXR0(3)
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: IPOSITIONSOBJECT
      INTEGER(4)            :: IBONDSOBJECT    
      INTEGER(4)            :: ICELLOBJECT    
      INTEGER(4)            :: IMOLECULEOBJECT
      INTEGER(4)            :: ISPHERECOLOROBJECT
      INTEGER(4)            :: ISPHERERADIUSOBJECT
      INTEGER(4)            :: I1,I2,I3
      INTEGER(4)            :: I,J,IB,IA
      REAL(8)               :: T1,T2,T3
      REAL(8)               :: X,Y,Z
!     ******************************************************************
      IF(.NOT.TSET) THEN
        CALL ERROR$MSG('CALL ..$START BEFORE ..$ADDFRAME')
        CALL ERROR$STOP('DXBALLSTICK$ADDFRAME')
      END IF
!     
!     ================================================================
!     ==  COLLECT DATA FROM LINKEDLIST                              ==
!     ================================================================
      CALL LINKEDLIST$GET($LIST,'NAT',4,NAT)
      IF(NAT.NE.NAT1) THEN
        CALL ERROR$MSG('NUMBER OF BALLS MUST NOT CHANGE DURING MOVIE')
        CALL ERROR$STOP('WRITEDX_BALLSTICK$ADDFRAME')
      END IF
      ALLOCATE(POS(3,NAT))
      CALL LINKEDLIST$GET($LIST,'R',8*3*NAT,POS)
      
      CALL LINKEDLIST$EXIST($LIST,'NBOND',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET($LIST,'NBOND',4,NBOND)
      ELSE
        NBOND=0
        CALL LINKEDLIST$SET($LIST,'NBOND',4,NBOND)
      END IF
      IF(NBOND.EQ.0) THEN
        NBOND=1
        CALL LINKEDLIST$SET($LIST,'NBOND',4,NBOND)
        ALLOCATE(IBOND(2,NBOND))
        IBOND(1,1)=1
        IBOND(2,1)=1
        CALL LINKEDLIST$SET($LIST,'IBOND',4*2*NBOND,IBOND)
      ELSE
        ALLOCATE(IBOND(2,NBOND))
      END IF
!     CALL DX_BOND(.FALSE.,NAT,POS,RAD,NBOND,IBOND)
!     CALL DX_BOND(.TRUE.,NAT,POS,RAD,NBOND,IBOND)
      CALL LINKEDLIST$GET($LIST,'IBOND',4*2*NBOND,IBOND)
!     CALL LINKEDLIST$REPORT($LIST,6)
!     
      CALL LINKEDLIST$EXIST($LIST,'BOXR0',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET($LIST,'BOXR0',8*3,BOXR0)
      ELSE
        BOXR0(1)=0.D0
        BOXR0(2)=0.D0
        BOXR0(3)=0.D0
        CALL LINKEDLIST$SET($LIST,'BOXR0',8*3,BOXR0)
      END IF
!     
      CALL LINKEDLIST$EXIST($LIST,'RBAS',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET($LIST,'RBAS',8*3*3,RBAS)
      ELSE
        DO I=1,3
          DO J=1,3
            RBAS(I,J)=0.D0
          ENDDO
        ENDDO
        CALL LINKEDLIST$SET($LIST,'RBAS',8*3*3,RBAS)
      END IF
!     
!     ================================================================
!     ==  POSITIONS                                                 ==
!     ================================================================
      NFRAME=NFRAME+1
      IPOSITIONSOBJECT=(NFRAME-1)*4+2+1
      IBONDSOBJECT    =(NFRAME-1)*4+2+2
      ICELLOBJECT     =(NFRAME-1)*4+2+3
      IMOLECULEOBJECT =(NFRAME-1)*4+2+4
!     
!     ================================================================
!     ==   POSITIONS ARRAY: ATOMIC POSITIONS                        ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"ATOMIC POSITIONS"/"#")')
      WRITE(NFIL,FMT='("OBJECT ",I10 &
     &            /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &            /"DATA FOLLOWS")')IPOSITIONSOBJECT,NAT
      DO IA=1,NAT
          WRITE(NFIL,FMT='(3F10.5)')(POS(I,IA),I=1,3)
      ENDDO  
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   CONNECTIONS ARRAY: BONDS                                 ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"BONDS"/"#")')
      WRITE(NFIL,FMT='("OBJECT ",I10 &
     &            /"CLASS ARRAY TYPE INT RANK 1 SHAPE 2 ITEMS ",I10 &
     &            /"DATA FOLLOWS")')IBONDSOBJECT,NBOND
      WRITE(NFIL,FMT='(2I5)')((IBOND(I,IB)-1,I=1,2),IB=1,NBOND)
      WRITE(NFIL,FMT='("ATTRIBUTE ""REF"" STRING ""POSITIONS""")')
      WRITE(NFIL &
     &     ,FMT='("ATTRIBUTE ""ELEMENT TYPE"" STRING ""LINES""")')
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   BOX: LATTICE VECTORS                                     ==
!     ================================================================
      WRITE(NFIL,FMT='("#"/"#",T10,"LATTICE VECTORS"/"#")')
      WRITE(NFIL,FMT='("OBJECT ",I10 &
     &            /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &            /"DATA FOLLOWS")')ICELLOBJECT,8
      DO I3=0,1 
         T3=DBLE(I3)-0.5D0
        DO I2=0,1
          T2=DBLE(I2)-0.5D0
          DO I1=0,1
            T1=DBLE(I1)-0.5D0
            IF (TLATTICE)THEN
               T1=T1+0.5D0
               T2=T2+0.5D0
               T3=T3+0.5D0
            ENDIF
            X=RBAS(1,1)*T1+RBAS(1,2)*T2+RBAS(1,3)*T3
            Y=RBAS(2,1)*T1+RBAS(2,2)*T2+RBAS(2,3)*T3
            Z=RBAS(3,1)*T1+RBAS(3,2)*T2+RBAS(3,3)*T3
            WRITE(NFIL,FMT='(3F10.5)')X,Y,Z
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,FMT='("#")')
!     
!     ================================================================
!     ==   OBJECT MOLECULE:                                         ==
!     ================================================================
      WRITE(NFIL,FMT='("OBJECT ",I10)')IMOLECULEOBJECT
      WRITE(NFIL,FMT='("CLASS FIELD")')
      WRITE(NFIL,FMT='("COMPONENT ""DATA"" VALUE ",I10)') &
     &      ISPHERERADIUSOBJECT
      WRITE(NFIL,FMT='("COMPONENT ""POSITIONS"" VALUE ",I10)') &
     &      IPOSITIONSOBJECT
      WRITE(NFIL,FMT='("COMPONENT ""CONNECTIONS"" VALUE ",I10)') &
     &      IBONDSOBJECT
      WRITE(NFIL,FMT='("COMPONENT ""BOX"" VALUE ",I10)') &
     &      ICELLOBJECT
      WRITE(NFIL,FMT='("COMPONENT ""COLORS"" VALUE ",I10)') &
     &     ISPHERECOLOROBJECT
      WRITE(NFIL,FMT='("ATTRIBUTE ""NAME"" STRING ""ATOMS""")')
      WRITE(NFIL,FMT='("#")')
      DEALLOCATE(POS)
      DEALLOCATE(IBOND)
      RETURN
      END
!     
!     .................................................................
      SUBROUTINE BALLSTICKEND(NFIL)
!     ******************************************************************
!     **  WRITE END                                                   **
!     ******************************************************************
      USE DX_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IFRAME
      INTEGER(4)            :: IMOLECULEOBJECT
!     ******************************************************************
      IF(.NOT.TSET) THEN
        CALL ERROR$MSG('CALL ..$START BEFORE ..$END')
        CALL ERROR$STOP('DXBALLSTICK$END')
      END IF
      TSET=.FALSE.
!     
!     ================================================================
!     ==                                                            ==
!     ================================================================
      IF (NFRAME.NE.1) THEN 
        WRITE(NFIL,FMT='("OBJECT ""SERIES""")')
        WRITE(NFIL,FMT='("CLASS SERIES")')
        DO IFRAME=1,NFRAME
          IMOLECULEOBJECT =(IFRAME-1)*4+2+4
          WRITE(NFIL &
     &         ,FMT='("MEMBER ",I10," VALUE ",I10," POSITION ",I10)') &
     &         IFRAME-1,IMOLECULEOBJECT,IFRAME
        ENDDO
        WRITE(NFIL,FMT='("#")')
      END IF
      WRITE(NFIL,FMT='("END")')
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE DX$DENSITY(NFIL,NR1,NR2,NR3,DENSITY_)
!     ******************************************************************
!     **  PLOT A DENSITY                                              **
!     **    REQUIRES: RBAS                                            **
!     ******************************************************************
      USE DX_MODULE, ONLY : $LIST
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      INTEGER(4),INTENT(IN) :: NR3
      REAL(8)   ,INTENT(IN) :: DENSITY_
      REAL(8)               :: RBAS(3,3)
      REAL(8)               :: BOXR0(3)
      REAL(8)               :: BOXVEC(3,3)
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: I,J
!     ******************************************************************
!
!     ================================================================
!     ==  GET DATA FROM LINKED LIST                                 ==
!     ================================================================
!     == CLIP BOX; DEFAULT IS ONE UNIT CELL WITH ONE CORNER AT 0,0,0
      CALL LINKEDLIST$GET($LIST,'RBAS',8*3*3,RBAS)
      CALL LINKEDLIST$EXIST($LIST,'BOXR0',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET($LIST,'BOXR0',8*3,BOXR0)
      ELSE
        DO I=1,3
          BOXR0(I)=0.D0+1.D-8
        ENDDO
        CALL LINKEDLIST$SET($LIST,'BOXR0',8*3,BOXR0)
      END IF
      CALL LINKEDLIST$EXIST($LIST,'BOXVEC',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET($LIST,'BOXVEC',8*3*3,BOXVEC)
      ELSE
        DO I=1,3
          DO J=1,3
            BOXVEC(I,J)=RBAS(I,J)*(1.D0-1.D-8)
          ENDDO
        ENDDO
        CALL LINKEDLIST$SET($LIST,'BOXVEC',8*3*3,BOXVEC)
      END IF
!     
!     ================================================================
!     ==  GET DATA FROM LINKED LIST                                 ==
!     ================================================================
      CALL DX_DENSITY(NFIL,NR1,NR2,NR3,RBAS,DENSITY_,BOXR0,BOXVEC)
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE DX$SPAGETTI(NFIL,NAT_,NSTEP,RTRA_,FTRA_)
!     ******************************************************************
!     **  SPAGETTI PLOT                                               **
!     ******************************************************************
      USE DX_MODULE, ONLY : $LIST
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT_
      INTEGER(4),INTENT(IN) :: NSTEP
      REAL(8)   ,INTENT(IN) :: RTRA_(3,NAT_,NSTEP)
      REAL(8)   ,INTENT(IN) :: FTRA_(3,NAT_,NSTEP)
!     ******************************************************************
      CALL DX_PATH(NFIL,NAT_,NSTEP,RTRA_,FTRA_)
      RETURN
      END
!     
!     ..................................................................
      SUBROUTINE DX$EPLOT(NFIL)
!     ******************************************************************
!     **  COULOMB POTENTIAL OF A POINT CHARGE DISTRIBUTION            **
!     ******************************************************************
      USE DX_MODULE, ONLY : $LIST
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NAT
      REAL(8)   ,ALLOCATABLE:: Q(:)    !(NAT)
      REAL(8)               :: RBAS(3,3)
      REAL(8)               :: BOXR0(3)
      REAL(8)               :: BOXVEC(3,3)
      REAL(8)   ,ALLOCATABLE:: POS(:,:)   !(3,NAT)
      INTEGER(4)            :: N1,N2,N3
      REAL(8)   ,ALLOCATABLE:: EPOT(:,:,:)   !(N1,N2,N3)
      INTEGER(4)            :: I,J
      LOGICAL(4)            :: TCHK
!     ******************************************************************
!
!     ================================================================
!     ==  GET DATA FROM LINKED LIST                                 ==
!     ================================================================
      CALL LINKEDLIST$GET($LIST,'NAT',4,NAT)
      ALLOCATE(Q(NAT))
      CALL LINKEDLIST$GET($LIST,'POINTCHARGE',8*NAT,Q)
!     == CLIP BOX; DEFAULT IS ONE UNIT CELL WITH ONE CORNER AT 0,0,0
      CALL LINKEDLIST$GET($LIST,'RBAS',8*3*3,RBAS)
      CALL LINKEDLIST$EXIST($LIST,'BOXR0',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET($LIST,'BOXR0',8*3,BOXR0)
      ELSE
        DO I=1,3
          BOXR0(I)=0.D0+1.D-8
        ENDDO
        CALL LINKEDLIST$SET($LIST,'BOXR0',8*3,BOXR0)
      END IF
      CALL LINKEDLIST$EXIST($LIST,'BOXVEC',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET($LIST,'BOXVEC',8*3*3,BOXVEC)
      ELSE
        DO I=1,3
          DO J=1,3
            BOXVEC(I,J)=RBAS(I,J)*(1.D0-1.D-8)
          ENDDO
        ENDDO
        CALL LINKEDLIST$SET($LIST,'BOXVEC',8*3*3,BOXVEC)
      END IF
      ALLOCATE(POS(3,NAT))
      CALL LINKEDLIST$GET($LIST,'R',8*3*NAT,POS)
!     
!     ================================================================
!     ==  GET DATA FROM LINKED LIST                                 ==
!     ================================================================
      N1=30
      N2=30
      N3=30
      ALLOCATE(EPOT(N1,N2,N3))
      CALL DX_EPOT(NAT,POS,Q,BOXR0,BOXVEC,N1,N2,N3,EPOT)
      CALL DX_NDENSITY(NFIL,N1,N2,N3,RBAS,EPOT,BOXR0,BOXVEC)
      DEALLOCATE(EPOT)
      DEALLOCATE(Q)
      DEALLOCATE(POS)
      RETURN
      END
!                                                                       
!     .....................................................BONDS .......
      SUBROUTINE DX_NDENSITY(NFIL,NR1,NR2,NR3,RBAS,DENSITY,BOXR0,BOXVEC)
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                         
      LOGICAL TCHK
      REAL(8) RBAS(3,3)
      REAL(8) DENSITY(NR1,NR2,NR3)
      DIMENSION BOXR0(3),BOXVEC(3,3)
      CALL BOXBOX(RBAS,BOXR0,BOXVEC,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
      XMIN1=XMIN1*NR1-1.D0
      XMAX1=XMAX1*NR1+1.D0
      XMIN2=XMIN2*NR2-1.D0
      XMAX2=XMAX2*NR2+1.D0
      XMIN3=XMIN3*NR3-1.D0
      XMAX3=XMAX3*NR3+1.D0
      MIN1=INT(ABS(XMIN1))*INT(SIGN(1.D0,XMIN1))
      MAX1=INT(ABS(XMAX1))*INT(SIGN(1.D0,XMAX1))
      MIN2=INT(ABS(XMIN2))*INT(SIGN(1.D0,XMIN2))
      MAX2=INT(ABS(XMAX2))*INT(SIGN(1.D0,XMAX2))
      MIN3=INT(ABS(XMIN3))*INT(SIGN(1.D0,XMIN3))
      MAX3=INT(ABS(XMAX3))*INT(SIGN(1.D0,XMAX3))
!     N1=MAX1-MIN1+1
!     N2=MAX2-MIN2+1
!     N3=MAX3-MIN3+1
!     PRINT*,'BOXR0 ',BOXR0
!     PRINT*,'BOXVEC ',BOXVEC
!     PRINT*,'RBAS ',RBAS
!     PRINT*,'N1,N2,N3 ',N1,N2,N3
!     PRINT*,'MIN1,MIN2,MIN3 ',MIN1,MIN2,MIN3
!     PRINT*,'MAX1,MAX2,MAX3 ',MAX1,MAX2,MAX3
!
!     ==================================================================
!     ==================================================================
!     ==   WRITE DX FILE                                              ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     ==   DATA ARRAY: DENSITY                                        ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 1" &
     &                /5X,"CLASS ARRAY" &
     &                ,5X,"TYPE FLOAT" &
     &                ,5X,"RANK 0" &
     &                ,5X,"ITEMS ",I10 &
     &                /"DATA FOLLOWS")')NR1*NR2*NR3
      WRITE(NFIL,FMT='(10F12.5)') &
     &     (((DENSITY(I,J,K),K=1,NR1),J=1,NR2),I=1,NR3)
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   POSITIONS ARRAY: ATOMIC POSITIONS                          ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 2" &
     &              /"CLASS GRIDPOSITIONS COUNTS",3I10)') &
     &              NR1,NR2,NR3
      WRITE(NFIL,FMT='("ORIGIN",3F10.5)')BOXR0
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(BOXVEC(I,1)/DBLE(NR1),I=1,3)
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(BOXVEC(I,2)/DBLE(NR2),I=1,3)
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(BOXVEC(I,3)/DBLE(NR3),I=1,3)
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   CONNECTIONS ARRAY: CONNECTIONS                             ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 3" &
     &              /"CLASS GRIDCONNECTIONS COUNTS",3I10)') &
     &              NR1,NR2,NR3
      WRITE(NFIL,FMT='("ATTRIBUTE ""ELEMENT TYPE"" STRING ""CUBES""")')
      WRITE(NFIL,FMT='("ATTRIBUTE ""REF"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ================================================================
!     ==   CLIPBOX                                                  ==
!     ================================================================
      WRITE(NFIL,FMT='("OBJECT ",I10 &
     &            /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &            /"DATA FOLLOWS")')4,8
      DO I3=0,1 
         T3=DBLE(I3)
        DO I2=0,1
          T2=DBLE(I2)
          DO I1=0,1
            T1=DBLE(I1)
            X=BOXR0(1)+BOXVEC(1,1)*T1+BOXVEC(1,2)*T2+BOXVEC(1,3)*T3
            Y=BOXR0(2)+BOXVEC(2,1)*T1+BOXVEC(2,2)*T2+BOXVEC(2,3)*T3
            Z=BOXR0(3)+BOXVEC(3,1)*T1+BOXVEC(3,2)*T2+BOXVEC(3,3)*T3
            WRITE(NFIL,FMT='(3F10.5)')X,Y,Z
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   OBJECT MOLECULE:                                           ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT ""DENSITY""")')
      WRITE(NFIL,FMT='("CLASS FIELD")')
      WRITE(NFIL,FMT='("COMPONENT ""DATA"" VALUE 1")')
      WRITE(NFIL,FMT='("COMPONENT ""POSITIONS"" VALUE 2")')
      WRITE(NFIL,FMT='("COMPONENT ""CONNECTIONS"" VALUE 3")')
      WRITE(NFIL,FMT='("COMPONENT ""BOX"" VALUE 4")')
      WRITE(NFIL,FMT='("ATTRIBUTE ""NAME"" STRING ""DENSITY""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   END                                                        ==
!     ==================================================================
      WRITE(NFIL,FMT='("END")')
      RETURN
      END
!                                                                       
!     .....................................................BONDS .......
      SUBROUTINE DX_DENSITY(NFIL,NR1,NR2,NR3,RBAS,DENSITY,BOXR0,BOXVEC)
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                         
      LOGICAL TCHK
      REAL(8) RBAS(3,3)
      REAL(8) DENSITY(NR1,NR2,NR3)
      DIMENSION BOXR0(3),BOXVEC(3,3)
      CALL BOXBOX(RBAS,BOXR0,BOXVEC,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
      XMIN1=XMIN1*NR1-1.D0
      XMAX1=XMAX1*NR1+1.D0
      XMIN2=XMIN2*NR2-1.D0
      XMAX2=XMAX2*NR2+1.D0
      XMIN3=XMIN3*NR3-1.D0
      XMAX3=XMAX3*NR3+1.D0
      MIN1=INT(ABS(XMIN1))*INT(SIGN(1.D0,XMIN1))
      MAX1=INT(ABS(XMAX1))*INT(SIGN(1.D0,XMAX1))
      MIN2=INT(ABS(XMIN2))*INT(SIGN(1.D0,XMIN2))
      MAX2=INT(ABS(XMAX2))*INT(SIGN(1.D0,XMAX2))
      MIN3=INT(ABS(XMIN3))*INT(SIGN(1.D0,XMIN3))
      MAX3=INT(ABS(XMAX3))*INT(SIGN(1.D0,XMAX3))
      N1=MAX1-MIN1+1
      N2=MAX2-MIN2+1
      N3=MAX3-MIN3+1
!     PRINT*,'BOXR0 ',BOXR0
!     PRINT*,'BOXVEC ',BOXVEC
!     PRINT*,'RBAS ',RBAS
!     PRINT*,'N1,N2,N3 ',N1,N2,N3
!     PRINT*,'MIN1,MIN2,MIN3 ',MIN1,MIN2,MIN3
!     PRINT*,'MAX1,MAX2,MAX3 ',MAX1,MAX2,MAX3
!
!     ==================================================================
!     ==================================================================
!     ==   WRITE DX FILE                                              ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     ==   DATA ARRAY: DENSITY                                        ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 1" &
     &                /5X,"CLASS ARRAY" &
     &                ,5X,"TYPE FLOAT" &
     &                ,5X,"RANK 0" &
     &                ,5X,"ITEMS ",I10 &
     &                /"DATA FOLLOWS")')N1*N2*N3
      WRITE(NFIL,FMT='(10F12.5)') &
     &     (((DENSITY(MOD(I+1000*NR1,NR1)+1 &
     &               ,MOD(J+1000*NR2,NR2)+1 &
     &               ,MOD(K+1000*NR3,NR3)+1) &
     &               ,K=MIN3,MAX3) &
     &               ,J=MIN2,MAX2) &
     &               ,I=MIN1,MAX1)
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   POSITIONS ARRAY: ATOMIC POSITIONS                          ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 2" &
     &              /"CLASS GRIDPOSITIONS COUNTS",3I10)') &
     &              N1,N2,N3
      SH1=RBAS(1,1)/DBLE(NR1)*MIN1 &
     &   +RBAS(1,2)/DBLE(NR2)*MIN2 &
     &   +RBAS(1,3)/DBLE(NR3)*MIN3
      SH2=RBAS(2,1)/DBLE(NR1)*MIN1 &
     &   +RBAS(2,2)/DBLE(NR2)*MIN2 &
     &   +RBAS(2,3)/DBLE(NR3)*MIN3
      SH3=RBAS(3,1)/DBLE(NR1)*MIN1 &
     &   +RBAS(3,2)/DBLE(NR2)*MIN2 &
     &   +RBAS(3,3)/DBLE(NR3)*MIN3
      WRITE(NFIL,FMT='("ORIGIN",3F10.5)')SH1,SH2,SH3
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(RBAS(I,1)/DBLE(NR1),I=1,3)
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(RBAS(I,2)/DBLE(NR2),I=1,3)
      WRITE(NFIL,FMT='("DELTA",3F10.5)')(RBAS(I,3)/DBLE(NR3),I=1,3)
      WRITE(NFIL,FMT='("ATTRIBUTE ""DEP"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   CONNECTIONS ARRAY: CONNECTIONS                             ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT 3" &
     &              /"CLASS GRIDCONNECTIONS COUNTS",3I10)') &
     &              N1,N2,N3
      WRITE(NFIL,FMT='("ATTRIBUTE ""ELEMENT TYPE"" STRING ""CUBES""")')
      WRITE(NFIL,FMT='("ATTRIBUTE ""REF"" STRING ""POSITIONS""")')
      WRITE(NFIL,FMT='("#")')
!
!     ================================================================
!     ==   CLIPBOX                                                  ==
!     ================================================================
      WRITE(NFIL,FMT='("OBJECT ",I10 &
     &            /"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",I10 &
     &            /"DATA FOLLOWS")')4,8
      DO I3=0,1 
         T3=DBLE(I3)
        DO I2=0,1
          T2=DBLE(I2)
          DO I1=0,1
            T1=DBLE(I1)
            X=BOXR0(1)+BOXVEC(1,1)*T1+BOXVEC(1,2)*T2+BOXVEC(1,3)*T3
            Y=BOXR0(2)+BOXVEC(2,1)*T1+BOXVEC(2,2)*T2+BOXVEC(2,3)*T3
            Z=BOXR0(3)+BOXVEC(3,1)*T1+BOXVEC(3,2)*T2+BOXVEC(3,3)*T3
            WRITE(NFIL,FMT='(3F10.5)')X,Y,Z
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   OBJECT MOLECULE:                                           ==
!     ==================================================================
      WRITE(NFIL,FMT='("OBJECT ""DENSITY""")')
      WRITE(NFIL,FMT='("CLASS FIELD")')
      WRITE(NFIL,FMT='("COMPONENT ""DATA"" VALUE 1")')
      WRITE(NFIL,FMT='("COMPONENT ""POSITIONS"" VALUE 2")')
      WRITE(NFIL,FMT='("COMPONENT ""CONNECTIONS"" VALUE 3")')
      WRITE(NFIL,FMT='("COMPONENT ""BOX"" VALUE 4")')
      WRITE(NFIL,FMT='("ATTRIBUTE ""NAME"" STRING ""DENSITY""")')
      WRITE(NFIL,FMT='("#")')
!
!     ==================================================================
!     ==   END                                                        ==
!     ==================================================================
      WRITE(NFIL,FMT='("END")')
      RETURN
      END
!                                                                       
!     .....................................................WRITE .......
      SUBROUTINE DX_PATH(NFIL,NAT,NSTEP,POSI,FORC)
!     **                                                              **
!     **  WRITE INPUT FOR THE DATAEXPLORER FOR PLOTTING TRAJECTORIES  **
!     **  AND FORCES                                                  **
!     **                                                              **
!     **  WARNING! POSI AND FORC ARE SINGLE PRECISION ARRAYS!!!       **
!     **                                             P. MARGL         **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      REAL(4)  POSI(3,NAT,NSTEP),FORC(3,NAT,NSTEP)
      REWIND NFIL
!
!     ================================================================
!     ==  WRITE LINES                                               ==
!     ================================================================
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'#  CONNECTIONS (COMMON FOR ALL ATOMS)'
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'OBJECT',0
      WRITE(NFIL,*)'CLASS GRIDCONNECTIONS COUNTS 1'
      WRITE(NFIL,*)'ATTRIBUTE "REF" STRING "POSITIONS"'
      WRITE(NFIL,*)'ATTRIBUTE "ELEMENT TYPE" STRING "LINES"'
      WRITE(NFIL,*)'#'
!
!     ================================================================
!     ==   LOPOP OVER ATOMS                                         ==
!     ================================================================
      DO IAT=1,NAT     
        IOBJECT1=(IAT-1)*3+1
        IOBJECT2=(IAT-1)*3+2
        IOBJECT3=(IAT-1)*3+3
!
!       ================================================================
!       ==   WRITE FORCES (DATA)                                      ==
!       ================================================================
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'#  FORCES OF ATOM ',IAT
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'OBJECT ',IOBJECT1
        WRITE(NFIL,*)'CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3', &
     & ' ITEMS',NSTEP
        WRITE(NFIL,*)'DATA FOLLOWS'
        DO ISTEP=1,NSTEP
          WRITE(NFIL,FMT='(1X,3F15.9)')(FORC(I,IAT,ISTEP),I=1,3)
        ENDDO
        WRITE(NFIL,*)'ATTRIBUTE "DEP" STRING "POSITIONS"'
        WRITE(NFIL,*)'#'
!
!       ================================================================
!       ==  WRITE THE PATHPOSITIONS                                   ==
!       ================================================================
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'#  POSITIONS OF ATOM ',IAT
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'OBJECT',IOBJECT2
        WRITE(NFIL,*)'CLASS ARRAY   TYPE FLOAT  RANK 1 SHAPE 3', &
     &  ' ITEMS',NSTEP
        WRITE(NFIL,*)'DATA FOLLOWS'
        DO ISTEP=1,NSTEP
          WRITE(NFIL,FMT='(1X,3F15.9)')(POSI(I,IAT,ISTEP),I=1,3)
        ENDDO
        WRITE(NFIL,*)'ATTRIBUTE "DEP" STRING "POSITIONS"'
        WRITE(NFIL,*)'#'
!
!       ================================================================
!       ==                                                            ==
!       ================================================================
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'#  FIELD OBJECT OF ATOM ',IAT
        WRITE(NFIL,*)'#'
        WRITE(NFIL,*)'OBJECT',IOBJECT3
        WRITE(NFIL,*)'CLASS FIELD'
        WRITE(NFIL,*)'COMPONENT "DATA" VALUE ',IOBJECT1
        WRITE(NFIL,*)'COMPONENT "POSITIONS" VALUE ',IOBJECT2
        WRITE(NFIL,*)'COMPONENT "CONNECTIONS" VALUE ',0
        WRITE(NFIL,*)'ATTRIBUTE "NAME" STRING "PATH"  '
        WRITE(NFIL,*)'#'
      ENDDO
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'#  COMPOSITE FIELD FOR PATH'
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'OBJECT ',3*NAT+1,' COMPOSITEFIELD'
      DO IAT=1,NAT
        WRITE(NFIL,*)'MEMBER ',IAT-1,' VALUE ',3*IAT
      ENDDO
      WRITE(NFIL,*)'#'
      WRITE(NFIL,*)'END'
!
      RETURN
      END
!
!     .....................................................
      SUBROUTINE DX_BONDS(TDO,NAT,R,RAD,NBOND,IBOND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TDO
      DIMENSION R(3,NAT),RAD(NAT)
      DIMENSION IBOND(2,NBOND)
      DATA SCALE/1.2D0/
      IB=0
      DO IAT1=1,NAT
        RAD1=RAD(IAT1)
        DO IAT2=IAT1+1,NAT
          DIS=0.D0
          DO I=1,3
            DIS=DIS+(R(I,IAT1)-R(I,IAT2))**2
          ENDDO
          DIS=DSQRT(DIS)-(RAD1+RAD(IAT2))*SCALE
          IF(DIS.LT.0.D0) THEN
            IB=IB+1
            IF(TDO) THEN
              IF(NBOND.GE.IB) THEN
                IBOND(1,IB)=IAT1
                IBOND(2,IB)=IAT2
              ELSE
                CALL ERROR$MSG('NUMBER OF BONDS EXCEEDS NBOND')
                I=0
                IF(TDO) I=1                
                CALL ERROR$I4VAL('TDO',I)
                CALL ERROR$I4VAL('IB',IB)
                CALL ERROR$I4VAL('NBOND',NBOND)
                CALL ERROR$STOP('DX_BONDS')
              END IF
            END IF
          END IF
        ENDDO
      ENDDO
      IF(TDO) THEN
        IF(NBOND.NE.IB) THEN
          CALL ERROR$MSG('NBOND INCONSISTENT')
          CALL ERROR$STOP('DX_BONDS')
        END IF
      ELSE
        NBOND=IB
      ENDIF
      RETURN
      END

!     ...........................................................EPLOT..
      SUBROUTINE DX_EPOT(NAT,R,Q,BOXO,BOXVEC,N1,N2,N3,POT)
!     **                                                              **
!     **  CALCULATES THE ELECTROSTATIC POTENTIAL OF A POINT CHARGE    **
!     **  DISTRIBUTION Q ON THE POSITIONS R ON A GRID SPECIFIED BY    **
!     **  AN ORIGIN BOXO, THREE DISPLACEMENT VECTORS DR AND THE       **
!     **  NUMBER OF DISPLACEMENTS IN THE THREE SPACIAL DIRECTIONS     **
!     **  N1,N2,N3.                                                   **
!     **                                                              **
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(3,NAT),Q(NAT)
      DIMENSION BOXO(3),BOXVEC(3,3)    
      DIMENSION POT(N1,N2,N3)
      PRINT*,'NAT',NAT
      PRINT*,'R  ',R
      PRINT*,'Q  ',Q
      PRINT*,'BOXO   ',BOXO
      PRINT*,'BOXVEC ',BOXVEC
      PRINT*,'N1,N2,N3 ',N1,N2,N3
      DO I1=1,N1
        T1=DBLE(I1-1)/DBLE(N1-1)
        DO I2=1,N2
          T2=DBLE(I2-1)/DBLE(N2-1)
          DO I3=1,N3
            T3=DBLE(I3-1)/DBLE(N3-1)
            XM=BOXO(1)+BOXVEC(1,1)*T1+BOXVEC(1,2)*T2+BOXVEC(1,3)*T3
            YM=BOXO(2)+BOXVEC(2,1)*T1+BOXVEC(2,2)*T2+BOXVEC(2,3)*T3
            ZM=BOXO(3)+BOXVEC(3,1)*T1+BOXVEC(3,2)*T2+BOXVEC(3,3)*T3
            SVAR=0.D0
            DO IAT=1,NAT
              DX=R(1,IAT)-XM
              DY=R(2,IAT)-YM
              DZ=R(3,IAT)-ZM
              SVAR=SVAR-Q(IAT)/DSQRT(DX*DX+DY*DY+DZ*DZ)
            ENDDO
            POT(I1,I2,I3)=SVAR
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

