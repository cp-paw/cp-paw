!************************************************************************
!  DEMONSTRATOR CODE FOR DMFT
!  CODE WRITTEN FOR A SINGLE SITE 
!************************************************************************
MODULE DMFT_MODULE
LOGICAL(4)             :: TON=.FALSE.
LOGICAL(4)             :: TINI=.FALSE.
INTEGER(4)             :: NOMEGA
REAL(8)   ,ALLOCATABLE :: OMEGA(:)
INTEGER(4)             :: NDIM=2
INTEGER(4),ALLOCATABLE :: NORB
INTEGER(4),ALLOCATABLE :: ORBPRO(:,:)     !(NPRO,NORB)
COMPLEX(8),ALLOCATABLE :: GREEN(:,:,:,:,:)  ! WEISS FUNCTION
COMPLEX(8),ALLOCATABLE :: HAMIL(:,:,:,:,:)
END MODULE DMFT_MODULE
!
!      .................................................................
       SUBROUTINE DMFT_INITIALIZE
!      *****************************************************************
!      **  EVALUATES THE COEFFICIENTS OF THE LOCAL ORBITALS           **
!      *****************************************************************
       USE DMFT_MODULE
       IMPLICIT NONE
       REAL(8)    :: OMEGAMIN,OMEGAMAX
       REAL(8)    :: SVAR
       INTEGER(4) :: I
!      *****************************************************************
       TINI=.TRUE.
       IF(NOMEGA.EQ.0) THEN        
         CALL ERROR$MSG('NOMEGA NOT DEFINED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('DMFT$SETI4')
       END IF
!
!      =================================================================
!      ==  EVALUATE OMEGA GRID                                        ==
!      =================================================================
       NOMEGA=100
       OMEGAMIN=-30.D0
       OMEGAMAX=10.D0
       ALLOCATE(OMEGA(NOMEGA))
       SVAR=(OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1,KIND=8)
       DO I=1,NOMEGA
         OMEGA(I)=OMEGAMIN+REAL(I-1,KIND=8)*SVAR
       ENDDO
!
!      =================================================================
!      ==  DEFINE ORBITAL PROJETORS                                   ==
!      =================================================================
       ALLOCATE(ORBPRO(NPRO,NORB))
       ORBPRO(:,:)=0.D0
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DMFT$SETI4(ID,VAL)
       USE DMFT_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(IN) :: VAL
!      *****************************************************************
       IF(ID.EQ.'NOMEGA') THEN
         NOMEGA=VAL
         CALL ERROR$MSG('NOMEGA IS STILL HARD WIRED; DO NOT SET')
         CALL ERROR$STOP('DMFT$SETI4')
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('DMFT$SETI4')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DMFT$SETL4(ID,VAL)
       USE DMFT_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       LOGICAL(4)  ,INTENT(IN) :: VAL
!      *****************************************************************
       IF(ID.EQ.'ON') THEN
         TON=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('DMFT$SETI4')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DMFT$SETPROJ(NDIM,NB,NPRO,PROJ,HAMILTON)
!      *****************************************************************
!      **  PROJECTS ONTO WAVE FUNCTIONS ONTO LOCAL ORBITALS           **
!      **  AND ADDS RESULT TO WEISS FUNCTION                          **
!      *****************************************************************
       USE DMFT_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NDIM       ! #(BANDS)
       INTEGER(4),INTENT(IN) :: NB         ! #(BANDS)
       INTEGER(4),INTENT(IN) :: NPRO       ! #(PARTIAL WAVES)
       INTEGER(4),INTENT(IN) :: NORB       ! #(LOCALIZED ORBITALS)
       COMPLEX(8),INTENT(IN) :: HAMILTON(NB,NB)
       COMPLEX(8),INTENT(IN) :: PROJ(NPRO,NB)
       COMPLEX(8)            :: ORBCOEFF(NDIM,NB,NORB)
       REAL(8)               :: EPSILON(NB)
       COMPLEX(8)            :: U(NB,NB)
!      *****************************************************************
       IF(.NOT.TON) RETURN
       IF(.NOT.TINI) CALL DMFT_INITIALIZE
       IF(.NOT.ALLOCATED(GREEN)) THEN
         ALLOCATE(GREEN(NOMEGA,NDIM,NORB,NDIM,NORB))
         ALLOCATE(GREEN(HAMIL,NDIM,NORB,NDIM,NORB))
         GREEN(:,:,:,:,:)=(0.D0,0.D0)
         HAMIL(:,:,:,:,:)=(0.D0,0.D0)
       END IF
!
!      =================================================================
!      ==  PROJECT WAVE FUNCTION ONTO LOCAL ORBITALS                  ==
!      =================================================================
       CALL DMFT_PROJECTIONS(NB,NPRO,NORB,PROJ,ORBPRO,ORBCOEFF)
!
!      =================================================================
!      ==  TRANSFORM TO EIGENSTATES                                   ==
!      =================================================================
       CALL LIB$DIAGC8(NB,HAMILTON,EPSILON,U)
       DO IORB=1,NORB
         DO IDIM=1,NDIM
           ORBCOFF(IDIM,:,IORB)=MATMUL(ORBCOFF(IDIM,:,IORB),U)
         ENDDO
       ENDDO
!
!      =================================================================
!      ==  ADD TO GREENS FUNCTION                                     ==
!      =================================================================
       CALL DMFT_ADDGREEN(NB,NORB,EPSILON,ORBCOEFF,NOMEGA,OMEGA,GREEN,HAMIL)
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DMFT$ETOT
!      *****************************************************************
!      **  PROJECTS ONTO WAVE FUNCTIONS ONTO LOCAL ORBITALS           **
!      **  AND ADDS RESULT TO WEISS FUNCTION                          **
!      *****************************************************************
       USE DMFT_MODULE
       IMPLICIT NONE
       COMPLEX(8),ALLOCATABLE :: WEISS(NDIM,NORB,NDIM,NORB,NOMEGA)
!      *****************************************************************
       IF(.NOT.TON) RETURN
       ALLOCATE(WEISS(NDIM,NORB,NDIM,NORB,NOMEGA)
       CALL DMFT_WEISS(NDIM*NORB,NOMEGA,GREEN,HAMIL,WEISS)
       DEALLOCATE(GREEN)
       DEALLOCATE(HAMIL)
       START=TRUE

       RETURN
       END
!
!      .................................................................
       SUBROUTINE DMFT_PROJECTIONS(NB,NPRO,NORB,PROJ,ORBPRO,ORBCOEFF)
!      *****************************************************************
!      **  EVALUATES THE COEFFICIENTS OF THE LOCAL ORBITALS           **
!      **  FOR EACH ONE-PARTICLE STATE                                **
!      **                                                             **
!      **      |PSI> = |CHI>*ORBPRO*<PTILDE|PSITILDE>                 **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB         ! #(BANDS)
       INTEGER(4),INTENT(IN) :: NPRO       ! #(PARTIAL WAVES)
       INTEGER(4),INTENT(IN) :: NORB       ! #(LOCALIZED ORBITALS)
       COMPLEX(8),INTENT(IN) :: PROJ(NDIM,NBH,NPRO)  ! <P|PSITILDE>
       REAL(8)   ,INTENT(IN) :: ORBPRO(NPRO,NORB)    !
       COMPLEX(8),INTENT(OUT):: ORBCOEFF(NDIM,NB,NORB)
       INTEGER(4)            :: IORB,IPRO,IB,IDIM
!      *****************************************************************
       ORBCOEFF(:,:,:)=(0.D0,0.D0)
       DO IORB=1,NORB
         DO IPRO=1,NPRO
           SVAR=ORBPRO(IPRO,IORB)
           IF(SVAR.EQ.0.D0) CYCLE
           DO IB=1,NB
             DO IDIM=1,NDIM
               ORBCOEFF(IDIM,IB,IORB)=ORBCOEFF(IDIM,IB,IORB) &
      &                              +PROJ(IDIM,IB,IPRO)*SVAR
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       RETURN
       END
!
!      .................................................................
       CALL DMFT_ADDGREEN(NB,NORB,EPSILON,ORBCOEFF,NOMEGA,OMEGA,GREEN,HAMIL)
!      *****************************************************************
!      **  EVALUATES THE MATRIX ELEMENTS OF TEH GREENS FUNCTION       **
!      **  AND ITS INVERSE BETWEEN THE LOCAL ORBITALS                 **
!      **                                                             **
!      *****************************************************************
       INTEGER(4),INTENT(IN) :: NB         ! #(BANDS)
       INTEGER(4),INTENT(IN) :: NORB       ! #(LOCALIZED ORBITALS)
       REAL(8)   ,INTENT(IN) :: EPSILON(NB)
       COMPLEX(8),INTENT(IN) :: ORBCOEFF(NDIM,NB,NORB)
       INTEGER(4),INTENT(IN) :: NOMEGA
       REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
       COMPLEX(8),INTENT(OUT):: GREEN(NORB,NORB,NOMEGA)
       COMPLEX(8),INTENT(OUT):: HAMILTON(NORB,NORB,NOMEGA)
       COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
       INTEGER(4)            :: IO1,IO2
       INTEGER(4)            :: IDIM1,IDIM2
       INTEGER(4)            :: IB,IOMEGA
       COMPLEX(8)            :: CSVAR,DE
!      *****************************************************************
       DO IO2=1,NORB
         DO IDIM2=1,NDIM
           DO IO1=1,NORB
             DO IDIM1=1,NDIM
               DO IB=1,NB
                 CSVAR=CONJG(ORBCOEFF(IDIM1,IB,IO1))*ORBCOEFF(IDIM2,IB,IO2)
                 DO IOMEGA=1,NOMEGA
                   DE=CI*OMEGA-EPS(IB)+MU         
                   GREEN(IOMEGA,IDIM1,IO1,IDIM2,IO2)=GREEN(IOMEGA,IDIM1,IO1,IDIM2,IO2)
      &                                             +CSVAR/DE
                   HAMIL(IOMEGA,IDIM1,IO1,IDIM2,IO2)=HAMIL(IOMEGA,IDIM1,IO1,IDIM2,IO2)
      &                                             +CSVAR*DE
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DMFT_WEISS(NORB,NOMEGA,GREEN,HAMIL,WEISS)
!      *****************************************************************
!      **  EVALUATES THE MATRIX ELEMENTS OF TEH GREENS FUNCTION       **
!      **  AND ITS INVERSE BETWEEN THE LOCAL ORBITALS                 **
!      **                                                             **
!      *****************************************************************
       INTEGER(4),INTENT(IN) :: NORB       ! #(LOCALIZED ORBITALS)
       INTEGER(4),INTENT(IN) :: NOMEGA
        COMPLEX(8),INTENT(IN) :: GREEN(NORB,NORB,NOMEGA)
       COMPLEX(8),INTENT(IN) :: HAMILTON(NORB,NORB,NOMEGA)
       COMPLEX(8),INTENT(OUT):: WEISS(NORB,NORB,NOMEGA)
       INTEGER(4)            :: IO
       INTEGER(4)            :: IOMEGA
       COMPLEX(8)            :: ONE(NORB,NORB)
       COMPLEX(8)            :: GH(NORB,NORB)
       COMPLEX(8)            :: GH2(NORB,NORB)
!      *****************************************************************
       ONE(:,:)=(0.D0,0.D0)
       DO IO=1,NORB
         ONE(IO,IO)=(1.D0,0.D0)
       ENDDO
       DO IOMEGA=1,NOMEGA
         GH(:,:)=MATMUL(GREEN(:,:,IOMEGA),HAMIL(:,:,IOMEGA))
         GH2(:,:)=MATMUL(GH,GH)
         WEISS(:,:,IOMEGA)=MATMUL(HAMIL(:,:,IOMEGA),3.D0*ONE-3.D0*GH+GH2)
       ENDDO
       RETURN
       END
