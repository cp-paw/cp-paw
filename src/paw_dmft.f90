!************************************************************************
!  DEMONSTRATOR CODE FOR dmft
!  CODE WRITTEN FOR A SINGLE SITE 
!************************************************************************
module dmft_module
logical(4)             :: TON=.FALSE.
logical(4)             :: TINI=.FALSE.
INTEGER(4)             :: NOMEGA
REAL(8)   ,ALLOCATABLE :: OMEGA(:)
INTEGER(4)             :: NDIM=2
INTEGER(4),ALLOCATABLE :: NORB
INTEGER(4),ALLOCATABLE :: ORBPRO(:,:)     !(NPRO,NORB)
COMPLEX(8),ALLOCATABLE :: GREEN(:,:,:,:,:)  ! WEISS FUNCTION
COMPLEX(8),ALLOCATABLE :: HAMIL(:,:,:,:,:)
end module dmft_module
!
!      .................................................................
       subroutine dmft_INITIALIZE
!      *****************************************************************
!      **  evaluates the coefficients of the local orbitals           **
!      *****************************************************************
       USE DMFT_MODULE
       IMPLICIT NONE
       REAL(8)    :: OMEGAMIN,OMEGAMAX
       REAL(8)    :: SVAR
       INTEGER(4) :: I
!      *****************************************************************
       TINI=.TRUE.
       IF(nOMEGA.eq.0) then        
         CALL ERROR$MSG('NOMEGA NOT DEFINED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$stop('dmft$seti4')
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
       subroutine dmft$setI4(ID,VAL)
       USE DMFT_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(IN) :: VAL
!      *****************************************************************
       IF(ID.EQ.'NOMEGA') THEN
         NOMEGA=VAL
         CALL ERROR$MSG('NOMEGA IS STILL HARD WIRED; DO NOT SET')
         CALL ERROR$stop('dmft$seti4')
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$stop('dmft$seti4')
       END IF
       RETURN
       END
!
!      .................................................................
       subroutine dmft$setL4(ID,VAL)
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
         CALL ERROR$stop('dmft$seti4')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DMFT$SETPROJ(NDIM,NB,NPRO,PROJ,Hamilton)
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
       COMPLEX(8),INTENT(in) :: hamilton(nb,nb)
       COMPLEX(8),INTENT(in) :: proj(npro,nb)
       COMPLEX(8)            :: ORBCOEFF(NDIM,NB,NORB)
       real(8)               :: epsilon(nb)
       complex(8)            :: u(nb,nb)
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
!      ==  project wave function onto local orbitals                  ==
!      =================================================================
       CALL DMFT_PROJECTIONS(NB,NPRO,NORB,PROJ,ORBPRO,ORBCOEFF)
!
!      =================================================================
!      ==  transform to eigenstates                                   ==
!      =================================================================
       call lib$diagc8(nb,hamilton,epsilon,u)
       do iorb=1,norb
         do idim=1,ndim
           orbcoff(idim,:,iorb)=matmul(orbcoff(idim,:,iorb),u)
         enddo
       enddo
!
!      =================================================================
!      ==  add to greens function                                     ==
!      =================================================================
       CALL DMFT_ADDGREEN(NB,NORB,epsilon,ORBCOEFF,nomega,omega,GREEN,HAMIL)
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
       complex(8),allocatable :: weiss(ndim,norb,ndim,norb,nomega)
!      *****************************************************************
       IF(.NOT.TON) RETURN
       allocate(weiss(ndim,norb,ndim,norb,nomega)
       CALL DMFT_WEISS(NDIM*NORB,NOMEGA,GREEN,HAMIL,WEISS)
       deallocate(green)
       deallocate(hamil)
       start=true

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
       CALL DMFT_ADDGREEN(NB,NORB,epsilon,ORBCOEFF,nomega,omega,GREEN,HAMIL)
!      *****************************************************************
!      **  EVALUATES THE matrix elements of teh greens function       **
!      **  and its inverse between the local orbitals                 **
!      **                                                             **
!      *****************************************************************
       INTEGER(4),INTENT(IN) :: NB         ! #(BANDS)
       INTEGER(4),INTENT(IN) :: NORB       ! #(LOCALIZED ORBITALS)
       real(8)   ,INTENT(in) :: epsilon(nb)
       COMPLEX(8),INTENT(in) :: ORBCOEFF(NDIM,NB,NORB)
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
       SUBROUTINE DMFT_weiss(NORB,nomega,GREEN,HAMIL,weiss)
!      *****************************************************************
!      **  EVALUATES THE matrix elements of teh greens function       **
!      **  and its inverse between the local orbitals                 **
!      **                                                             **
!      *****************************************************************
       INTEGER(4),INTENT(IN) :: NORB       ! #(LOCALIZED ORBITALS)
       INTEGER(4),INTENT(IN) :: NOmega
        COMPLEX(8),INTENT(in) :: GREEN(NORB,NORB,NOMEGA)
       COMPLEX(8),INTENT(in) :: HAMILTON(NORB,NORB,NOMEGA)
       COMPLEX(8),INTENT(OUT):: weiss(NORB,NORB,NOMEGA)
       INTEGER(4)            :: IO
       INTEGER(4)            :: iomega
       COMPLEX(8)            :: one(norb,norb)
       COMPLEX(8)            :: gh(norb,norb)
       COMPLEX(8)            :: gh2(norb,norb)
!      *****************************************************************
       ONE(:,:)=(0.D0,0.D0)
       DO IO=1,NORB
         ONE(IO,IO)=(1.D0,0.d0)
       ENDDO
       DO IOMEGA=1,NOMEGA
         GH(:,:)=MATMUL(GREEN(:,:,IOMEGA),HAMIL(:,:,IOMEGA))
         GH2(:,:)=MATMUL(GH,GH)
         WEISS(:,:,IOMEGA)=MATMUL(HAMIL(:,:,IOMEGA),3.D0*ONE-3.D0*GH+GH2)
       ENDDO
       RETURN
       END
