program main
call writeout
stop
end
!
!     ..................................................................
      SUBROUTINE WRITEOUT
!     **                                                              **
      IMPLICIT NONE
      real(8)              :: r1=1.056d-4
      real(8)              :: dex=0.05d0
      integer(4)           :: nr=250
      INTEGER(4)           :: NFIL=6
      REAL(8)              :: PSZ=0.D0   !DUMMY: NOT USED
      REAL(8)              :: aez=0.D0   
      REAL(8)              :: rcsmall=0.1
      integer(4)           :: nwave=9
      INTEGER(4)           :: IWAVE
      integer(4),allocatable :: lphi(:)
      real(8)   ,allocatable :: vhat(:)
      real(8)   ,allocatable :: aecore(:)
      real(8)   ,allocatable :: pscore(:)
      real(8)   ,allocatable :: dtkin(:,:)
      real(8)   ,allocatable :: do(:,:)
      real(8)   ,allocatable :: pro(:,:)
      real(8)   ,allocatable :: aephi(:,:,:)
      real(8)   ,allocatable :: psphi(:,:,:)
      real(8)              :: ri,xexp
      integer(4)           :: ir
!     ******************************************************************
      nfil=5
      read(NFIL,FMT='(F15.10,F10.5,2I4,2F5.2,F15.12)') &
     &                 R1,DEX,NR,NWAVE,PSZ,AEZ,RCSMALL
      allocate(lphi(nwave))
      allocate(vhat(nr))
      allocate(aecore(nr))
      allocate(pscore(nr))
      allocate(dtkin(nwave,nwave))
      allocate(do(nwave,nwave))
      allocate(pro(nr,nwave))
      allocate(psphi(nr,3,nwave))
      allocate(aephi(nr,3,nwave))
      READ(NFIL,FMT='(14I5)')LPHI(:)
      READ(NFIL,FMT='(SP,5E14.8)')VHAT(:)
!     ====  AECORE = CORE CHARGE DENSITY  ==============================
      READ(NFIL,FMT='(SP,5E14.8)')AECORE(:)
!     ====  PSCORE = PSEUDIZED CHARGE DENSITY===========================
      READ(NFIL,FMT='(SP,5E14.8)')PSCORE(:)
!     ====  DTKIN = <AEPHI|-DELTA/2|AEPHI> - <PSPHI|-DELTA/2|PSPHI> ====
      READ(NFIL,FMT='(SP,5E14.8)')DTKIN(:,:)
!     ====  DOVER = <AEPHI|AEPHI> - <PSPHI|PSPHI> ======================
      READ(NFIL,FMT='(SP,5E14.8)')DO(:,:)
      DO IWAVE=1,NWAVE
        READ(NFIL,FMT='(SP,5E14.8)')PRO(:,IWAVE)
        READ(NFIL,FMT='(SP,5E14.8)')AEPHI(:,1,IWAVE)
        READ(NFIL,FMT='(SP,5E14.8)')PSPHI(:,1,IWAVE)
      ENDDO
!=====================================================================
      aez=1.d0
      aecore=0.d0
      pscore=0.d0
      vhat=0.d0
      do=0.d0
      dtkin=0.d0
      psphi=0.d0
      aephi=psphi
!=====================================================================
       nfil=6
!      CALL FILEHANDLER$UNIT('SETUPO',NFIL)
      WRITE(NFIL,FMT='(F15.10,F10.5,2I4,2F5.2,F15.12)') &
     &                 R1,DEX,NR,NWAVE,PSZ,AEZ,RCSMALL
      WRITE(NFIL,FMT='(14I5)')LPHI(:)
      WRITE(NFIL,FMT='(SP,5E14.8)')VHAT(:)
!     ====  AECORE = CORE CHARGE DENSITY  ==============================
      WRITE(NFIL,FMT='(SP,5E14.8)')AECORE(:)
!     ====  PSCORE = PSEUDIZED CHARGE DENSITY===========================
      WRITE(NFIL,FMT='(SP,5E14.8)')PSCORE(:)
!     ====  DTKIN = <AEPHI|-DELTA/2|AEPHI> - <PSPHI|-DELTA/2|PSPHI> ====
      WRITE(NFIL,FMT='(SP,5E14.8)')DTKIN(:,:)
!     ====  DOVER = <AEPHI|AEPHI> - <PSPHI|PSPHI> ======================
      WRITE(NFIL,FMT='(SP,5E14.8)')DO(:,:)
      DO IWAVE=1,NWAVE
        WRITE(NFIL,FMT='(SP,5E14.8)')PRO(:,IWAVE)
        WRITE(NFIL,FMT='(SP,5E14.8)')AEPHI(:,1,IWAVE)
        WRITE(NFIL,FMT='(SP,5E14.8)')PSPHI(:,1,IWAVE)
      ENDDO
!     CALL FILEHANDLER$CLOSE('SETUPO')
      RETURN
      END
