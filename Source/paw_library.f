!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  INTERFACES TO SCIENTIFIC LIBRARY ROUTINES                        **
!**  USES THE ESSL LIBRARY                                            **
!**                                                                   **
!***********************************************************************
!***********************************************************************

#IF DEFINED(CPPVARIABLE_ABS)
!=============== ENVIRONMENT ABSOFT COMPILER =========================
!    IT ASSUMES THE ATLAS LAPACK AND BLAS ROUTINES TO BE LINKED
!    IT ASSUMES THE ABSOFT SUPPORT LIBRARY LIBU77.A  TO BE LINKED
!    IT ASSUMES THE FAST FOURIER TRANSFORM LIBRARY LIBFFTW.A  TO BE LINKED
#DEFINE CPPVAR_FFTW
#DEFINE CPPVAR_ATLAS
#DEFINE CPPVAR_U77

#ELIF DEFINED(CPPVARIABLE_XLF)
!=============== ENVIRONMENT IBM AIX ===================================
!    IT ASSUMES THE ESSL LIBRARY TO BE LINKED
!    IT USES XLF SPECIFIC SUPPORT ROUTINES AND LANGUAGE EXTENSIONS
#DEFINE CPPVAR_ESSL 
#DEFINE CPPVAR_XLFsup
#DEFINE CPPVAR_XLFext
#DEFINE CPPVAR_GETRUSAGE
!#UNDEFINE EXLICITERF

#ELIF DEFINED(CPPVARIABLE_DEC)
!=============== ENVIRONMENT DEC ALPHA =================================
#DEFINE CPPVAR_FFTW
#DEFINE CPPVAR_ATLAS
#DEFINE CPPVAR_GETRUSAGE

#ELSE
!=============== FREE ENVIRONMENT ===== =================================
! MAKES NO ASSUMPTIONS ABOUT THE ENVIRONMENT
#DEFINE CPPVAR_FFTpack

#ENDIF 

#IFNDEF CPPVARIABLE_XLF 
!!.......................................................................
!MODULE ARGS_MODULE
!!***********************************************************************
!!**                                                                   **
!!** ARGS EMULATES THE XLF SYSTEM ROUTINE TO GRAB ARGUMENTS            **
!!** FROM THE COMMAND LINE.                                            **
!!**                                                                   **
!!** THE ARGUMENTS ARE READ FROM STANDARD INPUT                        **
!!** THE ROUTINES ARE NOT COMPLETED                                    **
!!**                                                                   **
!!***********************************************************************
!CHARACTER(256) :: ARGS
!LOGICAL(4)     :: TINI=.FALSE.
!END MODULE ARGS_MODULE   
!!     ..................................................................
!      INTEGER FUNCTION IARGC()
!      USE ARGS_MODULE
!!     ******************************************************************
!      IF(.NOT.TINI)READ(*,FMT='(A)')ARGS
!      TINI=.TRUE.
!      IARGC=1
!      RETURN
!      END
!!     ..................................................................
!      SUBROUTINE GETARG(I,TEXT)
!      USE ARGS_MODULE
!      IMPLICIT NONE
!      INTEGER(4)  ,INTENT(IN)  :: I
!      CHARACTER(*),INTENT(OUT) :: TEXT
!!     ******************************************************************
!      IF(.NOT.TINI)READ(*,FMT='(A)')ARGS
!      TINI=.TRUE.
!      TEXT=ARGS
!      RETURN
!      END
#ENDIF
!
!     ..................................................................
      SUBROUTINE LIB$TIMES(NTIME)
!     ******************************************************************
!     **    TIMES RETURNS IN UNITS OF CLOCKTICS                       **
!     **    1)  USER TIME OF THE PARENT PROCESS                       **
!     **    2)  SYSTEM TIME OF PARENT PROCESS                         **
!     **    3)  USER TIME OF CHILD PROCESSES                          **
!     **    4)  SYSTEM TIME OF CHILD PROCESSES                        **
!     **                                                              **
!     **  A CLOCKTICK IS CONSIDERED 0.01 SEC                          **
!     **                                                              **
!     **  SUBROUTINES USED: TIMES (STANDARD C LIBRARY)                **
!     ******************************************************************
      INTEGER(4),INTENT(OUT) :: NTIME(4)
      REAL(4)   ,EXTERNAL    :: ETIME     ! TOTAL ELAPSED TIME
      REAL(4)                :: TARRAY(2) ! ELAPSED USER/SYSTEM TIME
      REAL(4)                :: TOTAL
!     ******************************************************************
#IF DEFINED(CPPVAR_XLFsup)
      CALL TIMES(NTIME)     !XLF
#ELIF DEFINED(CPPVAR_U77)
      TOTAL=ETIME(TARRAY)   ! FROM ABSOFT SUPPORT LIBRARY libU77.a
      NTIME=NINT(100*TOTAL)
#ELSE 
      NTIME=0
#ENDIF
      RETURN
      END
!     .................................................................
      SUBROUTINE LIB$GETUSAGE(ID,VALUE)
!     *****************************************************************
!     **                                                             **
!     **  PROVIDES INFORMATION  ON THE USAGE OF SYSTEM RESOURCES     **
!     **                                                             **
!     **  uses standard c-library routine getrusage                  **
!     **                                                             **
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      IMPLICIT NONE
      TYPE USG_TYPE 
        SEQUENCE
        INTEGER(4) :: UTIME(2)   ! USER TIME    (SECONDS,MICROSECONDS)
        INTEGER(4) :: STIME(2)   ! SYSTEM TIME  (SECONDS,MICROSECONDS)
        INTEGER(4) :: MAXRSS     ! MAX SIZE [KBYTE] 
        INTEGER(4) :: IXRSS      ! INTEGRAL SHARED MEMORY SIZE [KBYTE*SEC]         
        INTEGER(4) :: IDRSS      ! INTEGRAL UNSHARED DATA [KBYTE*SEC]    
        INTEGER(4) :: ISRSS      ! INTEGRAL UNSHARED STACK [KBYTE*SEC]    
        INTEGER(4) :: MINFLT     ! #(PAGE RECLAIMS)
        INTEGER(4) :: MAJFLT     ! #(PAGE FAULTS)
        INTEGER(4) :: NSWAP      ! #(SWAPPING PROCESS OUT OF MAIN)            
        INTEGER(4) :: INBLOCK    ! #(FILE INPUTS)                             
        INTEGER(4) :: OUBLOCK    ! #(FILE OUTPUTS)                            
        INTEGER(4) :: MSGSND     ! #(IPC MESSAGES SEND)                       
        INTEGER(4) :: MSGRCV     ! #(IPC MESSAGES RECEIVED)                   
        INTEGER(4) :: NSIGNALS   ! #(SIGNALS SENT)                            
        INTEGER(4) :: NVCSW      ! #(VOLUNTARY CONTEXT SWITCHES)
        INTEGER(4) :: NIVCSW     ! #(INVOLUNTARY CONTEXT SWITCHES)
      END TYPE USG_TYPE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VALUE
      TYPE(USG_TYPE)           :: USG
      REAL(8)                  :: KBYTE=2.D0**10
      INTEGER(4)               :: GETRUSAGE
      INTEGER(4)               :: RC
      REAL(8)                  :: CPUTIME
      EXTERNAL GETRUSAGE
!     ******************************************************************
      usg%UTIME=(/0,0/)
      usg%STIME=(/0,0/)
      usg%MAXRSS=0     
      usg%IXRSS=0      
      usg%IDRSS=0      
      usg%ISRSS=0      
      usg%MINFLT=0     
      usg%MAJFLT=0     
      usg%NSWAP=0      
      usg%INBLOCK=0    
      usg%OUBLOCK=0    
      usg%MSGSND=0     
      usg%MSGRCV=0     
      usg%NSIGNALS=0   
      usg%NVCSW=0      
      usg%NIVCSW=0     
!     ==================================================================
!     ==================================================================
#IF DEFined(CPPVAR_GETRUSAGE)
      RC=GETRUSAGE(%VAL(0),USG)    ! C-ROUTINE
      IF(RC.NE.0) THEN
        CALL ERROR$MSG('ERROR CALLING MYGETRUSAGE')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
#ENDIF
!
      CPUTIME=(REAL(USG%UTIME(1)+USG%STIME(1),KIND=8) &
     &        +REAL(USG%UTIME(2)+USG%STIME(2),KIND=8))*1.D-6
      IF(ID.EQ.'MAXMEM') THEN
        VALUE=REAL(USG%MAXRSS,KIND=8)*KBYTE
      ELSE IF(ID.EQ.'USRTIME') THEN
        VALUE=(REAL(USG%UTIME(1),KIND=8)+REAL(USG%UTIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'SYSTIME') THEN
        VALUE=(REAL(USG%STIME(1),KIND=8)+REAL(USG%STIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'CPUTIME') THEN
        VALUE=CPUTIME
      ELSE IF(ID.EQ.'PAGEFAULTRATE') THEN
        VALUE=REAL(USG%MAJFLT,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWAPRATE') THEN
        VALUE=REAL(USG%NSWAP,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWITCHRATE') THEN
        VALUE=REAL(USG%NVCSW+USG%NIVCSW,KIND=8)/CPUTIME
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
      RETURN
      END
#IF DEFined(Cppvar_xlfext)
!     ..................................................................
      SUBROUTINE LIB$ERFR8(X,Y)
!     ******************************************************************
!     **  ERROR FUNCTION                                              **
!     **    Y=2/SQRT(PI)  INT_0^X DZ EXP(-Z**2)                       **
!     **    Y=(INFINITY)=1                                            **
!     ******************************************************************
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
!     ******************************************************************
      Y=dERF(X) 
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$ERFCR8(X,Y)
!     ****************************************************************** 
!     **  COMPLEMENTARY ERROR FUNCTION                                **
!     **  Y=1-ERF(X)                                                  **
!     ******************************************************************
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
      REAL(8), EXTERNAL :: dERFC
!     ******************************************************************
      Y=dERFC(X)
      RETURN
      END
#ELSE
!
!     ..................................................................
      subroutine lib$erfr8(x,y)
!     ******************************************************************
!     **  Copyright(C) 1996 Takuya OOURA                              **
!     **  (email: ooura@mmm.t.u-tokyo.ac.jp).                         **
!     **  You may use, copy, modify this code for any purpose and     **
!     **  without fee. You may distribute this ORIGINAL package.        **
!     ******************************************************************
      implicit none
      real(8),intent(in) :: x
      real(8),intent(out):: y
      real(8)            :: w
      real(8)            :: t
      integer(4)         :: k,i
      real(8)            :: a(0:64)
      real(8)            :: b(0:64)
!     ******************************************************************
      data (a(i), i = 0, 12) / &
     &    0.00000000005958930743d0, -0.00000000113739022964d0, &
     &    0.00000001466005199839d0, -0.00000016350354461960d0, &
     &    0.00000164610044809620d0, -0.00001492559551950604d0, &
     &    0.00012055331122299265d0, -0.00085483269811296660d0, &
     &    0.00522397762482322257d0, -0.02686617064507733420d0, &
     &    0.11283791670954881569d0, -0.37612638903183748117d0, &
     &    1.12837916709551257377d0 / 
      data (a(i), i = 13, 25) /                                &
     &    0.00000000002372510631d0, -0.00000000045493253732d0, &
     &    0.00000000590362766598d0, -0.00000006642090827576d0, &
     &    0.00000067595634268133d0, -0.00000621188515924000d0, &
     &    0.00005103883009709690d0, -0.00037015410692956173d0, &
     &    0.00233307631218880978d0, -0.01254988477182192210d0, &
     &    0.05657061146827041994d0, -0.21379664776456006580d0, &
     &    0.84270079294971486929d0 / 
      data (a(i), i = 26, 38) /                                &
     &    0.00000000000949905026d0, -0.00000000018310229805d0, &
     &    0.00000000239463074000d0, -0.00000002721444369609d0, &
     &    0.00000028045522331686d0, -0.00000261830022482897d0, &
     &    0.00002195455056768781d0, -0.00016358986921372656d0, &
     &    0.00107052153564110318d0, -0.00608284718113590151d0, &
     &    0.02986978465246258244d0, -0.13055593046562267625d0, &
     &    0.67493323603965504676d0 / 
      data (a(i), i = 39, 51) /                                &
     &    0.00000000000382722073d0, -0.00000000007421598602d0, &
     &    0.00000000097930574080d0, -0.00000001126008898854d0, &
     &    0.00000011775134830784d0, -0.00000111992758382650d0, &
     &    0.00000962023443095201d0, -0.00007404402135070773d0, &
     &    0.00050689993654144881d0, -0.00307553051439272889d0, &
     &    0.01668977892553165586d0, -0.08548534594781312114d0, &
     &    0.56909076642393639985d0 / 
      data (a(i), i = 52, 64) /                                &
     &    0.00000000000155296588d0, -0.00000000003032205868d0, &
     &    0.00000000040424830707d0, -0.00000000471135111493d0, &
     &    0.00000005011915876293d0, -0.00000048722516178974d0, &
     &    0.00000430683284629395d0, -0.00003445026145385764d0, &
     &    0.00024879276133931664d0, -0.00162940941748079288d0, &
     &    0.00988786373932350462d0, -0.05962426839442303805d0, &
     &    0.49766113250947636708d0 / 
      data (b(i), i = 0, 12) /                                 &
     &   -0.00000000029734388465d0,  0.00000000269776334046d0, &
     &   -0.00000000640788827665d0, -0.00000001667820132100d0, & 
     &   -0.00000021854388148686d0,  0.00000266246030457984d0, &
     &    0.00001612722157047886d0, -0.00025616361025506629d0, &
     &    0.00015380842432375365d0,  0.00815533022524927908d0, &
     &   -0.01402283663896319337d0, -0.19746892495383021487d0, & 
     &    0.71511720328842845913d0 / 
      data (b(i), i = 13, 25) /                                 &
     &   -0.00000000001951073787d0, -0.00000000032302692214d0,  &
     &    0.00000000522461866919d0,  0.00000000342940918551d0,  &
     &   -0.00000035772874310272d0,  0.00000019999935792654d0,  &
     &    0.00002687044575042908d0, -0.00011843240273775776d0,  &
     &   -0.00080991728956032271d0,  0.00661062970502241174d0,  &
     &    0.00909530922354827295d0, -0.20160072778491013140d0,  &
     &    0.51169696718727644908d0 / 
      data (b(i), i = 26, 38) /                                 &
     &    0.00000000003147682272d0, -0.00000000048465972408d0,  &
     &    0.00000000063675740242d0,  0.00000003377623323271d0,  &
     &   -0.00000015451139637086d0, -0.00000203340624738438d0,  &
     &    0.00001947204525295057d0,  0.00002854147231653228d0,  &
     &   -0.00101565063152200272d0,  0.00271187003520095655d0,  &
     &    0.02328095035422810727d0, -0.16725021123116877197d0,  &
     &    0.32490054966649436974d0 / 
      data (b(i), i = 39, 51) /                                &
     &    0.00000000002319363370d0, -0.00000000006303206648d0, & 
     &   -0.00000000264888267434d0,  0.00000002050708040581d0, & 
     &    0.00000011371857327578d0, -0.00000211211337219663d0, & 
     &    0.00000368797328322935d0,  0.00009823686253424796d0, & 
     &   -0.00065860243990455368d0, -0.00075285814895230877d0, & 
     &    0.02585434424202960464d0, -0.11637092784486193258d0, & 
     &    0.18267336775296612024d0 / 
      data (b(i), i = 52, 64) /                                &
     &   -0.00000000000367789363d0,  0.00000000020876046746d0, & 
     &   -0.00000000193319027226d0, -0.00000000435953392472d0, & 
     &    0.00000018006992266137d0, -0.00000078441223763969d0, & 
     &   -0.00000675407647949153d0,  0.00008428418334440096d0, & 
     &   -0.00017604388937031815d0, -0.00239729611435071610d0, & 
     &    0.02064129023876022970d0, -0.06905562880005864105d0, & 
     &    0.09084526782065478489d0 / 
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - real(k,8)
          k = k * 13
          y = (((((((((((( a(k    )*t + a(k+ 1)) * t + &
     &        a(k+ 2))*t + a(k+ 3))*t + a(k+ 4)) * t + &
     &        a(k+ 5))*t + a(k+ 6))*t + a(k+ 7)) * t + &
     &        a(k+ 8))*t + a(k+ 9))*t + a(k+10)) * t + &
     &        a(k+11))*t + a(k+12))*w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - real(k,8)
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + &
     &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
     &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
     &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
     &        b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if(x.lt.0) y=-y
      return
      end
!     .......................................................................
      subroutine lib$erfcr8(x,y)
!     **  Copyright(C) 1996 Takuya OOURA                              **
!     **  (email: ooura@mmm.t.u-tokyo.ac.jp).                         **
!     **  You may use, copy, modify this code for any purpose and     **
!     **  without fee. You may distribute this ORIGINAL package.        **
      implicit none
      real(8),intent(in) :: x
      real(8),intent(out):: y
      real(8)     ,parameter :: pa  = 3.97886080735226000d+00
      real(8)     ,parameter :: p0  = 2.75374741597376782d-01
      real(8)     ,parameter :: p1  = 4.90165080585318424d-01 
      real(8)     ,parameter :: p2  = 7.74368199119538609d-01 
      real(8)     ,parameter :: p3  = 1.07925515155856677d+00 
      real(8)     ,parameter :: p4  = 1.31314653831023098d+00 
      real(8)     ,parameter :: p5  = 1.37040217682338167d+00 
      real(8)     ,parameter :: p6  = 1.18902982909273333d+00 
      real(8)     ,parameter :: p7  = 8.05276408752910567d-01 
      real(8)     ,parameter :: p8  = 3.57524274449531043d-01 
      real(8)     ,parameter :: p9  = 1.66207924969367356d-02 
      real(8)     ,parameter :: p10 =-1.19463959964325415d-01 
      real(8)     ,parameter :: p11 =-8.38864557023001992d-02
      real(8)     ,parameter :: p12 = 2.49367200053503304d-03 
      real(8)     ,parameter :: p13 = 3.90976845588484035d-02 
      real(8)     ,parameter :: p14 = 1.61315329733252248d-02 
      real(8)     ,parameter :: p15 =-1.33823644533460069d-02 
      real(8)     ,parameter :: p16 =-1.27223813782122755d-02 
      real(8)     ,parameter :: p17 = 3.83335126264887303d-03 
      real(8)     ,parameter :: p18 = 7.73672528313526668d-03 
      real(8)     ,parameter :: p19 =-8.70779635317295828d-04 
      real(8)     ,parameter :: p20 =-3.96385097360513500d-03 
      real(8)     ,parameter :: p21 = 1.19314022838340944d-04 
      real(8)     ,parameter :: p22 = 1.27109764952614092d-03
      real(8)                :: t,u
!     **********************************************************
      t = pa / (pa + abs(x))
      u = t - 0.5d0
      y = (((((((((p22 * u + p21) * u + p20) * u + &
     &    p19) * u + p18) * u + p17) * u + p16) * u + &
     &    p15) * u + p14) * u + p13) * u + p12 
      y = ((((((((((((y * u + p11) * u + p10) * u + &
     &    p9) * u + p8) * u + p7) * u + p6) * u + p5) * u + &
     &    p4) * u + p3) * u + p2) * u + p1) * u + p0) * t * &
     &    exp(-x * x)
      if (x .lt. 0) y = 2 - y
      return
      end
#ENDIF
!
!
!.......................................................................
MODULE RANDOM_MODULE
  INTEGER(4),SAVE        :: SEED=1
  INTEGER(4),PARAMETER   :: RANFAC1=22222221
  INTEGER(4),PARAMETER   :: RANFAC2=2**24
! CHOICE EXPLICIT CPPVARIABLE_XLF RNG
! INTEGER(8),SAVE        :: SEED=1_8
! INTEGER(8),PARAMETER   :: RANFAC1=44485709377909_8
! INTEGER(8),PARAMETER   :: RANFAC2=2_8**48
END MODULE RANDOM_MODULE
!
!     ..................................................................
      SUBROUTINE LIB$RANDOM(X)
!     ****************************************************************** 
!     **  returns a random number                                     **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: X
!     ******************************************************************
      SEED=MODULO(RANFAC1*SEED,RANFAC2)
      X=REAL(SEED,8)/REAL(RANFAC2,8)
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     CALL RANDOM_NUMBER(Y)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$RANDOMSEED
!     ****************************************************************** 
!     **  COMPLEMENTARY ERROR FUNCTION                                **
!     **  Y=1-ERF(X)                                                  **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
!     ******************************************************************
      SEED=1
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     SEED=1_8
!     == CHOICE ORIGINAL CPPVARIABLE_XLF RNG
!     CALL RANDOM_SEED(GENERATOR=2)
!     CALL RANDOM_SEED
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$INVERTR8(N,A,AINV)
!     ******************************************************************
!     **                                                              **
!     **  INVERTS THE REAL, SQUARE MATRIX A                           **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEICD                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: A(N,N)
      REAL(8)   ,INTENT(OUT):: AINV(N,N)
      INTEGER(4)            :: NAUX
      REAL(8)               :: AUX(100*N)
      INTEGER(4)            :: I,J
#IF DEFined(CPPVAR_ESSL)
      REAL(8)               :: RCOND
      REAL(8)               :: DET(2)
#ELSE 
      INTEGER(4)            :: IPIV(N)
      INTEGER(4)            :: INFO
#ENDIF
!     ******************************************************************
      NAUX=100*N
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        ENDDO
      ENDDO
#IF DEFined(CPPVAR_ESSL)
      CALL DGEICD(AINV,N,N,0,RCOND,DET,AUX,NAUX) !ESSL
#ELSE 
      CALL DGETRF(N,N,AINV,N,IPIV,INFO) !lapack
      CALL DGETRI(N,AINV,N,IPIV,AUX,NAUX,INFO) !lapack
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$DIAGR8(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **         U(K,I)*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                 **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL DSPEV  MATRIX DIAGONALIZATION P727                   **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
#if defined(cppvar_essl)
      INTERFACE
        SUBROUTINE EINFO(ICODE,INF1,INF2) !ESSL ERROR HANDLING ROUTINE
        INTEGER                       :: ICODE
        INTEGER ,INTENT(OUT),OPTIONAL :: INF1
        INTEGER ,INTENT(OUT),OPTIONAL :: INF2
        END SUBROUTINE EINFO
      END INTERFACE
#endif
      LOGICAL(4) ,PARAMETER :: TESSLERR=.FALSE.
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: WORK1((N*(N+1))/2)
#if defined(CPPVAR_ESSL)
      REAL(8)               :: WORK2(2*N)
#ELSE
      REAL(8)               :: WORK2(3*N)
#ENDIF
      INTEGER(4)            :: K,I,J
      CHARACTER(8)          :: SAV2101
      INTEGER(4)            :: I1,I2
      integer(4)            :: INFO
!     ******************************************************************
!
!     ==================================================================
!     ====  STORE IN LOWER PACKED STORAGE MODE FOR DIAGONALIZATION    ==
!     ==================================================================
      K=0
      DO J=1,N
        DO I=J,N
          K=K+1
          WORK1(K)=0.5D0*(H(I,J)+H(J,I))
        ENDDO
      ENDDO
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
#IF DEFined(CPPVAR_ESSL)
      IF(TESSLERR) THEN
        CALL EINFO(0)
        CALL ERRSAV(2101,SAV2101)
        CALL ERRSET(2101,255,0,0,0,2101)
!       CALL DSPEV(1,WORK1,E,U,N,N,WORK2,2*N,400) !->ESSL
        CALL DSPEV(1,WORK1,E,U,N,N,WORK2,2*N) !->ESSL
      ELSE
        CALL DSPEV(1,WORK1,E,U,N,N,WORK2,2*N) !->ESSL
      END IF
#ELSE
      CALL DSPEV('V','L',N,WORK1,E,U,N,WORK2,INFO) !->LAPACK
      IF(INFO.NE.0) THEN
        CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
        CALL ERROR$STOP('DIAG')
      END IF
#ENDIF
      RETURN
!
!     ==================================================================
!     == ESSL ERROR HANDLING                                          ==
!     ==================================================================
#IF DEFined(CPPVAR_ESSL)
 400  CONTINUE
      CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
      IF(TESSLERR) THEN
        CALL EINFO(2101,I1,I2)
        CALL ERROR$I4VAL('EIGENVALUE FAILED TO CONVERGE',I1)
        CALL ERROR$I4VAL('NUMBER OF ITERATIONS',I2)
      END IF
      DO I=1,N
        PRINT*,'H ',I,H(1:3,I)
      ENDDO
      CALL ERROR$STOP('LIB$DIAGR8')
      STOP
#ENDIF
      END
!
!     ..................................................................
      SUBROUTINE LIB$DIAGC8(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE HERMITEAN, SQUARE MATRIX H                 **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **      CONJG(U(K,I))*H(K,L)*U(L,J)=DELTA(I,J)*E(I)             **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: ZHPEV :  MATRIX DIAGONALIZATION  P727               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTERFACE
        SUBROUTINE EINFO(ICODE,INF1,INF2) !ESSL ERROR HANDLING ROUTINE
        INTEGER                       :: ICODE
        INTEGER ,INTENT(OUT),OPTIONAL :: INF1
        INTEGER ,INTENT(OUT),OPTIONAL :: INF2
        END SUBROUTINE EINFO
      END INTERFACE
      LOGICAL(4) ,PARAMETER :: TESSLERR=.FALSE.
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: U(N,N)
      COMPLEX(8)            :: WORK1((N*(N+1))/2)
#if defined(CPPVAR_ESSL)
      REAL(8)               :: RWORK(4*N)
#else
      COMPLEX(8)            :: CWORK(2*N)
      REAL(8)               :: RWORK(3*N)
#endif
      INTEGER(4)            :: K,I,J
      CHARACTER(8)          :: SAV2101
      INTEGER(4)            :: I1,I2
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8)            :: CSVAR
      LOGICAL(4)            :: TCHK
      COMPLEX(8),ALLOCATABLE:: WORK3(:,:)
      INTEGER(4)            :: INFO
!     ******************************************************************

!     ==================================================================
!     ====  ESSL ERROR HANDLING                                       ==
!     ==================================================================
!
!     ==================================================================
!     ====  STORE IN LOWER PACKED STORAGE MODE FOR DIAGONALIZATION    ==
!     ==================================================================
      K=0
      DO J=1,N
        DO I=J,N
          K=K+1
          WORK1(K)=0.5D0*(H(I,J)+CONJG(H(J,I)))
        ENDDO
      ENDDO
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
#IF DEFined(CPPVAR_ESSL)
      IF(TESSLERR) THEN
        CALL EINFO(0)
        CALL ERRSAV(2101,SAV2101)
        CALL ERRSET(2101,255,0,0,0,2101)
      END IF
      CALL ZHPEV(1,WORK1,E,U,N,N,RWORK,4*N)  !essl
#ELSE
      CALL ZHPEV('V','L',N,WORK1,E,U,N,CWORK,RWORK,INFO) !lapack
      IF(INFO.NE.0) THEN
        CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
        CALL ERROR$STOP('DIAG')
      END IF
#ENDIF
!
!     ==================================================================
!     ====  OPTIONAL TEST                                             ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(WORK3(N,N)) 
        TCHK=.TRUE.
        DO I=1,N
          DO J=1,N
            CSVAR=(0.D0,0.D0)
            DO K=1,N
              CSVAR=CSVAR+U(I,K)*E(K)*CONJG(U(J,K))
            ENDDO
            IF(ABS(CSVAR-H(I,J)).GT.1.D-8) THEN
              WRITE(*,FMT='(2I5,5F10.5)')I,J,CSVAR,H(I,J),ABS(CSVAR-H(I,J))
              TCHK=.FALSE.
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(WORK3)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$STOP('LIB$DIAGC8')
        END IF
      END IF
      RETURN
!
!     ==================================================================
!     == ESSL ERROR HANDLING                                          ==
!     ==================================================================
#IF DEFined(CPPVAR_ESSL)
 400  CONTINUE
      CALL ERROR$MSG('DIAGONALIZATION NOT CONVERGED')
      IF(TESSLERR) THEN
        CALL EINFO(2101,I1,I2)
        CALL ERROR$I4VAL('EIGENVALUE FAILED TO CONVERGE',I1)
        CALL ERROR$I4VAL('NUMBER OF ITERATIONS',I2)
      END IF
      DO I=1,N
        PRINT*,'H ',I,H(1:3,I)
      ENDDO
      CALL ERROR$STOP('LIB$DIAGC8')
      STOP
#ENDIF
      END
!
!     ..................................................................
      SUBROUTINE LIB$FFTC8(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     **  use fftw as standard                                        **
!     **  use fftessl if essl is installed                            **
!     **  use fftpack as backup if c-routines cannot be linked or     **
!     **      fftw is not available                                   **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
!     *******************************************************************
#IF DEFined(CPPVAR_ESSL)
      call LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
#ELIF defined(CPPVAR_FFTW)
      CALL LIB_FFTW(DIR,LEN,NFFT,X,Y)
#elif defined(cppvar_fftpack)
      CALL LIB_FFTpack(DIR,LEN,NFFT,X,Y)
#else
      call error$msg('no fft package selected during compilation')
      call error$stop('LIB$FFTC8')
#ENDIF
      RETURN
      END
!
#IF DEFined(CPPVAR_FFTW)
!     ..................................................................
      SUBROUTINE LIB_FFTW(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INCLUDE 'FFTW_F77.I'
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,I
      COMPLEX(8)              :: XDUMMY(LEN,NFFT)
      INTEGER(4),SAVE         :: NP=0
      INTEGER(4),PARAMETER    :: NPX=10 ! #(DIFFERENT FFT PLANS)
      INTEGER(4),SAVE         :: PLANS(NPX,2),PLAN=-1
      LOGICAL                 :: DEF
!     ******************************************************************
      XDUMMY=X
!
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('1D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST
        DEF=.FALSE.
        DO I=1,NP
          IF((LEN*ISIGN).EQ.PLANS(I,1)) THEN
            DEF=.TRUE.
            PLAN=PLANS(I,2)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==========================
        IF(.NOT.DEF) THEN
          WRITE(*,*) 'FFTW CREATE PLAN FOR: ', ISIGN,LEN,NP
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
          CALL FFTW_F77_CREATE_PLAN(PLANS(NP,2),LEN,ISIGN,FFTW_MEASURE)
          PLANS(NP,1)=ISIGN*LEN
          PLAN=PLANS(NP,2)
        END IF
        LENSAVE=LEN
        DIRSAVE=DIR
        SCALE=1.D0/REAL(LEN,KIND=8)
      END IF
!
!     ==================================================================
!     ==  NOW PERFORM FFT                                             ==
!     ==================================================================
      CALL FFTW_F77(PLAN,NFFT,XDUMMY(1,1),1,LEN,Y(1,1),1,LEN)
      IF (DIR.EQ.'RTOG') THEN
        DO IFFT=1,NFFT
          DO I=1,LEN
            Y(I,IFFT)=Y(I,IFFT)*SCALE
          ENDDO
        ENDDO
      END IF
      RETURN
      END
#ENDIF
!
!DCFT:  1-D FFT P765
#IF DEFined(CPPVAR_ESSL)
!     ..................................................................
      SUBROUTINE LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX1=20000
      INTEGER(4)  ,PARAMETER  :: NAUX2=20000
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      REAL(8)     ,SAVE       :: AUX1(NAUX1)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,i
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE.OR.NFFT.NE.NFFTSAVE) THEN
        IF(LEN.GT.8192) THEN
          CALL ERROR$MSG('FFT TOO LONG')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('FFT')
        END IF
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('FFT')
        END IF
        CALL DCFT(1,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
        DIRSAVE=DIR
        LENSAVE=LEN
        NFFTSAVE=NFFT
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      CALL DCFT(0,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
      RETURN
      END
#ENDIF
!
!DCFT:  1-D FFT P765
#if defined(cppvar_fftpack)
!     ..................................................................
      SUBROUTINE LIB_FFTpack(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX2=200
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      INTEGER(4)  ,SAVE       :: IAUX(30)
      real(8)                 :: sequence(2*len)
      CHARACTER(4),SAVE       :: DIRSAVE=''
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,i
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('lib_fftpack')
        END IF
        DIRSAVE=DIR
        LENSAVE=LEN
        IF(LEN.EQ.1) RETURN
        IF(NAUX2.LT.2*LEN) THEN
          CALL ERROR$MSG('AUXILIARY ARRAY TOO SMALL: INCREASE NAUX2')
          CALL ERROR$I4VAL('NAUX2',NAUX2)
          CALL ERROR$I4VAL('2*LEN',2*LEN)
          CALL ERROR$STOP('LIB_FFTPACK')
        END IF
        CALL CFFTI(LEN,AUX2,IAUX)
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      IF(LEN.EQ.1) RETURN
      IF(ISIGN.EQ.1) THEN
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT))
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTF(LEN,sequence,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I))*SCALE
          ENDDO
        ENDDO
      ELSE
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT))
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTB(LEN,sequence,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I))*SCALE
          ENDDO
        ENDDO
      END IF
      RETURN
      END
#endif
!
!     ...................................................FESSL..........
      SUBROUTINE LIB$FFTADJUSTGRD(NR)
!     ******************************************************************
!     **  THIS ROUTINE RETURNS THE ALLOWED FOURIER TRANSFORM LENGTH   ** 
!     **  THAT IS EQUAL OR LARGER THAN THE LENGTH SUPPLIED, BUT       ** 
!     **  BUT OTHERWISE AS SMALL AS POSSIBLE.                         **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT):: NR
      INTEGER(4),PARAMETER  :: MAXI=300
      LOGICAL(4),SAVE       :: TINIT=.TRUE.
      INTEGER(4),SAVE       :: COUNT
      INTEGER(4),SAVE       :: IFR(MAXI)
      INTEGER(4)            :: H,I,J,K,M
      INTEGER(4)            :: ISVAR
      REAL(8)               :: SVAR
!     ******************************************************************
      IF (TINIT) THEN
        TINIT=.FALSE.
#IF DEFined(CPPVAR_ESSL)
!       == allowed length for the ESSL FFT =============================
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,2
            DO J=0,1
              DO K=0,1
                DO M=0,1
                  IF(COUNT.GE.MAXI) EXIT OUTER
                  SVAR = 2.D0**H * 3.D0**I * 5.D0**J * 7.D0**K *11.D0**M
                  IF(SVAR.GT.37748736.D0) CYCLE
                  COUNT=COUNT+1
                  IFR(COUNT)=NINT(SVAR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO OUTER
#ELIF defined(CPPVAR_FFTW)
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,2
            DO J=0,1
              DO K=0,1
                DO M=0,1
                  IF(COUNT.GE.MAXI) EXIT OUTER
                  SVAR = 2.D0**H * 3.D0**I * 5.D0**J * 7.D0**K *11.D0**M
                  IF(SVAR.GT.37748736.D0) CYCLE
                  COUNT=COUNT+1
                  IFR(COUNT)=NINT(SVAR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO OUTER
#ELSE
!       == the general length has been removed because passf and passb
!       == do not work
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,4
            DO J=0,2
              IF(COUNT.GE.MAXI) EXIT OUTER
              SVAR = 2.D0**H * 3.D0**I * 5.D0**J 
              IF(SVAR.GT.37748736.D0) CYCLE
              COUNT=COUNT+1
              IFR(COUNT)=NINT(SVAR)
            ENDDO
          ENDDO
        ENDDO OUTER
#endif
        DO I=1,COUNT
          DO J=I+1,COUNT
            IF(IFR(I).GT.IFR(J)) THEN
              ISVAR=IFR(I)
              IFR(I)=IFR(J)
              IFR(J)=ISVAR
            END IF
          ENDDO
        ENDDO
      ENDIF
!
      DO I=2,COUNT
        IF(IFR(I).GE.NR) THEN
          NR=IFR(I)
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('REQUESTED GRIDPOINT OUT OF RANGE')
      CALL ERROR$I4VAL('NR',NR)
      CALL ERROR$STOP('PLANEWAVE$ADJUSTFFTGRD')
      STOP
      END
!
!     ..................................................................
      SUBROUTINE LIB$3DFFTC8(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFTW ROUTINES                                 **
!     **                                        Clemens Foerst, 2001  **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER(4)              :: plan
      INTEGER(4)  ,SAVE       :: np=0
      INTEGER(4),parameter    :: npx=10
      INTEGER(4)  ,SAVE       :: plans(npx,4)
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)  ,SAVE       :: DIMSAVE(3)=0
      CHARACTER(4),SAVE       :: DIRSAVE=''
      LOGICAL                 :: DEF
      INTEGER(4)              :: I
      INTEGER(4)  ,SAVE       :: ISIGN
      INCLUDE 'FFTW_F77.I'
!     ******************************************************************
      DIM(1)=N1
      DIM(2)=N2
      DIM(3)=N3
      IF(DIM(1).NE.DImSAVE(1).OR.DIm(2).NE.DImSAVE(2).OR. &
     &  DIm(3).NE.DImSAVE(3).OR.DIr.NE.DIrSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('3D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST
        DEF=.FALSE.
        DO I=1,NP
          IF((DIM(1)*ISIGN).EQ.PLANS(I,1).AND.(DIM(2)*ISIGN).EQ.PLANS(I,2)&
     &                     .AND.(DIM(3)*ISIGN).EQ.PLANS(I,3)) THEN
            DEF=.TRUE.
            PLAN=PLANS(I,4)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==========================
        IF(.NOT.DEF) THEN
          WRITE(*,*) '3D-FFTW CREATE PLAN FOR: ', ISIGN,DIM
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
          call fftwnd_f77_create_plan(PLAN,3,DIM,isign,FFTW_ESTIMATE)
          PLANS(NP,1:3)=ISIGN*DIM
          PLANS(NP,4)=PLAN
        END IF
        DIMSAVE=DIM
        DIRSAVE=DIR
        IF(DIR.EQ.'RTOG') SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
      END IF
!
!     ==================================================================
!     ==  NOW PERFORM FFT                                             ==
!     ==================================================================
      call fftwnd_f77_one(plan,X,Y)
      IF (DIR.EQ.'RTOG') Y=Y*SCALE
      RETURN
      END

!
!DGEMUL: MATRIX MULTIPLICATION P441
!     ..................................................................
      SUBROUTINE LIB$MATMULR8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(IN) :: B(M,L)
      REAL(8)   ,INTENT(OUT):: C(N,L)
      CHARACTER(8),PARAMETER:: LIB='ESSL'
      real(8)     ,PARAMETER:: ONE=1.D0
      real(8)     ,PARAMETER:: ZERO=0.D0
!     ******************************************************************
#IF DEFined(CPPVAR_ESSL)
      CALL DGEMUL(A,N,'N',B,M,'N',C,N,N,M,L)
#ELif defined(CPPVAR_ATLAS)
      CALL DGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#ELSE
      C=MATMUL(A,B)
#ENDIF
      RETURN
      END
!
!ZGEMUL: MATRIX MULTIPLICATION P441
!     . ................................................................
      SUBROUTINE LIB$MATMULC8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(IN) :: B(M,L)
      COMPLEX(8),INTENT(OUT):: C(N,L)
      CHARACTER(8),PARAMETER:: LIB='ESSL'
      COMPLEX(8)  ,PARAMETER:: ONE=(1.D0,0.D0)
      COMPLEX(8)  ,PARAMETER:: ZERO=(0.D0,0.D0)
!     ******************************************************************
#IF DEFined(CPPVAR_ESSL)
      CALL ZGEMUL(A,N,'N',B,M,'N',C,N,N,M,L)
#ELif defined(CPPVAR_ATLAS)
      CALL ZGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
#else
      C=MATMUL(A,B)
#ENDIF
      RETURN
      END
!
!ZGEMM(TRANSA,TRANSB,L,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC):
!C=BETA*C+ALPHA*A*B P456
!     ..................................................................
      SUBROUTINE LIB$ADDPRODUCTC8(TID,N,M,L,A,B,C)
!     ******************************************************************
!     **  ADD THE PRODUCT OF TWO MATRICES                             **
!     **  A=A+B*C                                                     **
!     **  OR A=A+A*C                                                  **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)   :: TID
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(IN)   :: M
      INTEGER(4),INTENT(IN)   :: L
      COMPLEX(8),INTENT(IN)   :: A(N,M)
      COMPLEX(8),INTENT(IN)   :: B(M,L)
      COMPLEX(8),INTENT(INOUT):: C(N,L)
      COMPLEX(8)              :: SUM
      COMPLEX(8),ALLOCATABLE  :: WORK(:,:)
      INTEGER(4)              :: I,J,K
!     ******************************************************************
#IF DEFined(CPPVAR_ESSL)
      IF(TID) THEN
        ALLOCATE(WORK(N,M))
        WORK=A
        DO I=1,N
          DO J=1,M
            SUM=(0.D0,0.D0)
            DO K=1,L
              SUM=SUM+WORK(I,K)*B(K,J)
            ENDDO
            C(I,J)=C(I,J)+SUM
          ENDDO
        ENDDO
        DEALLOCATE(WORK)
      ELSE
        CALL ZGEMM('N','N',N,L,M,(1.D0,0.D0),A,N,B,M,(1.D0,0.D0),C,N)
      END IF
#ELSE 
      IF(TID) THEN
        ALLOCATE(WORK(N,M))
        WORK=A
        DO I=1,N
          DO J=1,M
            SUM=(0.D0,0.D0)
            DO K=1,L
              SUM=SUM+WORK(I,K)*B(K,J)
            ENDDO
            C(I,J)=C(I,J)+SUM
          ENDDO
        ENDDO
        DEALLOCATE(WORK)
      ELSE
        DO I=1,N
          DO J=1,L
            SUM=(0.D0,0.D0)
            DO K=1,M
              SUM=SUM+A(I,K)*B(K,J)
            ENDDO
            C(I,J)=C(I,J)+SUM
          ENDDO
        ENDDO
      END IF
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMR8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: PSI1(LEN1,N)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN2,N)
      REAL(8)   ,INTENT(OUT):: OPERATOR(LEN1,LEN2)
      INTEGER(4)            :: I,J,K
      REAL(8)               :: SUM
!     ******************************************************************
#IF DEFINED(XLF)
      CALL DGEMUL(PSI1,LEN1,'N',PSI2,LEN2,'T',OPERATOR,LEN1,LEN1,N,LEN2)
#ELSE 
      DO I=1,LEN1
        DO J=1,LEN2
          SUM=0.D0
          DO K=1,N
            SUM=SUM+PSI1(I,K)*PSI2(J,K)
          ENDDO
          OPERATOR(I,J)=SUM
        ENDDO
      ENDDO
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMC8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: PSI1(LEN1,N)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN2,N)
      COMPLEX(8),INTENT(OUT):: OPERATOR(LEN1,LEN2)
      INTEGER(4)            :: I,J,K
      COMPLEX(8)            :: SUM
!     ******************************************************************
#IF DEFined(CPPVAR_ESSL)
      CALL ZGEMUL(PSI1,LEN1,'N',PSI2,LEN2,'C',OPERATOR,LEN1,LEN1,N,LEN2)
#ELSE 
      DO I=1,LEN1
        DO J=1,LEN2
          SUM=0.D0
          DO K=1,N
            SUM=SUM+PSI1(I,K)*CONJG(PSI2(J,K))
          ENDDO
          OPERATOR(I,J)=SUM
        ENDDO
      ENDDO
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTR8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      REAL(8)   ,INTENT(IN) :: PSI1(LEN,N1)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN,N2)
      REAL(8)   ,INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J,K
      REAL(8)               :: SUM
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTR8')
      END IF
#IF DEFined(CPPVAR_ESSL)
      IF(TID) THEN
!       ==  OVERLAP(I,J) = 0.D0*OVERLAP+1.D0*SUM_K:PSI1(K,I)*PSI1(K,J) =
        CALL DSYRK('U','T',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=OVERLAP(I,J)
          ENDDO
        ENDDO
      ELSE
        CALL DGEMUL(PSI1,LEN,'T',PSI2,LEN,'N',OVERLAP,N1,N1,LEN,N2)
      END IF
#ELSE 
      IF(TID) THEN
        DO I=1,N1
          DO J=I,N1
            SUM=0.D0
            DO K=1,LEN
              SUM=SUM+PSI1(K,I)*PSI1(K,J)
            ENDDO
            OVERLAP(I,J)=SUM
            OVERLAP(J,I)=SUM
          ENDDO
        ENDDO
      ELSE
        DO I=1,N1
          DO J=1,N2
            SUM=0.D0
            DO K=1,LEN
              SUM=SUM+PSI1(K,I)*PSI2(K,J)
            ENDDO
            OVERLAP(I,J)=SUM
          ENDDO
        ENDDO
      END IF
#ENDIF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTC8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      COMPLEX(8),INTENT(IN) :: PSI1(LEN,N1)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN,N2)
      COMPLEX(8),INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J,K
      COMPLEX(8)            :: SUM
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTC8')
      END IF
#IF DEFined(CPPVAR_ESSL)
      IF(TID) THEN
        CALL ZHERK('U','C',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=CONJG(OVERLAP(I,J))
          ENDDO
        ENDDO
      ELSE
        CALL ZGEMUL(PSI1,LEN,'C',PSI2,LEN,'N',OVERLAP,N1,N1,LEN,N2)
      END IF
#ELSE 
      IF(TID) THEN
        DO I=1,N1
          DO J=I,N1
            SUM=(0.D0,0.D0)
            DO K=1,LEN
              SUM=SUM+CONJG(PSI1(K,I))*PSI1(K,J)
            ENDDO
            OVERLAP(I,J)=SUM
            OVERLAP(J,I)=CONJG(SUM)
          ENDDO
        ENDDO
      ELSE
        DO I=1,N1
          DO J=1,N2
            SUM=(0.D0,0.D0)
            DO K=1,LEN
              SUM=SUM+CONJG(PSI1(K,I))*PSI2(K,J)
            ENDDO
            OVERLAP(I,J)=SUM
          ENDDO
        ENDDO
      END IF
#ENDIF
      RETURN
      END
!
!ZAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!     ..................................................................
      SUBROUTINE LIB$VECTORADDC8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      COMPLEX(8),INTENT(INOUT) :: X(N)
      COMPLEX(8),INTENT(IN)    :: FAC
      COMPLEX(8),INTENT(IN)    :: Y(N)
      INTEGER(4)               :: I
      CHARACTER(8),PARAMETER   :: LIB='ESSL'
!     ******************************************************************
#if defined(Cppvar_essl)
        CALL ZAXPY(N,FAC,X,1,Y,1)
#ELSE
        DO I=1,N
          X(I)=X(I)+FAC*Y(I)
        ENDDO
#endif
      RETURN
      END
!
!DAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!DVES(N,X,INCX,Y,INCY,X,INCZ)  Z=X-Y       P313
!     ..................................................................
      SUBROUTINE LIB$VECTORADDR8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      REAL(8)   ,INTENT(INOUT) :: X(N)
      REAL(8)   ,INTENT(IN)    :: FAC
      REAL(8)   ,INTENT(IN)    :: Y(N)
      INTEGER(4)               :: I
      CHARACTER(8),PARAMETER   :: LIB='ESSL'
!     ******************************************************************
#IF DEFined(CPPVAR_ESSL)
      IF(FAC.EQ.1.D0) THEN
        CALL DAXPY(N,FAC,X,1,Y,1)
      ELSE IF(FAC.EQ.-1.D0) THEN
        CALL DVES(N,X,1,Y,1,X,1)
      ELSE
        CALL DAXPY(N,FAC,X,1,Y,1)
      END IF
#ELSE
      IF(FAC.EQ.1.D0) THEN
        DO I=1,N
          X(I)=X(I)+Y(I)
        ENDDO
      ELSE IF(FAC.EQ.-1.D0) THEN
        DO I=1,N
          X(I)=X(I)-Y(I)
        ENDDO
      ELSE
        DO I=1,N
          X(I)=X(I)+FAC*Y(I)
        ENDDO
      END IF
#ENDIF
      RETURN
      END
!
!DGEF(A,LDA,N,IPVT)            MATRIX FACTORIZATION P507
!DGES(A,LDA,N,IPVT,BX)         BX=A^{-1}*BX (USES IPVT FROM DGEF)
!DGESVF                        SINGULAR VALUE DECOMPOSITION P696
!DGESVS  LEAST SQUARES SOLUTION USING SINGULAR VALUE DECOMPOSITION P703
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVE(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      REAL(8)   ,ALLOCATABLE:: AFACT(:,:)
      REAL(8)   ,ALLOCATABLE:: BFACT(:,:)
      INTEGER(4),ALLOCATABLE:: IPVT(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)   ,ALLOCATABLE:: S(:)
      INTEGER(4)            :: NAUX
      INTEGER(4)            :: NM
      INTEGER(4)            :: IRANK
      REAL(8)               :: TAU=1.D-6
      INTEGER(4)            :: INFO
      INTEGER(4)            :: I
!     ******************************************************************
      IF(N.EQ.M) THEN  !MATRIX FACTORIZATION
        ALLOCATE(AFACT(N,N))
        ALLOCATE(IPVT(N)) 
        AFACT=A
#IF DEFined(CPPVAR_ESSL)
        CALL DGEF(AFACT,N,N,IPVT)
        X=B
        CALL DGES(AFACT,N,N,IPVT,X,0)
#ELSE 
        CALL DGEFA(AFACT,N,N,IPVT,INFO)
        X=B
        CALL DGESL(AFACT,N,N,IPVT,X,0)
#ENDIF
        DEALLOCATE(AFACT)
        DEALLOCATE(IPVT)
      ELSE   !SINGULAR VALUE DECOMPOSITION
        NM=MAX(1,MAX(N,M))
#IF DEFined(CPPVAR_ESSL)
        ALLOCATE(AFACT(NM,M))
        ALLOCATE(BFACT(N,NEQ))
        AFACT(1:N,:)=A
        AFACT(N+1:,:)=0.D0
        BFACT=B
        NAUX=2*M+MAX(NM,NEQ)
        ALLOCATE(S(M))    !SINGULAR VALUES
        ALLOCATE(AUX(NAUX))
        CALL DGESVF(2,AFACT,N,BFACT,N,NEQ,S,N,M,AUX,NAUX) !->ESSL
        DEALLOCATE(AUX)
        CALL DGESVS(AFACT,NM,BFACT,N,NEQ,S,X,M,N,M,TAU) !->ESSL
        DEALLOCATE(S)
        DEALLOCATE(BFACT)
        DEALLOCATE(AFACT)
#ELSE
        ALLOCATE(AFACT(N,M))
        AFACT(:,:)=A(:,:)
        ALLOCATE(BFACT(NM,NEQ))
        BFACT(1:N,:)=B(:,:)
        ALLOCATE(S(M))    !SINGULAR VALUES
        NAUX=3*NM+MAX(2*MIN(M,N),NM,NEQ)
        ALLOCATE(AUX(NAUX))
        CALL DGELSS(N,M,NEQ,AFACT,NM,BFACT,N,S,TAU,IRANK,AUX,NAUX,INFO) !->LAPACK
        IF(INFO.NE.0) THEN
          IF(INFO.LT.0) THEN
            CALL ERROR$MSG('INCORRECT ARGUMENT FOR SINGULAR VALUE DECOMPOSITION')
            CALL ERROR$I4VAL('WRONG ARGUMENT NR.',-I)
          ELSE IF(INFO.GT.0) THEN
            CALL ERROR$MSG('SINGULAR VALUE DECOMPOSITION FAILED TO CONVERGE')
            CALL ERROR$I4VAL('#(SING. VALUES WITHOUT CONVERGENCE)',I)
          END IF
          CALL ERROR$STOP('LIB$MATRIXSOLVE')
        END IF
        DEALLOCATE(AUX)
        X(:,:)=BFACT(1:M,:)
        DEALLOCATE(S)
        DEALLOCATE(BFACT)
        DEALLOCATE(AFACT)
#ENDIF
      END IF
      RETURN
      END
!
!ISORT  SORTS ELEMENTS IN A SEQUENCE P904
!     .................................................................
      SUBROUTINE LIB$SORTI4(N,X)
!     ******************************************************************
!     **  SORTS THE ARRAY X IN ASCENDING ORDER                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(INOUT):: X(N)
      INTEGER(4)              :: IHELP
      INTEGER(4)              :: I,J
!     ******************************************************************
#IF DEFined(CPPVAR_ESSL)
      CALL ISORT(X,1,N)
#ELSE
      DO I=1,N-1
        DO J=I+1,N
          IF(X(I).GT.X(J)) THEN
            IHELP=X(I)
            X(I) =X(J)
            X(J) =IHELP
          ENDIF
        ENDDO
      ENDDO
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE lib$GETHOSTNAME(HOSTNAME)
!     *********************************************************************
!     **  COLLECTS THE HOST NAME OF THE EXECUTING MACHINE                **
!     *********************************************************************
      CHARACTER(*),INTENT(OUT)  :: HOSTNAME
      INTEGER(4)                :: RC
!     *********************************************************************
#IF DEFined(CPPVAR_XLFsup)
      RC=HOSTNM_(HOSTNAME)    ! xlf support library
      IF(RC.NE.0)HOSTNAME='UNKNOWN'
#ELSE
      HOSTNAME='UNKNOWN'
!     HOSTNM_=GETHOSTNAME(HOSTNAME,%VAL(LEN(HOSTNAME)))
#ENDIF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE lib$FLUSHFILE(N)
!     *********************************************************************
!     ** FLUSHES THE BUFFER FOR THE FILE CONNECTED TO FORTRAN UNIT N     **
!     *********************************************************************
      INTEGER(4),INTENT(IN) :: N
!     *********************************************************************
#IF DEFined(CPPVAR_XLFsup)
      CALL FLUSH_(N)  ! xlf uspport library
#ELif defined(CPPVAR_U77)
      CALL FLUSH(N) ! from absoft support library (underscore)
#ENDIF
      RETURN
      END 
!
#IFNDEF CPPVAR_ESSL
!
!     ..................................................................
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
!     ******************************************************************
!     ** DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN          **
!     ** ELIMINATION.                                                 **
!     **                                                               **
!     ** DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED       **
!     ** DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.     **
!     ** (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA).               **
!     **                                                              **
!     **  ON RETURN, A IS AN UPPER TRIANGULAR MATRIX                  **
!     **  AND THE MULTIPLIERS, WHICH WERE USED TO OBTAIN IT.          **
!     **  THE FACTORIZATION CAN BE WRITTEN  A = L*U,                  **
!     **  WHERE L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER         **
!     **  TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.            **
!     **                                                              **
!     **  INFO                                                        **
!     **     = 0  NORMAL VALUE.                                       **
!     **     = K  IF U(K,K).EQ.0.0 , THIS IS NOT AN ERROR CONDITION   **
!     **          FOR THIS SUBROUTINE, BUT IT DOES INDICATE           **
!     **          THAT DGESL OR DGEDI WILL DIVIDE BY ZERO IF CALLED.  **
!     **          USE RCOND IN DGECO FOR A RELIABLE INDICATION        **
!     **          OF SINGULARITY.                                     **
!     **                                                              **
!     **  SUBROUTINES AND FUNCTIONS CALLED:                           **
!     **     BLAS ROUTINES: DAXPY,DSCAL,IDAMAX                        **
!     **                                                              **
!     **  LINPACK. THIS VERSION DATED 08/14/78 .                      **
!     **  CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.**
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N        ! ORDER OF MATRIX A
      INTEGER(4),INTENT(IN)   :: LDA      ! LEADING DIM. OF MATRIX A
      REAL(8)   ,INTENT(INOUT):: A(LDA,N) ! MATRIX TO BE FACTORED
      INTEGER(4),INTENT(OUT)  :: IPVT(N)  ! PIVOT INDICES
      INTEGER(4),INTENT(OUT)  :: INFO     ! RETURN CODE
      REAL(8)                 :: T
      INTEGER(4)              :: IDAMAX,J,K,KP1,L,NM1
!     ******************************************************************
!     ==================================================================
!     ==  GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                  ==
!     ==================================================================
      INFO=0
      NM1 =N-1
      IF(NM1.LT.1) THEN
        IPVT(N)=N
         IF(A(N,N).EQ.0.0D0) INFO=N
        RETURN
      END IF
!
      DO K=1,NM1
        KP1=K+1
!       ==  FIND L = PIVOT INDEX ======================================
        L=IDAMAX(N-K+1,A(K,K),1) + K-1  !idamax<-blas library
        IPVT(K)=L
!       == ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED ======
        IF(A(L,K).NE.0.0D0) THEN
!         ==  INTERCHANGE IF NECESSARY ===============================
          IF(L.NE.K) THEN
            T=A(L,K)
            A(L,K)=A(K,K)
            A(K,K)=T
          END IF
!         == COMPUTE MULTIPLIERS ======================================
          T=-1.0D0/A(K,K)
          CALL DSCAL(N-K,T,A(K+1,K),1)
!         ==  ROW ELIMINATION WITH COLUMN INDEXING  ===================
          DO J=KP1,N
            T=A(L,J)
            IF(L.NE.K) THEN
              A(L,J)=A(K,J)
              A(K,J)=T
            END IF
            CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
          ENDDO
        ELSE
          INFO=K
        END IF
      ENDDO
      IPVT(N)=N
      IF(A(N,N).EQ.0.0D0) INFO=N
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
!     ******************************************************************
!     ** DGESL SOLVES THE DOUBLE PRECISION SYSTEM                     **
!     **       A*X=B  OR  TRANS(A)*X=B                                **
!     **  USING THE FACTORS COMPUTED BY DGECO OR DGEFA.               **
!     **                                                              **
!     **  JOB=0       TO SOLVE  A*X = B ,                             **
!     **     =NONZERO TO SOLVE  TRANS(A)*X = B,                       **
!     **                        ,WHERE TRANS(A)  IS THE TRANSPOSE.    **
!     **                                                              **
!     **  ERROR CONDITION: A DIVISION BY ZERO WILL OCCUR              **
!     **    IF THE INPUT FACTOR  CONTAINS A ZERO ON THE DIAGONAL.     **
!     **    TECHNICALLY THIS INDICATES SINGULARITY BUT IT IS OFTEN    **
!     **    CAUSED BY IMPROPER ARGUMENTS OR IMPROPER SETTING OF LDA.  **
!     **    IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY **
!     **    AND IF DGECO HAS SET RCOND>0. OR DGEFA HAS SET INFO=0.    **
!     **                                                              **
!     **  TO COMPUTE  INVERSE(A)*C WHERE C IS A MATRIX WITH P COLUMNS **
!     **    CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                          **
!     **    IF (RCOND IS TOO SMALL) THEN                              **
!     **      PRINT ERROR MESSAGE AND STOP...                         **
!     **    END IF                                                    **
!     **    DO J=1,P                                                  **
!     **      CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                       **
!     **    ENDDO                                                     **
!     **                                                              **
!     **   SUBROUTINES AND FUNCTIONS CALLED                           **
!     **      BLAS ROUTINES: DAXPY,DDOT                               **
!     **                                                              **
!     **  LINPACK. THIS VERSION DATED 08/14/78 .                      **
!     **  CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.**
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N        ! ORDER OF A
      INTEGER(4),INTENT(IN)   :: LDA      ! LEADING DIMENSION OF A
      REAL(8)   ,INTENT(IN)   :: A(LDA,N) ! OUTPUT OF DGECO OR DGEFA
      REAL(8)   ,INTENT(INOUT):: B(N)     ! INPUT:  RIGHT HAND SIDE B
                                          ! OUTPUT: SOLUTION VECTOR X
      INTEGER(4),INTENT(IN)   :: IPVT(N)  ! PIVOT VECTOR FROM DGECO OR DEGFA
      INTEGER(4),INTENT(IN)   :: JOB      ! A*X=B OR TRANS(A)X=B
      REAL(8)                 :: DDOT
      REAL(8)                 :: T
      INTEGER(4)              :: K,KB,L,NM1
!     ******************************************************************
      NM1=N-1
!     ==================================================================
!     == JOB=0: SOLVE A*X=B                                           ==
!     ==================================================================
      IF(JOB.EQ.0) THEN
!       == FIRST SOLVE  L*Y = B ========================================
        IF(NM1.GE.1) THEN
          DO K=1,NM1
            L = IPVT(K)
            T = B(L)
            IF(L.NE.K) THEN
              B(L) = B(K)
              B(K) = T
            END IF
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
          ENDDO
        END IF
!       ==  NOW SOLVE  U*X = Y ==========================================
        DO KB=1,N
          K =N+1-KB
          B(K)=B(K)/A(K,K)
          T =-B(K)
          CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
        ENDDO
!
!     ==================================================================
!     == JOB = NONZERO, SOLVE  TRANS(A) * X = B                       ==
!     ==================================================================
      ELSE
!       == FIRST SOLVE  TRANS(U)*Y=B ===================================
        DO K = 1, N
          T = DDOT(K-1,A(1,K),1,B(1),1)
          B(K) = (B(K) - T)/A(K,K)
        ENDDO
!       == NOW SOLVE TRANS(L)*X = Y =====================================
        IF (NM1.GE.1) THEN
          DO KB=1,NM1
            K = N-KB
            B(K)=B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF(L.NE.K) THEN
              T = B(L)
              B(L) = B(K)
              B(K) = T
            END IF
          ENDDO
        END IF
      END IF
      RETURN
      END
#endif
#if defined(cppvar_fftpack)
!
!     ..................................................................
      SUBROUTINE CFFTB(N,C,WA,IFAC)
!     ******************************************************************
!     ** CFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER         **
!     ** TRANSFORM (THE FOURIER SYNTHESIS).                           **
!     ** EQUIVALENTLY, CFFTB COMPUTES A COMPLEX PERIODIC SEQUENCE     **
!     ** FROM ITS FOURIER COEFFICIENTS.                               **
!     **                                                              **
!     **  THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED     **
!     **  TRANSFORM THE OUTPUT of CFFTF! MUST BE DIVIDED BY N.        **
!     **  OTHERWISE A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB       **
!     **  WILL MULTIPLY THE SEQUENCE BY N.                            **
!     **                                                              **
!     **  THE ARRAYS WA AND IFAC WHICH ARE USED BY SUBROUTINE CFFTF   **
!     **  MUST BE INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).   **
!     **  THE ARRAYS DEPEND ON THE LENGTH N OF THE TRANSFORM          **
!     **  THE SAME ARRAYS CAN BE USED FOR FORWARD AN BACK TRANSFORM   **
!     **  THE WORK ARRAYS ARE NOT CHANGED BY THE FOURIER TANSFORM     **
!     **  ROUTINES CFFTF AND CFFTB. AND CAN BE KEPT FOR SUBSEQUENT    **
!     **  TRANSFORMS WITH THE SAME LENGTH                             **
!     **                                                              **
!     **  THE METHOD IS MORE EFFICIENT IF THE LENGTH OF THE TRANSFORM **
!     **  IS A PRODUCT OF SMALL PRIME NUMBERS                         **
!     **                                                              **
!     **  THE COMPLEX SEQUENCE C OF LENGTH N IS TREATED INTERNALLY    **
!     **  AS A REAL ARRAY OF LENGTH 2*N                               **
!     **                                                              **
!     **  RESULT                                                      **
!     **     FOR K,J=1,...,N AND I=SQRT(-1)                           **
!     **     C(J)=SUM_{K=1}^N C(K)*EXP(+I*(J-1)*(K-1)*2*PI/N)         **
!     **                                                              **
!     **  THE ARRAYS WA,IFAC CONTAIN INITIALIZATION CALCULATIONS      **
!     **  WHICH MUST NOT BE DESTROYED BETWEEN CALLS                   **
!     **  OF SUBROUTINE CFFTF OR CFFTB                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      REAL(8)   ,INTENT(INOUT):: C(2*N)
      REAL(8)   ,INTENT(IN)   :: WA(2*N)
      INTEGER(4),INTENT(IN)   :: IFAC(30)
      REAL(8)                 :: CH(2*N)
      INTEGER(4)              :: NF,NA,L1,IW,L2,IDO,IDOT,IDL1,N2
      INTEGER(4)              :: IP
      INTEGER(4)              :: IX2,IX3,IX4
      INTEGER(4)              :: K1
      INTEGER(4)              :: I
      INTEGER(4)              :: NAC
!     ******************************************************************
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1=1,NF
        IP = IFAC(K1+2)
        L2 = IP*L1
        IDO = N/L2
        IDOT = IDO+IDO
        IDL1 = IDOT*L1
        IF(IP.EQ.4) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSB4(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
          ELSE
            CALL PASSB4(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.2) THEN
          IF(NA.EQ.0) THEN
            CALL PASSB2(IDOT,L1,C,CH,WA(IW))
          ELSE
            CALL PASSB2(IDOT,L1,CH,C,WA(IW))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.3) THEN
          IX2 = IW+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSB3(IDOT,L1,C,CH,WA(IW),WA(IX2))
          ELSE
            CALL PASSB3(IDOT,L1,CH,C,WA(IW),WA(IX2))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.5) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IX4 = IX3+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSB5(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          ELSE
            CALL PASSB5(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          END IF
          NA = 1-NA
        ELSE
          IF(NA.EQ.0) THEN
            CALL PASSB(NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
          ELSE
            CALL PASSB(NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
          END IF
          IF(NAC.NE.0) NA = 1-NA
        END IF
        L1 = L2
        IW = IW+(IP-1)*IDOT
      ENDDO
!
      IF(NA.NE.0) THEN
        N2 = N+N
        DO I=1,N2
          C(I) = CH(I)
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CFFTF(N,C,WA,IFAC)
!     ******************************************************************
!     **  CFFTF COMPUTES THE FORWARD COMPLEX DISCRETE FOURIER         **
!     **  TRANSFORM (THE FOURIER ANALYSIS).                           **
!     **  EQUIVALENTLY, CFFTF COMPUTES THE FOURIER COEFFICIENTS       **
!     **  OF A COMPLEX PERIODIC SEQUENCE.                             **
!     **                                                              **
!     **  THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED     **
!     **  TRANSFORM THE OUTPUT MUST BE DIVIDED BY N.                  **
!     **  OTHERWISE A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB       **
!     **  WILL MULTIPLY THE SEQUENCE BY N.                            **
!     **                                                              **
!     **  THE ARRAYS WA AND IFAC WHICH ARE USED BY SUBROUTINE CFFTF   **
!     **  MUST BE INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).   **
!     **  THE ARRAYS DEPEND ON THE LENGTH N OF THE TRANSFORM          **
!     **  THE SAME ARRAYS CAN BE USED FOR FORWARD AN BACK TRANSFORM   **
!     **  THE WORK ARRAYS ARE NOT CHANGED BY THE FOURIER TANSFORM     **
!     **  ROUTINES CFFTF AND CFFTB. AND CAN BE KEPT FOR SUBSEQUENT    **
!     **  TRANSFORMS WITH THE SAME LENGTH                             **
!     **                                                              **
!     **  THE METHOD IS MORE EFFICIENT IF THE LENGTH OF THE TRANSFORM **
!     **  IS A PRODUCT OF SMALL PRIME NUMBERS                         **
!     **                                                              **
!     **  THE COMPLEX SEQUENCE C OF LENGTH N IS TREATED INTERNALLY    **
!     **  AS A REAL ARRAY OF LENGTH 2*N                               **
!     **                                                              **
!     **  RESULT                                                      **
!     **     FOR K,J=1,...,N AND I=SQRT(-1)                           **
!     **     C(J)=SUM_{K=1}^N C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)         **
!     **                                                              **
!     **  THE ARRAYS WA,IFAC CONTAIN INITIALIZATION CALCULATIONS      **
!     **  WHICH MUST NOT BE DESTROYED BETWEEN CALLS                   **
!     **  OF SUBROUTINE CFFTF OR CFFTB                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N        ! LENGTH OF FFT
      REAL(8)   ,INTENT(INOUT) :: C(2*N)   ! SEQUENCE TO BE TRANSFORMED
      REAL(8)   ,INTENT(IN)    :: WA(2*N)  ! WORK ARRAY
      INTEGER(4),INTENT(IN)    :: IFAC(30) ! WORK ARRAY
      REAL(8)                  :: CH(2*N)  ! WORK ARRAY
      INTEGER(4)               :: NF,NA,L1,IW
      INTEGER(4)               :: IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4
      INTEGER(4)               :: N2
      INTEGER(4)               :: I,K1
      INTEGER(4)               :: NAC
!     ******************************************************************
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1=1,NF
        IP = IFAC(K1+2)
        L2 = IP*L1
        IDO = N/L2
        IDOT = IDO+IDO
        IDL1 = IDOT*L1
        IF(IP.EQ.4) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSF4(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
          ELSE
            CALL PASSF4(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.2) THEN
          IF(NA.EQ.0) THEN
            CALL PASSF2(IDOT,L1,C,CH,WA(IW))
          ELSE
            CALL PASSF2(IDOT,L1,CH,C,WA(IW))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.3) THEN
          IX2 = IW+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSF3(IDOT,L1,C,CH,WA(IW),WA(IX2))
          ELSE
            CALL PASSF3(IDOT,L1,CH,C,WA(IW),WA(IX2))
          END IF
          NA = 1-NA
        ELSE IF(IP.EQ.5) THEN
          IX2 = IW+IDOT
          IX3 = IX2+IDOT
          IX4 = IX3+IDOT
          IF(NA.EQ.0) THEN
            CALL PASSF5(IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          ELSE
            CALL PASSF5(IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
          ENDIF
          NA = 1-NA
        ELSE
          IF(NA.NE.0) THEN
            CALL PASSF(NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
          ELSE
            CALL PASSF(NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
          END IF
          IF (NAC .NE. 0) NA = 1-NA
        END IF
        L1 = L2
        IW = IW+(IP-1)*IDOT
      ENDDO
      IF(NA.NE.0) THEN
        N2 = N+N
        DO I=1,N2
          C(I) = CH(I)
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE CFFTI(N,WA,IFAC)
!     ******************************************************************
!     **  CFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN BOTH     **
!     **  CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH **
!     **  A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED    **
!     **  AND STORED IN wa and ifac.                                  **
!     **                                                              **
!     **  THE ARRAYS WA,IFAC CONTAIN INITIALIZATION CALCULATIONS      **
!     **  WHICH MUST NOT BE DESTROYED BETWEEN CALLS                   **
!     **  OF SUBROUTINE CFFTF OR CFFTB                                **
!     **                                                              **
!     ** ifac contains (1) the length of the transform,               **
!     **               (2) the number of factors and                  **
!     **           and (3,nf+2) nf factors                            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N        ! length of the fourier transform
      INTEGER(4),INTENT(OUT) :: IFAC(30) ! prime factorization
      REAL(8)   ,INTENT(OUT) :: WA(2*n)  ! trigonometric functions
      INTEGER(4),PARAMETER   :: NTRYH(4)=(/3,4,2,5/)
      INTEGER(4)             :: NL,NF,NQ,NR
      INTEGER(4)             :: NTRY,IB,L1,IP,LD,L2,IDO,IDOT,IPM,I1
      INTEGER(4)             :: I,J,II,K1
      REAL(8)                :: TPI,ARGH,FI,ARGLD,ARG
      REAL(8)                :: PI
!     ******************************************************************
!
!     ==================================================================
!     == PRIME FACTORIZATION                                          ==
!     ==================================================================
      NL = N    ! REMAINDER OF FACTORED LENGTH OF TRANSFORM
      NF = 0    !#(FACTORS)
      J=0
 101  CONTINUE
!     __ SELECT A NEW PRIME NUMBER NTRY FROM THE ARRAY 3,4,2,5,7,9,11,...
      J=J+1
      IF(J.LE.4) THEN
        NTRY=NTRYH(J)
      ELSE
        NTRY=NTRY+2
      END IF
104   CONTINUE
!     __  CHECK IF PRIME NUMBER IS A FACTOR_____________________________
      NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF(NR.NE.0) GOTO 101 ! if ntry is not a factor PICK A larger prime
!     __NEW FACTOR FOUND, ADD TO ARRAY IFAC_____________________________
      NF = NF+1
      IF(NF.GT.13) THEN
        CALL ERROR$MSG('ARRAY SIZE FOR FACTORIZATION TOO SMALL')
        CALL ERROR$MSG('INCREASE ARRAY SIZE FOR IFAC IN CFFTI,CFFTB,CFFTF')
        CALL ERROR$STOP('CFFTI')
      END IF
      IFAC(NF+2) = NTRY
      NL = NQ
!     __PLACE ANY PRIME FACTOR 2 UPFRONT________________________________
      IF(NTRY.EQ.2.AND.NF.NE.1) THEN  
        DO I=2,NF
          IB = NF-I+2
          IFAC(IB+2) = IFAC(IB+1)
        ENDDO
        IFAC(3) = 2
      END IF
      IF(NL.NE.1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
!
!     ==================================================================
!     == TABULATE TRIGONOMETRIC FUNCTIONS                             ==
!     ==================================================================
!     PI=3.14159265358979323846264338327950288419716939937510582097494D0
      PI=4.D0*DATAN(1.D0)
      TPI = 2.0D0*PI
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO K1=1,NF
        IP = IFAC(K1+2)
        LD = 0
        L2 = L1*IP
        IDO = N/L2
        IDOT = IDO+IDO+2
        IPM = IP-1
        DO J=1,IPM
          I1 = I
          WA(I-1) = 1.0D0
          WA(I)   = 0.0D0
          LD = LD+L1
          FI = 0.0D0
          ARGLD = FLOAT(LD)*ARGH
          DO II=4,IDOT,2
            I = I+2
            FI = FI+1.0D0
            ARG = FI*ARGLD
            WA(I-1)=COS(ARG)
            WA(I)  =SIN(ARG)
          ENDDO
          IF(IP.GT.5) THEN
            WA(I1-1) = WA(I-1)
            WA(I1)   = WA(I)
          END IF
        ENDDO
        L1 = L2
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)  :: NAC
      INTEGER(4),INTENT(IN)   :: IDO
      INTEGER(4),INTENT(IN)   :: IP
      INTEGER(4),INTENT(IN)   :: L1
      INTEGER(4),INTENT(IN)   :: IDL1
      REAL(8)   ,INTENT(OUT)  :: CH(IDO,L1,IP)
      REAL(8)   ,INTENT(IN)   :: CC(IDO,IP,L1)
      REAL(8)   ,INTENT(OUT)  :: C1(IDO,L1,IP)
      REAL(8)   ,INTENT(OUT)  :: C2(IDL1,IP)
      REAL(8)   ,INTENT(INOUT):: CH2(IDL1,IP)
      REAL(8)   ,INTENT(IN)   :: WA(*) 
!MAX(1+(IP-1)*IDO,2+((IP-1)/2)*IDO)??
      INTEGER(4)              :: IDOT,NT,IPP2,IPPH,IDP,IDL,INC,IDLJ,IDIJ,IDJ
      INTEGER(4)              :: JC,LC
      INTEGER(4)              :: I,J,K,L,IK
      REAL(8)                 :: WAR,WAI
!     ******************************************************************
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF(IDO.GE.L1) THEN
        DO J=2,IPPH
          JC = IPP2-J
          DO K=1,L1
            DO I=1,IDO
              CH(I,K,J)  = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO K=1,L1
          DO I=1,IDO
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      ELSE
        DO J=2,IPPH
          JC = IPP2-J
          DO I=1,IDO
            DO K=1,L1
              CH(I,K,J)  = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO I=1,IDO
          DO K=1,L1
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      ENDIF
!
      IDL = 2-IDO
      INC = 0
      DO L=2,IPPH
        LC = IPP2-L
        IDL = IDL+IDO
        DO IK=1,IDL1
          C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
          C2(IK,LC)=           WA(IDL)  *CH2(IK,IP)
        ENDDO
        IDLJ = IDL
        INC = INC+IDO
        DO J=3,IPPH
          JC = IPP2-J
          IDLJ = IDLJ+INC
          IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
          WAR = WA(IDLJ-1)
          WAI = WA(IDLJ)
          DO IK=1,IDL1
            C2(IK,L)  = C2(IK,L) +WAR*CH2(IK,J)
            C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
          ENDDO
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        DO IK=1,IDL1
          CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        JC = IPP2-J
        DO IK=2,IDL1,2
          CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
          CH2(IK-1,JC)= C2(IK-1,J)+C2(IK,JC)
          CH2(IK,J)   = C2(IK,J)  +C2(IK-1,JC)
          CH2(IK,JC)  = C2(IK,J)  -C2(IK-1,JC)
        ENDDO
      ENDDO
!
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO IK=1,IDL1
        C2(IK,1) = CH2(IK,1)
      ENDDO
      DO J=2,IP
        DO K=1,L1
          C1(1,K,J) = CH(1,K,J)
          C1(2,K,J) = CH(2,K,J)
        ENDDO
      ENDDO
!
      IF(IDOT.LE.L1) THEN
        IDIJ = 0
        DO J=2,IP
         IDIJ = IDIJ+2
         DO I=4,IDO,2
           IDIJ = IDIJ+2
           DO K=1,L1
             C1(I-1,K,J)=WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
             C1(I,K,J)  =WA(IDIJ-1)*CH(I,K,J)  +WA(IDIJ)*CH(I-1,K,J)
           ENDDO
         ENDDO
       ENDDO
      ELSE
        IDJ = 2-IDO
        DO J=2,IP
          IDJ = IDJ+IDO
          DO K=1,L1
            IDIJ = IDJ
            DO I=4,IDO,2
              IDIJ = IDIJ+2
              C1(I-1,K,J)= WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
              C1(I,K,J)  = WA(IDIJ-1)*CH(I,K,J)  +WA(IDIJ)*CH(I-1,K,J)
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB2(IDO,L1,CC,CH,WA1)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,2,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,2)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR2,TI2
!     ******************************************************************
      IF(IDO.LE.2) THEN
        DO K=1,L1
          CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
          CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
          CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
          CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,3,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,3)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)               :: TR2,TI2
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)   ,PARAMETER  :: TAUR=-0.5D0
      REAL(8)   ,PARAMETER  :: TAUI=0.866025403784439D0
      INTEGER(4)            :: K,I
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO K=1,L1
          TR2 = CC(1,2,K)+CC(1,3,K)
          CR2 = CC(1,1,K)+TAUR*TR2
          CH(1,K,1) = CC(1,1,K)+TR2
          TI2 = CC(2,2,K)+CC(2,3,K)
          CI2 = CC(2,1,K)+TAUR*TI2
          CH(2,K,1) = CC(2,1,K)+TI2
          CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
          CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
          CH(1,K,2) = CR2-CI3
          CH(1,K,3) = CR2+CI3
          CH(2,K,2) = CI2+CR3
          CH(2,K,3) = CI2-CR3
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB4(IDO,L1,CC,CH,WA1,WA2,WA3)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,4,L1)
      REAL(8)   ,INTENt(OUT):: CH(IDO,L1,4)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      REAL(8)               :: TR1,TI1
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
      INTEGER(4)            :: I,K
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO K=1,L1
          TI1 = CC(2,1,K)-CC(2,3,K)
          TI2 = CC(2,1,K)+CC(2,3,K)
          TR4 = CC(2,4,K)-CC(2,2,K)
          TI3 = CC(2,2,K)+CC(2,4,K)
          TR1 = CC(1,1,K)-CC(1,3,K)
          TR2 = CC(1,1,K)+CC(1,3,K)
          TI4 = CC(1,2,K)-CC(1,4,K)
          TR3 = CC(1,2,K)+CC(1,4,K)
          CH(1,K,1) = TR2+TR3
          CH(1,K,3) = TR2-TR3
          CH(2,K,1) = TI2+TI3
          CH(2,K,3) = TI2-TI3
          CH(1,K,2) = TR1+TR4
          CH(1,K,4) = TR1-TR4
          CH(2,K,2) = TI1+TI4
          CH(2,K,4) = TI1-TI4
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,5,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,5)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      REAL(8)   ,INTENT(IN) :: WA4(IDO)
      REAL(8)   ,PARAMETER  :: TR11= 0.309016994374947D0
      REAL(8)   ,PARAMETER  :: TI11= 0.951056516295154D0
      REAL(8)   ,PARAMETER  :: TR12=-0.809016994374947D0
      REAL(8)   ,PARAMETER  :: TI12= 0.587785252292473D0
      INTEGER(4)            :: I,K
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: TR5,TI5
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
      REAL(8)               :: CR5,CI5
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)               :: DR4,DI4
      REAL(8)               :: DR5,DI5
!     ******************************************************************
      IF(IDO.EQ.2) THEN
      DO K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1)   = CC(I,1,K)  +TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)  +TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)  +TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2)   = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3)   = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4)   = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5)   = WA4(I-1)*DI5+WA4(I)*DR5
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)  :: NAC
      INTEGER(4),INTENT(IN)   :: IDO
      INTEGER(4),INTENT(IN)   :: IP
      INTEGER(4),INTENT(IN)   :: L1
      INTEGER(4),INTENT(IN)   :: IDL1
      REAL(8)   ,INTENT(IN)   :: CC(IDO,IP,L1)
      REAL(8)   ,INTENT(OUT)  :: C1(IDO,L1,IP)
      REAL(8)   ,INTENT(OUT)  :: C2(IDL1,IP)
      REAL(8)   ,INTENT(OUT)  :: CH(IDO,L1,IP)
      REAL(8)   ,INTENT(INOUT):: CH2(IDL1,IP)
      REAL(8)   ,INTENT(IN)   :: WA(*)
      INTEGER(4)              :: IDOT,NT,IPP2,IPPH,IDP,INC
      INTEGER(4)              :: IDJ,IDIJ,IDL,IDLJ
      INTEGER(4)              :: JC,LC
      INTEGER(4)              :: I,J,K,L,IK
      REAL(8)                 :: WAR,WAI
!     ******************************************************************
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
!
      IF(IDO.GE.L1) THEN
        DO J=2,IPPH
          JC = IPP2-J
          DO K=1,L1
            DO I=1,IDO
              CH(I,K,J)  = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO K=1,L1
          DO I=1,IDO
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      ELSE
        DO J=2,IPPH
          JC = IPP2-J
          DO I=1,IDO
            DO K=1,L1
              CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
              CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
            ENDDO
          ENDDO
        ENDDO
        DO I=1,IDO
          DO K=1,L1
            CH(I,K,1) = CC(I,1,K)
          ENDDO
        ENDDO
      END IF
!
      IDL = 2-IDO
      INC = 0
      DO L=2,IPPH
        LC = IPP2-L
        IDL = IDL+IDO
        DO IK=1,IDL1
          C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
          C2(IK,LC)=          -WA(IDL)  *CH2(IK,IP)
        ENDDO
        IDLJ = IDL
        INC = INC+IDO
        DO J=3,IPPH
          JC = IPP2-J
          IDLJ = IDLJ+INC
          IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
          WAR = WA(IDLJ-1)
          WAI = WA(IDLJ)
          DO IK=1,IDL1
            C2(IK,L)  = C2(IK,L) +WAR*CH2(IK,J)
            C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
          ENDDO
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        DO IK=1,IDL1
          CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
        ENDDO
      ENDDO
!
      DO J=2,IPPH
        JC = IPP2-J
        DO IK=2,IDL1,2
          CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
          CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
          CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
          CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
        ENDDO
      ENDDO
!
      NAC = 1
      IF(IDO.EQ.2) RETURN
      NAC = 0
      DO IK=1,IDL1
        C2(IK,1) = CH2(IK,1)
      ENDDO
      DO J=2,IP
        DO K=1,L1
          C1(1,K,J) = CH(1,K,J)
          C1(2,K,J) = CH(2,K,J)
        ENDDO
      ENDDO
!
      IF(IDOT.LE.L1) THEN
        IDIJ = 0
        DO J=2,IP
          IDIJ = IDIJ+2
          DO I=4,IDO,2
            IDIJ = IDIJ+2
            DO K=1,L1
              C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
              C1(I,K,J)   = WA(IDIJ-1)*CH(I,K,J)  -WA(IDIJ)*CH(I-1,K,J)
            ENDDO
          ENDDO
        ENDDO
        RETURN
      END IF
!
      IDJ = 2-IDO
      DO J=2,IP
        IDJ = IDJ+IDO
        DO K=1,L1
          IDIJ = IDJ
          DO I=4,IDO,2
            IDIJ = IDIJ+2
            C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
            C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF2(IDO,L1,CC,CH,WA1)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,2,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,2)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR2,TI2
!     ******************************************************************
      IF(IDO.LE.2) THEN
        DO K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF3(IDO,L1,CC,CH,WA1,WA2)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,3,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,3)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR2,TI2
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)   ,PARAMETER  :: TAUR=-0.5D0
      REAL(8)   ,PARAMETER  :: TAUI=-0.866025403784439D0
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO  K=1,L1
          TR2 = CC(1,2,K)+CC(1,3,K)
          CR2 = CC(1,1,K)+TAUR*TR2
          CH(1,K,1) = CC(1,1,K)+TR2
          TI2 = CC(2,2,K)+CC(2,3,K)
          CI2 = CC(2,1,K)+TAUR*TI2
          CH(2,K,1) = CC(2,1,K)+TI2
          CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
          CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
          CH(1,K,2) = CR2-CI3
          CH(1,K,3) = CR2+CI3
          CH(2,K,2) = CI2+CR3
          CH(2,K,3) = CI2-CR3
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF4(IDO,L1,CC,CH,WA1,WA2,WA3)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,4,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,4)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR1,TI1
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
!     ******************************************************************
      IF(IDO.EQ.2) THEN
        DO K=1,L1
          TI1 = CC(2,1,K)-CC(2,3,K)
          TI2 = CC(2,1,K)+CC(2,3,K)
          TR4 = CC(2,2,K)-CC(2,4,K)
          TI3 = CC(2,2,K)+CC(2,4,K)
          TR1 = CC(1,1,K)-CC(1,3,K)
          TR2 = CC(1,1,K)+CC(1,3,K)
          TI4 = CC(1,4,K)-CC(1,2,K)
          TR3 = CC(1,2,K)+CC(1,4,K)
          CH(1,K,1) = TR2+TR3
          CH(1,K,3) = TR2-TR3
          CH(2,K,1) = TI2+TI3
          CH(2,K,3) = TI2-TI3
          CH(1,K,2) = TR1+TR4
          CH(1,K,4) = TR1-TR4
          CH(2,K,2) = TI1+TI4
          CH(2,K,4) = TI1-TI4
        ENDDO
        RETURN
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IDO
      INTEGER(4),INTENT(IN) :: L1
      REAL(8)   ,INTENT(IN) :: CC(IDO,5,L1)
      REAL(8)   ,INTENT(OUT):: CH(IDO,L1,5)
      REAL(8)   ,INTENT(IN) :: WA1(IDO)
      REAL(8)   ,INTENT(IN) :: WA2(IDO)
      REAL(8)   ,INTENT(IN) :: WA3(IDO)
      REAL(8)   ,INTENT(IN) :: WA4(IDO)
      INTEGER(4)            :: K,I
      REAL(8)               :: TR1,TI1
      REAL(8)               :: TR2,TI2
      REAL(8)               :: TR3,TI3
      REAL(8)               :: TR4,TI4
      REAL(8)               :: TR5,TI5
      REAL(8)               :: CR2,CI2
      REAL(8)               :: CR3,CI3
      REAL(8)               :: CR4,CI4
      REAL(8)               :: CR5,CI5
      REAL(8)               :: DR2,DI2
      REAL(8)               :: DR3,DI3
      REAL(8)               :: DR4,DI4
      REAL(8)               :: DR5,DI5
      REAL(8)   ,PARAMETER  :: TR11=+0.309016994374947D0
      REAL(8)   ,PARAMETER  :: TI11=-0.951056516295154D0
      REAL(8)   ,PARAMETER  :: TR12=-0.809016994374947D0
      REAL(8)   ,PARAMETER  :: TI12=-0.587785252292473D0
!     ******************************************************************
      IF(IDO.NE.2) THEN
        DO K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
        ENDDO
      ELSE
        DO K=1,L1
          DO I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
          ENDDO
        ENDDO
      END IF
      RETURN
      END
#endif











