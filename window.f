      SUBROUTINE WINDOW_TS(MDIM,MSAMP,MREC,MSRC, NDIM,NSAMP,NREC,NSRC,
     ;                     ITYPE,NPCT,DT,TWIN, TS,IERR) 
!
!     Windows the time series data base using Steanrds and David's
!     windowing routines
!
!     INPUT      MEANING
!     -----      -------
!     DT         sampling period (same units as twin)
!     ITYPE      (1) Rectangular window
!                (2) Tapered rectangle
!                (3) Triangle
!                (4) Hanning
!                (5) Hamming
!                (6) Blackman 
!     MDIM       leading dimension
!     MREC       leading dimension
!     MSAMP      leading dimension
!     MSRC       leading dimension
!     NDIM       number of components 
!     NREC       number of receivers
!     NSAMP      number of samples
!     NSRC       number of sources
!     TS         time series to window
!     TWIN       window times (start,stop)
!
!     OUTPUT     MEANING
!     ------     ------- 
!     TS         windowed time series 
!
!.... variable declarations
      IMPLICIT NONE
      REAL*4, INTENT(INOUT) :: TS(MDIM,MSAMP,MREC,*) 
      REAL*4, INTENT(IN) :: TWIN(MREC,MSRC,*), DT 
      INTEGER*4, INTENT(IN) :: MDIM,MSAMP,MREC,MSRC, NDIM,NSAMP,NREC,
     ;                         NSRC,NPCT,ITYPE
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      REAL*4, ALLOCATABLE :: TEMP(:)
      REAL*4 TSV, T0, T1 
      INTEGER*4 ISRC, IREC, I, I1,I2, LX, K, ION, IOFF 
      REAL*4, PARAMETER :: FLIN = 0.5
      REAL*8, PARAMETER :: FHIN = 0.0
      INTEGER*4, PARAMETER :: NLIN = 8
      INTEGER*4, PARAMETER :: NHIN = 0
!
!----------------------------------------------------------------------!
!
!.... window 
      IERR = 0 
      LX = NSAMP - 1 
      ALLOCATE(TEMP(NSAMP))
      DO 1 ISRC=1,NSRC
         DO 2 IREC=1,NREC
            T0 = TWIN(IREC,ISRC,1)
            T1 = TWIN(IREC,ISRC,2)
            I1 = MAX(    1,INT(T0/DT) + 1)
            I2 = MIN(NSAMP,INT(T1/DT) + 1)
            ION  = 1 + INT( FLOAT(I2 - I1 + 1)*FLOAT(NPCT)/100.0 + 0.5)
            IOFF = I2 - I1 + 1 
     ;           - INT( FLOAT(I2 - I1 + 1)*FLOAT(NPCT)/100.0 + 0.5) 
            LX = I2 - I1 
            DO 3 I=1,NDIM
               CALL TAPER4(TS(I,1:NSAMP,IREC,ISRC),I1,I2,NSAMP,0.10)
               !CALL SCOPY(NSAMP,TS(I,1:NSAMP,IREC,ISRC),1,TEMP,1)
               TEMP(1:NSAMP) = TS(I,1:NSAMP,IREC,ISRC)
               CALL IIR(TEMP, TS(I,1:NSAMP,IREC,ISRC), DT, NSAMP, 
     ;                  NLIN, FLIN, NHIN, FHIN)
               print *, minval(temp),maxval(temp), 
     ; minval(ts(i,1:nsamp,irec,isrc)),maxval(ts(i,1:nsamp,irec,isrc))
    3       CONTINUE
!           IF (LX.LT.0) THEN
!              WRITE(*,*) 'window: Window length is negative!'
!              IERR = 1
!              RETURN
!           ENDIF
!           DO 3 I=1,NDIM
!              CALL SPMASK_V2(2,TS(I,I1:I2,IREC,ISRC),LX,NPCT,ITYPE,TSV,
!    ;                        ION,IOFF, IERR)
!              !CALL SPMASK(TS(I,I1:I2,IREC,ISRC),LX,ITYPE,TSV,IERR)
!              IF (IERR.NE.0) THEN
!                 WRITE(*,*) 'window: Error calling SPMASK_V2',IERR
!                 RETURN
!              ENDIF
!              DO 4 K=1,I1-1
!                 TS(I,K,IREC,ISRC) = 0.0
!   4          CONTINUE
!              DO 5 K=I2+1,NSAMP
!                 TS(I,K,IREC,ISRC) = 0.0
!   5          CONTINUE 
!   3       CONTINUE !loop on components
    2    CONTINUE !loop on receivers
    1 CONTINUE !loop on sources
      DEALLOCATE(TEMP)
      RETURN 
      END
!                                                                      !
!======================================================================!
!                                                                      !
        subroutine taper4(a,ibeg, iend, n, tfrac)
c
c       Apply a ten percent sine squared filter to the
c       n point time series stored in array a.
c
c       ibeg   taper start
c       iend   taper end
c
      real*4, intent(in) :: tfrac
      integer*4, intent(in) :: n, ibeg, iend
      real*4, intent(inout) :: a(n)
      real*4 fac, w
      integer*4 n1, l, i, i1, i2, npt 
      real*4, parameter :: pi = 3.141592653589793

      if (ibeg.lt.0) then
         write(*,*) 'taper: ibeg < 0;, skipping!'
         return
      endif
      if (iend.gt.n) then
         write(*,*) 'taper: iend > 0; skipping!'
         return
      endif
      do 1 i = 1, ibeg
         a(i) = 0.0
    1 continue 
      do 2 i = iend, n
         a(i) = 0.0
    2 continue 
      npt = iend - ibeg + 1 
      n1 = int(float(npt)*tfrac)
      if(n1.le.1) return
      w = pi/(2.0*float(n1))
      l = iend
      i1 = ibeg + 1 
      i2 = i1 + n1
      do 3 i = i1, i2
         l = l - 1 
         fac = sin(w*float(i-i1))
         fac = fac*fac
         a(i) = a(i)*fac
         a(l) = a(l)*fac
    3 continue 
      return
      end 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE DFT_TS(MDIM,MFREQ,MSAMP,MREC, NDIM,NFREQ,NSAMP,NREC, 
     ;                  NSRC,IDIR, DT,STARTT, FREQ,TS, TFORM) 
!
!     Fourier transforms the time domain series
!
!     INPUT      MEANING
!     -----      -------
!     DT         sampling period (s)
!     FREQ       frequency list (Hz)
!     IDIR       -1 Time domain -> frequency standard definition
!     MDIM       leading dimension 
!     MFREQ      leading dimension
!     MREC       leading dimension
!     MSAMP      leading dimension
!     NDIM       number of components
!     NFREQ      number of frequencies
!     NREC       number of receivers
!     NSAMP      number of samples
!     NSRC       number of sources
!     STARTT     start time, probably 0 
!     TS         time series to DFT 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     TFORM      frequency domain components of ts 
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FREQ(NFREQ) 
      REAL*4, INTENT(IN) :: TS(MDIM,MSAMP,MREC,*), DT, STARTT
      INTEGER*4, INTENT(IN) :: MDIM,MFREQ,MSAMP,MREC, NDIM,NFREQ,NREC,
     ;                         NSAMP,NSRC,IDIR 
      COMPLEX*8, INTENT(OUT) :: TFORM(MDIM,MFREQ,MREC,*)
!.... local variables
      COMPLEX*8 CZERO, CFACT, CTS 
      REAL*4 TIME, DIR, OMEGA, ARG  
      INTEGER*4 IFREQ, ISRC, IREC, K, I
      PARAMETER(CZERO = CMPLX(0.0,0.0)) 
      REAL*8 PI /3.14159265358979323846D0/
!
!----------------------------------------------------------------------!
!
      DIR = FLOAT(IDIR) 
      DO 1 IFREQ=1,NFREQ
         OMEGA = SNGL( 2.D0*PI*FREQ(IFREQ) )
         DO 2 ISRC=1,NSRC
            DO 3 IREC=1,NREC
               TFORM(1:NDIM,IFREQ,IREC,ISRC) = CZERO  
               DO 4 K=1,NSAMP 
                  TIME = FLOAT(K - 1)*DT + STARTT 
                  ARG = DIR*OMEGA*TIME 
                  CFACT = CMPLX(CEXP(CMPLX(0.0,ARG)))
                  DO 5 I=1,NDIM 
                     CTS = TS(I,K,IREC,ISRC) 
                     TFORM(I,IFREQ,IREC,ISRC) = TFORM(I,IFREQ,IREC,ISRC)
     ;                                        + CFACT*CTS 
    5             CONTINUE
    4          CONTINUE 
    3       CONTINUE  
    2    CONTINUE
    1 CONTINUE
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SPMASK_V2(JOB,X,LX,NPCT,ITYPE,TSV,I1,I2, IERROR)
      IMPLICIT NONE
      REAL*4, INTENT(INOUT) :: X(0:LX)
      INTEGER*4, INTENT(IN) :: JOB, LX, NPCT, ITYPE, I1, I2
      REAL*4, INTENT(OUT) :: TSV 
      INTEGER*4, INTENT(OUT) :: IERROR 
!.... local variables
      REAL*8 DPWNDO
      REAL*4 W 
      INTEGER*4 LX2, IX2, K
!
!----------------------------------------------------------------------!
!
      IF (JOB.EQ.1) THEN
         CALL SPMASK(X,LX,ITYPE,TSV, IERROR)
      ELSE
         LX2 = 0 
         DO 1 K=0,LX
            IF (K.LT.I1) THEN 
               LX2 = LX2 + 1   
            ELSEIF (K.GT.I2) THEN
               LX2 = LX2 + 1 
            ENDIF
    1    CONTINUE
         IX2 =-1 
         DO 2 K=0,LX
            IF (K.LT.I1) THEN
               IX2 = IX2 + 1 
               W = SNGL( DPWNDO(ITYPE,LX2,NPCT,IX2) )
            ELSEIF (K.GT.I2) THEN
               IX2 = IX2 + 1 
               W = SNGL( DPWNDO(ITYPE,LX2,NPCT,IX2) )
            ELSE
               W = 1.0
            ENDIF
            X(K) = X(K)*W
    2    CONTINUE
      ENDIF
      RETURN
      END 

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SPMASK(X,LX,ITYPE,TSV, IERROR)
! 
!     APPLIES A DATA WINDOW TO THE DATA VECTOR X(0:LX).
!     ITYPE = 1(RECTANGULAR), 2(TAPERED RECTANGULAR), 3(TRIANGULAR)
!             4(HANNING), 5(HAMMING), 6(BLACKMAN)
!             (NOTE: TAPERED RECTANGULAR HAS COSINE-TAPERED 10% ENDS)
!     TSV = SUM OF SQUARED WINDOW VALUES
!     IERROR = 0 IF NO ERROR, 1 IF ITYPE OUT OF RANGE
      IMPLICIT NONE
      REAL*4, INTENT(INOUT) :: X(0:LX)
      INTEGER*4, INTENT(IN) :: LX, ITYPE
      REAL*4, INTENT(OUT) :: TSV
      INTEGER*4, INTENT(OUT) :: IERROR
      REAL*8 DPWNDO 
      REAL*4 W 
      INTEGER*4 K 
      IERROR = 1
      IF (ITYPE.LT.1 .OR. ITYPE.GT.6) RETURN
      TSV = 0.0
      DO 1 K=0,LX
         W = SNGL( DPWNDO(ITYPE,LX+1,10,K) )
         X(K) = X(K)*W
         TSV = TSV + W**2 
    1 CONTINUE
      IERROR = 0   
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      REAL*8 FUNCTION DPWNDO(ITYPE,N,NPCT,K)
!     
!     THIS GENERATES A FSINGLE SAMPLE OF A DATA WINDOW.
!     ITYPE = 1 (RECTANGULAR), 2 (TAPERED RECTANGULAR), 3 (TRIANGULAR),
!             4 (HANNING), 5 (HAMMING), OR 6 (BLACKMAN)
!             (NOTE: TAPERED RECTANGULAR HAS COSINE-TAPERED 10% ENDS)
!     N = SIZE (TOTAL NO. SAMPLES) OF WINDOW
!     K = SAMPLE NUMBER WITHIN WINDOW, FROM 0 THROUGH N-1
!         (IF K IS OUTSIDE THIS RANGE, SPWNDO IS SET TO 0.)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4, INTENT(IN) :: ITYPE,N,K,NPCT 
      PI = 4.D0*DATAN(1.D0)
      DPWNDO = 0.D0 
      IF (ITYPE.LT.1 .OR.ITYPE.GT.6) RETURN
      IF (ITYPE.EQ.2 .AND. NPCT.LT.0 .OR. NPCT.GT.100) THEN
         WRITE(*,*) 'dpwndo: npct must be between [0,100]',NPCT
         RETURN
      ENDIF
      IF (K.LT.0 .OR. K.GE.N) RETURN
      DPWNDO = 1.D0
      GOTO (1,2,3,4,5,6), ITYPE
    1 RETURN
    2 L = (N - 2)/NPCT !10
      IF (K.LE.L) DPWNDO = 0.5D0
     ;                   *(1.D0 - DCOS(DFLOAT(K)*PI/(DFLOAT(L + 1))))
      IF (K.GT.N-L-2) DPWNDO = 0.5D0*(1.D0
     ;                       - DCOS(DFLOAT(N - K - 1)*PI/DFLOAT(L + 1)))
      RETURN
    3 DPWNDO = 1.D0 - DABS(1.D0 - 2.D0*DFLOAT(K)/(DFLOAT(N) - 1.D0))
      RETURN
    4 DPWNDO = 0.5D0*(1.D0 - DCOS(2.D0*DFLOAT(K)*PI/(DFLOAT(N) - 1.D0)))
      RETURN
!.... alpha and beta optimized from wikipedia
    5 DPWNDO = 0.53836D0
     ;       - 0.46164D0*DCOS(2.D0*DFLOAT(K)*PI/(DFLOAT(N) - 1.D0))
c   5 DPWNDO = 0.54D0 
c    ;       - 0.46D0*DCOS(2.D0*DFLOAT(K)*PI/(DFLOAT(N) - 1.D0))
      RETURN
    6 DPWNDO = 0.42D0 - 0.5D0*DCOS(2.D0*DFLOAT(K)*PI/(DFLOAT(N) - 1.D0))
     ;       + 0.08D0*DCOS(4.D0*DFLOAT(K)*PI/(DFLOAT(N) - 1.D0))
      RETURN
      END

