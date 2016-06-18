      SUBROUTINE DFTWIG4(MDIM,MSAMP,MFREQ,MREC, NSAMP,NDIM,NFREQ,NREC,NSRC,IDIR,DT, &
                         FREQS,EST, WIGGLE)
      COMPLEX*8, INTENT(IN) :: EST(MDIM,MFREQ,MREC,*) 
      REAL*4, INTENT(IN) :: FREQS(NFREQ), DT
      INTEGER*4, INTENT(IN) :: MDIM,MSAMP,MFREQ,MREC, NSAMP,NDIM,NFREQ,NREC,NSRC,IDIR
      REAL*4, INTENT(OUT) :: WIGGLE(MDIM,MSAMP,MREC,*)
      COMPLEX*8 CFACT
      REAL*4 SCLR,SSAMP,EPS
      PARAMETER(EPS = 1.E-30) 
      REAL*4 PI /3.14159265358979323846E0/ 
! 
!----------------------------------------------------------------------------------------!
! 
!.... DFT each sample
      SCLR = REAL(IDIR) 
      SSAMP = REAL(NSAMP)
      DO 1 K=1,NSAMP 
         TIME = REAL(K - 1)*DT
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM 
                  WIGGLE(I,K,IREC,ISRC) = 0.0
                  DO 5 IFREQ=1,NFREQ 
                     OMEGA = 2.0*PI*FREQS(IFREQ)
                     ARG = OMEGA*TIME 
                     CFACT = CEXP(CMPLX(0.0,SCLR*ARG)) 
                     IF (ABS(OMEGA) <= EPS) THEN !0 frequency 
                        WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC) &
                                              + REAL(EST(I,IFREQ,IREC,ISRC))
                     ELSEIF (ABS(OMEGA - PI/DT) <= EPS) THEN !periodicity in 0 frequency
                        WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC) &
                                              + REAL(EST(I,IFREQ,IREC,ISRC)*CFACT)
                     ELSE !typical invsere-dft 
                        WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC) &
                                              + 2.0*REAL(EST(I,IFREQ,IREC,ISRC)*CFACT)
                     ENDIF
    5             CONTINUE !loop on frequencies 
                  WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC)/SSAMP
    4          CONTINUE !loop on spatial dimensions 
    3       CONTINUE !loop on sources
    2    CONTINUE !loop on recievers
    1 CONTINUE !loop on samples
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DFTWIG8(MDIM,MSAMP,MFREQ,MREC, NSAMP,NDIM,NFREQ,NREC,NSRC,IDIR,DT, &
                         FREQS,EST, WIGGLE)
      COMPLEX*16, INTENT(IN) :: EST(MDIM,MFREQ,MREC,*) 
      REAL*8, INTENT(IN) :: FREQS(NFREQ),DT
      INTEGER*4, INTENT(IN) :: MDIM,MSAMP,MFREQ,MREC, NSAMP,NDIM,NFREQ,NREC,NSRC,IDIR
      REAL*8, INTENT(OUT) :: WIGGLE(MDIM,MSAMP,MREC,*)
      COMPLEX*16 CFACT
      REAL*8 SCLR,DSAMP,EPS,OMEGA,ARG,TIME
      PARAMETER(EPS = 1.D-30)
      REAL*8 PI /3.14159265358979323846D0/ 
! 
!----------------------------------------------------------------------------------------!
! 
!.... DFT each sample
      SCLR = DFLOAT(IDIR) 
      DSAMP = DFLOAT(NSAMP)  
      DO 1 K=1,NSAMP 
         TIME = DFLOAT(K - 1)*DT
         DO 2 IREC=1,NREC 
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  WIGGLE(I,K,IREC,ISRC) = 0.D0 
                  DO 5 IFREQ=1,NFREQ
                     OMEGA = 2.D0*PI*FREQS(IFREQ)
                     ARG = OMEGA*TIME 
                     CFACT = CDEXP(DCMPLX(0.D0,SCLR*ARG)) 
                     IF (ABS(OMEGA) <= EPS) THEN !0 frequency 
                        WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC) &
                                              + DREAL(EST(I,IFREQ,IREC,ISRC))
                     ELSEIF (ABS(OMEGA - PI/DT) <= EPS) THEN !periodicity in 0 frequency
                        WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC) &
                                              + DREAL(EST(I,IFREQ,IREC,ISRC)*CFACT)
                     ELSE !typical invsere-dft 
                        WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC) &
                                              + 2.D0*DREAL(EST(I,IFREQ,IREC,ISRC)*CFACT)
                     ENDIF
    5             CONTINUE !Loop on frequencies
                  WIGGLE(I,K,IREC,ISRC) = WIGGLE(I,K,IREC,ISRC)/DSAMP
    4          CONTINUE !loop on spatial dimensions
    3       CONTINUE !loop on sources
    2    CONTINUE !loop on receivers
    1 CONTINUE !loop on smaples 
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DFTF4(NFREQ,NSAM,IDIR, DT,STARTT,FREQ,TRACE, TFORM)
! 
!     Discrete Fourier transform at selected frequencies 
      REAL*4 FREQ(NFREQ),TRACE(NSAM),DT,STARTT
      INTEGER*4 NFREQ,NSAM,IDIR
      COMPLEX*8 TFORM(NFREQ)
      COMPLEX*8 CZERO
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      COMPLEX*8 CFACT
      REAL*4 ARG, DIR, OMEGA, TIME
      INTEGER*4 IFREQ,J
      REAL*4 PI /3.14159265358979323846E0/ 
! 
!----------------------------------------------------------------------------------------!
! 
!.... discrete Fourier transform 
      DIR = REAL(IDIR) !transform direction
      DO 1 IFREQ=1,NFREQ
         OMEGA = 2.0*PI*FREQ(IFREQ)
         TFORM(IFREQ) = CZERO
         DO 2 J=1,NSAM
            TIME = REAL(J - 1)*DT + STARTT
            ARG = DIR*OMEGA*TIME
            CFACT = EXP(CMPLX(0.0,ARG))
            TFORM(IFREQ) = TFORM(IFREQ) + TRACE(J)*CFACT
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DFTF8(NFREQ,NSAM,IDIR, DT,STARTT,FREQ,TRACE, TFORM)
! 
!     Discrete Fourier transform at selected frequencies 
      REAL*8 FREQ(NFREQ),TRACE(NSAM),DT,STARTT
      INTEGER*4 NFREQ,NSAM,IDIR
      COMPLEX*16 TFORM(NFREQ)
      COMPLEX*16 CZERO
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      COMPLEX*16 CFACT
      REAL*8 ARG, DIR, OMEGA, TIME, PI
      INTEGER*4 IFREQ,J
      PARAMETER(PI = 3.1415926535897931D0)
! 
!----------------------------------------------------------------------------------------!
! 
!.... discrete Fourier transform 
      DIR = DFLOAT(IDIR) !transform direction
      DO 1 IFREQ=1,NFREQ
         OMEGA = 2.D0*PI*FREQ(IFREQ)
         TFORM(IFREQ) = CZERO
         DO 2 J=1,NSAM
            TIME = DFLOAT(J - 1)*DT + STARTT
            ARG = DIR*OMEGA*TIME
            CFACT = CDEXP(DCMPLX(0.D0,ARG))
            TFORM(IFREQ) = TFORM(IFREQ) + TRACE(J)*CFACT
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DFT_TF4(MDIM,MSAMP,MFREQ,MREC,  &
                         NDIM,NSAMP,NFREQ,NREC,NSRC, IDIR, & 
                         DT,STARTT, FREQ,TS, EST)
!
!     Performs the time to frequency DFT of a signal  
!
!     INPUT      MEANING
!     -----      -------
!     DT         sampling interval (s)
!     FREQ       frequency list to transform at (Hz)
!     IDIR       convention on FT sign; probably -1 
!     MDIM       leading dimension
!     MFREQ      leading dimension
!     MREC       leading dimension
!     MSAMP      leading dimension
!     NDIM       number of components
!     NFREQ      number of frequencies
!     NREC       number of receivers
!     NSAMP      number of samples
!     NSRC       number of sources
!     STARTT     start time (s)
!     TS         time series to transform
!
!     OUTPUT     MEANING
!     ------     -------
!     EST        frequency domain version of TS 
!
!.... variable declarations
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: TS(MDIM,MSAMP,MREC,*), FREQ(NFREQ), DT,STARTT
      INTEGER*4, INTENT(IN) :: MDIM,MSAMP,MFREQ,MREC, NDIM,NSAMP,NFREQ,NREC,NSRC, IDIR
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MFREQ,MREC,*)
!.... local variables
      COMPLEX*8 CFACT, CTRACE, CZERO 
      REAL*4 OMEGA, DIR, ARG, TWOPI, TIME 
      INTEGER*4 ISRC,IREC,IFREQ,K,I
      PARAMETER(TWOPI = 6.2831853071795862E0) !2*pi 
      PARAMETER(CZERO = CMPLX(0.0,0.0))
!
!----------------------------------------------------------------------------------------!
!
      DIR = FLOAT(IDIR)
      DO 1 ISRC=1,NSRC
         DO 2 IREC=1,NREC
            DO 3 IFREQ=1,NFREQ
               EST(1:NDIM,IFREQ,IREC,ISRC) = CZERO
               OMEGA = TWOPI*FREQ(IFREQ)
               DO 4 K=1,NSAMP
                  TIME = FLOAT(K - 1)*DT + STARTT
                  ARG = DIR*OMEGA*TIME
                  CFACT = CMPLX(CEXP(CMPLX(0.0,ARG))) 
                  DO 5 I=1,NDIM
                     CTRACE = CMPLX(TS(I,K,IREC,ISRC),0.0)
                     EST(I,IFREQ,IREC,ISRC) = EST(I,IFREQ,IREC,ISRC) + CFACT*CTRACE
    5             CONTINUE
    4          CONTINUE
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE 
      RETURN
      END
