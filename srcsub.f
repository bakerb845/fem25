
      SUBROUTINE SRCSUB_SB(PROJNM,MFREQ, NFREQ,NSRC,
     ;                     NFREQ_SRF,NFREQ_BDY,NSRC_SRF,NSRC_BDY, 
     ;                     CFTYPE,SRCTYP, FREQ, SOURCE,IERR) 
!
!     Reads the time domain STFs for surface and body waves and 
!     transforms to the frequency domain at the desired frequencies.  
!     The file extensions that must exist are projnm_srf.srcts if 
!     we are modeling surface waves and projnm_bdy.srcts if we are 
!     modeling body waves.  -  B. Baker May 2013
!
!     INPUT      MEANING
!     -----      -------
!     CFTYPE     frequency type 'B' -> body, 'S' -> surface wave
!     FREQ       frequency list
!     MFREQ      leading dimension for SOURCE
!     NFREQ      number of frequencies 
!     NFREQ_BDY  number of body wave frequencies
!     NFREQ_SRC  number of surface wave frequencies
!     NSRC       number of sources
!     NSRC_BDY   number of body wave sources
!     NSRC_SRF   number of surface wave sources
!     PROJNM     project name
!     SRCTYP     body wave or surface source
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     SOURCE     source time functions in frequency domain
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: PROJNM
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC) 
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ, NFREQ, NSRC, NFREQ_SRF, 
     ;                         NFREQ_BDY, NSRC_SRF, NSRC_BDY 
      COMPLEX*8, INTENT(OUT) :: SOURCE(MFREQ,NSRC) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*8, ALLOCATABLE :: SRCWRK(:,:) 
      REAL*8, ALLOCATABLE :: FREQW(:)  
      CHARACTER(80) FLNAME
      COMPLEX*8 CONE 
      INTEGER*4 IFREQ, ISRC, IFREQ_SRF, IFREQ_BDY, ISRC_SRF, ISRC_BDY 
      PARAMETER(CONE = CMPLX(1.0,0.0)) 
!
!----------------------------------------------------------------------!
!
      IERR = 0
      IF (NSRC_SRF.GT.0) THEN
         WRITE(*,*) 'srcsub_sb: Processing surface wave sources...'
         ALLOCATE(FREQW(NFREQ_SRF))
         IFREQ_SRF = 0
         DO 1 IFREQ=1,NFREQ
            IF (CFTYPE(IFREQ).EQ.'S') THEN
               IFREQ_SRF = IFREQ_SRF + 1
               FREQW(IFREQ_SRF) = FREQ(IFREQ)
            ENDIF
    1    CONTINUE
         IF (IFREQ_SRF.NE.NFREQ_SRF) THEN
            WRITE(*,*) 'srcsub_sb: Serious error ifreq_srf/=nfreq_srf!'
            IERR = 1 
            RETURN
         ENDIF
         ALLOCATE(SRCWRK(NFREQ_SRF,NSRC_SRF)) 
         FLNAME(1:80) = ' '
         FLNAME = TRIM(ADJUSTL(PROJNM))//'_srf'
         FLNAME = ADJUSTL(FLNAME)
         CALL SRCSUB(FLNAME,NFREQ_SRF,NSRC_SRF, NFREQ_SRF,FREQW, 
     ;               SRCWRK,IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'srcsub_sb: Setting surface wave STFs to unity!'
            IERR =-1
         ENDIF
         IFREQ_SRF = 0
         DO 2 IFREQ=1,NFREQ
            IF (CFTYPE(IFREQ).EQ.'S') THEN
               IFREQ_SRF = IFREQ_SRF + 1
               ISRC_SRF = 0 
               DO 3 ISRC=1,NSRC
                  IF (SRCTYP(ISRC)(1:1).EQ.'S') THEN
                     ISRC_SRF = ISRC_SRF + 1
                     IF (IERR.EQ.0) THEN
                        SOURCE(IFREQ,ISRC) = SRCWRK(IFREQ_SRF,ISRC_SRF)
                     ELSE
                        SOURCE(IFREQ,ISRC) = CONE
                     ENDIF
                  ENDIF
    3          CONTINUE  
            ENDIF 
    2    CONTINUE !loop on frequencies
         IERR = 0
         DEALLOCATE(FREQW)
         DEALLOCATE(SRCWRK)
      ENDIF
      IF (NSRC_BDY.GT.0) THEN
         WRITE(*,*) 'srcsub_sb: Processing body wave sources...'
         ALLOCATE(FREQW(NFREQ_BDY))
         IFREQ_BDY = 0
         DO 11 IFREQ=1,NFREQ
            IF (CFTYPE(IFREQ).EQ.'B') THEN
               IFREQ_BDY = IFREQ_BDY + 1
               FREQW(IFREQ_BDY) = FREQ(IFREQ)
            ENDIF
   11    CONTINUE
         IF (IFREQ_BDY.NE.NFREQ_BDY) THEN
            WRITE(*,*) 'srcsub_sb: Serious error ifreq_bdy/=nfreq_bdy!'
            IERR = 1
            RETURN
         ENDIF
         ALLOCATE(SRCWRK(NFREQ_BDY,NSRC_BDY))
         FLNAME(1:80) = ' '
         FLNAME = TRIM(ADJUSTL(PROJNM))//'_bdy'
         FLNAME = ADJUSTL(FLNAME)
         CALL SRCSUB(FLNAME,NFREQ_BDY,NSRC_BDY, NFREQ_BDY,FREQW, 
     ;               SRCWRK,IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'srcsub_sb: Setting body wave STFs to unity!'
            IERR =-1
         ENDIF
         IFREQ_BDY = 0 
         DO 12 IFREQ=1,NFREQ
            IF (CFTYPE(IFREQ).EQ.'B') THEN
               IFREQ_BDY = IFREQ_BDY + 1
               ISRC_BDY = 0
               DO 13 ISRC=1,NSRC
                  IF (SRCTYP(ISRC)(1:1).EQ.'P') THEN
                     ISRC_BDY = ISRC_BDY + 1
                     IF (IERR.EQ.0) THEN
                        SOURCE(IFREQ,ISRC) = SRCWRK(IFREQ_BDY,ISRC_BDY)
                     ELSE
                        SOURCE(IFREQ,ISRC) = CONE
                     ENDIF
                  ENDIF
   13          CONTINUE  
            ENDIF 
   12    CONTINUE !loop on frequencies
         IERR = 0
         DEALLOCATE(FREQW)
         DEALLOCATE(SRCWRK)
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SRCSUB(PROJNM,MFREQ,NSRC, NFREQ,FREQ, STF,IERR) 
! 
!     Reads and Fourier transforms the source file
! 
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       frequency list
!     MFREQ      leading dimension for STF
!     NFREQ      number of frequencies
!     NSRC       number of sources
!     PROJNM     project name 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     STF        source time function 
! 
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM 
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ,NSRC,NFREQ
      COMPLEX*8, INTENT(OUT) :: STF(MFREQ,NSRC)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local varaibles
      REAL*8, ALLOCATABLE :: TRACE(:,:) 
      COMPLEX*8 CFACT
      REAL*8 DT,STARTT,OMEGA,TIME,DIR,ARG 
      PARAMETER(IDIR =-1) !Time domain -> Frequency domain
      REAL*8 PI /3.14159265358979323846D0/
! 
!----------------------------------------------------------------------!
! 
!.... read source header for DT 
      CALL RDSRCT_HD(PROJNM, NSRCIN,NSAMP,DT,STARTT,IERR)
      IF (IERR.LT.0) THEN
         WRITE(*,*) 'srcsub: Warning setting Greens functions to unity!'
         DO 1 IFREQ=1,NFREQ
            DO 2 ISRC=1,NSRC
               STF(IFREQ,ISRC) = CMPLX(1.0,0.0)
    2       CONTINUE
    1    CONTINUE
         RETURN
      ENDIF 
      IF (IERR.NE.0) RETURN
      ALLOCATE(TRACE(NSAMP,NSRC)) 
!.... read source file and pack trace
      CALL RDSRCT(NSAMP,PROJNM, NSRC,TRACE,IERR) 
      IF (IERR.NE.0) RETURN
 
!.... discrete fourier transform
      DIR = DFLOAT(IDIR)
      DO 3 IFREQ=1,NFREQ 
         OMEGA = 2.D0*PI*FREQ(IFREQ)  
         DO 4 ISRC=1,NSRC
            STF(IFREQ,ISRC) = CMPLX(0.0,0.0) 
            DO 5 K=1,NSAMP 
               TIME = DFLOAT(K - 1)*DT + STARTT
               ARG = DIR*OMEGA*TIME
               CFACT = CMPLX(CDEXP(DCMPLX(0.D0,ARG)))
               STF(IFREQ,ISRC) = STF(IFREQ,ISRC) 
     ;                         + CFACT*CMPLX(SNGL(TRACE(K,ISRC)),0.0)
    5       CONTINUE
    4    CONTINUE
    3 CONTINUE 
      DEALLOCATE(TRACE)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDSRCT_HD(PROJNM, NSRC,NSAMP,DT,STARTT,IERR)
! 
!     Reads the source header 
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!     DT         time spacing in seconds
!     NSAMP      number of samples in trace
!     NSRC       number of sources in srcfl
!     STARTT     start time in seconds (likely 0.)
!
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM 
      REAL*8 DT,STARTT
!.... local variables
      CHARACTER(80) SRCFL
      CHARACTER(4) MODEL(4)    
      REAL*4 UNPACKR4
      LOGICAL*4 LSWAP, LEX
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(IUNIT = 73, LENHEAD = 4) 
! 
!----------------------------------------------------------------------!
!   
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
      IERR = 0 
      SRCFL(1:80) = ' '
      SRCFL = TRIM(PROJNM)//'.srcts'
      SRCFL = ADJUSTL(SRCFL)
      INQUIRE(FILE=TRIM(SRCFL),EXIST=LEX,SIZE=NBYTES) 
      IF (.NOT.LEX) THEN 
         WRITE(*,*) 'rdsrct_hd: Error cannot locate source time series'
         IERR =-1
         RETURN
      ENDIF
!     CALL STAT(TRIM(SRCFL),INFO)
!     NBYTES = INFO(8) 
      LRECS = LENHEAD
      NBYTES = 4*LRECS
      OPEN(UNIT=IUNIT,FILE=TRIM(SRCFL),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      READ(IUNIT,REC=1,IOSTAT=IERR) (MODEL(J),J=1,LENHEAD) 
      CLOSE(IUNIT)
      IF (IERR.NE.0) THEN 
         WRITE(*,*) 'rdsrc_hd: Error reading source file'
         RETURN
      ENDIF 
! 
!.... read the header
      NSRC   = UNPACKI4(LSWAP,MODEL(1))
      NSAMP  = UNPACKI4(LSWAP,MODEL(2))
      DT4    = UNPACKR4(LSWAP,MODEL(3))
      START4 = UNPACKR4(LSWAP,MODEL(4))
      DT = DBLE(DT4)
      STARTT = DBLE(START4)
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDSRCT(MSAMP,PROJNM, NSRC,TRACE,IERR) 
! 
!     Reads the source file
! 
!     INPUT      MEANING
!     -----      ------- 
!     MSAMP      max number of samples
!     NSRC       number of sources
!     SRCFL      source filename 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     TRACE      time domain source signatures
! 
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM 
      INTEGER*4, INTENT(IN) :: MSAMP,NSRC
      REAL*8 TRACE(MSAMP,NSRC)
!.... local variables
      CHARACTER(80) SRCFL
      CHARACTER(4), ALLOCATABLE :: MODEL(:)
      REAL*4 UNPACKR4, TEMP
      INTEGER*4 UNPACKI4, ENDIAN
      LOGICAL*4 LSWAP, LEX
      PARAMETER(IUNIT = 73)
! 
!----------------------------------------------------------------------!
! 
!.... determine endianness
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... read source file
      IERR = 0
      SRCFL(1:80) = ' '
      SRCFL = TRIM(PROJNM)//'.srcts'
      INQUIRE(FILE=TRIM(SRCFL),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdsrct: source file does not exist',TRIM(SRCFL)
         IERR = 1
         RETURN
      ENDIF
!     CALL STAT(TRIM(SRCFL),INFO)
!     NBYTES = INFO(8)
      LRECS = NBYTES/4
      ALLOCATE(MODEL(LRECS))
      OPEN(UNIT=IUNIT,FILE=TRIM(SRCFL),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      READ(IUNIT,REC=1,IOSTAT=IERR) (MODEL(J),J=1,LRECS)
      CLOSE(IUNIT)
      IF (IERR.NE.0) THEN
         WRITE(*,*) 'rdsrct: Error reading source file'
         DEALLOCATE(MODEL)
         RETURN
      ENDIF
! 
!.... read the header
      NSRCIN  = UNPACKI4(LSWAP,MODEL(1))
      NSAMPIN = UNPACKI4(LSWAP,MODEL(2))
      DTIN    = UNPACKR4(LSWAP,MODEL(3))
      STARTIN = UNPACKR4(LSWAP,MODEL(4))
! 
!.... unpack source
      NSL = MIN(NSRC,NSRCIN)
      MYLOC = 4
      DO 1 I=1,NSL
         DO 2 J=1,NSAMPIN
            MYLOC = MYLOC + 1
            TEMP = UNPACKR4(LSWAP,MODEL(MYLOC))
            TRACE(J,I) = DBLE(TEMP)
    2    CONTINUE
    1 CONTINUE
      IF (ALLOCATED(MODEL)) DEALLOCATE(MODEL)
! 
!.... possibly copy source
      IF (NSRCIN.LT.NSRC) THEN
         WRITE(*,*) 'rdsrct: Copying first source to others...'
         DO 3 J=1,NSAMPIN
            DO 4 I=NSRCIN+1,NSRC
               TRACE(J,I) = TRACE(J,1)
    4       CONTINUE
    3    CONTINUE
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE MODSRC(MFREQ, NFREQ,NSRC,IRESTP, SOURCE) 
!
!     In inversion may only want the source amplitude or source phase
!
!     INPUT      MEANING
!     -----      -------
!     IRESTP     =1 -> Phase only 
!                =2 -> Amplitude only 
!                =3 -> phase and amplitude (default)
!     MFREQ      leading dimension 
!     NFREQ      number of frequencies
!     NSRC       number of sources
!
!     OUTPUT     MEANING
!     ------     -------
!     SOURCE     modified source time function 
!
!.... variable declarations
      COMPLEX*8, INTENT(INOUT) :: SOURCE(MFREQ,*) 
      INTEGER*4, INTENT(IN) :: MFREQ, NFREQ, NSRC, IRESTP 
!.... local variables
      COMPLEX*8 CPHM2CM 
      REAL*4 SPHASE, PHASE, RMAG 
      INTEGER*4 IFREQ, ISRC 
!
!----------------------------------------------------------------------!
!
      IF (IRESTP.NE.1 .AND .IRESTP.NE.2) RETURN !Nothing to do
      DO 1 IFREQ=1,NFREQ
         DO 2 ISRC=1,NSRC
            RMAG  =   CABS(SOURCE(IFREQ,ISRC))  
            PHASE = SPHASE(SOURCE(IFREQ,ISRC))  
            IF (RMAG.GT.0.0) THEN
               IF (IRESTP.EQ.1) RMAG  = 1.0 !phase only inversion
               IF (IRESTP.EQ.2) PHASE = 0.0 !amplitude only inversion
               SOURCE(IFREQ,ISRC) = CPHM2CM(RMAG,PHASE)
            ELSE
               SOURCE(IFREQ,ISRC) = CPHM2CM(0.0,0.0)
            ENDIF
    2    CONTINUE
    1 CONTINUE  
      RETURN
      END
