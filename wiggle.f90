!
!     Fourier transforms the estimate wavefields for viewing with Gphase
      CHARACTER(80) PROJNM, PROJNM_WRK
      REAL*8 DT8,START8
      COMPLEX*8, ALLOCATABLE :: EST(:,:,:,:)
      REAL*4, ALLOCATABLE :: WIGGLE(:,:,:,:), FREQ(:)
      CHARACTER(1) C0, CIN
      LOGICAL*4 LPRATT, LEX 
      PARAMETER(LPRATT = .FALSE.) !do not use the pratt version of the F.T.
!
!----------------------------------------------------------------------------------------!
!
!.... read project name and options 
      PROJNM(1:80) = ' '
      WRITE(*,*) 'xwiggle25: Enter project name'
      READ(*,'(A)') PROJNM 
      PROJNM = ADJUSTL(PROJNM)
!
!.... determine if this is a mixed source simulation
      NTYPE = 1
      PROJNM_WRK(1:80) = ' '
      PROJNM_WRK = TRIM(PROJNM)//'.freqs'
      PROJNM_WRK = ADJUSTL(PROJNM_WRK)
      INQUIRE(FILE=TRIM(PROJNM_WRK),EXIST=LEX) 
      IF (LEX) THEN
         OPEN(UNIT=55,FILE=TRIM(PROJNM)//'.freqs')
         C0 = ' '
         READ(55,*,END=55)  
         READ(55,*,END=55) NFREQ_IN
         DO 1 IFREQ=1,NFREQ_IN
            READ(55,*,END=55) CIN
            IF (C0 == ' ' ) THEN
               C0 = CIN
            ELSE
               IF (C0 /= CIN) THEN
                  NTYPE = 2 
                  GOTO 55
               ENDIF
            ENDIF 
    1    CONTINUE
   55    CONTINUE 
         CLOSE(55) 
         IF (NTYPE == 1) THEN
            IF (CIN == 'S') THEN
               WRITE(*,*) 'xwiggle25: Will generate synthetics for surface waves...' 
            ELSE
               WRITE(*,*) 'xwiggle25: Will generate synthetics for body waves...'
            ENDIF
         ELSE
            WRITE(*,*) 'xwiggle25: Will generate synthetics for surface and body waves...'
         ENDIF
      ENDIF
!
!.... get number of samples
      WRITE(*,*) 'xwiggle25: Enter (1) to read the number of samples from a source file'
      !READ *, IEST
      iest = 1
      DO 100 ITYPE=1,NTYPE
         PROJNM_WRK(1:80) = ' '
         IF (NTYPE == 1) THEN
            PROJNM_WRK = PROJNM
         ELSE
            IF (ITYPE == 1) THEN
               PROJNM_WRK = TRIM(ADJUSTL(PROJNM))//'_srf'
            ELSE
               PROJNM_WRK = TRIM(ADJUSTL(PROJNM))//'_bdy'
            ENDIF
         ENDIF
         PROJNM_WRK = ADJUSTL(PROJNM_WRK)
         IF (IEST == 1) THEN
            WRITE(*,*) 'xwiggle25: Reading .srcts file'
            CALL RDSRCT_HD(PROJNM_WRK, NSRCIN,NSAMP,DT8,START8,IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xwiggle25: Enter number of samples'
               READ *, NSAMP
               WRITE(*,*) 'xwiggle25: Enter sampling interval (s)'
               READ *, DT
               STARTT = 0.0
            ELSE
               STARTT = SNGL(START8)
               DT = SNGL(DT8) 
            ENDIF
         ELSE
            WRITE(*,*) 'xwiggle25: Enter number of samples'
            READ *, NSAMP
            WRITE(*,*) 'xwiggle25: Enter sampling interval (s)'
            READ *, DT
            STARTT = 0.0 
         ENDIF
!
!....... read estimate file
         WRITE(*,*) 'xwiggle25: Reading estimate file...'
         CALL RDEEST_HD(PROJNM_WRK, NDIM,NFREQ,NREC,NSRC, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xwiggle25: Error could not read estimate file header'
            STOP
         ENDIF
         ALLOCATE(FREQ(NFREQ)) 
         CALL RDEEST_FR(PROJNM_WRK, NFREQ,FREQ,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xwiggle25: Error could not read frequencies'
            STOP
         ENDIF
         ALLOCATE(EST(NDIM,NFREQ,NREC,NSRC)) 
         CALL RDEEST25(PROJNM_WRK, NFREQ,NREC,NFREQ, NDIM,NREC,NSRC,          &
                       FREQ, EST(1,:,:,:),EST(2,:,:,:),EST(3,:,:,:),IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xwiggle25: Error could not read responses'
            STOP
         ENDIF 
!
!....... dft
         WRITE(*,*) 'xwiggle25: Performing dft...'
         ALLOCATE(WIGGLE(NDIM,NSAMP,NREC,NSRC)) 
         WIGGLE(:,:,:,:) = 0.0
         IDIR = 1 !inverse transform
         IF (LPRATT) IDIR =-1
         CALL DFTWIG4(NDIM,NSAMP,NFREQ,NREC, NSAMP,NDIM,NFREQ,NREC,NSRC,IDIR,DT, &
                      FREQ,EST, WIGGLE)
         WRITE(*,*) 'xwiggle25: Writing seismograms...'
         CALL WTSEISM(PROJNM_WRK, NDIM,NSAMP,NREC, NDIM,NSAMP,NREC,NSRC,    &
                      DT,STARTT, WIGGLE, IERR)
         IF (IERR /= 0) WRITE(*,*) 'xwiggle25: Error writing .pest file'
         DEALLOCATE(FREQ)
         DEALLOCATE(EST)
         DEALLOCATE(WIGGLE) 
  100 CONTINUE !Loop on types
      STOP
      END  
