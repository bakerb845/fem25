      SUBROUTINE RD_OBS_EST(PROJNM, MDIM,MFREQL,MREC,  
     ;                      NDIM,NFREQL,NREC,NSRC,  
     ;                      IRESTP,IMODSRC,IMODREC,  LUNWRAP, 
     ;                      AZMOD,FREQL,  LUPDREC,LUPDSRC, SOURCE,RECV,
     ;                      WGHT,OBS,EST,IERR)
!
!     Reads the observations and estimates.  This is useful for 
!     trying to explain windowed data.  If need be we can also estimate
!     the source time functions and receiver response functions.  
! 
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth (deg)
!     FREQL      frequency list for inversion
!     IMODREC    (1) phase, (2) amplitude, (3) both, 0 dont touch RRF
!     IMODSRC    (1) phase, (2) amplitude, (3) both, 0 dont touch STF
!     IRESTP     residual type (1) phase, (2) amplitude (3) both
!     LUNWRAP    True -> work with unwrapped phase
!     MDIM       leading dimension
!     MFREQL     leading dimension
!     MREC       leading dimension
!     NDIM       number of components
!     NFREQL     number of frequencies to invert at
!     NREC       number of receivers
!     NSRC       number of sources
!     PROJNM     project name
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     EST        estimate (U,V,W) format
!     LUPDREC    False -> RRFs have been updated 
!     LUPDSRC    False -> STFs have been updated
!     IERR       error flag
!     OBS        observations (N,E,Z) format
!     RECV       new estimate for receiver response functions
!     SOURCE     new estimate for source time function 
!     WGHT       observation weights
!  
!.... variable declarations
      IMPLICIT NONE 
      CHARACTER(80), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: FREQL(NFREQL), AZMOD 
      INTEGER*4, INTENT(IN) :: MDIM,MFREQL,MREC,  NDIM,NFREQL,NREC,NSRC,
     ;                         IMODSRC, IMODREC, IRESTP  
      LOGICAL*4, INTENT(IN) :: LUNWRAP 
 
      COMPLEX*8, INTENT(INOUT) :: SOURCE(MFREQL,*), RECV(MDIM,MFREQL,*) 
      COMPLEX*8, INTENT(OUT) :: OBS(MDIM,MFREQL,MREC,*), 
     ;                          EST(MDIM,MFREQL,MREC,*)
      REAL*4, INTENT(OUT) :: WGHT(MDIM,MFREQL,MREC,*)
      INTEGER*4, INTENT(OUT) :: IERR
      LOGICAL*4, INTENT(INOUT) :: LUPDSRC, LUPDREC
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FLNAME
      COMPLEX*8 CARG
      REAL*4 FREQIN, R1, Q1, XMAG, PHASE, TOL
      INTEGER*4 INFO(13), MYEND, NBYTES, NFREQIN, NSRCIN, NRECIN,
     ;          NDIMIN, NADV, ITYPE, ISGRNS, IUNIT, INDX, IFREQL,
     ;          I, ISRC, IREC, JFREQL
      LOGICAL*4 LSWAP, LEX
      REAL*4 UNPACKR4, WRAPPH
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(IUNIT = 44)
      PARAMETER(TOL = 1.11E-5)
!
!----------------------------------------------------------------------!
!
!.... get the observations 
      IERR = 0
      CALL RDTOBS25(PROJNM, MDIM,MFREQL,MREC,NFREQL, NDIM,NREC,
     ;              NSRC, LUNWRAP, FREQL, WGHT,OBS, IERR)
      IF (IERR.NE.0) THEN
         WRITE(*,*) 'rd_obs_est: Error reading observations!'
         RETURN
      ENDIF
!
!.... file detection and endianness
      FLNAME(1:80) = ' '
      FLNAME = './tobs/'//TRIM(PROJNM)//'.fwndo'
      FLNAME = ADJUSTL(FLNAME)
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... open and read
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN
         IERR = 1
         WRITE(*,*) 'rd_obs_est: Error file does not exist'
         RETURN
      ENDIF
!     CALL STAT(TRIM(FLNAME),INFO)
!     NBYTES = INFO(8)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(I),I=1,NBYTES)
      CLOSE(IUNIT)
      IF (IERR.NE.0) RETURN
!
!.... unpack the header
      NFREQIN = UNPACKI4(LSWAP,CDAT( 1:4))
      NSRCIN  = UNPACKI4(LSWAP,CDAT( 5:8))
      NRECIN  = UNPACKI4(LSWAP,CDAT( 9:12))
      NDIMIN  = UNPACKI4(LSWAP,CDAT(13:16))
      ITYPE   = UNPACKI4(LSWAP,CDAT(17:20))
      ISGRNS  = UNPACKI4(LSWAP,CDAT(21:24))
      IF (NRECIN.NE.NREC)
     ;WRITE(*,*) 'rd_obs_est: Warning sources do not match'
      IF (NSRCIN.NE.NSRC)
     ;WRITE(*,*) 'rd_obs_est: Warning sources do not match'
      IF (NDIMIN.NE.NDIM) THEN
         WRITE(*,*) 'rd_obs_est: Error dimensions do not match'
         IERR = 1
         RETURN
      ENDIF
      IF (ISGRNS.EQ.1) THEN
         WRITE(*,*) 'rd_obs_est: (u,v,w) windows are Greens functions'
         IF (IMODSRC.GT.0 .AND. LUPDSRC)
     ;   WRITE(*,*) 'rd_obs_est: I will estimate the STFs...'
         IF (IMODREC.GT.0 .AND. LUPDREC)
     ;   WRITE(*,*) 'rd_obs_est: I will estimate the RRFs...'
      ENDIF
      NADV = 4*2*NDIMIN*NSRCIN*NRECIN !real,imaginary, 3 components
!
!.... search for frequency in frequency list
      DO 1 IFREQL=1,NFREQL
         INDX = 24 !length of header
         DO 2 JFREQL=1,NFREQIN
            FREQIN = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
            INDX = INDX + 4
            IF (ABS(FREQIN - SNGL(FREQL(IFREQL))).LT.TOL) THEN
               GOTO 20
            ELSE
               INDX = INDX + NADV
            ENDIF
   2     CONTINUE
         WRITE(*,*) 'rd_obs_est: Error could not locate frequency'
         IERR = 1
         RETURN
  20     CONTINUE !got it
         DO 3 ISRC=1,NSRCIN
            DO 4 IREC=1,NRECIN
               DO 5 I=1,NDIMIN
                  R1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                  INDX = INDX + 4
                  Q1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                  INDX = INDX + 4
                  IF (.NOT.LUNWRAP) THEN !no phase unwrapping
                     IF (ITYPE.EQ.1) THEN !input stored mag,phase
                        CARG = CMPLX(R1,0.0)*CEXP(CMPLX(0.0,Q1))
                        R1 = REAL(CARG)
                        Q1 = IMAG(CARG)
                     ENDIF
                  ELSE !i am unwrapping the phase
                     IF (ITYPE.EQ.0) THEN !input stored complex
                        XMAG = SQRT(R1**2 + Q1**2)
                        PHASE = ATAN2(Q1,R1)
                        PHASE = WRAPPH(0,PHASE)
                        R1 = XMAG
                        Q1 = PHASE
                     ENDIF
                  ENDIF
                  EST(I,IFREQL,IREC,ISRC) = CMPLX(R1,Q1)
    5          CONTINUE
    4       CONTINUE !loop on recievers
    3    CONTINUE !looop on sources
    1 CONTINUE !loop on frequency list 
      DEALLOCATE(CDAT)
!
!.... update source
      IF (ISGRNS.EQ.1 .OR. LUPDSRC) THEN
         IF (IMODSRC.GT.0) THEN
            WRITE(*,*) 'rd_obs_est: Updating STF...'
            CALL SRCUPD(MDIM,MFREQL,MREC, NDIM,NFREQL,NREC,NSRC,
     ;                  LUNWRAP,IMODSRC, AZMOD,EST,OBS,
     ;                  SOURCE)
            IF (IRESTP.EQ.1) THEN
               WRITE(*,*) 'rd_obs_est: Setting STF amplitude to unity'
            ELSEIF (IRESTP.EQ.2) THEN
               WRITE(*,*) 'rd_obs_est: Setting RRF phase to 0'
            ENDIF
            CALL MODSRC(MFREQL, NFREQL,NSRC,IRESTP, SOURCE)
            WRITE(*,*) 'rd_obs_est: Updating estimates...'
            CALL CONVEST_STF(MDIM,MFREQL,MREC,
     ;                       NDIM,NFREQL,NREC,NSRC,
     ;                       LUNWRAP,SOURCE, EST)
            LUPDSRC = .FALSE. !no longer need to update STF
         ENDIF !end check on whether or not to modify STF
      ENDIF !end check on possibly updating STF 
      IF (ISGRNS.EQ.1 .OR. LUPDREC) THEN
         IF (IMODREC.GT.0) THEN
            WRITE(*,*) 'rd_obs_est: Updating RRF...'
            CALL RECUPD(MDIM,MFREQL,MREC, NDIM,NFREQL,NREC,NSRC,
     ;                  LUNWRAP,IMODREC, AZMOD, EST,OBS,
     ;                  RECV)
            IF (IRESTP.EQ.1) THEN
               WRITE(*,*) 'rd_obs_est: Setting RRF amplitude to unity'
            ELSEIF (IRESTP.EQ.2) THEN
               WRITE(*,*) 'rd_obs_est: Setting RRF phase to 0'
            ENDIF
            CALL MODREC(MDIM,MFREQL, NDIM,NFREQL,NREC, IRESTP,RECV)
            CALL CONVEST_RRF(MDIM,MFREQL,NREC,
     ;                       NDIM,NFREQL,NREC,NSRC,
     ;                       LUNWRAP,RECV, EST)
            LUPDREC = .FALSE. !no longer need to update RRF
         ENDIF !end check on whether or not to modify RRF
      ENDIF !end check on possibly updating RRF
      RETURN
      END
