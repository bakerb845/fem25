      SUBROUTINE SRCUPD(MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, 
     ;                  LUNWRAP,IRESTP, AZMOD, EST,OBS, SOURCE)
!
!     Updates the source time function.  Here we just keep one source
!     time function for all components.  This is a full Newton step and 
!     does not require iteration since the objective function is from 
!     an inner product, L2 norm, and amounts to finding the minima of
!     a parabaloid. - B. Baker February 2013 
!
!     INPUT      MEANING
!     -----      ------- 
!     AZMOD      model azimuth (degrees)
!     EST        estimates at receivers (already rotated!)
!     IRESTP     =1 -> Phase only 
!                =2 -> Amplitude only 
!                =3 -> phase and amplitude (default)
!     LUNWRAP    True -> data is unwrapped and stored (amp,phase)
!                False -> data is stored as a complex number 
!     MDIM       leading dimension for estimates/observations
!     MFREQ      leading dimension
!     NDIM       number of components for receivers
!     NFREQ      number of frequencies 
!     NREC       number of receivers
!     NSRC       number of sources
!     OBS        observations at receivers 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     SOURCE     updated source time function
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,MFREQ,MREC,*), 
     ;                         EST(MDIM,MFREQ,MREC,*)  
      REAL*8, INTENT(IN) :: AZMOD
      INTEGER*4, INTENT(IN) :: MFREQ, MDIM, MREC, NFREQ, NSRC, NREC, 
     ;                         NDIM, IRESTP
      LOGICAL*4 LUNWRAP 
      COMPLEX*8, INTENT(OUT) :: SOURCE(MFREQ,*)  
!.... local variables 
      COMPLEX*16 ESTR(3), ZNUM, ZZERO, U, V, W, D, COSAZ, SINAZ 
      COMPLEX*8 CZERO
      REAL*8 RDEN, ROBS, QOBS, REST, QEST, ADATA, PDATA, AEST, PEST,
     ;       PI180, OAMP, OPHS, EAMP, EPHS, DPH, 
     ;       UMAG, VMAG, WMAG, UPHASE, VPHASE, WPHASE
      REAL*4 SRCPHS, AMP
      INTEGER*4 IFREQ, ISRC, IREC, I
      PARAMETER(ZZERO = DCMPLX(0.D0,0.D0), CZERO = CMPLX(0.0,0.0))
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
!
!----------------------------------------------------------------------!
!
!.... initialize rotation angles
      COSAZ = DCMPLX(DCOS(AZMOD*PI180),0.D0)
      SINAZ = DCMPLX(DSIN(AZMOD*PI180),0.D0) 
      IF (IRESTP.EQ.1 .OR. IRESTP.EQ.3) THEN
         OPEN(31,FILE='src_phaseN.txt')
         OPEN(32,FILE='src_phaseE.txt')
         OPEN(33,FILE='src_phaseZ.txt')
      ENDIF
      IF (IRESTP.EQ.2 .OR. IRESTP.EQ.3) THEN
         OPEN(34,FILE='src_ampN.txt')
         OPEN(35,FILE='src_ampE.txt')
         OPEN(36,FILE='src_ampZ.txt')
      ENDIF
!.... loop on frequencies
      DO 100 IFREQ=1,NFREQ 
!
!....... loop on sources 
         DO 200 ISRC=1,NSRC
            ZNUM = ZZERO
            RDEN = 0.D0  
            DO 300 IREC=1,NREC
!
!............. rotate estimates counter-clockwise into obs frame
               IF (LUNWRAP) THEN
                  UMAG = DBLE(REAL(EST(1,IFREQ,IREC,ISRC)))
                  VMAG = DBLE(REAL(EST(2,IFREQ,IREC,ISRC)))
                  WMAG = DBLE(REAL(EST(3,IFREQ,IREC,ISRC)))
                  UPHASE = DBLE(IMAG(EST(1,IFREQ,IREC,ISRC)))
                  VPHASE = DBLE(IMAG(EST(2,IFREQ,IREC,ISRC)))
                  WPHASE = DBLE(IMAG(EST(3,IFREQ,IREC,ISRC)))
                  U = DCMPLX(UMAG,0.D0)*CDEXP(DCMPLX(0.D0,UPHASE))
                  V = DCMPLX(VMAG,0.D0)*CDEXP(DCMPLX(0.D0,VPHASE))
                  W = DCMPLX(WMAG,0.D0)*CDEXP(DCMPLX(0.D0,WPHASE))
               ELSE !data is complex
                  U = DCMPLX(EST(1,IFREQ,IREC,ISRC))
                  V = DCMPLX(EST(2,IFREQ,IREC,ISRC))
                  W = DCMPLX(EST(3,IFREQ,IREC,ISRC))
               ENDIF 
               ESTR(1) = U*COSAZ - V*SINAZ 
               ESTR(2) = U*SINAZ + V*COSAZ
               ESTR(3) = W 
!
!............. loop on components
               DO 400 I=1,NDIM
                  IF (OBS(I,IFREQ,IREC,ISRC).EQ.CZERO) THEN 
                     GOTO 450 !nothing here
                  ENDIF
!
!................ estimate phases are unwrapped 
                  IF (LUNWRAP) THEN 
                     OPHS = DBLE(IMAG(OBS(I,IFREQ,IREC,ISRC))) !phase
                     EPHS = DBLE(IMAG(EST(I,IFREQ,IREC,ISRC))) !phase
                     IF (IRESTP.EQ.1) THEN
                        OAMP = 1.D0
                        EAMP = 1.D0 
                     ELSE
                        OAMP = DBLE(REAL(OBS(I,IFREQ,IREC,ISRC)))
                        EAMP = DBLE(REAL(EST(I,IFREQ,IREC,ISRC))) 
                     ENDIF
                     DPH =-EPHS + OPHS !this is u*d 
                     ZNUM = ZNUM + DCMPLX(EAMP*OAMP,0.D0)
     ;                            *CDEXP(DCMPLX(0.D0,DPH))
                     RDEN = RDEN + EAMP**2
                  ELSE !estimates are complex
                     ROBS = DBLE(REAL(OBS(I,IFREQ,IREC,ISRC)))
                     QOBS = DBLE(IMAG(OBS(I,IFREQ,IREC,ISRC)))
                     REST = DREAL(ESTR(I))
                     QEST = DIMAG(ESTR(I))
                     IF (IRESTP.EQ.1 .OR. IRESTP.EQ.2) THEN
                        IF (IRESTP.EQ.2) THEN !amplitude only
                           ADATA = DSQRT(ROBS**2 + QOBS**2) 
                           PDATA = 0.D0
                           AEST  = DSQRT(REST**2 + QEST**2) 
                           PEST  = 0.D0
                        ELSE !phase is default here
                           ADATA = 1.D0
                           PDATA = DATAN2(QOBS,ROBS)
                           AEST  = 1.D0
                           PEST  = DATAN2(QEST,REST) 
                        ENDIF
                        D = DCMPLX(ADATA,0.D0)*CDEXP(DCMPLX(0.D0,PDATA))
                        U = DCMPLX(AEST ,0.D0)*CDEXP(DCMPLX(0.D0,PEST ))
                     ELSE !phase and amplitude 
                        D = DCMPLX(OBS(I,IFREQ,IREC,ISRC))
                        U = ESTR(I)  
                     ENDIF   
                     ZNUM = ZNUM + DCONJG(U)*D 
                     RDEN = RDEN + DREAL(DCONJG(U)*U) !|u|^2 = conj(u)*u
                  ENDIF
  450             CONTINUE !break ahead, no data 
  400          CONTINUE !loop on components 
  300       CONTINUE !loop on receivers
            IF (RDEN.GT.0.D0) THEN
               SOURCE(IFREQ,ISRC) = CMPLX(ZNUM/DCMPLX(RDEN,0.D0)) 
            ELSE
               WRITE(*,*) 'srcupd: Warning source',ISRC,'for frequency',
     ;                    IFREQ,'is dead'
               WRITE(*,*) '        Setting response to zero'
               SOURCE(IFREQ,ISRC) = CMPLX(0.0,0.0)
            ENDIF 
            IF (LUNWRAP) THEN
               DO I=1,NDIM
                  DO IREC=1,NREC
                     IF (CABS(OBS(I,IFREQ,IREC,ISRC)).GT.0.0) THEN
                        AMP = CABS(SOURCE(IFREQ,ISRC))
                        SRCPHS = ATAN2(IMAG(SOURCE(IFREQ,ISRC)),
     ;                                 REAL(SOURCE(IFREQ,ISRC)))
                        IF (IRESTP.EQ.1 .OR. IRESTP.EQ.3) THEN
                           WRITE(30+I,*)IFREQ,ISRC,IREC, SRCPHS, 
     ;                            IMAG(OBS(I,IFREQ,IREC,ISRC)),
     ;                            IMAG(EST(I,IFREQ,IREC,ISRC)), 
     ;                            IMAG(EST(I,IFREQ,IREC,ISRC))+SRCPHS
                        ENDIF
                        IF (IRESTP.EQ.2 .OR. IRESTP.EQ.3) THEN
                           WRITE(33+I,*)IFREQ,ISRC,IREC, AMP,
     ;                            REAL(OBS(I,IFREQ,IREC,ISRC)),
     ;                            REAL(EST(I,IFREQ,IREC,ISRC)),
     ;                            REAL(EST(I,IFREQ,IREC,ISRC))*AMP
                        ENDIF
                     ENDIF
                  ENDDO !loop on receivers
                  IF (IRESTP.EQ.1 .OR. IRESTP.EQ.3) WRITE(30+I,*) 
                  IF (IRESTP.EQ.2 .OR. IRESTP.EQ.3) WRITE(33+I,*)
               ENDDO 
            ENDIF 
  200    CONTINUE !loop on sources
  100 CONTINUE !loop on frequencies
      IF (IRESTP.EQ.1 .OR. IRESTP.EQ.3) THEN
         CLOSE(31)
         CLOSE(32) 
         CLOSE(33) 
      ENDIF
      IF (IRESTP.EQ.2 .OR. IRESTP.EQ.3) THEN
         CLOSE(34)
         CLOSE(35)
         CLOSE(36)
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDSRCINV_SB(PROJNM, MFREQ,NFREQ,NSRC, 
     ;                       NFREQ_SRF,NFREQ_BDY, NSRC_SRF,NSRC_BDY,
     ;                       CFTYPE,SRCTYP,FREQ, LEXS,LEXB,SOURCE,IERR)
!
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: PROJNM
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC)
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      REAL*8, INTENT(IN) :: FREQ(NFREQ) 
      INTEGER*4, INTENT(IN) :: MFREQ,NFREQ,NSRC, NFREQ_SRF,NFREQ_BDY,
     ;                         NSRC_SRF,NSRC_BDY 
      COMPLEX*8, INTENT(OUT) :: SOURCE(MFREQ,NSRC) 
      INTEGER*4, INTENT(OUT) :: IERR 
      LOGICAL*4, INTENT(OUT) :: LEXS, LEXB
!.... local variables
      CHARACTER(80) FILENM
      COMPLEX*8, ALLOCATABLE :: SRCWRK(:,:)  
      REAL*8, ALLOCATABLE :: FREQ_LOC(:) 
!
!----------------------------------------------------------------------!
!
      IERR = 0
      LEXS = .TRUE.
      LEXB = .TRUE.
      IF (NFREQ_SRF.GT.0 .AND. NSRC_SRF.GT.0) THEN
         ALLOCATE(SRCWRK(NFREQ_SRF,NSRC_SRF)) 
         ALLOCATE(FREQ_LOC(NFREQ_SRF))
         CALL GET_FREQ_MOD(NFREQ, NFREQ_SRF, .TRUE., CFTYPE,FREQ, 
     ;                     FREQ_LOC,IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'rdsrcinv_sb: Error calling get_freq_mod 1!'
            RETURN
         ENDIF
         FILENM(1:80) = ' '
         FILENM = TRIM(ADJUSTL(PROJNM))//'_srf'
         CALL RDSRCINV(NFREQ_SRF,FILENM, NFREQ_SRF,NSRC_SRF, FREQ_LOC, 
     ;                 LEXS,SRCWRK) 
         CALL SET_SRC_MOD(NFREQ,NFREQ_SRF, NFREQ,NSRC,
     ;                    NFREQ_SRF,NSRC_SRF, .TRUE., CFTYPE,SRCTYP, 
     ;                    SRCWRK, SOURCE,IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'rdsrcinv_sb: Error calling set_src_mod 1!'
            RETURN
         ENDIF 
         DEALLOCATE(SRCWRK) 
         DEALLOCATE(FREQ_LOC)
      ENDIF
      IF (NFREQ_BDY.GT.0 .AND. NSRC_BDY.GT.0) THEN
         ALLOCATE(SRCWRK(NFREQ_BDY,NSRC_BDY)) 
         ALLOCATE(FREQ_LOC(NFREQ_BDY))
         CALL GET_FREQ_MOD(NFREQ, NFREQ_BDY, .FALSE., CFTYPE,FREQ, 
     ;                     FREQ_LOC,IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'rdsrcinv_sb: Error calling get_freq_mod 2!'
            RETURN 
         ENDIF
         FILENM(1:80) = ' '
         FILENM = TRIM(ADJUSTL(PROJNM))//'_bdy'
         CALL RDSRCINV(NFREQ_BDY,FILENM, NFREQ_BDY,NSRC_BDY, FREQ_LOC,
     ;                 LEXB,SRCWRK) 
         CALL SET_SRC_MOD(NFREQ,NFREQ_BDY, NFREQ,NSRC, 
     ;                    NFREQ_BDY,NSRC_BDY, .FALSE., CFTYPE,SRCTYP, 
     ;                    SRCWRK, SOURCE, IERR)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'rdsrcinv_sb: Error calling set_src_mod 2!'
            RETURN
         ENDIF
         DEALLOCATE(SRCWRK) 
         DEALLOCATE(FREQ_LOC)
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDSRCINV(MFREQ,PROJNM, NOMINV,NSRC, FREQINV,
     ;                    LSRCEX,SOURCE)
!
!     Checks for a .srcinv file.  If it exists, reads it. 
!     If it doesn't then sets the source time function to unity
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQINV    frequencies in inversion block
!     MFREQ      leading dimension for SOURCE 
!     NOMINV     number of frequencies in inversion
!     NSRC       number of sources
!     PROJNM     project name   
!
!     OUTPUT     MEANING
!     ------     ------- 
!     LSRCEX     true -> source exists
!     SOURCE     source time function at frequencies
!
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: FREQINV(NOMINV)
      INTEGER*4, INTENT(IN) :: MFREQ,NOMINV
      COMPLEX*8, INTENT(OUT) :: SOURCE(MFREQ,*)
      LOGICAL*4, INTENT(OUT) :: LSRCEX
!.... local variables 
      CHARACTER(80) FLNAME
      CHARACTER(1), ALLOCATABLE :: CDAT(:) 
      LOGICAL*4 LEX, LONE,LSWAP
      REAL*4 UNPACKR4
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(TOL = 1.11E-4)
      PARAMETER(IUNIT = 55) 
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0 
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... detect if file exists
      LSRCEX = .FALSE.
      LONE = .FALSE.
      FLNAME(1:80) = ' ' 
      FLNAME = TRIM(PROJNM)//'.srcinv'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES)
      IF (LEX) THEN !file exists, try to read it
!        CALL STAT(TRIM(FLNAME),INFO)
!        NBYTES = INFO(8)
         OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;        ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
         ALLOCATE(CDAT(NBYTES))
         READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,NBYTES)
         CLOSE(IUNIT)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'srcinv: Error reading srcinv file'
            LONE = .TRUE.
            GOTO 50
         ENDIF
! 
!....... read the header
         NFREQIN = UNPACKI4(LSWAP,CDAT(1:4))
         NSRCIN  = UNPACKI4(LSWAP,CDAT(5:8))
         IF (NSRC.NE.NSRCIN) THEN
            WRITE(*,*) 'srcinv: Warning nsrcin != nsrc'
            WRITE(*,*) 'srcinv: Defaulting to unity STF'
            LONE = .TRUE.
            GOTO 50
         ENDIF
!
!....... attempt to parse file
         DO 1 IFREQ=1,NOMINV
            INDX = 8
            DO 2 JFREQ=1,NFREQIN
               FREQIN = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
               INDX = INDX + 4
               IF (ABS(FREQIN - SNGL(FREQINV(IFREQ))).LT.TOL) GOTO 20
               INDX = INDX + 4*2*NSRC
   2        CONTINUE
            WRITE(*,*) 'srcinv: Could not detect inversion frequency!'
            LONE = .TRUE.
            GOTO 50
  20        CONTINUE !source exists
!
!.......... got it, read the sources at this frequency
            DO 3 ISRC=1,NSRC
               R = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
               INDX = INDX + 4
               Q = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
               INDX = INDX + 4
               SOURCE(IFREQ,ISRC) = CMPLX(R,Q)
   3        CONTINUE
   1     CONTINUE !loop on inversion frequencies
      ELSE
         WRITE(*,*) 'srcinv: No srcinv file detected'
         WRITE(*,*) 'srcinv: Defaulting to Greens functions'
         LONE = .TRUE.
      ENDIF
!
!.... set to 1
   50 CONTINUE
      IF (LONE) THEN
         WRITE(*,*) 'srcinv: Setting source time function set to unity'
         DO 21 IFREQ=1,NOMINV
            DO 22 ISRC=1,NSRC
               SOURCE(IFREQ,ISRC) = CMPLX(1.0,0.0)
   22       CONTINUE
   21    CONTINUE
      ELSE
         LSRCEX = .TRUE.
      ENDIF
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE WTSTF(PROJNM, MFREQ, NFREQ,NSRC, IBLOCK,ITER, 
     ;                 FREQ,SOURCE, IERR)
!
!     Write the source time function for this inversion block  
!     and iteration 
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       frequency list
!     IBLOCK     block number in inversion
!     ITER       iteration in inversion
!     MFREQ      leading dimension for SOURCE
!     NFREQ      number of frequencies
!     NSRC       number of sources
!     PROJNM     project name
!     SOURCE     source time function for frequency/source
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: PROJNM
      COMPLEX*8, INTENT(IN) :: SOURCE(MFREQ,*) 
      REAL*8, INTENT(IN) :: FREQ(NFREQ) 
      INTEGER*4, INTENT(IN) :: MFREQ, NFREQ, NSRC, IBLOCK, ITER
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FLNAME 
      CHARACTER(5) CITER
      CHARACTER(3) CBLOCK
      CHARACTER(4), ALLOCATABLE :: CDAT(:)  
      REAL*4 R, Q
      INTEGER*4 NBYTES, LWORK, MYEND, INDX, IFREQ, ISRC, J, IUNIT  
      LOGICAL*4 LEX, LSWAP, LISDIR 
      CHARACTER(4) PACKI4, PACKR4
      INTEGER*4 ENDIAN 
      PARAMETER(IUNIT = 68) 
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0 
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE. 
!.... allocate space and pack header
      NBYTES = 8 + 4*NFREQ + 8*NFREQ*NSRC  
      LWORK = NBYTES/4
      ALLOCATE(CDAT(LWORK)) 
      CDAT(1) = PACKI4(LSWAP,NFREQ)
      CDAT(2) = PACKI4(LSWAP,NSRC)
      INDX = 3
      DO 1 IFREQ=1,NFREQ
         CDAT(INDX) = PACKR4(LSWAP,SNGL(FREQ(IFREQ)))
         INDX = INDX + 1
         DO 2 ISRC=1,NSRC
            R = REAL(SOURCE(IFREQ,ISRC))
            Q = IMAG(SOURCE(IFREQ,ISRC))
            CDAT(INDX) = PACKR4(LSWAP,R)
            INDX = INDX + 1
            CDAT(INDX) = PACKR4(LSWAP,Q)
            INDX = INDX + 1
    2    CONTINUE
    1 CONTINUE 
      INDX = INDX - 1 
      IF (INDX*4.NE.NBYTES) THEN
         WRITE(*,*) 'wtstf: Error this file is the wrong size!'
         IERR = 1
      ENDIF  
!
!.... file handling
      LEX = LISDIR('./stf')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./stf')
!.... set file name and write
      FLNAME(1:80) = ' ' 
      CITER(1:5) = ' '
      CBLOCK(1:3) = ' '
      WRITE(CITER,'(I5)') ITER
      WRITE(CBLOCK,'(I3)') IBLOCK
      CITER = ADJUSTL(CITER)
      CBLOCK = ADJUSTL(CBLOCK) 
      FLNAME = './stf/'//TRIM(PROJNM)//'-'//TRIM(CBLOCK)//'-'//
     ;         TRIM(CITER)//'.srcinv'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED', 
     ;     STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      WRITE(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,LWORK) 
      CLOSE(IUNIT)
      DEALLOCATE(CDAT) 
      RETURN
      END 
