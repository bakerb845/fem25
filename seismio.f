      SUBROUTINE WTEEST25(PROJNM, MFREQ,MREC,NFREQ, NREC,NSRC, 
     ;                    FREQ8, UEST,VEST,WEST, IERR)
!
!     Writes the estimates at receivers
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ8      frequency list
!     MFREQ      leading dimension for ?est
!     MREC       leading dimension for ?est
!     NFREQ      number of frequencies
!     NREC       number of recievers
!     PROJNM     project name
!     UEST       estimate wavefield in u 
!     VEST       estimate wavefield in v 
!     WEST       estimate waveifeld in w 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      COMPLEX*8, INTENT(IN) :: UEST(MFREQ,MREC,*), VEST(MFREQ,MREC,*),
     ;                         WEST(MFREQ,MREC,*)
      REAL*8, INTENT(IN) :: FREQ8(NFREQ) 
      INTEGER*4, INTENT(IN) :: MFREQ,MREC,NFREQ, NREC,NSRC
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FLNAME
      CHARACTER(4), ALLOCATABLE :: CDAT(:)
      LOGICAL*4 LEX,LSWAP,LISDIR
      CHARACTER(4) PACKI4, PACKR4
      INTEGER*4 ENDIAN
      PARAMETER(NDIM = 3, IUNIT = 45)
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE. 
!.... allocate space and pack header
      NBYTES = 16 + 4*NFREQ + 8*NDIM*NFREQ*NSRC*NREC 
      LWORK = NBYTES/4
      ALLOCATE(CDAT(LWORK)) 
      CDAT(1) = PACKI4(LSWAP,NFREQ)
      CDAT(2) = PACKI4(LSWAP,NSRC)
      CDAT(3) = PACKI4(LSWAP,NREC)
      CDAT(4) = PACKI4(LSWAP,NDIM)
!
!.... pack data
      INDX = 5
      DO 1 IFREQ=1,NFREQ
         CDAT(INDX) = PACKR4(LSWAP,SNGL(FREQ8(IFREQ)))
         INDX = INDX + 1
         DO 2 ISRC=1,NSRC
            DO 3 IREC=1,NREC
               CDAT(INDX) = PACKR4(LSWAP,REAL(UEST(IFREQ,IREC,ISRC)))
               INDX = INDX + 1
               CDAT(INDX) = PACKR4(LSWAP,IMAG(UEST(IFREQ,IREC,ISRC)))
               INDX = INDX + 1

               CDAT(INDX) = PACKR4(LSWAP,REAL(VEST(IFREQ,IREC,ISRC)))
               INDX = INDX + 1
               CDAT(INDX) = PACKR4(LSWAP,IMAG(VEST(IFREQ,IREC,ISRC)))
               INDX = INDX + 1

               CDAT(INDX) = PACKR4(LSWAP,REAL(WEST(IFREQ,IREC,ISRC)))
               INDX = INDX + 1
               CDAT(INDX) = PACKR4(LSWAP,IMAG(WEST(IFREQ,IREC,ISRC)))
               INDX = INDX + 1
    3       CONTINUE !loop on receivers
    2    CONTINUE !loop on sources
    1 CONTINUE !loop on frequencies
      INDX = INDX - 1
      IF (INDX*4.NE.NBYTES) THEN
         WRITE(*,*) 'wteest25: Error this file is the wrong size'
         IERR = 1
      ENDIF 
!
!.... file handling
      LEX = LISDIR('./dest')
!#ifdef INTEL
!      INQUIRE(DIRECTORY='./dest',EXIST=LEX)
!#else
!      INQUIRE(FILE='./dest',EXIST=LEX)
!#endif
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./dest')
      FLNAME(1:80) = ' '
      FLNAME = './dest/'//TRIM(PROJNM)//'.edest'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      WRITE(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,LWORK)
      CLOSE(IUNIT) 
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDTOBS_HD(PROJNM, NDIM,NFREQ,NREC,NSRC, IERR)
!
!     Reads the frequency domain observation headers 
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NDIM       number of componenets in file
!     NFREQ      number of frequencies in file
!     NREC       number of receivers in file
!     NSRC       number of sources in file
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(OUT) :: NDIM,NFREQ,NREC,NSRC, IERR 
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FLNAME
      LOGICAL*4 LEX, LSWAP
      PARAMETER(IUNIT = 44) 
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(TOL = 1.11E-5)
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... file detection
      FLNAME(1:80) = ' '
      FLNAME = './tobs/'//TRIM(PROJNM)//'.etobs'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdtobs25_hd: Error cannot detect observation file'
         IERR = 1
         RETURN
      ENDIF
!
!.... open and read header
      NBYTES = 16
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,NBYTES)
      CLOSE(IUNIT)
      IF (IERR.NE.0) THEN
         NFREQ = 0
         NSRC = 0
         NREC = 0
         NDIM = 0
         GOTO 500
      ENDIF
!
!.... unpack the header
      NFREQ = UNPACKI4(LSWAP,CDAT( 1:4))
      NSRC  = UNPACKI4(LSWAP,CDAT( 5:8))
      NREC  = UNPACKI4(LSWAP,CDAT( 9:12))
      NDIM  = UNPACKI4(LSWAP,CDAT(13:16))
  500 CONTINUE
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDTOBS_SB_HD(PROJNM, NDIM,NFREQ_SRF,NFREQ_BDY, 
     ;                        NREC,NSRC_SRF,NSRC_BDY, IERR)
!
!     Reads the frequency domain observation headers 
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NDIM       number of componenets in file
!     NFREQ_BDY  number of body wave frequencies in file
!     NFREQ_SRF  number of surface wave frequencies in file 
!     NREC       number of receivers in file
!     NSRC_BDY   number of body wave sources 
!     NSRC_SRF   number of surface wave sources
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(OUT) :: NDIM,NFREQ_SRF,NFREQ_BDY, 
     ;                          NREC,NSRC_SRF,NSRC_BDY, IERR 
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FLNAME
      LOGICAL*4 LEX, LSWAP
      PARAMETER(IUNIT = 44)  
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(TOL = 1.11E-5)
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... file detection
      FLNAME(1:80) = ' '
      FLNAME = './tobs/'//TRIM(PROJNM)//'.sb_etobs'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=FLNAME,EXIST=LEX)
      IF (.NOT.LEX) THEN 
         WRITE(*,*) 'rdtobs25_hd: Error cannot detect observation file'
         IERR = 1
         RETURN
      ENDIF
!
!.... open and read header
      NBYTES = 24
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,NBYTES)
      CLOSE(IUNIT)
      IF (IERR.NE.0) THEN
         NFREQ_SRF = 0
         NFREQ_BDY = 0 
         NSRC_SRF = 0
         NSRC_BDY = 0
         NREC = 0
         NDIM = 0
         GOTO 500
      ENDIF
!
!.... unpack the header
      NFREQ_SRF = UNPACKI4(LSWAP,CDAT( 1:4))
      NFREQ_BDY = UNPACKI4(LSWAP,CDAT( 5:8))
      NSRC_SRF  = UNPACKI4(LSWAP,CDAT( 9:12))
      NSRC_BDY  = UNPACKI4(LSWAP,CDAT(13:16))
      NREC      = UNPACKI4(LSWAP,CDAT(17:20))
      NDIM      = UNPACKI4(LSWAP,CDAT(21:24))
  500 CONTINUE 
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
!
      SUBROUTINE RDTOBS25(PROJNM, MDIM,MFREQL,MREC,NFREQL, NDIM,NREC,
     ;                    NSRC, LUNWRAP, FREQL, WGHT,TOBS, IERR)   
!
!     Reads the frequency domain observations frequencies FREQL 
!     and the corresponding data weights 
! 
!     INPUT      MEANING
!     -----      ------- 
!     FREQL      frequency list
!     LUNWRAP    True -> unwrapping phase so TOBS will be stored
!                        (magnitude,phase)
!                FALSE -> TOBS stored as a complex number
!     MFREQL     leading dimension for ?est
!     MREC       leading dimension for ?est
!     NDIM       number of spatial dimensions expected, should be 3
!     NFREQL     number of frequencies in list  
!     NREC       number of recievers we are expecting
!     NSRC       number of sources we are expecting
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     TOBS       (N,E,Vertical) observation at receivers
!     WGHT       (N,E,Vertical) weights on observations
! 
!.... variable declarations 
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: FREQL(NFREQL)
      INTEGER*4, INTENT(IN) :: MFREQL,MREC,NFREQL, NDIM,NREC,NSRC
      LOGICAL*4, INTENT(IN) :: LUNWRAP
      COMPLEX*8, INTENT(OUT) :: TOBS(MDIM,MFREQL,MREC,*)
      REAL*4, INTENT(OUT) :: WGHT(MDIM,MFREQL,MREC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FLNAME
      COMPLEX*8 CARG
      REAL*4 R1, Q1, XMAG, PHASE, WRAPPH
      LOGICAL*4 LSWAP, LEX
      PARAMETER(IUNIT = 44)
      REAL*4 UNPACKR4
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(TOL = 1.11E-5)
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... file detection
      FLNAME(1:80) = ' '
      FLNAME = './tobs/'//TRIM(PROJNM)//'.etobs'
      FLNAME = ADJUSTL(FLNAME)
!
!.... open and read
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN
         IERR = 1
         WRITE(*,*) 'rdtobs: Error file does not exist',TRIM(FLNAME)
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
      IF (NRECIN.NE.NREC)
     ;WRITE(*,*) 'rdtobs25: Warning sources do not match'
      IF (NSRCIN.NE.NSRC)
     ;WRITE(*,*) 'rdtobs25: Warning sources do not match'
      IF (NDIMIN.NE.NDIM) THEN
         WRITE(*,*) 'rdtobs25: Error dimensions do not match'
         IERR = 1
         RETURN
      ENDIF
      NADV = 36*NSRCIN*NRECIN
!
!.... search for frequency in frequency list
      DO 1 IFREQL=1,NFREQL
         INDX = 20 !length of header
         DO 2 JFREQL=1,NFREQIN 
            FREQIN = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
            INDX = INDX + 4
            IF (ABS(FREQIN - SNGL(FREQL(IFREQL))).LT.TOL) THEN
               GOTO 20
            ELSE
               INDX = INDX + NADV
            ENDIF
   2     CONTINUE
         WRITE(*,*) 'rdtobs25: Error could not locate frequency'
         IERR = 1
         RETURN
  20     CONTINUE !got it
         DO 3 ISRC=1,NSRCIN
            DO 4 IREC=1,NRECIN
               DO 5 I=1,NDIMIN
                  W1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                  INDX = INDX + 4
                  R1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                  INDX = INDX + 4
                  Q1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                  INDX = INDX + 4
                  IF (R1.EQ.0.0 .AND. Q1.EQ.0.0 .AND. W1.NE.0.0) THEN
                     WRITE(*,*) 'rdtobs25: Warning null obs has weight'
                     W1 = 0.0
                  ENDIF
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
                  TOBS(I,IFREQL,IREC,ISRC) = CMPLX(R1,Q1)
                  WGHT(I,IFREQL,IREC,ISRC) = W1
    5          CONTINUE
    4       CONTINUE !loop on recievers
    3    CONTINUE !looop on sources
    1 CONTINUE !loop on frequency list 
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDTOBS_SB25(PROJNM, RCV,SRC,FRQ,INV, IERR) 
!
!     Reads the observations for surface and body waves.  If we are 
!     windowing surface waves or body waves we need all the observations
!     for STF and RRF updates
!
!     INPUT          MEANING
!     -----          -------
!     LINVF          True -> frequency is an inversion frequency
!     LUNWRAP        True -> observations are (mag,phase)
!     MFREQ_INV      leading dimension 
!     NFREQ          number of frequencies to model (surface+body)
!     NFREQ_BDY      number of body wave frequencies to model
!     NFREQ_SRF      number of surface wave frequencies to model
!     NFREQ_BDY_INV  number of body wave inversion frequencies
!     NFREQ_SRF_INV  number of surface wave inversion frequencies
!     NREC           number of receivers
!     NSRC_BDY       number of body wave sources
!     NSRC_SRF       number of surface wave sources
!     PROJNM         project name
!     SRCTYP         source type surface or body wave 
! 
!     OUTPUT         MEANING
!     ------         -------
!     IERR           error flag
!     TOBS           observations 
!     WGHT           weights on observations
!   
!.... variable declarations 
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (RECV_INFO) RCV 
      TYPE (SRC_INFO) SRC
      TYPE (FRQ_INFO) FRQ
      TYPE (INV_INFO) INV
 
      CHARACTER(*), INTENT(IN) :: PROJNM
!     INTEGER*4, INTENT(IN) :: NFREQ_OBS
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:) 
      CHARACTER(80) FLNAME
      COMPLEX*8 CARG 
      REAL*4 UNPACKR4, XMAG, PHASE, W1, R1, Q1, TOL, FREQIN, WRAPPH  
      INTEGER*4 UNPACKI4, ENDIAN, NFREQ_SRF_IN, 
     ;          NFREQ_BDY_IN, NSRC_SRF_IN, NSRC_BDY_IN, NRECIN, NDIMIN,
     ;          ITYPE,  JFREQL, IFREQ, JFREQ, IOFF, ISRC, IREC, 
     ;          I, IUNIT, INDX, NADV, NBYTES, MYEND  
      LOGICAL*4 LSWAP, LEX 
      PARAMETER(IUNIT = 65)
      PARAMETER(TOL = 2.22E-5)
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... file detection
      FLNAME(1:80) = ' '
      FLNAME = './tobs/'//TRIM(PROJNM)//'.sb_etobs'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdtobs_sb25: Observation file does not exist!'
         IERR = 1
         RETURN
      ENDIF 
!
!.... open and read
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
      NFREQ_SRF_IN = UNPACKI4(LSWAP,CDAT( 1:4))
      NFREQ_BDY_IN = UNPACKI4(LSWAP,CDAT( 5:8))
      NSRC_SRF_IN  = UNPACKI4(LSWAP,CDAT( 9:12))
      NSRC_BDY_IN  = UNPACKI4(LSWAP,CDAT(13:16)) 
      NRECIN       = UNPACKI4(LSWAP,CDAT(17:20))
      NDIMIN       = UNPACKI4(LSWAP,CDAT(21:24))
      ITYPE        = UNPACKI4(LSWAP,CDAT(25:28))
      IF (NRECIN.NE.rcv%NREC) THEN
         WRITE(*,*) 'rdtobs_sb25: Error receiver mismatch'
         IERR = 1
      ENDIF
      IF (NSRC_SRF_IN.NE.src%NSRC_SRF) THEN
         WRITE(*,*) 'rdtobs_sb25: Error surface source mismatch'
         IERR = 1
      ENDIF
      IF (NSRC_BDY_IN.NE.src%NSRC_BDY) THEN
         WRITE(*,*) 'rdtobs_sb25: Error body source mismatch!'
         IERR = 1
      ENDIF
      IF (NDIMIN.NE.NDIM) THEN
         WRITE(*,*) 'rdtobs_sb25: Error dimensions do not match'
         IERR = 1 
      ENDIF
      IF (IERR.NE.0) RETURN

      IOFF = 28 
      NADV = 4*3*NDIM*NSRC_SRF_IN*NRECIN + 4*1 
!
!.... read the surface waves 
      DO 1 IFREQ=1,frq%NFREQ_SRF
         INDX = IOFF 
         DO 2 JFREQ=1,NFREQ_SRF_IN 
            JFREQL = JFREQ
            FREQIN = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
            IF (ABS(FREQIN - SNGL(frq%FREQ(IFREQ))).LT.TOL) THEN
               GOTO 20
            ELSE
               INDX = INDX + NADV
            ENDIF
    2    CONTINUE 
         WRITE(*,*) 'rdtobs_sb25: Error couldnt find surface wave freq'
         IERR = 1 
         RETURN
   20    CONTINUE 
         INDX = INDX + 4
         DO 3 ISRC=1,src%NSRC
            IF (src%SRCTYP(ISRC)(1:1).EQ.'S') THEN
               DO 4 IREC=1,rcv%NREC
                  DO 5 I=1,NDIM
                     W1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                     INDX = INDX + 4
                     R1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                     INDX = INDX + 4
                     Q1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                     INDX = INDX + 4
                     IF (.NOT.inv%LUNWRAP) THEN !no phase unwrapping
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
                     JFREQ = IFREQ
                     inv%OBS(I,JFREQ,IREC,ISRC) = CMPLX(R1,Q1)
                     inv%WGHTS(I,JFREQ,IREC,ISRC) = W1
    5             CONTINUE !Loop on components
    4          CONTINUE !Loop on receivers
            ENDIF !end check on source type
    3    CONTINUE 
    1 CONTINUE !loop on types
!
!.... repeat for body waves 
      IOFF = 28 
     ;     + 4*3*NDIM*NSRC_SRF_IN*NRECIN*NFREQ_SRF_IN
     ;     + 4*1*NFREQ_SRF_IN
      NADV = 4*3*NDIM*NSRC_BDY_IN*NRECIN + 4*1
      DO 11 IFREQ=1,frq%NFREQ_BDY
         INDX = IOFF
         DO 12 JFREQ=1,NFREQ_BDY_IN
            JFREQL = frq%NFREQ_SRF + JFREQ
            FREQIN = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
            IF (ABS(FREQIN 
     ;         -SNGL(frq%FREQ(frq%NFREQ_SRF+IFREQ))).LT.TOL) THEN
               GOTO 25
            ELSE
               INDX = INDX + NADV
            ENDIF
   12    CONTINUE 
         WRITE(*,*) 'rdtobs_sb25: Couldnt find body wave frequency!'
         IERR = 1
         RETURN 
   25    CONTINUE
         INDX = INDX + 4
         DO 13 ISRC=1,src%NSRC
            IF (src%SRCTYP(ISRC)(1:1).EQ.'P') THEN
               DO 14 IREC=1,rcv%NREC
                  DO 15 I=1,NDIM
                     W1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                     INDX = INDX + 4
                     R1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                     INDX = INDX + 4
                     Q1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
                     INDX = INDX + 4
                     IF (.NOT.inv%LUNWRAP) THEN !no phase unwrapping
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
                     JFREQ = IFREQ + frq%NFREQ_SRF 
                     inv%OBS(I,JFREQ,IREC,ISRC) = CMPLX(R1,Q1)
                     inv%WGHTS(I,JFREQ,IREC,ISRC) = W1
   15             CONTINUE
   14          CONTINUE 
            ENDIF !end check on source type
   13    CONTINUE
   11 CONTINUE
      DEALLOCATE(CDAT) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDWIN(PROJNM, NREC,NSRC, SRCTYP, WIN, IERR)
!
!     Reads the windows chopping down the synthetics
!
!     INPUT      MEANING
!     -----      -------
!     NREC       number of receivers
!     NSRC       number of sources
!     SRCTYP     source type
!     PROJNM     project name
!     
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     IWNDO      window type 
!     DT_BDY     sampling period (s) for body waves
!     DT_SRF     sampling period (s) for surfae waves
!     LWNDO_BDY  True -> windowing body wave synthetics
!     LWNDO_SRF  True -> windowing surface wave synthetics
!     NSAMP_BDY  number of samples in body wave synthetics
!     NSAMP_SRF  number of samples in surface wave synthetics
!     START_BDY  start time for body waves (s)
!     START_SRF  start time for surface waves (s)
!     TWIN       start stop windows (s) for windows 
! 
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (WIN_INFO) WIN 
      CHARACTER(*), INTENT(IN) :: PROJNM
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC) 
      INTEGER*4, INTENT(IN) :: NREC,NSRC
      INTEGER*4, INTENT(OUT) :: IERR
!.... variable declarations
      CHARACTER(80) FILENM
      REAL*4 TE, TB, TOL, STARTT, DT
      INTEGER*4 NREC_IN, NSRC_IN, IREC, ISRC, JREC, JSRC, IUNIT, NSAMP
      LOGICAL*4 LEX
      PARAMETER(IUNIT = 65) 
      PARAMETER(TOL = 1.11E-5)
!
!----------------------------------------------------------------------!
!
!.... initialize
      win%DT_SRF = 0.0
      win%DT_BDY = 0.0
      win%START_SRF = 0.0
      win%START_BDY = 0.0
      win%NSAMP_SRF = 0
      win%NSAMP_BDY = 0
      win%NPCT_SRF  = 0
      win%NPCT_BDY  = 0
      win%LWNDO_SRF = .FALSE.
      win%LWNDO_BDY = .FALSE. 
!.... check if the file exists
      IERR = 0
      FILENM(1:80) = ' '
      FILENM = './tobs/'//TRIM(PROJNM)//'.window' 
      FILENM = ADJUSTL(FILENM)
      INQUIRE(FILE=TRIM(FILENM),EXIST=LEX) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdwin: Window file does not exist!'
         IERR =-1
         RETURN
      ENDIF
      OPEN(UNIT=IUNIT,FILE=TRIM(FILENM)) 
      READ(IUNIT,*,END=60) NREC_IN,NSRC_IN,win%IWNDO,
     ;                     win%NPCT_SRF,win%NPCT_BDY
      IF (NSRC_IN.NE.NSRC) THEN
         WRITE(*,*) 'rdwin: Error source mismatch!'
         IERR = 1
         GOTO 60
      ENDIF
      IF (NREC_IN.NE.NREC) THEN
         WRITE(*,*) 'rdwin: Error receiver mismatch!'
         IERR = 1
         GOTO 60
      ENDIF 
!
!.... loop on sources
      DO 1 ISRC=1,NSRC
         READ(IUNIT,*,END=60) NSAMP,DT,STARTT
!
!....... get sampling rates
         IF (SRCTYP(ISRC)(1:1).EQ.'S') THEN 
            IF (win%DT_SRF.EQ.0.0) THEN 
               win%DT_SRF = DT 
               win%START_SRF = STARTT
               win%NSAMP_SRF = NSAMP
            ELSE 
               IF (ABS(win%DT_SRF - DT).GT.TOL) THEN 
                  WRITE(*,*) 'rdwin: Warning dt variation'
                  win%DT_SRF = MAX(DT,win%DT_SRF) !get away from nyq
               ENDIF
            ENDIF
         ELSE !body waves
            IF (win%DT_BDY.EQ.0.0) THEN 
               win%DT_BDY = DT 
               win%START_BDY = STARTT
               win%NSAMP_BDY = NSAMP
            ELSE 
               IF (ABS(win%DT_BDY - DT).GT.TOL) THEN 
                  WRITE(*,*) 'rdwin: Warning dt variation'
                  win%DT_BDY = MAX(DT,win%DT_BDY) !get away from nyq
               ENDIF
            ENDIF
         ENDIF
!
!....... read windows for receivers
         DO 2 IREC=1,NREC
            win%TWIN(IREC,ISRC,1) = 0.0
            win%TWIN(IREC,ISRC,2) = 0.0 
            READ(IUNIT,*,END=60) JREC,JSRC,TB,TE
            IF (JSRC.NE.ISRC) THEN
               WRITE(*,*) 'rdwin: Warning jsrc /= isrc'
               IERR = 1
               GOTO 60
            ENDIF
            IF (JREC.NE.IREC) THEN
               WRITE(*,*) 'rdwin: Warning jrec /= irec'
               IERR = 1
               GOTO 60
            ENDIF
            win%TWIN(IREC,ISRC,1) = TB
            win%TWIN(IREC,ISRC,2) = TE
            IF (ABS(TE - TB).LT.TOL) THEN 
               win%TWIN(IREC,ISRC,1) = 0.0 
               win%TWIN(IREC,ISRC,2) = 0.0
            ENDIF
            IF (win%TWIN(IREC,ISRC,1).GT.0.0 .AND. 
     ;          win%TWIN(IREC,ISRC,2).GT.0.0) THEN 
               IF (SRCTYP(ISRC)(1:1).EQ.'S') THEN
                  win%LWNDO_SRF = .TRUE.
               ELSE
                  win%LWNDO_BDY = .TRUE.
               ENDIF
            ENDIF
    2    CONTINUE
         READ(IUNIT,*,END=60) !blank line 
    1 CONTINUE !loop on frequencies 
      CLOSE(IUNIT)
!
!.... review
      IF (win%LWNDO_SRF) 
     ;WRITE(*,905) win%NSAMP_SRF,win%DT_SRF,win%START_SRF
      IF (win%LWNDO_BDY) 
     ;WRITE(*,906) win%NSAMP_BDY,win%DT_BDY,win%START_BDY
 905  FORMAT(' rdwin: Samples in surface wave synthetics:',I6,/, 
     ;       '        At sampling rate:',F12.4,/, 
     ;       '        With start time:',F8.3,/) 
 906  FORMAT(' rdwin: Samples in body wave synthetics:',I6,/,
     ;       '        At sampling rate:',F12.4,/,
     ;       '        With start time:',F8.3,/)
      RETURN
   60 CONTINUE
      CLOSE(IUNIT) 
      WRITE(*,*) 'rdwin: Premature end of window file!'
      IERR = 1
      CLOSE(60) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDEEST_HD(PROJNM, NDIM,NFREQ,NREC,NSRC, IERR) 
!
!     Reads the frequency response estimate seismogram headers
!
!     INPUT      MEANING
!     -----      ------- 
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     NDIM       number of componenets in file
!     NFREQ      number of frequencies in file
!     NREC       number of receivers in file
!     NSRC       number of sources in file
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(OUT) :: NDIM,NFREQ,NREC,NSRC, IERR 
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FLNAME
      LOGICAL*4 LEX, LSWAP
      PARAMETER(IUNIT = 44) 
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(TOL = 1.11E-5)
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0 
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... file detection
      FLNAME(1:80) = ' '
      FLNAME = './dest/'//TRIM(PROJNM)//'.edest'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=FLNAME,EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdeest25_hd: Error cannot detect estimate file'
         IERR = 1
         RETURN
      ENDIF
!
!.... open and read header
      NBYTES = 16
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,NBYTES)
      CLOSE(IUNIT)
      IF (IERR.NE.0) THEN
         NFREQ = 0
         NSRC = 0
         NREC = 0
         NDIM = 0
         GOTO 500
      ENDIF
!
!.... unpack the header
      NFREQ = UNPACKI4(LSWAP,CDAT( 1:4))
      NSRC  = UNPACKI4(LSWAP,CDAT( 5:8))
      NREC  = UNPACKI4(LSWAP,CDAT( 9:12))
      NDIM  = UNPACKI4(LSWAP,CDAT(13:16))
  500 CONTINUE
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDEEST_FR(PROJNM, NFREQ,FREQ4,IERR) 
!
!     Reads frequency list from frequency response estimate file
!
!     INPUT      MEANING
!     -----      -------
!     NFREQ      number of frequencies
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     FREQ4      frequency list from estimate file
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: NFREQ
      REAL*4, INTENT(OUT) :: FREQ4(NFREQ)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FLNAME
      LOGICAL*4 LSWAP, LEX
      PARAMETER(IUNIT = 44)
      REAL*4 UNPACKR4
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(TOL = 1.11E-5)

!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0 
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... file detection
      FLNAME(1:80) = ' ' 
      FLNAME = './dest/'//TRIM(PROJNM)//'.edest' 
      FLNAME = ADJUSTL(FLNAME)
!
!.... open and read
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN
         IERR = 1
         WRITE(*,*) 'rdeest_fr: Error file does not exist',TRIM(FLNAME)
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
      NADV = 2*4*NDIMIN*NSRCIN*NRECIN
!
!.... search for frequency in frequency list
      INDX = 16
      DO 1 IFREQL=1,NFREQIN 
         FREQ4(IFREQL) = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
         INDX = INDX + 4 + NADV
    1 CONTINUE 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDEEST25(PROJNM, MFREQL,MREC,NFREQL, NDIM,NREC,NSRC, 
     ;                    FREQL, UEST,VEST,WEST, IERR) 
!
!     Reads the frequency domain estimates at seismogram stations 
!     for plotting by a xwiggle25
! 
!     INPUT      MEANING
!     -----      ------- 
!     FREQL      frequency list
!     MFREQL     leading dimension for ?est
!     MREC       leading dimension for ?est
!     NDIM       number of spatial dimensions expected, should be 3
!     NFREQL     number of frequencies in list  
!     NREC       number of recievers we are expecting
!     NSRC       number of sources we are expecting
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     UEST       u estimate at receivers
!     VEST       v estimate at receivers
!     WEST       w estimate at receivers
! 
!.... variable declarations 
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*4, INTENT(IN) :: FREQL(NFREQL) 
      INTEGER*4, INTENT(IN) :: MFREQL,MREC,NFREQL, NDIM,NREC,NSRC
      COMPLEX*8, INTENT(OUT) :: UEST(MFREQL,MREC,*), 
     ;                          VEST(MFREQL,MREC,*), WEST(MFREQL,MREC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:) 
      CHARACTER(80) FLNAME
!     INTEGER*4 INFO(13)
      LOGICAL*4 LSWAP, LEX
      PARAMETER(IUNIT = 44) 
      REAL*4 UNPACKR4
      INTEGER*4 UNPACKI4, ENDIAN 
      PARAMETER(TOL = 1.11E-5) 
!
!----------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... file detection
      FLNAME(1:80) = ' '
      FLNAME = './dest/'//TRIM(PROJNM)//'.edest' 
      FLNAME = ADJUSTL(FLNAME)
!
!.... open and read
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN 
         IERR = 1
         WRITE(*,*) 'rdeest25: Error file does not exist',TRIM(FLNAME)
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
      IF (NRECIN.NE.NREC) 
     ;WRITE(*,*) 'rdeest25: Warning sources do not match'
      IF (NSRCIN.NE.NSRC)
     ;WRITE(*,*) 'rdeest25: Warning sources do not match'
      IF (NDIMIN.NE.NDIM)
     ;WRITE(*,*) 'rdeest25: Warning dimensions do not match'
      NADV = 2*4*NDIMIN*NSRCIN*NRECIN
!
!.... search for frequency in frequency list
      DO 1 IFREQL=1,NFREQL 
         INDX = 16
         DO 2 JFREQL=1,NFREQL
            FREQIN = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
            INDX = INDX + 4
            IF (ABS(FREQIN - FREQL(IFREQL)).LT.TOL) THEN
               GOTO 20 
            ELSE
               INDX = INDX + NADV
            ENDIF
   2     CONTINUE
         WRITE(*,*) 'rdeest25: Error could not locate frequency'
         IERR = 1
         RETURN
  20     CONTINUE !got it
         DO 3 ISRC=1,NSRCIN 
            DO 4 IREC=1,NRECIN
               R1 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
               Q1 = UNPACKR4(LSWAP,CDAT(INDX+5:INDX+8))
               INDX = INDX + 8
               R2 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
               Q2 = UNPACKR4(LSWAP,CDAT(INDX+5:INDX+8))
               INDX = INDX + 8 
               R3 = UNPACKR4(LSWAP,CDAT(INDX+1:INDX+4))
               Q3 = UNPACKR4(LSWAP,CDAT(INDX+5:INDX+8))
               INDX = INDX + 8
               UEST(IFREQL,IREC,ISRC) = CMPLX(R1,Q1)
               VEST(IFREQL,IREC,ISRC) = CMPLX(R2,Q2)
               WEST(IFREQL,IREC,ISRC) = CMPLX(R3,Q3)
    4      CONTINUE !loop on recievers
    3   CONTINUE !looop on sources
    1 CONTINUE !loop on frequency list 
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE WTSEISM(PROJNM, MDIM,MSAMP,MREC, NDIM,NSAMP,NREC,NSRC,
     ;                   DT,STARTT, WIGGLE, IERR) 
!
!     Writes seismogram data in a format Gphase can read
!
!     INPUT      MEANING
!     -----      ------- 
!     DT         sampling interval (seconds)
!     MDIM       leading dimension for wiggle
!     MREC       leading dimension for wiggle
!     MSAMP      leading dimension for wiggle
!     NDIM       number of components
!     NREC       number of receivers
!     NSAMP      number of samples
!     NSRC       number of sources
!     PROJNM     project name
!     STARTT     start time, probably 0
!     WIGGLE     2 or 3 component time domain seismic responses
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*4, INTENT(IN) :: WIGGLE(MDIM,MSAMP,MREC,*), DT,STARTT 
      INTEGER*4, INTENT(IN) :: MDIM,MSAMP,MREC, NDIM,NSAMP,NREC,NSRC
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(4), ALLOCATABLE :: CDAT(:) 
      CHARACTER(80) FLNAME
      LOGICAL*4 LSWAP
      CHARACTER(4) PACKI4,PACKR4
      INTEGER*4 ENDIAN
      PARAMETER(IUNIT = 45) 
!                                                                      !
!----------------------------------------------------------------------!
! 
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... pack header and model      
      LWORK = NDIM*NREC*NSRC*NSAMP + 5
      NBYTES = LWORK*4
      ALLOCATE(CDAT(LWORK))
      CDAT(1) = PACKI4(LSWAP,NSAMP)
      CDAT(2) = PACKI4(LSWAP,NSRC)
      CDAT(3) = PACKI4(LSWAP,NREC)
      DTMS = DT*1000.0 
      CDAT(4) = PACKR4(LSWAP,DTMS)
      CDAT(5) = PACKR4(LSWAP,STARTT)
      INDX = 5
      DO 1 ISRC=1,NSRC
         DO 2 IREC=1,NREC
            DO 3 I=1,NDIM
               DO 4 K=1,NSAMP
                  INDX = INDX + 1
                  CDAT(INDX) = PACKR4(LSWAP,WIGGLE(I,K,IREC,ISRC))
    4          CONTINUE
               IF (NDIM.EQ.2 .AND. I.EQ.1) THEN
                  DO 5 K=1,NSAMP
                     INDX = INDX + 1
                     CDAT(INDX) = PACKR4(LSWAP,0.)
    5             CONTINUE
               ENDIF
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE 
      IF (INDX*4.NE.NBYTES) THEN
         WRITE(*,*) 'wtseis: Error packing data'
         IERR = 1
      ENDIF
!
!.... write file
      FLNAME = TRIM(PROJNM)//'.pest'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED', 
     ;     STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR) 
      WRITE(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,LWORK)
      IF (IERR.NE.0) WRITE(*,*) 'wtseism: Error writing siesmograms'
      CLOSE(IUNIT) 
      DEALLOCATE(CDAT)
      RETURN
      END

