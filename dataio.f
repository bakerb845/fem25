      SUBROUTINE EWAVOUT2D(PROJNM,MDIM,MEN,MGNOD,NNPG,NDOF, NELEM,
     ;                     NDIM,NLXI,NLETA,NGNOD,ISRC,FREQ,
     ;                     IENG,LM, XLOCS,ZLOCS,SOL,IERR)
! 
!     Writes solution at anchor nodes.  Solution is single precision
! 
!     INPUT      MEANING
!     -----      ------- 
!     IENG       global anchor node pointer 
!     ISRC       current source
!     LM         location matrix 
!     FREQ       current frequency 
!     MDIM       max number of spatial dimensions
!     MEN        max number of element nodes 
!     MGNOD      max number of anchor nodes
!     NDIM       number of spatial dimensions
!     NDOF       number of degrees of freedom
!     NGNOD      number of anchor nodes on element
!     NLETA      number of lagrange interpolant points in eta 
!     NLXI       number of lagrange interpolant points in xi 
!     NNPG       number of global nodal points
!     PROJNM     project name
!     SOL        solution vector for this source
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!     
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM
      COMPLEX*8 SOL(NDOF)
      REAL*8 XLOCS(NNPG), ZLOCS(NNPG), FREQ
      INTEGER*4 LM(MDIM,MEN,*), IENG(MGNOD,*), ISRC
!.... local variables
      CHARACTER(80) FILENM
      CHARACTER(12) CFREQ
      CHARACTER(5) CSRC
      CHARACTER(4), ALLOCATABLE :: DAT(:)
      REAL*4 RESPR,RESPI
      INTEGER*4, ALLOCATABLE :: IDPTR(:,:), IENGNOD(:)
      LOGICAL*4 LSWAP,LEX,LISDIR
      CHARACTER(4) PACKI4,PACKR4
      INTEGER*4 ENDIAN
      PARAMETER(IUNIT = 31, LNULL =-5)
! 
!----------------------------------------------------------------------!
! 
!.... determine endianness, default write to little endian
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
      NIENGV = 0
      ALLOCATE(IENGNOD(NELEM))
      NIENGV = 0
      DO 1 IELEM=1,NELEM
         DO 2 IA=1,NGNOD
            IF (IENG(IA,IELEM).EQ.0) THEN !triangle
               IENGNOD(IELEM) = 3
               GOTO 20
            ENDIF
    2    CONTINUE
         IENGNOD(IELEM) = 4
   20    CONTINUE
    1 CONTINUE
! 
!.... create pointer to extract response
      ALLOCATE(IDPTR(NNPG,NDIM))
      CALL GENIDPTR(NNPG,MGNOD,MEN,MDIM, NNPG,NELEM,NDIM,
     ;              NLXI,NLETA,LNULL, IENGNOD,IENG,LM, IDPTR)
! 
!.... pack header 
      LRECS = 4 + 2*NNPG + NGNOD*NELEM + NDIM*NNPG*2
      ALLOCATE(DAT(LRECS))
      ILOC = 1
      DAT(ILOC) = PACKI4(LSWAP,NNPG)
      ILOC = ILOC + 1
      DAT(ILOC) = PACKI4(LSWAP,NDIM)
      ILOC = ILOC + 1
      DAT(ILOC) = PACKI4(LSWAP,NELEM)
      ILOC = ILOC + 1
      DAT(ILOC) = PACKI4(LSWAP,NGNOD)
! 
!.... x,z locations
      DO 6 INPG=1,NNPG
         ILOC = ILOC + 1
         DAT(ILOC) = PACKR4(LSWAP,SNGL(XLOCS(INPG)))
         ILOC = ILOC + 1
         DAT(ILOC) = PACKR4(LSWAP,SNGL(ZLOCS(INPG)))
    6 CONTINUE
! 
!.... mesh
      DO 7 IELEM=1,NELEM
         DO 8 IA=1,NGNOD
            ILOC = ILOC + 1
            DAT(ILOC) = PACKI4(LSWAP,IENG(IA,IELEM)-1)
    8    CONTINUE
    7 CONTINUE
! 
!.... response
      DO 9 INPG=1,NNPG
         DO 10 I=1,NDIM
            IF (IDPTR(INPG,I).EQ.LNULL) THEN
               WRITE(*,*) 'ewavout25: Error creating idptr'
               IERR = 1
               RETURN
            ENDIF
            IF (IDPTR(INPG,I).EQ.0) THEN
               RESPR = 0.
               RESPI = 0.
            ELSE
               RESPR = REAL(CMPLX(SOL(IDPTR(INPG,I))))
               RESPI = IMAG(CMPLX(SOL(IDPTR(INPG,I))))
            ENDIF
            ILOC = ILOC + 1
            DAT(ILOC) = PACKR4(LSWAP,RESPR)
            ILOC = ILOC + 1
            DAT(ILOC) = PACKR4(LSWAP,RESPI)
   10    CONTINUE
    9 CONTINUE
      DEALLOCATE(IDPTR)
      NBYTES = 4*LRECS
! 
!.... filename handling
      CALL NULLS(80,FILENM)
      CALL NULLS(12,CFREQ)
      CALL NULLS(5,CSRC)
      WRITE(CFREQ,'(F12.5)') FREQ
      CFREQ = ADJUSTL(CFREQ)
      WRITE(CSRC,'(I5)') ISRC
      CSRC = ADJUSTL(CSRC)
      LEX = LISDIR('./wave') 
!#ifdef INTEL
!      INQUIRE(DIRECTORY='./wave',EXIST=LEX)
!#else
!      INQUIRE(FILE='./wave',EXIST=LEX)
!#endif
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./wave')
      FILENM = './wave/'//TRIM(ADJUSTL(PROJNM))//'_wave-'//TRIM(CFREQ)//
     ;         '-'//TRIM(CSRC)//'.ewav'
      CALL DBLANK(80,FILENM,LENFL)
      OPEN(UNIT=IUNIT,FILE=FILENM(1:LENFL),FORM='UNFORMATTED',
     ;     STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      WRITE(IUNIT,REC=1) (DAT(J),J=1,LRECS)
      CLOSE(IUNIT)
      DEALLOCATE(DAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE EWAVHD25D(FLNAME, LSWAP,NNPG,NDIM,NELEM,NGNOD,IERR)
! 
!     Reads header off a .ewav file header
! 
!     INPUT      MEANING
!     -----      -------
!     FLNAME     name of file to read 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag 
!     LSWAP      True -> need to swap bytes in ewavin  
!     NDIM       number of spatial dimensions 
!     NELEM      number of elements 
!     NGNOD      number of anchor nodes 
!     NNPG       number of points in global geometry 
! 
!.... variable declarations
      CHARACTER(80) FLNAME
      CHARACTER(4), ALLOCATABLE :: DAT(:)
      LOGICAL*4 LSWAP
      INTEGER*4 UNPACKI4,ENDIAN
      PARAMETER(IUNIT = 80)
! 
!----------------------------------------------------------------------!
!
! 
!.... determine endianness, default write to little endian
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
! 
!.... open and read file
      IERR = 0
      LRECS = 4
      ALLOCATE(DAT(LRECS))
      NBYTES = LRECS*4
      CALL DBLANK(80,FLNAME,LENFL)
      OPEN(UNIT=IUNIT,FILE=FLNAME(1:LENFL),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      READ(IUNIT,REC=1,IOSTAT=IERR) (DAT(I),I=1,LRECS)
      CLOSE(IUNIT)
!   
!.... parse header
      LSWAP = .FALSE.
   60 CONTINUE
      ILOC = 1
      NNPG = UNPACKI4(LSWAP,DAT(ILOC))
      ILOC = ILOC + 1
      NDIM = UNPACKI4(LSWAP,DAT(ILOC))
      print *, nnpg,ndim
      IF (NDIM.NE.3) THEN
         IF (.NOT.LSWAP) THEN
            LSWAP = .TRUE.
         ELSE
            WRITE(*,*) 'ewavehd25d: cant parse header'
            IERR = 1
            RETURN
         ENDIF
         GOTO 60
      ENDIF
      ILOC = ILOC + 1
      NELEM = UNPACKI4(LSWAP,DAT(ILOC))
      ILOC = ILOC + 1
      NGNOD = UNPACKI4(LSWAP,DAT(ILOC))
      DEALLOCATE(DAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE EWAVIN2D(MNPG,MGNOD, FLNAME,LSWAP,NNPG,NDIM,NELEM,
     ;                    NGNOD, IENG,XLOCS4,ZLOCS4,WAVE8)
! 
!     reads the complex wavefield 
! 
!     INPUT      MEANING
!     -----      ------- 
!     FLNAME     name of input file 
!     LSWAP      True -> swap bytes 
!     MGNOD      leading dimension for ieng
!     MNPG       leading dimension for wave8
!     NDIM       number of spatial dimensions 
!     NELEM      number of elements 
!     NGNOD      number of anchor nodes 
!     NNPG       number of points in global geometry 
!     
!     OUTPUT     MEANING
!     ------     ------- 
!     IENG       global geometry pointer 
!     WAVE8      complex*8 complex response at anchor nodes 
!     XLOCS4     real*4 x locations at anchor nodes 
!     ZLOCS4     real*4 z locations at anchor nodes 
! 
!.... variable declarations
      CHARACTER(80) FLNAME
      COMPLEX*8 WAVE8(MNPG,NDIM)
      REAL*4 XLOCS4(NNPG), ZLOCS4(NNPG)
      INTEGER*4 IENG(MGNOD,NELEM)
!.... local variables
      CHARACTER(4), ALLOCATABLE :: DAT(:)
      REAL*4 UNPACKR4
      LOGICAL*4 LSWAP
      INTEGER*4 UNPACKI4,ENDIAN
      PARAMETER(IUNIT = 84)
! 
!----------------------------------------------------------------------!
! 
!.... determine endianness, default write to little endian
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... read file 
      LRECS = 4 + 2*NNPG + NGNOD*NELEM + NDIM*NNPG*2
      ALLOCATE(DAT(LRECS))
      NBYTES = LRECS*4
      CALL DBLANK(80,FLNAME,LENFL)
      OPEN(UNIT=IUNIT,FILE=FLNAME(1:LENFL),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      READ(IUNIT,REC=1,IOSTAT=IERR) (DAT(I),I=1,LRECS)
      CLOSE(IUNIT)
! 
!.... parse header
      ILOC = 1
      NNPGIN = UNPACKI4(LSWAP,DAT(ILOC))
      ILOC = ILOC + 1
      NDIMIN = UNPACKI4(LSWAP,DAT(ILOC))
      ILOC = ILOC + 1
      NELEMIN = UNPACKI4(LSWAP,DAT(ILOC))
      ILOC = ILOC + 1
      NGNODIN = UNPACKI4(LSWAP,DAT(ILOC))
! 
!.... x,z locations
      DO 6 INPG=1,NNPG
         ILOC = ILOC + 1
         XLOCS4(INPG) = UNPACKR4(LSWAP,DAT(ILOC))
         ILOC = ILOC + 1
         ZLOCS4(INPG) = UNPACKR4(LSWAP,DAT(ILOC))
    6 CONTINUE
! 
!.... mesh
      DO 7 IELEM=1,NELEM
         DO 8 IA=1,NGNOD
            ILOC = ILOC + 1
            IENG(IA,IELEM) = UNPACKI4(LSWAP,DAT(ILOC)) + 1
    8    CONTINUE
    7 CONTINUE
! 
!.... response
      DO 9 INPG=1,NNPG
         DO 10 I=1,NDIM
            ILOC = ILOC + 1
            RESPR = UNPACKR4(LSWAP,DAT(ILOC))
            ILOC = ILOC + 1
            RESPI = UNPACKR4(LSWAP,DAT(ILOC))
            WAVE8(INPG,I) = CMPLX(RESPR,RESPI)
   10    CONTINUE
    9 CONTINUE
      DEALLOCATE(DAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDREC_RESP(PROJNM, MDIM,MFREQ,NDIM,NFREQ,NREC, AZMOD, 
     ;                      FREQ, LRECEX,RECV)
!
!     Reads the receiver frequency response file
!
!     INPUT      MEANING
!     -----      -------
!     AZMOD      model azimuth (degrees)
!     FREQ       frequency list to read in
!     LRECEX     False -> .wrec file does not exist
!     MDIM       leading dimension for RECV
!     MFREQ      leading dimension for RECV
!     NDIM       number of components on receivers
!     NFREQ      number of frequencies to read in
!     NREC       number of receivers
!     PROJNM     project name
!     
!     OUTPUT     MEANING
!     ------     -------
!     LRECEX     False -> .wrec file does not exist or is invalid
!     RECV       receiver response function (u,v,w) frame
!
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: FREQ(NFREQ), AZMOD
      INTEGER*4, INTENT(IN) :: MDIM, NDIM, NFREQ, NREC
      COMPLEX*8, INTENT(OUT) :: RECV(MDIM,MFREQ,*)
      LOGICAL*4, INTENT(OUT) :: LRECEX 
!.... local variables
      CHARACTER(80) RECFL
      CHARACTER(1), ALLOCATABLE :: DAT(:)
      COMPLEX*8 RRFNEZ(3) 
      REAL*8 PI180 
      REAL*4 Q, R
      INTEGER*4 IERR
      REAL*4 UNPACKR4
      INTEGER*4 UNPACKI4, ENDIAN
      LOGICAL*4 LEX, LSWAP
      PARAMETER(TOL = 1.0E-4)
      PARAMETER(LENHEAD = 3) !number of header values
      PARAMETER(IUNIT = 44) 
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
!
!----------------------------------------------------------------------!
!
!.... check if file exists
      IERR = 0 
      RECFL(1:80) = ' ' 
      !RECFL = './recrsp'//TRIM(PROJNM)//'.wrec'
      RECFL = TRIM(PROJNM)//'.wrec'
      RECFL = ADJUSTL(RECFL)
      INQUIRE(FILE=TRIM(RECFL),EXIST=LEX,SIZE=NBYTES)
      IF (LEX) THEN
         LRECEX = .TRUE.
         MYEND = ENDIAN()
         LSWAP = .FALSE.
         IF (MYEND.NE.0) LSWAP = .TRUE.
         IERR = 0
         !CALL STAT(TRIM(RECFL),INFO)
         !NBYTES = INFO(8)
         LRECS = NBYTES/4
         ALLOCATE(DAT(NBYTES))
         OPEN(UNIT=IUNIT,FILE=TRIM(RECFL),FORM='UNFORMATTED',
     ;        ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
         READ(IUNIT,REC=1,IOSTAT=IERR) (DAT(I),I=1,NBYTES)
         CLOSE(IUNIT)
         IF (IERR.NE.0) THEN
            WRITE(*,*) 'rdrec_resp: Error reading wrec file'
            WRITE(*,*) 'rdrec_resp: Overriding receiver response'
            GOTO 100
         ENDIF
!
!....... read the header
         NFREQ_IN = UNPACKI4(LSWAP,DAT(1:4))
         NREC_IN  = UNPACKI4(LSWAP,DAT(5:8))
         NDIM_IN  = UNPACKI4(LSWAP,DAT(9:12))
         IF (NDIM.NE.NDIM_IN) THEN
            WRITE(*,*) 'rdrec_resp: Error component mismatch'
            WRITE(*,*) 'rdrec_resp: Overriding receiver response'
            IERR = 1
            GOTO 100
         ENDIF
         IF (NREC.NE.NREC_IN) THEN
            WRITE(*,*) 'rdrec_resp: Error receiver mismatch'
            WRITE(*,*) 'rdrec_resp: Overriding receiver response'
            IERR = 1
            GOTO 100
         ENDIF
         NADV = 2*NDIM_IN*NREC_IN*4 + 4 !2 for complex numbers
!
!....... loop on my input frequency list
         DO 1 IFREQ=1,NFREQ
!
!.......... locate the correct frequency
            INDX = 13 
            DO 2 JFREQ=1,NFREQ_IN
               FREQ_IN = UNPACKR4(LSWAP,DAT(INDX:INDX+3))
               IF (ABS(FREQ_IN - SNGL(FREQ(IFREQ))).LT.TOL) GOTO 20
               INDX = INDX + NADV
    2       CONTINUE
            WRITE(*,*) 'rdrec_resp: Error cannot locate frequency'
            WRITE(*,*) 'rdrec_resp: Overriding receiver resopnse'
            IERR = 1
            GOTO 100
   20       CONTINUE !got frequency, break ahead
            DO 3 IREC=1,NREC_IN
               DO 4 I=1,NDIM_IN
                  INDX = INDX + 4
                  R = UNPACKR4(LSWAP,DAT(INDX:INDX+3))
                  INDX = INDX + 4
                  Q = UNPACKR4(LSWAP,DAT(INDX:INDX+3))
                  RRFNEZ(I) = CMPLX(R,Q)
    4          CONTINUE !Loop on components 
               !(N,E) -> (u,v) clockwise rotation
               CALL CROTATE(-SNGL(AZMOD), RRFNEZ(1),RRFNEZ(2), 
     ;                      RECV(1,IFREQ,IREC),RECV(2,IFREQ,IREC))
               RECV(3,IFREQ,IREC) = RRFNEZ(3) 
    3       CONTINUE !loop on receivers
    1    CONTINUE !loop on frequencies
      ENDIF
!
!.... override receiver response or just set as unity
  100 CONTINUE
      IF (IERR.NE.0 .OR. .NOT.LEX) THEN
         LRECEX = .FALSE.
         WRITE(*,*) 'rdrec_resp: Setting receiver response to unity...'
         DO 11 IFREQ=1,NFREQ
            DO 12 IREC=1,NREC
               DO 13 I=1,NDIM
                  RECV(I,IFREQ,IREC) = CMPLX(1.0,0.0)
   13          CONTINUE
   12       CONTINUE
   11    CONTINUE
      ENDIF
!
!.... clean space
      IF (ALLOCATED(DAT)) DEALLOCATE(DAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDREC_RESP_HD(PROJNM, NDIM,NFREQ,NREC, IERR)
!
!     Reads the receiver frequency response header
!
!     INPUT      MEANING
!     -----      -------
!     PROJNM     project name
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     NDIM       number of components in wrec file 
!     NFREQ      number of frequencies in wrec file 
!     NREC       number of receivers in wrec file
!
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(OUT) :: NDIM, NFREQ, NREC, IERR
!.... local variables
      CHARACTER(80) RECFL
      INTEGER*4 UNPACKI4, ENDIAN
      LOGICAL*4 LEX, LSWAP
      PARAMETER(LENHEAD = 3) !number of header values
      CHARACTER(4) DAT(LENHEAD)
      PARAMETER(IUNIT = 44)
!
!----------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      NDIM = 0
      NFREQ = 0
      NREC = 0
!
!.... check if file exists
      IERR = 0
      RECFL(1:80) = ' '
      !RECFL = './recrsp'//TRIM(PROJNM)//'.wrec'
      RECFL = TRIM(PROJNM)//'.wrec'
      RECFL = ADJUSTL(RECFL)
      INQUIRE(FILE=TRIM(RECFL),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'rdrec_resp_hd: File does not exist!'
         IERR = 1
         RETURN
      ENDIF
!
!.... read the header
      LRECS = LENHEAD
      NBYTES = 4*LRECS
      OPEN(UNIT=IUNIT,FILE=TRIM(RECFL),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      READ(IUNIT,REC=1,IOSTAT=IERR) (DAT(I),I=1,LRECS)
      CLOSE(IUNIT)
      IF (IERR.NE.0) THEN
         WRITE(*,*) 'rdrec_resp_hd: Error reading file!'
         IERR = 1
         RETURN
      ENDIF
!
!.... read the header
      MYEND = ENDIAN()
      LSWAP = .FALSE.
      IF (MYEND.NE.0) LSWAP = .TRUE.
      NFREQ = UNPACKI4(LSWAP,DAT(1))
      NREC  = UNPACKI4(LSWAP,DAT(2))
      NDIM  = UNPACKI4(LSWAP,DAT(3))
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE RDSTF(PROJNM, MFREQ, NFREQ,NSRC, FREQ,SOURCE, IERR)
!
!     Reads the frequency domain STF function for the forward problem
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       frequency list to read
!     MFREQ      leading dimension for SOURCE
!     NFREQ      number of frequencies
!     NSRC       number of sources
!     PROJNM     project name
!
!     OUTPUT     MEANING 
!     ------     -------
!     IERR       error flag, /= 0 indicates STFs are unity
!     SOURCE     frequency domain source time functions  
!
!.... variable declarations
      CHARACTER(80), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ, NFREQ, NSRC
      COMPLEX*8, INTENT(OUT) :: SOURCE(MFREQ,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FILENM
      LOGICAL*4 LSWAP, LEX
      REAL*4 UNPACKR4 
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(IUNIT = 68)
!
!----------------------------------------------------------------------!
!
!.... set file name and read
      IERR =-1
      FILENM(1:80) = ' ' 
      FILENM = TRIM(PROJNM)//'.wstf'
      FILENM = ADJUSTL(FILENM)
      INQUIRE(FILE=TRIM(FILENM),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) GOTO 100
      !CALL STAT(TRIM(FILENM),INFO) 
      !NBYTES = INFO(8) 
      ALLOCATE(CDAT(NBYTES))
      OPEN(UNIT=IUNIT,FILE=TRIM(FILENM),FORM='UNFORMATTED',
     ;     STATUS='OLD',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,NBYTES)
      CLOSE(IUNIT) 
      IF (IERR.NE.0) THEN
         WRITE(*,*) 'rdstf: Error reading source time function'
         IERR =-1
         GOTO 100  
      ENDIF
!
!.... determine endianness
      MYEND = ENDIAN()
      LSWAP = .FALSE.
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... read the header
      NFREQ_IN = UNPACKI4(LSWAP,CDAT(1:4))
      NSRC_IN  = UNPACKI4(LSWAP,CDAT(5:8))
      IF (NFREQ_IN.NE.NFREQ) 
     ;WRITE(*,*) 'rdstf: Warning frequency number mismatch!'
      IF (NSRC_IN.NE.NSRC) THEN !this is a real error
         WRITE(*,*) 'rdstf: Source number mismatch!'
         IERR =-1
         GOTO 100
      ENDIF
!
!.... unpack the source time function by searching for frequencies
      NADV = 4 + 8*NSRC !frequency + nsrc complex numbers
      DO 1 IFREQ=1,NFREQ
         INDX = 9
         DO 2 IFREQ_IN=1,NFREQ_IN
            FREQ_IN = UNPACKR4(LSWAP,CDAT(INDX:INDX+3))
            IF (ABS(FREQ_IN - SNGL(FREQ(IFREQ))).LT.1.E-4) GOTO 20
            INDX = INDX + NADV 
    2    CONTINUE
         WRITE(*,*) 'rdstf: Could not locate frequency!'
         IERR =-1
         GOTO 100
   20    CONTINUE !got it
         DO 3 ISRC=1,NSRC
            INDX = INDX + 4
            R = UNPACKR4(LSWAP,CDAT(INDX:INDX+3))
            INDX = INDX + 4
            Q = UNPACKR4(LSWAP,CDAT(INDX:INDX+3))
            SOURCE(IFREQ,ISRC) = CMPLX(R,Q)
    3    CONTINUE
    1 CONTINUE !loop on frequencies
      IERR = 0 !success
  100 CONTINUE !break ahead for error
      IF (IERR.EQ.-1) THEN
         DO 30 IFREQ=1,NFREQ
            DO 31 ISRC=1,NSRC
               SOURCE(IFREQ,ISRC) = CMPLX(1.0,0.0)
   31       CONTINUE
   30    CONTINUE
      ENDIF
      IF (ALLOCATED(CDAT)) DEALLOCATE(CDAT)
      RETURN
      END

