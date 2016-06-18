      SUBROUTINE SRCINV(MFREQ,PROJNM, NOMINV,NSRC, FREQINV,
     ;                  LSRCEX,SOURCE)
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
      PARAMETER(TOL = 1.11E-5)
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
