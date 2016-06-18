!----------------------------------------------------------------------!
!                     Surface wave IO routines                         !
!----------------------------------------------------------------------!
      SUBROUTINE READSF_HD(PROJNM,ISRCIN, NNPE,ISRC,NFREQ)
!
!     Reads the header off the Green's function file
      CHARACTER(80), INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: ISRCIN
      INTEGER*4, INTENT(OUT) :: NNPE, ISRC, NFREQ
      CHARACTER(1), ALLOCATABLE :: CDAT(:) 
      CHARACTER(80) FLNAME
      CHARACTER(3) CSRC
      !INTEGER*4 INFO(13)
      LOGICAL*4 LSWAP, LEX
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(IUNIT = 31) 
      PARAMETER(LHEAD = 12) !lenght of header
!
!----------------------------------------------------------------------!
!
!.... detect endiannes
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... open and read file
      IERR = 0
      FLNAME(1:80) = ' '
      WRITE(CSRC,'(I3)') ISRCIN
      CSRC = ADJUSTL(CSRC)
      FLNAME = 'GRNS/'//TRIM(PROJNM)//'_'//TRIM(CSRC)//'.grns'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'readsf_hd: Cannot find Greens function file'
         IERR = 1
         RETURN
      ENDIF
      !CALL STAT(TRIM(FLNAME),INFO)
      !NBYTES = INFO(8)
      NBYTES = LHEAD 
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(I),I=1,NBYTES)
      CLOSE(IUNIT)
!
!.... unpack header
      NNPE = UNPACKI4(LSWAP,CDAT(1:4))
      ISRC = UNPACKI4(LSWAP,CDAT(5:8))
      NFREQ= UNPACKI4(LSWAP,CDAT(9:12))
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE READSF_XZ(PROJNM, NNPE,ISRC, XLOCSE,ZLOCSE)
!
!     Reads the header off the Green's function file
      CHARACTER*80, INTENT(IN) :: PROJNM
      INTEGER*4, INTENT(IN) :: NNPE, ISRC 
      REAL*8, INTENT(OUT) :: XLOCSE(NNPE), ZLOCSE(NNPE)
      CHARACTER*1, ALLOCATABLE :: CDAT(:)
      CHARACTER*80 FLNAME
      CHARACTER*3 CSRC
      INTEGER*4 INFO(13)
      LOGICAL*4 LSWAP, LEX
      REAL*8 UNPACKR8
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(IUNIT = 31) 
!
!----------------------------------------------------------------------!
!
!.... detect endiannes
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... open and read file
      IERR = 0
      FLNAME(1:80) = ' '
      WRITE(CSRC,'(I3)') ISRC
      CSRC = ADJUSTL(CSRC)
      FLNAME = 'GRNS/'//TRIM(PROJNM)//'_'//TRIM(CSRC)//'.grns'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES) 
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'readsf_xz: Cannot find Greens function file'
         IERR = 1
         RETURN
      ENDIF
!     CALL STAT(TRIM(FLNAME),INFO)
!     NBYTES = INFO(8)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(I),I=1,NBYTES)
      CLOSE(IUNIT)
!
!.... unpack header
      NNPEIN =  UNPACKI4(LSWAP,CDAT(1:4))
      IF (NNPEIN.NE.NNPE) THEN
         WRITE(*,*) 'readsf_xz: Warning NNPEIN != NNPE'
      ENDIF
!
!.... extract (x,z) locations
      INDX = 8 
      DO 1 INPE=1,NNPEIN
         XLOCSE(INPE) = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
         INDX = INDX + 8 
         ZLOCSE(INPE) = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
         INDX = INDX + 8 
    1 CONTINUE  
      DEALLOCATE(CDAT)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE READSF(PROJNM, NNPE,NFREQ, ISRC, FREQ,
     ;                  USPEC,VSPEC,WSPEC, IERR)
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       desired frequency (Hz)
!     PROJNM     project name 
!      
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     ?SPEC      (u,v,w) Green's functions at all nodes in E domain 
!                at desired frequency
!
!.... variable declaraitons 
      CHARACTER(80), INTENT(IN) :: PROJNM
      REAL*8, INTENT(IN) :: FREQ 
      INTEGER*4, INTENT(IN) :: ISRC
      COMPLEX*16, INTENT(OUT) :: USPEC(NNPE), VSPEC(NNPE), WSPEC(NNPE) 
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:)
      CHARACTER(80) FLNAME 
      CHARACTER(3) CSRC 
      REAL*8 FREQIN,R,Q,TOL 
      INTEGER*4 INFO(13) 
      LOGICAL*4 LSWAP, LKEEP, LEX 
      PARAMETER(IUNIT = 31) 
      REAL*8 UNPACKR8
      INTEGER*4 UNPACKI4, ENDIAN
      PARAMETER(TOL = 1.11D-7) 
!
!----------------------------------------------------------------------!
!
!.... estimate size
      LWORK = 4*3              !header
     ;      + 8*2*NNPE         !x and z locations
     ;      + 8*NFREQ          !frequencies
     ;      + 8*NFREQ*6*NNPE   !(u,v,w) responses at nodal points 
!
!.... detect endiannes
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND.NE.0) LSWAP = .TRUE.
!
!.... open and read file
      IERR = 0
      FLNAME(1:80) = ' ' 
      WRITE(CSRC,'(I3)') ISRC
      CSRC = ADJUSTL(CSRC)
      FLNAME = 'GRNS/'//TRIM(PROJNM)//'_'//TRIM(CSRC)//'.grns'
      FLNAME = ADJUSTL(FLNAME)
      INQUIRE(FILE=TRIM(FLNAME),EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'readsf: Cannot find Greens function file'
         IERR = 1
         RETURN
      ENDIF 
!     CALL STAT(TRIM(FLNAME),INFO) 
!     NBYTES = INFO(8) 
      IF (NBYTES.NE.LWORK) THEN
         WRITE(*,*) 'readsf: lworks != nbytes',LWORK,NBYTES
         IERR = 2
         RETURN
      ENDIF 
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED', 
     ;     ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(I),I=1,NBYTES)
      CLOSE(IUNIT) 
!
!.... unpack header
      NNPEIN  = UNPACKI4(LSWAP,CDAT(1:4)) 
      NFREQIN = UNPACKI4(LSWAP,CDAT(9:12))
      INDX = 12 + 8*2*NNPE  
!
!.... unpack responses
      DO 3 IFREQ=1,NFREQIN
         LKEEP = .FALSE.
         FREQIN = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8)) 
         INDX = INDX + 8
         IF (DABS(FREQIN - FREQ).LT.TOL) LKEEP = .TRUE.
         IF (LKEEP) THEN
            DO 4 INPE=1,NNPEIN 
               R = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
               INDX = INDX + 8
               Q = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
               USPEC(INPE) = DCMPLX(R,Q)
               INDX = INDX + 8

               R = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
               INDX = INDX + 8
               Q = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
               VSPEC(INPE) = DCMPLX(R,Q)
               INDX = INDX + 8

               R = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
               INDX = INDX + 8
               Q = UNPACKR8(LSWAP,CDAT(INDX+1:INDX+8))
               WSPEC(INPE) = DCMPLX(R,Q)
               INDX = INDX + 8
    4       CONTINUE
            GOTO 50
         ELSE
            INDX = INDX + 6*8*NNPE
         ENDIF
    3 CONTINUE  
      WRITE(*,*) 'readsf: Error could not locate frequency:',FREQ 
      IERR = 1
   50 CONTINUE
      DEALLOCATE(CDAT) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE WRITESF(MNPE,PROJNM, NNPE,NFREQ, ISRC, 
     ;                   XLOCSE,ZLOCSE, FREQ, USPEC,VSPEC,WSPEC, IERR) 
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       frequency list (Hz)
!     ISRC       source number
!     MNPE       leading dimension for uspec,vspec,wspec
!     NFREQ      number of frequencies
!     NNPE       number of nodal points in E layer
!     PROJNM     project name
!     ?SPEC      (u,v,w) Green's funs at all nodes for all frequencies
!     XLOCSE     x locations of nodes in E domain (m)
!     ZLOCSE     z locations of nodes in E domain (m)
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!
!.... variable declarations 
      CHARACTER(80), INTENT(IN) :: PROJNM
      COMPLEX*16, INTENT(IN) :: USPEC(MNPE,*), VSPEC(MNPE,*), 
     ;                          WSPEC(MNPE,*)
      REAL*8, INTENT(IN) :: XLOCSE(NNPE), ZLOCSE(NNPE), FREQ(NFREQ) 
      INTEGER*4, INTENT(IN) :: MNPE, NNPE, NFREQ, ISRC
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      CHARACTER(80) FLNAME
      CHARACTER(8) C8  
      CHARACTER(4) C4
      CHARACTER(1) CV8(8), CV4(4) 
      CHARACTER(3) CSRC
      CHARACTER(1), ALLOCATABLE :: CDAT(:) 
      LOGICAL*4 LSWAP, LEX 
      CHARACTER(8) PACKR8
      CHARACTER(4) PACKI4 
      INTEGER*4 ENDIAN
      PARAMETER(IUNIT = 31) 
      EQUIVALENCE(CV8,C8) 
      EQUIVALENCE(CV4,C4)
! 
!----------------------------------------------------------------------!
!
!.... calculate size and machine endianness
      LWORK = 4*3              !header
     ;      + 8*2*NNPE         !x and z locations
     ;      + 8*NFREQ          !frequencies
     ;      + 8*NFREQ*6*NNPE   !(u,v,w) responses at nodal points 
      ALLOCATE(CDAT(LWORK))
      CDAT(1:LWORK) = ' '
      MYEND = ENDIAN()
      LSWAP = .FALSE.
      IF (MYEND.NE.0) LSWAP = .TRUE.
!.... pack header
      C4 = PACKI4(LSWAP,NNPE)
      CDAT(1:4) = CV4(1:4)
      C4 = PACKI4(LSWAP,ISRC)
      CDAT(5:8) = CV4(1:4)
      C4 = PACKI4(LSWAP,NFREQ)
      CDAT(9:12)= CV4(1:4)
      INDX = 12
!.... pack (x,z) locations
      DO 1 INPE=1,NNPE
         C8 = PACKR8(LSWAP,XLOCSE(INPE))
         CDAT(INDX+1:INDX+8) = CV8(1:8)
         INDX = INDX + 8
         C8 = PACKR8(LSWAP,ZLOCSE(INPE))
         CDAT(INDX+1:INDX+8) = CV8(1:8)
         INDX = INDX + 8
    1 CONTINUE
!
!.... loop on frequencies and pack responses
      DO 2 IFREQ=1,NFREQ
         C8 = PACKR8(LSWAP,FREQ(IFREQ))
         CDAT(INDX+1:INDX+8) = CV8(1:8)
         INDX = INDX + 8
         DO 3 INPE=1,NNPE
            C8 = PACKR8(LSWAP,DREAL(USPEC(INPE,IFREQ)))
            CDAT(INDX+1:INDX+8) = CV8(1:8)
            INDX = INDX + 8
            C8 = PACKR8(LSWAP,DIMAG(USPEC(INPE,IFREQ)))
            CDAT(INDX+1:INDX+8) = CV8(1:8)
            INDX = INDX + 8

            C8 = PACKR8(LSWAP,DREAL(VSPEC(INPE,IFREQ)))
            CDAT(INDX+1:INDX+8) = CV8(1:8)
            INDX = INDX + 8
            C8 = PACKR8(LSWAP,DIMAG(VSPEC(INPE,IFREQ)))
            CDAT(INDX+1:INDX+8) = CV8(1:8)
            INDX = INDX + 8

            C8 = PACKR8(LSWAP,DREAL(WSPEC(INPE,IFREQ)))
            CDAT(INDX+1:INDX+8) = CV8(1:8)
            INDX = INDX + 8
            C8 = PACKR8(LSWAP,DIMAG(WSPEC(INPE,IFREQ)))
            CDAT(INDX+1:INDX+8) = CV8(1:8)
            INDX = INDX + 8
    3    CONTINUE
    2 CONTINUE !loop on frequencies
!
!.... open and write file 
      NBYTES = INDX
      IF (LWORK.NE.NBYTES) THEN
         WRITE(*,*) 'writesf: Error nbytes != lwork',NBYTES,LWORK
         IERR = 1
         GOTO 50
      ENDIF
      FLNAME(1:80) = ' '
      WRITE(CSRC,'(I3)') ISRC
      CSRC = ADJUSTL(CSRC)
      INQUIRE(FILE='GRNS',EXIST=LEX)
      IF (.NOT.LEX) CALL SYSTEM('mkdir GRNS')
      FLNAME = 'GRNS/'//TRIM(PROJNM)//'_'//TRIM(CSRC)//'.grns'
      FLNAME = ADJUSTL(FLNAME)
      OPEN(UNIT=IUNIT,FILE=TRIM(FLNAME),FORM='UNFORMATTED',
     ;     STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      WRITE(IUNIT,REC=1) (CDAT(J),J=1,NBYTES)
      CLOSE(IUNIT)
   50 CONTINUE
      DEALLOCATE(CDAT)
      RETURN
      END
