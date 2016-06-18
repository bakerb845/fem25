      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ
      CHARACTER(80) PROJNM,FREQFL,FILENM
      CHARACTER(2) SRCTYP 
      REAL*8 DT8, START8
      REAL*8, ALLOCATABLE :: FREQ_WORK(:)
      PARAMETER(MASTER = 0) 

      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, MPIERR)
      IF (MYID == MASTER) THEN
         CALL NULLS(80,PROJNM)
         CALL NULLS(80,FREQFL) 
         WRITE(*,*) 'xmovie: Enter name of project'
         READ(*,'(A)') PROJNM
         WRITE(*,*) 'xmovie: Enter number of slices'
         !nslices = 300
         READ *, NSLICES
         WRITE(*,*) 'xmovie: Enter source number:'
         READ *, ISRC
         PROJNM = ADJUSTL(PROJNM)
         CALL RDFREQ(PROJNM, FRQ,IERR) 
         IF (IERR /= 0) STOP 
         FILENM(1:80) = ' '
         FILENM = TRIM(PROJNM)//'.src'
         FILENM = ADJUSTL(FILENM)
         OPEN(UNIT=33,FILE=TRIM(FILENM))
         READ(33,*) !header
         READ(33,*) NSRC_IN 
         READ(33,*) !header
         IF (ISRC > NSRC_IN) THEN
            IERR = 1
            WRITE(*,*) 'xmovie: Error source numer too large'
         ENDIF
         SRCTYP(1:2) = ' '
         DO 1 JSRC=1,NSRC_IN
            IF (JSRC == ISRC) THEN
               READ(33,*,IOSTAT=IERR) SRCTYP
               GOTO 45
            ELSE
               READ(33,*,IOSTAT=IERR)
            ENDIF
   1     CONTINUE 
  45     CONTINUE 
         CLOSE(33) 
         FILENM(1:80) = ' '
         IF (SRCTYP(1:1) == 'S') THEN
            FILENM = TRIM(PROJNM)//'_srf' 
            WRITE(*,*) 'xmovie: This is a surface wave source'
         ELSEIF (SRCTYP(1:1) == 'P') THEN
            FILENM = TRIM(PROJNM)//'_bdy'
            WRITE(*,*) 'xmovie: This is a body wave source'
         ELSE
            WRITE(*,*) 'xmovie: Cant determine source type'
            IERR = 1
         ENDIF
         IF (IERR /= 0) CALL MPI_ABORT(MPI_COM_WORLD,30,MPIERR) 
         FILENM = ADJUSTL(FILENM)
         CALL RDSRCT_HD(FILENM, NSRC,NSAMP,DT8,START8,IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xmovie: Error reading source time funciton'
         ENDIF
         TMAX = SNGL( DT8*DFLOAT(NSAMP - 1) )
         DT = TMAX/FLOAT(NSLICES - 1)
         !CALL RDFREQ(FREQFL,2, NFREQT, FREQ)  
         !READ(30,*) TMAX
!
!....... modify the source list
         NFREQ = 0
         ALLOCATE(FREQ_WORK(frq%NFREQ))
         DO 2 IFREQ=1,frq%NFREQ
            IF (SRCTYP(1:1) == 'S' .AND. frq%CFTYPE(IFREQ)(1:1) == 'S' .OR. &
                SRCTYP(1:1) == 'P' .AND. frq%CFTYPE(IFREQ)(1:1) == 'B') THEN
               NFREQ = NFREQ + 1
               FREQ_WORK(NFREQ) = frq%FREQ(IFREQ)
            ENDIF
    2    CONTINUE
         IF (NFREQ == 0) THEN
            WRITE(*,*) 'xmovie: Error no frequencies!'
            CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
         ENDIF
         DEALLOCATE(frq%FREQ)
         frq%NFREQ = NFREQ
         ALLOCATE(frq%FREQ(frq%NFREQ))
         frq%FREQ(1:frq%NFREQ) = FREQ_WORK(1:frq%NFREQ)
         DEALLOCATE(FREQ_WORK) 
         print *, nslices, dt
         print *, dt*(nslices - 1),tmax
         IF (NSLICES*DT > TMAX) NSLICES = NSLICES - 1
      ENDIF
      CALL MPI_BCAST(frq%NFREQ,  1,MPI_INTEGER,MASTER, MPI_COMM_WORLD,MPIERR)
      IF (MYID /= MASTER) ALLOCATE(frq%FREQ(frq%NFREQ)) 
      CALL MPI_BCAST(NSLICES,1,MPI_INTEGER,MASTER, MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(ISRC   ,1,MPI_INTEGER,MASTER, MPI_COMM_WORLD,MPIERR) 
      CALL MPI_BCAST(DT  ,1,MPI_REAL,MASTER, MPI_COMM_WORLD,MPIERR)
      CALL MPI_BCAST(TMAX,1,MPI_REAL,MASTER, MPI_COMM_WORLD,MPIERR) 
      CALL MPI_BCAST(frq%FREQ,frq%NFREQ,MPI_DOUBLE_PRECISION,MASTER, MPI_COMM_WORLD,MPIERR)
      CALL VTKMOV(MYID,MASTER,MPI_COMM_WORLD,NPROCS, PROJNM,NSLICES,DT,ISRC, &
                  frq%NFREQ,frq%FREQ, IERR) 
      IF (IERR /= 0) THEN
         WRITE(*,*) 'xmovie25: An error occurred on process',MYID
         CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)
      ENDIF
      DEALLOCATE(frq%FREQ) 
      CALL MPI_FINALIZE(MPIERR)
      STOP
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE VTKMOV(MYID,MASTER,MYCOMM,NPROCS, PROJNM,NSLICE,DT,ISRC,NFREQ,FREQS, &
                        IERR) 
! 
!     Creates AVS/UCD movies 
      INCLUDE 'mpif.h'
      CHARACTER(80) PROJNM 
      REAL*8 FREQS(NFREQ) 
      CHARACTER(80) INFILE
      CHARACTER(12) CFREQ
      CHARACTER(5) CSRC
      ALLOCATABLE WAVE(:,:,:)
      COMPLEX*8 WAVE, CTWO, CFACT 
      ALLOCATABLE BUFF(:,:), BUFFW(:), UTSLICE(:), VTSLICE(:), WTSLICE(:),  &
                  XLOCS4(:), ZLOCS4(:)
      REAL*4 BUFF, UTSLICE, VTSLICE, WTSLICE, XLOCS4, ZLOCS4
      INTEGER*4 IENG(:,:)
      ALLOCATABLE IENG 
      LOGICAL*4 LEX,LSWAP 
      PARAMETER(CTWO = CMPLX(2.0,0.0))  
      PARAMETER(EPS = 1.E-30) 
      REAL*4 PI/3.1415926535897932384626434E0/
!
!----------------------------------------------------------------------------------------!
! 
!.... figure out file names and sizes 
      IERR = 0 
      CALL DBLANK(80,PROJNM,LENP)
      IF (MYID == MASTER) THEN
         INQUIRE(FILE='./movie_vtk',EXIST=LEX)
         IF (.NOT.LEX) CALL SYSTEM('mkdir movie_vtk')
         CALL NULLS(80,INFILE) 
         CALL NULLS(12,CFREQ) 
         CALL NULLS(5,CSRC) 
         WRITE(CSRC,'(I5)') ISRC
         WRITE(CFREQ,'(F12.5)') FREQS(1) 
         INFILE ='./wave/'//PROJNM(1:LENP)//'_wave-'//CFREQ//'-'//CSRC//'.ewav'
         CALL DBLANK(80,INFILE,LENIN) 
         CALL EWAVHD25D(INFILE, LSWAP,NNPG,NDIM,NELEM,NGNOD,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'vtkmov: Cannot parse header!'
            RETURN
         ENDIF 
      ENDIF 
      CALL MPI_BCAST(PROJNM,80,MPI_CHARACTER,MASTER, MYCOMM,MPIERR) 
      CALL MPI_BCAST(LSWAP,1,MPI_LOGICAL,MASTER, MYCOMM,MPIERR)
      CALL MPI_BCAST(NDIM ,1,MPI_INTEGER,MASTER, MYCOMM,MPIERR)
      CALL MPI_BCAST(NNPG ,1,MPI_INTEGER,MASTER, MYCOMM,MPIERR)
      CALL MPI_BCAST(NELEM,1,MPI_INTEGER,MASTER, MYCOMM,MPIERR)
      CALL MPI_BCAST(NGNOD,1,MPI_INTEGER,MASTER, MYCOMM,MPIERR) 
      CALL MPI_BCAST(LENP ,1,MPI_INTEGER,MASTER, MYCOMM,MPIERR) 
      MNPG = NNPG
      MGNOD = NGNOD
      ALLOCATE(XLOCS4(NNPG))
      ALLOCATE(ZLOCS4(NNPG))
      ALLOCATE(IENG(MGNOD,NELEM))
! 
!.... set up space
      MYFREQ = 0 
      DO 1 IFREQL=1,NFREQ,NPROCS
         IFREQ = IFREQL + MYID
         IF (IFREQ <= NFREQ) MYFREQ = MYFREQ + 1
    1 CONTINUE 
      MYFREQS = MYFREQ
      ALLOCATE(WAVE(MNPG,MYFREQS,NDIM)) 
      MYFREQ = 0
      DO 2 IFREQL=1,NFREQ,NPROCS
         IFREQ = IFREQL + MYID 
         CALL NULLS(80,INFILE)
         CALL NULLS(12,CFREQ)
         CALL NULLS(5,CSRC)
         WRITE(CSRC,'(I5)') ISRC
         IF (IFREQ <= NFREQ) WRITE(CFREQ,'(F12.5)') FREQS(IFREQ)
         INFILE ='./wave/'//PROJNM(1:LENP)//'_wave-'//CFREQ//'-'//CSRC//'.ewav'
         CALL DBLANK(80,INFILE,LENFL) 
         DO 3 IPROCS=1,NPROCS
            IF (MYID == IPROCS - 1) THEN
               IF (IFREQ <= NFREQ) THEN
                  WRITE(*,900) MYID,ISRC,FREQS(IFREQ)
  900             FORMAT(' ****************************************************** ',/, &
                         ' *   Process',I6,'                                    *',/, &
                         ' *   Reading source:',I3,' and frequency:',F12.5,'    * ',/,  &   
                         ' ****************************************************** ',/) 
               ENDIF
            ENDIF
            CALL MPI_BARRIER(MYCOMM,MPIERR)
    3    CONTINUE 
! 
!....... read file
         IF (IFREQ <= NFREQ) THEN
            MYFREQ = MYFREQ + 1  
            CALL EWAVIN2D(MNPG,MGNOD, INFILE,LSWAP,NNPG,NDIM,NELEM,NGNOD, &
                          IENG,XLOCS4,ZLOCS4,WAVE(1:MNPG,MYFREQ,1:NDIM))
         ENDIF 
    2 CONTINUE 
      ALLOCATE(UTSLICE(NNPG))
      ALLOCATE(VTSLICE(NNPG))
      ALLOCATE(WTSLICE(NNPG)) 
      IF (MYID == MASTER) THEN
         ALLOCATE(BUFF(MNPG,3))
         ALLOCATE(BUFFW(NNPG)) 
      ELSE
         ALLOCATE(BUFFW(1))
      ENDIF 
! 
!.... loop on slices
      DO 10 K=1,NSLICE 
! 
!....... fourier transform for this slice 
         TIME = FLOAT(K - 1)*DT
         UTSLICE(1:NNPG) = 0.0 
         VTSLICE(1:NNPG) = 0.0
         WTSLICE(1:NNPG) = 0.0
         MYFREQ = 0
         DO 11 IFREQL=1,NFREQ,NPROCS 
            IFREQ = IFREQL + MYID
            IF (IFREQ <= NFREQ) THEN 
               MYFREQ = MYFREQ + 1
               OMEGA = 2.0*PI*SNGL(FREQS(IFREQ)) 
               ARG = OMEGA*TIME 
               CFACT = CEXP(CMPLX(0.0,+ARG)) !e^{+iwt} 
               DO 12 INPG=1,NNPG
                  IF (ABS(OMEGA) <= EPS) THEN !0 frequency 
                     UTSLICE(INPG) = UTSLICE(INPG) + REAL(WAVE(INPG,MYFREQ,1))
                     VTSLICE(INPG) = VTSLICE(INPG) + REAL(WAVE(INPG,MYFREQ,2))
                     WTSLICE(INPG) = WTSLICE(INPG) + REAL(WAVE(INPG,MYFREQ,3))
                  ELSEIF (ABS(OMEGA - PI/DT) <= EPS) THEN !Nyquist 
                     UTSLICE(INPG) = UTSLICE(INPG) + REAL(WAVE(INPG,MYFREQ,1)*CFACT)
                     VTSLICE(INPG) = VTSLICE(INPG) + REAL(WAVE(INPG,MYFREQ,2)*CFACT)
                     WTSLICE(INPG) = WTSLICE(INPG) + REAL(WAVE(INPG,MYFREQ,3)*CFACT)
                  ELSE !typical invsere-dft 
                     UTSLICE(INPG) = UTSLICE(INPG) + 2.0*REAL(WAVE(INPG,MYFREQ,1)*CFACT)
                     VTSLICE(INPG) = VTSLICE(INPG) + 2.0*REAL(WAVE(INPG,MYFREQ,2)*CFACT)
                     WTSLICE(INPG) = WTSLICE(INPG) + 2.0*REAL(WAVE(INPG,MYFREQ,3)*CFACT)
                  ENDIF
   12          CONTINUE 
            ENDIF
   11    CONTINUE !loop on frequencies 
! 
!....... reduce onto head node 
         CALL MPI_REDUCE(UTSLICE,BUFFW,NNPG,MPI_REAL, MPI_SUM,MASTER, MYCOMM,MPIERR)
         IF (MYID == MASTER) BUFF(1:NNPG,1) = BUFFW(1:NNPG)
         CALL MPI_REDUCE(VTSLICE,BUFFW,NNPG,MPI_REAL, MPI_SUM,MASTER, MYCOMM,MPIERR)
         IF (MYID == MASTER) BUFF(1:NNPG,2) = BUFFW(1:NNPG)
         CALL MPI_REDUCE(WTSLICE,BUFFW,NNPG,MPI_REAL, MPI_SUM,MASTER, MYCOMM,MPIERR)
         IF (MYID == MASTER) BUFF(1:NNPG,3) = BUFFW(1:NNPG) 
          
! 
!....... write file
         IF (MYID == MASTER) THEN 
            UTSLICE(1:NNPG) = UTSLICE(1:NNPG)/REAL(NSLICE)
            VTSLICE(1:NNPG) = VTSLICE(1:NNPG)/REAL(NSLICE)
            WTSLICE(1:NNPG) = WTSLICE(1:NNPG)/REAL(NSLICE)
            CALL PLOT_ELRESP_VTK(MGNOD,MNPG,NNPG,PROJNM, NELEM,K,   &
                                 ISRC,3,IENG, XLOCS4,ZLOCS4,BUFF) 
         ENDIF
   10 CONTINUE !loop on slices 
      DEALLOCATE(WAVE) 
      DEALLOCATE(UTSLICE)
      DEALLOCATE(WTSLICE)
      IF (MYID == MASTER) THEN
         DEALLOCATE(BUFF)
      ENDIF
      RETURN
      END

