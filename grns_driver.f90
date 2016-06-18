      SUBROUTINE GEN_GRNS(MYID,MASTER,MGCOMM, TMPDIR, LCGRNS, LNSRF, &
                          SRC,M1D,FRQ,MSH, IERR) 
!
!     Generates the Greens functions for frequencies/sources to model.  If we were 
!     just doing body waves we could do these calculations on the fly however, 
!     surface waves cause a limitation in that their apparent slowness in y is 
!     frequency dependent.  The outer loop on frequency and inner loop on sources 
!     has been combined and parallelized.  The algorithm is blocking so it isn't the 
!     most efficient but I really don't care.  You are going to chew computer time 
!     with this code anyway, wait the extra minute.  
!     
!     A word on the source time function.  Here it is hardwired to unity so that 
!     we later can convolve our source time function.  For a forward simulation this 
!     of course is a somewhat silly strategy but in inversion our source time function 
!     may be updated at each it iteration and we may want to recycle our old Green's 
!     functions so convolution after the fact becomes an attractive approach.  
!
!     As discussed later I'm not thrilled with the surface wave calculation.  I wouldn't
!     use this code beyond the fundamental period.  
!     -  B Baker January 2013 
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       frequencies to model
!     LCGRNS     True -> calculate greens functions, 
!                False -> recycle greens functions and read read pytab
!     MASTER     master process ID
!     MGCOMM     global communicator
!     MYID       process ID  
!     NFREQ      number of frequencies
!     TMPDIR     temprorary directory for files 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     PYTAB      py at frequency and source (through common)  
!     VFAST      max velocity in 1D models (through common)
!
!.... variable declarations 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO)  MSH
      TYPE (SRC_INFO)   SRC 
      TYPE (MOD1D_INFO) M1D 
      TYPE (FRQ_INFO)   FRQ 
      CHARACTER(*), INTENT(IN) :: TMPDIR
      INTEGER*4, INTENT(IN) :: MYID,MASTER,MGCOMM
      LOGICAL*4 LCGRNS, LNSRF
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables 
      !COMPLEX*8, ALLOCATABLE :: UE(:) 
      COMPLEX*8 UE(1) 
      REAL*8, ALLOCATABLE :: CCRAY(:), CCLOV(:), PYTABB(:), &
                             BUFF(:)
      COMPLEX*16 STF
      REAL*8 PERIOD, XMLAT0,XMLON0, XMLAT1,XMLON1, VBASE, PI180, CVBASE, VFAST  
      INTEGER*4 LOOP, MYLOC, IPROC, ISRC, IFREQ, IPY, NPROCS, MAXMOD, NWORK, &
                INDX, MPIERR 
      LOGICAL*4 LEX, LSHIFT, LISDIR
      PARAMETER(STF = DCMPLX(1.D0,0.0)) 
      PARAMETER(PI180 = 0.017453292519943295D0)
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      IF (MYID == MASTER) WRITE(*,*) 
      CALL MPI_BCAST(LCGRNS,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      IF (.NOT.LCGRNS) THEN
         IF (MYID == MASTER) WRITE(*,*) 'gen_grns: Reading old py table'
         GOTO 800
      ENDIF
      CALL MPI_COMM_SIZE(MGCOMM,NPROCS,MPIERR)
!
!.... could just be reading 
      IF (MYID == MASTER) THEN
         LEX = LISDIR('./grns')
         IF (LEX) CALL SYSTEM('rm -rf ./grns')
         CALL SYSTEM('mkdir ./grns')
      ENDIF
      CALL MPI_BARRIER(MGCOMM,MPIERR) 

      XMLAT0 = msh%XLATMIN 
      XMLON0 = msh%XLONMIN
      XMLAT1 = msh%XLATMAX
      XMLON1 = msh%XLONMAX
      LSHIFT = .FALSE. !shift phase velocity velocity 
      MAXMOD = MAXVAL(src%MODE(1:src%NSRC),1) + 1 
      ALLOCATE(CCRAY(20)) 
      ALLOCATE(CCLOV(20))
      ALLOCATE(PYTABB(frq%NFREQ*src%NSRC))  
      IF (MAXMOD > 1 .AND. MYID == MASTER) THEN
         WRITE(*,*) 'gen_grns: Warning max mode > 0.  I dont know what will happen'
         WRITE(*,*) 'gen_grns: You should have little confidence it will work'
      ENDIF
      CCRAY(1:20) = 0.D0
      CCLOV(1:20) = 0.D0
      PERIOD = 0.D0
      PYTABB(1:frq%NFREQ*src%NSRC) = 0.D0
      VFAST = 0.D0
!
!.... loop on frequencies and sources
      DO 100 LOOP=1,frq%NFREQ*src%NSRC,NPROCS
         MYLOC = LOOP + MYID
         ISRC = MOD(MYLOC,src%NSRC) 
         IF (ISRC == 0) ISRC = src%NSRC !multiple is evenly divisible
         IFREQ = (MYLOC - ISRC)/src%NSRC + 1
         IPY = (IFREQ - 1)*src%NSRC + ISRC    !py index 
         DO 101 IPROC=0,NPROCS-1
            IF (MYID == IPROC .AND. IFREQ <= frq%NFREQ .AND. ISRC <= src%NSRC) THEN 
               WRITE(*,9500) MYID,ISRC,frq%FREQ(IFREQ)
 9500          FORMAT(' --------------------------------------------------------',/, &
                      ' -   gen_grns: Process:',I5,'                            -',/, &
                      ' -             Calculating 1D solution for source:',I5,' -',/, &
                      ' -             At frequency:',F12.5,' Hz             -',/,  &
                      ' --------------------------------------------------------',/)
            ENDIF
            CALL MPI_BARRIER(MGCOMM,MPIERR)
  101    CONTINUE 
         INDX = (ISRC - 1)*frq%NFREQ + IFREQ 
         IF (ISRC > src%NSRC .OR. IFREQ > frq%NFREQ) GOTO 1001 !nothing for you to do 
         IF (src%SRCTYP(ISRC)(1:1) == 'P') THEN !body waves
            VBASE = CVBASE(src%SRCTYP(ISRC),src%CSIDE(ISRC), m1d%NL1D_LT,m1d%NL1D_RT,   &
                           m1d%VP1D_LT,m1d%VS1D_LT, m1d%VP1D_RT,m1d%VS1D_RT)
            PYTABB(INDX) = DSIN(src%BAZN(ISRC)*PI180)*DSIN(src%AOI(ISRC)*PI180)/VBASE 
            IF (DABS(PYTABB(INDX)) < 1.D-10) PYTABB(INDX) = 0.D0 
            CALL GRNS_BDY(src%SRCTYP(ISRC), NDIM,msh%NDOF,msh%NNPE, .TRUE.,           &
                           m1d%NL1D_LT,m1d%NL1D_RT, src%CSIDE(ISRC), NDIM, ISRC,       &
                           frq%FREQ(IFREQ),msh%FREQ0,src%AOI(ISRC),src%BAZN(ISRC),     &
                           msh%XBLKL,msh%XBLKR, msh%ZBASE_INT,  &
                           msh%XMOD0,msh%XMOD1, STF, msh%IDOFSE,                       & 
                           m1d%VP1D_LT,m1d%VS1D_LT,m1d%RH1D_LT,                        &
                           m1d%QP1D_LT,m1d%QS1D_LT, m1d%Z1D_LT,             &
                           m1d%VP1D_RT,m1d%VS1D_RT,m1d%RH1D_RT,             & 
                           m1d%QP1D_RT,m1d%QS1D_RT, m1d%Z1D_RT,              &
                           msh%XLOCSE,msh%ZLOCSE, UE,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'gen_grns: An error occurred in grns_bdy on process:',MYID
               GOTO 1000
            ENDIF
            VFAST = DMAX1(VFAST,VBASE)
         ELSEIF (src%SRCTYP(ISRC)(1:1) == 'S') THEN !surface waves
            !print *, myid,srctyp(isrc), ndim,ndof,nnpe, nl1d_lt,nl1d_rt, cside(isrc)
            !print *, ifreq,isrc
            !print *, mode(ifreq),freq(ifreq), slat(isrc),slon(isrc),sdep(isrc),smag(isrc)
            CALL GRNS_SRF(TMPDIR, LNSRF, &
                          MYID,src%SRCTYP(ISRC), NDIM,msh%NDOF,msh%NNPE, .TRUE.,LSHIFT,  &
                          m1d%NL1D_LT,m1d%NL1D_RT, src%CSIDE(ISRC), NDIM,                &
                          src%MODE(ISRC), ISRC,msh%FREQ0,frq%FREQ(IFREQ), src%BAZN(ISRC),&
                          src%SLAT(ISRC),src%SLON(ISRC),src%SDEP(ISRC),src%SMAG(ISRC),   &
                          src%STRIKE(ISRC),src%DIP(ISRC),src%RAKE(ISRC),                 &
                          msh%XBLKL,msh%XBLKR, msh%XMOD0,msh%XMOD1,                      &
                          XMLAT0,XMLON0,XMLAT1,XMLON1,  &
                          STF,msh%IDOFSE,  &
                          m1d%VPD_RLLT,m1d%VSD_RLLT,m1d%ROD_RLLT,m1d%HDD_RLLT,           &
                          m1d%VPD_LVLT,m1d%VSD_LVLT,m1d%ROD_LVLT,m1d%HDD_LVLT,           &
                          m1d%VPD_RLRT,m1d%VSD_RLRT,m1d%ROD_RLRT,m1d%HDD_RLRT,           &
                          m1d%VPD_LVRT,m1d%VSD_LVRT,m1d%ROD_LVRT,m1d%HDD_LVRT,           &
                          m1d%QP1D_LT,m1d%QP1D_RT, m1d%QS1D_LT,m1d%QS1D_RT,              & 
                          m1d%Z1D_LT, m1d%Z1D_RT, msh%XLOCSE,msh%ZLOCSE,                 &
                          PERIOD,CCRAY,CCLOV, UE,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'gen_grns: An error occurred in grns_srf on process:',MYID
               GOTO 1000
            ENDIF
!
!.......... average the phase velocity for higher modes?
!           NAVG = 0
!           IF (SRCTYP(ISRC)(2:2) == 'L') THEN
!              DO IMODE=0,MODE(ISRC)-1
!                 IF (CCLOV(IMODE) > 0.D0) 
!                    NAVG = NAVG + 1
!                    PYTABB(INDX) = CCLOV(IMODE+1)
!                 ENDIF
!              ENDDO
!           ELSE
!              DO IMODE=0,MODE(ISRC)-1 
!                 IF (CCRAY(IMODE) > 0.D0) 
!                    NAVG = NAVG + 1
!                    PYTABB(INDX) = CCRAY(IMODE+1)
!                 ENDIF   
!              ENDDO   
!           ENDIF
!           IF (NAVG > 0) &
!           PYTABB(INDX) = DSIN(BAZN(ISRC)*PI180)/PYTABB(INDX)/DFLOAT(NAVG)
            IF (src%SRCTYP(ISRC)(2:2) == 'L') THEN
               PYTABB(INDX) = DSIN(src%BAZN(ISRC)*PI180)/(CCLOV(1)*1.D3)
               VFAST = DMAX1(VFAST,CCLOV(1)*1.D3)
            ELSE
               PYTABB(INDX) = DSIN(src%BAZN(ISRC)*PI180)/(CCRAY(1)*1.D3)
               VFAST = DMAX1(VFAST,CCRAY(1)*1.D3) 
            ENDIF
            IF (DABS(PYTABB(INDX)) < 1.D-10) PYTABB(INDX) = 0.D0
            !LNEW = .FALSE.
         ELSE !i don't know what to do 
            WRITE(*,*) 'gen_grns: Invalid source type!',MYID
            IERR = 1 
            GOTO 1000
         ENDIF
 1001    CONTINUE 
  100 CONTINUE 
 1000 CONTINUE !error break
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_grns: An error occurred on process:',MYID
         RETURN
      ENDIF 
!
!.... have head node figure out py table
      NWORK = frq%NFREQ*src%NSRC
      ALLOCATE(BUFF(NWORK)) 
      CALL MPI_ALLREDUCE(PYTABB,BUFF,NWORK, MPI_DOUBLE_PRECISION,MPI_SUM, MGCOMM,MPIERR)
      CALL MPI_ALLREDUCE(VFAST,m1d%VFAST,1, MPI_DOUBLE_PRECISION,MPI_MAX, MGCOMM,MPIERR)
      m1d%VFAST_BDY = m1d%VFAST
      m1d%VFAST_SRF = m1d%VFAST
      INDX = 0
      DO 10 ISRC=1,src%NSRC
         DO 11 IFREQ=1,frq%NFREQ
            INDX = INDX + 1
            src%PYTAB(IFREQ,ISRC) = BUFF(INDX) 
   11    CONTINUE
   10 CONTINUE 
      DEALLOCATE(BUFF) 
      DEALLOCATE(PYTABB) 
!
!.... also have head node dump the py table
      IF (MYID == MASTER) &
      CALL WRITE_PYTAB(frq%NFREQ,frq%NFREQ,src%NSRC, m1d%VFAST,m1d%VFAST_BDY,m1d%VFAST_SRF, &
                       frq%FREQ,src%PYTAB, IERR)
!
!.... done with the 1D models
      IF (MYID == MASTER) THEN
         WRITE(*,*) 'gen_grns: Deallocating 1D models...'
         WRITE(*,*) 
      ENDIF
!
!.... break ahead for read only
  800 CONTINUE 
      IF (MYID == MASTER .AND. .NOT.LCGRNS) THEN 
         CALL READ_PYTAB(frq%NFREQ,frq%NFREQ,src%NSRC, frq%FREQ,  &
                         m1d%VFAST,m1d%VFAST_BDY,m1d%VFAST_SRF,src%PYTAB,IERR) 
         IF (IERR /= 0) WRITE(*,*) 'gen_grns: Error reading py table'
      ENDIF
      IF (.NOT.LCGRNS) THEN
         CALL MPI_BCAST(m1d%VFAST,1,MPI_DOUBLE_PRECISION, MASTER,MGCOMM,MPIERR)
         DO 12 ISRC=1,src%NSRC
            CALL MPI_BCAST(src%PYTAB(1:frq%NFREQ,ISRC),frq%NFREQ,MPI_DOUBLE_PRECISION, &
                           MASTER,MGCOMM,MPIERR)
   12    CONTINUE
      ENDIF

      IF (ASSOCIATED(m1d%VP1D_LT))  DEALLOCATE(m1d%VP1D_LT)
      IF (ASSOCIATED(m1d%VS1D_LT))  DEALLOCATE(m1d%VS1D_LT)
      IF (ASSOCIATED(m1d%RH1D_LT))  DEALLOCATE(m1d%RH1D_LT)
      IF (ASSOCIATED(m1d%QP1D_LT))  DEALLOCATE(m1d%QP1D_LT)
      IF (ASSOCIATED(m1d%QS1D_LT))  DEALLOCATE(m1d%QS1D_LT)
      IF (ASSOCIATED(m1d% Z1D_LT))  DEALLOCATE(m1d% Z1D_LT)
      IF (ASSOCIATED(m1d%VPD_RLLT)) DEALLOCATE(m1d%VPD_RLLT)
      IF (ASSOCIATED(m1d%VSD_RLLT)) DEALLOCATE(m1d%VSD_RLLT)
      IF (ASSOCIATED(m1d%ROD_RLLT)) DEALLOCATE(m1d%ROD_RLLT)
      IF (ASSOCIATED(m1d%HDD_RLLT)) DEALLOCATE(m1d%HDD_RLLT)
      IF (ASSOCIATED(m1d%VPD_LVLT)) DEALLOCATE(m1d%VPD_LVLT)
      IF (ASSOCIATED(m1d%VSD_LVLT)) DEALLOCATE(m1d%VSD_LVLT)
      IF (ASSOCIATED(m1d%ROD_LVLT)) DEALLOCATE(m1d%ROD_LVLT)
      IF (ASSOCIATED(m1d%HDD_LVLT)) DEALLOCATE(m1d%HDD_LVLT)

      IF (ASSOCIATED(m1d%VP1D_RT))  DEALLOCATE(m1d%VP1D_RT)
      IF (ASSOCIATED(m1d%VS1D_RT))  DEALLOCATE(m1d%VS1D_RT)
      IF (ASSOCIATED(m1d%RH1D_RT))  DEALLOCATE(m1d%RH1D_RT)
      IF (ASSOCIATED(m1d%QP1D_RT))  DEALLOCATE(m1d%QP1D_RT)
      IF (ASSOCIATED(m1d%QS1D_RT))  DEALLOCATE(m1d%QS1D_RT)
      IF (ASSOCIATED(m1d% Z1D_RT))  DEALLOCATE(m1d% Z1D_RT)
      IF (ASSOCIATED(m1d%VPD_RLRT)) DEALLOCATE(m1d%VPD_RLRT)
      IF (ASSOCIATED(m1d%VSD_RLRT)) DEALLOCATE(m1d%VSD_RLRT)
      IF (ASSOCIATED(m1d%ROD_RLRT)) DEALLOCATE(m1d%ROD_RLRT)
      IF (ASSOCIATED(m1d%HDD_RLRT)) DEALLOCATE(m1d%HDD_RLRT)
      IF (ASSOCIATED(m1d%VPD_LVRT)) DEALLOCATE(m1d%VPD_LVRT)
      IF (ASSOCIATED(m1d%VSD_LVRT)) DEALLOCATE(m1d%VSD_LVRT)
      IF (ASSOCIATED(m1d%ROD_LVRT)) DEALLOCATE(m1d%ROD_LVRT)
      IF (ASSOCIATED(m1d%HDD_LVRT)) DEALLOCATE(m1d%HDD_LVRT)
 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GEN_GRNS_SB(MYID,MASTER,MGCOMM, TMPDIR, LCGRNS,LNSRF, &
                             SRC,M1D,FRQ,MSH, IERR) 
!
!     Generates the Greens functions for frequencies/sources to model.  If we were 
!     just doing body waves we could do these calculations on the fly however, 
!     surface waves cause a limitation in that their apparent slowness in y is 
!     frequency dependent.  The outer loop on frequency and inner loop on sources 
!     has been combined and parallelized.  The algorithm is blocking so it isn't the 
!     most efficient but I really don't care.  You are going to chew computer time 
!     with this code anyway, wait the extra minute.  
!     
!     A word on the source time function.  Here it is hardwired to unity so that 
!     we later can convolve our source time function.  For a forward simulation this 
!     of course is a somewhat silly strategy but in inversion our source time function 
!     may be updated at each it iteration and we may want to recycle our old Green's 
!     functions so convolution after the fact becomes an attractive approach.  
!
!     As discussed later I'm not thrilled with the surface wave calculation.  I wouldn't
!     use this code beyond the fundamental period.  
!
!     The modification is this handles mixed frequency source lists
!     -  B Baker May 2013 
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       frequencies to model
!     LCGRNS     True -> calculate greens functions, 
!                False -> recycle greens functions and read read pytab
!     MASTER     master process ID
!     MGCOMM     global communicator
!     MYID       process ID  
!     NFREQ      number of frequencies
!     TMPDIR     temporary directory for files
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     PYTAB      py at frequency and source (through common)  
!     VFAST      max velocity in 1D models (through common)
!
!.... variable declarations 
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO)  MSH  
      TYPE (SRC_INFO)   SRC  
      TYPE (MOD1D_INFO) M1D  
      TYPE (FRQ_INFO)   FRQ  
      CHARACTER(*), INTENT(IN) :: TMPDIR 
      INTEGER*4, INTENT(IN) :: MYID,MASTER,MGCOMM
      LOGICAL*4 LCGRNS, LNSRF
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables 
      COMPLEX*8 UE(1)
      REAL*8, ALLOCATABLE :: CCRAY(:), CCLOV(:), PYTABB(:), &
                             BUFF(:)
      INTEGER*4, ALLOCATABLE :: IFREQ_TAB(:), ISRC_TAB(:) 
      COMPLEX*16 STF
      REAL*8 PERIOD, XMLAT0,XMLON0, XMLAT1,XMLON1, VBASE, PI180, CVBASE, &
             VFAST_SRF,  VFAST_BDY 
      INTEGER*4 LOOP, MYLOC, IPROC, ISRC, IFREQ, IPY, NPROCS, MAXMOD, NWORK, &
                JFREQ, JSRC, INDX, MPIERR
      LOGICAL*4 LEX, LISDIR, LSHIFT
      PARAMETER(STF = DCMPLX(1.D0,0.0))
      PARAMETER(PI180 = 0.017453292519943295D0)
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      IF (MYID == MASTER) WRITE(*,*)
      CALL MPI_BCAST(LCGRNS,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      IF (.NOT.LCGRNS) THEN
         IF (MYID == MASTER) WRITE(*,*) 'gen_grns_sb: Reading old py table'
         GOTO 800
      ENDIF
      CALL MPI_COMM_SIZE(MGCOMM,NPROCS,MPIERR)
!
!.... could just be reading 
      IF (MYID == MASTER) THEN
         LEX = LISDIR('./grns')
         IF (LEX) CALL SYSTEM('rm -rf ./grns')
         CALL SYSTEM('mkdir ./grns')
      ENDIF
      CALL MPI_BARRIER(MGCOMM,MPIERR)

      XMLAT0 = msh%XLATMIN
      XMLON0 = msh%XLONMIN
      XMLAT1 = msh%XLATMAX
      XMLON1 = msh%XLONMAX
      LSHIFT = .FALSE. !shift phase velocity 
      MAXMOD = MAXVAL(src%MODE(1:src%NSRC),1) + 1
      ALLOCATE(CCRAY(20))
      ALLOCATE(CCLOV(20))
      ALLOCATE(PYTABB(frq%NFREQ*src%NSRC))
      IF (MAXMOD > 1 .AND. MYID == MASTER) THEN
         WRITE(*,*) 'gen_grns_sb: Warning max mode > 0.  I dont know what will happen'
         WRITE(*,*) 'gen_grns_sb: You should have little confidence it will work'
      ENDIF
      CCRAY(1:20) = 0.D0
      CCLOV(1:20) = 0.D0
      PERIOD = 0.D0
      PYTABB(1:frq%NFREQ*src%NSRC) = 0.D0
      VFAST_BDY = 0.D0
      VFAST_SRF = 0.D0
      NWORK = frq%NFREQ_SRF*src%NSRC_SRF + frq%NFREQ_BDY*src%NSRC_BDY
      ALLOCATE(ISRC_TAB(NWORK)) 
      ALLOCATE(IFREQ_TAB(NWORK)) 
      IF (MYID == MASTER) THEN
         CALL GEN_FRSRC_TAB(NWORK,frq%NFREQ,src%NSRC, src%SRCTYP,frq%CFTYPE, &
                            IFREQ_TAB,ISRC_TAB)
      ENDIF
      CALL MPI_BCAST(IFREQ_TAB,NWORK,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(ISRC_TAB ,NWORK,MPI_INTEGER, MASTER,MGCOMM,MPIERR)
!
!.... loop on frequencies and sources
      DO 100 LOOP=1,NWORK,NPROCS !frq%NFREQ*src%NSRC,NPROCS
         MYLOC = LOOP + MYID
!        ISRC = MOD(MYLOC,src%NSRC)
!        IF (ISRC == 0) ISRC = src%NSRC !multiple is evenly divisible
!        IFREQ = (MYLOC - ISRC)/src%NSRC + 1
         IFREQ = frq%NFREQ + 1
         ISRC  = src%NSRC + 1 
         IF (MYLOC <= NWORK) THEN
            IFREQ = IFREQ_TAB(MYLOC)
            ISRC  = ISRC_TAB(MYLOC) 
         ENDIF
         IPY = (IFREQ - 1)*src%NSRC + ISRC    !py index 
         DO 101 IPROC=0,NPROCS-1
!           IF (MYID == IPROC .AND. IFREQ <= frq%NFREQ .AND. ISRC <= src%NSRC) THEN
            IF (MYID == MASTER) THEN
               JFREQ = frq%NFREQ + 1 
               JSRC = src%NSRC + 1
               IF (MYLOC + IPROC <= NWORK) THEN
                  JFREQ = IFREQ_TAB(MYLOC+IPROC)
                  JSRC  = ISRC_TAB(MYLOC+IPROC) 
               ENDIF
               IF (JFREQ <= frq%NFREQ .AND. JSRC <= src%NSRC) &
               WRITE(*,9500) MYID+IPROC,JSRC,frq%FREQ(JFREQ)
 9500          FORMAT(' -----------------------------------------------------------',/, &
                      ' -   gen_grns_sb: Process:',I5,'                            -',/, &
                      ' -                Calculating 1D solution for source:',I5,' -',/, &
                      ' -                At frequency:',F12.5,' Hz             -',/,  &
                      ' -----------------------------------------------------------',/)
            ENDIF
            CALL MPI_BARRIER(MGCOMM,MPIERR)
  101    CONTINUE
         INDX = (ISRC - 1)*frq%NFREQ + IFREQ
         IF (ISRC > src%NSRC .OR. IFREQ > frq%NFREQ) GOTO 1001 !nothing for you to do 
!
!....... body waves
         IF (src%SRCTYP(ISRC)(1:1) == 'P' .AND. frq%CFTYPE(IFREQ) == 'B') THEN !body waves
            VBASE = CVBASE(src%SRCTYP(ISRC),src%CSIDE(ISRC), m1d%NL1D_LT,m1d%NL1D_RT,   &
                           m1d%VP1D_LT,m1d%VS1D_LT, m1d%VP1D_RT,m1d%VS1D_RT)
            PYTABB(INDX) = DSIN(src%BAZN(ISRC)*PI180)*DSIN(src%AOI(ISRC)*PI180)/VBASE
            IF (DABS(PYTABB(INDX)) < 1.D-10) PYTABB(INDX) = 0.D0
            CALL GRNS_BDY(src%SRCTYP(ISRC), NDIM,msh%NDOF,msh%NNPE, .TRUE.,           &
                           m1d%NL1D_LT,m1d%NL1D_RT, src%CSIDE(ISRC), NDIM, ISRC,       &
                           frq%FREQ(IFREQ),msh%FREQ0,src%AOI(ISRC),src%BAZN(ISRC),     &
                           msh%XBLKL,msh%XBLKR,msh%ZBASE_INT,  &
                           msh%XMOD0,msh%XMOD1, STF, msh%IDOFSE,                       &
                           m1d%VP1D_LT,m1d%VS1D_LT,m1d%RH1D_LT,                        &
                           m1d%QP1D_LT,m1d%QS1D_LT, m1d%Z1D_LT,             &
                           m1d%VP1D_RT,m1d%VS1D_RT,m1d%RH1D_RT,             &
                           m1d%QP1D_RT,m1d%QS1D_RT, m1d%Z1D_RT,              &
                           msh%XLOCSE,msh%ZLOCSE, UE,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'gen_grns_sb: An error occurred in grns_bdy on process:',MYID
               GOTO 1000
            ENDIF
            VFAST_BDY = DMAX1(VFAST_BDY,VBASE)
!
!....... surface waves
         ELSEIF (src%SRCTYP(ISRC)(1:1) == 'S' .AND. frq%CFTYPE(IFREQ) == 'S') THEN 
            CALL GRNS_SRF(TMPDIR,  LNSRF, &
                          MYID,src%SRCTYP(ISRC), NDIM,msh%NDOF,msh%NNPE, .TRUE.,LSHIFT,  &
                          m1d%NL1D_LT,m1d%NL1D_RT, src%CSIDE(ISRC), NDIM,                &
                          src%MODE(ISRC), ISRC,msh%FREQ0,frq%FREQ(IFREQ), src%BAZN(ISRC),&
                          src%SLAT(ISRC),src%SLON(ISRC),src%SDEP(ISRC),src%SMAG(ISRC),   &
                          src%STRIKE(ISRC),src%DIP(ISRC),src%RAKE(ISRC),                 &
                          msh%XBLKL,msh%XBLKR, msh%XMOD0,msh%XMOD1,                      &
                          XMLAT0,XMLON0,XMLAT1,XMLON1,  &
                          STF,msh%IDOFSE,  &
                          m1d%VPD_RLLT,m1d%VSD_RLLT,m1d%ROD_RLLT,m1d%HDD_RLLT,           &
                          m1d%VPD_LVLT,m1d%VSD_LVLT,m1d%ROD_LVLT,m1d%HDD_LVLT,           &
                          m1d%VPD_RLRT,m1d%VSD_RLRT,m1d%ROD_RLRT,m1d%HDD_RLRT,           &
                          m1d%VPD_LVRT,m1d%VSD_LVRT,m1d%ROD_LVRT,m1d%HDD_LVRT,           &
                          m1d%QP1D_LT,m1d%QP1D_RT, m1d%QS1D_LT,m1d%QS1D_RT,              &
                          m1d%Z1D_LT, m1d%Z1D_RT, msh%XLOCSE,msh%ZLOCSE,                 &
                          PERIOD,CCRAY,CCLOV,  UE,IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'gen_grns_sb: An error occurred in grns_srf on process:',MYID
               GOTO 1000
            ENDIF
            IF (src%SRCTYP(ISRC)(2:2) == 'L') THEN
               PYTABB(INDX) = DSIN(src%BAZN(ISRC)*PI180)/(CCLOV(1)*1.D3)
               VFAST_SRF = DMAX1(VFAST_SRF,CCLOV(1)*1.D3)
            ELSE
               PYTABB(INDX) = DSIN(src%BAZN(ISRC)*PI180)/(CCRAY(1)*1.D3)
               VFAST_SRF = DMAX1(VFAST_SRF,CCRAY(1)*1.D3)
            ENDIF
            IF (DABS(PYTABB(INDX)) < 1.D-10) PYTABB(INDX) = 0.D0
            !LNEW = .FALSE.
         ELSE !i may not know what to do 
            IF (src%SRCTYP(ISRC)(1:1) /= 'S' .AND. src%SRCTYP(ISRC)(1:1) /= 'P') THEN
               WRITE(*,*) 'gen_grns_sb: Invalid source type!',MYID
               IERR = 1
               GOTO 1000
            ENDIF
            IF (frq%CFTYPE(IFREQ) /= 'S' .AND. frq%CFTYPE(IFREQ) /= 'B') THEN
               WRITE(*,*) 'gen_grns_sb: Invalid frequency type!',MYID
               IERR = 1
               GOTO 1000
            ENDIF
         ENDIF
 1001    CONTINUE
  100 CONTINUE
 1000 CONTINUE !error break
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_grns: An error occurred on process:',MYID
         RETURN
      ENDIF
!
!.... have head node figure out py table
      NWORK = frq%NFREQ*src%NSRC
      ALLOCATE(BUFF(NWORK))
      CALL MPI_ALLREDUCE(PYTABB,BUFF,NWORK, MPI_DOUBLE_PRECISION,MPI_SUM, MGCOMM,MPIERR)
      CALL MPI_ALLREDUCE(VFAST_BDY,m1d%VFAST_BDY,1, MPI_DOUBLE_PRECISION,MPI_MAX, &
                         MGCOMM,MPIERR)
      CALL MPI_ALLREDUCE(VFAST_SRF,m1d%VFAST_SRF,1, MPI_DOUBLE_PRECISION,MPI_MAX, &
                         MGCOMM,MPIERR) 
      IF (m1d%VFAST_BDY > 0.D0) THEN
         m1d%VFAST = m1d%VFAST_BDY
      ELSE
         m1d%VFAST = m1d%VFAST_SRF
      ENDIF
      INDX = 0
      DO 10 ISRC=1,src%NSRC
         DO 11 IFREQ=1,frq%NFREQ
            INDX = INDX + 1
            src%PYTAB(IFREQ,ISRC) = BUFF(INDX)
   11    CONTINUE
   10 CONTINUE
      DEALLOCATE(BUFF)
      DEALLOCATE(PYTABB)
      DEALLOCATE(ISRC_TAB)
      DEALLOCATE(IFREQ_TAB) 
!
!.... also have head node dump the py table
      IF (MYID == MASTER) &
      CALL WRITE_PYTAB(frq%NFREQ,frq%NFREQ,src%NSRC, m1d%VFAST,m1d%VFAST_BDY,m1d%VFAST_SRF, &
                       frq%FREQ,src%PYTAB, IERR)
!
!.... done with the 1D models
      IF (MYID == MASTER) THEN
         WRITE(*,*) 'gen_grns: Deallocating 1D models...'
         WRITE(*,*)
      ENDIF
!
!.... break ahead for read only
  800 CONTINUE
      IF (MYID == MASTER .AND. .NOT.LCGRNS) THEN
         CALL READ_PYTAB(frq%NFREQ,frq%NFREQ,src%NSRC, frq%FREQ, &
                         m1d%VFAST,m1d%VFAST_BDY,m1d%VFAST_SRF,src%PYTAB,IERR)
         IF (IERR /= 0) WRITE(*,*) 'gen_grns: Error reading py table'
      ENDIF
      IF (.NOT.LCGRNS) THEN
         CALL MPI_BCAST(m1d%VFAST    ,1,MPI_DOUBLE_PRECISION, MASTER,MGCOMM,MPIERR)
         CALL MPI_BCAST(m1d%VFAST_SRF,1,MPI_DOUBLE_PRECISION, MASTER,MGCOMM,MPIERR)
         CALL MPI_BCAST(m1d%VFAST_BDY,1,MPI_DOUBLE_PRECISION, MASTER,MGCOMM,MPIERR)
         DO 12 ISRC=1,src%NSRC
            CALL MPI_BCAST(src%PYTAB(1:frq%NFREQ,ISRC),frq%NFREQ,MPI_DOUBLE_PRECISION, &
                           MASTER,MGCOMM,MPIERR)
   12    CONTINUE
      ENDIF
      IF (ASSOCIATED(m1d%VP1D_LT))  DEALLOCATE(m1d%VP1D_LT)
      IF (ASSOCIATED(m1d%VS1D_LT))  DEALLOCATE(m1d%VS1D_LT)
      IF (ASSOCIATED(m1d%RH1D_LT))  DEALLOCATE(m1d%RH1D_LT)
      IF (ASSOCIATED(m1d%QP1D_LT))  DEALLOCATE(m1d%QP1D_LT)
      IF (ASSOCIATED(m1d%QS1D_LT))  DEALLOCATE(m1d%QS1D_LT)
      IF (ASSOCIATED(m1d% Z1D_LT))  DEALLOCATE(m1d% Z1D_LT)
      IF (ASSOCIATED(m1d%VPD_RLLT)) DEALLOCATE(m1d%VPD_RLLT)
      IF (ASSOCIATED(m1d%VSD_RLLT)) DEALLOCATE(m1d%VSD_RLLT)
      IF (ASSOCIATED(m1d%ROD_RLLT)) DEALLOCATE(m1d%ROD_RLLT)
      IF (ASSOCIATED(m1d%HDD_RLLT)) DEALLOCATE(m1d%HDD_RLLT)
      IF (ASSOCIATED(m1d%VPD_LVLT)) DEALLOCATE(m1d%VPD_LVLT)
      IF (ASSOCIATED(m1d%VSD_LVLT)) DEALLOCATE(m1d%VSD_LVLT)
      IF (ASSOCIATED(m1d%ROD_LVLT)) DEALLOCATE(m1d%ROD_LVLT)
      IF (ASSOCIATED(m1d%HDD_LVLT)) DEALLOCATE(m1d%HDD_LVLT)

      IF (ASSOCIATED(m1d%VP1D_RT))  DEALLOCATE(m1d%VP1D_RT)
      IF (ASSOCIATED(m1d%VS1D_RT))  DEALLOCATE(m1d%VS1D_RT)
      IF (ASSOCIATED(m1d%RH1D_RT))  DEALLOCATE(m1d%RH1D_RT)
      IF (ASSOCIATED(m1d%QP1D_RT))  DEALLOCATE(m1d%QP1D_RT)
      IF (ASSOCIATED(m1d%QS1D_RT))  DEALLOCATE(m1d%QS1D_RT)
      IF (ASSOCIATED(m1d% Z1D_RT))  DEALLOCATE(m1d% Z1D_RT)
      IF (ASSOCIATED(m1d%VPD_RLRT)) DEALLOCATE(m1d%VPD_RLRT)
      IF (ASSOCIATED(m1d%VSD_RLRT)) DEALLOCATE(m1d%VSD_RLRT)
      IF (ASSOCIATED(m1d%ROD_RLRT)) DEALLOCATE(m1d%ROD_RLRT)
      IF (ASSOCIATED(m1d%HDD_RLRT)) DEALLOCATE(m1d%HDD_RLRT)
      IF (ASSOCIATED(m1d%VPD_LVRT)) DEALLOCATE(m1d%VPD_LVRT)
      IF (ASSOCIATED(m1d%VSD_LVRT)) DEALLOCATE(m1d%VSD_LVRT)
      IF (ASSOCIATED(m1d%ROD_LVRT)) DEALLOCATE(m1d%ROD_LVRT)
      IF (ASSOCIATED(m1d%HDD_LVRT)) DEALLOCATE(m1d%HDD_LVRT)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GEN_FRSRC_TAB(NWORK,NFREQ,NSRC, SRCTYP,CFTYPE, IFREQ_TAB,ISRC_TAB) 
!
!     Generates a frequency source pair table 
      IMPLICIT NONE
      CHARACTER(2), INTENT(IN) :: SRCTYP(NSRC) 
      CHARACTER(1), INTENT(IN) :: CFTYPE(NFREQ)
      INTEGER*4, INTENT(IN) :: NWORK, NFREQ, NSRC 
      INTEGER*4, INTENT(OUT) :: IFREQ_TAB(NWORK), ISRC_TAB(NWORK) 
!.... local variables
      INTEGER*4 ITYPE, IFREQ, ISRC, INDX 
!
!----------------------------------------------------------------------------------------!
!
      INDX = 0
      IFREQ_TAB(1:NWORK) = 0
      ISRC_TAB(1:NWORK)  = 0
      DO 1 ITYPE=1,2
         DO 2 IFREQ=1,NFREQ
            IF (ITYPE == 1 .AND. CFTYPE(IFREQ) == 'S' .OR. &
                ITYPE == 2 .AND. CFTYPE(IFREQ) == 'B') THEN
               DO 3 ISRC=1,NSRC
                  IF (ITYPE == 1 .AND. SRCTYP(ISRC)(1:1) == 'S' .OR. &
                      ITYPE == 2 .AND. SRCTYP(ISRC)(1:1) == 'P') THEN
                     INDX = INDX + 1
                     IFREQ_TAB(INDX) = IFREQ
                     ISRC_TAB(INDX)  = ISRC 
                  ENDIF
    3          CONTINUE !loop on sources
            ENDIF
    2    CONTINUE !Loop on frequencies
    1 CONTINUE !loop on types 
      IF (INDX /= NWORK) WRITE(*,*) 'gen_fsrc_tab: Warning indx/= nwork'
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LOAD_GRNS25(MDIM,NDOF,NNPE,NDIM, SRCTYP,ISRC,    &
                             FREQ,STF,IDOFSE,                     &
                             UE,IERR) 
!
!     Loads the Greens functions onto ue vector for convolution with impedance matrix
!
!     INPUT      MEANING
!     -----      ------- 
!     IDOFSE     maps bielak points to global DOFs 
!     ISRC       current source number
!     FREQ       current frequency
!     MDIM       leading dimension for IDOFSE 
!     NDIM       number of components in solution
!     NDOF       number of degrees of freedom
!     NNPE       number of nodal points in bielak layer
!     SRCTYP     source type; surface or body wave
!     STF        source time function
!  
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     UE         1D solution on vector of length ndof
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(2), INTENT(IN) :: SRCTYP
      COMPLEX*8, INTENT(IN) :: STF
      REAL*8, INTENT(IN) :: FREQ
      INTEGER*4, INTENT(IN) :: IDOFSE(MDIM,*), MDIM, NNPE, NDOF, NDIM, ISRC
      COMPLEX*8, INTENT(OUT) :: UE(NDOF) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*8, ALLOCATABLE :: UGRNS(:), VGRNS(:), WGRNS(:) 
      LOGICAL*4, ALLOCATABLE :: LNPINI(:)
      INTEGER*4 JOB
!
!----------------------------------------------------------------------------------------!
!
!.... load greens functions
      ALLOCATE(UGRNS(NNPE))
      ALLOCATE(VGRNS(NNPE))
      ALLOCATE(WGRNS(NNPE)) 
      IF (SRCTYP(1:1) == 'P') THEN
         JOB = 1 
      ELSE
         JOB = 2
      ENDIF 
      CALL READ_GRNS(MDIM, NNPE,NDIM,JOB,ISRC,.FALSE., FREQ, IDOFSE, &
                     UGRNS,VGRNS,WGRNS,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_grns25: There was an error reading the Greens functions!'
         RETURN
      ENDIF 

!.... stick onto ue
      ALLOCATE(LNPINI(1:NNPE)) 
      LNPINI(1:NNPE) = .TRUE.
      CALL FILL_GRNS(MDIM,NDOF,NNPE, NDIM,LNPINI,IDOFSE, UGRNS,VGRNS,WGRNS, UE,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_grns25: There was a warning when filling ue!'
         IERR = 0
      ENDIF
      DEALLOCATE(LNPINI)
      DEALLOCATE(UGRNS)
      DEALLOCATE(VGRNS)
      DEALLOCATE(WGRNS) 
!.... convolve source time function
      IF (STF /= CMPLX(1.0,0.0)) CALL CSCAL(NDOF,STF,UE,1)  
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GRNS_BDY(SRCTYP, MDIM,NDOF,NNPE, LSAVE,                    &
                          NL1D_LT,NL1D_RT, CSIDE, NDIM, ISRC,               &
                          FREQ,FREQ0,AOI,BAZN,                              &
                          XBLKL,XBLKR,ZBASE_INT,                   &
                          XMOD0,XMOD1, STF, IDOFSE,                         &
                          VP1D_LT,VS1D_LT,RH1D_LT,  &
                          QP1D_LT,QS1D_LT,Z1D_LT,   &
                          VP1D_RT,VS1D_RT,RH1D_RT,  &
                          QP1D_RT,QS1D_RT,Z1D_RT,   &
                          XLOCSE,ZLOCSE, UE,IERR)
!
!     Calculates the Green's functions at nodal points in the 'E' 
!     domain.  
!
!     Also note that the convention is Z (or w) positive up in finite 
!     elements.  If the model is striking east-west, (i.e., azmod = 0),
!     then  X is positive North, and Y is positive East.  Hence if the 
!     back azimuth (baz) is 180 degrees the incoming wave will be 
!     advancing in the +X (pure north) direction. In general, the 
!     slowness vector (px, py, pz) will be:
!        (sin(ai)cos(baz-pi), sin(ai)sin(baz-pi), cos(ai))/v
!     The azimuth of the model is defined by the direction of the +X 
!     axis relative to north.  Hence for a non-zero azmod we modify 
!     above form for (px,py,pz) to be
!        (sin(ai)cos(baz-pi-azmod), sin(ai)sin(baz-pi-azmod), cos(ai))/v
!     However, bazn has been corrected when read in 
!
!     INPUT      MEANING
!     -----      ------- 
!     AOI        angle of incidence degrees
!     BAZN       corrected back azimuth, degrees (see above) 
!     CSIDE      model side to take as 1D base model
!     IDOFSE     DOF pointer for nodes elements in Bielak domain
!     ISRC       source number
!     FREQ       frequency of interest (Hz)
!     FREQ0      reference frequency for dispersion
!     LSAVE      True -> save Greens fns in bielak boundary
!                False -> return on ue vector
!     MDIM       leading dimension for IDOFSE
!     NDIM       number of spatial dimensions
!     NDOF       number of degrees of freedom
!     NL1D_LT    number of points in left 1D model
!     NL1D_RT    number of points in right 1D model
!     NNPE       number of nodal points in Bielak domain
!     RH1D_LT    density 1d model on left
!     RH1D_RT    density 1d model on right
!     SRCTYP     source type
!     STF        source time function to convolve
!     QP1D_LT    1D P quality factor left
!     QS1D_LT    1D S quality factor left
!     QP1D_RT    1D P quality factor right
!     QS1D_RT    1D S quality factor right
!     VP1D_LT    vp 1d model on left
!     VP1D_RT    vp 1d model on right
!     VS1D_LT    vs 1d model on left
!     VS1D_RT    vs 1d model on right
!     XBLKL      left bielak/internal boundary
!     XBLKR      right bielak/internal boundary
!     XLOCSE     x locations of points in Bielak domain
!     XMOD0      left Bielak/absorbing boundary position in x
!     XMOD1      right Bielak/absorbing boundary position in x
!     ZBAES_INT  bottom of interior model
!     ZLOCSE     z locations of points in Bielak domain
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     UE         Greens functions at nodal points in Bielak domain (LSAVE = FALSE) 
! 
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(2), INTENT(IN) :: SRCTYP
      CHARACTER(1), INTENT(IN) :: CSIDE
      COMPLEX*16, INTENT(IN) :: STF
      REAL*8, INTENT(IN) :: XLOCSE(NNPE), ZLOCSE(NNPE),                                  &
              VP1D_LT(NL1D_LT), VS1D_LT(NL1D_LT), RH1D_LT(NL1D_LT), Z1D_LT(NL1D_LT),     &
              VP1D_RT(NL1D_RT), VS1D_RT(NL1D_RT), RH1D_RT(NL1D_RT), Z1D_RT(NL1D_RT),     &
              QP1D_RT(NL1D_RT), QP1D_LT(NL1D_LT), QS1D_RT(NL1D_RT), QS1D_LT(NL1D_LT),    &
              FREQ,FREQ0,AOI,BAZN, XBLKL,XBLKR, XMOD0,XMOD1, ZBASE_INT
      INTEGER*4, INTENT(IN) :: IDOFSE(MDIM,*), MDIM,NDOF,NNPE,    &
                               NL1D_LT,NL1D_RT, NDIM, ISRC
      LOGICAL*4, INTENT(IN) :: LSAVE 
      COMPLEX*8, INTENT(OUT) :: UE(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CLP(:)
      COMPLEX*16, ALLOCATABLE :: UGRN1D(:), WGRN1D(:)
      COMPLEX*8, ALLOCATABLE :: UGRNS(:), VGRNS(:), WGRNS(:) 
      REAL*8, ALLOCATABLE :: ZPTS1D(:) 
      LOGICAL*4, ALLOCATABLE :: LNPINI(:)
      COMPLEX*16 CCAZ, CSAZ, U,V,W
      REAL*8 POFF, OMEGA, PX, PY, CAZ, SAZ, XOFF, YOFF, ARG,    &
             TOL, TWOPI, PI180, VPBASE,VSBASE,     &
             SIGNPX, SIGNPY
      INTEGER*4 NNP1D, INP1D,INPE
      INTEGER*4 LOCATE8
      LOGICAL*4 LFLIP, LDISP
      PARAMETER(TOL = 1.11D-7)
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
      PARAMETER(LFLIP = .TRUE.)
      PARAMETER(YOFF = 0.D0) !we say that  y is in the plane
!
!----------------------------------------------------------------------------------------!
!
!.... input errors
      IERR = 0
      IF (SRCTYP(2:2) == 'H') THEN
         WRITE(*,*) 'grns_bdy: Error srctyp:',SRCTYP,' not yet programmed'
         IERR = 1
         RETURN
      ENDIF
!.... initialize 
      ALLOCATE(UGRNS(NNPE))
      ALLOCATE(VGRNS(NNPE))
      ALLOCATE(WGRNS(NNPE)) 
      UGRNS(1:NNPE) = CMPLX(0.,0.) 
      VGRNS(1:NNPE) = CMPLX(0.,0.)
      WGRNS(1:NNPE) = CMPLX(0.,0.)
      ALLOCATE(LNPINI(NNPE))
      LNPINI(1:NNPE) = .FALSE. 
      ALLOCATE(ZPTS1D(NNPE+1))
      ALLOCATE(CLP(NNPE)) 
!
!.... set phase shift offset based on direction of propagation
      IF (CSIDE == 'L') THEN !wave moving in +x, left is zero time 
         POFF = XMOD0
         SIGNPX = 1.D0
      ELSE !wave moving in -x, right is zero time
         POFF = XMOD1
         SIGNPX =-1.D0
      ENDIF
      SIGNPY = 1.D0 
!     IF (BAZN >= 0.D0 .AND. BAZN <= 180.D0) THEN !lower half plane +y
      IF (BAZN >= 180.D0 .AND. BAZN <= 360.D0) THEN !Lower half plane; +y
!        SIGNPY = 1.D0 !moving +y in finite elements 
      ELSE !upper half plane
!        SIGNPY =-1.D0 !moving -y in finite elements
      ENDIF 
!     CAZ = SIGNPX*DABS(DCOS(BAZN*PI180)) !>0 moving left to right
!     SAZ = SIGNPY*DABS(DSIN(BAZN*PI180)) !>0 moving towards observer
      CAZ = SIGNPX*DABS(DCOS(BAZN*PI180)) !>0 moving left to right
      SAZ = DSIN(BAZN*PI180)
      IF (DABS(CAZ) < 1.11D-15) CAZ = 0.D0
      IF (DABS(SAZ) < 1.11D-15) SAZ = 0.D0
      CCAZ = DCMPLX(CAZ,0.D0)
      CSAZ = DCMPLX(SAZ,0.D0)
      OMEGA = TWOPI*FREQ
!
!.... extract baesment velocity
      IF (CSIDE == 'L') THEN
         VPBASE = VP1D_LT(NL1D_LT)
         VSBASE = VS1D_LT(NL1D_LT)
      ELSE
         VPBASE = VP1D_RT(NL1D_RT)
         VSBASE = VS1D_RT(NL1D_RT)
      ENDIF
!
!.... calculate py
      IF (SRCTYP(2:2) == 'P') THEN
         PX = DSIN(AOI*PI180)*CAZ/VPBASE
         PY = DSIN(AOI*PI180)*SAZ/VPBASE
      ELSE
         PX = DSIN(AOI*PI180)*CAZ/VSBASE
         PY = DSIN(AOI*PI180)*SAZ/VSBASE
      ENDIF

      LDISP = .TRUE. !assume dispersion
      IF (FREQ0 == 0.D0) LDISP = .FALSE.
!
!.... generate the 1D solutions for the appropriate model
      IF (CSIDE == 'L') THEN
         CLP(1:NNPE) = 'L'
         CALL SET_GRNS1D(NNPE,1,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D)
         ALLOCATE(UGRN1D(NNP1D))
         ALLOCATE(WGRN1D(NNP1D))
         IF (LDISP) THEN
            CALL HASKATTN(NL1D_LT,NNP1D,SRCTYP,LFLIP, FREQ,FREQ0,AOI,      &
                          Z1D_LT,VP1D_LT,VS1D_LT,RH1D_LT,QP1D_LT,QS1D_LT,  &
                          ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ELSE
            CALL HASKGRN(NL1D_LT,NNP1D,SRCTYP,LFLIP, FREQ,AOI, Z1D_LT,VP1D_LT, &
                         VS1D_LT,RH1D_LT,ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ENDIF 
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'grns_bdy: Error calling haskgrn 1'
            RETURN
         ENDIF
      ELSE
         CLP(1:NNPE) = 'R'
         CALL SET_GRNS1D(NNPE,2,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D)
         ALLOCATE(UGRN1D(NNP1D))
         ALLOCATE(WGRN1D(NNP1D))
         IF (LDISP) THEN
            CALL HASKATTN(NL1D_RT,NNP1D,SRCTYP,LFLIP, FREQ,FREQ0,AOI,      &    
                          Z1D_RT,VP1D_RT,VS1D_RT,RH1D_RT,QP1D_RT,QS1D_RT,  &
                          ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ELSE
            CALL HASKGRN(NL1D_RT,NNP1D,SRCTYP,LFLIP, FREQ,AOI, Z1D_RT,VP1D_RT, &
                         VS1D_RT,RH1D_RT,ZPTS1D(1:NNP1D), UGRN1D,WGRN1D, IERR)
         ENDIF
         IF (IERR /= 0) THEN
            WRITE(*,*) 'grns_bdy: Error calling haskgrn 1'
            RETURN
         ENDIF
      ENDIF
      DEALLOCATE(CLP) 
!
!.... set points
      DO 100 INPE=1,NNPE
         INP1D = LOCATE8(NNP1D,TOL,ZLOCSE(INPE),ZPTS1D)
         IF (INP1D < 1) THEN
            WRITE(*,*) 'grns_bdy: Could not locate point!'
            IERR = 1
            GOTO 55
         ENDIF
         LNPINI(INPE) = .TRUE.
!
!....... calculate offset and phase shifter
         XOFF = XLOCSE(INPE) - POFF
         ARG =-OMEGA*(XOFF*PX + YOFF*PY) !-omega*p.x
!
!....... rotate for (u,v,w) response
         U = UGRN1D(INP1D)*CCAZ !u component, pure radial
         V = UGRN1D(INP1D)*CSAZ !v component, pure transverse
         W = WGRN1D(INP1D)      !w component, pure vertical
!
!....... convolve STF and shift
         U = U*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
         V = V*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
         W = W*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve 
         IF (SRCTYP(2:2) == 'P' .OR. SRCTYP(2:2) == 'S') THEN !incoming P/SV
            UGRNS(INPE) = CMPLX(U)
            VGRNS(INPE) = CMPLX(V)
         ELSE !incoming SH polarization
            UGRNS(INPE) =-CMPLX(V)
            VGRNS(INPE) = CMPLX(U)
         ENDIF
         WGRNS(INPE) =-CMPLX(W)!point up now 
  100 CONTINUE
   55 CONTINUE
!
!.... fill in the DOFs  or save ?
      IF (.NOT.LSAVE) THEN 
         CALL FILL_GRNS(MDIM,NDOF,NNPE, NDIM,LNPINI,IDOFSE, UGRNS,VGRNS,WGRNS, UE,IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'grns_bdy: Serious warning in fill_grns'
            IERR = 0
         ENDIF
      ELSE 
         DO 45 INPE=1,NNPE
            IF (.NOT.LNPINI(INPE)) & 
            WRITE(*,*) 'grns_bdy: Serious warning in grns_bdy, unitialized point',INPE
   45    CONTINUE 
         CALL SAVE_GRNS(NDIM, NNPE,NDIM,1,ISRC,.FALSE., FREQ,IDOFSE,UGRNS,VGRNS,WGRNS) 
      ENDIF
      DEALLOCATE(UGRNS)
      DEALLOCATE(WGRNS) 
      DEALLOCATE(UGRN1D)
      DEALLOCATE(WGRN1D)
      DEALLOCATE(ZPTS1D)
      RETURN
      END
!
!.... classify the points as left or right for this source
!     ALLOCATE(CLP(NNPE))
!     CALL CLASS_NNPE(NNPE,CSIDE,TOL,XBLKL,XBLKR, XLOCSE, CLP,IERR)
!     IF (IERR /= 0) THEN
!        WRITE(*,*) 'grns_bdy: Error calling class_nnpe!'
!        RETURN
!     ENDIF
!
!.... calculate the distance between the base model bases
!     IF (CSIDE == 'L') THEN 
!        XOFF = DABS(XBLKR - POFF)
!     ELSE !wave moving right to left
!        XOFF = XBLKL - POFF
!     ENDIF
!     IF (SRCTYP(2:2) == 'P') THEN
!        PX = DSIN(AOI*PI180)*CAZ/VPBASE
!        PY = DSIN(AOI*PI180)*SAZ/VPBASE
!     ELSE
!        PX = DSIN(AOI*PI180)*CAZ/VSBASE
!        PY = DSIN(AOI*PI180)*SAZ/VSBASE
!     ENDIF
!     ZBASE_BLK = MINVAL(ZLOCSE) !+z up
!     CALL HASKGRN(NL1D_LT,1,SRCTYP,LFLIP, FREQ,AOI, Z1D_LT,VP1D_LT,      &
!                  VS1D_LT,RH1D_LT,ZBASE_BLK, UL,WL, IERR)
!     IF (IERR /= 0) THEN
!        WRITE(*,*) 'grns_driver: Error calling haskgrn 1'
!        RETURN
!     ENDIF
!     CALL HASKGRN(NL1D_RT,1,SRCTYP,LFLIP, FREQ,AOI, Z1D_RT,VP1D_RT,      &
!                  VS1D_RT,RH1D_RT,ZBASE_BLK, UR,WR, IERR)
!     IF (IERR /= 0) THEN
!        WRITE(*,*) 'grns_driver: Error calling haskgrn 2'
!        RETURN
!     ENDIF
!     print *, xoff, px
!     ARG =-OMEGA*(XOFF*PX + YOFF*PY)
!     print *, px,py
!     IF (CSIDE == 'L') THEN 
!        US = UL*CDEXP(DCMPLX(0.D0,ARG)) !shift the left model
!        WS = WL*CDEXP(DCMPLX(0.D0,ARG)) !shift, but keep + down 
!        U  = UR
!        W  = WR !no negative, keep + down
!     ELSE
!        US = UR*CDEXP(DCMPLX(0.D0,ARG)) !shift the right model
!        WS = WR*CDEXP(DCMPLX(0.D0,ARG)) !shift, but keep + down
!        U  = UL
!        W  = WL !no negative, keep + down 
!     ENDIF
!
!.... need base phase at phase shifted signals so shift phase at base dphi
!     DPHIU = DATAN2(DIMAG(US),DREAL(US)) - DATAN2(DIMAG(U),DREAL(U))
!     DPHIW = DATAN2(DIMAG(WS),DREAL(WS)) - DATAN2(DIMAG(W),DREAL(W))
!
!.... make the 1D models
!     CALL SET_GRNS1D(NNPE,1,CLP,TOL,ZLOCSE, NNP1D_LT,ZPTS1D)
!     ALLOCATE(U1DL(NNP1D_LT))
!     ALLOCATE(W1DL(NNP1D_LT))
!     CALL HASKGRN(NL1D_LT,NNP1D_LT,SRCTYP,LFLIP, FREQ,AOI, Z1D_LT,VP1D_LT, &
!                  VS1D_LT,RH1D_LT,ZPTS1D(1:NNP1D_LT), U1DL,W1DL, IERR)
!     IF (IERR /= 0) THEN
!        WRITE(*,*) 'grns_bdy: Error calling haskgrn 3'
!        RETURN
!     ENDIF
!
!.... make the 1D models
!     CALL SET_GRNS1D(NNPE,2,CLP,TOL,ZLOCSE, NNP1D_RT,ZPTS1D)
!     ALLOCATE(U1DR(NNP1D_RT))
!     ALLOCATE(W1DR(NNP1D_RT))
!     CALL HASKGRN(NL1D_RT,NNP1D_RT,SRCTYP,LFLIP, FREQ,AOI, Z1D_RT,VP1D_RT, &
!                  VS1D_RT,RH1D_RT,ZPTS1D(1:NNP1D_RT), U1DR,W1DR, IERR)
!     IF (IERR /= 0) THEN
!        WRITE(*,*) 'grns_bdy: Error calling haskgrn 4'
!        RETURN
!     ENDIF
!
!.... shift the model on the other side
!     IF (CSIDE == 'L') THEN
!        ZSHIFT = CDEXP(DCMPLX(0.D0,DPHIU))
!        CALL ZSCAL(NNP1D_RT,ZSHIFT,U1DR,1)
!        ZSHIFT = CDEXP(DCMPLX(0.D0,DPHIW)) 
!        CALL ZSCAL(NNP1D_RT,ZSHIFT,W1DR,1)
!     ELSE
!        ZSHIFT = CDEXP(DCMPLX(0.D0,DPHIU))
!        CALL ZSCAL(NNP1D_LT,ZSHIFT,U1DL,1)
!        ZSHIFT = CDEXP(DCMPLX(0.D0,DPHIW))
!        CALL ZSCAL(NNP1D_LT,ZSHIFT,W1DL,1)
!     ENDIF
!     if (cside == 'L') then
!     print *, 'dp',dphiu,dphiw
!     print *, 'a',datan2(dimag(us),dreal(us)),datan2(dimag(u1dr(1)),dreal(u1dr(1)))
!     print *, 'a',datan2(dimag(ws),dreal(ws)),datan2(dimag(w1dr(1)),dreal(w1dr(1)))
!     else
!     print *, 'dp',dphiu,dphiw
!     print *, 'b',datan2(dimag(us),dreal(us)),datan2(dimag(u1dl(1)),dreal(u1dl(1)))
!     print *, 'b',datan2(dimag(ws),dreal(ws)),datan2(dimag(w1dl(1)),dreal(w1dl(1)))
!     endif
!     !pause
!
!.... fill the greens fns at each point
!     DO 2 ISIDE=1,2
!        CALL SET_GRNS1D(NNPE,ISIDE,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D) 
!        ALLOCATE(UGRN1D(NNP1D))
!        ALLOCATE(WGRN1D(NNP1D))
!        IF (ISIDE == 1) THEN !left, copy left 1D greens fns 
!           CALL ZCOPY(NNP1D,U1DL,1,UGRN1D,1)
!           CALL ZCOPY(NNP1D,W1DL,1,WGRN1D,1)
!        ELSE !right, copy 1D greens fns
!           CALL ZCOPY(NNP1D,U1DR,1,UGRN1D,1)
!           CALL ZCOPY(NNP1D,W1DR,1,WGRN1D,1)
!        ENDIF 
!        INP1D = 1
!        DO 11 INPE=1,NNPE
!           LCHECK = .FALSE.
!           IF (ISIDE == 1 .AND. CLP(INPE) == 'L') LCHECK = .TRUE.
!           IF (ISIDE == 2 .AND. CLP(INPE) == 'R') LCHECK = .TRUE.
!           IF (LCHECK) THEN
!              INP1D = LOCATE8(NNP1D,TOL,ZLOCSE(INPE),ZPTS1D) 
!              IF (INP1D < 1) THEN
!                 WRITE(*,*) 'gengrns: Could not locate point!'
!                 IERR = 1
!                 GOTO 50
!              ENDIF
!              IF (LNPINI(INPE)) &
!              WRITE(*,*) 'gengrns: Warning point already initialized'
!              LNPINI(INPE) = .TRUE. 
!              IF (CSIDE == 'L' .AND. ISIDE == 1 .OR. &
!                  CSIDE == 'R' .AND. ISIDE == 2) THEN
!                 IF (CSIDE == 'L') THEN
!                    XOFF = DABS(XLOCSE(INPE) - POFF) !+(x - x_0) =-(x_0 - x) 
!                 ELSE
!                    XOFF = XLOCSE(INPE) - POFF
!                 ENDIF
!              ELSE
!                 IF (CSIDE == 'L') THEN
!                    XOFF = DABS(XLOCSE(INPE) - XBLKR)
!                 ELSE
!                    XOFF = XLOCSE(INPE) - XBLKL
!                 ENDIF
!              ENDIF
!              !this is a wave number vector where omega*p = k
!              !haskell has given us the response at depths  
!              !so we don't need z
!              ARG =-OMEGA*(XOFF*PX + YOFF*PY) !-omega*p.x
!              U = UGRN1D(INP1D)*CCAZ !u component, pure radial
!              V = UGRN1D(INP1D)*CSAZ !v component, pure transverse
!              W = WGRN1D(INP1D)      !w component, pure vertical
!
!............. match phase at model base
!              IF (CSIDE == 'L' .AND. ISIDE == 2 .OR. &
!                  CSIDE == 'R' .AND. ISIDE == 1) THEN
!                 U = U*CDEXP(DCMPLX(0.D0,ARG+DPHIU))*STF 
!                 V = V*CDEXP(DCMPLX(0.D0,ARG))*STF
!                 W = W*CDEXP(DCMPLX(0.D0,ARG+DPHIW))*STF
!              ELSE
!                 U = U*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
!                 V = V*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve
!                 W = W*CDEXP(DCMPLX(0.D0,ARG))*STF!phase shift/convolve 
!              ENDIF
!              IF (SRCTYP(2:2) == 'P' .OR. SRCTYP(2:2) == 'S') THEN !incoming P/SV
!                 UGRNS(INPE) = CMPLX(U) 
!                 VGRNS(INPE) = CMPLX(V)
!              ELSE !incoming SH polarization
!                 UGRNS(INPE) =-CMPLX(V)
!                 VGRNS(INPE) = CMPLX(U)
!              ENDIF
!              WGRNS(INPE) =-CMPLX(W)!point up now 
!           ENDIF
!  11    CONTINUE !loop on nodal points
!        IF (ALLOCATED(UGRN1D)) DEALLOCATE(UGRN1D)
!        IF (ALLOCATED(WGRN1D)) DEALLOCATE(WGRN1D) 
!   2 CONTINUE !loop on sides 
!
!.... fill in the DOFs  or save ?
!     IF (.NOT.LSAVE) THEN 
!        CALL FILL_GRNS(MDIM,NDOF,NNPE, NDIM,LNPINI,IDOFSE, UGRNS,VGRNS,WGRNS, UE,IERR)
!        IF (IERR /= 0) THEN
!           WRITE(*,*) 'gengrns: Serious warning in fill_grns'
!           IERR = 0
!        ENDIF
!     ELSE 
!        DO 30 INPE=1,NNPE
!           IF (.NOT.LNPINI(INPE)) & 
!           WRITE(*,*) 'gengrns: Serious warning in grns_bdy, unitialized point',INPE
!  30    CONTINUE 
!        CALL SAVE_GRNS(NDIM, NNPE,NDIM,1,ISRC,.FALSE., FREQ, IDOFSE,UGRNS,VGRNS,WGRNS) 
!     ENDIF
!  50 CONTINUE
!     DEALLOCATE(U1DL)
!     DEALLOCATE(W1DL)
!     DEALLOCATE(U1DR)
!     DEALLOCATE(W1DR)
!     DEALLOCATE(UGRNS)
!     DEALLOCATE(VGRNS)
!     DEALLOCATE(WGRNS) 
!     DEALLOCATE(CLP) 
!     DEALLOCATE(LNPINI)
!     RETURN
!     END
!                                                                                        ! 
!========================================================================================!
!
      SUBROUTINE GRNS_SRF(TMPDIR, LNSRF, &
                          MYID,SRCTYP, MDIM,NDOF,NNPE, LSAVE,LSHIFT,     &
                          NL1D_LT,NL1D_RT, CSIDE, NDIM,                          &
                          MODE, ISRC,FREQ0,FREQ,BAZN,                            &
                          SLAT,SLON,SDEP,SMAG,                                   &
                          STRIKE,DIP,RAKE,                                       &
                          XBLKL,XBLKR, XMOD0,XMOD1,                              &
                          XMLAT0,XMLON0,XMLAT1,XMLON1,  &
                          STF, IDOFSE, &
                          VPD_RLLT,VSD_RLLT,ROD_RLLT,HDD_RLLT,                   &
                          VPD_LVLT,VSD_LVLT,ROD_LVLT,HDD_LVLT,                   &
                          VPD_RLRT,VSD_RLRT,ROD_RLRT,HDD_RLRT,                   &
                          VPD_LVRT,VSD_LVRT,ROD_LVRT,HDD_LVRT,                   &
                          QP1D_LT,QP1D_RT, QS1D_LT,QS1D_RT,                      &
                          Z1D_LT,Z1D_RT, XLOCSE,ZLOCSE,                          &
                          PERIOD,CCRAY,CCLOV, UE,IERR)
!
!     Generates the surface wave Greens functions at model DOFs in the Bielak layer.
!
!     The surface waves are a mess.  Seriously, this really is not well programmed. 
!     Instead of piecemealing information from subprogram to subprogram via file IO 
!     I should estimate space requirements and return the results of calculation.  
!     I would fix this but I literally have no idea how things are being written so 
!     I can't predict the space needed.  Moreover, I really don't know what in the 
!     common blocks is important.  I'm just trying to get this work for Rayleigh 
!     fundamental modes.  You're on your own after. - B Baker January 2012
!
!     INPUT      MEANING
!     -----      ------- 
!     BAZN       corrected back-azimuth (degrees) 
!     CSIDE      'L' wave approaching from left, 'R' wave aproaching from right
!     DIP        source dip angle (degrees)
!     ISRC       current source number
!     FREQ       current frequency (Hz)
!     FREQ0      reference frequency
!     HDD_???T   Depth interfaces in Rayleigh/Love model on left/right
!     LSAVE      True -> save the Greens functions to disk
!     LSHIFT     True -> shift phase of argument
!     MDIM       leading dimension for IDOFSE
!     MODE       surface wave mode to model
!     MYID       process ID for file IO
!     NDIM       number of components in solution
!     NDOF       number of degrees of freedom
!     NL1D_?T    number of layers in 1D model on left or right
!     NNPE       number of nodes in bielak layer
!     RAKE       source rake angle (degrees)
!     ROD_???T   Density 1D model fo Rayleigh/Love waves on left/right
!     SDEP       source depth (km)
!     SLAT       source latitude (degrees)
!     SLON       source longitude (degrees)
!     SMAG       source magnitude (dyne-cm) 
!     SRCTYP     type of source, SV(ertical), SR(ayleigh), SB(oth) or SL(ove)
!     STF        source time function to convolve 
!     STRIKE     source strike angle (degrees)
!     TMPDIR     temporary directory for files
!     VSD_???T   Vs 1D model for Rayleigh/Love waves on left/right
!     VPD_???T   Vp 1D model for Rayleigh/Love waves on left/right
!     XBLKL      left x position of bielak/internal boundary
!     XBLKR      right x position of bielak/internal boundary
!     XLOCSE     x locations of bielak nodes
!     XMLAT0     latitude of XBLKL
!     XMLAT1     latitude of XBLKR
!     XMLON0     longitude of XBLKL 
!     XMLAT1     longitude of XBLKR 
!     XMOD0      left bielak/absorbing position in x
!     XMOD1      right bielak/absorbing position in x
!     Z1D_?T     z locations in 1D model on left/right
!     ZLOCSE     z locations of bielak nodes 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     CCLOV      phase velocity (km/s) for love waves at modes
!     CCRAY      phase velocity (km/s) for rayleigh waves at modes
!     IERR       error flag
!     PERIOD     1/frequency
!     UE         if LSAVE = False the 1D solution in bielak layer
!  
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: TMPDIR
      CHARACTER(2), INTENT(IN) :: SRCTYP
      CHARACTER(1), INTENT(IN) :: CSIDE
      COMPLEX*16, INTENT(IN) :: STF 
      REAL*8, INTENT(IN) :: XLOCSE(NNPE), ZLOCSE(NNPE), &
              VPD_RLLT(NL1D_LT),VSD_RLLT(NL1D_LT),ROD_RLLT(NL1D_LT),HDD_RLLT(NL1D_LT),  &
              VPD_LVLT(NL1D_LT),VSD_LVLT(NL1D_LT),ROD_LVLT(NL1D_LT),HDD_LVLT(NL1D_LT),  &
              VPD_RLRT(NL1D_RT),VSD_RLRT(NL1D_RT),ROD_RLRT(NL1D_RT),HDD_RLRT(NL1D_RT),  &
              VPD_LVRT(NL1D_RT),VSD_LVRT(NL1D_RT),ROD_LVRT(NL1D_RT),HDD_LVRT(NL1D_RT),  &
              QP1D_LT(NL1D_LT), QP1D_RT(NL1D_RT), QS1D_LT(NL1D_LT), QS1D_RT(NL1D_RT),   &
              Z1D_LT(NL1D_LT), Z1D_RT(NL1D_RT), & 
              FREQ, XBLKL,XBLKR,      &
              XMOD0,XMOD1, XMLAT0,XMLON0, XMLAT1,XMLON1,                   & 
              SLAT,SLON,SDEP,SMAG, STRIKE,DIP,RAKE, BAZN, FREQ0 
      INTEGER*4, INTENT(IN) :: IDOFSE(MDIM,*), MODE,MDIM,NDOF,NNPE,    &   
                               NL1D_LT,NL1D_RT, NDIM, MYID, ISRC
      LOGICAL*4, INTENT(IN) :: LSAVE, LSHIFT, LNSRF
      COMPLEX*8, INTENT(OUT) :: UE(*) 
      REAL*8, INTENT(OUT) :: CCRAY(*), CCLOV(*), PERIOD 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CLP(:)
      COMPLEX*16, ALLOCATABLE :: UGRN1D(:), VGRN1D(:), WGRN1D(:)
      COMPLEX*8, ALLOCATABLE :: UGRNS(:), VGRNS(:), WGRNS(:)
      REAL*8, ALLOCATABLE :: VPD_RL(:), VSD_RL(:), ROD_RL(:), HDD_RL(:), &
                             VPD_LV(:), VSD_LV(:), ROD_LV(:), HDD_LV(:), &
                             QA(:), QB(:), ZPTS1D(:), XLAYR(:), XDIST(:), Z1D(:), &
                             AZIMS(:)   
      INTEGER*4, ALLOCATABLE :: MYNNPE(:) 
      LOGICAL*4, ALLOCATABLE :: LNPINI(:)
      CHARACTER(80) PHASEFL, DISPFL, RAYFL, LOVFL, PHASE_RAY, PHASE_LOV  
      CHARACTER(10) CIL 
      CHARACTER(5) CID
      COMPLEX*16 CCBAZ, CSBAZ, U,V,W
      REAL*8 PERRAYW(20), PERLOVW(20), CCRAYW(20), CCLOVW(20),                       &
             POFF, OMEGA, CBAZ, SBAZ,                                                &
             VPBASE_RL,VSBASE_RL,RHBASE_RL, VPBASE_LV,VSBASE_LV,RHBASE_LV,           &
             ZBASE_RL, ZBASE_LV, ZTOP_RL, ZTOP_LV, ZKM, PX, PY, XOFF, YOFF, DSTDEG,  &
             DELT0,AZ0, BAZ0, ARG, cc_save
      INTEGER*4 IEVTIME(4), NL, NNP1D, NSRCIN, NFREQIN, NBRAN2P, JCOM, &
                INPE,JNPE,ISIDE,IL,IZ1,IZ2,IXPTS,IEFL, KNPE,IZINTER,NXPTS, &
                IORD1,IORD, IBSECT8
      LOGICAL*4 LCHECK,LEX,LDEL
      REAL*8 SSEC,QPP,QSS,C1,C2,TREF,TOL,TWOPI,PI180,T0
      INTEGER*4 MDMIN,MDMAX,NBRAN1,NBRAN2
      LOGICAL*4 LPRATT,LBIN,LFLIP 
      PARAMETER(YOFF = 0.D0) !working in plane
      PARAMETER(QPP = 600.D0, QSS = 300.D0) !Qp and Qs quality factors 
      PARAMETER(C1  = 0.D0  , C2 = 0.D0) !min and max phase velocities;earthsr 
                                         !c1 = c2 = 0 -> earthsr picks phase velocities
      PARAMETER(NBRAN1 = 0  ) 
      !PARAMETER(NBRAN2 = 4)  !min and max modes;earthsr should be an input
      PARAMETER(IEFL = 0)     !1 -> apply earth flattening correction in earthsr
                              !0 -> do not.  pry use this b/c our mesh is in a
                              !flat model
      !PARAMETER(TREF = 10.D0) !reference period for material dispersion correction;earhtsr
      PARAMETER(T0   = 0.D0)  !start time if debugging srgramf 
      PARAMETER(LPRATT = .FALSE.) !don't use pratt's convention on Fourier transform
      PARAMETER(LFLIP = .TRUE.) !flip the coordinate system
      PARAMETER(LBIN = .TRUE.) !want binary files in surface wave computation
      PARAMETER(TOL = 1.11D-7) !closer than this means its the same point
      PARAMETER(MDMIN = 0) !min number of modes
      !PARAMETER(MDMAX = 10) !min/max modes, but not equal to nbran?;earthsr
      PARAMETER(TWOPI = 6.2831853071795862D0) !2*pi
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180 
      PARAMETER(NSRCIN = 1) !only do one source at a time here
      PARAMETER(NFREQIN = 1) !only one frequency at a time
      DATA IEVTIME/2008,1,0,0/ !year, julian day, hour minute
      DATA SSEC/0.D0/  
!
!----------------------------------------------------------------------------------------!
!
!.... initialization
      IERR = 0
      CID(1:5) = ' '
      WRITE(CID,'(I5)') MYID
      CID = ADJUSTL(CID) 
      ALLOCATE(UGRN1D(NNPE)) !this is really a worst case
      ALLOCATE(VGRN1D(NNPE)) 
      ALLOCATE(WGRN1D(NNPE))
      IF (.NOT.LSAVE) UE(1:NDOF) = CMPLX(0.0,0.0)
      ALLOCATE(UGRNS(NNPE))
      ALLOCATE(VGRNS(NNPE))
      ALLOCATE(WGRNS(NNPE)) 
      ALLOCATE(XLAYR(NNPE)) 
      ALLOCATE(MYNNPE(NNPE))
      UGRNS(1:NNPE) = CMPLX(0.,0.) 
      VGRNS(1:NNPE) = CMPLX(0.,0.)
      WGRNS(1:NNPE) = CMPLX(0.,0.) 
      ALLOCATE(LNPINI(NNPE))
      LNPINI(1:NNPE) = .FALSE. 
      ALLOCATE(ZPTS1D(NNPE))
!
!.... reference period
      IF (FREQ0 > 0.D0) THEN
         TREF = 1.D0/FREQ0
      ELSE
         TREF = 10.D0
      ENDIF
!
!.... adjust from absolute distance from bielak side 
      IF (CSIDE == 'L') THEN !wave moving in +x, left is zero time 
         POFF = XMOD0
      ELSE !wave moving in -x, right is zero time
         POFF = XMOD1
      ENDIF
      CBAZ = DCOS(BAZN*PI180) !> 0 moving left to right
      SBAZ = DSIN(BAZN*PI180)
      IF (DABS(CBAZ) < 1.11D-15) CBAZ = 0.D0
      IF (DABS(SBAZ) < 1.11D-15) SBAZ = 0.D0
      CCBAZ = DCMPLX(CBAZ,0.D0)
      CSBAZ = DCMPLX(SBAZ,0.D0)
      OMEGA = TWOPI*FREQ
!
!.... force base layers to match on sides 
      IF (CSIDE == 'L') THEN !wave moving +x
         VPBASE_RL = VPD_RLLT(NL1D_LT)
         VSBASE_RL = VSD_RLLT(NL1D_LT)
         RHBASE_RL = ROD_RLLT(NL1D_LT)
         ZBASE_RL  = HDD_RLLT(NL1D_LT)
         ZTOP_RL   = HDD_RLLT(1)
         VPBASE_LV = VPD_LVLT(NL1D_LT)
         VSBASE_LV = VSD_LVLT(NL1D_LT)
         RHBASE_LV = ROD_LVLT(NL1D_LT)
         ZBASE_LV  = HDD_LVLT(NL1D_LT)
         ZTOP_LV   = HDD_LVLT(1) 
         CALL VINCENTY(.FALSE., SLAT,SLON, XMLAT0,XMLON0,   &
                       DSTDEG,DELT0,AZ0,BAZ0,IERR) 
      ELSE !wave moving -x
         VPBASE_RL = VPD_RLRT(NL1D_RT)
         VSBASE_RL = VSD_RLRT(NL1D_RT)
         RHBASE_RL = ROD_RLRT(NL1D_RT)
         ZBASE_RL  = HDD_RLRT(NL1D_RT)
         ZTOP_RL   = HDD_RLRT(1) 
         VPBASE_LV = VPD_LVRT(NL1D_RT)
         VSBASE_LV = VSD_LVRT(NL1D_RT)
         RHBASE_LV = ROD_LVRT(NL1D_RT)
         ZBASE_LV  = HDD_LVRT(NL1D_RT)
         ZTOP_LV   = HDD_LVRT(1)
         CALL VINCENTY(.FALSE., SLAT,SLON, XMLAT1,XMLON1,   &
                       DSTDEG,DELT0,AZ0,BAZ0,IERR) 
      ENDIF
!
!.... set space for phase info
      NBRAN2 = MAX(MODE,0) 
      MDMAX = NBRAN2 
!
!.... classify the points as left, bottom, or right
      ALLOCATE(CLP(NNPE))
      CALL CLASS_NNPE(NNPE,CSIDE,TOL,XBLKL,XBLKR, XLOCSE, CLP,IERR)
      IF (CSIDE == 'L') THEN
         CALL SET_GRNS1D(NNPE,1,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D)
      ELSE
         CALL SET_GRNS1D(NNPE,2,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D)
      ENDIF
!
!.... fool the initialization check
      DO 1 INPE=1,NNPE
         IF (CSIDE /= CLP(INPE)) THEN 
            LNPINI(INPE) = .TRUE.
            UGRNS(INPE) = CMPLX(0.0,0.0)
            VGRNS(INPE) = CMPLX(0.0,0.0)
            WGRNS(INPE) = CMPLX(0.0,0.0)
         ENDIF
    1 CONTINUE 
!
!.... set the 1D model
      IF (CSIDE == 'L') THEN !set model to left 1D model
         NL = NL1D_LT
         ALLOCATE(VPD_RL(NL))
         ALLOCATE(VSD_RL(NL))
         ALLOCATE(ROD_RL(NL))
         ALLOCATE(HDD_RL(NL)) 
         ALLOCATE(VPD_LV(NL))
         ALLOCATE(VSD_LV(NL))
         ALLOCATE(ROD_LV(NL))
         ALLOCATE(HDD_LV(NL)) 
         ALLOCATE(QA(NL))
         ALLOCATE(QB(NL))
         ALLOCATE(Z1D(NL)) 
         VPD_RL(1:NL) = VPD_RLLT(1:NL) 
         VSD_RL(1:NL) = VSD_RLLT(1:NL)
         ROD_RL(1:NL) = ROD_RLLT(1:NL)
         HDD_RL(1:NL) = HDD_RLLT(1:NL)
         VPD_LV(1:NL) = VPD_LVLT(1:NL)
         VSD_LV(1:NL) = VSD_LVLT(1:NL)
         ROD_LV(1:NL) = ROD_LVLT(1:NL)
         HDD_LV(1:NL) = HDD_LVLT(1:NL)
         Z1D(1:NL) = Z1D_LT(1:NL)
         QA(1:NL) = QP1D_LT(1:NL) !P quality factor
         QB(1:NL) = QS1D_LT(1:NL) !S quality factor
      ELSE !set model to right 1D model
         NL = NL1D_RT
         ALLOCATE(VPD_RL(NL))
         ALLOCATE(VSD_RL(NL))
         ALLOCATE(ROD_RL(NL))
         ALLOCATE(HDD_RL(NL)) 
         ALLOCATE(VPD_LV(NL))
         ALLOCATE(VSD_LV(NL))
         ALLOCATE(ROD_LV(NL))
         ALLOCATE(HDD_LV(NL))
         ALLOCATE(QA(NL))
         ALLOCATE(QB(NL))
         ALLOCATE(Z1D(NL)) 
         VPD_RL(1:NL) = VPD_RLRT(1:NL)
         VSD_RL(1:NL) = VSD_RLRT(1:NL)
         ROD_RL(1:NL) = ROD_RLRT(1:NL)
         HDD_RL(1:NL) = HDD_RLRT(1:NL)
         VPD_LV(1:NL) = VPD_LVRT(1:NL)
         VSD_LV(1:NL) = VSD_LVRT(1:NL)
         ROD_LV(1:NL) = ROD_LVRT(1:NL)
         HDD_LV(1:NL) = HDD_LVRT(1:NL)
         Z1D(1:NL) = Z1D_RT(1:NL)
         QA(1:NL) = QP1D_RT(1:NL) !P quality factor
         QB(1:NL) = QS1D_RT(1:NL) !S quality factor
      ENDIF !end check on side
      VPD_RL(NL) = VPBASE_RL !force base to be consistent
      VSD_RL(NL) = VSBASE_RL 
      ROD_RL(NL) = RHBASE_RL  
      VPD_LV(NL) = VPBASE_LV !force base to be consistent
      VSD_LV(NL) = VSBASE_LV
      ROD_LV(NL) = RHBASE_LV
!
!.... generate a list for depths
      DO 105 IL=1,NNP1D !loop on layers
!
!....... generate an x location list
         INPE = IBSECT8(NNPE,1,TOL,ZPTS1D(IL),ZLOCSE)
         IF (INPE <= 0) THEN
            WRITE(*,*) 'gengrns: Error could not locate inpe!'
            IERR = 1
            GOTO 500
         ENDIF
         KNPE = INPE
         DO 106 JNPE=INPE,NNPE
            IF (ABS(ZLOCSE(JNPE) - ZPTS1D(IL)) > TOL) GOTO 160
            KNPE = JNPE
   106   CONTINUE
   160   CONTINUE 
         IZ1 = 1
         IZ2 = 0 
         DO 107 JNPE=INPE,KNPE
            LCHECK = .FALSE.
            IF (CSIDE == CLP(JNPE)) LCHECK = .TRUE.
            IF (LCHECK) THEN
               IZ2 = IZ2 + 1 
               MYNNPE(IZ2) = JNPE 
               XLAYR(IZ2) = XLOCSE(JNPE)
            ENDIF !end check on working
  107    CONTINUE !loop on nodal points in bielak layer
!        IF (IZ2 < 1) GOTO 150 !nothing to do
         IF (IZ2 < 1) THEN
            WRITE(*,*) 'grns_srf: Error couldnt locate any points!'
            IERR = 1
            RETURN
         ENDIF
         IF (IZ2 > NNPE) THEN
            WRITE(*,*) 'grns_srf: Warning iz2 > nnpe'
            IZ2 = NNPE
         ENDIF
         IF (IZ1 > IZ2) THEN
            WRITE(*,*) 'grns_srf: Error iz2 < iz1'
            IERR = 1 
            GOTO 500
         ENDIF
         NXPTS = IZ2 - IZ1 + 1
!
!....... set phase velocities, source depths, receivier depths for srgramf 
         NBRAN2P = NBRAN2
         JCOM = 1 !rayleigh
         IF (LFLIP) THEN
            ZKM = (Z1D(1) - ZPTS1D(IL))*1.D-3 !depth to interfaces 
         ELSE
            ZKM = ZPTS1D(IL)*1.D-3 
         ENDIF
         !print *, zkm,zpts1d(il)
         DISPFL(1:80) = ' '
         PHASEFL(1:80) = ' '
         RAYFL(1:80) = ' '
         CIL(1:10) = ' '
         WRITE(CIL,'(I10)') IL 
         CIL = ADJUSTL(CIL) 
         DISPFL = TRIM(ADJUSTL(TMPDIR))//'surf_ray_'//TRIM(CIL)//'-'//TRIM(CID)//'.out'
         IF (CSIDE == 'L') THEN
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'srf_ray_phvelL-'//TRIM(CIL)//'-'// &
                      TRIM(CID)//'.out'
         ELSE
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'srf_ray_phvelR-'//TRIM(CIL)//'-'// &
                      TRIM(CID)//'.out'
         ENDIF
         DISPFL = ADJUSTL(DISPFL)
         PHASEFL = ADJUSTL(PHASEFL)
         PHASE_RAY(1:80) = ' '
         PHASE_RAY = TRIM(ADJUSTL(PHASEFL))
         PHASE_RAY = ADJUSTL(PHASE_RAY) 
         RAYFL = TRIM(ADJUSTL(TMPDIR))//'ray-'//TRIM(CID)
         RAYFL = ADJUSTL(RAYFL) 
         LDEL = .TRUE. !delete the phase file?
         IF (IL == NNP1D) LDEL = .FALSE. 
         !IF (.NOT.LDEL) THEN
         !   PHASEFL(1:80) = ' '
         !   PHASEFL = 'surf_ray_'//TRIM(CID)//'.out'
         !   PHASEFL = ADJUSTL(PHASEFL) 
         !ENDIF
         !IF (LNEW) THEN
             INQUIRE(FILE=TRIM(PHASEFL),EXIST=LEX)
             IF (LEX) THEN
                OPEN(UNIT=66,FILE=TRIM(PHASEFL),STATUS='OLD')
                CLOSE(66,STATUS='DELETE')
             ENDIF
         !ENDIF
         CALL EARTHSR(NL,JCOM, NSRCIN,NFREQIN,IEFL,NBRAN1,NBRAN2P,C1,C2,          &
                      TREF,ZKM,FREQ,SDEP, HDD_RL,VPD_RL,VSD_RL,ROD_RL,QA,QB,      &
                      LBIN,LDEL,DISPFL,PHASEFL,RAYFL,                             &
                      PERRAYW,PERLOVW,CCRAYW,CCLOVW, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'gengrns: Error calling earthsr for Rayleigh waves' 
            GOTO 500
         ENDIF
!
!....... save phase velocity at surface?
         IF (IL == 1) THEN
!           IF (CSIDE == 'L' .AND. ISIDE == 1 .OR. &
!               CSIDE == 'R' .AND. ISIDE == 2) THEN
             DO 15 IORD=0,NBRAN2P
                IORD1 = IORD + 1 
                 IF (IORD == 0) PERIOD = PERRAYW(IORD1)
                 CCRAY(IORD1) = CCRAYW(IORD1)
   15        CONTINUE 
!           ENDIF
         ENDIF 
         cc_save = ccrayw(1)
!
!....... repeat for love waves
         NBRAN2P = NBRAN2
         JCOM = 2 !love
         DISPFL(1:80) = ' '
         PHASEFL(1:80) = ' '
         DISPFL = TRIM(ADJUSTL(TMPDIR))//'surf_lov_'//TRIM(CIL)//'-'//TRIM(CID)//'.out'
         IF (CSIDE == 'L') THEN
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'srf_lov_phvelL-'//TRIM(CIL)//'-'// &
                      TRIM(CID)//'.out'
         ELSE
            PHASEFL = TRIM(ADJUSTL(TMPDIR))//'srf_lov_phvelR-'//TRIM(CIL)//'-'// &
                      TRIM(CID)//'.out'
         ENDIF
         DISPFL = ADJUSTL(DISPFL)
         PHASEFL = ADJUSTL(PHASEFL)
         PHASE_LOV(1:80) = ' '
         PHASE_LOV = TRIM(ADJUSTL(PHASEFL))
         PHASE_LOV = ADJUSTL(PHASE_LOV)
         LOVFL = TRIM(ADJUSTL(TMPDIR))//'lov-'//TRIM(CID)
         LOVFL = ADJUSTL(LOVFL)
         !IF (.NOT.LDEL) THEN 
         !   PHASEFL(1:80) = ' '
         !   PHASEFL = 'surf_lov_'//TRIM(CID)//'.out'
         !   PHASEFL = ADJUSTL(PHASEFL) 
         !ENDIF
         !IF (LNEW) THEN 
             INQUIRE(FILE=TRIM(PHASEFL),EXIST=LEX)
             IF (LEX) THEN 
                OPEN(UNIT=66,FILE=TRIM(PHASEFL),STATUS='OLD')
                CLOSE(66,STATUS='DELETE')
             ENDIF
         !ENDIF
         CALL EARTHSR(NL,JCOM, NSRCIN,NFREQIN,IEFL,NBRAN1,NBRAN2P,C1,C2,            &
                      TREF,ZKM,FREQ,SDEP, HDD_LV,VPD_LV,VSD_LV,ROD_LV,QA,QB,        &
                      LBIN,LDEL,DISPFL,PHASEFL,LOVFL,                               &
                      PERRAYW,PERLOVW,CCRAYW,CCLOVW, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'grns_srf: Error calling earthsr for Love waves'
            GOTO 500
         ENDIF
         IF (IL == 1) THEN
!           IF (CSIDE == 'L' .AND. ISIDE == 1 .OR. &
!               CSIDE == 'R' .AND. ISIDE == 2) THEN
            DO 16 IORD=0,NBRAN2P
               IORD1 = IORD + 1
               CCLOV(IORD1) = CCLOVW(IORD1)
   16       CONTINUE
!           ENDIF
         ENDIF 
!
!....... calculate offset by projecting points along model azimuth then
!....... calculating source to point via BJDAZ2
         ALLOCATE(XDIST(NXPTS))
         ALLOCATE(AZIMS(NXPTS))
         CALL XPTS2DIST3(NXPTS, CSIDE,SLAT,SLON, XMOD0,XMOD1,           &
                         XMLAT0,XMLON0, XMLAT1,XMLON1, XLAYR(1:NXPTS),  &
                         XDIST,AZIMS,IERR) 
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'grns_srf: Error calling xpts2dist3!'
            WRITE(*,*) 'grns_srf: Ben you have to program a backup!'
            RETURN
         ENDIF
!        IF (IERR /= 0) THEN !switch to linear interpolation
!           CALL XPTS2DIST(NXPTS,.FALSE.,POFF, SLAT,SLON, XMLAT0,XMLON0,   &
!                          XMLAT1,XMLON1, XMOD0,XMOD1, XLAYR(1:NXPTS),XDIST)
!        ENDIF
!
!....... now technically you would call srgramf at different modes and get the 
!....... greens functions for that mode, but that is for someone else to finish 
!        DO 20 IMODE=0,MODE 
            !IMODE1 = IMODE + 1
            !MDMAX = MODE 
            UGRN1D(1:NXPTS) = DCMPLX(0.D0,0.D0)
            VGRN1D(1:NXPTS) = DCMPLX(0.D0,0.D0)
            WGRN1D(1:NXPTS) = DCMPLX(0.D0,0.D0)
            IZINTER = IL 
            CALL SRGRAMF(NFREQIN,LPRATT,LBIN, IZINTER,NXPTS, IEVTIME(1:4),SSEC,  &
                         TMPDIR,RAYFL,LOVFL, MDMIN,MDMAX, .TRUE.,      &
                         SLAT,SLON,SDEP,                        &
                         STRIKE,DIP,RAKE,SMAG,T0,               &
                         ISRC,XMLAT0,XMLON0,XMLAT1,XMLON1,      &
                         MYID,XDIST,AZIMS,FREQ,                 &
                         NXPTS,UGRN1D(1:NXPTS),VGRN1D(1:NXPTS),WGRN1D(1:NXPTS), IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'grns_srf: Error calling srgramf'
               GOTO 500
            ENDIF
!  20    CONTINUE
!
!....... file clean up
         !IF (LNEW) THEN
            INQUIRE(FILE=TRIM(PHASE_RAY),EXIST=LEX) 
            IF (LEX) THEN 
               OPEN(UNIT=66,FILE=TRIM(PHASE_RAY),STATUS='OLD')
               CLOSE(66,STATUS='DELETE')
            ENDIF
            INQUIRE(FILE=TRIM(PHASE_LOV),EXIST=LEX) 
            IF (LEX) THEN 
               OPEN(UNIT=66,FILE=TRIM(PHASE_LOV),STATUS='OLD')
               CLOSE(66,STATUS='DELETE')
            ENDIF
         !ENDIF
!
!....... save the greens fns to points in bielak boundary
         DO 115 IXPTS=1,NXPTS 
!
!.......... modify greens functions
            IF (SRCTYP(2:2) == 'R') THEN !rayleigh only
               VGRN1D(IXPTS) = DCMPLX(0.D0,0.D0)
            ELSEIF (SRCTYP(2:2) == 'L') THEN !love only 
               UGRN1D(IXPTS) = DCMPLX(0.D0,0.D0)
               WGRN1D(IXPTS) = DCMPLX(0.D0,0.D0)
            ELSEIF (SRCTYP(2:2) == 'V') THEN !vertical only, not smart
               UGRN1D(IXPTS) = DCMPLX(0.D0,0.D0) !vertical only
               VGRN1D(IXPTS) = DCMPLX(0.D0,0.D0) !vertical only
            ELSE
               IF (SRCTYP(2:2) /= 'B') THEN
                  WRITE(*,*) 'grns_srf: Unknown source type!'
                  IERR = 1
                  GOTO 500
               ENDIF
            ENDIF
!
!.......... rotate into plane
            INPE = MYNNPE(IXPTS) 
            IF (LNPINI(INPE)) & 
            WRITE(*,*) 'grns_srf: Warning point already initialized!',CSIDE,INPE,IXPTS
            LNPINI(INPE) = .TRUE. 
            U = UGRN1D(IXPTS)*CCBAZ - VGRN1D(IXPTS)*CSBAZ
            V = UGRN1D(IXPTS)*CSBAZ + VGRN1D(IXPTS)*CCBAZ
            W = WGRN1D(IXPTS) !already point up 
            IF (LSHIFT) THEN
               !ccray is in km/s
               IF (CC_SAVE > 0.D0) THEN
                  !WRITE(*,*) 'grns_srf: Error zero phase velocity!'
!                 IERR = 1 
!                 RETURN 
                  !apparent slowness in s/km, yoff = 0.0 by default b/c we are in plane
                  PX = DCOS(AZIMS(IXPTS)*PI180)/CC_SAVE !CCRAYW(1) !slowness in x, note sin(90) = 1
                  PY = DSIN(AZIMS(IXPTS)*PI180)/CC_SAVE !CCRAYW(1) !slowness in y, note sin(90) = 1
                  XOFF = DABS(XDIST(IXPTS) - DELT0)      !km
                  ARG =-OMEGA*(XOFF*PX + YOFF*PY) !-omega p.x 
                  U = U*CDEXP(DCMPLX(0.D0,ARG))
                  V = V*CDEXP(DCMPLX(0.D0,ARG))
                  W = W*CDEXP(DCMPLX(0.D0,ARG))
               ELSE
                  print *, 'weird'
                  WRITE(*,*) 'grns_srf: Error zero phase velocity!'
                  IERR = 1
                  RETURN
               ENDIF
            ENDIF
!           IF (SRCTYP(2:2) == 'R') THEN
!              V = DCMPLX(0.D0,0.D0) 
!           ELSEIF (SRCTYP(2:2) == 'L') THEN
!              U = DCMPLX(0.D0,0.D0)
!              W = DCMPLX(0.D0,0.D0) 
!           ELSEIF (SRCTYP(2:2) == 'V') THEN
!              U = DCMPLX(0.D0,0.D0)
!              V = DCMPLX(0.D0,0.D0)
!           ENDIF
            UGRNS(INPE) = CMPLX(U*STF)
            VGRNS(INPE) = CMPLX(V*STF)
            WGRNS(INPE) = CMPLX(W*STF)
  115    CONTINUE !loop on x locations
         DEALLOCATE(XDIST)
         DEALLOCATE(AZIMS)
! 150    CONTINUE !break ahead, no points on this side
  105 CONTINUE !loop on depths
!
!.... clean for next pass
      DEALLOCATE(VPD_RL)
      DEALLOCATE(VSD_RL)
      DEALLOCATE(ROD_RL)
      DEALLOCATE(HDD_RL) 
      DEALLOCATE(VPD_LV)
      DEALLOCATE(VSD_LV)
      DEALLOCATE(ROD_LV)
      DEALLOCATE(HDD_LV) 
      DEALLOCATE(QA)
      DEALLOCATE(QB)
      DEALLOCATE(Z1D) 
!   1 CONTINUE !loop on sides
!
!.... map the bielak solution to the global vector
      IF (.NOT.LSAVE) THEN
         IF (LNSRF) THEN 
            DO 29 INPE=1,NNPE
               IF (CABS(UGRNS(INPE)) > 0.0) UGRNS(INPE) = UGRNS(INPE)/CABS(UGRNS(INPE))
               IF (CABS(VGRNS(INPE)) > 0.0) VGRNS(INPE) = VGRNS(INPE)/CABS(VGRNS(INPE))
               IF (CABS(WGRNS(INPE)) > 0.0) WGRNS(INPE) = WGRNS(INPE)/CABS(WGRNS(INPE))
  29        CONTINUE
         ENDIF
         CALL FILL_GRNS(MDIM,NDOF,NNPE, NDIM,LNPINI,IDOFSE, UGRNS,VGRNS,WGRNS, UE,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'grns_srf: Serious warning in fill_grns'
            IERR = 0
         ENDIF
      ELSE 
         DO 30 INPE=1,NNPE
            IF (.NOT.LNPINI(INPE)) & 
            WRITE(*,*) 'grns_srf: Serious warning in grns_bdy, unitialized point',INPE
            IF (LNSRF) THEN
               IF (CABS(UGRNS(INPE)) > 0.0) UGRNS(INPE) = UGRNS(INPE)/CABS(UGRNS(INPE))
               IF (CABS(VGRNS(INPE)) > 0.0) VGRNS(INPE) = VGRNS(INPE)/CABS(VGRNS(INPE))
               IF (CABS(WGRNS(INPE)) > 0.0) WGRNS(INPE) = WGRNS(INPE)/CABS(WGRNS(INPE))
            ENDIF
   30    CONTINUE 
         CALL SAVE_GRNS(NDIM, NNPE,NDIM,2,ISRC,.FALSE., FREQ, IDOFSE,UGRNS,VGRNS,WGRNS) 
      ENDIF
  500 CONTINUE !break ahead for errors
      DEALLOCATE(UGRN1D)
      DEALLOCATE(VGRN1D) 
      DEALLOCATE(WGRN1D)  
      DEALLOCATE(UGRNS)
      DEALLOCATE(VGRNS)
      DEALLOCATE(WGRNS)
      DEALLOCATE(CLP) 
      DEALLOCATE(LNPINI) 
      DEALLOCATE(MYNNPE)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE UPD_PYTAB(PHSFL, MODE, FREQ,BAZN, PYTAB, IERR) 
!
!     Updates the py table.  Based on SWR's initiosr.f  - B Baker Jan 2013 
!
!     INPUT      MEANING
!     -----      ------- 
!     BAZN       back azimuth
!     FREQ       current frequency
!     MODE       desired mode to model 
!     PHSFL      holds phase velocities at period for modes
! 
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     PYTAB      py table for this frequency and source
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: PHSFL 
      REAL*8, INTENT(IN) :: BAZN,FREQ
      INTEGER*4, INTENT(IN) :: MODE 
      REAL*8, INTENT(OUT) :: PYTAB
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: PERIOD(:), CPHASE(:) 
      REAL*8 PERIN,CCIN,X1,X2,Y1,Y2,X,A,B,PHCZ,PI180
      INTEGER*4 J,NPERIOD,MODEIN,IUNIT  
      PARAMETER(IUNIT = 65) 
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
!
!----------------------------------------------------------------------------------------!
!
!.... open file 
      IERR = 0
      OPEN(UNIT=IUNIT,FILE=TRIM(PHSFL),STATUS='OLD',IOSTAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'upd_pytab: There was an error opening:',TRIM(PHSFL)
         RETURN
      ENDIF
!
!.... count the number of periods for this mode
      J = 0 
    1 READ(IUNIT,*,END=60) 
         READ(IUNIT,*,END=60) MODEIN,PERIN,CCIN
         IF (MODEIN == MODE) J = J + 1
      GOTO 1
   60 REWIND(IUNIT) 
      NPERIOD = J
!
!.... now read the modes
      J = 0 
      ALLOCATE(PERIOD(NPERIOD))
      ALLOCATE(CPHASE(NPERIOD)) 
    2 READ(IUNIT,*,END=65) 
         READ(IUNIT,*,END=65) MODEIN,PERIN,CCIN
         IF (MODEIN == MODE) THEN
            J = J + 1
            PERIOD(J) = PERIN !period
            CPHASE(J) = CCIN  !phase velocity km/s
         ENDIF
      GOTO 2 
   65 CONTINUE 
      CLOSE(IUNIT) 
!
!.... set limits to the period 
      PHCZ = 0.D0
      IF (1.D0/FREQ <= PERIOD(1)) THEN
         PHCZ = CPHASE(1)
      ELSEIF (1.D0/FREQ >= PERIOD(NPERIOD)) THEN
         PHCZ = CPHASE(NPERIOD)
      ELSE !linearly interoplate via two point
         DO 3 J=1,NPERIOD-1
            IF (1.D0/FREQ >= PERIOD(J) .AND. 1.D0/FREQ <= PERIOD(J+1)) THEN 
               Y2 = CPHASE(J+1)
               Y1 = CPHASE(J) 
               X2 = PERIOD(J+1)
               X1 = PERIOD(J) 
               X = 1.D0/FREQ - X1 
               A = (Y2 - Y1)/(X2 - X1)
               B = Y1  
               PHCZ = A*X + B 
            ENDIF
    3    CONTINUE 
      ENDIF 
      IF (PHCZ == 0.D0) THEN
         WRITE(*,*) 'upd_pytab: Error phcz = 0.'
         IERR = 1
         RETURN
      ENDIF
      PHCZ = PHCZ*1000.D0 !convert to m/s
!     there is no negative sign b/c this is a forward transform w/ e{-iwt}
!     of course for modeling we will take the absolute value, we just want to correct
!     the stations
      PYTAB = SIN(BAZN*PI180)/PHCZ
      DEALLOCATE(PERIOD)
      DEALLOCATE(CPHASE)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE CLASS_NNPE(NNPE,CSIDE,TOL,XBLKL,XBLKR, XLOCSE, CLP,IERR)
!
!     Classifies the Bielak locations as either "Left" or "Right" for this source.
!     For example if a wave is coming from the "Left" side then everything up to the 
!     right boundary is in the "Left".  After that boundary it is on the "Right"
!
!     INPUT      MEANING
!     -----      ------- 
!     CSIDE      L -> wave approaching from left, R -> wave approaching from right
!     NNPE       number of nodal points in Bielak boundary
!     XBLKL      x bielak/internal left boundary 
!     XBLKR      x bielak/internal right boundary
!     XLOCSE     x locations of points in Bielak layer 
!     
!     OUTPUT     MEANING
!     ------     ------- 
!     CLP        determins whether a node is on the left or right
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CSIDE
      REAL*8, INTENT(IN) :: XLOCSE(NNPE), TOL,XBLKL,XBLKR 
      INTEGER*4, INTENT(IN) :: NNPE
      CHARACTER(1), INTENT(OUT) :: CLP(NNPE) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local varaibles
      INTEGER*4 INPE
!
!----------------------------------------------------------------------------------------!
!
!.... loop on points
      IERR = 0
      CLP(1:NNPE) = ' '
      DO 1 INPE=1,NNPE
         IF (CSIDE == 'L') THEN
            IF (XLOCSE(INPE) >= XBLKR - TOL) THEN
               CLP(INPE) = 'R'
            ELSE
               CLP(INPE) = 'L'
            ENDIF
         ELSE
            IF (XLOCSE(INPE) <= XBLKL + TOL) THEN
               CLP(INPE) = 'L'
            ELSE
               CLP(INPE) = 'R'
            ENDIF
         ENDIF
    1 CONTINUE
!
!.... check my work
      DO 2 INPE=1,NNPE
         IF (CLP(INPE) == ' ') THEN
            WRITE(*,*) 'class_nnpe: Error locating point!',INPE
            IERR = 1
         ENDIF 
    2 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE CLASS_NNPE_BDY(NNPE,CSIDE,TOL,XBLKL,XBLKR,ZBASE, XLOCSE,ZLOCSE,  &
                                CLP,IERR)
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CSIDE
      REAL*8, INTENT(IN) :: XLOCSE(NNPE),ZLOCSE(NNPE), TOL,XBLKL,XBLKR, ZBASE
      INTEGER*4, INTENT(IN) :: NNPE 
      CHARACTER(1), INTENT(OUT) :: CLP(NNPE)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 INPE 
!
!----------------------------------------------------------------------------------------!
!
!.... classes the nodal points for body waves
      IERR = 0
      CLP(1:NNPE) = ' '
      DO 1 INPE=1,NNPE
         IF (CSIDE == 'L') THEN !wave coming from left
            IF (XLOCSE(INPE) >= XBLKR - TOL) THEN
               CLP(INPE) = 'O'
            ELSE !side of model
               CLP(INPE) = 'S'
            ENDIF
         ELSE !come from right
            IF (XLOCSE(INPE) >= XBLKL + TOL) THEN
               CLP(INPE) = 'S'
            ELSE
               CLP(INPE) = 'O'
            ENDIF
         ENDIF
         IF (ZLOCSE(INPE) <= ZBASE + TOL) CLP(INPE) = 'B'
    1 CONTINUE
!
!.... fidelity check
      DO 2 INPE=1,NNPE
         IF (CLP(INPE) == ' ') THEN
            WRITE(*,*) 'class_nnpe_bdy: Error locating point!',INPE
            IERR = 1
         ENDIF
    2 CONTINUE
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE CLASS_NNPE_SRF(NNPE,CSIDE,TOL, XBLKL,XBLKR, XLOCSE, CLP,IERR) 
      CHARACTER(1), INTENT(IN) :: CSIDE 
      REAL*8, INTENT(IN) :: XLOCSE(NNPE), XBLKL, XBLKR 
      INTEGER*4, INTENT(IN) :: NNPE 
      CHARACTER(1), INTENT(OUT) :: CLP(NNPE) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4 INPE 
!
!----------------------------------------------------------------------------------------!
!
      IERR = 0
      CLP(1:NNPE) = ' '
      DO 1 INPE=1,NNPE
         IF (CSIDE == 'L') THEN !wave coming from left
            IF (XLOCSE(INPE) <= XBLKL + TOL) THEN
               CLP(INPE) = 'L'
            ELSE
               IF (XLOCSE(INPE) <= XBLKR - TOL) THEN
                  CLP(INPE) = 'B'
               ELSE
                  CLP(INPE) = 'R'
               ENDIF
            ENDIF
         ELSE !wave coming from right
            IF (XLOCSE(INPE) >= XBLKR - TOL) THEN
               CLP(INPE) = 'R'
            ELSE
               IF (XLOCSE(INPE) >= XBLKL + TOL) THEN
                  CLP(INPE) = 'B'
               ELSE
                  CLP(INPE) = 'L'
               ENDIF
            ENDIF
        ENDIF !end check on direction
    1 CONTINUE
!
!.... double check we got everything
      DO 2 INPE=1,NNPE
         IF (CLP(INPE) == ' ') THEN
            WRITE(*,*) 'class_nnpe_srf: Error locating point!',INPE
            IERR = 1
         ENDIF
    2 CONTINUE
      RETURN
      END
 
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE SET_GRNS1D(NNPE,ISIDE,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D) 
!
!     Sets all the 1D points to expand for this side 
!
!     INPUT      MEANING
!     -----      ------- 
!     CLP        holds if z location is on left or right side
!     ISIDE      = 1 -> want left side, = 2 -> want right side
!     NNPE       number of nodal points in Bielak layer
!     TOL        tolerance for when i call a point the same point
!     ZLOCSE     z locations of Bielak nodal points 
! 
!     OUTPUT     MEANING
!     ------     -------
!     ZPTS1D     1D points to set 
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CLP(NNPE) 
      REAL*8, INTENT(IN) :: ZLOCSE(NNPE), TOL  
      INTEGER*4, INTENT(IN) :: NNPE, ISIDE 
      REAL*8, INTENT(OUT) :: ZPTS1D(NNPE) 
      INTEGER*4, INTENT(OUT) :: NNP1D 
!.... local variables
      REAL*8 ZMAX 
      INTEGER*4 INP1D, JNP1D, INPE 
      LOGICAL*4 LCHECK 
!
!----------------------------------------------------------------------------------------!
!
      ZMAX =-HUGE(1.D0)
      ZPTS1D(1:NNPE) = 0.D0
      INP1D = 0
      DO 1 INPE=1,NNPE
         LCHECK = .FALSE.
         IF (ISIDE == 1 .AND. CLP(INPE) == 'L') LCHECK = .TRUE.
         IF (ISIDE == 2 .AND. CLP(INPE) == 'R') LCHECK = .TRUE.
         IF (LCHECK) THEN
            IF (INPE < NNPE) THEN
               IF (ZLOCSE(INPE) /= ZLOCSE(INPE+1) .OR.  &
                   CLP(INPE)    /= CLP   (INPE+1)) THEN
                  DO 2 JNP1D=1,INP1D !no repeats
                     IF (DABS(ZLOCSE(INPE) - ZPTS1D(INP1D)) < TOL) GOTO 10 
    2             CONTINUE
                  INP1D = INP1D + 1
                  ZPTS1D(INP1D) = ZLOCSE(INPE)
               ENDIF
            ELSE
               IF (ZLOCSE(NNPE) /= ZLOCSE(NNPE-1) .OR.  &
                   CLP(NNPE)    /= CLP   (NNPE-1)) THEN
                  DO 3 JNP1D=1,INP1D !no repeats
                     IF (DABS(ZLOCSE(NNPE) - ZPTS1D(INP1D)) < TOL) GOTO 10
    3             CONTINUE
                  INP1D = INP1D + 1
                  ZPTS1D(INP1D) = ZLOCSE(NNPE)
               ENDIF
            ENDIF
            ZMAX = MAX(ZMAX,ZLOCSE(INPE))
         ENDIF
   10    CONTINUE !break ahead for repeat
    1 CONTINUE
!
!.... i can miss the last point
      IF (DABS(ZPTS1D(INP1D) - ZMAX) > TOL .AND. ZMAX /= -HUGE(1.D0)) THEN
         INP1D = INP1D + 1
         ZPTS1D(INP1D) = ZMAX
      ENDIF
      NNP1D = INP1D
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SET_GRNS1D_BDY(NNPE,JOB,CLP,TOL,ZLOCSE, NNP1D,ZPTS1D)
      IMPLICIT NONE 
      CHARACTER(1), INTENT(IN) :: CLP(NNPE) 
      REAL*8, INTENT(IN) :: ZLOCSE(NNPE), TOL  
      INTEGER*4, INTENT(IN) :: NNPE, JOB 
      REAL*8, INTENT(OUT) :: ZPTS1D(NNPE) 
      INTEGER*4, INTENT(OUT) :: NNP1D 
!.... local variables
      REAL*8 ZMAX
      INTEGER*4 INP1D, JNP1D, INPE
      LOGICAL*4 LCHECK
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      ZMAX =-HUGE(1.D0)
      ZPTS1D(1:NNPE) = 0.D0
!
!.... extract unique z locations at model base or model side
      INP1D = 0
      DO 1 INPE=1,NNPE
         LCHECK = .FALSE.
         IF (JOB == 1 .AND. CLP(INPE) == 'S') LCHECK = .TRUE.
         IF (JOB == 2 .AND. CLP(INPE) == 'B') LCHECK = .TRUE. 
         IF (LCHECK) THEN
            IF (INPE < NNPE) THEN
               IF (ZLOCSE(INPE) /= ZLOCSE(INPE+1) .OR. &
                   CLP(INPE)    /= CLP   (INPE+1)) THEN
                  DO 2 JNP1D=1,INP1D !no repeats
                     IF (DABS(ZLOCSE(INPE) - ZPTS1D(INP1D)) < TOL) GOTO 10
    2             CONTINUE
                  INP1D = INP1D + 1
                  ZPTS1D(INP1D) = ZLOCSE(INPE)
               ENDIF
            ELSE
               IF (ZLOCSE(NNPE) /= ZLOCSE(NNPE-1) .OR. &
                   CLP   (NNPE) /= CLP   (NNPE-1)) THEN
                  DO 3 JNP1D=1,INP1D !No repeats
                     IF (DABS(ZLOCSE(NNPE) - ZPTS1D(INP1D)) < TOL) GOTO 10
    3             CONTINUE
                  INP1D = INP1D + 1
                  ZPTS1D(INP1D) = ZLOCSE(NNPE)
              ENDIF
            ENDIF
            ZMAX = MAX(ZMAX,ZLOCSE(INPE))
         ENDIF
   10    CONTINUE !break ahead for repeat
    1 CONTINUE ! loop
!
!.... i can miss the last point
      IF (DABS(ZPTS1D(INP1D) - ZMAX) > TOL .AND. ZMAX /= -HUGE(1.D0)) THEN
         INP1D = INP1D + 1
         ZPTS1D(INP1D) = ZMAX
      ENDIF
      NNP1D = INP1D
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FILL_GRNS(MDIM,NDOF,NNPE, NDIM,LNPINI,IDOFSE, UGRNS,VGRNS,WGRNS, UE,IERR)
!
!     Fills in the Greens functions on a global vector
!
!     INPUT      MEANING
!     -----      ------- 
!     IDOFSE     maps Bielak nodes to global DOF numbers
!     LNPINI     determines if point was set 
!     MDIM       leading dimension for IDOFSE
!     NDIM       number of components in solution
!     NDOF       number of degrees of freedom in solution
!     NNPE       number of nodal poits in bielak layer
!     UGRNS      1D Greens functions in u component
!     VGRNS      1D Greens functions in v component 
!     WGRNS      1D Greens functions in w component 
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     UE         analytic solution in Bielak boundary 
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: UGRNS(NNPE), VGRNS(*), WGRNS(NNPE) 
      INTEGER*4, INTENT(IN) :: IDOFSE(MDIM,*), MDIM, NNPE, NDIM, NDOF 
      LOGICAL*4, INTENT(IN) :: LNPINI(NNPE) 
      COMPLEX*8, INTENT(OUT) :: UE(NDOF) 
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variable
      COMPLEX*8 CZERO
      INTEGER*4 INPE, I, IDOF 
      PARAMETER(CZERO = CMPLX(0.0,0.0)) 
!
!----------------------------------------------------------------------------------------!
!
!.... null out ue and fill in points
      IERR = 0
      UE(1:NDOF) = CZERO 
      DO 1 INPE=1,NNPE
         IF (.NOT.LNPINI(INPE)) THEN
            WRITE(*,*) 'fill_grns: Error uninitialized nodal point:',INPE
            IERR = 1
         ENDIF
         DO 2 I=1,NDIM
            IDOF = IDOFSE(I,INPE)
            IF (IDOF > 0) THEN
               IF (NDIM == 3) THEN
                  IF (I == 1) THEN
                     UE(IDOF) = UGRNS(INPE)
                  ELSEIF (I == 2) THEN
                     UE(IDOF) = VGRNS(INPE)
                  ELSEIF (I == 3) THEN
                     UE(IDOF) = WGRNS(INPE)
                  ENDIF
               ELSE
                  IF (I == 1) THEN
                     UE(IDOF) = UGRNS(INPE)
                  ELSE
                     UE(IDOF) = WGRNS(INPE)
                  ENDIF
               ENDIF 
            ENDIF !end check on DOF
    2    CONTINUE !loop on components
    1 CONTINUE !loop on nodal points
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE SAVE_GRNS(MDIM, NNPE,NDIM,JOB,ISRC,LDOFS, FREQ, IDOFSE,UGRNS,VGRNS,WGRNS)
!
!     Saves the Greens functions calculated at global DOFs in bielak boundary
!  
!     INPUT      MEANING
!     -----      -------
!     FREQ       current frequency
!     IDOFSE     maps bielak nodes to global DOFs 
!     JOB        1 -> Body waves; 2 -> Surface waves
!     LDOFS      True -> write the IDOFSE map
!     FREQ       frequency (Hz)
!     MDIM       leading dimension for IDOFSE
!     NDIM       number of components in solution
!     NNPE       number of nodal points in bielak layer
!     UGRNS      Greens functions on component 1
!     VGRNS      Greens functions on component 2
!     WGRNS      Greens functions on component 3 
! 
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: UGRNS(NNPE), VGRNS(*), WGRNS(NNPE)
      REAL*8, INTENT(IN) :: FREQ
      INTEGER*4, INTENT(IN) :: IDOFSE(MDIM,*),MDIM,NNPE,NDIM,JOB,ISRC
      LOGICAL*4, INTENT(IN) :: LDOFS 
!.... local variables
      CHARACTER(80) FILENM
      CHARACTER(12) CFREQ 
      CHARACTER(5) CSRC
      INTEGER*4 IUNIT, I, INPE
      PARAMETER(IUNIT = 66)
!
!----------------------------------------------------------------------------------------!
!
!.... set file name and open file
      FILENM(1:80) = ' '
      CFREQ(1:12) = ' '
      CSRC(1:5) = ' '
      WRITE(CFREQ,'(F12.5)') FREQ
      WRITE(CSRC ,'(I5)') ISRC
      CFREQ = ADJUSTL(CFREQ)
      CSRC = ADJUSTL(CSRC) 
      IF (JOB == 1) THEN !body wave
         FILENM = './grns/'//'body-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
      ELSE !surface wave
         FILENM = './grns/'//'surf-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
      ENDIF
      FILENM = ADJUSTL(FILENM) 
      OPEN(UNIT=IUNIT,FILE=TRIM(FILENM),STATUS='NEW',FORM='UNFORMATTED')
!.... write file
      WRITE(IUNIT) NDIM
      WRITE(IUNIT) NNPE
      WRITE(IUNIT) ISRC
      WRITE(IUNIT) FREQ
      WRITE(IUNIT) (UGRNS(INPE),INPE=1,NNPE)
      IF (NDIM > 2) WRITE(IUNIT) (VGRNS(INPE),INPE=1,NNPE)
      WRITE(IUNIT) (WGRNS(INPE),INPE=1,NNPE)
      IF (LDOFS) WRITE(IUNIT) ((IDOFSE(I,INPE),I=1,NDIM),INPE=1,NNPE) 
      CLOSE(IUNIT)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_GRNS(MDIM, NNPE,NDIM,JOB,ISRC,LDOFS, FREQ, IDOFSE, & 
                           UGRNS,VGRNS,WGRNS,IERR)
!
!     Loads the Greens functions calculated at global DOFs in bielak boundary
!  
!     INPUT      MEANING
!     -----      -------
!     FREQ       current frequency
!     JOB        1 -> Body waves; 2 -> Surface waves
!     LDOFS      True -> write the IDOFSE map
!     FREQ       frequency (Hz)
!     MDIM       leading dimension for IDOFSE
!     NDIM       number of components in solution
!     NNPE       number of nodal points in bielak layer
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     IDOFSE     LDOFS True -> read IDOFSE map and output
!     UGRNS      Greens functions on component 1
!     VGRNS      Greens functions on component 2
!     WGRNS      Greens functions on component 3 
! 
!.... variable declarations
      IMPLICIT NONE 
      REAL*8, INTENT(IN) :: FREQ 
      INTEGER*4, INTENT(IN) :: MDIM,NNPE,NDIM,JOB,ISRC
      LOGICAL*4, INTENT(IN) :: LDOFS 
      INTEGER*4 IDOFSE(MDIM,*) 
      COMPLEX*8, INTENT(OUT) :: UGRNS(NNPE), VGRNS(*), WGRNS(NNPE)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      CHARACTER(80) FILENM
      CHARACTER(12) CFREQ 
      CHARACTER(5) CSRC 
      REAL*8 FREQ_IN
      INTEGER*4 IUNIT, NDIM_IN, NNPE_IN, ISRC_IN, I, INPE 
      PARAMETER(IUNIT = 66)
!
!----------------------------------------------------------------------------------------!
!
!.... set file name and open file
      FILENM(1:80) = ' '
      CFREQ(1:12) = ' '
      CSRC(1:5) = ' '
      WRITE(CFREQ,'(F12.5)') FREQ 
      WRITE(CSRC ,'(I5)') ISRC 
      CFREQ = ADJUSTL(CFREQ)
      CSRC = ADJUSTL(CSRC) 
      IF (JOB == 1) THEN !body wave
         FILENM = './grns/'//'body-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
      ELSE !surface wave
         FILENM = './grns/'//'surf-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
      ENDIF
      FILENM = ADJUSTL(FILENM) 
      OPEN(UNIT=IUNIT,FILE=TRIM(FILENM),STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'load_grns: Error reading file:',TRIM(FILENM)
         RETURN
      ENDIF
!.... read 
      READ(IUNIT,END=60,IOSTAT=IERR) NDIM_IN
      READ(IUNIT,END=60,IOSTAT=IERR) NNPE_IN
      READ(IUNIT,END=60,IOSTAT=IERR) ISRC_IN
      READ(IUNIT,END=60,IOSTAT=IERR) FREQ_IN 
!
!.... check headers
      IF (ISRC_IN /= ISRC) WRITE(*,*) 'load_grns: Warning isrc_in /= isrc'
      IF (ABS(FREQ_IN - FREQ) > 1.1E-7) WRITE(*,*) 'load_grns: Warning freq_in /= freq'
      IF (NDIM_IN /= NDIM) WRITE(*,*) 'load_grns: Serious warning ndim_in /= ndim'
      IF (NNPE_IN /= NNPE) WRITE(*,*) 'load_grns: Serious warning nnpe_in /= nnpe' 
      READ(IUNIT,END=60,IOSTAT=IERR) (UGRNS(INPE),INPE=1,NNPE)
      IF (NDIM > 2) READ(IUNIT,END=60,IOSTAT=IERR) (VGRNS(INPE),INPE=1,NNPE)
      READ(IUNIT) (WGRNS(INPE),INPE=1,NNPE)
      IF (LDOFS) READ(IUNIT,END=60,IOSTAT=IERR) ((IDOFSE(I,INPE),I=1,NDIM),INPE=1,NNPE)
      CLOSE(IUNIT)
      RETURN
   60 CONTINUE
      WRITE(*,*) 'load_grns: Premature end of Greens function file!'
      IERR = 1
      CLOSE(IUNIT) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE WRITE_PYTAB(MFREQ,NFREQ,NSRC, VFAST,VFB,VFS, FREQ,PYTAB, IERR)
!
!     Writes the py table
!
!     INPUT      MEANING
!     -----      -------
!     FREQ       frequency list
!     MFREQ      leading dimension for pytable 
!     NFREQ      number of frequencies
!     NSRC       number of sources
!     PYTAB      py table 
!     VFAST      fast velocity for calculating source groups
!     VFB        fastest body wave velocity 
!     VFS        fastest surface wave velocity
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       =1 indicates an error occurred when writing the py table
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: PYTAB(MFREQ,*), FREQ(NFREQ), VFAST, VFB, VFS 
      INTEGER*4, INTENT(IN) :: MFREQ, NFREQ, NSRC
      INTEGER*4, INTENT(OUT) :: IERR 
!.... variable declarations
      CHARACTER(4), ALLOCATABLE :: CDAT(:)
      CHARACTER(4) PACKR4
      CHARACTER(4) PACKI4 
      INTEGER*4 NBYTES, LWORK, INDX, IFREQ, ISRC, J, MYEND, ENDIAN, IUNIT
      LOGICAL*4 LSWAP, LEX, LISDIR 
      PARAMETER(IUNIT = 65) 
!
!----------------------------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0 
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND /= 0) LSWAP = .TRUE. 
!.... allocate space and pack header 
      NBYTES = 20 + 4*NFREQ + 4*NFREQ*NSRC  
      LWORK = NBYTES/4 
      ALLOCATE(CDAT(LWORK)) 
      CDAT(1) = PACKI4(LSWAP,NFREQ) 
      CDAT(2) = PACKI4(LSWAP,NSRC) 
      CDAT(3) = PACKR4(LSWAP,SNGL(VFAST))
      CDAT(4) = PACKR4(LSWAP,SNGL(VFB))
      CDAT(5) = PACKR4(LSWAP,SNGL(VFS))  
      INDX = 6
!.... pack data
      DO 1 IFREQ=1,NFREQ  
         CDAT(INDX) = PACKR4(LSWAP,SNGL(FREQ(IFREQ)))
         INDX = INDX + 1
         DO 2 ISRC=1,NSRC
            CDAT(INDX) = PACKR4(LSWAP,SNGL(PYTAB(IFREQ,ISRC)))
            INDX = INDX + 1
    2    CONTINUE 
    1 CONTINUE
      LEX = LISDIR('./grns')
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./grns')
      OPEN(UNIT=IUNIT,FILE=TRIM('./grns/pytab.bin'),FORM='UNFORMATTED',  &
           STATUS='REPLACE',ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR) 
      IF (IERR /= 0) THEN
         WRITE(*,*) 'write_pytab: Error opening py table!'
         GOTO 500
      ENDIF
      WRITE(IUNIT,REC=1,IOSTAT=IERR) (CDAT(J),J=1,LWORK)
      IF (IERR /= 0) WRITE(*,*) 'write_pytab: Error writing py table!'
      CLOSE(IUNIT)
  500 CONTINUE !break ahead for error
      DEALLOCATE(CDAT)
      CLOSE(IUNIT)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_PYTAB(MFREQ,NFREQ,NSRC, FREQ, VFAST,VFB,VFS,PYTAB,IERR) 
!
!     Reads the py table 
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       frequency list
!     MFREQ      leading dimension for py table
!     NFREQ      number of frequencies to read
!     NSRC       number of sources
!    
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     PYTAB      py table
!     VFAST      max velocity (m/s)
!     VFB        fastest body wave velocity (m/s)
!     VFS        fastest surface wave velocity (m/s)
!
!.... variable declarations 
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FREQ(NFREQ)
      INTEGER*4, INTENT(IN) :: MFREQ, NFREQ, NSRC
      REAL*8, INTENT(OUT) :: PYTAB(MFREQ,NSRC), VFAST, VFB, VFS
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: CDAT(:) 
      REAL*8 FREQ_IN
      REAL*4 R4, UNPACKR4
      INTEGER*4 NBYTES, NADV, IUNIT, ISRC, IFREQ, IFREQ_IN, NSRC_IN, NFREQ_IN, INDX, I, &
                MYEND, UNPACKI4, ENDIAN
      LOGICAL*4 LSWAP, LEX 
      PARAMETER(IUNIT = 65)
!
!----------------------------------------------------------------------------------------!
!
!.... detect endianness
      IERR = 0
      LSWAP = .FALSE.
      MYEND = ENDIAN()
      IF (MYEND /= 0) LSWAP = .TRUE.
!.... get file size and set space
      INQUIRE(FILE='./grns/pytab.bin',EXIST=LEX,SIZE=NBYTES)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'read_pytab: py table binary file does not exist!'
         IERR = 1
         RETURN
      ENDIF
      !CALL STAT(TRIM('./grns/pytab.bin'),INFO)
      !NBYTES = INFO(8) 
      OPEN(UNIT=IUNIT,FILE=TRIM('./grns/pytab.bin'),FORM='UNFORMATTED',  &
           ACCESS='DIRECT',RECL=NBYTES,IOSTAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'read_pytab: Error opening py table!'
         RETURN
      ENDIF
      ALLOCATE(CDAT(NBYTES))
      READ(IUNIT,REC=1,IOSTAT=IERR) (CDAT(I),I=1,NBYTES)
      CLOSE(IUNIT)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'read_pytab: Error reading py table!'
         GOTO 500
      ENDIF 
!.... header
      NFREQ_IN = UNPACKI4(LSWAP,CDAT(1:4))
      NSRC_IN  = UNPACKI4(LSWAP,CDAT(5:8))
      R4       = UNPACKR4(LSWAP,CDAT(9:12))
      VFAST = DBLE(R4) 
      R4       = UNPACKR4(LSWAP,CDAT(13:16))
      VFB = DBLE(R4)
      R4       = UNPACKR4(LSWAP,CDAT(17:20))
      VFS = DBLE(R4)
      IF (NSRC_IN /= NSRC) THEN
         WRITE(*,*) 'read_pytab: Sources do not match!'
         IERR = 1
         GOTO 500
      ENDIF
!.... scan file for frequencies we need
      NADV = 4*NSRC_IN + 4 
      DO 1 IFREQ=1,NFREQ 
         INDX = 21 !offset after header
         DO 2 IFREQ_IN=1,NFREQ_IN
            R4 = UNPACKR4(LSWAP,CDAT(INDX:INDX+3))
            FREQ_IN = DBLE(R4) 
            IF (ABS(FREQ_IN - FREQ(IFREQ)) < 1.D-4) GOTO 20 
            INDX = INDX + NADV 
    2    CONTINUE  
         WRITE(*,*) 'read_pytab: Error could not locate frequency:',FREQ(IFREQ)
         IERR = 1
         GOTO 500
   20    CONTINUE
         DO 3 ISRC=1,NSRC
            INDX = INDX + 4
            R4 = UNPACKR4(LSWAP,CDAT(INDX:INDX+3))
            PYTAB(IFREQ,ISRC) = DBLE(R4) 
    3    CONTINUE
    1 CONTINUE 
  500 CONTINUE 
      DEALLOCATE(CDAT) 
      RETURN
      END
