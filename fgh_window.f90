      SUBROUTINE FGH_WINDOW(MASTER,MYID,MYNID,IPGROUP,NPGROUPS, NPPGRP,           &
                            MGCOMM,MYSLV_COMM,MYHD_COMM,                          &
                            LUPDREC,LUPDSRC,LSRCEX,LRECEX, LFUNC,LGRAD,      & 
                            IBLOCK,K,IALPHA,                                      &
                            MID,INV,MSH,SRC,RCV,FRQ,M1D,WIN, IERR) 
!
!     This is an expert driver routine for generating the objective function and 
!     Jacobian for inversion.  We have a few jobs to perform: 
!
!     If we are windowing data we have two frequency lists, a frequency modeling
!     list and as a subset an inversion list.  So part of the routine is modeling 
!     and, when a modeling frequency corresponds to an inversion frequency, we also 
!     calculate the Jacobian
!
!     After the modeled waveforms and Jacobians are claculated we then Fourier Transform
!     to the time domain, window, and inverse Fourier transform back to the frequency 
!     domain.
!
!     Next, if required, we update the STFs and RRFs.  
!
!     Finally, we multiply adj(J) delta d and calculate the objective function.
!
!     - B. Baker May 2013 
!
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      INCLUDE 'cmumps_struc.h'
      TYPE (CMUMPS_STRUC) MID 
      TYPE (INV_INFO)     INV 
      TYPE (MESH_INFO)    MSH
      TYPE (SRC_INFO)     SRC
      TYPE (RECV_INFO)    RCV
      TYPE (FRQ_INFO)     FRQ 
      TYPE (MOD1D_INFO)   M1D  
      TYPE (WIN_INFO)     WIN 
      INTEGER*4, INTENT(IN) :: MASTER,MYID,MYNID, MGCOMM,MYSLV_COMM,MYHD_COMM, &
                               NPGROUPS,NPPGRP,IPGROUP, IBLOCK,K,IALPHA 
      LOGICAL*4 LUPDREC,LUPDSRC,LSRCEX,LRECEX,LFUNC,LGRAD
      INTEGER*4, INTENT(OUT) :: IERR

!.... local variables
      COMPLEX*8, ALLOCATABLE :: UE(:), WAVE(:), FMAT_DIST(:) 
      REAL*4, ALLOCATABLE :: GRADL_SRF(:), GRADL_BDY(:), WORK(:)
      INTEGER*4, ALLOCATABLE :: NSG_LSTB(:), NSG_LST(:), ISRCLST(:) 
      LOGICAL*4, ALLOCATABLE :: LOBS(:,:) 
      COMPLEX*8 CZERO
      REAL*8 PY, PYAVG
      REAL*4 HESSB(1), EPSS, THRESHS, OBJ 
      INTEGER*4 STAT(MPI_STATUS_SIZE), NOBSLST, NOBS, NFREQ_OBS, IBPHASE,  &
                IFREQL,IFREQ,JFREQ, ISRC,JSRC,KSRC,  I,IREC, &
                ISG, JPGROUP, NSGMAX, NSG_GRP, MYSRC, MYDEST, MYTAG, MPIERR
      LOGICAL*4 LSTF, LPGNEW, LEXIST, LISDIR, LSURF, LBODY, LCJAC, LOAD_JAC, &
                LJACOB_SRF, LJACOB_BDY, LCOVD
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      PARAMETER(EPSS = 0.2, THRESHS = 0.2)
      PARAMETER(LPGNEW = .FALSE.) !we dont need a gradient preconditioner
      PARAMETER(IBPHASE = 0) !no geometric spreading correction
      character(80) filenm
!
!----------------------------------------------------------------------------------------!
!
!.... initial checks 
      IF (.NOT.inv%LSURF .AND. .NOT.inv%LBODY) THEN
         WRITE(*,*) 'fgh_window: Error no suface or body waves to propagate!'
         IERR = 1
      ENDIF
!
!.... send logicals
      CALL MPI_BCAST(   LGRAD,  1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(   LFUNC,  1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(   LRECEX, 1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(   LSRCEX, 1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(   LUPDSRC,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(   LUPDREC,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(inv%LGNEWT,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
!
!.... do we need a jacobian?
      IF (MYID == MASTER) THEN
         LJACOB_SRF = .FALSE.
         LJACOB_BDY = .FALSE.
         IF (LUPDSRC .OR. inv%LGNEWT .OR. win%LWNDO_SRF) LJACOB_SRF = .TRUE.
         IF (LUPDSRC .OR. inv%LGNEWT .OR. win%LWNDO_BDY) LJACOB_BDY = .TRUE.
         IF (.NOT.inv%LSURF) LJACOB_SRF = .FALSE. !no surface wave inversion
         IF (.NOT.inv%LBODY) LJACOB_BDY = .FALSE. !no body wave inversion 
      ENDIF
      CALL MPI_BCAST(LJACOB_SRF,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      CALL MPI_BCAST(LJACOB_BDY,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
      IF (.NOT.LJACOB_SRF) THEN
         IF (MYNID == MASTER) THEN
            ALLOCATE(GRADL_SRF(inv%NA35))
         ELSE
            ALLOCATE(GRADL_SRF(1))
         ENDIF
         GRADL_SRF(:) = 0.0
      ENDIF
      IF (.NOT.LJACOB_BDY) THEN
         IF (MYNID == MASTER) THEN
            ALLOCATE(GRADL_BDY(inv%NA35))
         ELSE
            ALLOCATE(GRADL_BDY(1))
         ENDIF
         GRADL_BDY(:) = 0.0
      ENDIF
 
!
!.... set space on wavefields and observation list
      ALLOCATE(UE(msh%NDOF))
      ALLOCATE(WAVE(msh%NDOF))
      UE(1:msh%NDOF) = CZERO
      WAVE(1:msh%NDOF) = CZERO
      
      ALLOCATE(LOBS(NDIM,rcv%NREC))
      IF (LGRAD) ALLOCATE(FMAT_DIST(inv%NZ_FDIST))
!
!.... first propagate the surface waves
      IF (MYID == MASTER .AND. inv%LSURF) &
      WRITE(*,*) 'fgh_window: Will propagate surface waves...' 
      IF (MYID == MASTER .AND. inv%LBODY) &
      WRITE(*,*) 'fgh_window: Will propagate body waves...'
!
!.... set space for estimates
      IF (MYNID == MASTER) THEN
         !IF (win%LWNDO_SRF) THEN
         !   NFREQ_OBS = frq%NFREQ_SRF  
         !ELSE
         !   NFREQ_OBS = frq%NFREQ_SRF_INV
         !ENDIF
         !IF (win%LWNDO_BDY) THEN
         !   NFREQ_OBS = NFREQ_OBS + frq%NFREQ_BDY 
         !ELSE
         !   NFREQ_OBS = NFREQ_OBS + frq%NFREQ_BDY_INV
         !ENDIF
         !inv%EST(1:NDIM,1:NFREQ_OBS,1:rcv%NREC,1:src%NSRC) = CZERO
         mid%RHS(1:msh%NDOF) = CZERO
      ENDIF
      CALL MPI_BCAST(NFREQ_OBS,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR)
!
!.... file handling
      load_jac = .false.
      IF (MYID == MASTER .AND. inv%LGNEWT) THEN
         IF (LJACOB_SRF .OR. LJACOB_BDY) THEN
            LEXIST = LISDIR('./scratch')
            IF (load_jac) then
               IF (.NOT.LEXIST) load_jac = .true.
            else 
               IF (LEXIST) CALL SYSTEM('rm -r ./scratch')
               CALL SYSTEM('mkdir ./scratch')
            Endif
         ENDIF
      ENDIF
      IF (MYNID == MASTER) THEN
         ALLOCATE(NSG_LSTB(NPGROUPS)) 
         ALLOCATE(NSG_LST(NPGROUPS))
         ALLOCATE(ISRCLST(src%NSRC))
      ENDIF
      CALL MPI_BARRIER(MGCOMM,MPIERR)
!
!.... loop on frequencies
      DO 1000 IFREQL=1,frq%NFREQ
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1
         IF ( (IFREQL - 1)*NPGROUPS + 1 > frq%NFREQ) GOTO 1005 !we are done
!
!....... set and send source group info  
         src%NSG = 0
         IF (MYNID == MASTER) THEN
            src%NSG = 0
            IF (IFREQ <= frq%NFREQ) THEN
               IF (frq%CFTYPE(IFREQ) == 'S') THEN !surface wave
                  IF (win%LWNDO_SRF .OR. frq%LINVF(IFREQ)) THEN
                     CALL FILL_SRCPRM_SB(IFREQ,.TRUE., m1d%VFAST_SRF,msh%AZTOL, &
                                          msh%AOITOL, SRC)
                  ENDIF
               ELSEIF (frq%CFTYPE(IFREQ) == 'B') THEN !body wave
                  IF (win%LWNDO_BDY .OR. frq%LINVF(IFREQ)) THEN
                     CALL FILL_SRCPRM_SB(IFREQ,.FALSE.,m1d%VFAST_BDY,msh%AZTOL, &
                                         msh%AOITOL, SRC)
                  ENDIF
               ELSE !no clue
                  WRITE(*,*) 'fgh_window: Cannot determine frequency type!',MYID
                  IERR = 1
                  RETURN
               ENDIF
            ENDIF !end check on frequency
            NSG_LSTB(1:NPGROUPS) = 0 
            NSG_LSTB(IPGROUP+1) = src%NSG
            CALL MPI_ALLREDUCE(src%NSG,NSGMAX,1,MPI_INTEGER, MPI_MAX,MYHD_COMM,MPIERR)
            CALL MPI_ALLREDUCE(NSG_LSTB,NSG_LST,NPGROUPS,MPI_INTEGER,MPI_SUM, &
                               MYHD_COMM,MPIERR)
         ENDIF 
         CALL MPI_BCAST(NSGMAX,1,MPI_INTEGER,MASTER, MYSLV_COMM,MPIERR)
         IF (IFREQ <= frq%NFREQ) CALL BCAST_SRC_PTRS(MYNID,MYSLV_COMM,MASTER, SRC)
!
!....... loop on source groups
         DO 2000 ISG=1,NSGMAX 
!
!.......... have master flash stats 
            IF (MYID == MASTER) THEN 
               MYSRC = 0
               PYAVG = 0.D0
               DO 201 JPGROUP=0,NPGROUPS-1
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1)
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN
                     IF (frq%CFTYPE(JFREQ) == 'S') THEN !surface wave frequency
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_SRF) THEN
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPPGRP
                              CALL MPI_RECV(PYAVG,1,MPI_DOUBLE_PRECISION,MYSRC,  &
                                            MPI_ANY_TAG, MGCOMM,STAT,MPIERR)
                           ELSE
                              PYAVG = src%PYAVG(ISG) 
                           ENDIF
                           IF (.NOT.frq%LINVF(JFREQ)) THEN  !strictly modeling
                              WRITE(*,9408) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ) 
                           ELSE !inversion
                              WRITE(*,9409) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                           ENDIF
                        ENDIF
                     ELSE !body wave frequency 
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_BDY) THEN
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPPGRP
                              CALL MPI_RECV(PYAVG,1,MPI_DOUBLE_PRECISION,MYSRC, &
                                            MPI_ANY_TAG, MGCOMM,STAT,MPIERR)
                           ELSE
                              PYAVG = src%PYAVG(ISG)
                           ENDIF 
                           IF (.NOT.frq%LINVF(JFREQ)) THEN !strictly modeling
                              WRITE(*,9410) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                           ELSE !inversion
                              WRITE(*,9411) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                           ENDIF
                        ENDIF !end check on whether to use freq
                     ENDIF !end check on frequency type 
                  ENDIF !end check on frequency range
  201          CONTINUE !loop on groups
            ELSE
               IF (MYNID == MASTER) THEN 
                  MYDEST = MASTER
                  MYTAG  = MYID
                  IF (IFREQ <= frq%NFREQ .AND. src%NSG > 0 .AND. ISG <= src%NSG) THEN
                     IF (frq%CFTYPE(IFREQ) == 'S') THEN !surface wave frequency
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_SRF) THEN
                           CALL MPI_SEND(src%PYAVG(ISG),1,MPI_DOUBLE_PRECISION, MASTER, &
                                         MYTAG, MGCOMM,MPIERR)
                        ENDIF !end check if using freq
                     ELSE !body wave frequency
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_BDY) THEN
                           CALL MPI_SEND(src%PYAVG(ISG),1,MPI_DOUBLE_PRECISION, MASTER, &
                                         MYTAG, MGCOMM,MPIERR)
                        ENDIF !end check if using freq 
                     ENDIF !end check on cftype
                  ENDIF !end check on range
               ENDIF !end check on mynid
            ENDIF
            !IF (ISG > src%NSG) GOTO 2001
            !IF (IFREQ > frq%NFREQ) GOTO 2001 !this group is done w/ things to do
!
!.......... assemble and factor
            IF (MYID == MASTER) WRITE(*,*) 'fgh_window: Assembling impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN 
               CALL ASMBLE_DRIVER(msh%NZLOC,frq%CFTYPE(IFREQ),frq%FREQ(IFREQ), &
                                  src%PYAVG(ISG), MSH, mid%A_LOC,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'fgh_window: Error calling asmble_driver!',MYID
                  GOTO 500 
               ENDIF
            ENDIF
            IF (MYID == MASTER) WRITE(*,*) 'fgh_window: Factoring impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN
               mid%JOB = 2 
               CALL CMUMPS(MID) 
               IF (mid%INFO(1) < 0) THEN
                  WRITE(*,*) 'fgh_window: An error occurred in the factorization!',MYID
                  IERR = 1 
                  GOTO 500 
               ENDIF
            ENDIF
!
!.......... have master tell which sources are being used
            LCJAC = .FALSE.
            IF (MYID == MASTER) THEN
               MYSRC = 0
               ISRCLST(1:src%NSRC) = 0
               DO 290 JPGROUP=0,NPGROUPS-1
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1)
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN
                     IF (frq%CFTYPE(JFREQ) == 'S') THEN
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_SRF) THEN
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPPGRP
                              CALL MPI_RECV(ISRCLST,src%NSRC,MPI_INTEGER,  &
                                            MYSRC,MPI_ANY_TAG, MGCOMM,STAT,MPIERR)
                           ELSE !mine
                              KSRC = 0
                              DO 291 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1) - 1
                                 KSRC = KSRC + 1
                                 ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  291                         CONTINUE
                           ENDIF
                           DO 292 KSRC=1,src%NSRC
                              ISRC = ISRCLST(KSRC)
                              IF (ISRC == 0) GOTO 282
                              WRITE(*,9421) ISRC, src%AOI(ISRC), src%BAZN(ISRC), &
                                            src%PYTAB(JFREQ,ISRC)
                              IF (frq%LINVF(JFREQ)) LCJAC = .TRUE. 
  292                      CONTINUE
  282                      CONTINUE 
                        ENDIF !end check on working
                     ELSE !body wave frequency
                        IF (frq%LINVF(JFREQ) .OR. win%LWNDO_BDY) THEN
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPPGRP
                              CALL MPI_RECV(ISRCLST,src%NSRC,MPI_INTEGER,  &
                                            MYSRC,MPI_ANY_TAG, MGCOMM,STAT,MPIERR)
                           ELSE !mine
                              KSRC = 0
                              DO 293 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1) - 1
                                 KSRC = KSRC + 1
                                 ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  293                         CONTINUE
                           ENDIF
                           DO 294 KSRC=1,src%NSRC
                              ISRC = ISRCLST(KSRC)
                              IF (ISRC == 0) GOTO 284
                              WRITE(*,9423) ISRC, src%AOI(ISRC), src%BAZN(ISRC), &
                                            src%PYTAB(JFREQ,ISRC)
                              IF (frq%LINVF(JFREQ)) LCJAC = .TRUE. 
  294                      CONTINUE 
  284                      CONTINUE
                        ENDIF !end check on working 
                     ENDIF !end check on source type
                  ENDIF !end check on bounds
  290          CONTINUE !loop on groups
            ELSE
               IF (MYNID == MASTER) THEN
                  MYDEST = MASTER
                  MYTAG = MYID
                  ISRCLST(1:src%NSRC) = 0
                  IF (IFREQ <= frq%NFREQ .AND. src%NSG > 0 .AND. ISG <= src%NSG) THEN
                     IF (frq%CFTYPE(IFREQ) == 'S') THEN !surface wave
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_SRF) THEN
                           KSRC = 0 
                           DO 297 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
                              KSRC = KSRC + 1 
                              ISRCLST(KSRC) = src%ISRCPRM(JSRC) 
  297                      CONTINUE  
                           CALL MPI_SEND(ISRCLST,src%NSRC,MPI_INTEGER, MASTER,MYTAG, &
                                         MGCOMM,MPIERR)
                        ENDIF  !end check on working
                     ELSE !body wave
                        IF (frq%LINVF(IFREQ) .OR. win%LWNDO_BDY) THEN
                           KSRC = 0
                           DO 298 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
                              KSRC = KSRC + 1
                              ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  298                      CONTINUE
                           CALL MPI_SEND(ISRCLST,src%NSRC,MPI_INTEGER, MASTER,MYTAG, &
                                         MGCOMM,MPIERR)
                        ENDIF !end check on working
                     ENDIF !end check on frequency type
                  ENDIF !end check on bounds
               ENDIF !end check on mynid
            ENDIF !end check on myid
            IF (MYID == MASTER) THEN
               LSURF = .FALSE. 
               IF (frq%CFTYPE(IFREQ) == 'S') LSURF = .TRUE.
               IF (LCJAC .AND. LGRAD) THEN
                  IF (LSURF) THEN
                     IF (LJACOB_SRF) THEN
                        WRITE(*,*)'fgh_window: Will calculate surface wave Jacobian...'
                     ELSE
                        WRITE(*,*) 'fgh_window: Will use adjgrad for surface waves...' 
                     ENDIF
                  ELSE
                     IF (LJACOB_BDY) THEN
                        WRITE(*,*) 'fgh_window: Will calculate body wave Jacobian...'
                     ELSE
                        WRITE(*,*) 'fgh_window: Will use adjgrad for body waves...'
                     ENDIF
                  ENDIF
               ENDIF 
            ENDIF
            IF (ISG > src%NSG)     GOTO 2001 !no solves, leave
            IF (IFREQ > frq%NFREQ) GOTO 2001 !no solves, leave
!----------------------------------------------------------------------------------------!
!                                 Loop on sources in source group                        !
!----------------------------------------------------------------------------------------!
!
!.......... loop on sources in source group
            DO 3000 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
               ISRC = src%ISRCPRM(JSRC) !extract original source ID
!
!............. load the 1D solution, this is where we convolve the STF
               IF (MYNID == MASTER) THEN
!
!................ this is where we convolve the source time function 
                  CALL LOAD_GRNS25(NDIM,msh%NDOF,msh%NNPE,NDIM, src%SRCTYP(ISRC),ISRC,   &
                                   frq%FREQ(IFREQ), src%SOURCE(IFREQ,ISRC),msh%IDOFSE,   &
                                   UE,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'fgh_window: Error calling load_grns25!'
                     GOTO 500
                  ENDIF
!
!................ do not want scaling in [dS/dm u] if inversion is phase only
!................ haskell will give me unit scaling; hence only applies to surface waves.
!................ problem is PML
                  !IF (inv%IRESTP == 1 .AND. frq%CFTYPE(IFREQ) == 'S') THEN 
                  !   DO IDOF=1,msh%NDOF
                  !      IF (CABS(UE(IDOF)) > 0.0) UE(IDOF) = UE(IDOF)/CABS(UE(IDOF))
                  !   ENDDO 
                  !ENDIF
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'fgh_window: Calculating effective forces...'
               ENDIF !end check on mynid
!
!............. calculate the force distribution
               CALL MPI_BCAST(UE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
               CALL ASMB_PEFF(MASTER,MYSLV_COMM, msh%NDOF,msh%NDOFL,msh%NZLOC, msh%CNP,  &
                              msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC, mid%A_LOC,UE, WAVE)
               IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,WAVE,1,mid%RHS,1)
               mid%ICNTL(9) = 1 !always solve Ax = b
               mid%JOB = 3 !solve phase
               !IF (MYID == MASTER) WRITE(*,*) 'fgh_window: Solving Ax=b...'
               CALL CMUMPS(MID)
               IF (mid%INFO(1) < 0) THEN
                  WRITE(*,*) 'fgh_window: An error occurred in the solution!',MYID
                  IERR = 1 
                  RETURN
               ENDIF
!
!............. have master processes extract the responses
               IF (MYNID == MASTER) THEN
                  CALL ADDBLK(msh%NDOF,msh%CNP,UE, mid%RHS) !add background field in
                  IF (rcv%NREC <= 0) THEN
                     WRITE(*,*) 'fgh_window: Severe error nrec <= 0',rcv%NREC,MYID
                     IERR = 1 
                     RETURN
                  ENDIF
!
!................ do not want scaling in [dS/dm u] if inversion is phase only
!................ haskell will give me unit scaling; hence only applies to surface waves.
!................ problem is PML
!                 IF (inv%IRESTP == 1 .AND. frq%CFTYPE(IFREQ) == 'S') THEN 
!                    DO IDOF=1,msh%NDOF
!                       IF (CABS(mid%RHS(IDOF)) > 0.0) &
!                       mid%RHS(IDOF) = mid%RHS(IDOF)/CABS(mid%RHS(IDOF))
!                    ENDDO 
!                 ENDIF

                  IF (inv%LUNWRAP) THEN
                     CALL EXRESP_MP_2(NDIM,rcv%NREC,msh%NDOF,                            &
                                      1,rcv%NREC,NDIM, frq%FREQ(IFREQ),                  &
                                      src%PYTAB(IFREQ,ISRC),rcv%YREC,mid%RHS, rcv%MRDOF, &
                                      rcv%RECV(1:NDIM,IFREQ,1:rcv%NREC),                 &
                                      inv%EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC))
                  ELSE
                     CALL EXRESP25D(NDIM,rcv%NREC,msh%NDOF,                              &
                                    1,rcv%NREC,NDIM, frq%FREQ(IFREQ),                    &
                                    src%PYTAB(IFREQ,ISRC),rcv%YREC, mid%RHS, rcv%MRDOF,  &
                                    rcv%RECV(1:NDIM,IFREQ,1:rcv%NREC),                   &
                                    inv%EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC))
                  ENDIF
                     !CALL PLOT_ELCRESP_VTK(NGNOD,NDIM,msh%NEN,NDIM,msh%NNPG,msh%NDOF,   &   
                     !                    msh%NLXI,msh%NLETA,'fgh', msh%NELEM,3,ISRC,   &   
                     !                    3,frq%FREQ(IFREQ), msh%LM,msh%IENG,            &   
                     !                    msh%XLOCS,msh%ZLOCS, mid%RHS) 
               ENDIF !done extracting response
!
!............. determine if this is an inversion frequency and we have observations 
               NOBS = 0 !default to not working 
               IF (MYNID == MASTER .AND. frq%LINVF(IFREQ)) THEN
                  NOBS = NOBSLST(NDIM, rcv%NREC,NDIM,       &
                                 inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC)) 
               ENDIF !end check on myid and if there are observations
               CALL MPI_BCAST(NOBS,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR)
!
!............. if there are observations at this inversion frequency calc jacobian
               IF (LGRAD .AND. NOBS > 0 .AND. frq%LINVF(IFREQ)) THEN
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'fgh_window: Calculating distributed F matrix...'
                  IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,mid%RHS,1,WAVE,1)
                  CALL MPI_BCAST(WAVE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
                  PY = src%PYTAB(IFREQ,ISRC) 
                  IF (.not.load_jac) &
                  CALL CFMAT_DIST(inv%NZ_FDIST,msh%NDOF, frq%FREQ(IFREQ),PY, &
                                  WAVE, MSH,INV, FMAT_DIST,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'fgh_window: Error calling cfmat_dist!'
                     RETURN
                  ENDIF
                  IF (MYNID == MASTER) THEN
                     CALL LOBSLST(NDIM, rcv%NREC,NDIM,                            &
                                  inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC), NOBS,LOBS)
                  ENDIF
                  CALL MPI_BCAST(NOBS,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR)
                  CALL BCAST_LOBS(MYSLV_COMM,MASTER, NDIM,NDIM,rcv%NREC, LOBS)
                  LSURF = .TRUE.
                  LBODY = .FALSE.
                  IF (frq%CFTYPE(IFREQ) == 'B') LSURF = .FALSE. 
                  IF (.NOT.LSURF) LBODY = .TRUE.
!
!................ say what we're doing
                  IF (MYID == MASTER) THEN
                     IF (LSURF.AND.LJACOB_SRF .OR. LBODY.AND.LJACOB_BDY) THEN 
                        WRITE(*,*) 'fgh_window: Calculating Jacobian...' 
                     ELSE
                        WRITE(*,*) 'fgh_window: Calculating gradient...'
                     ENDIF
                  ENDIF
!
!................ need a jacobian? 
                  IF (LSURF.AND.LJACOB_SRF .OR. LBODY.AND.LJACOB_BDY) THEN
                     IF (.not.load_jac) &
                     CALL PCHADJ(MYID,MYNID,MASTER,MYSLV_COMM, inv%NZ_FDIST,            &
                                 .TRUE.,LSURF,LPGNEW,IFREQ,ISRC,frq%FREQ(IFREQ),  LOBS, & 
                                 FMAT_DIST, RCV,INV,MID, HESSB,IERR) 
                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'fgh_window: Error calling pchadj!'
                        RETURN
                     ENDIF
                  ELSE !just a gradient/fn evaluation
                     IF (LSURF) THEN
                        LCOVD = inv%LCOVD_SRF
                        CALL ADJGRAD(MASTER,MYID,MYNID,MYSLV_COMM, inv%NZ_FDIST,msh%NDOF,&
                                     .FALSE.,LCOVD,IFREQ,ISRC, msh%AZMOD,frq%FREQ(IFREQ),&
                                     src%PYTAB(IFREQ,ISRC), &
                                     EPSS,inv%DTMAX_SRF, FMAT_DIST, RCV,MSH,INV,MID,     &
                                     GRADL_SRF,IERR) 
                     ELSE
                        LCOVD = inv%LCOVD_BDY
                        CALL ADJGRAD(MASTER,MYID,MYNID,MYSLV_COMM,inv%NZ_FDIST,msh%NDOF, &
                                     .FALSE.,LCOVD,IFREQ,ISRC, msh%AZMOD,frq%FREQ(IFREQ),&
                                     src%PYTAB(IFREQ,ISRC), &
                                     EPSS,inv%DTMAX_BDY, FMAT_DIST, RCV,MSH,INV,MID,     &
                                     GRADL_BDY,IERR) 
                     ENDIF
                     IF (IERR /= 0) THEN
                        WRITE(*,*) 'fgh_window: Error calling adjgrad!'
                        RETURN
                     ENDIF 
                  ENDIF
               ENDIF !end check 
 3000       CONTINUE !loop on sources   
 2001       CONTINUE !break ahead for loop on source group
 2000    CONTINUE !Loop on source groups
         IF (src%NSG > 0 .AND. IFREQ <= frq%NFREQ) THEN
            IF (ASSOCIATED(src%ISGPTR))   DEALLOCATE(src%ISGPTR)
            IF (ASSOCIATED(src%ISRCPRM))  DEALLOCATE(src%ISRCPRM)
            IF (ASSOCIATED(src%PYAVG))    DEALLOCATE(src%PYAVG)
         ENDIF
         CALL MPI_BARRIER(MGCOMM,MPIERR) 
 1000 CONTINUE !loop on frequencies 
 1005 CONTINUE 
      IF (ALLOCATED(NSG_LSTB)) DEALLOCATE(NSG_LSTB)
      IF (ALLOCATED(NSG_LST))  DEALLOCATE(NSG_LST)
      IF (ALLOCATED(ISRCLST))  DEALLOCATE(ISRCLST)
      CALL MPI_BARRIER(MGCOMM,MPIERR)
!
!.... take estimates back to time domain and window
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) WRITE(*,*) 'fgh_window: Reducing estimates...'
         CALL REDEST4(MYID,MASTER,MYHD_COMM,NPGROUPS,IPGROUP,                         &
                      NDIM,frq%NFREQ,rcv%NREC,                                        &
                      NDIM,frq%NFREQ,rcv%NREC,src%NSRC,2, inv%EST)
!
!....... null obsrevations shouldn't have estimates 
         DO IFREQ=1,frq%NFREQ
            DO ISRC=1,src%NSRC
               DO IREC=1,rcv%NREC
                  DO I=1,NDIM
                     IF (CABS(inv%OBS(I,IFREQ,IREC,ISRC)) == 0.0) &
                     inv%EST(I,IFREQ,IREC,ISRC) = CMPLX(0.0,0.0)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF (MYID == MASTER) THEN
            filenm(1:80) = ' '
            filenm = 'debug_surf'
            filenm = adjustl(filenm)
            CALL WRITE_EST_SB(FILENM, .TRUE., INV,FRQ,SRC,RCV, IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'fgh_window: Error calling write_est_sb 1'
               RETURN
            ENDIF
            filenm(1:80) = ' '
            filenm = 'debug_body'
            filenm = adjustl(filenm)
            CALL WRITE_EST_SB(FILENM,.FALSE., INV,FRQ,SRC,RCV, IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'fgh_window: Error calling write_est_sb 2'
               RETURN
            ENDIF
            IF (win%LWNDO_SRF) THEN
               IF (.NOT.LSRCEX .OR. .NOT.LRECEX) THEN
                  IF (LUPDSRC) &
                  WRITE(*,*) 'fgh_window: Setting initial surface wave STFs...'
                  IF (LUPDREC) &
                  WRITE(*,*) 'fgh_window: Setting initial surface wave RRFs...'
                  CALL UPD_SRCREC_WNDO(.TRUE. ,LUPDSRC,LUPDREC, msh%AZMOD,  &
                                       FRQ,INV,SRC,WIN,RCV)
               ENDIF
            ENDIF
            IF (win%LWNDO_BDY) THEN
               IF (LUPDSRC .OR. LUPDREC) THEN 
                  IF (LUPDSRC) &
                  WRITE(*,*) 'fgh_window: Setting initial body wave STFs...'
                  IF (LUPDREC) &
                  WRITE(*,*) 'fgh_window: Setting initial body wave RRFs...'
                  CALL UPD_SRCREC_WNDO(.FALSE.,LUPDSRC,LUPDREC, msh%AZMOD,  &
                                       FRQ,INV,SRC,WIN,RCV)
               ENDIF
            ENDIF
!
!.......... window the estimates  
            IF (win%LWNDO_SRF) THEN 
               WRITE(*,*) 'fgh_window: Windowing time domain surface wave synthetics...'
               CALL WINDOW_DRIVER(IBLOCK,K,IALPHA, .TRUE., INV,FRQ,SRC,RCV,WIN, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'fgh_window: Error calling window_driver 1'
                  RETURN
               ENDIF
            ENDIF
            IF (win%LWNDO_BDY) THEN
               WRITE(*,*) 'fgh_window: Windowing time domain body wave synthetics...'
               CALL WINDOW_DRIVER(IBLOCK,K,IALPHA,.FALSE., INV,FRQ,SRC,RCV,WIN, IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'fgh_window: Error calling window_driver 2'
                  RETURN
               ENDIF 
            ENDIF
            filenm(1:80) = ' '
            filenm = 'debug_surf2'
            filenm = adjustl(filenm)
            CALL WRITE_EST_SB(FILENM, .TRUE., INV,FRQ,SRC,RCV, IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'fgh_window: Error calling write_est_sb 3!'
               RETURN
            ENDIF
            filenm(1:80) = ' '
            filenm = 'debug_body2'
            filenm = adjustl(filenm)
            CALL WRITE_EST_SB(FILENM,.FALSE., INV,FRQ,SRC,RCV, IERR)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'fgh_window: Error calling write_est_sb 4!'
               RETURN
            ENDIF
!
!.......... update the STFs and RRFs in the windows
            IF (LUPDSRC .OR. LUPDREC) THEN
               CALL UPD_SRCREC_WNDO(.TRUE. ,LUPDSRC,LUPDREC, msh%AZMOD, FRQ,INV,SRC,WIN,RCV)
               CALL UPD_SRCREC_WNDO(.FALSE.,LUPDSRC,LUPDREC, msh%AZMOD, FRQ,INV,SRC,WIN,RCV)
            ENDIF
            LSTF = .FALSE. !assume we havent had to do anything
            IF (LUPDSRC) LSTF = .TRUE.
!
!.......... calculate the objective function
            IF (LFUNC) THEN
               inv%FOBJ_SRF = 0.0
               inv%FOBJ_BDY = 0.0
               IF (inv%LSURF) THEN
                  WRITE(*,*) 'fgh_window: Calculating surface wave objective function...'
                  CALL COBJ_SB( .TRUE.,inv%DTMAX_SRF,msh%AZMOD, INV,FRQ,SRC,RCV,  &
                               OBJ, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'fgh_window: Error calling cobj_sb 1!'
                     RETURN
                  ENDIF
                  inv%FOBJ_SRF = OBJ
               ENDIF
               IF (inv%LBODY) THEN
                  WRITE(*,*) 'fgh_window: Calculating body wave objective function...'
                  CALL COBJ_SB(.FALSE.,inv%DTMAX_BDY,msh%AZMOD, INV,FRQ,SRC,RCV,  &
                               OBJ, IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'fgh_window: Error calling cobj_sb 2!'
                     RETURN
                  ENDIF
                  inv%FOBJ_BDY = OBJ
               ENDIF 
               inv%FOBJ = inv%FOBJ_SRF + inv%SCLBDY*inv%FOBJ_BDY 
               WRITE(*,*) 'fgh_window: Total objective function:',inv%FOBJ 
            ENDIF
         ENDIF
!
!....... calculate gradient?
         CALL MPI_BARRIER(MYHD_COMM,MPIERR)
         IF (LGRAD) THEN
            CALL MPI_BCAST(LSTF,1,MPI_LOGICAL, MASTER,MYHD_COMM,MPIERR) 
            IF (inv%LSURF .AND. LJACOB_SRF) THEN !adj(J)*res 
               CALL MULJAC_RESID_MPI(MYID,MASTER,MYHD_COMM, NPPGRP, &
                                     .TRUE. ,LSTF, msh%AZMOD, INV,FRQ,SRC,RCV, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'fgh_window: Error calling muljac_resid_mpi1',MYID
                  RETURN
               ENDIF 
            ELSE !can do gradient with adjiont
               IF (inv%LSURF) THEN
                  IF (MYID == MASTER) THEN
                     ALLOCATE(WORK(inv%NA35))
                  ELSE
                     ALLOCATE(WORK(1))
                  ENDIF
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'fgh_window: Reducing surface wave gradient...'
                  CALL MPI_REDUCE(GRADL_SRF,WORK,inv%NA35,MPI_REAL,MPI_SUM, &
                                  MASTER,MYHD_COMM,MPIERR)
                  IF (MYID == MASTER) CALL SCOPY(inv%NA35,WORK,1,inv%GRAD_SRF,1)
                  DEALLOCATE(WORK)
               ENDIF
            ENDIF
            IF (inv%LBODY .AND. LJACOB_BDY) THEN !adj(J)*res 
               CALL MULJAC_RESID_MPI(MYID,MASTER,MYHD_COMM, NPPGRP,     &
                                     .FALSE.,LSTF, msh%AZMOD, INV,FRQ,SRC,RCV, IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'fgh_window: Error calling muljac_resid_mpi',MYID
                  RETURN
               ENDIF
            ELSE !can do gradient with adjoint method
               IF (inv%LBODY) THEN 
                  IF (MYID == MASTER) THEN 
                     ALLOCATE(WORK(inv%NA35))
                  ELSE 
                     ALLOCATE(WORK(1))
                  ENDIF
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'fgh_window: Reducing body wave gradient...'
                  CALL MPI_REDUCE(GRADL_BDY,WORK,inv%NA35,MPI_REAL,MPI_SUM, &
                                  MASTER,MYHD_COMM,MPIERR)
                  IF (MYID == MASTER) CALL SCOPY(inv%NA35,WORK,1,inv%GRAD_BDY,1)
                  DEALLOCATE(WORK)
               ENDIF 
            ENDIF !end check on type of gradient calculation
         ENDIF !end check on if we need a gradient calculation
         IF (MYID == MASTER .AND. inv%LSURF .AND. inv%LBODY) THEN 
            inv%GRAD(1:inv%NA35) = inv%GRAD_SRF(1:inv%NA35) &
                                 + SNGL(inv%SCLBDY)*inv%GRAD_BDY(1:inv%NA35)
         ENDIF
      ENDIF !end check on mynid
!
!.... clean
  500 CONTINUE
      IF (ALLOCATED(WAVE))      DEALLOCATE(WAVE)
      IF (ALLOCATED(UE))        DEALLOCATE(UE) 
      IF (ALLOCATED(LOBS))      DEALLOCATE(LOBS) 
      IF (ALLOCATED(FMAT_DIST)) DEALLOCATE(FMAT_DIST)
      IF (ALLOCATED(GRADL_SRF)) DEALLOCATE(GRADL_SRF)
      IF (ALLOCATED(GRADL_BDY)) DEALLOCATE(GRADL_BDY) 
!
!.... format statements
 9408 FORMAT(' ----------------------------------------------------------',/,&
             ' - fgh_window: Group:',I4,'                                 -'   ,/,& 
             ' -             Propagating surface wave source group ',I3,' -'  ,/,& 
             ' -             With average slowness in y',E12.4,'     - '       ,/,& 
             ' -             At frequency',F12.5,' Hz                -'        ,/,& 
             ' ----------------------------------------------------------',/)
 9409 FORMAT(' ----------------------------------------------------------',/,&
            ' - fgh_window: Group:',I4,'                                 -'   ,/,& 
            ' -             Propagating surface wave source group ',I3,'  -'  ,/,& 
            ' -             With average slowness in y',E12.4,'     - '       ,/,& 
            ' -             At frequency',F12.5,' Hz                -'        ,/,& 
            ' -             This is an inversion frequency             -',/,&
            ' ----------------------------------------------------------',/)

 9410 FORMAT(' ----------------------------------------------------------',/,&
             ' - fgh_window: Group:',I4,'                                 -'   ,/,& 
             ' -             Propagating body wave source group ',I5,'   -'    ,/,& 
             ' -             With average slowness in y',E12.4,'     - '       ,/,& 
             ' -             At frequency',F12.5,' Hz                -'        ,/,& 
             ' ----------------------------------------------------------',/)
 9411 FORMAT(' ----------------------------------------------------------',/,&
             ' - fgh_window: Group:',I4,'                                 -'   ,/,& 
             ' -             Propagating body wave source group ',I5,'   -'    ,/,& 
             ' -             With average slowness in y',E12.4,'     - '       ,/,& 
             ' -             At frequency',F12.5,' Hz                -'        ,/,& 
             ' -             This is an inversion frequency             -',/,&
             ' ----------------------------------------------------------',/)
 9421 FORMAT(/,' fgh_window: Will process surface wave source,',I3,/   &
               '             With angle of incidence ',F8.3,' degrees',/, &
               '             Adjusted azimuth ',F8.3,' degrees',/,    &    
               '             And y slowness ',G10.3,' s/m',/)
 9423 FORMAT(/,' fgh_window: Will process body wave source,',I3,/   &    
               '             With angle of incidence ',F8.3,' degrees',/, &
               '             Adjusted azimuth ',F8.3,' degrees',/,    &    
               '             And y slowness ',G10.3,' s/m',/)

      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE WRITE_EST_SB(FILENM,LSURF, INV,FRQ,SRC,RCV, IERR) 
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV  
      TYPE (FRQ_INFO) FRQ
      TYPE (SRC_INFO) SRC
      TYPE (RECV_INFO) RCV
      CHARACTER(80), INTENT(IN) :: FILENM
      LOGICAL*4, INTENT(IN) :: LSURF
      INTEGER*4, INTENT(OUT) :: IERR
!.... variable declarations
      COMPLEX*8, ALLOCATABLE :: EST_LOC(:,:,:,:) 
      REAL*8, ALLOCATABLE :: FREQ_LST(:)
      COMPLEX*8 CZERO
      INTEGER*4 NFREQ, NSRC, NREC
      LOGICAL*4 LBODY
      PARAMETER(CZERO = CMPLX(0.0,0.0))
!
!----------------------------------------------------------------------------------------!
!
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
      NREC = rcv%NREC
      IF (LSURF) THEN
         NFREQ = frq%NFREQ_SRF
         NSRC  = src%NSRC_SRF
      ELSE
         NFREQ = frq%NFREQ_BDY
         NSRC = src%NSRC_BDY
      ENDIF
      IF (NFREQ == 0 .OR. NSRC == 0) RETURN
      NREC = rcv%NREC
      ALLOCATE(EST_LOC(NDIM,NFREQ,NREC,NSRC))
      ALLOCATE(FREQ_LST(NFREQ))
      EST_LOC(1:NDIM,1:NFREQ,1:NREC,1:NSRC) = CZERO
      IF (LSURF) THEN
         CALL GET_EST_SB( .TRUE., NDIM,NFREQ,NREC, INV,FRQ,SRC,RCV, EST_LOC,IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'write_est_sb: Error calling get_est_sb 1!'
            RETURN
         ENDIF
         CALL GET_FREQ_MOD(frq%NFREQ, NFREQ, .TRUE., frq%CFTYPE,frq%FREQ, &
                           FREQ_LST,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'write_est_sb: Error calling get_freq_mod 1!' 
            RETURN
         ENDIF
      ELSE
         CALL GET_EST_SB(.FALSE., NDIM,NFREQ,NREC, INV,FRQ,SRC,RCV, EST_LOC,IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'write_est_sb: Error calling get_est_sb 2!' 
            RETURN
         ENDIF
         CALL GET_FREQ_MOD(frq%NFREQ, NFREQ, .FALSE., frq%CFTYPE,frq%FREQ, &
                           FREQ_LST,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'write_est_sb: Error calling get_freq_mod 2!'
            RETURN
         ENDIF
      ENDIF
!
!.... write estimates
      CALL WTEEST25(FILENM, NFREQ,NREC,NFREQ, NREC,NSRC,      &    
                    FREQ_LST,EST_LOC(1,:,:,:),EST_LOC(2,:,:,:),EST_LOC(3,:,:,:), IERR)
      IF (IERR /= 0) WRITE(*,*) 'write_est_sb: Error calling wteest25'
      DEALLOCATE(FREQ_LST)
      DEALLOCATE(EST_LOC) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE UPD_SRCREC_WNDO(LSURF,LUPDSRC,LUPDREC, AZMOD, FRQ,INV,SRC,WIN,RCV) 
!
!     Updates the source time and receiver respone functions for surface or body waves
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ 
      TYPE (INV_INFO) INV 
      TYPE (SRC_INFO) SRC 
      TYPE (WIN_INFO) WIN
      TYPE (RECV_INFO) RCV
      REAL*8, INTENT(IN) :: AZMOD
      LOGICAL*4, INTENT(IN) :: LSURF, LUPDSRC, LUPDREC 
      COMPLEX*8, ALLOCATABLE :: EST_WRK(:,:,:,:), OBS_WRK(:,:,:,:), RCV_WRK(:,:,:), &
                                SRC_WRK(:,:) 
      COMPLEX*8 CFACT 
      INTEGER*4 NFREQ, NSRC, IFREQ,JFREQ,IREC,ISRC,JSRC,I
      LOGICAL*4 LBODY, LSKIP 
!
!.... quick return
      IF (.NOT.LUPDSRC .AND. .NOT.LUPDREC) THEN
         WRITE(*,*) 'upd_srcrec_wndo: Nothing to update, returning'
         RETURN
      ENDIF
!
!.... body or surface waves?
      LBODY = .FALSE.  
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... determine size and set size
      IF (LSURF) THEN
         IF (win%LWNDO_SRF) THEN
            NFREQ = frq%NFREQ_SRF
         ELSE
            NFREQ = frq%NFREQ_SRF_INV
         ENDIF
         NSRC = src%NSRC_SRF
      ELSE
         IF (win%LWNDO_BDY) THEN
            NFREQ = frq%NFREQ_BDY
         ELSE 
            NFREQ = frq%NFREQ_BDY_INV
         ENDIF
         NSRC = src%NSRC_BDY
      ENDIF 
      ALLOCATE(OBS_WRK(NDIM,NFREQ,rcv%NREC,NSRC))
      ALLOCATE(EST_WRK(NDIM,NFREQ,rcv%NREC,NSRC))
      ALLOCATE(RCV_WRK(NDIM,NFREQ,rcv%NREC))
      ALLOCATE(SRC_WRK(NFREQ,NSRC))
!
!.... copy vectors
      JFREQ = 0       
      DO 1 IFREQ=1,frq%NFREQ 
         LSKIP = .FALSE. 
         IF (LSURF) THEN
            IF (.NOT.win%LWNDO_SRF .AND. .NOT.frq%LINVF(IFREQ)) LSKIP = .TRUE.
            IF (frq%CFTYPE(IFREQ) == 'B') LSKIP = .TRUE.
         ENDIF
         IF (LBODY) THEN
            IF (.NOT.win%LWNDO_BDY .AND. .NOT.frq%LINVF(IFREQ)) LSKIP = .TRUE.
            IF (frq%CFTYPE(IFREQ) == 'S') LSKIP = .TRUE. 
         ENDIF
         IF (LSKIP) GOTO 100
         JFREQ = JFREQ + 1
!
!....... loop on sources
         JSRC = 0 
         DO 2 ISRC=1,NSRC  
            IF (LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
                LBODY .AND. src%SRCTYP(ISRC)(1:1) == 'B') THEN
               JSRC = JSRC + 1
               DO 3 IREC=1,rcv%NREC 
                  DO 4 I=1,NDIM
                     EST_WRK(I,JFREQ,IREC,JSRC) = inv%EST(I,IFREQ,IREC,ISRC)
                     OBS_WRK(I,JFREQ,IREC,JSRC) = inv%OBS(I,IFREQ,IREC,ISRC) 
                     RCV_WRK(I,JFREQ,IREC) = CMPLX(1.0,0.0)
    4             CONTINUE
    3          CONTINUE
               SRC_WRK(JFREQ,JSRC) = CMPLX(1.0,0.0)
            ENDIF
    2    CONTINUE
  100    CONTINUE !break ahead nothing to do 
    1 CONTINUE 
!
!.... perform the updates
      IF (LUPDSRC) THEN
         WRITE(*,*) 'upd_srcrec_wndo: Solving new STFs...'
         CALL SRCUPD(NDIM,NFREQ,rcv%NREC, NDIM,NFREQ,rcv%NREC,NSRC, &
                     inv%LUNWRAP,inv%IMODSRC, AZMOD, EST_WRK,OBS_WRK, SRC_WRK)
         !IF (inv%IRESTP == 1) THEN
         !   WRITE(*,*) 'upd_srcrec_wndo: Setting source amplitude to unity...'
         !ELSEIF (inv%IRESTP == 2) THEN
         !   WRITE(*,*) 'upd_srcrec_wndo: Setting source phase to zero...'
         !ENDIF
         !CALL MODSRC(NFREQ, NFREQ,NSRC,inv%IRESTP, SRC_WRK)
      ENDIF
      IF (LUPDREC) THEN
         WRITE(*,*) 'upd_srcrec_wndo: Solving new RRFs...'
         CALL RECUPD(NDIM,NFREQ,rcv%NREC, NDIM,NFREQ,rcv%NREC,NSRC, &
                     inv%LUNWRAP,inv%IMODREC, AZMOD, EST_WRK,OBS_WRK, RCV_WRK)
         !IF (inv%IRESTP == 1) THEN
         !   WRITE(*,*) 'upd_srcrec_wndo: Setting receiver amplitude to unity...'
         !ELSEIF (inv%IRESTP == 2) THEN
         !   WRITE(*,*) 'upd_srcrec_wndo: Setting receiver phase to zero...'
         !ENDIF
         !CALL MODREC(NDIM,NFREQ, NDIM,NFREQ,rcv%NREC, inv%IRESTP, RCV_WRK)
      ENDIF 
!
!.... convolve back 
      JFREQ = 0
      DO 11 IFREQ=1,frq%NFREQ
         LSKIP = .FALSE.
         IF (LSURF) THEN
            IF (.NOT.win%LWNDO_SRF .AND. .NOT.frq%LINVF(IFREQ)) LSKIP = .TRUE.
            IF (frq%CFTYPE(IFREQ) == 'B') LSKIP = .TRUE.
         ENDIF
         IF (LBODY) THEN
            IF (.NOT.win%LWNDO_BDY .AND. .NOT.frq%LINVF(IFREQ)) LSKIP = .TRUE.
            IF (frq%CFTYPE(IFREQ) == 'S') LSKIP = .TRUE. 
         ENDIF
         IF (LSKIP) GOTO 110 
         JFREQ = JFREQ + 1 
!
!....... loop on sources
         JSRC = 0 
         DO 12 ISRC=1,NSRC  
            IF (LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
                LBODY .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
               JSRC = JSRC + 1 
               DO 13 IREC=1,rcv%NREC 
                  DO 14 I=1,NDIM
                     CFACT = RCV_WRK(I,JFREQ,IREC)*SRC_WRK(JFREQ,JSRC)
                     inv%EST(I,IFREQ,IREC,ISRC) = EST_WRK(I,JFREQ,IREC,JSRC)*CFACT
   14             CONTINUE
   13          CONTINUE
            ENDIF
   12    CONTINUE
  110    CONTINUE !break ahead nothing to do 
   11 CONTINUE 
      DEALLOCATE(OBS_WRK)
      DEALLOCATE(EST_WRK)  
      DEALLOCATE(SRC_WRK)
      DEALLOCATE(RCV_WRK) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FILL_RESID_SB(NWORK,IBPHASE, LSURF, AZMOD, INV,FRQ,RCV,SRC,WIN,  &
                               OBJ,RESID, IERR) 
!
!     Fills the surface wave or body wave residual vector 
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO)  INV
      TYPE (FRQ_INFO)  FRQ
      TYPE (RECV_INFO) RCV
      TYPE (SRC_INFO)  SRC
      TYPE (WIN_INFO)  WIN
      REAL*8, INTENT(IN) :: AZMOD
      INTEGER*4, INTENT(IN) :: NWORK, IBPHASE 
      LOGICAL*4, INTENT(IN) :: LSURF 
      REAL*4, INTENT(OUT) :: RESID(NWORK), OBJ 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      COMPLEX*8 OBSR(3), EST(3), CZERO, DIF, CRESID, HFACT, HDIFFER, QN, QE, QZ 
      REAL*4 W2, RMAG, PHASE, XNEG 
      LOGICAL*4 LBODY 
      INTEGER*4 JNDX, JNDX0, IFREQ,JFREQ,KFREQ, ISRC,IREC,I,INDX 
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      COMPLEX*8 CPHM2CM  
!
!----------------------------------------------------------------------------------------!
!
!.... body or surface waves?
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
!
!.... offset s.t. the first set of observations are the reals, second imaginary
      INDX = 0 
      IF (LSURF) THEN
         JNDX = NDIM*rcv%NREC*src%NSRC_SRF*frq%NFREQ_SRF_INV  
      ELSE
         JNDX = NDIM*rcv%NREC*src%NSRC_BDY*frq%NFREQ_BDY_INV
      ENDIF
      IF (JNDX == 0) THEN
         IF (LSURF) THEN
            WRITE(*,*) 'fill_resid_sb: No surface wave residuals to calculate'
            RETURN
         ELSE
            WRITE(*,*) 'fill_resid_sb: No body wave residuals to calculate'
            RETURN
         ENDIF
      ENDIF 
      JNDX0 = JNDX  
!
!.... loop on frequencies, sources, receivers, components
      OBJ = 0.0 
      XNEG =-HUGE(1.0) 
      RESID(1:NWORK) = XNEG 
      DO 1 IFREQ=1,frq%NFREQ
         IF (LSURF .AND. frq%CFTYPE(IFREQ) == 'B') GOTO 100
         IF (LBODY .AND. frq%CFTYPE(IFREQ) == 'S') GOTO 100  
!
!....... is this an inversion frequency
         IF (frq%LINVF(IFREQ)) THEN
!
!.......... locate the frequency
            KFREQ = 0
            DO 21 JFREQ=1,frq%NFREQ 
               IF (frq%CFTYPE(JFREQ) == 'S') THEN 
                  IF (win%LWNDO_SRF .OR. frq%LINVF(JFREQ)) KFREQ = KFREQ + 1
               ELSE 
                  IF (win%LWNDO_BDY .OR. frq%LINVF(JFREQ)) KFREQ = KFREQ + 1
               ENDIF
               IF (JFREQ == IFREQ) GOTO 180
   21       CONTINUE  
  180       CONTINUE 
            HFACT = HDIFFER(IBPHASE,frq%FREQ(IFREQ)) !option for half-differentiator
            DO 2 ISRC=1,src%NSRC
               IF (LSURF .AND. src%SRCTYP(ISRC)(1:1) == 'S' .OR.  &
                   LBODY .AND. src%SRCTYP(ISRC)(1:1) == 'P') THEN
                  DO 3 IREC=1,rcv%NREC
!
!................... extract residuals and rotate
                     QN = CZERO
                     QE = CZERO
                     QZ = CZERO 
                     IF (inv%LUNWRAP) THEN
                        RMAG  = REAL(inv%EST(1,KFREQ,IREC,ISRC))
                        PHASE = IMAG(inv%EST(1,KFREQ,IREC,ISRC))
                        EST(1) = CPHM2CM(RMAG,PHASE)
                        RMAG  = REAL(inv%EST(2,KFREQ,IREC,ISRC))
                        PHASE = IMAG(inv%EST(2,KFREQ,IREC,ISRC))
                        EST(2) = CPHM2CM(RMAG,PHASE)
                        RMAG  = REAL(inv%EST(3,KFREQ,IREC,ISRC))
                        PHASE = IMAG(inv%EST(3,KFREQ,IREC,ISRC))
                        EST(3) = CPHM2CM(RMAG,PHASE)
                        IF (CABS(inv%OBS(1,KFREQ,IREC,ISRC)) > 0.0) THEN 
                           RMAG  = REAL(inv%OBS(1,KFREQ,IREC,ISRC))
                           PHASE = IMAG(inv%OBS(1,KFREQ,IREC,ISRC))
                           QN = CPHM2CM(RMAG,PHASE)  
                        ENDIF
                        IF (CABS(inv%OBS(2,KFREQ,IREC,ISRC)) > 0.0) THEN 
                           RMAG  = REAL(inv%OBS(2,KFREQ,IREC,ISRC))
                           PHASE = IMAG(inv%OBS(2,KFREQ,IREC,ISRC))
                           QE = CPHM2CM(RMAG,PHASE) 
                        ENDIF
                        IF (CABS(inv%OBS(3,KFREQ,IREC,ISRC)) > 0.0) THEN 
                           RMAG  = REAL(inv%OBS(3,KFREQ,IREC,ISRC))
                           PHASE = IMAG(inv%OBS(3,KFREQ,IREC,ISRC))
                           QZ = CPHM2CM(RMAG,PHASE)  
                        ENDIF
                     ELSE
                        IF (CABS(inv%OBS(1,KFREQ,IREC,ISRC)) > 0.0) &
                        QN = inv%OBS(1,KFREQ,IREC,ISRC) 
                        IF (CABS(inv%OBS(2,KFREQ,IREC,ISRC)) > 0.0) & 
                        QE = inv%OBS(2,KFREQ,IREC,ISRC) 
                        IF (CABS(inv%OBS(3,KFREQ,IREC,ISRC)) > 0.0) &
                        QZ = inv%OBS(3,KFREQ,IREC,ISRC)
                        EST(1) = inv%EST(1,KFREQ,IREC,ISRC) 
                        EST(2) = inv%EST(2,KFREQ,IREC,ISRC)
                        EST(3) = inv%EST(3,KFREQ,IREC,ISRC)
                     ENDIF
                     CALL CROTATE(-SNGL(AZMOD),QN,QE, OBSR(1),OBSR(2)) 
                     OBSR(3) = QZ
                     DO 4 I=1,NDIM
                        DIF = CZERO 
                        W2 = inv%WGHTS(I,KFREQ,IREC,ISRC)**2  
                        IF (CABS(OBSR(I)) > 0.0) THEN
                           DIF = CRESID(inv%IRESTP,OBSR(I), EST(I)) 
                        ENDIF
                        INDX = INDX + 1
                        JNDX = JNDX + 1
                        OBJ = OBJ + W2*CABS(DIF)**2
                        DIF = DIF*HFACT
                        RESID(INDX) = W2*REAL(DIF)
                        RESID(JNDX) = W2*IMAG(DIF)
                        write(48,*) kfreq,irec,isrc,i,indx,jndx,resid(indx),resid(jndx)
    4                CONTINUE !loop on components 
    3             CONTINUE ! loop on receiver
               ENDIF !end check on body wave surface source match
    2       CONTINUE !loop on sources
         ENDIF !end check on inversion frequency
  100    CONTINUE !break ahead, not in right spot
    1 CONTINUE !loop on frequencies 
      print *, nwork,indx,jndx, minval(resid),maxval(resid)
      IF (INDX /= JNDX0) THEN
         WRITE(*,*) 'fill_resid_sb: You have an indexing error!'
         IERR = 1
         RETURN
      ENDIF
!
!.... check on initialization
      IF (MINVAL(RESID(1:JNDX)) == XNEG) THEN
         WRITE(*,*) 'fill_resid_sb: You have missed an index!'
         IERR = 1
         RETURN
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE WINDOW_DRIVER(IBLOCK,K,IALPHA, LSURF, INV,FRQ,SRC,RCV,WIN, IERR)
!
!     Handles the windowing
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h' 
      TYPE (INV_INFO)   INV
      TYPE (FRQ_INFO)   FRQ
      TYPE (SRC_INFO)   SRC
      TYPE (WIN_INFO)   WIN 
      TYPE(RECV_INFO)   RCV 
      INTEGER*4, INTENT(IN) :: IBLOCK,K,IALPHA
      LOGICAL*4, INTENT(IN) :: LSURF
      INTEGER*4, INTENT(OUT) :: IERR

      COMPLEX*8, ALLOCATABLE :: EST_LOC(:,:,:,:)
      REAL*8, ALLOCATABLE :: FREQ8(:)
      REAL*4, ALLOCATABLE :: WIGGLE(:,:,:,:), TWIN(:,:,:), FREQ_LOC(:)
      CHARACTER(80) PROJNM 
      CHARACTER(5) CBLOCK, CK, CALPHA
      COMPLEX*8 CPHM2CM 
      REAL*4 STARTT, DT, RMAG, PHASE, SPHASE
      INTEGER*4 NREC, NFREQ_WRK,NSRC_WRK, NSAMP, IFREQ,ISRC,IREC,I,IDIR,NPCT  
      LOGICAL*4 LBODY, LEX, LISDIR  
!
!----------------------------------------------------------------------------------------!
!
!.... quick return 
      IERR = 0
      LBODY = .FALSE.
      IF (.NOT.LSURF) LBODY = .TRUE.
      IF (.NOT.win%LWNDO_SRF .AND. LSURF) THEN
         WRITE(*,*) 'window_driver: Windowing turned off for surface waves, returning...'
         RETURN
      ENDIF
      IF (.NOT.win%LWNDO_BDY .AND. LBODY) THEN
         WRITE(*,*) 'window_driver: Windowing turned off for body waves, returning...'
         RETURN
      ENDIF
!
!.... initialize
      IF (LSURF) THEN
         NFREQ_WRK = frq%NFREQ_SRF
         NSRC_WRK  = src%NSRC_SRF
         NSAMP = win%NSAMP_SRF
         STARTT = win%START_SRF
         DT = win%DT_SRF
         NPCT = win%NPCT_SRF
      ELSE
         NFREQ_WRK = frq%NFREQ_BDY
         NSRC_WRK  = src%NSRC_BDY
         NSAMP = win%NSAMP_BDY
         STARTT = win%START_BDY
         DT = win%DT_BDY
         NPCT = win%NPCT_BDY
      ENDIF
      IF (NFREQ_WRK <= 0 .OR. NSRC_WRK <= 0) THEN
         IF (LSURF) THEN
            WRITE(*,*) 'window_driver: Skipping windowing for surface waves...'
            RETURN
         ELSE
            WRITE(*,*) 'window_driver: Skipping windowing for body waves...'
            RETURN
         ENDIF
      ENDIF
!
!.... extract the estimates
      NREC = rcv%NREC
      ALLOCATE(EST_LOC(NDIM,NFREQ_WRK,NREC,NSRC_WRK))
      CALL GET_EST_SB(LSURF, NDIM,NFREQ_WRK,NREC, INV,FRQ,SRC,RCV, EST_LOC,IERR) 
      IF (IERR /= 0) THEN
         WRITE(*,*) 'xwindow_driver: Error calling get_est_sb!'
         RETURN
      ENDIF
!
!.... get the windows
      ALLOCATE(TWIN(NREC,NSRC_WRK,2))
      CALL GET_TWIN(NREC,NSRC_WRK, NREC, LSURF,NSRC_WRK, WIN,SRC, TWIN,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'window_driver: Error calling get_twin!'
         RETURN
      ENDIF
!
!.... get the frequency list
      ALLOCATE(FREQ8(NFREQ_WRK))
      ALLOCATE(FREQ_LOC(NFREQ_WRK))
      CALL GET_FREQ_MOD(frq%NFREQ, NFREQ_WRK, LSURF, frq%CFTYPE,frq%FREQ, &
                        FREQ8,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'window_driver: Error calling get_freq_mod!'
         RETURN
      ENDIF
      FREQ_LOC(1:NFREQ_WRK) = SNGL(FREQ8(1:NFREQ_WRK))
!
!.... unpack phase
      IF (inv%LUNWRAP) THEN
         DO 11 IFREQ=1,NFREQ_WRK
            DO 12 ISRC=1,NSRC_WRK  
               DO 13 IREC=1,NREC
                  DO 14 I=1,NDIM
                     RMAG  = REAL(EST_LOC(I,IFREQ,IREC,ISRC))
                     PHASE = IMAG(EST_LOC(I,IFREQ,IREC,ISRC))
                     EST_LOC(I,IFREQ,IREC,ISRC) = CPHM2CM(RMAG,PHASE)
   14             CONTINUE
   13          CONTINUE
   12       CONTINUE
   11    CONTINUE
      ENDIF
!
!.... window
      WRITE(*,*) 'window_driver: Inverse Fourier transforming...'
      ALLOCATE(WIGGLE(NDIM,NSAMP,NREC,NSRC_WRK))
      IDIR = 1 !inverse transform 
      CALL DFTWIG4(NDIM,NSAMP,NFREQ_WRK, NREC,                     &
                   NSAMP,NDIM,NFREQ_WRK, NREC,NSRC_WRK, IDIR,DT,   &
                   FREQ_LOC,EST_LOC, WIGGLE) 
!
!.... directory handling
      LEX = LISDIR('./pest_inv') 
      IF (.NOT.LEX) CALL SYSTEM('mkdir ./pest_inv')
!
!.... dump the time domain windowed estimates 
      PROJNM(1:80) = ' '
      CALPHA(1:5) = ' '
      CK(1:5) = ' '
      CBLOCK(1:5) = ' '
      WRITE(CALPHA,'(I5)') IALPHA
      WRITE(CK,'(I5)') K
      WRITE(CBLOCK,'(I5)') IBLOCK
      CALPHA = ADJUSTL(CALPHA)
      CK = ADJUSTL(CK) 
      CBLOCK = ADJUSTL(CBLOCK) 
      IF (IALPHA > 0) THEN 
         PROJNM = './pest_inv/body-'//TRIM(CBLOCK)//'-'//TRIM(CK)//'-'//TRIM(CALPHA)
      ELSE 
         PROJNM = './pest_inv/body-'//TRIM(CBLOCK)//'-'//TRIM(CK)
      ENDIF
      PROJNM = ADJUSTL(PROJNM)
      CALL WTSEISM(PROJNM, NDIM,NSAMP,NREC, NDIM,NSAMP,NREC,NSRC_WRK,  &
                   DT,STARTT, WIGGLE, IERR)

      WRITE(*,*) 'window_driver: Windowing data...'
      CALL WINDOW_TS(NDIM,NSAMP,NREC,NSRC_WRK, NDIM,NSAMP,NREC,NSRC_WRK, &
                     win%IWNDO,NPCT, DT,TWIN, WIGGLE,IERR) 
      IF (IERR /= 0) THEN  
         WRITE(*,*) 'window_driver: Fourier transforming data...'  
         RETURN
      ENDIF
      WRITE(*,*) 'window_driver: Transforming back to frequency domain...'
      IDIR =-1 !forward transform
      CALL DFT_TF4(NDIM,NSAMP,NFREQ_WRK,NREC,                 &
                   NDIM,NSAMP,NFREQ_WRK,NREC, NSRC_WRK, IDIR, &
                   DT,STARTT, FREQ_LOC,WIGGLE, EST_LOC) 
!
!.... rewrap phase?
      IF (inv%LUNWRAP) THEN
         DO 21 IFREQ=1,NFREQ_WRK
            DO 22 ISRC=1,NSRC_WRK  
               DO 23 IREC=1,NREC
                  DO 24 I=1,NDIM
                     RMAG  = CABS(EST_LOC(I,IFREQ,IREC,ISRC))
                     PHASE = SPHASE(EST_LOC(I,IFREQ,IREC,ISRC))
                     EST_LOC(I,IFREQ,IREC,ISRC) = CMPLX(RMAG,PHASE) 
   24             CONTINUE
   23          CONTINUE
   22       CONTINUE
   21    CONTINUE
      ENDIF
!
!.... put windowed estimates back into inv%EST
      CALL SET_EST_SB(LSURF, NDIM,NFREQ_WRK,NREC, EST_LOC, INV,FRQ,SRC,RCV, IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'window_driver: Error calling set_est_sb!'
         RETURN
      ENDIF
!
!.... file name 
      PROJNM(1:80) = ' '
      IF (IALPHA > 0) THEN 
         PROJNM = './pest_inv/body_win-'//TRIM(CBLOCK)//'-'//TRIM(CK)//'-'//TRIM(CALPHA)
      ELSE
         PROJNM = './pest_inv/body_win-'//TRIM(CBLOCK)//'-'//TRIM(CK)
      ENDIF
      PROJNM = ADJUSTL(PROJNM) 
      CALL WTSEISM(PROJNM, NDIM,NSAMP,NREC, NDIM,NSAMP,NREC,NSRC_WRK,  &
                   DT,STARTT, WIGGLE, IERR)
!
!.... clean 
      DEALLOCATE(WIGGLE) 
      DEALLOCATE(EST_LOC)  
      DEALLOCATE(FREQ_LOC) 
      DEALLOCATE(FREQ8) 
      DEALLOCATE(TWIN) 
      RETURN
      END
