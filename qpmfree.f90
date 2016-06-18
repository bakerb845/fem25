       SUBROUTINE QPMFREE(MASTER,MYID,MYNID,IPGROUP,NPGROUPS,NPPGRP,          & 
                          MGCOMM,MYSLV_COMM,MYHD_COMM,    PROJNM,             &
                          P1, WIN,INV,MSH,SRC,RCV,FRQ,M1D, MID,               &
                          P2SRF,P2BDY, IERR) 
!
!     Handles the matrix free matrix vector multiply via adjoint state 
!     of Q p.  This follows the papaer Full Waveform Inversion and the Truncated
!     Newton Method Method - L. Metivier, R. Brossier, J. Virieux, and S. Operto 
!     SIAM 2012 with a correction by B. Baker to remove the adj(J) J term
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'fwd_struc.h'
      TYPE (CMUMPS_STRUC) MID  !MUMPS 
      TYPE (INV_INFO)     INV  !inversion information 
      TYPE (MESH_INFO)    MSH  !mesh, model parameters
      TYPE (SRC_INFO)     SRC  !source information
      TYPE (RECV_INFO)    RCV  !receiver information
      TYPE (FRQ_INFO)     FRQ  !frequency information
      TYPE (MOD1D_INFO)   M1D  !1D model
      TYPE (WIN_INFO)     WIN  !window information 
      CHARACTER(*), INTENT(IN) :: PROJNM
      REAL*4, INTENT(IN) :: P1(*)
      INTEGER*4, INTENT(IN) :: MASTER,MYID,MYNID, MGCOMM,MYSLV_COMM, MYHD_COMM, &
                               NPGROUPS,NPPGRP, IPGROUP 
      REAL*4, INTENT(OUT) :: P2SRF(*), P2BDY(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FILENM
      COMPLEX*8, ALLOCATABLE :: ADJLAM(:), ADJNU(:), ALPHA(:), UE(:), WAVE(:)
      REAL*4, ALLOCATABLE :: V(:), PBUF1_SRF(:), PBUF2_SRF(:), PBUF3_SRF(:),  &
                             PBUF1_BDY(:), PBUF2_BDY(:), PBUF3_BDY(:),  &
                             P1BUF(:), P2BUF(:), P3BUF(:), PWORK(:)   
      INTEGER*4, ALLOCATABLE :: NSG_LSTB(:), NSG_LST(:), ISRCLST(:)
      COMPLEX*8 CZERO
      REAL*8 PYAVG 
      REAL*4 EPSS, THRESHS 
      INTEGER*4 STAT(MPI_STATUS_SIZE), MYDEST, MYTAG, MYSRC, NSG_GRP, & 
                JFREQ, IFREQ, IFREQL, ISG, JPGROUP, KSRC, JSRC, ISRC, NSGMAX, NOBS, &
                MPIERR
      LOGICAL*4 LSKIP !, LPSURF, LPBODY
      INTEGER*4 NOBSLST
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      PARAMETER(EPSS = 0.2)    !scales thresholding criteria in Huber norm, untested
      PARAMETER(THRESHS = 0.2) !scales thresholding criteria for L1/L2 norm, untested
      REAL*8, PARAMETER :: PI180 = 0.017453292519943295D0
!
!
!----------------------------------------------------------------------------------------!
!
!.... set space and initialize
      IERR = 0
      ALLOCATE(V(inv%NA35))              !complex search direction vector
      IF (MYID == MASTER) CALL SCOPY(inv%NA35,P1,1,V,1) 
      CALL MPI_BCAST(V,inv%NA35,MPI_REAL,MASTER,MGCOMM,MPIERR)
      print *, minval(v),maxval(v)
      ALLOCATE(UE(msh%NDOF))             !background solution
      ALLOCATE(WAVE(msh%NDOF))           !total wavefield 
      ALLOCATE(ADJLAM(msh%NDOF))         !solution of 2.24
      ALLOCATE(ADJNU (msh%NDOF))         !solution of 2.23; modified to remove adj(J) J
      ALLOCATE(ALPHA (msh%NDOF))         !solution of 2.15
      IF (MYNID == MASTER) THEN
         ALLOCATE(P3BUF(inv%NA35))
         ALLOCATE(P2BUF(inv%NA35))
         ALLOCATE(P1BUF(inv%NA35))
         ALLOCATE(PBUF1_SRF(inv%NA35)) 
         ALLOCATE(PBUF2_SRF(inv%NA35))
         ALLOCATE(PBUF3_SRF(inv%NA35)) 
         PBUF1_SRF(1:inv%NA35) = 0.0 
         PBUF2_SRF(1:inv%NA35) = 0.0
         PBUF3_SRF(1:inv%NA35) = 0.0
         ALLOCATE(PBUF1_BDY(inv%NA35))
         ALLOCATE(PBUF2_BDY(inv%NA35))
         ALLOCATE(PBUF3_BDY(inv%NA35))
         PBUF1_BDY(1:inv%NA35) = 0.0
         PBUF2_BDY(1:inv%NA35) = 0.0
         PBUF3_BDY(1:inv%NA35) = 0.0
      ELSE
         ALLOCATE(P3BUF(1))
         ALLOCATE(P2BUF(1))
         ALLOCATE(P1BUF(1)) 
         ALLOCATE(PBUF1_SRF(1)) 
         ALLOCATE(PBUF2_SRF(1))
         ALLOCATE(PBUF3_SRF(1)) 
         ALLOCATE(PBUF1_BDY(1))
         ALLOCATE(PBUF2_BDY(1))
         ALLOCATE(PBUF3_BDY(1))
         PBUF1_SRF(1) = 0.0
         PBUF2_SRF(1) = 0.0
         PBUF3_SRF(1) = 0.0
         PBUF1_BDY(1) = 0.0
         PBUF2_BDY(1) = 0.0
         PBUF3_BDY(1) = 0.0
      ENDIF
      IF (MYNID == MASTER) THEN
         ALLOCATE(NSG_LSTB(NPGROUPS))
         ALLOCATE(NSG_LST(NPGROUPS))
         ALLOCATE(ISRCLST(src%NSRC))
      ENDIF
      CALL MPI_BARRIER(MGCOMM,MPIERR)
!
!.... loop on frequencies
      !LPSURF = .FALSE. 
      !LPBODY = .FALSE.
      DO 1000 IFREQL=1,frq%NFREQ
!
!....... get frequency
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
            LSKIP = .TRUE.
            IF (MYID == MASTER) THEN 
               MYSRC = 0
               PYAVG = 0.D0 
               DO 201 JPGROUP=0,NPGROUPS-1
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1)
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN 
                     IF (frq%CFTYPE(JFREQ) == 'S') THEN !surface wave frequency
                        IF (frq%LINVF(JFREQ)) THEN 
                           LSKIP = .FALSE.
                           IF (JPGROUP > 0) THEN 
                              MYSRC = JPGROUP*NPPGRP
                              CALL MPI_RECV(PYAVG,1,MPI_DOUBLE_PRECISION,MYSRC, &
                                            MPI_ANY_TAG, MGCOMM,STAT,MPIERR)
                           ELSE 
                              PYAVG = src%PYAVG(ISG) 
                           ENDIF
                           !LPSURF = .TRUE.
                           WRITE(*,9409) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                        ENDIF !end check on inversion frequency
                     ELSE 
                        IF (frq%LINVF(JFREQ)) THEN
                           LSKIP = .FALSE.
                           IF (JPGROUP > 0) THEN
                              MYSRC = JPGROUP*NPPGRP
                              CALL MPI_RECV(PYAVG,1,MPI_DOUBLE_PRECISION,MYSRC, &
                                            MPI_ANY_TAG, MGCOMM,STAT,MPIERR)
                           ELSE
                              PYAVG = src%PYAVG(ISG)
                           ENDIF
                           !LPBODY = .TRUE.
                           WRITE(*,9411) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                        ENDIF !end check on inversion frequency
                     ENDIF !end check on frequency type 
                  ENDIF !end check on frequency range
  201          CONTINUE !loop on groups
            ELSE 
               IF (MYNID == MASTER) THEN
                  MYDEST = MASTER
                  MYTAG  = MYID
                  IF (IFREQ <= frq%NFREQ .AND. src%NSG > 0 .AND. ISG <= src%NSG) THEN
                     IF (frq%CFTYPE(IFREQ) == 'S') THEN !surface wave frequency
                        IF (frq%LINVF(IFREQ)) THEN
                           CALL MPI_SEND(src%PYAVG(ISG),1,MPI_DOUBLE_PRECISION, MASTER, &
                                         MYTAG, MGCOMM,MPIERR)
                        ENDIF !end check if using freq
                     ELSE !body wave frequency
                        IF (frq%LINVF(IFREQ)) THEN
                           CALL MPI_SEND(src%PYAVG(ISG),1,MPI_DOUBLE_PRECISION, MASTER, &
                                         MYTAG, MGCOMM,MPIERR)
                        ENDIF !end check if using freq 
                     ENDIF !end check on cftype
                  ENDIF !end check on range
               ENDIF !end check on mynid
            ENDIF
            CALL MPI_BCAST(LSKIP,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
            IF (LSKIP) THEN
               IF (MYID == MASTER) &
               WRITE(*,*) 'qpmfree: Skipping frequency group; no inversion frequencies...'
               GOTO 2001
            ENDIF
!
!.......... assemble and factor
            IF (MYID == MASTER) WRITE(*,*) 'qpmfree: Assembling impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN !check bounds
               IF (frq%LINVF(IFREQ)) THEN !only interested in inversion frequencies
                  CALL ASMBLE_DRIVER(msh%NZLOC,frq%CFTYPE(IFREQ),frq%FREQ(IFREQ), &
                                     src%PYAVG(ISG), MSH, mid%A_LOC,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'qpmfree: Error calling asmble_driver!'
                     GOTO 500
                  ENDIF
               ENDIF !end check on inversion frequency
            ENDIF !end check on bounds
            IF (MYID == MASTER) WRITE(*,*) 'qpmfree: Factoring impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN !check bounds
               IF (frq%LINVF(IFREQ)) THEN !only interested in inv. freqs
                  MID%JOB = 2
                  CALL CMUMPS(MID)
                  IF (MID%INFO(1) < 0) THEN
                     WRITE(*,*) 'qpmfree: An error occurred in the factorization!'
                     IERR = 1
                     GOTO 500
                  ENDIF
               ENDIF !end check on inversion frequency
            ENDIF !end check on bounds
!
!.......... have master tell which sources are being used
            IF (MYID == MASTER) THEN
               MYSRC = 0
               ISRCLST(1:src%NSRC) = 0
               DO 290 JPGROUP=0,NPGROUPS-1
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1)
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN
                     IF (frq%CFTYPE(JFREQ) == 'S') THEN
                        IF (frq%LINVF(JFREQ)) THEN
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
  292                      CONTINUE
  282                      CONTINUE
                        ENDIF !end check on working
                     ELSE !body wave frequency
                        IF (frq%LINVF(JFREQ)) THEN
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
                        IF (frq%LINVF(IFREQ)) THEN
                           KSRC = 0
                           DO 297 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
                              KSRC = KSRC + 1
                              ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  297                      CONTINUE
                           CALL MPI_SEND(ISRCLST,src%NSRC,MPI_INTEGER, MASTER,MYTAG, &
                                         MGCOMM,MPIERR)
                        ENDIF  !end check on working
                     ELSE !body wave
                        IF (frq%LINVF(IFREQ)) THEN
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
            IF (ISG > src%NSG)         GOTO 2001 !no solves, leave
            IF (IFREQ > frq%NFREQ)     GOTO 2001 !no solves, leave
            IF (.NOT.frq%LINVF(IFREQ)) GOTO 2001 !not an inversion frequency
!
!.......... loop on sources in source group
            DO 3000 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
               ISRC = src%ISRCPRM(JSRC) !extract original source ID 

               IF (MYNID == MASTER) THEN
                  NOBS = NOBSLST(NDIM, rcv%NREC,NDIM,                            &
                                 inv%OBS(1:NDIM,IFREQ,1:rcv%NREC,ISRC))
               ENDIF
               CALL MPI_BCAST(NOBS,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR)
               IF (NOBS == 0) GOTO 3005

!
!............. load the 1D solution
               IF (MYNID == MASTER) THEN
                  IF (MYID == MASTER) WRITE(*,9410) ISRC, src%AOI(ISRC), &
                                                    src%BAZN(ISRC), src%PYTAB(IFREQ,ISRC)
 9410             FORMAT(/,' qpmfree: Calculating ue field for source ',I3,/,          &
                           '          With angle of incidence ',F8.3,' degrees',/,     &
                           '          And adjusted azimuth ',F8.3,' degrees',/,        &
                           '          And y slowness ',G10.3,' s/m',/)
!
!................ this is where we convolve the source time function 
                  CALL LOAD_GRNS25(NDIM,msh%NDOF,msh%NNPE,NDIM, src%SRCTYP(ISRC),ISRC,   &
                                   frq%FREQ(IFREQ), src%SOURCE(IFREQ,ISRC),msh%IDOFSE,   &
                                   UE,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'qpmfree: Error calling load_grns25!'
                     GOTO 500
                  ENDIF
!
!................ do not want scaling in [dS/dm u] if inversion is phase only
!................ haskell will give me unit scaling; hence only applies to surface waves
!                 IF (inv%IRESTP == 1 .AND. frq%CFTYPE(IFREQ) == 'S') THEN
!                    DO IDOF=1,msh%NDOF
!                       IF (CABS(UE(IDOF)) > 0.0) UE(IDOF) = UE(IDOF)/CABS(UE(IDOF))
!                    ENDDO
!                 ENDIF
                  IF (MYID == MASTER) &
                  WRITE(*,*) 'qpmfree: Calculating effective forces...'
               ENDIF !end check on mynid
!
!............. calculate the force distribution
               CALL MPI_BCAST(UE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
               CALL ASMB_PEFF(MASTER,MYSLV_COMM, msh%NDOF,msh%NDOFL,msh%NZLOC, msh%CNP,  &
                              msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC, mid%A_LOC,UE, WAVE)
               IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,WAVE,1,mid%RHS,1)
               mid%ICNTL(9) = 1 !always solve Ax = b 
               MID%JOB = 3 !solve phase
               !IF (MYID == MASTER) WRITE(*,*) 'qpmfree: Solving Ax=b...'
               CALL CMUMPS(MID)
               IF (mid%INFO(1) < 0) THEN
                  WRITE(*,*) 'qpmfree: An error occurred in the solution!',MYID
                  IERR = 1
                  RETURN
               ENDIF
!
!............. dont extract response; because we already calculated them; but do get U  
               IF (MYNID == MASTER) THEN
                  CALL ADDBLK(msh%NDOF,msh%CNP,UE, mid%RHS) !add background field in
                  CALL CCOPY(msh%NDOF,mid%RHS,1,WAVE,1) 
               ENDIF !end check on myid
               IF (MYID == MASTER) WRITE(*,*) 'qpmfree: Calculating lambda field...'
               CALL CADJLAM(MASTER,MYNID,MYSLV_COMM, msh%NDOF, .FALSE.,IFREQ,ISRC, 0.0, &   
                            msh%AZMOD,EPSS,frq%FREQ(IFREQ), RCV,INV,MID, ADJLAM,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'qpmfree: Error calling cadjlam!'
                  GOTO 500
               ENDIF
               IF (MYID == MASTER) WRITE(*,*) 'qpmfree: Calculating nu and alpha field...'
               CALL MPI_BCAST(WAVE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
               CALL CADJ_NU_ALPHA(MASTER,MYNID,MYSLV_COMM, msh%NDOF,inv%NA35,            &
                                  frq%FREQ(IFREQ),src%PYTAB(IFREQ,ISRC),  V,ADJLAM,WAVE, &
                                  MSH,INV,MID,  P3BUF,ADJNU,ALPHA, IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'qpmfree: Error calling cadj_nu_alpha!'
                  GOTO 500
               ENDIF  
!              if (mynid == master) then
!                 call plot_adjgrad_vtk(ngnod,ndim,msh%nen, ndim,msh%nnpg,msh%ndof,   &
!                                       msh%nlxi,msh%nleta, msh%nelem, ndim, isrc,2,  &
!                                       frq%freq(ifreq), msh%lm,msh%ieng, msh%xlocs,msh%zlocs, adjnu)
!                 call plot_adjgrad_vtk(ngnod,ndim,msh%nen, ndim,msh%nnpg,msh%ndof,   &
!                                       msh%nlxi,msh%nleta, msh%nelem, ndim, isrc,3,  &
!                                       frq%freq(ifreq), msh%lm,msh%ieng, msh%xlocs,msh%zlocs, alpha)
!              endif
               IF (MYID == MASTER) WRITE(*,*) 'qpmfree: Calling met225_mod...'
               CALL MET225_MOD(MASTER,MYNID,MYSLV_COMM, msh%NDOF,                  &
                               frq%FREQ(IFREQ),src%PYTAB(IFREQ,ISRC),              &
                               ADJLAM,ADJNU,ALPHA,WAVE, MSH,INV, P1BUF,P2BUF, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'qpmfree: Error calling met225_mod!'
                  GOTO 500
               ENDIF
               IF (MYID == MASTER) WRITE(*,*) 'qpmfree: Updating pbuf...'
               IF (MYNID == MASTER) THEN
                  IF (frq%CFTYPE(IFREQ) == 'S') THEN
                     CALL SAXPY(inv%NA35,1.0,P1BUF,1,PBUF1_SRF,1)
                     CALL SAXPY(inv%NA35,1.0,P2BUF,1,PBUF2_SRF,1)
                     CALL SAXPY(inv%NA35,1.0,P3BUF,1,PBUF3_SRF,1)
                     print *, minval(pbuf1_srf),maxval(pbuf1_srf)
                     print *, minval(pbuf2_srf),maxval(pbuf2_srf)
                     print *, minval(pbuf3_srf),maxval(pbuf3_srf)
                  ELSE !body wave update
                     CALL SAXPY(inv%NA35,1.0,P1BUF,1,PBUF1_BDY,1)
                     CALL SAXPY(inv%NA35,1.0,P2BUF,1,PBUF2_BDY,1)
                     CALL SAXPY(inv%NA35,1.0,P3BUF,1,PBUF3_BDY,1)
                     print *, minval(pbuf1_bdy),maxval(pbuf1_bdy)
                     print *, minval(pbuf2_bdy),maxval(pbuf2_bdy)
                     print *, minval(pbuf3_bdy),maxval(pbuf3_bdy)
                  ENDIF
               ENDIF !end check onmynid
 3005          CONTINUE !nothing to do
 3000       CONTINUE !loop on sources in group 
 2001       CONTINUE !break ahead for loop on source group
 2000    CONTINUE !loop on sources in source group
         IF (src%NSG > 0 .AND. IFREQ <= frq%NFREQ) THEN
            IF (ASSOCIATED(src%ISGPTR))   DEALLOCATE(src%ISGPTR)
            IF (ASSOCIATED(src%ISRCPRM))  DEALLOCATE(src%ISRCPRM)
            IF (ASSOCIATED(src%PYAVG))    DEALLOCATE(src%PYAVG)
         ENDIF
         CALL MPI_BARRIER(MGCOMM,MPIERR)
 1000 CONTINUE !loop on frequencies in inversion 
 1005 CONTINUE !we are out of frequencies
  500 CONTINUE !break ahead for errors
      IF (IERR /= 0) THEN
         WRITE(*,*) 'qpmfree: An error was detected on process:',MYID
         RETURN
      ENDIF
!
!.... reduce onto head nodes
      IF (MYNID == MASTER) THEN
         !CALL MPI_BCAST(LPSURF,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR) 
         !CALL MPI_BCAST(LPBODY,1,MPI_LOGICAL, MASTER,MGCOMM,MPIERR)
         ALLOCATE(PWORK(inv%NA35)) 
         PWORK(:) = 0.0 
         !IF (LPSURF) THEN
         IF (inv%LSURF) THEN
            CALL MPI_REDUCE(PBUF1_SRF,PWORK,inv%NA35, MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
            IF (MYID == MASTER) THEN
               FILENM(1:80) = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_srf_t1'
               FILENM = ADJUSTL(FILENM)
               !P2(1:inv%NA35) = PWORK(1:inv%NA35) 
               CALL SCOPY(inv%NA35,PWORK,1,P2SRF,1) 
               CALL PLOT_SHGRAD_VTK(FILENM,NGNOD,msh%NNPG, msh%NELEM,          &
                                    inv%CINVTYPE, 5,2,1,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,            &
                                    msh%XLOCS,msh%ZLOCS, pwork)
            ENDIF
            CALL MPI_REDUCE(PBUF2_SRF,PWORK,inv%NA35, MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
            IF (MYID == MASTER) THEN
               FILENM(1:80) = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_srf_t2'
               FILENM = ADJUSTL(FILENM)
               !P2(1:inv%NA35) = P2(1:inv%NA35) + PWORK(1:inv%NA35) 
               CALL SAXPY(inv%NA35,1.0,PWORK,1,P2SRF,1)
               CALL PLOT_SHGRAD_VTK(FILENM,NGNOD,msh%NNPG, msh%NELEM,          &
                                    inv%CINVTYPE, 5,2,1,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,            &
                                    msh%XLOCS,msh%ZLOCS, pwork)
            ENDIF
            CALL MPI_REDUCE(PBUF3_SRF,PWORK,inv%NA35, MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
            IF (MYID == MASTER) THEN
               FILENM(1:80) = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_srf_t3'
               FILENM = ADJUSTL(FILENM)
               !P2(1:inv%NA35) = P2(1:inv%NA35) + PWORK(1:inv%NA35) 
               CALL SAXPY(inv%NA35,1.0,PWORK,1,P2SRF,1) 
               CALL PLOT_SHGRAD_VTK(FILENM,NGNOD,msh%NNPG, msh%NELEM,          &
                                    inv%CINVTYPE, 5,2,1,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,            &
                                    msh%XLOCS,msh%ZLOCS, pwork)
            ENDIF
         ENDIF
         !IF (LPBODY) THEN
         IF (inv%LBODY) THEN
            CALL MPI_REDUCE(PBUF1_BDY,PWORK,inv%NA35, MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
            IF (MYID == MASTER) THEN 
               FILENM(1:80) = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_bdy_t1'
               FILENM = ADJUSTL(FILENM)
               !P2(1:inv%NA35) = PWORK(1:inv%NA35) 
               CALL SCOPY(inv%NA35,PWORK,1,P2BDY,1) 
               CALL PLOT_SHGRAD_VTK(FILENM,NGNOD,msh%NNPG, msh%NELEM,          &
                                    inv%CINVTYPE, 5,2,1,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,            &
                                    msh%XLOCS,msh%ZLOCS, pwork)
            ENDIF
            CALL MPI_REDUCE(PBUF2_BDY,PWORK,inv%NA35, MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
            IF (MYID == MASTER) THEN 
               FILENM(1:80) = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_bdy_t2'
               FILENM = ADJUSTL(FILENM)
               !P2(1:inv%NA35) = P2(1:inv%NA35) + PWORK(1:inv%NA35) 
               CALL SAXPY(inv%NA35,1.0,PWORK,1,P2BDY,1)
               CALL PLOT_SHGRAD_VTK(FILENM,NGNOD,msh%NNPG, msh%NELEM,          &
                                    inv%CINVTYPE, 5,2,1,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,            &
                                    msh%XLOCS,msh%ZLOCS, pwork)
            ENDIF
            CALL MPI_REDUCE(PBUF3_BDY,PWORK,inv%NA35, MPI_REAL,MPI_SUM, &
                            MASTER,MYHD_COMM,MPIERR)
            IF (MYID == MASTER) THEN 
               FILENM(1:80) = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_bdy_t3'
               FILENM = ADJUSTL(FILENM)
               !P2(1:inv%NA35) = P2(1:inv%NA35) + PWORK(1:inv%NA35) 
               CALL SAXPY(inv%NA35,1.0,PWORK,1,P2BDY,1) 
               CALL PLOT_SHGRAD_VTK(FILENM,NGNOD,msh%NNPG, msh%NELEM,          &
                                    inv%CINVTYPE, 5,2,1,0,                     &
                                    msh%CDOMAIN,inv%MASKG,msh%IENG,            &
                                    msh%XLOCS,msh%ZLOCS, pwork)
            ENDIF

         ENDIF
         IF (ALLOCATED(PWORK)) DEALLOCATE(PWORK) 
      ENDIF
!
!.... clean space
      IF (ALLOCATED(V))            DEALLOCATE(V) 
      IF (ALLOCATED(UE))           DEALLOCATE(UE)
      IF (ALLOCATED(WAVE))         DEALLOCATE(WAVE)
      IF (ALLOCATED(ADJLAM))       DEALLOCATE(ADJLAM)
      IF (ALLOCATED(ADJNU))        DEALLOCATE(ADJNU)
      IF (ALLOCATED(ALPHA))        DEALLOCATE(ALPHA) 
      IF (ALLOCATED(PBUF1_SRF))    DEALLOCATE(PBUF1_SRF) 
      IF (ALLOCATED(PBUF2_SRF))    DEALLOCATE(PBUF2_SRF)
      IF (ALLOCATED(PBUF3_SRF))    DEALLOCATE(PBUF3_SRF) 
      IF (ALLOCATED(PBUF1_BDY))    DEALLOCATE(PBUF1_BDY)
      IF (ALLOCATED(PBUF2_BDY))    DEALLOCATE(PBUF2_BDY)
      IF (ALLOCATED(PBUF3_BDY))    DEALLOCATE(PBUF3_BDY)
      IF (ALLOCATED(P1BUF))        DEALLOCATE(P1BUF) 
      IF (ALLOCATED(P2BUF))        DEALLOCATE(P2BUF)
      IF (ALLOCATED(P3BUF))        DEALLOCATE(P3BUF)
      IF (ALLOCATED(NSG_LSTB))     DEALLOCATE(NSG_LSTB)
      IF (ALLOCATED(NSG_LST))      DEALLOCATE(NSG_LST)
      IF (ALLOCATED(ISRCLST))      DEALLOCATE(ISRCLST)
!
!.... format statements
 9409 FORMAT(' ----------------------------------------------------------',/,&
            ' - qpmfree: Group:',I4,'                                 -',/,&
            ' -          Propagating surface wave source group ',I3,'  -',/,&
            ' -          With average slowness in y',E12.4,'     - ',/,&
            ' -          At frequency',F12.5,' Hz                -',/,&
            ' ----------------------------------------------------------',/)
 9411 FORMAT(' ----------------------------------------------------------',/,&
             ' - qpmfree: Group:',I4,'                                 -',/,& 
             ' -          Propagating body wave source group ',I5,'   -',/,& 
             ' -          With average slowness in y',E12.4,'     - ',/,& 
             ' -          At frequency',F12.5,' Hz                -',/,& 
             ' ----------------------------------------------------------',/)
 9421 FORMAT(/,' qpmfree: Will process surface wave source,',I3,/   &
               '          With angle of incidence ',F8.3,' degrees',/, &
               '          Adjusted azimuth ',F8.3,' degrees',/,    &
               '          And y slowness ',G10.3,' s/m',/)
 9423 FORMAT(/,' qpmfree: Will process body wave source,',I3,/   &
               '          With angle of incidence ',F8.3,' degrees',/, &
               '          Adjusted azimuth ',F8.3,' degrees',/,    &
               '          And y slowness ',G10.3,' s/m',/)

      RETURN
      END

      REAL*8 FUNCTION ELEM_AREA(MGNOD,NNPG,NGNOD, IELEM,IENG, XLOCS,ZLOCS)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG) 
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MGNOD,NNPG,NGNOD, IELEM 
      REAL*8 RLOC, X1, X2, Z1, Z2 
      INTEGER*4 IA1, IA2, JNPG1, JNPG2, IMOD  
      RLOC = 0.D0 
      IMOD = 5
      IF (MINVAL(IENG(1:4,IELEM)) == 0) IMOD = 4
      DO 1 IA1=1,NGNOD  
         JNPG1 = IENG(IA1,IELEM)
         IF (JNPG1 == 0) GOTO 30 !break for triangle 
         IA2 = MOD(IA1+1,IMOD)
         IF (IA2 == 0) IA2 = 1
         JNPG2 = IENG(IA2,IELEM)
         X1 = XLOCS(JNPG1)
         X2 = XLOCS(JNPG2)
         Z1 = ZLOCS(JNPG1)
         Z2 = ZLOCS(JNPG2)
         RLOC = RLOC + X1*Z2 - X2*Z1
    1 CONTINUE  
   30 CONTINUE !break ahead for triangle
      IF (RLOC <= 0.D0) THEN 
         WRITE(*,*) 'elem_wts: Warning your mesh is wrong'
         RLOC =-RLOC
      ENDIF 
      RLOC = RLOC*0.5D0 !only going one way here
      ELEM_AREA = RLOC
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE CADJ_NU_ALPHA(MASTER,MYNID,MYCOMM, NDOF,NA35, &
                               FREQ,PY, V,ADJLAM,WAVE,         &
                               MSH,INV,MID,  P3,NU,ALPHA, IERR) 
!
!     Calculates the nu adjoint wavefield -[dS/dm1 v1 + ... + dS/dmm vm]^dagger lambda
!     Note that we perform the transpose matrix vector multiply by mirroring sgemv.f
!     whereupon 
!        y_i <- y_i + a_{ij} x_j is instead performed by
!        y_j <- y_j + a_{ij} x_i 
! 
!     And, we can calculate the alpha wavefield -[dS/dm1 v1 + ... + dS/dmm vm] u
!     since we already have dS/dm1 and the wavefield u by this point
!
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'fwd_struc.h'
      TYPE (CMUMPS_STRUC) MID !MUMPS 
      TYPE (INV_INFO)     INV !inversion information
      TYPE (MESH_INFO)    MSH !mesh information
      COMPLEX*8, INTENT(IN) :: ADJLAM(NDOF), WAVE(NDOF)
      REAL*8, INTENT(IN) :: FREQ, PY
      REAL*4, INTENT(IN) :: V(NA35)  
      INTEGER*4, INTENT(IN) :: MASTER,MYNID,MYCOMM, NDOF, NA35 
      COMPLEX*8, INTENT(OUT) :: NU(NDOF), ALPHA(NDOF)  
      REAL*4, INTENT(OUT) :: P3(*) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*16, ALLOCATABLE :: DKDA(:,:), DKDB(:,:), DKDA2(:,:), DKDB2(:,:)
      COMPLEX*8, ALLOCATABLE :: TEMP(:,:), RHS_DENSE(:), RHSBNU(:), RHSALP(:), CARG(:),  &
                                CARG2(:), RHSBUF(:)    
      REAL*8, ALLOCATABLE :: SHL(:,:,:,:), SHG(:,:,:,:), DAB(:,:,:), DAB2(:,:,:),     &
                             XIGLL(:,:), ETAGLL(:,:), DMDA(:,:), DMDB(:,:), DMDA2(:,:), &
                             DMDB2(:,:), DET(:,:) 
      REAL*4, ALLOCATABLE :: P3BUF(:) 
      COMPLEX*8 CWGHT, CZERO
      REAL*8 OMEGA, TWOPI 
      INTEGER*4 NWORK, NINTX, NINTZ, NRHSS, NEE, NNZROW, ISAVE9, ISAVE20, INPGL, INPG, &
                ICON, IZERO, INPINV, INZROW, IVINV, I1,I2, IAE,IBE,I,J, IQE,IPE,       &
                IDOF,JDOF, IELEM, IA35, MPIERR   
      LOGICAL*4 LISP, LISS
      INTEGER*4 IBSECT
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      PARAMETER(TWOPI = 6.2831853071795862D0)
       real*8 elem_area, area
      !complex*8 cdotc 
!
!----------------------------------------------------------------------------------------!
! 
!.... first error check
      IERR = 0
      IF (inv%CINVTYPE /= 'PP' .AND. inv%CINVTYPE /= 'pp' .AND. &
          inv%CINVTYPE /= 'SS' .AND. inv%CINVTYPE /= 'ss' .AND. &
          inv%CINVTYPE /= 'PS' .AND. inv%CINVTYPE /= 'ps') THEN
         WRITE(*,*) 'cadj_nu_alpha2: Error anisotropy not programmed!'
         IERR = 1 
         RETURN
      ENDIF
!
!.... set inversion type
      LISP = .FALSE.
      LISS = .FALSE.
      IF (inv%CINVTYPE == 'PP' .OR. inv%CINVTYPE == 'pp') LISP = .TRUE.
      IF (inv%CINVTYPE == 'SS' .OR. inv%CINVTYPE == 'ss') LISS = .TRUE.
      IF (inv%CINVTYPE == 'PS' .OR. inv%CINVTYPE == 'ps') THEN
         LISP = .TRUE.
         LISS = .TRUE.
      ENDIF
!
!.... save RHSs
      IF (MYNID == MASTER) THEN
         NRHSS   = mid%NRHS      !number of RHSs
         ISAVE9  = mid%ICNTL(9)  !save transpose problem 
         ISAVE20 = mid%ICNTL(20) !sparse RHS?
         IF (ISAVE20 == 0) THEN 
            NWORK = NRHSS*MID%LRHS
            ALLOCATE(RHS_DENSE(NWORK))
            !RHS_DENSE(1:NWORK) = mid%RHS(1:NWORK)
            CALL CCOPY(NWORK,mid%RHS,1,RHS_DENSE,1)
            DEALLOCATE(mid%RHS)
            ALLOCATE(mid%RHS(mid%N))
         ENDIF
         mid%ICNTL(20) = 0 !dense RHS
         mid%NRHS = 1      !only one RHS
      ENDIF
!
!.... initialize element shape fns and zero out gradient buffer
      NINTX = 2*msh%NORD + 2         !Lobatto quad, k = 2n - 3, k->2k for mass
      NINTX = NINTX + 3          !add a little extra for deformed elements
      NINTZ = NINTX
      ALLOCATE(SHL(3,msh%NEN,NINTX,NINTZ))
      ALLOCATE(XIGLL (NINTX,2))
      ALLOCATE(ETAGLL(NINTZ,2))
      CALL DSHL(msh%NEN,NINTX,NINTZ, msh%IITYPE,NINTX,NINTZ, msh%NLXI,msh%NLETA,    &   
                msh%XIPTS,msh%ETAPTS, XIGLL,ETAGLL,SHL)
      OMEGA = TWOPI*FREQ
      NEE = NDIM*msh%NEN
      ALLOCATE(SHG(3,msh%NEN,NINTX,NINTZ)) 
      ALLOCATE(DAB2(NINTX,NINTZ,2))
      ALLOCATE(DAB (NINTX,NINTZ,2))
      ALLOCATE(DET(NINTX,NINTZ))
      ALLOCATE(DMDA(NEE,NEE))
      ALLOCATE(DKDA(NEE,NEE))
      ALLOCATE(DMDB(NEE,NEE))
      ALLOCATE(DKDB(NEE,NEE))
      ALLOCATE(DMDA2(NEE,NEE))
      ALLOCATE(DKDA2(NEE,NEE)) 
      ALLOCATE(DMDB2(NEE,NEE))
      ALLOCATE(DKDB2(NEE,NEE))  
      ALLOCATE(RHSBNU(msh%NDOF))
      ALLOCATE(RHSALP(msh%NDOF)) 
      RHSBNU(1:msh%NDOF) = CZERO 
      RHSALP(1:msh%NDOF) = CZERO
      ALLOCATE(TEMP(inv%MBUFRHS,inv%NVINV))
      TEMP(1:inv%MBUFRHS,1:inv%NVINV) = CMPLX(0.0,0.0)
!     ALLOCATE(TEMP(NDOF,inv%NVINV)) 
      ALLOCATE(CARG(inv%NA35))
      ALLOCATE(CARG2(inv%NA35)) 
      ALLOCATE(P3BUF(inv%NA35)) 
      P3BUF(1:inv%NA35) = 0.0 
!
!.... loop on local Hessian points 
      DO 100 INPGL=1,inv%NNPGL,inv%NVINV 
!
!....... locate anchor node to invert at
         DO 101 INPG=1,msh%NNPG
            INPINV = inv%MASKG(INPG)
            IF ((INPINV - 1)*inv%NVINV + 1 == inv%MYGRAD(INPGL)) GOTO 110 
  101    CONTINUE  
         WRITE(*,*) 'cadj_nu_alpha: Cant find anchor node!'
         IERR = 1 
         RETURN 
  110    CONTINUE 
!
!....... loop on connectivity 
         DO 200 ICON=1,inv%NCON
            IELEM = inv%MCONN(INPG,ICON) 
            IF (IELEM <= 0) GOTO 250 !done with connections
            !IF (msh%CDOMAIN(IELEM) /= 'I') GOTO 200 !only use inversion elements
!
!.......... integrate element
            IF (msh%LISISO) THEN !P and S gradients 
               CALL DJACAB2(msh%NEN,NINTX,NINTZ,msh%NNPG, msh%LISISO,             &
                            NINTX,NINTZ,msh%NEN,NGNOD,msh%NNPG, INPG,             &
                            msh%IENG(1:NGNOD,IELEM),                              &
                            XIGLL(:,1),ETAGLL(:,1),msh%XLOCS,msh%ZLOCS,           &
                            msh%DENS,msh%ECOEFF,SHL,                              &
                            DET,DAB,DAB2,SHG,IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) ''
                  RETURN
               ENDIF
!
!............. numerically integrate element matrices
               !CALL DABMKISOB25(msh%NEN,NEE,NINTX,NINTZ, msh%NEN,NINTX,NINTZ,   & 
               !                 LISP,LISS,OMEGA,PY,                             &
               !                 XIGLL(:,2),ETAGLL(:,2),               &
               !                 DET,DAB, SHG, DMDA, DMDB, DKDA, DKDB)  
               !print *, 'a',minval(cdabs(dkda)),minval(cdabs(dkdb)),maxval(cdabs(dkda)),maxval(cdabs(dkdb))
               !print *, 'a',minval(dmda),minval(dmdb),maxval(dmda),maxval(dmdb)
               !CALL DABMKISOB25(msh%NEN,NEE,NINTX,NINTZ, msh%NEN,NINTX,NINTZ,   &
               !                 LISP,LISS,OMEGA,PY,                   &
               !                 XIGLL(:,2),ETAGLL(:,2),               &
               !                 DET,DAB2,SHG, DMDA2,DMDB2,DKDA2,DKDB2) 
               !print *, 'a',minval(cdabs(dkda2)),minval(cdabs(dkdb2)),maxval(cdabs(dkda2)),maxval(cdabs(dkdb2))
               !print *, 'a',minval(dmda2),minval(dmdb2),maxval(dmda2),maxval(dmdb2)
               CALL DABM2KISOB25(msh%NEN,NEE,NINTX,NINTZ, msh%NEN,NINTX,NINTZ,   &
                                 LISP,LISS,OMEGA,PY,                   &
                                 XIGLL(:,2),ETAGLL(:,2),               &
                                 DET,DAB,DAB2, SHG, DMDA, DMDB, DKDA, DKDB, &
                                 DMDA2,DMDB2,DKDA2,DKDB2)
               !print *, minval(dmda),minval(dmdb), maxval(dmda),maxval(dmdb)
               !print *, 'b',minval(cdabs(dkda)),minval(cdabs(dkdb)),maxval(cdabs(dkda)),maxval(cdabs(dkdb))
               !print *, 'b',minval(cdabs(dkda2)),minval(cdabs(dkdb2)),maxval(cdabs(dkda2)),maxval(cdabs(dkdb2))
               !print *, 'b',minval(dmda),minval(dmdb), maxval(dmda),maxval(dmdb)
               !print *, 'b',minval(dmda2),minval(dmdb2), maxval(dmda2),maxval(dmdb2)
            ELSE  
               WRITE(*,*) 'cadj_nu_alpha: No anisotropy yet!'
               IERR = 1 
               RETURN 
            ENDIF 
!
!.......... null out these column vector of 2nd derivative matrix
            I1 = inv%JCSC_FDIST(INPGL)
            I2 = inv%JCSC_FDIST(INPGL+1) - 1 
            NNZROW = I2 - I1 + 1 
            DO 150 IVINV=1,inv%NVINV
               TEMP(1:NNZROW,IVINV) = CZERO
  150       CONTINUE
!           TEMP(1:NDOF,1:inv%NVINV) = CZERO
!           AREA = ELEM_AREA(NGNOD,msh%NNPG,NGNOD, IELEM,msh%IENG, msh%XLOCS,msh%ZLOCS)

!
!.......... perform the scaled transpose matrix vector multiply; loop on columns 
            IPE = 0
            DO 201 IAE=1,msh%NEN
               DO 202 I=1,NDIM
                  IPE = IPE + 1
                  IDOF = msh%LM(I,IAE,IELEM)
                  IF (IDOF <= 0) GOTO 220 !not a dof
                  INZROW = IBSECT(NNZROW,IDOF,inv%ICSC_FDIST(I1:I2))
                  IF (INZROW <= 0 .OR. INZROW > NNZROW) THEN 
                     WRITE(*,*) 'cadj_nu_alpha: Error cant find inzrow!'
                     IERR = 1
                     RETURN
                  ENDIF
!
!................ loop on columns of S to perform column vector multiply
                  IQE = 0
                  DO 203 IBE=1,msh%NEN
                     DO 204 J=1,NDIM
                        IQE = IQE + 1
                        JDOF = msh%LM(J,IBE,IELEM)
                        IF (JDOF <= 0) GOTO 240 !not a dof 
                        CARG(1:inv%NVINV) = CZERO
                        IF (msh%LISISO) THEN
                           IVINV = 0
                           IF (LISP) THEN !dS/dalpha v(alpha) 
                              IVINV = IVINV + 1
                              CARG(IVINV) = CMPLX(DCMPLX(DMDA (IPE,IQE),0.D0)            &
                                                       + DKDA (IPE,IQE))
                              CARG2(IVINV)= CMPLX(DCMPLX(DMDA2(IPE,IQE),0.D0)            &
                                                       + DKDA2(IPE,IQE))
                           ENDIF
                           IF (LISS) THEN !dS/dbeta v(beta)
                              IVINV = IVINV + 1
                              CARG(IVINV) = CMPLX(DCMPLX(DMDB (IPE,IQE),0.D0)            &
                                                       + DKDB (IPE,IQE))
                              CARG2(IVINV)= CMPLX(DCMPLX(DMDB2(IPE,IQE),0.D0)            &
                                                       + DKDB2(IPE,IQE))
                                                      
                           ENDIF
                        ELSE
                           WRITE(*,*) 'cadj_nu_alpha: No anisotropy!'
                           IERR = 1
                           RETURN
                        ENDIF
!----------------------------------------------------------------------------------------!
!        We calculate three terms here:                                                  !
!          (1) nu:              -[v_i dS/dm_i]^dagger lambda                             !
!          (2) alpha:           -[v_i dS/dm_i] u                                         !
!          (3) 2nd derivatives: e_j [v_j d^S/dm_j^2 u]                                   ! 
!----------------------------------------------------------------------------------------!
                        DO 205 IVINV=1,inv%NVINV
                           IA35 = (INPINV - 1)*inv%NVINV + IVINV 
                           CWGHT = CMPLX(SNGL(inv%ELEM_WTS(IELEM))*V(IA35),0.0)
                           !CWGHT = CMPLX(V(IA35)*SNGL(1.D0/AREA),0.0)
                           RHSBNU(JDOF) = RHSBNU(JDOF) &
                                        + ADJLAM(IDOF)*CONJG(CWGHT*CARG(IVINV))
                           RHSALP(IDOF) = RHSALP(IDOF) &
                                        + CWGHT*CARG(IVINV)*WAVE(JDOF)
                           !CWGHT = CMPLX(V(IA35)*SNGL(1.D0/(AREA**2)),0.0) 
                           !CWGHT = CMPLX(SNGL(inv%ELEM_WTS(IELEM))**2*V(IA35),0.0)
                           CWGHT = CMPLX(SNGL(inv%ELEM_WTS(IELEM))*V(IA35),0.0)
                           TEMP(INZROW,IVINV) = TEMP(INZROW,IVINV) &
                                              + CWGHT*CARG(IVINV)*WAVE(JDOF)
!                          TEMP(IDOF,IVINV) = TEMP(IDOF,IVINV) &
!                                           + CWGHT*CARG(IVINV)*WAVE(JDOF) 
  205                   CONTINUE !loop on inversion variables 
  204                CONTINUE !loop on components
  240                CONTINUE !not a DOF
  203             CONTINUE !loop on element nodes 
  220             CONTINUE !not a DOF
  202          CONTINUE !loop on components
  201       CONTINUE !loop on element nodes
  200    CONTINUE !loop on connectivity 
  250    CONTINUE !break ahead, out of connections
!
!....... finish the vector vector multiply for 2nd derivatives
         DO 210 IVINV=1,inv%NVINV
            IA35 = (INPINV - 1)*inv%NVINV + IVINV 
!           P3BUF(IA35) = P3BUF(IA35) - REAL( CDOTC(NDOF,TEMP(1:NDOF,IVINV),1,ADJLAM,1) )
            I1 = inv%JCSC_FDIST(INPGL+IVINV-1)
            I2 = inv%JCSC_FDIST(INPGL+IVINV) - 1 
            INZROW = 0
            DO 211 IZERO=I1,I2 
               INZROW = INZROW + 1
               IDOF = inv%ICSC_FDIST(IZERO)
               P3BUF(IA35) = P3BUF(IA35) - REAL( CONJG(TEMP(INZROW,IVINV))*ADJLAM(IDOF) )
  211       CONTINUE
  210    CONTINUE
  100 CONTINUE !loop on inversion variables
!
!.... sum contributions for nu 
      IF (MYNID == MASTER) THEN
         ALLOCATE(RHSBUF(msh%NDOF)) 
      ELSE  
         ALLOCATE(RHSBUF(1)) 
      ENDIF
      CALL MPI_REDUCE(RHSBNU,RHSBUF,msh%NDOF, MPI_COMPLEX,MPI_SUM, MASTER,MYCOMM,MPIERR)
!
!.... solution phase for S^T nu* =-[dS/dm1 v1 + ... + dS/dm_m vm]^T lam^*
      IF (MYNID == MASTER) THEN
         DO 405 IDOF=1,msh%NDOF
            mid%RHS(IDOF) = CONJG(RHSBUF(IDOF)) 
  405    CONTINUE 
      ENDIF
!
!.... solution phase for S^T nu* =-[dS/dm1 v1 + ... + dS/dm_m vm]^T lam^* 
      IF (MYNID == MASTER) mid%ICNTL(9)  = 0 !want to solve transpose problem
      mid%JOB = 3
      CALL CMUMPS(MID)
      IF (MID%INFO(1) < 0) THEN
         IF (MYNID == MASTER) WRITE(*,*) 'cadj_nu_alpha: Error in solution phase 1'
         IERR = 1
         RETURN
      ENDIF
!
!.... extract nu and send to all processes 
      IF (MYNID == MASTER) THEN 
         DO 500 IDOF=1,msh%NDOF
            NU(IDOF) = CONJG(mid%RHS(IDOF))
  500    CONTINUE
      ENDIF 
      CALL MPI_BCAST(NU,msh%NDOF,MPI_COMPLEX, MASTER,MYCOMM,MPIERR)
!
!.... sum contributions for alpha 
      CALL MPI_REDUCE(RHSALP,RHSBUF,msh%NDOF, MPI_COMPLEX,MPI_SUM, MASTER,MYCOMM,MPIERR)
      IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,RHSBUF,1,mid%RHS,1)
      DEALLOCATE(RHSBUF) 
!
!.... solution phase for S alpha =-[dS/dm1 v1 + ... + dS/dm_m vm] u
      IF (MYNID == MASTER) mid%ICNTL(9)  = 1 !want to solve regular problem
      mid%JOB = 3
      CALL CMUMPS(MID)
      IF (MID%INFO(1) < 0) THEN 
         IF (MYNID == MASTER) WRITE(*,*) 'cadj_nu_alpha: Error in solution phase 2'
         IERR = 1
         RETURN
      ENDIF
!
!.... extract alpha and send to all processes 
      IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,mid%RHS,1,ALPHA,1)
      CALL MPI_BCAST(ALPHA,msh%NDOF,MPI_COMPLEX, MASTER,MYCOMM,MPIERR) 
! 
!.... restore memory 
      IF (MYNID == MASTER) THEN
         mid%NRHS = NRHSS
         mid%ICNTL(9)  = ISAVE9
         mid%ICNTL(20) = ISAVE20
         DEALLOCATE(mid%RHS)
         IF (MID%ICNTL(20) == 0) THEN
            ALLOCATE(mid%RHS(NWORK))
            !mid%RHS(1:NWORK) = RHS_DENSE(1:NWORK)
            CALL CCOPY(NWORK,RHS_DENSE,1,mid%RHS,1) 
            DEALLOCATE(RHS_DENSE)
         ENDIF
      ENDIF
!
!.... reduce second derivatives on head node
      CALL MPI_REDUCE(P3BUF,P3,inv%NA35, MPI_REAL,MPI_SUM, MASTER,MYCOMM,MPIERR)
!
!.... free space
      DEALLOCATE(SHL)
      DEALLOCATE(XIGLL)
      DEALLOCATE(ETAGLL) 
      DEALLOCATE(SHG) 
      DEALLOCATE(DAB2)
      DEALLOCATE(DAB)
      DEALLOCATE(DET)
      DEALLOCATE(DMDA)
      DEALLOCATE(DKDA)
      DEALLOCATE(DMDB)
      DEALLOCATE(DKDB)
      DEALLOCATE(DMDA2)
      DEALLOCATE(DKDA2) 
      DEALLOCATE(DMDB2)
      DEALLOCATE(DKDB2)  
      DEALLOCATE(TEMP)
      DEALLOCATE(CARG)
      DEALLOCATE(CARG2) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE MET225_MOD(MASTER,MYNID,MYCOMM, NDOF,  &
                            FREQ,PY,                    &
                            ADJLAM,NU,ALPHA,WAVE,  MSH,INV, P1BUF,P2BUF, IERR) 
!
!     Modified matrix vector multiply of Metivier's equation 2.25 where I've removed 
!     the adj(J) J v contribution
!
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO)     INV !inversion information
      TYPE (MESH_INFO)    MSH !mesh information
      COMPLEX*8, INTENT(IN) :: ADJLAM(NDOF), ALPHA(NDOF), NU(NDOF), WAVE(NDOF)
      REAL*8, INTENT(IN) :: FREQ, PY
      INTEGER*4, INTENT(IN) :: MASTER,MYNID,MYCOMM, NDOF
      REAL*4, INTENT(OUT) :: P1BUF(*), P2BUF(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*16, ALLOCATABLE :: DKDA(:,:), DKDB(:,:)
      COMPLEX*8, ALLOCATABLE :: T1(:,:), T2(:,:), CARG(:)
      REAL*8, ALLOCATABLE :: SHL(:,:,:,:), SHG(:,:,:,:), DAB(:,:,:),                   &
                             XIGLL(:,:), ETAGLL(:,:), DMDA(:,:), DMDB(:,:), DET(:,:)  
      REAL*4, ALLOCATABLE :: PBUF1(:), PBUF2(:), PWORK(:) 
      COMPLEX*8 CWGHT, CZERO
      REAL*8 OMEGA, TWOPI
      INTEGER*4 NINTX, NINTZ, NEE, NNZROW, INPGL, INPG, &
                ICON, IZERO, INPINV, INZROW, IVINV, I1,I2, IAE,IBE,I,J, IQE,IPE,       &
                IDOF,JDOF, IELEM, IA35, MPIERR
      LOGICAL*4 LISP, LISS
      INTEGER*4 IBSECT
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      PARAMETER(TWOPI = 6.2831853071795862D0)
      real*8 area, elem_area
      !complex*8 cdotc 
!
!----------------------------------------------------------------------------------------!
!
!.... first error check
      IERR = 0
      IF (inv%CINVTYPE /= 'PP' .AND. inv%CINVTYPE /= 'pp' .AND. &
          inv%CINVTYPE /= 'SS' .AND. inv%CINVTYPE /= 'ss' .AND. &
          inv%CINVTYPE /= 'PS' .AND. inv%CINVTYPE /= 'ps') THEN
         WRITE(*,*) 'met225_mod: Error anisotropy not programmed!'
         IERR = 1
         RETURN
      ENDIF
!
!.... set inversion type
      LISP = .FALSE.
      LISS = .FALSE.
      IF (inv%CINVTYPE == 'PP' .OR. inv%CINVTYPE == 'pp') LISP = .TRUE.
      IF (inv%CINVTYPE == 'SS' .OR. inv%CINVTYPE == 'ss') LISS = .TRUE.
      IF (inv%CINVTYPE == 'PS' .OR. inv%CINVTYPE == 'ps') THEN
         LISP = .TRUE.
         LISS = .TRUE.
      ENDIF
!
!.... initialize element shape fns and zero out gradient buffer
      NINTX = 2*msh%NORD + 2         !Lobatto quad, k = 2n - 3, k->2k for mass
      NINTX = NINTX + 3          !add a little extra for deformed elements
      NINTZ = NINTX
      ALLOCATE(SHL(3,msh%NEN,NINTX,NINTZ))
      ALLOCATE(XIGLL (NINTX,2))
      ALLOCATE(ETAGLL(NINTZ,2))
      CALL DSHL(msh%NEN,NINTX,NINTZ, msh%IITYPE,NINTX,NINTZ, msh%NLXI,msh%NLETA,    &   
                msh%XIPTS,msh%ETAPTS, XIGLL,ETAGLL,SHL)
      OMEGA = TWOPI*FREQ
      NEE = NDIM*msh%NEN
      ALLOCATE(SHG(3,msh%NEN,NINTX,NINTZ)) 
      ALLOCATE(DAB (NINTX,NINTZ,2))
      ALLOCATE(DET(NINTX,NINTZ))
      ALLOCATE(DMDA(NEE,NEE))
      ALLOCATE(DKDA(NEE,NEE))
      ALLOCATE(DMDB(NEE,NEE))
      ALLOCATE(DKDB(NEE,NEE))
      ALLOCATE(PBUF1(inv%NA35))
      ALLOCATE(PBUF2(inv%NA35))
      PBUF1(1:inv%NA35) = 0.0 
      PBUF2(1:inv%NA35) = 0.0 
      ALLOCATE(CARG(inv%NVINV)) 
      ALLOCATE(T1(inv%MBUFRHS,inv%NVINV))
      ALLOCATE(T2(Inv%MBUFRHS,inv%NVINV)) 
!     ALLOCATE(T1(NDOF,inv%NVINV))
!     ALLOCATE(T2(NDOF,inv%NVINV)) 
!
!.... loop on local Hessian points 
      DO 100 INPGL=1,inv%NNPGL,inv%NVINV
!
!....... locate anchor node to invert at
         DO 101 INPG=1,msh%NNPG
            INPINV = inv%MASKG(INPG)
            IF ((INPINV - 1)*inv%NVINV + 1 == inv%MYGRAD(INPGL)) GOTO 110
  101    CONTINUE
         WRITE(*,*) 'met225_mod: Cant find anchor node!'
         IERR = 1
         RETURN
  110    CONTINUE
!
!....... loop on connectivity 
         DO 200 ICON=1,inv%NCON
            IELEM = inv%MCONN(INPG,ICON)
            IF (IELEM <= 0) GOTO 250 !done with connections
            CWGHT = CMPLX(SNGL(inv%ELEM_WTS(IELEM)),0.0)
            AREA = ELEM_AREA(NGNOD,msh%NNPG,NGNOD, IELEM,msh%IENG, msh%XLOCS,msh%ZLOCS)
            !CWGHT = CMPLX(SNGL(1.D0/AREA),0.0) 
!
!.......... integrate element
            IF (msh%LISISO) THEN !P and S gradients 
               CALL DJACAB(msh%NEN,NINTX,NINTZ,msh%NNPG, msh%LISISO,             &
                           NINTX,NINTZ,msh%NEN,NGNOD,msh%NNPG, INPG,             &
                           msh%IENG(1:NGNOD,IELEM),                              &
                           XIGLL(:,1),ETAGLL(:,1),msh%XLOCS,msh%ZLOCS,           &
                           msh%DENS,msh%ECOEFF,SHL,                              &
                           DET,DAB,SHG,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) ''
                  RETURN
               ENDIF
!
!............. numerically integrate element matrices
               CALL DABMKISOB25(msh%NEN,NEE,NINTX,NINTZ, msh%NEN,NINTX,NINTZ,   &
                                LISP,LISS,OMEGA,PY,                             &
                                XIGLL(:,2),ETAGLL(:,2),               &
                                DET,DAB, SHG, DMDA, DMDB, DKDA, DKDB)
            ELSE
               WRITE(*,*) 'met225_mod: No anisotropy yet!'
               IERR = 1
               RETURN
            ENDIF
!
!.......... null out these column vector of 2nd derivative matrix
            I1 = inv%JCSC_FDIST(INPGL)
            I2 = inv%JCSC_FDIST(INPGL+1) - 1
            NNZROW = I2 - I1 + 1
            DO 150 IVINV=1,inv%NVINV
               T1(1:NNZROW,IVINV) = CZERO
               T2(1:NNZROW,IVINV) = CZERO
  150       CONTINUE
!           T1(1:NDOF,1:inv%NVINV) = CZERO
!           T2(1:NDOF,1:inv%NVINV) = CZERO

!
!.......... perform the scaled transpose matrix vector multiply; loop on columns 
            IPE = 0
            DO 201 IAE=1,msh%NEN
               DO 202 I=1,NDIM
                  IPE = IPE + 1
                  IDOF = msh%LM(I,IAE,IELEM)
                  IF (IDOF <= 0) GOTO 220 !not a dof
                  INZROW = IBSECT(NNZROW,IDOF,inv%ICSC_FDIST(I1:I2))
                  IF (INZROW <= 0 .OR. INZROW > NNZROW) THEN
                     WRITE(*,*) 'met225_mod: Error cant find inzrow!'
                     IERR = 1
                     RETURN
                  ENDIF
!
!................ loop on columns of S to perform column vector multiply
                  IQE = 0
                  DO 203 IBE=1,msh%NEN
                     DO 204 J=1,NDIM
                        IQE = IQE + 1
                        JDOF = msh%LM(J,IBE,IELEM)
                        IF (JDOF <= 0) GOTO 240 !not a dof 
                        CARG(1:inv%NVINV) = CZERO
                        IF (msh%LISISO) THEN
                           IVINV = 0
                           IF (LISP) THEN !dS/dalpha v(alpha) 
                              IVINV = IVINV + 1
                              CARG(IVINV) = CMPLX(DCMPLX(DMDA(IPE,IQE),0.D0)            &
                                                       + DKDA(IPE,IQE))
                           ENDIF
                           IF (LISS) THEN !dS/dbeta v(beta)
                              IVINV = IVINV + 1
                              CARG(IVINV) = CMPLX(DCMPLX(DMDB(IPE,IQE),0.D0)            &
                                                       + DKDB(IPE,IQE))
                           ENDIF
                        ELSE
                           WRITE(*,*) 'met225_mod: No anisotropy!'
                           IERR = 1
                           RETURN
                        ENDIF
                        DO 205 IVINV=1,inv%NVINV
                           T1(INZROW,IVINV) = T1(INZROW,IVINV) &
                                            + CWGHT*CARG(IVINV)*WAVE(JDOF)
                           T2(INZROW,IVINV) = T2(INZROW,IVINV) &
                                            + CWGHT*CARG(IVINV)*ALPHA(JDOF)
!                          T1(IDOF,IVINV) = T1(IDOF,IVINV) &
!                                         + CWGHT*CARG(IVINV)*WAVE(JDOF)
!                          T2(IDOF,IVINV) = T2(IDOF,IVINV) &
!                                         + CWGHT*CARG(IVINV)*ALPHA(JDOF) 
  205                   CONTINUE !loop on inversion variables
  204                CONTINUE !loop on components
  240                CONTINUE !not a DOF
  203             CONTINUE !loop on element nodes 
  220             CONTINUE !not a DOF
  202          CONTINUE !loop on components
  201       CONTINUE !loop on element nodes
  200    CONTINUE !loop on connectivity 
  250    CONTINUE !break ahead, out of connections
!
!....... finish the vector vector multiply 
         DO 210 IVINV=1,inv%NVINV
            IA35 = (INPINV - 1)*inv%NVINV + IVINV 
!           PBUF1(IA35) = PBUF1(IA35) &
!                       - REAL(CDOTC(NDOF,T1(1:NDOF,IVINV),1,NU,1))
!           PBUF2(IA35) = PBUF2(IA35) &
!                       - REAL(CDOTC(NDOF,T2(1:NDOF,IVINV),1,ADJLAM,1))
            I1 = inv%JCSC_FDIST(INPGL+IVINV-1)
            I2 = inv%JCSC_FDIST(INPGL+IVINV) - 1
            INZROW = 0
            DO 211 IZERO=I1,I2
               INZROW = INZROW + 1
               IDOF = inv%ICSC_FDIST(IZERO)
               PBUF1(IA35) = PBUF1(IA35) - REAL(CONJG(T1(INZROW,IVINV))*NU(IDOF))
               PBUF2(IA35) = PBUF2(IA35) - REAL(CONJG(T2(INZROW,IVINV))*ADJLAM(IDOF)) 
  211       CONTINUE
  210    CONTINUE
  100 CONTINUE !loop on inversion variables
!
!.... reduce onto buffer
      IF (MYNID == MASTER) THEN
         ALLOCATE(PWORK(inv%NA35))
      ELSE
         ALLOCATE(PWORK(1))
      ENDIF
      CALL MPI_REDUCE(PBUF1,PWORK,inv%NA35, MPI_REAL,MPI_SUM, MASTER,MYCOMM,MPIERR) 
      IF (MYNID == MASTER) CALL SCOPY(inv%NA35,PWORK,1,P1BUF,1)
      CALL MPI_REDUCE(PBUF2,PWORK,inv%NA35, MPI_REAL,MPI_SUM, MASTER,MYCOMM,MPIERR)
      IF (MYNID == MASTER) CALL SCOPY(inv%NA35,PWORK,1,P2BUF,1)
          
!
!.... clean space
      DEALLOCATE(PWORK) 
      DEALLOCATE(SHL)
      DEALLOCATE(XIGLL)
      DEALLOCATE(ETAGLL)
      DEALLOCATE(SHG)
      DEALLOCATE(DAB)
      DEALLOCATE(DET)
      DEALLOCATE(DMDA)
      DEALLOCATE(DKDA)
      DEALLOCATE(DMDB)
      DEALLOCATE(DKDB)
      DEALLOCATE(T1)
      DEALLOCATE(T2) 
      DEALLOCATE(CARG)
      DEALLOCATE(PBUF1)
      DEALLOCATE(PBUF2) 
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!     SUBROUTINE MATV_QV( ) 
!
!     Performs the matrix vector multiply Qv which B. Baker's modified form for 
!     Metivier's equation 2.25

!
!----------------------------------------------------------------------------------------!
! 
!.... first error check
!     IERR = 0
!     IF (inv%CINVTYPE /= 'PP' .AND. inv%CINVTYPE /= 'pp' .AND. &
!         inv%CINVTYPE /= 'SS' .AND. inv%CINVTYPE /= 'ss' .AND. &
!         inv%CINVTYPE /= 'PS' .AND. inv%CINVTYPE /= 'ps') THEN
!        WRITE(*,*) 'cadjnu: Error anisotropy not programmed!'
!        IERR = 1
!        RETURN
!     ENDIF
!
!.... save RHSs
!     IF (MYNID == MASTER) THEN
!        NRHSS   = mid%NRHS      !number of RHSs
!        ISAVE9  = mid%ICNTL(9)  !save transpose problem 
!        ISAVE20 = mid%ICNTL(20) !sparse RHS?
!        IF (ISAVE20 == 0) THEN
!           NWORK = NRHSS*MID%LRHS
!           ALLOCATE(RHS_DENSE(NWORK))
!           RHS_DENSE(1:NWORK) = mid%RHS(1:NWORK)
!           DEALLOCATE(mid%RHS)
!           ALLOCATE(mid%RHS(mid%N))
!        ENDIF
!        mid%ICNTL(20) = 0 !dense RHS
!        mid%NRHS = 1      !only one RHS
!     ENDIF
!
!.... initialize element shape fns and zero out gradient buffer
!     NINTX = 2*msh%NORD + 2         !Lobatto quad, k = 2n - 3, k->2k for mass
!     NINTX = NINTX + 3          !add a little extra for deformed elements
!     NINTZ = NINTX
!     ALLOCATE(SHL(3,msh%NEN,NINTX,NINTZ))
!     ALLOCATE(XIGLL (NINTX,2))
!     ALLOCATE(ETAGLL(NINTZ,2))
!     CALL DSHL(msh%NEN,NINTX,NINTZ, msh%IITYPE,NINTX,NINTZ, msh%NLXI,msh%NLETA,    &
!               msh%XIPTS,msh%ETAPTS, XIGLL,ETAGLL,SHL)
!     OMEGA = TWOPI*FREQ
!     NEE = NDIM*msh%NEN
!
!.... loop on local Hessian points 
!     DO 100 INPGL=1,inv%NNPGL,inv%NVINV
!
!....... locate anchor node to invert at
!        DO 101 INPG=1,msh%NNPG
!           INPINV = inv%MASKG(INPG)
!           IF ((INPINV - 1)*inv%NVINV + 1 == inv%MYGRAD(INPGL)) GOTO 110
! 101    CONTINUE
!        WRITE(*,*) 'cadjnu: Cant find anchor node!'
!        IERR = 1
!        GOTO 730
! 110    CONTINUE
!
!....... loop on connectivity 
!        DO 200 ICON=1,NCON
!           IELEM = inv%MCONN(INPG,ICON)
!           IF (IELEM <= 0) GOTO 250 !done with connections
!
!.......... integrate element
!           IF (msh%LISISO) THEN !P and S gradients 
!
!............. evaluate model derivatives at integration points
!              CALL DJACAB(msh%NEN,NINTX,NINTZ,msh%NNPG, msh%LISISO,                  &
!                          NINTX,NINTZ,msh%NEN,NGNOD,msh%NNPG, INPG,                  &
!                          msh%IENG(1:NGNOD,IELEM),                                   &
!                          XIGLL(1:NINTX,1),ETAGLL(1:NINTZ,1),                        &
!                          msh%XLOCS,msh%ZLOCS,msh%DENS,msh%ECOEFF(1:msh%NNPG,1:2),   &
!                          SHL(1:3,1:NEN,1:NINTX,1:NINTZ),                            &
!                          DET,DAB,SHG,IERR)
!              IF (IERR /= 0) THEN
!                 WRITE(*,*) 'cadjnu: Error calling djacab'
!                 GOTO 730
!              ENDIF
!
!............. numerically integrate element matrices
!              CALL DABMKISOB25(msh%NEN,NEE,NINTX,NINTZ, msh%NEN,NINTX,NINTZ,   &
!                               OMEGA,PY, XIGLL(1:NINTX,2),ETAGLL(1:NINTZ,2),   &
!                               DET,DAB,SHG, DMDA,DMDB,DKDA,DKDB)
!           ELSE
!              WRITE(*,*) 'cadjnu: No anisotropy yet!'
!              IERR = 1
!              GOTO 730
!           ENDIF
!
!.......... perform the scaled transpose matrix vector multiply; loop on columns 
!           DSMU(1:NEE) = CZERO
!           DSMA(1:NEE) = CZERO
!           IDOFS(1:NEE) = 0
!           IPE = 0
!           DO 201 IAE=1,NEN
!              DO 202 I=1,NDIM
!                 IPE = IPE + 1
!                 IDOF = LM(I,IAE,IELEM)
!                 IF (IDOF <= 0) GOTO 220 !not a dof
!                 IDOFS(IPE) = IDOF
!
!................ loop on columns of S to perform column vector multiply
!                 IQE = 0
!                 DO 203 IBE=1,NEN
!                    DO 204 J=1,NDIM
!                       IQE = IQE + 1
!                       JDOF = LM(J,IBE,IELEM)
!                       IF (JDOF <= 0) GOTO 240 !not a dof 
!                       CARG(1:NVINV) = CZERO
!                       IF (LISISO) THEN
!                          IVINV = 0
!                          IF (LISP) THEN !dS/dalphau
!                             IVINV = IVINV + 1
!                             CARG(IVINV) = CMPLX(DCMPLX(DMDA(IPE,IQE),0.D0)  &
!                                                      + DKDA(IPE,IQE)
!                             CARG2(IVINV)= CMPLX(DCMPLX(DMDA2(IPE,IQE),0.D0) &
!                                                      + DKDA2(IPE,IQE)
!                          ENDIF
!                          IF (LISS) THEN !dS/dbeta)
!                             IVINV = IVINV + 1
!                             CARG(IVINV) = CMPLX(DCMPLX(DMDB(IPE,IQE),0.D0)
!                                                      + DKDB(IPE,IQE))
!                             CARG2(IVINV)= CMPLX(DCMPLX(DMDB2(IPE,IQE),0.D0) &
!                                                      + DKDB2(IPE,IQE))  
!                          ENDIF
!                       ELSE
!                          WRITE(*,*) 'cadjnu: No anisotropy!'
!                          IERR = 1
!                          RETURN
!                       ENDIF
!
!...................... update dS/dm u and dS/dm alpha 
!                       DO 206 IVINV=1,NVINV
!                          DSMU(IPE)  = DSMU(IPE)  + CARG(IVINV)*U(JDOF) 
!                          DSMA(IPE)  = DSMA(IPE)  + CARG(IVINV)*ALPHA(JDOF) 
!                          DSMU2(IPE) = DSMU2(IPE) + CARG2(IVINV)*U(JDOF)  
! 206                   CONTINUE
! 240                   CONTINUE !not a DOF
! 204                CONTINUE !loop on components
! 203             CONTINUE !loop on element nodes 
! 220             CONTINUE !not a DOF
! 202          CONTINUE !loop on components
! 201       CONTINUE !loop on element nodes
!
!.......... multiply i'th column of [dS/dm_i u] by nu, [dS/dm_i alpha] lambda, and 
!.......... v_i [d^2 S/ dm_i^2 u]^dagger lambda
!           IPE = 0 
!           DO 210 IAE=1,NEN
!              DO 211 I=1,NDIM
!                 IPE = IPE + 1 
!                 IDOF = IDOFC(IPE) 
!                 IF (IDOF > 0) THEN
!                    DO 212 IVINV=1,NVINV 
!                       IA35 = (INPINV - 1)*inv%NVINV + IVINV
!                       QVB(IA35) = QVB(IA35) + CONJG(DSMU(IPE))*NU(IDOF)
!                       QVB(IA35) = QVB(IA35) + CONJG(DSMA(IPE))*ADJLAM(IDOF)
!                       QVB(IA35) = QVB(IA35) + CONJG(DSMU2(IPE))*ADJLAM(IDOF)*V(IA35) 
! 212                CONTINUE
!                 ENDIF
! 211          CONTINUE !loop on components
! 210       CONTINUE !loop on elements nodes
! 200    CONTINUE !loop on connectivity 
! 110    CONTINUE !node not in my partition
! 100 CONTINUE !loop on inversion variables
!
!.... free space
!     DEALLOCATE(DSMU)  
!     DEALLOCATE(DSMA)
!     DEALLOCATE(DSMU2) 
!     RETURN
!     END
