      IMPLICIT NONE
!.... MUMPS
      INCLUDE 'cmumps_struc.h'
      TYPE (CMUMPS_STRUC) MID
!.... MPI
      INCLUDE 'mpif.h'
      REAL*8 TSSIM,TESIM
      INTEGER*4 MASTER,MYID,MYNID, MYSLV_COMM,MYHD_COMM,IPGROUP,JPGROUP,MPIERR,     &
                NPGROUPS,NPROCS,NPARTS, STAT(MPI_STATUS_SIZE), MYTAG, MYDEST, MYSRC 
      PARAMETER(MASTER = 0)
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO)  MSH  !mesh, model parameters 
      TYPE (MOD1D_INFO) M1D  !1D models
      TYPE (SRC_INFO)   SRC  !source information
      TYPE (RECV_INFO)  RCV  !receiver information
      TYPE (FRQ_INFO)   FRQ  !frequency information 
!.... local varaiables
      CHARACTER(80) PROJNM, FILENM, TMPDIR
      COMPLEX*8, ALLOCATABLE :: EST(:,:,:,:), EST_LOC(:,:,:,:), UE(:), WAVE(:) 
      REAL*8, ALLOCATABLE :: XREC(:), ZREC(:), FREQ_LOC(:) 
      INTEGER*4, ALLOCATABLE :: NSG_LST(:), NSG_LSTB(:), ISRCLST(:) 
      LOGICAL*4, ALLOCATABLE :: LFSURF(:)
      REAL*8 PYAVG
      INTEGER*4 NABS, NELEME, NSGMAX, IFREQ, JFREQ, IFREQL, ISRC, JSRC, KSRC,  &
                ISG, IREC, I, ICGRNS, NSG_GRP, IERR  
      LOGICAL*4 LFILES, LNSRF, LEX, LISDIR, LRECEX, LCGRNS 
!.... functions
      INTEGER*4 ICELEML
!
!----------------------------------------------------------------------------------------!
!
!.... initialize MPI 
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, MPIERR)
!----------------------------------------------------------------------------------------!
!                    Read information on mesh and model                                  !
!----------------------------------------------------------------------------------------!
!
!.... have head node read model information
      lfiles = .true.
      LCGRNS = .TRUE. !.false.   !always calculate 1D greens functions in forward modeling
      TSSIM = MPI_WTIME() 
      IERR = 0
      IF (MYID == MASTER) THEN
         WRITE(*,9407) 
 9407    FORMAT(' ---------------------------------------------------------------',/, &
                ' -   xbielak25: A massively parallel program for simulating    -',/, &
                ' -              plane wave and surface wave propragation in    -',/, &
                ' -              locally hetergeneous media                     -',/, &
                ' ---------------------------------------------------------------',/) 
         PROJNM(1:80) = ' '
         WRITE(*,*) 'xbielak25: Enter project name:'
         READ(*,'(A)') PROJNM
         !projnm = 'surf'
         !projnm = 'half'
         !!projnm = 'simple'
         !projnm = 'tester'
         !!projnm = 'tsmigr'
         PROJNM = ADJUSTL(PROJNM)         
         WRITE(*,*) 'xbielak25: Enter number of process groups:'
         READ *, NPGROUPS
         WRITE(*,*) 'xbielak25: Enter 1 to load the Greens functions from disk'
         WRITE(*,*) '           Any other number will calculate Greens functions'
         READ *, ICGRNS
         IF (ICGRNS == 1) THEN
            LCGRNS = .FALSE.
         ELSE
            LCGRNS = .TRUE.
         ENDIF
         !npgroups = 1
         !npgroups = 1
!
!....... get variables from the .ini files
         CALL READ_FWD_INI(PROJNM,.TRUE., TMPDIR,LFILES,LNSRF, MSH,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error I cannot locate your spec file!'
            GOTO 500
         ENDIF
         !msh%freq0 = 20.d0
         !msh%aztol = 10.d0
         !msh%aoitol = 5.d0
         !lfiles = .true.
!....... initial warnings on NPGROUPS
         IF (MOD(NPROCS,NPGROUPS) /= 0) THEN 
            WRITE(*,*) 'xbielak25: Error I cant divide nprocs by npgroups evenly!'
            IERR = 1
            GOTO 500
         ENDIF
         NPARTS = NPROCS/NPGROUPS
         WRITE(*,*) 'xbielak25: Reading mesh...'
         CALL RDMESHBK(PROJNM, NELEME,NABS, MSH, IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error reading mesh files!'
            GOTO 500 
         ENDIF 

         WRITE(*,*) 'xbielak25: Setting interpolation points...'
         ALLOCATE(msh%XIPTS (msh%NLXI))
         ALLOCATE(msh%ETAPTS(msh%NLETA)) 
         CALL GENPTS(msh%IITYPE, msh%NLXI,msh%NLETA, msh%XIPTS,msh%ETAPTS) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error reading mesh files!'
            GOTO 500 
         ENDIF 
         !print *, 'hacking pmls'
         !do inpg=1,nnpg
         !   xtemp = xd(inpg)
         !   ztemp = zd(inpg)
         !   if (xtemp > 0.d0 .and. ztemp == 0.d0) zd(inpg) = xtemp!/2.d0
         !   if (ztemp > 0.d0 .and. xtemp == 0.d0) xd(inpg) = ztemp!/2.d0
         !enddo 
         CALL PLOT_ELMOD_VTK(PROJNM,NGNOD,msh%NNPG,msh%NNPG, msh%NELEM, msh%LISISO, &
                             msh%IENG,msh%XLOCS,msh%ZLOCS, msh%DENS,msh%ECOEFF)
!
!....... read the frequency list
         CALL RDFREQ(PROJNM, FRQ,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error reading frequency list'
            GOTO 500 
         ENDIF
!
!....... read the source list
         CALL RDSRC_EQ(PROJNM, msh%XLATMIN,msh%XLONMIN,msh%XLATMAX,msh%XLONMAX, & 
                       msh%AZMOD, SRC,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'xbielak25: Error calling rdsrc_eq'
            GOTO 500
         ENDIF
!
!....... check for a receiver file
         CALL RDREC_HD(PROJNM, rcv%NREC,IERR)
         IF (IERR /= 0) THEN
            IF (IERR > 0) THEN
               WRITE(*,*) 'xbielak25: Error reading recv header'
               GOTO 500
            ELSE
               WRITE(*,*) 'xbielak25: Warning there is no receiver file!'
               IERR = 0 
            ENDIF 
            ALLOCATE(LFSURF(1)) 
            ALLOCATE(XREC(1))
            ALLOCATE(rcv%YREC(1))
            ALLOCATE(ZREC(1))
         ELSE !read the receiver file
            ALLOCATE(LFSURF  (rcv%NREC))
            ALLOCATE(    XREC(rcv%NREC))
            ALLOCATE(rcv%YREC(rcv%NREC))
            ALLOCATE(    ZREC(rcv%NREC))
            CALL RDREC(PROJNM, rcv%NREC, LFSURF,XREC,rcv%YREC,ZREC, IERR)
         ENDIF
!
!....... check for a source file
         ALLOCATE(src%SOURCE(frq%NFREQ,src%NSRC))
         CALL SRCSUB_SB(PROJNM,frq%NFREQ, frq%NFREQ,src%NSRC,                   &
                        frq%NFREQ_SRF,frq%NFREQ_BDY,src%NSRC_SRF,src%NSRC_BDY,  &
                        frq%CFTYPE,src%SRCTYP, frq%FREQ, src%SOURCE,IERR)
         IF (IERR /= 0) THEN
!           CALL RDSTF(PROJNM, frq%NFREQ, frq%NFREQ,src%NSRC, frq%FREQ,src%SOURCE, IERR)
            IF (IERR < 0) THEN
               WRITE(*,*) 'xbielak25: Warning Greens functions are unity'
               WRITE(*,*) '           You can convolve an appropriate STF later'
               IERR = 0
            ENDIF
         ENDIF
!
!....... read the frequency response list
         ALLOCATE(rcv%RECV(NDIM,frq%NFREQ,rcv%NREC))
         CALL RDREC_RESP(PROJNM, NDIM,frq%NFREQ,NDIM,frq%NFREQ,rcv%NREC, msh%AZMOD, &
                         frq%FREQ, LRECEX,rcv%RECV)
         IF (.NOT.LRECEX) THEN
            WRITE(*,*) 'xbielak25: Receiver response functions are set to unity' 
         ENDIF
!
!....... generate the 1D models
         IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Generating 1D models...'
         CALL FILL1D(MSH,M1D) 
!
!....... make directories (well load balanced programs will see an error otherwise)
         IF (.NOT.LFILES) THEN
            LEX = LISDIR('./blk_vtk')
            IF (.NOT.LEX) CALL SYSTEM('mkdir ./blk_vtk') !bielak response
            LEX = LISDIR('./frc_vtk')
            IF (.NOT.LEX) CALL SYSTEM('mkdir ./frc_vtk') !virtual forces on boundary
            LEX = LISDIR('./wav_vtk')
            IF (.NOT.LEX) CALL SYSTEM('mkdir ./wav_vtk') !total reponse
         ENDIF 
         LEX = LISDIR('./wave')
         IF (.NOT.LEX) CALL SYSTEM('mkdir ./wave')
         LEX = LISDIR('./dest')
         IF (.NOT.LEX) CALL SYSTEM('mkdir ./dest')
!
!....... split the process groups
         WRITE(*,*) 'xbielak25: Splitting process groups...'
      ENDIF 
      CALL MKMPI_GROUPS(MASTER,MYID,NPGROUPS,NPROCS,MPI_COMM_WORLD, &
                        IPGROUP,MYNID,MYHD_COMM,MYSLV_COMM)
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Initializing MUMPS...'
      mid%COMM = MYSLV_COMM !set communicator
      mid%SYM = 0 !unsymmetric
      mid%JOB =-1 !initialize
      mid%PAR = 1 !host host working 
      CALL CMUMPS(MID) 
      mid%ICNTL(3) = 0 !suppress output
      mid%ICNTL(4) = 1 !only error messages
!----------------------------------------------------------------------------------------!
!                        Graph generation                                                !
!----------------------------------------------------------------------------------------!
!
!.... generate a graph
      IF (MYNID == MASTER) THEN
         CALL MPI_BCAST(LFILES,1,MPI_LOGICAL, MASTER,MYHD_COMM,MPIERR) 
         IF (MYID == MASTER) THEN
            WRITE(*,*) 'xbielak25: Generating graph...'
            CALL GEN_GRAPH25(.TRUE.,NPARTS, LFSURF,XREC,ZREC, MSH,RCV,IERR) 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xbielak25: Error generating graph!'
               GOTO 500
            ENDIF
            mid%N  = msh%NDOF
            mid%NZ = msh%NZERO
            CALL PLOT_MESHP_VTK(PROJNM,NGNOD,NDIM,msh%NEN, msh%NNPG,msh%NELEM,msh%NDOF, &
                                NDIM,msh%NEN, msh%PART,msh%IENG,msh%LM, &
                                msh%XLOCS,msh%ZLOCS) 
            !CALL DUMP_CSRFULL(PROJNM, msh%ndof,msh%nzero,msh%PART,msh%IRPTR,msh%JCPTR)
            WRITE(*,*)
            WRITE(*,9408) msh%NORD,msh%NELEM,NABS,NELEME, &
                          msh%NNPG,msh%NNPE,msh%NDOF,msh%NZERO  
 9408       FORMAT(' xbielak25: Polynomial order:'                    ,I4 ,/,         &
                   '            Number of elements in mesh:'          ,I10,/,         &
                   '            Number of absorbing elements:'        ,I8 ,/,         &
                   '            Number of Bielak elements:'           ,I8 ,/,         &
                   '            Number of anchor nodes in mesh:'      ,I10,/,         & 
                   '            Number of nodes in Bielak boundary:'  ,I10,/,         &
                   '            Number of degrees of freedom:'        ,I14,/,         &
                   '            Number of non-zeros in global matrix:',I16,/)
!
!.......... free some unwatend variables 
            IF (ALLOCATED(LFSURF)) DEALLOCATE(LFSURF) 
            IF (ALLOCATED(XREC))   DEALLOCATE(XREC)
            IF (ALLOCATED(ZREC))   DEALLOCATE(ZREC)
         ENDIF
         mid%ICNTL(18) = 3 !i will tell MUMPS how to split the mesh
      ENDIF
!----------------------------------------------------------------------------------------!
!                          Broadcast Forward Problem Parameters                          !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Broadcasting 2D models...'
      CALL MPI_BCAST(NPARTS,1,MPI_INTEGER, MASTER,MPI_COMM_WORLD,MPIERR) 
      CALL BCAST_MESH_INFO(MYID,MPI_COMM_WORLD,MASTER, &
                           LNSRF, PROJNM,TMPDIR, MSH) 
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Broadcasting 1D models...'
      CALL BCAST_MOD1D_INFO(MYID,MPI_COMM_WORLD,MASTER, MSH,M1D) 
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Broadcasting source details...' 
      CALL BCAST_SRC_INFO(MYID,MPI_COMM_WORLD,MASTER, SRC) 
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Broadcasting frequency information...'
      CALL BCAST_FRQ_INFO(MYID,MPI_COMM_WORLD,MASTER, FRQ)
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Broadcasting receiver information...'
         CALL BCAST_RCV_INFO(MYID,MYHD_COMM,MASTER, RCV) 
         CALL BCAST_RCV_RESP(MYID,MYHD_COMM,MASTER, .TRUE.,frq%NFREQ, RCV)
         IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Broadcasting source time function...'
         CALL BCAST_SRC_STF(MYID,MYHD_COMM,MASTER,  .TRUE.,frq%NFREQ, SRC) 
      ENDIF
!----------------------------------------------------------------------------------------!
!               Calculate 1D analytic solutions and free superfluous space               !
!----------------------------------------------------------------------------------------!
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Calculating 1D solutions...'
      ALLOCATE(src%PYTAB(frq%NFREQ,src%NSRC))
      src%PYTAB(1:frq%NFREQ,1:src%NSRC) = 0.D0
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR) 
      !CALL GEN_GRNS(MYID,MASTER,MPI_COMM_WORLD, TMPDIR,LCGRNS,LNSRF, &
      !              SRC,M1D,FRQ,MSH, IERR) 
      CALL GEN_GRNS_SB(MYID,MASTER,MPI_COMM_WORLD, TMPDIR,LCGRNS,LNSRF,  &
                       SRC,M1D,FRQ,MSH,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'xbielak25: Error calling gen_grns_sb on process:',MYID
         GOTO 500
      ENDIF
!.... source
      IF (MYNID /= MASTER) THEN
         IF (ASSOCIATED(src%SRCTYP))  DEALLOCATE(src%SRCTYP)
         IF (ASSOCIATED(src%BAZN))    DEALLOCATE(src%BAZN)
         IF (ASSOCIATED(src%AOI))     DEALLOCATE(src%AOI)
      ENDIF
      IF (ASSOCIATED(src%SLAT))    DEALLOCATE(src%SLAT)
      IF (ASSOCIATED(src%SLON))    DEALLOCATE(src%SLON)
      IF (ASSOCIATED(src%SDEP))    DEALLOCATE(src%SDEP)
      IF (ASSOCIATED(src%STRIKE))  DEALLOCATE(src%STRIKE)
      IF (ASSOCIATED(src%DIP))     DEALLOCATE(src%DIP)
      IF (ASSOCIATED(src%RAKE))    DEALLOCATE(src%RAKE)
      IF (ASSOCIATED(src%SMAG))    DEALLOCATE(src%SMAG)
      IF (ASSOCIATED(src%MODE))    DEALLOCATE(src%MODE)
!----------------------------------------------------------------------------------------!
!                 Graph distribution.  Note we do not reorder the matrix as MUMPS        !
!                 calls METIS_NodeWND which seems to be a bit more clever than           !
!                 METIS_NodeND when the matrix is distributed                            !
!----------------------------------------------------------------------------------------!
!
!.... have masters save matrix sizes 
      IF (MYNID == MASTER) THEN
         mid%N = MSH%NDOF
         mid%NZ = MSH%NZERO
         mid%ICNTL(18) = 3 !matrix distributed by user 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
!.... count number of local non-zeros
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Generating local graphs...'
      CALL ICNDZ_LOC(MYNID,1,msh%NDOF, msh%PART,msh%IRPTR, msh%NDOFL,msh%NZLOC)
!.... generate distributed CSR structure
      ALLOCATE(msh%MYDOFS(msh%NDOFL))
      ALLOCATE(msh%IRPTR_LOC(msh%NDOFL+1))
      ALLOCATE(msh%JCPTR_LOC(msh%NZLOC))
      CALL PART2CRSL(MSH%NZERO,MSH%NZLOC, MSH%NDOF,MSH%NDOFL, MYNID,1,             &
                     MSH%IRPTR,MSH%JCPTR,MSH%PART, MSH%MYDOFS,MSH%IRPTR_LOC,MSH%JCPTR_LOC)
!.... counter number of local elements
      msh%NELEML = ICELEML(MYNID, NDIM,msh%NEN,msh%NDOF, msh%NELEM, 1,             &
                           msh%NEN,NDIM, msh%PART,msh%LM)
!.... have process ave its elements 
      ALLOCATE(msh%MYELEM(msh%NELEML))
      CALL GENELEML(MYNID, NDIM,msh%NEN,msh%NDOF, msh%NELEM,msh%NELEML, 1,         &
                    msh%NEN,NDIM, msh%PART,MSH%LM, msh%MYELEM)
      IF (ASSOCIATED(MSH%IRPTR)) DEALLOCATE(MSH%IRPTR)
      IF (ASSOCIATED(MSH%JCPTR)) DEALLOCATE(MSH%JCPTR)
      IF (ASSOCIATED(MSH%PART))  DEALLOCATE(MSH%PART)
!.... create the distributed COO storage for MUMPS
      mid%NZ_LOC = MSH%NZLOC
      ALLOCATE(mid%IRN_LOC(mid%NZ_LOC))
      ALLOCATE(mid%JCN_LOC(mid%NZ_LOC)) 
      CALL CRS2COOLOC(MSH%NZLOC,MSH%NDOFL, MSH%MYDOFS,MSH%IRPTR_LOC,MSH%JCPTR_LOC,      &
                      mid%IRN_LOC,mid%JCN_LOC)
!.... reorder the matrix
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Reordering matrix...'
      mid%JOB = 1
      CALL CMUMPS(MID) 
      ALLOCATE(mid%A_LOC(mid%NZ_LOC))
      IF (MYNID == MASTER) THEN
         mid%LRHS = mid%N
         ALLOCATE(mid%RHS(mid%N))
         ALLOCATE(EST(NDIM,frq%NFREQ,rcv%NREC,src%NSRC))
         EST(1:NDIM,1:frq%NFREQ,1:rcv%NREC,1:src%NSRC) = CMPLX(0.0,0.0) 
      ENDIF
      IF (MYID == MASTER) WRITE(*,*) 
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
!----------------------------------------------------------------------------------------!
!                 Loop on frequencies                                                    !
!----------------------------------------------------------------------------------------!
      ALLOCATE(UE(msh%NDOF))
      IF (MYNID == MASTER) THEN
         ALLOCATE(WAVE(msh%NDOF))
         ALLOCATE(NSG_LSTB(NPGROUPS))
         ALLOCATE(NSG_LST(NPGROUPS)) 
         ALLOCATE(ISRCLST(src%NSRC))
      ELSE
         ALLOCATE(WAVE(1))
      ENDIF
      DO 1000 IFREQL=1,frq%NFREQ
!
!....... get frequency
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1
         IF ( (IFREQL - 1)*NPGROUPS + 1 > frq%NFREQ) GOTO 1000 !done
!
!....... set and send source group info
         src%NSG = 0
         IF (MYNID == MASTER) THEN
            IF (IFREQ <= frq%NFREQ) THEN
               IF (frq%CFTYPE(IFREQ) == 'S') THEN
                  CALL FILL_SRCPRM_SB(IFREQ, .TRUE., m1d%VFAST_SRF,msh%AZTOL,  &
                                      msh%AOITOL, SRC)
               ELSEIF (frq%CFTYPE(IFREQ) == 'B') THEN
                  CALL FILL_SRCPRM_SB(IFREQ,.FALSE., m1d%VFAST_BDY,msh%AZTOL,  &
                                      msh%AOITOL, SRC)
               ELSE
                  WRITE(*,*) 'xbielak25: Cant determine frequency type!'
                  IERR = 1 
                  GOTO 500
               ENDIF
            ELSE
               src%NSG = 0
            ENDIF
            NSG_LSTB(1:NPGROUPS) = 0
            NSG_LSTB(IPGROUP+1) = src%NSG
            CALL MPI_ALLREDUCE(src%NSG ,NSGMAX,        1,MPI_INTEGER,MPI_MAX, &
                               MYHD_COMM,MPIERR)
            CALL MPI_ALLREDUCE(NSG_LSTB,NSG_LST,NPGROUPS,MPI_INTEGER,MPI_SUM, &
                               MYHD_COMM,MPIERR)
         ENDIF
         CALL MPI_BCAST(NSGMAX,1,MPI_INTEGER, MASTER,MYSLV_COMM,MPIERR)
         IF (IFREQ <= frq%NFREQ) CALL BCAST_SRC_PTRS(MYNID,MYSLV_COMM,MASTER, SRC)  
!
!....... loop on source groups
         DO 2000 ISG=1,NSGMAX
!
!.......... have master flash states
            IF (MYID == MASTER) THEN
               WRITE(*,*) 
               MYSRC = 0
               PYAVG = 0.D0
               DO 201 JPGROUP=0,NPGROUPS-1
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1)
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN
                     IF (JPGROUP > 0) THEN
                        MYSRC = JPGROUP*NPROCS/NPGROUPS
                        CALL MPI_RECV(PYAVG,1,MPI_DOUBLE_PRECISION, MYSRC,MPI_ANY_TAG, &
                                      MPI_COMM_WORLD,STAT,MPIERR)
                     ELSE
                        PYAVG = src%PYAVG(ISG)  
                     ENDIF !end check on recieving or me
!
!................... write surface wave or body wave
                     IF (frq%CFTYPE(JFREQ) == 'S') THEN  
                        WRITE(*,9410) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                     ELSE
                        WRITE(*,9411) JPGROUP+1,ISG,PYAVG,frq%FREQ(JFREQ)
                     ENDIF 
                  ENDIF !end check on working frequency
  201          CONTINUE !loop on processs groups
            ELSE
               IF (MYNID == MASTER) THEN
                  MYDEST = MASTER
                  MYTAG  = MYID
                  IF (IFREQ <= frq%NFREQ .AND. src%NSG > 0 .AND. ISG <= src%NSG) THEN 
                     CALL MPI_SEND(src%PYAVG(ISG),1,MPI_DOUBLE_PRECISION, MASTER,MYTAG, & 
                                   MPI_COMM_WORLD,MPIERR)
                  ENDIF
               ENDIF
            ENDIF !end check on myid
            CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
!
!.......... assemble impedance matrix for average py value in group 
            IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Assembling impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN
               CALL ASMBLE_DRIVER(msh%nzloc,frq%CFTYPE(IFREQ),frq%FREQ(IFREQ), &
                                  src%PYAVG(ISG), MSH, mid%A_LOC,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'xbielak25: An error occurred calling asmble_driver',MYID
                  GOTO 500
               ENDIF
            ENDIF
            IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Factoring impedance matrix...'
            IF (ISG <= src%NSG .AND. IFREQ <= frq%NFREQ) THEN
               mid%JOB = 2 !factorization
               CALL CMUMPS(MID) 
               IF (mid%INFO(1) < 0) THEN
                  WRITE(*,*) 'xbielak25: An error occurred in the factorization!',MYID
                  IERR = 1
                  GOTO 500 
               ENDIF
            ENDIF !end check 
!
!.......... have master tell which sources are being solved
            IF (MYID == MASTER) THEN
               MYSRC = 0 
               ISRCLST(1:src%NSRC) = 0
               DO 290 JPGROUP=0,NPGROUPS-1 
                  JFREQ = (IFREQL - 1)*NPGROUPS + JPGROUP + 1
                  NSG_GRP = NSG_LST(JPGROUP+1) 
                  IF (JFREQ <= frq%NFREQ .AND. NSG_GRP > 0 .AND. ISG <= NSG_GRP) THEN
                     IF (JPGROUP > 0) THEN !otherwise srclist already set
                        MYSRC = JPGROUP*NPROCS/NPGROUPS
                        IF (NSG_LST(JPGROUP+1) > 0 .AND. JFREQ <= frq%NFREQ) THEN
                           CALL MPI_RECV(ISRCLST,src%NSRC,MPI_INTEGER,  &
                                         MYSRC,MPI_ANY_TAG,             &
                                         MPI_COMM_WORLD,STAT,MPIERR)
                        ENDIF
                     ELSE
                        KSRC = 0
                        DO 294 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1) - 1 
                           KSRC = KSRC + 1
                           ISRCLST(KSRC) = src%ISRCPRM(JSRC)
  294                   CONTINUE
                     ENDIF !end check on jpgroup
                     DO 296 KSRC=1,src%NSRC
                        ISRC = ISRCLST(KSRC)
                        IF (ISRC == 0) GOTO 295
                        IF (frq%CFTYPE(JFREQ) == 'S') THEN
                           WRITE(*,9412) JPGROUP+1,ISRC,src%AOI(ISRC),src%BAZN(ISRC), &
                                         src%PYTAB(JFREQ,ISRC)
                        ELSE
                           WRITE(*,9413) JPGROUP+1,ISRC,src%AOI(ISRC),src%BAZN(ISRC), &
                                         src%PYTAB(JFREQ,ISRC)
                        ENDIF !end check on surface/body wave source
  296                CONTINUE !out of sources
  295                CONTINUE !break ahead, out of sources
                  ENDIF !end check on working
  290          CONTINUE !loop on sources groups
            ELSE !slaves on master communicator
               IF (MYNID == MASTER) THEN
                  MYDEST = MASTER
                  MYTAG =  MYID
                  ISRCLST(1:src%NSRC) = 0
                  IF (IFREQ <= frq%NFREQ .AND. src%NSG > 0 .AND. ISG <= src%NSG) THEN
                     KSRC = 0
                     DO 297 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
                        KSRC = KSRC + 1
                        ISRCLST(KSRC) = src%ISRCPRM(JSRC) 
  297                CONTINUE  
                     CALL MPI_SEND(ISRCLST,src%NSRC,MPI_INTEGER, MASTER,MYTAG, & 
                                   MPI_COMM_WORLD,MPIERR)
                  ENDIF !end check on working
               ENDIF !end check on myid
            ENDIF !end check on myid
            IF (MYID == MASTER) THEN 
               WRITE(*,*) 'xbielak25: Calculating effective forces'
               WRITE(*,*) 'xbielak25: Solving Ax = b' 
               IF (rcv%NREC > 0) WRITE(*,*) 'xbielak25: Extracting responses'
            ENDIF
            IF (ISG > src%NSG)     GOTO 2001 !no solves, leave
            IF (IFREQ > frq%NFREQ) GOTO 2001 !no solves, leave
!----------------------------------------------------------------------------------------!
!                                 Loop on sources in source group                        !
!----------------------------------------------------------------------------------------!
            DO 3000 JSRC=src%ISGPTR(ISG),src%ISGPTR(ISG+1)-1
               ISRC = src%ISRCPRM(JSRC) !extract original source ID
!
!............. generate 1D solution
               IF (MYNID == MASTER) THEN
                  CALL LOAD_GRNS25(NDIM,msh%NDOF,msh%NNPE,NDIM, src%SRCTYP(ISRC),ISRC,   &
                                   frq%FREQ(IFREQ), src%SOURCE(IFREQ,ISRC),msh%IDOFSE,   &
                                   UE,IERR)
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xbielak25: Error in load_grns25!'
                     GOTO 500
                  ENDIF
                  IF (.NOT.LFILES) THEN
                     CALL PLOT_ELCRESP_VTK(NGNOD,NDIM,msh%NEN, NDIM,msh%NNPG,msh%NDOF,   &
                                           msh%NLXI,msh%NLETA,PROJNM, msh%NELEM,3, ISRC, &
                                           1,frq%FREQ(IFREQ), msh%LM,msh%IENG,           &
                                           msh%XLOCS,msh%ZLOCS, UE)
                  ENDIF
!                 IF (MYID == MASTER) &
!                 WRITE(*,*) 'xbielak25: Calculating effective forces...'
               ENDIF  !end check on master
               CALL MPI_BCAST(UE,msh%NDOF,MPI_COMPLEX, MASTER,MYSLV_COMM,MPIERR)
               CALL ASMB_PEFF(MASTER,MYSLV_COMM, msh%NDOF,msh%NDOFL,msh%NZLOC, msh%CNP,  &
                              msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC, mid%A_LOC,UE, WAVE)
               IF (MYNID == MASTER) CALL CCOPY(msh%NDOF,WAVE,1,mid%RHS,1)
               IF (MYNID == MASTER .AND. .NOT.LFILES) THEN
                  CALL PLOT_ELCRESP_VTK(NGNOD,NDIM,msh%NEN, NDIM,msh%NNPG,msh%NDOF,      &
                                        msh%NLXI,msh%NLETA,PROJNM, msh%NELEM,3, ISRC,    &
                                        2,frq%FREQ(IFREQ), msh%LM,msh%IENG,              &
                                        msh%XLOCS,msh%ZLOCS, mid%RHS)
               ENDIF
               mid%JOB = 3 !solve phase
               CALL CMUMPS(MID) 
               IF (mid%INFO(1) < 0) THEN
                  WRITE(*,*) 'xbielak25: An error occurred in the solve phase!'
                  IERR = 1 
                  GOTO 500 
              ENDIF
!
!............. head processes extract solution
               IF (MYNID == MASTER) THEN 
                  CALL ADDBLK(msh%NDOF,msh%CNP,UE, mid%RHS) !add background field in
                  IF (rcv%NREC > 0) THEN
!                    IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Extracting responses...'
                     CALL EXRESP25D(NDIM,rcv%NREC,msh%NDOF,                              &
                                    1,rcv%NREC,NDIM, frq%FREQ(IFREQ),                    &
                                    src%PYTAB(IFREQ,ISRC),rcv%YREC, mid%RHS, rcv%MRDOF,  &
                                    rcv%RECV(1:NDIM,IFREQ,1:rcv%NREC),                   &
                                    EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC))
                  ENDIF
                  IF (.NOT.LFILES) THEN
                     CALL PLOT_ELCRESP_VTK(NGNOD,NDIM,msh%NEN, NDIM,msh%NNPG,msh%NDOF,   &
                                         msh%NLXI,msh%NLETA,PROJNM, msh%NELEM,3, ISRC,   &
                                         3,frq%FREQ(IFREQ), msh%LM,msh%IENG,             &
                                         msh%XLOCS,msh%ZLOCS, mid%RHS) 
                  ENDIF
                  CALL EWAVOUT2D(PROJNM,NDIM,msh%NEN,NGNOD,msh%NNPG,msh%NDOF, msh%NELEM, &
                                 NDIM,msh%NLXI,msh%NLETA,NGNOD,ISRC,frq%FREQ(IFREQ),     &
                                 msh%IENG,msh%LM, msh%XLOCS,msh%ZLOCS,mid%RHS,IERR) 
                  IF (IERR /= 0) THEN
                     WRITE(*,*) 'xbielak25: Error writing 2.5D wavefield'
                     GOTO 500
                  ENDIF
               ENDIF !end check on master
               !if (mynid == master) work(:) = mid%rhs(:)
               !call verify_steve(mynid,master,myslv_comm, msh%ndof,msh%ndofl,msh%nzloc, msh%cnp, &
               !                  msh%mydofs,msh%irptr_loc,msh%jcptr_loc, mid%a_loc,ue,work) 
 3000       CONTINUE !loop on sources
 2001       CONTINUE !break ahead for loop on source group
 2000    CONTINUE !loop on source groups
!
!....... clean source info
         IF (src%NSG > 0 .AND. IFREQ <= frq%NFREQ) THEN
            IF (ASSOCIATED(src%ISRCPRM)) DEALLOCATE(src%ISRCPRM)
            IF (ASSOCIATED(src%ISGPTR))  DEALLOCATE(src%ISGPTR)
            IF (ASSOCIATED(src%PYAVG))   DEALLOCATE(src%PYAVG)
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
 1000 CONTINUE !Loop on frequencies 
!
!.... reduce estimates
      IF (MYNID == MASTER .AND. rcv%NREC > 0) THEN
         IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Reducing estimates...'
         CALL REDEST4(MYID,MASTER,MYHD_COMM,NPGROUPS,IPGROUP,                           &
                      NDIM,frq%NFREQ,rcv%NREC, NDIM,frq%NFREQ,rcv%NREC,src%NSRC,1, EST)
         IF (MYID == MASTER) THEN
            IF (frq%NFREQ_SRF == 0 .OR. frq%NFREQ_BDY == 0) THEN
               WRITE(*,*) 'xbielak25: Writing estimates...'
               CALL WTEEST25(PROJNM, frq%NFREQ,rcv%NREC,frq%NFREQ, rcv%NREC,src%NSRC,   &
                             frq%FREQ, EST(1,:,:,:),EST(2,:,:,:),EST(3,:,:,:), IERR) 
            ELSE
               WRITE(*,*) 'xbielak25: Writing surface wave estimates...'
               ALLOCATE(EST_LOC(NDIM,frq%NFREQ_SRF,rcv%NREC,src%NSRC_SRF))
               ALLOCATE(FREQ_LOC(frq%NFREQ_SRF))
               JFREQ = 0
               DO 50 IFREQ=1,frq%NFREQ
                  IF (frq%CFTYPE(IFREQ) == 'S') THEN
                     JFREQ = JFREQ + 1
                     FREQ_LOC(JFREQ) = frq%FREQ(IFREQ)
                     JSRC = 0
                     DO 51 ISRC=1,src%NSRC
                        IF (src%SRCTYP(ISRC)(1:1) == 'S') THEN
                           JSRC = JSRC + 1
                           DO 52 IREC=1,rcv%NREC
                              DO 53 I=1,NDIM
                                 EST_LOC(I,JFREQ,IREC,JSRC) = EST(I,IFREQ,IREC,ISRC)
   53                         CONTINUE
   52                      CONTINUE
                        ENDIF
   51                CONTINUE
                  ENDIF
   50          CONTINUE
               FILENM = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_srf'
               FILENM = ADJUSTL(FILENM)
               CALL WTEEST25(FILENM, frq%NFREQ_SRF,rcv%NREC,frq%NFREQ_SRF,             &
                             rcv%NREC,src%NSRC_SRF, FREQ_LOC,                          &
                             EST_LOC(1,:,:,:),EST_LOC(2,:,:,:),EST_LOC(3,:,:,:), IERR)
               DEALLOCATE(EST_LOC)
               DEALLOCATE(FREQ_LOC)
               WRITE(*,*) 'xbielak25: Writing body wave estimates...'
               ALLOCATE(EST_LOC(NDIM,frq%NFREQ_BDY,rcv%NREC,src%NSRC_BDY))
               ALLOCATE(FREQ_LOC(frq%NFREQ_BDY))
               JFREQ = 0
               DO 54 IFREQ=1,frq%NFREQ
                  IF (frq%CFTYPE(IFREQ) == 'B') THEN
                     JFREQ = JFREQ + 1
                     FREQ_LOC(JFREQ) = frq%FREQ(IFREQ)
                     JSRC = 0
                     DO 55 ISRC=1,src%NSRC
                        IF (src%SRCTYP(ISRC)(1:1) == 'P') THEN
                           JSRC = JSRC + 1
                           DO 56 IREC=1,rcv%NREC
                              DO 57 I=1,NDIM
                                 EST_LOC(I,JFREQ,IREC,JSRC) = EST(I,IFREQ,IREC,ISRC)
   57                         CONTINUE
   56                      CONTINUE
                        ENDIF
   55                CONTINUE
                  ENDIF
   54          CONTINUE
               FILENM = ' '
               FILENM = TRIM(ADJUSTL(PROJNM))//'_bdy'
               FILENM = ADJUSTL(FILENM)
               CALL WTEEST25(FILENM, frq%NFREQ_BDY,rcv%NREC,frq%NFREQ_BDY,             &
                             rcv%NREC,src%NSRC_BDY, FREQ_LOC,                          &
                             EST_LOC(1,:,:,:),EST_LOC(2,:,:,:),EST_LOC(3,:,:,:), IERR)
               DEALLOCATE(EST_LOC)
               DEALLOCATE(FREQ_LOC)
            ENDIF
            IF (IERR /= 0) THEN
               WRITE(*,*) 'xbielak25: Error writing estimate field...'
               GOTO 500
            ENDIF
         ENDIF
      ENDIF
!
!.... clean Greens functions files
!     IF (MYID == MASTER .AND. LFILES) THEN
!        LEX = LISDIR('./grns') 
!        IF (LEX .AND. .NOT.LCGRNS) CALL SYSTEM('rm -rf ./grns')
!     ENDIF
!----------------------------------------------------------------------------------------!
!                 Error handling and memory cleaning                                     !
!----------------------------------------------------------------------------------------!
  500 CONTINUE !break ahead for errors or out of frequencies
      IF (IERR /= 0) THEN
         WRITE(*,9406) MYID
 9406    FORMAT(' -------------------------------------------------------------',/, &
                ' -               An Error was Detected on Process ',I6,'     -',/, &
                ' -                  I will now abort the program             -',/, &
                ' -------------------------------------------------------------',/)
         CALL MPI_ABORT(MPI_COMM_WORLD,30,MPIERR)  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF (MYID == MASTER) WRITE(*,*) 'xbielak25: Freeing memory...'
!.... background field
      IF (ALLOCATED(UE)) DEALLOCATE(UE) 
!.... communication list
      IF (ALLOCATED(NSG_LSTB)) DEALLOCATE(NSG_LSTB)
      IF (ALLOCATED(NSG_LST))  DEALLOCATE(NSG_LST) 
      IF (ALLOCATED(ISRCLST))  DEALLOCATE(ISRCLST) 
!.... model
      IF (ASSOCIATED(msh%ECOEFF))    DEALLOCATE(msh%ECOEFF)
      IF (ASSOCIATED(msh%DENS))      DEALLOCATE(msh%DENS)
      IF (ASSOCIATED(msh%XD))        DEALLOCATE(msh%XD)
      IF (ASSOCIATED(msh%ZD))        DEALLOCATE(msh%ZD)
      IF (ASSOCIATED(msh%QP))        DEALLOCATE(msh%QP)
      IF (ASSOCIATED(msh%QS))        DEALLOCATE(msh%QS) 
      IF (ASSOCIATED(msh%XLOCS))     DEALLOCATE(msh%XLOCS)
      IF (ASSOCIATED(msh%ZLOCS))     DEALLOCATE(msh%ZLOCS)
      IF (ASSOCIATED(msh%XLOCSE))    DEALLOCATE(msh%XLOCSE)
      IF (ASSOCIATED(msh%ZLOCSE))    DEALLOCATE(msh%ZLOCSE)
      IF (ASSOCIATED(msh%XIPTS))     DEALLOCATE(msh%XIPTS)
      IF (ASSOCIATED(msh%ETAPTS))    DEALLOCATE(msh%ETAPTS) 
      IF (ASSOCIATED(msh%CNNPG))     DEALLOCATE(msh%CNNPG)
      IF (ASSOCIATED(msh%CDOMAIN))   DEALLOCATE(msh%CDOMAIN) 
      IF (ASSOCIATED(msh%CNP))       DEALLOCATE(msh%CNP) 
!.... graph
      IF (ASSOCIATED(msh%LM))        DEALLOCATE(msh%LM)
      IF (ASSOCIATED(msh%IENG))      DEALLOCATE(msh%IENG)
      IF (ASSOCIATED(msh%IDOFSE))    DEALLOCATE(msh%IDOFSE)
      IF (ASSOCIATED(msh%MYDOFS))    DEALLOCATE(msh%MYDOFS)
      IF (ASSOCIATED(msh%MYELEM))    DEALLOCATE(msh%MYELEM)       
      IF (ASSOCIATED(msh%IRPTR_LOC)) DEALLOCATE(msh%IRPTR_LOC)
      IF (ASSOCIATED(msh%JCPTR_LOC)) DEALLOCATE(msh%JCPTR_LOC)
!.... source
      IF (ASSOCIATED(src%PYTAB))   DEALLOCATE(src%PYTAB) 
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(src%SOURCE))  DEALLOCATE(src%SOURCE)
         IF (ASSOCIATED(src%SRCTYP))  DEALLOCATE(src%SRCTYP)
         IF (ASSOCIATED(src%BAZN))    DEALLOCATE(src%BAZN)
         IF (ASSOCIATED(src%AOI))     DEALLOCATE(src%AOI)
      ENDIF
!.... receiver
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(rcv%YREC))  DEALLOCATE(rcv%YREC)
         IF (ASSOCIATED(rcv%MRDOF)) DEALLOCATE(rcv%MRDOF)
      ENDIF
!.... frequency
      IF (ASSOCIATED(frq%FREQ))    DEALLOCATE(frq%FREQ)
!.... mumps
      IF (ASSOCIATED(mid%A_LOC))   DEALLOCATE(mid%A_LOC)
      IF (ASSOCIATED(mid%IRN_LOC)) DEALLOCATE(mid%IRN_LOC)
      IF (ASSOCIATED(mid%JCN_LOC)) DEALLOCATE(mid%JCN_LOC)
      IF (ALLOCATED(WAVE))         DEALLOCATE(WAVE)
      IF (MYNID == MASTER) THEN
         IF (ASSOCIATED(mid%RHS))   DEALLOCATE(mid%RHS)
         IF (ALLOCATED(EST))        DEALLOCATE(EST)
         IF (ASSOCIATED(rcv%RECV))  DEALLOCATE(rcv%RECV)
         IF (ASSOCIATED(rcv%MRDOF)) DEALLOCATE(rcv%MRDOF) 
      ENDIF
      mid%JOB =-2 
      CALL CMUMPS(MID) 
!
!.... all done 
      TESIM = MPI_WTIME()
      IF (MYID == MASTER) THEN
         WRITE(*,9405) (TESIM - TSSIM)/60.D0
 9405    FORMAT(' xbielak25: Simulation time:',F8.2,' minutes')
      ENDIF 
      CALL MPI_FINALIZE(MPIERR)
!----------------------------------------------------------------------------------------!
!                                     Format statements                                  !
!----------------------------------------------------------------------------------------!
 9410 FORMAT(' ---------------------------------------------------------'   ,/, &
             ' -  xbielak25: Group:',I4,'                                -' ,/, &
             ' -             Propagating surface wave source group ',I3,' -',/, &
             ' -             With average slowness in y',E12.4,'    -'      ,/, &
             ' -             At frequency',F12.5,' Hz               -'      ,/, &
             ' ---------------------------------------------------------',/)
 9411 FORMAT(' ---------------------------------------------------------'   ,/, &
             ' -  xbielak25: Group:',I4,'                                -' ,/, &
             ' -             Propagating body wave source group',I3,'     -',/, &
             ' -             With average slowness in y',E12.4,'    -'      ,/, &
             ' -             At frequency',F12.5,' Hz               -'      ,/, &
             ' ---------------------------------------------------------',/) 
 9412 FORMAT(/,' xbielak25: Process group: ',I3,/, &
               '            Will process surface wave source ',I3,/,           &
               '            With angle of incidence ',F8.3,' degrees',/,  &
               '            Corrected azimuth ',F8.3,' degrees',/,    &
               '            And y slowness ',G10.3,' s/m',/)
 9413 FORMAT(/,' xbielak25: Process group: ',I3,/, &
               '            Will process body wave source ',I3,/,           &
               '            With angle of incidence ',F8.3,' degrees',/,  &
               '            Corrected azimuth ',F8.3,' degrees',/,    &
               '            And y slowness ',G10.3,' s/m',/)

      STOP
      END
