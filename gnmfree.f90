      SUBROUTINE GNMFREE(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM, &
                         IPGROUP,NPGROUPS,  FRQ,SRC,INV, P,      &
                         Z,IERR)
!
!     Loads the Jacobian matrices vectors from disk and multplies by the search 
!     direction.  The idea is to calculate, 
!        Re{ adj(J)J }p  
!       =Re{ [Re{J^T} - iIm{J^T}]*[Re{J} + iIm{J}] }p 
!       ={ Re{J^T}Re{J} + Im{J^T}Im{J} } p 
!       =Re{J^T}Re{J}p + Im{J^T}Im{J}p = z 
!
!       [adj(J)_1, adj(J)_2, ... adj(J)_nf*nsrc][ J_1 ] {p}  
!                                               [ J_2 ]
!                                               [  .  ]
!                                               [  .  ]
!                                               [  .  ]
!                                               [ J_n ] 
!     = [adj(J)_1 J_1 p + adj(J)_2 J_2 p + ...]
!     = z
!    
!     Of course, we set the calculation so that each process frequency group has part
!     of the Jacobians, and for each frequency group the Jacobian is itself distributed
!     amongst all process in that frequency group.  This form of matrix vector 
!     multiplication should be useful if solving the Gauss-Newton approximation:
!        B_a p =-grad 
!        where B_a = \sum_{nf} \sum_{nsrc} adj(J_i) J_i
! 
!    - B. Baker January 1 2013.   
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ 
      TYPE (SRC_INFO) SRC 
      TYPE (INV_INFO) INV 
      REAL*8, INTENT(IN) :: P(*)
      INTEGER*4, INTENT(IN) :: MYID,MYNID,MASTER,MYHD_COMM, MYSLV_COMM, &
                               IPGROUP,NPGROUPS
      REAL*8, INTENT(OUT) :: Z(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*16, ALLOCATABLE :: JP1(:), JP2(:)
      REAL*8, ALLOCATABLE :: Z1(:), Z2(:), RWORK(:), PW(:)
      REAL*4, ALLOCATABLE :: RCJLOC(:,:), CCJLOC(:,:)
      CHARACTER(80) FLNAME
      CHARACTER(12) CFREQ
      CHARACTER(5) CID, CSRC
      COMPLEX*16 CARG, CZERO
      REAL*8 ALPHAM
      INTEGER*4 NOBS_IN, MOBS, IFREQL, IFREQ, ISRC, INPGL, IUNIT, &
                IOBS, ILOC, I, IREC, NA35, NNPGL_IN, NDIM_IN, NREC_IN, MPIERR
      LOGICAL*4 LOBS
      PARAMETER(IUNIT = 65) 
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
!
!----------------------------------------------------------------------------------------!
!
!.... send search direction to all frequency group masters
      IERR = 0
      NA35 = inv%NA35
      ALLOCATE(RWORK(NA35)) !model buffer; slave to master
      CALL DSCAL(NA35,0.D0,RWORK,1) !zeros won't change since NNPGL is constant
      ALLOCATE(PW(NA35))    !current search direction p 
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) CALL DCOPY(NA35,P,1,PW,1)
         CALL MPI_BCAST(PW,NA35,MPI_DOUBLE_PRECISION, MASTER,MYHD_COMM,MPIERR) 
         ALLOCATE(Z1(NA35))
         ALLOCATE(Z2(NA35))
         CALL DSCAL(NA35,0.D0,Z2,1) !null out running sum 
      ELSE
         ALLOCATE(Z1(1))
         ALLOCATE(Z2(1))
      ENDIF
      CALL MPI_BCAST(PW,NA35,MPI_DOUBLE_PRECISION, MASTER,MYSLV_COMM,MPIERR) 
!
!.... calculate max space required
      CID(1:5) = ' '
      WRITE(CID,'(I5)') MYID
      CID = ADJUSTL(CID)
      MOBS = 0
      DO 50 IFREQL=1,frq%NFREQ
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1
         IF ( (IFREQL - 1)*NPGROUPS + 1 > frq%NFREQ) GOTO 85
         IF (IFREQ > frq%NFREQ) GOTO 91
         DO 51 ISRC=1,src%NSRC 
!
!.......... read the jacobian header 
            FLNAME(1:80) = ' ' 
            CFREQ(1:12) = ' ' 
            CSRC(1:5) = ' ' 
            WRITE(CFREQ,'(F12.5)') frq%FREQ(IFREQ)
            CFREQ = ADJUSTL(CFREQ) 
            WRITE(CSRC ,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC) 
            FLNAME ='./scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
            FLNAME = ADJUSTL(FLNAME)
            OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
            READ(IUNIT,IOSTAT=IERR) NOBS_IN,NNPGL_IN, NDIM_IN, NREC_IN
            MOBS = MAX(NOBS_IN,MOBS) 
            CLOSE(IUNIT)
            IF (IERR /= 0) THEN
               WRITE(*,*) 'gnmfree: Error reading Jacobian header:',MYID
               RETURN
            ENDIF
            IF (NDIM_IN /= NDIM) THEN
               WRITE(*,*) 'gnmfree: Dimension mismatch!:',MYID
               IERR = 1
               RETURN
            ENDIF
            IF (NNPGL_IN /= inv%NNPGL) THEN
               WRITE(*,*) 'gnmfree: NNPGL inconsistent!:',MYID
               IERR = 1
               RETURN
            ENDIF
   51    CONTINUE !Loop on sources
   91    CONTINUE !break ahead for new frequency
   50 CONTINUE !Loop on frequencies
   85 CONTINUE !done with frequencies
!
!.... set max buffers 
!     ALLOCATE(R(MOBS))
!     ALLOCATE(Q(MOBS))
!     ALLOCATE(RW(MOBS))
!     ALLOCATE(QW(MOBS))
      ALLOCATE(JP1(MOBS))
      ALLOCATE(JP2(MOBS)) 
      ALLOCATE(RCJLOC(MOBS,inv%NNPGL),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gnmfree: Error setting space for RCJLOC on process',MYID
         RETURN
      ENDIF
      ALLOCATE(CCJLOC(MOBS,inv%NNPGL)) 
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gnmfree: Error setting space for CCJLOC on process',MYID
         RETURN
      ENDIF
!
!.... loop on frequencies, sources for matrix matrix vector multiplies 
      DO 100 IFREQL=1,frq%NFREQ
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1
         IF ( (IFREQL - 1)*NPGROUPS + 1 > frq%NFREQ) GOTO 105 
         IF (IFREQ > frq%NFREQ) GOTO 101
!
!....... loop on sources
         DO 200 ISRC=1,src%NSRC
!
!.......... read the jacobian 
            FLNAME(1:80) = ' '
            CFREQ(1:12) = ' '
            CSRC(1:5) = ' '
            WRITE(CFREQ,'(F12.5)') frq%FREQ(IFREQ)
            CFREQ = ADJUSTL(CFREQ) 
            WRITE(CSRC ,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC) 
            FLNAME ='./scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
            FLNAME = ADJUSTL(FLNAME)
            OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
            READ(IUNIT,IOSTAT=IERR) NOBS_IN,NNPGL_IN, NDIM_IN, NREC_IN
            IF (IERR /= 0) GOTO 105
            READ(IUNIT,IOSTAT=IERR)((LOBS,I=1,NDIM),IREC=1,NREC_IN)
            IF (IERR /= 0) GOTO 105
            READ(IUNIT,IOSTAT=IERR)((RCJLOC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,inv%NNPGL)
            IF (IERR /= 0) GOTO 105
            READ(IUNIT,IOSTAT=IERR)((CCJLOC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,inv%NNPGL)
            IF (IERR /= 0) GOTO 105
            IF (IERR /= 0) THEN
               WRITE(*,*) 'gnmfree: Error reading cjloc!'
               CLOSE(IUNIT)
               IERR = 3
               GOTO 105
            ENDIF
            CLOSE(IUNIT) 
!
!.......... migration: Jacobian corresponds built w/ Greens fns, need an updated STF
            IF (MYNID == MASTER) THEN
               IF (inv%LMIGR) THEN !we need to add in source time functions
                  ALPHAM = CABS(src%SOURCE(IFREQ,ISRC))**2
               ELSE
                  ALPHAM = 1.0
               ENDIF 
            ENDIF
!
!.......... matrix vector multiply; Re{J}m, Im{J}m 
            DO 500 IOBS=1,NOBS_IN
               JP1(IOBS) = CZERO 
!              R(IOBS) = 0.D0
!              Q(IOBS) = 0.D0
               DO 501 INPGL=1,inv%NNPGL 
                  ILOC = inv%MYGRAD(INPGL) 
                  CARG = DCMPLX(DBLE(RCJLOC(IOBS,INPGL)), DBLE(CCJLOC(IOBS,INPGL)))
                  JP1(IOBS) = JP1(IOBS) + CARG*DCMPLX(PW(ILOC),0.D0) 
!                 R(IOBS) = R(IOBS) + DBLE(RCJLOC(IOBS,INPGL))*PW(ILOC) 
!                 Q(IOBS) = Q(IOBS) + DBLE(CCJLOC(IOBS,INPGL))*PW(ILOC) 
 501           CONTINUE !loop on rows
 500        CONTINUE !loop on columns
            CALL MPI_ALLREDUCE(JP1,JP2,NOBS_IN,MPI_DOUBLE_COMPLEX, MPI_SUM, &
                               MYSLV_COMM,MPIERR) 
!
!.......... transpose matrix vector multiply
            RWORK(1:inv%NA35) = 0.D0
            DO 502 INPGL=1,inv%NNPGL
               ILOC = inv%MYGRAD(INPGL) 
!              RWORK(ILOC) = RWORK(ILOC) &
!                          + SDOT(NOBS_IN,RCJLOC(1:NOBS_IN,INPGL),1,RW,1)
!              RWORK(ILOC) = RWORK(ILOC) &
!                          + SDOT(NOBS_IN,CCJLOC(1:NOBS_IN,INPGL),1,QW,1)
               DO 503 IOBS=1,NOBS_IN 
                  CARG = DCMPLX(DBLE(RCJLOC(IOBS,INPGL)),-DBLE(CCJLOC(IOBS,INPGL)))
                  RWORK(ILOC) = RWORK(ILOC) + DREAL(CARG*JP2(IOBS))
 503           CONTINUE !loop on observations read in
 502        CONTINUE !loop on local points in gradient
!
!.......... synchronize loops by reducing onto head node
            CALL MPI_REDUCE(RWORK,Z1,NA35,MPI_DOUBLE_PRECISION, MPI_SUM, &
                            MASTER, MYSLV_COMM,MPIERR)   
            IF (MYNID == MASTER) CALL DAXPY(NA35,ALPHAM,Z1,1,Z2,1) !z2 = z1 + alpha*z2 
  200    CONTINUE 
  101    CONTINUE !break ahead, frequency group done
  100 CONTINUE !loop on frequencies
  105 CONTINUE !break ahead, done with frequencies or had an errror
!
!.... reduce onto head node
      IF (MYNID == MASTER) &
      CALL MPI_REDUCE(Z2,Z,NA35,MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MASTER, MYHD_COMM,MPIERR)  
!
!.... free space
      IF (ALLOCATED(RCJLOC)) DEALLOCATE(RCJLOC) 
      IF (ALLOCATED(CCJLOC)) DEALLOCATE(CCJLOC) 
      IF (ALLOCATED(RWORK))  DEALLOCATE(RWORK) 
      IF (ALLOCATED(JP1))    DEALLOCATE(JP1)
      IF (ALLOCATED(JP2))    DEALLOCATE(JP2)
!     IF (ALLOCATED(R))      DEALLOCATE(R)
!     IF (ALLOCATED(Q))      DEALLOCATE(Q) 
!     IF (ALLOCATED(RW))     DEALLOCATE(RW)
!     IF (ALLOCATED(QW))     DEALLOCATE(QW)
      IF (ALLOCATED(PW))     DEALLOCATE(PW) 
      IF (ALLOCATED(Z1))     DEALLOCATE(Z1)
      IF (ALLOCATED(Z2))     DEALLOCATE(Z2)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GNGRAD(MYID,MYNID,MASTER,MYHD_COMM,MYSLV_COMM, &
                        IPGROUP,NPGROUPS,  MSH,RCV,FRQ,SRC,INV,IERR) 
!
!     This is a debugging routine to generate a gradient
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ
      TYPE (SRC_INFO) SRC
      TYPE (INV_INFO) INV
      TYPE (RECV_INFO) RCV
      TYPE (MESH_INFO) MSH
      INTEGER*4, INTENT(IN) :: MYID, MYNID, MASTER, MYHD_COMM, MYSLV_COMM, &
                               IPGROUP, NPGROUPS 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*8, ALLOCATABLE :: RESID(:,:), OBSR(:,:) 
      REAL*4, ALLOCATABLE :: RCJLOC(:,:), CCJLOC(:,:), GRADB(:), &
                             GRADL(:), TEMP(:)  
      LOGICAL*4, ALLOCATABLE :: LOBS(:,:)
      INTEGER*4 NOBS_IN, NNPGL_IN, NDIM_IN, NREC_IN, IFREQL, IFREQ, ISRC, I,  &
                IREC, IOBS, ILOC, INPGL, IUNIT, MPIERR 
      COMPLEX*8 CARG, SAZ, CAZ 
      REAL*8 PI180 
      REAL*4 WMIN, WMAX
      CHARACTER(80) FLNAME
      CHARACTER(12) CFREQ
      CHARACTER(5) CID, CSRC
      PARAMETER(IUNIT = 65) 
      PARAMETER(PI180 = 0.017453292519943295D0)
!
!----------------------------------------------------------------------------------------!
!
!.... set space
      IERR = 0
      ALLOCATE(RESID(NDIM,rcv%NREC))
      ALLOCATE(GRADB(inv%NA35))
      ALLOCATE(TEMP (inv%NA35)) 
      ALLOCATE(CCJLOC(NDIM*rcv%NREC,inv%NNPGL))
      ALLOCATE(RCJLOC(NDIM*rcv%NREC,inv%NNPGL))
      ALLOCATE(LOBS(NDIM,rcv%NREC)) 
!
!.... rotation information
      IF (MYNID == MASTER) THEN
         CAZ = CMPLX(SNGL(DCOS(msh%AZMOD*PI180)),0.0)
         SAZ = CMPLX(SNGL(DSIN(msh%AZMOD*PI180)),0.0)
         ALLOCATE(OBSR(NDIM,rcv%NREC))
         ALLOCATE(GRADL(inv%NA35))
         GRADL(1:inv%NA35) = 0.0
      ENDIF
!.... set character process ID for file IO
      CID(1:5) = ' '
      WRITE(CID,'(I5)') MYID
      CID = ADJUSTL(CID)
      DO 100 IFREQL=1,frq%NFREQ
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1 
         IF ( (IFREQL - 1)*NPGROUPS + 1 > frq%NFREQ) GOTO 105 
         IF (IFREQ > frq%NFREQ) GOTO 101 
!
!....... loop on sources
         DO 200 ISRC=1,src%NSRC
!
!.......... read the jacobian 
            FLNAME(1:80) = ' ' 
            CFREQ(1:12) = ' ' 
            CSRC(1:5) = ' ' 
            WRITE(CFREQ,'(F12.5)') frq%FREQ(IFREQ)
            CFREQ = ADJUSTL(CFREQ) 
            WRITE(CSRC ,'(I5)') ISRC
            CSRC = ADJUSTL(CSRC) 
            FLNAME ='./scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'//TRIM(CSRC)//'.dat'
            FLNAME = ADJUSTL(FLNAME)
            OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
            READ(IUNIT,IOSTAT=IERR) NOBS_IN,NNPGL_IN, NDIM_IN, NREC_IN
            IF (IERR /= 0) GOTO 105 
            READ(IUNIT,IOSTAT=IERR)((LOBS(I,IREC),I=1,NDIM),IREC=1,NREC_IN)
            IF (IERR /= 0) GOTO 105 
            READ(IUNIT,IOSTAT=IERR)((RCJLOC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,inv%NNPGL)
            IF (IERR /= 0) GOTO 105 
            READ(IUNIT,IOSTAT=IERR)((CCJLOC(IOBS,INPGL),IOBS=1,NOBS_IN),INPGL=1,inv%NNPGL)
            IF (IERR /= 0) GOTO 105 
            IF (IERR /= 0) THEN
               WRITE(*,*) 'gnmfree: Error reading cjloc!'
               CLOSE(IUNIT)
               IERR = 3 
               GOTO 105 
            ENDIF
            CLOSE(IUNIT) 
!
!.......... calculate the residuals in (u,v,w) frame
            IF (MYNID == MASTER) THEN
               DO 201 IREC=1,rcv%NREC
                  OBSR(1,IREC) = CAZ*inv%OBS(1,IFREQ,IREC,ISRC) &
                               + SAZ*inv%OBS(2,IFREQ,IREC,ISRC)
                  OBSR(2,IREC) =-SAZ*inv%OBS(1,IFREQ,IREC,ISRC) &
                               + CAZ*inv%OBS(2,IFREQ,IREC,ISRC)
                  OBSR(3,IREC) = inv%OBS(3,IFREQ,IREC,ISRC) 
  201          CONTINUE
               CALL BPRESID4(NDIM,rcv%NREC,NDIM, inv%NORM,inv%IRESTP,  0.2,  & 
                             inv%WGHTS(1:NDIM,IFREQ,1:rcv%NREC,ISRC),         &   
                             OBSR, inv%EST(1:NDIM,IFREQ,1:rcv%NREC,ISRC),     &   
                             RESID)
            ENDIF
!
!.......... broadcast residuals
            DO 203 I=1,NDIM
               CALL MPI_BCAST(RESID(I,1:rcv%NREC),rcv%NREC,MPI_COMPLEX, & 
                              MASTER,MYSLV_COMM,MPIERR)  
  203       CONTINUE
!
!.......... dstirbuted adjoint multiply
            IOBS = 0 
            GRADB(1:inv%NA35) = 0.0
            DO 211 IREC=1,rcv%NREC 
               DO 212 I=1,NDIM
                  IF (LOBS(I,IREC)) THEN
                     IOBS = IOBS + 1
                     IF (IOBS > NOBS_IN) THEN
                        WRITE(*,*) 'Error, nobs > nobs_in'
                     ENDIF
                     DO 213 INPGL=1,inv%NNPGL
                        ILOC = inv%MYGRAD(INPGL) 
                        CARG = CMPLX(RCJLOC(IOBS,INPGL),-CCJLOC(IOBS,INPGL))
                        GRADB(ILOC) = GRADB(ILOC) + REAL(CARG*RESID(I,IREC))  
  213                CONTINUE
                  ENDIF
  212          CONTINUE
  211       CONTINUE 
!           print *, 'reducing',maxval(gradb) 
            CALL MPI_REDUCE(GRADB,TEMP,inv%NA35,MPI_REAL, MPI_SUM,MASTER, &
                            MYSLV_COMM,MPIERR)  
            IF (MYNID == MASTER) CALL SAXPY(inv%NA35,1.0,TEMP,1,GRADL,1) !stack
  200    CONTINUE !loop on sources
  101    CONTINUE !nothing for this frequency group to do
  100 CONTINUE !loop on frequencies
  105 CONTINUE !out of frequencies 
!
!.... reduce the radient
      IF (MYNID == MASTER) THEN
         CALL MPI_REDUCE(GRADL,inv%GRAD,inv%NA35,MPI_REAL,MPI_SUM, &
                         MASTER,MYHD_COMM,MPIERR)
         IF (MYID == MASTER) THEN
            WMIN = MINVAL(inv%WMASK)
            WMAX = MAXVAL(inv%WMASK)
            IF (WMIN /= 1.0 .AND. WMAX /= 1.0) THEN
               WRITE(*,*) 'funcgradh: Multiplying gradient mask...'
               CALL WTGRAD(inv%NA35,inv%WMASK, inv%GRAD)
            ENDIF
            IF (inv%NCASC > 0 .and. .not.inv%lpgrad) THEN
               WRITE(*,*) 'funcgradh: Smoothing gradient...'
               CALL AVGRAD(NGNOD,msh%NNPG, inv%NA35,msh%NNPG, &
                           NGNOD,inv%NCON,inv%NVINV,inv%NCASC, &
                           inv%MASKG,inv%MCONN,msh%IENG, inv%GRAD)  
            ENDIF
         ENDIF
      ENDIF
      IF (ALLOCATED(RCJLOC)) DEALLOCATE(RCJLOC) 
      IF (ALLOCATED(CCJLOC)) DEALLOCATE(CCJLOC) 
      IF (ALLOCATED(RESID))  DEALLOCATE(RESID) 
      IF (ALLOCATED(GRADL))  DEALLOCATE(GRADL) 
      IF (ALLOCATED(GRADB))  DEALLOCATE(GRADB) 
      IF (ALLOCATED(TEMP))   DEALLOCATE(TEMP) 
      IF (ALLOCATED(LOBS))   DEALLOCATE(LOBS) 
      IF (ALLOCATED(OBSR))   DEALLOCATE(OBSR) 
      RETURN
      END
