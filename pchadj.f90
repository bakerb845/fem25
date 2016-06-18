      SUBROUTINE PCHADJ(MYID,MYNID,MASTER,MYCOMM,  NSPACE,           &
                        LJOINT,LSURF, LPGRAD,IFREQ,ISRC, FREQ, LOBS, &
                        FMAT_DIST, RCV,INV,MID, HESSB,IERR) 
!
!     Calculates the pre-conditioner adj(J) trans(A) A J for the inversion. 
!     Here, we operate on J =-A S^{-1} [dS/dm_1 u, ... dS/dm_m u]
!     To approach this we extract the Green's functions: 
!       G =-A S^{-1}
!       G S =-A 
!       S^T G^T =-A^T
!       (LU)^T G^T =-A^T
!     and solve -(LU)^T G^T =-A^T where G^T is [ndof x n_{obs}] and contains 
!     the Greens functions for a receiver to every point in the medium. 
!     -  Ben Baker December 2012
!
!     INPUT      MEANING
!     -----      ------- 
!     IFREQ      frequency number
!     ISRC       current source for file IO
!     FREQ       current frequency for file IO
!     LJOINT     True -> this is a joint inversion, need lsurf
!     LOBS       True -> theres an observation for component,receiver 
!     LPGRAD     True -> set the block diagonal pre-conditioner
!     LSURF      if ljoint is true then this specifies surface or body wave
!     MASTER     master ID
!     MYCOMM     communicator for MYNID 
!     MYID       ID on MPI_COMM_WORLD
!     MYNID      ID in frequency group
!     NDIM       number of components in solution
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     HESSB      updated local hessian on head process
!     IERR       error flag
! 
!.... variable declarations  
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_struc.h'
      INCLUDE 'fwd_struc.h'
      TYPE (CMUMPS_STRUC) MID
      TYPE (RECV_INFO) RCV
      TYPE (INV_INFO) INV
      COMPLEX*8, INTENT(IN) :: FMAT_DIST(NSPACE)
      REAL*8, INTENT(IN) :: FREQ
      INTEGER*4, INTENT(IN) :: NSPACE, MYID,MYNID,MASTER,MYCOMM, IFREQ, ISRC
      LOGICAL*4, INTENT(IN) :: LOBS(NDIM,*), LPGRAD, LJOINT, LSURF 
      COMPLEX*8, INTENT(INOUT) :: HESSB(*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*8, ALLOCATABLE :: JACL(:,:), HMATLOC(:,:), RHS_SAVE(:),  GRNSFNS(:)
      REAL*4, ALLOCATABLE :: RCJLOC(:,:), CCJLOC(:,:), HBUFF(:), HBUFF2(:) 
      CHARACTER(80) FLNAME 
      CHARACTER(12) CFREQ
      CHARACTER(5) CID, CSRC 
      INTEGER*4  ISAVE9, ISAVE21, NRHSS,          &
                 NRHS, NRHSF, I, K, IVINV, JVINV, IRHS, KRHS, IREC, I1,I2,         &
                 INZ, INPINV, INPGL,  IDOF, IHDEST, JDOF, JLOC, LRHS,                &
                 NGFN, NWORK, MPIERR, INDX, IUNIT 
      LOGICAL*4 LRTRN
      COMPLEX*8, PARAMETER :: CZERO = CMPLX(0.0,0.0) 
      COMPLEX*8, PARAMETER :: CONE  = CMPLX(1.0,0.0)
      PARAMETER(IUNIT = 68) 
!
!----------------------------------------------------------------------------------------!
!
!.... check for early break
      IF (MYNID == MASTER) THEN
         DO 1 IREC=1,rcv%NREC
            DO 2 I=1,NDIM
               IF (LOBS(I,IREC)) THEN
                  LRTRN = .FALSE.
                  GOTO 100
               ENDIF
    2       CONTINUE
    1    CONTINUE
         LRTRN = .TRUE.
  100    CONTINUE
      ENDIF
      CALL MPI_BCAST(LRTRN,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      IF (LRTRN) THEN
         IF (MYNID == MASTER) WRITE(*,*) 'pchadj: Early return; no observations!'
         RETURN
      ENDIF 
!.... set the sparse RHS of trans(A)_[ndof,ncomp*nrec] 
      IERR = 0
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) THEN
            WRITE(*,*) 
            WRITE(*,*) 'pchadj: Setting sparse RHS...'
         ENDIF
         ISAVE9  = MID%ICNTL(9)  !save transpose problem
         ISAVE21 = mid%ICNTL(21) !RHS dense (0) or destributed (1)
         NRHSS = mid%NRHS
         IF (ISAVE21 == 0) THEN  
            NWORK = NRHSS*MID%LRHS 
            ALLOCATE(RHS_SAVE(NWORK))
            RHS_SAVE(1:NWORK) = MID%RHS(1:NWORK)
            DEALLOCATE(MID%RHS)
         ENDIF
         IF (MYID == MASTER) &
         WRITE(*,*) 'pchadj: Max number of RHSs:',rcv%NREC*NDIM 
!
!....... number of RHSs
         NRHS = rcv%NREC*NDIM 
!        K = 0 
!        DO 1 IREC=1,rcv%NREC
!           DO 2 I=1,NDIM
!              IF (LOBS(I,IREC)) K = K + 1
!   2       CONTINUE !loop on components
!   1    CONTINUE  
!        NRHS = K 
         mid%NRHS = NRHS !number of RHSs
         ALLOCATE(mid%RHS(mid%LRHS*mid%NRHS))
         mid%RHS(1:mid%LRHS*mid%NRHS) = CMPLX(0.0,0.0)
         K = 0 
         DO 3 IREC=1,rcv%NREC
            DO 4 I=1,NDIM
!              IF (LOBS(I,IREC)) THEN
                  K = K + 1
                  IDOF = rcv%MRDOF(I,IREC)
                  IF (IDOF > 0) THEN
                     INDX = (K - 1)*mid%LRHS + IDOF
                     IF (INDX < 0 .OR. INDX > mid%LRHS*mid%NRHS) THEN
                        WRITE(*,*) 'pchadj: Location error!',INDX
                        IERR = 1
                        RETURN
                     ENDIF
                     !mid%RHS(INDX) =-CONE !remember the negative
                     IF (     LSURF .AND. inv%LCOVD_SRF .OR. &
                         .NOT.LSURF .AND. inv%LCOVD_BDY) THEN  
                        IF (inv%WGHTS(I,IFREQ,IREC,ISRC) > 0.0) THEN
                           mid%RHS(INDX) =-CMPLX(1.0,0.0)
                        ELSE
                           mid%RHS(INDX) =-CMPLX(0.0,0.0)
                        ENDIF
                     ELSE
                        mid%RHS(INDX) =-CMPLX(inv%WGHTS(I,IFREQ,IREC,ISRC),0.0) 
                     ENDIF
                  ELSE
                     WRITE(*,*) 'pchadj: Receiver not collocated to DOF!'
                  ENDIF
!              ENDIF !end check on observation
    4       CONTINUE 
    3    CONTINUE 
         NRHS = K 
         IF (MYID == MASTER) WRITE(*,*) 'pchadj: Number of RHSs:',NRHS
         mid%NRHS = NRHS   !new number of RHSs
!        mid%NZ_RHS = K    !new number of non-zeros in RHS
         mid%ICNTL(9)  = 0 !want to solve transpose problem
!        mid%ICNTL(20) = 1 !sparse RHS
!        mid%ICNTL(21) = 0 !want centralized RHS
!        ALLOCATE(MID%RHS(mid%LRHS*mid%NRHS))
         LRHS = mid%LRHS 
      ENDIF
      CALL MPI_BCAST(LRHS   ,1,MPI_INTEGER,MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NRHS   ,1,MPI_INTEGER,MASTER,MYCOMM,MPIERR)
!
!.... extract Green's functions in m_{rec} columns, now and [ndof x m_{rec}]; G^T 
      IF (MYID == MASTER) WRITE(*,*) 'pchadj: Calculating receiver Greens functions...'
      mid%JOB = 3 
      CALL CMUMPS(MID) 
      IF (mid%INFO(1) < 0) THEN
         IF (MYID == MASTER) WRITE(*,*) 'pchadj: An error occurred in cmumps'
         IERR = 1 
         RETURN
      ENDIF 
      IF (LPGRAD) THEN
         ALLOCATE(HBUFF(inv%NHSIZE))
         HBUFF(1:inv%NHSIZE) = 0.0 !CZERO 
         ALLOCATE(HMATLOC(inv%NVINV,inv%NVINV))
         HMATLOC(1:inv%NVINV,1:inv%NVINV) = CZERO
      ENDIF
      NGFN = NRHS*LRHS 
      ALLOCATE(GRNSFNS(NGFN),STAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'pchadj: Error allocating space for Greens function buffer!'
         WRITE(*,*) 'pchadj: The error occured on proess:',MYID
         RETURN
      ENDIF
      IF (MYNID == MASTER) THEN
         IF (MYID == MASTER) WRITE(*,*) 'pchadj: Copying Greens functions...'
         !GRNSFNS(1:NGFN) = CONJG(mid%RHS(1:NGFN)) !copy solution back to buffer
         CALL CCOPY(NGFN,mid%RHS,1,GRNSFNS,1) !copy G^T onto buffer 
         !CALL CCOPY(NGFN,CONJG(mid%RHS),1,GRNSFNS,1) !copy solution to buffer 
         IF (MYID == MASTER) WRITE(*,*) 'pchadj: Restoring data structures...' 
         mid%NRHS = NRHSS
         mid%ICNTL(9)  = ISAVE9
         mid%ICNTL(21) = ISAVE21
         IF (ISAVE21 == 0) THEN  
            NWORK = mid%NRHS*mid%LRHS 
            ALLOCATE(mid%RHS(NWORK))
            CALL CCOPY(NWORK,RHS_SAVE,1,mid%RHS,1) !copy saved RHS back to MUMPS RHS
            DEALLOCATE(RHS_SAVE)
         ENDIF
         IF (MYID == MASTER) WRITE(*,*) 'pchadj: Broadcasting Greens functions...'
      ENDIF 
      CALL MPI_BCAST(GRNSFNS,NGFN,MPI_COMPLEX, MASTER,MYCOMM,MPIERR)
      IF (MYID == MASTER) THEN
         IF (LPGRAD.AND.inv%LGNEWT) THEN 
            WRITE(*,*) 'pchadj: Updating preconditioner and calculating Jacobian...'
         ELSEIF (inv%LGNEWT .AND..NOT.LPGRAD) THEN
            WRITE(*,*) 'pchadj: Calculating Jacobian...'
         ELSEIF (LPGRAD .AND..NOT.inv%LGNEWT) THEN
            WRITE(*,*) 'pchadj: Updating preconditioner...'
         ENDIF
      ENDIF
!
!.... all process have Greens functions, now do matrix matrix multiply (columns)
      !ALLOCATE(ADJAC(NRHS,inv%NVINV)) !holds adj(J) 
      ALLOCATE(JACL(NRHS,inv%NVINV))  !holds Jacobian (locally)
      IF (inv%LGNEWT) THEN
         NRHSF = rcv%NREC*NDIM
         ALLOCATE(RCJLOC(NRHSF,inv%NNPGL),STAT=IERR) ![n_{obs} x m] 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error setting space in rcjloc on process',MYID
            RETURN
         ENDIF
         ALLOCATE(CCJLOC(NRHSF,inv%NNPGL),STAT=IERR) ![n_{obs} x m]
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error setting space in ccjloc on process',MYID
            RETURN
         ENDIF
      ENDIF
      DO 500 INPGL=1,inv%NNPGL,inv%NVINV
!
!....... loop on ndim*nrec rows of J
         INPINV = inv%MYGRAD(INPGL) 
         IRHS = 0
         KRHS = 0 
         DO 501 IREC=1,rcv%NREC
            DO 502 I=1,NDIM  
               KRHS = KRHS + 1
!              IF (LOBS(I,IREC)) THEN 
                  IRHS = IRHS + 1
!
!................ loop on inversion variables (columns of J)
                  DO 503 IVINV=1,inv%NVINV
                     INDX = INPGL + IVINV - 1
                     I1 = inv%JCSC_FDIST(INDX)
                     I2 = inv%JCSC_FDIST(INDX+1) - 1
                     IF (I1 < 0 .OR. I2 > NSPACE) THEN
                        WRITE(*,*) 'pchadj: Weird problem',I1,I2
                        IERR = 1
                        RETURN
                     ENDIF
                     !ADJAC(IVINV,IRHS) = CZERO
                     JACL(IRHS,IVINV) = CZERO
                     DO 504 INZ=I1,I2 !multiply row of G^T with column of F
                        JDOF = inv%ICSC_FDIST(INZ) !points to the row of F
                        IF (JDOF == 0) GOTO 504 !this would be weird
                        JLOC = (IRHS - 1)*LRHS + JDOF 
                        IF (JLOC < 1 .OR. JLOC > NGFN) THEN
                           WRITE(*,*) 'pchadj: Sizing error!', JLOC,NGFN
                           IERR = 1
                           RETURN 
                        ENDIF
!                       ADJAC(IVINV,IRHS) = ADJAC(IVINV,IRHS)                        &
!                                         + CONJG(FMAT_DIST(INZ))*GRNSFNS(JLOC)
                        JACL(IRHS,IVINV) = JACL(IRHS,IVINV)  & 
                                         + GRNSFNS(JLOC)*FMAT_DIST(INZ) 
  504                CONTINUE !loop on CSC rows
  503             CONTINUE !loop on inversion variables
!
!................ save Jacobian?
                  IF (inv%LGNEWT) THEN
                     DO 505 IVINV=1,inv%NVINV !copy inversion variables in block
                        INDX = INPGL + IVINV - 1
                        RCJLOC(KRHS,INDX) = REAL(JACL(IRHS,IVINV)) !REAL(ADJAC(IVINV,IRHS)) !transpose
                        CCJLOC(KRHS,INDX) = IMAG(JACL(IRHS,IVINV)) !-IMAG(ADJAC(IVINV,IRHS)) !conjuate transpose
  505                CONTINUE
                  ENDIF !end check on jacobian save
!              ELSE
!                 IF (inv%LGNEWT) THEN
!                    DO 506 IVINV=1,inv%NVINV !zero row
!                       INDX = INPGL + IVINV - 1
!                       RCJLOC(KRHS,INDX) = 0.0
!                       CCJLOC(KRHS,INDX) = 0.0 
! 506                CONTINUE 
!                 ENDIF
!              ENDIF !end check on  if observation is defined
  502       CONTINUE !loop on components
  501    CONTINUE !loop on receivers
!
!....... update preconditioner?
         IF (LPGRAD) THEN
!.......... computation time of both is indistinguishable
!           CALL CGEMM('N','C',inv%NVINV,inv%NVINV,NRHS,CONE, ADJAC,inv%NVINV, ADJAC,inv%NVINV, &
!                      CZERO,HMATLOC,inv%NVINV)
            CALL CGEMM('C','N',inv%NVINV,inv%NVINV,NRHS,CONE, JACL, NRHS, JACL,NRHS, &
                       CZERO,HMATLOC,inv%NVINV) 
            !HMATLOC = MATMUL(ADJAC,TRANSPOSE(CONJG(ADJAC))) 
            !HMATLOC = MATMUL(TRANSPOSE(CONJG(JACL)),JACL) 
            INPINV = inv%MYGRAD(INPGL)       !get the right block 
            IHDEST = (INPINV - 1)*inv%NVINV  !now offset the banded ordering
            DO 520 IVINV=1,inv%NVINV
               DO 521 JVINV=1,inv%NVINV
                  IHDEST = IHDEST + 1 
                  if (ihdest < 0) print *, 'err'
                  if (ihdest > inv%nhsize) print *, 'err', ihdest,inv%nhsize
!                 IF (IVINV == JVINV) THEN 
!                    HBUFF(IHDEST) = HBUFF(IHDEST) + CMPLX(REAL(HMATLOC(IVINV,JVINV)),0.0)
!                 ELSE
!                    HBUFF(IHDEST) = HBUFF(IHDEST) + HMATLOC(IVINV,JVINV)
!                 ENDIF
                  HBUFF(IHDEST) = HBUFF(IHDEST) + REAL(HMATLOC(IVINV,JVINV)) 
  521          CONTINUE
  520       CONTINUE
         ENDIF !end check on pre-conditioner update
  500 CONTINUE
!
!.... free space
      IF (ALLOCATED(GRNSFNS)) DEALLOCATE(GRNSFNS)
      !IF (ALLOCATED(ADJAC))   DEALLOCATE(ADJAC) 
      IF (ALLOCATED(JACL))    DEALLOCATE(JACL)
      IF (ALLOCATED(HMATLOC)) DEALLOCATE(HMATLOC)
!
!.... update preconditioner 
      IF (LPGRAD) THEN
         IF (MYNID == MASTER) THEN
            ALLOCATE(HBUFF2(inv%NHSIZE))
         ELSE
            ALLOCATE(HBUFF2(1))
         ENDIF
         NWORK = inv%NHSIZE
!        CALL MPI_REDUCE(HBUFF,HBUFF2,NWORK,MPI_COMPLEX,MPI_SUM, MASTER,MYCOMM,MPIERR) 
         CALL MPI_REDUCE(HBUFF,HBUFF2,NWORK,MPI_REAL,MPI_SUM, MASTER,MYCOMM,MPIERR) 
         DEALLOCATE(HBUFF)
         IF (MYNID == MASTER) THEN
!           CALL CAXPY(inv%NHSIZE,CONE,HBUFF2,1,HESSB,1) !hessl = hessl + hbuff2
            CALL SAXPY(inv%NHSIZE,1.0,HBUFF2,1,HESSB,1) !hessl = hessl + hbuff2
            !print *, maxval(cabs(hbuff2)),minval(cabs(hbuff2))
            !print *, maxloc(cabs(hbuff2))
            !print *, minloc(cabs(hbuff2))
         ENDIF
         IF (ALLOCATED(HBUFF2)) DEALLOCATE(HBUFF2)
      ENDIF
!
!.... dump the local jacobian
      IF (inv%LGNEWT .OR. LJOINT) THEN
         IF (MYID == MASTER) WRITE(*,*) 'pchadj: Writing local Jacobians to disk...'
         FLNAME(1:80) = ' '
         CFREQ(1:12) = ' ' 
         CID(1:5) = ' '
         CSRC(1:5) = ' '
         WRITE(CFREQ,'(F12.5)') FREQ
         WRITE(CID ,'(I5)') MYNID
         WRITE(CSRC,'(I5)') ISRC
         CFREQ = ADJUSTL(CFREQ)
         CID   = ADJUSTL(CID) 
         CSRC  = ADJUSTL(CSRC)
         IF (LJOINT) THEN
            IF (LSURF) THEN 
               FLNAME = './scratch/jac_surf_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                        TRIM(CSRC)//'.dat'
            ELSE
               FLNAME = './scratch/jac_body_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                        TRIM(CSRC)//'.dat'
            ENDIF
         ELSE
            FLNAME = './scratch/jac_'//TRIM(CID)//'-'//TRIM(CFREQ)//'-'// &
                     TRIM(CSRC)//'.dat'
         ENDIF
         FLNAME = ADJUSTL(FLNAME)
         OPEN(FILE=TRIM(FLNAME),UNIT=IUNIT,STATUS='REPLACE',   &
              FORM='UNFORMATTED',IOSTAT=IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error opening file:',TRIM(FLNAME),IERR
            RETURN
         ENDIF
         WRITE(IUNIT,IOSTAT=IERR) NRHSF,inv%NNPGL,NDIM,rcv%NREC
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error 1',IERR
            RETURN
         ENDIF
         WRITE(IUNIT,IOSTAT=IERR) ((LOBS(I,IREC),I=1,NDIM),IREC=1,rcv%NREC)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error 2',IERR
            RETURN
         ENDIF
         WRITE(IUNIT,IOSTAT=IERR) ((RCJLOC(IRHS,INPGL),IRHS=1,NRHSF),INPGL=1,inv%NNPGL) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error 3',IERR 
            RETURN
         ENDIF
         WRITE(IUNIT,IOSTAT=IERR) ((CCJLOC(IRHS,INPGL),IRHS=1,NRHSF),INPGL=1,inv%NNPGL)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error 4',IERR
            RETURN
         ENDIF
         WRITE(IUNIT,IOSTAT=IERR) (inv%MYGRAD(INPGL),INPGL=1,inv%NNPGL)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'pchadj: Error 5',IERR
            RETURN
         ENDIF
         CLOSE(IUNIT) 
      ENDIF  
      IF (ALLOCATED(RCJLOC)) DEALLOCATE(RCJLOC) 
      IF (ALLOCATED(CCJLOC)) DEALLOCATE(CCJLOC)
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE CADJLOC(MVINV, NVINV,NRHS, NZ_FDIST,NNPGL,NNZCOL, LCONJG, INPGL, &
                         JCSC_FDIST, FMAT_DIST,GRNSLOC, ADJJLOC) 
!
!     This routine handles the calculation of adj(F) adj(S^{-1}) A^T where 
!     adj(S^{-1}) A^T has already been solved for and stored in GRNSLOC.  Note, since
!     adj(F) = dS/dm_i u is sparse we truncate GRNSLOC to match this.   
!
!     INPUT      MEANING
!     -----      ------- 
!     INPGL      local inversion variable we are working on 
!     FMAT_DIST  distributed dS/dm_i 
!     GRNSLOC    pertinent Greens functions (likely not conjugated) for each receiver 
!     JCSC_FDIST CSC column pointer for distributed F matrix 
!     LCONJG     True -> slave process conjugates the Greens functions (likely True)
!     MVINV      leading dimension for  
!     NNZCOL     number of non-zeros in the columns of adj(F) for this node 
!     NNPGL      number of local inversion variables 
!     NZ_FDIST   number of non-zeros in FMAT_DIST
!     NRHS       number of RHSs held in GRNSLOC (likely nrec*ndim) 
!     NVINV      number of inversion variables at this node
!
!     OUTPUT     MEANING
!     ------     ------- 
!     ADJJLOC    for this nodal point the [nvinv x ndim*nrec] matrix adj(J)
!
!.... variable declarations 
      COMPLEX*8 FMAT_DIST(NZ_FDIST), GRNSLOC(NRHS*NNZCOL)  
      INTEGER*4, INTENT(IN) :: JCSC_FDIST(NNPGL+1), MVINV,NRHS,NVINV,  &
                               NZ_FDIST,NNPGL,NNZCOL, INPGL 
      LOGICAL*4, INTENT(IN) :: LCONJG
      COMPLEX*8, INTENT(OUT) :: ADJJLOC(MVINV,*) 
!.... local variables
      COMPLEX*8 CZERO
      INTEGER*4 IRHS,IVINV,IDOFL,INZHL, I1,I2 
      PARAMETER(CZERO = CMPLX(0.0,0.0))  
!
!----------------------------------------------------------------------------------------!
!
!.... calculate the ndof x ndof inner produts for each (ivinv,irhs) pair 
      DO 1 IRHS=1,NRHS !loop on columns
         DO 2 IVINV=1,NVINV !Loop on rows
            ADJJLOC(IVINV,IRHS) = CZERO 
            I1 = JCSC_FDIST(INPGL+IVINV-1)
            I2 = JCSC_FDIST(INPGL+IVINV) - 1 
            IDOFL = (IRHS - 1)*NRHS 
!
!.......... inner product on DOFs, conjg is b/c B was never conjugated
            DO 3 INZHL=I1,I2
               IDOFL = IDOFL + 1 
               IF (LCONJG) THEN
                  ADJJLOC(IVINV,IRHS) = ADJJLOC(IVINV,IRHS) &
                                      + FMAT_DIST(INZHL)*CONJG(GRNSLOC(IDOFL))
               ELSE 
                  ADJJLOC(IVINV,IRHS) = ADJJLOC(IVINV,IRHS) &
                                      + FMAT_DIST(INZHL)*CONJG(GRNSLOC(IDOFL))
               ENDIF
    3       CONTINUE  
    2    CONTINUE
    1 CONTINUE !loop on columns
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GETHBUFF(MVINV,NVINV,HMATLOC, HBUFF) 
!     Simple utility to convert the HMATLOC -> column major ordering.  Master process
!     can later unpack it.
      COMPLEX*8, INTENT(IN) :: HMATLOC(MVINV,*) 
      INTEGER*4, INTENT(IN) :: MVINV,NVINV 
      REAL*4, INTENT(OUT) :: HBUFF(NVINV**2) 
! 
!----------------------------------------------------------------------------------------!
!
      K = 0 
      DO 1 IVINV=1,NVINV
         DO 2 JVINV=1,NVINV  
            K = K + 1 
            HBUFF(K) = REAL(HMATLOC(IVINV,JVINV)) 
    2    CONTINUE !loop on rows 
    1 CONTINUE !loop on columns 
      RETURN
      END
!
!.... adj(J)_{loc} J_{loc}
!     HMATLOC = MATMUL(ADJJLOC,TRANSPOSE(CONJG(ADJJLOC)))
!
!.... reshape into row major buffer, done w/ imaginary part
!     K = 0 
!     DO 304 IVINV=1,NVINV
!        DO 305 JVINV=1,NVINV  
!           K = K + 1 
!                       HBUFF(K) = REAL(HMATLOC(IVINV,JVINV)) 
! 305                CONTINUE !loop on rows 
! 304             CONTINUE !loop on columns 

