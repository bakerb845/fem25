      SUBROUTINE CFMAT_DIST(NWORK,N, FREQ,PY, WAVE, MSH,INV, FMAT_DIST,IERR) 
!
!     Calculates each processes columns of the matrix [dS/dm1 u, dS/dm2 u,...,dS/dmn u] 
!     Realize that the columns are distributed and sparse
!
!     INPUT      MEANING
!     -----      -------
!     FREQ       current frequency
!     N          size of WAVE = NDOF
!     NWORK      size of FMAT_DIST which is local number of non-zeros NZ_FDIST
!     PY         apparent slowness
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     FMAT_DIST  processes' contributed to FMAT matrix
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH
      TYPE (INV_INFO)  INV
      COMPLEX*8, INTENT(IN) :: WAVE(N)
      REAL*8, INTENT(IN) :: FREQ,PY
      INTEGER*4, INTENT(IN) :: NWORK,N
      COMPLEX*8, INTENT(OUT) :: FMAT_DIST(NWORK) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*8, ALLOCATABLE :: FORCEV(:,:)
      REAL*8, ALLOCATABLE :: SHL(:,:,:,:), XIGLL(:,:), ETAGLL(:,:) 
      COMPLEX*8 CZERO
      REAL*8 OMEGA, TWOPI
      INTEGER*4 NINTX,NINTZ,NEE,NNZROW, INPG,INPGL, I1,I2,K1,K2, INPINV,IVINV,  &
                INDX
      PARAMETER(CZERO = CMPLX(0.0,0.0))
      PARAMETER(TWOPI = 6.2831853071795862D0)
!
!----------------------------------------------------------------------------------------!
!
!.... first error check
      IERR = 0 
      IF (inv%CINVTYPE /= 'PP' .AND. inv%CINVTYPE /= 'pp' .AND. &
          inv%CINVTYPE /= 'SS' .AND. inv%CINVTYPE /= 'ss' .AND. &
          inv%CINVTYPE /= 'PS' .AND. inv%CINVTYPE /= 'ps') THEN
         WRITE(*,*) 'cfmat_dist: Error anisotropy not programmed!'
         IERR = 1 
         RETURN
      ENDIF
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
      ALLOCATE(FORCEV(inv%MBUFRHS,inv%NVINV)) 
!
!.... loop on local Hessian points 
      DO 100 INPGL=1,inv%NNPGL,inv%NVINV 
!
!....... locate anchor node to invert at
         DO 101 INPG=1,msh%NNPG
            INPINV = inv%MASKG(INPG)
            IF ((INPINV - 1)*inv%NVINV + 1 == inv%MYGRAD(INPGL)) GOTO 110
  101    CONTINUE  
         WRITE(*,*) 'cfmat_dist: Cant find anchor node!'
         IERR = 1
         GOTO 730
  110    CONTINUE 
!
!....... null out fmat_dist
         I1 = inv%JCSC_FDIST(INPGL)
         I2 = inv%JCSC_FDIST(INPGL+1) - 1
         FMAT_DIST(I1:I2) = CZERO 
         NNZROW = I2 - I1 + 1
         DO 30 IVINV=2,inv%NVINV
            INDX = INPGL + IVINV - 1
            IF (INDX > inv%NNPGL) THEN
               WRITE(*,*) 'cfmat_dist: INDX out of bounds!'
               IERR = 1
               GOTO 730
            ENDIF 
            K1 = inv%JCSC_FDIST(INDX) 
            K2 = inv%JCSC_FDIST(INDX+1) - 1
            IF (K2 - K1 + 1 /= NNZROW) THEN
               WRITE(*,*) 'cfmat_dist: Rows are inconsistent size!'
               IERR = 1
               GOTO 730
            ENDIF
            FMAT_DIST(K1:K2) = CZERO 
   30    CONTINUE 
!
!....... calculate [dS/dm_i u] for this grid point 
         CALL VFORCE25(NNZROW,msh%NEN,NINTX,NINTZ,msh%NNPG,NDIM,NGNOD,                   &
                              msh%NEN,NINTX,NINTZ,msh%NNPG,     NGNOD,                   &
                       inv%NNPGL,inv%NCON,inv%NVINV, inv%NZ_FDIST,msh%NDOF,msh%NELEM,    &
                       inv%CINVTYPE,msh%LISISO, INPG,INPGL,                              &
                       inv%JCSC_FDIST,inv%ICSC_FDIST, msh%IENG,inv%MCONN,msh%LM,         &
                       OMEGA,PY, inv%VPVS, msh%XLOCS,msh%ZLOCS,msh%DENS,                 &
                       msh%ECOEFF(1:msh%NNPG,1:2), XIGLL,ETAGLL,SHL, inv%ELEM_WTS,WAVE,  &
                       FORCEV(1:NNZROW,1:inv%NVINV), IERR) 
         IF (IERR /= 0) THEN
            WRITE(*,*) 'cfmat_dist: Error calling vforce25'
            GOTO 730
         ENDIF
!
!....... copy
         DO 31 IVINV=1,inv%NVINV
            INDX = INPGL + IVINV - 1
            K1 = inv%JCSC_FDIST(INDX)
            K2 = inv%JCSC_FDIST(INDX+1) - 1
            NNZROW = K2 - K1 + 1
            IF (NNZROW > inv%MBUFRHS) THEN
               WRITE(*,*) 'cfmat_dist: Sizing error'
               IERR = 1
               RETURN
            ENDIF
            CALL CCOPY(NNZROW,FORCEV(1:NNZROW,IVINV),1,FMAT_DIST(K1:K2),1) 
   31    CONTINUE
  100 CONTINUE !loop on anchor nodes
  730 CONTINUE !break ahead for an error
!
!.... free space
      DEALLOCATE(XIGLL)
      DEALLOCATE(ETAGLL)
      DEALLOCATE(SHL)
      DEALLOCATE(FORCEV)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      REAL*8 FUNCTION AREA_GRNS(MGNOD,NNPG,NGNOD,NCON,  MCONN,IENG, &
                                XLOCS,ZLOCS) 
!
!     Returns the area of the patch with which we calculate each virtual force.  
!     Note each element is simple, a triangle or quadrilateral, so the total area 
!     is the sum of all areas of the quadrilaterals and simplexes - B. Baker
!
!     INPUT     MEANING
!     -----     ------- 
!     IENG      global anchor node IEN matrix
!     MGNOD     leading dimension
!     MCONN     holds current anchor node's element connections
!     NCON      max number of connections 
!     NGNOD     number of anchor nodes on element
!     NNPG      number of anchor nodes in mesh
!     XLOCS     x locations of anchor nodes
!     ZLOCS     z locations of anchor nodes 
!
!     OUTPUT    MEANING
!     ------    ------- 
!     AREA_GRNS area of patch
!
!.... variable declarations
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG) 
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MCONN(NCON), MGNOD, NGNOD, NNPG, NCON 
!.... local variables
      REAL*8 RLOC, X1,X2, Z1,Z2 
      INTEGER*4 IMOD, ICON, IELEM, IA1, IA2, JNPG1, JNPG2   
      LOGICAL*4 LQUAD 
!
!----------------------------------------------------------------------------------------!
!
!.... loop on anchor node connectivity 
      AREA_GRNS = 0.D0
      DO 1 ICON=1,NCON 
         IELEM = MCONN(ICON)
         IF (IELEM <= 0) GOTO 100
!
!....... calculate area of each neighboring element
         RLOC = 0.D0
         IMOD = 5 
         LQUAD = .TRUE.  !likely a quad
         IF (MINVAL(IENG(1:NGNOD,IELEM)) == 0) THEN
            LQUAD = .FALSE.
            IMOD = 4
         ENDIF
         DO 2 IA1=1,NGNOD !loop on anchor nodes
            JNPG1 = IENG(IA1,IELEM)
            IF (JNPG1 == 0) GOTO 200 !break for triangle
            IA2 = MOD(IA1+1,IMOD) 
            IF (IA2 == 0) IA2 = 1  
            JNPG2 = IENG(IA2,IELEM)
            X1 = XLOCS(JNPG1)
            X2 = XLOCS(JNPG2)
            Z1 = ZLOCS(JNPG1) 
            Z2 = ZLOCS(JNPG2) 
            RLOC = RLOC + X1*Z2 - X2*Z1 
    2    CONTINUE
  200    CONTINUE !break ahead, out of anchor nodes
         IF (RLOC < 0.D0) THEN
            WRITE(*,*) 'area_grns: Warning your mesh is wrong'
            RLOC =-RLOC
         ENDIF 
         RLOC = RLOC*0.5D0 !only going one way here
         AREA_GRNS = AREA_GRNS + RLOC 
    1 CONTINUE 
  100 CONTINUE !break ahead, out of elements 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
!     SUBROUTINE GHESSTRCT(MYID,MASTER,MYCOMM, MDIM,MEN, NNPG,   &
!                          NEN,NDIM, LM, IERR)
!
!     Sets the FMAT_DIST structure.  The idea is that when calculating a preconditioner 
!     that can calculate adj(F) adj(S^{-1}) A^T = adj(F) B^T where the B^T is the matrix 
!     of Green's functions for each receiver.  B^T will be centralized on the host.  Then,
!     each process will request a selection of B^T corresponding the sparse adj(F) 
!     structure for the inversion point and create a subset of adj(J)_inpinv which 
!     is of dimension [nvinv,ndim*nrec].  The process can then calculate the 
!     [nvinv x nvinv] adj(Jloc) Jloc matrix matrix multiply to obtain a block diagonal 
!     segment of the precodntioner.  -  Ben Baker December 2012
!
!     INPUT      MEANING
!     -----      ------- 
!     LM         maps element -> global DOF
!     MASTER     master process ID
!     MDIM       leading dimension for LM
!     MEN        leading dimension for LM
!     MYCOMM     MPI communicator
!     MYID       process ID on group
!     NDIM       number of componenets in solution
!     NEN        number of element nodes
!     NNPG       number of anchor noes in mesh
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!
!.... variable declarations
!     IMPLICIT NONE
!     INCLUDE 'mpif.h'
!     INCLUDE 'mesh_inv.inc'
!     INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), MYID,MASTER, MYCOMM,    &   
!                              MDIM,MEN,NNPG, NEN,NDIM 
!     INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
!     INTEGER*4, ALLOCATABLE :: IROW(:), IWORK(:)  
!     LOGICAL*4, ALLOCATABLE :: LINIT(:) 
!     INTEGER*4 STAT(MPI_STATUS_SIZE), INPINV, ICON, INZCOL, IELEM, IR, IROW1, IROW2,    &
!               IDOF, IAE, I, INZFL, IVINV, INPG, INPGL, ICOMBIG, INPG0, ILOC,         &
!               NEE, NNZCOL, NROW, NWORK, NPROCS, MYTAG, MYDEST, MYSRC, MPIERR 
!     INTEGER*4 ICNPGL
!     LOGICAL*4 LDEBUG 
!     PARAMETER(LDEBUG = .FALSE.) 
!
!----------------------------------------------------------------------------------------!
!
!.... may not be much for a poorly partition mesh  
!     IF (NNPGL <= 0) THEN
!        NNPGL = 0
!        NZ_FDIST = 0
!        GOTO 80 !really nothing to do
!     ENDIF
!
!.... set space
!     NEE = NDIM*NEN
!     NWORK = NNPG*NCON*NEE
!     ALLOCATE(IROW(NCON*NEE + 1)) !holds non-zeros in row
!     ALLOCATE(JCSC_FDIST(NNPGL+1)) !CSC column pointer
!     ALLOCATE(ICSC_FDIST(NWORK))   !CSC row pointer
!     JCSC_FDIST(1:NNPGL+1) = 0
!     ICSC_FDIST(1:NWORK) = 0
!     ICOMBIG = 0
!
!.... loop on local points in inversion and generate a CSC structure
!     INZFL = 0
!     INPGL = 1 
!     JCSC_FDIST(1) = 1 
!     DO 101 INPG=1,NNPG
!        INPINV = MASKG(INPG)
!        IF (INPINV > 0) THEN !node is in inversion
!           IF (IGPART(INPINV) == MYID + 1) THEN !in my partition 
!
!............. investigate anchor nodes connectivity in column 
!              INZCOL = 0
!              DO 102 ICON=1,NCON
!                 IELEM = MCONN(INPG,ICON)
!                 IF (IELEM > 0) THEN
!                    DO 103 IAE=1,NEN
!                       DO 104 I=1,NDIM
!                          IDOF = LM(I,IAE,IELEM)
!                          IF (IDOF > 0) THEN !non-zero row? 
!                             INZCOL = INZCOL + 1
!                             IROW(INZCOL) = IDOF
!                          ENDIF !end check on dof
! 104                   CONTINUE !loop on components
! 103                CONTINUE !loop on element nodes
!                 ELSE
!                    GOTO 110
!                 ENDIF !end check on connected element
! 102          CONTINUE !loop on connectivity
! 110          CONTINUE !break ahead 
!              NNZCOL = INZCOL
!              CALL ISHELL1(NNZCOL,IROW) !sort 
!              IROW(INZCOL+1) = 0 !cheapy way to make next step work 
!
!............. remove repeats and dump to row pointer
!              DO 105 INZCOL=1,NNZCOL
!                 IF (IROW(INZCOL+1) /= IROW(INZCOL)) THEN
!                    INZFL = INZFL + 1
!                    IF (INZFL > NWORK) THEN
!                       WRITE(*,*) 'ghesstrct: ICSC_FDIST too small!'
!                       IERR = 1
!                       RETURN
!                    ENDIF
!                    ICSC_FDIST(INZFL) = IROW(INZCOL)
!                 ENDIF
! 105          CONTINUE
!              JCSC_FDIST(INPGL+1) = INZFL + 1
!
!............. copy structure forward in inversion variables
!              IROW1 = JCSC_FDIST(INPGL)
!              IROW2 = JCSC_FDIST(INPGL+1) - 1
!              ICOMBIG = MAX(ICOMBIG,IROW2 - IROW1 + 1)
!              DO 106 IVINV=2,NVINV
!                 DO 107 IR=IROW1,IROW2
!                    INZFL = INZFL + 1
!                    IF (INZFL > NWORK) THEN
!                       WRITE(*,*) 'ghesstrct: ICSC_FDIST too small!'
!                       IERR = 1
!                       RETURN
!                    ENDIF
!                    ICSC_FDIST(INZFL) = ICSC_FDIST(IR)
! 107            CONTINUE !loop on rows 
!                INPGL = INPGL + 1
!                JCSC_FDIST(INPGL+1) = INZFL + 1
! 106          CONTINUE !loop on inversion 
!              INPGL = INPGL + 1
!           ENDIF !end check on partition 
!        ENDIF !end check on node in inversion
! 101 CONTINUE
!     NZ_FDIST = INZFL
!     DEALLOCATE(IROW) 
!
!.... fidelity check
!     INPGL = INPGL - 1 
!     IF (INPGL /= NNPGL) THEN 
!        WRITE(*,*) 'ghesstrct: Error inphl != nnphl!' 
!        IERR = 1 
!        RETURN
!     ENDIF
!
!.... resize
!     ALLOCATE(IWORK(NZ_FDIST))
!     IWORK(1:NZ_FDIST) = ICSC_FDIST(1:NZ_FDIST)
!     DEALLOCATE(ICSC_FDIST)
!     ALLOCATE(ICSC_FDIST(NZ_FDIST))
!     ICSC_FDIST(1:NZ_FDIST) = IWORK(1:NZ_FDIST)
!     DEALLOCATE(IWORK)

!  80 CONTINUE
!     CALL MPI_ALLREDUCE(ICOMBIG,MBUFRHS,1,MPI_INTEGER,MPI_MAX, &
!                        MYCOMM,MPIERR)
!
!.... this is a template for communicating
!     IF (LDEBUG) THEN
!        INPGL = 0
!        INPG0 = 0
!        ALLOCATE(LINIT(NA35)) 
!        LINIT(1:NA35) = .FALSE.
!        CALL MPI_COMM_SIZE(MYCOMM, NPROCS, MPIERR)
!        IF (MYID == MASTER) THEN
!           ALLOCATE(IROW(MBUFRHS))
!           DO 301 INPG=1,NNPG
!              INPINV = MASKG(INPG) !mask h better be in order!
!              IF (INPINV > 0) THEN
!                 IF (INPINV < INPG0) THEN
!                    WRITE(*,*) 'ghesstrct: maskh out or order!' 
!                    IERR = 1 
!                 ENDIF
!                 INPG0 = INPINV
!                 MYSRC = IGPART(INPINV) - 1 
!                 DO 302 IVINV=1,NVINV
!                    IF (MYSRC == MYID) THEN !copy my own column 
!                       INPGL = INPGL + 1 
!                       IROW1 = JCSC_FDIST(INPGL)
!                       IROW2 = JCSC_FDIST(INPGL+1) - 1 
!                       NROW = IROW2 - IROW1 + 1 
!                       IROW(1:NROW) = ICSC_FDIST(IROW1:IROW2)
!                    ELSE !receive a slave's row
!                       CALL MPI_RECV(NROW,           1,MPI_INTEGER, MYSRC,  &
!                                     MPI_ANY_TAG, MYCOMM,STAT,MPIERR) 
!                       CALL MPI_RECV(IROW(1:NROW),NROW,MPI_INTEGER, MYSRC,  &
!                                     MPI_ANY_TAG, MYCOMM,STAT,MPIERR)
!                    ENDIF !end check on source
!                    ILOC = (INPINV - 1)*NVINV + IVINV
!                    IF (LINIT(ILOC)) THEN
!                       IERR = 1 
!                       WRITE(*,*) 'ghesstrct: Wrong double init!',ILOC
!                    ENDIF
!                    LINIT(ILOC) = .TRUE.
!                    DO 303 IR=1,NROW
!                       WRITE(45+NPROCS,*) IROW(IR)
! 303                CONTINUE 
!                    !LINIT(ILOC) = .true.
! 302             CONTINUE !loop on inversion points
!              ENDIF !end check if node is in inversion 
! 301       CONTINUE !loop on anchor nodes 
!           DEALLOCATE(IROW) 
!           DO 305 ILOC=1,NA35
!              IF (.NOT.LINIT(ILOC)) THEN
!                 WRITE(*,*) 'ghesstrct: Did not initialize inversion variable:',ILOC
!                 IERR = 1
!              ENDIF
! 305       CONTINUE
!           DEALLOCATE(LINIT) 
!           CLOSE(45+NPROCS)
!        ELSE !slave process sends
!           MYDEST = MASTER
!           DO 401 INPG=1,NNPG
!              INPINV = MASKG(INPG)
!              IF (INPINV > 0) THEN
!                 IF (IGPART(INPINV) == MYID+1) THEN !send
!                    DO 402 IVINV=1,NVINV
!                       INPGL = INPGL + 1 
!                       IROW1 = JCSC_FDIST(INPGL)
!                       IROW2 = JCSC_FDIST(INPGL+1) - 1 
!                       NROW = IROW2 - IROW1 + 1 
!                       CALL MPI_SEND(NROW,                     1,MPI_INTEGER,   &
!                                     MYDEST,MYTAG,MYCOMM,MPIERR)
!                       CALL MPI_SEND(ICSC_FDIST(IROW1:IROW2),NROW,MPI_INTEGER,  &
!                                     MYDEST,MYTAG,MYCOMM,MPIERR)
! 402                CONTINUE !loop on variables in inversion
!                 ENDIF
!              ENDIF !end check if INPINV belongs to process 
! 401       CONTINUE !loop on anchor nodes
!        ENDIF !end check on ID
!     ENDIF !end checkon debugging
!     RETURN
!     END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      INTEGER*4 FUNCTION ICNPGL(MYID,NNPG,NNPINV,NVINV,IGPART,MASKG)
!
!     Counts the number of nodal points in Hessian partition
!
!     INPUT        MEANING
!     -----        ------- 
!     IGPART       holds hessian node partition
!     MASKG        gradient nodal point mask -> inpinv
!     MYID         process ID
!     NNPINV       number of nodal points in inversion
!     NNPG         number of anchor nodes
!     NVINV        number of variables in inversion
!     
!     OUTPUT       MEANING
!     ------       ------- 
!     ICNHPL       holds number of points in Hessian partition
!
!.... variable declarations
      INTEGER*4, INTENT(IN) :: MASKG(NNPG), IGPART(NNPINV), MYID, NNPG, NNPINV, NVINV 
!
!----------------------------------------------------------------------------------------!
!
      ICNPGL = 0 
      DO 100 INPG=1,NNPG
         INPINV = MASKG(INPG)
         IF (INPINV > 0) THEN !in inversion
            IF (IGPART(INPINV) == MYID + 1) ICNPGL = ICNPGL + 1 
         ENDIF !end check on node in inversion
  100 CONTINUE !loop on anchor nodes
      RETURN
      END 
