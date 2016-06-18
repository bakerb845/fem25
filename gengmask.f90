      sUBROUTINE GENGMASK(PROJNM, MDIM,MEN,MGNOD, NDOF,NNPG,NELEM,                &
                          NLXI,NLETA,LINREC,LKBDRY,                               &
                          CDOMAIN,CNNPG,PART,                                     &
                          LM,IENG, RCV,INV, IERR) 
!
!     Sets mask pointer for the inverse problem.  At each anchor node determines 
!     if point is active, value > 0, and then if it is, what its number in the 
!     inverse problem is.  Note that  maxval(maskg) = nnpinv.  Additionally, we 
!     read the weighting mask for the search direction.  wmask has values between 0 
!     and 1.  wmask would be helpful if we really want to focus on a region in the
!     inversion  or apply a cheapy gain control - B. Baker Dec. 2012
!
!     INPUT      MEANING
!     -----      -------
!     CDOMAIN    element domains 'I' interior, 'E' exterior, or 'A' absorbing
!     CNNPG      character description of anchor node.  only invert at 
!                interior nodes 'II' and free surface 'FS'
!     IENG       element -> anchor node number
!     LINREC     True -> let receivers be included in inversion
!     LKBDRY     True -> keep elements on bielak/interior layer
!     LM         element -> global DOF
!     MDIM       leading dimension for LM
!     MEN        max number of element nodes
!     MGNOD      max number of anchor nodes
!     NELEM      number of elements
!     NLETA      number of lagrange interpolation points in eta
!     NLXI       number of lagrange interpolation points in xi
!     NNPG       number of anchor nodes 
!     PROJNM     project name
!    
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     MASKG      inversion point -> anchor node pointer (output through common)
!     NNPINV     number of nodal points in inversion (output through common)
!     WMASK      weights to apply to search direction (output through common)
!
!.... variable declarations  
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV
      TYPE (RECV_INFO) RCV
      CHARACTER(*), INTENT(IN) :: PROJNM 
      CHARACTER(2), INTENT(IN) :: CNNPG(NNPG) 
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*),IENG(MGNOD,*),PART(NDOF), MDIM,MEN,MGNOD, &
                               NDOF,NNPG, NELEM, NLXI,NLETA
      LOGICAL*4, INTENT(IN) :: LINREC, LKBDRY 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) FILENM
      REAL*4, ALLOCATABLE :: TEMP(:) 
      LOGICAL*4 LMASK, LKEEP, LSHIFT, LBDR
      INTEGER*4 NNPG_IN, N0, NNPSV, NPROCS, IN0, INPINV, INPG, IREC, I, IA, IAE,  &
                IDX_PART, IELEM, IFAIL, INDX, IVINV, IUNIT
      PARAMETER(IUNIT = 44) 
      INTEGER*4 ICMCON, IA2IAE  
!
!----------------------------------------------------------------------------------------!
!
!.... first set the anchor node element connectivit, count connectivity of elements
      inv%NCON = ICMCON(NGNOD,NNPG, NELEM,NGNOD,IENG) 
!
!.... and connect them
      ALLOCATE(inv%MCONN(NNPG,inv%NCON))
      CALL CONNECT(NGNOD,NNPG, NNPG,NELEM,NGNOD,inv%NCON, IENG,inv%MCONN)

!.... check if a mask file exists
      IERR = 0
      FILENM(1:80) = ' ' 
      FILENM = TRIM(ADJUSTL(PROJNM))//'.mask'
      FILENM = ADJUSTL(FILENM)
      INQUIRE(FILE=TRIM(FILENM),EXIST=LMASK) 
      ALLOCATE(TEMP(NNPG))  
      IF (LMASK) THEN 
         IFAIL = 2
         WRITE(*,*) 'invgraph: Reading mask file...'
         OPEN(UNIT=IUNIT,FILE=TRIM(FILENM))
         READ(IUNIT,*,END=50) NNPG_IN 
         IF (NNPG_IN /= NNPG) THEN
            WRITE(*,*) 'invgraph: nnpg_in != nnpg'
            IFAIL = 1
            GOTO 50
         ENDIF
         DO 1 INPG=1,NNPG
            READ(IUNIT,*,END=50) TEMP(INPG)  
            IF (TEMP(INPG) > 1.0) THEN
               WRITE(*,*) 'invgraph: Warning wmask > 1, resetting to 1'
               TEMP(INPG) = 1.0
            ELSEIF (TEMP(INPG) < 0.0) THEN
               WRITE(*,*) 'invgraph: Warning wmask < 0., resseint to 0'
               TEMP(INPG) = 0.0
            ENDIF
    1    CONTINUE
         IFAIL = 0  
   50    CONTINUE 
         CLOSE(IUNIT) 
         IF (IFAIL == 2) & 
         WRITE(*,*) 'invgraph: Premature end of mask file.  Ignoring'
         IF (IFAIL > 0) TEMP(1:NNPG) = 1.0 !set mask weights to unity
      ELSE
         TEMP(1:NNPG) = 1.0
      ENDIF
!
!.... count the inversion points in the mesh 
      INPINV = 0 
      ALLOCATE(inv%MASKG(NNPG)) 
      inv%MASKG(1:NNPG) = 0 
      DO 10 INPG=1,NNPG
!
!....... can get FS points in bielak and absorbing elements
         LKEEP  = .TRUE.
         IELEM = inv%MCONN(INPG,1) 
         IF (IELEM > 0) THEN
            IF (CDOMAIN(IELEM) /= 'I') LKEEP = .FALSE. 
         ENDIF 
!
!....... check point is interior
         IF (LKEEP .AND. (CNNPG(INPG) == 'II' .OR. &
                          CNNPG(INPG) == 'BI' .OR. &
                          CNNPG(INPG) == 'FS')) THEN
            IF (TEMP(INPG) > 0.0) THEN !we aren't externally masking it
               INPINV = INPINV + 1
               inv%MASKG(INPG) = INPINV
            ELSE
               inv%MASKG(INPG) = 0
            ENDIF
         ELSE 
            inv%MASKG(INPG) = 0
         ENDIF
   10 CONTINUE 
      LSHIFT = .FALSE.
      inv%NNPINV = INPINV 
      NNPSV  = inv%NNPINV
!
!.... now we need to make sure inversion pionts are attached to DOFs
      DO 11 IELEM=1,NELEM
         DO 12 IA=1,NGNOD
            INPG = IENG(IA,IELEM)
            IF (INPG == 0) THEN
               WRITE(*,*) 'invgraph: Triangles not yet programmed!'
               IERR = 1
               RETURN
            ENDIF
            IF (inv%MASKG(INPG) > 0) THEN !nodal point is active
               IAE = IA2IAE(NLXI,NLETA,IA)
               IF (IAE > 0) THEN
                  DO 13 I=1,NDIM
                     IF (LM(I,IAE,IELEM) <= 0) THEN
                        inv%NNPINV = inv%NNPINV - 1 
                        inv%MASKG(INPG) = 0 
                        LSHIFT = .TRUE.
                        GOTO 130
                     ENDIF
   13             CONTINUE
!
!................ it is possible we may also want to discount receivers
                  IF (.NOT.LINREC .AND. IAE > 0) THEN
                     DO 14 IREC=1,rcv%NREC
                        IF (LM(1,IAE,IELEM) == rcv%MRDOF(1,IREC)) THEN 
                           inv%NNPINV = inv%NNPINV - 1
                           inv%MASKG(INPG) = 0
                           LSHIFT = .TRUE.
                           GOTO 130
                        ENDIF
   14                CONTINUE 
                  ENDIF
  130             CONTINUE 
               ELSE
                  WRITE(*,*) 'invgraph: Error in ia2iae'
                  IERR = 1
                  RETURN
               ENDIF
            ENDIF
   12    CONTINUE !loop on anchor nodes
   11 CONTINUE !loop on elements
!
!.... kick out elements bordering bielak layer
      IF (.NOT.LKBDRY) THEN
         DO 15 IELEM=1,NELEM
            LBDR = .FALSE.
            DO 16 IA=1,NGNOD  
               INPG = IENG(IA,IELEM)   
               IF (INPG > 0) THEN
                  IF (CNNPG(INPG) /= 'II' .AND. CNNPG(INPG) /= 'FS') LBDR = .TRUE.
               ENDIF 
   16       CONTINUE
            IF (LBDR) THEN
               DO 17 IA=1,NGNOD
                  INPG = IENG(IA,IELEM)
                  IF (inv%MASKG(INPG) > 0) THEN
                     inv%MASKG(INPG) = 0 
                     inv%NNPINV = inv%NNPINV - 1
                     LSHIFT = .TRUE.
                  ENDIF
   17          CONTINUE
            ENDIF
   15    CONTINUE 
      ENDIF
  
!
!.... double check we have point to invert at
      IF (inv%NNPINV < 1) THEN
         WRITE(*,*) 'invgraph: Error no points to invert at!'
         IERR = 1
      ELSEIF (inv%NNPINV > NNPG) THEN
         WRITE(*,*) 'invgraph: Error number of inversion points > than anchor nodes!'
         IERR = 1
      ELSE
         WRITE(*,*) 'invgraph: Number of inversion points:',inv%NNPINV
      ENDIF
!
!.... need to pull gradient inversion points backwards
      IF (LSHIFT) THEN
         INPINV = 0
         DO 18 INPG=1,NNPG
            IF (inv%MASKG(INPG) > 0) THEN
               INPINV = INPINV + 1
               inv%MASKG(INPG) = INPINV
            ENDIF
   18    CONTINUE
      ENDIF 
!
!.... reorder to encourage diagonal dominance
!     CALL REORDER_GRAD(MGNOD,NNPG,NNPG, NGNOD,inv%NCON, IENG,inv%MCONN, inv%MASKG,IERR)
!     call  REORDER_GRAD(MGNOD,NNPG,NNPG, NGNOD,inv%NCON, IENG,inv%MCONN, inv%MASKG,IERR)
!
!.... load balance partition by approximately equally distributing
      NPROCS = MAXVAL(PART) 
      N0 = MINVAL(PART)
      IF (N0 /= 1) WRITE(*,*) 'invgraph: WARNING! process 0 is not ID 0'
      ALLOCATE(inv%IGPART(inv%NNPINV)) 
      IDX_PART = inv%NNPINV/NPROCS 
      IN0 = 1
      DO 21 INPINV=1,inv%NNPINV
         IF (INPINV > IN0*IDX_PART) THEN
            IN0 = IN0 + 1
            N0 = N0 + 1
            IF (N0 > NPROCS) N0 = NPROCS 
         ENDIF
         inv%IGPART(INPINV) = N0
   21 CONTINUE 
      IF (N0 > NPROCS) THEN
         WRITE(*,*) 'invgraph: N0 is too large!'
         IERR = 1
      ENDIF
!
!.... resize mask 
      ALLOCATE(inv%WMASK(inv%NNPINV*inv%NVINV))  
      INPINV = 0
      DO 22 INPG=1,NNPG
         INPINV = inv%MASKG(INPG) 
         IF (INPINV > 0) THEN
            DO 23 IVINV=1,inv%NVINV
               INDX = (IVINV - 1)*inv%NNPINV + INPINV 
               inv%WMASK(INDX)  =  TEMP(INPG) 
   23       CONTINUE
         ENDIF
   22 CONTINUE
      DEALLOCATE(TEMP) 
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE REORDER_GRAD(MGNOD,MNPG,NNPG, NGNOD,NCON, IENG,MCONN, MASKG,IERR)
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MCONN(MNPG,*), MGNOD,MNPG,NNPG, &
                               NGNOD,NCON  
      INTEGER*4, INTENT(INOUT) :: MASKG(NNPG)  
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4, ALLOCATABLE :: IRPTR(:), JCPTR(:), XADJ(:), ADJNCY(:), PERM(:), &
                                IPERM(:), JWORK(:) 
      INTEGER*4 IOPT8(8) 
!
!----------------------------------------------------------------------------------------!
!
!.... null out distance
      IERR = 0
      WRITE(*,*) 'reorder_grad: Generating a list...'
      NNPINV = MAXVAL(MASKG) 
      N     = NNPINV              !numbers of rows/columsn in graph
      NWORK = NNPINV*NCON*NGNOD   !over-estimate number of non-zeros in graph
      ALLOCATE(IRPTR(N+1))
      ALLOCATE(JCPTR(NWORK)) 
      NWORK2 = NCON*NGNOD + 1
      ALLOCATE(JWORK(NWORK2)) 
!
!.... loop on anchor nodes
      IZERO = 0
      IRPTR(1) = 1
      I = 0 
      DO 1 INPG=1,NNPG
         INPINV = MASKG(INPG)
         IF (INPINV > 0) THEN
            I = I + 1
            IF (INPINV /= I) THEN
               WRITE(*,*) 'reorder_grad: Warning need to think harder'
               IERR = 1
            ENDIF
            IWORK = 0
            DO 2 ICON=1,NCON
               IELEM = MCONN(INPG,ICON)
               IF (IELEM == 0) GOTO 20
               DO 3 IA=1,NGNOD
                  JNPG = IENG(IA,IELEM)
                  IF (JNPG == 0) GOTO 30
                  JNPINV = MASKG(JNPG) 
                  IF (JNPINV > 0) THEN !is in inversion mesh
                     IWORK = IWORK + 1
                     JWORK(IWORK) = JNPINV
                  ENDIF !end check on inversion mesh
    3          CONTINUE !loop on anchor nodes
   30          CONTINUE !break ahead, out of anchor nodes
    2       CONTINUE !loop on connections
   20       CONTINUE !break ahead, out of connections
!
!.......... now sort the column 
            NSORT = IWORK
            CALL ISHELL1(NSORT,JWORK(1:NSORT)) 
            JWORK(NSORT+1) =-1 
!
!.......... and add the non-repeats to list
            DO 4 ISORT=1,NSORT
               IF (JWORK(ISORT) /= JWORK(ISORT+1)) THEN
                  IZERO = IZERO + 1 
                  JCPTR(IZERO) = JWORK(ISORT) 
               ENDIF  
    4       CONTINUE  
            IRPTR(I+1) = IZERO + 1  
         ENDIF !end check on in inversion mesh
    1 CONTINUE !loop on anchor nodes
      NZERO = IZERO  
      do i=1,n
         j1 = irptr(i)
         j2 = irptr(i+1) - 1
         do j=j1,j2
            write(44,*) n + 1 - i,jcptr(j) 
         enddo
      enddo
      close(44)

      NADJ = NZERO - N 
      ALLOCATE(XADJ(N+1))
      ALLOCATE(ADJNCY(NADJ))
      CALL CRS2ADJ(NZERO,N, IRPTR,JCPTR, XADJ,ADJNCY)
!
!.... minimize the bandwidth  
      WRITE(*,*) 'reorder_grad: Generating RCM ordering...'
      ALLOCATE(PERM(N)) 
      ALLOCATE(IPERM(N)) 
      IOPT8(1) = 0
      NUMFLAG = 1
      CALL METIS_NodeND(N,XADJ,ADJNCY,NUMFLAG,IOPT8,IPERM,PERM)
!     CALL GENRCM(N,NADJ,XADJ,ADJNCY,PERM)
!
!.... reorder
      DO 12 INPG=1,NNPG
         INPINV = MASKG(INPG)
         IF (INPINV > 0) then
            print *, inpinv,perm(inpinv)
            MASKG(INPINV) = PERM(INPINV) 
         endif
   12 CONTINUE 
      do i=1,n
         ip = i !perm(i)
         j1 = irptr(ip)
         j2 = irptr(ip+1) - 1
         do j=j1,j2
            write(45,*) n + 1 - ip,iperm(jcptr(j)) 
         enddo
      enddo
      DEALLOCATE(IRPTR)
      DEALLOCATE(JCPTR) 
      DEALLOCATE(XADJ)
      DEALLOCATE(ADJNCY) 
      DEALLOCATE(PERM) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE WTGRAD(NA35,WGRAD, GRAD) 
!
!     Applies a scaling of (0,1] to inversion points in the gradient. 
! 
!     INPUT      MEANING
!     -----      ------- 
!     GRAD       gradient to multiply by coefficent (0,1]
!     NA35       size of gradient
!     WGRAD      gradient weight (0,1]
!
!     OUTPUT     MEANING
!     ------     ------- 
!     GRAD       rescaled gradient
!
!.... variable declarations 
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: WGRAD(NA35) 
      INTEGER*4, INTENT(IN) :: NA35 
      REAL*4, INTENT(INOUT) :: GRAD(NA35)
      INTEGER*4 IA35
! 
!----------------------------------------------------------------------------------------!
! 
      DO 1 IA35=1,NA35  
         GRAD(IA35) = GRAD(IA35)*WGRAD(IA35)
    1 CONTINUE
!     DO 1 INPINV=1,NNPINV
!        DO 2 IVINV=1,NVINV 
!           INDX = (IVINV - 1)*NNPINV + INPINV  
!           GRAD(INDX) = GRAD(INDX)*WGRAD(INPINV) 
!   2    CONTINUE
!   1 CONTINUE 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE WTSRCH8(NA35,WGRAD, SRCH)
!
!     Applies a scaling of (0,1] to search direction. 
! 
!     INPUT      MEANING
!     -----      ------- 
!     GRAD       gradient to multiply by coefficent (0,1]
!     NA35       size of gradient
!     WGRAD      gradient weight (0,1]
!
!     OUTPUT     MEANING
!     ------     ------- 
!     SRCH       rescaled search direction
!
!.... variable declarations
      IMPLICIT NONE 
      REAL*4, INTENT(IN) :: WGRAD(NA35)
      INTEGER*4, INTENT(IN) :: NA35
      REAL*8, INTENT(INOUT) :: SRCH(NA35)
      INTEGER*4 IA35
      DO 1 IA35=1,NA35
         SRCH(IA35) = SRCH(IA35)*WGRAD(IA35)
    1 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE CWTGRAD(NA35,WGRAD, GRAD) 
!
!     Applies a scaling of (0,1] to inversion points in the gradient. 
! 
!     INPUT      MEANING
!     -----      ------- 
!     GRAD       gradient to multiply by coefficent (0,1]
!     NA35       size of gradient
!     WGRAD      gradient weight (0,1]
!
!     OUTPUT     MEANING
!     ------     ------- 
!     GRAD       rescaled gradient
!
!.... variable declarations 
      IMPLICIT NONE
      REAL*4, INTENT(IN) :: WGRAD(NA35) 
      INTEGER*4, INTENT(IN) :: NA35
      COMPLEX*8, INTENT(INOUT) :: GRAD(NA35)
      INTEGER*4 IA35
! 
!----------------------------------------------------------------------------------------!
! 
      DO 1 IA35=1,NA35
         GRAD(IA35) = GRAD(IA35)*CMPLX(WGRAD(IA35),0.0)
    1 CONTINUE
!     DO 1 INPINV=1,NNPINV
!        DO 2 IVINV=1,NVINV 
!           INDX = (IVINV - 1)*NNPINV + INPINV  
!           GRAD(INDX) = GRAD(INDX)*CMPLX(WGRAD(INPINV),0.0)
!   2    CONTINUE
!   1 CONTINUE 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      REAL*8 FUNCTION GET_STP0(NA35,DP0, XMOD,P)
!
!     Returns steplength on desired max perturbation
!
!     INPUT      MEANING
!     DP0        max perturbation (percentage)
!     P          search direction
!     XMOD       model
! 
!     OUTPUT     MEANING
!     ------     -------
!     GET_STP0   step length corresponding to max perturbation
!
!.... local variables
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XMOD(NA35), P(NA35), DP0
      INTEGER*4, INTENT(IN) :: NA35
!.... local variables
      REAL*8 DPMAX, PCT
      INTEGER*4 IMAX, IA35
!
!----------------------------------------------------------------------------------------!
!
      IF (DP0 <= 0.D0) THEN
         GET_STP0 = 1.D0 
      ELSE
         IMAX = 0
         DPMAX = 0.D0
         DO 1 IA35=1,NA35
            PCT = DABS(P(IA35)/XMOD(IA35)) !( x1 - x0 )/x0 = (x0 + dp - x0)/x0 = p/x0 
            IF (PCT > DPMAX) THEN
               DPMAX = PCT
               IMAX = IA35
            ENDIF
    1    CONTINUE
         IF (DPMAX == 0.D0) THEN
            WRITE(*,*) 'get_stp0: Error dpmax = 0.0; returning 1'
            GET_STP0 = 1.D0
            RETURN
         ENDIF
         IF (P(IMAX) == 0.D0) THEN
            WRITE(*,*) 'get_stp0: pmax == 0.0; returning 1'
            GET_STP0 = 1.D0
         ENDIF
         !dp = (x0 - (x0 + alpha p_{max}))/x0 = alpha p_{max}/x0
         !alpha = x0 dp/p_{max}
         GET_STP0 = DP0/100.D0/DPMAX !1./dpmax = X0/P(IMAX) 
     ENDIF
     RETURN
     END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE UPDMOD8(PROJNM,LFILES,NA35, IBLOCK,K,IALPHA,    &
                         VPMIN_INV,VPMAX_INV, XMOD, INV,MSH,IERR) 
!
!     Updates the model and acts as an interface
      IMPLICIT NONE 
      INCLUDE 'fwd_struc.h'
      CHARACTER(*), INTENT(IN) :: PROJNM 
      REAL*8, INTENT(IN) :: XMOD(NA35), VPMIN_INV, VPMAX_INV 
      TYPE (MESH_INFO) MSH
      TYPE (INV_INFO)  INV  
      INTEGER*4, INTENT(IN) :: NA35, IBLOCK,K,IALPHA 
      LOGICAL*4, INTENT(IN) :: LFILES 
      INTEGER*4, INTENT(OUT) :: IERR 
      REAL*8 ALPHA, BETA, SQRT3, POISSON, AMIN, BMIN, AMAX, BMAX, PMIN, PMAX, DMFRHO, &
             VPVS2
      INTEGER*4 INPG, INPINV, ILOC1, ILOC2, JOB  
      PARAMETER(SQRT3 = 1.7320508075688772D0)
!
!----------------------------------------------------------------------------------------!
!
!.... classify inverse problem
      IERR = 0  
      IF (inv%CINVTYPE == 'PP' .OR. inv%CINVTYPE == 'pp') THEN 
         JOB = 1  
      ELSEIF (inv%CINVTYPE == 'SS' .OR. inv%CINVTYPE == 'ss') THEN 
         JOB = 2  
      ELSEIF (inv%CINVTYPE == 'PS' .OR. inv%CINVTYPE == 'ps') THEN 
         JOB = 3  
      ELSE 
         WRITE(*,*) 'upmod: Cannot do anisotropy yet!'
         IERR = 1  
         RETURN
      ENDIF
      IERR = 0
      VPVS2 = inv%VPVS**2 
      AMIN = HUGE(1.D0)
      BMIN = HUGE(1.D0)
      AMAX = 0.D0
      BMAX = 0.D0
      DO 1 INPG=1,msh%NNPG
         INPINV = inv%MASKG(INPG)
         IF (INPINV > 0) THEN
            ILOC1 = (INPINV - 1)*inv%NVINV + 1
            ILOC2 = (INPINV - 1)*inv%NVINV + 2
            IF (JOB == 1) THEN
               ALPHA = XMOD(ILOC1)
               BETA  = ALPHA/inv%VPVS !VPVS2 !SQRT3
            ELSEIF (JOB == 2) THEN
               BETA  = XMOD(ILOC1)
               ALPHA = BETA*inv%VPVS !VPVS2 !SQRT3
            ELSEIF (JOB == 3) THEN
               ALPHA = XMOD(ILOC1)
               BETA  = XMOD(ILOC2)
            ELSE
               WRITE(*,*) 'updmod8: Cannot do anisotropy yet!'
            ENDIF
!
!.......... warning
            IF (ALPHA < 0.D0) THEN
               WRITE(*,*) 'upmod8: Error Vp < 0.',ALPHA
               IERR = 1
               RETURN
            ENDIF
            IF (BETA < 0.D0) THEN
               WRITE(*,*) 'upmod8: Error Vs < 0.',BETA
               IERR = 1
               RETURN
            ENDIF
!
!.......... clipping? 
            IF (VPMIN_INV >= 0.D0) THEN
               IF (ALPHA < VPMIN_INV) THEN
                  ALPHA = VPMIN_INV
                  BETA  = ALPHA/inv%VPVS !VPVS2 !DSQRT(3.D0)
               ENDIF
            ENDIF
            IF (VPMAX_INV > 0.D0) THEN
               IF (ALPHA > VPMAX_INV) THEN
                  ALPHA = VPMAX_INV
                  BETA  = ALPHA/inv%VPVS !VPVS2 !DSQRT(3.D0)
               ENDIF
            ENDIF
!
!.......... check poisson's ratio for fluid like layers
            POISSON = (ALPHA**2 - 2.D0*BETA**2)/(2.D0*(ALPHA**2 - BETA**2))
            IF (POISSON < 0.22D0) &
            WRITE(*,*) 'updmod8: Warning on Poissons ratio:',POISSON
            IF (POISSON > 0.35D0) &
            WRITE(*,* )'updmod8: Warning on Poissons ratio:',POISSON

            msh%DENS(INPG) = DMFRHO(ALPHA) !update density 
            msh%ECOEFF(INPG,2) = msh%DENS(INPG)*BETA**2
            msh%ECOEFF(INPG,1) = msh%DENS(INPG)*ALPHA**2 - 2.D0*msh%ECOEFF(INPG,2)

            AMIN = DMIN1(AMIN,ALPHA)
            BMIN = DMIN1(BMIN,BETA)
            AMAX = DMAX1(AMAX,ALPHA)
            BMAX = DMAX1(BMAX,BETA)
         ENDIF
    1 CONTINUE
      PMIN = (AMIN**2 - 2.D0*BMIN**2)/(2.D0*(AMIN**2 - BMIN**2))
      PMAX = (AMAX**2 - 2.D0*BMAX**2)/(2.D0*(AMAX**2 - BMAX**2))
      WRITE(*,900) AMIN,BMIN, AMAX,BMAX, PMIN,PMAX
  900 FORMAT(' updmod8: Minimum P and S wave velocity: ',F12.3,F12.3,' m/s',/, &
             '          Maximum P and S wave velocity: ',F12.3,F12.3,' m/s',/, &
             '          Minimum and max Poissons ratio:',F12.3,F12.3,/)
!
!.... write updated models
      CALL PLOT_ELMOD_UPD(PROJNM,NGNOD,msh%NNPG,msh%NNPG, msh%NELEM, msh%LISISO, &
                          IBLOCK,K,IALPHA, msh%CDOMAIN, &
                          msh%IENG,msh%XLOCS,msh%ZLOCS, msh%DENS,msh%ECOEFF)
      RETURN
      END
!             

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE UPDMOD(PROJNM,LFILES, IBLOCK,K,IALPHA, XMOD, INV,MSH,IERR) 
!
!     Updates the model and acts as an interface
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      CHARACTER(80), INTENT(IN) :: PROJNM 
      REAL*4, INTENT(IN) :: XMOD(*) 
      TYPE (MESH_INFO) MSH
      TYPE (INV_INFO)  INV
      INTEGER*4, INTENT(IN) :: IBLOCK,K,IALPHA 
      LOGICAL*4, INTENT(IN) :: LFILES 
      INTEGER*4, INTENT(OUT) :: IERR
      REAL*8 ALPHA, BETA, SQRT3, POISSON, AMIN, BMIN, AMAX, BMAX, PMIN, PMAX, DMFRHO, &
             VPVS2
      INTEGER*4 INPG, INPINV, ILOC1, ILOC2, JOB  
      PARAMETER(SQRT3 = 1.7320508075688772D0)
!
!----------------------------------------------------------------------------------------!
!
!.... classify inverse problem
      IERR = 0 
      IF (inv%CINVTYPE == 'PP' .OR. inv%CINVTYPE == 'pp') THEN
         JOB = 1 
      ELSEIF (inv%CINVTYPE == 'SS' .OR. inv%CINVTYPE == 'ss') THEN
         JOB = 2 
      ELSEIF (inv%CINVTYPE == 'PS' .OR. inv%CINVTYPE == 'ps') THEN
         JOB = 3 
      ELSE
         WRITE(*,*) 'upmod: Cannot do anisotropy yet!'
         IERR = 1 
         RETURN
      ENDIF
      VPVS2 = inv%VPVS**2

      IERR = 0
      AMIN = HUGE(1.D0)
      BMIN = HUGE(1.D0)
      AMAX = 0.D0
      BMAX = 0.D0
      DO 1 INPG=1,msh%NNPG
         INPINV = inv%MASKG(INPG)
         IF (INPINV > 0) THEN
            ILOC1 = (INPINV - 1)*inv%NVINV + 1
            ILOC2 = (INPINV - 1)*inv%NVINV + 2
            IF (JOB == 1) THEN
               ALPHA = DBLE(XMOD(ILOC1))
               BETA  = ALPHA/inv%VPVS !VPVS2 !SQRT3
            ELSEIF (JOB == 2) THEN
               BETA  = DBLE(XMOD(ILOC1))
               ALPHA = BETA*inv%VPVS !VPVS2 !SQRT3 
            ELSEIF (JOB == 3) THEN 
               ALPHA = DBLE(XMOD(ILOC1))
               BETA  = DBLE(XMOD(ILOC2))
            ELSE 
               WRITE(*,*) 'updmod: Cannot do anisotropy yet!' 
            ENDIF
!
!.......... warning
            IF (ALPHA < 0.D0) THEN
               WRITE(*,*) 'upmod: Error Vp < 0.',ALPHA
               IERR = 1
               RETURN
            ENDIF
            IF (BETA < 0.D0) THEN
               WRITE(*,*) 'upmod: Error Vs < 0.',BETA
               IERR = 1
               RETURN
            ENDIF
            POISSON = (ALPHA**2 - 2.D0*BETA**2)/(2.D0*(ALPHA**2 - BETA**2)) 
            IF (POISSON < 0.22D0) &
            WRITE(*,*) 'updmod: Warning on Poissons ratio:',POISSON
            IF (POISSON > 0.35D0) &
            WRITE(*,* )'updmod: Warning on Poissons ratio:',POISSON
 
            msh%DENS(INPG) = DMFRHO(ALPHA) !update density 
            msh%ECOEFF(INPG,2) = msh%DENS(INPG)*BETA**2
            msh%ECOEFF(INPG,1) = msh%DENS(INPG)*ALPHA**2 - 2.D0*msh%ECOEFF(INPG,2)

            AMIN = DMIN1(AMIN,ALPHA)
            BMIN = DMIN1(BMIN,BETA)
            AMAX = DMAX1(AMAX,ALPHA)
            BMAX = DMAX1(BMAX,BETA) 
         ENDIF
    1 CONTINUE 

      PMIN = (AMIN**2 - 2.D0*BMIN**2)/(2.D0*(AMIN**2 - BMIN**2))
      PMAX = (AMAX**2 - 2.D0*BMAX**2)/(2.D0*(AMAX**2 - BMAX**2)) 
      WRITE(*,900) AMIN,BMIN, AMAX,BMAX, PMIN,PMAX 
  900 FORMAT(' updmod: Minimum P and S wave velocity: ',F12.3,F12.3,' m/s',/, &
             '         Maximum P and S wave velocity: ',F12.3,F12.3,' m/s',/, &
             '         Minimum and max Poissons ratio:',F12.3,F12.3,/)
!
!.... write updated models
      CALL PLOT_ELMOD_UPD(PROJNM,NGNOD,msh%NNPG,msh%NNPG, msh%NELEM, msh%LISISO,  &
                          IBLOCK,K,IALPHA, msh%CDOMAIN, &
                          msh%IENG,msh%XLOCS,msh%ZLOCS, msh%DENS,msh%ECOEFF)  
      RETURN 
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE ELEM_WTS(MSH,INV,IERR) 
!
!     For the virtual forces we calculate the virtual force matrix 
!     [dS/dm_1 u, ..., dS/dm_m u] however, where dS/dm_i are the weights for the 
!     i'th virtual force.  Since, some elements are larger than others each the 
!     numerical integration over an element artificially pays attention to larger 
!     elements and ignores smaller elements.  Consequently, we create a series of 
!     weights to renormalize the dS/dm u calculation to w dS/dm u. 

      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH
      TYPE (INV_INFO) INV
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      LOGICAL*4, ALLOCATABLE :: LINIT(:) 
      REAL*8 WTMAX, X1,X2,Z1,Z2,RLOC  
      LOGICAL*4 LQUAD 
      INTEGER*4 INPG, ICON, IELEM, IMOD, JNPG1,JNPG2, IA1,IA2 
!
!----------------------------------------------------------------------------------------!
!
!.... initialize 
      IERR = 0
      ALLOCATE(inv%ELEM_WTS(msh%NELEM)) 
      ALLOCATE(LINIT(msh%NELEM)) 
      LINIT(1:msh%NELEM) = .FALSE. 
      WTMAX = 0.D0 
      inv%ELEM_WTS(1:msh%NELEM) = 0.D0
!
!.... loop on anchor nodes
      DO 1 INPG=1,msh%NNPG
         IF (inv%MASKG(INPG) > 0) THEN !in inversion
            DO 2 ICON=1,inv%NCON
               IELEM = inv%MCONN(INPG,ICON) 
               IF (IELEM <= 0) GOTO 20
               IF (.NOT.LINIT(IELEM)) THEN
                  LQUAD = .TRUE.
                  IMOD = 5
                  IF (MINVAL(msh%IENG(1:NGNOD,IELEM)) == 0) THEN
                     LQUAD = .FALSE.
                     IMOD = 4
                  ENDIF
!
!................ calculate the element area with Green's theorem 
                  RLOC = 0.D0 
                  DO 3 IA1=1,NGNOD  
                     JNPG1 = msh%IENG(IA1,IELEM)
                     IF (JNPG1 == 0) GOTO 30 !break for triangle 
                     IA2 = MOD(IA1+1,IMOD)
                     IF (IA2 == 0) IA2 = 1
                     JNPG2 = msh%IENG(IA2,IELEM)
                     X1 = msh%XLOCS(JNPG1)
                     X2 = msh%XLOCS(JNPG2)
                     Z1 = msh%ZLOCS(JNPG1)
                     Z2 = msh%ZLOCS(JNPG2)
                     RLOC = RLOC + X1*Z2 - X2*Z1
    3             CONTINUE  
   30             CONTINUE !break ahead for triangle
                  IF (RLOC <= 0.D0) THEN
                     WRITE(*,*) 'elem_wts: Warning your mesh is wrong'
                     RLOC =-RLOC
                  ENDIF 
                  RLOC = RLOC*0.5D0 !only going one way here
                  inv%ELEM_WTS(IELEM) = RLOC
                  LINIT(IELEM) = .TRUE.
                  WTMAX = MAX(inv%ELEM_WTS(IELEM),WTMAX)
               ENDIF !end check on intialized 
    2       CONTINUE !loop on connectivity
   20       CONTINUE !break ahead  
         ENDIF !end check on node in inversion 
    1 CONTINUE 
      DEALLOCATE(LINIT) 
      IF (WTMAX <= 0.D0) THEN
         WRITE(*,*) 'elem_wts: Error weighting is <= 0!',WTMAX
      ENDIF
!
!.... normalize the elements
      DO 4 IELEM=1,msh%NELEM
         IF (inv%ELEM_WTS(IELEM) > 0.D0) inv%ELEM_WTS(IELEM) = WTMAX/inv%ELEM_WTS(IELEM)
         !if (inv%elem_wts(ielem) > 0.d0) inv%elem_wts(ielem) = 1.d0
    4 CONTINUE 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GENXMOD8(MNPG,NA35,NNPG, CINVTYPE, MASKG,DENS,ECOEFF, XMOD,IERR) 
      
      CHARACTER(2), INTENT(IN) :: CINVTYPE
      REAL*8, INTENT(IN) :: ECOEFF(MNPG,*), DENS(NNPG)
      INTEGER*4, INTENT(IN) :: MASKG(NNPG), MNPG,NA35,NNPG
      REAL*8, INTENT(OUT) :: XMOD(NA35) 
      INTEGER*4, INTENT(OUT) :: IERR 
      REAL*8 ALPHA,BETA
!
!----------------------------------------------------------------------------------------!
!
!.... classify inverse problem
      IERR = 0
      NVINV = 0
      IF (CINVTYPE == 'PP' .OR. CINVTYPE == 'pp') THEN
         JOB = 1
         NVINV = 1
      ELSEIF (CINVTYPE == 'SS' .OR. CINVTYPE == 'ss') THEN
         JOB = 2
         NVINV = 1
      ELSEIF (CINVTYPE == 'PS' .OR. CINVTYPE == 'ps') THEN
         JOB = 3
         NVINV = 2
      ELSE
         WRITE(*,*) 'genxmod: Cannot do anisotropy yet!'
         IERR = 1
         RETURN
      ENDIF
      IMARK = 0
      DO 1 INPG=1,NNPG
         INPINV = MASKG(INPG)
         IF (INPINV > 0) THEN
            IMARK = IMARK + 1
            ALPHA = DSQRT( (ECOEFF(INPG,1) + 2.D0*ECOEFF(INPG,2))/DENS(INPG) )
            BETA  = DSQRT( ECOEFF(INPG,2)/DENS(INPG) )
            ILOC1 = (INPINV - 1)*NVINV + 1
            ILOC2 = (INPINV - 1)*NVINV + 2
            IF (JOB == 1) THEN
               XMOD(ILOC1) = ALPHA
            ELSEIF (JOB == 2) THEN
               XMOD(ILOC1) = BETA
            ELSEIF (JOB == 3) THEN
               XMOD(ILOC1) = ALPHA
               XMOD(ILOC2) = BETA
            ELSE
               WRITE(*,*) 'genxmod: No anisotropy yet'
               IERR = 1
            ENDIF
         ENDIF
    1 CONTINUE
      IF (IMARK*NVINV /= NA35) THEN
         WRITE(*,*) 'genxmod: Error inpvin /= na35'
         IERR = 2
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GENXMOD(MNPG,NA35,NNPG, CINVTYPE, MASKG,DENS,ECOEFF, XMOD,IERR) 
      CHARACTER(2), INTENT(IN) :: CINVTYPE
      REAL*8, INTENT(IN) :: ECOEFF(MNPG,*), DENS(NNPG)
      INTEGER*4, INTENT(IN) :: MASKG(NNPG), MNPG,NA35,NNPG
      REAL*4, INTENT(OUT) :: XMOD(NA35) 
      INTEGER*4, INTENT(OUT) :: IERR
      REAL*8 ALPHA,BETA
!
!----------------------------------------------------------------------------------------!
!
!.... classify inverse problem
      IERR = 0
      NVINV = 0
      IF (CINVTYPE == 'PP' .OR. CINVTYPE == 'pp') THEN
         JOB = 1 
         NVINV = 1
      ELSEIF (CINVTYPE == 'SS' .OR. CINVTYPE == 'ss') THEN
         JOB = 2 
         NVINV = 1
      ELSEIF (CINVTYPE == 'PS' .OR. CINVTYPE == 'ps') THEN
         JOB = 3 
         NVINV = 2
      ELSE
         WRITE(*,*) 'genxmod: Cannot do anisotropy yet!'
         IERR = 1
         RETURN
      ENDIF
      IMARK = 0
      DO 1 INPG=1,NNPG
         INPINV = MASKG(INPG) 
         IF (INPINV > 0) THEN
            IMARK = IMARK + 1
            ALPHA = DSQRT( (ECOEFF(INPG,1) + 2.D0*ECOEFF(INPG,2))/DENS(INPG) )
            BETA  = DSQRT( ECOEFF(INPG,2)/DENS(INPG) )
            ILOC1 = (INPINV - 1)*NVINV + 1
            ILOC2 = (INPINV - 1)*NVINV + 2 
            IF (JOB == 1) THEN
               XMOD(ILOC1) = SNGL(ALPHA)
            ELSEIF (JOB == 2) THEN
               XMOD(ILOC1) = SNGL(BETA) 
            ELSEIF (JOB == 3) THEN
               XMOD(ILOC1) = SNGL(ALPHA)
               XMOD(ILOC2) = SNGL(BETA)
            ELSE
               WRITE(*,*) 'genxmod: No anisotropy yet'
               IERR = 1
            ENDIF
         ENDIF
    1 CONTINUE
      IF (IMARK*NVINV /= NA35) THEN
         WRITE(*,*) 'genxmod: Error inpvin /= na35'
         IERR = 2
      ENDIF 
      RETURN
      END  
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GRADPTRS(MYID,MASTER,MYCOMM, MDIM,MEN, NEN,NNPG, LM,     &
                          INV,IERR)
!
!     Here we generate the local CRS gradient pointers such that we can efficiently 
!     perform the sparse matrix vector multiplication grad = adj(F) v.  For more 
!     details look at adjgrad.f90 
!
!     INPUT      MEANING
!     -----      ------- 
!     LM         element node -> global DOF 
!     MASTER     master process ID 
!     MDIM       leading dimension for LM
!     MEN        leading dimension for LM
!     MYCOMM     MPI communicator
!     MYID       process ID in frequency group 
!     NEN        number of element nodes 
!     NNPG       number of anchor nodes 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     ICSC_FDIST local CSC distributed F matrix row pointer
!     IERR       error flag
!     JCSC_FDIST local CSC distributed F matrix column pointer
!     MYGRAD     holds the global gradient inversion nodes
!     NNPGL      number of local gradient inversion points
!     NZ_FDIST   number of non-zeros in distributed F matrix
!  
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), MASTER,MYID,MYCOMM,   &
                               MDIM,MEN, NEN,NNPG
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4, ALLOCATABLE :: IROW(:), ICSC_FDIST(:)
      LOGICAL*4, ALLOCATABLE :: LINIT(:)
      INTEGER*4 STAT(MPI_STATUS_SIZE), INPINV, ICON, INZCOL, IELEM, IR, IROW1, IROW2,    &
                IDOF, IAE, I, INZFL, IVINV, INPG, INPGL, ICOMBIG, INPG0, ILOC,         &
                NEE, NNZCOL, NROW, NWORK, NNPGEST, NPROCS, MYTAG, MYDEST, MYSRC, MPIERR
      LOGICAL*4 LDEBUG
      PARAMETER(LDEBUG = .FALSE.)
!
!----------------------------------------------------------------------------------------!
!
!.... count number of local points in gradeint
      IERR = 0
      inv%NNPGL = 0 
      DO 1 INPG=1,NNPG
         INPINV = inv%MASKG(INPG)
         IF (INPINV > 0) THEN
            IF (inv%IGPART(INPINV) == MYID + 1) inv%NNPGL = inv%NNPGL + inv%NVINV 
         ENDIF
    1 CONTINUE  
      CALL MPI_ALLREDUCE(inv%NNPGL,NNPGEST,1,MPI_INTEGER,MPI_SUM,  &
                         MYCOMM,MPIERR)
      IF (NNPGEST /= inv%NA35) THEN
         IF (MYID == MASTER) WRITE(*,*) 'gradptrs: Error distributed gradient'
      ENDIF
      IF (inv%NNPGL == 0) GOTO 780 !nothing to do here 
!
!.... set space
      NEE = NDIM*NEN
      NWORK = NNPG*inv%NCON*NEE
      ALLOCATE(IROW(inv%NCON*NEE + 1))  !holds non-zeros in row
      ALLOCATE(inv%JCSC_FDIST(inv%NNPGL+1)) !CSC column pointer
      ALLOCATE(ICSC_FDIST(NWORK))   !CSC row pointer
      ALLOCATE(inv%MYGRAD(inv%NNPGL))       !holds gradient points
      inv%JCSC_FDIST(1:inv%NNPGL+1) = 0 
      ICSC_FDIST(1:NWORK) = 0 
      ICOMBIG = 0 
!
!.... loop on local points in inversion and generate a CSC structure
      INZFL = 0 
      INPGL = 1 
      inv%JCSC_FDIST(1) = 1 
      DO 101 INPG=1,NNPG
         INPINV = inv%MASKG(INPG)
         IF (INPINV > 0) THEN !node is in inversion
            IF (inv%IGPART(INPINV) == MYID + 1) THEN !in my partition 
!
!............. investigate anchor nodes connectivity in column 
               INZCOL = 0 
               DO 102 ICON=1,inv%NCON
                  IELEM = inv%MCONN(INPG,ICON)
                  IF (IELEM > 0) THEN
                     DO 103 IAE=1,NEN
                        DO 104 I=1,NDIM
                           IDOF = LM(I,IAE,IELEM)
                           IF (IDOF > 0) THEN !non-zero row? 
                              INZCOL = INZCOL + 1 
                              IROW(INZCOL) = IDOF
                           ENDIF !end check on dof
  104                   CONTINUE !loop on components
  103                CONTINUE !loop on element nodes
                  ELSE
                     GOTO 110 
                  ENDIF !end check on connected element
  102          CONTINUE !loop on connectivity
  110          CONTINUE !break ahead 
               NNZCOL = INZCOL
               CALL ISHELL1(NNZCOL,IROW) !sort 
               IROW(INZCOL+1) = 0 !cheapy way to make next step work 
!
!............. remove repeats and dump to row pointer
               DO 105 INZCOL=1,NNZCOL
                  IF (IROW(INZCOL+1) /= IROW(INZCOL)) THEN
                     INZFL = INZFL + 1
                     IF (INZFL > NWORK) THEN
                        WRITE(*,*) 'ggradptrs: ICSC_FDIST too small!'
                        IERR = 1
                        RETURN
                     ENDIF
                     ICSC_FDIST(INZFL) = IROW(INZCOL)
                  ENDIF
  105          CONTINUE
               inv%MYGRAD(INPGL) = (INPINV - 1)*inv%NVINV + 1
               inv%JCSC_FDIST(INPGL+1) = INZFL + 1
!
!............. copy structure forward in inversion variables
               IROW1 = inv%JCSC_FDIST(INPGL)
               IROW2 = inv%JCSC_FDIST(INPGL+1) - 1
               ICOMBIG = MAX(ICOMBIG,IROW2 - IROW1 + 1)
               DO 106 IVINV=2,inv%NVINV
                  DO 107 IR=IROW1,IROW2
                     INZFL = INZFL + 1
                     IF (INZFL > NWORK) THEN
                        WRITE(*,*) 'gradptrs: ICSC_FDIST too small!'
                        IERR = 1
                        RETURN
                     ENDIF
                     ICSC_FDIST(INZFL) = ICSC_FDIST(IR)
  107            CONTINUE !loop on rows 
                 INPGL = INPGL + 1
                 inv%MYGRAD(INPGL) = (INPINV - 1)*inv%NVINV + IVINV
                 inv%JCSC_FDIST(INPGL+1) = INZFL + 1
  106          CONTINUE !loop on inversion 
               INPGL = INPGL + 1
            ENDIF !end check on partition 
         ENDIF !end check on node in inversion
  101 CONTINUE
      inv%NZ_FDIST = INZFL
      DEALLOCATE(IROW)
!
!.... fidelity check
      INPGL = INPGL - 1
      IF (INPGL /= inv%NNPGL) THEN
         WRITE(*,*) 'gradptrs: Error inpgl != nnpgl'
         IERR = 1
         RETURN
      ENDIF
!
!.... resize
      ALLOCATE(inv%ICSC_FDIST(inv%NZ_FDIST))
      inv%ICSC_FDIST(1:inv%NZ_FDIST) = ICSC_FDIST(1:inv%NZ_FDIST)
      DEALLOCATE(ICSC_FDIST)
  780 CONTINUE !break ahead, nothing to do 
      CALL MPI_ALLREDUCE(ICOMBIG,inv%MBUFRHS,1,MPI_INTEGER,MPI_MAX, &
                         MYCOMM,MPIERR)
!
!.... this is a template for communicating
      IF (LDEBUG) THEN
         INPGL = 0
         INPG0 = 0
         ALLOCATE(LINIT(inv%NA35))
         LINIT(1:inv%NA35) = .FALSE.
         CALL MPI_COMM_SIZE(MYCOMM, NPROCS, MPIERR)
         IF (MYID == MASTER) THEN
            ALLOCATE(IROW(inv%MBUFRHS))
            DO 301 INPG=1,NNPG
               INPINV = inv%MASKG(INPG) !mask h better be in order!
               IF (INPINV > 0) THEN
                  IF (INPINV < INPG0) THEN
                     WRITE(*,*) 'gradptrs: maskg out or order!'
                     IERR = 1
                  ENDIF
                  INPG0 = INPINV
                  MYSRC = inv%IGPART(INPINV) - 1
                  DO 302 IVINV=1,inv%NVINV
                     IF (MYSRC == MYID) THEN !copy my own column 
                        INPGL = INPGL + 1
                        IROW1 = inv%JCSC_FDIST(INPGL)
                        IROW2 = inv%JCSC_FDIST(INPGL+1) - 1
                        NROW = IROW2 - IROW1 + 1
                        IROW(1:NROW) = inv%ICSC_FDIST(IROW1:IROW2)
                     ELSE !receive a slave's row
                        CALL MPI_RECV(NROW,           1,MPI_INTEGER, MYSRC,  &
                                      MPI_ANY_TAG, MYCOMM,STAT,MPIERR)
                        CALL MPI_RECV(IROW(1:NROW),NROW,MPI_INTEGER, MYSRC,  &
                                      MPI_ANY_TAG, MYCOMM,STAT,MPIERR)
                     ENDIF !end check on source
                     ILOC = (INPINV - 1)*inv%NVINV + IVINV
                     IF (LINIT(ILOC)) THEN
                        IERR = 1
                        WRITE(*,*) 'gradptrs: Wrong double init!',ILOC
                     ENDIF
                     LINIT(ILOC) = .TRUE.
                     DO 303 IR=1,NROW
                        WRITE(45+NPROCS,*) IROW(IR)
  303                CONTINUE
                     !LINIT(ILOC) = .true.
  302             CONTINUE !loop on inversion points
               ENDIF !end check if node is in inversion 
  301       CONTINUE !loop on anchor nodes 
            DEALLOCATE(IROW)
            DO 305 ILOC=1,inv%NA35
               IF (.NOT.LINIT(ILOC)) THEN
                  WRITE(*,*) 'gradptrs: Did not initialize inversion variable:',ILOC
                  IERR = 1
               ENDIF
  305       CONTINUE
            DEALLOCATE(LINIT)
            CLOSE(45+NPROCS)
         ELSE !slave process sends
            MYDEST = MASTER
            DO 401 INPG=1,NNPG
               INPINV = inv%MASKG(INPG)
               IF (INPINV > 0) THEN
                  IF (inv%IGPART(INPINV) == MYID+1) THEN !send
                     DO 402 IVINV=1,inv%NVINV
                        INPGL = INPGL + 1
                        IROW1 = inv%JCSC_FDIST(INPGL)
                        IROW2 = inv%JCSC_FDIST(INPGL+1) - 1
                        NROW = IROW2 - IROW1 + 1
                        CALL MPI_SEND(NROW,                         1,MPI_INTEGER,   &
                                      MYDEST,MYTAG,MYCOMM,MPIERR)
                        CALL MPI_SEND(inv%ICSC_FDIST(IROW1:IROW2),NROW,MPI_INTEGER,  &
                                      MYDEST,MYTAG,MYCOMM,MPIERR)
  402                CONTINUE !loop on variables in inversion
                  ENDIF
               ENDIF !end check if INPINV belongs to process 
  401       CONTINUE !loop on anchor nodes
         ENDIF !end check on ID
      ENDIF !end checkon debugging
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE AVGRAD(MGNOD,MNPG, NA35,NNPG, NGNOD,NCON,NVINV,NCASC, &
                         MASKG,MCONN,IENG, GRAD) 
!
!     Averages 
      IMPLICIT NONE
      REAL*4, INTENT(INOUT) :: GRAD(NA35) 
      INTEGER*4, INTENT(IN) :: MCONN(MNPG,*), IENG(MGNOD,*), MASKG(NNPG), &
                               MNPG, MGNOD, NA35, NNPG, NGNOD, NVINV,     &
                               NCON, NCASC 
      REAL*4, ALLOCATABLE :: WORK(:), XAVG(:)
      INTEGER*4, ALLOCATABLE :: NAVG(:) 
      INTEGER*4 ICASC, INPG, JNPG, INPINV, JNPINV, ICON, IA35, IVINV,  &
                IELEM, IA 

      IF (NCASC <= 0) RETURN !early return?
      ALLOCATE(WORK(NA35)) 
      ALLOCATE(XAVG(NVINV)) 
      ALLOCATE(NAVG(NVINV)) 
!
!.... loop on cascades
      DO 100 ICASC=1,NCASC
         CALL SCOPY(NA35,GRAD,1,WORK,1) !copy gradient to workspace
         DO 1 INPG=1,NNPG
            INPINV = MASKG(INPG) 
!
!.......... check anchor node is in gradient
            IF (INPINV > 0) THEN 
               XAVG(1:NVINV) = 0.0 
               NAVG(1:NVINV) = 0
!
!............. loop on element connectivity and average element
               DO 2 ICON=1,NCON
                  IELEM = MCONN(INPG,ICON) 
                  IF (IELEM == 0) GOTO 20 !out of connections 
!
!................ average 
                  DO 3 IA=1,NGNOD
                     JNPG = IENG(IA,IELEM) 
                     IF (JNPG == 0) GOTO 30 !out of anchor nodes 
                     JNPINV = MASKG(JNPG)
                     IF (JNPINV > 0) THEN
                        DO 4 IVINV=1,NVINV
                           IA35 = (JNPINV - 1)*NVINV + IVINV  
                           XAVG(IVINV) = XAVG(IVINV) + GRAD(IA35) 
                           NAVG(IVINV) = NAVG(IVINV) + 1 
    4                   CONTINUE
                     ENDIF !end check if connection is inversion
    3             CONTINUE !loop on anchor nodes
   30             CONTINUE !out of anchor nodes
    2          CONTINUE !Loop on element connectivity
   20          CONTINUE !out of connections
!
!............. average
               DO 5 IVINV=1,NVINV
                  IA35 = (INPINV - 1)*NVINV + IVINV 
                  IF (NAVG(IVINV) > 0) &
                  WORK(IA35) = XAVG(IVINV)/FLOAT(NAVG(IVINV)) 
    5          CONTINUE 
            ENDIF !end check point is in gradient 
    1    CONTINUE
         CALL SCOPY(NA35,WORK,1,GRAD,1) !copy averaged gradient back onto grad
  100 CONTINUE !loop on cascades
      IF (ALLOCATED(WORK)) DEALLOCATE(WORK) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE AVGRAD8(MGNOD,MNPG, NA35,NNPG, NGNOD,NCON,NVINV,NCASC, &
                         MASKG,MCONN,IENG, GRAD) 
!
!     Averages gradient (double precision version)  
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: GRAD(NA35) 
      INTEGER*4, INTENT(IN) :: MCONN(MNPG,*), IENG(MGNOD,*), MASKG(NNPG), &
                               MNPG, MGNOD, NA35, NNPG, NGNOD, NVINV,     &   
                               NCON, NCASC 
      REAL*8, ALLOCATABLE :: WORK(:), XAVG(:)
      INTEGER*4, ALLOCATABLE :: NAVG(:) 
      INTEGER*4 ICASC, INPG, JNPG, INPINV, JNPINV, ICON, IA35, IVINV,  &
                IELEM, IA  

      IF (NCASC <= 0) RETURN !early return?
      ALLOCATE(WORK(NA35)) 
      ALLOCATE(XAVG(NVINV)) 
      ALLOCATE(NAVG(NVINV)) 
!
!.... loop on cascades
      DO 100 ICASC=1,NCASC
         CALL DCOPY(NA35,GRAD,1,WORK,1) !copy gradient to workspace
         DO 1 INPG=1,NNPG
            INPINV = MASKG(INPG)
!
!.......... check anchor node is in gradient
            IF (INPINV > 0) THEN
               XAVG(1:NVINV) = 0.D0
               NAVG(1:NVINV) = 0
!
!............. loop on element connectivity and average element
               DO 2 ICON=1,NCON
                  IELEM = MCONN(INPG,ICON)
                  IF (IELEM == 0) GOTO 20 !out of connections 
!
!................ average 
                  DO 3 IA=1,NGNOD
                     JNPG = IENG(IA,IELEM)
                     IF (JNPG == 0) GOTO 30 !out of anchor nodes 
                     JNPINV = MASKG(JNPG)
                     IF (JNPINV > 0) THEN
                        DO 4 IVINV=1,NVINV
                           IA35 = (JNPINV - 1)*NVINV + IVINV
                           XAVG(IVINV) = XAVG(IVINV) + GRAD(IA35)
                           NAVG(IVINV) = NAVG(IVINV) + 1
    4                   CONTINUE
                     ENDIF !end check if connection is inversion
    3             CONTINUE !loop on anchor nodes
   30             CONTINUE !out of anchor nodes
    2          CONTINUE !Loop on element connectivity
   20          CONTINUE !out of connections
!
!............. average
               DO 5 IVINV=1,NVINV
                  IA35 = (INPINV - 1)*NVINV + IVINV
                  IF (NAVG(IVINV) > 0) &
                  WORK(IA35) = XAVG(IVINV)/DFLOAT(NAVG(IVINV))
    5          CONTINUE
            ENDIF !end check point is in gradient 
    1    CONTINUE
         CALL DCOPY(NA35,WORK,1,GRAD,1) !copy averaged gradient back onto grad
  100 CONTINUE !loop on cascades
      IF (ALLOCATED(WORK)) DEALLOCATE(WORK)
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE AVGMOD4(MGNOD,NELEM,NNPG,NGNOD, NCASC, IENG, XMOD)  
!
!     A simple averaging for for REAL*4 model vectors.  Here we 
!     calculate the average element property, take that to the nodes,
!     then average that with the adjacent element properties.  If 
!     this simple filter is unsatisfactory we can cascade it.  
!     - B. Baker February 2013 
!
!     INPUT     MEANING
!     -----     ------- 
!     IENG      global IEN pointer (as a vector for python) 
!     NCASC     number of cascades in averaging filter
!     NELEM     number of elements
!     NNPG      number of anchor nodes
!     NGNOD     number of anchor nodes on element
!     MGNOD     leading dimension for IENG 
! 
!     OUTPUT    MEANING
!     ------    ------- 
!     XMOD      model to average 
! 
!.... variable declarations
      REAL*4, INTENT(INOUT) :: XMOD(NNPG) 
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*), MGNOD, NGNOD, NNPG, NELEM, NCASC
!.... local variables
      REAL*4, ALLOCATABLE :: XSUM(:), XHIT(:) 
      REAL*4 XAVG  
!
!----------------------------------------------------------------------------------------!
!
!.... loop on the cascade
      ALLOCATE(XSUM(NNPG))
      ALLOCATE(XHIT(NNPG))
      DO 1 ICASC=1,NCASC
!
!....... calculate average element property, sum on nodes
         XSUM(1:NNPG) = 0.0 
         XHIT(1:NNPG) = 0.0 
         DO 2 IELEM=1,NELEM
!.......... calculate average element property
            XAVG = 0.0
            JA = 0
            DO 3 IA=1,NGNOD
               INPG = IENG(IA,IELEM)
               IF (INPG > 0) THEN
                  XAVG = XAVG + XMOD(INPG)
               ELSE
                  GOTO 30
               ENDIF
               JA = JA + 1
    3       CONTINUE !loop on anchor nodes
   30       CONTINUE !break ahead
            XAVG = XAVG/FLOAT(JA)
!.......... set element nodal values to average element value
            DO 4 IA=1,NGNOD
               INDX = (IELEM - 1)*MGNOD + IA
               INPG = IENG(IA,IELEM)
               IF (INPG > 0) THEN
                  XSUM(INPG) = XSUM(INPG) + XAVG
                  XHIT(INPG) = XHIT(INPG) + 1.0
               ELSE
                  GOTO 40
               ENDIF
    4       CONTINUE !loop on anchor nodes
   40       CONTINUE !break ahead
    2    CONTINUE !loop on elements
!....... average nodes
         DO 5 INPG=1,NNPG
            IF (XHIT(INPG) > 0.0) XMOD(INPG) = XSUM(INPG)/XHIT(INPG)
    5    CONTINUE
    1 CONTINUE !loop on cascade 
      DEALLOCATE(XSUM)
      DEALLOCATE(XHIT)
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      REAL*8 FUNCTION DMFRHO(VPR)
!
!     Estimates density (kg/m**3) from Darcy McPhees formula given 
!     P wave velocity in m/s
      REAL*8, INTENT(IN) :: VPR
      REAL*8  PGRCM3, VKMPS
      REAL*8 DCOEF(9) 
      DATA DCOEF /-21893558144627.d-7 ,  30290041149786.d-7,   &
                  -18300316791235.d-7 ,  63062641966165.d-8,   &
                  -13556725156168.d-8 ,  18616741015779.d-9,   &
                  -15948394116079.d-10,  77924797412933.d-12,  &
                  -16626306058716.d-13 /
!
!----------------------------------------------------------------------------------------!
!
      VKMPS = VPR/1000.D0
      IF (VKMPS >= 5.93D0) THEN
         PGRCM3 = 0.9893D0 + 0.2891D0*VKMPS
      ELSEIF (VKMPS >= 5.5D0.AND.VKMPS < 5.93D0) THEN
         PGRCM3 = 0.D0
         DO 1 IJ=1,8
            PGRCM3 = (PGRCM3 + DCOEF(10-IJ))*VKMPS
    1    CONTINUE
         PGRCM3 = PGRCM3 + DCOEF(1)
      ELSE
         PGRCM3 = 1.7407D0*(VKMPS**0.25D0)
      ENDIF
      DMFRHO = PGRCM3*1000.D0
      RETURN
      END
