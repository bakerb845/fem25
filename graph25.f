      SUBROUTINE GENLM(MDIM,MEN,NEN,NNP, NDIM,NELEM,ID,IEN, LM)
! 
!     Generates the LM matrix from the ID matrix.  Hughes 2.10.1 pg 92
! 
!     INPUT      MEANING
!     -----      ------- 
!     ID         determines whether a node is a DOF or not 
!     IEN        holds element node numbers
!     LM         location matrix 
!     MDIM       max spatial dimensions
!     MEN        max element nodes
!     NDIM       number of spatial dimensions
!     NELEM      number of elements
!     NEN        number of element nodes 
!     NNP        number of nodal points in mesh 
!   
!     OUTPUT     MEANING
!     ------     ------- 
!     LM         location  matrix 
! 
!.... variable declarations
      INTEGER*4, INTENT(IN) :: IEN(MEN,NELEM), ID(MDIM,NNP), 
     ;                         MDIM,MEN,NEN,NNP, NDIM,NELEM
      INTEGER*4, INTENT(OUT) :: LM(MDIM,MEN,NELEM)
!
!----------------------------------------------------------------------!
!
!.... loop on elements, nodes, components
      DO 1 IELEM=1,NELEM
         DO 2 IAE=1,NEN
            DO 3 I=1,NDIM
               LM(I,IAE,IELEM) = ID(I,IEN(IAE,IELEM))
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CONNECT(MGNOD,MNP, NNP,NELEM,NGNOD,NCON, IEN,MCONN)
!
!     Determines which elements are sharing anchor nodes 
! 
!     INPUT      MEANING 
!     -----      ------- 
!     IEN        global anchor element node numbers 
!     MGNOD      max number of anchor nodes 
!     MNP        max number of nodal points 
!     NCON       max number of connections 
!     NELEM      number of elements 
!     NGNOD      number of anchor nodes, probably 4 
!  
!     OUTPUT     MEANING
!     ------     ------- 
!     MCONN      holds each nodes element connections 
! 
!.... variable declarations 
      INTEGER*4 IEN(MGNOD,*), MGNOD,MNP, NNP,NELEM,NGNOD,NCON
      INTEGER*4 MCONN(MNP,*)
      INTEGER*4 IELEM,IA,ICON,INP
! 
!----------------------------------------------------------------------!
! 
!.... null out mconn
      DO 1 INP=1,NNP
         DO 2 ICON=1,NCON 
            MCONN(INP,ICON) = 0 
    2    CONTINUE
    1 CONTINUE 
! 
!.... determine who shares global anchor nodes
      DO 3 IELEM=1,NELEM
         DO 4 IA=1,NGNOD
            INP = IEN(IA,IELEM) 
            DO 5 ICON=1,NCON 
               IF (MCONN(INP,ICON).EQ.0) THEN
                  MCONN(INP,ICON) = IELEM
                  GOTO 60
               ENDIF
    5       CONTINUE
   60       CONTINUE
    4    CONTINUE
    3 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION ICMCON(MGNOD,NNPG, NELEM,NGNOD,IENG) 
! 
!     Finds the the max number of connectivities.  Useful for 
!     allocating space
!
!     INPUT     MEANING
!     -----     -------
!     IENG      element node to anchor node pointer
!     MGNOD     leading dimension for IENG
!     NELEM     number of elements
!     NGNOD     number of anchor nodes on an element
!     NNPG      number of anchor nodes
!    
!     RESULT    MEANING
!     ------    ------- 
!     ICMCON    the max number of connectivities 
!
!.... variable declarations 
      INTEGER*4 IENG(MGNOD,*), MGNOD,NNPG, NELEM,NGNOD
!.... local variables
      ALLOCATABLE ICOUNT(:)
      INTEGER*4 ICOUNT
!
!----------------------------------------------------------------------!
!
      ALLOCATE(ICOUNT(NNPG))
      DO 1 I=1,NNPG 
         ICOUNT(I) = 0 
    1 CONTINUE
      DO 2 IELEM=1,NELEM
         DO 3 IA=1,NGNOD
            INODE = IENG(IA,IELEM)
            ICOUNT(INODE) = ICOUNT(INODE) + 1 
    3    CONTINUE
    2 CONTINUE
      ICMCON = ICOUNT(1)
      DO 4 I=2,NNPG
         ICMCON = MAX0(ICMCON,ICOUNT(I))
    4 CONTINUE
      DEALLOCATE(ICOUNT) 
      RETURN
      END 

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GENCRS_QUICK(MDIM,MEN,MNP,MSPACE, NCONF,NELEM,NEN,NDIM,
     ;                   NDOF,MCONN,IEN,LM, NZERO,IRPTR,JCPTR, IERR) 
! 
!     Quick and dirty way to generate a CRS ordering
!
!     INPUT      MEANING
!     -----      ------- 
!     IEN        element node to global node pointer
!     LM         element node to dof pointer
!     MCONN      holds the connectivity for each node
!     MDIM       leading dimension for LM matrix
!     MEN        leading dimension for LM and IEN matrices
!     MNP        leading dimension for MCONN matrix 
!     MSPACE     workspace size for JCPTR
!     NCONF      number of connections per node to search
!     NDIM       number of spatial dimensions
!     NELEM      number of elements 
!     NEN        number of element nodes
!     NDOF       number of degrees of freedom
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag, likely a resize issue
!     IRPTR      CRS row pointer
!     JCPTR      CRS column pointer
!     NZERO      number of non-zeros
!
!.... variable declarations
      INTEGER*4 LM(MDIM,MEN,*), IEN(MEN,*), MCONN(MNP,*) 
      INTEGER*4 IRPTR(NDOF+1), JCPTR(MSPACE)
!.... local variables
      INTEGER*4, ALLOCATABLE :: JZWORK(:), IRPTRW(:), JCPTRW(:), 
     ;                          IPERM(:), PERM(:) 
!
!----------------------------------------------------------------------!
!
!.... generate a connectivity list
      NWORK = NDIM*NCONF*NEN + 1
      ALLOCATE(JZWORK(NWORK))
      ALLOCATE(IRPTRW(NDOF+1))
      ALLOCATE(JCPTRW(MSPACE))
      ALLOCATE(IPERM(NDOF))  
      ALLOCATE(PERM(NDOF)) 
      PERM(1:NDOF) = 0 
!
!.... loop on elements 
      WRITE(*,*) 'gencrs_quick: Generating initial ordering...'
      IERR = 0
      IFOUND = 0
      IZERO = 0  
      IRPTRW(1) = 1
      DO 1 IELEM=1,NELEM
         DO 2 IAE=1,NEN
            INP = IEN(IAE,IELEM)
            DO 3 I=1,NDIM
!
!............. check DOFs
               JRR = 0
               IDOF = LM(I,IAE,IELEM)
               IF (IDOF.GT.0) THEN
                  IF (PERM(IDOF).EQ.0) THEN
                     DO 4 ICONF=1,NCONF
                        JELEM = MCONN(INP,ICONF)
                        IF (JELEM.EQ.0) GOTO 40
                        DO 5 IBE=1,NEN
                           DO 6 J=1,NDIM 
                              JDOF = LM(J,IBE,JELEM)
                              IF (JDOF.GT.0) THEN
                                 JRR = JRR + 1
                                 IF (JRR.GT.NWORK) THEN
                                    WRITE(*,900) 
                                    IERR = 1
                                    GOTO 730
                                 ENDIF 
                                 JZWORK(JRR) = JDOF
                              ENDIF
    6                      CONTINUE !Loop on dimensions
    5                   CONTINUE !Loop on element nodes
    4                CONTINUE !loop on connectivity
   40                CONTINUE !out of elements
                     NRR = JRR
                     JZWORK(NRR+1) = 0
                     CALL ISHELL1(NRR,JZWORK) !sort
                     IZR = 0
                     DO 7 JRR=1,NRR
                        IF (JZWORK(JRR).NE.JZWORK(JRR+1)) THEN
                           IZR = IZR + 1 
                           IZERO = IZERO + 1 
                           IF (IZERO.GT.MSPACE) THEN
                              WRITE(*,901) 
                              IERR = 2
                              GOTO 730
                           ENDIF
                           JCPTRW(IZERO) = JZWORK(JRR)
                        ENDIF
    7                CONTINUE 
                     IFOUND = IFOUND + 1
                     PERM(IDOF) = IFOUND
                     IPERM(IFOUND) = IDOF
                     IRPTRW(IFOUND+1) = IZERO + 1 
                  ENDIF !end check on initialization
               ENDIF  !end check on DOF 
    3       CONTINUE 
    2    CONTINUE
    1 CONTINUE
  730 CONTINUE !break ahead for an error
      DEALLOCATE(JZWORK)
      IF (IERR.NE.0) RETURN
!
!.... sort the permutation matrix in case things arent in order
      IZERO = 0
      IRPTR(1) = 1 
      DO 10 IDOF=1,NDOF
         IROW = PERM(IDOF)
         DO 11 JZERO=IRPTRW(IROW),IRPTRW(IROW+1)-1 
            IZERO = IZERO + 1 
            JCPTR(IZERO) = JCPTRW(JZERO) 
   11    CONTINUE 
         IRPTR(IDOF+1) = IZERO + 1  
   10 CONTINUE
      NZERO = IZERO 
      DEALLOCATE(JCPTRW)
      DEALLOCATE(IRPTRW)
      DEALLOCATE(IPERM)
      DEALLOCATE(PERM) 
!
!.... error messages
  900 FORMAT(' gencrs_quick: Error NWORK is too small')
  901 FORMAT(' gencrs_quick: Error MSPACE is too small')
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GENCRS_DUMB(MDIM,MEN,NDOF,NWORK, NDIM,NEN,NELEM, 
     ;                       LDEBUG,LM, NZERO,IRPTR,JCPTR,IERR) 
      INTEGER*4 LM(MDIM,MEN,*)
      LOGICAL*4 LDEBUG
      INTEGER*4 IRPTR(NDOF+1), JCPTR(NWORK) 
      INTEGER*4, ALLOCATABLE :: IMAT(:,:) 
      ALLOCATE(IMAT(NDOF,NDOF))
      IMAT(1:NDOF,1:NDOF) = 0 
      DO 1 IELEM=1,NELEM
         DO 2 IAE=1,NEN
            DO 3 I=1,NDIM
               IDOF = LM(I,IAE,IELEM)
               IF (IDOF.GT.0) THEN
                  DO 4 IBE=1,NEN
                     DO 5 J=1,NDIM
                        JDOF = LM(J,IBE,IELEM)
                        IF (JDOF.GT.0) IMAT(IDOF,JDOF) = 1 
    5                CONTINUE 
    4             CONTINUE
               ENDIF
    3       CONTINUE 
    2    CONTINUE
    1 CONTINUE
      IF (LDEBUG) OPEN(UNIT=34,FILE='crs_dumb_mat.txt')
      IERR = 0
      IZ = 0 
      IRPTR(1) = 1
      DO 12 I=1,NDOF
         DO 13 J=1,NDOF
            IF (IMAT(I,J).EQ.1) THEN
               IZ = IZ + 1 
               IF (LDEBUG) WRITE(34,*) NDOF+1-I,J
               JCPTR(IZ) = J
             ENDIF
   13    CONTINUE
         IRPTR(I+1) = J + 1
   12 CONTINUE
      NZERO = IZ 
      IF (LDEBUG) CLOSE(34)
      DEALLOCATE(IMAT) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GENID(MDIM,MGNOD,MEN, NNPG,NNP,NELEM,NLXI,NLETA, 
     ;                 NGNOD,NDIM, CDOMAIN,CNNPG, IENG,IEN, 
     ;                 XIPTS,ETAPTS,XLOCS,ZLOCS, NDOF,ID,IERR) 
!
!     Generates the ID matrix.  
!
!     Changelog: S. Roecker discovered a bug where I can miss nodes 
!     on elements in the corners which may have three fixed nodes 
!     versus 2.  This was corrected in February 2014. 
!
!     INPUT      MEANING
!     -----      ------- 
!     CDOMAIN    element domain 'A', 'E', or 'B'
!     CNNPG      anchor node domain
!     ETAPTS     eta interpolation points
!     IEN        fine mesh IEN matrix
!     IENG       global IEN matrix
!     MDIM       max number of spatial dimensions 
!     MEN        max number of element nodes
!     MGNOD      max number of anchor nodes
!     NDIM       number of spatial dimensions
!     NELEM      number of elements in mesh
!     NGNOD      number of anchor nodes
!     NLETA      number of lagrange interpolant points in eta
!     NLXI       number of lagrange inteprolant points in xi 
!     NNP        number of nodal points in mesh 
!     NNPG       number of anchor nodes in mesh 
!     XIPTS      xi interpolation points
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!
!     OUTPUT     MEANING
!     ------     ------- 
!     ID         ID matrix
!     IERR       error flag
!     NDOF       number of degrees of freedom 
!
!.... variable declarations
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER(2), INTENT(IN) :: CNNPG(NNPG)
      CHARACTER(1), INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), 
     ;                      XIPTS(NLXI), ETAPTS(NLETA) 
      INTEGER*4, INTENT(IN) :: IEN(MEN,nelem), IENG(MGNOD,nelem)
      INTEGER*4, INTENT(OUT) :: ID(MDIM,nnp), NDOF, IERR
!.... local variables
      REAL*8, ALLOCATABLE :: SF(:)
      REAL*8 V1(3), V2(3)
      INTEGER*4 IAV(3)
      PARAMETER(TOL = 1.D-10)  
      PARAMETER(LNULL =-5) 
!
!----------------------------------------------------------------------!
!
!.... null out ID array
      IERR = 0
      DO 1 INP=1,NNP
         ID(1:NDIM,INP) = LNULL
    1 CONTINUE
!
!.... loop on anchor nodes
      ALLOCATE(SF(NGNOD)) 
      DO 2 IELEM=1,NELEM
         IF (CDOMAIN(IELEM).EQ.'A') THEN !only fixed B.C.s in PML
            ICHECK = 0
            IA1 = 0
            IA2 = 0
            IAV(1:3) = 0
            DO 3 IA=1,NGNOD
               ILOC = IENG(IA,IELEM)
               IF (CNNPG(ILOC).EQ.'FL' .OR. CNNPG(ILOC).EQ.'FB' .OR.
     ;             CNNPG(ILOC).EQ.'FR') THEN
                  ICHECK = ICHECK + 1
                  IA1 = ILOC 
                  IAV(ICHECK) = IA1
                  DO 4 JA=1,NGNOD
                     JLOC = IENG(JA,IELEM)
                     IF (CNNPG(JLOC).EQ.'FL'.OR.CNNPG(JLOC).EQ.'FB'.OR.
     ;                   CNNPG(JLOC).EQ.'FR') THEN
                        IF (ILOC.NE.JLOC) THEN
                           ICHECK = ICHECK + 1
                           IA2 = JLOC 
                           IAV(ICHECK) = JLOC 
                           IF (ICHECK == 3) GOTO 20 !break ahead
                        ENDIF
                     ENDIF
    4             CONTINUE
                  GOTO 20 !break ahead
               ENDIF !end check on fixed condition
    3       CONTINUE !loop on anchor nodes 
   20       CONTINUE !break ahead
            IF (ICHECK.EQ.0) GOTO 500 !element not on boundary
            !print *, xlocs(ia1),zlocs(ia1) 
!
!.......... investigate points along each potential side
            V1(1:3) = 0.D0
            IA1 = IAV(1)
            DO 5 LOOP=1,2
               IF (ICHECK.GT.1) THEN
                  IA1 = IAV(LOOP)
                  IA2 = IAV(LOOP+1)
                  IF (IA2 == 0) GOTO 500 !only one side to check
                  V1(1) = XLOCS(IA2) - XLOCS(IA1)
                  V1(2) = 0.D0
                  V1(3) = ZLOCS(IA2) - ZLOCS(IA1)
               ELSE !just 1 node on boundary
                  IF (LOOP == 2) GOTO 500 !done
               ENDIF
!
!............. loop on element nodes
               DO 6 ILETA=1,NLETA
                  ETA = ETAPTS(ILETA)
                  DO 7 ILXI=1,NLXI
                     IAE = (ILETA - 1)*NLXI + ILXI 
                     IF (ILXI .EQ.1 .OR. ILXI .EQ.NLXI .OR.
     ;                   ILETA.EQ.1 .OR. ILETA.EQ.NLETA) THEN
                        INP = IEN(IAE,IELEM)
                        XI = XIPTS(ILXI)
                        CALL CSF2DND(NGNOD, XI,ETA, SF,IERR) 
                        IF (IERR.NE.0) THEN
                           WRITE(*,*) 'genid: Error calling csf2dnd'
                           IERR = 1 
                           GOTO 70  
                        ENDIF
                        X = 0.D0
                        Z = 0.D0
                        DO 8 IA=1,NGNOD
                           ILOC = IENG(IA,IELEM)
                           X = X + SF(IA)*XLOCS(ILOC)
                           Z = Z + SF(IA)*ZLOCS(ILOC)
    8                   CONTINUE !loop on anchor nodes
!
!...................... just one anchor node attached to boundary
                        IF (ICHECK.EQ.1) THEN
                           IF (DABS(X - XLOCS(IA1)).LT.TOL .AND.
     ;                         DABS(Z - ZLOCS(IA1)).LT.TOL) THEN
                              ID(1:NDIM,INP) = 0
                              GOTO 400
                           ENDIF
                        ELSE !along line?
                           V2(1) = X - XLOCS(IA1) !line cant be parallel
                           V2(2) = 0.D0
                           V2(3) = Z - ZLOCS(IA1)
                           DET = V1(1)*V2(3) - V1(3)*V2(1)
                           IF (DABS(DET).LT.TOL) THEN
                              ID(1:NDIM,INP) = 0
                              GOTO 400
                           ENDIF 
                        ENDIF !end check on whether to check
                     ENDIF !end check on element side
  400                CONTINUE !break ahead for new xi point
    7             CONTINUE !loop on xi points
    6          CONTINUE !loop on eta points
    5       CONTINUE !loop on points to check
  500       CONTINUE !break ahead for new element 
         ENDIF !end check on domain
    2 CONTINUE
   70 CONTINUE !break ahead
      DEALLOCATE(SF) 
      IF (IERR.NE.0) RETURN
!
!.... easy now, we just go and count
      IDOF = 0
      DO 11 IELEM=1,NELEM
         DO 12 ILETA=1,NLETA
            DO 13 ILXI=1,NLXI
               IAE = (ILETA - 1)*NLXI + ILXI
               INP = IEN(IAE,IELEM)
               DO 14 I=1,NDIM
                  IF (ID(I,INP).EQ.LNULL) THEN
                     IDOF = IDOF + 1
                     ID(I,INP) = IDOF
                  ENDIF
   14          CONTINUE !loop on spatial dimensions
   13       CONTINUE !loop on xi
   12    CONTINUE !loop on eta 
   11 CONTINUE
      NDOF = IDOF
!
!.... fidelity check
      IERR = 0
      DO 21 INP=1,NNP
         DO 22 I=1,NDIM
            IF (ID(I,INP).EQ.LNULL) THEN
               WRITE(*,*) 'genid: Error generating ID',I,INP
               IERR = 1
            ENDIF
   22    CONTINUE
   21 CONTINUE 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GENIEN(MGNOD,MEN,NNPG,NELEM,NLXI,NLETA,NGNOD, LVERB,
     ;                  NEN,IENG,XIPTS,ETAPTS,XLOCS,ZLOCS, 
     ;                  NNP,IEN,IERR)
! 
!     Generates the IEN vector using the geoemtry.  First we 
!     locate the physical locations of points on the side of an 
!     element and number them.  Next we number the points inside the 
!     element 
! 
!.... variable declarations 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8, INTENT(IN) :: XLOCS(NNPG),ZLOCS(NNPG), 
     ;                      XIPTS(NLXI),ETAPTS(NLETA)
      INTEGER*4, INTENT(IN) :: IENG(MGNOD,*)
      LOGICAL*4, INTENT(IN) :: LVERB
      INTEGER*4, INTENT(OUT) :: IEN(MEN,*),NNP,IERR
      REAL*8, ALLOCATABLE :: SF(:),XNODES(:),ZNODES(:)
      LOGICAL*4 LSORTED
      PARAMETER(TOL = 1.D-8)
      PARAMETER(LNULL =-5)
! 
!----------------------------------------------------------------------!
! 
!.... calculate locations at every node on element side 
      MNODES = NELEM*(2*(NLXI - 1) + 2*NLETA)
      ALLOCATE(XNODES(MNODES))
      ALLOCATE(ZNODES(MNODES))
      ALLOCATE(SF(NGNOD))
      XMAX = MAXVAL(XLOCS)
      ZMAX = MAXVAL(ZLOCS)
      TOOBIG = DMAX1(XMAX,ZMAX)*100.D0
      TOOBIG = HUGE(1.D0) 
      XNODES(1:MNODES) = TOOBIG
      ZNODES(1:MNODES) = TOOBIG
      IINT = 0
      INODES = 0 
      DO 1 IELEM=1,NELEM
         DO 2 ILETA=1,NLETA
            ETA = ETAPTS(ILETA)
            DO 3 ILXI=1,NLXI
               IAE = (ILETA - 1)*NLXI + ILXI 
               IF (ILXI .EQ.1 .OR. ILXI .EQ.NLXI  .OR. 
     ;             ILETA.EQ.1 .OR. ILETA.EQ.NLETA) THEN 
                  XI = XIPTS(ILXI)
                  CALL CSF2DND(NGNOD, XI,ETA, SF,IERR)
                  IF (IERR.NE.0) THEN
                     WRITE(*,*) 'genien: Error calling csf2dnd 1'
                     IERR = 1
                     GOTO 730
                  ENDIF
                  X = 0.D0
                  Z = 0.D0
                  DO 4 IA=1,NGNOD
                     ILOC = IENG(IA,IELEM)
                     X = X + SF(IA)*XLOCS(ILOC)
                     Z = Z + SF(IA)*ZLOCS(ILOC)
    4             CONTINUE
!
!................ kick out the repeats
                  DO 5 JNODES=1,INODES
                     IF (DABS(XNODES(JNODES) - X).LT.TOL .AND.
     ;                   DABS(ZNODES(JNODES) - Z).LT.TOL) GOTO 45 
    5             CONTINUE 
                  INODES = INODES + 1
                  XNODES(INODES) = X
                  ZNODES(INODES) = Z
   45             CONTINUE !break ahead
               ENDIF !end check on side
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      NNODES = INODES
! 
!.... heap sort on x locations
      IF (LVERB) WRITE(*,*) 'genien: Heapsorting x locations...'
      CALL DHPSRT2(NNODES,XNODES,ZNODES)
! 
!.... shell sort z locations 
      IF (LVERB) WRITE(*,*) 'genien: Shell sorting z locations...'
      LSORTED = .FALSE.
      JBEG = 1
      DO 10 INODES=1,NNODES
! 
!....... special case
         IF (INODES+1.EQ.NNODES) THEN
            IF (JBEG.EQ.NNODES) THEN
               GOTO 50 !we're done
            ELSE
               JEND = NNODES
               JVAR = JEND - JBEG + 1
               CALL SHELL2(JVAR,ZNODES(JBEG:JEND),XNODES(JBEG:JEND))
            ENDIF
            GOTO 50
         ENDIF
         !IF (XNODES(INODES).NE.XNODES(INODES+1)) THEN
         IF (DABS(XNODES(INODES+1) - XNODES(INODES)).GT.TOL) THEN
            JEND = INODES
            JVAR = JEND - JBEG + 1
            IF (JVAR.GT.1) THEN
               CALL SHELL2(JVAR,ZNODES(JBEG:JEND),XNODES(JBEG:JEND))
            ENDIF
            JBEG = JEND + 1
            LSORTED = .TRUE.
            GOTO 100
         ENDIF
         LSORTED = .FALSE.
  100    CONTINUE
   10 CONTINUE
   50 CONTINUE
      NNPS = NNODES
      XNODES(NNPS+1) = TOOBIG
      ZNODES(NNPS+1) = TOOBIG
      IF (LVERB) WRITE(*,*) 'genien: Assigning nodal numbers...'
! 
!.... create IEN vector
      IINT = NNPS !start interior counter at nnp sides
      IEN(1:NEN,1:NELEM) = LNULL 
      DO 30 IELEM=1,NELEM
         DO 31 ILETA=1,NLETA
            ETA = ETAPTS(ILETA)
            DO 32 ILXI=1,NLXI
               IAE = (ILETA - 1)*NLXI + ILXI 
               IF (ILXI .EQ.1 .OR. ILXI .EQ.NLXI  .OR. 
     ;             ILETA.EQ.1 .OR. ILETA.EQ.NLETA) THEN 
                  XI = XIPTS(ILXI)
                  CALL CSF2DND(NGNOD,XI,ETA, SF,IERR)
                  IF (IERR.NE.0) THEN
                     WRITE(*,*) 'genien: Error calling csf2dnd 2'
                     IERR = 2
                     GOTO 730
                  ENDIF
                  X = 0.D0
                  Z = 0.D0
                  DO 33 IA=1,NGNOD
                     ILOC = IENG(IA,IELEM)
                     X = X + SF(IA)*XLOCS(ILOC)
                     Z = Z + SF(IA)*ZLOCS(ILOC)
   33             CONTINUE
!................ search for x location
                  JBEG = IBSECT8(NNPS+1,1,TOL,X,XNODES)
                  IF (JBEG.LT.0) THEN
                     WRITE(*,*) 'genien: Error calling ibsect8'
                     IERR = 3 
                     GOTO 730
                  ENDIF
                  IF (DABS(XNODES(JBEG) - X).GT.TOL) THEN
                     WRITE(*,*) 'genien: Could not locate x point'
                     IERR = 4
                     GOTO 730
                  ENDIF
!
!................ get z through linear search
                  DO 34 J=JBEG,NNPS
                     IF (DABS(XNODES(J) - X).LT.TOL .AND. 
     ;                   DABS(ZNODES(J) - Z).LT.TOL) THEN
                        IEN(IAE,IELEM) = J
                        GOTO 60 
                     ENDIF
   34            CONTINUE
                 IERR = 4 
                 WRITE(*,*) 'genien: Could not locate z point'
                 GOTO 730 
   60            CONTINUE !break forward
               ELSE !element interior
                  IINT = IINT + 1
                  IEN(IAE,IELEM) = IINT
               ENDIF
   32       CONTINUE
   31    CONTINUE
   30 CONTINUE
  730 CONTINUE
      DEALLOCATE(XNODES)
      DEALLOCATE(ZNODES)
      DEALLOCATE(SF)
      IF (IERR.GT.0) RETURN
      NNINT = IINT - NNPS
      NNP = NNINT + NNPS !#number of nodal points = interior + sides 
! 
!.... fidelity check
      IERR = 0
      DO 40 IELEM=1,NELEM
         DO 41 IAE=1,NEN
            IF (IEN(IAE,IELEM).EQ.LNULL) THEN
               WRITE(*,*) 'genien: Error ien not classified'
               IERR = 1
            ENDIF
   41    CONTINUE
   40 CONTINUE
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CRS2COO(NZERO,N, IRPTR,JCPTR, IRN,JCN)
      INTEGER*4 IRPTR(N+1),JCPTR(NZERO)
      INTEGER*4 IRN(NZERO),JCN(NZERO)
      DO 1 I=1,N
         JSTRT = IRPTR(I)
         JSTOP = IRPTR(I+1) - 1   
         DO 2 J=JSTRT,JSTOP
            IRN(J) = I   
            JCN(J) = JCPTR(J)
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CRS2ADJ(NZERO,N, IRPTR,JCPTR, XADJ,ADJNCY)
      INTEGER*4 IRPTR(N+1),JCPTR(NZERO),NZERO,N
      INTEGER*4 XADJ(N+1), ADJNCY(NZERO-N)
      XADJ(1) = 1
      ITER = 0
      DO 1 I=1,N
         JBEG = IRPTR(I)
         JEND = IRPTR(I+1) - 1
         DO 2 J=JBEG,JEND
            IF (JCPTR(J).NE.I) THEN
               ITER = ITER + 1
               ADJNCY(ITER) = JCPTR(J)
            ENDIF
    2    CONTINUE
         XADJ(I+1) = ITER + 1
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE ICNDZ_LOC(MYID,NFLAG,NDOF, PART,IRPTR, NDOFL,NZLOC)
!
!     Counts number of local DOFs and non-zeros
!
!     INPUT     MEANING 
!     -----     ------- 
!     IRPTR     CRS row pointer
!     MYID      process ID to match to part
!     NDOF      number of degrees of freedom
!     NFLAG     numbering flag used in Metis (0) C (1) Fortran (default)
!     PART      holds partition numbers for each DOF
!
!     RESULT    MEANING
!     ------    ------- 
!     NDOFL     number of local DOFs
!     NZLOC     number of local non-zeros
!
      INTEGER*4, INTENT(IN) :: IRPTR(NDOF+1),PART(NDOF), MYID,NFLAG,NDOF
      INTEGER*4, INTENT(OUT) :: NDOFL, NZLOC
!
!----------------------------------------------------------------------!
!
      IF (NFLAG.EQ.0) THEN !C process IDs match C partition from Metis
         MYTID = MYID
      ELSE !offset C numbered process IDs w/ Fortran part from Metis
         MYTID = MYID + 1
      ENDIF
      NDOFL = 0 
      NZLOC = 0
      DO 1 IDOF=1,NDOF
         IF (PART(IDOF).EQ.MYTID) THEN
            NDOFL = NDOFL + 1
            NZROW = IRPTR(IDOF+1) - IRPTR(IDOF)
            NZLOC = NZLOC + NZROW
         ENDIF
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PART2CRSL(NZERO,NZLOC, NDOF,NDOFL, MYID,NFLAG, 
     ;                     IRPTR,JCPTR,PART, MYDOFS,IRPTR_LOC,JCPTR_LOC)
!
!     Converts the global CRS storage to the local CRS storage
!
!     INPUT      MEANING
!     -----      ------- 
!     IRPTR      CRS global row pointer
!     JCPTR      CRS global column pointer
!     MYID       process ID
!     NDOF       number of degrees of freedom
!     NDDOFL     number of local DOFs
!     NFLAG      numbering flag, (0) C, (1) Fortran used in Metis
!     NZERO      number of non-zeros in global matrix
!     NZLOC      number of local non-zeros
!     PART       holds partition number of each DOF
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IRPTR_LOC  local CRS row pointer
!     JCPTR_LOC  local CRS column pointer
!     MYDOFS     global DOF numbers 
!  
!.... variable declarations
      INTEGER*4, INTENT(IN) :: IRPTR(NDOF+1),JCPTR(NZERO), PART(NDOF),
     ;                         NZERO,NZLOC,NDOF,NDOFL, MYID,NFLAG
      INTEGER*4, INTENT(OUT) :: IRPTR_LOC(NDOFL+1), JCPTR_LOC(NZLOC), 
     ;                          MYDOFS(NDOFL)
!
!----------------------------------------------------------------------!
!
      IF (NFLAG.EQ.0) THEN !C process IDs match C partition from Metis
         MYTID = MYID
      ELSE !offset C numbered processes IDs w Fortran part from Metis
         MYTID = MYID + 1
      ENDIF
      IDOFL = 0
      IZLOC = 0
      IRPTR_LOC(1) = 1
      DO 1 IDOF=1,NDOF
         IF (PART(IDOF).EQ.MYTID) THEN 
            IDOFL = IDOFL + 1
            MYDOFS(IDOFL) = IDOF
            J1 = IRPTR(IDOF)
            J2 = IRPTR(IDOF+1) - 1
            NZROW = J2 - J1 + 1
            L1 = IZLOC + 1
            L2 = L1 + NZROW - 1  
            JCPTR_LOC(L1:L2) = JCPTR(J1:J2)
            IZLOC = IZLOC + NZROW
            IRPTR_LOC(IDOFL+1) = IRPTR_LOC(IDOFL) + NZROW
         ENDIF
    1 CONTINUE
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GRAPHREORD(NZERO,N,IRPTR,JCPTR, PERM) 
! 
!     Reorders numbering on a mesh to reduce fill in 
      INTEGER*4 IRPTR(N+1),JCPTR(NZERO) 
      INTEGER*4 PERM(N)
      ALLOCATABLE XADJ(:),ADJNCY(:),IPERM(:)
      INTEGER*4 XADJ,ADJNCY,IOPT(5),IOPT8(8),KWTF,NUMFLAG
! 
!----------------------------------------------------------------------!
! 
!.... generate compressed storage format of george and liu
      ALLOCATE(XADJ(N+1))
      ALLOCATE(ADJNCY(NZERO-N))
      CALL CRS2ADJ(NZERO,N, IRPTR,JCPTR, XADJ,ADJNCY)
      KWTF = 0   
      IOPT(1) = 0 !defaults
      NUMFLAG = 1 !Fortran numbering
! 
!     NB: IPERM and PERM are reversed to be consistent w/ George and Liu
      ALLOCATE(IPERM(N))
      IOPT8(1) = 0 !defaults
      CALL METIS_NodeND(N,XADJ,ADJNCY,NUMFLAG,IOPT8,IPERM,PERM)
      DEALLOCATE(XADJ)
      DEALLOCATE(ADJNCY)
      DEALLOCATE(IPERM) 
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE CRS2COOLOC(NZLOC,NDOFL, MYDOFS,IRPTR_LOC,JCPTR_LOC, 
     ;                      IRN_LOC,JCN_LOC) 
! 
!     Converts the local CRS to the local COO structure
! 
!     INPUT      MEANING 
!     -----      ------- 
!     IRPTR_LOC  local non-zero pointers in CRS format
!     JCPTR_LOC  global column numbers in CRS format
!     MYDOFS     global degree of freedom numbers this process holds 
!     NDOFL      number of local degrees of freedom 
!     NZLOC      number of local non-zeros
!    
!     OUTPUT     MEANING 
!     ------     -------
!     IRN_LOC    points to global rows 
!     JCN_LOC    points to global columns 
! 
!.... variable delcarations
      INTEGER*4, INTENT(IN) :: IRPTR_LOC(NDOFL+1), JCPTR_LOC(NZLOC), 
     ;                         MYDOFS(NDOFL), NZLOC,NDOFL
      INTEGER*4, INTENT(OUT) :: IRN_LOC(NZLOC), JCN_LOC(NZLOC)
! 
!---------------------------------------------------------------------!
! 
      DO 1 IDOFL=1,NDOFL
         JBEG = IRPTR_LOC(IDOFL) 
         JEND = IRPTR_LOC(IDOFL+1) - 1   
         IROW = MYDOFS(IDOFL) 
         DO 2 J=JBEG,JEND
            IRN_LOC(J) = IROW 
            JCN_LOC(J) = JCPTR_LOC(J) 
    2    CONTINUE
    1 CONTINUE 
      RETURN 
      END  
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GENCRS_LOC(MZLOC,NZERO,NDOF, NDOFL,MYDOFS, IRPTR,JCPTR,
     ;                      NZLOC,IRPTR_LOC,JCPTR_LOC)
! 
!     Generates the local compressed row storage 
! 
!     INPUT      MEANING
!     -----      ------- 
!     IRPTR      global CSR row pointer 
!     JCPTR      global CSR column pointer 
!     MYDOFS     vector holding dofs in partition
!     MZLOC      max number of local non-zeros
!     NDOF       number of degrees of freedom 
!     NDOFL      number of local degrees of freedom 
!     NZERO      number of non-zeros in global matrix 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     NZ_LOC     number of local non-zeros
!     IRPTR_LOC  local row pointer 
!     JCPTR_LOC  local column pointer 
! 
      INTEGER*4, INTENT(IN) :: IRPTR(NDOF+1),JCPTR(NZERO),MYDOFS(NDOFL),
     ;                         MZLOC,NZERO,NDOF,NDOFL
      INTEGER*4, INTENT(OUT) :: IRPTR_LOC(NDOFL+1), JCPTR_LOC(MZLOC), 
     ;                          NZLOC 
! 
!----------------------------------------------------------------------!
!
      NZLOC = 0 
      JZERO = 0 
      IRPTR_LOC(1) = 1 
      DO 1 IDOFL=1,NDOFL
         IDOF = MYDOFS(IDOFL)
         JBEG = IRPTR(IDOF)  
         JEND = IRPTR(IDOF+1) - 1 
         NZ1 = JEND - JBEG + 1 
         NZLOC = NZLOC + NZ1 
         IRPTR_LOC(IDOFL+1) = IRPTR_LOC(IDOFL) + NZ1 
         JBEGL = IRPTR_LOC(IDOFL)
         JENDL = IRPTR_LOC(IDOFL+1) - 1 
         DO 2 J=JBEG,JEND 
            JZERO = JZERO + 1 
            JCPTR_LOC(JZERO) = JCPTR(J) 
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION ICELEML(MYID, MDIM,MEN,NDOF, NELEM, NFLAG, 
     ;                           NEN,NDIM, PART,LM) 
!
!     Counts local elements 
!
!     INPUT     MEANING
!     -----     -------
!     LM        location matrix
!     MDIM      leading dimnesion for LM
!     MEN       leading dimension for LM
!     MYID      process ID 
!     NDIM      number of spatial dimensions
!     NDOF      number of degrees of freedom
!     NELEM     number of elements 
!     NEN       number of element nodes
!     NFLAG     numbering flag used in Metis (0) C (1) Fortran (default)
!     PART      holds the partition IDs
!
!     RESULT    MEANING
!     ------    ------- 
!     ICELEML   number of elements local to process MYID 
!
!.... variable declarations
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), PART(NDOF), MDIM,MEN, 
     ;           NDOF,NELEM,NFLAG,NEN,NDIM, MYID 
!
!----------------------------------------------------------------------!
!
      IF (NFLAG.EQ.0) THEN 
         MYTID = MYID 
      ELSE 
         MYTID = MYID + 1
      ENDIF
      ICELEML = 0  
      DO 1 IELEM=1,NELEM
         DO 2 IAE=1,NEN
            DO 3 I=1,NDIM
               IDOF = LM(I,IAE,IELEM)
               IF (IDOF.GT.0) THEN
                  IF (PART(IDOF).EQ.MYTID) THEN
                     ICELEML = ICELEML + 1
                     GOTO 20
                  ENDIF !end check
                ENDIF !end check on DOF 
    3       CONTINUE !loop on dimensions 
    2   CONTINUE !loop on element nodes
   20   CONTINUE ! 
    1 CONTINUE !loop on elements
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GENELEML(MYID, MDIM,MEN,NDOF, NELEM,NELEML, NFLAG,
     ;                    NEN,NDIM, PART,LM, MYELEM) 
!
!     Saves local elements into MYELEM
!
!     INPUT     MEANING
!     -----     -------
!     LM        location matrix
!     MDIM      leading dimnesion for LM
!     MEN       leading dimension for LM
!     MYID      process ID 
!     NDIM      number of spatial dimensions
!     NDOF      number of degrees of freedom
!     NELEM     number of elements 
!     NELEML    number of local elements
!     NEN       number of element nodes
!     NFLAG     numbering flag used in Metis (0) C (1) Fortran (default)
!     PART      holds the partition IDs
!
!     OUTPUT    MEANING
!     ------    ------- 
!     MYELEM    holds process' elements
!
!.... variable declarations
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), PART(NDOF), MDIM,MEN, 
     ;           NDOF,NELEM,NELEML,NFLAG,NEN,NDIM, MYID
      INTEGER*4, INTENT(OUT) :: MYELEM(NELEML)
!
!----------------------------------------------------------------------!
!
      IF (NFLAG.EQ.0) THEN
         MYTID = MYID
      ELSE
         MYTID = MYID + 1
      ENDIF
      IELEML = 0 
      DO 1 IELEM=1,NELEM
         DO 2 IAE=1,NEN
            DO 3 I=1,NDIM
               IDOF = LM(I,IAE,IELEM)
               IF (IDOF.GT.0) THEN
                  IF (PART(IDOF).EQ.MYTID) THEN
                     IELEML = IELEML + 1 
                     MYELEM(IELEML) = IELEM
                     GOTO 20
                  ENDIF !end check
                ENDIF !end check on DOF 
    3       CONTINUE !loop on dimensions 
    2   CONTINUE !loop on element nodes
   20   CONTINUE ! 
    1 CONTINUE !loop on elements
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE COLOR_MESH(MDIM,MEN,NDOF,NZERO, NDIM,NEN,NELEM, 
     ;                      NPARTS, LWGT, CDOMAIN,IRPTR,JCPTR,LM, PART) 
! 
!     Colors the mesh.   
! 
!     INPUT      MEANING
!     -----      ------- 
!     CDOMAIN    hodls the element domain 'A' -> absorbing
!     IRPTR      CRS row poitner
!     JCPTR      CRS column pointer
!     LM         location matrix
!     LWGT       True -> dofs in PML will be weighted
!     MDIM       leading dimension for LM
!     MEN        leading dimension for LM
!     NDIM       number of spatial dimensions
!     NDOF       number of degrees of freedom
!     NELEM      number of element in mesh
!     NEN        number of element nodes
!     NPARTS     number of parts to subdivide mesh
!     NZERO      number of non-zeros
!
!     OUTPUT     MEANING
!     ------     ------- 
!     PART       holds each DOFs partition number 
! 
!.... variable declarations
      CHARACTER*1, INTENT(IN) :: CDOMAIN(NELEM)
      INTEGER*4, INTENT(IN) ::  LM(MDIM,MEN,*),IRPTR(NDOF+1), 
     ;           JCPTR(NZERO), MDIM,MEN,NDOF,NZERO, NDIM,NEN,NELEM, 
     ;           NPARTS 
      LOGICAL*4, INTENT(IN) :: LWGT
      INTEGER*4, INTENT(OUT) :: PART(NDOF) 
!.... local variables
      INTEGER*4, ALLOCATABLE :: XADJ(:), ADJNCY(:), VWGT(:) 
      INTEGER*4 OPTIONS(5), WGTFLAG 
!
!----------------------------------------------------------------------!
!
!.... 1 part is pretty easy 
      IF (NPARTS.EQ.1) THEN
         PART(1:NDOF) = 1   
      ELSE
         IF (LWGT) THEN
            ALLOCATE(VWGT(NDOF)) 
            VWGT(1:NDOF) = 1 
            DO 11 IELEM=1,NELEM 
               IF (CDOMAIN(IELEM).EQ.'A') THEN !abosrbing domain
                  DO 12 IAE=1,NEN
                     DO 13 I=1,NDIM
                        IDOF = LM(I,IAE,IELEM)
                        IF (IDOF.GT.0) VWGT(IDOF) = 2
   13                CONTINUE
   12             CONTINUE
               ENDIF !end check on domain 
   11       CONTINUE !loop on elements
            WGTFLAG = 2 !weights on verticles only
         ELSE
            WGTFLAG = 0
         ENDIF
         NSADJ = NZERO - NDOF
         ALLOCATE(XADJ(NDOF+1))
         ALLOCATE(ADJNCY(NSADJ))
         CALL CRS2ADJ(NZERO,NDOF, IRPTR,JCPTR, XADJ,ADJNCY)
         N = NDOF
         OPTIONS(1) = 0 !defaults
         NUMFLAG = 1 !Fortran numbering
         NP = NPARTS !number of parts
         IF (LWGT) THEN
            WRITE(*,*) 'color_mesh: Coloring PML weighted mesh...'
            CALL METIS_PartGraphKway(N,XADJ,ADJNCY, VWGT,NULL,WGTFLAG,
     ;                               NUMFLAG,NP,OPTIONS, NEDGE,PART)
         ELSE
            WRITE(*,*) 'color_mesh: Coloring unweighted mesh...'
            CALL METIS_PartGraphKway(N,XADJ,ADJNCY, NULL,NULL,WGTFLAG,
     ;                               NUMFLAG,NP,OPTIONS, NEDGE,PART)
         ENDIF
         DEALLOCATE(XADJ)
         DEALLOCATE(ADJNCY)
         IF (ALLOCATED(VWGT)) DEALLOCATE(VWGT)
      ENDIF !end check on nparts
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      INTEGER*4 FUNCTION IENZERO(MNP,MDIM, NNP,NDIM,NCONF,NEN, ID,MCONN)
!
!     A crude over estimate of the number of non-zeros 
!
!     INPUT      MEANING
!     -----      ------- 
!     ID         holds whether each node and its components are DOFs
!     MCONN      holds nodal points connectivities
!     MDIM       leading dimension for ID
!     MNP        leading dimension for MCONN
!     NCONF      number of connetions
!     NDIM       number of components in solution
!     NEN        number of element nodes
!     NNP        number of nodal points in mesh
!  
!     OUTPUT     MEANING
!     ------     ------- 
!     IENZERO    over estimate of number of non-zeros
!
!.... variable declarations 
      INTEGER*4, INTENT(IN) :: MCONN(MNP,NCONF), ID(MDIM,NNP), 
     ;                         MNP,MDIM, NNP,NDIM,NCONF,NEN
!
!----------------------------------------------------------------------!
!
!.... loop on nodal points 
      NEE = NDIM*NEN !length of a subrow
      IENZERO = 0
      DO 1 INP=1,NNP
         DO 2 I=1,NDIM !spatial dimensions in solution
            IF (ID(I,INP).GT.0) THEN !check its a DOF
               DO 3 ICONF=1,NCONF
                  IF (MCONN(INP,ICONF).GT.0) THEN
                     IF (ICONF.EQ.1) THEN
                        IENZERO = IENZERO + NEE
                     ELSE
                        IENZERO = IENZERO + NEE - 1
                     ENDIF
                  ELSE
                     GOTO 30
                  ENDIF
    3          CONTINUE
   30          CONTINUE !break ahead
            ENDIF
    2    CONTINUE !Loop on spatial dimensions 
    1 CONTINUE
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GIDOFFS(MDIM,MEN,MGNOD, NNPG,NLXI,NLETA, 
     ;                   NELEM,NGNOD,NDIM, CNNPG,IENG,LM,
     ;                   XIPTS,ETAPTS,XLOCS,ZLOCS, NFS,IDOF_FS,IERR)
!
!     Calculates the x locations at the free surface and their x 
!     locations.  The x locations are sorted. 
!
!     INPUT      MEANING
!     -----      ------- 
!     CNNPG      anchor node type 
!     ETAPTS     eta interpolationpoints
!     IENG       anchor node IEN pointer
!     LM         maps element nodes to global DOFs
!     MDIM       leading dimension
!     MEN        max number of element nodes
!     MGNOD      max number of nachor nodes
!     NDIM       number of components
!     NELEM      number of elements in mesh
!     NGNOD      number of anchor nodes on an element
!     NLETA      number eta interpolant points
!     NLXI       number of xi interpolant points
!     NNPG       number of anchor nodes in mesh
!     XIPTS      xi interpolation points
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IDOF_FS    DOFs at free surface fine nodes
!     IERR       error flag
!     NFS        number of nodes at free surface
!    
!.... variable declarations
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER(2), INTENT(IN) :: CNNPG(NNPG)
      REAL*8, INTENT(IN) :: XLOCS(NNPG),ZLOCS(NNPG), 
     ;                      XIPTS(NLXI),ETAPTS(NLETA)
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*),
     ;           MDIM,MEN,MGNOD, NNPG,NLXI,NLETA,
     ;           NELEM,NGNOD,NDIM
      INTEGER*4, INTENT(OUT) :: IERR 
      INTEGER*4, INTENT(OUT) :: IDOF_FS(MDIM,*) 
!.... local variables
      REAL*8, ALLOCATABLE :: SF(:), XFS(:) 
      INTEGER*4, ALLOCATABLE :: IDOF_FSW(:) 
      REAL*8 V1(3), V2(3), TOL 
      LOGICAL*4 LSFREE, LFS  
      PARAMETER(TOL = 2.22D-7)
!
!----------------------------------------------------------------------!
!
!.... loop on geometry
      IERR = 0
      JFS = 0
      ALLOCATE(SF(NGNOD)) 
      ALLOCATE(XFS(NLXI*NNPG))      !over-estimate
      ALLOCATE(IDOF_FSW(NLXI*NNPG)) !over-estimate 
      DO 11 IELEM=1,NELEM
         KEEP = 0 
         ILOC1 = 0 
         ILOC2 = 0 
         LSFREE = .FALSE.
         DO 101 IA=1,NGNOD
            ILOC = IENG(IA,IELEM)
            IF (CNNPG(ILOC).EQ.'FS') THEN !free surface?
               ILOC1 = ILOC
               LSFREE = .TRUE.
               KEEP = KEEP + 1 
               DO 102 JA=1,NGNOD
                  JLOC = IENG(JA,IELEM)
                  IF (JLOC.NE.ILOC .AND. CNNPG(JLOC).EQ.'FS') THEN
                     ILOC2 = JLOC
                     KEEP = KEEP + 1 
                     GOTO 40  
                  ENDIF !two nodes on free surface
  102          CONTINUE !loop on anchor nodes
               GOTO 40
            ENDIF !end check if node is on free surface
  101    CONTINUE !loop on anchor nodes
   40    CONTINUE !break ahead
         V1(1:3) = 0.D0
         IF (KEEP.EQ.2) THEN  
            V1(1) = XLOCS(ILOC2) - XLOCS(ILOC1)
            V1(2) = 0.D0
            V1(3) = ZLOCS(ILOC2) - ZLOCS(ILOC1)
         ENDIF
!
!....... begin loop on geometry of element 
         DO 12 ILETA=1,NLETA
            ETA = ETAPTS(ILETA)
            DO 13 ILXI=1,NLXI
               IAE = (ILETA - 1)*NLXI + ILXI
               XI = XIPTS(ILXI)
               CALL CSF2DND(NGNOD,XI,ETA, SF, IERR)
               IF (IERR.NE.0) THEN
                  WRITE(*,*) 'locrec2d: Error calling csf2dnd'
                  GOTO 100
               ENDIF
               X = 0.D0
               Z = 0.D0
               DO 14 IA=1,NGNOD
                  ILOC = IENG(IA,IELEM)
                  X = X + SF(IA)*XLOCS(ILOC)
                  Z = Z + SF(IA)*ZLOCS(ILOC)
   14          CONTINUE !loop on anchor ndoes
               LFS = .FALSE.
               IF (ILXI .EQ.1 .OR. ILXI .EQ.NLXI .OR.
     ;             ILETA.EQ.1 .OR. ILETA.EQ.NLETA) THEN!element side?
                  IF (KEEP.EQ.1) THEN !check same point
                     IF (DABS(X - XLOCS(ILOC1)).LT.TOL .AND.
     ;                   DABS(Z - ZLOCS(ILOC1)).LT.TOL) LFS = .TRUE.
                  ELSEIF (KEEP.EQ.2) THEN !check along line
                     V2(1) = X - XLOCS(ILOC1) !line cant be parallel
                     V2(2) = 0.D0
                     V2(3) = Z - ZLOCS(ILOC1)
                     DET = V1(1)*V2(3) - V1(3)*V2(1)
                     IF (DABS(DET).LT.TOL) LFS = .TRUE.
                  ENDIF !end check on keep
               ENDIF !end check on element side
!
!............. double check free surface node not already claimed
               IF (LFS) THEN
                  DO 15 IFS=1,JFS 
                     IDOF = LM(1,IAE,IELEM) 
                     IF (IDOF.GT.0) THEN
                        IF (IDOF_FSW(JFS).EQ.IDOF) GOTO 150
                     ENDIF
   15             CONTINUE
                  JFS = JFS + 1
                  IDOF_FSW(JFS) = LM(1,IAE,IELEM) 
                  XFS(JFS) = X 
  150             CONTINUE
               ENDIF
   13       CONTINUE !loop on xi anchor nodes
   12    CONTINUE !loop on eta anchor nodes
   11 CONTINUE !loop on elements
      NFS = JFS 
      DEALLOCATE(SF) 
!
!.... sort 
      CALL SHELLR8I2(NFS,XFS,IDOF_FSW) 
!
!.... unpack 
      DO 20 IFS=1,NFS
         DO 21 I=1,NDIM
            IDOF_FS(I,IFS) = IDOF_FSW(IFS) + I - 1 !triplets
   21    CONTINUE
   20 CONTINUE
  100 CONTINUE 
      DEALLOCATE(IDOF_FSW) 
      DEALLOCATE(XFS)
      RETURN
      END 
