      SUBROUTINE LOCREC2D(MDIM,MEN,MGNOD, NNPG,NREC,NLXI,NLETA,
     ;                    LVERB, NELEM,NGNOD,NDIM, CNNPG,LFSURF,IENG,LM,
     ;                    XIPTS,ETAPTS,XREC,ZREC,XLOCS,ZLOCS, 
     ;                    MRDOF,IERR) 
!
!     This pairs the receivers with a degree of freedom based on 
!     collocation to the nearest degree of freedom and the receivers
!     (x,z) position
!
!.... variable declarations
      IMPLICIT REAL*8 (A-H,O-Z) 
      CHARACTER(2), INTENT(IN) :: CNNPG(NNPG) 
      REAL*8, INTENT(IN) :: XREC(NREC), ZREC(NREC), 
     ;        XLOCS(NNPG),ZLOCS(NNPG), XIPTS(NLXI),ETAPTS(NLETA)
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*),
     ;           MDIM,MEN,MGNOD, NNPG,NREC,NLXI,NLETA,  
     ;           NELEM,NGNOD,NDIM  
      LOGICAL*4, INTENT(IN) :: LFSURF(NREC), LVERB 
      INTEGER*4, INTENT(OUT) :: MRDOF(MDIM,*), IERR
!.... local variables
      REAL*8, ALLOCATABLE :: DMIN(:), SF(:)
      REAL*8 V1(3),V2(3), BIG,X,Z,XI,ETA,DIST,DET, TOL
      LOGICAL*4 LCFREE, LSFREE, LFS
      PARAMETER(LNULL =-5)
      PARAMETER(TOL = 2.22D-7) 
!
!----------------------------------------------------------------------!
!
!.... null out the DOFs and intiialize distance to something big
      IERR = 0
      BIG = HUGE(1.D0)
      LCFREE = .FALSE.
      ALLOCATE(DMIN(NREC))
      DO 1 IREC=1,NREC
         DMIN(IREC) = BIG
         DO 2 I=1,NDIM
            MRDOF(I,IREC) = LNULL 
    2    CONTINUE
         IF (LFSURF(IREC)) LCFREE = .TRUE.
    1 CONTINUE 
      ALLOCATE(SF(NGNOD))
!
!.... loop on geometry
      DO 11 IELEM=1,NELEM
         KEEP = 0
         ILOC1 = 0
         ILOC2 = 0
         LSFREE = .FALSE.
         IF (LCFREE) THEN !need to check if element is on free surface 
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
  102             CONTINUE !loop on anchor nodes
                  GOTO 40
               ENDIF !end check if node is on free surface
  101       CONTINUE !loop on anchor nodes
   40       CONTINUE !break ahead
         ENDIF
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
               IF (LCFREE) THEN !check if node is at free surface
                  IF (ILXI .EQ.1 .OR. ILXI .EQ.NLXI .OR. 
     ;                ILETA.EQ.1 .OR. ILETA.EQ.NLETA) THEN!element side?
                     IF (KEEP.EQ.1) THEN !check same point
                        IF (DABS(X - XLOCS(ILOC1)).LT.TOL .AND. 
     ;                      DABS(Z - ZLOCS(ILOC1)).LT.TOL) LFS = .TRUE. 
                     ELSEIF (KEEP.EQ.2) THEN !check along line
                        V2(1) = X - XLOCS(ILOC1) !line cant be parallel
                        V2(2) = 0.D0 
                        V2(3) = Z - ZLOCS(ILOC1) 
                        DET = V1(1)*V2(3) - V1(3)*V2(1)
                        IF (DABS(DET).LT.TOL) LFS = .TRUE. 
                     ENDIF !end check on keep
                  ENDIF !end check on element side
               ENDIF !end check on whether or not to calculate f.s.
!
!............. check we have all DOFs at nodal point
               IF (MINVAL(LM(1:NDIM,IAE,IELEM)).GT.0) THEN 
                  DO 15 IREC=1,NREC
                     DIST = (X - XREC(IREC))**2 + (Z - ZREC(IREC))**2
                     IF (LFSURF(IREC)) THEN  !take closest pt at FS
                        IF (DIST.LT.DMIN(IREC).AND.LFS) THEN
                           DO 16 I=1,NDIM
                              MRDOF(I,IREC) = LM(I,IAE,IELEM)
   16                      CONTINUE
                           DMIN(IREC) = DIST
                        ENDIF
                     ELSE !just take closest point
                        IF (DIST.LT.DMIN(IREC)) THEN
                           DO 17 I=1,NDIM
                              MRDOF(I,IREC) = LM(I,IAE,IELEM)
   17                      CONTINUE
                           DMIN(IREC) = DIST
                        ENDIF
                     ENDIF
   15             CONTINUE !loop on receivers
               ENDIF !end check on point being a complete DOF
   13       CONTINUE !loop on xi points
   12    CONTINUE !loop on eta points
   11 CONTINUE !loop on elements
!
!.... quality control, check if receivers were initialized
      DO 21 IREC=1,NREC
         IF (MINVAL(MRDOF(1:NDIM,IREC)).EQ.LNULL) THEN
            WRITE(*,*) 'locrec2d: Error locating receiver',IREC
            IERR = 2
         ENDIF
         IF (LVERB) THEN 
            IF (IREC.EQ.1) WRITE(*,*)
            IF (IREC.EQ.1) WRITE(*,905)
            WRITE(*,906) IREC,DSQRT(DMIN(IREC)) 
            IF (IREC.EQ.NREC) WRITE(*,*) 
         ENDIF 
   21 CONTINUE !loop on receivers
  100 CONTINUE !break ahead for error
      DEALLOCATE(SF)
      DEALLOCATE(DMIN) 
  905 FORMAT(' locrec: Distance of recievers to collocation points:')
  906 FORMAT('         Receiver:',I4,' distance ',F12.2, '(m)')
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE REC_SPACE1D(NREC,XREC, DX,IERR)
!
!     Calculates the receiver spacing along a line.  It is assumed
!     the receivers are in ascending or descending order
!
!     INPUT      MEANING
!     -----      -------
!     NREC       number of receivers
!     XREC       receiver x locations (in ascending/descending order)
!  
!     OUTPUT     MEANING
!     ------     -------
!     DX         normalized distance between receivers 
!     IERR       error flag
!
!.... variable declarations
      implicit none
      REAL*8, INTENT(IN) :: XREC(NREC)
      INTEGER*4, INTENT(IN) :: NREC
      REAL*8, INTENT(OUT) :: DX(NREC-1)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8 DXMIN, XDIR
      INTEGER*4 IREC
!
!----------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      DX(1:NREC-1) = 0.D0
      IF (NREC.EQ.1) THEN
         WRITE(*,*) 'offset1d: Error only one receiver!'
         IERR = 1
         RETURN
      ENDIF
!
!.... calculate offsets
      XDIR = XREC(2) - XREC(1) 
      DXMIN = HUGE(1.D0) 
      DO 1 IREC=1,NREC-1
         DX(IREC) = DABS(XREC(IREC+1) - XREC(IREC))
         DXMIN = DMIN1(DXMIN,DX(IREC))
         IF (DX(IREC).LT.1.D-10) THEN
            WRITE(*,*) 'offset1d: Error duplicate receiver locations!'
            DX(1:NREC-1) = 0.D0
            IERR = 1
            RETURN
         ENDIF
         IF (IREC.GT.1) THEN
            IF (XREC(IREC).GT.XREC(IREC-1) .AND. XDIR.LT.0.D0 .OR.
     ;          XREC(IREC).LT.XREC(IREC-1) .AND. XDIR.GT.0.D0) THEN
               WRITE(*,*) 'offset1d: Error xrec not in order!'
               DX(1:NREC-1) = 0.D0
               IERR = 1
               RETURN 
            ENDIF
          ENDIF
    1 CONTINUE
!
!.... normalize
      DO 2 IREC=1,NREC-1
         DX(IREC) = DX(IREC)/DXMIN
    2 CONTINUE
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LAPLACE_1D(LDC,NREC,DX, COVD,IERR)
!
!     This applies a 1D laplacian differencing scheme to each 
!     component of motion.  So that we may respect varying offsets 
!     we use a finite element scheme which in the case of a constnat 
!     dx reduces to a second order differencing matrix [-1 2 -1].  
!     The matrix can then multiply the Jacobian or residuals
!
!     INPUT      MEANING
!     -----      ------- 
!     DX         offset (likely normalized) 
!     LDC        leading dimension
!     NREC       number of receivers
!
!     OUTPUT     MEANING
!     ------     -------
!     COVD       data covariance matrix (1d FEM lapalace discretization)
!     IERR       error flag
!
!.... variable declarations
      REAL*8, INTENT(IN) :: DX(NREC-1)
      INTEGER*4, INTENT(IN) :: LDC,NREC
      REAL*8, INTENT(OUT) :: COVD(LDC,ldc) 
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      INTEGER*4, ALLOCATABLE :: LM(:,:,:), ID(:,:), IEN(:,:)
      REAL*8 K(2,2)
      INTEGER*4 NNPG,NELEM,NDOF, IELEM,INPG,IA,JA,IDOF,JDOF,I,J
      REAL*8, PARAMETER :: TWO  = 2.D0 
      REAL*8, PARAMETER :: HALF = 0.5D0
      INTEGER*4, PARAMETER :: NDIM = 3  
      INTEGER*4, PARAMETER :: NGNOD = 2  
!
!----------------------------------------------------------------------!
!
!.... generate a graph
      NELEM = NREC - 1  
      ALLOCATE(IEN(NGNOD,NELEM))
      INPG = 0  
      DO 1 IELEM=1,NELEM
         DO 2 IA=1,NGNOD
            INPG = INPG + 1  
            IEN(IA,IELEM) = INPG 
    2    CONTINUE
         INPG = INPG - 1
    1 CONTINUE
      NNPG = INPG + 1
      IF (NNPG.NE.NREC) THEN
         WRITE(*,*) 'laplace_1d: nnpg initialization error'
         IERR = 1
         RETURN
      ENDIF
      NDOF = NDIM*NNPG
      ALLOCATE(ID(NDIM,NNPG))
      IDOF = 0
      DO 3 INPG=1,NNPG
         DO 4 I=1,NDIM
            IDOF = IDOF + 1
            ID(I,INPG) = IDOF
    4    CONTINUE
    3 CONTINUE
      IF (IDOF.NE.NDOF) THEN
         WRITE(*,*) 'laplace_1d: idof initialization error'
         IERR = 1
         RETURN
      ENDIF
      ALLOCATE(LM(3,NGNOD,NELEM))
      DO 5 IELEM=1,NELEM
         DO 6 IA=1,NGNOD
            DO 7 I=1,NDIM
               LM(I,IA,IELEM) = ID(I,IEN(IA,IELEM))
    7       CONTINUE
    6    CONTINUE
    5 CONTINUE
      DEALLOCATE(ID)
      DEALLOCATE(IEN)
!
!.... loop on elements
      K(1,1) = HALF
      K(1,2) =-HALF
      K(2,1) =-HALF
      K(2,2) = HALF
      COVD(1:LDC,1:2*NREC*NDIM) = 0.D0
      DO 8 IELEM=1,NELEM
         DO 9 IA=1,NGNOD
            DO 10 I=1,NDIM
               IDOF = LM(I,IA,IELEM)
               DO 11 JA=1,NGNOD
                  DO 12 J=1,NDIM
                     JDOF = LM(J,JA,IELEM)
                     IF (I.EQ.J) THEN
                        COVD(IDOF,JDOF) = COVD(IDOF,JDOF) 
     ;                                  + TWO/DX(IELEM)*K(IA,JA)
                        COVD(NDOF+IDOF,NDOF+JDOF)  
     ;                = COVD(NDOF+IDOF,NDOF+JDOF)+TWO/DX(IELEM)*K(IA,JA)
                     ENDIF
   12             CONTINUE
   11          CONTINUE
   10       CONTINUE
    9    CONTINUE
    8 CONTINUE
      DEALLOCATE(LM)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE LAPLACE_1DV2(LDC,MDIM,NREC,LDWGHT, DX, OBS, COVD,IERR)
!
!     Builds the data roughening matrix  
!          [ Re{L_w}    0       0         0        0       0     ]
!          [    0    Re{L_v}    0         0        0       0     ]
!          [    0       0    Re{L_u}      0        0       0     ] 
!          [    0       0       0      Im{L_w}     0       0     ]
!          [    0       0       0         0     Im{L_v}    0     ]
!          [    0       0       0         0        0     Im{L_u} ]
! 
!     where each submatrix is ordered so that only receivers 
!     with observations are listed first, then, back filled with
!     0's 
!
!     INPUT     MEANING
!     -----     ------- 
!     DX        normalized distances between receivers in x 
!     LDC       leading dimension
!     LDWGHT    True -> weight by receiver distance
!     MDIM      leading dimension
!     NREC      number of receviers
!     OBS       observations for frequency source pair 
!
!     OUTPUT    MEANING 
!     ------    -------
!     COVD      data covariance matrix
!     IERR      error flag
!
!.... variable declarations
      IMPLICIT NONE 
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,*) 
      REAL*8, INTENT(IN) :: DX(NREC-1)
      INTEGER*4, INTENT(IN) :: LDC,MDIM,NREC 
      LOGICAL*4, INTENT(IN) :: LDWGHT
      REAL*8, INTENT(OUT) :: COVD(LDC,*)
      INTEGER*4, INTENT(OUT) :: IERR 
!.... local variables
      REAL*8, ALLOCATABLE :: DXLOC(:) 
      INTEGER*4, ALLOCATABLE :: LM(:,:), IEN(:,:), ID(:) 
      LOGICAL*4, ALLOCATABLE :: LACTIVE(:)
      REAL*8 K(2,2), XSUM 
      INTEGER*4 NDOF,NELEM,NNPG, IREC,JREC,INPG,IELEM,JELEM, IA,JA,
     ;          IDOF,JDOF, IROW1,IROW2,JCOL1,JCOL2, I 
      REAL*8, PARAMETER :: TWO = 2.D0 
      REAL*8, PARAMETER :: ONE = 1.D0 
      REAL*8, PARAMETER :: HALF = 0.5D0
      INTEGER*4, PARAMETER :: NDIM = 3  
      INTEGER*4, PARAMETER :: NGNOD = 2  
!
!----------------------------------------------------------------------!
!
!.... initialize and check for early return
      IERR = 0  
      COVD(1:LDC,1:2*NREC*NDIM) = 0.D0 
      IF (MAXVAL(CABS(OBS(1:NDIM,1:NREC))) == 0.0) THEN 
         WRITE(*,*) 'laplace_1dv2: This is a dead frequency!'
         RETURN
      ENDIF
      K(1,1) = HALF 
      K(1,2) =-HALF
      K(2,1) =-HALF
      K(2,2) = HALF 
      IDOF = 0
!
!.... loop on components 
      DO 100 I=NDIM,1,-1
         ALLOCATE(LACTIVE(NREC))
         NELEM =-1
         DO 1 IREC=1,NREC
            LACTIVE(IREC) = .FALSE.
            IF (CABS(OBS(I,IREC)) > 0.0) THEN
               NELEM = NELEM + 1
               LACTIVE(IREC) = .TRUE.
            ENDIF
    1    CONTINUE
         IF (NELEM ==-1) GOTO 105 !dead component 
!
!....... generate local dx vector
         ALLOCATE(DXLOC(NELEM))
         DXLOC(1:NELEM) = 0.D0
         IF (LDWGHT) THEN
            IELEM = 0
            DO 2 IREC=1,NREC !IELEM=1,NELEM 
               IF (LACTIVE(IREC) .AND. IELEM + 1 <= NELEM) THEN
                  IELEM = IELEM + 1
                  XSUM = 0.D0
                  JELEM = IREC - 1
                  DO 3 JREC=IREC,NREC-1
                     JELEM = JELEM + 1
                     XSUM = XSUM + DX(JELEM)
                     IF (LACTIVE(JREC+1)) GOTO 30
    3             CONTINUE
   30             CONTINUE
                  DXLOC(IELEM) = XSUM
                  IF (XSUM == 0.D0) THEN
                     WRITE(*,*) 'laplace_1d: xsum error',IREC,JREC
                     IERR = 1
                     RETURN
                  ENDIF
               ENDIF
    2       CONTINUE
         ELSE
            DXLOC(1:NELEM) = ONE
         ENDIF
         DEALLOCATE(LACTIVE)
         IF (MINVAL(DXLOC) == 0.D0) THEN
            WRITE(*,*) 'laplace_1d: dxloc init error'
            IERR = 1
            RETURN
         ENDIF
!
!....... generate a graph 
         ALLOCATE(IEN(NGNOD,NELEM))
         INPG = 0
         DO 4 IELEM=1,NELEM
            DO 5 IA=1,NGNOD
               INPG = INPG + 1
               IEN(IA,IELEM) = INPG
    5       CONTINUE
            INPG = INPG - 1
    4    CONTINUE
         NNPG = INPG + 1
         IF (NNPG /= NELEM + 1) THEN
            WRITE(*,*) 'laplace_1d: nnpg initialization error',NNPG
            IERR = 1
            RETURN
         ENDIF
         NDOF = NNPG
         ALLOCATE(ID(NNPG))
         INPG = 0
         DO 6 IREC=1,NREC !INPG=1,NNPG
            IF (CABS(OBS(I,IREC)) > 0.0) THEN
               INPG = INPG + 1
               IDOF = IDOF + 1
               ID(INPG) = IDOF
            ENDIF
    6    CONTINUE
         ALLOCATE(LM(NGNOD,NELEM))
         DO 7 IELEM=1,NELEM
            DO 8 IA=1,NGNOD
               LM(IA,IELEM) = ID(IEN(IA,IELEM))
    8       CONTINUE
    7    CONTINUE
         DEALLOCATE(ID)
         DEALLOCATE(IEN)
!
!....... loop on elements
         DO 21 IELEM=1,NELEM
            DO 22 IA=1,NGNOD
               IDOF = LM(IA,IELEM)
               IF (IDOF == 0) GOTO 220
               DO 23 JA=1,NGNOD
                  JDOF = LM(JA,IELEM)
                  IF (JDOF == 0) GOTO 230
                  IROW1 = IDOF
                  IROW2 = NDIM*NREC + IDOF
                  JCOL1 = JDOF
                  JCOL2 = NDIM*NREC + JDOF
                  COVD(IROW1,JCOL1) = COVD(IROW1,JCOL1) 
     ;                              + TWO/DXLOC(IELEM)*K(IA,JA)
                  COVD(IROW2,JCOL2) = COVD(IROW2,JCOL2) 
     ;                              + TWO/DXLOC(IELEM)*K(IA,JA)
  230             CONTINUE !not a dof
   23          CONTINUE
  220          CONTINUE !not a dof
   22       CONTINUE
   21    CONTINUE
  105    CONTINUE !dead component
         IF (ALLOCATED(LACTIVE)) DEALLOCATE(LACTIVE)
         IF (ALLOCATED(DXLOC))   DEALLOCATE(DXLOC)
         IF (ALLOCATED(LM))      DEALLOCATE(LM)
  100 CONTINUE !loop on components
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE PERM_OBS(MDIM,NDIM,NREC, OBS,  
     ;                    IPERM_REC,IERR)
!
!     Generates a permutation list for the observations.  We want the
!     data stored stored active observations for the vertical component
!     then null vertical receivers, then active receivers for 
!     the v component then null receivers, and finally the active
!     u receivers then null receivers 
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,*)
      INTEGER*4, INTENT(IN) :: MDIM,NDIM,NREC 
      INTEGER*4, INTENT(OUT) :: IPERM_REC(NDIM*NREC), IERR 
!
!----------------------------------------------------------------------!
!
!.... set the observations
      IPERM_REC(1:NDIM*NREC) = 0
      IOBS = 0
      DO 1 I=NDIM,1,-1
         DO 2 IREC=1,NREC 
            IF (CABS(OBS(I,IREC)) > 0.0) THEN 
               INDX = (IREC - 1)*NDIM + I
               IOBS = IOBS + 1
               IPERM_REC(INDX) = IOBS 
            ENDIF
    2    CONTINUE
         DO 3 IREC=1,NREC
            IF (CABS(OBS(I,IREC)) == 0.0) THEN 
               INDX = (IREC - 1)*NDIM + I
               IOBS = IOBS + 1
               IPERM_REC(INDX) = IOBS 
            ENDIF
    3    CONTINUE 
    1 CONTINUE
      IF (MINVAL(IPERM_REC) == 0) THEN 
         WRITE(*,*) 'perm_obs: Error setting permutation vector!'
         IERR = 1
      ELSE 
         IERR = 0
      ENDIF
      RETURN
      END  
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE LAPLACE_DIST1D(ICTXT,NPROW,NPCOL, LDC,NREC,MCL,NCL, 
     ;                          DESCC, DX, COVD,IERR)
!
!     This applies a 1D laplacian differencing scheme to each component
!     of motion.  So that we may respect varying offsets we use a 
!     finite element scheme which in the case of a constnat dx reduces 
!     to a second order differencing matrix [-1 2 -1].  The matrix can 
!     then multiply the Jacobian or residuals
!
!     INPUT      MEANING
!     -----      ------- 
!     DX         offset (likely normalized) 
!     LDC        leading dimension
!     MCL        number of local rows in distributed data covariance
!     matrix
!     NCLO       number of local columns in distributed data covariance
!     matrix
!     NREC       number of receivers
!
!     OUTPUT     MEANING
!     ------     -------
!     COVD       data covariance matrix (1d FEM lapalace discretization)
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DX(NREC-1)
      INTEGER*4, INTENT(IN) :: DESCC(9), LDC,MCL,NCL,NREC, 
     ;                         NPROW,NPCOL,ICTXT
      REAL*8, INTENT(OUT) :: COVD(LDC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      INTEGER*4, ALLOCATABLE :: LM(:,:,:), ID(:,:), IEN(:,:)
      REAL*8 K(2,2)
      INTEGER*4 NNPG,NELEM,NDOF, IELEM,INPG,IA,JA,IDOF,JDOF,I,J
      INTEGER*4 MYROW,MYCOL, IAROW,IACOL, IIA,JJA, IROW,JCOL
      REAL*8, PARAMETER :: TWO  = 2.D0
      REAL*8, PARAMETER :: HALF = 0.5D0
      INTEGER*4, PARAMETER :: NDIM = 3
      INTEGER*4, PARAMETER :: NGNOD = 2
!
!----------------------------------------------------------------------!
!
!.... recall BLACS information
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!.... generate a graph
      NELEM = NREC - 1
      ALLOCATE(IEN(NGNOD,NELEM))
      INPG = 0
      DO 1 IELEM=1,NELEM
         DO 2 IA=1,NGNOD
            INPG = INPG + 1
            IEN(IA,IELEM) = INPG
    2    CONTINUE
         INPG = INPG - 1
    1 CONTINUE
      NNPG = INPG + 1
      IF (NNPG.NE.NREC) THEN
         WRITE(*,*) 'laplace_dist1d: nnpg initialization error'
         IERR = 1
         RETURN
      ENDIF
      NDOF = NDIM*NNPG
      ALLOCATE(ID(NDIM,NNPG))
      IDOF = 0
      DO 3 INPG=1,NNPG
         DO 4 I=1,NDIM
            IDOF = IDOF + 1
            ID(I,INPG) = IDOF
    4    CONTINUE
    3 CONTINUE
      IF (IDOF.NE.NDOF) THEN
         WRITE(*,*) 'laplace_dist1d: idof initialization error'
         IERR = 1
         RETURN
      ENDIF
      ALLOCATE(LM(3,NGNOD,NELEM))
      DO 5 IELEM=1,NELEM
         DO 6 IA=1,NGNOD
            DO 7 I=1,NDIM
               LM(I,IA,IELEM) = ID(I,IEN(IA,IELEM))
    7       CONTINUE
    6    CONTINUE
    5 CONTINUE
      DEALLOCATE(ID)
      DEALLOCATE(IEN)
!
!.... loop on elements
      K(1,1) = HALF
      K(1,2) =-HALF
      K(2,1) =-HALF
      K(2,2) = HALF
      COVD(1:MCL,1:NCL) = 0.D0
      DO 8 IELEM=1,NELEM
         DO 9 IA=1,NGNOD
            DO 10 I=1,NDIM
               IDOF = LM(I,IA,IELEM)
               DO 11 JA=1,NGNOD
                  DO 12 J=1,NDIM
                     JDOF = LM(J,JA,IELEM)
                     IF (I == J) THEN
                        IROW = IDOF
                        JCOL = JDOF
                        CALL INFOG2L(IROW,JCOL,DESCC, NPROW,NPCOL,
     ;                               MYROW,MYCOL, IIA,JJA, IAROW,IACOL)
                        IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN 
                           COVD(IIA,JJA) = COVD(IIA,JJA)
     ;                                   + TWO/DX(IELEM)*K(IA,JA)
                        ENDIF
                        IROW = NDOF + IDOF
                        JCOL = NDOF + JDOF
                        CALL INFOG2L(IROW,JCOL,DESCC, NPROW,NPCOL,
     ;                               MYROW,MYCOL, IIA,JJA, IAROW,IACOL)
                        IF (MYROW == IAROW .AND. MYCOL == IACOL) THEN 
                           COVD(IIA,JJA) = COVD(IIA,JJA)
     ;                                   + TWO/DX(IELEM)*K(IA,JA)
                        ENDIF
                     ENDIF
   12             CONTINUE
   11          CONTINUE
   10       CONTINUE
    9    CONTINUE
    8 CONTINUE
      DEALLOCATE(LM)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE LAPLACE_DIST1DV2(ICTXT,NPROW,NPCOL,                    
     ;                        LDC,MDIM,NREC,MCL,NCL,DESCC,  LDWGHT, DX, 
     ;                        OBS,COVD,IERR)
!
!     This applies a 1D laplacian differencing scheme to each component of motion.  
!     So that we may respect varying offsets we use a finite element scheme which in 
!     the case of a constnat dx reduces to a second order differencing matrix [-1 2 -1].  
!     The matrix can then multiply the Jacobian or residuals
!
!     INPUT      MEANING
!     -----      ------- 
!     DESCC      covariance matrix descriptor
!     DX         offset (likely normalized) 
!     LDC        leading dimension
!     LDWGHT     true use distance weighting
!     MCL        number of local rows in distributed data covariance matrix
!     MDIM       leading dimension
!     NCL        number of local columns in distributed data covariance matrix
!     NREC       number of receivers
!     OBS        observations for frequency source pair
!
!     OUTPUT     MEANING
!     ------     -------
!     COVD       data covariance matrix (1d FEM lapalace discretization)
!     IERR       error flag
!
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: OBS(MDIM,*)
      REAL*8, INTENT(IN) :: DX(NREC-1)
      INTEGER*4, INTENT(IN) :: DESCC(9), LDC,MDIM,MCL,NCL,NREC, 
     ;                         NPROW,NPCOL,ICTXT
      LOGICAL*4, INTENT(IN) :: LDWGHT
      REAL*8, INTENT(OUT) :: COVD(LDC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*8, ALLOCATABLE :: DXLOC(:)
      INTEGER*4, ALLOCATABLE :: LM(:,:), ID(:), IEN(:,:)
      LOGICAL*4, ALLOCATABLE :: LACTIVE(:)
      REAL*8 K(2,2), XSUM
      INTEGER*4 NNPG,NELEM,NDOF, IELEM,JELEM,INPG,IA,JA,IDOF,JDOF,I,
     ;          IREC,JREC
      INTEGER*4 MYROW,MYCOL, IAROW,IACOL, IIA,JJA, IROW,JCOL
      REAL*8, PARAMETER :: TWO  = 2.D0
      REAL*8, PARAMETER :: ONE  = 1.D0
      REAL*8, PARAMETER :: HALF = 0.5D0
      INTEGER*4, PARAMETER :: NDIM = 3
      INTEGER*4, PARAMETER :: NGNOD = 2
!
!----------------------------------------------------------------------------------------!
!
!.... recall BLACS information
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!
!.... initialize and check for early return
      IERR = 0
      COVD(1:MCL,1:NCL) = 0.D0
      IF (MAXVAL(CABS(OBS(1:NDIM,1:NREC))) == 0.0) THEN
         WRITE(*,*) 'laplace_1dv2: This is a dead frequency!'
         RETURN
      ENDIF
      K(1,1) = HALF
      K(1,2) =-HALF
      K(2,1) =-HALF
      K(2,2) = HALF
      IDOF = 0
!
!.... loop on components 
      DO 100 I=NDIM,1,-1
         ALLOCATE(LACTIVE(NREC))
         NELEM =-1
         DO 1 IREC=1,NREC
            LACTIVE(IREC) = .FALSE.
            IF (CABS(OBS(I,IREC)) > 0.0) THEN
               NELEM = NELEM + 1
               LACTIVE(IREC) = .TRUE.
            ENDIF
    1    CONTINUE
         IF (NELEM ==-1) GOTO 105 !dead component 
!
!....... generate local dx vector
         ALLOCATE(DXLOC(NELEM))
         DXLOC(1:NELEM) = 0.D0
         IF (LDWGHT) THEN
            IELEM = 0
            DO 2 IREC=1,NREC !IELEM=1,NELEM 
               IF (LACTIVE(IREC) .AND. IELEM + 1 <= NELEM) THEN
                  IELEM = IELEM + 1
                  XSUM = 0.D0
                  JELEM = IREC - 1
                  DO 3 JREC=IREC,NREC-1
                     JELEM = JELEM + 1
                     XSUM = XSUM + DX(JELEM)
                     IF (LACTIVE(JREC+1)) GOTO 30
    3             CONTINUE
   30             CONTINUE
                  DXLOC(IELEM) = XSUM
                  IF (XSUM == 0.D0) THEN
                     WRITE(*,*) 'laplace_dist1d: xsum error',IREC,JREC
                     IERR = 1
                     RETURN
                  ENDIF
               ENDIF
    2       CONTINUE
         ELSE
            DXLOC(1:NELEM) = ONE
         ENDIF
         DEALLOCATE(LACTIVE)
         IF (MINVAL(DXLOC) == 0.D0) THEN
            WRITE(*,*) 'laplace_dist1d: dxloc init error'
            IERR = 1
            RETURN
         ENDIF
!
!....... generate a graph 
         ALLOCATE(IEN(NGNOD,NELEM))
         INPG = 0
         DO 4 IELEM=1,NELEM
            DO 5 IA=1,NGNOD
               INPG = INPG + 1
               IEN(IA,IELEM) = INPG
    5       CONTINUE
            INPG = INPG - 1
    4    CONTINUE
         NNPG = INPG + 1
         IF (NNPG /= NELEM + 1) THEN
            WRITE(*,*) 'laplace_dist1d: nnpg initialization error',NNPG
            IERR = 1
            RETURN
         ENDIF
         NDOF = NNPG
         ALLOCATE(ID(NNPG))
         INPG = 0
         DO 6 IREC=1,NREC !INPG=1,NNPG
            IF (CABS(OBS(I,IREC)) > 0.0) THEN
               INPG = INPG + 1
               IDOF = IDOF + 1
               ID(INPG) = IDOF
            ENDIF
    6    CONTINUE
         ALLOCATE(LM(NGNOD,NELEM))
         DO 7 IELEM=1,NELEM
            DO 8 IA=1,NGNOD
               LM(IA,IELEM) = ID(IEN(IA,IELEM))
    8       CONTINUE
    7    CONTINUE
         DEALLOCATE(ID)
         DEALLOCATE(IEN)
         DO 21 IELEM=1,NELEM
            DO 22 IA=1,NGNOD
               IDOF = LM(IA,IELEM)
               IF (IDOF == 0) GOTO 220
               DO 23 JA=1,NGNOD
                  JDOF = LM(JA,IELEM)
                  IF (JDOF == 0) GOTO 230
                  IROW = IDOF
                  JCOL = JDOF
                  CALL INFOG2L(IROW,JCOL,DESCC, NPROW,NPCOL,MYROW,MYCOL,
     ;                         IIA,JJA, IAROW,IACOL)
                  IF (MYROW == IAROW .AND. MYCOL == IACOL)  
     ;            COVD(IIA,JJA) = COVD(IIA,JJA) + TWO/DX(IELEM)*K(IA,JA)
                  IROW = NDIM*NREC + IDOF
                  JCOL = NDIM*NREC + JDOF
                  CALL INFOG2L(IROW,JCOL,DESCC, NPROW,NPCOL,MYROW,MYCOL,
     ;                         IIA,JJA, IAROW,IACOL)
                  IF (MYROW == IAROW .AND. MYCOL == IACOL)  
     ;            COVD(IIA,JJA) = COVD(IIA,JJA) + TWO/DX(IELEM)*K(IA,JA)
  230             CONTINUE !not a dof; break ahead 
   23          CONTINUE !loop on anchor nodes
  220          CONTINUE
   22       CONTINUE
   21    CONTINUE
  105    CONTINUE !dead component
         IF (ALLOCATED(LACTIVE)) DEALLOCATE(LACTIVE)
         IF (ALLOCATED(DXLOC))   DEALLOCATE(DXLOC)
         IF (ALLOCATED(LM))      DEALLOCATE(LM)
  100 CONTINUE !loop on components
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE EXRESP25D(MDIM,MREC,NDOF, NSRC,NREC,NDIM, FREQ,
     ;                     PY, YREC, SOL, MRDOF, 
     ;                     RECV,EST)
! 
!     Extract the response at all receivers for all sources for this 
!     frequency.  Additionally we apply the phase shifter ala 
!     Kennett e^{i(omega t - kx x - ky y - kz z)} is wave propagating
!     forward in time, with local conventions +x right, +z up, then
!     we'd like our y to be +90 degrees from x, so we phase shift 
!     at e^{+i omega py y} e^{i(omega t - kx x - kz z)} 
!
!     INPUT      MEANING
!     -----      ------- 
!     FREQ       current frequency in Hz
!     MDIM       max spatial dimension; leading dimension
!     MRDOF      holds the receivers DOF numbers
!     MREC       max number of recievers; leading dimension
!     NDIM       number of components in solution 
!     NDOF       number of degrees of freedom
!     NREC       number of receivers to extract at
!     NSRC       number of sources to extract
!     PY         apparent slownesses at sources
!     RECV       receiver reseponse function
!     SOL        solution of AX = B 
!     YREC       y position of recievers
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     EST        estimated response from wavefield at receiver location
!
!.... varibale declarations
      COMPLEX*8, INTENT(IN) :: RECV(MDIM,*), SOL(NDOF*NSRC)
      REAL*8, INTENT(IN) :: YREC(NREC), PY(NSRC),FREQ
      INTEGER*4, INTENT(IN) :: MRDOF(MDIM,*), MDIM,MREC,NDOF,  
     ;                         NSRC,NREC,NDIM
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MREC,*) 
      COMPLEX*8 CFACT
      REAL*8 TWOPI, OMEGA
      PARAMETER(TWOPI = 6.2831853071795862D0)
! 
!----------------------------------------------------------------------!
!
      OMEGA = TWOPI*FREQ
      DO 1 ISRC=1,NSRC
         LOFF = (ISRC - 1)*NDOF 
         DO 2 IREC=1,NREC
            CFACT = CEXP(CMPLX(0.D0,+OMEGA*PY(ISRC)*YREC(IREC)))
            DO 3 I=1,NDIM
               IF (MRDOF(I,IREC).NE.0) THEN
                  ILOC = LOFF + MRDOF(I,IREC)
                  EST(I,IREC,ISRC) = RECV(I,IREC)*SOL(ILOC)*CFACT 
               ELSE
                  EST(I,IREC,ISRC) = CMPLX(0.0,0.0) 
               ENDIF
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE EXRESP_MP_2(MDIM,MREC,NDOF, NSRC,NREC,NDIM, FREQ,
     ;                     PY, YREC, SOL, MRDOF,
     ;                     RECV,EST) 
!
!     Extracts the magnitude and phase in terms of [Mag,phase] 
      COMPLEX*8, INTENT(IN) :: RECV(MDIM,*), SOL(NDOF*NSRC)
      REAL*8, INTENT(IN) :: YREC(NREC), PY(NSRC),FREQ
      INTEGER*4, INTENT(IN) :: MRDOF(MDIM,*), MDIM,MREC,NDOF,
     ;                         NSRC,NREC,NDIM
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MREC,*)
      COMPLEX*8 CFACT, ESTT
      REAL*8 TWOPI, OMEGA
      REAL*4 SPHASE 
      PARAMETER(TWOPI = 6.2831853071795862D0)

      OMEGA = TWOPI*FREQ
      DO 1 ISRC=1,NSRC
         LOFF = (ISRC - 1)*NDOF
         DO 2 IREC=1,NREC
            CFACT = CEXP(CMPLX(0.D0,+OMEGA*PY(ISRC)*YREC(IREC)))
            DO 3 I=1,NDIM
               IF (MRDOF(I,IREC).NE.0) THEN
                  ILOC = LOFF + MRDOF(I,IREC)
                  ESTT = RECV(I,IREC)*SOL(ILOC)*CFACT
                  EST(I,IREC,ISRC) = CMPLX(CABS(ESTT),SPHASE(ESTT))
               ELSE
                  EST(I,IREC,ISRC) = CMPLX(0.0,0.0)
               ENDIF
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE

      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE EXRESP_MP(MDIM,MREC,NDOF, NSRC,NREC,NFS, FREQ,
     ;                     PY, YREC, USOL, CSIDE,
     ;                     MRDOF,IDOF_FS,RECV,
     ;                     EST, IERR)
!
!     Extracts the magntude and unwrapped phase response.
!     Additionally we apply the phase shifter ala 
!     Kennett e^{i(omega t - kx x - ky y - kz z)} is wave propagating
!     forward in time, with local conventions +x right, +z up, then
!     we'd like our y to be +90 degrees from x, so we phase shift 
!     at e^{+i py y} e^{i(omega t - kx x - kz z)} 
!
!     INPUT      MEANING
!     -----      ------- 
!     CSIDE      side source is approaching from
!     FREQ       current frequency in Hz
!     IDOF_FS    holds DOF numbers at fine nodes at free surface
!     MDIM       max spatial dimension; leading dimension
!     MRDOF      holds the receivers DOF numbers
!     MREC       max number of recievers; leading dimension
!     NDOF       number of degrees of freedom
!     NFS        number of fine nodes at free surface
!     NREC       number of receivers to extract at
!     NSRC       number of sources to extract
!     PY         apparent slownesses at sources
!     RECV       receiver reseponse function
!     SOL        solution of AX = B 
!     YREC       y position of recievers
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     EST        estimated response from wavefield at receiver location
!     IERR       error flag, couldn't find receiver
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(1), INTENT(IN) :: CSIDE(NSRC)
      COMPLEX*8, INTENT(IN) :: RECV(MDIM,*), USOL(NDOF*NSRC)
      REAL*8, INTENT(IN) :: YREC(NREC), PY(NSRC), FREQ 
      INTEGER*4, INTENT(IN) :: MRDOF(MDIM,*), IDOF_FS(MDIM,*), 
     ;                         MDIM,MREC,NDOF, NSRC,NREC,NFS
      COMPLEX*8, INTENT(OUT) :: EST(MDIM,MREC,*)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      REAL*4, ALLOCATABLE :: PHASE(:,:) 
      INTEGER*4, ALLOCATABLE :: LFOUND(:) 
      COMPLEX*8 U, V, W 
      REAL*4 XU, XV, XW, PU, PV, PW, PCU, PCV, PCW, XMU, XMV, XMW 
      INTEGER*4 ISRC, IREC, LOFF, IFS_BEG, IFS_END, IFS_DIR , IFS, I
      REAL*8 TWOPI, OMEGA 
      PARAMETER(TWOPI = 6.2831853071795862D0)
!
!----------------------------------------------------------------------!
!
!.... initialize
      OMEGA = TWOPI*FREQ
      ALLOCATE(PHASE(3,NFS))
      ALLOCATE(LFOUND(NREC)) 
      IERR = 0
!
!.... loop on sources
      DO 1 ISRC=1,NSRC
         LOFF = (ISRC - 1)*NDOF 
         IF (CSIDE(ISRC).EQ.'R') THEN
            IFS_BEG = NFS
            IFS_END = 1 
            IFS_DIR =-1
         ELSE
            IFS_BEG = 1
            IFS_END = NFS
            IFS_DIR = 1
         ENDIF
         DO 2 IFS=IFS_BEG,IFS_END,IFS_DIR
            U = USOL(LOFF + IDOF_FS(1,IFS))
            V = USOL(LOFF + IDOF_FS(2,IFS))
            W = USOL(LOFF + IDOF_FS(3,IFS))
            PHASE(1,IFS) = ATAN2(IMAG(U),REAL(U))
            PHASE(2,IFS) = ATAN2(IMAG(V),REAL(V))
            PHASE(3,IFS) = ATAN2(IMAG(W),REAL(W))
    2    CONTINUE 
!
!....... unwrap
         DO 3 I=1,3
            CALL SPUNWR(PHASE(I,1:NFS),NFS-1,1)
    3    CONTINUE
!
!....... extract and phase shift
         LFOUND(1:NREC) = 0
         DO 4 IREC=1,NREC
            DO 5 IFS=1,NFS
               IF (MRDOF(1,IREC).EQ.IDOF_FS(1,IFS)) THEN 
                  LFOUND(IREC) = 1
                  PCU = ATAN2(IMAG(RECV(1,IREC)),REAL(RECV(1,IREC)))
                  PCV = ATAN2(IMAG(RECV(2,IREC)),REAL(RECV(2,IREC)))
                  PCW = ATAN2(IMAG(RECV(3,IREC)),REAL(RECV(3,IREC)))
                  XMU = CABS(RECV(1,IREC))
                  XMV = CABS(RECV(2,IREC))
                  XMW = CABS(RECV(3,IREC))
!................ phase + phase corrections
                  PU = PHASE(1,IFS) + SNGL(OMEGA*PY(ISRC)*YREC(IREC)) 
     ;                              + PCU 
                  PV = PHASE(2,IFS) + SNGL(OMEGA*PY(ISRC)*YREC(IREC)) 
     ;                              + PCV
                  PW = PHASE(3,IFS) + SNGL(OMEGA*PY(ISRC)*YREC(IREC)) 
     ;                              + PCW
!................ amplitude * amplitude corrections
                  XU = CABS(USOL(LOFF + IDOF_FS(1,IFS)))*XMU
                  XV = CABS(USOL(LOFF + IDOF_FS(2,IFS)))*XMV
                  XW = CABS(USOL(LOFF + IDOF_FS(3,IFS)))*XMW
                  EST(1,IREC,ISRC) = CMPLX(XU,PU) 
                  EST(2,IREC,ISRC) = CMPLX(XV,PV)
                  EST(3,IREC,ISRC) = CMPLX(XW,PW)
               ENDIF
    5       CONTINUE 
    4    CONTINUE
!
!....... check we got everyone
         IF (MINVAL(LFOUND).EQ.0) THEN
            WRITE(*,*) 'exresp_mp: Error missed receivers',LFOUND
            IERR = 1
         ENDIF 
    1 CONTINUE !Loop on sources
      DEALLOCATE(LFOUND) 
      DEALLOCATE(PHASE)
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE REDEST4(MYID,MASTER,MYCOMM,NPGROUPS,IPGROUP,  
     ;                   MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC,JOB, EST)
! 
!     Reduce the frequency responses onto all or head process
!
!     INPUT      MEANING
!     -----      ------- 
!     EST        estimate response
!     JOB        = 1 (default) for reduction on master
!                = 2 for reduction on all processes
!     IPGROUP    process group ID, 0,1,2,...,npgroups-1
!     MASTER     head process ID
!     MDIM       leading dimension for est
!     MFREQ      leading dimension for est
!     MREC       leading dimension for est
!     MYCOMM     MPI communicator
!     MYID       process ID on communicator
!     NDIM       number of components
!     NFREQ      number of frequencies
!     NPGROUPS   number of process groups
!     NREC       number of receivers
!     NSC        number of sources
!
!     OUTPUT     MEANING
!     ------     ------- 
!     EST        estimate response on head node or all processes, se job
!
!.... variable declarations 
      INCLUDE 'mpif.h'
      COMPLEX*8 EST(MDIM,MFREQ,MREC,*) 
      INTEGER*4, INTENT(IN) :: MYID,MASTER,MYCOMM,NPGROUPS,IPGROUP, 
     ;           MDIM,MFREQ,MREC, NDIM,NFREQ,NREC,NSRC, JOB
!.... local variables
      ALLOCATABLE BUFF(:),QBUFF(:)
      COMPLEX*8 BUFF,QBUFF,ZERO
      PARAMETER(ZERO = CMPLX(0.0,0.0))
! 
!----------------------------------------------------------------------!
!
!.... nothing to do
      IF (NPGROUPS.EQ.1) RETURN
!
!.... set space
      NVARS = NDIM*NFREQ*NREC*NSRC
      ALLOCATE(BUFF(NVARS))
      BUFF(1:NVARS) = ZERO 
      IF (MYID.EQ.MASTER .OR. JOB.EQ.2) THEN
         ALLOCATE(QBUFF(NVARS))
      ELSE
         ALLOCATE(QBUFF(1))
      ENDIF
      DO 1 IFREQL=1,NFREQ
         IFREQ = (IFREQL - 1)*NPGROUPS + IPGROUP + 1 
         IF ( (IFREQL - 1)*NPGROUPS + 1.GT.NFREQ) GOTO 2500!out of freqs
         IF (IFREQ.GT.NFREQ) GOTO 1 !no frequency for this group
         ITER = (IFREQ - 1)*NDIM*NREC*NSRC
         DO 2 IREC=1,NREC
            DO 3 ISRC=1,NSRC
               DO 4 I=1,NDIM
                  ITER = ITER + 1 
                  BUFF(ITER) = EST(I,IFREQ,IREC,ISRC)
    4          CONTINUE
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE
 2500 CONTINUE !reduce
!
!.... reduce onto all
      IF (JOB.EQ.2) THEN
         CALL MPI_ALLREDUCE(BUFF,QBUFF,NVARS,MPI_COMPLEX,MPI_SUM,
     ;                      MYCOMM,MPIERR)
      ELSE
         CALL MPI_REDUCE(BUFF,QBUFF,NVARS,MPI_COMPLEX,MPI_SUM, 
     ;                   MASTER,MYCOMM,MPIERR)   
      ENDIF
      IF (MYID.EQ.MASTER .OR. JOB.EQ.2) THEN 
         ITER = 0 
         DO 10 IFREQ=1,NFREQ
            DO 11 IREC=1,NREC
               DO 12 ISRC=1,NSRC
                  DO 13 I=1,NDIM
                     ITER = ITER + 1 
                     EST(I,IFREQ,IREC,ISRC) = QBUFF(ITER)
   13             CONTINUE
   12          CONTINUE
   11       CONTINUE
   10    CONTINUE
      ENDIF
      DEALLOCATE(BUFF)
      DEALLOCATE(QBUFF) 
      RETURN
      END
