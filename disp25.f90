      implicit real*8 (a-h,o-z)
      real*8 gmat(6,6), rm(3,3), eigs(3), grpv(3,3), gvec(3), u(3,3), s(3)
!.... material arrays 
      CHARACTER(2), ALLOCATABLE :: CNNPG(:)      !anchor node domain, everything interior
      CHARACTER(1), ALLOCATABLE :: CDOMAIN(:)    !element domain, everything is interior
      COMPLEX*16, ALLOCATABLE :: UWAVE(:)        !analytic wavefield
      REAL*8, ALLOCATABLE :: XLOCS(:), ZLOCS(:)  !(x,z) locations of anchor nodes 
      REAL*8, ALLOCATABLE :: XIPTS(:), ETAPTS(:) !interpolation points
      REAL*8, ALLOCATABLE :: ECOEFF(:,:)         !lambda/mu at anchor nodes
      REAL*8, ALLOCATABLE :: DENS(:)             !density at anchor nodes
      INTEGER*4, ALLOCATABLE :: IENG(:,:)        !global IEN vector
!.... assembly matrices
      INTEGER*4, ALLOCATABLE :: LM(:,:,:)        !element to global DOF
      INTEGER*4, ALLOCATABLE :: IRPTR(:)         !CRS row pointer
      INTEGER*4, ALLOCATABLE :: JCPTR(:)         !CRS column pointer
      INTEGER*4, ALLOCATABLE :: MYDOFS(:)        !holds DOFs local to process, but since 
                                                 !we won't use distributed assembly this
                                                 !holds all the DOFs
      COMMON /CRS_BLOCK/ IRPTR_PTR, JCPTR_PTR    
      INTEGER*4, POINTER :: IRPTR_PTR(:), JCPTR_PTR(:) 
      REAL*8 KAPPA
!.... parameters
      LOGICAL*4 LISISO 
      PARAMETER(PI = 3.1415926535897931D0) !pi
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
      PARAMETER(NGNOD  = 4)      !quads only 
      PARAMETER(NDIM = 3)        !will have 3 components of motion
      PARAMETER(NELEM = 16)       !minimum number of elements for analysis
      PARAMETER(LISISO = .TRUE.) !isotropy only
      

      dpois = 0.d0
      dtheta = 0.d0
      dbaz = 0.d0
      dg = 0.1d0
      nordm = 8

      IITYPE = 1 !interpolation type
      G = 4.D0 !grid points per wavelength
      C = 4100.D0 !P wave velcotiy
      RHO = 2700.D0 !density 
      FREQ = 1000.D0   !characteristic frequency of 1 Hz
      POIS0 = 0.27D0 
      POIS1 = 0.4D0 
      BAZ0 = 0.D0 
      BAZE = 45.D0 
!
!.... things that drop out from parameters 
      NNPG = INT(SQRT(REAL(NELEM)) + 1)**2
      OMEGA = 2.D0*PI*FREQ
      DELTA = .5D0 !length of element
!
!.... set the anchor node locations
      ALLOCATE(ECOEFF(NNPG,2)) !elasticity only
      ALLOCATE(DENS(NNPG)) 
      ALLOCATE(XLOCS(NNPG)) 
      ALLOCATE(ZLOCS(NNPG)) 
      ALLOCATE(IENG(NGNOD,NELEM))
      ALLOCATE(CDOMAIN(NELEM)) 
      ALLOCATE(CNNPG(NNPG))
      CALL GEN_QMESH_LOCS(NGNOD,NNPG, NELEM,NGNOD, CDOMAIN,CNNPG,IENG,XLOCS,ZLOCS)
!
!.... an outer loop on polynomial order
      DO 100 IORD=2,2 !NORDM 
         NORD = IORD      !polynomial order
         NLXI = NORD + 1  !number of lagrange interpolant points in xi
         NLETA = NLXI     !conforming
         NEN = NLXI*NLETA !number of element nodes
         ALLOCATE(XIPTS(NLXI)) 
         ALLOCATE(ETAPTS(NLETA)) 
         CALL GENPTS(IITYPE, NLXI,NLETA, XIPTS,ETAPTS)
         print *, xipts
         print *, etapts
!
!....... generate the graph
         ALLOCATE(LM(NDIM,NEN,NELEM))
         CALL GRAPH_ONLY(NDIM,NEN,NNPG,NELEM,NLXI,NLETA, NGNOD, CDOMAIN,CNNPG, &
                         IENG(1:NGNOD,1:NELEM),XIPTS,ETAPTS,XLOCS,ZLOCS, &
                         NDOF,NZERO,LM(1:NDIM,1:NEN,1:NELEM), IERR) 
         ALLOCATE(IRPTR(NDOF+1))
         ALLOCATE(JCPTR(NZERO)) 
         IRPTR(1:NDOF+1) = IRPTR_PTR(1:NDOF+1)
         JCPTR(1:NZERO)  = JCPTR_PTR(1:NZERO) 
         DEALLOCATE(IRPTR_PTR)
         DEALLOCATE(JCPTR_PTR)
         ALLOCATE(MYDOFS(NDOF))
         DO 1 IDOF=1,NDOF
            MYDOFS(IDOF) = IDOF
    1    CONTINUE 
         ALLOCATE(UWAVE(NDOF)) 
         CALL ZSCAL(NDOF,DCMPLX(0.D0,0.D0),UWAVE,1)
!
!....... second loop on poissons ratios
         DO 200 IPOIS=1,1 !NPOIS
            POISON = POIS0 + DFLOAT(IPOIS - 1)*DPOIS 
            VP = C 
            CALL FILL_HOMO_ECD(NNPG,NNPG, VP,RHO,POISON, DENS,ECOEFF)
            WRITE(*,*) 'xdisp25: Vp, Poissons Ratio, Density:', &
            '',REAL(VP),REAL(POISON),REAL(RHO)
!
!.......... third loop on back-azimuth
            DO 300 IBAZ=1,1 !NBAZ
               BAZ = BAZ0 + DFLOAT(IBAZ - 1)*DBAZ 
               PHI = BAZ
!
!............. loop on angle of incidence
               DO 400 ITHETA=1,1 !NTHETA
                  THETA = 22.D0 + DFLOAT(ITHETA - 1)*DTHETA
                  AOI = THETA
                  PY = SIN(BAZ*PI180)*SIN(AOI*PI180)/VP
                  KAPPA = 2.D0*PI/(G*DELTA)  
                  RKAPX = KAPPA*SIN(AOI*PI180)
                  RKAPZ = KAPPA*COS(AOI*PI180) 
                  CALL FILL_UWAVE_QMESH(NDIM,NEN,NGNOD, NDOF,NNPG,NLXI,NLETA,NELEM,        &
                                        NGNOD,NDIM, .TRUE., RKAPX,RKAPZ,BAZ,AOI, IENG,LM, &
                                        XIPTS,ETAPTS,XLOCS,ZLOCS, UWAVE,IERR)

!
!................ we are now ready to assemble and calculate the dispersion 
                  CALL DISP25(NNPG,NGNOD,NEN,NDIM, NDOF,NZERO,NNPG,NGNOD,NEN, &
                              NLXI,NLETA,NELEM,NDIM, IITYPE, PY,THETA,KAPPA, &
                              MYDOFS,IRPTR,JCPTR, IENG,LM,  &
                              XIPTS,ETAPTS,XLOCS,ZLOCS,DENS,ECOEFF, Uwave)
!
!................ loop on grid points per wavelength
                  DO 500 IGRD=1,1 !NGRD
                     G = 1.D0 !grid points per wavelength

  500             CONTINUE !loop on grid points per wavelength
  400          CONTINUE !loop on theta
  300       CONTINUE !loop on back-azimuth
  200    CONTINUE  !loop on poissons ratios
         DEALLOCATE(LM) 
         DEALLOCATE(IRPTR)
         DEALLOCATE(JCPTR) 
         DEALLOCATE(UWAVE) 
         DEALLOCATE(XIPTS)
         DEALLOCATE(ETAPTS) 
  100 CONTINUE !loop on polynomial order
!
!.... clean remaining stuff
 5500 CONTINUE !error break
      DEALLOCATE(ECOEFF)
      DEALLOCATE(DENS) 
      DEALLOCATE(XLOCS)
      DEALLOCATE(ZLOCS)
      DEALLOCATE(IENG) 
      DEALLOCATE(CNNPG) 
      DEALLOCATE(CDOMAIN) 
 
      STOP
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE DISP25(MNPG,MGNOD,MEN,MDIM, NDOF,NZERO,NNPG,NGNOD,NEN,        &
                        NLXI,NLETA,NELEM,NDIM,  IITYPE, PY,THETA,KAPPA,             &
                        MYDOFS,IRPTR,JCPTR, IENG,LM, &
                        XIPTS,ETAPTS,XLOCS,ZLOCS,DENS,ECOEFF, UE)
!
!     This program calculates the dispersion error for a 2.5D medium

      COMPLEX*16, INTENT(IN) :: UE(NDOF) 
      REAL*8, INTENT(IN) :: ECOEFF(MNPG,*), XLOCS(NNPG), ZLOCS(NNPG), DENS(NNPG),  &
                            XIPTS(NLXI),ETAPTS(NLETA),PY,THETA,KAPPA 
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), MYDOFS(NDOF), &
                               IRPTR(NDOF+1), JCPTR(NZERO), &
                               MNPG,MGNOD,MEN,MDIM, NDOF,NZERO,NGNOD,NEN,NLXI,NLETA,NDIM,&
                               IITYPE   

!.... local variables
      REAL*8, ALLOCATABLE :: SHG(:,:,:,:), SHL(:,:,:,:), DMAT(:,:,:), &
                             XIGLL(:,:), ETAGLL(:,:), RHO(:,:), DET(:,:)   
      COMPLEX*16, ALLOCATABLE :: ZMAT(:,:,:), CE(:,:), U(:), KVEC(:), EIGS(:) 
      REAL*8, ALLOCATABLE :: ME(:,:), MVEC(:) 
      REAL*8 KAPPAN, KAPPAX, KAPPAZ, PI180, PERT  
      COMPLEX*16 CZERO
      LOGICAL*4 LISISO 
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0)) 
      PARAMETER(LISISO = .TRUE.) !must be isotropic 
      PARAMETER(NPERT = 3) !forward + backward difference and analytic 
      PARAMETER(PERT = 1.D-7) 
      PARAMETER(PI180 = 0.017453292519943295D0) !pi/180
!
!----------------------------------------------------------------------------------------!
!
!.... set parameters
      NORD = NLXI - 1
      NINTX = 2*NORD + 2 + 2           !Lobatto quadruature, k = 2n - 3, k -> 2k for mass
      NINTZ = NINTX                    !conforming
      THETAR = THETA*PI180 
!
!.... set the grid and pointers
      ALLOCATE(SHL(3,NEN,NINTX,NINTZ))
      ALLOCATE(XIGLL(NINTX,2))
      ALLOCATE(ETAGLL(NINTZ,2)) 
      CALL DSHL(NEN,NINTX,NINTZ, IITYPE,NINTX,NINTZ, NLXI,NLETA, &
                XIPTS,ETAPTS, XIGLL,ETAGLL,SHL)
      ALLOCATE(KVEC(NZERO))
      ALLOCATE(MVEC(NZERO))
      CALL DSCAL(NZERO, 0.D0,MVEC,1)
      CALL ZSCAL(NZERO,CZERO,KVEC,1) 
      ALLOCATE(SHG(3,NEN,NINTX,NINTZ))
      ALLOCATE(RHO (NINTX,NINTZ))
      ALLOCATE(DET (NINTX,NINTZ))
      ALLOCATE(DMAT(NINTX,NINTZ,2))
      NEE = NDIM*NEN
      ALLOCATE(CE  (NEE,NEE))
      ALLOCATE(ME  (NEE,NEE)) 
!
!.... loop on elements and assemble
      DO 50 IELEM=1,NELEM
         CALL CJMK_TVI(NEN,NINTX,NINTZ,MNPG, LISISO,                              &   
                       NINTX,NINTZ,NEN,NGNOD,NNPG, IENG(1:NGNOD,IELEM),           &   
                       XIGLL(:,1),ETAGLL(:,1),XLOCS,ZLOCS,DENS,ECOEFF,SHL,        &   
                       DET,RHO,DMAT,SHG, IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'disp25: Error calling cjmk_tvi!'
            RETURN
         ENDIF
         CALL KM25ISOQ(NEN,NEE,NINTX, NEN,NINTX,NINTZ, &
                       1.D0,PY, XIGLL(:,2),ETAGLL(:,2), DMAT(:,:,1),DMAT(:,:,2), & 
                       RHO,DET,SHG, ME,CE)  
         CALL  ASM_DISP(NEE,MDIM,NZERO,NDOF, NEN,NDIM,  &
                        MYDOFS,IRPTR,JCPTR,LM(1:NDIM,1:NEN,IELEM), ME,CE,  &
                        MVEC,KVEC,IERR)
         IF (IERR /= 0) THEN 
            WRITE(*,*) 'disp25: Error calling asm_disp!'
            RETURN
         ENDIF

   50 CONTINUE 
!
!.... set the matrices in the generalized eigenvalue problem  
      IELEM = 11 !Only interested in the top right element
      NEQN = (NORD*NDIM)**2
      ALLOCATE(ZMAT(3,NEQN,NEQN))
      DO 100 IPERT=1,1 !NPERT
         IF (IPERT == 1) THEN !standard evaluation
            KAPPAN = KAPPA  
         ELSE !set up directional derivative
            KAPPAN = KAPPA + PERT 
         ENDIF
         KAPPAX = KAPPAN*COS(THETAR)
         KAPPAZ = KAPPAN*SIN(THETAR)
!
!....... group my complex matrices
         DO 13 I=1,NEQN
            DO 14 J=1,NEQN 
               ZMAT(1,I,J) = CZERO 
               ZMAT(2,I,J) = CZERO
               ZMAT(3,I,J) = CZERO
   14       CONTINUE
   13    CONTINUE
!
!....... loop on rows
         IEQN = 0 
         DO 21 ILETA=1,NLETA-1
            DO 22 ILXI=1,NLXI-1 
               IAE = (ILETA - 1)*NLXI + ILXI
               DO 23 I=1,NDIM
                  IEQN = IEQN + 1
                  IDOF = LM(I,IAE,IELEM)
                  IF (IDOF == 0) THEN
                     WRITE(*,*) 'disp25: Quite peculiar'
                     GOTO 170
                  ENDIF
                  IDOFL = IBSECT2(NDOF,IDOF,MYDOFS)
                  IF (IDOFL > 0) THEN
                     JBEG = IRPTR(IDOFL)
                     JEND = IRPTR(IDOFL+1) - 1
                     NVAR = JEND - JBEG + 1 
                  ELSE
                     WRITE(*,*) 'disp25: You are wrong'
                     IERR = 1
                     RETURN
                  ENDIF
!
!................ and now multiply the columns 
                  JEQN = 0 
                  DO 25 JLETA=1,NLETA-1
                     DO 26 JLXI=1,NLXI-1
                        IBE = (JLETA - 1)*NLXI + JLXI
                        DO 27 J=1,NDIM
                           JEQN = JEQN + 1 
                           JDOF = LM(J,IBE,IELEM)
                           INDX = IBSECT(NVAR,JDOF,JCPTR(JBEG:JEND))
                           IF (INDX < 0) THEN
                              WRITE(*,*) 'disp25: You are definitely wrong'
                              IERR = 1 
                              RETURN
                           ENDIF
                           IZLOC = JBEG - 1 + INDX
                           IF (I == 1 .AND. J == 2 .OR. &
                               I == 2 .AND. J == 1 .OR. &
                               I == 2 .AND. J == 3 .OR. &
                               I == 3 .AND. J == 2) THEN
                              ZMAT(2,IEQN,JEQN) = KVEC(IZLOC)*UE(JDOF) 
                           ELSE
                              ZMAT(1,IEQN,JEQN) = KVEC(IZLOC)*UE(JDOF)
                           ENDIF
                           IF (ABS(MVEC(IZLOC)) > 2.22D-15) &
                           ZMAT(3,IEQN,JEQN) = DCMPLX(MVEC(IZLOC),0.D0)*UE(JDOF) 
   27                   CONTINUE !loop on components
   26                CONTINUE !loop on xi
   25             CONTINUE !loop on eta
  
!                 !IORD = MYTAG(IDOF) 
!                 DO 18 J=JSTRT,JSTOP
!                    JDOF = JCPTR(J) 
!                    !JORD = MYTAG(JDOF) 
!                    ZME(IORD,JORD) = ZME(IORD,JORD) + DCMPLX(MVEC(J),0.D0)*ZSHIFT*U(JDOF)
!                    ZKE(IORD,JORD) = ZKE(IORD,JORD) + KVEC(J)*ZSHIFT*U(JDOF)
!  18             CONTINUE 
  170             CONTINUE !not a DOF, really?
   23          CONTINUE
   22       CONTINUE
   21    CONTINUE
  100 CONTINUE 
      ALLOCATE(EIGS(2*NEQN)) 
      CALL QPEEIG2(3,NEQN,NEQN,2, ZMAT, EIGS, IERR)
      do i=1,neqn
         print *, eigs(i)/kappa
      enddo
!
!.... clean 
      DEALLOCATE(SHL) 
      DEALLOCATE(XIGLL) 
      DEALLOCATE(ETAGLL) 
      DEALLOCATE(SHG) 
      DEALLOCATE(RHO)
      DEALLOCATE(DET)
      DEALLOCATE(DMAT) 
      DEALLOCATE(CE)
      DEALLOCATE(ME) 
      DEALLOCATE(MVEC)
      DEALLOCATE(KVEC) 
      DEALLOCATE(ZMAT) 
      DEALLOCATE(EIGS) 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GEN_QMESH_LOCS(MGNOD,NNPG, NELEM,NGNOD, CDOMAIN,CNNPG,IENG,XLOCS,ZLOCS)
!
!     Generates the IEN anchor node pointer with the x z locations
!
!     INPUT      MEANING
!     -----      ------- 
!     MGNOD      leading dimension for IENG
!     NELEM      number of elements (a perfect square)
!     NNPG       number of anchor nodes
!     
!     OUTPUT     MEANING
!     ------     ------- 
!     CDOMAIN    element domain; everything is in interior
!     CNNPG      anchor node domain; everything is inteiror
!     IENG       anchor node IEN pointer
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!
!.... variable declarations
      INTEGER*4, INTENT(IN) :: MGNOD,NNPG, NELEM,NGNOD 
      CHARACTER*2, INTENT(OUT) :: CNNPG(NNPG) 
      CHARACTER*1, INTENT(OUT) :: CDOMAIN(NELEM)
      REAL*8, INTENT(OUT) :: XLOCS(NNPG), ZLOCS(NNPG) 
      INTEGER*4, INTENT(OUT) :: IENG(MGNOD,*) 
!.... local variables
      REAL*8 X0, Z0, XP, ZP, DX, DZ 
      PARAMETER(X0 =-1.D0, Z0 =-1.D0, XP = 1.D0, ZP = 1.D0)
!
!----------------------------------------------------------------------------------------!
!
      NELEMX = INT(SQRT(REAL(NELEM))) 
      NELEMZ = NELEMX 
      DX = (XP - X0)/DFLOAT(NELEMX)
      DZ = DX 
      INPG = 0
      DO 1 IELEMZ=1,NELEMZ+1
         DO 2 IELEMX=1,NELEMX+1
            INPG = INPG + 1
            XLOCS(INPG) = DFLOAT(IELEMX - 1)*DX + X0  
            ZLOCS(INPG) = DFLOAT(IELEMZ - 1)*DZ + Z0 
            CNNPG(INPG) = 'II'
    2    CONTINUE
    1 CONTINUE 
!
!.... set the anchor node pointer
      IENG(1:NGNOD,1:NELEM) = 0
      NNPGX = NELEMX + 1
      IELEM = 0
      DO 3 IELEMZ=1,NELEMZ
         DO 4 IELEMX=1,NELEMX
            IELEM = IELEM + 1
            IA1 = (IELEMZ - 1)*NNPGX + IELEMX
            IA2 = IA1 + 1
            IA4 = IA1 + NNPGX
            IA3 = IA4 + 1
            IENG(1,IELEM) = IA1
            IENG(2,IELEM) = IA2
            IENG(3,IELEM) = IA3
            IENG(4,IELEM) = IA4
            CDOMAIN(IELEM) = 'I'
    4    CONTINUE
    3 CONTINUE
      IF (MAXVAL(IENG(1:NGNOD,1:NELEM)) < 0) THEN
         WRITE(*,*) 'gen_disp_qmesh: Failed to initialize a point!'
         IERR = 1
      ENDIF
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE GRAPH_ONLY(NDIM,NEN,NNPG,NELEM,NLXI,NLETA, NGNOD, CDOMAIN,CNNPG, &
                            IENG,XIPTS,ETAPTS,XLOCS,ZLOCS, &
                            NDOF,NZERO,LM, IERR) 
!
!     Driver routine that generates the global graph in CRS format
!
!     INPUT      MEANING
!     -----      ------- 
!     CDOMAIN    element domain (absorbing, bielak, interior)
!     CNNPG      anchor node domain 
!     ETAPTS     eta lagrange interpolation points
!     IENG       anchor node IEN pointer
!     NDIM       number of components in solution
!     NELEM      number of elements
!     NEN        number of element nodes
!     NGNOD      number of anchor nodes 
!     NLETA      number of lagrange interpolant points in eta
!     NLXI       number of lagrange interpolant points in xi
!     NNPG       number of anchor nodes
!     XIPTS      xi lagrange interpolation points
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes
!      
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     IRPTR      CRS row pointer 
!     JCPTR      CRS column pointer
!     LM         element to global DOF
!     NDOF       number of degrees of freedom
!     NZERO      number of non-zeros in global matrix
! 
!.... variable declarations
      IMPLICIT NONE
      CHARACTER*2, INTENT(IN) :: CNNPG(NNPG)
      CHARACTER*1, INTENT(IN) :: CDOMAIN(NELEM)
      REAL*8, INTENT(IN) :: XLOCS(NNPG), ZLOCS(NNPG), XIPTS(NLXI), ETAPTS(NLETA) 
      INTEGER*4, INTENT(IN) :: IENG(NGNOD,*),NDIM,NEN,NNPG,NELEM,NLXI,NLETA, NGNOD
      INTEGER*4, INTENT(OUT) :: LM(NDIM,NEN,*), NZERO, NDOF, IERR
      COMMON /CRS_BLOCK/ IRPTR,JCPTR
      INTEGER*4, POINTER :: IRPTR(:), JCPTR(:)
!.... local variables
      INTEGER*4, ALLOCATABLE :: IEN(:,:), ID(:,:), MCONN(:,:)
      INTEGER*4 NNP, NCONF, MSPACE, IDOF
      INTEGER*4 ICMCON, IENZERO
!
!----------------------------------------------------------------------------------------!
!
!.... generate the graph
      IERR = 0
      ALLOCATE(IEN(NEN,NELEM))
      WRITE(*,*) 'graph_only: Generating IEN matrix...'
      CALL GENIEN(NGNOD,NEN,NNPG, NELEM,NLXI,NLETA,NGNOD, .FALSE.,  &
                  NEN,IENG,XIPTS,ETAPTS,XLOCS,ZLOCS, NNP,IEN,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'graph_only: Error generating IEN matrix...' 
         RETURN 
      ENDIF 
      NCONF = ICMCON(NEN,NNP, NELEM,NEN,IEN)
      ALLOCATE(MCONN(NNP,NCONF))
      CALL CONNECT(NEN,NNP, NNP,NELEM,NEN,NCONF, IEN,MCONN)
      WRITE(*,*) 'graph_only: Generating ID vector...'
      ALLOCATE(ID(NDIM,NNP))
      CALL GENID(NDIM,NGNOD,NEN, NNPG,NNP,NELEM,NLXI,NLETA,       &   
                 NGNOD,NDIM, CDOMAIN,CNNPG, IENG,IEN,             &   
                 XIPTS,ETAPTS,XLOCS,ZLOCS, NDOF,ID,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'graph_only: Error calling genid'
         RETURN 
      ENDIF
      WRITE(*,*) 'graph_only: Generating LM matrix...'
      CALL GENLM(NDIM,NEN,NEN,NNP, NDIM,NELEM,ID,IEN, LM) 
      MSPACE = IENZERO(NNP,NDIM, NNP,NDIM,NCONF,NEN, ID,MCONN)
      DEALLOCATE(ID)
      WRITE(*,*) 'graph_only: Generating compressed row storage format...'
      ALLOCATE(IRPTR(NDOF+1))
      ALLOCATE(JCPTR(MSPACE))
      CALL GENCRS_QUICK(NDIM,NEN,NNP,MSPACE, NCONF,NELEM,NEN,NDIM,  &
                        NDOF,MCONN,IEN,LM, NZERO,IRPTR,JCPTR, IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'graph_only: Error calling gencrs_quick'
         RETURN
      ENDIF
      DEALLOCATE(IEN)
      DEALLOCATE(MCONN)
      RETURN
      END

!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FILL_HOMO_ECD(MNPG,NNPG, VP,RHO,POISON, DENS,ECOEFF) 
!
!     Fills the density and isotropic stiffness parameters at anchor nodes
!
!     INPUT      MEANING
!     -----      ------- 
!     MNPG       leading dimension for ECOEFF
!     NNPG       number of anchor nodes 
!     POISON     Poisson's ratio 
!     RHO        constant density (kg/m**3)
!     VP         constant compressional velocity (m/s)
!
!     OUTPUT     MEANING
!     ------     ------- 
!     DENS       constant density at achor nodes (kg/m**3) 
!     ECOEFF     (lambda,mu) at anchor nodes (Pa)
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: VP,RHO,POISON
      INTEGER*4, INTENT(IN) :: MNPG,NNPG
      REAL*8, INTENT(OUT) :: ECOEFF(MNPG,*), DENS(NNPG) 
!.... local variables
      REAL*8 ALPHA, BETA, RLAM, RMU 
      INTEGER*4 INPG
!
!----------------------------------------------------------------------------------------!
!
!.... set the materials
      ALPHA = VP
      BETA =  ALPHA*SQRT( (0.5D0 - POISON)/(1.D0 - POISON) ) 
      DO 1 INPG=1,NNPG
         RMU = RHO*BETA**2
         RLAM = RHO*ALPHA**2 - 2.D0*RMU 
         ECOEFF(INPG,1) = RLAM 
         ECOEFF(INPG,2) = RMU
         DENS(INPG) = RHO
    1 CONTINUE !loop on anchor nodes 
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE QPEEIG2(LDA1,LDA2,N,NP, AIN, EIGS, IERR)
!
!     Fortran implementation of polynomial eigenvalue solver:
!     http://mind.cog.jhu.edu/courses/680/octave/Installers/
!     Octave/Octave.OSX10.6/Applications/MATLAB_R2009b.app/
!     toolbox/matlab/matfun/polyeig.m
!  
!     For (A_0 + lambda A_1 + lambda^2 A_2)x = 0
!
!     Step 1: Build the [n*p x n*p] matrices
!
!     A = [A0  0  0]   B = [-A1  -A2  -A3]
!         [ 0  I  0]       [  I    0    0]
!         [ 0  0  I]       [  0    I    0]
!
!     I only need the eigenvalues so this is all I bother to 
!     return - B. Baker January 2013
!
!     INPUT      MEANING
!     -----      ------- 
!     AIN        A matrices in order A0, A1, A2... as seen above
!     LDA1       leading dimension 1 for AIN; >= np + 1
!     LDA2       leading dimension 2 for AIN; >= N
!     N          number of rows/columns of each matrix
!     NP         order of polynomial: ex np is givne above
!
!     OUTPUT     MEANING
!     ------     ------- 
!     EIGS       eigenvalues 
!     IERR       error flag
!
!.... variable declarations 
      COMPLEX*16, INTENT(IN) :: AIN(LDA1,LDA2,n)
      COMPLEX*16, INTENT(OUT) :: EIGS(*)
      INTEGER*4, INTENT(IN) :: LDA1, LDA2, N, NP  
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*16, ALLOCATABLE :: A(:,:), B(:,:), VL(:,:), VR(:,:), & 
                                 ALPHA(:), BETA(:), WORK(:)
      REAL*8, ALLOCATABLE :: RWORK(:)  
      LOGICAL*4, ALLOCATABLE :: LKEEP(:)
      CHARACTER*1 JOBVL, JOBVR 
      COMPLEX*16 CZERO, CONE, Z
      REAL*8 ZR, ZI  
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0))
      PARAMETER(CONE  = DCMPLX(1.D0,0.D0)) 
!
!----------------------------------------------------------------------------------------!
!
!.... initialize, set space 
      IERR = 0
      NN = NP*N
      ALLOCATE(A(NN,NN))
      ALLOCATE(B(NN,NN))
!
!.... build the two [n*p x n*p] matrices
      I1 = 0
      DO 1 I=1,NN
         I1 = I1 + 1
         JTRIP = 1
         IP = 0
         J1 = 0
         DO 2 J=1,NN
            IPF = IP + 1 !fortran indexing
            J1 = J1 + 1
            A(I,J) = CZERO
            IF (I <= N .AND. J <= N) THEN
               A(I,J) = AIN(1,I1,J1)
            ELSE
               IF (I == J) THEN !identity on diagonal
                  A(I,J) = CONE
               ELSE
                  A(I,J) = CZERO
               ENDIF
            ENDIF
            B(I,J) = CZERO
            IF (I <= N) THEN
               B(I,J) =-AIN(IPF+1,I1,J1)
            ELSE
               IF (I > (IP+1)*N .AND. I <= (IP+2)*N) THEN
                  IF (I1 == J1) B(I,J) = CONE
               ENDIF
            ENDIF
            IF (J1 == N) THEN
               IP = IP + 1
               J1 = 0
            ENDIF
    2    CONTINUE
         IF (I1 == N) I1 = 0
    1 CONTINUE
!
!.... solve the right generalized eigenvalue problem
      JOBVL = 'N'  !'V'
      JOBVR = 'N'  !'V'
      LDA  = NN
      LDB  = NN
      LDVL = 1  !NN
      LDVR = 1  !NN
      ALLOCATE(VL(LDVL,NN))
      ALLOCATE(VR(LDVR,NN))
      ALLOCATE(ALPHA(NN))
      ALLOCATE(BETA(NN))
      LWORK = MAX(1,2*NN) + 10 !add some extra
      ALLOCATE(WORK(LWORK))
      ALLOCATE(RWORK(8*NN))
      CALL ZGGEV(JOBVL,JOBVR,NN,A,LDA,B,LDB, ALPHA,BETA,   &
                 VL, LDVL, VR, LDVR, WORK,LWORK, RWORK, INFO)
!
!.... clean
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(VR)
      DEALLOCATE(VL)
      DEALLOCATE(A)
      DEALLOCATE(B)
!
!.... error check
      IF (INFO < 0) THEN
         WRITE(*,*) 'qpeeig2: Error illegal argument:',INFO
         IERR = 1
         RETURN
      ELSEIF (INFO > 0) THEN
         IF (INFO >= 1 .AND.INFO <= NN) THEN
            WRITE(*,*) 'qpeeig2: No eigenvectors, eigenvalues up to',INFO
            IERR = 1
            RETURN
         ELSE
            IF (INFO == NN+1) THEN
               WRITE(*,*) 'qpeeig2: zggev failed in dhgeqz'
            ELSE
               WRITE(*,*) 'qpeeig2: zggev failed in dtgevc',INFO
            ENDIF
            IERR = 1
            RETURN
         ENDIF
      ENDIF
      ALLOCATE(LKEEP(NN))
      DO 10 I=1,NN
         IF (CDABS(BETA(I)) > 2.22D-16) THEN
            Z = ALPHA(I)/BETA(I)
            ZR = DREAL(Z)
            ZI = DIMAG(Z)
            IF (ABS(ZR) < 2.22D-15) ZR = 0.D0
            IF (ABS(ZI) < 2.22D-15) ZI = 0.D0
            EIGS(I) = DCMPLX(ZR,ZI)
            LKEEP(I) = .TRUE.
            print *, cmplx(eigs(i)) 
         ELSE
            EIGS(I) = CZERO
            LKEEP(I) = .FALSE.
         ENDIF
   10 CONTINUE
      DEALLOCATE(ALPHA)
      DEALLOCATE(BETA)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE FILL_UWAVE_QMESH(MDIM,MEN,MGNOD, NDOF,NNPG,NLXI,NLETA,NELEM,           &
                                  NGNOD,NDIM, LCHECK, KAPPAX,KAPPAZ,BAZ,AOI, IENG,LM,   &
                                  XIPTS,ETAPTS,XLOCS,ZLOCS, UE,IERR)
!
!     Fills the analytic solution vector for a plane wave in homogeneous unbounded media
!
!     INPUT      MEANING
!     -----      ------- 
!     AOI        angle of indidence (degrees)
!     BAZ        back-azimuth (degrees)
!     ETAPTS     eta interpolation points
!     IENG       IEN pointer for anchor nodes
!     KAPPAX     x wavenumber rad/m
!     KAPPAZ     z wavenumber rad/m
!     LM         maps element DOFs to global DOFs
!     MDIM       leading dimension for LM
!     MEN        max number of element nodes
!     MGNOD      max number of anchor nodes
!     NDIM       number of components in solution
!     NDOF       number of degrees of freedom
!     NELEM      number of elements
!     NGNOD      number of anchor nodes
!     NLETA      number of lagrange interpolant points in eta
!     NLXI       number of lagrange interpolant points in xi
!     NNPG       number of anchor nodes in geometry
!     XIPTS      xi interpolation points 
!     XLOCS      x locations of anchor nodes
!     ZLOCS      z locations of anchor nodes 
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     UE         exact solution 
!
!.... variable declarations
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XIPTS(NLXI), ETAPTS(NLETA), XLOCS(NNPG), ZLOCS(NNPG), &
                            BAZ, AOI, KAPPAX, KAPPAZ  
      INTEGER*4, INTENT(IN) :: LM(MDIM,MEN,*), IENG(MGNOD,*), MDIM, MEN, MGNOD, &
                               NDOF, NNPG, NLXI, NLETA, NELEM, NGNOD,NDIM
      LOGICAL*4, INTENT(IN) :: LCHECK
      COMPLEX*16, INTENT(OUT) :: UE(NDOF)
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables 
      REAL*8, ALLOCATABLE :: SF(:), XOUT(:), ZOUT(:), WORK(:)   
      INTEGER*4, ALLOCATABLE :: IPERM(:), JPERM(:) 
      LOGICAL*4, ALLOCATABLE :: LINIT(:) 
      COMPLEX*16 CCOSAZ, CSINAZ, CCOSI, CSINI, CARG, U, V, W, CZERO
      REAL*8 ETA,XI,X,Z, COSAZ,SINAZ,COSI,SINI,PI180
      INTEGER*4 IDOFS(3), IELEM, ILETA, ILXI, IAE, I, IDOF, ILOC, IA, NNP, INP, JNP, KNP 
      LOGICAL*4 LBLANK 
      PARAMETER(CZERO = DCMPLX(0.D0,0.D0)) 
      PARAMETER(PI180 = 0.017453292519943295D0)
!
!----------------------------------------------------------------------------------------!
!
!.... initialize
      IERR = 0
      COSAZ = COS(BAZ*PI180)
      SINAZ = SIN(BAZ*PI180)
      COSI = COS(AOI*PI180)
      SINI = SIN(AOI*PI180)
      CCOSAZ = DCMPLX(COSAZ) 
      CSINAZ = DCMPLX(SINAZ)
      CCOSI = DCMPLX(COSI)
      CSINI = DCMPLX(SINI)
      ALLOCATE(SF(NGNOD)) 
      ALLOCATE(LINIT(NDOF)) 
      LINIT(1:NDOF) = .FALSE. 
      CALL ZSCAL(NDOF,CZERO,UE,1) 
      IF (LCHECK) THEN
         NNP = NDOF/NDIM
         ALLOCATE(XOUT(NNP))
         ALLOCATE(ZOUT(NNP))
         ALLOCATE(IPERM(NNP))
      ENDIF
         
!
!.... loop on elements 
      DO 1 IELEM=1,NELEM
         DO 2 ILETA=1,NLETA
            ETA = ETAPTS(ILETA)
            DO 3 ILXI=1,NLXI
               IAE = (ILETA - 1)*NLXI + ILXI
!
!............. check if point has been initialized
               DO 4 I=1,NDIM
                  IDOF = LM(I,IAE,IELEM)
                  IF (IDOF > 0) THEN
                     IF (LINIT(IDOF)) GOTO 300 
                  ENDIF
    4          CONTINUE 

               XI  = XIPTS(ILXI) 
               CALL CSF2DND(NGNOD,XI,ETA, SF, IERR) 
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'fill_uwave_qmesh: Error calling csf2dnd'
                  RETURN
               ENDIF 
               X = 0.D0
               Z = 0.D0 
               DO 5 IA=1,NGNOD
                  ILOC = IENG(IA,IELEM) 
                  X = X + SF(IA)*XLOCS(ILOC)
                  Z = Z + SF(IA)*ZLOCS(ILOC)
    5          CONTINUE !loop on anchor nodes
!
!............. plane wave evaluation
               CARG = CDEXP(DCMPLX(0.D0,-KAPPAX*X - KAPPAZ*Z))
               U = CCOSAZ*CSINI*CARG
               V = CSINAZ*CSINI*CARG
               W =        CCOSI*CARG
               DO 6 I=1,NDIM
                  IDOF = LM(I,IAE,IELEM) 
                  IF (IDOF > 0) THEN
                     LINIT(IDOF) = .TRUE.
                     IF (LCHECK .AND. I == 1) THEN
                        INP = (IDOF + (NDIM - 1))/NDIM
                        XOUT(INP) = X
                        ZOUT(INP) = Z
                        IPERM(INP) = INP
                     ENDIF
                     IF (NDIM == 2) THEN
                        IF (I == 1) THEN
                           UE(IDOF) = U
                        ELSE
                           UE(IDOF) = W
                        ENDIF
                     ELSE
                        IF (I == 1) THEN
                           UE(IDOF) = U
                        ELSEIF (I == 2) THEN
                           UE(IDOF) = V
                        ELSE
                           UE(IDOF) = W 
                        ENDIF
                     ENDIF 
                  ENDIF
    6          CONTINUE !loop on components
  300          CONTINUE !break ahead, point done
    3       CONTINUE !loop on xi points
    2    CONTINUE !loop on eta points
    1 CONTINUE !loop on elements
      DEALLOCATE(SF) 
      DO 10 IDOF=1,NDOF
         IF (.NOT.LINIT(IDOF)) THEN
            WRITE(*,*) 'fill_uwave_qmesh: Did not initialize dof:',IDOF
            IERR = 1
         ENDIF
   10 CONTINUE
      DEALLOCATE(LINIT)
!
!.... write wavefield
      IF (LCHECK) THEN
         ALLOCATE(WORK(NNP))
         ALLOCATE(JPERM(NNP))
         CALL DCOPY(NNP,ZOUT,1,WORK,1) 
         JPERM(1:NNP) = IPERM(1:NNP) 
         CALL SHELLR8I2(NNP, WORK,IPERM)
         WORK(1:NNP) = 0.D0  
         JNP = 0
         KNP = 1
         DO 20 INP=1,NNP
            LBLANK = .FALSE. 
            JNP  = JNP + 1
            WORK(JNP) = XOUT(IPERM(INP))  
            JPERM(INP) = IPERM(INP) 
            IF (INP < NNP) THEN
               IF (ABS(ZOUT(IPERM(INP)) - ZOUT(IPERM(INP+1))) > 1.11D-7) THEN 
                  CALL SHELLR8I2(JNP, WORK,JPERM(KNP:INP)) 
                  KNP = INP + 1
                  LBLANK = .TRUE. 
                  JNP = 0
               ENDIF
            ELSE
               CALL SHELLR8I2(JNP, WORK,JPERM(KNP:NNP))
            ENDIF
   20    CONTINUE
         DEALLOCATE(IPERM) 
         OPEN(UNIT=44,FILE='ue_wave_gnu.dat')
         DO 21 INP=1,NNP
            JNP = JPERM(INP) 
            DO 22 I=1,NDIM
               IDOFS(I) = (JNP - 1)*NDIM + I 
   22       CONTINUE
            WRITE(44,*) SNGL(XOUT(JNP)),SNGL(ZOUT(JNP)), &
                        SNGL(REAL(UE(IDOFS(1:3)))),SNGL(IMAG(UE(IDOFS(1:3))))
            IF (INP < NNP) THEN
               IF (ABS(ZOUT(JNP) - ZOUT(JPERM(INP+1))) > 1.11D-7) WRITE(44,*)
            ENDIF
   21    CONTINUE
         CLOSE(44) 
         DEALLOCATE(JPERM)
         DEALLOCATE(XOUT) 
         DEALLOCATE(ZOUT)
         DEALLOCATE(WORK)
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE ASM_DISP(MEE,MDIM,NZLOC,NDOFL, NEN,NDIM,  &
                          MYDOFS,IAOLOC,JAOLOC,LM, ME,KE,  &
                          MVEC,KVEC,IERR)
! 
!     Assembles the sparse matrix completely on one process.  You 
!     should use ASMBLE_DIST instead and put all dofs on MYDOFS 
! 
!     INPUT      MEANING
!     -----      ------- 
!     CE         element damping matrix
!     IRPTR      CRS row pointer 
!     JCPTR      CRS column pointer 
!     LM         location matrix 
!     MDIM       max number of spatial dimensions
!     ME         element mass
!     MEE        max number of element equations
!     NDIM       number of spatial dimensions
!     NDOF       number of degrees of freedom 
!     NEN        number of element nodes
!     NZERO      number of non-zeros in global matrices 
! 
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error flag
!     KVEC       global stiffness matrix
!     MVEC       global mass matrix
!  
!.... variable declarations
      IMPLICIT NONE
      COMPLEX*16, INTENT(IN) :: KE(MEE,*)
      REAL*8, INTENT(IN) :: ME(MEE,*)
      INTEGER*4, INTENT(IN) :: LM(MDIM,*), MYDOFS(NDOFL), & 
                 IAOLOC(NDOFL+1), JAOLOC(NZLOC), MEE,MDIM,NZLOC,NDOFL, & 
                 NEN,NDIM 
      COMPLEX*16, INTENT(INOUT) :: KVEC(NZLOC)
      REAL*8, INTENT(OUT) :: MVEC(NZLOC) 
      INTEGER*4, INTENT(OUT) :: IERR
      INTEGER*4 IBSECT2,IBSECT, IAE,IDOF,IDOFL,IBE,JDOF,JBEG,JEND,  &
                NVAR,IZLOC, IPE,IQE,I,J,INDX 
! 
!----------------------------------------------------------------------------------------!
! 
!.... loop on element nodes
      IERR = 0 
      IPE = 0 
      DO 1 IAE=1,NEN
         DO 2 I=1,NDIM
            IPE = IPE + 1 
            IDOF = LM(I,IAE)
            IF (IDOF /= 0) THEN
               IDOFL = IBSECT2(NDOFL,IDOF,MYDOFS)
               IF (IDOFL > 0) THEN
                  JBEG = IAOLOC(IDOFL)
                  JEND = IAOLOC(IDOFL+1) - 1
                  NVAR = JEND - JBEG + 1
                  IQE = 0
                  DO 3 IBE=1,NEN
                     DO 4 J=1,NDIM
                        IQE = IQE + 1
                        JDOF = LM(J,IBE)
                        IF (JDOF /= 0) THEN
                           INDX = IBSECT(NVAR,JDOF,JAOLOC(JBEG:JEND))
                           IF (INDX <= 0) THEN
                              IERR = 1
                              WRITE(*,*) 'asmble_cdist: Error INDX < 0'
                              RETURN
                           ENDIF
                           IZLOC = JBEG - 1 + INDX
                           KVEC(IZLOC) = KVEC(IZLOC) + KE(IPE,IQE) 
                           MVEC(IZLOC) = MVEC(IZLOC) + ME(IPE,IQE) 
                        ENDIF !endif on jdof
    4                CONTINUE !loop on ndim 
    3             CONTINUE !loop on ibe
               ENDIF !end check if DOF belongs to this process
            ENDIF !endif on idof
    2    CONTINUE !loop on ndim
    1 CONTINUE !loop on iae 
      RETURN
      END
