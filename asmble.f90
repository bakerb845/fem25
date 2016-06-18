      SUBROUTINE ASMBLE_DRIVER(NZLOC,CFTYPE,FREQ, &
                               PY, MSH, CDIST,IERR)
!
!     Assembles the global impedance matrix 
!
!     INPUT      MEANING
!     -----      ------- 
!     CDOMAIN    holds the element domain 
!     CFTYPE     frequency type
!     DENS       density at anchor nodes 
!     ECOEFF     elastic coefficients at anchor nodes
!     ETAPTS     eta interpolation points
!     FREQ       current frequency (Hz)
!     IENG       global IEN vector
!     IITYPE     interpolation type
!     IRPTR_LOC  local compressed row storage row pointer
!     JCPTR_LOC  loal compressed row storage column pointer
!     LISISO     True -> isotropic simulation
!     LM         location matrix for assembly
!     MDIM       max number of spatial dimensions
!     MEN        max number of element nodes
!     MGNOD      max number of anchor nodes
!     MNPG       max number of anchor noes in mesh 
!     MYDOFS     holds processes DOFS
!     MYELEM     holds processes elements
!     NDIM       number of spatial dimensions
!     NDOF       number of degrees of freedom
!     NDOFL      number of local degrees of freedom
!     NELEM      number of elements 
!     NELEML     number of local elements
!     NEN        number of element nodes
!     NGNOD      number of anchor nodes
!     NLETA      number of lagrange interpolant points in eta
!     NLXI       number of lagrange interpolant points in xi
!     NNPG       number of anchor nodes in mesh
!     XIPTS      xi interoplation points
!     XD         x distance into PML at anchor nodes
!     XLOCS      x locations of anchor nodes
!     XWIDTH     x width of PML
!     ZD         z distance into PML at anchor nodes
!     ZLOCS      z locations of anchor nodes
!     ZWIDTH     z width of PML
!
!     OUTPUT     MEANING
!     ------     ------- 
!     CDIST      distributed global impedance matrix
!     IERR       error flag
! 
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH
      CHARACTER(1), INTENT(IN) :: CFTYPE
      REAL*8, INTENT(IN) :: FREQ, PY 
      INTEGER*4, INTENT(IN) :: NZLOC
      COMPLEX*8, INTENT(OUT) :: CDIST(NZLOC) 
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      COMPLEX*16, ALLOCATABLE :: CE(:,:), DAMPX(:,:), DAMPZ(:,:)
      REAL*8, ALLOCATABLE :: SHG(:,:,:,:), SHL(:,:,:,:), DMAT(:,:,:),                    &
                             XIGLL(:,:),ETAGLL(:,:), RHO(:,:),DET(:,:), ME(:,:),         &
                             VP(:,:), VS(:,:), QPP(:,:), QSS(:,:) 
      REAL*8 XWIDTH, ZWIDTH, RCOEFF,DMAX_IN,RKAPPA,FCENT 
      INTEGER*4 MDIM,MEN,MGNOD,MNPG, NNPG,NLXI,NLETA,  &
                NELEM,NEN,NDOFL,NELEML, IITYPE,              &
                NINTX,NINTZ,IELEM,IELEML, NEE,NCOEFF,NORD 
      LOGICAL*4 LISISO 
      REAL*8 OMEGA, TWOPI
      LOGICAL*4 LCONV
      PARAMETER(TWOPI = 6.2831853071795862D0)
      PARAMETER(LCONV = .TRUE.) !convolutional PML
!
!----------------------------------------------------------------------------------------!
! 
!.... copy over sizes and PML info
      XWIDTH = msh%XWIDTH
      ZWIDTH = msh%ZWIDTH
      IF (CFTYPE == 'S') THEN
         RCOEFF  = msh%RCOEFF_SRF
         DMAX_IN = msh%DMAX_SRF
         RKAPPA  = msh%KAPPA_SRF
         FCENT   = msh%FCENT_SRF
      ELSE !body wave pml
         RCOEFF  = msh%RCOEFF_BDY
         DMAX_IN = msh%DMAX_BDY
         RKAPPA  = msh%KAPPA_BDY
         FCENT   = msh%FCENT_BDY
      ENDIF

      MDIM  = NDIM 
      MEN   = msh%NEN 
      MGNOD = NGNOD 
      MNPG  = msh%NNPG 
      NNPG  = msh%NNPG
      NLXI  = msh%NLXI 
      NLETA = msh%NLETA
      NELEM = msh%NELEM 
      NEN   = msh%NEN 
      NDOFL = msh%NDOFL 
      NELEML = msh%NELEML 
      NCOEFF = msh%NCOEFF
      NORD   = msh%NORD
      IITYPE = msh%IITYPE  

      LISISO = msh%LISISO 
!       
!.... define element shape functions
      IERR = 0 
      OMEGA = TWOPI*FREQ
      NINTX = 2*NORD + 2         !Lobatto quad, k = 2n - 3, k->2k for mass 
      NINTX = NINTX + 3          !add a little extra for deformed elements
      NINTZ = NINTX
      ALLOCATE(SHL(3,NEN,NINTX,NINTZ)) 
      ALLOCATE(XIGLL(NINTX,2))
      ALLOCATE(ETAGLL(NINTZ,2))
      CALL DSHL(NEN,NINTX,NINTZ, IITYPE,NINTX,NINTZ, NLXI,NLETA,    &   
                msh%XIPTS,msh%ETAPTS, XIGLL,ETAGLL,SHL)
!
!.... matrices for assembly
      ALLOCATE(SHG (3,NEN,NINTX,NINTZ))
      IF (.NOT.msh%LDISP) THEN
         ALLOCATE(DMAT(NINTX,NINTZ,NCOEFF))
         ALLOCATE(RHO (NINTX,NINTZ))
      ELSE
         ALLOCATE(VP  (NINTX,NINTZ))
         ALLOCATE(VS  (NINTX,NINTZ))
         ALLOCATE(QPP (NINTX,NINTZ))
         ALLOCATE(QSS (NINTX,NINTZ))
      ENDIF
      ALLOCATE(DET (NINTX,NINTZ))
      ALLOCATE(DAMPX(NINTX,NINTZ))
      ALLOCATE(DAMPZ(NINTX,NINTZ))
      NEE = NDIM*NEN
      ALLOCATE(CE(NEE,NEE))
      ALLOCATE(ME(NEE,NEE))
!.... null matrices
      CDIST(1:NZLOC) = CMPLX(0.0,0.0)
!
!.... loop on local elements
      DO 100 IELEML=1,NELEML
         IELEM = msh%MYELEM(IELEML)
         IF (msh%CDOMAIN(IELEM) /= 'A') THEN 
            IF (.NOT.msh%LDISP) THEN !standard
               CALL CJMK_TVI(NEN,NINTX,NINTZ,MNPG, LISISO,                       &
                             NINTX,NINTZ,NEN,NGNOD,NNPG,                         &
                             msh%IENG(1:NGNOD,IELEM),XIGLL(:,1),ETAGLL(:,1),     &
                             msh%XLOCS,msh%ZLOCS,msh%DENS,msh%ECOEFF,SHL,        &
                             DET,RHO,DMAT,SHG, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'asmble_driver: Error calling cjmk_tvi'
                  RETURN
               ENDIF 
               IF (LISISO) THEN
                  CALL KM25ISOQ(NEN,NEE,NINTX, NEN,NINTX,NINTZ,                          &
                              OMEGA,PY, XIGLL(:,2),ETAGLL(:,2), DMAT(:,:,1),DMAT(:,:,2), &
                              RHO,DET,SHG, ME,CE)
               ENDIF
            ELSE !medium is dispersive
               CALL CJKM_QUAL(NEN,NINTX,NNPG, LISISO,                                    &
                              NINTX,NINTZ,NEN,NGNOD,NNPG,                                &
                              msh%IENG(1:NGNOD,IELEM),XIGLL(:,1),ETAGLL(:,1),            & 
                              msh%XLOCS,msh%ZLOCS,msh%DENS,msh%QP,msh%QS,msh%ECOEFF,SHL, &
                              DET,VP,VS,QPP,QSS,SHG,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'asmble_driver: Error calling cjkm_qual'
                  RETURN
               ENDIF
               CALL KMDAMP(NEN,NEE,NINTX, NEN,NINTX,NINTZ,                      &
                           FREQ,msh%FREQ0,PY, XIGLL(:,2),ETAGLL(:,2), VP,VS,    &
                           QPP,QSS, DET,SHG, ME,CE)  
            ENDIF
         ELSE !absorbing domain
            IF (.NOT.msh%LDISP) THEN !standard
               CALL CJACCQ(NEN,NINTX,NINTZ,MNPG, LISISO,                                 &
                           LCONV,NINTX,NINTZ,NEN,NGNOD,NNPG,                             &
                           OMEGA,XWIDTH,ZWIDTH,                                          &
                           RCOEFF,DMAX_IN,RKAPPA,FCENT,                                  &
                           msh%IENG(1:NGNOD,IELEM), XIGLL(:,1),ETAGLL(:,1),              &
                           msh%XD,msh%ZD,msh%XLOCS,msh%ZLOCS,  msh%DENS,msh%ECOEFF,      &
                           SHL, DET,RHO,DMAT,SHG, DAMPX,DAMPZ,IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'asmble_driver: Error calling cjaccq'
                  RETURN
               ENDIF
               CALL CM25ISOQ(NEN,NEE,NINTX, NEN,NINTX,NINTZ, &
                             OMEGA,PY, XIGLL(:,2),ETAGLL(:,2), DMAT(:,:,1),DMAT(:,:,2),  &
                             RHO,DET, SHG, DAMPX,DAMPZ, ME,CE)
            ELSE !medium is dispersive
               CALL CJCM_QUAL(NEN,NINTX,NNPG, LISISO,                                    &
                              NINTX,NINTZ,NEN,NGNOD,NNPG,                                &
                              msh%IENG(1:NGNOD,IELEM), XIGLL(:,1),ETAGLL(:,1),           &
                              OMEGA,XWIDTH,ZWIDTH,                                       &
                              RCOEFF,DMAX_IN,RKAPPA,FCENT, msh%XD,msh%ZD,                &
                              msh%XLOCS,msh%ZLOCS,msh%DENS,msh%QP,msh%QS,msh%ECOEFF,SHL, &
                              DET,VP,VS,QPP,QSS,SHG, DAMPX,DAMPZ, IERR)
               IF (IERR /= 0) THEN
                  WRITE(*,*) 'asmble_driver: ERror calling cjaccq'
                  RETURN
               ENDIF
               CALL CMDAMP(NEN,NEE,NINTX, NEN,NINTX,NINTZ,                          &
                           FREQ,msh%FREQ0,PY,  XIGLL(:,2),ETAGLL(:,2), VP,VS,       &
                           QPP,QSS, DET,SHG, DAMPX,DAMPZ, ME,CE)
            ENDIF 
         ENDIF
         CALL ASMBLE_CDIST(NEE,NDIM,NZLOC,NDOFL, NEN,NDIM,                               &
                      msh%MYDOFS,msh%IRPTR_LOC,msh%JCPTR_LOC,msh%LM(1:NDIM,1:NEN,IELEM), &
                      ME,CE, CDIST,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'asmble_driver: Error assembling sparse matrix'
            RETURN
         ENDIF
  100 CONTINUE !loop on elements
!
!.... clean
      DEALLOCATE(SHL)
      DEALLOCATE(XIGLL)
      DEALLOCATE(ETAGLL)
      DEALLOCATE(SHG)
      IF (ALLOCATED(DMAT)) DEALLOCATE(DMAT)
      IF (ALLOCATED(RHO))  DEALLOCATE(RHO)
      IF (ALLOCATED(VP))   DEALLOCATE(VP)
      IF (ALLOCATED(VS))   DEALLOCATE(VS) 
      IF (ALLOCATED(QPP))  DEALLOCATE(QPP)
      IF (ALLOCATED(QSS))  DEALLOCATE(QSS) 
      DEALLOCATE(DET)
      DEALLOCATE(DAMPX)
      DEALLOCATE(DAMPZ) 
      DEALLOCATE(CE)
      DEALLOCATE(ME)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE ASMBLE_CDIST(MEE,MDIM,NZLOC,NDOFL, NEN,NDIM,  &
                              MYDOFS,IAOLOC,JAOLOC,LM,         &
                              ME,CE, CDIST,IERR)
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
!     CDIST      distributed global impedance matrix 
!     IERR       error flag
!  
!.... variable declarations
      COMPLEX*16, INTENT(IN) :: CE(MEE,*)
      REAL*8, INTENT(IN) :: ME(MEE,*)
      INTEGER*4, INTENT(IN) :: LM(MDIM,*), MYDOFS(NDOFL),                & 
                 IAOLOC(NDOFL+1), JAOLOC(NZLOC), MEE,MDIM,NZLOC,NDOFL,   & 
                 NEN,NDIM 
      COMPLEX*8, INTENT(INOUT) :: CDIST(NZLOC)
      INTEGER*4, INTENT(OUT) :: IERR
      REAL*8 R1, Q1
      COMPLEX*16 CARG
      INTEGER*4 IBSECT2,IBSECT, IAE,IDOF,IBE,JDOF,JBEG,JEND,NVAR,IZLOC,  &
                IPE,IQE,I,J,INDX 
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
            IF (IDOF.NE.0) THEN
               IDOFL = IBSECT2(NDOFL,IDOF,MYDOFS)
               IF (IDOFL.GT.0) THEN
                  JBEG = IAOLOC(IDOFL)
                  JEND = IAOLOC(IDOFL+1) - 1
                  NVAR = JEND - JBEG + 1
                  IQE = 0
                  DO 3 IBE=1,NEN
                     DO 4 J=1,NDIM
                        IQE = IQE + 1
                        JDOF = LM(J,IBE)
                        IF (JDOF.NE.0) THEN
                           INDX = IBSECT(NVAR,JDOF,JAOLOC(JBEG:JEND))
                           IF (INDX.LE.0) THEN
                              IERR = 1
                              WRITE(*,*) 'asmble_cdist: Error INDX < 0'
                              RETURN
                           ENDIF
                           IZLOC = JBEG - 1 + INDX
                           R1 = ME(IPE,IQE) + DREAL(CE(IPE,IQE))
                           Q1 =             + DIMAG(CE(IPE,IQE))
                           CARG = DCMPLX(R1,Q1)
                           CDIST(IZLOC) = CDIST(IZLOC) + CMPLX(CARG)
                        ENDIF !endif on jdof
    4                CONTINUE !loop on ndim 
    3             CONTINUE !loop on ibe
               ENDIF !end check if DOF belongs to this process
            ENDIF !endif on idof
    2    CONTINUE !loop on ndim
    1 CONTINUE !loop on iae 
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE ASMB_PEFF(MASTER,MYCOMM, NDOF,NDOFL,NZLOC, CNP,                       &
                           MYDOFS,IRPTR_LOC,JCPTR_LOC, CDIST,UE, PEFF)
!
!     INPUT      MEANING
!     -----      ------- 
!     CDIST      distributed impedance matrix
!     CNP        DOFs position in mesh
!     IRPTR_LOC  local CRS row pointer
!     JCPTR_LOC  local CRS column pointer
!     MASTER     master process ID
!     MYCOMM     MPI communicator
!     MYDOFS     global DOFs local to process
!     NDOF       number of degrees of freedom (global)
!     NDOFL      number of degrees of freedom local to process
!     NZLOC      number of non-zeros locally
!     UE         background wavefield in bielak layer
!
!     OUTPUT     MEANING
!     ------     ------- 
!     PEFF       an effective force distribution
!
!.... variable declarations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      CHARACTER(1), INTENT(IN) :: CNP(NDOF) 
      COMPLEX*8, INTENT(IN) :: UE(NDOF), CDIST(NZLOC)  
      INTEGER*4, INTENT(IN) :: MYDOFS(NDOFL), IRPTR_LOC(NDOFL+1), JCPTR_LOC(NZLOC),  &
                               NDOF,NDOFL,NZLOC, MASTER,MYCOMM  
      COMPLEX*8, INTENT(OUT) :: PEFF(*)
!.... local variables
      COMPLEX*8, ALLOCATABLE :: BUFF(:) 
      INTEGER*4 IDOFL, IDOF, JBEG, JEND, J, JDOF, MPIERR
      CHARACTER(2) CTYPE
! 
!========================================================================================!
!
!.... perform a distributed sparse matrix vector multiply 
      ALLOCATE(BUFF(NDOF)) 
      BUFF(1:NDOF) = CMPLX(0.0,0.0) 
      DO 1 IDOFL=1,NDOFL
         IDOF = MYDOFS(IDOFL) 
         CTYPE(1:1) = CNP(IDOF) 
         IF (CTYPE(1:1) == 'B' .OR. CTYPE(1:1) == 'E') THEN !check we may be on be bdry
            JBEG = IRPTR_LOC(IDOFL)
            JEND = IRPTR_LOC(IDOFL+1) - 1 
            DO 2 J=JBEG,JEND
               JDOF = JCPTR_LOC(J) 
               CTYPE(2:2) = CNP(JDOF)
               IF (CTYPE == 'BE') THEN
                  BUFF(IDOF) = BUFF(IDOF) - CDIST(J)*UE(JDOF)
               ELSEIF (CTYPE == 'EB') THEN
                  BUFF(IDOF) = BUFF(IDOF) + CDIST(J)*UE(JDOF) 
               ENDIF
    2       CONTINUE !Loop on row
         ENDIF !end check on location
    1 CONTINUE !loop on local DOFs
      CALL MPI_REDUCE(BUFF,PEFF,NDOF,MPI_COMPLEX,MPI_SUM, MASTER,MYCOMM,MPIERR)  
      DEALLOCATE(BUFF) 
      RETURN
      END 
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE VERIFY_STEVE(MYID,MASTER,MYCOMM, NDOF,NDOFL,NZLOC, CNP, &   
                              MYDOFS,IRPTR_LOC,JCPTR_LOC, CDIST,U0,SOL)
!
!     Steve's claim is P_e = S_{eb} u_b^0 + S_{ee} u_e^0 
!                          = S_{eb} u_e   + S_{ee} u_b 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      CHARACTER(1), INTENT(IN) :: CNP(NDOF)
      COMPLEX*8, INTENT(IN) :: SOL(*), U0(NDOF), CDIST(NZLOC)
      INTEGER*4, INTENT(IN) :: MYDOFS(NDOFL), IRPTR_LOC(NDOFL+1), JCPTR_LOC(NZLOC), & 
                               NDOF, NDOFL, NZLOC, MASTER, MYID, MYCOMM
!.... local variables
      COMPLEX*8, ALLOCATABLE :: BUFF1(:), BUFF2(:), UE(:)
      REAL*4 DIF, DIF1
      INTEGER*4 IDOF, IDOFL, JBEG, JEND, J, JDOF, MPIERR 
      CHARACTER(2) CTYPE 
!
!----------------------------------------------------------------------------------------!
!
      ALLOCATE(UE(NDOF))
      IF (MYID == MASTER) CALL CCOPY(NDOF,SOL,1,UE,1)
      CALL MPI_BCAST(UE,NDOF,MPI_COMPLEX, MASTER,MYCOMM,MPIERR)
      ALLOCATE(BUFF1(NDOF))
      ALLOCATE(BUFF2(NDOF))
      BUFF1(1:NDOF) = CMPLX(0.0,0.0)
      BUFF2(1:NDOF) = CMPLX(0.0,0.0)
      DO 1 IDOFL=1,NDOFL
         IDOF = MYDOFS(IDOFL)
         CTYPE(1:1) = CNP(IDOF)
         IF (CTYPE(1:1) == 'E') THEN 
            JBEG = IRPTR_LOC(IDOFL)
            JEND = IRPTR_LOC(IDOFL+1) - 1
            DO 2 J=JBEG,JEND
               JDOF = JCPTR_LOC(J)
               CTYPE(2:2) = CNP(JDOF)
               IF (CTYPE == 'EB') THEN !could be .OR. CTYPE == 'EE') THEN
                  BUFF1(IDOF) = BUFF1(IDOF) + CDIST(J)*U0(JDOF)  
                  BUFF2(IDOF) = BUFF2(IDOF) + CDIST(J)*UE(JDOF)
               ELSEIF (CTYPE == 'EE') THEN
                  BUFF1(IDOF) = BUFF1(IDOF) + CDIST(J)*U0(JDOF)  
                  BUFF2(IDOF) = BUFF2(IDOF) + CDIST(J)*UE(JDOF)
               ENDIF
    2       CONTINUE
         ENDIF
    1 CONTINUE
      DIF1 = MAXVAL(ABS(BUFF1 - BUFF2))
      CALL MPI_REDUCE(DIF1,DIF,1,MPI_REAL, MPI_MAX, MASTER,MYCOMM,MPIERR)
      IF (MYID == MASTER) WRITE(*,*) 'verify_steve: Max difference:',DIF1
      DEALLOCATE(BUFF1)
      DEALLOCATE(BUFF2)
      DEALLOCATE(UE)
      RETURN
      END
