!     Generates the graph for the forward problem
      SUBROUTINE GEN_GRAPH25(LVERB,NPARTS,  LFSURF,XREC,ZREC, MSH,RCV, IERR) 
!
!     Generates the graph for the 2.5D simulation
!
!.... variable declarations
      IMPLICIT NONE 
      REAL*8, INTENT(IN) :: XREC(*), ZREC(*)
      INTEGER*4, INTENT(IN) :: NPARTS
      LOGICAL*4, INTENT(IN) :: LFSURF(*), LVERB
      INTEGER*4, INTENT(OUT) :: IERR 
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH
      TYPE (RECV_INFO) RCV 
!.... local variables
      REAL*8 DZBASE 
      REAL*8, ALLOCATABLE :: XLOCSE(:), ZLOCSE(:)  
      INTEGER*4, ALLOCATABLE :: IEN(:,:), MCONN(:,:), ID(:,:), JCPTR(:) 
      INTEGER*4 NCONF, NNP, NWORK 
      INTEGER*4 ICMCON, IENZERO 
!
!----------------------------------------------------------------------------------------!
!
!.... generate the IEN matrix, ID matrix, LM matrix, then CRS pointers
      IF (LVERB) THEN
         WRITE(*,*) 
         WRITE(*,*) 'gen_graph25: Creating IEN matrix...'
      ENDIF 
      ALLOCATE(IEN(msh%NEN,msh%NELEM)) 
      CALL GENIEN(NGNOD,msh%NEN,msh%NNPG,msh%NELEM,msh%NLXI,msh%NLETA,NGNOD, LVERB,      &
                  msh%NEN,msh%IENG,msh%XIPTS,msh%ETAPTS,msh%XLOCS,msh%ZLOCS,             &
                  NNP,IEN,IERR) 
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_graph25: Error calling genien'
         RETURN
      ENDIF
      NCONF = ICMCON(msh%NEN,NNP, msh%NELEM,msh%NEN,IEN)
      ALLOCATE(MCONN(NNP,NCONF))
      CALL CONNECT(msh%NEN,NNP, NNP,msh%NELEM,msh%NEN,NCONF, IEN,MCONN)
      IF (LVERB) WRITE(*,*) 'gen_graph25: Generating ID vector...'
      ALLOCATE(ID(NDIM,NNP))
      CALL GENID(NDIM,NGNOD,msh%NEN, msh%NNPG,NNP,msh%NELEM,msh%NLXI,msh%NLETA,       &
                 NGNOD,NDIM, msh%CDOMAIN,msh%CNNPG, msh%IENG,IEN,                     &
                 msh%XIPTS,msh%ETAPTS,msh%XLOCS,msh%ZLOCS, msh%NDOF,ID,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_graph25: Error calling genid'
         RETURN 
      ENDIF
      IF (LVERB) WRITE(*,*) 'gen_graph25: Generating LM matrix...'
      ALLOCATE(msh%LM(NDIM,msh%NEN,msh%NELEM))
      CALL GENLM(NDIM,msh%NEN,msh%NEN,NNP, NDIM,msh%NELEM,ID,IEN, msh%LM)
      NWORK = IENZERO(NNP,NDIM, NNP,NDIM,NCONF,msh%NEN, ID,MCONN)
      DEALLOCATE(ID) 
      IF (LVERB) WRITE(*,*) 'gen_graph25: Generating compressed row storage format...'
      ALLOCATE(msh%IRPTR(msh%NDOF+1))
      ALLOCATE(JCPTR(NWORK))
      CALL GENCRS_QUICK(NDIM,msh%NEN,NNP,NWORK, NCONF,msh%NELEM,msh%NEN,NDIM,  &
                        msh%NDOF,MCONN,IEN,msh%LM, msh%NZERO,msh%IRPTR,JCPTR, IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_graph25: Error calling gencrs_quick'
         RETURN
      ENDIF
      DEALLOCATE(IEN)
      DEALLOCATE(MCONN)
      ALLOCATE(msh%JCPTR(msh%NZERO))
      msh%JCPTR(1:msh%NZERO) = JCPTR(1:msh%NZERO)
!
!.... classify DOFs 
      IF (LVERB) WRITE(*,*) 'gen_graph25: Locating Bielak locations...'
      ALLOCATE(msh%CNP(msh%NDOF))
      CALL CNMESH(NDIM,msh%NEN, msh%NDOF,msh%NELEM, NDIM,msh%NEN, msh%CDOMAIN,msh%LM,   &
                  msh%CNP,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_graph25: Error calling cnmesh'
         RETURN 
      ENDIF
!
!.... get number of GLL nodes in bielak layer and their (x,z) locations
      NWORK = MSH%NEN*MSH%NELEM
      ALLOCATE(XLOCSE(NWORK))
      ALLOCATE(ZLOCSE(NWORK))
      CALL GXZPTS_BLK2(NDIM,msh%NEN,NGNOD,msh%NDOF,msh%NNPG,msh%NELEM,msh%NLXI,msh%NLETA,&
                       NWORK, NGNOD,NDIM, msh%CDOMAIN,msh%CNP, msh%IENG,msh%LM,          &
                       msh%XIPTS,msh%ETAPTS, msh%XLOCS,msh%ZLOCS,                        &
                       msh%NNPE,XLOCSE,ZLOCSE,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_graph25: Error calling gxzpts_blk2'
         RETURN
      ENDIF
      ALLOCATE(msh%XLOCSE(msh%NNPE))
      ALLOCATE(msh%ZLOCSE(msh%NNPE))
      msh%XLOCSE(1:msh%NNPE) = XLOCSE(1:msh%NNPE)
      msh%ZLOCSE(1:msh%NNPE) = ZLOCSE(1:msh%NNPE) 
      DEALLOCATE(XLOCSE)
      DEALLOCATE(ZLOCSE)
!
!.... connect Bielak points to DOFS and determine nodal point type
      IF (LVERB) WRITE(*,*) 'gen_graph25: Connecting Bielak points to DOF numbers...'
      ALLOCATE(msh%IDOFSE(NDIM,msh%NNPE))
      CALL CDOF2XZE2(NDIM,msh%NEN,NGNOD, msh%NDOF,msh%NNPG,msh%NELEM,msh%NLXI,msh%NLETA, &
                     msh%NNPE, NGNOD,NDIM, msh%CDOMAIN,msh%CNP, msh%IENG,msh%LM,         &
                     msh%XIPTS,msh%ETAPTS, msh%XLOCS,msh%ZLOCS, msh%XLOCSE,msh%ZLOCSE,   &
                     msh%IDOFSE,IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'gen_graph25: Error calling cdof2xze'
         RETURN 
      ENDIF
      msh%ZBASE_INT = DZBASE(msh%NNPG,1,msh%CNNPG,msh%ZLOCS) 
      msh%ZBASE_BLK = DZBASE(msh%NNPG,0,msh%CNNPG,msh%ZLOCS)
!
!     useful for checking mesh
!     do ielem=1,nelem
!        iae = 0
!        do ileta=1,nleta
!           do ilxi=1,nlxi
!              iae = iae + 1
!              i = 1
!              idof = lm(i,iae,ielem)
!              if (idof > 0) then
!                 call csf2dnd(ngnod,xipts(ilxi),etapts(ileta), sf, ierr) 
!                 x = 0.d0
!                 z = 0.d0
!                 do ia=1,ngnod
!                    iloc = ieng(ia,ielem)
!                    x = x + sf(ia)*xlocs(iloc)
!                    z = z + sf(ia)*zlocs(iloc)
!                 enddo
!                 if (cnp(idof) == 'A') write(70,*) x,z
!                 if (cnp(idof) == 'E') write(71,*) x,z
!                 if (cnp(idof) == 'B') write(72,*) x,z
!                 if (cnp(idof) == 'I') write(73,*) x,z  
!              else
!                 call csf2dnd(ngnod,xipts(ilxi),etapts(ileta), sf, ierr)
!                 x = 0.d0
!                 z = 0.d0
!                 do ia=1,ngnod
!                    iloc = ieng(ia,ielem)
!                    x = x + sf(ia)*xlocs(iloc)
!                    z = z + sf(ia)*zlocs(iloc)
!                 enddo
!                 write(74,*) x,z
!              endif
!           enddo
!       enddo
!       iloc1 = ieng(1,ielem)
!       if (cdomain(ielem) == 'I') then
!          do ia=1,ngnod
!             iloc = ieng(ia,ielem)
!             write(80,*) xlocs(iloc),zlocs(iloc)
!          enddo
!          write(80,*) xlocs(iloc1),zlocs(iloc1)
!          write(80,*)
!       elseif (cdomain(ielem) == 'A') then
!          do ia=1,ngnod
!             iloc = ieng(ia,ielem)
!             write(81,*) xlocs(iloc),zlocs(iloc)
!          enddo
!          write(81,*) xlocs(iloc1),zlocs(iloc1)
!          write(81,*)
!       elseif (cdomain(ielem) == 'E') then
!          do ia=1,ngnod
!             iloc = ieng(ia,ielem)
!             write(82,*) xlocs(iloc),zlocs(iloc)
!          enddo
!          iloc1 = ieng(1,ielem)
!          write(82,*) xlocs(iloc1),zlocs(iloc1) 
!          write(82,*)
!       else
!         write(*,*) 'fuck',cdomain(ielem)
!       endif
!     enddo
!     stop
!
!.... color the mesh
      IF (LVERB) WRITE(*,*) 'gen_graph25: Coloring mesh...'
      ALLOCATE(msh%PART(msh%NDOF))
      CALL COLOR_MESH(NDIM,msh%NEN,msh%NDOF,msh%NZERO, NDIM,msh%NEN,msh%NELEM,        &
                      NPARTS, .FALSE.,msh%CDOMAIN,msh%IRPTR,msh%JCPTR,msh%LM, msh%PART)
!
!.... if need be locate the receivers
      IF (rcv%NREC > 0) THEN
         IF (LVERB) WRITE(*,*) 'gen_graph25: Locating receivers...'
         ALLOCATE(rcv%MRDOF(NDIM,rcv%NREC))
         CALL LOCREC2D(NDIM,msh%NEN,NGNOD, msh%NNPG,rcv%NREC,msh%NLXI,msh%NLETA,       &
                       LVERB, msh%NELEM,NGNOD,NDIM, msh%CNNPG,LFSURF,msh%IENG,msh%LM,  &
                       msh%XIPTS,msh%ETAPTS,XREC,ZREC,msh%XLOCS,msh%ZLOCS,             &
                       rcv%MRDOF,IERR)
         IF (IERR /= 0) THEN
            WRITE(*,*) 'gen_graph25: Error calling locrec2d'
            RETURN
         ENDIF
      ENDIF !end check on locating receivers
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE GEN_GIDOFFS(MSH,IERR) 
!
!     For inversion we may want to use the unwrapped phase estimates at the free 
!     surface.  This routine extracts the DOFs at the free surface so we can unwrap
!     the phase and collocate to a receiver
!
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH
      INTEGER*4, INTENT(OUT) :: IERR
      INTEGER*4, ALLOCATABLE :: IDOF_FS(:,:)
      INTEGER*4 I  
!
!----------------------------------------------------------------------------------------!
!
      ALLOCATE(IDOF_FS(NDIM,msh%NLXI*msh%NNPG)) !more than enough spcae
      CALL GIDOFFS(NDIM,msh%NEN,NGNOD, msh%NNPG,msh%NLXI,msh%NLETA,                 &
                   msh%NELEM,NGNOD,NDIM, msh%CNNPG,msh%IENG,msh%LM,                 &
                   msh%XIPTS,msh%ETAPTS,msh%XLOCS,msh%ZLOCS, msh%NFS,IDOF_FS,IERR) 
      IF (IERR /= 0) THEN 
         WRITE(*,*) 'gen_idoffs: Error generating free surface DOF numbers!'
         DEALLOCATE(IDOF_FS)
         RETURN
      ENDIF
      ALLOCATE(msh%IDOF_FS(NDIM,msh%NFS)) 
      DO 1 I=1,NDIM
         msh%IDOF_FS(I,1:msh%NFS) = IDOF_FS(I,1:msh%NFS)
    1 CONTINUE
      DEALLOCATE(IDOF_FS) 
      RETURN
      END
