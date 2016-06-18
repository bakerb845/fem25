      SUBROUTINE BCAST_MESH_INFO(MYID,MYCOMM,MASTER, 
     ;                           LNSRF, PROJNM,TMPDIR, MSH) 
!
!     Broadcasts forward problem.  Memory handling is performed here 
!     to clarify the driver program.
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      CHARACTER(80) PROJNM, TMPDIR 
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER 
      LOGICAL*4 LNSRF
      TYPE (MESH_INFO) MSH
!.... local variables
      INTEGER*4 IA, IAE, I, ICOEFF, INPG, MPIERR
!
!----------------------------------------------------------------------!
!
!.... model seizes
      CALL MPI_BCAST(PROJNM,   80,MPI_CHARACTER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(TMPDIR,   80,MPI_CHARACTER, MASTER,MYCOMM,MPIERR)

      CALL MPI_BCAST(msh%LISISO,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%LDISP ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR) 
      CALL MPI_BCAST(LNSRF     ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
 
      CALL MPI_BCAST(msh%NELEM ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NDOF  ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NZERO ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NNPG  ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NNPE  ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NORD  ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NLXI  ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NLETA ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%NEN   ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%IITYPE,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%NCOEFF,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR) 
!
!.... model REAL*8 parameters 
      CALL MPI_BCAST(msh%XWIDTH   ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%ZWIDTH   ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%FCENT_SRF ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%FCENT_BDY ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%KAPPA_SRF ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%KAPPA_BDY ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%RCOEFF_SRF,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%RCOEFF_BDY,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%DMAX_SRF  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%DMAX_BDY  ,1,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XMOD0    ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%XMOD1    ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%XBLKL    ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XBLKR    ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XLATMIN  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XLATMAX  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XLONMIN  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%XLONMAX  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%AZMOD    ,1,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%AZTOL    ,1,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%AOITOL   ,1,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%FREQ0    ,1,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%ZBASE_INT,1,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%ZBASE_BLK,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
!
!.... allocate forward problem
      IF (MYID.NE.MASTER) THEN
         ALLOCATE(msh%CDOMAIN(msh%NELEM))
         ALLOCATE(msh%CNP(msh%NDOF))
         ALLOCATE(msh%CNNPG(msh%NNPG))

         ALLOCATE(msh%LM(NDIM,msh%NEN,msh%NELEM))
         ALLOCATE(msh%IENG(NGNOD,msh%NELEM))
         ALLOCATE(msh%IDOFSE(NDIM,msh%NNPE)) 
         ALLOCATE(msh%IRPTR(msh%NDOF+1))
         ALLOCATE(msh%JCPTR(msh%NZERO))
         ALLOCATE(msh%PART(msh%NDOF))

         ALLOCATE(msh%ECOEFF(msh%NNPG,msh%NCOEFF))
         ALLOCATE(msh%QP    (msh%NNPG))
         ALLOCATE(msh%QS    (msh%NNPG))
         ALLOCATE(msh%DENS  (msh%NNPG))
         ALLOCATE(msh%XD    (msh%NNPG))
         ALLOCATE(msh%ZD    (msh%NNPG))
         ALLOCATE(msh%XLOCS (msh%NNPG))
         ALLOCATE(msh%ZLOCS (msh%NNPG))
         ALLOCATE(msh%XLOCSE(msh%NNPE))
         ALLOCATE(msh%ZLOCSE(msh%NNPE))
         ALLOCATE(msh%XIPTS (msh%NLXI))
         ALLOCATE(msh%ETAPTS(msh%NLETA)) 
      ENDIF
!
!.... characters first
      CALL MPI_BCAST(msh%CDOMAIN, msh%NELEM, MPI_CHARACTER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%CNP,      msh%NDOF, MPI_CHARACTER, MASTER,
     ;               MYCOMM,MPIERR)
      DO 1 INPG=1,msh%NNPG
         CALL MPI_BCAST(msh%CNNPG(INPG),2, MPI_CHARACTER,
     ;                  MASTER,MYCOMM,MPIERR)
    1 CONTINUE
!
!.... integers
      CALL MPI_BCAST(msh%PART, msh%NDOF  ,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%IRPTR,msh%NDOF+1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%JCPTR,msh%NZERO ,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      DO 2 I=1,NDIM
         CALL MPI_BCAST(msh%IDOFSE(I,1:msh%NNPE),msh%NNPE ,MPI_INTEGER, 
     ;                  MASTER,MYCOMM,MPIERR)
    2 CONTINUE
      DO 3 IA=1,NGNOD
         CALL MPI_BCAST(msh%IENG(IA,1:msh%NELEM),msh%NELEM,MPI_INTEGER, 
     ;                  MASTER,MYCOMM,MPIERR)
    3 CONTINUE
      DO 4 IAE=1,msh%NEN
         DO 5 I=1,NDIM
            CALL MPI_BCAST(msh%LM(I,IAE,1:msh%NELEM),msh%NELEM, 
     ;                     MPI_INTEGER, MASTER,MYCOMM,MPIERR)
    5    CONTINUE
    4 CONTINUE

      CALL MPI_BCAST(msh%XD   ,msh%NNPG,  MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%ZD   ,msh%NNPG,  MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XLOCS,msh%NNPG,  MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%ZLOCS,msh%NNPG,  MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%DENS ,msh%NNPG,  MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%QP   ,msh%NNPG,  MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%QS   ,msh%NNPG,  MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)

      CALL MPI_BCAST(msh%XLOCSE,msh%NNPE, MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%ZLOCSE,msh%NNPE, MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR) 

      CALL MPI_BCAST(msh%XIPTS ,msh%NLXI ,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%ETAPTS,msh%NLETA,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      DO 6 ICOEFF=1,msh%NCOEFF
         CALL MPI_BCAST(msh%ECOEFF(1:msh%NNPG,ICOEFF),msh%NNPG,
     ;                  MPI_DOUBLE_PRECISION, MASTER,MYCOMM,MPIERR)
    6 CONTINUE 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_FS_INFO(MYID,MYCOMM,MASTER, MSH) 
!     Broadcasts info on the free surface
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER
      TYPE (MESH_INFO) MSH
      INTEGER*4 I, MPIERR
!----------------------------------------------------------------------!
      CALL MPI_BCAST(msh%NFS,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      IF (MYID /= MASTER) ALLOCATE(msh%IDOF_FS(NDIM,msh%NFS))
      DO 1 I=1,NDIM
        CALL MPI_BCAST(msh%IDOF_FS(I,1:msh%NFS),msh%NFS,MPI_INTEGER, 
     ;                 MASTER,MYCOMM,MPIERR)
    1 CONTINUE 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_MOD1D_INFO(MYID,MYCOMM,MASTER, MSH,M1D)
!     Broadcasts the 1D models.  Also, since xmod0,xmod1, xblkl, and 
!     xlkbr are calculated with the 1d models we remove the potential
!     for miscalling sequences by broadcasting here
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER 
      TYPE (MESH_INFO) MSH
      TYPE (MOD1D_INFO) M1D  
      INTEGER*4 MPIERR 
!----------------------------------------------------------------------!
      CALL MPI_BCAST(m1d%NL1D_LT,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%NL1D_RT,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XMOD0  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%XMOD1  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(msh%XBLKL  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(msh%XBLKR  ,1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      IF (MYID.NE.MASTER) THEN
         ALLOCATE(m1d%VP1D_LT (m1d%NL1D_LT))
         ALLOCATE(m1d%VS1D_LT (m1d%NL1D_LT))
         ALLOCATE(m1d%RH1D_LT (m1d%NL1D_LT))
         ALLOCATE(m1d%QP1D_LT (m1d%NL1D_LT))
         ALLOCATE(m1d%QS1D_LT (m1d%NL1D_LT))
         ALLOCATE(m1d% Z1D_LT (m1d%NL1D_LT))
         ALLOCATE(m1d%VPD_RLLT(m1d%NL1D_LT))
         ALLOCATE(m1d%VSD_RLLT(m1d%NL1D_LT))
         ALLOCATE(m1d%ROD_RLLT(m1d%NL1D_LT))
         ALLOCATE(m1d%HDD_RLLT(m1d%NL1D_LT))
         ALLOCATE(m1d%VPD_LVLT(m1d%NL1D_LT))
         ALLOCATE(m1d%VSD_LVLT(m1d%NL1D_LT))
         ALLOCATE(m1d%ROD_LVLT(m1d%NL1D_LT))
         ALLOCATE(m1d%HDD_LVLT(m1d%NL1D_LT))

         ALLOCATE(m1d%VP1D_RT (m1d%NL1D_RT))
         ALLOCATE(m1d%VS1D_RT (m1d%NL1D_RT))
         ALLOCATE(m1d%RH1D_RT (m1d%NL1D_RT))
         ALLOCATE(m1d%QP1D_RT (m1d%NL1D_RT))
         ALLOCATE(m1d%QS1D_RT (m1d%NL1D_RT))
         ALLOCATE(m1d% Z1D_RT (m1d%NL1D_RT))
         ALLOCATE(m1d%VPD_RLRT(m1d%NL1D_RT))
         ALLOCATE(m1d%VSD_RLRT(m1d%NL1D_RT))
         ALLOCATE(m1d%ROD_RLRT(m1d%NL1D_RT))
         ALLOCATE(m1d%HDD_RLRT(m1d%NL1D_RT))
         ALLOCATE(m1d%VPD_LVRT(m1d%NL1D_RT))
         ALLOCATE(m1d%VSD_LVRT(m1d%NL1D_RT))
         ALLOCATE(m1d%ROD_LVRT(m1d%NL1D_RT))
         ALLOCATE(m1d%HDD_LVRT(m1d%NL1D_RT))
      ENDIF
!
!.... left
      CALL MPI_BCAST(m1d%VP1D_LT,  m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%VS1D_LT,  m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%RH1D_LT,  m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%QP1D_LT,  m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%QS1D_LT,  m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d% Z1D_LT,  m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
!.... left rayleigh
      CALL MPI_BCAST(m1d%VPD_RLLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%VSD_RLLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%ROD_RLLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%HDD_RLLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
!.... left love
      CALL MPI_BCAST(m1d%VPD_LVLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%VSD_LVLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%ROD_LVLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%HDD_LVLT,m1d%NL1D_LT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
!
!.... right
      CALL MPI_BCAST(m1d%VP1D_RT,  m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%VS1D_RT,  m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%RH1D_RT,  m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%QP1D_RT,  m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR) 
      CALL MPI_BCAST(m1d%QS1D_RT,  m1d%NL1D_RT,MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR) 
      CALL MPI_BCAST(m1d% Z1D_RT,  m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
!.... right rayleigh
      CALL MPI_BCAST(m1d%VPD_RLRT,m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%VSD_RLRT,m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%ROD_RLRT, m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%HDD_RLRT,m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
!.... right love
      CALL MPI_BCAST(m1d%VPD_LVRT,m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%VSD_LVRT,m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%ROD_LVRT,m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(m1d%HDD_LVRT,m1d%NL1D_RT,MPI_DOUBLE_PRECISION,
     ;               MASTER,MYCOMM,MPIERR)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_SRC_INFO(MYID,MYCOMM,MASTER, SRC)
!     Broadcasts source information
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (SRC_INFO) SRC 
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER
      INTEGER*4 ISRC, MPIERR 
!----------------------------------------------------------------------!
      CALL MPI_BCAST(src%NSRC    ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(src%NSRC_SRF,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(src%NSRC_BDY,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR) 
      IF (MYID.NE.MASTER) THEN
         ALLOCATE(src%SRCTYP(src%NSRC))
         ALLOCATE(src%CSIDE(src%NSRC))

         ALLOCATE(src%MODE(src%NSRC))

         ALLOCATE(src%AOI(src%NSRC))
         ALLOCATE(src%BAZN(src%NSRC))
         ALLOCATE(src%SLAT(src%NSRC))
         ALLOCATE(src%SLON(src%NSRC))
         ALLOCATE(src%SDEP(src%NSRC))
         ALLOCATE(src%STRIKE(src%NSRC))
         ALLOCATE(src%DIP(src%NSRC))
         ALLOCATE(src%RAKE(src%NSRC))
         ALLOCATE(src%SMAG(src%NSRC))
      ENDIF
!.... characters
      DO 1 ISRC=1,src%NSRC 
         CALL MPI_BCAST(src%CSIDE(ISRC),1,MPI_CHARACTER,MASTER, 
     ;                  MYCOMM,MPIERR)
    1 CONTINUE
      DO 2 ISRC=1,src%NSRC
         CALL MPI_BCAST(src%SRCTYP(ISRC),2,MPI_CHARACTER,
     ;                  MASTER,MYCOMM,MPIERR)
    2 CONTINUE
!.... integers
      CALL MPI_BCAST(src%MODE,src%NSRC,MPI_INTEGER,MASTER, 
     ;               MYCOMM,MPIERR)
!.... real*8
      CALL MPI_BCAST(src%AOI   ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%BAZN  ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%SLAT  ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%SLON  ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(src%SDEP  ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%STRIKE,src%NSRC,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(src%DIP   ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%RAKE  ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%SMAG  ,src%NSRC,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)

      RETURN
      END
!                                                                      !

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_SRC_STF(MYID,MYCOMM,MASTER, LALLOC,NFREQ, SRC) 
!     Broadcasts the source time functions 
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER, NFREQ
      LOGICAL*4, INTENT(IN) :: LALLOC
      TYPE (SRC_INFO) SRC 
      INTEGER*4 IFREQ, MPIERR
!----------------------------------------------------------------------!
      IF (MYID.NE.MASTER .AND. LALLOC) 
     ;ALLOCATE(src%SOURCE(NFREQ,src%NSRC)) 
      DO 1 IFREQ=1,NFREQ
         CALL MPI_BCAST(SRC%SOURCE(IFREQ,1:src%NSRC),src%NSRC, 
     ;                  MPI_COMPLEX,MASTER,MYCOMM, MPIERR) 
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_SRC_PTRS(MYID,MYCOMM,MASTER, SRC) 
!     Broadcasts the source pointers
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER
      TYPE (SRC_INFO) SRC
      INTEGER*4 MPIERR
!----------------------------------------------------------------------!
      CALL MPI_BCAST(src%NSG,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)  
      IF (src%NSG <= 0) RETURN !group may have no sources
      IF (MYID.NE.MASTER) THEN
         ALLOCATE(src%ISGPTR(src%NSG+1))
         ALLOCATE(src%ISRCPRM(src%NSRC))
         ALLOCATE(src%PYAVG(src%NSG))
      ENDIF
      CALL MPI_BCAST(src%ISGPTR,src%NSG+1,MPI_INTEGER, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%ISRCPRM,src%NSRC,MPI_INTEGER, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(src%PYAVG,src%NSG,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_RCV_INFO(MYID,MYCOMM,MASTER, RCV)
!     Broadcasts the receiver information
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER
      TYPE (RECV_INFO) RCV
      INTEGER*4 I, MPIERR
!----------------------------------------------------------------------!
      CALL MPI_BCAST(rcv%NREC,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      IF (rcv%NREC <= 1) RETURN !don't necessarily need recievers
      IF (MYID.NE.MASTER) THEN
         ALLOCATE(rcv%YREC(rcv%NREC))
         ALLOCATE(rcv%MRDOF(NDIM,rcv%NREC))
      ENDIF
      DO 1 I=1,NDIM
         CALL MPI_BCAST(rcv%MRDOF(I,1:rcv%NREC),rcv%NREC,MPI_INTEGER, 
     ;                  MASTER,MYCOMM,MPIERR)
    1 CONTINUE 
      CALL MPI_BCAST(rcv%YREC,rcv%NREC,MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_RCV_RESP(MYID,MYCOMM,MASTER, LALLOC,NFREQ, RCV) 
!     Broadcasts the frequency domain receiver response functions
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER, NFREQ
      LOGICAL*4, INTENT(IN) :: LALLOC
      TYPE (RECV_INFO) RCV
      INTEGER*4 IFREQ, I, MPIERR
!----------------------------------------------------------------------!
      IF (MYID.NE.MASTER .AND. LALLOC) 
     ;ALLOCATE(rcv%RECV(NDIM,NFREQ,rcv%NREC)) 
      DO 1 IFREQ=1,NFREQ
         DO 2 I=1,NDIM
            CALL MPI_BCAST(rcv%RECV(I,IFREQ,1:rcv%NREC),rcv%NREC, 
     ;                     MPI_COMPLEX, MASTER,MYCOMM,MPIERR) 
    2    CONTINUE
    1 CONTINUE
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_FRQ_INFO(MYID,MYCOMM,MASTER, FRQ) 
!     Broadcast the frequencies to model
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h' 
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER
      TYPE (FRQ_INFO) FRQ 
      INTEGER*4 I, MPIERR
!
!----------------------------------------------------------------------!
!
      CALL MPI_BCAST(frq%NFREQ,    1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_SRF,1,MPI_INTEGER, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_BDY,1,MPI_INTEGER, MASTER, 
     ;               MYCOMM,MPIERR) 
      IF (MYID.NE.MASTER) THEN
         ALLOCATE(frq%CFTYPE(frq%NFREQ))
         ALLOCATE(frq%FREQ(frq%NFREQ))
      ENDIF
      DO 1 I=1,frq%NFREQ
         CALL MPI_BCAST(frq%CFTYPE(I),1,MPI_CHARACTER, MASTER, 
     ;                  MYCOMM,MPIERR)
    1 CONTINUE 
      CALL MPI_BCAST(frq%FREQ  ,frq%NFREQ,MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR) 
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_FRQ_JOINT_INFO(MYID,MYCOMM,MASTER, FRQ) 
!     Broadcast the frequencies to model/invert
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h' 
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER
      TYPE (FRQ_INFO) FRQ
      INTEGER*4 I, MPIERR
!
!----------------------------------------------------------------------!
!
      CALL MPI_BCAST(frq%NFREQ,        1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_SRF,    1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_BDY,    1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_SRF_INV,1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_BDY_INV,1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_INV    ,1,MPI_INTEGER, MASTER, 
     ;               MYCOMM,MPIERR)
      IF (MYID.NE.MASTER) THEN
         ALLOCATE(frq%CFTYPE(frq%NFREQ))
         ALLOCATE(frq%LINVF(frq%NFREQ))
         ALLOCATE(frq%FREQ(frq%NFREQ))
      ENDIF
      DO 1 I=1,frq%NFREQ
         CALL MPI_BCAST(frq%CFTYPE(I),1,MPI_CHARACTER, 
     ;                  MASTER,MYCOMM,MPIERR)
    1 CONTINUE 
      CALL MPI_BCAST(frq%LINVF,frq%NFREQ,MPI_LOGICAL, 
     ;                  MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%FREQ ,frq%NFREQ,MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR) 
      RETURN
      END
