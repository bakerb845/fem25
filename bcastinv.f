      SUBROUTINE BCAST_INV_PARMS(MYID,MYCOMM,MASTER, NNPG,NELEM,NBLOCKS,
     ;                           NREC, ISTOP_PT, INV)
!     Broadcasts global inversion parameters
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV
      INTEGER*4, INTENT(IN) :: MYID,MYCOMM,MASTER, NNPG, NELEM, NREC 
      INTEGER*4 NBLOCKS, ISTOP_PT
      INTEGER*4 ICON,MPIERR

      CALL MPI_BCAST(inv%CINVTYPE,2,MPI_CHARACTER, MASTER,MYCOMM,MPIERR)

      CALL MPI_BCAST(ISTOP_PT   ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NBLOCKS    ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%MCJLOC ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NHSIZE ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NVINV  ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NNPINV ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NA35   ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NCON   ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      
      CALL MPI_BCAST(inv%MAXIT   ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NORM    ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%IRESTP  ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%IMODSRC ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%IMODREC ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%IBPHASE ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NCASC   ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%NREC_MIN,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)

      CALL MPI_BCAST(inv%PTHRESH,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)

      CALL MPI_BCAST(inv%LMIGR     ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%LUNWRAP   ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%LBODY     ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%LSURF     ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%LAHESS_HOT,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%LCOVD_SRF ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%LCOVD_BDY ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%LDWGHT    ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)

      CALL MPI_BCAST(inv%DTMAX_SRF,1,MPI_REAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%DTMAX_BDY,1,MPI_REAL, MASTER,MYCOMM,MPIERR)

      CALL MPI_BCAST(inv%SCLBDY,1,MPI_DOUBLE_PRECISION,MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%VPVS  ,1,MPI_DOUBLE_PRECISION,MASTER,
     ;               MYCOMM,MPIERR)
      IF (MYID /= MASTER) THEN
         ALLOCATE(inv%MASKG(NNPG))
         ALLOCATE(inv%IGPART(inv%NNPINV)) 
         ALLOCATE(inv%MCONN(NNPG,inv%NCON))
         ALLOCATE(inv%WMASK(NNPG))
         ALLOCATE(inv%ELEM_WTS(NELEM))
         ALLOCATE(inv%DX(NREC-1))
      ENDIF 
      CALL MPI_BCAST(inv%MASKG,       NNPG,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%IGPART,inv%NNPINV,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      DO 1 ICON=1,inv%NCON
         CALL MPI_BCAST(inv%MCONN(1:NNPG,ICON),NNPG,MPI_INTEGER, 
     ;                  MASTER,MYCOMM,MPIERR)
    1 CONTINUE
      CALL MPI_BCAST(inv%WMASK,NNPG,MPI_REAL   , MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%ELEM_WTS,NELEM,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(inv%DX,NREC-1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_LINVF(MYID,MYCOMM,MASTER, FRQ)
!
!     Broadcasts the inversion frequency list and number of inversion
!     frequencies 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (FRQ_INFO) FRQ
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER
      INTEGER*4 MPIERR

      CALL MPI_BCAST(frq%NFREQ_INV    ,1,MPI_INTEGER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_SRF_INV,1,MPI_INTEGER, MASTER, 
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(frq%NFREQ_BDY_INV,1,MPI_INTEGER, MASTER, 
     ;               MYCOMM,MPIERR)
      IF (MYID.NE.MASTER) ALLOCATE(frq%LINVF(frq%NFREQ))
      CALL MPI_BCAST(frq%NFREQ,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_OBS_INFO(MYID,MYCOMM,MASTER, NFREQ,NREC,NSRC, 
     ;                          INV) 
!     Broadcasts observations
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV 
      INTEGER*4, INTENT(IN) :: MYID,MYCOMM,MASTER, NFREQ, NREC, NSRC 
      INTEGER*4 IFREQ,ISRC,I,MPIERR

      IF (MYID.NE.MASTER) THEN
         ALLOCATE(inv%OBS(NDIM,NFREQ,NREC,NSRC))  
         ALLOCATE(inv%WGHTS(NDIM,NFREQ,NREC,NSRC))
      ELSE
         WRITE(*,*) 'bcast_obs: Broadcasting observations...'
      ENDIF
      DO 1 IFREQ=1,NFREQ
         DO 2 ISRC=1,NSRC
            DO 3 I=1,NDIM
               CALL MPI_BCAST(inv%  OBS(I,IFREQ,1:NREC,ISRC),NREC, 
     ;                        MPI_COMPLEX,MASTER, MYCOMM,MPIERR)
               CALL MPI_BCAST(inv%WGHTS(I,IFREQ,1:NREC,ISRC),NREC,
     ;                        MPI_REAL,MASTER, MYCOMM,MPIERR)
    3       CONTINUE !loop on components
    2    CONTINUE
    1 CONTINUE !loop on frequencies
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_MOD_UPD(MYCOMM,MASTER, MSH) 
!     Just broadcasts the material updates
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH
      INTEGER*4, INTENT(IN) :: MYCOMM, MASTER
      INTEGER*4 ICOEFF, MPIERR

      CALL MPI_BCAST(msh%DENS ,msh%NNPG,  MPI_DOUBLE_PRECISION, 
     ;               MASTER,MYCOMM,MPIERR)
      DO 1 ICOEFF=1,msh%NCOEFF
         CALL MPI_BCAST(msh%ECOEFF(1:msh%NNPG,ICOEFF),msh%NNPG,
     ;                  MPI_DOUBLE_PRECISION, MASTER,MYCOMM,MPIERR)
    1 CONTINUE

      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_WIN(MYID,MYCOMM,MASTER, MREC,MSRC,NREC,NSRC, WIN)
!     Broadcasts the window information
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'fwd_struc.h'
      TYPE (WIN_INFO) WIN 
      INTEGER*4, INTENT(IN) :: MYID, MYCOMM, MASTER, MREC,MSRC,NREC,NSRC
      INTEGER*4 IREC, I, MPIERR
!
!----------------------------------------------------------------------!
!
      CALL MPI_BCAST(win%IWNDO    ,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(win%NSAMP_SRF,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(win%NSAMP_BDY,1,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(win%LWNDO_SRF,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(win%LWNDO_BDY,1,MPI_LOGICAL, MASTER,MYCOMM,MPIERR) 
      CALL MPI_BCAST(win%DT_SRF,   1,MPI_REAL,MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(win%DT_BDY,   1,MPI_REAL,MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(win%START_SRF,1,MPI_REAL,MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(win%START_BDY,1,MPI_REAL,MASTER,MYCOMM,MPIERR)
      IF (MYID.NE.MASTER) ALLOCATE(win%TWIN(MREC,MSRC,2))
      DO 1 IREC=1,NREC
         DO 2 I=1,2 
            CALL MPI_BCAST(win%TWIN(IREC,1:NSRC,I),NSRC, MPI_REAL,
     ;                     MASTER,MYCOMM,MPIERR)
    2    CONTINUE
    1 CONTINUE
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      ! 
      SUBROUTINE BCAST_INV_P1(MYCOMM,MASTER, PROJNM,CINVTYPE, IITYPE, 
     ;                        LISISO,LPHESS,LGNEWT, 
     ;                        NPARTS,NNPINV,NVINV,NA35,NHSIZE,
     ;                        MAXIT,NSRC,NSG,NREC, NORD,NLXI,NLETA,NEN, 
     ;                        NNPG,NNPE,NELEM,NBLOCKS, NL1D_LT,NL1D_RT, 
     ;                        NDOF,NZERO,NCON) 
!     First part of inverse problem broadcasts.  Allows us to 
!     determine sizes for vectors
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      CHARACTER(80) PROJNM
      CHARACTER(2)  CINVTYPE
      INTEGER*4, INTENT(IN) :: MYCOMM,MASTER 
      INTEGER*4 NPARTS,NVINV, NSRC,NSG,NREC, NORD,NLXI,NLETA, NEN, 
     ;          NNPG,NNPE,NELEM,NBLOCKS, NL1D_LT,NL1D_RT, NDOF,NZERO,
     ;          MAXIT,NNPINV,NHSIZE, IITYPE,NCON,NA35
      LOGICAL*4 LISISO, LPHESS, LGNEWT 
      INTEGER*4 MPIERR
!.... characters
      CALL MPI_BCAST(PROJNM,  80,MPI_CHARACTER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(CINVTYPE, 2,MPI_CHARACTER, MASTER,MYCOMM,MPIERR)
!.... logicals
      CALL MPI_BCAST(LISISO, 1,MPI_LOGICAL,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(LPHESS, 1,MPI_LOGICAL,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(LGNEWT, 1,MPI_LOGICAL,     MASTER,MYCOMM,MPIERR) 
!.... general sizes
      CALL MPI_BCAST(IITYPE, 1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NPARTS, 1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NVINV,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NNPINV, 1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NA35,   1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NHSIZE, 1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR) 
      CALL MPI_BCAST(MAXIT,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR) 
      CALL MPI_BCAST(NSRC ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NSG  ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NREC ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NORD ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NLXI ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NLETA,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NEN  ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NNPG ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NNPE ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NELEM,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NBLOCKS,1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NDOF ,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NZERO,  1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NL1D_LT,1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NL1D_RT,1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(NCON   ,1,MPI_INTEGER,     MASTER,MYCOMM,MPIERR)
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE BCAST_INV_MD(MYCOMM,MASTER, LMOD,LISISO, 
     ;                MDIM,MEN,MGNOD,MNPG,  
     ;                NDIM,NEN,NGNOD,NNPG, NNPE,NELEM,NDOF, NLXI,NLETA, 
     ;                NNPINV,NCON, AZMOD,XWIDTH,ZWIDTH, 
     ;                CDOMAIN,CNP, CNNPG, IGPART,MASKG,
     ;                IDOFSE,IENG,LM, 
     ;                MCONN, WMASK, XIPTS,ETAPTS, XLOCS,ZLOCS, XD,ZD, 
     ;                XLOCSE,ZLOCSE, DENS,ECOEFF)
!     Broadcasts model parameters with mesh and model 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER*4, INTENT(IN) :: MYCOMM,MASTER, 
     ;           MDIM,MEN,MGNOD,MNPG,
     ;           NDIM,NEN,NGNOD,NNPG, NNPE, NELEM,NDOF, NLXI,NLETA,
     ;           NCON,NNPINV 
      LOGICAL*4, INTENT(IN) :: LMOD, LISISO
      CHARACTER(2) CNNPG(NNPG)
      CHARACTER(1) CDOMAIN(NELEM), CNP(NDOF) 
      REAL*8 ECOEFF(MNPG,*), DENS(NNPG), XLOCS(NNPG),ZLOCS(NNPG), 
     ;       XD(NNPG),ZD(NNPG), XLOCSE(NNPE),ZLOCSE(NNPE), 
     ;       XIPTS(NLXI),ETAPTS(NLETA), 
     ;       XWIDTH,ZWIDTH,AZMOD 
      REAL*4 WMASK(NNPINV) 
      INTEGER*4 LM(MDIM,MEN,*), IENG(MGNOD,*), IDOFSE(MDIM,*), 
     ;          MCONN(MNPG,*), MASKG(NNPG), IGPART(NNPINV) 
      INTEGER*4 MPIERR,NCOEFF,IC,IA,IAE,I
!
!----------------------------------------------------------------------!
!
!.... possibly just updating model
      IF (LMOD) GOTO 100 
!
!.... characters
      CALL MPI_BCAST(CDOMAIN,NELEM,             MPI_CHARACTER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(CNP    , NDOF,             MPI_CHARACTER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(CNNPG(1:NNPG)(1:1)  , NNPG,MPI_CHARACTER, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(CNNPG(1:NNPG)(2:2)  , NNPG,MPI_CHARACTER, MASTER, 
     ;               MYCOMM,MPIERR)
!.... integers
      CALL MPI_BCAST(MASKG ,NNPG  ,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      CALL MPI_BCAST(IGPART,NNPINV,MPI_INTEGER, MASTER,MYCOMM,MPIERR)
      DO 1 IA=1,NGNOD
         CALL MPI_BCAST(IENG(IA,1:NELEM),    NELEM,MPI_INTEGER, MASTER,
     ;                  MYCOMM,MPIERR)
    1 CONTINUE
      DO 2 I=1,NDIM
         CALL MPI_BCAST(IDOFSE(I,1:NNPE),     NNPE,MPI_INTEGER, MASTER,
     ;                  MYCOMM,MPIERR) 
         DO 3 IAE=1,NEN
            CALL MPI_BCAST(LM(I,IAE,1:NELEM),NELEM,MPI_INTEGER, MASTER,
     ;                     MYCOMM,MPIERR)
    3    CONTINUE
    2 CONTINUE
      DO 4 I=1,NCON
         CALL MPI_BCAST(MCONN(1:NNPG,I),     NNPG,MPI_INTEGER, MASTER,
     ;                  MYCOMM,MPIERR)
    4 CONTINUE  
!.... reals
      CALL MPI_BCAST(WMASK,NNPINV,MPI_REAL,             MASTER,
     ;               MYCOMM,MPIERR)

!.... double precision reals
      CALL MPI_BCAST(AZMOD,     1,MPI_DOUBLE_PRECISION, MASTER, 
     ;               MYCOMM,MPIERR) 
      CALL MPI_BCAST(XWIDTH,    1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(ZWIDTH,    1,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(XIPTS , NLXI,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(ETAPTS,NLETA,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(XLOCS , NNPG,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(ZLOCS , NNPG,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(XLOCSE, NNPE,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(ZLOCSE, NNPE,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(XD    , NNPG,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      CALL MPI_BCAST(ZD    , NNPG,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
  100 CONTINUE !model parameters only; useful for updates
      CALL MPI_BCAST(DENS  , NNPG,MPI_DOUBLE_PRECISION, MASTER,
     ;               MYCOMM,MPIERR)
      IF (LISISO) THEN
         NCOEFF = 2
      ELSE
         NCOEFF = 5
      ENDIF
      DO 5 IC=1,NCOEFF
         CALL MPI_BCAST(ECOEFF(1:NNPG,IC),NNPG,MPI_DOUBLE_PRECISION,
     ;                  MASTER, MYCOMM,MPIERR)
    5 CONTINUE
      RETURN
      END

      SUBROUTINE BCAST_LOBS(MYCOMM,MASTER, MDIM,NDIM,NREC, LOBS) 
!     Broadcasts the list active observations for this frequency and 
!     source
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER*4, INTENT(IN) :: MDIM, NDIM, NREC, MYCOMM,MASTER 
      LOGICAL*4 LOBS(MDIM,*) 
      INTEGER*4 I,MPIERR 
      DO 1 I=1,NDIM
         CALL MPI_BCAST(LOBS(I,1:NREC),NREC,MPI_LOGICAL, MASTER,
     ;                  MYCOMM,MPIERR)
    1 CONTINUE
      RETURN
      END

      SUBROUTINE BCAST_OBS(MYCOMM,MASTER, MDIM,MFREQ,MREC,  
     ;                     NDIM,NFREQ,NREC,NSRC, WTS,OBS) 
!     Broadcasts the observations and data weights 
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      COMPLEX*8 OBS(MDIM,MFREQ,MREC,*)  
      REAL*4    WTS(MDIM,MFREQ,MREC,*) 
      INTEGER*4, INTENT(IN) :: MYCOMM,MASTER, MDIM,MFREQ,MREC, 
     ;                         NDIM,NFREQ,NREC,NSRC
      INTEGER*4 IRC,ISRC,I,MPIERR
      DO 1 IRC=1,NREC
         DO 2 ISRC=1,NSRC
            DO 3 I=1,NDIM
               CALL MPI_BCAST(OBS(I,1:NFREQ,IRC,ISRC),NFREQ,MPI_COMPLEX,
     ;                        MASTER,MYCOMM,MPIERR)
               CALL MPI_BCAST(WTS(I,1:NFREQ,IRC,ISRC),NFREQ,MPI_REAL   ,
     ;                        MASTER,MYCOMM,MPIERR)
    3       CONTINUE
    2    CONTINUE
    1 CONTINUE 
      RETURN
      END

      SUBROUTINE BCAST_STF(MYCOMM,MASTER, MFREQ, NFREQ,NSRC, SOURCE)
!     Broadcasts the source time function
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      COMPLEX*8 SOURCE(MFREQ,*) 
      INTEGER*4, INTENT(IN) :: MYCOMM,MASTER, MFREQ, NFREQ,NSRC  
      INTEGER*4 ISRC,MPIERR
      DO 1 ISRC=1,NSRC
         CALL MPI_BCAST(SOURCE(1:NFREQ,ISRC),NFREQ,MPI_COMPLEX, 
     ;                  MASTER,MYCOMM,MPIERR)
    1 CONTINUE
      RETURN
      END
