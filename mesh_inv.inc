!
!.... observations and corresponding weights (ndim,nfreq,nrec,nsrc)
      COMPLEX*8, POINTER :: OBS(:,:,:,:)  !Observations at receivers
      COMPLEX*8, POINTER :: EST(:,:,:,:)  !Estimate responses at recievers
      REAL*4, POINTER ::  WGHTS(:,:,:,:)  !Data weights
!
!.... gradient, hessian, objective function, and mask 
      REAL*4, POINTER :: GRAD(:)       !Gradient length nnpinv*nvinv 
      REAL*4, POINTER :: PCHESS(:)     !Gradient preconditioner length nnpinv*nvinv**2 
      REAL*4, POINTER :: WMASK(:) 
      INTEGER*4, POINTER :: MCONN(:,:) !holds anchor node connections 
      INTEGER*4, POINTER :: MASKG(:)   !Length nnpg points to inpinv in gradient
      INTEGER*4, POINTER :: IGPART(:)  !Holds gradient partition IDs, length nnpinv
      INTEGER*4, POINTER :: IPIVH(:)   !Pivots for PCHESS factorization, length nnpinv*nvinv
      REAL*4 FOBJ                      !Objective function
!
!.... gradient calculation should save a graph structure
      INTEGER*4, POINTER :: MYGRAD(:)     !points to global gradient nodes, like MYDOFS
      INTEGER*4, POINTER :: ICRS_GLOC(:)  !local CRS gradient row pointer
      INTEGER*4, POINTER :: JCRS_GLOC(:)  !local CRS gradient column pointer
      INTEGER*4 NNPGL   !number of rows in distributed gradient CRS pointer
      INTEGER*4 NNZGL   !number of non-zeros in distributed gradient CRS pointer
!.... local on slaves
      INTEGER*4, POINTER :: ICSC_FDIST(:)  !local virtual force matrix CSC row pointer
      INTEGER*4, POINTER :: JCSC_FDIST(:)  !local virtual force matrix CSC column pointer 
      INTEGER*4 NZ_FDIST          !Number of non-zeros in local virtual force F matrix 
      INTEGER*4 NHSIZE            !size of PCHESS
      INTEGER*4 MBUFRHS           !max buffer size for RHS communication
!
!.... inversion type, number of variables in inversion, number of nodal pts
      CHARACTER*2 CINVTYPE         !inversion type
      INTEGER*4 NVINV              !number of variables to invert for: P, S...
      INTEGER*4 NNPINV             !number of nodal points to invert at <= nnpg
      INTEGER*4 NOMINV             !number of frequencies to invert at
      INTEGER*4 NA35               !number of variables in grad = nnpinv*nvinv
      INTEGER*4 NCON               !max number of connected elements
!
!.... standard gardient inversion info
      COMMON /INVBLK/ OBS,EST, WGHTS, WMASK, MCONN,  &
                      FOBJ, NA35,NVINV,NOMINV,NNPINV,NCON, CINVTYPE 
!
!.... gradient bock
      COMMON /GRDBLK/ GRAD, MASKG,IGPART, ICRS_GLOC,JCRS_GLOC,MYGRAD, NNPGL,NNZGL
!
!.... hessian block, needs F = dS/dm u which is distributed
      COMMON /FMATBLK/ PCHESS, ICSC_FDIST,JCSC_FDIST, NHSIZE, NZ_FDIST, MBUFRHS  
!
!.... source information
      CHARACTER*2, POINTER :: SRCTYP(:)
      CHARACTER*1, POINTER :: CSIDE(:) 
      COMPLEX*16, POINTER :: SOURCE(:,:)
      REAL*8, POINTER :: AOI(:), BAZN(:), PYAVG(:) 
      INTEGER*4, POINTER :: ISRCPRM(:), ISGPTR(:)
      INTEGER*4 NSG  !number of source groups
      INTEGER*4 NSRC !Number of sources
      COMMON /SRCINFO/ SOURCE, AOI,BAZN,PYAVG, SRCTYP,CSIDE,  &
                       ISRCPRM,ISGPTR, NSRC,NSG
!
!.... receiver information
      REAL*8, POINTER :: YREC(:)         !Y receiver locations
      INTEGER*4, POINTER :: MRDOF(:,:)   !receiver DOF numbers
      INTEGER*4 NREC                     !Number of receivers
      COMMON /RECV/ MRDOF,YREC,NREC
