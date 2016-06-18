!
!     Mesh parameters
      CHARACTER*2, POINTER :: CNNPG(:)   !anchor node location
      CHARACTER*1, POINTER :: CDOMAIN(:) !elements domain `A', `E', or `I' 
      CHARACTER*1, POINTER :: CNP(:)     !fine node location 

      COMPLEX*16, POINTER :: U0(:)       !1D wavefield at all DOFs

      REAL*8, POINTER :: ECOEFF(:,:)     !elastic coefficients at anchor nodes
      REAL*8, POINTER :: DENS(:)         !density at anchor nodes 
      REAL*8, POINTER :: XD(:)           !x depths into PML at anchor nodes
      REAL*8, POINTER :: ZD(:)           !z depths into PML at anchor nodes
      REAL*8, POINTER :: XLOCS(:)        !x locations of anchor nodes
      REAL*8, POINTER :: ZLOCS(:)        !z locations of anchor nodes
      REAL*8, POINTER :: XLOCSE(:)       !x locations of bielak nodes
      REAL*8, POINTER :: ZLOCSE(:)       !z locations of bielak nodes
      REAL*8, POINTER :: XIPTS(:)        !xi interpolation points
      REAL*8, POINTER :: ETAPTS(:)       !eta interpolation points
      REAL*8 XWIDTH                      !x width of PML
      REAL*8 ZWIDTH                      !z width of PML

      INTEGER*4, POINTER :: LM(:,:,:)    !location matrix
      INTEGER*4, POINTER :: IENG(:,:)    !global anchor node pointer
      INTEGER*4, POINTER :: IDOFSE(:,:)  !holds DOFs in Bielak domain
      INTEGER*4, POINTER :: MYDOFS(:)    !DOFs that belong to a process
      INTEGER*4, POINTER :: MYELEM(:)    !Holds my elements
      INTEGER*4, POINTER :: IRPTR_LOC(:) !Local CRS row pointer
      INTEGER*4, POINTER :: JCPTR_LOC(:) !Local CRS column pointer
      INTEGER*4, POINTER :: IRPTR(:)     !CRS row pointer
      INTEGER*4, POINTER :: JCPTR(:)     !CRS column pointer
      INTEGER*4, POINTER :: PART(:)      !Partition number of each element

      INTEGER*4 NELEM                    !Number of elements in mesh
      INTEGER*4 NELEML                   !Number of local elements
      INTEGER*4 NDOF                     !Number of degrees of freedom
      INTEGER*4 NDOFL                    !Number of local degrees of freedom
      INTEGER*4 NZERO                    !Number of non-zeros in global matrix
      INTEGER*4 NZLOC                    !Number of local non-zeros
      INTEGER*4 NNPG                     !Number of anchor nodes
      INTEGER*4 NNPE                     !Number of nodal points in Bielak domain
      INTEGER*4 NLXI                     !Number of lagrange interpolation points in xi
      INTEGER*4 NLETA                    !Number of lagrange interpolation points in eta 
      INTEGER*4 NEN                      !Number of element nodes
      INTEGER*4 IITYPE                   !Interpolation type (1) GLL, default (2) GLC (3) equal

      LOGICAL*4 LISISO                   !True then simulaiton is isotropic 

      REAL*8, POINTER :: VP1D_LT(:)      !Vp velocity of 1D model on left side
      REAL*8, POINTER :: VS1D_LT(:)      !Vs velocity of 1D model on left side
      REAL*8, POINTER :: RH1D_LT(:)      !density of 1D model on left side
      REAL*8, POINTER ::  Z1D_LT(:)      !z depths of 1D model on left side
      REAL*8, POINTER :: VP1D_RT(:)      !Vp velocity of 1D model on right side
      REAL*8, POINTER :: VS1D_RT(:)      !Vs velocity of 1D model on right side
      REAL*8, POINTER :: RH1D_RT(:)      !density of 1D model on right side
      REAL*8, POINTER ::  Z1D_RT(:)      !z depths of 1D model on right side
      REAL*8 XMOD0                       !left Bielak/Abosrbing boundary x location
      REAL*8 XMOD1                       !right Bielak/Absorbing boundary x location
      REAL*8 XBLKL                       !left Bielak/Internal interface x location
      REAL*8 XBLKR                       !right Bielak/Internal interface x location
      INTEGER*4 NL1D_LT                  !number of 1D points on left model
      INTEGER*4 NL1D_RT                  !number of 1D points in right model

      INTEGER*4 NDIM,NGNOD
      PARAMETER(NDIM = 3)                !Number of dimensions in solution
      PARAMETER(NGNOD = 4)               !Number of anchor nodes on element

      COMMON /MATERIAL/ ECOEFF,DENS,XD,ZD, XLOCS,ZLOCS,XLOCSE,ZLOCSE, XIPTS,ETAPTS,      &
                        XWIDTH,ZWIDTH, CNNPG,CDOMAIN, NNPG,NNPE, NLXI,NLETA,NEN,         &
                        NELEM,NELEML, IITYPE, LISISO
      COMMON /GRAPH/ LM,IENG,IDOFSE,MYDOFS,MYELEM, IRPTR_LOC,JCPTR_LOC,PART,IRPTR,JCPTR, &
                     CNP, NDOF,NDOFL, NZERO,NZLOC 

      COMMON /MOD1D/ VP1D_LT,VS1D_LT,RH1D_LT,Z1D_LT, VP1D_RT,VS1D_RT,RH1D_RT,Z1D_RT,     &
                     XBLKL,XBLKR, XMOD0,XMOD1, NL1D_LT,NL1D_RT

