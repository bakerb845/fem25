      INTEGER*4 NDIM, NGNOD
      PARAMETER(NDIM = 3)   !2.5D simulation which generates 3 components of motion 
      PARAMETER(NGNOD = 4)  !at most we will use linear quadrilaterals
!----------------------------------------------------------------------------------------!
!     This structure contains all parameters for the mesh and inversion as to provide    !
!     easy modifications when backtracking and adding water, anisotropy, and viscous     !
!     damping.  To those who modify this, please be generous with comments since I can   !
!     never remember what these things do  - Ben Baker March 2013                        !
!----------------------------------------------------------------------------------------!
      TYPE MESH_INFO
         !SEQUENCE !force memory to be continugous
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Character descriptors of mesh.  CDOMAIN will need to be CHARACTER(2) when I     !
!        include water, i.e., 'AA', 'AE', 'AI', and 'EA', 'EE', and 'EI' for acoustic    !
!        and elastic distinction.  Funding pending of course                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CHARACTER(2), POINTER :: CNNPG(:)     !anchor node location 
         CHARACTER(1), POINTER :: CDOMAIN(:)   !element domain 'A' - absorbing 
                                               !               'E' - bielak 1D, 
                                               !               'I' - interior
         CHARACTER(1), POINTER :: CNP(:)       !GLL node locations 'A' - absorbing 
                                               !                   'E' - bielak 1D
                                               !                   'I' - interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Mesh parameters: These are the physical properties of the model including       !
!        the boundaries in (lat,lon) and locally.                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         REAL*8, POINTER :: ECOEFF(:,:)       !elastic coefficients (Pa), at present
                                              !lambda and mu, but later may be 
                                              !c_{11},c_{13},c_{33},c_{44},c_{66} 
         REAL*8, POINTER :: DENS(:)           !density at anchor nodes (kg/m**3)
         REAL*8, POINTER :: QP(:)             !Qp quality factor (not yet programmed)
         REAL*8, POINTER :: QS(:)             !Qs quality factor (not yet programmed)
         REAL*8, POINTER :: XD(:)             !x depths into PML at anchor nodes (m)
         REAL*8, POINTER :: ZD(:)             !z depths into PML at anchor nodes (m)
         REAL*8, POINTER :: XLOCS(:)          !x locations of anchor nodes (m)
         REAL*8, POINTER :: ZLOCS(:)          !z locations of anchor nodes (m)
         REAL*8, POINTER :: XLOCSE(:)         !x locations of GLL nodes in Bielak (m)
         REAL*8, POINTER :: ZLOCSE(:)         !z locations of GLL nodes in Bielak (m)
         REAL*8, POINTER :: XIPTS(:)          !xi interpolation points [-1,1]
         REAL*8, POINTER :: ETAPTS(:)         !eta interpolation points [-1,1] 
                                              !b/c mesh is conforming = xipts

         REAL*8 XWIDTH                        !x width of PML (m)
         REAL*8 ZWIDTH                        !z width of PLM (m)   

         REAL*8 XMOD0                         !left Bielak/Abosrbing boundary x location
         REAL*8 XMOD1                         !right Bielak/Absorbing boundary x location
         REAL*8 XBLKL                         !left Bielak/Internal interface x location
         REAL*8 XBLKR                         !right Bielak/Internal interface x location
         REAL*8 XLATMIN                       !left Bielak/Absorbing latitude
         REAL*8 XLONMIN                       !left Bielak/Absorbing longitude
         REAL*8 XLATMAX                       !right Bielak/Absorbing latitude
         REAL*8 XLONMAX                       !right Bielak/Absorbing longitude
         REAL*8 AZMOD                         !model azimuth degrees
         REAL*8 AZTOL                         !azimuth tolerance (deg) for py groups
         REAL*8 AOITOL                        !angle of inc. tol (deg) for py groups
         REAL*8 FREQ0                         !reference frequency for dispersion (Hz)
         REAL*8 FCENT_SRF                     !PML surface wave center frequency (Hz)
         REAL*8 FCENT_BDY                     !PML body wave center frequency (Hz)
         REAL*8 KAPPA_SRF                     !PML kappa number for surface waves
         REAL*8 KAPPA_BDY                     !PML kappa number for body waves
         REAL*8 RCOEFF_SRF                    !PML reflection coefficient surface waves
         REAL*8 RCOEFF_BDY                    !PML reflection coefficient body waves
         REAL*8 DMAX_SRF                      !PML dmax override
         REAL*8 DMAX_BDY                      !PML dmax override
         REAL*8 ZBASE_INT                     !z base of interior domain
         REAL*8 ZBASE_BLK                     !z base of bielak domain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Assembly pointers: These are the pointers critical for assembly.  The assembly  ! 
!        is distributed, i.e., each process holds a subset of the global matrix.  The    !
!        FORTRAN friendly storage format I use is compressed sparse row and it's         !
!        distributed compressed sparse row counterpart                                   !
!        http://netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000 !
!        Moreover, the mesh partition results in distributed matrices that are disjoint  !
!        so the global matrix is not shared along an edgecut thereby reducing memory     !
!        consumption                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         INTEGER*4, POINTER :: LM(:,:,:)           !location matrix [ndim,nen,nelem]
                                                   !element node -> global DOF # 
         INTEGER*4, POINTER :: IENG(:,:)           !anchor node IEN matrix [ngnod,nelem]
                                                   !element anchor node -> anchor node #
         INTEGER*4, POINTER :: IDOFSE(:,:)         !bielak DOF numbers [ndim,nnpe]
         INTEGER*4, POINTER :: IDOF_FS(:,:)        !DOF numbers of nodes at free surface
                                                   ![ndim,nfs] 
         INTEGER*4, POINTER :: IRPTR(:)            !global row CRS pointer 
                                                   !freed after distributed graph made
                                                   ![ndof+1]
         INTEGER*4, POINTER :: JCPTR(:)            !global column CRS pointer
                                                   !freed after distributed graph made
                                                   ![nzero]
         INTEGER*4, POINTER :: MYELEM(:)           !elements local to process for assembly
         INTEGER*4, POINTER :: MYDOFS(:)           !this is a processes global DOF #s 
                                                   ![ndofl]
         INTEGER*4, POINTER :: IRPTR_LOC(:)        !local CRS row pointer, each index 
                                                   !corresponds to global DOF in MYDOFS
                                                   ![ndofl+1]
         INTEGER*4, POINTER :: JCPTR_LOC(:)        !local CRS column pointer, each entry
                                                   !corresponds to global DOF # [nzloc]
         INTEGER*4, POINTER :: PART(:)             !determines which DOF belongs to which
                                                   !process.  for example, 3 processes 
                                                   !means PART only has values (1,2,3) 
                                                   ![ndof] 
         !LOGICAL*4, POINTER :: LQUAD(:)            !in the future i'd lke triangles in 
                                                    !mesh as well
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Integer sizes                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         INTEGER*4 NELEM           !number of elements in mesh
         INTEGER*4 NELEML          !number of elements involved in processes' asssembly
         INTEGER*4 NDOF            !number of rows of global impedance matrix 
         INTEGER*4 NDOFL           !number of rows to local impedance matrix
                                   !sum ndofl_{i} = ndof
         INTEGER*4 NZERO           !number of non-zeros in global impedance matrix 
         INTEGER*4 NZLOC           !number of non-zeros in local impedance matrix 
                                   !sum nzloc_{i} = nzero
         INTEGER*4 NNPG            !number of anchor nodes in mesh
         INTEGER*4 NNPE            !number of bielak nodes in Bielak domain (GLL grid)
         INTEGER*4 NORD            !polynomial order for interpolation
         INTEGER*4 NLXI            !number of lagrange interpolant points in xi 
         INTEGER*4 NLETA           !lagrange interp points in eta, = nlxi b/c mesh is 
                                   !conforming
         INTEGER*4 NEN             !number of element nodes = nlxi*nleta 
         INTEGER*4 IITYPE          !interpolation type 
                                   ! (1) GLL -> Gauss-Legendre-Lobatto default
                                   ! (2) GLC -> Gauss-Legendre-Chebyshev 
                                   ! (3) Equal spacing (strong discouraged)
         INTEGER*4 NCOEFF          !number of coefficients in ECOEFF
         INTEGER*4 NFS             !number of fine nodes at free surface

         LOGICAL*4 LDISP           !True -> medium is dispersive (Qp/Qs,freq0 active)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        A lazy descriptor on anisotropy level.  What I should do is classify the        !
!        elastic coefficients with character descriptors then assemble based on the      !
!        number of aniostropic coefficients.  As for now, all we have is isotropy but    !
!        someone can rewrite the IO routines, update CJAC, piggy back off the heuristics !
!        in the element matrix integration to finish the job                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         LOGICAL*4 LISISO          !True -> simulation is isotropic 
      END TYPE MESH_INFO
!----------------------------------------------------------------------------------------!
!     This concludes the mesh and model parameters.  From this we can now derived the    !
!     1D models used for generating background solutions in the Bielak layers.           !
!----------------------------------------------------------------------------------------!
      TYPE MOD1D_INFO  
         !SEQUENCE
         REAL*8, POINTER :: VP1D_LT(:)    !Vp velocity of 1D model on left side
         REAL*8, POINTER :: VS1D_LT(:)    !Vs velocity of 1D model on left side
         REAL*8, POINTER :: RH1D_LT(:)    !density of 1D model on left side
         REAL*8, POINTER :: QP1D_LT(:)    !P quality factor left side
         REAL*8, POINTER :: QS1D_LT(:)    !S quality factor left side
         REAL*8, POINTER ::  Z1D_LT(:)    !z depths of 1D model on left side
         REAL*8, POINTER :: VP1D_RT(:)    !Vp velocity of 1D model on right side
         REAL*8, POINTER :: VS1D_RT(:)    !Vs velocity of 1D model on right side
         REAL*8, POINTER :: RH1D_RT(:)    !density of 1D model on right side
         REAL*8, POINTER :: QP1D_RT(:)    !P quality factor right side
         REAL*8, POINTER :: QS1D_RT(:)    !S quality factor right side
         REAL*8, POINTER ::  Z1D_RT(:)    !z depths of 1D model on right side
         REAL*8, POINTER :: VPD_RLLT(:)   !Vp vel of 1D unflattened rayleigh model, left
         REAL*8, POINTER :: VSD_RLLT(:)   !Vs vel of 1D unflattened rayleigh model, left
         REAL*8, POINTER :: ROD_RLLT(:)   !Density of 1D unflattened rayleigh model, left
         REAL*8, POINTER :: HDD_RLLT(:)   !Depths of 1D unflattened rayleigh model, left
         REAL*8, POINTER :: VPD_RLRT(:)   !Vp vel of 1D unflattened rayleigh model, right
         REAL*8, POINTER :: VSD_RLRT(:)   !Vs vel of 1D unflattened rayleigh model, right
         REAL*8, POINTER :: ROD_RLRT(:)   !Density of 1D unflattened rayleigh model, right
         REAL*8, POINTER :: HDD_RLRT(:)   !Depths of 1D unflattened rayleigh model, right
         REAL*8, POINTER :: VPD_LVLT(:)   !Vp velocity of 1D unflattened love model, left
         REAL*8, POINTER :: VSD_LVLT(:)   !Vs velocity of 1D unflattened love model, left
         REAL*8, POINTER :: ROD_LVLT(:)   !Density of 1D unflattened love model, left
         REAL*8, POINTER :: HDD_LVLT(:)   !Depths of 1D unflattened love model, left
         REAL*8, POINTER :: VPD_LVRT(:)   !Vp velocity of 1D unflattened love model, right
         REAL*8, POINTER :: VSD_LVRT(:)   !Vs velocity of 1D unflattened love model, right
         REAL*8, POINTER :: ROD_LVRT(:)   !Density of 1D unflattened love model, right
         REAL*8, POINTER :: HDD_LVRT(:)   !Depths of 1D unflattened love model, right

         REAL*8 VFAST                     !a max velocity for calculating py tolerances
         REAL*8 VFAST_BDY                 !max body wave velocity
         REAL*8 VFAST_SRF                 !max surface wave velocity

         INTEGER*4 NL1D_LT                !number of 1D points on left model
         INTEGER*4 NL1D_RT                !number of 1D points in right model
      END TYPE MOD1D_INFO
!----------------------------------------------------------------------------------------!
!     This concludes the 1D models.  Now we define the receiver information.             !
!----------------------------------------------------------------------------------------!
      TYPE RECV_INFO
         !SEQUENCE
         COMPLEX*8, POINTER :: RECV(:,:,:)    !receiver response functions
                                              !if it cant be read from file it is set 
                                              !to one [ndim,nfreq,nrec]
         REAL*8, POINTER :: YREC(:)           !y receiver locations (m) [nrec]
         INTEGER*4, POINTER :: MRDOF(:,:)     !maps receiver/component to global DOF #
                                              ![ndim,nrec]
         INTEGER*4 NREC                       !number of receivers
      END TYPE RECV_INFO
!----------------------------------------------------------------------------------------!
!     This concludes the receiver information.  Now we store the source information.     !
!----------------------------------------------------------------------------------------!
      TYPE SRC_INFO 
         !SEQUENCE
         CHARACTER(2), POINTER :: SRCTYP(:)   !source type body/surface  
         CHARACTER(1), POINTER :: CSIDE(:)    !side of model sources approaches from 
                                              !'L' or 'R' for left or right 
         COMPLEX*8, POINTER :: SOURCE(:,:)    !source time function [nfreq,nsrc]
         REAL*8, POINTER :: PYTAB(:,:)        !apparent slowness in y table [nfreq,nsrc] 
         REAL*8, POINTER :: AOI(:)            !angle of incididene (degrees)
         REAL*8, POINTER :: BAZN(:)           !corrected source to model azimuth (degrees)
                                              !i apologize for the name.  this is NOT
                                              !back-azimuth or new back-azimuth [nsrc]
         REAL*8, POINTER :: SLAT(:)           !source latitude (degrees)
         REAL*8, POINTER :: SLON(:)           !source longitude (degrees)
         REAL*8, POINTER :: SDEP(:)           !source depth (km)
         REAL*8, POINTER :: STRIKE(:)         !source strike angle (degrees) 
         REAL*8, POINTER :: DIP(:)            !source dip angle (degrees)
         REAL*8, POINTER :: RAKE(:)           !source rake angle (degrees)
         REAL*8, POINTER :: SMAG(:)           !source magnitude (dyne-cm)

         REAL*8, POINTER :: PYAVG(:)          !average py value in group 
         INTEGER*4, POINTER :: MODE(:)        !desired mode 
         INTEGER*4, POINTER :: ISRCPRM(:)     !source permutation pointer [nsrc]
         INTEGER*4, POINTER :: ISGPTR(:)      !source group pointer [nsg]

         INTEGER*4 NSG                        !number of source groups
         INTEGER*4 NSRC                       !Number of sources = nsrc_srf + nsrc_bdy
         INTEGER*4 NSRC_SRF                   !Number of surface wave sources
         INTEGER*4 NSRC_BDY                   !Number of body wave sources
      END TYPE SRC_INFO 
!----------------------------------------------------------------------------------------!
!     This concludes the source information.  Now the frequency information              !
!----------------------------------------------------------------------------------------!
      TYPE FRQ_INFO
         !SEQUENCE
         CHARACTER(1), POINTER :: CFTYPE(:) !frequency type S - surface wave
                                            !               B - body wave
         REAL*8, POINTER :: FREQ(:)         !frequency list (Hz)
         LOGICAL*4, POINTER :: LINVF(:)     !True -> frequency is an inversion frequency
         INTEGER*4 NFREQ                    !number of frequencies in list
         INTEGER*4 NFREQ_SRF                !number of surface wave modeling frequencies
         INTEGER*4 NFREQ_BDY                !number of body wave modeling frequencies
         INTEGER*4 NFREQ_INV                !number of inversion frequencies
         INTEGER*4 NFREQ_SRF_INV            !Number of surface wave inversion frequencies
         INTEGER*4 NFREQ_BDY_INV            !Number of body wave inversion frequencies
      END TYPE FRQ_INFO
!----------------------------------------------------------------------------------------! 
!                      THIS IS THE END OF THE FORWARD MODELING                           !
!                        THIS BEGINS THE INVERSION VARIABLES                             !
!----------------------------------------------------------------------------------------!
      TYPE INV_INFO
         !SEQUENCE

         COMPLEX*8, POINTER :: EST(:,:,:,:)  !estimates    [ndim,nfreq,nrec,nsrc]
         COMPLEX*8, POINTER :: OBS(:,:,:,:)  !observations [ndim,nfreq,nrec,nsrc]
                                             !ndim format -> [N,E,Z] coordinates
         REAL*8, POINTER :: ELEM_WTS(:)      !normalization for each element in [dS/dm u]
                                             !virtual force calculation [nelem]
         REAL*8, POINTER :: DX(:)            !receiver spacing for data covar matrix
         REAL*4, POINTER :: WGHTS(:,:,:,:)   !data weights [ndim,nfreq,nrec,nsrc]
         REAL*4, POINTER :: GRADPC(:)        !gradient preconditioner [nhsize]
         REAL*4, POINTER :: GRAD(:)          !gradient [na35]
         REAL*4, POINTER :: GRAD_SRF(:)      !gradient [na35] for surface waves 
         REAL*4, POINTER :: GRAD_BDY(:)      !gradient [na35] for body waves
         REAL*4, POINTER :: WMASK(:)         !gradient mask [na35]
         REAL*4, POINTER :: RHS(:)           !RHS for Gauss-Newton (residuals) or 
                                             !migration (data) [nobs] 

         INTEGER*4, POINTER :: MCONN(:,:)    !anchor nodes connection [nnpg,ncon]
         INTEGER*4, POINTER :: MASKG(:)      !anchor node -> gradient nodal point [nnpg]
         INTEGER*4, POINTER :: IGPART(:)     !each inpinv partition number [nnpinv]
         INTEGER*4, POINTER :: IPIVH(:)      !pivots for PCHESS [na35]
         INTEGER*4, POINTER :: MYGRAD(:)     !points to global anchor nodes in gradient
                                             !that belong to process [nnpgl]
         INTEGER*4, POINTER :: ICSC_FDIST(:) !local CSC row pointer for gradient 
                                             ![nz_fdist] 
         INTEGER*4, POINTER :: JCSC_FDIST(:) !Local CSC column pointer for gradient  
                                             ![na35 + 1]  

         REAL*8 VPVS                         !Vp/Vs ratio 
         REAL*8 SCLBDY                       !grad = grad_srf + sclbdy*grad_bdy

         REAL*4 FOBJ                         !objective function
         REAL*4 FOBJ_SRF                     !objective function from surface wave misfit 
         REAL*4 FOBJ_BDY                     !objective function from body wave misfit
         REAL*4 PTHRESH                      !singular value threshold percentage [0,100] 
                                             !if pthresh = 2., we would extract the 
                                             !corresponding singular value at the pthresh
                                             !,s2, so the block diagonal pre-conditioner
                                             !is inflated D = U S V^H = U (S + s2 I) V^H 
         REAL*4 DTMAX_SRF                    !max residual (s) for surface waves
         REAL*4 DTMAX_BDY                    !max residual (s) for body waves
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     Integer Sizes                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         INTEGER*4 MCJLOC                    !leading dimension for Jacobian [ndim*nrec]
         INTEGER*4 NHSIZE                    !hessian preconditioner size 
                                             !if cinvtype == 'PS' [nnpinv*nvinv**2]
                                             !which is the nvinv x nvinv block diagonal  
                                             !for nnpinv nodes structure 
                                             !if cinvtype == 'PP' or 'SS' [na35]
         INTEGER*4 NVINV                     !number of variables to invert for 
                                             !e.x.: 1 -> P or S, 2 -> P and S
         INTEGER*4 NNPINV                    !number of nodal points to invert at  
         INTEGER*4 NA35                      !number of variables in inversion; 
         INTEGER*4 NCON                      !max number of connections for an node

         INTEGER*4 MBUFRHS                   !max buffer communication size in Jacobian 
                                             !calculation
         INTEGER*4 NNPGL                     !number of nodal points in local gradient 
                                             ! sum nnpgl_i = na35 
         INTEGER*4 NZ_FDIST                  !number of non-zeros in distributed 
                                             !virtual force matrix for process 
         INTEGER*4 NOBS                      !number of observations (residuals) 
         INTEGER*4 NREC_MIN                  !min number of recs when making data covar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   Controls on inversion                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         INTEGER*4 MAXIT                     !max number of iterations in inversion
         INTEGER*4 NORM                      !norm in inversion 
                                             !   1 -> L1 
                                             !   2 -> L2 (default) 
         INTEGER*4 IRESTP                    !Residual type for inversion
                                             !   1 -> Phase only
                                             !   2 -> Amplitude only
                                             !   3 -> Phase and amplitude (default)
         INTEGER*4 IMODSRC                   !Residual type for source correction  
                                             !Same conventions as IRESTP except 0 
                                             !toggles STF update off 
         INTEGER*4 IMODREC                   !Residual type for receiver correction
                                             !Same convetnions as IRESTP except 0 
                                             !toggles receiver update off
         INTEGER*4 IBPHASE                   !Geometric spreading correction in backprop
                                             !   0 -> 1 (default)
                                             !   1 -> sqrt(-i omega)
                                             !   2 ->-sqrt(-i omega)
                                             !   3 -> sqrt( i omega)
                                             !   4 ->-sqrt( i omega)
                                             !   5 -> i (90 degrees) 
         INTEGER*4 NCASC                     !Number of cascades to run a nearest 
                                             !neighbor averaging filter on gradient
         LOGICAL*4 LPGRAD                    !calculate gradient pre-conditioner
                                             ! defined as diag[adj(J) J]
                                             !where the gradient is given by 
                                             ! g <- inv[diag(adj(J) J) + Pad] g
                                             !For padding, see pthresh
         LOGICAL*4 LMIGR                     !True -> only perform migration
         LOGICAL*4 LGNEWT                    !True -> performing Gauss-Newton step
                                             !     -> This will save the Jacobians
         LOGICAL*4 LFUNC                     !calculate an objective function
         LOGICAL*4 LGRAD                     !calculate a gradient
         LOGICAL*4 LUNWRAP                   !True OBS is held as (mag,phase) and we 
                                             !calculate only phase residuals 
                                             !False OBS is held as a complex number and 
                                             !we can invert, phase, amplitude, or both 
         LOGICAL*4 LSURF                     !Inverting surface waves 
         LOGICAL*4 LBODY                     !Inverting body waves
         LOGICAL*4 LCOVD_SRF                 !True -> surface FD data covariance matrix
         LOGICAL*4 LCOVD_BDY                 !True -> body wave FD data covariance matrix 
         LOGICAL*4 LAHESS_HOT                !True -> approximate higher order terms in 
                                             !Hessian (default false)
         LOGICAL*4 LDWGHT                    !True -> data covariance distance weighting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   Inversion Descriptors                                ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CHARACTER(2) CINVTYPE               !inversion type: 
                                             !'PP' -> P velocity inversion
                                             !'SS' -> S velocity inversion
                                             !'PS' -> P and S velocity inversion
      END TYPE INV_INFO
!----------------------------------------------------------------------------------------!
!                        THIS IS THE END OF THE INVERSION INFO                           !
!                            BEGIN THE WINDOW INFORMATION                                !
!----------------------------------------------------------------------------------------!
      TYPE WIN_INFO 
         REAL*4, POINTER :: TWIN(:,:,:)   !Window start and stop times [nrec,nsrc,2] 

         REAL*4 DT_SRF      !sampling period (s) for surface wave seismograms 
         REAL*4 DT_BDY      !sampling period (s) for body wave seismograms  
         REAL*4 START_SRF   !start time (s) for surface wave seismograms 
         REAL*4 START_BDY   !start time (s) for body wave seismograms

         INTEGER*4 NSAMP_SRF !number of samples in surface wave seismograms
         INTEGER*4 NSAMP_BDY !number of samples in body wave seismograms
         INTEGER*4 IWNDO     !(1) Rectangle window
                             !(2) Tapered rectangle 10 pct length rolloff (default)
                             !(3) Triangle window 
                             !(4) Hanning window
                             !(5) Hamming window
                             !(6) Blackman (optimized coeffs ala Wikipedia) 
         INTEGER*4 NPCT_SRF  !percent of surface wave synthetic to apply window to
         INTEGER*4 NPCT_BDY  !percent of body wave synthetic to apply window to
 
         LOGICAL*4 LWNDO_SRF !True windowing surface wave synthetics
         LOGICAL*4 LWNDO_BDY !True windowing body wave syntheitcs

      END TYPE WIN_INFO 
