      SUBROUTINE ORDKRG(MPTS,NPTS, NVAR, DECSCL,DECFAC,X0,Y0,
     ;                  XLOCS,YLOCS, OBS, RKRG,IERR) 
!
!     Performs an ordinary kriging on a data set.  Here we use a 
!     scaling factor based on distance for calculating the covariance
!     DECFAC*exp(-h/DECSCL).  Algorithm 4.2 of Geostatistics for 
!     Engineers and Earth Scientsits - Olea.   
!     -  Ben Baker 
!
!     INPUT      MEANING
!     -----      ------- 
!     DECFAC     distance scalaing factor in exponential > 0
!     DECSCL     distance scaling factor on exponential > 0
!     MPTS       leading dimension for observation list
!     NPTS       number of observation points
!     NVAR       number of variables to interpolate at each point 
!     OBS        observations at xlocs, ylocs
!     X0         x point to estimate
!     XLOCS      x locations of observations
!     Y0         y point to estimate
!     YLOCS      y locations of observations
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error in factorization
!     RKRG       interpolated observations 
!
!.... variable declarations
      REAL*8, INTENT(IN) :: XLOCS(NPTS), YLOCS(NPTS), X0,Y0, 
     ;                      DECFAC,DECSCL 
      REAL*4, INTENT(IN) :: OBS(MPTS,*) 
      INTEGER*4 MPTS,NPTS,NVAR 
      REAL*4, INTENT(OUT) :: RKRG(NVAR)
!.... local variables
      REAL*8, ALLOCATABLE :: VMAT(:,:), V(:), OBS8(:) 
      INTEGER*4, ALLOCATABLE :: IPIV(:) 
!
!----------------------------------------------------------------------!
!
!.... calculate the SPD covariance matrix and RHS 
      IERR = 0
      LDV = NPTS + 1
      ALLOCATE(VMAT(LDV,LDV))
      N = NPTS + 1
      ALLOCATE(V(N))
      CALL GVMAT_KRIGE(LDV, NPTS,1,DECSCL,DECFAC,X0,Y0,
     ;                 XLOCS,YLOCS, V,VMAT) 
!
!.... LU factor matrix 
      ALLOCATE(IPIV(N))
      CALL DGETRF(N,N,VMAT,LDV,IPIV,INFO)
      IF (INFO.LT.0) THEN
         WRITE(*,*) 'ordkrg: DGETRF illegal variable:',INFO
         IERR = 1
         GOTO 730
      ENDIF
      IF (INFO.GT.0) THEN
         WRITE(*,*) 'ordkrg: DGETRF Singular matrix:',INFO
         IERR = 1
         GOTO 730
      ENDIF  
      CALL DGETRS('N',N,1,VMAT,LDV,IPIV,V,N,INFO)
      IF (INFO.LT.0) THEN
         WRITE(*,*) 'ordkrg: DGETRS illegal variable:',INFO
         IERR = 1
         GOTO 730
      ENDIF
!
!.... dot the weights with the observations to estimate value
      ALLOCATE(OBS8(N))
      OBS8(N) = 0.D0
      DO 1 IVAR=1,NVAR
         DO 2 IPTS=1,NPTS
            OBS8(IPTS) = DBLE(OBS(IPTS,IVAR))
    2    CONTINUE 
         RKRG(IVAR) = SNGL(DOT_PRODUCT(OBS8,V))
    1 CONTINUE  
!
!.... clear memory
  730 CONTINUE !break ahead for errors
      DEALLOCATE(VMAT)
      DEALLOCATE(V) 
      DEALLOCATE(IPIV) 
      IF (ALLOCATED(OBS8)) DEALLOCATE(OBS8) 
      RETURN
      END
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE SIMKRG(MPTS,NPTS, NVAR, DECSCL,DECFAC,X0,Y0,
     ;                  XLOCS,YLOCS, OBS, RKRG,IERR) 
!
!     Performs a normal kriging on a data set.  Here we use a 
!     scaling factor based on distance for calculating the covariance
!     DECFAC*exp(-h/DECSCL).  Algorithm 2.1 of Geostatistics for 
!     Engineers and Earth Scientsits - Olea.   
!     -  Ben Baker 
!
!     INPUT      MEANING
!     -----      ------- 
!     DECFAC     distance scalaing factor in exponential > 0
!     DECSCL     distance scaling factor on exponential > 0
!     MPTS       leading dimension for observation list
!     NPTS       number of observation points
!     NVAR       number of variables to interpolate at each point 
!     OBS        observations at xlocs, ylocs
!     X0         x point to estimate
!     XLOCS      x locations of observations
!     Y0         y point to estimate
!     YLOCS      y locations of observations
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       error in factorization
!     RKRG       interpolated observations 
!
!.... variable declarations
      REAL*8, INTENT(IN) :: XLOCS(NPTS), YLOCS(NPTS), X0,Y0, 
     ;                      DECFAC,DECSCL 
      REAL*4, INTENT(IN) :: OBS(MPTS,*) 
      INTEGER*4 MPTS,NPTS,NVAR 
      REAL*4, INTENT(OUT) :: RKRG(NVAR)
!.... local variables
      REAL*8, ALLOCATABLE :: VMAT(:,:), V(:), OBS8(:) 
      REAL*8 DMEAN
!
!----------------------------------------------------------------------!
!
!.... calculate the SPD covariance matrix and RHS 
      IERR = 0
      LDV = NPTS 
      ALLOCATE(VMAT(LDV,LDV))
      N = NPTS 
      ALLOCATE(V(N))
      CALL GVMAT_KRIGE(LDV, NPTS,0,DECSCL,DECFAC,X0,Y0,
     ;                 XLOCS,YLOCS, V,VMAT)
!
!.... Cholesky factor matrix 
      CALL DPOTRF('U',N,VMAT,LDV,INFO)
      IF (INFO.LT.0) THEN
         WRITE(*,*) 'ordkrg: DPOTRF illegal variable:',INFO
         IERR = 1
         GOTO 730
      ENDIF
      IF (INFO.GT.0) THEN
         WRITE(*,*) 'ordkrg: DPOTRF Singular matrix:',INFO
         IERR = 1
         GOTO 730
      ENDIF  
!
!.... V now becomes the krige weights
      CALL DPOTRS('U',N,1,VMAT,LDV,V,N,INFO)
      IF (INFO.LT.0) THEN
         WRITE(*,*) 'ordkrg: DPOTRS illegal variable:',INFO
         IERR = 1
         GOTO 730
      ENDIF
!
!.... dot the weights with the observations to estimate value
      ALLOCATE(OBS8(N))
      DO 1 IVAR=1,NVAR
!
!....... copy observations and demean
         DMEAN = 0.D0
         DO 2 IPTS=1,NPTS
            OBS8(IPTS) = DBLE(OBS(IPTS,IVAR))
            DMEAN = DMEAN + OBS8(IPTS)
    2    CONTINUE
         DMEAN = DMEAN/DFLOAT(NPTS)
         DO 3 IPTS=1,NPTS
            OBS8(IPTS) = OBS8(IPTS) - DMEAN
    3    CONTINUE 
         RKRG(IVAR) = SNGL(DMEAN + DOT_PRODUCT(OBS8,V))
    1 CONTINUE
!
!.... clear memory
  730 CONTINUE !break ahead for errors
      DEALLOCATE(VMAT)
      DEALLOCATE(V)
      IF (ALLOCATED(OBS8)) DEALLOCATE(OBS8)
      RETURN
      END

!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GDMAT_KRIGE(MPTS, NPTS,X0,Y0,XLOCS,YLOCS, DIST,DMAT)
!
!     Generates the distance matrix and distance from known points 
!     to krige point 
!
!     INPUT      MEANING
!     -----      ------- 
!     MPTS       leading dimension for DMAT
!     NPTS       number of points in x and z locations
!     X0         x location of point to krige 
!     XLOCS      x locations of observations
!     Y0         y location of point to krige
!     YLOCS      y locations of observations
!
!     OUTPUT     MEANING
!     ------     ------- 
!     DIST       distance of krige point to known points 
!     DMAT       distance matrix between known points
!
!.... variable declarations
      REAL*8, INTENT(IN) :: XLOCS(NPTS), YLOCS(NPTS), X0,Y0
      INTEGER*4, INTENT(IN) :: MPTS,NPTS 
      REAL*8, INTENT(OUT) :: DMAT(MPTS,*), DIST(NPTS)
!.... local variables
      REAL*8 DIST2
!
!----------------------------------------------------------------------!
!
!.... calculate the SPD distance matrix
      DO 1 IPTS=1,NPTS
         DO 2 JPTS=IPTS,NPTS
            DIST2 = (XLOCS(IPTS) - XLOCS(JPTS))**2
     ;            + (YLOCS(IPTS) - YLOCS(JPTS))**2 
            DMAT(IPTS,JPTS) = DSQRT(DIST2) 
            IF (IPTS.NE.JPTS) DMAT(JPTS,IPTS) = DMAT(IPTS,JPTS)
    2    CONTINUE
         DIST(IPTS) = DSQRT( (X0 - XLOCS(IPTS))**2 
     ;                     + (Y0 - YLOCS(IPTS))**2 )
    1 CONTINUE  
      RETURN
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE GVMAT_KRIGE(LDV, NPTS,JOB,DECSCL,DECFAC,X0,Y0,
     ;                       XLOCS,YLOCS, V,VMAT)
!
!     Generates the V matrix for a normal or simple krige (default)
!
!     INPUT      MEANING
!     -----      ------- 
!     DECFAC     expontial decay factor; > 0.0 in neighbors importance
!     DECSCL     scale factor on exponetial decay; > 0.0
!     DMAT       matrix of distances to points
!     JOB        = 1 ordinary krige, otherwise simple krige (default) 
!     LDD        leading dimension of DMAT
!     LDV        leading dimension for VMAT
!     NPTS       number of observations
!     X0         x location of point to krige
!     XLOCS      x locations of observation points
!     Y0         y location of point to krige
!     YLOCS      y locations of observation points
!     
!     OUTPUT     MEANING
!     ------     ------- 
!     VMAT       covariance matrix for an ordinary or simple krige
!     V          RHS for krige
!
!.... variable declarations
      REAL*8, INTENT(IN) :: XLOCS(NPTS), YLOCS(NPTS), X0,Y0, 
     ;                      DECSCL,DECFAC
      INTEGER*4, INTENT(IN) :: LDV,NPTS,JOB
      REAL*8, INTENT(OUT) :: VMAT(LDV,*), V(*) 
!.... local variables
      REAL*8 ARG, DIST
!
!---------------------------------------------------------------------!
!
!.... calculate covariance off diagonals
      IROW = 1
      DO 1 IPTS=1,NPTS
         DIST = DSQRT( (XLOCS(IPTS) - X0)**2 + (YLOCS(IPTS) - Y0)**2 )
         ARG =-DIST/DECFAC 
         V(IPTS) = DECSCL*DEXP(ARG)
         VMAT(IROW,IROW) = DECSCL
         IROW = IROW + 1 
         JCOL = IROW
         KPTS = 0 
         DO 2 JPTS=IPTS,NPTS
            JCOL = JCOL + 1
            DIST = DSQRT( (XLOCS(IPTS) - XLOCS(JPTS))**2 
     ;                  + (YLOCS(IPTS) - YLOCS(JPTS))**2 )
            ARG =-DIST/DECFAC
            VMAT(IPTS,JPTS) = DECSCL*DEXP(ARG)
            VMAT(JPTS,IPTS) = VMAT(IPTS,JPTS)
    2    CONTINUE
    1 CONTINUE 
!
!.... constraints for ordinary krige 
      IF (JOB.EQ.1) THEN
         V(NPTS+1) = 1.D0 
         VMAT(NPTS+1,NPTS+1) = 0.D0
         DO 3 IPTS=1,NPTS
            VMAT(NPTS+1,IPTS) = 1.D0
            VMAT(IPTS,NPTS+1) = 1.D0
    3    CONTINUE
      ENDIF
      RETURN
      END
