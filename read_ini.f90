!     implicit real*8 (a-h,o-z)
!     include 'fwd_struc.h'
!     type (mesh_info) msh
!     type (inv_info) inv
!     logical*4 lsame_igcaps, lfiles, LPRIOR,LJSRCH,LINREC,LKBDRY
!     character(80) projnm
!     projnm(1:80) = ' '
!     projnm = 'tswhl'
!     projnm = adjustl(projnm)
!     call read_fwd_ini(projnm,.true., lfiles, msh, ierr)
!     if (ierr /= 0) then
!        write(*,*) 'error calling read_fwd_ini!'
!        stop
!     endif 
!     call read_inv_ini(projnm,.true., &
!                       LPRIOR,LJSRCH,LINREC,LKBDRY, &
!                       MAXFEV,KERNEL,                                 &   
!                       GTOL8,XTOL8,FTOL8, STPMIN,STPMAX, DS0,         &   
!                       SIGMA_BDYX,SIGMA_BDYZ, SIGMA_SRFX,SIGMA_SRFZ,  &
!                       VPVAR_BDY,VPVAR_SRF, INV,ierr) 
!     stop
!     end

      SUBROUTINE READ_INV_INI(PROJNM,LVERB,  &
                              LPRIOR,LJSRCH,LINREC,LKBDRY, &
                              ISTOP_PT,MAXFEV,KERNEL,                        &   
                              GTOL8,XTOL8,FTOL8, STPMIN,STPMAX, DS0,DP0,     &   
                              SIGMA_BDYX,SIGMA_BDYZ, SIGMA_SRFX,SIGMA_SRFZ,  &
                              VPVAR_BDY,VPVAR_SRF, VPMAX_INV,VPMIN_INV,      &
                              INV,IERR) 
!     Reads the initialization file for the inversion 
      implicit none
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV
      CHARACTER(*), INTENT(IN) :: PROJNM
      LOGICAL*4, INTENT(IN) :: LVERB
      REAL*8, INTENT(OUT) :: GTOL8, XTOL8, FTOL8, STPMIN, STPMAX, DS0,DP0, &
                             SIGMA_BDYX, SIGMA_BDYZ, SIGMA_SRFX, SIGMA_SRFZ, &
                             VPVAR_BDY, VPVAR_SRF, VPMAX_INV, VPMIN_INV 
      INTEGER*4, INTENT(OUT) :: MAXFEV, KERNEL, ISTOP_PT, IERR
      LOGICAL*4, INTENT(OUT) :: LPRIOR,LJSRCH,LINREC,LKBDRY
!.... local variables
      CHARACTER(80) SPECFL
      LOGICAL*4 LEX
      INTEGER*4, PARAMETER :: IUNIT = 65
!
!----------------------------------------------------------------------------------------!
!
      IF (LVERB) WRITE(*,*) 'read_inv_ini: Setting defaults...'
      CALL INV_DEFAULTS(LPRIOR,LJSRCH,LINREC,LKBDRY,                   &   
                        MAXFEV,KERNEL,ISTOP_PT,                        &   
                        GTOL8,XTOL8,FTOL8, STPMIN,STPMAX, DS0,DP0,     &   
                        SIGMA_BDYX,SIGMA_BDYZ, SIGMA_SRFX,SIGMA_SRFZ,  &
                        VPVAR_BDY,VPVAR_SRF, VPMAX_INV,VPMIN_INV, INV) 
      IF (LVERB) WRITE(*,*) 'read_inv_ini: Reading spec file...'
!
!.... now try reading
      SPECFL(1:80) = ' ' 
      SPECFL = TRIM(ADJUSTL(PROJNM))//'.spec' 
      SPECFL = ADJUSTL(SPECFL)
      INQUIRE(FILE=TRIM(SPECFL),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'read_ini: Error spec file does not exist!' 
         IERR = 1 
         RETURN
      ENDIF 
!.... open file
      OPEN(UNIT=IUNIT,FILE=TRIM(SPECFL),STATUS='OLD',ACTION='READ') 
      !standard inversion parameters
      CALL READ_VALUE_STRING('CINVTYPE',LVERB, IUNIT,.TRUE., inv%CINVTYPE,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LJSRCH'   ,LVERB, IUNIT,.FALSE., LJSRCH,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LAHESS'   ,LVERB, IUNIT,.FALSE., inv%LAHESS_HOT,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LCOVD_SRF',LVERB, IUNIT,.FALSE., inv%LCOVD_SRF,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LCOVD_BDY',LVERB, IUNIT,.FALSE., inv%LCOVD_BDY,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LDWGHT',   LVERB, IUNIT,.FALSE., inv%LDWGHT,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('IBPHASE', LVERB, IUNIT,.FALSE., inv%IBPHASE,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('IRESTP',  LVERB, IUNIT,.FALSE., inv%IRESTP,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('IMODSRC', LVERB, IUNIT,.FALSE., inv%IMODSRC,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('IMODREC', LVERB, IUNIT,.FALSE., inv%IMODREC,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LINREC',  LVERB, IUNIT,.FALSE., LINREC,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LKBDRY',  LVERB, IUNIT,.FALSE., LKBDRY,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('NCASC',   LVERB, IUNIT,.FALSE., inv%NCASC,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('NORM',    LVERB, IUNIT,.FALSE., inv%NORM, IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('NREC_MIN',LVERB, IUNIT,.FALSE., inv%NREC_MIN,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_INTEGER('MAXIT',   LVERB, IUNIT,.FALSE., inv%MAXIT,IERR)
      IF (IERR /= 0) RETURN
      ! regularization
      CALL READ_VALUE_INTEGER('KERNEL',LVERB, IUNIT,.FALSE., KERNEL,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('SIGMA_BDYX',LVERB, IUNIT,.FALSE., SIGMA_BDYX,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('SIGMA_BDYZ',LVERB, IUNIT,.FALSE., SIGMA_BDYZ,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('SIGMA_SRFX',LVERB, IUNIT,.FALSE., SIGMA_SRFX,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('SIGMA_SRFZ',LVERB, IUNIT,.FALSE., SIGMA_SRFZ,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('VPVAR_BDY',LVERB,  IUNIT,.FALSE., VPVAR_BDY ,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('VPVAR_SRF',LVERB,  IUNIT,.FALSE., VPVAR_SRF, IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('SCLBDY',   LVERB, IUNIT,.FALSE., inv%SCLBDY, IERR)
      ! line search
      CALL READ_VALUE_INTEGER('MAXFEV',LVERB, IUNIT,.FALSE., MAXFEV,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('GTOL8'  ,LVERB, IUNIT,.FALSE., GTOL8,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('XTOL8'  ,LVERB, IUNIT,.FALSE., XTOL8,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('FTOL8'  ,LVERB, IUNIT,.FALSE., FTOL8,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('DS0'    ,LVERB, IUNIT,.FALSE.,   DS0,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('DP0'    ,LVERB, IUNIT,.FALSE.,   DP0,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('STPMIN' ,LVERB, IUNIT,.FALSE.,STPMIN,IERR) 
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('STPMAX' ,LVERB, IUNIT,.FALSE.,STPMAX,IERR)
      IF (IERR /= 0) RETURN
      ! controls in inversion on max and min velocity
      CALL READ_VALUE_DOUBLE('VPMAX_INV' ,LVERB, IUNIT,.FALSE.,VPMAX_INV,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE('VPMIN_INV' ,LVERB, IUNIT,.FALSE.,VPMIN_INV,IERR)
      IF (IERR /= 0) RETURN
      ! maybe an early end to program
      CALL READ_VALUE_INTEGER('STOP_PT',LVERB, IUNIT,.FALSE.,ISTOP_PT,IERR)
      IF (IERR /= 0) RETURN
      ! max residuals
      CALL READ_VALUE_SINGLE('DTMAX_SRF',LVERB, IUNIT,.FALSE.,inv%DTMAX_SRF,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_SINGLE('DTMAX_BDY',LVERB, IUNIT,.FALSE.,inv%DTMAX_BDY,IERR)
      IF (IERR /= 0) RETURN
      ! vp/vs ratio
      CALL READ_VALUE_DOUBLE('VPVS',LVERB, IUNIT,.FALSE.,inv%VPVS,IERR)
!
!.... all done
      CLOSE(IUNIT)
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_FWD_INI(PROJNM,LVERB, TMPDIR,LFILES,LNSRF, MSH,IERR)  
!     Reads the initialization file for the forward simulation
      IMPLICIT NONE
      INCLUDE 'fwd_struc.h' 
      TYPE (MESH_INFO) MSH
      CHARACTER(*), INTENT(IN) :: PROJNM
      LOGICAL*4, INTENT(IN) :: LVERB
      CHARACTER(*), INTENT(OUT) :: TMPDIR
      LOGICAL*4, INTENT(OUT) :: LFILES, LNSRF
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(80) SPECFL
      LOGICAL*4 LEX
      INTEGER*4 INDX 
      INTEGER*4, PARAMETER :: IUNIT = 65
!
!----------------------------------------------------------------------------------------!
!
!.... get the forward problem defaults
      IF (LVERB) WRITE(*,*) 'read_fwd_ini: Setting defaults...'
      CALL FWD_DEFAULTS(TMPDIR,LFILES,LNSRF,MSH)
      IF (LVERB) WRITE(*,*) 'read_fwd_ini: Reading spec file...'
!
!.... now try reading
      SPECFL(1:80) = ' ' 
      SPECFL = TRIM(ADJUSTL(PROJNM))//'.spec' 
      SPECFL = ADJUSTL(SPECFL)
      INQUIRE(FILE=TRIM(SPECFL),EXIST=LEX)
      IF (.NOT.LEX) THEN
         WRITE(*,*) 'read_ini: Error spec file does not exist!' 
         IERR = 1
         RETURN
      ENDIF 
!.... open file
      OPEN(UNIT=IUNIT,FILE=TRIM(SPECFL),STATUS='OLD',ACTION='READ') 
!.... read forward problem parameters; none are necessary
      CALL READ_VALUE_LOGICAL('LFILES',LVERB, IUNIT,.FALSE.,     LFILES,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_LOGICAL('LNSRF', LVERB, IUNIT,.FALSE.,      LNSRF,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('FREQ0', LVERB, IUNIT,.FALSE.,  msh%FREQ0,IERR)  
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('AZTOL', LVERB, IUNIT,.FALSE.,  msh%AZTOL,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('AOITOL',LVERB, IUNIT,.FALSE., msh%AOITOL,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('FCENT_SRF', LVERB, IUNIT,.FALSE., msh%FCENT_SRF, IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('FCENT_BDY', LVERB, IUNIT,.FALSE., msh%FCENT_BDY, IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('KAPPA_SRF', LVERB, IUNIT,.FALSE., msh%KAPPA_SRF, IERR) 
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('KAPPA_BDY', LVERB, IUNIT,.FALSE., msh%KAPPA_BDY, IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('RCOEFF_SRF',LVERB, IUNIT,.FALSE., msh%RCOEFF_SRF,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('RCOEFF_BDY',LVERB, IUNIT,.FALSE., msh%RCOEFF_BDY,IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('DMAX_SRF',  LVERB, IUNIT,.FALSE., msh%DMAX_SRF,  IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_DOUBLE ('DMAX_BDY',  LVERB, IUNIT,.FALSE., msh%DMAX_BDY,  IERR)
      IF (IERR /= 0) RETURN
      CALL READ_VALUE_STRING ('TMPDIR',    LVERB, IUNIT,.FALSE., TMPDIR, IERR)
      INDX = LEN_TRIM(TMPDIR)
      IF (TMPDIR(INDX:INDX) /= '/') TMPDIR(INDX+1:INDX+1) = '/'
      CLOSE(IUNIT) !finished; close file
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !

      SUBROUTINE FWD_DEFAULTS(TMPDIR,LFILES,LNSRF,MSH) 
      INCLUDE 'fwd_struc.h'
      TYPE (MESH_INFO) MSH
      CHARACTER(*), INTENT(OUT) :: TMPDIR
      LOGICAL*4, INTENT(OUT) :: LFILES, LNSRF
      LFILES = .TRUE.  !True; limit output
      LNSRF  = .FALSE. !True; normalize surface wave Greens fns to unity
      msh%FREQ0 = 20.D0    !reference frequency for dispersion
      msh%AZTOL = 10.D0    !azimuth tolerance for py groups
      msh%AOITOL = 5.D0    !angle of incidence tolerance for py groups
      !PML stuff
      msh%FCENT_SRF = 0.D0   !center frequency (hz) 
      msh%FCENT_BDY = 0.D0   !center frequency (hz)
      msh%KAPPA_SRF = 1.D0   !seems to work with Qp and Qs
      msh%KAPPA_BDY = 1.D0   !seems to work with Qp and Qs
      msh%RCOEFF_SRF = 1.D-4 !reflection coefficient
      msh%RCOEFF_BDY = 1.D-4 !reflection coefficient
      msh%DMAX_SRF = 0.D0    !override on dmax function
      msh%DMAX_BDY = 0.D0    !override on dmax function
      !Greens fns
      TMPDIR(1:80) = ' ' 
      TMPDIR = './'          !workspace directory for surface wave calculations
      TMPDIR = ADJUSTL(TMPDIR) 
      RETURN
      END

      SUBROUTINE INV_DEFAULTS(LPRIOR,LJSRCH,LINREC,LKBDRY,                   &
                              MAXFEV,KERNEL,ISTOP_PT,                        &
                              GTOL8,XTOL8,FTOL8, STPMIN,STPMAX, DS0,DP0,     &
                              SIGMA_BDYX,SIGMA_BDYZ, SIGMA_SRFX,SIGMA_SRFZ,  &
                              VPVAR_BDY,VPVAR_SRF, VPMAX_INV,VPMIN_INV, INV) 
      IMPLICIT NONE                        
      INCLUDE 'fwd_struc.h'
      TYPE (INV_INFO) INV
      REAL*8, INTENT(OUT) :: GTOL8, XTOL8, FTOL8, STPMIN, STPMAX, DS0,DP0, &
                             SIGMA_BDYX, SIGMA_BDYZ, SIGMA_SRFX, SIGMA_SRFZ, &
                             VPVAR_BDY, VPVAR_SRF, VPMAX_INV, VPMIN_INV 
      INTEGER*4, INTENT(OUT) :: MAXFEV, KERNEL, ISTOP_PT
      LOGICAL*4, INTENT(OUT) :: LPRIOR,LINREC,LJSRCH,LKBDRY

!.... very important
      inv%CINVTYPE = 'PP'
!.... logical variables
      LPRIOR = .FALSE.  !prior model
      LJSRCH = .FALSE.  !dump output for a jacobian grid search
      inv%LUNWRAP = .FALSE. !don't try to work with unwrapped phase
      inv%LAHESS_HOT = .FALSE. !don't try estimating higher order terms in Hessian
      inv%LCOVD_SRF = .FALSE.  !True -> use surface wave FD data covariance matrix 
      inv%LCOVD_BDY = .FALSE.  !True -> use body wave FD data covariance matrix
      inv%LDWGHT    = .TRUE.   !True -> use normalized distances in data covariance matrix
!.... inversion control 
      inv%IBPHASE = 0   !geometric spreading correction in residual backpropagation
      inv%IRESTP  = 3   !1 -> phase only; 2 -> amplitude only; 3 -> phase/amplitude
      inv%IMODSRC = 3   !0 -> no STF update; 1 -> phase only; 2 -> amplitude only; 3 -> phase/amplitude
      inv%IMODREC = 3   !0 -> no RRF update; 1 -> phase only; 2 -> amplitude only; 3 -> phase/amplitude
      LINREC  = .TRUE.  !include receiver locations in inversion
      LKBDRY  = .TRUE.  !include bielak/interior boundry in inversion 
      inv%NCASC   = 0       !can smooth gradient by nearest neighbor averaging; 
                            !this is the amount of times to cascade that procedure
      inv%NORM    = 2       !norm to calculate residual; (1) L1, (2) L2 
      inv%MAXIT   = 3       !max number of iterations
      inv%NREC_MIN = 5      !minimum number of receivers that must be active
                            !when creating the data covariance matrix
!.... line search
      MAXFEV = 4        !max number of function evaluations in line search
                        !if this fails we take the minimum obj. fn evaluation
      GTOL8 = 9.D-1     !gradient tolerance < 1; smaller -> more exact  
      XTOL8 = EPSILON(1.D0)*100.D0  !x tolerance in line search  
      FTOL8 = 1.D-8                 !necesary decrease; 1.d-4 is a good place to start
      STPMIN = 1.D-20               !min step
      STPMAX = 1.D+20               !max step
      DS0 = 150.D0                  !force model to change this much  
      DP0 = 5.D0                    !force model to change by a percentage
!.... regularization  
      KERNEL = 1               !(1) Gaussian kernel; (2) exponential kernel
      SIGMA_BDYX = 3.D0        !body wave correlation length in x
      SIGMA_BDYZ = 3.D0        !body wave correlation length in z
      SIGMA_SRFX = 10.D3       !surface wave correlation length in x 
      SIGMA_SRFZ = 10.D3       !surface wave correlation length in z
      VPVAR_BDY  = 124500.D0   !scales body wave covariance matrix
      VPVAR_SRF  = 7.07D0      !scales surface wave covariance matrix
      inv%SCLBDY = 1.D0        !grad = grad_srf + sclbdy*grad_bdy
!.... early termination for testing
      ISTOP_PT = 0             !default; let the program run
                               !1 -> end after initialization
                               !2 -> end after initial search/gradient calculation
!.... cutoffs on max and min Vp velocities
      VPMAX_INV = 0.D0         !max velocity to clip at in inversion; <= 0 -> toggle off
      VPMIN_INV = 0.D0         !min velocity to clip at in inversion; <= 0 -> toggle off
!.... max residuals
      inv%DTMAX_SRF = 0.0      !max body wave residual
      inv%DTMAX_BDY = 0.0      !max surface wave residual
!.... vp/vs ratio
      inv%VPVS = DSQRT(3.D0)   !Vp/Vs ratio; hardwire Poisson solid
      RETURN
      END
     
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_VALUE_LOGICAL(VNAME,LVERB, IUNIT,LREQ, LVAR,IERR)
!
!     Reads a logical variable VNAME from file. 
!
!     INPUT      MEANING
!     -----      ------- 
!     IUNIT      unit number of file
!     LREQ       True -> the variable must be located
!     LVERB      True -> will communicate what was read
!     VNAME      variable name in file
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     IVAR       if LREQ is true integer value from file, otherwise
!                value from file if present or default
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: VNAME !variable name i'm looking for
      INTEGER*4, INTENT(IN) :: IUNIT    !unit number of input file
      LOGICAL*4, INTENT(IN) :: LREQ     !True -> this variable is a must have
      LOGICAL*4, INTENT(IN) :: LVERB    !True -> routine says what it does
      INTEGER*4, INTENT(OUT) :: IERR    !error flag
      LOGICAL*4, INTENT(INOUT) :: LVAR  !value of variable
      CHARACTER(256) STRING
      CALL READ_VAR_LINE(VNAME, IUNIT,LREQ,STRING,IERR)
      IF (IERR == 0) THEN 
         READ(STRING,*) LVAR 
         IF (LVERB) &
         WRITE(*,*) 'read_value_logical: ',TRIM(VNAME),' value from file:',LVAR
      ELSE
         IF (IERR == 1) THEN
            WRITE(*,*) 'read_value_logical: Error calling READ_VAR_LINE!'
         ELSE
            IF (LVERB) &
            WRITE(*,*) 'read_value_logical: ',TRIM(VNAME),' set to default:',LVAR
            IERR = 0
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_VALUE_INTEGER(VNAME,LVERB, IUNIT,LREQ, IVAR,IERR)
!
!     Reads an integer value VNAME from file.
!
!     INPUT      MEANING
!     -----      ------- 
!     IUNIT      unit number of file
!     LREQ       True -> the variable must be located
!     LVERB      True -> echo variable value
!     VNAME      variable name in file
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     IVAR       if LREQ is true integer value from file, otherwise
!                value from file if present or default
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: VNAME
      INTEGER*4, INTENT(IN) :: IUNIT 
      LOGICAL*4, INTENT(IN) :: LREQ, LVERB
      INTEGER*4, INTENT(OUT) :: IERR
      INTEGER*4, INTENT(INOUT) :: IVAR
!.... local varaibles
      CHARACTER(256) STRING
!
!----------------------------------------------------------------------------------------!
!
      CALL READ_VAR_LINE(VNAME, IUNIT,LREQ, STRING,IERR)
      IF (IERR == 0) THEN
         READ(STRING,*) IVAR
         IF (LVERB) &
         WRITE(*,*) 'read_value_integer: ',TRIM(VNAME),' value from file:',IVAR
      ELSE
         IF (IERR == 1) THEN
            WRITE(*,*) 'read_value_integer: Error in READ_VAR_LINE!'
         ELSE
            IF (LVERB) &
            WRITE(*,*) 'read_value_integer: ',TRIM(VNAME),' set to default:',IVAR
            IERR = 0
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_VALUE_SINGLE(VNAME,LVERB, IUNIT,LREQ, VAR,IERR)
!     
!     Reads a double precision value VNAME from file.
!
!     INPUT      MEANING
!     -----      ------- 
!     IUNIT      unit number of file
!     LREQ       True -> the variable must be located
!     LVERB      True -> echo variable value
!     VNAME      variable name in file
!        
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     VAR        if LREQ is true integer value from file, otherwise
!                value from file if present or default
!     
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: VNAME
      INTEGER*4, INTENT(IN) :: IUNIT
      LOGICAL*4, INTENT(IN) :: LREQ, LVERB
      INTEGER*4, INTENT(OUT) :: IERR
      REAL*4, INTENT(INOUT) :: VAR
!.... local varaibles
      CHARACTER(256) STRING
!
!----------------------------------------------------------------------------------------!
!
      CALL READ_VAR_LINE(VNAME, IUNIT,LREQ, STRING,IERR)
      IF (IERR == 0) THEN
         READ(STRING,*) VAR
         IF (LVERB) &
         WRITE(*,*) 'read_value_single: ',TRIM(VNAME),' value from file:',VAR
      ELSE
         IF (IERR == 1) THEN
            WRITE(*,*) 'read_value_single: Error in READ_VAR_LINE!'
         ELSE
            IF (LVERB) &
            WRITE(*,*) 'read_value_single: ',TRIM(VNAME),' set to default:',VAR
            IERR = 0
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_VALUE_DOUBLE(VNAME,LVERB, IUNIT,LREQ, VAR,IERR)
!
!     Reads a double precision value VNAME from file.
!
!     INPUT      MEANING
!     -----      ------- 
!     IUNIT      unit number of file
!     LREQ       True -> the variable must be located
!     LVERB      True -> echo variable value
!     VNAME      variable name in file
!
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     VAR        if LREQ is true integer value from file, otherwise
!                value from file if present or default
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: VNAME
      INTEGER*4, INTENT(IN) :: IUNIT
      LOGICAL*4, INTENT(IN) :: LREQ, LVERB
      INTEGER*4, INTENT(OUT) :: IERR
      REAL*8, INTENT(INOUT) :: VAR
!.... local varaibles
      CHARACTER(256) STRING
!
!----------------------------------------------------------------------------------------!
!
      CALL READ_VAR_LINE(VNAME, IUNIT,LREQ, STRING,IERR)
      IF (IERR == 0) THEN
         READ(STRING,*) VAR
         IF (LVERB) &
         WRITE(*,*) 'read_value_double: ',TRIM(VNAME),' value from file:',VAR
      ELSE
         IF (IERR == 1) THEN
            WRITE(*,*) 'read_value_double: Error in READ_VAR_LINE!'
         ELSE
            IF (LVERB) &
            WRITE(*,*) 'read_value_double: ',TRIM(VNAME),' set to default:',VAR
            IERR = 0
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE READ_VALUE_STRING(VNAME, LVERB,IUNIT,LREQ, VAR,IERR)
!     
!     Reads a double precision value VNAME from file.
!
!     INPUT      MEANING
!     -----      ------- 
!     IUNIT      unit number of file
!     LREQ       True -> the variable must be located
!     LVERB      True -> echo variable value
!     VNAME      variable name in file
!        
!     OUTPUT     MEANING
!     ------     -------
!     IERR       error flag
!     VAR        if LREQ is true integer value from file, otherwise
!                value from file if present or default
!     
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: VNAME
      INTEGER*4, INTENT(IN) :: IUNIT
      LOGICAL*4, INTENT(IN) :: LREQ, LVERB
      INTEGER*4, INTENT(OUT) :: IERR
      CHARACTER(*), INTENT(INOUT) :: VAR
!.... local varaibles
      CHARACTER(256) STRING
!
!----------------------------------------------------------------------------------------!
!
      CALL READ_VAR_LINE(VNAME, IUNIT,LREQ, STRING,IERR)
      IF (IERR == 0) THEN
         READ(STRING,'(A)') VAR
         IF (LVERB) &
         WRITE(*,*) 'read_value_string: ',TRIM(VNAME),' value from file: ',TRIM(VAR)
      ELSE
         IF (IERR == 1) THEN
            WRITE(*,*) 'read_value_string: Error in READ_VAR_LINE!'
         ELSE
            IF (LVERB) &
            WRITE(*,*) 'read_value_string: ',TRIM(VNAME),' set to default: ',TRIM(VAR)
            IERR = 0
         ENDIF
      ENDIF
      RETURN
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
    
      SUBROUTINE READ_VAR_LINE(VNAME, IUNIT,LREQ, STRING,IERR)
!
!     Searches a file for a variable.  The variable name must be followed by an equality
!     sign.  Comments are denoted with a # sign
!
!     INPUT      MEANING
!     -----      ------- 
!     IUNIT      unit number of file to search
!     LREQ       True -> variable is required; if not found will ierr = 1
!     VNAME      variable name to search file for
!
!     OUTPUT     MEANING
!     ------     ------- 
!     IERR       (0) no error; -1 variable does not exist; +1 an error occurred
!     STRING     character holding variable if ierr = 0
!
!.... variable declarations
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: VNAME
      INTEGER*4, INTENT(IN) :: IUNIT !unit number of file 
      LOGICAL*4, INTENT(IN) :: LREQ
      CHARACTER(256), INTENT(INOUT) :: STRING
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(256) CWORK
      INTEGER*4 INDX
      LOGICAL*4 LSAME_IGCAPS
!
!----------------------------------------------------------------------------------------!
!
!.... rewind back to the beginning
      IERR = 0 
      STRING(1:256) = ' '
      CWORK(1:256)  = ' '
      REWIND(IUNIT,IOSTAT=IERR)
      IF (IERR /= 0) THEN
         WRITE(*,*) 'read_var_line: There was an error rewinding the file!'
         RETURN
      ENDIF
!
!.... read the file
    1 CONTINUE ! loop on input until we get a non-blank line 
         READ(IUNIT,FMT="(A256)",IOSTAT=IERR,END=60) STRING 
         IF (IERR.NE.0) THEN 
            WRITE(*,*) 'read_var_line: Warning there was an error reading the spec file!'
            IF (.NOT.LREQ) IERR =-1 !this is okay; variable not required
            RETURN
         ENDIF
         STRING = ADJUSTL(STRING) !suppress leading blank spaces
         !Windows/DOS may leave carriage returns; remove that 
         IF (INDEX(STRING,ACHAR(13)).GT.0) & 
         STRING = STRING(1:INDEX(STRING,ACHAR(13))-1) !remove carriage return 
         IF (LEN_TRIM(STRING) == 0) GOTO 100 !line is blank, continue 
         IF (STRING(1:1) == '#') GOTO 100    !line is a comment 
!
!....... get the variable name
         CWORK = STRING(1:LEN_TRIM(STRING))  !suppress trailing blank spaces 
         CWORK = ADJUSTL(CWORK)              !suppress leading blank spaces
         INDX = INDEX(STRING,'=')            !trim anything after equals sign
         IF (INDX <= 1 .OR. INDX == LEN(STRING)) THEN
            WRITE(*,*) 'read_var_line: No equal sign with variable!'
            IERR = 1
            RETURN
         ENDIF
         CWORK = CWORK(1:INDX-1)          !copy everything up to equals sign
         CWORK = ADJUSTL(CWORK)           !remove blank spaces before variable name
         CWORK = CWORK(1:LEN_TRIM(CWORK)) !remove any blank spaces before '=' sign
         IF (LSAME_IGCAPS(CWORK,VNAME,IERR)) THEN
            IF (IERR /= 0) THEN
               WRITE(*,*) 'read_var_line: Error caling LSAME_IGCAPS!'
               RETURN
            ENDIF
            GOTO 200 !done; have the right line
         ENDIF
  100    CONTINUE
      GOTO 1 
   60 CONTINUE !break ahead; end of file
      IF (LREQ) THEN
         WRITE(*,*) 'read_var_line: I couldnt locate required variable: ', &
         TRIM(ADJUSTL(VNAME))
         IERR = 1
      ELSE
         IERR =-1
      ENDIF
      RETURN
  200 CONTINUE !this line likely contains a variable 
!
!.... pull out the string containing the variable (after = sign, before # sign)
      STRING = STRING(1:LEN_TRIM(STRING)) !suppress trailing blank spaces 
      IF (INDEX(STRING,'#') > 0)  &
      STRING = STRING(1:INDEX(STRING,'#')-1) !remove trailing comments
      INDX = INDEX(STRING,'=') !suppress variables leading up to = sign
      IF (INDX <= 1 .OR. INDX == LEN(STRING)) THEN
         WRITE(*,*) 'read_next_line: The format of your input file is incorrect!'
         IERR = 1
         RETURN
      ENDIF
      STRING = STRING(INDX+1:LEN_TRIM(STRING)) !trim anything before '=' sign
      STRING = ADJUSTL(STRING)                 !remove left blank spaces
      STRING = STRING(1:LEN_TRIM(STRING))      !finally remove any trailing blank spaces
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      SUBROUTINE UNUSED_STRING(S)
      CHARACTER(*), INTENT(IN) :: S
      IF (LEN(S) == 1) CONTINUE
      RETURN 
      END
!                                                                                        !
!========================================================================================!
!                                                                                        !
      LOGICAL*4 FUNCTION LSAME_IGCAPS(VNAME,VTARG,IERR) 
!
!     Determines if vname and vtarg are the same whilst ignoring capitlization 
!
!     INPUT     MEANING
!     -----     ------- 
!     VNAME     variable name 
!     VTARG     variable name VNAME is matched to
!
!     OUTPUT    MEANING
!     ------    -------
!     IERR      error raised
!
!.... variable declarations
      CHARACTER(*), INTENT(IN) :: VNAME, VTARG
      INTEGER*4, INTENT(OUT) :: IERR
!.... local variables
      CHARACTER(1), ALLOCATABLE :: VW1(:), VW2(:)
      CHARACTER(1) CU1, CU2
      CHARACTER(37), PARAMETER :: CASEU = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
      CHARACTER(37), PARAMETER :: CASEL = 'abcdefghijklmnopqrstuvwxyz0123456789_'
!
!----------------------------------------------------------------------------------------!
!
!.... get the string length
      IERR = 0
      LSAME_IGCAPS = .FALSE.
      LENS1 = LEN_TRIM(ADJUSTL(VNAME))
      LENS2 = LEN_TRIM(ADJUSTL(VTARG))
      IF (LENS1 /= LENS2) RETURN !different lengths -> different words 
      LENS = LENS1 
!
!.... copy names 
      ALLOCATE(VW1(LENS))
      ALLOCATE(VW2(LENS))
      LCOPY = LEN(VNAME)
      I1 = 0
      DO 1 I=1,LCOPY
         IF (VNAME(I:I) /= ' ') THEN
            I1 = I1 + 1
            VW1(I1) = VNAME(I:I) 
         ENDIF
    1 CONTINUE 
      LCOPY = LEN(VTARG)
      I1 = 0
      DO 2 I=1,LCOPY
         IF (VTARG(I:I) /= ' ') THEN
            I1 = I1 + 1
            VW2(I1) = VTARG(I:I) 
         ENDIF
    2 CONTINUE
!
!.... loop on string length 
      DO 3 ILEN=1,LENS
         DO 4 I=1,37 
            IF (VW1(ILEN) == CASEL(I:I)) THEN
               CU1 = CASEU(I:I)
               GOTO 40
            ENDIF
            IF (VW1(ILEN) == CASEU(I:I)) THEN
               CU1 = CASEU(I:I)
               GOTO 40
            ENDIF
    4    CONTINUE
         WRITE(*,*) 'lsame_caps: I dont understand this symbol 1!'
         IERR = 1
         RETURN
   40    CONTINUE !break ahead; letter found

         DO 5 I=1,37 
            IF (VW2(ILEN) == CASEL(I:I)) THEN
               CU2 = CASEU(I:I)
               GOTO 50
            ENDIF
            IF (VW2(ILEN) == CASEU(I:I)) THEN
               CU2 = CASEU(I:I)
               GOTO 50
            ENDIF
    5    CONTINUE
         WRITE(*,*) 'lsame_caps: I dont understand this symbol 2!'
         IERR = 1 
         RETURN
   50    CONTINUE !break ahead; letter found
         IF (CU1 /= CU2) GOTO 100 
    3 CONTINUE
      LSAME_IGCAPS = .TRUE.
  100 CONTINUE
      DEALLOCATE(VW1)
      DEALLOCATE(VW2)
      RETURN
      END
